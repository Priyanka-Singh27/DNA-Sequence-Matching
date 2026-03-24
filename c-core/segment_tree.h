#ifndef SEGMENT_TREE_H
#define SEGMENT_TREE_H

/*
 * segment_tree.h
 * ==============
 * Segment Tree with coordinate compression for genome SNP position queries.
 *
 * YOUR JOB (segment_tree.c):
 *   Implement every function declared in this header.
 *   Do NOT modify this header file.
 *   Do NOT use any external DS libraries.
 *   Do NOT printf inside library functions.
 *
 * DEPENDENCIES:
 *   #include "dna_utils.h"   (already written — use freely)
 *
 * BUILD & TEST:
 *   gcc -O2 -Wall -std=c11 -DSGT_TEST dna_utils.c segment_tree.c -o test_sgt -lm
 *   ./test_sgt
 *
 * BUILD OBJECT FILE:
 *   gcc -O2 -Wall -std=c11 -fPIC -c segment_tree.c -o segment_tree.o
 *
 * USED BY FEATURES:
 *   Feature 6 — "Find all known mutations between position 1000 and 5000"
 *   Feature 9 — Map matched SNP positions to disease labels and severity
 *
 * KEY CONCEPT — WHY COORDINATE COMPRESSION?
 *   The human genome has ~3,000,000,000 positions.
 *   We cannot allocate an array of 3 billion entries.
 *   BUT ClinVar only gives us ~10,000 known SNP positions.
 *   Coordinate compression maps those 10,000 real positions
 *   to dense indices 0..9999.
 *   All segment tree operations work on compressed indices.
 *   Results are converted back to real positions before returning.
 *
 * EXAMPLE:
 *   Real positions:   [17, 1234, 45678, 1000000]
 *   Compressed index: [ 0,    1,     2,        3]
 *
 *   query(real_left=1, real_right=50000)
 *     → compress to [0, 2]
 *     → segment tree range query on [0, 2]
 *     → found at compressed indices 0 and 1
 *     → return results with real positions 17 and 1234
 */

#include <stddef.h>
#include "dna_utils.h"

/* ============================================================
   CONSTANTS
   ============================================================ */

/* Maximum number of SNP positions the segment tree can hold */
#define SGT_MAX_POSITIONS   100000

/* Return codes */
#define SGT_OK    1
#define SGT_ERR   0

/* ============================================================
   STRUCTS
   ============================================================ */

/*
 * SegmentQueryResult
 * ------------------
 * One entry in a range query result.
 * disease is a heap-allocated string — caller must free it.
 */
typedef struct {
    int   position;   /* REAL genome position (not compressed) */
    char* disease;    /* heap-allocated disease name            */
    int   severity;   /* 0=low, 1=medium, 2=high               */
} SegmentQueryResult;

/*
 * SegmentTree
 * -----------
 * The complete segment tree structure.
 *
 * INTERNAL ARRAYS:
 *
 *   tree[]
 *     The segment tree itself. Array-based, 1-indexed.
 *     Size = 4 * n (standard allocation for segment tree on n elements).
 *     tree[node] = maximum severity value in the range covered by node.
 *     (We store max severity so we can quickly skip ranges with no mutations.)
 *
 *   lazy[]
 *     Lazy propagation array. Same size as tree[].
 *     Used for bulk updates when loading the full ClinVar dataset.
 *
 *   real_positions[]
 *     The SORTED array of real genome positions.
 *     real_positions[i] = actual genome coordinate of compressed index i.
 *     Size = n.
 *
 *   diseases[]
 *     diseases[i] = heap-allocated disease name at compressed index i.
 *     Size = n.
 *
 *   severities[]
 *     severities[i] = severity at compressed index i.
 *     Size = n.
 *
 *   n = number of SNP positions stored.
 */
typedef struct {
    int*   tree;             /* segment tree array, size = 4 * n  */
    int*   lazy;             /* lazy propagation, size = 4 * n    */
    int*   real_positions;   /* sorted real genome positions       */
    char** diseases;         /* disease name at each index         */
    int*   severities;       /* severity at each index             */
    int    n;                /* number of SNP positions            */
} SegmentTree;

/* ============================================================
   CORE API — YOU MUST IMPLEMENT ALL OF THESE
   ============================================================ */

/* ------------------------------------------------------------
 * segment_tree_build
 * ------------------
 * Build a segment tree from arrays of SNP positions and metadata.
 *
 * INPUT:
 *   positions  — array of REAL genome positions (unsorted is OK)
 *   diseases   — disease name at each position (parallel array)
 *   severities — severity at each position (parallel array)
 *   count      — number of entries (length of all three arrays)
 *
 * OUTPUT:
 *   Heap-allocated SegmentTree* on success.
 *   NULL on failure (NULL input, count == 0, alloc failure).
 *
 * ALGORITHM:
 *   Step 1 — Sort positions[] ascending, keeping diseases[] and
 *             severities[] aligned (sort all three arrays together
 *             by positions[] key).
 *
 *   Step 2 — Remove duplicates: if two entries have the same position,
 *             keep the one with higher severity.
 *
 *   Step 3 — Store sorted positions in real_positions[].
 *             Copy disease strings (heap-allocate each one).
 *             Copy severities[].
 *             Set tree->n = deduplicated count.
 *
 *   Step 4 — Allocate tree[] and lazy[] of size 4 * n.
 *             Initialise all to 0.
 *
 *   Step 5 — For each index i in [0, n-1]:
 *             Call the internal point-update on compressed index i
 *             with value = severities[i].
 *             (This fills the tree[] array correctly.)
 *
 * TIME:  O(n log n) for sorting + O(n log n) for building
 * ------------------------------------------------------------ */
SegmentTree* segment_tree_build(int*   positions,
                                 char** diseases,
                                 int*   severities,
                                 int    count);


/* ------------------------------------------------------------
 * segment_tree_update
 * -------------------
 * Add or update a single SNP position in the tree.
 *
 * INPUT:
 *   st           — the segment tree
 *   real_pos     — real genome position to update
 *   disease      — disease name (heap-copied internally)
 *   severity     — new severity value
 *
 * ALGORITHM:
 *   Step 1 — Find compressed index: idx = compress(st, real_pos)
 *             If real_pos is not in real_positions[]:
 *               This is a new position we haven't seen.
 *               For simplicity, only update existing positions.
 *               (Full dynamic insert would require rebuilding.)
 *               Return SGT_ERR.
 *
 *   Step 2 — Update diseases[idx] (free old string, strdup new one).
 *             Update severities[idx] = severity.
 *
 *   Step 3 — Point update on tree[]:
 *             Standard segment tree point update at compressed index idx.
 *             Update all ancestors by propagating max value upward.
 *
 * Returns SGT_OK on success, SGT_ERR if position not found.
 * ------------------------------------------------------------ */
int segment_tree_update(SegmentTree* st,
                         int          real_pos,
                         const char*  disease,
                         int          severity);


/* ------------------------------------------------------------
 * segment_tree_query
 * ------------------
 * Find all SNPs with known mutations in the real position range
 * [left_real, right_real] (inclusive).
 *
 * INPUT:
 *   st         — the segment tree
 *   left_real  — real genome position, left bound (inclusive)
 *   right_real — real genome position, right bound (inclusive)
 *   count      — OUTPUT: set to number of results found
 *
 * OUTPUT:
 *   Heap-allocated SegmentQueryResult[] of length *count.
 *   Each entry has a heap-allocated disease string.
 *   Returns NULL and *count=0 if nothing found or on error.
 *   CALLER must call segment_tree_results_free() on the result.
 *
 * ALGORITHM:
 *   Step 1 — Coordinate compression:
 *     cl = lower_bound(real_positions, n, left_real)
 *          (first index where real_positions[cl] >= left_real)
 *     cr = upper_bound(real_positions, n, right_real) - 1
 *          (last index where real_positions[cr] <= right_real)
 *     If cl > cr → no positions in range → return NULL, *count=0
 *
 *   Step 2 — Standard segment tree range query on [cl, cr]:
 *     Walk the tree, collect all LEAF nodes in range [cl, cr]
 *     where severity > 0 (i.e. has a known mutation).
 *     Do NOT just check the aggregate — collect individual leaves.
 *
 *   Step 3 — For each found compressed index i:
 *     Build SegmentQueryResult:
 *       position = real_positions[i]
 *       disease  = strdup(diseases[i])
 *       severity = severities[i]
 *
 *   Step 4 — Sort results by position ascending before returning.
 *
 * TIME:  O(log n + k) where k = number of results
 * ------------------------------------------------------------ */
SegmentQueryResult* segment_tree_query(SegmentTree* st,
                                        int          left_real,
                                        int          right_real,
                                        int*         count);


/* ------------------------------------------------------------
 * segment_tree_point_query
 * ------------------------
 * Find the SNP at exactly one genome position.
 *
 * INPUT:
 *   st       — the segment tree
 *   real_pos — exact real genome position to look up
 *
 * OUTPUT:
 *   Heap-allocated SegmentQueryResult* with one entry, or NULL if
 *   no SNP exists at that position or on error.
 *   CALLER must call segment_tree_results_free() with count=1.
 *
 * IMPLEMENTATION:
 *   Calls segment_tree_query(st, real_pos, real_pos, &count).
 *   Returns results[0] if count==1, else NULL.
 * ------------------------------------------------------------ */
SegmentQueryResult* segment_tree_point_query(SegmentTree* st,
                                              int          real_pos);


/* ------------------------------------------------------------
 * segment_tree_load_clinvar
 * -------------------------
 * Build a SegmentTree by loading all SNP positions from clinvar_snps.tsv.
 *
 * FILE FORMAT (tab-separated, no header):
 *   snp_id \t sequence \t disease \t severity \t genome_position
 *
 * WAIT — our clinvar_snps.tsv has 4 columns, not 5.
 * For the segment tree we need genome positions.
 * Use this EXTENDED format with 5 columns:
 *   rs334 \t GTGCACCTGACTCCTGTG \t Sickle Cell Anemia \t 2 \t 5248232
 *
 * ALGORITHM:
 *   Parse all lines → collect positions[], diseases[], severities[]
 *   Call segment_tree_build(positions, diseases, severities, count)
 *   Set *st_out = result
 *   Return count of SNPs loaded, or -1 on file error.
 *
 * NOTE: The trie_errors module reads the same TSV but uses the
 * sequence column. The segment tree uses the position column.
 * Both can read from the same file — they just use different columns.
 * ------------------------------------------------------------ */
int segment_tree_load_clinvar(SegmentTree** st_out,
                               const char*   filepath);


/* ------------------------------------------------------------
 * segment_tree_results_free
 * -------------------------
 * Free a SegmentQueryResult array.
 *
 * For each entry: free(results[i].disease)
 * Then: free(results)
 *
 * Safe to call with NULL results or count=0.
 * ------------------------------------------------------------ */
void segment_tree_results_free(SegmentQueryResult* results, int count);


/* ------------------------------------------------------------
 * segment_tree_free
 * -----------------
 * Free all memory used by the segment tree.
 *
 * Steps:
 *   1. For each i in [0, n-1]: free(st->diseases[i])
 *   2. free(st->diseases)
 *   3. free(st->real_positions)
 *   4. free(st->severities)
 *   5. free(st->tree)
 *   6. free(st->lazy)
 *   7. free(st)
 *
 * Safe to call with NULL.
 * ------------------------------------------------------------ */
void segment_tree_free(SegmentTree* st);


/* ============================================================
   JSON SERIALISATION — implement for Flask integration
   ============================================================ */

/* ------------------------------------------------------------
 * segment_tree_results_to_json
 * ----------------------------
 * Serialise a SegmentQueryResult array to heap-allocated JSON.
 *
 * Output format:
 * {
 *   "count": 2,
 *   "mutations": [
 *     { "position": 17,   "disease": "Sickle Cell Anemia", "severity": 2 },
 *     { "position": 1234, "disease": "Beta-Thalassemia",   "severity": 1 }
 *   ]
 * }
 *
 * Returns NULL on allocation failure or NULL input.
 * CALLER must free() the returned string.
 * ------------------------------------------------------------ */
char* segment_tree_results_to_json(SegmentQueryResult* results, int count);


/* ============================================================
   INTERNAL HELPERS — declare here so main_test.c can test them
   ============================================================ */

/* ------------------------------------------------------------
 * sgt_lower_bound
 * ---------------
 * Binary search: returns first index i where arr[i] >= val.
 * Returns n if all elements < val.
 * arr must be sorted ascending.
 * ------------------------------------------------------------ */
int sgt_lower_bound(int* arr, int n, int val);


/* ------------------------------------------------------------
 * sgt_upper_bound
 * ---------------
 * Binary search: returns first index i where arr[i] > val.
 * Equivalently, last index where arr[i] <= val is (result - 1).
 * Returns 0 if all elements > val.
 * arr must be sorted ascending.
 * ------------------------------------------------------------ */
int sgt_upper_bound(int* arr, int n, int val);


/* ============================================================
   OPTIONAL DEBUG HELPERS
   ============================================================ */
#ifdef SGT_DEBUG
/* Print the tree array (for debugging) */
void segment_tree_print(SegmentTree* st);

/* Print all stored positions and their diseases */
void segment_tree_print_positions(SegmentTree* st);
#endif

#endif /* SEGMENT_TREE_H */