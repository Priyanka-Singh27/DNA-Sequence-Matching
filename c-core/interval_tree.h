#ifndef INTERVAL_TREE_H
#define INTERVAL_TREE_H

/*
 * interval_tree.h
 * ===============
 * Augmented BST (Interval Tree) for gene region overlap queries.
 * Each node stores a genomic interval [low, high] representing a gene.
 * The tree answers: "which genes overlap this position or range?"
 *
 * YOUR JOB (interval_tree.c):
 *   Implement every function declared in this header.
 *   Do NOT modify this header file.
 *   Do NOT use any external DS libraries.
 *   Do NOT printf inside library functions.
 *
 * DEPENDENCIES:
 *   #include "dna_utils.h"   (already written — use freely)
 *
 * BUILD & TEST:
 *   gcc -O2 -Wall -std=c11 -DIT_TEST dna_utils.c interval_tree.c -o test_it -lm
 *   ./test_it
 *
 * BUILD OBJECT FILE:
 *   gcc -O2 -Wall -std=c11 -fPIC -c interval_tree.c -o interval_tree.o
 *
 * USED BY FEATURES:
 *   Feature 6 — Mutation detection: "which gene does position X fall in?"
 *
 * KEY CONCEPT — THE max FIELD (AUGMENTATION):
 *
 *   A normal BST would require O(n) per overlap query (visit all nodes).
 *   Augmenting each node with the maximum 'high' value in its subtree
 *   allows us to PRUNE entire subtrees during search:
 *
 *   If node->left->max < query_low:
 *     NO interval in the left subtree can overlap [query_low, query_high]
 *     → Skip entire left subtree
 *
 *   This reduces query time to O(log n + k) where k = results found.
 *
 * EXAMPLE:
 *   Genes:
 *     HBB:   [5246695, 5250625]    chr11
 *     BRCA1: [43044295, 43125483]  chr17
 *     TP53:  [7661779, 7687550]    chr17
 *
 *   Query: point_query(position=5248000)
 *     → finds HBB (5246695 <= 5248000 <= 5250625)
 *     → does NOT visit BRCA1 subtree (max pruning)
 */

#include <stddef.h>
#include "dna_utils.h"

/* ============================================================
   CONSTANTS
   ============================================================ */

/* Maximum number of gene intervals the tree can hold */
#define IT_MAX_GENES     100000

/* Return codes */
#define IT_OK    1
#define IT_ERR   0

/* ============================================================
   STRUCTS
   ============================================================ */

/*
 * IntervalNode
 * ------------
 * One node in the interval tree (BST node).
 *
 * BST PROPERTY: sorted by 'low' value.
 *   node->left  contains intervals with low < node->low
 *   node->right contains intervals with low >= node->low
 *
 * AUGMENTATION: max field.
 *   max = maximum 'high' value among ALL nodes in this subtree.
 *   Updated on every insert, on the path back up to the root.
 *
 *   RULE: node->max = max(node->high, left->max, right->max)
 *         where left->max = 0 if left child is NULL
 *               right->max = 0 if right child is NULL
 *
 * STRINGS:
 *   gene_name and chromosome are heap-allocated.
 *   Must be freed in interval_tree_free().
 */
typedef struct IntervalNode {
    int    low;                    /* gene start position (BST key)      */
    int    high;                   /* gene end position                  */
    int    max;                    /* max high in subtree (AUGMENTATION) */
    char*  gene_name;              /* heap-allocated e.g. "HBB"         */
    char*  chromosome;             /* heap-allocated e.g. "chr11"       */
    struct IntervalNode* left;
    struct IntervalNode* right;
} IntervalNode;

/*
 * IntervalTree
 * ------------
 * The complete interval tree.
 * root = NULL for an empty tree.
 * size = number of gene intervals inserted.
 */
typedef struct {
    IntervalNode* root;
    int           size;
} IntervalTree;

/*
 * IntervalResult
 * --------------
 * One entry in the result of an overlap query.
 * gene_name and chromosome are heap-allocated strings.
 * CALLER must free them and the result array with
 * interval_tree_results_free().
 */
typedef struct {
    int   low;          /* gene start position */
    int   high;         /* gene end position   */
    char* gene_name;    /* heap-allocated      */
    char* chromosome;   /* heap-allocated      */
} IntervalResult;

/* ============================================================
   CORE API — YOU MUST IMPLEMENT ALL OF THESE
   ============================================================ */

/* ------------------------------------------------------------
 * interval_tree_create
 * --------------------
 * Allocate and initialise an empty IntervalTree.
 *
 * Sets root = NULL, size = 0.
 * Returns heap-allocated IntervalTree* on success, NULL on failure.
 * ------------------------------------------------------------ */
IntervalTree* interval_tree_create(void);


/* ------------------------------------------------------------
 * interval_tree_free
 * ------------------
 * Free all memory used by the interval tree.
 *
 * TRAVERSAL ORDER: post-order DFS (children before parent).
 * For each node:
 *   1. Recursively free left subtree.
 *   2. Recursively free right subtree.
 *   3. free(node->gene_name)
 *   4. free(node->chromosome)
 *   5. free(node)
 * After all nodes: free(tree)
 *
 * Safe to call with NULL.
 * ------------------------------------------------------------ */
void interval_tree_free(IntervalTree* it);


/* ------------------------------------------------------------
 * interval_tree_insert
 * --------------------
 * Insert a gene interval into the tree.
 *
 * INPUT:
 *   it         — the interval tree
 *   low        — gene start position (must be <= high)
 *   high       — gene end position
 *   gene_name  — gene name string (will be heap-copied)
 *   chromosome — chromosome name e.g. "chr11" (will be heap-copied)
 *
 * ALGORITHM:
 *   Step 1 — Create new node:
 *     node->low        = low
 *     node->high       = high
 *     node->max        = high     (no children yet, max = own high)
 *     node->gene_name  = strdup(gene_name)
 *     node->chromosome = strdup(chromosome)
 *     node->left = node->right = NULL
 *
 *   Step 2 — BST insert by 'low':
 *     Walk from root:
 *       if new.low < current.low  → go left
 *       if new.low >= current.low → go right
 *     Insert at the NULL position found.
 *
 *   Step 3 — UPDATE max ON THE WAY BACK UP:
 *     This is the AUGMENTATION step.
 *     After inserting, walk back up the path to root.
 *     At each ancestor node, recalculate:
 *       node->max = max(node->high,
 *                       node->left  ? node->left->max  : 0,
 *                       node->right ? node->right->max : 0)
 *
 *     IMPLEMENTATION TIP: Use a parent pointer stack or implement
 *     insert recursively (recursive insert naturally updates max
 *     on the way back up via the return path).
 *
 *   Step 4 — it->size++
 *
 * Does nothing if it or gene_name or chromosome is NULL.
 * Does nothing if low > high.
 * ------------------------------------------------------------ */
void interval_tree_insert(IntervalTree* it,
                           int           low,
                           int           high,
                           const char*   gene_name,
                           const char*   chromosome);


/* ------------------------------------------------------------
 * interval_tree_query
 * -------------------
 * Find all genes whose interval overlaps [query_low, query_high].
 *
 * Two intervals [a,b] and [c,d] OVERLAP if: a <= d AND b >= c
 * In our case:  node->low <= query_high AND node->high >= query_low
 *
 * INPUT:
 *   it          — the interval tree
 *   query_low   — left bound of query range (inclusive)
 *   query_high  — right bound of query range (inclusive)
 *   count       — OUTPUT: set to number of overlapping genes found
 *
 * OUTPUT:
 *   Heap-allocated IntervalResult[] of length *count.
 *   Each entry has heap-allocated gene_name and chromosome strings.
 *   Returns NULL and *count=0 if nothing found or on error.
 *   CALLER must call interval_tree_results_free() on the result.
 *
 * ALGORITHM (recursive DFS with pruning):
 *
 *   visit(node, query_low, query_high, results):
 *
 *   Base case: if node == NULL → return
 *
 *   PRUNING: if node->max < query_low → return
 *     (No interval in this subtree can overlap — all highs < query_low)
 *
 *   Recurse left: visit(node->left, ...)
 *
 *   CHECK THIS NODE:
 *     if node->low <= query_high AND node->high >= query_low:
 *       → OVERLAP FOUND
 *       → append IntervalResult to results array:
 *            low        = node->low
 *            high       = node->high
 *            gene_name  = strdup(node->gene_name)
 *            chromosome = strdup(node->chromosome)
 *
 *   Recurse right (only if node->low <= query_high):
 *     visit(node->right, ...)
 *     (If node->low > query_high, all right subtree nodes have
 *      low > query_high too, so they cannot overlap.)
 *
 * RESULT COLLECTION:
 *   Use a dynamic array (start with capacity 16, realloc as needed).
 *   After DFS, copy to exactly-sized array for return.
 *
 * TIME: O(log n + k) due to max pruning
 * ------------------------------------------------------------ */
IntervalResult* interval_tree_query(IntervalTree* it,
                                     int           query_low,
                                     int           query_high,
                                     int*          count);


/* ------------------------------------------------------------
 * interval_tree_point_query
 * -------------------------
 * Find all genes that contain a single genome position.
 *
 * A gene [low, high] contains position p if: low <= p <= high
 *
 * INPUT:
 *   it       — the interval tree
 *   position — genome position to query
 *   count    — OUTPUT: set to number of genes containing position
 *
 * OUTPUT:
 *   Same as interval_tree_query — IntervalResult[] or NULL.
 *   CALLER must call interval_tree_results_free() on result.
 *
 * IMPLEMENTATION:
 *   Calls interval_tree_query(it, position, position, count).
 * ------------------------------------------------------------ */
IntervalResult* interval_tree_point_query(IntervalTree* it,
                                           int           position,
                                           int*          count);


/* ------------------------------------------------------------
 * interval_tree_load_genes
 * ------------------------
 * Load gene regions from a TSV file into the tree.
 *
 * FILE FORMAT (tab-separated, no header):
 *   gene_name \t chromosome \t start \t end
 *
 * EXAMPLE LINES:
 *   HBB	chr11	5246695	5250625
 *   BRCA1	chr17	43044295	43125483
 *   TP53	chr17	7661779	7687550
 *   OPN1LW	chrX	154144243	154159032
 *   CFTR	chr7	117480025	117668665
 *
 * ALGORITHM:
 *   fopen(filepath, "r")
 *   For each line:
 *     Parse 4 fields by splitting on '\t'
 *     start = atoi(start_field)
 *     end   = atoi(end_field)
 *     interval_tree_insert(it, start, end, gene_name, chromosome)
 *     count++
 *   fclose
 *
 * Returns number of genes successfully loaded.
 * Returns -1 if file cannot be opened.
 * Skips malformed lines (not exactly 4 fields) silently.
 * ------------------------------------------------------------ */
int interval_tree_load_genes(IntervalTree* it, const char* filepath);


/* ------------------------------------------------------------
 * interval_tree_results_free
 * --------------------------
 * Free an IntervalResult array returned by any query function.
 *
 * For each entry:
 *   free(results[i].gene_name)
 *   free(results[i].chromosome)
 * Then: free(results)
 *
 * Safe to call with NULL results or count=0.
 * ------------------------------------------------------------ */
void interval_tree_results_free(IntervalResult* results, int count);


/* ============================================================
   JSON SERIALISATION — implement for Flask integration
   ============================================================ */

/* ------------------------------------------------------------
 * interval_tree_results_to_json
 * -----------------------------
 * Serialise an IntervalResult array to heap-allocated JSON.
 *
 * Output format:
 * {
 *   "count": 2,
 *   "genes": [
 *     {
 *       "gene":       "HBB",
 *       "chromosome": "chr11",
 *       "start":      5246695,
 *       "end":        5250625
 *     },
 *     {
 *       "gene":       "TP53",
 *       "chromosome": "chr17",
 *       "start":      7661779,
 *       "end":        7687550
 *     }
 *   ]
 * }
 *
 * Returns NULL on allocation failure or NULL input.
 * CALLER must free() the returned string.
 * ------------------------------------------------------------ */
char* interval_tree_results_to_json(IntervalResult* results, int count);


/* ============================================================
   OPTIONAL DEBUG HELPERS
   ============================================================ */
#ifdef IT_DEBUG
/* Print the tree structure (in-order traversal) */
void interval_tree_print(IntervalTree* it);

/* Print max values at each node */
void interval_tree_print_max(IntervalTree* it);

/* Verify max field is correct at every node (for testing) */
int interval_tree_verify_max(IntervalTree* it);
#endif

#endif /* INTERVAL_TREE_H */