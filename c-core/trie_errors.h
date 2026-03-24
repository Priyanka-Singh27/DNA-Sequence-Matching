#ifndef TRIE_ERRORS_H
#define TRIE_ERRORS_H

/*
 * trie_errors.h
 * =============
 * Trie (prefix tree) with fuzzy/approximate string matching.
 * Stores known SNP (mutation) marker sequences.
 * Supports exact search AND fuzzy search allowing k mismatches.
 *
 * YOUR JOB (trie_errors.c):
 *   Implement every function declared in this header.
 *   Do NOT modify this header file.
 *   Do NOT use any external DS libraries.
 *   Do NOT printf inside library functions.
 *
 * DEPENDENCIES:
 *   #include "dna_utils.h"   (already written — use freely)
 *
 * BUILD & TEST:
 *   gcc -O2 -Wall -std=c11 -DTE_TEST dna_utils.c trie_errors.c -o test_te -lm
 *   ./test_te
 *
 * BUILD OBJECT FILE:
 *   gcc -O2 -Wall -std=c11 -fPIC -c trie_errors.c -o trie_errors.o
 *
 * USED BY FEATURES:
 *   Feature 6 — Mutation detection: match input DNA against known SNPs
 *   Feature 9 — Disease risk: fuzzy-match user snippet → disease labels
 *
 * KEY INSIGHT — WHY FUZZY MATCHING?
 *   A patient's DNA may have 1-2 base differences from the reference SNP
 *   due to sequencing errors or nearby polymorphisms.
 *   Exact matching would miss these — fuzzy matching catches them.
 *   max_errors = 0 → exact only
 *   max_errors = 1 → allow 1 base mismatch (most common clinical use)
 *   max_errors = 2 → allow 2 base mismatches
 */

#include <stddef.h>
#include "dna_utils.h"

/* ============================================================
   CONSTANTS
   ============================================================ */

/*
 * TE_ALPHA_SIZE — alphabet size for trie children.
 * Only A, T, G, C — NO sentinel '$' (SNP sequences don't need it).
 * children[base_to_index(c)] for c in {A,T,G,C}.
 */
#define TE_ALPHA_SIZE    4

/* Maximum allowed errors in fuzzy search */
#define TE_MAX_ERRORS    2

/* Maximum SNP sequence length we store */
#define TE_MAX_SNP_LEN   100

/* Maximum TSV line length when loading from file */
#define TE_MAX_LINE      1024

/* Severity levels */
#define TE_SEV_LOW       0
#define TE_SEV_MEDIUM    1
#define TE_SEV_HIGH      2

/* Return codes */
#define TE_OK            1
#define TE_ERR           0

/* ============================================================
   STRUCTS
   ============================================================ */

/*
 * TrieNode
 * --------
 * One node in the trie.
 *
 * children[i] = child node for base with index i (A=0,T=1,G=2,C=3).
 * NULL means no child for that base.
 *
 * is_end == 1 means this node is the last character of a stored SNP.
 * Only end nodes have non-NULL disease_label and snp_id.
 *
 * MEMORY:
 *   disease_label and snp_id are heap-allocated strings.
 *   trie_free() must free these at every end node.
 */
typedef struct TrieNode {
    struct TrieNode* children[TE_ALPHA_SIZE];
    int              is_end;           /* 1 if this ends a stored SNP        */
    char*            disease_label;    /* heap-allocated, e.g. "Sickle Cell" */
    char*            snp_id;           /* heap-allocated, e.g. "rs334"       */
    int              severity;         /* TE_SEV_LOW / MEDIUM / HIGH         */
} TrieNode;

/*
 * Trie
 * ----
 * The complete trie structure.
 * root is always an empty node (no label, is_end=0).
 * size = number of SNPs successfully inserted.
 */
typedef struct {
    TrieNode* root;
    int       size;
} Trie;

/*
 * TrieMatch
 * ---------
 * One result from trie_search_fuzzy().
 * Represents one SNP that matched the query within max_errors mismatches.
 *
 * snp_id and disease_label are pointers INTO the trie node —
 * do NOT free them individually. They are owned by the trie.
 * Free only the TrieMatch array itself.
 */
typedef struct {
    const char* snp_id;         /* points to trie node string — do not free */
    const char* disease_label;  /* points to trie node string — do not free */
    int         severity;       /* TE_SEV_LOW / MEDIUM / HIGH               */
    int         mismatches;     /* how many bases differed (0, 1, or 2)     */
    int         match_position; /* position in query where match started    */
} TrieMatch;

/*
 * TrieMatchList
 * -------------
 * Return type of trie_search_fuzzy().
 * matches[] is heap-allocated. Free with trie_match_list_free().
 */
typedef struct {
    TrieMatch* matches;   /* heap-allocated array of matches */
    int        count;     /* number of matches found         */
} TrieMatchList;

/* ============================================================
   STACK FRAME — used internally by fuzzy DFS
   ============================================================ */

/*
 * FuzzyFrame
 * ----------
 * One entry on the explicit DFS stack used by trie_search_fuzzy().
 *
 * WHY EXPLICIT STACK?
 *   Recursive DFS risks stack overflow for large tries with max_errors=2.
 *   An explicit stack on the heap is safer and gives us full control.
 *
 * node           — current trie node
 * query_pos      — how far along the query we have consumed
 * errors_used    — how many mismatches used so far on this path
 * match_start    — position in query where this match attempt started
 */
typedef struct {
    TrieNode* node;
    int       query_pos;
    int       errors_used;
    int       match_start;
} FuzzyFrame;

/* ============================================================
   CORE API — YOU MUST IMPLEMENT ALL OF THESE
   ============================================================ */

/* ------------------------------------------------------------
 * trie_create
 * -----------
 * Allocate and initialise an empty Trie.
 *
 * Creates root node with all children NULL, is_end=0.
 * Sets size = 0.
 *
 * Returns heap-allocated Trie* on success, NULL on alloc failure.
 * ------------------------------------------------------------ */
Trie* trie_create(void);


/* ------------------------------------------------------------
 * trie_free
 * ---------
 * Free all memory used by the trie.
 *
 * TRAVERSAL ORDER: post-order DFS (children before parent).
 * For each node:
 *   1. Recursively free all non-NULL children first.
 *   2. If node->is_end: free(node->disease_label), free(node->snp_id).
 *   3. free(node).
 * After all nodes freed: free(trie).
 *
 * Safe to call with NULL (no-op).
 * ------------------------------------------------------------ */
void trie_free(Trie* t);


/* ------------------------------------------------------------
 * trie_insert
 * -----------
 * Insert one SNP marker sequence into the trie.
 *
 * INPUT:
 *   t            — the trie
 *   snp_sequence — DNA sequence of the SNP marker (A,T,G,C only)
 *   disease      — disease name string (will be heap-copied)
 *   snp_id       — SNP identifier string e.g. "rs334" (will be heap-copied)
 *   severity     — TE_SEV_LOW, TE_SEV_MEDIUM, or TE_SEV_HIGH
 *
 * ALGORITHM (standard trie insert):
 *   Start at root.
 *   For each base c in snp_sequence:
 *     idx = base_to_index(c)          ← from dna_utils.h
 *     if node->children[idx] == NULL:
 *       create new TrieNode, set all children NULL, is_end=0
 *       node->children[idx] = new_node
 *     node = node->children[idx]
 *   At final node:
 *     node->is_end = 1
 *     node->disease_label = strdup(disease)
 *     node->snp_id        = strdup(snp_id)
 *     node->severity      = severity
 *   t->size++
 *
 * Does nothing if t or snp_sequence is NULL.
 * Does nothing if snp_sequence fails dna_validate().
 * ------------------------------------------------------------ */
void trie_insert(Trie*       t,
                 const char* snp_sequence,
                 const char* disease,
                 const char* snp_id,
                 int         severity);


/* ------------------------------------------------------------
 * trie_search_exact
 * -----------------
 * Check if an exact SNP sequence is in the trie.
 *
 * INPUT:
 *   t     — the trie
 *   query — DNA sequence to look up
 *
 * OUTPUT:
 *   TE_OK  (1) if query is in the trie (exact match, is_end == 1)
 *   TE_ERR (0) if not found, or on NULL/invalid input
 *
 * ALGORITHM:
 *   Start at root.
 *   For each base c in query:
 *     idx = base_to_index(c)
 *     if node->children[idx] == NULL → return TE_ERR
 *     node = node->children[idx]
 *   return node->is_end ? TE_OK : TE_ERR
 * ------------------------------------------------------------ */
int trie_search_exact(Trie* t, const char* query);


/* ------------------------------------------------------------
 * trie_search_fuzzy
 * -----------------
 * Find all SNPs in the trie that match query within max_errors mismatches.
 *
 * INPUT:
 *   t          — the trie (loaded with SNPs)
 *   query      — DNA sequence to match against (A,T,G,C)
 *   max_errors — maximum mismatches allowed (0, 1, or 2)
 *                clamped to TE_MAX_ERRORS if larger
 *
 * OUTPUT:
 *   Heap-allocated TrieMatchList* with all matches found.
 *   count == 0 if nothing matched.
 *   Returns NULL on error (NULL inputs, invalid query).
 *   CALLER must free with trie_match_list_free().
 *
 * ALGORITHM — Iterative DFS with explicit stack:
 *
 *   Initialise stack with one frame:
 *     { node=root, query_pos=0, errors_used=0, match_start=0 }
 *
 *   While stack not empty:
 *     Pop frame { node, query_pos, errors_used, match_start }
 *
 *     If query_pos == strlen(query):
 *       If node->is_end:
 *         → MATCH FOUND
 *         → add TrieMatch to results:
 *              snp_id        = node->snp_id
 *              disease_label = node->disease_label
 *              severity      = node->severity
 *              mismatches    = errors_used
 *              match_position = match_start
 *       Continue to next frame (don't push more)
 *
 *     Else (still characters to consume):
 *       current_base_idx = base_to_index(query[query_pos])
 *       For each child index c in [0, 1, 2, 3]:
 *         If node->children[c] == NULL: skip
 *         If c == current_base_idx:
 *           → EXACT match: push { children[c], query_pos+1,
 *                                  errors_used, match_start }
 *         Else if errors_used < max_errors:
 *           → MISMATCH: push { children[c], query_pos+1,
 *                               errors_used+1, match_start }
 *
 * STACK SIZE:
 *   Worst case: 4^max_errors * query_length frames.
 *   For max_errors=2, query_length=50 → ~800 frames max.
 *   Allocate stack with initial capacity 256, realloc if needed.
 *
 * NOTE:
 *   match_start is always 0 in the basic implementation
 *   (we search the full query against each SNP).
 *   If you want to find SNPs as substrings of a longer sequence,
 *   call this function for each window of the sequence separately
 *   (the Flask route handles windowing).
 *
 * TIME:  O(4^k * m) where k = max_errors, m = query length
 * ------------------------------------------------------------ */
TrieMatchList* trie_search_fuzzy(Trie*       t,
                                 const char* query,
                                 int         max_errors);


/* ------------------------------------------------------------
 * trie_search_in_sequence
 * -----------------------
 * Find all SNP matches anywhere within a longer DNA sequence.
 * Slides a window of each SNP length along the sequence and
 * calls trie_search_fuzzy() on each window.
 *
 * INPUT:
 *   t          — the trie
 *   sequence   — full DNA sequence to scan (e.g. HBB gene, 626 bp)
 *   max_errors — mismatches allowed per window
 *
 * OUTPUT:
 *   Heap-allocated TrieMatchList* with all matches found.
 *   match_position in each TrieMatch = position in sequence.
 *   CALLER must free with trie_match_list_free().
 *
 * IMPLEMENTATION NOTE:
 *   Since all SNPs in our TSV are the same length (18 bp),
 *   you can use a fixed window size = TE_MAX_SNP_LEN.
 *   For variable-length SNPs, iterate over all possible window sizes.
 * ------------------------------------------------------------ */
TrieMatchList* trie_search_in_sequence(Trie*       t,
                                       const char* sequence,
                                       int         max_errors);


/* ------------------------------------------------------------
 * trie_load_from_file
 * -------------------
 * Load SNP markers from a TSV file into the trie.
 *
 * FILE FORMAT (tab-separated, one SNP per line, no header):
 *   snp_id \t sequence \t disease_name \t severity
 *
 * EXAMPLE LINES:
 *   rs334	GTGCACCTGACTCCTGTG	Sickle Cell Anemia	2
 *   rs28897672	ATGGATTTATCTGCTCTT	Breast Cancer (BRCA1)	2
 *   rs11549407	ATGGTGCACCTGACTCCT	Beta-Thalassemia	1
 *
 * ALGORITHM:
 *   fopen(filepath, "r")
 *   For each line:
 *     Parse 4 fields by splitting on '\t'
 *     severity = atoi(severity_field)   [0, 1, or 2]
 *     trie_insert(t, sequence, disease, snp_id, severity)
 *     count++
 *   fclose
 *
 * Returns number of SNPs successfully loaded.
 * Returns -1 if file cannot be opened.
 * Skips malformed lines silently (lines without exactly 4 fields).
 * ------------------------------------------------------------ */
int trie_load_from_file(Trie* t, const char* filepath);


/* ------------------------------------------------------------
 * trie_match_list_free
 * --------------------
 * Free a TrieMatchList returned by trie_search_fuzzy() or
 * trie_search_in_sequence().
 *
 * IMPORTANT: Do NOT free matches[i].snp_id or matches[i].disease_label.
 * These are pointers into the trie nodes — the trie owns them.
 * Only free the matches[] array and the TrieMatchList struct.
 *
 * Safe to call with NULL.
 * ------------------------------------------------------------ */
void trie_match_list_free(TrieMatchList* list);


/* ============================================================
   JSON SERIALISATION — implement for Flask integration
   ============================================================ */

/* ------------------------------------------------------------
 * trie_match_list_to_json
 * -----------------------
 * Serialise a TrieMatchList to a heap-allocated JSON string.
 *
 * Output format:
 * {
 *   "count": 2,
 *   "matches": [
 *     {
 *       "snp_id":    "rs334",
 *       "disease":   "Sickle Cell Anemia",
 *       "severity":  2,
 *       "mismatches": 0,
 *       "position":  0
 *     },
 *     {
 *       "snp_id":    "rs11549407",
 *       "disease":   "Beta-Thalassemia",
 *       "severity":  1,
 *       "mismatches": 1,
 *       "position":  0
 *     }
 *   ]
 * }
 *
 * Returns NULL on allocation failure or NULL input.
 * CALLER must free() the returned string.
 * ------------------------------------------------------------ */
char* trie_match_list_to_json(TrieMatchList* list);


/* ============================================================
   OPTIONAL DEBUG HELPERS
   Only compiled when TE_DEBUG is defined.
   ============================================================ */
#ifdef TE_DEBUG
/* Print all SNPs stored in the trie */
void trie_print_all(Trie* t);

/* Print info about one TrieMatch */
void trie_match_print(TrieMatch* m);
#endif

#endif /* TRIE_ERRORS_H */