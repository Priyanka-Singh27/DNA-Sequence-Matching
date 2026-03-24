#ifndef SUFFIX_TREE_H
#define SUFFIX_TREE_H

/*
 * suffix_tree.h
 * =============
 * Suffix Tree built using Ukkonen's O(n) online algorithm.
 *
 * YOUR JOB (suffix_tree.c):
 *   Implement every function declared in this header.
 *   Do NOT modify this header file.
 *   Do NOT use any external DS libraries.
 *   Do NOT printf inside library functions.
 *
 * ALGORITHM: Ukkonen's O(n) suffix tree construction.
 *   Reference: https://cp-algorithms.com/string/suffix-tree-ukkonen.html
 *
 * DEPENDENCIES:
 *   #include "dna_utils.h"   (already written — use freely)
 *
 * BUILD & TEST:
 *   gcc -O2 -Wall -std=c11 -DST_TEST dna_utils.c suffix_tree.c -o test_st -lm
 *   ./test_st
 *
 * BUILD OBJECT FILE (for libdna):
 *   gcc -O2 -Wall -std=c11 -fPIC -c suffix_tree.c -o suffix_tree.o
 *
 * USED BY FEATURES:
 *   Feature 1  — LCS between normal and mutated sequences
 *   Feature 5  — Pathogen identification (search unknown seq vs database)
 *   Feature 7  — Find conserved/matching blocks in alignment viewer
 */

#include <stddef.h>
#include "dna_utils.h"

/* ============================================================
   CONSTANTS
   ============================================================ */

/*
 * ST_ALPHA_SIZE — number of children per node.
 * We support: A=0, T=1, G=2, C=3, $=4
 * $ is the sentinel character appended to every sequence before building.
 * It makes all suffixes explicit (terminates every suffix at a leaf).
 */
#define ST_ALPHA_SIZE    5

/*
 * ST_INF — used as "infinity" for the end of a leaf edge.
 * Leaf edges in Ukkonen's algorithm extend to the end of the string.
 * We represent this as ST_INF — updated globally via leaf_end pointer.
 */
#define ST_INF           1000000000

/*
 * ST_MAX_SEQ_LEN — maximum input sequence length we support.
 * Sequences longer than this are rejected by suffix_tree_build().
 */
#define ST_MAX_SEQ_LEN   1000000    /* 1 million bases */

/* Return codes */
#define ST_OK            1
#define ST_ERR           0

/* ============================================================
   STRUCTS
   ============================================================ */

/*
 * SuffixTreeNode
 * --------------
 * One node in the suffix tree.
 *
 * EDGE LABEL:
 *   The edge from parent → this node is labelled by text[start .. *end].
 *   For internal nodes, end is a fixed value set at creation.
 *   For leaf nodes, end points to the global leaf_end integer in SuffixTree,
 *   so all leaves automatically extend when a new character is processed.
 *
 * SUFFIX INDEX:
 *   -1 for internal nodes.
 *   >= 0 for leaf nodes: the start position of the suffix in text[].
 *
 * SUFFIX LINK:
 *   Only used by internal nodes during Ukkonen's construction.
 *   Points to the node representing the same string minus its first character.
 *   Makes the algorithm O(n) instead of O(n^2).
 *
 * CHILDREN:
 *   children[base_to_index(c)] = child node reached by character c.
 *   children[4] = child reached by '$' sentinel.
 *   NULL means no child for that character.
 */
typedef struct SuffixTreeNode {
    int   start;                           /* start index of edge label in text[] */
    int*  end;                             /* pointer to end index of edge label   */
    int   suffix_index;                    /* -1=internal node, >=0=leaf node      */
    struct SuffixTreeNode* children[ST_ALPHA_SIZE];
    struct SuffixTreeNode* suffix_link;    /* Ukkonen suffix link (internal nodes) */
} SuffixTreeNode;

/*
 * SuffixTree
 * ----------
 * The complete suffix tree structure.
 *
 * text[]     : the input DNA string with '$' appended.
 *              text is heap-allocated — free it in suffix_tree_free().
 * length     : strlen(text) — includes the '$' sentinel.
 * leaf_end   : global end pointer shared by ALL leaf nodes.
 *              Ukkonen's algorithm increments this each phase.
 * size       : total node count (useful for debugging/stats).
 *
 * HOW leaf_end WORKS:
 *   All leaf nodes store end = &(tree->leaf_end).
 *   When leaf_end is incremented, ALL leaves extend simultaneously.
 *   This is what makes Ukkonen's algorithm O(n) — Rule 1 extensions
 *   happen implicitly just by incrementing leaf_end.
 */
typedef struct {
    SuffixTreeNode* root;     /* root node — always an internal node */
    char*           text;     /* input string + '$', heap-allocated  */
    int             length;   /* length of text[] including '$'      */
    int             leaf_end; /* global end for all leaf nodes       */
    int             size;     /* total number of nodes in tree       */
} SuffixTree;

/*
 * STSearchResult
 * --------------
 * Return type of suffix_tree_search().
 * Contains all positions where a pattern was found.
 */
typedef struct {
    int* positions;   /* heap-allocated array of match start positions */
    int  count;       /* number of matches found                       */
} STSearchResult;

/*
 * STRepeat
 * --------
 * One repeated substring found by suffix_tree_find_repeats().
 */
typedef struct {
    char* sequence;   /* heap-allocated repeat string                  */
    int   length;     /* length of the repeat                          */
    int   count;      /* how many times it appears in the text         */
} STRepeat;

/* ============================================================
   CORE API — YOU MUST IMPLEMENT ALL OF THESE
   ============================================================ */

/* ------------------------------------------------------------
 * suffix_tree_build
 * -----------------
 * Build a suffix tree from a DNA sequence using Ukkonen's algorithm.
 *
 * INPUT:
 *   dna    — valid DNA string (A, T, G, C only). Must be null-terminated.
 *            dna_validate() is called internally — invalid chars → NULL.
 *   length — strlen(dna). Must be > 0 and <= ST_MAX_SEQ_LEN.
 *
 * OUTPUT:
 *   Heap-allocated SuffixTree* on success.
 *   NULL on failure (invalid input, length == 0, alloc failure).
 *
 * ALGORITHM OUTLINE (Ukkonen's):
 *
 *   Step 0 — Append '$' to dna → store as tree->text.
 *
 *   Step 1 — Create root node. root->suffix_link = root.
 *
 *   Step 2 — Initialise active point:
 *              active_node   = root
 *              active_edge   = 0     (index into text[])
 *              active_length = 0
 *
 *   Step 3 — For each character text[i] (i = 0 .. length):
 *     a. Increment leaf_end (extends all leaves — Rule 1)
 *     b. Set last_new_internal = NULL
 *     c. remaining++
 *     d. While remaining > 0:
 *          - Look at active_edge character
 *          - RULE 3: if text[active_node->children[active_edge]->start
 *                             + active_length] == text[i]
 *                     → active_length++, break   (suffix already in tree)
 *          - RULE 2: otherwise, split or create leaf:
 *                     → if active_length == 0: add leaf from active_node
 *                     → else: split edge, create internal node + new leaf
 *                     → if last_new_internal != NULL:
 *                          last_new_internal->suffix_link = new_internal
 *                     → last_new_internal = new_internal
 *                     → follow suffix link or move active_node to root
 *                     → remaining--
 *
 *   Step 4 — DFS to set suffix_index on all leaf nodes.
 *
 * IMPORTANT DETAILS:
 *   - All leaf nodes share &(tree->leaf_end) as their end pointer.
 *   - Internal node end pointers are heap-allocated ints (not shared).
 *   - After construction, call a DFS to fill suffix_index on leaves.
 *     suffix_index = length - depth_of_leaf (where depth = edge label sum).
 *   - root->suffix_link = root (self-link, used as sentinel).
 *
 * TIME:  O(n)
 * SPACE: O(n)
 * ------------------------------------------------------------ */
SuffixTree* suffix_tree_build(const char* dna, int length);


/* ------------------------------------------------------------
 * suffix_tree_search
 * ------------------
 * Find all positions in the original sequence where pattern occurs.
 *
 * INPUT:
 *   st      — suffix tree built by suffix_tree_build().
 *   pattern — DNA pattern to search for. Must be valid DNA.
 *
 * OUTPUT:
 *   Heap-allocated STSearchResult* with positions[] and count.
 *   positions[] contains 0-based indices into the ORIGINAL dna
 *   (not including the '$' sentinel).
 *   count == 0 if pattern not found (positions == NULL).
 *   Returns NULL on error (NULL inputs, invalid pattern).
 *
 * ALGORITHM:
 *   Phase 1 — Traverse:
 *     Start at root. For each character in pattern:
 *       - Get child: node = current->children[base_to_index(pattern[i])]
 *       - If child is NULL → pattern not found → return count=0
 *       - Walk along the edge, matching pattern characters
 *       - If mismatch on edge → pattern not found → return count=0
 *       - If edge exhausted → move to child node, continue
 *     After consuming all pattern characters, we are at a node (or mid-edge).
 *
 *   Phase 2 — Collect:
 *     DFS from the node we landed on.
 *     Collect suffix_index from every leaf node reached.
 *     Those suffix_index values are all positions where pattern occurs.
 *
 *   Phase 3 — Sort positions[] (ascending) before returning.
 *
 * CALLER must free result->positions and the result struct itself.
 *
 * TIME:  O(m + k) where m = pattern length, k = number of matches
 * ------------------------------------------------------------ */
STSearchResult* suffix_tree_search(SuffixTree*  st,
                                   const char*  pattern);


/* ------------------------------------------------------------
 * suffix_tree_lcs
 * ---------------
 * Find the Longest Common Substring of two DNA sequences.
 *
 * INPUT:
 *   seq1 — first DNA sequence (valid A,T,G,C)
 *   seq2 — second DNA sequence (valid A,T,G,C)
 *
 * OUTPUT:
 *   Heap-allocated char* containing the LCS string.
 *   Returns "" (empty heap string, not NULL) if no common substring.
 *   Returns NULL on error (NULL input, invalid bases, alloc failure).
 *   CALLER must free() the returned string.
 *
 * ALGORITHM (Generalised Suffix Tree):
 *   Step 1 — Build combined string:
 *              combined = seq1 + '#' + seq2 + '$'
 *              Use dna_concat_with_separator(seq1, seq2, '#') from dna_utils.h
 *              '#' is not A/T/G/C so it acts as a separator.
 *
 *   Step 2 — Build one suffix tree on combined string.
 *
 *   Step 3 — DFS, marking each node as:
 *              LEFT  : subtree has a leaf with suffix_index < len(seq1)
 *              RIGHT : subtree has a leaf with suffix_index > len(seq1)+1
 *              BOTH  : subtree has leaves from BOTH seq1 and seq2
 *                      → this node represents a common substring candidate
 *
 *   Step 4 — Among all BOTH nodes, find the one with greatest
 *              string depth (sum of edge lengths from root to node).
 *              That depth = length of the longest common substring.
 *
 *   Step 5 — Extract the substring from combined[node_start .. node_start+depth].
 *              Make sure to stop before '#' or '$'.
 *
 * EXAMPLE:
 *   seq1 = "ATGCATGC"
 *   seq2 = "ATGCTGCA"
 *   combined = "ATGCATGC#ATGCTGCA$"
 *   LCS = "ATGC" (length 4)
 *
 * TIME:  O(n + m)
 * SPACE: O(n + m)
 * ------------------------------------------------------------ */
char* suffix_tree_lcs(const char* seq1, const char* seq2);


/* ------------------------------------------------------------
 * suffix_tree_find_repeats
 * ------------------------
 * Find all substrings of length >= min_length that appear
 * more than once in the sequence.
 *
 * INPUT:
 *   st         — suffix tree
 *   min_length — minimum repeat length to report (must be >= 1)
 *   count      — OUTPUT: set to number of repeats found
 *
 * OUTPUT:
 *   Heap-allocated STRepeat* array of length *count.
 *   Each entry has sequence (heap-allocated), length, count.
 *   Returns NULL and *count=0 if none found or on error.
 *   CALLER must free each repeat->sequence and the array itself.
 *
 * ALGORITHM:
 *   DFS through suffix tree.
 *   Every INTERNAL node represents a substring that occurs >= 2 times
 *   (because it has >= 2 children).
 *   If the string depth of an internal node >= min_length:
 *     → extract the substring it represents
 *     → count = number of leaf nodes in its subtree
 *     → add to results
 *
 * ------------------------------------------------------------ */
STRepeat* suffix_tree_find_repeats(SuffixTree* st,
                                   int         min_length,
                                   int*        count);


/* ------------------------------------------------------------
 * suffix_tree_free
 * ----------------
 * Free all memory used by the suffix tree.
 *
 * Traversal order: post-order DFS (children before parent).
 * For each node:
 *   - If internal node: free node->end (heap-allocated int)
 *   - If leaf node: do NOT free node->end (it points to tree->leaf_end)
 *   - Free the node struct itself
 * After all nodes freed:
 *   - free(tree->text)
 *   - free(tree)
 *
 * Safe to call with NULL (no-op).
 * ------------------------------------------------------------ */
void suffix_tree_free(SuffixTree* st);


/* ------------------------------------------------------------
 * suffix_tree_search_result_free
 * --------------------------------
 * Free an STSearchResult returned by suffix_tree_search().
 * Safe to call with NULL.
 * ------------------------------------------------------------ */
void suffix_tree_search_result_free(STSearchResult* result);


/* ------------------------------------------------------------
 * suffix_tree_repeats_free
 * ------------------------
 * Free an STRepeat array returned by suffix_tree_find_repeats().
 * count = number of entries in the array.
 * Safe to call with NULL.
 * ------------------------------------------------------------ */
void suffix_tree_repeats_free(STRepeat* repeats, int count);


/* ============================================================
   JSON SERIALISATION — implement these for Flask integration
   ============================================================ */

/* ------------------------------------------------------------
 * suffix_tree_search_to_json
 * --------------------------
 * Serialise an STSearchResult to a heap-allocated JSON string.
 *
 * Output format:
 * {
 *   "pattern":   "ATG",
 *   "count":     3,
 *   "positions": [0, 45, 203]
 * }
 *
 * Returns NULL on allocation failure or NULL input.
 * CALLER must free() the returned string.
 * ------------------------------------------------------------ */
char* suffix_tree_search_to_json(STSearchResult* result,
                                 const char*     pattern);


/* ------------------------------------------------------------
 * suffix_tree_lcs_to_json
 * -----------------------
 * Serialise an LCS result to a heap-allocated JSON string.
 *
 * Output format:
 * {
 *   "lcs":    "ATGCATGC",
 *   "length": 8
 * }
 *
 * Returns NULL on allocation failure or NULL input.
 * CALLER must free() the returned string.
 * ------------------------------------------------------------ */
char* suffix_tree_lcs_to_json(const char* lcs);


/* ============================================================
   OPTIONAL DEBUG HELPERS
   Only compiled when ST_DEBUG is defined.
   Do NOT use in production / Flask code.
   ============================================================ */
#ifdef ST_DEBUG
/* Print the tree structure to stdout (for debugging) */
void suffix_tree_print(SuffixTree* st);

/* Print all suffixes stored in the tree */
void suffix_tree_print_suffixes(SuffixTree* st);

/* Return total node count */
int suffix_tree_node_count(SuffixTree* st);
#endif

#endif /* SUFFIX_TREE_H */