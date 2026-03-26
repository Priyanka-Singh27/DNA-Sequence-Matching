#ifndef SUFFIX_TREE_H
#define SUFFIX_TREE_H

#include <stddef.h>
#include "dna_utils.h"

/* =========================
   CONSTANTS
   ========================= */

#define ST_ALPHA_SIZE 5      /* A T G C $ */
#define ST_MAX_SEQ_LEN 1000000

/* =========================
   NODE STRUCTURE
   ========================= */

/*
 * Each node represents an edge label using:
 * text[start ... *end]
 */
typedef struct SuffixTreeNode {
    int start;                          /* start index of substring */
    int* end;                           /* end index of substring   */
    int suffix_index;                   /* leaf: >=0, internal: -1  */
    struct SuffixTreeNode* children[ST_ALPHA_SIZE];
} SuffixTreeNode;

/* =========================
   TREE STRUCTURE
   ========================= */

typedef struct {
    SuffixTreeNode* root;   /* root node */
    char* text;             /* original DNA + '$' */
    int length;             /* length including '$' */
} SuffixTree;

/* =========================
   SEARCH RESULT
   ========================= */

typedef struct {
    int* positions;   /* array of match positions */
    int  count;       /* number of matches */
} STSearchResult;

/* =========================
   CORE FUNCTIONS
   ========================= */

/* Build suffix tree (NAIVE O(n^2)) */
SuffixTree* suffix_tree_build(const char* dna, int length);

/* Search pattern */
STSearchResult* suffix_tree_search(SuffixTree* st,
                                   const char* pattern);

/* Free tree */
void suffix_tree_free(SuffixTree* st);

/* Free search result */
void suffix_tree_search_result_free(STSearchResult* result);

#endif