#include "suffix_tree.h"
#include <stdlib.h>

SuffixTree* suffix_tree_build(const char* dna, int length) { return NULL; }
STSearchResult* suffix_tree_search(SuffixTree* st, const char* pattern) { return NULL; }
char* suffix_tree_lcs(const char* seq1, const char* seq2) { return NULL; }
STRepeat* suffix_tree_find_repeats(SuffixTree* st, int min_length, int* count) { *count = 0; return NULL; }
void suffix_tree_free(SuffixTree* st) {}
void suffix_tree_search_result_free(STSearchResult* result) {}
void suffix_tree_repeats_free(STRepeat* repeats, int count) {}
char* suffix_tree_search_to_json(STSearchResult* result, const char* pattern) { return NULL; }
char* suffix_tree_lcs_to_json(const char* lcs) { return NULL; }
