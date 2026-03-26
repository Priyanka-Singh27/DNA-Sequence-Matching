#include "interval_tree.h"
#include <stdlib.h>

IntervalTree* interval_tree_create(void) { return NULL; }
void interval_tree_free(IntervalTree* it) {}
void interval_tree_insert(IntervalTree* it, int low, int high, const char* gene_name, const char* chromosome) {}
IntervalResult* interval_tree_query(IntervalTree* it, int query_low, int query_high, int* count) { *count = 0; return NULL; }
IntervalResult* interval_tree_point_query(IntervalTree* it, int position, int* count) { *count = 0; return NULL; }
int interval_tree_load_genes(IntervalTree* it, const char* filepath) { return 0; }
void interval_tree_results_free(IntervalResult* results, int count) {}
char* interval_tree_results_to_json(IntervalResult* results, int count) { return NULL; }
