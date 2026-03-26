#include "segment_tree.h"
#include <stdlib.h>

SegmentTree* segment_tree_build(int* positions, char** diseases, int* severities, int count) { return NULL; }
int segment_tree_update(SegmentTree* st, int real_pos, const char* disease, int severity) { return 0; }
SegmentQueryResult* segment_tree_query(SegmentTree* st, int left_real, int right_real, int* count) { *count = 0; return NULL; }
SegmentQueryResult* segment_tree_point_query(SegmentTree* st, int real_pos) { return NULL; }
int segment_tree_load_clinvar(SegmentTree** st_out, const char* filepath) { *st_out = NULL; return 0; }
void segment_tree_results_free(SegmentQueryResult* results, int count) {}
void segment_tree_free(SegmentTree* st) {}
char* segment_tree_results_to_json(SegmentQueryResult* results, int count) { return NULL; }
int sgt_lower_bound(int* arr, int n, int val) { return 0; }
int sgt_upper_bound(int* arr, int n, int val) { return 0; }
