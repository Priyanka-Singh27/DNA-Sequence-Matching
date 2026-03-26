#include "union_find.h"
#include <stdlib.h>

UnionFind* union_find_create(int n, char** species_names) { return NULL; }
void union_find_free(UnionFind* uf) {}
int union_find_find(UnionFind* uf, int x) { return 0; }
int union_find_union(UnionFind* uf, int x, int y, double similarity) { return 0; }
int union_find_connected(UnionFind* uf, int x, int y) { return 0; }
ClusterList* union_find_get_clusters(UnionFind* uf) { return NULL; }
MergeEvent* union_find_get_merge_log(UnionFind* uf, int* count) { *count = 0; return NULL; }
void union_find_cluster_species(UnionFind* uf, double** similarity_matrix, double threshold) {}
void cluster_list_free(ClusterList* list) {}
char* union_find_clusters_to_json(UnionFind* uf) { return NULL; }
char* union_find_merge_log_to_json(UnionFind* uf) { return NULL; }
