#ifndef UNION_FIND_H
#define UNION_FIND_H

/*
 * union_find.h
 * ============
 * Union-Find (Disjoint Set Union) with path compression and union by rank.
 * Used to progressively cluster species by evolutionary similarity.
 * Records every merge event so D3.js can animate the phylogenetic tree.
 *
 * YOUR JOB (union_find.c):
 *   Implement every function declared in this header.
 *   Do NOT modify this header file.
 *   Do NOT use any external DS libraries.
 *   Do NOT printf inside library functions.
 *
 * DEPENDENCIES:
 *   #include "dna_utils.h"   (already written — use freely)
 *
 * BUILD & TEST:
 *   gcc -O2 -Wall -std=c11 -DUF_TEST dna_utils.c union_find.c -o test_uf -lm
 *   ./test_uf
 *
 * BUILD OBJECT FILE:
 *   gcc -O2 -Wall -std=c11 -fPIC -c union_find.c -o union_find.o
 *
 * USED BY FEATURES:
 *   Feature 3 — Human evolution: cluster human + Neanderthal + Chimp + Gorilla
 *   Feature 4 — Cross-species: pairwise similarity matrix clusters
 *   Feature 8 — Phylogenetic tree: merge log drives D3.js tree animation
 *
 * HOW IT BUILDS THE PHYLOGENETIC TREE:
 *
 *   Input: similarity matrix (from NW scores)
 *     Human-Neanderthal : 99.7%
 *     Human-Chimp       : 98.7%
 *     Human-Gorilla     : 98.3%
 *     Human-Mouse       : 85.0%
 *     Human-Zebrafish   : 70.0%
 *
 *   Process (union_find_cluster_species):
 *     Sort all pairs by similarity descending
 *     Union highest pairs first:
 *       Step 1: Union(Human, Neanderthal, 99.7) → MergeEvent recorded
 *       Step 2: Union(Human, Chimp, 98.7)       → MergeEvent recorded
 *       Step 3: Union(Human, Gorilla, 98.3)     → MergeEvent recorded
 *       Step 4: Union(Primates, Mouse, 85.0)    → MergeEvent recorded
 *       Step 5: Union(All, Zebrafish, 70.0)     → MergeEvent recorded
 *
 *   Output: merge_log[] array
 *     D3.js reads this log → draws nodes and edges progressively
 *     Each step becomes one branch in the phylogenetic tree
 */

#include <stddef.h>
#include "dna_utils.h"

/* ============================================================
   CONSTANTS
   ============================================================ */

/* Maximum number of species the union-find can hold */
#define UF_MAX_SPECIES    100

/* Return codes */
#define UF_OK     1
#define UF_ERR    0

/* ============================================================
   STRUCTS
   ============================================================ */

/*
 * MergeEvent
 * ----------
 * Records one Union operation — used by D3.js to draw the tree.
 *
 * species_a, species_b — indices of the two species that were merged.
 *                         NOT their roots — the original species indices.
 * similarity           — NW identity score at time of merge (0.0-100.0).
 * step                 — sequential step number (1, 2, 3, ...).
 *
 * The frontend uses this to animate:
 *   "At step 1, Human and Neanderthal merged at 99.7% similarity"
 *   → draw an edge between Human and Neanderthal nodes
 */
typedef struct {
    int    species_a;    /* index of first species in the merge  */
    int    species_b;    /* index of second species in the merge */
    double similarity;   /* NW identity % at time of merge       */
    int    step;         /* merge step number (1-indexed)        */
} MergeEvent;

/*
 * SpeciesPair
 * -----------
 * Used internally by union_find_cluster_species() to sort
 * all pairwise similarities before merging.
 */
typedef struct {
    int    a;             /* index of first species  */
    int    b;             /* index of second species */
    double similarity;    /* NW identity score       */
} SpeciesPair;

/*
 * UnionFind
 * ---------
 * The complete Union-Find structure.
 *
 * parent[i]            — parent of species i. If parent[i] == i, i is a root.
 * rank[i]              — rank of species i (used for union by rank).
 *                        Rank is an upper bound on subtree height.
 * species_names[i]     — heap-allocated name of species i.
 * merge_similarity[i]  — similarity score at which species i was merged
 *                        into another set. 0.0 if not yet merged.
 * n                    — total number of species.
 * merge_log            — heap-allocated array of MergeEvent, size n-1 max.
 * merge_count          — number of merges that have happened so far.
 */
typedef struct {
    int*        parent;             /* parent[i] = parent of i    */
    int*        rank;               /* rank[i] = tree height est  */
    char**      species_names;      /* names[i] = species i name  */
    double*     merge_similarity;   /* similarity when i merged   */
    int         n;                  /* total species count        */
    MergeEvent* merge_log;          /* log of all merges          */
    int         merge_count;        /* how many merges done       */
} UnionFind;

/*
 * ClusterList
 * -----------
 * Return type of union_find_get_clusters().
 * Contains all current clusters as groups of species indices.
 */
typedef struct {
    int** groups;         /* groups[i] = array of species indices in cluster i */
    int*  group_sizes;    /* group_sizes[i] = number of species in cluster i   */
    int   cluster_count;  /* total number of clusters                           */
} ClusterList;

/* ============================================================
   CORE API — YOU MUST IMPLEMENT ALL OF THESE
   ============================================================ */

/* ------------------------------------------------------------
 * union_find_create
 * -----------------
 * Allocate and initialise a UnionFind for n species.
 *
 * INPUT:
 *   n             — number of species (must be > 0, <= UF_MAX_SPECIES)
 *   species_names — array of n strings (will be heap-copied)
 *
 * INITIALISATION:
 *   parent[i]           = i     (each species is its own root)
 *   rank[i]             = 0     (all ranks start at zero)
 *   merge_similarity[i] = 0.0
 *   merge_log           = malloc(n * sizeof(MergeEvent))
 *   merge_count         = 0
 *   species_names[i]    = strdup(species_names[i])
 *
 * Returns heap-allocated UnionFind* on success, NULL on failure.
 * ------------------------------------------------------------ */
UnionFind* union_find_create(int n, char** species_names);


/* ------------------------------------------------------------
 * union_find_free
 * ---------------
 * Free all memory used by the UnionFind.
 *
 * Steps:
 *   1. For each i: free(uf->species_names[i])
 *   2. free(uf->species_names)
 *   3. free(uf->parent)
 *   4. free(uf->rank)
 *   5. free(uf->merge_similarity)
 *   6. free(uf->merge_log)
 *   7. free(uf)
 *
 * Safe to call with NULL.
 * ------------------------------------------------------------ */
void union_find_free(UnionFind* uf);


/* ------------------------------------------------------------
 * union_find_find
 * ---------------
 * Find the root (representative) of the set containing species x.
 *
 * INPUT:
 *   uf — the UnionFind
 *   x  — species index (0-indexed, must be < uf->n)
 *
 * OUTPUT:
 *   Root index of x's set.
 *   Returns -1 on invalid input.
 *
 * ALGORITHM — PATH COMPRESSION (two-step / halving):
 *
 *   While parent[x] != x:
 *     parent[x] = parent[parent[x]]   ← skip one level (path halving)
 *     x = parent[x]
 *   return x
 *
 * WHY PATH COMPRESSION?
 *   Without it, a chain of n unions takes O(n) per find.
 *   With path compression, amortised time per operation = O(α(n)) ≈ O(1).
 *   α = inverse Ackermann function — essentially constant for all n.
 *
 * ALTERNATIVE (full path compression — also acceptable):
 *   int root = x;
 *   while (parent[root] != root) root = parent[root];
 *   while (parent[x] != root) { int next = parent[x]; parent[x] = root; x = next; }
 *   return root;
 * ------------------------------------------------------------ */
int union_find_find(UnionFind* uf, int x);


/* ------------------------------------------------------------
 * union_find_union
 * ----------------
 * Merge the sets containing species x and species y.
 *
 * INPUT:
 *   uf         — the UnionFind
 *   x          — first species index
 *   y          — second species index
 *   similarity — NW identity score between x and y (0.0-100.0)
 *
 * OUTPUT:
 *   UF_OK  (1) if merge happened (x and y were in different sets)
 *   UF_ERR (0) if already in same set OR invalid input
 *
 * ALGORITHM — UNION BY RANK:
 *
 *   Step 1 — Find roots:
 *     rx = union_find_find(uf, x)
 *     ry = union_find_find(uf, y)
 *     if rx == ry → already same set → return UF_ERR
 *
 *   Step 2 — Union by rank (attach smaller tree under larger):
 *     if rank[rx] < rank[ry]:
 *       parent[rx] = ry           (rx becomes child of ry)
 *     else if rank[rx] > rank[ry]:
 *       parent[ry] = rx           (ry becomes child of rx)
 *     else:
 *       parent[ry] = rx           (arbitrary: rx becomes root)
 *       rank[rx]++                (rank increases only on tie)
 *
 *   WHY UNION BY RANK?
 *     Keeps tree height O(log n) in the worst case.
 *     Combined with path compression → O(α(n)) per operation.
 *
 *   Step 3 — Record merge event:
 *     merge_log[merge_count].species_a  = x
 *     merge_log[merge_count].species_b  = y
 *     merge_log[merge_count].similarity = similarity
 *     merge_log[merge_count].step       = merge_count + 1
 *     merge_count++
 *
 *   Step 4 — Store merge similarity on the new root:
 *     merge_similarity[new_root] = similarity
 *
 *   Return UF_OK.
 * ------------------------------------------------------------ */
int union_find_union(UnionFind* uf, int x, int y, double similarity);


/* ------------------------------------------------------------
 * union_find_connected
 * --------------------
 * Check if two species are in the same cluster.
 *
 * Returns 1 if find(x) == find(y), 0 otherwise.
 * Returns 0 on invalid input.
 * ------------------------------------------------------------ */
int union_find_connected(UnionFind* uf, int x, int y);


/* ------------------------------------------------------------
 * union_find_get_clusters
 * -----------------------
 * Get all current clusters as groups of species indices.
 *
 * OUTPUT:
 *   Heap-allocated ClusterList* with all current groups.
 *   Each group is a heap-allocated int[] of species indices.
 *   CALLER must call cluster_list_free() on result.
 *
 * ALGORITHM:
 *   Step 1 — For each species i: root = find(i)
 *             Group species by root.
 *
 *   Step 2 — Build ClusterList:
 *             cluster_count = number of distinct roots
 *             groups[k]     = array of species indices with the same root
 *             group_sizes[k]= length of groups[k]
 *
 * NOTE: The cluster IDs are not stable — they change as merges happen.
 * ------------------------------------------------------------ */
ClusterList* union_find_get_clusters(UnionFind* uf);


/* ------------------------------------------------------------
 * union_find_get_merge_log
 * ------------------------
 * Get a pointer to the internal merge log array.
 *
 * OUTPUT:
 *   Pointer to uf->merge_log (NOT a copy — owned by UnionFind).
 *   Sets *count = uf->merge_count.
 *   Returns NULL if uf is NULL or no merges have happened.
 *
 * CALLER must NOT free the returned pointer.
 * The log is valid as long as the UnionFind exists.
 * ------------------------------------------------------------ */
MergeEvent* union_find_get_merge_log(UnionFind* uf, int* count);


/* ------------------------------------------------------------
 * union_find_cluster_species
 * --------------------------
 * Main driver: cluster all species using a similarity matrix.
 * This is what Feature 8 (phylogenetic tree) calls.
 *
 * INPUT:
 *   uf                — the UnionFind (already created with n species)
 *   similarity_matrix — n×n matrix where matrix[i][j] = NW identity
 *                       between species i and species j (0.0-100.0).
 *                       matrix[i][i] = 100.0 (diagonal).
 *                       matrix[i][j] == matrix[j][i] (symmetric).
 *   threshold         — minimum similarity to merge two species (0.0-100.0).
 *                       Pairs below threshold are not merged.
 *                       Typical value: 70.0 (merge anything above 70% similar)
 *
 * ALGORITHM:
 *   Step 1 — Collect all unique pairs (i < j) into SpeciesPair[]:
 *              pair.a = i, pair.b = j, pair.similarity = matrix[i][j]
 *
 *   Step 2 — Sort pairs by similarity DESCENDING.
 *              (Highest similarity pairs merged first — builds tree correctly)
 *
 *   Step 3 — For each pair in sorted order:
 *              if pair.similarity >= threshold:
 *                union_find_union(uf, pair.a, pair.b, pair.similarity)
 *
 *   After this, uf->merge_log contains the complete merge history
 *   that D3.js uses to draw the phylogenetic tree.
 *
 * NOTE: similarity_matrix is a 2D array passed as double**.
 *       Access as similarity_matrix[i][j].
 * ------------------------------------------------------------ */
void union_find_cluster_species(UnionFind* uf,
                                 double**   similarity_matrix,
                                 double     threshold);


/* ------------------------------------------------------------
 * cluster_list_free
 * -----------------
 * Free a ClusterList returned by union_find_get_clusters().
 *
 * For each i: free(list->groups[i])
 * free(list->groups)
 * free(list->group_sizes)
 * free(list)
 *
 * Safe to call with NULL.
 * ------------------------------------------------------------ */
void cluster_list_free(ClusterList* list);


/* ============================================================
   JSON SERIALISATION — implement for Flask integration
   ============================================================ */

/* ------------------------------------------------------------
 * union_find_clusters_to_json
 * ---------------------------
 * Serialise the current clusters AND merge log to JSON.
 * This is the main output consumed by the frontend.
 *
 * Output format:
 * {
 *   "merge_log": [
 *     { "species_a": "Human",   "species_b": "Neanderthal",
 *       "similarity": 99.7,     "step": 1 },
 *     { "species_a": "Human",   "species_b": "Chimp",
 *       "similarity": 98.7,     "step": 2 },
 *     { "species_a": "Human",   "species_b": "Gorilla",
 *       "similarity": 98.3,     "step": 3 }
 *   ],
 *   "clusters": [
 *     { "id": 0, "members": ["Human", "Neanderthal", "Chimp", "Gorilla"] },
 *     { "id": 1, "members": ["Mouse"] },
 *     { "id": 2, "members": ["Zebrafish"] }
 *   ]
 * }
 *
 * Uses uf->merge_log and union_find_get_clusters() internally.
 * Returns NULL on allocation failure or NULL input.
 * CALLER must free() the returned string.
 * ------------------------------------------------------------ */
char* union_find_clusters_to_json(UnionFind* uf);


/* ------------------------------------------------------------
 * union_find_merge_log_to_json
 * ----------------------------
 * Serialise only the merge log to JSON.
 * Used by the phylogenetic tree page for step-by-step animation.
 *
 * Output format:
 * {
 *   "count": 3,
 *   "merges": [
 *     { "species_a": "Human", "species_b": "Neanderthal",
 *       "similarity": 99.7, "step": 1 },
 *     ...
 *   ]
 * }
 *
 * Returns NULL on allocation failure or NULL input.
 * CALLER must free() the returned string.
 * ------------------------------------------------------------ */
char* union_find_merge_log_to_json(UnionFind* uf);


/* ============================================================
   OPTIONAL DEBUG HELPERS
   ============================================================ */
#ifdef UF_DEBUG
/* Print parent[] and rank[] arrays */
void union_find_print(UnionFind* uf);

/* Print all clusters with species names */
void union_find_print_clusters(UnionFind* uf);

/* Print the full merge log */
void union_find_print_merge_log(UnionFind* uf);
#endif

#endif /* UNION_FIND_H */