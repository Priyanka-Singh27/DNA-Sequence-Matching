#include <stdlib.h>
#include <string.h>
#include "suffix_tree.h"

/* =========================
   CREATE NODE
   ========================= */
static SuffixTreeNode* create_node(int start, int* end) {
    SuffixTreeNode* node = malloc(sizeof(SuffixTreeNode));

    node->start = start;
    node->end = end;
    node->suffix_index = -1;

    for (int i = 0; i < ST_ALPHA_SIZE; i++)
        node->children[i] = NULL;

    return node;
}

/* =========================
   INSERT SUFFIX (NAIVE)
   ========================= */
static void insert_suffix(SuffixTree* st, int pos) {
    SuffixTreeNode* current = st->root;
    int i = pos;

    while (i < st->length) {
        int idx = base_to_index(st->text[i]);

        /* No edge → create leaf */
        if (!current->children[idx]) {
            int* end = malloc(sizeof(int));
            *end = st->length - 1;

            current->children[idx] = create_node(i, end);
            current->children[idx]->suffix_index = pos;
            return;
        }

        /* Traverse existing edge */
        SuffixTreeNode* next = current->children[idx];
        int k = next->start;

        while (k <= *(next->end) &&
               i < st->length &&
               st->text[k] == st->text[i]) {
            k++;
            i++;
        }

        /* Mismatch → split edge */
        if (k <= *(next->end)) {
            int* split_end = malloc(sizeof(int));
            *split_end = k - 1;

            SuffixTreeNode* split = create_node(next->start, split_end);
            current->children[idx] = split;

            /* Adjust old node */
            next->start = k;
            split->children[base_to_index(st->text[k])] = next;

            /* New leaf */
            int* new_end = malloc(sizeof(int));
            *new_end = st->length - 1;

            split->children[base_to_index(st->text[i])] =
                create_node(i, new_end);

            return;
        }

        current = next;
    }
}

/* =========================
   BUILD TREE
   ========================= */
SuffixTree* suffix_tree_build(const char* dna, int length) {
    if (!dna || length <= 0 || length > ST_MAX_SEQ_LEN)
        return NULL;

    if (!dna_validate(dna))
        return NULL;

    SuffixTree* st = malloc(sizeof(SuffixTree));

    st->length = length + 1;
    st->text = malloc(st->length + 1);

    strcpy(st->text, dna);
    st->text[length] = '$';
    st->text[length + 1] = '\0';

    int* root_end = malloc(sizeof(int));
    *root_end = -1;

    st->root = create_node(-1, root_end);

    /* Insert all suffixes */
    for (int i = 0; i < st->length; i++)
        insert_suffix(st, i);

    return st;
}

/* =========================
   DFS COLLECT LEAVES
   ========================= */
static void collect_leaves(SuffixTreeNode* node, int* res, int* idx) {
    if (!node) return;

    if (node->suffix_index >= 0) {
        res[(*idx)++] = node->suffix_index;
        return;
    }

    for (int i = 0; i < ST_ALPHA_SIZE; i++)
        collect_leaves(node->children[i], res, idx);
}

/* =========================
   SEARCH
   ========================= */
STSearchResult* suffix_tree_search(SuffixTree* st,
                                   const char* pattern) {

    if (!st || !pattern || !dna_validate(pattern))
        return NULL;

    SuffixTreeNode* current = st->root;
    int i = 0;

    while (pattern[i]) {
        int idx = base_to_index(pattern[i]);

        if (!current->children[idx]) {
            STSearchResult* res = malloc(sizeof(STSearchResult));
            res->positions = NULL;
            res->count = 0;
            return res;
        }

        SuffixTreeNode* next = current->children[idx];
        int k = next->start;

        while (k <= *(next->end) && pattern[i]) {
            if (st->text[k] != pattern[i]) {
                STSearchResult* res = malloc(sizeof(STSearchResult));
                res->positions = NULL;
                res->count = 0;
                return res;
            }
            k++;
            i++;
        }

        current = next;
    }

    int* positions = malloc(sizeof(int) * st->length);
    int count = 0;

    collect_leaves(current, positions, &count);

    STSearchResult* res = malloc(sizeof(STSearchResult));
    res->positions = positions;
    res->count = count;

    return res;
}

/* =========================
   FREE TREE
   ========================= */
static void free_node(SuffixTreeNode* node) {
    if (!node) return;

    for (int i = 0; i < ST_ALPHA_SIZE; i++)
        free_node(node->children[i]);

    free(node->end);
    free(node);
}

void suffix_tree_free(SuffixTree* st) {
    if (!st) return;

    free_node(st->root);
    free(st->text);
    free(st);
}

/* =========================
   FREE SEARCH RESULT
   ========================= */
void suffix_tree_search_result_free(STSearchResult* result) {
    if (!result) return;
    free(result->positions);
    free(result);
}