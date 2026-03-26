/*
 * trie_errors.c
 * =============
 * Trie with fuzzy string matching. Simplified structure.
 */

#define _POSIX_C_SOURCE 200809L
#include "trie_errors.h"
#include "dna_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAX_RESULTS 1024

/* ============================================================
   CREATE & FREE
   ============================================================ */

Trie* trie_create(void) {
    Trie* t = (Trie*)calloc(1, sizeof(Trie));
    t->root = (TrieNode*)calloc(1, sizeof(TrieNode));
    return t;
}

static void node_free(TrieNode* node) {
    if (!node) return;
    for (int i = 0; i < TE_ALPHA_SIZE; i++) node_free(node->children[i]);
    if (node->is_end) {
        free(node->disease_label);
        free(node->snp_id);
    }
    free(node);
}

void trie_free(Trie* t) {
    if (!t) return;
    node_free(t->root);
    free(t);
}

/* ============================================================
   INSERT & EXACT SEARCH
   ============================================================ */

void trie_insert(Trie* t, const char* seq, const char* disease, const char* snp_id, int severity) {
    if (!t || !seq || !dna_validate(seq)) return;
    
    TrieNode* node = t->root;
    for (const char* p = seq; *p; p++) {
        int idx = base_to_index(toupper((unsigned char)*p));
        if (idx < 0 || idx >= TE_ALPHA_SIZE) return;
        if (!node->children[idx]) node->children[idx] = (TrieNode*)calloc(1, sizeof(TrieNode));
        node = node->children[idx];
    }
    
    if (node->is_end) {
        free(node->disease_label);
        free(node->snp_id);
    } else {
        node->is_end = 1;
        t->size++;
    }
    
    node->disease_label = strdup(disease);
    node->snp_id        = strdup(snp_id);
    node->severity      = severity;
}

int trie_search_exact(Trie* t, const char* query) {
    if (!t || !query || !dna_validate(query)) return TE_ERR;
    TrieNode* node = t->root;
    for (const char* p = query; *p; p++) {
        int idx = base_to_index(toupper((unsigned char)*p));
        if (idx < 0 || idx >= TE_ALPHA_SIZE || !node->children[idx]) return TE_ERR;
        node = node->children[idx];
    }
    return node->is_end ? TE_OK : TE_ERR;
}

/* ============================================================
   RECURSIVE FUZZY SEARCH
   ============================================================ */

static void dfs_recursive(TrieNode* node, const char* q, int pos, int len, 
                          int errs_used, int max_errs, TrieMatchList* list) {
    if (!node || list->count >= MAX_RESULTS) return;

    if (pos == len) {
        if (node->is_end) {
            TrieMatch m;
            m.disease_label = node->disease_label;
            m.snp_id = node->snp_id;
            m.severity = node->severity;
            m.mismatches = errs_used;
            m.match_position = 0;
            list->matches[list->count++] = m;
        }
        return;
    }

    int cur_idx = base_to_index(toupper((unsigned char)q[pos]));
    for (int c = 0; c < TE_ALPHA_SIZE; c++) {
        if (!node->children[c]) continue;
        dfs_recursive(node->children[c], q, pos + 1, len, errs_used + (c == cur_idx ? 0 : 1), max_errs, list);
    }
}

TrieMatchList* trie_search_fuzzy(Trie* t, const char* query, int max_errors) {
    if (!t || !query || !dna_validate(query)) return NULL;

    TrieMatchList* list = (TrieMatchList*)malloc(sizeof(TrieMatchList));
    list->matches = (TrieMatch*)malloc(MAX_RESULTS * sizeof(TrieMatch));
    list->count = 0;

    dfs_recursive(t->root, query, 0, strlen(query), 0, max_errors < 0 ? 0 : max_errors, list);
    return list;
}

/* ============================================================
   SEQUENCE SCANNER
   ============================================================ */

TrieMatchList* trie_search_in_sequence(Trie* t, const char* seq, int max_errors) {
    if (!t || !seq || !dna_validate(seq)) return NULL;

    TrieMatchList* all = (TrieMatchList*)malloc(sizeof(TrieMatchList));
    all->matches = (TrieMatch*)malloc(MAX_RESULTS * sizeof(TrieMatch));
    all->count = 0;

    int seq_len = strlen(seq);
    for (int wlen = 1; wlen <= TE_MAX_SNP_LEN; wlen++) {
        for (int start = 0; start + wlen <= seq_len; start++) {
            char window[64];
            memcpy(window, seq + start, wlen);
            window[wlen] = '\0';

            TrieMatchList* wr = trie_search_fuzzy(t, window, max_errors);
            if (wr) {
                for (int i = 0; i < wr->count && all->count < MAX_RESULTS; i++) {
                    TrieMatch m = wr->matches[i];
                    m.match_position = start;
                    
                    int dup = 0;
                    for (int j = 0; j < all->count; j++) {
                        if (all->matches[j].match_position == start && all->matches[j].snp_id == m.snp_id) {
                            if (m.mismatches < all->matches[j].mismatches) all->matches[j] = m;
                            dup = 1; break;
                        }
                    }
                    if (!dup) all->matches[all->count++] = m;
                }
                trie_match_list_free(wr);
            }
        }
    }
    return all;
}

/* ============================================================
   FILE IO & CLEANUP
   ============================================================ */

int trie_load_from_file(Trie* t, const char* filepath) {
    FILE* fp = fopen(filepath, "r");
    if (!fp) return -1;

    char line[TE_MAX_LINE];
    int count = 0;
    while (fgets(line, sizeof(line), fp)) {
        size_t len = strlen(line);
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) line[--len] = '\0';
        if (len == 0) continue;

        char* s_id = strtok(line, "\t");
        char* seq  = strtok(NULL, "\t");
        char* dis  = strtok(NULL, "\t");
        char* sev  = strtok(NULL, "\t");

        if (s_id && seq && dis && sev) {
            trie_insert(t, seq, dis, s_id, atoi(sev));
            count++;
        }
    }
    fclose(fp);
    return count;
}

void trie_match_list_free(TrieMatchList* list) {
    if (!list) return;
    free(list->matches);
    free(list);
}

/* ============================================================
   JSON EXPORT
   ============================================================ */

char* trie_match_list_to_json(TrieMatchList* list) {
    if (!list) return NULL;

    int capacity = 512 * 1024; /* Static 512 KB buffer is massive overkill */
    char* buf = (char*)malloc(capacity);
    int written = snprintf(buf, capacity, "{\"count\":%d,\"matches\":[", list->count);

    for (int i = 0; i < list->count; i++) {
        TrieMatch* m = &list->matches[i];
        written += snprintf(buf + written, capacity - written, "%s{\"snp_id\":\"%s\",\"disease\":\"%s\",\"severity\":%d,\"mismatches\":%d,\"position\":%d}",
            i > 0 ? "," : "", m->snp_id, m->disease_label, m->severity, m->mismatches, m->match_position);
    }
    
    snprintf(buf + written, capacity - written, "]}");
    return buf;
}
