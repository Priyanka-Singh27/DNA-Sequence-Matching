/*
 * needleman_wunsch.c
 * ------------------
 * Needleman-Wunsch global sequence alignment.
 * Simplified variant.
 */

#define _POSIX_C_SOURCE 200809L
#include "needleman_wunsch.h"
#include "dna_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#define DIR_DIAG 1
#define DIR_UP   2
#define DIR_LEFT 3

static inline int max3(int a, int b, int c) {
    int m = a > b ? a : b;
    return m > c ? m : c;
}

NwParams nw_default_params(void) {
    NwParams p = { NW_DEFAULT_MATCH, NW_DEFAULT_MISMATCH, NW_DEFAULT_GAP };
    return p;
}

void nw_result_free(AlignmentResult* r) {
    if (r) { free(r->aligned_a); free(r->aligned_b); free(r); }
}

AlignmentResult* nw_align(const char* seq_a, const char* seq_b, NwParams p) {
    if (!seq_a || !seq_b || !dna_validate(seq_a) || !dna_validate(seq_b)) return NULL;
    
    size_t n = strlen(seq_a), m = strlen(seq_b);
    if (n == 0 || m == 0 || n > NW_MAX_FULL_LEN || m > NW_MAX_FULL_LEN) return NULL;

    int* dp = (int*)malloc((n+1) * (m+1) * sizeof(int));
    char* dir = (char*)malloc((n+1) * (m+1));
    if (!dp || !dir) { free(dp); free(dir); return NULL; }

    for (size_t i = 0; i <= n; i++) { dp[i*(m+1)] = i * p.gap; dir[i*(m+1)] = DIR_UP; }
    for (size_t j = 0; j <= m; j++) { dp[j] = j * p.gap; dir[j] = DIR_LEFT; }
    dp[0] = 0; dir[0] = 0;

    for (size_t i = 1; i <= n; i++) {
        for (size_t j = 1; j <= m; j++) {
            int diag = dp[(i-1)*(m+1)+(j-1)] + (seq_a[i-1] == seq_b[j-1] ? p.match : p.mismatch);
            int up   = dp[(i-1)*(m+1)+j] + p.gap;
            int left = dp[i*(m+1)+(j-1)] + p.gap;
            int best = max3(diag, up, left);
            dp[i*(m+1)+j] = best;
            dir[i*(m+1)+j] = (best == diag) ? DIR_DIAG : (best == up ? DIR_UP : DIR_LEFT);
        }
    }

    size_t max_len = n + m + 1, pos = max_len - 1, i = n, j = m;
    char* tmp_a = (char*)malloc(max_len);
    char* tmp_b = (char*)malloc(max_len);
    tmp_a[pos] = '\0'; tmp_b[pos] = '\0';

    int matches = 0, mismatches = 0, gaps = 0;
    while (i > 0 || j > 0) {
        char d = dir[i*(m+1)+j];
        pos--;
        if (d == DIR_DIAG && i > 0 && j > 0) {
            tmp_a[pos] = seq_a[--i]; tmp_b[pos] = seq_b[--j];
            if (tmp_a[pos] == tmp_b[pos]) matches++; else mismatches++;
        } else if (d == DIR_UP && i > 0) {
            tmp_a[pos] = seq_a[--i]; tmp_b[pos] = NW_GAP_CHAR; gaps++;
        } else {
            tmp_a[pos] = NW_GAP_CHAR; tmp_b[pos] = seq_b[--j]; gaps++;
        }
    }

    AlignmentResult* res = (AlignmentResult*)calloc(1, sizeof(AlignmentResult));
    res->aligned_a = strdup(tmp_a + pos);
    res->aligned_b = strdup(tmp_b + pos);
    res->score = dp[n*(m+1)+m];
    res->matches = matches;
    res->mismatches = mismatches;
    res->gaps = gaps;
    res->length = max_len - 1 - pos;
    res->identity = res->length > 0 ? ((double)matches / res->length) * 100.0 : 0;
    res->similarity = res->length > 0 ? ((double)(matches + mismatches) / res->length) * 100.0 : 0;
    res->params = p;

    free(dp); free(dir); free(tmp_a); free(tmp_b);
    return res;
}

int nw_score_only(const char* seq_a, const char* seq_b, NwParams p) {
    if (!seq_a || !seq_b || !dna_validate(seq_a) || !dna_validate(seq_b)) return INT_MIN;
    size_t n = strlen(seq_a), m = strlen(seq_b);
    if (n == 0 || m == 0) return INT_MIN;

    const char *longer = (m <= n) ? seq_a : seq_b, *shorter = (m <= n) ? seq_b : seq_a;
    size_t long_len = (m <= n) ? n : m, short_len = (m <= n) ? m : n;

    int *prev = (int*)malloc((short_len + 1) * sizeof(int)), *curr = (int*)malloc((short_len + 1) * sizeof(int));
    for (size_t j = 0; j <= short_len; j++) prev[j] = (int)j * p.gap;

    for (size_t i = 1; i <= long_len; i++) {
        curr[0] = (int)i * p.gap;
        for (size_t j = 1; j <= short_len; j++) {
            int diag = prev[j-1] + (longer[i-1] == shorter[j-1] ? p.match : p.mismatch);
            int up = prev[j] + p.gap, left = curr[j-1] + p.gap;
            curr[j] = max3(diag, up, left);
        }
        int* tmp = prev; prev = curr; curr = tmp;
    }
    
    int score = prev[short_len];
    free(prev); free(curr);
    return score;
}

AlignmentResult* nw_align_default(const char* seq_a, const char* seq_b) {
    return nw_align(seq_a, seq_b, nw_default_params());
}

double nw_identity(const char* seq_a, const char* seq_b) {
    AlignmentResult* r = nw_align_default(seq_a, seq_b);
    if (!r) return -1.0;
    double id = r->identity;
    nw_result_free(r);
    return id;
}

char* nw_result_to_json(const AlignmentResult* r) {
    if (!r) return NULL;
    size_t sz = r->length * 2 + 512;
    char* buf = (char*)malloc(sz);
    snprintf(buf, sz, "{\"aligned_a\":\"%s\",\"aligned_b\":\"%s\",\"score\":%d,\"matches\":%d,\"mismatches\":%d,\"gaps\":%d,\"length\":%d,\"identity\":%.2f,\"similarity\":%.2f}",
        r->aligned_a ? r->aligned_a : "", r->aligned_b ? r->aligned_b : "", r->score, r->matches, r->mismatches, r->gaps, r->length, r->identity, r->similarity);
    return buf;
}