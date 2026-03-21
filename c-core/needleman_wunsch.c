/*
 * needleman_wunsch.c
 * ------------------
 * Needleman-Wunsch global sequence alignment.
 *
 * Build object file:
 *   gcc -O2 -Wall -std=c11 -fPIC -c needleman_wunsch.c -o needleman_wunsch.o
 *
 * Self-test:
 *   gcc -O2 -Wall -std=c11 -DNW_TEST dna_utils.c needleman_wunsch.c -o test_nw -lm
 *   ./test_nw
 */

#include "needleman_wunsch.h"
#include "dna_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

/* ============================================================
   INTERNAL HELPERS
   ============================================================ */

/* Traceback direction codes */
#define DIR_DIAG  1   /* diagonal: match or mismatch */
#define DIR_UP    2   /* up:       gap in seq_b      */
#define DIR_LEFT  3   /* left:     gap in seq_a      */

static inline int max3(int a, int b, int c) {
    int m = a > b ? a : b;
    return m > c ? m : c;
}

/* ============================================================
   DEFAULT PARAMS
   ============================================================ */

NwParams nw_default_params(void) {
    NwParams p;
    p.match    = NW_DEFAULT_MATCH;
    p.mismatch = NW_DEFAULT_MISMATCH;
    p.gap      = NW_DEFAULT_GAP;
    return p;
}

/* ============================================================
   nw_result_free
   ============================================================ */

void nw_result_free(AlignmentResult* r) {
    if (r == NULL) return;
    free(r->aligned_a);
    free(r->aligned_b);
    free(r);
}

/* ============================================================
   nw_align  (full DP + traceback)
   ============================================================ */

AlignmentResult* nw_align(const char* seq_a,
                           const char* seq_b,
                           NwParams    params)
{
    /* input guards */
    if (!seq_a || !seq_b)     return NULL;
    if (!dna_validate(seq_a)) return NULL;
    if (!dna_validate(seq_b)) return NULL;

    size_t n = strlen(seq_a);
    size_t m = strlen(seq_b);

    if (n == 0 || m == 0)                            return NULL;
    if (n > NW_MAX_FULL_LEN || m > NW_MAX_FULL_LEN) return NULL;

    /*
     * DP matrix flat layout: dp[i*(m+1)+j]
     * Direction matrix same layout, stored as char (1 byte each).
     */
    size_t cells = (n + 1) * (m + 1);
    int*  dp  = (int*)  malloc(cells * sizeof(int));
    char* dir = (char*) malloc(cells * sizeof(char));

    if (!dp || !dir) { free(dp); free(dir); return NULL; }

    /* initialise borders */
    for (size_t i = 0; i <= n; i++) {
        dp [i*(m+1)]  = (int)i * params.gap;
        dir[i*(m+1)]  = DIR_UP;
    }
    for (size_t j = 0; j <= m; j++) {
        dp [j]  = (int)j * params.gap;
        dir[j]  = DIR_LEFT;
    }
    dp[0]  = 0;
    dir[0] = 0;

    /* fill */
    for (size_t i = 1; i <= n; i++) {
        for (size_t j = 1; j <= m; j++) {
            int s    = (seq_a[i-1] == seq_b[j-1]) ? params.match : params.mismatch;
            int diag = dp[(i-1)*(m+1)+(j-1)] + s;
            int up   = dp[(i-1)*(m+1)+ j   ] + params.gap;
            int left = dp[ i   *(m+1)+(j-1)] + params.gap;
            int best = max3(diag, up, left);

            dp [i*(m+1)+j] = best;
            dir[i*(m+1)+j] = (best == diag) ? DIR_DIAG
                            : (best == up)   ? DIR_UP
                                             : DIR_LEFT;
        }
    }

    int final_score = dp[n*(m+1)+m];

    /* traceback — build reversed alignment */
    size_t max_len = n + m + 1;
    char* tmp_a = (char*)malloc(max_len);
    char* tmp_b = (char*)malloc(max_len);
    if (!tmp_a || !tmp_b) {
        free(dp); free(dir); free(tmp_a); free(tmp_b);
        return NULL;
    }

    size_t pos = 0;
    size_t i = n, j = m;

    while (i > 0 || j > 0) {
        char d = dir[i*(m+1)+j];
        if (d == DIR_DIAG && i > 0 && j > 0) {
            tmp_a[pos] = seq_a[i-1];
            tmp_b[pos] = seq_b[j-1];
            i--; j--;
        } else if (d == DIR_UP && i > 0) {
            tmp_a[pos] = seq_a[i-1];
            tmp_b[pos] = NW_GAP_CHAR;
            i--;
        } else {
            tmp_a[pos] = NW_GAP_CHAR;
            tmp_b[pos] = seq_b[j-1];
            j--;
        }
        pos++;
    }

    free(dp);
    free(dir);

    size_t aln_len = pos;

    /* allocate result */
    AlignmentResult* res = (AlignmentResult*)calloc(1, sizeof(AlignmentResult));
    if (!res) { free(tmp_a); free(tmp_b); return NULL; }

    res->aligned_a = (char*)malloc(aln_len + 1);
    res->aligned_b = (char*)malloc(aln_len + 1);
    if (!res->aligned_a || !res->aligned_b) {
        nw_result_free(res); free(tmp_a); free(tmp_b);
        return NULL;
    }

    /* reverse into result */
    for (size_t k = 0; k < aln_len; k++) {
        res->aligned_a[k] = tmp_a[aln_len - 1 - k];
        res->aligned_b[k] = tmp_b[aln_len - 1 - k];
    }
    res->aligned_a[aln_len] = '\0';
    res->aligned_b[aln_len] = '\0';

    free(tmp_a);
    free(tmp_b);

    /* statistics */
    int matches = 0, mismatches = 0, gaps = 0;
    for (size_t k = 0; k < aln_len; k++) {
        char ca = res->aligned_a[k], cb = res->aligned_b[k];
        if      (ca == NW_GAP_CHAR || cb == NW_GAP_CHAR) gaps++;
        else if (ca == cb)                                matches++;
        else                                              mismatches++;
    }

    res->score      = final_score;
    res->matches    = matches;
    res->mismatches = mismatches;
    res->gaps       = gaps;
    res->length     = (int)aln_len;
    res->identity   = aln_len > 0 ? ((double)matches  / aln_len) * 100.0 : 0.0;
    res->similarity = aln_len > 0 ? ((double)(matches + mismatches) / aln_len) * 100.0 : 0.0;
    res->params     = params;

    return res;
}

/* ============================================================
   nw_score_only  (two-row space-optimised DP)
   ============================================================ */

int nw_score_only(const char* seq_a,
                  const char* seq_b,
                  NwParams    params)
{
    if (!seq_a || !seq_b)     return INT_MIN;
    if (!dna_validate(seq_a)) return INT_MIN;
    if (!dna_validate(seq_b)) return INT_MIN;

    size_t n = strlen(seq_a);
    size_t m = strlen(seq_b);
    if (n == 0 || m == 0)     return INT_MIN;

    /* keep shorter sequence as columns to minimise memory */
    const char* longer;
    const char* shorter;
    size_t long_len, short_len;

    if (m <= n) {
        shorter = seq_b; short_len = m;
        longer  = seq_a; long_len  = n;
    } else {
        shorter = seq_a; short_len = n;
        longer  = seq_b; long_len  = m;
    }

    int* prev = (int*)malloc((short_len + 1) * sizeof(int));
    int* curr = (int*)malloc((short_len + 1) * sizeof(int));
    if (!prev || !curr) { free(prev); free(curr); return INT_MIN; }

    /* initialise first row */
    for (size_t j = 0; j <= short_len; j++)
        prev[j] = (int)j * params.gap;

    /* fill row by row */
    for (size_t i = 1; i <= long_len; i++) {
        curr[0] = (int)i * params.gap;
        for (size_t j = 1; j <= short_len; j++) {
            int s    = (longer[i-1] == shorter[j-1]) ? params.match : params.mismatch;
            int diag = prev[j-1] + s;
            int up   = prev[j  ] + params.gap;
            int left = curr[j-1] + params.gap;
            curr[j]  = max3(diag, up, left);
        }
        int* tmp = prev; prev = curr; curr = tmp;
    }

    int score = prev[short_len];
    free(prev);
    free(curr);
    return score;
}

/* ============================================================
   CONVENIENCE WRAPPERS
   ============================================================ */

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

/* ============================================================
   JSON SERIALISATION
   ============================================================ */

char* nw_result_to_json(const AlignmentResult* r) {
    if (!r) return NULL;

    /* buf_size: aligned strings + overhead for all numeric fields */
    size_t buf_size = (size_t)(r->length) * 2 + 1024;
    char*  buf      = (char*)malloc(buf_size);
    if (!buf) return NULL;

    int written = snprintf(buf, buf_size,
        "{"
          "\"aligned_a\":\"%s\","
          "\"aligned_b\":\"%s\","
          "\"score\":%d,"
          "\"matches\":%d,"
          "\"mismatches\":%d,"
          "\"gaps\":%d,"
          "\"length\":%d,"
          "\"identity\":%.2f,"
          "\"similarity\":%.2f"
        "}",
        r->aligned_a  ? r->aligned_a  : "",
        r->aligned_b  ? r->aligned_b  : "",
        r->score,
        r->matches,
        r->mismatches,
        r->gaps,
        r->length,
        r->identity,
        r->similarity
    );

    if (written < 0 || (size_t)written >= buf_size) {
        free(buf);
        return NULL;
    }
    return buf;
}


/* ============================================================
   SELF-TEST
   ============================================================ */

#ifdef NW_TEST

#include <assert.h>
#include <math.h>

static int tests_run = 0, tests_passed = 0;

#define TEST(name) do { tests_run++; printf("  %-58s ", name); fflush(stdout); } while(0)
#define PASS()     do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(msg)  do { printf("FAIL -- %s\n", msg); } while(0)

static void test_params(void) {
    printf("\n[default params]\n");
    NwParams p = nw_default_params();
    TEST("match = 2");    assert(p.match == 2);    PASS();
    TEST("mismatch = -1");assert(p.mismatch == -1);PASS();
    TEST("gap = -2");     assert(p.gap == -2);     PASS();
}

static void test_identical(void) {
    printf("\n[identical sequences]\n");
    NwParams p = nw_default_params();

    TEST("score = len * match_score");
    AlignmentResult* r = nw_align("ATGC", "ATGC", p);
    assert(r && r->score == 8); nw_result_free(r); PASS();

    TEST("identity = 100.0");
    r = nw_align("ATGC", "ATGC", p);
    assert(r && r->identity == 100.0); nw_result_free(r); PASS();

    TEST("no gaps in output");
    r = nw_align("ATGC", "ATGC", p);
    assert(r && strchr(r->aligned_a, '-') == NULL);
    nw_result_free(r); PASS();
}

static void test_sickle_cell(void) {
    printf("\n[sickle cell SNP (A->T at HBB codon 6)]\n");
    const char* normal = "ATGGTGCACCTGACTCCTGAGGAGAAGTCT";
    const char* sickle = "ATGGTGCACCTGACTCCTGTGGAGAAGTCT";
    NwParams p = nw_default_params();

    TEST("exactly 1 mismatch");
    AlignmentResult* r = nw_align(normal, sickle, p);
    assert(r && r->mismatches == 1 && r->matches == 29 && r->gaps == 0);
    nw_result_free(r); PASS();

    TEST("identity = 96.67%");
    r = nw_align(normal, sickle, p);
    assert(r && fabs(r->identity - 96.67) < 0.01);
    nw_result_free(r); PASS();

    TEST("score_only matches full score");
    r = nw_align(normal, sickle, p);
    int fs = r->score; nw_result_free(r);
    assert(nw_score_only(normal, sickle, p) == fs); PASS();
}

static void test_gaps(void) {
    printf("\n[gap insertion]\n");
    NwParams p = nw_default_params();

    TEST("single deletion introduces a gap");
    AlignmentResult* r = nw_align("ATGCAT", "ATCAT", p);
    assert(r && r->gaps > 0 && r->length == 6);
    nw_result_free(r); PASS();

    TEST("aligned strings always same length");
    r = nw_align("ATGCAT", "ATCAT", p);
    assert(r && strlen(r->aligned_a) == strlen(r->aligned_b));
    nw_result_free(r); PASS();

    TEST("no position has gap in both strands simultaneously");
    r = nw_align("ATGCAT", "ATCAT", p);
    for (int k = 0; k < r->length; k++)
        assert(!(r->aligned_a[k] == '-' && r->aligned_b[k] == '-'));
    nw_result_free(r); PASS();
}

static void test_score_only(void) {
    printf("\n[nw_score_only]\n");
    NwParams p = nw_default_params();

    TEST("score_only == full score (identical)");
    AlignmentResult* r = nw_align("ATGCATGC", "ATGCATGC", p);
    int s1 = r->score; nw_result_free(r);
    assert(nw_score_only("ATGCATGC", "ATGCATGC", p) == s1); PASS();

    TEST("score_only handles asymmetric lengths");
    int s = nw_score_only("ATGCATGCATGC", "ATGC", p);
    assert(s != INT_MIN); PASS();

    TEST("score_only NULL -> INT_MIN");
    assert(nw_score_only(NULL, "ATGC", p) == INT_MIN); PASS();

    TEST("score_only invalid base -> INT_MIN");
    assert(nw_score_only("ATXC", "ATGC", p) == INT_MIN); PASS();
}

static void test_edge_cases(void) {
    printf("\n[edge cases]\n");
    NwParams p = nw_default_params();

    TEST("NULL seq_a -> NULL");
    assert(nw_align(NULL, "ATGC", p) == NULL); PASS();

    TEST("NULL seq_b -> NULL");
    assert(nw_align("ATGC", NULL, p) == NULL); PASS();

    TEST("invalid base -> NULL");
    assert(nw_align("ATXC", "ATGC", p) == NULL); PASS();

    TEST("single base identical: score = match");
    AlignmentResult* r = nw_align("A", "A", p);
    assert(r && r->score == p.match); nw_result_free(r); PASS();

    TEST("single base mismatch: score = mismatch");
    r = nw_align("A", "T", p);
    assert(r && r->score == p.mismatch); nw_result_free(r); PASS();

    TEST("custom params applied correctly");
    NwParams c = {5, -3, -4};
    r = nw_align("AAAA", "AAAA", c);
    assert(r && r->score == 20); nw_result_free(r); PASS();

    TEST("nw_result_free(NULL) safe");
    nw_result_free(NULL); PASS();
}

static void test_json(void) {
    printf("\n[JSON serialisation]\n");

    TEST("output non-NULL");
    AlignmentResult* r = nw_align_default("ATGC", "ATGC");
    char* j = nw_result_to_json(r);
    assert(j != NULL); nw_result_free(r); free(j); PASS();

    TEST("contains aligned_a key");
    r = nw_align_default("ATGC", "ATGC");
    j = nw_result_to_json(r);
    assert(strstr(j, "\"aligned_a\"")); nw_result_free(r); free(j); PASS();

    TEST("contains score key");
    r = nw_align_default("ATGC", "ATGC");
    j = nw_result_to_json(r);
    assert(strstr(j, "\"score\"")); nw_result_free(r); free(j); PASS();

    TEST("contains identity key");
    r = nw_align_default("ATGC", "ATGC");
    j = nw_result_to_json(r);
    assert(strstr(j, "\"identity\"")); nw_result_free(r); free(j); PASS();

    TEST("NULL result -> NULL JSON");
    assert(nw_result_to_json(NULL) == NULL); PASS();
}

static void test_convenience(void) {
    printf("\n[convenience wrappers]\n");

    TEST("nw_identity identical = 100.0");
    assert(nw_identity("ATGC", "ATGC") == 100.0); PASS();

    TEST("nw_identity NULL -> -1.0");
    assert(nw_identity(NULL, "ATGC") == -1.0); PASS();

    TEST("nw_identity HBB ~96.67");
    double id = nw_identity("ATGGTGCACCTGACTCCTGAGGAGAAGTCT",
                             "ATGGTGCACCTGACTCCTGTGGAGAAGTCT");
    assert(fabs(id - 96.67) < 0.01); PASS();
}

int main(void) {
    printf("=======================================================\n");
    printf("  GenomeX -- needleman_wunsch.c self-test\n");
    printf("=======================================================\n");

    test_params();
    test_identical();
    test_sickle_cell();
    test_gaps();
    test_score_only();
    test_edge_cases();
    test_json();
    test_convenience();

    printf("\n=======================================================\n");
    printf("  Results: %d / %d tests passed\n", tests_passed, tests_run);
    printf("=======================================================\n");
    return (tests_passed == tests_run) ? 0 : 1;
}

#endif /* NW_TEST */