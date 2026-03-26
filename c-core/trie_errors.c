/*
 * trie_errors.c
 * =============
 * Trie with fuzzy/approximate string matching for SNP detection.
 *
 * Build object file:
 *   gcc -O2 -Wall -std=c11 -fPIC -c trie_errors.c -o trie_errors.o
 *
 * Self-test:
 *   gcc -O2 -Wall -std=c11 -DTE_TEST dna_utils.c trie_errors.c -o test_te -lm
 *   ./test_te
 */

#define _POSIX_C_SOURCE 200809L  /* enables strdup in strict C11 mode */

#include "trie_errors.h"
#include "dna_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* ============================================================
   INTERNAL HELPERS
   ============================================================ */

/* Allocate and zero-initialise a single TrieNode */
static TrieNode* node_create(void) {
    TrieNode* n = (TrieNode*)calloc(1, sizeof(TrieNode));
    /* calloc zeroes all fields:
       children[] = {NULL,NULL,NULL,NULL}
       is_end = 0, disease_label = NULL, snp_id = NULL, severity = 0 */
    return n;
}

/* ============================================================
   trie_create
   ============================================================ */

Trie* trie_create(void) {
    Trie* t = (Trie*)malloc(sizeof(Trie));
    if (!t) return NULL;

    t->root = node_create();
    if (!t->root) { free(t); return NULL; }

    t->size = 0;
    return t;
}

/* ============================================================
   trie_free  (post-order recursive helper)
   ============================================================ */

static void node_free(TrieNode* node) {
    if (!node) return;

    /* Free all children first (post-order) */
    for (int i = 0; i < TE_ALPHA_SIZE; i++) {
        node_free(node->children[i]);
    }

    /* Free heap strings at end nodes */
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
   trie_insert
   ============================================================ */

void trie_insert(Trie*       t,
                 const char* snp_sequence,
                 const char* disease,
                 const char* snp_id,
                 int         severity)
{
    /* Guard: reject NULL or invalid inputs */
    if (!t || !snp_sequence || !disease || !snp_id) return;
    if (!dna_validate(snp_sequence)) return;

    TrieNode* node = t->root;

    /* Walk / create path for each base in the sequence */
    for (const char* p = snp_sequence; *p != '\0'; p++) {
        int idx = base_to_index((char)toupper((unsigned char)*p));
        if (idx < 0 || idx >= TE_ALPHA_SIZE) return; /* safety: skip $ */

        if (!node->children[idx]) {
            node->children[idx] = node_create();
            if (!node->children[idx]) return; /* alloc failure */
        }
        node = node->children[idx];
    }

    /* Mark end node with metadata.
       If sequence already existed, free old strings before overwriting. */
    if (node->is_end) {
        if (node->disease_label) { free(node->disease_label); node->disease_label = NULL; }
        if (node->snp_id)        { free(node->snp_id);        node->snp_id        = NULL; }
    } else {
        node->is_end = 1;
        t->size++;
    }

    node->disease_label = strdup(disease);
    node->snp_id        = strdup(snp_id);
    node->severity      = severity;
}

/* ============================================================
   trie_search_exact
   ============================================================ */

int trie_search_exact(Trie* t, const char* query) {
    if (!t || !query || !t->root) return TE_ERR;
    if (!dna_validate(query))     return TE_ERR;

    TrieNode* node = t->root;

    for (const char* p = query; *p != '\0'; p++) {
        int idx = base_to_index((char)toupper((unsigned char)*p));
        if (idx < 0 || idx >= TE_ALPHA_SIZE) return TE_ERR;
        if (!node->children[idx])             return TE_ERR;
        node = node->children[idx];
    }

    return node->is_end ? TE_OK : TE_ERR;
}

/* ============================================================
   DYNAMIC ARRAY — for collecting TrieMatch results
   ============================================================ */

typedef struct {
    TrieMatch* data;
    int        size;
    int        capacity;
} MatchArray;

static MatchArray match_array_create(void) {
    MatchArray a;
    a.capacity = 16;
    a.size     = 0;
    a.data     = (TrieMatch*)malloc(a.capacity * sizeof(TrieMatch));
    return a;
}

static int match_array_push(MatchArray* a, TrieMatch m) {
    if (a->size >= a->capacity) {
        int new_cap = a->capacity * 2;
        TrieMatch* new_data = (TrieMatch*)realloc(a->data,
                                new_cap * sizeof(TrieMatch));
        if (!new_data) return 0; /* alloc failure */
        a->data     = new_data;
        a->capacity = new_cap;
    }
    a->data[a->size++] = m;
    return 1;
}

/* ============================================================
   EXPLICIT DFS STACK — for trie_search_fuzzy
   ============================================================ */

typedef struct {
    FuzzyFrame* data;
    int         top;
    int         capacity;
} FuzzyStack;

static FuzzyStack stack_create(int initial_cap) {
    FuzzyStack s;
    s.capacity = initial_cap;
    s.top      = 0;
    s.data     = (FuzzyFrame*)malloc(s.capacity * sizeof(FuzzyFrame));
    return s;
}

static int stack_push(FuzzyStack* s, FuzzyFrame f) {
    if (s->top >= s->capacity) {
        int new_cap = s->capacity * 2;
        FuzzyFrame* nd = (FuzzyFrame*)realloc(s->data,
                            new_cap * sizeof(FuzzyFrame));
        if (!nd) return 0;
        s->data     = nd;
        s->capacity = new_cap;
    }
    s->data[s->top++] = f;
    return 1;
}

static FuzzyFrame stack_pop(FuzzyStack* s) {
    return s->data[--s->top];
}

/* ============================================================
   trie_search_fuzzy
   ============================================================ */

TrieMatchList* trie_search_fuzzy(Trie*       t,
                                  const char* query,
                                  int         max_errors)
{
    if (!t || !query || !t->root) return NULL;
    if (!dna_validate(query))     return NULL;

    /* Clamp max_errors */
    if (max_errors > TE_MAX_ERRORS) max_errors = TE_MAX_ERRORS;
    if (max_errors < 0)             max_errors = 0;

    int query_len = (int)strlen(query);

    /* Initialise result collector and DFS stack */
    MatchArray  results = match_array_create();
    FuzzyStack  stack   = stack_create(256);

    if (!results.data || !stack.data) {
        free(results.data);
        free(stack.data);
        return NULL;
    }

    /* Push initial frame */
    FuzzyFrame init = { t->root, 0, 0, 0 };
    stack_push(&stack, init);

    /* ---- Iterative DFS ---- */
    while (stack.top > 0) {
        FuzzyFrame frame = stack_pop(&stack);
        TrieNode*  node  = frame.node;

        if (frame.query_pos == query_len) {
            /* Consumed entire query — check if this is a complete SNP */
            if (node->is_end) {
                TrieMatch m;
                m.snp_id         = node->snp_id;
                m.disease_label  = node->disease_label;
                m.severity       = node->severity;
                m.mismatches     = frame.errors_used;
                m.match_position = frame.match_start;
                match_array_push(&results, m);
            }
            continue;
        }

        /* Still have query characters to consume */
        int cur_idx = base_to_index(
                        (char)toupper((unsigned char)query[frame.query_pos]));

        for (int c = 0; c < TE_ALPHA_SIZE; c++) {
            if (!node->children[c]) continue;

            if (c == cur_idx) {
                /* Exact match — no error consumed */
                FuzzyFrame next = {
                    node->children[c],
                    frame.query_pos + 1,
                    frame.errors_used,
                    frame.match_start
                };
                stack_push(&stack, next);
            } else if (frame.errors_used < max_errors) {
                /* Mismatch — consume one error */
                FuzzyFrame next = {
                    node->children[c],
                    frame.query_pos + 1,
                    frame.errors_used + 1,
                    frame.match_start
                };
                stack_push(&stack, next);
            }
            /* else: would exceed max_errors — skip this child */
        }
    }

    free(stack.data);

    /* Build TrieMatchList from results */
    TrieMatchList* list = (TrieMatchList*)malloc(sizeof(TrieMatchList));
    if (!list) { free(results.data); return NULL; }

    list->count   = results.size;
    list->matches = results.data;   /* transfer ownership */

    return list;
}

/* ============================================================
   trie_search_in_sequence
   ============================================================ */

TrieMatchList* trie_search_in_sequence(Trie*       t,
                                        const char* sequence,
                                        int         max_errors)
{
    if (!t || !sequence) return NULL;
    if (!dna_validate(sequence)) return NULL;

    int seq_len = (int)strlen(sequence);

    MatchArray all = match_array_create();
    if (!all.data) return NULL;

    /*
     * Slide windows of every possible SNP length (1..TE_MAX_SNP_LEN)
     * across the sequence.
     *
     * KEY INSIGHT: trie_search_fuzzy() only matches queries whose
     * length equals the length of a stored SNP.  So for each start
     * position we must try EVERY window length, not just one fixed
     * length.  In practice our SNPs are all 18 bp, so only the
     * window of length 18 will ever find anything — but this loop
     * handles variable-length SNPs correctly too.
     */
    for (int wlen = 1; wlen <= TE_MAX_SNP_LEN; wlen++) {
        /* Slide this window length across the whole sequence */
        for (int start = 0; start + wlen <= seq_len; start++) {
            char window[TE_MAX_SNP_LEN + 1];
            memcpy(window, sequence + start, (size_t)wlen);
            window[wlen] = '\0';

            TrieMatchList* wr = trie_search_fuzzy(t, window, max_errors);

            if (wr && wr->count > 0) {
                for (int i = 0; i < wr->count; i++) {
                    TrieMatch m       = wr->matches[i];
                    m.match_position  = start;

                    /* Deduplicate same snp at same position */
                    int dup = 0;
                    for (int j = 0; j < all.size; j++) {
                        if (all.data[j].match_position == start &&
                            all.data[j].snp_id         == m.snp_id) {
                            if (m.mismatches < all.data[j].mismatches)
                                all.data[j] = m;
                            dup = 1;
                            break;
                        }
                    }
                    if (!dup) match_array_push(&all, m);
                }
            }
            trie_match_list_free(wr);
        }
    }

    TrieMatchList* result = (TrieMatchList*)malloc(sizeof(TrieMatchList));
    if (!result) { free(all.data); return NULL; }

    result->count   = all.size;
    result->matches = all.data;
    return result;
}

/* ============================================================
   trie_load_from_file
   ============================================================ */

int trie_load_from_file(Trie* t, const char* filepath) {
    if (!t || !filepath) return -1;

    FILE* fp = fopen(filepath, "r");
    if (!fp) return -1;

    char line[TE_MAX_LINE];
    int  count = 0;

    while (fgets(line, sizeof(line), fp)) {
        /* Strip trailing newline */
        size_t len = strlen(line);
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r'))
            line[--len] = '\0';

        if (len == 0) continue; /* skip empty lines */

        /*
         * Parse tab-separated fields:
         *   snp_id \t sequence \t disease \t severity
         */
        char* fields[4];
        int   field_count = 0;
        char* token = strtok(line, "\t");

        while (token && field_count < 4) {
            fields[field_count++] = token;
            token = strtok(NULL, "\t");
        }

        if (field_count != 4) continue; /* malformed line — skip */

        char* snp_id   = fields[0];
        char* sequence = fields[1];
        char* disease  = fields[2];
        int   severity = atoi(fields[3]);

        /* Clamp severity to valid range */
        if (severity < TE_SEV_LOW)  severity = TE_SEV_LOW;
        if (severity > TE_SEV_HIGH) severity = TE_SEV_HIGH;

        trie_insert(t, sequence, disease, snp_id, severity);
        count++;
    }

    fclose(fp);
    return count;
}

/* ============================================================
   trie_match_list_free
   ============================================================ */

void trie_match_list_free(TrieMatchList* list) {
    if (!list) return;
    /*
     * Do NOT free list->matches[i].snp_id or disease_label —
     * they are owned by the trie nodes.
     */
    free(list->matches);
    free(list);
}

/* ============================================================
   trie_match_list_to_json
   ============================================================ */

char* trie_match_list_to_json(TrieMatchList* list) {
    if (!list) return NULL;

    /*
     * Estimate buffer size:
     *   Per match: ~200 chars (snp_id + disease + numbers + keys)
     *   Overhead:  ~64 chars
     */
    size_t buf_size = (size_t)(list->count * 256) + 128;
    if (buf_size < 256) buf_size = 256;

    char* buf = (char*)malloc(buf_size);
    if (!buf) return NULL;

    int written = 0;

    written += snprintf(buf + written, buf_size - (size_t)written,
                        "{\"count\":%d,\"matches\":[", list->count);

    for (int i = 0; i < list->count; i++) {
        TrieMatch* m = &list->matches[i];

        /* Grow buffer if needed */
        size_t needed = (size_t)written + 300;
        if (needed > buf_size) {
            buf_size = needed * 2;
            char* new_buf = (char*)realloc(buf, buf_size);
            if (!new_buf) { free(buf); return NULL; }
            buf = new_buf;
        }

        written += snprintf(buf + written, buf_size - (size_t)written,
            "%s{"
              "\"snp_id\":\"%s\","
              "\"disease\":\"%s\","
              "\"severity\":%d,"
              "\"mismatches\":%d,"
              "\"position\":%d"
            "}",
            i > 0 ? "," : "",
            m->snp_id        ? m->snp_id        : "",
            m->disease_label ? m->disease_label : "",
            m->severity,
            m->mismatches,
            m->match_position
        );
    }

    /* Close array and object */
    size_t needed = (size_t)written + 4;
    if (needed > buf_size) {
        char* new_buf = (char*)realloc(buf, needed);
        if (!new_buf) { free(buf); return NULL; }
        buf = new_buf;
    }
    snprintf(buf + written, needed - (size_t)written, "]}");

    return buf;
}

/* ============================================================
   DEBUG HELPERS
   ============================================================ */

#ifdef TE_DEBUG

static void node_print_all(TrieNode* node, char* path, int depth) {
    if (!node) return;
    if (node->is_end) {
        path[depth] = '\0';
        printf("  SNP: %-20s  disease: %-30s  id: %-12s  sev: %d\n",
               path, node->disease_label, node->snp_id, node->severity);
    }
    for (int i = 0; i < TE_ALPHA_SIZE; i++) {
        if (node->children[i]) {
            path[depth] = index_to_base(i);
            node_print_all(node->children[i], path, depth + 1);
        }
    }
}

void trie_print_all(Trie* t) {
    if (!t) { printf("[Trie: NULL]\n"); return; }
    printf("[Trie: %d SNPs]\n", t->size);
    char path[TE_MAX_SNP_LEN + 1];
    node_print_all(t->root, path, 0);
}

void trie_match_print(TrieMatch* m) {
    if (!m) { printf("[TrieMatch: NULL]\n"); return; }
    printf("[TrieMatch] snp=%-12s  disease=%-30s  sev=%d  mm=%d  pos=%d\n",
           m->snp_id        ? m->snp_id        : "?",
           m->disease_label ? m->disease_label : "?",
           m->severity, m->mismatches, m->match_position);
}

#endif /* TE_DEBUG */


/* ============================================================
   SELF-TEST
   ============================================================ */

#ifdef TE_TEST

#include <assert.h>
#include <math.h>

static int tests_run = 0, tests_passed = 0;

#define TEST(name) do { \
    tests_run++; \
    printf("  %-60s ", name); \
    fflush(stdout); \
} while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(msg) do { printf("FAIL -- %s\n", msg); } while(0)

/* ---- sample SNP data ---- */
static const char* SNP_HBB_SEQ  = "GTGCACCTGACTCCTGTG";
static const char* SNP_HBB_DIS  = "Sickle Cell Anemia";
static const char* SNP_HBB_ID   = "rs334";

static const char* SNP_BRCA_SEQ = "ATGGATTTATCTGCTCTT";
static const char* SNP_BRCA_DIS = "Breast Cancer (BRCA1)";
static const char* SNP_BRCA_ID  = "rs28897672";

static const char* SNP_THAL_SEQ = "ATGGTGCACCTGACTCCT";
static const char* SNP_THAL_DIS = "Beta-Thalassemia";
static const char* SNP_THAL_ID  = "rs11549407";

/* Helper: load 3 test SNPs */
static Trie* make_test_trie(void) {
    Trie* t = trie_create();
    trie_insert(t, SNP_HBB_SEQ,  SNP_HBB_DIS,  SNP_HBB_ID,  TE_SEV_HIGH);
    trie_insert(t, SNP_BRCA_SEQ, SNP_BRCA_DIS, SNP_BRCA_ID, TE_SEV_HIGH);
    trie_insert(t, SNP_THAL_SEQ, SNP_THAL_DIS, SNP_THAL_ID, TE_SEV_MEDIUM);
    return t;
}

/* ---- test groups ---- */

static void test_create_free(void) {
    printf("\n[create / free]\n");

    TEST("trie_create returns non-NULL");
    Trie* t = trie_create();
    assert(t != NULL); PASS();

    TEST("initial size is 0");
    assert(t->size == 0); PASS();

    TEST("root node exists");
    assert(t->root != NULL); PASS();

    TEST("root has no children");
    for (int i = 0; i < TE_ALPHA_SIZE; i++)
        assert(t->root->children[i] == NULL);
    PASS();

    TEST("trie_free(NULL) is safe");
    trie_free(NULL); PASS();

    TEST("trie_free normal trie");
    trie_free(t); PASS();   /* valgrind verifies no leak */
}

static void test_insert(void) {
    printf("\n[insert]\n");
    Trie* t = trie_create();

    TEST("insert valid SNP increments size");
    trie_insert(t, SNP_HBB_SEQ, SNP_HBB_DIS, SNP_HBB_ID, TE_SEV_HIGH);
    assert(t->size == 1); PASS();

    TEST("insert second SNP: size = 2");
    trie_insert(t, SNP_BRCA_SEQ, SNP_BRCA_DIS, SNP_BRCA_ID, TE_SEV_HIGH);
    assert(t->size == 2); PASS();

    TEST("insert duplicate sequence: size stays 2");
    trie_insert(t, SNP_HBB_SEQ, "Other Disease", "rs999", TE_SEV_LOW);
    assert(t->size == 2); PASS();

    TEST("insert NULL sequence: no crash");
    trie_insert(t, NULL, SNP_HBB_DIS, SNP_HBB_ID, 0); PASS();

    TEST("insert invalid base: not added");
    int sz_before = t->size;
    trie_insert(t, "ATGX", "Test", "rs1", 0);
    assert(t->size == sz_before); PASS();

    trie_free(t);
}

static void test_exact_search(void) {
    printf("\n[exact search]\n");
    Trie* t = make_test_trie();

    TEST("exact search: inserted SNP found");
    assert(trie_search_exact(t, SNP_HBB_SEQ) == TE_OK); PASS();

    TEST("exact search: second SNP found");
    assert(trie_search_exact(t, SNP_BRCA_SEQ) == TE_OK); PASS();

    TEST("exact search: not-inserted sequence → TE_ERR");
    assert(trie_search_exact(t, "AAAAAAAAAAAAAAAAAAA") == TE_ERR); PASS();

    TEST("exact search: NULL query → TE_ERR");
    assert(trie_search_exact(t, NULL) == TE_ERR); PASS();

    TEST("exact search: prefix only → TE_ERR");
    /* First 5 chars of HBB SNP — not a complete SNP */
    assert(trie_search_exact(t, "GTGCA") == TE_ERR); PASS();

    TEST("exact search: empty string → TE_ERR");
    assert(trie_search_exact(t, "") == TE_ERR); PASS();

    trie_free(t);
}

static void test_fuzzy_exact(void) {
    printf("\n[fuzzy search — 0 errors (exact mode)]\n");
    Trie* t = make_test_trie();

    TEST("fuzzy max_errors=0: exact match found");
    TrieMatchList* r = trie_search_fuzzy(t, SNP_HBB_SEQ, 0);
    assert(r && r->count == 1);
    assert(r->matches[0].mismatches == 0);
    assert(strcmp(r->matches[0].snp_id, SNP_HBB_ID) == 0);
    trie_match_list_free(r); PASS();

    TEST("fuzzy max_errors=0: no match → count=0");
    r = trie_search_fuzzy(t, "AAAAAAAAAAAAAAAAAAA", 0);
    assert(r && r->count == 0);
    trie_match_list_free(r); PASS();

    TEST("fuzzy max_errors=0: 1 mismatch → not found");
    /* Mutate position 0: G→A */
    char mutated[64];
    strcpy(mutated, SNP_HBB_SEQ);
    mutated[0] = 'A';
    r = trie_search_fuzzy(t, mutated, 0);
    assert(r && r->count == 0);
    trie_match_list_free(r); PASS();

    trie_free(t);
}

static void test_fuzzy_one_error(void) {
    printf("\n[fuzzy search — 1 error]\n");
    Trie* t = make_test_trie();

    TEST("fuzzy max_errors=1: exact match still found");
    TrieMatchList* r = trie_search_fuzzy(t, SNP_HBB_SEQ, 1);
    assert(r && r->count >= 1);
    /* Find the exact match in results */
    int exact_found = 0;
    for (int i = 0; i < r->count; i++)
        if (r->matches[i].mismatches == 0 &&
            strcmp(r->matches[i].snp_id, SNP_HBB_ID) == 0)
            exact_found = 1;
    assert(exact_found);
    trie_match_list_free(r); PASS();

    TEST("fuzzy max_errors=1: 1 mismatch found with mismatches=1");
    char one_mm[64];
    strcpy(one_mm, SNP_HBB_SEQ);
    one_mm[5] = (one_mm[5] == 'A') ? 'T' : 'A'; /* flip position 5 */
    r = trie_search_fuzzy(t, one_mm, 1);
    int found_mm1 = 0;
    for (int i = 0; r && i < r->count; i++)
        if (r->matches[i].mismatches == 1 &&
            strcmp(r->matches[i].snp_id, SNP_HBB_ID) == 0)
            found_mm1 = 1;
    assert(found_mm1);
    trie_match_list_free(r); PASS();

    TEST("fuzzy max_errors=1: 2 mismatches NOT found");
    char two_mm[64];
    strcpy(two_mm, SNP_HBB_SEQ);
    two_mm[0] = (two_mm[0] == 'A') ? 'T' : 'A';
    two_mm[1] = (two_mm[1] == 'A') ? 'T' : 'A';
    r = trie_search_fuzzy(t, two_mm, 1);
    int found_hbb = 0;
    for (int i = 0; r && i < r->count; i++)
        if (strcmp(r->matches[i].snp_id, SNP_HBB_ID) == 0)
            found_hbb = 1;
    assert(!found_hbb);
    trie_match_list_free(r); PASS();

    trie_free(t);
}

static void test_fuzzy_two_errors(void) {
    printf("\n[fuzzy search — 2 errors]\n");
    Trie* t = make_test_trie();

    TEST("fuzzy max_errors=2: 2 mismatches found");
    char two_mm[64];
    strcpy(two_mm, SNP_HBB_SEQ);
    two_mm[0] = (two_mm[0] == 'A') ? 'T' : 'A';
    two_mm[1] = (two_mm[1] == 'A') ? 'T' : 'A';
    TrieMatchList* r = trie_search_fuzzy(t, two_mm, 2);
    int found = 0;
    for (int i = 0; r && i < r->count; i++)
        if (r->matches[i].mismatches == 2 &&
            strcmp(r->matches[i].snp_id, SNP_HBB_ID) == 0)
            found = 1;
    assert(found);
    trie_match_list_free(r); PASS();

    TEST("fuzzy max_errors=2: 3 mismatches NOT found");
    char three_mm[64];
    strcpy(three_mm, SNP_HBB_SEQ);
    three_mm[0] = (three_mm[0]=='A') ? 'T' : 'A';
    three_mm[1] = (three_mm[1]=='A') ? 'T' : 'A';
    three_mm[2] = (three_mm[2]=='A') ? 'T' : 'A';
    r = trie_search_fuzzy(t, three_mm, 2);
    found = 0;
    for (int i = 0; r && i < r->count; i++)
        if (strcmp(r->matches[i].snp_id, SNP_HBB_ID) == 0)
            found = 1;
    assert(!found);
    trie_match_list_free(r); PASS();

    TEST("fuzzy clamps max_errors > TE_MAX_ERRORS");
    r = trie_search_fuzzy(t, SNP_HBB_SEQ, 99);
    assert(r != NULL); /* should not crash */
    trie_match_list_free(r); PASS();

    trie_free(t);
}

static void test_search_in_sequence(void) {
    printf("\n[search in sequence]\n");
    Trie* t = make_test_trie();

    TEST("scan sequence containing HBB SNP: found");
    /* Build a longer sequence with SNP_HBB_SEQ embedded at position 10 */
    char seq[128] = "ATGCATGCAT";
    strcat(seq, SNP_HBB_SEQ);
    strcat(seq, "GCTAGCTAGC");

    TrieMatchList* r = trie_search_in_sequence(t, seq, 0);
    int found = 0;
    for (int i = 0; r && i < r->count; i++)
        if (strcmp(r->matches[i].snp_id, SNP_HBB_ID) == 0 &&
            r->matches[i].match_position == 10)
            found = 1;
    assert(found);
    trie_match_list_free(r); PASS();

    TEST("scan sequence with no SNP: count=0");
    r = trie_search_in_sequence(t, "AAAAAAAAAAAAAAAAAAAAAAAAAAAA", 0);
    assert(r && r->count == 0);
    trie_match_list_free(r); PASS();

    TEST("scan NULL sequence: returns NULL");
    assert(trie_search_in_sequence(t, NULL, 0) == NULL); PASS();

    trie_free(t);
}

static void test_load_from_file(void) {
    printf("\n[load from file]\n");

    /* Write a temp TSV file */
    const char* tmpfile = "test_snps.tsv";
    FILE* fp = fopen(tmpfile, "w");
    assert(fp != NULL);
    fprintf(fp, "rs334\tGTGCACCTGACTCCTGTG\tSickle Cell Anemia\t2\n");
    fprintf(fp, "rs28897672\tATGGATTTATCTGCTCTT\tBreast Cancer (BRCA1)\t2\n");
    fprintf(fp, "rs11549407\tATGGTGCACCTGACTCCT\tBeta-Thalassemia\t1\n");
    fprintf(fp, "bad_line_only_3_fields\tATGC\tmissing\n"); /* malformed */
    fclose(fp);

    TEST("load_from_file returns correct count (3 valid)");
    Trie* t = trie_create();
    int count = trie_load_from_file(t, tmpfile);
    assert(count == 3);
    PASS();

    TEST("loaded SNPs are searchable");
    assert(trie_search_exact(t, "GTGCACCTGACTCCTGTG") == TE_OK); PASS();

    TEST("load_from_file bad path returns -1");
    assert(trie_load_from_file(t, "/nonexistent/path.tsv") == -1); PASS();

    trie_free(t);
    remove(tmpfile);
}

static void test_match_list_free(void) {
    printf("\n[match list free]\n");

    TEST("trie_match_list_free(NULL) is safe");
    trie_match_list_free(NULL); PASS();

    TEST("free result from empty search");
    Trie* t = make_test_trie();
    TrieMatchList* r = trie_search_fuzzy(t, "AAAAAAAAAAAAAAAAAAA", 0);
    trie_match_list_free(r); /* valgrind checks no leak */
    trie_free(t); PASS();
}

static void test_json(void) {
    printf("\n[JSON serialisation]\n");
    Trie* t = make_test_trie();

    TEST("JSON output non-NULL");
    TrieMatchList* r = trie_search_fuzzy(t, SNP_HBB_SEQ, 0);
    char* json = trie_match_list_to_json(r);
    assert(json != NULL); PASS();

    TEST("JSON contains count key");
    assert(strstr(json, "\"count\"") != NULL); PASS();

    TEST("JSON contains matches key");
    assert(strstr(json, "\"matches\"") != NULL); PASS();

    TEST("JSON contains snp_id value");
    assert(strstr(json, SNP_HBB_ID) != NULL); PASS();

    TEST("JSON contains disease value");
    assert(strstr(json, "Sickle") != NULL); PASS();

    free(json);
    trie_match_list_free(r);

    TEST("JSON for empty list");
    r = trie_search_fuzzy(t, "AAAAAAAAAAAAAAAAAAA", 0);
    json = trie_match_list_to_json(r);
    assert(json != NULL);
    assert(strstr(json, "\"count\":0") != NULL);
    free(json);
    trie_match_list_free(r); PASS();

    TEST("trie_match_list_to_json(NULL) returns NULL");
    assert(trie_match_list_to_json(NULL) == NULL); PASS();

    trie_free(t);
}

static void test_severity(void) {
    printf("\n[severity levels]\n");
    Trie* t = trie_create();

    trie_insert(t, "ATGCATGCATGCATGCAT", "Low Risk",  "rs001", TE_SEV_LOW);
    trie_insert(t, "GCTAGCTAGCTAGCTAGC", "Med Risk",  "rs002", TE_SEV_MEDIUM);
    trie_insert(t, "TACGTACGTACGTACGTA", "High Risk", "rs003", TE_SEV_HIGH);

    TEST("severity LOW stored correctly");
    TrieMatchList* r = trie_search_fuzzy(t, "ATGCATGCATGCATGCAT", 0);
    assert(r && r->count == 1 && r->matches[0].severity == TE_SEV_LOW);
    trie_match_list_free(r); PASS();

    TEST("severity HIGH stored correctly");
    r = trie_search_fuzzy(t, "TACGTACGTACGTACGTA", 0);
    assert(r && r->count == 1 && r->matches[0].severity == TE_SEV_HIGH);
    trie_match_list_free(r); PASS();

    trie_free(t);
}

int main(void) {
    printf("=======================================================\n");
    printf("  GenomeX -- trie_errors.c self-test\n");
    printf("=======================================================\n");

    test_create_free();
    test_insert();
    test_exact_search();
    test_fuzzy_exact();
    test_fuzzy_one_error();
    test_fuzzy_two_errors();
    test_search_in_sequence();
    test_load_from_file();
    test_match_list_free();
    test_json();
    test_severity();

    printf("\n=======================================================\n");
    printf("  Results: %d / %d tests passed\n", tests_passed, tests_run);
    printf("=======================================================\n");

    return (tests_passed == tests_run) ? 0 : 1;
}

#endif /* TE_TEST */