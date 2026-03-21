/*
 * dna_utils.c
 * -----------
 * Implementation of shared DNA utility functions for GenomeX.
 *
 * Build:
 *   gcc -O2 -Wall -std=c11 -c dna_utils.c -o dna_utils.o
 *
 * Test:
 *   gcc -O2 -Wall -std=c11 -DDNA_DEBUG -DDNA_UTILS_TEST \
 *       dna_utils.c -o test_utils && ./test_utils
 */

#include "dna_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

/* ============================================================
   THREAD-LOCAL LAST ERROR
   ============================================================ */

/* _Thread_local requires C11 — falls back to static in single-threaded use */
#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
    static _Thread_local DnaError g_last_error = DNA_SUCCESS;
#else
    static DnaError g_last_error = DNA_SUCCESS;
#endif

/* Internal setter — only used within this file */
static void set_error(DnaError err) {
    g_last_error = err;
}

DnaError dna_last_error(void) {
    return g_last_error;
}

/* ============================================================
   LOOKUP TABLES (built once at startup via initializer)
   ============================================================ */

/*
 * BASE_TO_IDX_TABLE[256]
 * ----------------------
 * Maps every ASCII character to its base index.
 * Valid: A/a=0, T/t=1, G/g=2, C/c=3, $=4
 * All others: -1
 *
 * Using a 256-entry table makes base_to_index() a single array
 * lookup — O(1) with no branching. Critical for hot paths.
 */
static int BASE_TO_IDX_TABLE[256];
static char IDX_TO_BASE_TABLE[5] = { 'A', 'T', 'G', 'C', '$' };

/* COMPLEMENT_TABLE[256]: maps each base char to its complement char */
static char COMPLEMENT_TABLE[256];

/* Called once before any function runs */
static void init_tables(void) {
    /* Fill base→index table with -1 (invalid) */
    for (int i = 0; i < 256; i++) BASE_TO_IDX_TABLE[i] = -1;

    /* Valid uppercase bases */
    BASE_TO_IDX_TABLE[(unsigned char)'A'] = DNA_BASE_A;
    BASE_TO_IDX_TABLE[(unsigned char)'T'] = DNA_BASE_T;
    BASE_TO_IDX_TABLE[(unsigned char)'G'] = DNA_BASE_G;
    BASE_TO_IDX_TABLE[(unsigned char)'C'] = DNA_BASE_C;
    BASE_TO_IDX_TABLE[(unsigned char)'$'] = DNA_BASE_END;

    /* Accept lowercase too (internally we always uppercase) */
    BASE_TO_IDX_TABLE[(unsigned char)'a'] = DNA_BASE_A;
    BASE_TO_IDX_TABLE[(unsigned char)'t'] = DNA_BASE_T;
    BASE_TO_IDX_TABLE[(unsigned char)'g'] = DNA_BASE_G;
    BASE_TO_IDX_TABLE[(unsigned char)'c'] = DNA_BASE_C;

    /* Complement table */
    for (int i = 0; i < 256; i++) COMPLEMENT_TABLE[i] = '\0';
    COMPLEMENT_TABLE[(unsigned char)'A'] = 'T';
    COMPLEMENT_TABLE[(unsigned char)'T'] = 'A';
    COMPLEMENT_TABLE[(unsigned char)'G'] = 'C';
    COMPLEMENT_TABLE[(unsigned char)'C'] = 'G';
    COMPLEMENT_TABLE[(unsigned char)'a'] = 'T';
    COMPLEMENT_TABLE[(unsigned char)'t'] = 'A';
    COMPLEMENT_TABLE[(unsigned char)'g'] = 'C';
    COMPLEMENT_TABLE[(unsigned char)'c'] = 'G';
}

/*
 * GCC/Clang constructor: runs init_tables() before main().
 * MSVC alternative uses #pragma init_seg.
 */
#if defined(__GNUC__) || defined(__clang__)
__attribute__((constructor))
static void dna_utils_init(void) {
    init_tables();
}
#else
/* Fallback: caller must invoke dna_utils_init() explicitly */
void dna_utils_init(void) {
    init_tables();
}
#endif

/* ============================================================
   BASE ENCODING
   ============================================================ */

int base_to_index(char base) {
    return BASE_TO_IDX_TABLE[(unsigned char)base];
}

char index_to_base(int index) {
    if (index < 0 || index > 4) return '\0';
    return IDX_TO_BASE_TABLE[index];
}

char dna_complement_base(char base) {
    return COMPLEMENT_TABLE[(unsigned char)base];
}

/* ============================================================
   VALIDATION
   ============================================================ */

int dna_validate(const char* seq) {
    if (seq == NULL) {
        set_error(DNA_ERR_NULL_INPUT);
        return DNA_ERR;
    }
    if (*seq == '\0') {
        set_error(DNA_ERR_EMPTY_SEQ);
        return DNA_ERR;
    }

    size_t len = 0;
    for (const char* p = seq; *p != '\0'; p++, len++) {
        /* Allow only A, T, G, C (upper or lower) */
        unsigned char uc = (unsigned char)*p;
        if (BASE_TO_IDX_TABLE[uc] == -1) {
            set_error(DNA_ERR_INVALID_BASE);
            return DNA_ERR;
        }
        if (len > DNA_MAX_SEQ_LEN) {
            set_error(DNA_ERR_TOO_LONG);
            return DNA_ERR;
        }
    }

    set_error(DNA_SUCCESS);
    return DNA_OK;
}

int dna_validate_n(const char* seq, size_t n) {
    if (seq == NULL) { set_error(DNA_ERR_NULL_INPUT); return DNA_ERR; }
    if (n == 0)      { set_error(DNA_ERR_EMPTY_SEQ);  return DNA_ERR; }
    if (n > DNA_MAX_SEQ_LEN) { set_error(DNA_ERR_TOO_LONG); return DNA_ERR; }

    for (size_t i = 0; i < n; i++) {
        if (BASE_TO_IDX_TABLE[(unsigned char)seq[i]] == -1) {
            /* $ is valid internally but not from user input */
            if (seq[i] != '$') {
                set_error(DNA_ERR_INVALID_BASE);
                return DNA_ERR;
            }
        }
    }
    set_error(DNA_SUCCESS);
    return DNA_OK;
}

/* ============================================================
   STRING HELPERS
   ============================================================ */

char* dna_uppercase(char* seq) {
    if (seq == NULL) return NULL;
    for (char* p = seq; *p != '\0'; p++) {
        *p = (char)toupper((unsigned char)*p);
    }
    return seq;
}

/* ============================================================
   SEQUENCE MANAGEMENT
   ============================================================ */

DnaSequence* dna_sequence_create(const char* seq) {
    if (seq == NULL) { set_error(DNA_ERR_NULL_INPUT); return NULL; }

    size_t len = strlen(seq);
    if (len == 0)              { set_error(DNA_ERR_EMPTY_SEQ); return NULL; }
    if (len > DNA_MAX_SEQ_LEN) { set_error(DNA_ERR_TOO_LONG);  return NULL; }

    DnaSequence* ds = (DnaSequence*)malloc(sizeof(DnaSequence));
    if (ds == NULL) { set_error(DNA_ERR_ALLOC); return NULL; }

    ds->data = (char*)malloc(len + 1);
    if (ds->data == NULL) {
        free(ds);
        set_error(DNA_ERR_ALLOC);
        return NULL;
    }

    memcpy(ds->data, seq, len + 1);
    dna_uppercase(ds->data);

    /* Validate after uppercasing */
    if (!dna_validate(ds->data)) {
        free(ds->data);
        free(ds);
        return NULL;   /* error already set by dna_validate */
    }

    ds->length = len;
    ds->header[0] = '\0';

    set_error(DNA_SUCCESS);
    return ds;
}

void dna_sequence_free(DnaSequence* ds) {
    if (ds == NULL) return;
    if (ds->data != NULL) {
        free(ds->data);
        ds->data = NULL;
    }
    free(ds);
}

DnaSequence* dna_sequence_copy(const DnaSequence* ds) {
    if (ds == NULL) { set_error(DNA_ERR_NULL_INPUT); return NULL; }

    DnaSequence* copy = (DnaSequence*)malloc(sizeof(DnaSequence));
    if (copy == NULL) { set_error(DNA_ERR_ALLOC); return NULL; }

    copy->data = (char*)malloc(ds->length + 1);
    if (copy->data == NULL) {
        free(copy);
        set_error(DNA_ERR_ALLOC);
        return NULL;
    }

    memcpy(copy->data, ds->data, ds->length + 1);
    copy->length = ds->length;
    strncpy(copy->header, ds->header, DNA_MAX_HEADER - 1);
    copy->header[DNA_MAX_HEADER - 1] = '\0';

    set_error(DNA_SUCCESS);
    return copy;
}

/* ============================================================
   FASTA PARSING
   ============================================================ */

/*
 * Internal helper: strips whitespace and newlines from a DNA line.
 * Writes result into buf. Returns number of characters written.
 */
static size_t strip_dna_line(const char* line, char* buf, size_t buf_size, size_t written) {
    for (const char* p = line; *p != '\0' && *p != '\n' && *p != '\r'; p++) {
        char uc = (char)toupper((unsigned char)*p);
        /* Skip whitespace silently */
        if (uc == ' ' || uc == '\t') continue;
        /* Only copy valid DNA bases */
        if (BASE_TO_IDX_TABLE[(unsigned char)uc] != -1 && uc != '$') {
            if (written < buf_size - 1) {
                buf[written++] = uc;
            }
        }
    }
    return written;
}

DnaSequence* fasta_parse_file(const char* filepath, DnaError* err_out) {
    if (filepath == NULL) {
        if (err_out) *err_out = DNA_ERR_NULL_INPUT;
        set_error(DNA_ERR_NULL_INPUT);
        return NULL;
    }

    FILE* fp = fopen(filepath, "r");
    if (fp == NULL) {
        if (err_out) *err_out = DNA_ERR_FILE;
        set_error(DNA_ERR_FILE);
        return NULL;
    }

    DnaSequence* ds = (DnaSequence*)calloc(1, sizeof(DnaSequence));
    if (ds == NULL) {
        fclose(fp);
        if (err_out) *err_out = DNA_ERR_ALLOC;
        set_error(DNA_ERR_ALLOC);
        return NULL;
    }
    ds->header[0] = '\0';

    /* Allocate a large working buffer (grows if needed) */
    size_t buf_capacity = 65536;   /* Start with 64KB */
    char* buf = (char*)malloc(buf_capacity);
    if (buf == NULL) {
        fclose(fp);
        free(ds);
        if (err_out) *err_out = DNA_ERR_ALLOC;
        set_error(DNA_ERR_ALLOC);
        return NULL;
    }

    char line[4096];
    int header_found = 0;
    size_t written = 0;

    while (fgets(line, sizeof(line), fp) != NULL) {
        if (line[0] == '>') {
            if (header_found) {
                /* Second sequence starts — stop (we only read the first) */
                break;
            }
            header_found = 1;
            /* Copy header (strip '>' and newline) */
            size_t hlen = strlen(line + 1);
            if (hlen > 0 && (line[hlen] == '\n' || line[hlen] == '\r')) hlen--;
            size_t copy_len = hlen < (DNA_MAX_HEADER - 1) ? hlen : (DNA_MAX_HEADER - 1);
            strncpy(ds->header, line + 1, copy_len);
            ds->header[copy_len] = '\0';
            continue;
        }

        if (line[0] == ';') continue;  /* FASTA comment line */
        if (!header_found)  continue;  /* Skip lines before first '>' */

        /* Strip and accumulate DNA from this line */
        /* Grow buffer if needed */
        if (written + sizeof(line) + 2 > buf_capacity) {
            buf_capacity *= 2;
            if (buf_capacity > DNA_MAX_SEQ_LEN + 2) buf_capacity = DNA_MAX_SEQ_LEN + 2;
            char* new_buf = (char*)realloc(buf, buf_capacity);
            if (new_buf == NULL) {
                fclose(fp);
                free(buf);
                free(ds);
                if (err_out) *err_out = DNA_ERR_ALLOC;
                set_error(DNA_ERR_ALLOC);
                return NULL;
            }
            buf = new_buf;
        }

        written = strip_dna_line(line, buf, buf_capacity, written);
    }

    fclose(fp);

    if (!header_found || written == 0) {
        free(buf);
        free(ds);
        if (err_out) *err_out = DNA_ERR_BAD_FASTA;
        set_error(DNA_ERR_BAD_FASTA);
        return NULL;
    }

    buf[written] = '\0';

    /* Transfer buffer ownership to DnaSequence */
    ds->data   = buf;
    ds->length = written;

    if (err_out) *err_out = DNA_SUCCESS;
    set_error(DNA_SUCCESS);
    return ds;
}

DnaSequence* fasta_parse_string(const char* fasta_text, DnaError* err_out) {
    if (fasta_text == NULL) {
        if (err_out) *err_out = DNA_ERR_NULL_INPUT;
        set_error(DNA_ERR_NULL_INPUT);
        return NULL;
    }

    DnaSequence* ds = (DnaSequence*)calloc(1, sizeof(DnaSequence));
    if (!ds) { if (err_out) *err_out = DNA_ERR_ALLOC; return NULL; }

    size_t text_len = strlen(fasta_text);
    size_t buf_capacity = text_len + 2;
    char*  buf = (char*)malloc(buf_capacity);
    if (!buf) {
        free(ds);
        if (err_out) *err_out = DNA_ERR_ALLOC;
        return NULL;
    }

    const char* p    = fasta_text;
    int header_found = 0;
    size_t written   = 0;

    while (*p != '\0') {
        /* Find end of current line */
        const char* line_end = p;
        while (*line_end != '\0' && *line_end != '\n') line_end++;

        size_t line_len = (size_t)(line_end - p);

        if (line_len > 0 && p[0] == '>') {
            if (header_found) break;
            header_found = 1;
            size_t hlen = line_len - 1;  /* skip '>' */
            if (hlen >= DNA_MAX_HEADER) hlen = DNA_MAX_HEADER - 1;
            strncpy(ds->header, p + 1, hlen);
            ds->header[hlen] = '\0';
        } else if (header_found && line_len > 0 && p[0] != ';') {
            /* Temporary line buffer */
            char linebuf[4096];
            size_t copy = line_len < sizeof(linebuf)-1 ? line_len : sizeof(linebuf)-1;
            memcpy(linebuf, p, copy);
            linebuf[copy] = '\0';
            written = strip_dna_line(linebuf, buf, buf_capacity, written);
        }

        p = line_end;
        if (*p == '\n') p++;
    }

    if (!header_found || written == 0) {
        /* Try treating entire string as raw sequence */
        free(buf);
        free(ds);
        return fasta_parse_raw(fasta_text, err_out);
    }

    buf[written] = '\0';
    ds->data   = buf;
    ds->length = written;

    if (err_out) *err_out = DNA_SUCCESS;
    set_error(DNA_SUCCESS);
    return ds;
}

DnaSequence* fasta_parse_raw(const char* raw_seq, DnaError* err_out) {
    if (raw_seq == NULL) {
        if (err_out) *err_out = DNA_ERR_NULL_INPUT;
        set_error(DNA_ERR_NULL_INPUT);
        return NULL;
    }

    size_t raw_len = strlen(raw_seq);
    if (raw_len == 0) {
        if (err_out) *err_out = DNA_ERR_EMPTY_SEQ;
        set_error(DNA_ERR_EMPTY_SEQ);
        return NULL;
    }

    char* buf = (char*)malloc(raw_len + 1);
    if (!buf) {
        if (err_out) *err_out = DNA_ERR_ALLOC;
        set_error(DNA_ERR_ALLOC);
        return NULL;
    }

    size_t written = 0;
    for (size_t i = 0; i < raw_len; i++) {
        char uc = (char)toupper((unsigned char)raw_seq[i]);
        if (uc == ' ' || uc == '\t' || uc == '\n' || uc == '\r') continue;
        if (BASE_TO_IDX_TABLE[(unsigned char)uc] != -1 && uc != '$') {
            buf[written++] = uc;
        }
        /* Silently skip invalid characters to be lenient with user input */
    }

    if (written == 0) {
        free(buf);
        if (err_out) *err_out = DNA_ERR_INVALID_BASE;
        set_error(DNA_ERR_INVALID_BASE);
        return NULL;
    }

    buf[written] = '\0';

    DnaSequence* ds = (DnaSequence*)calloc(1, sizeof(DnaSequence));
    if (!ds) {
        free(buf);
        if (err_out) *err_out = DNA_ERR_ALLOC;
        set_error(DNA_ERR_ALLOC);
        return NULL;
    }

    ds->data   = buf;
    ds->length = written;
    ds->header[0] = '\0';

    if (err_out) *err_out = DNA_SUCCESS;
    set_error(DNA_SUCCESS);
    return ds;
}

/* ============================================================
   STRING OPERATIONS
   ============================================================ */

char* dna_reverse_complement(const char* seq) {
    if (seq == NULL) { set_error(DNA_ERR_NULL_INPUT); return NULL; }

    size_t len = strlen(seq);
    if (len == 0) { set_error(DNA_ERR_EMPTY_SEQ); return NULL; }

    char* result = (char*)malloc(len + 1);
    if (result == NULL) { set_error(DNA_ERR_ALLOC); return NULL; }

    for (size_t i = 0; i < len; i++) {
        char base = (char)toupper((unsigned char)seq[len - 1 - i]);
        char comp = COMPLEMENT_TABLE[(unsigned char)base];
        if (comp == '\0') {
            /* Invalid base in input */
            free(result);
            set_error(DNA_ERR_INVALID_BASE);
            return NULL;
        }
        result[i] = comp;
    }
    result[len] = '\0';

    set_error(DNA_SUCCESS);
    return result;
}

char* dna_concat_with_separator(const char* seq1, const char* seq2, char sep) {
    if (seq1 == NULL || seq2 == NULL) {
        set_error(DNA_ERR_NULL_INPUT);
        return NULL;
    }

    size_t len1   = strlen(seq1);
    size_t len2   = strlen(seq2);
    /* Layout: seq1 + sep + seq2 + $ + \0 */
    size_t total  = len1 + 1 + len2 + 1 + 1;

    char* result = (char*)malloc(total);
    if (result == NULL) { set_error(DNA_ERR_ALLOC); return NULL; }

    memcpy(result, seq1, len1);
    result[len1] = sep;
    memcpy(result + len1 + 1, seq2, len2);
    result[len1 + 1 + len2]     = '$';
    result[len1 + 1 + len2 + 1] = '\0';

    set_error(DNA_SUCCESS);
    return result;
}

char* dna_subsequence(const char* seq, size_t start, size_t len) {
    if (seq == NULL) { set_error(DNA_ERR_NULL_INPUT); return NULL; }

    size_t seq_len = strlen(seq);
    if (start >= seq_len || start + len > seq_len) {
        set_error(DNA_ERR_INVALID_BASE);   /* re-using as "out of bounds" */
        return NULL;
    }

    char* sub = (char*)malloc(len + 1);
    if (sub == NULL) { set_error(DNA_ERR_ALLOC); return NULL; }

    memcpy(sub, seq + start, len);
    sub[len] = '\0';

    set_error(DNA_SUCCESS);
    return sub;
}

/* ============================================================
   SIMILARITY METRICS
   ============================================================ */

double dna_similarity_simple(const char* a, const char* b) {
    if (a == NULL || b == NULL) { set_error(DNA_ERR_NULL_INPUT); return 0.0; }

    size_t len_a = strlen(a);
    size_t len_b = strlen(b);
    size_t min_len = len_a < len_b ? len_a : len_b;

    if (min_len == 0) { set_error(DNA_ERR_EMPTY_SEQ); return 0.0; }

    size_t matches = 0;
    for (size_t i = 0; i < min_len; i++) {
        char ca = (char)toupper((unsigned char)a[i]);
        char cb = (char)toupper((unsigned char)b[i]);
        if (ca == cb) matches++;
    }

    set_error(DNA_SUCCESS);
    return ((double)matches / (double)min_len) * 100.0;
}

int dna_hamming_distance(const char* a, const char* b) {
    if (a == NULL || b == NULL) { set_error(DNA_ERR_NULL_INPUT); return -1; }

    size_t len_a = strlen(a);
    size_t len_b = strlen(b);

    if (len_a != len_b) {
        set_error(DNA_ERR_INVALID_BASE);
        return -1;   /* Hamming undefined for different lengths */
    }

    int dist = 0;
    for (size_t i = 0; i < len_a; i++) {
        char ca = (char)toupper((unsigned char)a[i]);
        char cb = (char)toupper((unsigned char)b[i]);
        if (ca != cb) dist++;
    }

    set_error(DNA_SUCCESS);
    return dist;
}

double dna_gc_content(const char* seq) {
    if (seq == NULL) { set_error(DNA_ERR_NULL_INPUT); return 0.0; }

    size_t total = 0, gc = 0;
    for (const char* p = seq; *p != '\0'; p++) {
        char uc = (char)toupper((unsigned char)*p);
        if (uc == 'A' || uc == 'T' || uc == 'G' || uc == 'C') {
            total++;
            if (uc == 'G' || uc == 'C') gc++;
        }
    }

    if (total == 0) { set_error(DNA_ERR_EMPTY_SEQ); return 0.0; }

    set_error(DNA_SUCCESS);
    return ((double)gc / (double)total) * 100.0;
}

/* ============================================================
   K-MER UTILITIES
   ============================================================ */

int dna_kmer_index(const char* kmer, int k) {
    if (kmer == NULL || k <= 0 || k > 10) return -1;

    int index = 0;
    for (int i = 0; i < k; i++) {
        int b = base_to_index(kmer[i]);
        if (b < 0 || b > 3) return -1;  /* invalid base ($ not allowed) */
        index = index * 4 + b;
    }
    return index;
}

int* dna_count_kmers(const char* seq, int k) {
    if (seq == NULL || k <= 0 || k > 10) {
        set_error(k > 10 ? DNA_ERR_TOO_LONG : DNA_ERR_NULL_INPUT);
        return NULL;
    }

    /* Number of possible k-mers = 4^k */
    int table_size = 1;
    for (int i = 0; i < k; i++) {
        table_size *= 4;
        /* Overflow guard */
        if (table_size > 1048576) { /* 4^10 = 1,048,576 */
            set_error(DNA_ERR_TOO_LONG);
            return NULL;
        }
    }

    int* counts = (int*)calloc((size_t)table_size, sizeof(int));
    if (counts == NULL) { set_error(DNA_ERR_ALLOC); return NULL; }

    size_t seq_len = strlen(seq);
    if ((size_t)k > seq_len) {
        /* Sequence shorter than k — return all-zero counts, not an error */
        set_error(DNA_SUCCESS);
        return counts;
    }

    for (size_t i = 0; i <= seq_len - (size_t)k; i++) {
        int idx = dna_kmer_index(seq + i, k);
        if (idx >= 0) {
            counts[idx]++;
        }
        /* Skip k-mers containing invalid bases (e.g. N) silently */
    }

    set_error(DNA_SUCCESS);
    return counts;
}

/* ============================================================
   ERROR REPORTING
   ============================================================ */

const char* dna_error_string(DnaError err) {
    switch (err) {
        case DNA_SUCCESS:          return "Success";
        case DNA_ERR_NULL_INPUT:   return "NULL input pointer";
        case DNA_ERR_INVALID_BASE: return "Invalid base character (expected A, T, G, C)";
        case DNA_ERR_EMPTY_SEQ:    return "Empty sequence";
        case DNA_ERR_TOO_LONG:     return "Sequence exceeds maximum allowed length";
        case DNA_ERR_ALLOC:        return "Memory allocation failed";
        case DNA_ERR_FILE:         return "File open/read failed";
        case DNA_ERR_BAD_FASTA:    return "Malformed FASTA format (missing header or sequence)";
        default:                   return "Unknown error";
    }
}

/* ============================================================
   DEBUG HELPERS
   ============================================================ */

#ifdef DNA_DEBUG
void dna_print_sequence(const DnaSequence* ds) {
    if (ds == NULL) { printf("[DnaSequence: NULL]\n"); return; }
    printf("[DnaSequence]\n");
    printf("  header : %s\n", ds->header[0] ? ds->header : "(none)");
    printf("  length : %zu\n", ds->length);
    printf("  gc     : %.1f%%\n", dna_gc_content(ds->data));
    printf("  data   : %.80s%s\n", ds->data, ds->length > 80 ? "..." : "");
}

void dna_print_bases(const char* seq, size_t max_len) {
    if (seq == NULL) { printf("(NULL)\n"); return; }
    size_t len = strlen(seq);
    size_t show = len < max_len ? len : max_len;
    for (size_t i = 0; i < show; i++) {
        char c = seq[i];
        /* Color-code bases in terminal */
        switch (c) {
            case 'A': printf("\033[32mA\033[0m"); break;  /* green */
            case 'T': printf("\033[36mT\033[0m"); break;  /* cyan  */
            case 'G': printf("\033[35mG\033[0m"); break;  /* magenta */
            case 'C': printf("\033[33mC\033[0m"); break;  /* yellow */
            default:  printf("%c", c);
        }
    }
    if (len > max_len) printf("...[%zu more]", len - max_len);
    printf("\n");
}
#endif  /* DNA_DEBUG */


/* ============================================================
   SELF-TEST (compile with -DDNA_UTILS_TEST to run)
   ============================================================ */

#ifdef DNA_UTILS_TEST

#include <assert.h>

static int tests_run    = 0;
static int tests_passed = 0;

#define TEST(name) do { \
    tests_run++; \
    printf("  %-50s ", name); \
    fflush(stdout); \
} while(0)

#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(msg) do { printf("FAIL — %s\n", msg); } while(0)

static void test_base_encoding(void) {
    printf("\n[base encoding]\n");

    TEST("base_to_index A=0");
    assert(base_to_index('A') == 0); PASS();

    TEST("base_to_index T=1");
    assert(base_to_index('T') == 1); PASS();

    TEST("base_to_index G=2");
    assert(base_to_index('G') == 2); PASS();

    TEST("base_to_index C=3");
    assert(base_to_index('C') == 3); PASS();

    TEST("base_to_index $=4");
    assert(base_to_index('$') == 4); PASS();

    TEST("base_to_index lowercase a=0");
    assert(base_to_index('a') == 0); PASS();

    TEST("base_to_index invalid X=-1");
    assert(base_to_index('X') == -1); PASS();

    TEST("index_to_base 0='A'");
    assert(index_to_base(0) == 'A'); PASS();

    TEST("index_to_base 4='$'");
    assert(index_to_base(4) == '$'); PASS();

    TEST("index_to_base 5='\\0' (invalid)");
    assert(index_to_base(5) == '\0'); PASS();
}

static void test_validation(void) {
    printf("\n[validation]\n");

    TEST("validate valid uppercase");
    assert(dna_validate("ATGC") == DNA_OK); PASS();

    TEST("validate valid lowercase");
    assert(dna_validate("atgc") == DNA_OK); PASS();

    TEST("validate invalid character X");
    assert(dna_validate("ATGX") == DNA_ERR); PASS();

    TEST("validate NULL");
    assert(dna_validate(NULL) == DNA_ERR); PASS();

    TEST("validate empty string");
    assert(dna_validate("") == DNA_ERR); PASS();

    TEST("validate long valid sequence");
    char* long_seq = (char*)malloc(1001);
    memset(long_seq, 'A', 1000);
    long_seq[1000] = '\0';
    assert(dna_validate(long_seq) == DNA_OK);
    free(long_seq); PASS();
}

static void test_sequence_create(void) {
    printf("\n[DnaSequence create/free]\n");

    TEST("create valid sequence");
    DnaSequence* ds = dna_sequence_create("ATGC");
    assert(ds != NULL);
    assert(ds->length == 4);
    assert(strcmp(ds->data, "ATGC") == 0);
    dna_sequence_free(ds); PASS();

    TEST("create uppercases input");
    ds = dna_sequence_create("atgc");
    assert(ds != NULL);
    assert(strcmp(ds->data, "ATGC") == 0);
    dna_sequence_free(ds); PASS();

    TEST("create NULL returns NULL");
    assert(dna_sequence_create(NULL) == NULL); PASS();

    TEST("create invalid base returns NULL");
    assert(dna_sequence_create("ATGX") == NULL); PASS();

    TEST("free NULL is safe");
    dna_sequence_free(NULL); PASS();   /* should not crash */

    TEST("copy sequence");
    ds = dna_sequence_create("GCTA");
    DnaSequence* copy = dna_sequence_copy(ds);
    assert(copy != NULL);
    assert(strcmp(copy->data, ds->data) == 0);
    assert(copy->data != ds->data);   /* deep copy */
    dna_sequence_free(ds);
    dna_sequence_free(copy); PASS();
}

static void test_fasta_parsing(void) {
    printf("\n[FASTA parsing]\n");

    TEST("parse FASTA string");
    const char* fasta = ">test_header\nATGCATGC\nATGC\n";
    DnaError err;
    DnaSequence* ds = fasta_parse_string(fasta, &err);
    assert(ds != NULL);
    assert(err == DNA_SUCCESS);
    assert(strcmp(ds->data, "ATGCATGCATGC") == 0);
    assert(ds->length == 12);
    dna_sequence_free(ds); PASS();

    TEST("parse FASTA header stored correctly");
    ds = fasta_parse_string(fasta, &err);
    assert(strncmp(ds->header, "test_header", 11) == 0);
    dna_sequence_free(ds); PASS();

    TEST("parse raw sequence (no header)");
    ds = fasta_parse_raw("ATGCATGC", &err);
    assert(ds != NULL);
    assert(ds->length == 8);
    dna_sequence_free(ds); PASS();

    TEST("parse raw with spaces/newlines stripped");
    ds = fasta_parse_raw("ATG CAT\nGC", &err);
    assert(ds != NULL);
    assert(strcmp(ds->data, "ATGCATGC") == 0);
    dna_sequence_free(ds); PASS();

    TEST("parse NULL returns NULL");
    assert(fasta_parse_string(NULL, &err) == NULL);
    assert(err == DNA_ERR_NULL_INPUT); PASS();
}

static void test_string_operations(void) {
    printf("\n[string operations]\n");

    TEST("reverse complement ATGC -> GCAT");
    char* rc = dna_reverse_complement("ATGC");
    assert(rc != NULL);
    assert(strcmp(rc, "GCAT") == 0);
    free(rc); PASS();

    TEST("reverse complement single base A -> T");
    rc = dna_reverse_complement("A");
    assert(strcmp(rc, "T") == 0);
    free(rc); PASS();

    TEST("reverse complement palindrome ATAT -> ATAT");
    rc = dna_reverse_complement("ATAT");
    assert(strcmp(rc, "ATAT") == 0);
    free(rc); PASS();

    TEST("reverse complement NULL returns NULL");
    assert(dna_reverse_complement(NULL) == NULL); PASS();

    TEST("concat with separator");
    char* concat = dna_concat_with_separator("ATGC", "TACG", '#');
    assert(concat != NULL);
    assert(strcmp(concat, "ATGC#TACG$") == 0);
    free(concat); PASS();

    TEST("subsequence middle");
    char* sub = dna_subsequence("ATGCATGC", 2, 4);
    assert(sub != NULL);
    assert(strcmp(sub, "GCAT") == 0);
    free(sub); PASS();

    TEST("subsequence out of bounds returns NULL");
    assert(dna_subsequence("ATGC", 3, 5) == NULL); PASS();
}

static void test_metrics(void) {
    printf("\n[similarity metrics]\n");

    TEST("simple similarity identical = 100%");
    double sim = dna_similarity_simple("ATGC", "ATGC");
    assert(sim == 100.0); PASS();

    TEST("simple similarity zero overlap");
    sim = dna_similarity_simple("AAAA", "TTTT");
    assert(sim == 0.0); PASS();

    TEST("simple similarity 50%");
    sim = dna_similarity_simple("AATT", "AAGC");
    assert(sim == 50.0); PASS();

    TEST("hamming distance 0");
    assert(dna_hamming_distance("ATGC", "ATGC") == 0); PASS();

    TEST("hamming distance 1");
    assert(dna_hamming_distance("ATGC", "TTGC") == 1); PASS();

    TEST("hamming distance different lengths = -1");
    assert(dna_hamming_distance("ATGC", "ATG") == -1); PASS();

    TEST("GC content 50% (ATGC)");
    double gc = dna_gc_content("ATGC");
    assert(gc == 50.0); PASS();

    TEST("GC content 0% (AAAA)");
    gc = dna_gc_content("AAAA");
    assert(gc == 0.0); PASS();

    TEST("GC content 100% (GCGC)");
    gc = dna_gc_content("GCGC");
    assert(gc == 100.0); PASS();
}

static void test_kmers(void) {
    printf("\n[k-mer utilities]\n");

    TEST("kmer index AA = 0");
    assert(dna_kmer_index("AA", 2) == 0); PASS();  /* A=0, 0*4+0=0 */

    TEST("kmer index AT = 1");
    assert(dna_kmer_index("AT", 2) == 1); PASS();  /* A=0, T=1, 0*4+1=1 */

    TEST("kmer index TC = 7");
    assert(dna_kmer_index("TC", 2) == 7); PASS();  /* T=1, C=3, 1*4+3=7 */

    TEST("count kmers k=1 AAATGC");
    int* counts = dna_count_kmers("AAATGC", 1);
    assert(counts != NULL);
    assert(counts[0] == 3);  /* A count = 3 */
    assert(counts[1] == 1);  /* T count = 1 */
    assert(counts[2] == 1);  /* G count = 1 */
    assert(counts[3] == 1);  /* C count = 1 */
    free(counts); PASS();

    TEST("count kmers k=2 ATAT");
    counts = dna_count_kmers("ATAT", 2);
    assert(counts != NULL);
    /* AT appears at pos 0 and 2 → count=2; TA appears at pos 1 → count=1 */
    assert(counts[dna_kmer_index("AT", 2)] == 2);
    assert(counts[dna_kmer_index("TA", 2)] == 1);
    free(counts); PASS();

    TEST("count kmers k>10 returns NULL");
    assert(dna_count_kmers("ATGC", 11) == NULL); PASS();
}

static void test_error_reporting(void) {
    printf("\n[error reporting]\n");

    TEST("error string SUCCESS");
    const char* s = dna_error_string(DNA_SUCCESS);
    assert(s != NULL && strlen(s) > 0); PASS();

    TEST("error string NULL_INPUT");
    s = dna_error_string(DNA_ERR_NULL_INPUT);
    assert(s != NULL && strlen(s) > 0); PASS();

    TEST("last error updates after invalid operation");
    dna_validate(NULL);
    assert(dna_last_error() == DNA_ERR_NULL_INPUT); PASS();

    TEST("last error clears after valid operation");
    dna_validate("ATGC");
    assert(dna_last_error() == DNA_SUCCESS); PASS();
}

int main(void) {
    printf("=======================================================\n");
    printf("  GenomeX — dna_utils.c self-test\n");
    printf("=======================================================\n");

    test_base_encoding();
    test_validation();
    test_sequence_create();
    test_fasta_parsing();
    test_string_operations();
    test_metrics();
    test_kmers();
    test_error_reporting();

    printf("\n=======================================================\n");
    printf("  Results: %d / %d tests passed\n", tests_passed, tests_run);
    printf("=======================================================\n");

    return (tests_passed == tests_run) ? 0 : 1;
}

#endif  /* DNA_UTILS_TEST */