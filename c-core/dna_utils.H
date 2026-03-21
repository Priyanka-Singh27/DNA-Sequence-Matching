#ifndef DNA_UTILS_H
#define DNA_UTILS_H

/*
 * dna_utils.h
 * -----------
 * Shared utility functions for the GenomeX DNA analysis engine.
 * All other modules (#include "dna_utils.h") depend on this.
 *
 * Conventions:
 *   - All functions return error codes where applicable (0=fail, 1=success)
 *   - All heap-allocated return values must be freed by the caller
 *   - DNA strings must contain only: A, T, G, C (uppercase)
 *   - FASTA sequences are stripped of headers and newlines
 */

#include <stddef.h>   /* size_t */

/* ============================================================
   CONSTANTS
   ============================================================ */

#define DNA_ALPHA_SIZE    5     /* A, T, G, C, $ (end marker) */
#define DNA_BASE_A        0
#define DNA_BASE_T        1
#define DNA_BASE_G        2
#define DNA_BASE_C        3
#define DNA_BASE_END      4     /* $ sentinel used by suffix tree */

#define DNA_MAX_SEQ_LEN   100000000   /* 100 MB — hard cap on sequence length */
#define DNA_MAX_HEADER    512          /* max FASTA header line length */

#define DNA_OK            1
#define DNA_ERR           0

/* ============================================================
   ERROR CODES
   ============================================================ */

typedef enum {
    DNA_SUCCESS         =  0,
    DNA_ERR_NULL_INPUT  = -1,   /* NULL pointer passed */
    DNA_ERR_INVALID_BASE= -2,   /* character not in {A,T,G,C} */
    DNA_ERR_EMPTY_SEQ   = -3,   /* zero-length sequence */
    DNA_ERR_TOO_LONG    = -4,   /* exceeds DNA_MAX_SEQ_LEN */
    DNA_ERR_ALLOC       = -5,   /* malloc/calloc returned NULL */
    DNA_ERR_FILE        = -6,   /* file open/read failed */
    DNA_ERR_BAD_FASTA   = -7,   /* malformed FASTA format */
} DnaError;

/* ============================================================
   STRUCTS
   ============================================================ */

/*
 * DnaSequence — a validated, heap-allocated DNA string with metadata.
 * Use dna_sequence_create() / dna_sequence_free() to manage.
 */
typedef struct {
    char*  data;        /* null-terminated DNA string (uppercase) */
    size_t length;      /* length NOT including null terminator */
    char   header[DNA_MAX_HEADER]; /* FASTA header (empty if not from FASTA) */
} DnaSequence;

/* ============================================================
   VALIDATION
   ============================================================ */

/*
 * dna_validate
 * ------------
 * Returns DNA_OK (1) if every character in seq is A, T, G, or C.
 * Returns DNA_ERR (0) otherwise.
 * seq must be null-terminated.
 *
 * Call this at the entry point of every DS operation.
 */
int dna_validate(const char* seq);

/*
 * dna_validate_n
 * --------------
 * Same as dna_validate but checks exactly n characters.
 * Use when seq is NOT null-terminated.
 */
int dna_validate_n(const char* seq, size_t n);

/* ============================================================
   BASE ENCODING
   ============================================================ */

/*
 * base_to_index
 * -------------
 * Maps a DNA base character to its array index.
 *   A → 0
 *   T → 1
 *   G → 2
 *   C → 3
 *   $ → 4  (sentinel — used internally by suffix tree)
 *
 * Returns -1 for any unrecognized character.
 * Hot path — called millions of times. Implemented as a lookup table.
 */
int base_to_index(char base);

/*
 * index_to_base
 * -------------
 * Inverse of base_to_index.
 *   0 → 'A', 1 → 'T', 2 → 'G', 3 → 'C', 4 → '$'
 * Returns '\0' for invalid index.
 */
char index_to_base(int index);

/* ============================================================
   FASTA PARSING
   ============================================================ */

/*
 * fasta_parse_file
 * ----------------
 * Opens a FASTA file and reads the FIRST sequence it contains.
 *
 * FASTA format:
 *   >header line (ignored beyond storing in DnaSequence.header)
 *   ATGCATGC...   (may span multiple lines)
 *   ATGCATGC...
 *
 * Returns: heap-allocated DnaSequence* on success, NULL on failure.
 * Caller must call dna_sequence_free().
 *
 * Sets *err_out to a DnaError code. Pass NULL to ignore.
 */
DnaSequence* fasta_parse_file(const char* filepath, DnaError* err_out);

/*
 * fasta_parse_string
 * ------------------
 * Same as fasta_parse_file but reads from an in-memory string.
 * Useful for testing without disk I/O.
 */
DnaSequence* fasta_parse_string(const char* fasta_text, DnaError* err_out);

/*
 * fasta_parse_raw
 * ---------------
 * Parses a raw DNA string (no FASTA header).
 * Strips whitespace/newlines, uppercases, validates.
 * Returns heap-allocated DnaSequence* or NULL.
 */
DnaSequence* fasta_parse_raw(const char* raw_seq, DnaError* err_out);

/* ============================================================
   SEQUENCE MANAGEMENT
   ============================================================ */

/*
 * dna_sequence_create
 * -------------------
 * Allocates a new DnaSequence from a plain C string.
 * Validates, uppercases, and copies the data.
 * Returns NULL on failure.
 */
DnaSequence* dna_sequence_create(const char* seq);

/*
 * dna_sequence_free
 * -----------------
 * Frees a DnaSequence and its internal data buffer.
 * Safe to call with NULL (no-op).
 */
void dna_sequence_free(DnaSequence* ds);

/*
 * dna_sequence_copy
 * -----------------
 * Returns a deep copy of a DnaSequence.
 * Returns NULL if ds is NULL or allocation fails.
 */
DnaSequence* dna_sequence_copy(const DnaSequence* ds);

/* ============================================================
   STRING OPERATIONS
   ============================================================ */

/*
 * dna_reverse_complement
 * ----------------------
 * Returns the reverse complement of seq.
 *   A ↔ T, G ↔ C, then reverse the whole string.
 *
 * Example: ATGC → complement: TACG → reverse: GCAT
 *
 * Returns heap-allocated string. Caller must free().
 * Returns NULL if seq is NULL or invalid.
 */
char* dna_reverse_complement(const char* seq);

/*
 * dna_uppercase
 * -------------
 * Converts seq to uppercase in-place.
 * Only affects lowercase a,t,g,c — other chars untouched.
 * Returns seq for convenience (same pointer).
 */
char* dna_uppercase(char* seq);

/*
 * dna_concat_with_separator
 * -------------------------
 * Concatenates two DNA sequences with a separator character.
 * Used by suffix_tree_lcs: "ATGC" + "#" + "TACG" → "ATGC#TACG$"
 *
 * sep is typically '#' (must not be A,T,G,C).
 * The result is always terminated with '$' for suffix tree.
 *
 * Returns heap-allocated string. Caller must free().
 */
char* dna_concat_with_separator(const char* seq1, const char* seq2, char sep);

/*
 * dna_subsequence
 * ---------------
 * Returns a heap-allocated substring of seq from index start, length len.
 * Returns NULL if out of bounds or allocation fails.
 * Caller must free().
 */
char* dna_subsequence(const char* seq, size_t start, size_t len);

/* ============================================================
   SIMILARITY METRICS
   ============================================================ */

/*
 * dna_similarity_simple
 * ---------------------
 * Fast pairwise similarity: counts matching bases at same positions.
 * Does NOT align — just compares character by character up to min(|a|,|b|).
 *
 * Returns percentage 0.0 – 100.0.
 * Used by union_find for rough species clustering.
 *
 * For accurate similarity, use nw_score_only() from needleman_wunsch.h.
 */
double dna_similarity_simple(const char* a, const char* b);

/*
 * dna_hamming_distance
 * --------------------
 * Counts positions where a and b differ (must be same length).
 * Returns -1 if lengths differ.
 * Used by trie_errors for mismatch counting.
 */
int dna_hamming_distance(const char* a, const char* b);

/*
 * dna_gc_content
 * --------------
 * Returns GC content as a percentage (0.0 – 100.0).
 * GC content = (G + C) / total_bases × 100
 * Used as a quick descriptor in the dataset explorer.
 */
double dna_gc_content(const char* seq);

/* ============================================================
   COMPLEMENT HELPERS
   ============================================================ */

/*
 * dna_complement_base
 * -------------------
 * Returns the complement of a single base.
 *   A↔T, G↔C
 * Returns '\0' for invalid input.
 */
char dna_complement_base(char base);

/* ============================================================
   K-MER UTILITIES
   ============================================================ */

/*
 * dna_count_kmers
 * ---------------
 * Counts all k-mers of length k in seq.
 * Returns a heap-allocated array of (4^k) counts indexed by base_to_index.
 *
 * For k=2, index = base_to_index(b1)*4 + base_to_index(b2)
 * For k=3, index = base_to_index(b1)*16 + base_to_index(b2)*4 + base_to_index(b3)
 *
 * Caller must free() the returned array.
 * Returns NULL if k > 10 (too large) or allocation fails.
 */
int* dna_count_kmers(const char* seq, int k);

/*
 * dna_kmer_index
 * --------------
 * Computes the integer index for a k-mer string.
 * kmer must be exactly k characters, all valid bases.
 * Returns -1 on error.
 */
int dna_kmer_index(const char* kmer, int k);

/* ============================================================
   ERROR REPORTING
   ============================================================ */

/*
 * dna_error_string
 * ----------------
 * Returns a human-readable string for a DnaError code.
 * Return value is a string literal — do NOT free.
 */
const char* dna_error_string(DnaError err);

/*
 * dna_last_error
 * --------------
 * Returns the last DnaError set by any utility function.
 * Thread-local in a thread-safe build.
 */
DnaError dna_last_error(void);

/* ============================================================
   DEBUG / PRINT HELPERS
   (Only active when DNA_DEBUG is defined — not in release build)
   ============================================================ */

#ifdef DNA_DEBUG
void dna_print_sequence(const DnaSequence* ds);
void dna_print_bases(const char* seq, size_t max_len);
#endif

#endif /* DNA_UTILS_H */