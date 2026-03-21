#ifndef NEEDLEMAN_WUNSCH_H
#define NEEDLEMAN_WUNSCH_H

/*
 * needleman_wunsch.h
 * ------------------
 * Global sequence alignment using the Needleman-Wunsch dynamic programming
 * algorithm. Produces a gapped alignment string and alignment statistics.
 *
 * Two variants:
 *   nw_align()      — full DP with traceback  (sequences up to ~10,000 bp)
 *   nw_score_only() — two-row DP, score only  (large sequences, species scoring)
 *
 * Dependencies: dna_utils.h
 *
 * Build (object file):
 *   gcc -O2 -Wall -std=c11 -fPIC -c needleman_wunsch.c -o needleman_wunsch.o
 *
 * Self-test:
 *   gcc -O2 -Wall -std=c11 -DNW_TEST dna_utils.c needleman_wunsch.c -o test_nw && ./test_nw
 */

#include <stddef.h>
#include <limits.h>
#include "dna_utils.h"

/* ============================================================
   SCORING DEFAULTS
   ============================================================ */

#define NW_DEFAULT_MATCH      2
#define NW_DEFAULT_MISMATCH  -1
#define NW_DEFAULT_GAP       -2

/* Maximum sequence length for full DP with traceback.
   Above this, nw_align() returns NULL — use nw_score_only(). */
#define NW_MAX_FULL_LEN      10000

/* Gap character in aligned output strings */
#define NW_GAP_CHAR          '-'

/* ============================================================
   STRUCTS
   ============================================================ */

/*
 * NwParams — scoring parameters for a single alignment run.
 */
typedef struct {
    int match;       /* reward for matching bases    (positive) */
    int mismatch;    /* penalty for mismatched bases (negative) */
    int gap;         /* penalty per gap character   (negative) */
} NwParams;

/*
 * AlignmentResult — complete output of nw_align().
 * All pointer fields are heap-allocated.
 * Free with nw_result_free().
 */
typedef struct {
    char*    aligned_a;    /* seq A with NW_GAP_CHAR inserted  */
    char*    aligned_b;    /* seq B with NW_GAP_CHAR inserted  */
    int      score;        /* final NW alignment score         */
    int      matches;      /* identical base pair count        */
    int      mismatches;   /* substitution count               */
    int      gaps;         /* total gap characters (both seqs) */
    int      length;       /* alignment length (incl. gaps)    */
    double   identity;     /* matches / length * 100.0         */
    double   similarity;   /* (matches + mismatches) / length  */
    NwParams params;       /* params used for this alignment   */
} AlignmentResult;

/* ============================================================
   CORE API
   ============================================================ */

/*
 * nw_default_params — returns standard NW scoring parameters.
 */
NwParams nw_default_params(void);

/*
 * nw_align
 * --------
 * Full global alignment with traceback.
 * Both sequences must be valid DNA, length <= NW_MAX_FULL_LEN.
 *
 * Returns heap-allocated AlignmentResult* or NULL on error.
 * Caller must free with nw_result_free().
 *
 * Time:  O(n x m)
 * Space: O(n x m)  (full DP matrix stored for traceback)
 */
AlignmentResult* nw_align(const char* seq_a,
                           const char* seq_b,
                           NwParams    params);

/*
 * nw_score_only
 * -------------
 * Alignment score without traceback. Uses only two DP rows.
 * Use for large sequences or bulk species comparisons.
 *
 * Returns alignment score, or INT_MIN on error.
 *
 * Time:  O(n x m)
 * Space: O(min(n, m))
 */
int nw_score_only(const char* seq_a,
                  const char* seq_b,
                  NwParams    params);

/*
 * nw_result_free — frees AlignmentResult and all its fields.
 * Safe to call with NULL.
 */
void nw_result_free(AlignmentResult* result);

/* ============================================================
   CONVENIENCE WRAPPERS
   ============================================================ */

/* nw_align with NW_DEFAULT_MATCH / MISMATCH / GAP */
AlignmentResult* nw_align_default(const char* seq_a, const char* seq_b);

/* Returns percent identity after default alignment, or -1.0 on error */
double nw_identity(const char* seq_a, const char* seq_b);

/* ============================================================
   JSON SERIALISATION  (consumed by Flask / c_bridge.py)
   ============================================================ */

/*
 * nw_result_to_json
 * -----------------
 * Converts AlignmentResult to a heap-allocated JSON string.
 *
 * Output format:
 * {
 *   "aligned_a":  "ATGC-ATGC",
 *   "aligned_b":  "ATGCTATGC",
 *   "score":      14,
 *   "matches":    8,
 *   "mismatches": 0,
 *   "gaps":       1,
 *   "length":     9,
 *   "identity":   88.89,
 *   "similarity": 100.0
 * }
 *
 * Returns heap-allocated string. Caller must free().
 */
char* nw_result_to_json(const AlignmentResult* result);

#endif /* NEEDLEMAN_WUNSCH_H */