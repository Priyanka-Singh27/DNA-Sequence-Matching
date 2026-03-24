"""
c_bridge.py
-----------
ctypes bridge between Python/Flask and the C shared library (libdna.so).

Usage:
    from c_bridge import align, score_only

All functions return plain Python dicts/values — no ctypes types leak out.
The bridge handles all type marshalling, error checking, and memory freeing.

Library path resolution order:
  1. LIBDNA_PATH environment variable
  2. ./c-core/libdna.so  (relative to this file)
  3. ../c-core/libdna.so
"""

import ctypes
import os
import json
from pathlib import Path


# ============================================================
# LIBRARY LOADING
# ============================================================

def _find_lib() -> Path:
    """Locate libdna.so — checks env var then relative paths."""
    env = os.environ.get("LIBDNA_PATH")
    if env and Path(env).exists():
        return Path(env)

    here = Path(__file__).parent
    candidates = [
        here / "c-core" / "libdna.dll",   # Windows
        here / "c-core" / "libdna.so",    # Linux/Mac
        here / "../c-core" / "libdna.dll",
        here / "../c-core" / "libdna.so",
        here / "libdna.dll",
        here / "libdna.so",
    ]
    for p in candidates:
        if p.resolve().exists():
            return p.resolve()

    raise FileNotFoundError(
        "libdna.so not found. Build it with:\n"
        "  cd c-core && gcc -O2 -fPIC -shared dna_utils.c needleman_wunsch.c "
        "-o libdna.so -lm"
    )


_lib_path = _find_lib()
_lib      = ctypes.CDLL(str(_lib_path))


# ============================================================
# STRUCT DEFINITIONS  (must mirror C exactly)
# ============================================================

class NwParams(ctypes.Structure):
    """Maps to NwParams in needleman_wunsch.h"""
    _fields_ = [
        ("match",    ctypes.c_int),
        ("mismatch", ctypes.c_int),
        ("gap",      ctypes.c_int),
    ]


class AlignmentResult(ctypes.Structure):
    """Maps to AlignmentResult in needleman_wunsch.h"""
    _fields_ = [
        ("aligned_a",   ctypes.c_char_p),
        ("aligned_b",   ctypes.c_char_p),
        ("score",       ctypes.c_int),
        ("matches",     ctypes.c_int),
        ("mismatches",  ctypes.c_int),
        ("gaps",        ctypes.c_int),
        ("length",      ctypes.c_int),
        ("identity",    ctypes.c_double),
        ("similarity",  ctypes.c_double),
        ("params",      NwParams),
    ]


# ============================================================
# FUNCTION SIGNATURES
# ============================================================

# --- dna_utils ---

_lib.dna_validate.argtypes  = [ctypes.c_char_p]
_lib.dna_validate.restype   = ctypes.c_int

_lib.dna_gc_content.argtypes = [ctypes.c_char_p]
_lib.dna_gc_content.restype  = ctypes.c_double

_lib.dna_similarity_simple.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
_lib.dna_similarity_simple.restype  = ctypes.c_double

_lib.dna_hamming_distance.argtypes  = [ctypes.c_char_p, ctypes.c_char_p]
_lib.dna_hamming_distance.restype   = ctypes.c_int

_lib.dna_reverse_complement.argtypes = [ctypes.c_char_p]
_lib.dna_reverse_complement.restype  = ctypes.c_char_p

_lib.base_to_index.argtypes = [ctypes.c_char]
_lib.base_to_index.restype  = ctypes.c_int

# --- needleman_wunsch ---

_lib.nw_default_params.argtypes = []
_lib.nw_default_params.restype  = NwParams

_lib.nw_align.argtypes = [ctypes.c_char_p, ctypes.c_char_p, NwParams]
_lib.nw_align.restype  = ctypes.POINTER(AlignmentResult)

_lib.nw_score_only.argtypes = [ctypes.c_char_p, ctypes.c_char_p, NwParams]
_lib.nw_score_only.restype  = ctypes.c_int

_lib.nw_result_free.argtypes = [ctypes.POINTER(AlignmentResult)]
_lib.nw_result_free.restype  = None

_lib.nw_result_to_json.argtypes = [ctypes.POINTER(AlignmentResult)]
_lib.nw_result_to_json.restype  = ctypes.c_char_p

_lib.nw_identity.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
_lib.nw_identity.restype  = ctypes.c_double


# ============================================================
# PYTHON-FRIENDLY WRAPPERS
# ============================================================

# Default NW scoring constants (mirror C defaults)
DEFAULT_MATCH    =  2
DEFAULT_MISMATCH = -1
DEFAULT_GAP      = -2


def _encode(s: str) -> bytes:
    """UTF-8 encode and uppercase a Python str for passing to C."""
    return s.upper().encode("utf-8")


def validate(seq: str) -> bool:
    """
    Returns True if every character in seq is A, T, G, or C.
    Calls dna_validate() in C — fast lookup-table implementation.
    """
    if not seq:
        return False
    return bool(_lib.dna_validate(_encode(seq)))


def gc_content(seq: str) -> float:
    """
    Returns GC content as a percentage (0.0 – 100.0).
    Returns 0.0 on invalid input.
    """
    if not validate(seq):
        return 0.0
    return _lib.dna_gc_content(_encode(seq))


def similarity_simple(seq_a: str, seq_b: str) -> float:
    """
    Fast pairwise similarity by position-wise comparison (no alignment).
    Returns 0.0 – 100.0.
    """
    if not validate(seq_a) or not validate(seq_b):
        return 0.0
    return _lib.dna_similarity_simple(_encode(seq_a), _encode(seq_b))


def hamming_distance(seq_a: str, seq_b: str) -> int:
    """
    Hamming distance between two equal-length sequences.
    Returns -1 if lengths differ or input invalid.
    """
    return _lib.dna_hamming_distance(_encode(seq_a), _encode(seq_b))


def reverse_complement(seq: str) -> str | None:
    """
    Returns the reverse complement of seq.
    Returns None on invalid input.
    """
    if not validate(seq):
        return None
    result = _lib.dna_reverse_complement(_encode(seq))
    return result.decode("utf-8") if result else None


def align(
    seq_a: str,
    seq_b: str,
    match: int    = DEFAULT_MATCH,
    mismatch: int = DEFAULT_MISMATCH,
    gap: int      = DEFAULT_GAP,
) -> dict | None:
    """
    Full Needleman-Wunsch global alignment with traceback.

    Parameters
    ----------
    seq_a, seq_b : DNA strings (A/T/G/C, case-insensitive)
    match        : reward for matching bases       (default +2)
    mismatch     : penalty for mismatched bases    (default -1)
    gap          : penalty per gap character       (default -2)

    Returns
    -------
    dict with keys:
        aligned_a   str   seq_a with gap characters inserted
        aligned_b   str   seq_b with gap characters inserted
        score       int   NW alignment score
        matches     int   number of identical base pairs
        mismatches  int   number of substitutions
        gaps        int   total gap characters
        length      int   alignment length
        identity    float percent identity (0-100)
        similarity  float percent similarity (matches+mismatches / length)

    Returns None if either sequence is invalid or too long (> 10 000 bp).
    """
    if not validate(seq_a) or not validate(seq_b):
        return None

    params          = NwParams()
    params.match    = match
    params.mismatch = mismatch
    params.gap      = gap

    ptr = _lib.nw_align(_encode(seq_a), _encode(seq_b), params)
    if not ptr:
        return None

    r = ptr.contents
    result = {
        "aligned_a":  r.aligned_a.decode("utf-8")  if r.aligned_a  else "",
        "aligned_b":  r.aligned_b.decode("utf-8")  if r.aligned_b  else "",
        "score":      r.score,
        "matches":    r.matches,
        "mismatches": r.mismatches,
        "gaps":       r.gaps,
        "length":     r.length,
        "identity":   round(r.identity,   2),
        "similarity": round(r.similarity, 2),
    }

    _lib.nw_result_free(ptr)   # free C heap memory
    return result


def score_only(
    seq_a: str,
    seq_b: str,
    match: int    = DEFAULT_MATCH,
    mismatch: int = DEFAULT_MISMATCH,
    gap: int      = DEFAULT_GAP,
) -> int | None:
    """
    Alignment score without traceback.
    Uses two-row DP — suitable for large sequences.

    Returns integer score, or None on error.
    """
    if not validate(seq_a) or not validate(seq_b):
        return None

    params          = NwParams()
    params.match    = match
    params.mismatch = mismatch
    params.gap      = gap

    score = _lib.nw_score_only(_encode(seq_a), _encode(seq_b), params)

    # INT_MIN signals error from C
    if score == -2147483648:
        return None
    return score


def identity(seq_a: str, seq_b: str) -> float | None:
    """
    Returns percent identity between two sequences after alignment.
    Convenience wrapper — aligns internally, extracts identity, frees.

    Returns None on error.
    """
    if not validate(seq_a) or not validate(seq_b):
        return None
    val = _lib.nw_identity(_encode(seq_a), _encode(seq_b))
    return round(val, 2) if val >= 0 else None


# ============================================================
# BRIDGE SELF-TEST
# ============================================================

if __name__ == "__main__":
    print("=" * 55)
    print("  c_bridge.py self-test")
    print(f"  Library: {_lib_path}")
    print("=" * 55)

    _counts = [0, 0]  # [run, passed]

    def test(name, condition):
        _counts[0] += 1
        status = "PASS" if condition else "FAIL"
        if condition:
            _counts[1] += 1
        print(f"  {'%-52s' % name} {status}")

    # --- validate ---
    print("\n[validate]")
    test("valid sequence",        validate("ATGC"))
    test("lowercase accepted",    validate("atgc"))
    test("invalid base rejected", not validate("ATGX"))
    test("empty rejected",        not validate(""))
    test("None rejected",         not validate(None))

    # --- gc_content ---
    print("\n[gc_content]")
    test("ATGC = 50%",  gc_content("ATGC") == 50.0)
    test("AAAA = 0%",   gc_content("AAAA") == 0.0)
    test("GCGC = 100%", gc_content("GCGC") == 100.0)

    # --- similarity_simple ---
    print("\n[similarity_simple]")
    test("identical = 100%", similarity_simple("ATGC", "ATGC") == 100.0)
    test("all diff  = 0%",   similarity_simple("AAAA", "TTTT") == 0.0)

    # --- hamming_distance ---
    print("\n[hamming_distance]")
    test("identical = 0",         hamming_distance("ATGC", "ATGC") == 0)
    test("one mismatch = 1",      hamming_distance("ATGC", "TTGC") == 1)
    test("diff lengths = -1",     hamming_distance("ATG",  "ATGC") == -1)

    # --- reverse_complement ---
    print("\n[reverse_complement]")
    test("ATGC -> GCAT", reverse_complement("ATGC") == "GCAT")
    test("A -> T",       reverse_complement("A")    == "T")

    # --- align ---
    print("\n[align]")
    r = align("ATGC", "ATGC")
    test("identical: result not None",    r is not None)
    test("identical: score = 8",          r["score"] == 8)
    test("identical: identity = 100.0",   r["identity"] == 100.0)
    test("identical: no gaps",            r["gaps"] == 0)

    r = align("ATGGTGCACCTGACTCCTGAGGAGAAGTCT",
              "ATGGTGCACCTGACTCCTGTGGAGAAGTCT")
    test("HBB: 1 mismatch",              r is not None and r["mismatches"] == 1)
    test("HBB: identity ~96.67",         r is not None and abs(r["identity"] - 96.67) < 0.01)

    r = align("ATGCAT", "ATCAT")
    test("deletion: gaps > 0",           r is not None and r["gaps"] > 0)

    test("invalid seq -> None",           align("ATXC", "ATGC") is None)
    test("custom params",                 align("AAAA", "AAAA", match=5) is not None)

    # --- score_only ---
    print("\n[score_only]")
    s_full  = align("ATGCATGC", "ATGCATGC")["score"]
    s_fast  = score_only("ATGCATGC", "ATGCATGC")
    test("score_only == full score",     s_full == s_fast)
    test("invalid -> None",              score_only("ATXC", "ATGC") is None)

    # --- identity ---
    print("\n[identity]")
    test("identical = 100.0",            identity("ATGC", "ATGC") == 100.0)
    test("None on invalid",              identity("ATXC", "ATGC") is None)

    print(f"\n{'=' * 55}")
    print(f"  Results: {_counts[1]} / {_counts[0]} tests passed")
    print(f"{'=' * 55}")