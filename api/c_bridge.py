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
        here / "c-core" / "libdna.so",
        here / "../c-core" / "libdna.so",
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
# TRIE STRUCTS  (mirror trie_errors.h exactly)
# ============================================================

class TrieMatch(ctypes.Structure):
    """Maps to TrieMatch in trie_errors.h"""
    _fields_ = [
        ("snp_id",         ctypes.c_char_p),
        ("disease_label",  ctypes.c_char_p),
        ("severity",       ctypes.c_int),
        ("mismatches",     ctypes.c_int),
        ("match_position", ctypes.c_int),
    ]


class TrieMatchList(ctypes.Structure):
    """Maps to TrieMatchList in trie_errors.h"""
    _fields_ = [
        ("matches", ctypes.POINTER(TrieMatch)),
        ("count",   ctypes.c_int),
    ]


# ============================================================
# TRIE FUNCTION SIGNATURES
# ============================================================

_lib.trie_create.argtypes = []
_lib.trie_create.restype  = ctypes.c_void_p

_lib.trie_free.argtypes = [ctypes.c_void_p]
_lib.trie_free.restype  = None

_lib.trie_insert.argtypes = [
    ctypes.c_void_p,   # Trie*
    ctypes.c_char_p,   # snp_sequence
    ctypes.c_char_p,   # disease
    ctypes.c_char_p,   # snp_id
    ctypes.c_int,      # severity
]
_lib.trie_insert.restype = None

_lib.trie_search_exact.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
_lib.trie_search_exact.restype  = ctypes.c_int

_lib.trie_search_fuzzy.argtypes = [
    ctypes.c_void_p,   # Trie*
    ctypes.c_char_p,   # query
    ctypes.c_int,      # max_errors
]
_lib.trie_search_fuzzy.restype = ctypes.POINTER(TrieMatchList)

_lib.trie_search_in_sequence.argtypes = [
    ctypes.c_void_p,   # Trie*
    ctypes.c_char_p,   # sequence
    ctypes.c_int,      # max_errors
]
_lib.trie_search_in_sequence.restype = ctypes.POINTER(TrieMatchList)

_lib.trie_load_from_file.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
_lib.trie_load_from_file.restype  = ctypes.c_int

_lib.trie_match_list_free.argtypes = [ctypes.POINTER(TrieMatchList)]
_lib.trie_match_list_free.restype  = None

_lib.trie_match_list_to_json.argtypes = [ctypes.POINTER(TrieMatchList)]
_lib.trie_match_list_to_json.restype  = ctypes.c_char_p


# ============================================================
# TRIE — module-level singleton
# ============================================================
# The trie is loaded ONCE at startup from clinvar_snps.tsv.
# All route calls share this single in-memory trie.
# Call trie_init() from app.py after Flask starts.

_trie_handle = None   # raw C void* pointer


def trie_init(tsv_filepath: str) -> int:
    """
    Load the SNP trie from a TSV file.
    Call once at server startup.

    Returns number of SNPs loaded, or -1 on failure.
    """
    global _trie_handle

    # Free old trie if exists
    if _trie_handle is not None:
        _lib.trie_free(_trie_handle)

    _trie_handle = _lib.trie_create()
    if not _trie_handle:
        return -1

    count = _lib.trie_load_from_file(
        _trie_handle,
        tsv_filepath.encode("utf-8")
    )
    return count


def trie_insert_snp(snp_sequence: str,
                    disease: str,
                    snp_id: str,
                    severity: int) -> None:
    """
    Insert a single SNP into the loaded trie.
    Requires trie_init() to have been called first.
    """
    if _trie_handle is None:
        raise RuntimeError("Trie not initialised — call trie_init() first")
    _lib.trie_insert(
        _trie_handle,
        _encode(snp_sequence),
        disease.encode("utf-8"),
        snp_id.encode("utf-8"),
        ctypes.c_int(severity),
    )


def trie_exact(query: str) -> bool:
    """
    Exact SNP lookup.
    Returns True if query exactly matches a stored SNP sequence.
    """
    if _trie_handle is None or not validate(query):
        return False
    return bool(_lib.trie_search_exact(_trie_handle, _encode(query)))


def trie_fuzzy(query: str, max_errors: int = 1) -> list[dict]:
    """
    Fuzzy SNP search — finds all stored SNPs within max_errors mismatches.

    Returns list of dicts:
        [{ "snp_id": "rs334", "disease": "Sickle Cell Anemia",
           "severity": 2, "mismatches": 0, "position": 0 }, ...]

    Returns [] on no match or error.
    """
    if _trie_handle is None or not validate(query):
        return []

    ptr = _lib.trie_search_fuzzy(
        _trie_handle,
        _encode(query),
        ctypes.c_int(max_errors)
    )
    if not ptr:
        return []

    result = _match_list_to_python(ptr)
    _lib.trie_match_list_free(ptr)
    return result


def trie_scan_sequence(sequence: str, max_errors: int = 1) -> list[dict]:
    """
    Scan a full DNA sequence for any embedded SNP markers.
    Slides a window across the sequence and checks each position.

    Returns list of dicts (same format as trie_fuzzy).
    """
    if _trie_handle is None or not validate(sequence):
        return []

    ptr = _lib.trie_search_in_sequence(
        _trie_handle,
        _encode(sequence),
        ctypes.c_int(max_errors)
    )
    if not ptr:
        return []

    result = _match_list_to_python(ptr)
    _lib.trie_match_list_free(ptr)
    return result


def _match_list_to_python(ptr) -> list[dict]:
    """
    Convert a C TrieMatchList* into a Python list of dicts.
    Internal helper — do not call directly.
    """
    if not ptr:
        return []

    ml = ptr.contents
    results = []

    for i in range(ml.count):
        m = ml.matches[i]
        results.append({
            "snp_id":     m.snp_id.decode("utf-8")        if m.snp_id        else "",
            "disease":    m.disease_label.decode("utf-8")  if m.disease_label else "",
            "severity":   m.severity,
            "mismatches": m.mismatches,
            "position":   m.match_position,
        })

    return results


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

# ============================================================
# C-CORE: Suffix Tree
# ============================================================

_lib.suffix_tree_build.argtypes = [ctypes.c_char_p, ctypes.c_int]
_lib.suffix_tree_build.restype = ctypes.c_void_p

_lib.suffix_tree_search_to_json.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
_lib.suffix_tree_search_to_json.restype = ctypes.c_void_p

_lib.suffix_tree_search.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
_lib.suffix_tree_search.restype = ctypes.c_void_p

_lib.suffix_tree_free.argtypes = [ctypes.c_void_p]
_lib.suffix_tree_free.restype = None

_lib.suffix_tree_search_result_free.argtypes = [ctypes.c_void_p]
_lib.suffix_tree_search_result_free.restype = None

def suffix_tree_build(sequence: str):
    b_seq = sequence.encode('utf-8')
    return _lib.suffix_tree_build(b_seq, len(sequence))

def suffix_tree_search(st_ptr, pattern: str):
    b_pat = pattern.encode('utf-8')
    res_ptr = _lib.suffix_tree_search(st_ptr, b_pat)
    if not res_ptr: return None
    
    json_ptr = _lib.suffix_tree_search_to_json(res_ptr, b_pat)
    if not json_ptr:
        _lib.suffix_tree_search_result_free(res_ptr)
        return None
        
    s = ctypes.cast(json_ptr, ctypes.c_char_p).value.decode('utf-8')
    _libc.free(json_ptr)
    _lib.suffix_tree_search_result_free(res_ptr)
    
    import json
    try:
        return json.loads(s)
    except:
        return s

def suffix_tree_free(st_ptr):
    if st_ptr:
        _lib.suffix_tree_free(st_ptr)


# ============================================================
# C-CORE: Segment Tree
# ============================================================
# (Since the segment tree relies on clinvar loading)

_lib.segment_tree_load_clinvar.argtypes = [ctypes.POINTER(ctypes.c_void_p), ctypes.c_char_p]
_lib.segment_tree_load_clinvar.restype = ctypes.c_int

_lib.segment_tree_query.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_int)]
_lib.segment_tree_query.restype = ctypes.c_void_p

_lib.segment_tree_results_to_json.argtypes = [ctypes.c_void_p, ctypes.c_int]
_lib.segment_tree_results_to_json.restype = ctypes.c_void_p

_lib.segment_tree_results_free.argtypes = [ctypes.c_void_p, ctypes.c_int]
_lib.segment_tree_results_free.restype = None

_st_handle = None

def segment_tree_init(filepath: str) -> int:
    global _st_handle
    if _st_handle:
        _lib.segment_tree_free(_st_handle)
        _st_handle = None
    
    ptr = ctypes.c_void_p()
    b_path = filepath.encode('utf-8')
    res = _lib.segment_tree_load_clinvar(ctypes.byref(ptr), b_path)
    if res > 0:
        _st_handle = ptr
    return res

def segment_tree_query(low: int, high: int):
    if not _st_handle: return None
    
    count = ctypes.c_int(0)
    res_ptr = _lib.segment_tree_query(_st_handle, low, high, ctypes.byref(count))
    if not res_ptr: return []
    
    json_ptr = _lib.segment_tree_results_to_json(res_ptr, count.value)
    if not json_ptr:
        _lib.segment_tree_results_free(res_ptr, count.value)
        return []
        
    s = ctypes.cast(json_ptr, ctypes.c_char_p).value.decode('utf-8')
    _libc.free(json_ptr)
    _lib.segment_tree_results_free(res_ptr, count.value)
    
    import json
    try:
        return json.loads(s)
    except:
        return s


# ============================================================
# C-CORE: Interval Tree
# ============================================================

_lib.interval_tree_create.argtypes = []
_lib.interval_tree_create.restype = ctypes.c_void_p

_lib.interval_tree_load_genes.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
_lib.interval_tree_load_genes.restype = ctypes.c_int

_lib.interval_tree_query.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_int)]
_lib.interval_tree_query.restype = ctypes.c_void_p

_lib.interval_tree_results_to_json.argtypes = [ctypes.c_void_p, ctypes.c_int]
_lib.interval_tree_results_to_json.restype = ctypes.c_void_p

_lib.interval_tree_results_free.argtypes = [ctypes.c_void_p, ctypes.c_int]
_lib.interval_tree_results_free.restype = None

_it_handle = None

def interval_tree_init(filepath: str) -> int:
    global _it_handle
    if not _it_handle:
        _it_handle = _lib.interval_tree_create()
    
    b_path = filepath.encode('utf-8')
    return _lib.interval_tree_load_genes(_it_handle, b_path)

def interval_tree_query(low: int, high: int):
    if not _it_handle: return None
    
    count = ctypes.c_int(0)
    res_ptr = _lib.interval_tree_query(_it_handle, low, high, ctypes.byref(count))
    if not res_ptr: return []
    
    json_ptr = _lib.interval_tree_results_to_json(res_ptr, count.value)
    if not json_ptr:
        _lib.interval_tree_results_free(res_ptr, count.value)
        return []
        
    s = ctypes.cast(json_ptr, ctypes.c_char_p).value.decode('utf-8')
    _libc.free(json_ptr)
    _lib.interval_tree_results_free(res_ptr, count.value)
    
    import json
    try:
        return json.loads(s)
    except:
        return s


# ============================================================
# C-CORE: Union Find
# ============================================================

_lib.union_find_create.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_char_p)]
_lib.union_find_create.restype = ctypes.c_void_p

_lib.union_find_cluster_species.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.c_double]
_lib.union_find_cluster_species.restype = None

_lib.union_find_clusters_to_json.argtypes = [ctypes.c_void_p]
_lib.union_find_clusters_to_json.restype = ctypes.c_void_p

_lib.union_find_merge_log_to_json.argtypes = [ctypes.c_void_p]
_lib.union_find_merge_log_to_json.restype = ctypes.c_void_p

_lib.union_find_free.argtypes = [ctypes.c_void_p]
_lib.union_find_free.restype = None

def union_find_process(species_list, similarity_matrix, threshold=95.0):
    n = len(species_list)
    
    # Create array of char pointers for species
    c_species = (ctypes.c_char_p * n)()
    for i, s in enumerate(species_list):
        c_species[i] = s.encode('utf-8')
        
    uf_ptr = _lib.union_find_create(n, c_species)
    if not uf_ptr: return None
    
    # Create 2D array of doubles
    # Python list of lists to C double**
    c_matrix = (ctypes.POINTER(ctypes.c_double) * n)()
    # Keep references to the rows to avoid garbage collection!!
    rows = []
    for i in range(n):
        row = (ctypes.c_double * n)(*similarity_matrix[i])
        c_matrix[i] = ctypes.cast(row, ctypes.POINTER(ctypes.c_double))
        rows.append(row)
        
    _lib.union_find_cluster_species(uf_ptr, c_matrix, threshold)
    
    clusters_json_ptr = _lib.union_find_clusters_to_json(uf_ptr)
    log_json_ptr = _lib.union_find_merge_log_to_json(uf_ptr)
    
    import json
    out = {}
    
    if clusters_json_ptr:
        cstr = ctypes.cast(clusters_json_ptr, ctypes.c_char_p).value.decode('utf-8')
        _libc.free(clusters_json_ptr)
        try: out["clusters"] = json.loads(cstr)
        except: out["clusters"] = cstr
        
    if log_json_ptr:
        lstr = ctypes.cast(log_json_ptr, ctypes.c_char_p).value.decode('utf-8')
        _libc.free(log_json_ptr)
        try: out["merge_log"] = json.loads(lstr)
        except: out["merge_log"] = lstr
        
    _lib.union_find_free(uf_ptr)
    return out
