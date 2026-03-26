"""
routes/trie.py
--------------
Flask blueprint for SNP fuzzy matching and mutation scanning.
All computation delegated to trie_errors.c via c_bridge.py.

Endpoints
---------
POST /api/snp/fuzzy
    Fuzzy-match a short DNA snippet against all stored SNP markers.
    Allows 0, 1, or 2 base mismatches (Trie with Errors).
    Used by: Feature 9 (Disease Risk Demo)

    Body:   { "query": "GTGCACCTGACTCCTGTG", "max_errors": 1 }
    Returns:
    {
      "query":      "GTGCACCTGACTCCTGTG",
      "max_errors": 1,
      "count":      2,
      "matches": [
        { "snp_id": "rs334", "disease": "Sickle Cell Anemia",
          "severity": 2, "mismatches": 0, "position": 0 },
        { "snp_id": "rs11549407", "disease": "Beta-Thalassemia",
          "severity": 1, "mismatches": 1, "position": 0 }
      ]
    }

POST /api/snp/exact
    Exact lookup of a single SNP sequence.
    Used by: internal validation, Feature 6

    Body:   { "query": "GTGCACCTGACTCCTGTG" }
    Returns: { "found": true, "query": "GTGCACCTGACTCCTGTG" }

POST /api/mutations/scan
    Scan a full gene sequence for all embedded SNP markers.
    Slides a window across the sequence at every position.
    Used by: Feature 6 (Mutation & Disease Detection)

    Body:   { "sequence": "ATGGTGCACCTGACT...", "max_errors": 1 }
    Returns:
    {
      "sequence_length": 626,
      "max_errors":      1,
      "count":           2,
      "matches": [
        { "snp_id": "rs334", "disease": "Sickle Cell Anemia",
          "severity": 2, "mismatches": 0, "position": 17 },
        ...
      ]
    }

GET /api/snp/stats
    Returns stats about the loaded trie (SNP count, status).
    Used by: Dataset Explorer, health dashboard

    Returns: { "loaded": true, "snp_count": 8, "source": "clinvar_snps.tsv" }
"""

from flask import Blueprint, request, jsonify
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import c_bridge as bridge

trie_bp = Blueprint("trie", __name__)


# ============================================================
# HELPERS  (shared with alignment.py style)
# ============================================================

def _bad_request(msg: str):
    return jsonify({"error": msg}), 400


def _require_json(*fields):
    """Returns (data, error_response). error_response is None on success."""
    data = request.get_json(silent=True)
    if data is None:
        return None, _bad_request("Request body must be JSON")
    for f in fields:
        if f not in data or data[f] is None:
            return None, _bad_request(f"Missing field: '{f}'")
    return data, None


def _trie_ready():
    """Check trie is loaded. Returns (True, None) or (False, error_response)."""
    if bridge._trie_handle is None:
        return False, (
            jsonify({
                "error": "SNP database not loaded. "
                         "Call trie_init() in app.py startup, "
                         "or POST /api/snp/load first."
            }), 503
        )
    return True, None


def _clamp_errors(val) -> int:
    """Parse and clamp max_errors to [0, 2]."""
    try:
        v = int(val)
    except (TypeError, ValueError):
        v = 1
    return max(0, min(2, v))


# ============================================================
# POST /api/snp/fuzzy
# ============================================================

@trie_bp.route("/api/snp/fuzzy", methods=["POST"])
def api_snp_fuzzy():
    """
    Fuzzy-match a short DNA snippet against all stored SNP markers.

    Request JSON:
        query       string  required  DNA snippet to match (A,T,G,C)
        max_errors  int     optional  0-2 mismatches allowed (default 1)

    Response JSON:
        query       string  echoed back
        max_errors  int     echoed back
        count       int     number of SNP matches found
        matches     list    array of match objects
    """
    ready, err = _trie_ready()
    if not ready:
        return err

    data, err = _require_json("query")
    if err:
        return err

    query      = data["query"].strip().upper()
    max_errors = _clamp_errors(data.get("max_errors", 1))

    if not bridge.validate(query):
        return _bad_request(
            "query contains invalid characters — only A, T, G, C allowed"
        )

    if len(query) < 4:
        return _bad_request("query too short — minimum 4 bases required")

    matches = bridge.trie_fuzzy(query, max_errors=max_errors)

    return jsonify({
        "query":      query,
        "max_errors": max_errors,
        "count":      len(matches),
        "matches":    matches,
    }), 200


# ============================================================
# POST /api/snp/exact
# ============================================================

@trie_bp.route("/api/snp/exact", methods=["POST"])
def api_snp_exact():
    """
    Exact SNP sequence lookup.

    Request JSON:
        query  string  required  DNA sequence to check

    Response JSON:
        query  string  echoed back
        found  bool    true if exact match exists in trie
    """
    ready, err = _trie_ready()
    if not ready:
        return err

    data, err = _require_json("query")
    if err:
        return err

    query = data["query"].strip().upper()

    if not bridge.validate(query):
        return _bad_request("Invalid DNA sequence")

    found = bridge.trie_exact(query)

    return jsonify({
        "query": query,
        "found": found,
    }), 200


# ============================================================
# POST /api/mutations/scan
# ============================================================

@trie_bp.route("/api/mutations/scan", methods=["POST"])
def api_mutations_scan():
    """
    Scan a full gene sequence for all embedded SNP markers.
    Slides a window of each SNP length across every position.

    Request JSON:
        sequence    string  required  Full DNA sequence to scan
        max_errors  int     optional  0-2 mismatches allowed (default 1)

    Response JSON:
        sequence_length  int   length of input sequence
        max_errors       int   echoed back
        count            int   number of mutations found
        matches          list  array of match objects with position
    """
    ready, err = _trie_ready()
    if not ready:
        return err

    data, err = _require_json("sequence")
    if err:
        return err

    sequence   = data["sequence"].strip().upper()
    max_errors = _clamp_errors(data.get("max_errors", 1))

    if not bridge.validate(sequence):
        return _bad_request(
            "sequence contains invalid characters — only A, T, G, C allowed"
        )

    if len(sequence) < 10:
        return _bad_request("sequence too short — minimum 10 bases required")

    matches = bridge.trie_scan_sequence(sequence, max_errors=max_errors)

    # Sort by position ascending for consistent frontend rendering
    matches.sort(key=lambda m: m["position"])

    return jsonify({
        "sequence_length": len(sequence),
        "max_errors":      max_errors,
        "count":           len(matches),
        "matches":         matches,
    }), 200


# ============================================================
# GET /api/snp/stats
# ============================================================

@trie_bp.route("/api/snp/stats", methods=["GET"])
def api_snp_stats():
    """
    Returns info about the currently loaded SNP trie.

    Response JSON:
        loaded     bool   true if trie has been initialised
        snp_count  int    number of SNPs in the trie (0 if not loaded)
        source     string filepath of the loaded TSV (or null)
    """
    loaded = bridge._trie_handle is not None

    # trie->size is not directly exposed yet —
    # we track it at the Python level via a module variable set by trie_init()
    snp_count = _trie_snp_count if loaded else 0

    return jsonify({
        "loaded":    loaded,
        "snp_count": snp_count,
        "source":    _trie_source if loaded else None,
    }), 200


# ============================================================
# POST /api/snp/load
# ============================================================

@trie_bp.route("/api/snp/load", methods=["POST"])
def api_snp_load():
    """
    (Re)load the SNP trie from a TSV file.
    Normally called once at server startup from app.py.
    This endpoint lets you reload without restarting the server.

    Request JSON:
        filepath  string  required  path to clinvar_snps.tsv

    Response JSON:
        loaded     bool  true on success
        snp_count  int   number of SNPs loaded
        filepath   string echoed back
    """
    data, err = _require_json("filepath")
    if err:
        return err

    filepath = data["filepath"].strip()

    if not os.path.exists(filepath):
        return _bad_request(f"File not found: {filepath}")

    count = bridge.trie_init(filepath)

    if count < 0:
        return _bad_request(f"Failed to open file: {filepath}")

    # Update module-level tracking vars
    global _trie_snp_count, _trie_source
    _trie_snp_count = count
    _trie_source    = filepath

    return jsonify({
        "loaded":    True,
        "snp_count": count,
        "filepath":  filepath,
    }), 200


# ============================================================
# MODULE-LEVEL STATE
# Tracks SNP count and source path after trie_init() is called.
# ============================================================

_trie_snp_count: int       = 0
_trie_source:    str | None = None


def init_trie_from_file(filepath: str) -> int:
    """
    Called by app.py at startup to pre-load the SNP trie.

    Example in app.py:
        from routes.trie import init_trie_from_file
        init_trie_from_file("data/clinvar_snps.tsv")

    Returns number of SNPs loaded.
    """
    global _trie_snp_count, _trie_source

    if not os.path.exists(filepath):
        print(f"[trie] WARNING: {filepath} not found — trie not loaded")
        return 0

    count = bridge.trie_init(filepath)

    if count < 0:
        print(f"[trie] ERROR: failed to load {filepath}")
        return 0

    _trie_snp_count = count
    _trie_source    = filepath
    print(f"[trie] Loaded {count} SNPs from {filepath}")
    return count