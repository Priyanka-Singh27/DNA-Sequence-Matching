"""
routes/alignment.py
-------------------
Flask blueprint for all alignment-related API endpoints.

Endpoints
---------
POST /api/align
    Full NW alignment with traceback.
    Body:  { "seq_a": "ATGC...", "seq_b": "ATGC...",
             "match": 2, "mismatch": -1, "gap": -2 }
    Returns: AlignmentResult JSON

POST /api/score
    Score-only alignment (no traceback, fast for large seqs).
    Body:  { "seq_a": "...", "seq_b": "..." }
    Returns: { "score": 42 }

POST /api/identity
    Returns percent identity between two sequences.
    Body:  { "seq_a": "...", "seq_b": "..." }
    Returns: { "identity": 96.67 }

POST /api/gc
    Returns GC content of a sequence.
    Body:  { "seq": "ATGC..." }
    Returns: { "gc_content": 50.0 }

POST /api/validate
    Validates a DNA sequence.
    Body:  { "seq": "ATGC..." }
    Returns: { "valid": true }

POST /api/reverse_complement
    Returns reverse complement.
    Body:  { "seq": "ATGC..." }
    Returns: { "result": "GCAT" }
"""

from flask import Blueprint, request, jsonify
import sys
import os

# Ensure parent dir is on path so c_bridge can be found
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import c_bridge as bridge

alignment_bp = Blueprint("alignment", __name__)


# ============================================================
# HELPERS
# ============================================================

def _bad_request(msg: str):
    return jsonify({"error": msg}), 400


def _require_json(*fields):
    """Returns (data, error_response). error_response is None on success."""
    data = request.get_json(silent=True)
    if data is None:
        return None, _bad_request("Request body must be JSON")
    for f in fields:
        if f not in data or not data[f]:
            return None, _bad_request(f"Missing or empty field: '{f}'")
    return data, None


# ============================================================
# POST /api/align
# ============================================================

@alignment_bp.route("/api/align", methods=["POST"])
def api_align():
    """
    Full NW alignment with traceback.

    Request JSON:
        seq_a     string  required  DNA sequence A
        seq_b     string  required  DNA sequence B
        match     int     optional  default 2
        mismatch  int     optional  default -1
        gap       int     optional  default -2

    Response JSON:
        aligned_a   string
        aligned_b   string
        score       int
        matches     int
        mismatches  int
        gaps        int
        length      int
        identity    float
        similarity  float
    """
    data, err = _require_json("seq_a", "seq_b")
    if err:
        return err

    seq_a    = data["seq_a"].strip().upper()
    seq_b    = data["seq_b"].strip().upper()
    match    = int(data.get("match",    bridge.DEFAULT_MATCH))
    mismatch = int(data.get("mismatch", bridge.DEFAULT_MISMATCH))
    gap      = int(data.get("gap",      bridge.DEFAULT_GAP))

    # Validate before calling C
    if not bridge.validate(seq_a):
        return _bad_request("seq_a contains invalid characters (only A,T,G,C allowed)")
    if not bridge.validate(seq_b):
        return _bad_request("seq_b contains invalid characters (only A,T,G,C allowed)")

    result = bridge.align(seq_a, seq_b, match=match, mismatch=mismatch, gap=gap)
    if result is None:
        return _bad_request(
            "Alignment failed — sequences may be too long (max 10,000 bp for full alignment). "
            "Use /api/score for longer sequences."
        )

    return jsonify(result), 200


# ============================================================
# POST /api/score
# ============================================================

@alignment_bp.route("/api/score", methods=["POST"])
def api_score():
    """
    Score-only alignment — fast, no traceback.
    Suitable for sequences up to the full genome length.

    Request JSON:
        seq_a     string  required
        seq_b     string  required
        match     int     optional
        mismatch  int     optional
        gap       int     optional

    Response JSON:
        score  int
    """
    data, err = _require_json("seq_a", "seq_b")
    if err:
        return err

    seq_a    = data["seq_a"].strip().upper()
    seq_b    = data["seq_b"].strip().upper()
    match    = int(data.get("match",    bridge.DEFAULT_MATCH))
    mismatch = int(data.get("mismatch", bridge.DEFAULT_MISMATCH))
    gap      = int(data.get("gap",      bridge.DEFAULT_GAP))

    if not bridge.validate(seq_a):
        return _bad_request("seq_a contains invalid characters")
    if not bridge.validate(seq_b):
        return _bad_request("seq_b contains invalid characters")

    score = bridge.score_only(seq_a, seq_b, match=match, mismatch=mismatch, gap=gap)
    if score is None:
        return _bad_request("Scoring failed")

    return jsonify({"score": score}), 200


# ============================================================
# POST /api/identity
# ============================================================

@alignment_bp.route("/api/identity", methods=["POST"])
def api_identity():
    """
    Percent identity after full NW alignment.

    Request JSON:
        seq_a  string  required
        seq_b  string  required

    Response JSON:
        identity  float  0.0 – 100.0
    """
    data, err = _require_json("seq_a", "seq_b")
    if err:
        return err

    seq_a = data["seq_a"].strip().upper()
    seq_b = data["seq_b"].strip().upper()

    val = bridge.identity(seq_a, seq_b)
    if val is None:
        return _bad_request("Could not compute identity — check sequence length and characters")

    return jsonify({"identity": val}), 200


# ============================================================
# POST /api/gc
# ============================================================

@alignment_bp.route("/api/gc", methods=["POST"])
def api_gc():
    """
    GC content of a sequence.

    Request JSON:
        seq  string  required

    Response JSON:
        gc_content  float  0.0 – 100.0
    """
    data, err = _require_json("seq")
    if err:
        return err

    seq = data["seq"].strip().upper()
    if not bridge.validate(seq):
        return _bad_request("Invalid DNA sequence")

    return jsonify({"gc_content": round(bridge.gc_content(seq), 2)}), 200


# ============================================================
# POST /api/validate
# ============================================================

@alignment_bp.route("/api/validate", methods=["POST"])
def api_validate():
    """
    Validate a DNA sequence.

    Request JSON:
        seq  string  required

    Response JSON:
        valid   bool
        length  int
    """
    data = request.get_json(silent=True)
    if data is None or "seq" not in data:
        return _bad_request("Missing field: 'seq'")

    seq   = (data["seq"] or "").strip().upper()
    valid = bridge.validate(seq)

    return jsonify({"valid": valid, "length": len(seq)}), 200


# ============================================================
# POST /api/reverse_complement
# ============================================================

@alignment_bp.route("/api/reverse_complement", methods=["POST"])
def api_reverse_complement():
    """
    Reverse complement of a sequence.

    Request JSON:
        seq  string  required

    Response JSON:
        result  string
    """
    data, err = _require_json("seq")
    if err:
        return err

    seq    = data["seq"].strip().upper()
    result = bridge.reverse_complement(seq)

    if result is None:
        return _bad_request("Invalid DNA sequence")

    return jsonify({"result": result}), 200