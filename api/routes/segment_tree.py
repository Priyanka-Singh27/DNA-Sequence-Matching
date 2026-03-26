from flask import Blueprint, request, jsonify
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import c_bridge as bridge

segment_tree_bp = Blueprint("segment_tree", __name__)

_loaded = False

def init_segment_tree(filepath: str):
    global _loaded
    if os.path.exists(filepath):
        res = bridge.segment_tree_init(filepath)
        if res > 0:
            _loaded = True
            print(f"[SegmentTree] Loaded {res} SNPs from {filepath}")
            return res
    return 0

@segment_tree_bp.route("/api/segment/query", methods=["POST"])
def api_segment_query():
    if not _loaded:
        return jsonify({"error": "Segment tree not loaded"}), 503
        
    data = request.get_json()
    if not data or "start" not in data or "end" not in data:
        return jsonify({"error": "Missing start or end"}), 400
        
    try:
        start = int(data["start"])
        end = int(data["end"])
    except ValueError:
        return jsonify({"error": "start and end must be integers"}), 400
        
    res = bridge.segment_tree_query(start, end)
    return jsonify({"results": res}), 200
