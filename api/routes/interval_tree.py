from flask import Blueprint, request, jsonify
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import c_bridge as bridge

interval_tree_bp = Blueprint("interval_tree", __name__)

_loaded = False

def init_interval_tree(filepath: str):
    global _loaded
    if os.path.exists(filepath):
        res = bridge.interval_tree_init(filepath)
        if res > 0:
            _loaded = True
            print(f"[IntervalTree] Loaded {res} genes from {filepath}")
            return res
    return 0

@interval_tree_bp.route("/api/interval/query", methods=["POST"])
def api_interval_query():
    if not _loaded:
        return jsonify({"error": "Interval tree not loaded"}), 503
        
    data = request.get_json()
    if not data or "start" not in data or "end" not in data:
        return jsonify({"error": "Missing start or end"}), 400
        
    try:
        start = int(data["start"])
        end = int(data["end"])
    except ValueError:
        return jsonify({"error": "start and end must be integers"}), 400
        
    res = bridge.interval_tree_query(start, end)
    return jsonify({"results": res}), 200
