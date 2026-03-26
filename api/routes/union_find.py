from flask import Blueprint, request, jsonify
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import c_bridge as bridge

union_find_bp = Blueprint("union_find", __name__)

@union_find_bp.route("/api/evolution/cluster", methods=["POST"])
def api_evolution_cluster():
    data = request.get_json()
    if not data or "species" not in data or "matrix" not in data:
        return jsonify({"error": "Missing species list or similarity matrix"}), 400
        
    species = data["species"]
    matrix = data["matrix"]
    threshold = float(data.get("threshold", 95.0))
    
    if len(species) != len(matrix):
        return jsonify({"error": "Dimension mismatch"}), 400
        
    res = bridge.union_find_process(species, matrix, threshold)
    if not res:
        return jsonify({"error": "Failed to map dependencies in C-core backend"}), 500
        
    return jsonify({"results": res}), 200
