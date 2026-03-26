from flask import Blueprint, request, jsonify
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import c_bridge as bridge

suffix_tree_bp = Blueprint("suffix_tree", __name__)

pathogen_trees = {}

def init_pathogens():
    # Load the 3 pathogens if available
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    paths = {
        "Influenza A (H1N1)": "data/pathogens/influenza_h1n1.fasta",
        "SARS-CoV-2": "data/pathogens/sars_cov2.fasta",
        "E. coli K12": "data/pathogens/ecoli_k12.fasta"
    }
    
    for name, rel_path in paths.items():
        fp = os.path.join(base_dir, rel_path.replace("/", os.sep))
        if os.path.exists(fp):
            with open(fp, "r") as f:
                lines = f.readlines()
                seq = "".join(l.strip() for l in lines if not l.startswith(">")).upper()
                tree_ptr = bridge.suffix_tree_build(seq)
                if tree_ptr:
                    pathogen_trees[name] = {"ptr": tree_ptr, "length": len(seq), "seq": seq}
                    print(f"[SuffixTree] Loaded {name} (length: {len(seq)})")

@suffix_tree_bp.route("/api/pathogen/identify", methods=["POST"])
def api_pathogen_identify():
    data = request.get_json()
    if not data or "sequence" not in data:
        return jsonify({"error": "Missing sequence"}), 400
        
    query = data["sequence"].strip().upper()
    if not bridge.validate(query):
        return jsonify({"error": "Invalid sequence"}), 400
        
    results = []
    
    for name, pdata in pathogen_trees.items():
        res = bridge.suffix_tree_search(pdata["ptr"], query)
        
        # Determine match score based on LCS
        # Pathogen matching: if we just search exactly it might fail, 
        # normally we should just return if there is a hit.
        # But we only compiled `search` into our bridge. 
        # If `suffix_tree_search` returns matches, it's a 100% exact substring match.
        
        if res and isinstance(res, dict) and res.get("count", 0) > 0:
            results.append({
                "pathogen": name,
                "confidence": 100.0,
                "matches": res["matches"]
            })
            
    # Sort by confidence
    results.sort(key=lambda x: x["confidence"], reverse=True)
    
    return jsonify({"results": results}), 200

@suffix_tree_bp.route("/api/suffix/search", methods=["POST"])
def api_suffix_search():
    data = request.get_json()
    if not data or "sequence" not in data or "pattern" not in data:
        return jsonify({"error": "Missing sequence or pattern"}), 400
        
    seq = data["sequence"].strip().upper()
    pat = data["pattern"].strip().upper()
    
    st_ptr = bridge.suffix_tree_build(seq)
    res = bridge.suffix_tree_search(st_ptr, pat)
    bridge.suffix_tree_free(st_ptr)
    
    return jsonify({"results": res}), 200
