"""
app.py
------
GenomeX Flask API server.

Registers all route blueprints and starts the development server.

Run:
    cd api
    LIBDNA_PATH=../c-core/libdna.so python app.py

Production:
    gunicorn -w 4 -b 0.0.0.0:5000 app:app
"""

from flask import Flask, jsonify
from flask_cors import CORS
import os

# ── blueprints ────────────────────────────────────────────
from routes.alignment import alignment_bp
from routes.trie      import trie_bp, init_trie_from_file
from routes.suffix_tree   import suffix_tree_bp, init_pathogens
from routes.segment_tree  import segment_tree_bp, init_segment_tree
from routes.interval_tree import interval_tree_bp, init_interval_tree
from routes.union_find    import union_find_bp

# ============================================================
# APP FACTORY
# ============================================================

def create_app() -> Flask:
    app = Flask(__name__)

    # Allow requests from React dev server (localhost:3000)
    CORS(app, resources={r"/api/*": {"origins": [
        "http://localhost:3000",
        "http://127.0.0.1:3000",
        "http://localhost:5500",   # Live Server (VS Code)
        "null",                    # file:// origin for direct HTML open
    ]}})

    # ── register blueprints ──────────────────────────────
    app.register_blueprint(alignment_bp)
    app.register_blueprint(trie_bp)
    app.register_blueprint(suffix_tree_bp)
    app.register_blueprint(segment_tree_bp)
    app.register_blueprint(interval_tree_bp)
    app.register_blueprint(union_find_bp)

    # ── pre-load SNP trie at startup ─────────────────────
    _snp_tsv = os.path.normpath(
        os.path.join(os.path.dirname(__file__), "..", "data", "clinvar_snps.tsv")
    )
    init_trie_from_file(_snp_tsv)
    
    _gene_tsv = os.path.normpath(
        os.path.join(os.path.dirname(__file__), "..", "data", "gene_regions.tsv")
    )
    init_segment_tree(_snp_tsv)
    init_interval_tree(_gene_tsv)
    init_pathogens()

    # ── health check ─────────────────────────────────────
    @app.route("/api/health", methods=["GET"])
    def health():
        return jsonify({
            "status": "ok",
            "service": "GenomeX API",
            "version": "0.2.0",
            "endpoints_active": [
                "POST /api/align",
                "POST /api/score",
                "POST /api/identity",
                "POST /api/gc",
                "POST /api/validate",
                "POST /api/reverse_complement",
                "POST /api/snp/fuzzy",
                "POST /api/snp/exact",
                "POST /api/snp/load",
                "POST /api/mutations/scan",
                "GET  /api/snp/stats",
            ]
        }), 200

    # ── 404 handler ──────────────────────────────────────
    @app.errorhandler(404)
    def not_found(e):
        return jsonify({"error": "Endpoint not found"}), 404

    # ── 500 handler ──────────────────────────────────────
    @app.errorhandler(500)
    def server_error(e):
        return jsonify({"error": "Internal server error", "detail": str(e)}), 500

    return app


app = create_app()

if __name__ == "__main__":
    print("=" * 50)
    print("  GenomeX API  —  http://localhost:5000")
    print("  Health:  GET  /api/health")
    print("  Align:   POST /api/align")
    print("=" * 50)
    app.run(host="0.0.0.0", port=5000, debug=True)