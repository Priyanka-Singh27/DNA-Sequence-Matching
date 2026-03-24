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

# ── blueprints ────────────────────────────────────────────
from routes.alignment import alignment_bp
# Future blueprints (uncomment as implemented):
# from routes.suffix_tree  import suffix_tree_bp
# from routes.trie         import trie_bp
# from routes.segment_tree import segment_tree_bp
# from routes.interval_tree import interval_tree_bp
# from routes.union_find   import union_find_bp

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
    # app.register_blueprint(suffix_tree_bp)
    # app.register_blueprint(trie_bp)
    # app.register_blueprint(segment_tree_bp)
    # app.register_blueprint(interval_tree_bp)
    # app.register_blueprint(union_find_bp)

    # ── health check ─────────────────────────────────────
    @app.route("/api/health", methods=["GET"])
    def health():
        return jsonify({
            "status": "ok",
            "service": "GenomeX API",
            "version": "0.1.0",
            "endpoints_active": [
                "POST /api/align",
                "POST /api/score",
                "POST /api/identity",
                "POST /api/gc",
                "POST /api/validate",
                "POST /api/reverse_complement",
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