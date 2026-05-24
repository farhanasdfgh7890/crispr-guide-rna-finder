"""
app.py — Flask Web Server  (Phase 4: ML + deployment-ready)
"""

import os
import json
from flask import Flask, request, jsonify, render_template
from flask_cors import CORS

from validator       import validate, ValidationError
from pam_scanner     import scan
from guide_extractor import extract
from analyzer        import analyze, filter_guides, quality_label
from utils           import sequence_stats

app = Flask(__name__)
CORS(app)

# Pre-warm ML model at startup so first request isn't slow
try:
    from ml_predictor import predict_efficiency, efficiency_label
    predict_efficiency("ATCGATCGATCGATCGATCG")
    ML_AVAILABLE = True
    print("  ✓  ML efficiency model loaded")
except Exception as e:
    ML_AVAILABLE = False
    print(f"  ⚠  ML model unavailable: {e}")


# ── Routes ────────────────────────────────────────────────────────────────────

@app.route("/")
def index():
    return render_template("index.html")


@app.route("/health")
def health():
    return jsonify({"status": "ok", "version": "4.0", "ml_available": ML_AVAILABLE})


@app.route("/analyze", methods=["POST"])
def analyze_sequence():
    """
    Full analysis pipeline.
    Body JSON: { sequence, guide_len, strand, min_gc, max_gc,
                 min_score, exclude_high_offtarget, use_ml }
    """
    data = request.get_json(silent=True) or {}
    raw  = data.get("sequence", "").strip()
    if not raw:
        return jsonify({"success": False, "error": "No sequence provided."}), 400

    try:
        seq = validate(raw)
    except ValidationError as e:
        return jsonify({"success": False, "error": str(e)}), 400

    guide_len              = int(data.get("guide_len", 20))
    strand                 = data.get("strand", "both")
    min_gc                 = float(data.get("min_gc", 0))
    max_gc                 = float(data.get("max_gc", 100))
    min_score              = int(data.get("min_score", 0))
    exclude_high_offtarget = bool(data.get("exclude_high_offtarget", False))
    use_ml                 = bool(data.get("use_ml", True)) and ML_AVAILABLE

    pam_sites = scan(seq, guide_len=guide_len)
    if strand == "forward": pam_sites = [s for s in pam_sites if s.strand == "+"]
    elif strand == "reverse": pam_sites = [s for s in pam_sites if s.strand == "-"]

    guides = extract(seq, pam_sites, guide_len=guide_len)
    guides = analyze(guides, seq=seq, use_ml=use_ml)
    guides = filter_guides(guides, min_gc=min_gc, max_gc=max_gc,
                           min_score=min_score,
                           exclude_high_offtarget=exclude_high_offtarget)

    return jsonify({
        "success":      True,
        "ml_available": ML_AVAILABLE,
        "stats":        sequence_stats(seq),
        "guides":       [_guide_dict(i+1, g) for i, g in enumerate(guides)],
    })


@app.route("/predict", methods=["POST"])
def predict_single():
    """
    Quick ML prediction for a single guide.
    Body JSON: { "guide": "ATCGATCGATCGATCGATCG" }
    """
    if not ML_AVAILABLE:
        return jsonify({"success": False, "error": "ML model not available."}), 503

    data  = request.get_json(silent=True) or {}
    guide = data.get("guide", "").strip().upper()

    if len(guide) < 20:
        return jsonify({"success": False, "error": "Guide must be at least 20 nt."}), 400

    eff = predict_efficiency(guide[:20])
    return jsonify({
        "success":    True,
        "guide":      guide[:20],
        "efficiency": eff,
        "label":      efficiency_label(eff),
    })


@app.route("/upload-fasta", methods=["POST"])
def upload_fasta():
    if "file" not in request.files:
        return jsonify({"success": False, "error": "No file uploaded."}), 400

    f = request.files["file"]
    if not f.filename.endswith((".fa", ".fasta", ".fas", ".txt")):
        return jsonify({"success": False, "error": "Upload a .fa / .fasta file."}), 400

    content = f.read().decode("utf-8", errors="ignore")
    records = _parse_fasta_string(content)
    if not records:
        return jsonify({"success": False, "error": "No sequences found in file."}), 400

    guide_len = int(request.form.get("guide_len", 20))
    strand    = request.form.get("strand", "both")
    min_score = int(request.form.get("min_score", 0))
    use_ml    = request.form.get("use_ml", "true").lower() == "true" and ML_AVAILABLE

    results = []
    for header, raw_seq in records.items():
        try:
            seq    = validate(raw_seq)
            sites  = scan(seq, guide_len=guide_len)
            if strand == "forward": sites = [s for s in sites if s.strand == "+"]
            elif strand == "reverse": sites = [s for s in sites if s.strand == "-"]
            guides = extract(seq, sites, guide_len=guide_len)
            guides = analyze(guides, seq=seq, use_ml=use_ml)
            guides = filter_guides(guides, min_score=min_score)
            results.append({
                "header": header,
                "stats":  sequence_stats(seq),
                "guides": [_guide_dict(i+1, g) for i, g in enumerate(guides)],
            })
        except ValidationError as e:
            results.append({"header": header, "error": str(e), "guides": []})

    return jsonify({"success": True, "ml_available": ML_AVAILABLE, "records": results})


# ── helpers ───────────────────────────────────────────────────────────────────

def _guide_dict(rank, g):
    try:
        from ml_predictor import efficiency_label
        ml_label = efficiency_label(g.ml_efficiency)
    except Exception:
        ml_label = "N/A"

    return {
        "rank":             rank,
        "guide":            g.guide,
        "pam":              g.pam,
        "position":         g.position,
        "strand":           g.strand,
        "gc_percent":       g.gc_percent,
        "score":            g.score,
        "quality":          quality_label(g.score),
        "off_target_hits":  g.off_target_hits,
        "off_target_risk":  g.off_target_risk,
        "ml_efficiency":    g.ml_efficiency,
        "ml_label":         ml_label,
    }


def _parse_fasta_string(text: str) -> dict:
    records, current_header, current_seq = {}, None, []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if current_header:
                records[current_header] = "".join(current_seq)
            current_header, current_seq = line[1:], []
        else:
            current_seq.append(line)
    if current_header:
        records[current_header] = "".join(current_seq)
    return records


# ── run ───────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    port  = int(os.environ.get("PORT", 5000))
    debug = os.environ.get("FLASK_DEBUG", "true").lower() == "true"
    print(f"\n  CRISPR Guide RNA Finder v4.0")
    print(f"  Open: http://localhost:{port}\n")
    app.run(host="0.0.0.0", port=port, debug=debug)
