"""
utils.py — Shared helpers used across modules.
"""

import os
import json
from datetime import datetime
from guide_extractor import GuideRNA
from analyzer import quality_label


# ── FASTA I/O ─────────────────────────────────────────────────────────────────

def read_fasta(path: str) -> dict[str, str]:
    """
    Parse a FASTA file and return a dict of {header: sequence}.
    Sequences are cleaned (uppercase, no whitespace).
    """
    records: dict[str, str] = {}
    current_header = None
    current_seq: list[str] = []

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if current_header is not None:
                    records[current_header] = "".join(current_seq).upper()
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line.replace(" ", ""))

    if current_header is not None:
        records[current_header] = "".join(current_seq).upper()

    return records


def write_fasta(path: str, records: dict[str, str], line_width: int = 60) -> None:
    """Write a dict of {header: sequence} to a FASTA file."""
    with open(path, "w") as fh:
        for header, seq in records.items():
            fh.write(f">{header}\n")
            for i in range(0, len(seq), line_width):
                fh.write(seq[i:i + line_width] + "\n")


# ── Result export ─────────────────────────────────────────────────────────────

def export_tsv(guides: list[GuideRNA], path: str) -> None:
    """Export guide RNA results to a TSV file."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as fh:
        fh.write("rank\tguide\tpam\tposition\tstrand\tgc_percent\tscore\tquality\n")
        for i, g in enumerate(guides, 1):
            fh.write(
                f"{i}\t{g.guide}\t{g.pam}\t{g.position}\t"
                f"{g.strand}\t{g.gc_percent}\t{g.score}\t{quality_label(g.score)}\n"
            )
    print(f"  Results saved → {path}")


def export_json(guides: list[GuideRNA], path: str, meta: dict | None = None) -> None:
    """Export guide RNA results to a JSON file."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    data = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "meta": meta or {},
        "guides": [
            {
                "rank": i,
                "guide": g.guide,
                "pam": g.pam,
                "position": g.position,
                "strand": g.strand,
                "gc_percent": g.gc_percent,
                "score": g.score,
                "quality": quality_label(g.score),
            }
            for i, g in enumerate(guides, 1)
        ],
    }
    with open(path, "w") as fh:
        json.dump(data, fh, indent=2)
    print(f"  Results saved → {path}")


# ── Sequence stats ────────────────────────────────────────────────────────────

def sequence_stats(seq: str) -> dict:
    """Return a dict of basic statistics for a DNA sequence."""
    counts = {nt: seq.count(nt) for nt in "ATCG"}
    total = len(seq)
    gc = counts["G"] + counts["C"]
    return {
        "length": total,
        "A": counts["A"],
        "T": counts["T"],
        "C": counts["C"],
        "G": counts["G"],
        "gc_percent": round((gc / total) * 100, 1) if total else 0.0,
        "at_percent": round(((counts["A"] + counts["T"]) / total) * 100, 1) if total else 0.0,
    }
