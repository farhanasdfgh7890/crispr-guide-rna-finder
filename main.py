"""
main.py — CRISPR Guide RNA Finder  (Phase 2)
Entry point that wires together all modules including off-target prediction,
batch FASTA processing, and TSV/JSON export.

Usage:
    python main.py                            # interactive prompt
    python main.py -s ATCGATCG...            # inline sequence
    python main.py -f data/sample.fa         # single or multi-record FASTA
    python main.py -f data/sample.fa --tsv results/out.tsv --json results/out.json
    python main.py -f data/sample.fa --min-gc 40 --max-gc 70 --min-score 50 --top 10
"""

import argparse
import sys
import time

from validator      import validate, ValidationError
from pam_scanner    import scan
from guide_extractor import extract
from analyzer       import analyze, filter_guides
from visualizer     import print_results, visualize_sequence, print_position_map
from utils          import read_fasta, export_tsv, export_json, sequence_stats


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="CRISPR Guide RNA Finder — SpCas9 / NGG PAM",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python main.py -s ATCGATCGGGGATCGATCGGGG
  python main.py -f data/sample.fa --tsv results/out.tsv
  python main.py -f data/sample.fa --min-gc 40 --max-gc 70 --min-score 50
  python main.py -f data/sample.fa --exclude-high-offtarget --top 5
        """,
    )
    src = p.add_mutually_exclusive_group()
    src.add_argument("-s", "--sequence",  help="DNA sequence string")
    src.add_argument("-f", "--fasta",     help="Path to FASTA file (single or multi-record)")

    p.add_argument("--guide-len",  type=int,   default=20, choices=[17,18,19,20,21])
    p.add_argument("--strand",     default="both", choices=["both","forward","reverse"])
    p.add_argument("--min-gc",     type=float, default=0.0)
    p.add_argument("--max-gc",     type=float, default=100.0)
    p.add_argument("--min-score",  type=int,   default=0)
    p.add_argument("--top",        type=int,   default=None, help="Show only top N guides")
    p.add_argument("--exclude-high-offtarget", action="store_true",
                   help="Remove guides with High off-target risk")
    p.add_argument("--no-vis",     action="store_true", help="Skip sequence visualization")
    p.add_argument("--tsv",        help="Export results to TSV  (e.g. results/out.tsv)")
    p.add_argument("--json",       help="Export results to JSON (e.g. results/out.json)")
    return p.parse_args()


# ── pipeline ──────────────────────────────────────────────────────────────────

def run(raw_seq: str, args: argparse.Namespace, label: str = "sequence") -> list:
    """
    Full analysis pipeline for one sequence.
    Returns the final guide list (useful when called from other code).
    """
    t0 = time.perf_counter()

    # 1 ─ Validate
    print(f"\n  {'─'*54}")
    print(f"  Sequence : {label}")
    try:
        seq = validate(raw_seq)
    except ValidationError as e:
        print(f"  ✗  Validation failed: {e}")
        return []

    stats = sequence_stats(seq)
    print(f"  Length   : {stats['length']} bp   "
          f"GC {stats['gc_percent']}%   "
          f"A:{stats['A']}  T:{stats['T']}  C:{stats['C']}  G:{stats['G']}")

    # 2 ─ Scan PAM sites
    pam_sites = scan(seq, guide_len=args.guide_len)
    if args.strand == "forward":
        pam_sites = [s for s in pam_sites if s.strand == "+"]
    elif args.strand == "reverse":
        pam_sites = [s for s in pam_sites if s.strand == "-"]
    print(f"  PAM sites: {len(pam_sites)}  ({args.strand} strand{'s' if args.strand=='both' else ''})")

    # 3 ─ Extract guides
    guides = extract(seq, pam_sites, guide_len=args.guide_len)
    print(f"  Guides   : {len(guides)} unique guide RNAs extracted")
    if not guides:
        print("  No guides found — try a longer sequence.")
        return []

    # 4 ─ Analyze: GC + score + off-target
    print(f"  Scoring + off-target prediction …")
    guides = analyze(guides, seq=seq)          # pass seq for off-target scan
    guides = filter_guides(
        guides,
        min_gc=args.min_gc,
        max_gc=args.max_gc,
        min_score=args.min_score,
        exclude_high_offtarget=args.exclude_high_offtarget,
    )
    if args.top:
        guides = guides[:args.top]

    elapsed = time.perf_counter() - t0
    print(f"  Done in  : {elapsed:.2f}s   →  {len(guides)} guides after filters")

    # 5 ─ Report + visualize
    print_results(guides, seq_length=len(seq))
    if not args.no_vis:
        visualize_sequence(seq, guides)
        print_position_map(guides, seq_length=len(seq))

    # 6 ─ Export
    if args.tsv:
        export_tsv(guides, args.tsv)
    if args.json:
        meta = {
            "sequence_label": label,
            "length": len(seq),
            "guide_len": args.guide_len,
            "strand": args.strand,
        }
        export_json(guides, args.json, meta=meta)

    return guides


# ── entry point ───────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()

    if args.fasta:
        # ── BATCH FASTA MODE ─────────────────────────────────────────────────
        records = read_fasta(args.fasta)
        if not records:
            print("No sequences found in FASTA file.")
            sys.exit(1)

        print(f"\n  FASTA file : {args.fasta}")
        print(f"  Records    : {len(records)}")

        all_guides = []
        for header, seq in records.items():
            guides = run(seq, args, label=header)
            all_guides.extend(guides)

        if len(records) > 1:
            print(f"\n  ══  Batch summary: {len(all_guides)} total guides across {len(records)} sequences  ══")

    elif args.sequence:
        run(args.sequence, args)

    else:
        # ── INTERACTIVE MODE ──────────────────────────────────────────────────
        print("\n  ╔══════════════════════════════════════════╗")
        print("  ║   CRISPR Guide RNA Finder  v2.0          ║")
        print("  ║   SpCas9 · NGG PAM · off-target scoring  ║")
        print("  ╚══════════════════════════════════════════╝")
        print("\n  Paste a DNA sequence or FASTA block, then press Enter twice.")
        print("  (type 'quit' to exit)\n")

        lines = []
        while True:
            try:
                line = input()
            except (EOFError, KeyboardInterrupt):
                break
            if line.strip().lower() == "quit":
                break
            if line == "" and lines:
                break
            lines.append(line)

        raw = "\n".join(lines).strip()
        if raw:
            run(raw, args)
        else:
            print("  No input provided.")


if __name__ == "__main__":
    main()
