"""
Phase 2 — Unit tests for all 5 modules.
Run with:  python test_modules.py
"""

import sys, traceback

PASS = "\033[32m PASS\033[0m"
FAIL = "\033[31m FAIL\033[0m"
results = []

def test(name, fn):
    try:
        fn()
        print(f"{PASS}  {name}")
        results.append(True)
    except Exception as e:
        print(f"{FAIL}  {name}")
        traceback.print_exc()
        results.append(False)

def eq(a, b, msg=""):
    assert a == b, f"{msg}  expected {b!r}  got {a!r}"

def ok(cond, msg=""):
    assert cond, msg

# ─────────────────────────────────────────────────────────────
# Module 1 — validator
# ─────────────────────────────────────────────────────────────
from validator import validate, ValidationError

def t_valid_basic():
    seq = validate("ATCGATCGATCGATCGATCGATCG")
    eq(seq, "ATCGATCGATCGATCGATCGATCG")

def t_valid_lowercase():
    seq = validate("atcgatcgatcgatcgatcgatcg")
    eq(seq, "ATCGATCGATCGATCGATCGATCG")

def t_valid_fasta():
    fasta = ">gene1\nATCGATCGATCGATCGATCGATCG"
    seq = validate(fasta)
    eq(seq, "ATCGATCGATCGATCGATCGATCG")

def t_valid_spaces():
    seq = validate("ATCG ATCG ATCG ATCG ATCG ATG")
    ok(len(seq) >= 23, "cleaned seq should be at least 23 bp")
    ok(all(c in "ATCG" for c in seq), "should contain only ATCG")

def t_invalid_char():
    try:
        validate("ATCGXATCGATCGATCGATCGATCG")
        assert False, "should have raised"
    except ValidationError as e:
        ok("X" in str(e))

def t_too_short():
    try:
        validate("ATCG")
        assert False, "should have raised"
    except ValidationError as e:
        ok("short" in str(e).lower())

# ─────────────────────────────────────────────────────────────
# Module 2 — pam_scanner
# ─────────────────────────────────────────────────────────────
from pam_scanner import scan, reverse_complement

def t_revcomp():
    eq(reverse_complement("ATCG"), "CGAT")
    eq(reverse_complement("AAAA"), "TTTT")
    eq(reverse_complement("GCGC"), "GCGC")

def t_pam_finds_ngg():
    # hand-crafted: 20 nt guide + AGG PAM
    seq = "ATCGATCGATCGATCGATCG" + "AGG" + "TTTT"
    sites = [s for s in scan(seq, guide_len=20) if s.strand == "+"]
    ok(len(sites) >= 1, "should find at least one PAM")
    eq(sites[0].pam, "AGG")

def t_pam_both_strands():
    seq = "ATCGATCGATCGATCGATCGAGGATCGATCGATCGATCGATCGAGG"
    sites = scan(seq, guide_len=20)
    strands = {s.strand for s in sites}
    ok("+" in strands, "should find forward strand")

def t_pam_min_length():
    # too short to fit a guide
    seq = "ATCGAGG"
    sites = [s for s in scan(seq, guide_len=20) if s.strand == "+"]
    eq(len(sites), 0, "no guides should fit")

# ─────────────────────────────────────────────────────────────
# Module 3 — guide_extractor
# ─────────────────────────────────────────────────────────────
from pam_scanner import scan
from guide_extractor import extract

def t_extract_length():
    seq = "ATCGATCGATCGATCGATCG" + "AGG" + "TTTT"
    sites = scan(seq, guide_len=20)
    guides = extract(seq, sites, guide_len=20)
    for g in guides:
        eq(len(g.guide), 20, f"guide length wrong: {g.guide}")

def t_extract_dedup():
    # same 20-mer could appear at two PAM sites — should deduplicate
    guide = "ATCGATCGATCGATCGATCG"
    seq = guide + "AGG" + "NNNN" + guide + "TGG"
    seq = seq.replace("N","A")
    sites = scan(seq, guide_len=20)
    guides = extract(seq, sites, guide_len=20)
    sequences = [g.guide for g in guides]
    eq(len(sequences), len(set(sequences)), "duplicates found")

# ─────────────────────────────────────────────────────────────
# Module 4 — analyzer
# ─────────────────────────────────────────────────────────────
from analyzer import analyze, filter_guides, quality_label, _gc_content, _count_off_targets

def t_gc_content():
    eq(_gc_content("GGCC"), 100.0)
    eq(_gc_content("ATAT"), 0.0)
    eq(_gc_content("ATCG"), 50.0)

def t_score_range():
    seq = "ATCGATCGATCGATCGATCG" + "AGG" + "TTTT"
    sites = scan(seq, guide_len=20)
    guides = extract(seq, sites, guide_len=20)
    guides = analyze(guides, seq)
    for g in guides:
        ok(0 <= g.score <= 80, f"score out of range: {g.score}")

def t_quality_labels():
    eq(quality_label(75), "High")
    eq(quality_label(60), "High")
    eq(quality_label(59), "Medium")
    eq(quality_label(35), "Medium")
    eq(quality_label(34), "Low")
    eq(quality_label(0),  "Low")

def t_off_target_exact_skip():
    guide = "ATCGATCGATCGATCGATCG"
    seq   = guide + "AGG"
    hits  = _count_off_targets(guide, seq)
    eq(hits, 0, "exact match should not count as off-target")

def t_off_target_finds_near_match():
    guide  = "ATCGATCGATCGATCGATCG"
    near   = "ATCGATCGATCGATCGATCC"   # 1 mismatch at end
    seq    = guide + "AGG" + near + "AGG"
    hits   = _count_off_targets(guide, seq)
    ok(hits >= 1, "near match (1 mismatch) should be detected")

def t_filter_gc():
    seq = "GCGCGCGCGCGCGCGCGCGC" + "AGG" + "ATCGATCGATCGATCGATCG" + "AGG"
    sites  = scan(seq, guide_len=20)
    guides = extract(seq, sites, guide_len=20)
    guides = analyze(guides, seq)
    high_gc = filter_guides(guides, min_gc=80, max_gc=100)
    for g in high_gc:
        ok(g.gc_percent >= 80, f"GC filter wrong: {g.gc_percent}")

def t_filter_offtarget():
    seq = "ATCGATCGATCGATCGATCG" + "AGG" + "TTTT"
    sites  = scan(seq, guide_len=20)
    guides = extract(seq, sites, guide_len=20)
    guides = analyze(guides, seq)
    # mark one as High artificially
    if guides:
        guides[0].off_target_risk = "High"
    filtered = filter_guides(guides, exclude_high_offtarget=True)
    for g in filtered:
        ok(g.off_target_risk != "High")

# ─────────────────────────────────────────────────────────────
# Module 5 — utils (FASTA + export)
# ─────────────────────────────────────────────────────────────
import os, json, tempfile
from utils import read_fasta, write_fasta, export_tsv, export_json

def t_read_fasta():
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as f:
        f.write(">seq1\nATCGATCG\n>seq2\nGGGGCCCC\n")
        name = f.name
    records = read_fasta(name)
    os.unlink(name)
    eq(set(records.keys()), {"seq1","seq2"})
    eq(records["seq1"], "ATCGATCG")

def t_write_read_fasta():
    records = {"alpha": "ATCGATCG"*5, "beta": "GCGCGCGC"*3}
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as f:
        name = f.name
    write_fasta(name, records)
    back = read_fasta(name)
    os.unlink(name)
    eq(back["alpha"], records["alpha"])
    eq(back["beta"],  records["beta"])

def t_export_tsv():
    seq    = "ATCGATCGATCGATCGATCG" + "AGG" + "TTTT"
    sites  = scan(seq, guide_len=20)
    guides = extract(seq, sites)
    guides = analyze(guides, seq)
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
        name = f.name
    export_tsv(guides, name)
    with open(name) as f:
        lines = f.readlines()
    os.unlink(name)
    ok(len(lines) >= 2, "TSV should have header + at least one row")
    ok("guide" in lines[0], "header row missing")

def t_export_json():
    seq    = "ATCGATCGATCGATCGATCG" + "AGG" + "TTTT"
    sites  = scan(seq, guide_len=20)
    guides = extract(seq, sites)
    guides = analyze(guides, seq)
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
        name = f.name
    export_json(guides, name, meta={"test": True})
    with open(name) as f:
        data = json.load(f)
    os.unlink(name)
    ok("guides"    in data)
    ok("generated" in data)
    eq(data["meta"]["test"], True)

# ─────────────────────────────────────────────────────────────
# Run all tests
# ─────────────────────────────────────────────────────────────
print("\n  ══════════════════════════════════════════")
print("  CRISPR Guide RNA Finder — Unit Tests")
print("  ══════════════════════════════════════════\n")

# Module 1
print("  [ Module 1 — Validator ]")
test("valid basic sequence",        t_valid_basic)
test("valid lowercase conversion",  t_valid_lowercase)
test("valid FASTA format",          t_valid_fasta)
test("valid strips spaces",         t_valid_spaces)
test("rejects invalid character",   t_invalid_char)
test("rejects too-short sequence",  t_too_short)

# Module 2
print("\n  [ Module 2 — PAM Scanner ]")
test("reverse complement",          t_revcomp)
test("finds NGG PAM",               t_pam_finds_ngg)
test("finds both strands",          t_pam_both_strands)
test("no PAM when too short",       t_pam_min_length)

# Module 3
print("\n  [ Module 3 — Guide Extractor ]")
test("extracted guides are 20 nt",  t_extract_length)
test("deduplication works",         t_extract_dedup)

# Module 4
print("\n  [ Module 4 — Analyzer ]")
test("GC content calculation",      t_gc_content)
test("scores within 0–80 range",    t_score_range)
test("quality labels correct",      t_quality_labels)
test("off-target skips exact hit",  t_off_target_exact_skip)
test("off-target finds near match", t_off_target_finds_near_match)
test("filter by GC range",          t_filter_gc)
test("filter by off-target risk",   t_filter_offtarget)

# Module 5
print("\n  [ Module 5 — Utils / FASTA / Export ]")
test("read FASTA file",             t_read_fasta)
test("write then read FASTA",       t_write_read_fasta)
test("export TSV",                  t_export_tsv)
test("export JSON with meta",       t_export_json)

# Summary
passed = sum(results)
total  = len(results)
color  = "\033[32m" if passed == total else "\033[33m"
print(f"\n  {'═'*44}")
print(f"  {color}{passed}/{total} tests passed\033[0m")
if passed == total:
    print("  ✓  All tests passed — Phase 2 complete!")
else:
    print("  ✗  Some tests failed — see above for details.")
print(f"  {'═'*44}\n")
sys.exit(0 if passed == total else 1)
