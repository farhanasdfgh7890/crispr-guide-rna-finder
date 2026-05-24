# CRISPR Guide RNA Finder

A bioinformatics tool that analyzes DNA sequences and identifies SpCas9 guide RNA
target sites using PAM (NGG) detection and scoring.

---

## Project Structure

```
crispr-guide-rna-finder/
├── main.py            # Entry point — wires all modules together
├── validator.py       # Module 1 — DNA sequence validation
├── pam_scanner.py     # Module 2 — NGG PAM motif detection
├── guide_extractor.py # Module 3 — Guide RNA extraction & deduplication
├── analyzer.py        # Module 4 — GC content, scoring, filtering
├── visualizer.py      # Module 5 — Terminal visualization
├── utils.py           # FASTA I/O, TSV/JSON export, sequence stats
│
├── data/
│   └── sample.fa      # Example FASTA sequences
└── results/           # Output directory (TSV / JSON exports)
```

---

## Requirements

- Python 3.10+ (uses `str | None` union syntax)
- No external libraries required — pure Python standard library only

---

## Usage

### Interactive mode
```bash
python main.py
```
Paste a DNA sequence or FASTA block, press Enter twice to analyze.

### Inline sequence
```bash
python main.py -s ATCGATCGTGCAGGCTACGGTAGCTATCGATCGATCGGGG
```

### FASTA file
```bash
python main.py -f data/sample.fa
```

### With filters and export
```bash
python main.py -f data/sample.fa \
  --guide-len 20 \
  --strand both \
  --min-gc 40 \
  --max-gc 70 \
  --min-score 50 \
  --top 10 \
  --tsv results/guides.tsv \
  --json results/guides.json
```

---

## CLI Options

| Option | Default | Description |
|---|---|---|
| `-s`, `--sequence` | — | Inline DNA string |
| `-f`, `--fasta` | — | Path to FASTA file |
| `--guide-len` | 20 | Guide RNA length (17–21) |
| `--strand` | both | `both` / `forward` / `reverse` |
| `--min-gc` | 0 | Minimum GC% filter |
| `--max-gc` | 100 | Maximum GC% filter |
| `--min-score` | 0 | Minimum quality score (0–80) |
| `--top` | all | Show only top N guides |
| `--no-vis` | off | Skip sequence visualization |
| `--tsv` | — | Export results to TSV |
| `--json` | — | Export results to JSON |

---

## Modules

### Module 1 — validator.py
- Strips FASTA headers
- Removes whitespace and digits
- Converts to uppercase
- Rejects non-ACGT characters
- Enforces minimum length (23 bp)

### Module 2 — pam_scanner.py
- Scans forward strand for NGG motifs
- Scans reverse strand via NCC detection on the forward strand (reverse complement of NGG)
- Returns `PAMSite` objects with position, PAM sequence, and strand

### Module 3 — guide_extractor.py
- Extracts 20 nt upstream of each PAM site
- Handles both strands correctly
- Deduplicates by guide sequence

### Module 4 — analyzer.py
- Calculates GC% using the formula: `(G + C) / Total × 100`
- Scores each guide out of 80 points:
  - +30 GC 40–70% (optimal)
  - +15 GC 30–40% or 70–80%
  - +15 No poly-T run of 4+
  - +10 No poly-A run of 5+
  - +10 No repeated 4-mer
  - +10 Terminal G (position 20)
  - +5 Starts with G
- Filters by GC range and minimum score

### Module 5 — visualizer.py
- ANSI-colored terminal output
- Highlights guide regions (teal) and PAM sites (purple)
- ASCII position map across the full sequence
- Quality labels: High (≥60) / Medium (≥35) / Low (<35)

---

## Output Format

```
Guide RNA : TCGTAGGCTAACCGGATCGA
PAM       : AGG
Position  : 120
GC%       : 55.0%
Score     : 65/80  [High]
Strand    : (+)
```

---

## Learning Concepts

**Biology**: DNA structure, CRISPR-Cas9 mechanism, guide RNA, PAM motifs, genomic targeting

**Computer Science**: String pattern matching, sequence algorithms, modular software design,
data filtering and ranking, file I/O

**Bioinformatics**: Computational genomics, FASTA format, reverse complement, GC content analysis
