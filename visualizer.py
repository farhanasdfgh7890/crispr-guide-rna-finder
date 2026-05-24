"""
Module 5 — Visualizer  (Phase 2 upgraded)
Terminal visualization — now includes off-target risk in output.
"""

from guide_extractor import GuideRNA
from analyzer import quality_label

RESET  = "\033[0m";  BOLD   = "\033[1m";  DIM = "\033[2m"
BG_TEAL   = "\033[48;5;30m";  BG_PURPLE = "\033[48;5;55m"
FG_WHITE  = "\033[97m";       FG_BLACK  = "\033[30m"
FG_GREEN  = "\033[32m";  FG_YELLOW = "\033[33m";  FG_RED = "\033[31m"
FG_CYAN   = "\033[36m";  FG_BLUE   = "\033[34m"


def print_results(guides: list[GuideRNA], seq_length: int) -> None:
    print()
    _print_header("  CRISPR Guide RNA Analysis Results  ")
    print(f"  Sequence length : {seq_length} bp")
    print(f"  Guides found    : {len(guides)}")
    fwd = sum(1 for g in guides if g.strand == "+")
    print(f"  Forward / Rev   : {fwd} / {len(guides)-fwd}")
    print()
    for i, g in enumerate(guides, 1):
        _print_guide(i, g)
    _print_legend()


def visualize_sequence(seq: str, guides: list[GuideRNA], max_display: int = 120) -> None:
    display = seq[:max_display]
    n = len(display)
    marker = [""] * n
    for g in [g for g in guides if g.strand == "+"][:5]:
        start = g.position - 1
        for i in range(start, min(start + len(g.guide), n)): marker[i] = "guide"
        for i in range(start + len(g.guide), min(start + len(g.guide) + 3, n)): marker[i] = "pam"

    print()
    _print_header(f"  Sequence visualization (first {min(max_display,len(seq))} bp)  ")
    for row in range(0, n, 60):
        print(f"{row+1:>5} ", end="")
        for i in range(row, min(row+60, n)):
            nt = display[i]
            if   marker[i] == "guide": print(f"{BG_TEAL}{FG_WHITE}{nt}{RESET}", end="")
            elif marker[i] == "pam":   print(f"{BG_PURPLE}{FG_WHITE}{nt}{RESET}", end="")
            else:                      print(f"{DIM}{nt}{RESET}", end="")
        print()
    print()
    _print_legend()


def print_position_map(guides: list[GuideRNA], seq_length: int, width: int = 60) -> None:
    print()
    _print_header("  Guide RNA position map  ")
    print(f"  0{'':<{width-4}}{seq_length}")
    print(f"  {'─'*width}")
    for i, g in enumerate(guides[:10], 1):
        frac   = (g.position - 1) / max(seq_length, 1)
        bar_pos = int(frac * width)
        bar_len = max(int(len(g.guide) / seq_length * width), 1)
        color   = FG_CYAN if g.strand == "+" else FG_YELLOW
        line    = [" "] * width
        for j in range(bar_pos, min(bar_pos+bar_len, width)): line[j] = "█"
        risk_col = {"Low": FG_GREEN, "Medium": FG_YELLOW, "High": FG_RED}.get(g.off_target_risk, "")
        print(f"  {color}{''.join(line)}{RESET}  "
              f"#{i} pos {g.position}  "
              f"OT: {risk_col}{g.off_target_risk}{RESET}")
    print()


# ── internal helpers ──────────────────────────────────────────────────────────

def _print_guide(idx: int, g: GuideRNA) -> None:
    ql    = quality_label(g.score)
    qcol  = FG_GREEN if ql=="High" else (FG_YELLOW if ql=="Medium" else FG_RED)
    otcol = {"Low": FG_GREEN, "Medium": FG_YELLOW, "High": FG_RED}.get(g.off_target_risk, "")
    strand = "(+)" if g.strand == "+" else "(−)"

    print(f"  {BOLD}#{idx:>3}{RESET}  {FG_CYAN}{g.guide}{RESET}  "
          f"{BG_PURPLE}{FG_WHITE} {g.pam} {RESET}")
    print(f"       Position    : {g.position:<8}  Strand : {strand}")
    print(f"       GC%%         : {g.gc_percent:<6.1f}    Score  : "
          f"{qcol}{g.score:>2}/80  [{ql}]{RESET}")
    print(f"       Off-target  : {otcol}{g.off_target_hits} hit(s)  [{g.off_target_risk} risk]{RESET}")
    print()


def _print_header(title: str) -> None:
    w = len(title) + 2
    print(f"  {'═'*w}")
    print(f"  {BOLD}{title}{RESET}")
    print(f"  {'═'*w}")
    print()


def _print_legend() -> None:
    print(f"  Legend: "
          f"{BG_TEAL}{FG_WHITE} guide {RESET} fwd guide  "
          f"{BG_PURPLE}{FG_WHITE} PAM {RESET} NGG site  "
          f"{FG_GREEN}■{RESET} High  {FG_YELLOW}■{RESET} Medium  {FG_RED}■{RESET} Low  "
          f"OT = off-target risk")
    print()
