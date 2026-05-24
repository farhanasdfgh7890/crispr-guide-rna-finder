"""
Module 2 — PAM Scanner
Scans a DNA sequence for SpCas9 NGG PAM motifs on forward and reverse strands.
"""

from dataclasses import dataclass, field


@dataclass
class PAMSite:
    position: int       # 1-based position of the PAM in the original sequence
    pam: str            # 3-nt PAM string (e.g. 'AGG')
    strand: str         # '+' or '-'
    guide_start: int    # 1-based start of the 20-nt guide window


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA string."""
    comp = str.maketrans("ATCG", "TAGC")
    return seq.translate(comp)[::-1]


def scan(seq: str, guide_len: int = 20) -> list[PAMSite]:
    """
    Find all NGG PAM sites in *seq* on both strands.

    Returns a list of PAMSite objects sorted by position.
    Only sites where a full guide window fits upstream are returned.
    """
    sites = []
    sites.extend(_scan_forward(seq, guide_len))
    sites.extend(_scan_reverse(seq, guide_len))
    sites.sort(key=lambda s: (s.position, s.strand))
    return sites


# ── internal scanners ─────────────────────────────────────────────────────────

def _scan_forward(seq: str, guide_len: int) -> list[PAMSite]:
    """Scan the forward (+) strand for NGG."""
    sites = []
    for i in range(guide_len, len(seq) - 2):
        # PAM is at seq[i], seq[i+1], seq[i+2]
        if seq[i + 1] == "G" and seq[i + 2] == "G":
            pam = seq[i:i + 3]
            guide_start = i - guide_len + 1   # 1-based
            sites.append(PAMSite(
                position=i + 1,               # 1-based PAM position
                pam=pam,
                strand="+",
                guide_start=guide_start,
            ))
    return sites


def _scan_reverse(seq: str, guide_len: int) -> list[PAMSite]:
    """
    Scan the reverse (–) strand for NGG by looking for NCC on the forward
    strand (the complement of NGG read 3'→5').
    """
    sites = []
    rc = reverse_complement(seq)
    n = len(seq)

    for i in range(guide_len, len(rc) - 2):
        if rc[i + 1] == "G" and rc[i + 2] == "G":
            pam_rc = rc[i:i + 3]
            # Map rc position back to original sequence coordinate (1-based)
            orig_pam_pos = n - (i + 3) + 1
            guide_start_rc = i - guide_len + 1
            orig_guide_start = n - (guide_start_rc + guide_len - 1)

            sites.append(PAMSite(
                position=orig_pam_pos,
                pam=reverse_complement(pam_rc),   # PAM as it appears on rev strand
                strand="-",
                guide_start=orig_guide_start,
            ))
    return sites
