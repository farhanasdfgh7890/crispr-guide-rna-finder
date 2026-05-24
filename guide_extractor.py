from dataclasses import dataclass

@dataclass
class GuideRNA:
    guide: str
    pam: str
    position: int
    strand: str
    gc_percent: float = 0.0
    score: int = 0
    off_target_hits: int = 0
    off_target_risk: str = "Unknown"
    ml_efficiency: float = 0.0

def reverse_complement(seq):
    comp = str.maketrans("ATCG", "TAGC")
    return seq.translate(comp)[::-1]

def extract(seq, pam_sites, guide_len=20):
    guides = []
    seen = set()
    for site in pam_sites:
        guide_seq = _extract_guide(seq, site, guide_len)
        if guide_seq is None:
            continue
        if guide_seq in seen:
            continue
        seen.add(guide_seq)
        guides.append(GuideRNA(
            guide=guide_seq,
            pam=site.pam,
            position=site.guide_start,
            strand=site.strand,
        ))
    return guides

def _extract_guide(seq, site, guide_len):
    if site.strand == "+":
        pam_idx = site.position - 1
        start = pam_idx - guide_len
        if start < 0:
            return None
        return seq[start:pam_idx]
    else:
        rc = reverse_complement(seq)
        n = len(seq)
        rc_pam_idx = n - site.position
        guide_start_rc = rc_pam_idx - guide_len
        if guide_start_rc < 0:
            return None
        return rc[guide_start_rc:rc_pam_idx]