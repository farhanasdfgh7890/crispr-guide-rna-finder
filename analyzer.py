from guide_extractor import GuideRNA

def analyze(guides, seq="", use_ml=True):
    if use_ml and guides:
        try:
            from ml_predictor import predict_batch
            efficiencies = predict_batch([g.guide for g in guides])
            for g, eff in zip(guides, efficiencies):
                g.ml_efficiency = eff
        except Exception:
            pass
    for g in guides:
        g.gc_percent = _gc_content(g.guide)
        g.score = _score(g.guide)
        if seq:
            hits = _count_off_targets(g.guide, seq)
            g.off_target_hits = hits
            g.off_target_risk = _risk_label(hits)
    guides.sort(key=lambda g: 0.6 * g.score + 0.4 * g.ml_efficiency * 80, reverse=True)
    return guides

def filter_guides(guides, min_gc=0.0, max_gc=100.0, min_score=0, exclude_high_offtarget=False):
    result = [g for g in guides if min_gc <= g.gc_percent <= max_gc and g.score >= min_score]
    if exclude_high_offtarget:
        result = [g for g in result if g.off_target_risk != "High"]
    return result

def quality_label(score):
    if score >= 60: return "High"
    if score >= 35: return "Medium"
    return "Low"

def _gc_content(seq):
    if not seq: return 0.0
    gc = sum(1 for nt in seq if nt in "GC")
    return round((gc / len(seq)) * 100, 1)

def _score(guide):
    score = 0
    gc = _gc_content(guide)
    if 40 <= gc <= 70:   score += 30
    elif 30 <= gc <= 80: score += 15
    if "TTTT"  not in guide: score += 15
    if "AAAAA" not in guide: score += 10
    if not _has_repeated_kmer(guide): score += 10
    if guide[-1] == "G": score += 10
    if guide[0]  == "G": score += 5
    return score

def _has_repeated_kmer(seq, k=4):
    seen = set()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if kmer in seen: return True
        seen.add(kmer)
    return False

def _count_off_targets(guide, seq, max_mismatches=3):
    guide_len = len(guide)
    seed_start = guide_len - 12
    hits = 0
    for i in range(len(seq) - guide_len + 1):
        window = seq[i:i+guide_len]
        if window == guide: continue
        cost = sum(2 if j >= seed_start else 1 for j,(a,b) in enumerate(zip(guide,window)) if a != b)
        if cost <= max_mismatches: hits += 1
    return hits

def _risk_label(hits):
    if hits == 0:  return "Low"
    if hits <= 3:  return "Medium"
    return "High"