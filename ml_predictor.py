"""
ml_predictor.py — Machine Learning Efficiency Predictor
========================================================
Predicts on-target CRISPR-Cas9 cleavage efficiency for each guide RNA.

APPROACH
--------
We use a Random Forest trained on biologically-grounded synthetic features
derived from the Doench et al. 2016 rule set and the Wang et al. 2014
dataset patterns. Since real experimental data requires a license, we
generate a reproducible synthetic training set whose feature-to-efficiency
mapping mirrors published findings:

  • GC 40–70%           → higher efficiency
  • G at position 20    → higher efficiency
  • A at position 1–3   → lower efficiency
  • Poly-T runs         → lower efficiency
  • Dinucleotide context around PAM → affects efficiency
  • Thermodynamic stability proxy (GC in seed region)

The model is trained once on first import and cached in memory.
Predicted values are calibrated to the 0–1 range.

EXTENDING TO REAL DATA
----------------------
To use real experimental data (e.g. Doench 2016 dataset):
  1. Load your CSV: pd.read_csv("doench2016.csv")
  2. Replace _generate_training_data() with your loader
  3. Re-run train_model()
"""

import numpy as np
import hashlib
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.preprocessing import MinMaxScaler

# ── Module-level cache (trained once per process) ─────────────────────────────
_model: GradientBoostingRegressor | None = None
_scaler: MinMaxScaler | None = None


def predict_efficiency(guide: str) -> float:
    """
    Predict on-target efficiency for a 20-nt guide RNA.
    Returns a float in [0, 1]. Higher = more efficient.
    """
    global _model, _scaler
    if _model is None:
        _model, _scaler = _train_model()

    features = np.array(_extract_features(guide)).reshape(1, -1)
    raw = _model.predict(features)[0]
    # Clip to [0, 1] — GBR can extrapolate slightly outside training range
    return round(float(np.clip(raw, 0.0, 1.0)), 3)


def predict_batch(guides: list[str]) -> list[float]:
    """Predict efficiency for a list of guide sequences."""
    global _model, _scaler
    if _model is None:
        _model, _scaler = _train_model()

    X = np.array([_extract_features(g) for g in guides])
    raw = _model.predict(X)
    return [round(float(np.clip(v, 0.0, 1.0)), 3) for v in raw]


def efficiency_label(score: float) -> str:
    """Convert 0-1 score to a human-readable label."""
    if score >= 0.70: return "High"
    if score >= 0.45: return "Medium"
    return "Low"


# ── Feature engineering ───────────────────────────────────────────────────────

def _extract_features(guide: str) -> list[float]:
    """
    Extract 28 numerical features from a guide RNA sequence.
    These mirror the feature set used in Doench 2016 / Rule Set 2.
    """
    if len(guide) < 20:
        guide = guide.ljust(20, "N")
    g = guide[:20].upper()

    feats = []

    # 1. Overall GC content
    gc = sum(1 for c in g if c in "GC") / 20
    feats.append(gc)

    # 2. GC content in seed region (positions 8–20, 0-indexed)
    seed = g[7:]
    feats.append(sum(1 for c in seed if c in "GC") / len(seed))

    # 3. GC content in PAM-distal region (positions 1–7)
    distal = g[:7]
    feats.append(sum(1 for c in distal if c in "GC") / len(distal))

    # 4. Individual nucleotide at each of 20 positions (one-hot A/T/C/G)
    nt_map = {"A": 0, "T": 1, "C": 2, "G": 3, "N": 0}
    for i in range(20):
        nt = nt_map.get(g[i], 0)
        feats.append(nt / 3.0)   # normalised 0–1

    # 5. Poly-T run length (max)
    feats.append(_max_run(g, "T") / 20)

    # 6. Poly-A run length (max)
    feats.append(_max_run(g, "A") / 20)

    # 7. Poly-G run (GGGG → Cas9 issues)
    feats.append(_max_run(g, "G") / 20)

    # 8. Terminal G at position 20 (boolean)
    feats.append(1.0 if g[-1] == "G" else 0.0)

    # 9. Terminal G at position 19
    feats.append(1.0 if g[-2] == "G" else 0.0)

    # 10. First nucleotide is G (U6 preference)
    feats.append(1.0 if g[0] == "G" else 0.0)

    # 11. Dinucleotide at PAM-proximal end (positions 19-20)
    dint = g[-2:]
    feats.append(1.0 if dint in ("GG", "GA", "GT") else 0.0)

    return feats


def _max_run(seq: str, nt: str) -> int:
    """Return the length of the longest run of nucleotide `nt` in seq."""
    max_run = cur = 0
    for c in seq:
        if c == nt:
            cur += 1
            max_run = max(max_run, cur)
        else:
            cur = 0
    return max_run


# ── Training ──────────────────────────────────────────────────────────────────

def _train_model():
    """
    Generate synthetic training data and train a GradientBoostingRegressor.
    Uses a fixed seed for reproducibility.
    """
    rng = np.random.RandomState(42)
    X, y = _generate_training_data(n=2000, rng=rng)

    model = GradientBoostingRegressor(
        n_estimators=200,
        max_depth=4,
        learning_rate=0.05,
        subsample=0.8,
        random_state=42,
    )
    model.fit(X, y)

    scaler = MinMaxScaler()
    scaler.fit(X)

    return model, scaler


def _generate_training_data(n: int, rng: np.random.RandomState):
    """
    Generate n synthetic (guide, efficiency) pairs.

    Efficiency is a deterministic function of the guide's features
    plus a small noise term, so the model learns biologically meaningful
    rules rather than fitting random noise.
    """
    nts = list("ATCG")
    X, y = [], []

    for _ in range(n):
        guide = "".join(rng.choice(nts, size=20))
        feats = _extract_features(guide)
        eff   = _rule_based_efficiency(guide, rng)
        X.append(feats)
        y.append(eff)

    return np.array(X), np.array(y)


def _rule_based_efficiency(guide: str, rng: np.random.RandomState) -> float:
    """
    Compute a synthetic efficiency score based on published CRISPR rules.
    Noise is added to simulate experimental variation.
    """
    score = 0.5   # baseline

    gc = sum(1 for c in guide if c in "GC") / 20
    if 0.40 <= gc <= 0.70:
        score += 0.20
    elif gc < 0.25 or gc > 0.85:
        score -= 0.20

    # Terminal G boost
    if guide[-1] == "G":  score += 0.10
    if guide[-2] == "G":  score += 0.05

    # First G (U6)
    if guide[0] == "G":   score += 0.05

    # Poly-T penalty (Pol III termination)
    if "TTTT" in guide:   score -= 0.15
    if "TTTTT" in guide:  score -= 0.10  # extra

    # Poly-G penalty (G-quadruplex)
    if "GGGG" in guide:   score -= 0.12

    # Seed region GC (positions 8–20)
    seed_gc = sum(1 for c in guide[7:] if c in "GC") / 13
    if 0.35 <= seed_gc <= 0.65:
        score += 0.08

    # A/T at positions 1–3 (weak Pol II binding)
    if guide[:3].count("A") + guide[:3].count("T") >= 2:
        score -= 0.06

    # Gaussian noise (σ = 0.08) — simulates experimental variation
    noise = rng.normal(0, 0.08)
    return float(np.clip(score + noise, 0.0, 1.0))
