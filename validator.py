"""
Module 1 — DNA Validator
Validates and cleans input DNA sequences, including FASTA format.
"""


class ValidationError(Exception):
    pass


def validate(raw: str, min_length: int = 23) -> str:
    """
    Validate and clean a raw DNA string (or FASTA block).

    Returns the cleaned uppercase sequence string on success.
    Raises ValidationError with a descriptive message on failure.
    """
    if not raw or not raw.strip():
        raise ValidationError("Empty input — please provide a DNA sequence.")

    seq = _strip_fasta(raw)
    seq = _clean(seq)
    _check_characters(seq)
    _check_length(seq, min_length)

    return seq


# ── helpers ──────────────────────────────────────────────────────────────────

def _strip_fasta(text: str) -> str:
    """Remove FASTA header lines (lines starting with '>')."""
    lines = text.splitlines()
    seq_lines = [l for l in lines if not l.startswith(">")]
    return "".join(seq_lines)


def _clean(seq: str) -> str:
    """Remove whitespace, digits (position numbers), and convert to uppercase."""
    import re
    seq = re.sub(r"[\s\d]", "", seq)
    return seq.upper()


def _check_characters(seq: str) -> None:
    """Raise ValidationError if any character is not A/T/C/G."""
    import re
    bad = re.findall(r"[^ATCG]", seq)
    if bad:
        unique = sorted(set(bad))
        raise ValidationError(
            f"Invalid character(s) found: {', '.join(repr(c) for c in unique)}. "
            "Only A, T, C, G are allowed."
        )


def _check_length(seq: str, min_length: int) -> None:
    """Raise ValidationError if sequence is too short."""
    if len(seq) < min_length:
        raise ValidationError(
            f"Sequence is too short ({len(seq)} bp). "
            f"Need at least {min_length} bp to find a guide + PAM."
        )
