"""Edit-distance helpers (Myers bit-vector)."""
from __future__ import annotations


def myers(pattern: str, text: str, band: int | None = None) -> int:
    """Return the Levenshtein distance using Myers' bit-vector algorithm."""
    pattern = pattern.upper()
    text = text.upper()

    m = len(pattern)
    if m == 0:
        return len(text)
    if m > 64:
        raise ValueError("Reference Myers implementation supports pattern length <= 64.")

    peq: dict[str, int] = {}
    for ch in set(pattern):
        peq[ch] = 0
    for i, ch in enumerate(pattern):
        peq[ch] = peq.get(ch, 0) | (1 << i)

    pv = ~0
    mv = 0
    score = m
    high_bit = 1 << (m - 1)

    for idx, ch in enumerate(text):
        eq = peq.get(ch, 0)
        xv = eq | mv
        xh = (((eq & pv) + pv) ^ pv) | eq
        ph = mv | ~(xh | pv)
        mh = pv & xh

        if ph & high_bit:
            score += 1
        elif mh & high_bit:
            score -= 1

        pv = (mh << 1) | ~(xh | (ph << 1))
        mv = (ph << 1) & xh

        if band is not None:
            lower = score - max(0, idx - m + 1)
            if lower > band:
                return band + 1

    return score


def myers_search(pattern: str, text: str, max_distance: int) -> list[dict[str, int]]:
    """Return approximate-match hits (start/end/score) using Myers' algorithm."""
    pattern = pattern.upper()
    text = text.upper()

    m = len(pattern)
    if m == 0:
        return [{"start": i, "end": i, "score": 0} for i in range(len(text) + 1)]
    if m > 64:
        raise ValueError("Reference Myers implementation supports pattern length <= 64.")

    peq: dict[str, int] = {}
    for ch in set(pattern):
        peq[ch] = 0
    for i, ch in enumerate(pattern):
        peq[ch] = peq.get(ch, 0) | (1 << i)

    pv = ~0
    mv = 0
    score = m
    high_bit = 1 << (m - 1)
    hits: list[dict[str, int]] = []

    for idx, ch in enumerate(text):
        eq = peq.get(ch, 0)
        xv = eq | mv
        xh = (((eq & pv) + pv) ^ pv) | eq
        ph = mv | ~(xh | pv)
        mh = pv & xh

        if ph & high_bit:
            score += 1
        elif mh & high_bit:
            score -= 1

        pv = (mh << 1) | ~(xh | (ph << 1))
        mv = (ph << 1) & xh

        if score <= max_distance:
            start = idx - m + 1
            if start >= 0:
                hits.append({"start": start, "end": idx + 1, "score": score})

    return hits
