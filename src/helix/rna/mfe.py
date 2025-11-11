"""Zuker-style minimal MFE folding with traceback."""
from __future__ import annotations

from typing import Dict, List, Tuple

from .params import DEFAULTS, PAIRS

INF = 1e9


def mfe_dotbracket(seq: str, params: Dict[str, object] | None = None) -> Dict[str, object]:
    p = {**DEFAULTS, **(params or {})}
    s = seq.upper().replace("T", "U")
    n = len(s)
    if n == 0:
        return {"dotbracket": "", "energy": 0.0, "pairs": [], "trace": None}

    V = [[INF] * n for _ in range(n)]
    W = [[INF] * n for _ in range(n)]
    bt_V: List[List[Tuple]] = [[None] * n for _ in range(n)]
    bt_W: List[List[Tuple]] = [[None] * n for _ in range(n)]

    for i in range(n):
        W[i][i] = 0.0

    hairpin_min = p["hairpin_min"]

    def can_pair(i: int, j: int) -> bool:
        return (s[i], s[j]) in PAIRS

    def hairpin_energy(i: int, j: int) -> float:
        loop = j - i - 1
        if loop < hairpin_min:
            return INF
        return p["hairpin_penalty"](loop)

    def stack_energy(i: int, j: int) -> float:
        if i + 1 >= j:
            return INF
        return p["stack_energy"](s[i], s[j], s[i + 1], s[j - 1])

    def bulge_energy(size: int) -> float:
        return p["bulge_penalty"](size)

    def internal_energy(l1: int, l2: int) -> float:
        return p["internal_penalty"](l1, l2)

    span_range = range(1, n)
    for span in span_range:
        for i in range(0, n - span):
            j = i + span

            best = INF
            choice = None
            # unpaired i
            if W[i + 1][j] < best:
                best = W[i + 1][j]
                choice = ("unpaired_i", i + 1, j)
            # unpaired j
            if W[i][j - 1] < best:
                best = W[i][j - 1]
                choice = ("unpaired_j", i, j - 1)

            if can_pair(i, j):
                e = hairpin_energy(i, j)
                pair_choice = ("hairpin",)
                if i + 1 < j and V[i + 1][j - 1] < INF:
                    e_stack = V[i + 1][j - 1] + stack_energy(i, j)
                    if e_stack < e:
                        e = e_stack
                        pair_choice = ("stack", i + 1, j - 1)

                for bulge in range(1, 3):
                    if i + bulge < j and V[i + bulge][j - 1] < INF:
                        e_b = V[i + bulge][j - 1] + bulge_energy(bulge)
                        if e_b < e:
                            e = e_b
                            pair_choice = ("bulge_left", i + bulge, j - 1)
                    if i + 1 < j - bulge and V[i + 1][j - bulge] < INF:
                        e_b = V[i + 1][j - bulge] + bulge_energy(bulge)
                        if e_b < e:
                            e = e_b
                            pair_choice = ("bulge_right", i + 1, j - bulge)

                if i + 2 < j - 1 and V[i + 2][j - 2] < INF:
                    e_internal = V[i + 2][j - 2] + internal_energy(1, 1)
                    if e_internal < e:
                        e = e_internal
                        pair_choice = ("internal_11", i + 2, j - 2)

                V[i][j] = e
                bt_V[i][j] = pair_choice
                if e < best:
                    best = e
                    choice = ("pair", i, j)

            for k in range(i, j):
                e_split = W[i][k] + W[k + 1][j]
                if e_split < best:
                    best = e_split
                    choice = ("split", i, k, k + 1, j)

            W[i][j] = best
            bt_W[i][j] = choice

    pairs: List[Tuple[int, int]] = []

    def traceback_W(i: int, j: int) -> None:
        if i >= j or bt_W[i][j] is None:
            return
        tag = bt_W[i][j][0]
        if tag == "unpaired_i":
            _, ni, nj = bt_W[i][j]
            traceback_W(ni, nj)
        elif tag == "unpaired_j":
            _, ni, nj = bt_W[i][j]
            traceback_W(ni, nj)
        elif tag == "pair":
            _, a, b = bt_W[i][j]
            pairs.append((a, b))
            traceback_V(a, b)
        elif tag == "split":
            _, i1, k, k1, j1 = bt_W[i][j]
            traceback_W(i1, k)
            traceback_W(k1, j1)

    def traceback_V(i: int, j: int) -> None:
        if bt_V[i][j] is None:
            return
        tag = bt_V[i][j][0]
        if tag in {"stack", "bulge_left", "bulge_right", "internal_11"}:
            _, a, b = bt_V[i][j]
            traceback_V(a, b)
        # hairpin terminates

    traceback_W(0, n - 1)

    dotbracket = ["."] * n
    for i, j in pairs:
        dotbracket[i] = "("
        dotbracket[j] = ")"

    return {
        "dotbracket": "".join(dotbracket),
        "energy": W[0][n - 1],
        "pairs": sorted(pairs),
        "trace": None,
    }
