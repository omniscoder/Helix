"""McCaskill-style partition function and ensemble helpers."""
from __future__ import annotations

import math
from typing import Dict, List, Tuple

from .params import DEFAULTS, PAIRS


def partition_posteriors(seq: str, params: Dict[str, object] | None = None, beta: float | None = None) -> Dict[str, object]:
    cfg = {**DEFAULTS, **(params or {})}
    beta = beta if beta is not None else cfg["beta"]
    s = seq.upper().replace("T", "U")
    n = len(s)
    if n == 0:
        return {"Q": 1.0, "P": [], "p_unpaired": [], "entropy": []}

    Q = [[0.0 for _ in range(n)] for _ in range(n)]
    Qb = [[0.0 for _ in range(n)] for _ in range(n)]

    for i in range(n):
        Q[i][i] = 1.0

    hairpin_min = cfg["hairpin_min"]

    def hairpin_energy(i: int, j: int) -> float:
        loop = j - i - 1
        if loop < hairpin_min:
            return math.inf
        return cfg["hairpin_penalty"](loop)

    def stack_energy(i: int, j: int) -> float:
        if i + 1 >= j:
            return math.inf
        return cfg["stack_energy"](s[i], s[j], s[i + 1], s[j - 1])

    def boltz(E: float) -> float:
        if math.isinf(E):
            return 0.0
        return math.exp(-beta * E)

    for span in range(1, n):
        for i in range(0, n - span):
            j = i + span

            if (s[i], s[j]) in PAIRS:
                total = boltz(hairpin_energy(i, j))
                if i + 1 < j and Qb[i + 1][j - 1] > 0:
                    total += Qb[i + 1][j - 1] * boltz(stack_energy(i, j))
                for bulge in range(1, 3):
                    if i + bulge < j and Qb[i + bulge][j - 1] > 0:
                        total += Qb[i + bulge][j - 1] * boltz(cfg["bulge_penalty"](bulge))
                    if i + 1 < j - bulge and Qb[i + 1][j - bulge] > 0:
                        total += Qb[i + 1][j - bulge] * boltz(cfg["bulge_penalty"](bulge))
                if i + 2 < j - 1 and Qb[i + 2][j - 2] > 0:
                    total += Qb[i + 2][j - 2] * boltz(cfg["internal_penalty"](1, 1))
                Qb[i][j] = total
            else:
                Qb[i][j] = 0.0

            total_Q = Q[i + 1][j] if i + 1 <= j else 1.0
            total_Q += Q[i][j - 1] if i <= j - 1 else 1.0
            total_Q -= Q[i + 1][j - 1] if i + 1 <= j - 1 else 1.0

            pair_sum = 0.0
            for k in range(i, j):
                if Qb[k][j] == 0:
                    continue
                left = Q[i][k - 1] if i <= k - 1 else 1.0
                pair_sum += left * Qb[k][j]
            total_Q += pair_sum
            Q[i][j] = max(total_Q, 1e-12)

    Z = Q[0][n - 1]
    posterior = [[0.0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if Qb[i][j] == 0:
                continue
            left = Q[0][i - 1] if i - 1 >= 0 else 1.0
            right = Q[j + 1][n - 1] if j + 1 < n else 1.0
            posterior[i][j] = (left * Qb[i][j] * right) / Z
            posterior[j][i] = posterior[i][j]

    p_unpaired = []
    entropy = []
    for i in range(n):
        pair_sum = sum(posterior[i][j] for j in range(n))
        pu = max(0.0, 1.0 - pair_sum)
        p_unpaired.append(pu)
        H = 0.0
        if pu > 0:
            H -= pu * math.log(pu)
        for j in range(n):
            pij = posterior[i][j]
            if pij > 0:
                H -= pij * math.log(pij)
        entropy.append(H)

    return {
        "Q": Z,
        "P": posterior,
        "p_unpaired": p_unpaired,
        "entropy": entropy,
    }


def mea_structure(seq: str, posterior: List[List[float]], gamma: float = 1.0) -> Dict[str, object]:
    s = seq.upper().replace("T", "U")
    n = len(s)
    if n == 0:
        return {"dotbracket": "", "score": 0.0, "pairs": []}
    pu = [max(0.0, 1.0 - sum(posterior[i][j] for j in range(n))) for i in range(n)]
    S = [[0.0 for _ in range(n)] for _ in range(n)]
    B = [[None for _ in range(n)] for _ in range(n)]

    for span in range(1, n):
        for i in range(0, n - span):
            j = i + span
            best = S[i + 1][j] + pu[i]
            choice = ("i_unpaired", i + 1, j)
            if S[i][j - 1] + pu[j] > best:
                best = S[i][j - 1] + pu[j]
                choice = ("j_unpaired", i, j - 1)
            gain = 2 * gamma * posterior[i][j]
            if gain + (S[i + 1][j - 1] if i + 1 <= j - 1 else 0.0) > best:
                best = gain + (S[i + 1][j - 1] if i + 1 <= j - 1 else 0.0)
                choice = ("pair", i + 1, j - 1)
            for k in range(i, j):
                sc = S[i][k] + S[k + 1][j]
                if sc > best:
                    best = sc
                    choice = ("split", i, k, k + 1, j)
            S[i][j] = best
            B[i][j] = choice

    pairs: List[Tuple[int, int]] = []

    def traceback(i: int, j: int) -> None:
        if i >= j or B[i][j] is None:
            return
        tag = B[i][j][0]
        if tag in {"i_unpaired", "j_unpaired"}:
            _, a, b = B[i][j]
            traceback(a, b)
        elif tag == "pair":
            _, a, b = B[i][j]
            pairs.append((i, j))
            traceback(a, b)
        elif tag == "split":
            _, i1, k, k1, j1 = B[i][j]
            traceback(i1, k)
            traceback(k1, j1)

    traceback(0, n - 1)
    dotbracket = ["."] * n
    for i, j in pairs:
        dotbracket[i] = "("
        dotbracket[j] = ")"
    return {"dotbracket": "".join(dotbracket), "score": S[0][n - 1], "pairs": pairs}


def centroid_structure(seq: str, posterior: List[List[float]]) -> Dict[str, object]:
    return mea_structure(seq, posterior, gamma=1.0)
