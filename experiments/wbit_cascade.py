"""§126 -- two W-bit cells in series: does CRNL's near-miss mechanism hold in W-bit's own dynamics?

Cell 1: W-bit R = 2 cell, self-coupling a = 1.5, T = 0.5, dt = 0.5, no bias (bistable, saddle at p = 1/2).
Cell 2: the same element driven by cell 1 through Theta_21 = c I and a bias b2: drive d = b2 + c (2 p1 - 1)
(W-bit updates synchronously, so cell 2 sees cell 1's previous state). Both cells get their own N(0, 2 sigma^2)
logit noise, exactly as `AnalogWbitNetwork.step` (validated in §125, W1).

Cell 2 is bistable while d > -d_c (d_c = 0.20755). Its EDGE e2 is the least p1 that keeps a high state:
b2 + c (2 e2 - 1) = -d_c. Built like CRNL's Hill-threshold lever (§124.9): cell 2's drive at a HEALTHY input
(p1 at its high well 0.9293) is fixed at h = 0.1 in both pairs -- identical cell-2 barrier there -- and only
the coupling c moves the edge:
    Z   c = 0.6713, edge 0.70  (above cell 1's saddle: a "zone" exists)
    N   c = 0.2445, edge 0.30  (below it: no zone)

Observables, exact (2-D transfer operator on an M x M grid, product kernel, direct solve):
started at both cells' high wells, exits are a DIP (p1 < 1/2) or a committed FALL (p2 reaches its low well,
p2 < p2_low + 0.02, with p1 never below 1/2). q_fall, q_dip = exit probabilities; R = q_fall / q_dip, the
W-bit analogue of CRNL's never-dipped falls per dip. Large parameter: 1/sigma^2.

PREDICTIONS, WRITTEN BEFORE RUNNING.
  C1  GRID. R at sigma = 0.30 changes by < 5% from M = 120 to M = 160 in both pairs (else the trend is not read).
  C2  THE MECHANISM'S SIGN, the transfer of §124.8 / §124.9: as sigma decreases over 0.34, 0.30, 0.27, 0.25,
      R INCREASES in Z and DECREASES in N.
  C3  PER-CELL RELIABILITY DOES NOT COMPOSE. In Z at the smallest sigma, cell 2's fall probability per cell-1 dip
      (R) exceeds what cell 2's own healthy-input escape gives -- its 1-D flip rate with input frozen at p1 = 0.9293,
      divided by cell 1's dip rate -- by more than 10x.

OUTCOMES (appended after running; predictions above left as written). C1 HELD: <= 0.6% from M = 120 to 160
(M = 200 exceeds SuperLU's memory here). C2 HELD: Z's R 0.523 / 0.576 / 0.604 / 0.615 rises, N's 0.119 /
0.080 / 0.048 / 0.029 falls, over sigma = 0.34 / 0.30 / 0.27 / 0.25. C3 HELD: in Z, R / per-cell estimate =
2.4 / 4.6 / 9.0 / 16.1 (N: 0.54 / 0.64 / 0.72 / 0.77). FINDINGS §126.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from experiments.wbit_cell import fixed_points, flip_rate, kernel

A, T, DT = 1.5, 0.5, 0.5
DC = 0.207546
H = 0.1
HI1 = 0.929284


def pair(edge):
    c = (DC + H) / (2 * HI1 - 1 - (2 * edge - 1))     # fixed healthy drive H, edge moved by c alone
    b2 = H - c * (2 * HI1 - 1)
    return c, b2


PAIRS = {"Z": pair(0.70), "N": pair(0.30)}


def joint_kernel(M, sigma, c, b2):
    """P[(i,j),(i',j')] = K1[i,i'] * K2^(i)[j,j'], cell 2 driven by cell 1's current state p1_i."""
    K1, x = kernel(M, A, T, DT, sigma, 0.0)
    K1 = K1.tocsr()
    blocks_r, blocks_c, blocks_v = [], [], []
    for i in range(M):
        K2, _ = kernel(M, A, T, DT, sigma, b2 + c * (2 * x[i] - 1))
        K2 = K2.tocoo()
        r1 = K1.indices[K1.indptr[i]:K1.indptr[i + 1]]
        v1 = K1.data[K1.indptr[i]:K1.indptr[i + 1]]
        for ip, w in zip(r1, v1):
            blocks_r.append(i * M + K2.row); blocks_c.append(ip * M + K2.col); blocks_v.append(w * K2.data)
    P = sp.csr_matrix((np.concatenate(blocks_v), (np.concatenate(blocks_r), np.concatenate(blocks_c))),
                      shape=(M * M, M * M))
    return P, x


def exits(M, sigma, c, b2):
    P, x = joint_kernel(M, sigma, c, b2)
    lo2 = fixed_points(A, T, DT, H)[0]
    I1, I2 = np.meshgrid(np.arange(M), np.arange(M), indexing="ij")
    p1, p2 = x[I1.ravel()], x[I2.ravel()]
    dip = p1 < 0.5
    fall = (~dip) & (p2 < lo2 + 0.02)
    alive = ~(dip | fall)
    ia = np.nonzero(alive)[0]
    Paa = P[ia][:, ia]
    Pa = P[ia]
    b_fall = np.asarray(Pa[:, np.nonzero(fall)[0]].sum(axis=1)).ravel()
    b_dip = np.asarray(Pa[:, np.nonzero(dip)[0]].sum(axis=1)).ravel()
    lu = spla.splu((sp.identity(len(ia)) - Paa).tocsc())
    qf, qd = lu.solve(b_fall), lu.solve(b_dip)
    hi2 = fixed_points(A, T, DT, H)[-1]
    k = int(np.argmin((p1[ia] - HI1) ** 2 + (p2[ia] - hi2) ** 2))
    return {"q_fall": float(qf[k]), "q_dip": float(qd[k]), "R": float(qf[k] / qd[k]),
            "mass_check": float(qf[k] + qd[k])}


def main():
    sigmas = [float(s) for s in sys.argv[1].split(",")] if len(sys.argv) > 1 else [0.34, 0.30, 0.27, 0.25]
    M = int(sys.argv[2]) if len(sys.argv) > 2 else 160
    out = Path("results/wbit_cascade.json")
    rep = json.loads(out.read_text()) if out.exists() else {}
    rep["pairs"] = {k: {"c": v[0], "b2": v[1]} for k, v in PAIRS.items()}
    for name, (c, b2) in PAIRS.items():
        for s in sigmas:
            key = f"{name}:{s}:{M}"
            if key in rep:
                continue
            r = exits(M, s, c, b2)
            r["cell2_own_rate"] = flip_rate(1600, A, T, DT, s, H)
            r["dip_rate"] = flip_rate(1600, A, T, DT, s, 0.0)
            rep[key] = r
            print(key, json.dumps(r), flush=True)
            out.write_text(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main()
