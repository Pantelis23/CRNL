"""Python side of krsolve (§134.4): write a lattice system, run the KernRift solver, read the solution.

kr_solve(M, rhs, coords) solves M x = rhs. coords = (x, y): each unknown's lattice coordinates; the pattern
must connect only lattice neighbours (it is symmetrised here with explicit zeros, as the solver requires).
The binary is experiments/krsolve/solver (build: krc --arch=x86_64 solver.kr -o solver), or KRSOLVE_BIN.
KRSOLVE_CMD may prefix the command (e.g. "wsl -d NVIDIA-Workbench --" with KRSOLVE_BIN a Linux path), and
KRSOLVE_DIR sets where the exchange files go (default: a temp directory); KRSOLVE_FACTOR, if set, is the file the
factor is streamed to (out-of-core, in the solver's own path space).
"""

from __future__ import annotations

import os
import shlex
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import scipy.sparse as sp

MAGIC = 0x564C4F53524B
HERE = Path(__file__).resolve().parent
LAST = {}


def _symmetrised(M):
    C = sp.coo_matrix(M)
    rows = np.concatenate([C.row, C.col])
    cols = np.concatenate([C.col, C.row])
    vals = np.concatenate([C.data, np.zeros(len(C.data))])
    S = sp.csr_matrix((vals, (rows, cols)), shape=M.shape)    # duplicates summed; T-only slots stay 0
    S.sort_indices()
    return S


def write_system(path, M, rhs, coords):
    S = _symmetrised(M)
    n = S.shape[0]
    with open(path, "wb") as fh:
        np.array([MAGIC, n, S.nnz], dtype="<u8").tofile(fh)
        np.asarray(coords[0], dtype="<u4").tofile(fh)
        np.asarray(coords[1], dtype="<u4").tofile(fh)
        S.indptr.astype("<u8").tofile(fh)
        S.indices.astype("<u4").tofile(fh)
        S.data.astype("<f8").tofile(fh)
        np.asarray(rhs, dtype="<f8").tofile(fh)
    return n, S.nnz


def kr_solve(M, rhs, coords, threads=None):
    rhs = np.asarray(rhs, dtype=float)
    if rhs.ndim != 1:
        return np.column_stack([kr_solve(M, rhs[:, j], coords, threads) for j in range(rhs.shape[1])])
    threads = threads or int(os.environ.get("KRSOLVE_THREADS", str(os.cpu_count() or 1)))
    d = os.environ.get("KRSOLVE_DIR") or tempfile.mkdtemp(prefix="krsolve_")
    Path(d).mkdir(parents=True, exist_ok=True)
    fin, fout = str(Path(d) / "system.bin"), str(Path(d) / "x.bin")
    n, nnz = write_system(fin, M, rhs, coords)
    del M
    binary = os.environ.get("KRSOLVE_BIN", str(HERE / "solver"))
    pre = shlex.split(os.environ.get("KRSOLVE_CMD", ""))
    to_bin = os.environ.get("KRSOLVE_PATHMAP")          # e.g. "C:/=/mnt/c/" for WSL
    a, b = (to_bin.split("=", 1) if to_bin else ("", ""))
    mp = (lambda p: p.replace("\\", "/").replace(a, b, 1)) if to_bin else (lambda p: p)
    cmd = pre + [binary, mp(fin), mp(fout), str(threads)]
    if os.environ.get("KRSOLVE_FACTOR"):                 # out-of-core: factor blocks streamed to this file
        cmd.append(os.environ["KRSOLVE_FACTOR"])
    r = subprocess.run(cmd, capture_output=True, text=True)
    LAST.clear()
    LAST.update({"stdout": r.stdout, "returncode": r.returncode, "n": n, "nnz": nnz})
    if r.returncode != 0:
        raise RuntimeError(f"krsolve failed ({r.returncode}): {r.stdout}{r.stderr}")
    x = np.fromfile(fout, dtype="<f8")
    if len(x) != n:
        raise RuntimeError("krsolve: wrong output length")
    for line in r.stdout.splitlines():
        if "factor_entries=" in line:
            LAST["factor_entries"] = int(line.split("factor_entries=")[1].split()[0])
        if "residual" in line:
            LAST["residual"] = float(line.split("=")[-1])
    if not os.environ.get("KRSOLVE_DIR"):
        for f in (fin, fout):
            os.remove(f)
        os.rmdir(d)
    return x
