"""§134.4: the KernRift lattice solver agrees with SuperLU, across thread counts, in and out of core."""

import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp
import scipy.sparse.linalg as spla

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "experiments" / "krsolve"

pytestmark = pytest.mark.skipif(shutil.which("krc") is None, reason="krc (KernRift) not installed")


@pytest.fixture(scope="module")
def binary(tmp_path_factory):
    out = tmp_path_factory.mktemp("kr") / "solver"
    subprocess.run(["krc", "--arch=x86_64", "solver.kr", "-o", str(out)], cwd=SRC, check=True,
                   capture_output=True)
    return out


def lattice_system(seed=0):
    rng = np.random.default_rng(seed)
    M1, M2 = 37, 53
    I1, I2 = np.meshgrid(np.arange(M1), np.arange(M2), indexing="ij")
    I1, I2 = I1.ravel(), I2.ravel()
    keep = (I2 < 40 - (I1 // 3)) | (I1 > 20)                 # an irregular region, like the region below Sigma
    I1, I2 = I1[keep], I2[keep]
    n = len(I1)
    idx = -np.ones((M1, M2), int)
    idx[I1, I2] = np.arange(n)
    rows, cols, vals = [], [], []
    for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        J1, J2 = I1 + di, I2 + dj
        ok = (J1 >= 0) & (J1 < M1) & (J2 >= 0) & (J2 < M2)
        t = np.where(ok, idx[np.clip(J1, 0, M1 - 1), np.clip(J2, 0, M2 - 1)], -1)
        ok &= t >= 0
        r = rng.uniform(0.5, 2.0, n)
        rows.append(np.nonzero(ok)[0]); cols.append(t[ok]); vals.append(r[ok])
    Q = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n)).tocsr()
    Q = Q - sp.diags(np.asarray(Q.sum(axis=1)).ravel() + (rng.random(n) < 0.05) * 1.0)
    At = Q.T.tocsr().tolil()
    At[7, :] = 0
    At[7, 7] = 1.0                                           # pinned row, as qprocess_J
    b = np.zeros(n)
    b[7] = 1.0
    return At.tocsr(), b, (I1, I2)


def run(binary, tmp_path, M, b, coords, threads, ooc=False):
    from experiments.krsolve.krsolve import write_system
    fin, fout = tmp_path / "s.bin", tmp_path / "x.bin"
    write_system(fin, M, b, coords)
    cmd = [str(binary), str(fin), str(fout), str(threads)] + ([str(tmp_path / "f.bin")] if ooc else [])
    subprocess.run(cmd, check=True, capture_output=True)
    return np.fromfile(fout)


def test_matches_superlu_threads_and_out_of_core(binary, tmp_path):
    M, b, coords = lattice_system()
    x0 = spla.spsolve(M.tocsc(), b)
    x1 = run(binary, tmp_path, M, b, coords, 1)
    x4 = run(binary, tmp_path, M, b, coords, 4)
    xo = run(binary, tmp_path, M, b, coords, 4, ooc=True)
    assert np.max(np.abs(x1 - x0)) / np.max(np.abs(x0)) < 1e-13
    assert np.array_equal(x1, x4) and np.array_equal(x4, xo)


def test_symbolic_grid_mode(binary):
    out = subprocess.run([str(binary), "--grid", "100", "120"], check=True, capture_output=True, text=True).stdout
    assert "factor_entries=" in out and "flops~" in out
