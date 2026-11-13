"""§134.3: the CSR row pin that replaced the LIL route in qprocess_J is identical to it."""

import numpy as np
import scipy.sparse as sp

import experiments.episode_committor as ec


def test_pin_row_matches_lil_route():
    M = sp.random(300, 300, density=0.03, random_state=2, format="csr") + sp.identity(300)
    for k in (0, 5, 150, 299):
        old = M.T.tocsr().tolil()
        old[k, :] = 0.0
        old[k, k] = 1.0
        old = old.tocsr()
        new = ec.pin_row(M.T.tocsr(), k)
        assert old.nnz == new.nnz
        assert abs(old - new).max() == 0


def test_pin_row_without_a_diagonal_entry():
    M = sp.csr_matrix(np.array([[0.0, 2.0], [3.0, 4.0]]))
    P = ec.pin_row(M, 0).toarray()
    assert np.array_equal(P, np.array([[1.0, 0.0], [3.0, 4.0]]))
