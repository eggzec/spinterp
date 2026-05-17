"""Tests for SPINTERP_CC and SPCMPVALS_CC (full interpolation workflow)."""

import numpy as np
import pytest
import spinterp


def _spvals(f, n, d):
    """Compute hierarchical surpluses for f:[0,1]^d -> R at level n.

    Returns (list_of_levelseqs, list_of_surplus_arrays), one per level 0..n.
    """
    all_seq = []
    all_surpluses = []

    for k in range(n + 1):
        nlevels_k = spinterp.spnlevels(k, d)
        seq_k = spinterp.spgetseq(k, d, nlevels_k)  # shape (nlevels_k, d)
        tp_k = spinterp.spdim_cc(seq_k)
        x_k = spinterp.spgrid_cc(seq_k, tp_k)  # shape (tp_k, d)

        fvals = np.array([f(*x_k[i, :]) for i in range(tp_k)])

        if k == 0:
            surpluses_k = fvals.copy()
        else:
            z_prev = np.concatenate(all_surpluses)
            # Combine all previous level sequences into one array
            seq_prev = np.vstack(all_seq)

            # Interpolant at new grid points using previous-level surpluses
            interp_at_new = spinterp.spcmpvals_cc(z_prev, x_k, seq_k, seq_prev)
            surpluses_k = fvals - interp_at_new

        all_seq.append(seq_k)
        all_surpluses.append(surpluses_k)

    return all_seq, all_surpluses


def _spinterp(all_seq, all_surpluses, y):
    """Evaluate the assembled sparse grid interpolant at query points y.

    y: ndarray of shape (npoints, d).
    """
    ninterp = y.shape[0]
    result = np.zeros(ninterp)
    for seq_k, surp_k in zip(all_seq, all_surpluses):
        result += spinterp.spinterp_cc(surp_k, y, seq_k)
    return result


# ------------------------------------------------------------------
# Test: linear function  f(x) = x  in 1-D
# A level-1 CC grid reproduces linear functions exactly.
# ------------------------------------------------------------------
def test_linear_1d_exact():
    f = lambda x: x
    d, n = 1, 1
    all_seq, all_surp = _spvals(f, n, d)
    y = np.array([[0.0], [0.25], [0.5], [0.75], [1.0]])
    ip = _spinterp(all_seq, all_surp, y)
    assert ip == pytest.approx(y[:, 0], abs=1e-12)


# ------------------------------------------------------------------
# Test: constant function  f(x) = 3.14
# ------------------------------------------------------------------
def test_constant_1d():
    f = lambda x: 3.14
    d, n = 1, 2
    all_seq, all_surp = _spvals(f, n, d)
    y = np.linspace(0, 1, 20).reshape(-1, 1)
    ip = _spinterp(all_seq, all_surp, y)
    assert ip == pytest.approx(np.full(20, 3.14), abs=1e-12)


# ------------------------------------------------------------------
# Test: quadratic f(x) = x^2 error decreases with level in 1-D
# ------------------------------------------------------------------
def test_quadratic_1d_convergence():
    f = lambda x: x * x
    y = np.linspace(0.05, 0.95, 30).reshape(-1, 1)
    errors = []
    for n in range(1, 5):
        all_seq, all_surp = _spvals(f, n, 1)
        ip = _spinterp(all_seq, all_surp, y)
        errors.append(np.max(np.abs(ip - y[:, 0] ** 2)))
    for i in range(len(errors) - 1):
        assert errors[i + 1] <= errors[i] + 1e-14


# ------------------------------------------------------------------
# Test: 2-D linear  f(x,y) = x + y  at level 2
# ------------------------------------------------------------------
def test_bilinear_2d():
    f = lambda x, y: x + y
    d, n = 2, 2
    all_seq, all_surp = _spvals(f, n, d)
    rng = np.random.default_rng(42)
    pts = rng.random((20, 2))
    ip = _spinterp(all_seq, all_surp, pts)
    assert ip == pytest.approx(pts[:, 0] + pts[:, 1], abs=1e-10)


# ------------------------------------------------------------------
# Test: interpolant matches function exactly at the grid nodes
# ------------------------------------------------------------------
def test_interpolant_at_nodes():
    f = lambda x, y: x**2 + y**2 - 2 * x * y
    d, n = 2, 3
    all_seq, all_surp = _spvals(f, n, d)

    pts = np.vstack([
        spinterp.spgrid_cc(seq_k, spinterp.spdim_cc(seq_k)) for seq_k in all_seq
    ])

    ip = _spinterp(all_seq, all_surp, pts)
    exact = np.array([f(pts[i, 0], pts[i, 1]) for i in range(len(pts))])
    assert ip == pytest.approx(exact, abs=1e-10)
