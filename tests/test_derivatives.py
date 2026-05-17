"""Tests for derivative routines, DCT surpluses, heap, and pp_deriv."""

import numpy as np
import pytest
import spinterp

from tests.test_spinterpcc import _spvals


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _build_z_seq(f, n, d):
    """Return (z_flat, seq_flat) for all levels 0..n."""
    all_seq, all_surp = _spvals(f, n, d)
    return np.concatenate(all_surp), np.vstack(all_seq)


# ---------------------------------------------------------------------------
# SPDERIV_CC (dense, full levelseq) — already partially tested via grids tests
# ---------------------------------------------------------------------------


def test_spderiv_cc_constant_zero_gradient():
    f = lambda x, y: 5.0
    z, seq = _build_z_seq(f, 2, 2)
    pts = np.array([[0.3, 0.4], [0.7, 0.2]])
    ip, grad = spinterp.spderiv_cc(z, pts, seq)
    assert ip == pytest.approx([5.0, 5.0], abs=1e-10)
    assert np.allclose(grad, 0.0, atol=1e-10)


def test_spderiv_cc_linear_exact():
    f = lambda x, y: 3 * x - 2 * y
    z, seq = _build_z_seq(f, 2, 2)
    pts = np.array([[0.25, 0.75], [0.5, 0.5], [0.1, 0.9]])
    ip, grad = spinterp.spderiv_cc(z, pts, seq)
    exact_ip = 3 * pts[:, 0] - 2 * pts[:, 1]
    assert ip == pytest.approx(exact_ip, abs=1e-9)
    assert grad[:, 0] == pytest.approx(np.full(3, 3.0), abs=1e-9)
    assert grad[:, 1] == pytest.approx(np.full(3, -2.0), abs=1e-9)


def test_spderiv_cc_1d_quadratic_convergence():
    f = lambda x: x**2
    errors = []
    pts = np.linspace(0.05, 0.95, 20).reshape(-1, 1)
    exact_grad = 2 * pts[:, 0]
    for n in range(1, 5):
        z, seq = _build_z_seq(f, n, 1)
        _, grad = spinterp.spderiv_cc(z, pts, seq)
        errors.append(np.max(np.abs(grad[:, 0] - exact_grad)))
    for i in range(len(errors) - 1):
        assert errors[i + 1] <= errors[i] + 1e-13


# ---------------------------------------------------------------------------
# SPCONT_DERIV_CC + PP_DERIV
# ---------------------------------------------------------------------------


def test_spcont_deriv_cc_constant():
    f = lambda x, y: 4.2
    z, seq = _build_z_seq(f, 3, 2)
    pts = np.array([[0.3, 0.4], [0.6, 0.7]])
    maxlev = 3
    ip, ipder, ipder2 = spinterp.spcont_deriv_cc(z, pts, seq, maxlev)
    assert ip == pytest.approx([4.2, 4.2], abs=1e-10)
    assert np.allclose(ipder, 0.0, atol=1e-10)
    assert np.allclose(ipder2, 0.0, atol=1e-10)


def test_pp_deriv_no_crash():
    """PP_DERIV should run without error and return modified ipder."""
    pts = np.linspace(0.1, 0.9, 5).reshape(-1, 1)
    y = np.asfortranarray(np.hstack([pts, pts]))
    ipder = np.asfortranarray(np.ones((5, 2)))
    ipder2 = np.asfortranarray(np.ones((5, 2)) * 1.1)
    maxlevvec = np.array([3, 3], dtype=np.int32)
    spinterp.pp_deriv(ipder, ipder2, maxlevvec, y)
    assert np.all(np.isfinite(ipder))


def test_pp_deriv_same_sign_linear_interp():
    """When ipd1 and ipd2 have same sign, pp_deriv linearly interpolates."""
    # Single point, single dimension
    y = np.array([[0.3]])
    ipder = np.array([[2.0]])
    ipder2 = np.array([[4.0]])
    maxlevvec = np.array([1], dtype=np.int32)  # maxlev=1: ytd1=0.25, step=0.5
    spinterp.pp_deriv(ipder, ipder2, maxlevvec, y)
    # Both positive → linear interp: 2 + (4-2)/0.5 * (0.3-0.25) = 2 + 4*0.05 = 2.2
    assert ipder[0, 0] == pytest.approx(2.2, abs=1e-10)


# ---------------------------------------------------------------------------
# SPDERIV_CB (Chebyshev derivative)
# ---------------------------------------------------------------------------


def test_spderiv_cb_constant_zero_gradient():
    f = lambda x: 7.0
    z, seq = _build_z_seq(f, 2, 1)
    pts = np.linspace(0.1, 0.9, 8).reshape(-1, 1)
    ip, grad = spinterp.spderiv_cb(z, pts, seq)
    assert ip == pytest.approx(np.full(8, 7.0), abs=1e-10)
    assert np.allclose(grad, 0.0, atol=1e-10)


def test_spderiv_cb_linear_1d():
    f = lambda x: 2 * x - 1
    z, seq = _build_z_seq_cb(f, 2, 1)
    pts = np.linspace(0.1, 0.9, 10).reshape(-1, 1)
    ip, grad = spinterp.spderiv_cb(z, pts, seq)
    exact_ip = 2 * pts[:, 0] - 1
    # Values should be accurate
    assert ip == pytest.approx(exact_ip, abs=1e-9)
    # Gradient should equal 2 (derivative of 2x-1)
    assert grad[:, 0] == pytest.approx(np.full(10, 2.0), abs=1e-7)


def _build_z_seq_cb(f, n, d):
    """Build CB surpluses (CB grid, CB interpolation)."""
    all_seq, all_surp = [], []
    for k in range(n + 1):
        nl = spinterp.spnlevels(k, d)
        seq_k = spinterp.spgetseq(k, d, nl)
        tp_k = spinterp.spdim_cc(seq_k)  # same point count as CC
        x_k = spinterp.spgrid_cb(seq_k, tp_k)
        fvals = np.array([f(*x_k[i, :]) for i in range(tp_k)])
        if k == 0:
            surp_k = fvals.copy()
        else:
            z_prev = np.concatenate(all_surp)
            seq_prev = np.vstack(all_seq)
            interp = spinterp.spcmpvals_cb(z_prev, x_k, seq_k, seq_prev)
            surp_k = fvals - interp
        all_seq.append(seq_k)
        all_surp.append(surp_k)
    return np.concatenate(all_surp), np.vstack(all_seq)


def test_spderiv_cb_gradient_convergence():
    f = lambda x: x * x
    pts = np.linspace(0.05, 0.95, 20).reshape(-1, 1)
    errors = []
    for n in range(2, 6):
        z, seq = _build_z_seq_cb(f, n, 1)
        _, grad = spinterp.spderiv_cb(z, pts, seq)
        errors.append(np.max(np.abs(grad[:, 0] - 2 * pts[:, 0])))
    for i in range(len(errors) - 1):
        assert errors[i + 1] <= errors[i] + 1e-12


# ---------------------------------------------------------------------------
# SPCMPVALS_CB_DCT vs SPCMPVALS_CB (should give same surpluses)
# ---------------------------------------------------------------------------


def test_spcmpvals_cb_dct_matches_barycentric():
    """DCT and barycentric surplus computation should agree."""
    f = lambda x: x**2 - 0.5
    d, n = 1, 3
    # Build first level surpluses (level 0 only)
    nlevels_0 = spinterp.spnlevels(0, d)
    seq_0 = spinterp.spgetseq(0, d, nlevels_0)
    tp_0 = spinterp.spdim_cc(seq_0)
    x_0 = spinterp.spgrid_cc(seq_0, tp_0)
    surp_0 = np.array([f(*x_0[i, :]) for i in range(tp_0)])

    # Compute level 1 surpluses with both methods
    nlevels_1 = spinterp.spnlevels(1, d)
    seq_1 = spinterp.spgetseq(1, d, nlevels_1)
    tp_1 = spinterp.spdim_cc(seq_1)
    x_1 = spinterp.spgrid_cc(seq_1, tp_1)
    fvals_1 = np.array([f(*x_1[i, :]) for i in range(tp_1)])

    # Barycentric method
    interp_bary = spinterp.spcmpvals_cb(surp_0, x_1, seq_1, seq_0)
    surp_1_bary = fvals_1 - interp_bary

    # DCT method
    interp_dct = spinterp.spcmpvals_cb_dct(surp_0, x_1, seq_1, seq_0)
    surp_1_dct = fvals_1 - interp_dct

    np.testing.assert_allclose(
        surp_1_dct,
        surp_1_bary,
        atol=1e-10,
        err_msg="DCT and barycentric surpluses disagree",
    )


# ---------------------------------------------------------------------------
# SORT_HEAP and POP_HEAP
# ---------------------------------------------------------------------------


def test_sort_heap_basic():
    """After inserting element at end, sort_heap should restore max-heap."""
    g = np.array([10.0, 5.0, 8.0, 3.0, 7.0])
    # Build a valid heap manually (1-indexed): [8,5,7,3,10] → after sort
    # We simulate: existing heap [8,5,7,3], add 10 at position 5
    a = np.array([1, 2, 3, 4, 5], dtype=np.int32)  # indices into g
    # After sorting: a[0] should hold index with max g value
    spinterp.sort_heap(a, g)
    # Index stored at a[0] should have the maximum g value
    best_idx = a[0] - 1  # convert to 0-based
    assert g[best_idx] == g.max()


def test_pop_heap_returns_max():
    """POP_HEAP should return the index of the maximum element."""
    g = np.array([0.0, 4.0, 9.0, 2.0, 7.0])  # 0-based values
    # Build a max-heap: create the heap array manually
    # Max heap for g[1..5] (1-based): index 3 has g=9 → should be at top
    # Build by repeated sort_heap
    a = np.array([1, 2, 3, 4, 5], dtype=np.int32)
    # Insert elements one by one
    for i in range(1, 6):
        spinterp.sort_heap(a[:i], g)

    # Now pop: should give index 3 (0-based g[2]=9 via 1-based a index=3)
    idx = spinterp.pop_heap(a, g)
    # idx is 1-based into g
    assert g[idx - 1] == g.max()


def test_sort_then_pop_heap():
    """Sort + pop should extract maximum repeatedly."""
    n = 6
    g = np.array([0.0, 3.0, 8.0, 1.0, 5.0, 4.0, 7.0])  # 0-based, 1..6 used
    a = np.arange(1, n + 1, dtype=np.int32)

    # Build heap
    for i in range(1, n + 1):
        spinterp.sort_heap(a[:i], g)

    # Pop all in order — should give descending g values
    vals = []
    na = n
    while na > 0:
        idx = spinterp.pop_heap(a[:na], g)
        vals.append(g[idx - 1])
        na -= 1
        if na > 0:
            spinterp.sort_heap(a[:na], g)

    # Should be sorted descending
    for i in range(len(vals) - 1):
        assert vals[i] >= vals[i + 1]
