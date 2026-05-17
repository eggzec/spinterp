"""Tests for NoBoundary, Maximum, Chebyshev, and GP grid types."""

import numpy as np
import pytest
import spinterp


# ---------------------------------------------------------------------------
# Generic helper: build surpluses for a given grid type
# ---------------------------------------------------------------------------


def _spvals(f, n, d, gridtype="cc"):
    """Build hierarchical surpluses for f:[0,1]^d -> R at level n.

    gridtype: 'cc', 'nb', 'm', 'cb', 'gp'
    """
    dim_fns = {
        "cc": (
            spinterp.spdim_cc,
            spinterp.spgrid_cc,
            spinterp.spcmpvals_cc,
            spinterp.spinterp_cc,
        ),
        "nb": (
            lambda s: spinterp.spdim_m(s, boundary=0),
            spinterp.spgrid_nb,
            spinterp.spcmpvals_nb,
            spinterp.spinterp_nb,
        ),
        "m": (
            lambda s: spinterp.spdim_m(s, boundary=1),
            spinterp.spgrid_m,
            spinterp.spcmpvals_m,
            spinterp.spinterp_m,
        ),
        "cb": (
            spinterp.spdim_cc,
            spinterp.spgrid_cb,
            spinterp.spcmpvals_cb,
            spinterp.spinterp_cb,
        ),
        "gp": (
            lambda s: spinterp.spdim_m(s, boundary=0),
            spinterp.spgrid_gp,
            spinterp.spcmpvals_gp,
            spinterp.spinterp_gp,
        ),
    }
    spdim, spgrid, spcmpvals, _ = dim_fns[gridtype]

    all_seq, all_surp = [], []

    for k in range(n + 1):
        nlevels_k = spinterp.spnlevels(k, d)
        seq_k = spinterp.spgetseq(k, d, nlevels_k)
        tp_k = spdim(seq_k)
        x_k = spgrid(seq_k, tp_k)

        fvals = np.array([f(*x_k[i, :]) for i in range(tp_k)])

        if k == 0:
            surp_k = fvals.copy()
        else:
            z_prev = np.concatenate(all_surp)
            seq_prev = np.vstack(all_seq)
            interp = spcmpvals(z_prev, x_k, seq_k, seq_prev)
            surp_k = fvals - interp

        all_seq.append(seq_k)
        all_surp.append(surp_k)

    return all_seq, all_surp


def _spinterp(all_seq, all_surp, y, gridtype="cc"):
    interp_fns = {
        "cc": spinterp.spinterp_cc,
        "nb": spinterp.spinterp_nb,
        "m": spinterp.spinterp_m,
        "cb": spinterp.spinterp_cb,
        "gp": spinterp.spinterp_gp,
    }
    fn = interp_fns[gridtype]
    result = np.zeros(y.shape[0])
    for seq_k, surp_k in zip(all_seq, all_surp):
        result += fn(surp_k, y, seq_k)
    return result


# ---------------------------------------------------------------------------
# NoBoundary grid tests
# ---------------------------------------------------------------------------


def test_nb_grid_in_unit_cube():
    for d in range(1, 4):
        for n in range(4):
            nlevels = spinterp.spnlevels(n, d)
            seq = spinterp.spgetseq(n, d, nlevels)
            tp = spinterp.spdim_m(seq, boundary=0)
            x = spinterp.spgrid_nb(seq, tp)
            assert np.all(x > 0.0) and np.all(x < 1.0), (
                f"NB grid d={d} n={n}: nodes outside (0,1)"
            )


def test_nb_grid_d1_n0():
    nlevels = spinterp.spnlevels(0, 1)
    seq = spinterp.spgetseq(0, 1, nlevels)
    tp = spinterp.spdim_m(seq, boundary=0)
    x = spinterp.spgrid_nb(seq, tp)
    assert tp == 1
    assert x[0, 0] == pytest.approx(0.5)


def test_nb_constant_interpolation():
    f = lambda x: 2.71828
    all_seq, all_surp = _spvals(f, 3, 1, "nb")
    y = np.linspace(0.05, 0.95, 20).reshape(-1, 1)
    ip = _spinterp(all_seq, all_surp, y, "nb")
    assert ip == pytest.approx(np.full(20, 2.71828), abs=1e-10)


def test_nb_linear_convergence():
    f = lambda x: x
    errors = []
    y = np.linspace(0.05, 0.95, 30).reshape(-1, 1)
    for n in range(1, 5):
        all_seq, all_surp = _spvals(f, n, 1, "nb")
        ip = _spinterp(all_seq, all_surp, y, "nb")
        errors.append(np.max(np.abs(ip - y[:, 0])))
    for i in range(len(errors) - 1):
        assert errors[i + 1] <= errors[i] + 1e-13


# ---------------------------------------------------------------------------
# Maximum grid tests
# ---------------------------------------------------------------------------


def test_m_grid_d1_n0():
    nlevels = spinterp.spnlevels(0, 1)
    seq = spinterp.spgetseq(0, 1, nlevels)
    tp = spinterp.spdim_m(seq, boundary=1)
    x = spinterp.spgrid_m(seq, tp)
    assert tp == 3
    coords = sorted(x[:, 0].tolist())
    assert coords == pytest.approx([0.0, 0.5, 1.0])


def test_m_grid_d1_n1():
    nlevels = spinterp.spnlevels(1, 1)
    seq = spinterp.spgetseq(1, 1, nlevels)
    tp = spinterp.spdim_m(seq, boundary=1)
    x = spinterp.spgrid_m(seq, tp)
    assert tp == 2
    coords = sorted(x[:, 0].tolist())
    assert coords == pytest.approx([0.25, 0.75])


def test_m_constant_interpolation():
    f = lambda x: 1.618
    all_seq, all_surp = _spvals(f, 3, 1, "m")
    y = np.linspace(0.0, 1.0, 20).reshape(-1, 1)
    ip = _spinterp(all_seq, all_surp, y, "m")
    assert ip == pytest.approx(np.full(20, 1.618), abs=1e-10)


def test_m_linear_convergence():
    f = lambda x: x
    errors = []
    y = np.linspace(0.0, 1.0, 30).reshape(-1, 1)
    for n in range(1, 5):
        all_seq, all_surp = _spvals(f, n, 1, "m")
        ip = _spinterp(all_seq, all_surp, y, "m")
        errors.append(np.max(np.abs(ip - y[:, 0])))
    for i in range(len(errors) - 1):
        assert errors[i + 1] <= errors[i] + 1e-13


# ---------------------------------------------------------------------------
# Chebyshev grid tests
# ---------------------------------------------------------------------------


def test_cb_grid_same_size_as_cc():
    for d in range(1, 4):
        for n in range(4):
            nlevels = spinterp.spnlevels(n, d)
            seq = spinterp.spgetseq(n, d, nlevels)
            tp_cc = spinterp.spdim_cc(seq)
            tp_cb = spinterp.spdim_cc(seq)
            assert tp_cc == tp_cb


def test_cb_grid_in_unit_cube():
    for d in range(1, 3):
        for n in range(4):
            nlevels = spinterp.spnlevels(n, d)
            seq = spinterp.spgetseq(n, d, nlevels)
            tp = spinterp.spdim_cc(seq)
            x = spinterp.spgrid_cb(seq, tp)
            assert np.all(x >= 0.0) and np.all(x <= 1.0)


def test_cb_constant_interpolation():
    f = lambda x: 3.14
    all_seq, all_surp = _spvals(f, 2, 1, "cb")
    y = np.linspace(0, 1, 20).reshape(-1, 1)
    ip = _spinterp(all_seq, all_surp, y, "cb")
    assert ip == pytest.approx(np.full(20, 3.14), abs=1e-10)


def test_cb_linear_convergence():
    f = lambda x: x
    errors = []
    y = np.linspace(0.05, 0.95, 30).reshape(-1, 1)
    for n in range(1, 5):
        all_seq, all_surp = _spvals(f, n, 1, "cb")
        ip = _spinterp(all_seq, all_surp, y, "cb")
        errors.append(np.max(np.abs(ip - y[:, 0])))
    for i in range(len(errors) - 1):
        assert errors[i + 1] <= errors[i] + 1e-12


# ---------------------------------------------------------------------------
# GP abscissae and weights
# ---------------------------------------------------------------------------


def test_gp_absc_level0():
    x = spinterp.gp_absc(0, 1)
    assert x[0] == pytest.approx(0.5)


def test_gp_absc_level1():
    x = spinterp.gp_absc(1, 3)
    assert sorted(x.tolist()) == pytest.approx(
        [
            0.5 - 0.77459666924148337704 / 2,
            0.5,
            0.5 + 0.77459666924148337704 / 2,
        ],
        abs=1e-12,
    )


def test_gp_absc_symmetric():
    for level in range(1, 7):
        nx = 2 ** (level + 1) - 1
        x = spinterp.gp_absc(level, nx)
        mid = nx // 2
        assert x[mid] == pytest.approx(0.5, abs=1e-12)
        for i in range(mid):
            assert x[i] + x[nx - 1 - i] == pytest.approx(1.0, abs=1e-12)


def test_gp_grid_in_unit_cube():
    for d in range(1, 3):
        for n in range(4):
            nlevels = spinterp.spnlevels(n, d)
            seq = spinterp.spgetseq(n, d, nlevels)
            tp = spinterp.spdim_m(seq, boundary=0)
            x = spinterp.spgrid_gp(seq, tp)
            assert np.all(x >= 0.0) and np.all(x <= 1.0)


# ---------------------------------------------------------------------------
# Quadrature weights
# ---------------------------------------------------------------------------


def test_quadw_cc_all_positive():
    for d in range(1, 4):
        for n in range(5):
            nlevels = spinterp.spnlevels(n, d)
            seq = spinterp.spgetseq(n, d, nlevels)
            tp = spinterp.spdim_cc(seq)
            w = spinterp.spquadw_cc(seq, tp)
            assert len(w) == tp
            assert np.all(w > 0)


def test_quadw_nb_all_positive():
    for d in range(1, 3):
        for n in range(4):
            nlevels = spinterp.spnlevels(n, d)
            seq = spinterp.spgetseq(n, d, nlevels)
            tp = spinterp.spdim_m(seq, boundary=0)
            w = spinterp.spquadw_nb(seq, tp)
            assert np.all(w > 0)


def test_quadw_m_all_positive():
    for d in range(1, 3):
        for n in range(4):
            nlevels = spinterp.spnlevels(n, d)
            seq = spinterp.spgetseq(n, d, nlevels)
            tp = spinterp.spdim_m(seq, boundary=1)
            w = spinterp.spquadw_m(seq, tp)
            assert np.all(w > 0)


# ---------------------------------------------------------------------------
# Derivative
# ---------------------------------------------------------------------------


def test_deriv_cc_constant():
    f = lambda x, y: 5.0
    d, n = 2, 2
    all_seq, all_surp = _spvals(f, n, d)
    z = np.concatenate(all_surp)
    seq = np.vstack(all_seq)
    y = np.array([[0.3, 0.7], [0.5, 0.5]])
    ip, ipder = spinterp.spderiv_cc(z, y, seq)
    assert ip == pytest.approx([5.0, 5.0], abs=1e-10)
    assert np.allclose(ipder, 0.0, atol=1e-10)


def test_deriv_cc_linear():
    f = lambda x, y: x + 2 * y
    d, n = 2, 2
    all_seq, all_surp = _spvals(f, n, d)
    z = np.concatenate(all_surp)
    seq = np.vstack(all_seq)
    y = np.array([[0.3, 0.4], [0.6, 0.2]])
    ip, ipder = spinterp.spderiv_cc(z, y, seq)
    exact = y[:, 0] + 2 * y[:, 1]
    assert ip == pytest.approx(exact, abs=1e-10)
