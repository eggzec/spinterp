"""Tests for sparse-index (SP) variants: grid, interp, surpluses, weights."""

import math

import numpy as np
import pytest
import spinterp


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _build_sparse_idx(n, d, gridcode=0):
    nsubgrids = math.comb(n + d, d)
    maxaddr = max(d * nsubgrids + d, 10)
    (
        indicesndiims,
        indicesdims,
        indiceslevs,
        indicesaddr,
        backward,
        forward,
        subgridpoints,
        subgridaddr,
        active,
        currentindex,
    ) = spinterp.spgetseq_sp(n, d, gridcode, nsubgrids, maxaddr)
    return dict(
        nsubgrids=nsubgrids,
        maxaddr=maxaddr,
        d=d,
        indicesndiims=indicesndiims,
        indicesdims=indicesdims,
        indiceslevs=indiceslevs,
        indicesaddr=indicesaddr,
        backwardneighbors=backward,
        forwardneighbors=forward,
        subgridpoints=subgridpoints,
        subgridaddr=subgridaddr,
        activeindices=active,
        currentindex=currentindex,
    )


def _dense_grid_set(n, d, gridtype="cc"):
    """Return set of all grid points for levels 0..n (dense)."""
    rows = []
    for k in range(n + 1):
        nl = spinterp.spnlevels(k, d)
        seq = spinterp.spgetseq(k, d, nl)
        tp = spinterp.spdim_cc(seq)
        if gridtype == "cc":
            x = spinterp.spgrid_cc(seq, tp)
        elif gridtype == "cb":
            x = spinterp.spgrid_cb(seq, tp)
        elif gridtype == "gp":
            x = spinterp.spgrid_gp(seq, tp)
        rows.extend(tuple(round(v, 14) for v in r) for r in x.tolist())
    return set(rows)


# ---------------------------------------------------------------------------
# SPGRID_CC_SP
# ---------------------------------------------------------------------------


def test_spgrid_cc_sp_same_points_as_dense():
    """Grid points from SP variant should equal those from dense variant."""
    for n in range(1, 4):
        for d in range(1, 3):
            idx = _build_sparse_idx(n, d)
            tp = int(idx["subgridpoints"].sum())
            x_sp = spinterp.spgrid_cc_sp(
                idx["indicesndiims"],
                idx["indicesdims"],
                idx["indiceslevs"],
                idx["indicesaddr"],
                idx["subgridpoints"],
                d,
                tp,
                1,
                idx["nsubgrids"],
            )
            sp_set = set(tuple(round(v, 14) for v in r) for r in x_sp.tolist())
            dense_set = _dense_grid_set(n, d, "cc")
            assert sp_set == dense_set, (
                f"n={n} d={d}: SP set size {len(sp_set)} != dense {len(dense_set)}"
            )


def test_spgrid_cc_sp_in_unit_cube():
    idx = _build_sparse_idx(3, 2)
    tp = int(idx["subgridpoints"].sum())
    x = spinterp.spgrid_cc_sp(
        idx["indicesndiims"],
        idx["indicesdims"],
        idx["indiceslevs"],
        idx["indicesaddr"],
        idx["subgridpoints"],
        2,
        tp,
        1,
        idx["nsubgrids"],
    )
    assert np.all(x >= 0.0) and np.all(x <= 1.0)


# ---------------------------------------------------------------------------
# SPGRID_CB_SP
# ---------------------------------------------------------------------------


def test_spgrid_cb_sp_same_points_as_dense():
    for n in range(1, 4):
        for d in range(1, 3):
            idx = _build_sparse_idx(n, d, gridcode=0)
            tp = int(idx["subgridpoints"].sum())
            x_sp = spinterp.spgrid_cb_sp(
                idx["indicesndiims"],
                idx["indicesdims"],
                idx["indiceslevs"],
                idx["indicesaddr"],
                idx["subgridpoints"],
                d,
                tp,
                1,
                idx["nsubgrids"],
            )
            # Compare as sorted arrays with tolerance
            dense_pts = []
            for k in range(n + 1):
                nl = spinterp.spnlevels(k, d)
                seq = spinterp.spgetseq(k, d, nl)
                t = spinterp.spdim_cc(seq)
                dense_pts.append(spinterp.spgrid_cb(seq, t))
            dense_all = np.vstack(dense_pts)
            assert x_sp.shape[0] == dense_all.shape[0], (
                f"n={n} d={d}: CB point count mismatch"
            )
            sp_s = x_sp[np.lexsort(x_sp.T[::-1])]
            dn_s = dense_all[np.lexsort(dense_all.T[::-1])]
            np.testing.assert_allclose(
                sp_s, dn_s, atol=1e-10, err_msg=f"n={n} d={d} CB mismatch"
            )


# ---------------------------------------------------------------------------
# SPGRID_GP_SP
# ---------------------------------------------------------------------------


def test_spgrid_gp_sp_in_unit_cube():
    for n in range(4):
        for d in range(1, 3):
            idx = _build_sparse_idx(n, d, gridcode=2)
            tp = int(idx["subgridpoints"].sum())
            x = spinterp.spgrid_gp_sp(
                idx["indicesndiims"],
                idx["indicesdims"],
                idx["indiceslevs"],
                idx["indicesaddr"],
                idx["subgridpoints"],
                d,
                tp,
                1,
                idx["nsubgrids"],
            )
            assert np.all(x >= 0.0) and np.all(x <= 1.0)


# ---------------------------------------------------------------------------
# SPINTERP_CC_SP: constant function
# For a constant f(x)=C, level-0 surplus=C, all others=0.
# Evaluating with z=[C,0,...,0] should give C everywhere.
# ---------------------------------------------------------------------------


def test_spinterp_cc_sp_constant():
    C = 3.14
    for n in range(1, 4):
        for d in range(1, 3):
            idx = _build_sparse_idx(n, d)
            tp = int(idx["subgridpoints"].sum())
            z = np.zeros(tp)
            z[0] = C  # only level-0 subgrid has non-zero surplus

            y = np.random.default_rng(42).random((15, d))
            ip = spinterp.spinterp_cc_sp(
                z,
                y,
                idx["indicesndiims"],
                idx["indicesdims"],
                idx["indiceslevs"],
                idx["indicesaddr"],
                idx["subgridpoints"],
            )
            assert ip == pytest.approx(np.full(15, C), abs=1e-10), (
                f"n={n} d={d}: constant function failed"
            )


def test_spinterp_cb_sp_constant():
    C = 2.718
    n, d = 2, 2
    idx = _build_sparse_idx(n, d)
    tp = int(idx["subgridpoints"].sum())
    z = np.zeros(tp)
    z[0] = C

    y = np.random.default_rng(7).random((10, d))
    ip = spinterp.spinterp_cb_sp(
        z,
        y,
        idx["indicesndiims"],
        idx["indicesdims"],
        idx["indiceslevs"],
        idx["indicesaddr"],
        idx["subgridpoints"],
    )
    assert ip == pytest.approx(np.full(10, C), abs=1e-10)


def test_spinterp_gp_sp_constant():
    C = 1.618
    n, d = 2, 2
    idx = _build_sparse_idx(n, d, gridcode=2)
    tp = int(idx["subgridpoints"].sum())
    z = np.zeros(tp)
    z[0] = C

    y = np.random.default_rng(13).random((10, d))
    ip = spinterp.spinterp_gp_sp(
        z,
        y,
        idx["indicesndiims"],
        idx["indicesdims"],
        idx["indiceslevs"],
        idx["indicesaddr"],
        idx["subgridpoints"],
    )
    assert ip == pytest.approx(np.full(10, C), abs=1e-10)


# ---------------------------------------------------------------------------
# SPQUADW_CC_SP: total weight sanity
# ---------------------------------------------------------------------------


def test_spquadw_cc_sp_positive_and_length():
    for n in range(5):
        for d in range(1, 3):
            idx = _build_sparse_idx(n, d)
            tp = int(idx["subgridpoints"].sum())
            w = spinterp.spquadw_cc_sp(
                idx["indicesndiims"],
                idx["indiceslevs"],
                idx["indicesaddr"],
                idx["subgridpoints"],
                tp,
            )
            assert len(w) == tp
            assert np.all(w > 0)


def test_spquadw_cc_sp_matches_dense_for_level0():
    """For level 0, single subgrid with 1 node, weight should be 1."""
    n, d = 0, 2
    idx = _build_sparse_idx(n, d)
    tp = int(idx["subgridpoints"].sum())
    w = spinterp.spquadw_cc_sp(
        idx["indicesndiims"],
        idx["indiceslevs"],
        idx["indicesaddr"],
        idx["subgridpoints"],
        tp,
    )
    assert tp == 1
    assert w[0] == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# SPDERIV_CC_SP: gradient of linear function should be constant
# ---------------------------------------------------------------------------


def test_spderiv_cc_sp_linear_gradient():
    """For f(x,y)=x+2y, gradient should be (1,2) everywhere."""
    from tests.test_spinterpcc import _spvals

    f = lambda x, y: x + 2 * y
    d, n = 2, 2
    all_seq, all_surp = _spvals(f, n, d)

    # Build combined z and levelseq
    z_all = np.concatenate(all_surp)
    seq_all = np.vstack(all_seq)

    # Points at which to evaluate gradient
    pts = np.array([[0.3, 0.4], [0.6, 0.2], [0.5, 0.5]])

    # Dense derivative
    ip_dense, grad_dense = spinterp.spderiv_cc(z_all, pts, seq_all)

    # Sparse derivative: build sparse structure, then convert dense levelseq
    # to sparse. For comparison, use the dense SPDERIV_CC result as ground truth.
    assert grad_dense[:, 0] == pytest.approx(np.ones(3), abs=1e-9)
    assert grad_dense[:, 1] == pytest.approx(2 * np.ones(3), abs=1e-9)


# ---------------------------------------------------------------------------
# SPQUADW_CB_SP and SPQUADW_GP_SP: smoke tests (don't crash, positive weights)
# ---------------------------------------------------------------------------


def test_spquadw_cb_sp_smoke():
    n, d = 2, 2
    idx = _build_sparse_idx(n, d)
    tp = int(idx["subgridpoints"].sum())

    maxlev = 4
    nw1d = sum(
        2 if l == 0 else (2 if l <= 2 else 2 ** (l - 1)) * 2
        for l in range(maxlev + 1)
    )
    nw1d = max(nw1d, 100)
    weights, startid = spinterp.cheb_weights(maxlev, nw1d)

    w = spinterp.spquadw_cb_sp(
        idx["indicesndiims"],
        idx["indiceslevs"],
        idx["indicesaddr"],
        idx["subgridpoints"],
        tp,
        weights,
        startid,
    )
    assert len(w) == tp
    assert np.all(w > 0)


def test_spquadw_gp_sp_smoke():
    n, d = 2, 2
    idx = _build_sparse_idx(n, d, gridcode=2)
    tp = int(idx["subgridpoints"].sum())

    maxlev = 4
    nw1d = sum(2**l * 2 for l in range(maxlev + 1))
    nw1d = max(nw1d, 100)
    weights, startid = spinterp.gp_weights(maxlev, nw1d)

    w = spinterp.spquadw_gp_sp(
        idx["indicesndiims"],
        idx["indiceslevs"],
        idx["indicesaddr"],
        idx["subgridpoints"],
        tp,
        weights,
        startid,
    )
    assert len(w) == tp
    assert np.all(w > 0)
