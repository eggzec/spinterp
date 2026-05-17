"""Tests for SPGET_NPOINTS_*, SPSEQ2FULL, REORDER_VALS."""

import math

import numpy as np
import spinterp


# ---------------------------------------------------------------------------
# SPGET_NPOINTS_CC
# ---------------------------------------------------------------------------


def test_spget_npoints_cc_matches_spdim():
    """SPGET_NPOINTS_CC total should equal SPDIM_CC."""
    for d in range(1, 4):
        for n in range(5):
            nlevels = spinterp.spnlevels(n, d)
            seq = spinterp.spgetseq(n, d, nlevels)
            tp_from_spdim = spinterp.spdim_cc(seq)
            tp, npoints = spinterp.spget_npoints_cc(seq)
            assert tp == tp_from_spdim
            assert npoints.sum() == tp
            assert len(npoints) == nlevels


def test_spget_npoints_cc_per_subgrid():
    # D=1, N=3: seq=[[3]], lev=3>=3 → 2^(3-1)=4 points
    seq = spinterp.spgetseq(3, 1, 1)
    tp, npoints = spinterp.spget_npoints_cc(seq)
    assert tp == 4
    assert npoints[0] == 4


# ---------------------------------------------------------------------------
# SPGET_NPOINTS_M
# ---------------------------------------------------------------------------


def test_spget_npoints_m_matches_spdim():
    """SPGET_NPOINTS_M total should equal SPDIM_M(..., boundary=1)."""
    for d in range(1, 4):
        for n in range(4):
            nlevels = spinterp.spnlevels(n, d)
            seq = spinterp.spgetseq(n, d, nlevels)
            tp_spdim = spinterp.spdim_m(seq, boundary=1)
            tp, npoints = spinterp.spget_npoints_m(seq)
            assert tp == tp_spdim
            assert npoints.sum() == tp


def test_spget_npoints_m_lev0():
    # D=1, N=0: seq=[[0]], lev=0 → 3 points
    seq = spinterp.spgetseq(0, 1, 1)
    tp, npoints = spinterp.spget_npoints_m(seq)
    assert tp == 3
    assert npoints[0] == 3


# ---------------------------------------------------------------------------
# SPGET_NPOINTS_NB
# ---------------------------------------------------------------------------


def test_spget_npoints_nb_matches_spdim():
    """SPGET_NPOINTS_NB total should equal SPDIM_M(..., boundary=0)."""
    for d in range(1, 4):
        for n in range(4):
            nlevels = spinterp.spnlevels(n, d)
            seq = spinterp.spgetseq(n, d, nlevels)
            tp_spdim = spinterp.spdim_m(seq, boundary=0)
            tp, npoints = spinterp.spget_npoints_nb(seq)
            assert tp == tp_spdim
            assert npoints.sum() == tp


def test_spget_npoints_nb_lev0():
    # D=1, N=0: seq=[[0]], lev=0 → 2^0=1 point
    seq = spinterp.spgetseq(0, 1, 1)
    tp, npoints = spinterp.spget_npoints_nb(seq)
    assert tp == 1
    assert npoints[0] == 1


# ---------------------------------------------------------------------------
# SPSEQ2FULL (requires SPGETSEQ_SP)
# ---------------------------------------------------------------------------


def _build_sparse_idx(n, d, gridcode=0):
    """Return the sparse index arrays from SPGETSEQ_SP."""
    nsubgrids = math.comb(n + d, d)
    maxaddr = max(d * nsubgrids + d, 10)
    return spinterp.spgetseq_sp(n, d, gridcode, nsubgrids, maxaddr)


def test_spseq2full_d2_n2():
    """SPSEQ2FULL should reconstruct the same multi-index set as the dense approach."""
    n, d = 2, 2
    nsubgrids = math.comb(n + d, d)

    # Build sparse structure
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
    ) = _build_sparse_idx(n, d)

    # Convert sparse → full
    fullseq = spinterp.spseq2full(
        indicesndiims, indicesdims, indiceslevs, indicesaddr, d
    )
    assert fullseq.shape == (nsubgrids, d)
    assert np.all(fullseq >= 0)

    # Build dense multi-index set (all levels 0..n combined)
    dense_rows = []
    for k in range(n + 1):
        nl = spinterp.spnlevels(k, d)
        seq_k = spinterp.spgetseq(k, d, nl)
        for row in seq_k:
            dense_rows.append(tuple(row.tolist()))
    dense_set = set(dense_rows)

    # Sparse full set
    sparse_set = set(tuple(row.tolist()) for row in fullseq)

    assert sparse_set == dense_set, (
        f"Sparse multi-index set {sparse_set} != dense {dense_set}"
    )


def test_spseq2full_sum_constraints():
    """Each row of SPSEQ2FULL should sum to at most n."""
    for n in range(1, 4):
        for d in range(1, 4):
            (indicesndiims, indicesdims, indiceslevs, indicesaddr, *_) = (
                _build_sparse_idx(n, d)
            )
            fullseq = spinterp.spseq2full(
                indicesndiims, indicesdims, indiceslevs, indicesaddr, d
            )
            row_sums = fullseq.sum(axis=1)
            assert np.all(row_sums <= n)
            assert np.all(row_sums >= 0)


# ---------------------------------------------------------------------------
# REORDER_VALS
# ---------------------------------------------------------------------------


def test_reorder_vals_noop_for_d1():
    """For D=1 there is only one dimension — reordering should be identity."""
    n, d = 3, 1
    nlevels = spinterp.spnlevels(n, d)
    seq = spinterp.spgetseq(n, d, nlevels)
    tp = spinterp.spdim_cc(seq)
    z_orig = np.arange(1.0, tp + 1)
    z_out = z_orig.copy()
    spinterp.reorder_vals(z_out, seq)
    np.testing.assert_array_equal(z_out, z_orig)


def test_reorder_vals_preserves_content():
    """REORDER_VALS should preserve the set of values, only reorder."""
    n, d = 2, 2
    nlevels = spinterp.spnlevels(n, d)
    seq = spinterp.spgetseq(n, d, nlevels)
    tp = spinterp.spdim_cc(seq)
    rng = np.random.default_rng(0)
    z = rng.random(tp)
    z_out = z.copy()
    spinterp.reorder_vals(z_out, seq)
    # Content preserved (sorted values equal)
    np.testing.assert_allclose(np.sort(z_out), np.sort(z))
