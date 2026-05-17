"""Tests for SPDIM_CC and SPGRID_CC."""

import numpy as np
import pytest
import spinterp


def _get_grid(n, d):
    """Helper: build CC sparse grid at level n, dimension d."""
    nlevels = spinterp.spnlevels(n, d)
    seq = spinterp.spgetseq(n, d, nlevels)
    totalpoints = spinterp.spdim_cc(seq)
    x = spinterp.spgrid_cc(seq, totalpoints)
    return x, seq, nlevels, totalpoints


def test_spdim_cc_d1_n0():
    # D=1, N=0: one subgrid [0], one midpoint 0.5
    nlevels = spinterp.spnlevels(0, 1)
    seq = spinterp.spgetseq(0, 1, nlevels)
    tp = spinterp.spdim_cc(seq)
    assert tp == 1


def test_spdim_cc_d1_n1():
    # D=1, N=1: subgrid [1] -> 2 boundary points
    nlevels = spinterp.spnlevels(1, 1)
    seq = spinterp.spgetseq(1, 1, nlevels)
    tp = spinterp.spdim_cc(seq)
    assert tp == 2


def test_spdim_cc_d1_n2():
    # D=1, N=2: subgrid [2] -> 2 interior points
    nlevels = spinterp.spnlevels(2, 1)
    seq = spinterp.spgetseq(2, 1, nlevels)
    tp = spinterp.spdim_cc(seq)
    assert tp == 2


def test_grid_in_unit_cube():
    for d in range(1, 4):
        for n in range(4):
            x, seq, nlevels, tp = _get_grid(n, d)
            assert x.shape == (tp, d)
            assert np.all(x >= 0.0)
            assert np.all(x <= 1.0)


def test_grid_d1_n0_coords():
    # Level 0, dim 1: single midpoint at 0.5
    x, *_ = _get_grid(0, 1)
    assert x.shape == (1, 1)
    assert x[0, 0] == pytest.approx(0.5)


def test_grid_d1_n1_coords():
    # Level 1, dim 1: boundary nodes at 0 and 1
    x, *_ = _get_grid(1, 1)
    coords = sorted(x[:, 0].tolist())
    assert coords == pytest.approx([0.0, 1.0])


def test_grid_d1_n2_coords():
    # Level 2, dim 1: interior nodes at 0.25 and 0.75
    x, *_ = _get_grid(2, 1)
    coords = sorted(x[:, 0].tolist())
    assert coords == pytest.approx([0.25, 0.75])


def test_grid_no_duplicates_d2_n2():
    x, *_ = _get_grid(2, 2)
    rows = [tuple(r) for r in x.tolist()]
    assert len(rows) == len(set(rows)), "duplicate grid points found"
