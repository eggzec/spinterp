"""Tests for SPGETSEQ and SPNLEVELS."""

import numpy as np
import spinterp


def test_spnlevels_d1():
    # C(N+0, 0) = 1 for any N when D=1
    assert spinterp.spnlevels(0, 1) == 1
    assert spinterp.spnlevels(3, 1) == 1
    assert spinterp.spnlevels(5, 1) == 1


def test_spnlevels_d2():
    # C(N+1, 1) = N+1
    for n in range(6):
        assert spinterp.spnlevels(n, 2) == n + 1


def test_spnlevels_d3():
    # C(N+2, 2) = (N+1)(N+2)/2
    for n in range(6):
        expected = (n + 1) * (n + 2) // 2
        assert spinterp.spnlevels(n, 3) == expected


def test_spgetseq_d1():
    # Only one multi-index for D=1: [N]
    seq = spinterp.spgetseq(3, 1, 1)
    assert seq.shape == (1, 1)
    assert seq[0, 0] == 3


def test_spgetseq_d2_n2():
    # D=2, N=2: indices summing to 2
    # Expected: [2,0], [1,1], [0,2]
    nlevels = spinterp.spnlevels(2, 2)
    assert nlevels == 3
    seq = spinterp.spgetseq(2, 2, nlevels)
    assert seq.shape == (3, 2)
    row_sums = seq.sum(axis=1)
    assert np.all(row_sums == 2)
    rows = set(map(tuple, seq.tolist()))
    assert (2, 0) in rows
    assert (1, 1) in rows
    assert (0, 2) in rows


def test_spgetseq_d3_n3():
    # D=3, N=3: C(5,2)=10 multi-indices
    nlevels = spinterp.spnlevels(3, 3)
    assert nlevels == 10
    seq = spinterp.spgetseq(3, 3, nlevels)
    assert seq.shape == (10, 3)
    row_sums = seq.sum(axis=1)
    assert np.all(row_sums == 3)


def test_spgetseq_first_row():
    # First row always has full level in dimension 1
    for d in range(1, 5):
        for n in range(1, 5):
            nlevels = spinterp.spnlevels(n, d)
            seq = spinterp.spgetseq(n, d, nlevels)
            assert seq[0, 0] == n
            for l in range(1, d):
                assert seq[0, l] == 0


def test_spgetseq_all_nonneg():
    seq = spinterp.spgetseq(4, 4, spinterp.spnlevels(4, 4))
    assert np.all(seq >= 0)
