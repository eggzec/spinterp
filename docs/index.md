# spinterp

**Sparse Grid Interpolation Toolbox**

`spinterp` is a high-performance toolkit for **multivariate interpolation, quadrature, and
gradient computation on sparse grids**.  It supports five grid types — from piecewise-linear
Clenshaw-Curtis to polynomial Chebyshev and Gauss-Patterson grids — and scales to hundreds
of dimensions through hierarchical surpluses and dimension-adaptive refinement.

---

## Overview

Sparse grid methods overcome the *curse of dimensionality*: the exponential growth of grid
points in a full tensor-product discretisation.  Using Smolyak's construction, one selects
only those grid points whose hierarchical surplus exceeds a tolerance, achieving accuracy
comparable to a full grid at a fraction of the cost — making interpolation, integration, and
optimisation in $d \gg 10$ dimensions practical.

The error of the sparse grid interpolant $A_{q,d}(f)$ satisfies

\[
\|f - A_{q,d}(f)\|_\infty = O\!\left(N^{-r}\,(\log N)^{(d-1)(r+1)}\right)
\]

where $N$ is the number of support nodes and $r$ depends on the smoothness of $f$ and the
chosen basis.

---

## Grid types

| Grid type | Basis | Max depth |
|---|---|---|
| **Clenshaw-Curtis** | Piecewise linear | 8 |
| **Chebyshev** | Polynomial (DCT) | 10 |
| **Gauss-Patterson** | Nested Gaussian | 6 |
| **Maximum** | Piecewise linear | 8 |
| **NoBoundary** | Piecewise linear | 8 |

---

## Core routines

| Function | Description |
|---|---|
| `spgetseq` | Generate multi-index level sequences |
| `spgrid_cc` / `spgrid_cb` / … | Sparse grid node coordinates per grid type |
| `spinterp_cc` / `spinterp_cb` / … | Evaluate the interpolant at query points |
| `spcmpvals_cc` / `spcmpvals_cb` / … | Compute hierarchical surpluses |
| `spderiv_cc` / `spderiv_cb` | Interpolant values and exact gradient vectors |
| `spquadw_cc` / `spquadw_cb` / … | Quadrature weight vectors for integration |

---

## Attribution

Original MATLAB toolbox by **W. Andreas Klimke**, Universität Stuttgart.

> Klimke, A., Wohlmuth, B. (2005).
> *Algorithm 847: spinterp — Piecewise Multilinear Hierarchical Sparse Grid Interpolation in MATLAB.*
> ACM Transactions on Mathematical Software, **31**(4), 561–579.
