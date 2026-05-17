# Theory

## What is the Sparse Grid Interpolation Toolbox?

### Introduction

The interpolation problem addressed by sparse grid methods is an *optimal recovery* problem:
selecting a set of evaluation points such that a smooth multivariate function can be
approximated with a suitable interpolation formula.

Depending on the characteristics of the function (degree of smoothness, periodicity),
various interpolation techniques based on sparse grids exist.  All of them employ
**Smolyak's construction**, which forms the basis of all sparse grid methods.

With Smolyak's method, well-known univariate interpolation formulas are extended to the
multivariate case by using tensor products in a special, selective way.  The resulting
method requires significantly fewer support nodes than conventional interpolation on a
full grid — the difference can be several orders of magnitude as the problem dimension $d$
increases.

More formally, let $U_l^{(i)}$ be a univariate interpolation operator at level $l$ in
dimension $i$.  The Smolyak formula at level $n$ in dimension $d$ is

\[
A_{n,d}(f) = \sum_{\substack{|\mathbf{l}|_1 \leq n+d \\ \mathbf{l} \geq \mathbf{1}}}
(-1)^{n+d-|\mathbf{l}|_1} \binom{d-1}{n+d-|\mathbf{l}|_1}
\left(U_{l_1}^{(1)} \otimes \cdots \otimes U_{l_d}^{(d)}\right)(f)
\]

where $|\mathbf{l}|_1 = l_1 + \cdots + l_d$.

The most important property of the method is that the asymptotic error decay of full-grid
interpolation is preserved **up to a logarithmic factor**.  An additional benefit is its
**hierarchical structure**, which provides an estimate of the current approximation error
and enables automatic termination when a desired accuracy is reached.

### Major Features

`spinterp` includes hierarchical sparse grid interpolation algorithms based on both
piecewise multilinear and polynomial basis functions.  Special emphasis is placed on
efficient implementation that performs well even for very large dimensions $d > 10$.

Features include:

- **Five grid types** — choose the basis (piecewise linear or polynomial) and node
  distribution (Clenshaw-Curtis, Chebyshev, Gauss-Patterson, Maximum, NoBoundary)
  best suited to your function.
- **Gradients on demand** — exact derivatives of the interpolant are computed on-the-fly
  with no additional memory overhead.
- **Numerical integration** — integrate the sparse grid interpolant over its domain.
- **Dimension-adaptive grids** — automatically detect which dimensions carry more
  information and allocate points accordingly; essential for $d > 30$.
- **Sparse index sets** — optimised data structures reduce memory and runtime for
  high-dimensional problems.

---

## Pages in this section

| Page | Contents |
|---|---|
| [Linear Basis Functions](linear-basis.md) | Piecewise multilinear grids, error bounds, grid comparison |
| [Polynomial Basis Functions](polynomial-basis.md) | Chebyshev and Gauss-Patterson grids, polynomial error bounds |
| [Dimensional Adaptivity](dimension-adaptive.md) | Dimension-adaptive algorithm, Smolyak index sets |
| [Bibliography](bibliography.md) | Selected references |
