# Dimensional Adaptivity

## Motivation

With the standard (isotropic) sparse grid approach all dimensions are treated equally: the
number of grid points in each coordinate direction is the same.  In practice, however,
many high-dimensional functions are *anisotropic* — some input variables carry far more
information than others.  Unfortunately, which dimensions are important is generally not
known in advance.

**Dimension-adaptive sparse grids** address this by automatically detecting which dimensions
matter and allocating more nodes there, without wasting function evaluations.

---

## The Gerstner-Griebel algorithm

The dimension-adaptive algorithm implemented in `spinterp` is based on Gerstner and
Griebel [[8]](bibliography.md), enhanced with the performance improvements from Klimke
[[3, ch. 3]](bibliography.md).

The key idea is to maintain a **profit indicator** $g(\mathbf{l})$ for each admissible
multi-index $\mathbf{l}$ that measures the expected gain from adding the corresponding
subgrid to the interpolant.  In the greedy version,

\[
g(\mathbf{l}) = \frac{|\alpha_{\mathbf{l}}|}{n(\mathbf{l})}
\]

where $|\alpha_{\mathbf{l}}|$ is the $\ell^\infty$ norm of the hierarchical surpluses on
subgrid $\mathbf{l}$ and $n(\mathbf{l})$ is the number of new points this subgrid adds.

A multi-index $\mathbf{l}$ is **admissible** (i.e. may be added to the active set) if all
its backward neighbours $\mathbf{l} - \mathbf{e}_k$ (for every $k$ with $l_k > 0$) are
already present in the interpolant.

The degree of dimensional adaptivity is controlled by the scalar parameter
$\delta \in [0, 1]$:

- $\delta = 1$ — fully adaptive (greedy): always refine the subgrid with the highest
  profit.
- $\delta = 0$ — no adaptivity: reduces to the standard isotropic sparse grid.

In practice, a purely greedy strategy can occasionally underestimate the true error in some
dimensions (because indicators are based on a finite set of surpluses), which may stop
refinement too early in those directions.  To avoid this, the original MATLAB toolbox also
introduces an explicit **degree of dimensional adaptivity** control that interpolates between
greedy and conservative refinement.

### Degree balancing

Version 5.1 of the MATLAB toolbox defines the *actual* adaptivity degree as the ratio

\[
r = \frac{n_{\text{adapt}}}{n_{\text{total}}}
\]

where $n_{\text{adapt}}$ counts points added by the adaptive (greedy) rule and
$n_{\text{total}}$ is the total number of sparse-grid points.

Given a target degree (configured by `DimadaptDegree` in MATLAB), the algorithm maintains
both adaptive and regular candidate index sets and chooses the next refinement source so that
$r$ stays close to the target value.  This **degree-balancing** strategy preserves the
benefits of adaptivity while injecting enough conservative refinement to improve robustness.

---

## Example: the Trid function

Consider the quadratic Trid function

\[
f(\mathbf{x}) = \sum_{i=1}^{d} (x_i - 1)^2 - \sum_{i=2}^{d} x_i\, x_{i-1}
\]

defined on the domain $[-d^2, d^2]^d$.  The function has a tridiagonal Hessian and clearly
exhibits **additive structure** — the coupling is only through nearest-neighbour terms
$x_i x_{i-1}$.

For $d = 100$, a traditional full tensor-product approach would require at least $2^{100}$
nodes.  The dimension-adaptive sparse grid detects the near-additive structure automatically
and recovers the function with $O(d^2)$ points:

```python
# Piecewise-linear CC grid, relative tolerance 0.1%
# (illustrative Python pseudocode)
d = 100
# ... dimension-adaptive loop using spgetseq_sp + SPGETSEQ_SP ...
# Result: ~27000 evaluations, estimated relative error < 0.03%
```

With the polynomial CGL grid, the quadratic function is recovered to **floating-point
accuracy** ($\approx 10^{-14}$ relative error) using only $\sim 20000$ evaluations.

---

## Sparse index sets

For high-dimensional problems, the sparse index data structure (sparse multi-index format)
is essential.  Instead of storing a dense $N_{\text{subgrids}} \times d$ matrix, only the
non-zero entries are kept:

| Array | Description |
|---|---|
| `indicesndiims(k)` | Number of active dimensions in subgrid $k$ |
| `indicesdims(a)` | Dimension index at packed address $a$ |
| `indiceslevs(a)` | Level at packed address $a$ |
| `indicesaddr(k)` | Start address in packed arrays for subgrid $k$ |
| `backwardneighbors(a)` | Index of backward-neighbour subgrid |
| `forwardneighbors(k,i)` | Forward-neighbour in dimension $i$ |

The `SPGETSEQ_SP` Fortran subroutine builds this structure; `SPSEQ2FULL` converts it back
to the dense representation for verification or export.

---

## Separability detection

A function $f(\mathbf{x}) = g_1(x_1) + g_2(x_2) + \cdots + g_d(x_d)$ is **fully
separable**.  The hierarchical surpluses for any subgrid $\mathbf{l}$ with two or more
active dimensions will be exactly zero for such a function, so the dimension-adaptive
algorithm concentrates all nodes on the 1-D subgrids.

More generally, for **nearly additive** functions, the surpluses decay rapidly for subgrids
with many active dimensions, and the algorithm automatically discovers this structure.
