# Numerical Integration (Quadrature)

Once a sparse grid interpolant has been constructed, its integral over the domain can be
computed exactly (with respect to the interpolant) by taking the dot product of the
**quadrature weights** with the **hierarchical surpluses**.

For the unit domain $[0,1]^d$:

\[
\int_{[0,1]^d} A_{n,d}(f)(\mathbf{x})\,\mathrm{d}\mathbf{x}
= \sum_{\mathbf{l}} \sum_j w_j^{(\mathbf{l})}\, \alpha_j^{(\mathbf{l})}
= \mathbf{w}^\top \mathbf{z}
\]

where $\mathbf{w}$ is the vector of quadrature weights and $\mathbf{z}$ is the flat
surplus array.

```python
import spinterp, numpy as np

# Build quadrature weights for the CC grid
# (seq and z come from the standard surplus construction)
tp_all = sum(spinterp.spdim_cc(s) for s in all_seq)
w = spinterp.spquadw_cc(seq, tp_all)
integral = float(np.dot(w, z))
```

---

## Integration of regular sparse grid interpolants

### Test: product-type function in $d = 5$

Consider the integration test function

\[
f(\mathbf{x}) = \left(1 + \frac{1}{d}\right)^d \prod_{i=1}^{d} x_i^{1/d}
\]

over $[0,1]^5$.  The exact integral is $1$.  The table below reproduces results from
Table 1 of *Gerstner & Griebel, Numerical Algorithms 18(3–4), 1998*:

| Depth | CC points | CC error | Cheby points | Cheby error | GP points | GP error |
|---|---|---|---|---|---|---|
| 0 | 1 | 2.44e-01 | 1 | 2.44e-01 | 1 | 2.44e-01 |
| 1 | 11 | 1.08e+00 | 11 | 6.39e-01 | 11 | 8.94e-03 |
| 2 | 61 | 7.58e-02 | 61 | 1.44e-01 | 71 | 8.07e-04 |
| 3 | 241 | 2.86e-01 | 241 | 1.24e-01 | 351 | 2.07e-04 |
| 4 | 801 | 1.08e-01 | 801 | 6.65e-03 | 1471 | 2.26e-05 |
| 5 | 2433 | 8.00e-02 | 2433 | 1.06e-02 | 5503 | 1.42e-06 |
| 6 | 6993 | 5.03e-02 | 6993 | 1.74e-03 | 18943 | 3.44e-09 |

!!! note
    The grid types *Clenshaw-Curtis* and *Chebyshev* in this toolbox correspond to the
    "Trapez rule" and "Clenshaw-Curtis rule" sparse grids in Gerstner-Griebel (1998),
    respectively.

---

## Integration of dimension-adaptive interpolants

For high-dimensional problems, the dimension-adaptive approach can achieve dramatically
better accuracy with the same number of function evaluations compared to Monte Carlo.

### Absorption problem ($d = 20$)

Consider the absorption integral from Morokoff and Caflisch (1995).  Using the smooth
representation of the integrand with parameters $\gamma = 0.5$, $x = 0$, $d = 20$, the
exact solution is

\[
I = \frac{1}{\gamma} - \frac{1-\gamma}{\gamma}\,e^{\gamma(1-x)}
  \approx 0.3513.
\]

The dimension-adaptive Chebyshev (CGL) sparse grid converges orders of magnitude faster
than crude Monte Carlo with the same number of points:

| Points | CGL error | MC error (avg.) |
|---|---|---|
| 41 | 4.62e-04 | 1.30e-02 |
| 87 | 5.61e-06 | 6.01e-03 |
| 177 | 6.01e-07 | 7.79e-03 |
| 367 | 1.57e-07 | 3.44e-03 |
| 739 | 3.89e-08 | 4.76e-03 |
| 1531 | 2.46e-08 | 1.90e-03 |
| 3085 | 1.06e-09 | 1.04e-03 |
| 6181 | 2.75e-09 | 6.15e-04 |
| 24795 | 3.01e-10 | 5.42e-04 |
| 49739 | 1.79e-10 | 3.31e-04 |

The improvement is particularly pronounced here because the smooth integrand (second
representation) is well-suited to polynomial approximation, and its near-additive structure
is exploited by the dimension-adaptive algorithm.

---

## Quadrature weight functions

| Grid type | Weight function | Notes |
|---|---|---|
| Clenshaw-Curtis | `spquadw_cc(seq, nw)` | Simple uniform weights $1/\text{wval}$ |
| Chebyshev | `spquadw_cb(seq, nw, w1d, startid)` | Needs 1-D table from `cheb_weights` |
| Gauss-Patterson | `spquadw_gp(seq, nw, w1d, startid)` | Needs 1-D table from `gp_weights` |
| Maximum | `spquadw_m(seq, nw)` | Explicit 1-D weight array per dim |
| NoBoundary | `spquadw_nb(seq, nw)` | Endpoint weights doubled |

All weight functions also have sparse-index variants (`_sp` suffix) for high-dimensional
problems.
