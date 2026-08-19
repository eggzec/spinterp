![spinterp](https://raw.githubusercontent.com/eggzec/spinterp/master/docs/assets/spinterp-banner.png)

# spinterp

**Sparse Grid Interpolation Toolbox for Python**

[![Tests](https://github.com/eggzec/spinterp/actions/workflows/test.yml/badge.svg)](https://github.com/eggzec/spinterp/actions/workflows/test.yml)
[![Documentation](https://github.com/eggzec/spinterp/actions/workflows/docs.yml/badge.svg)](https://github.com/eggzec/spinterp/actions/workflows/docs.yml)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

[![codecov](https://codecov.io/github/eggzec/spinterp/graph/badge.svg)](https://codecov.io/github/eggzec/spinterp)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=eggzec_spinterp&metric=alert_status)](https://sonarcloud.io/project/overview?id=eggzec_spinterp)
[![License](https://img.shields.io/badge/license-GPL%203.0-blue.svg)](./LICENSE)

[![PyPI Downloads](https://img.shields.io/pypi/dm/spinterp.svg?label=PyPI%20downloads)](https://pypi.org/project/spinterp/)
[![Python versions](https://img.shields.io/pypi/pyversions/spinterp.svg)](https://pypi.org/project/spinterp/)

Sparse Grid Interpolation Toolbox for Python.

## Grid types

| Grid type        | Basis              | Max depth |
|------------------|--------------------|-----------|
| Clenshaw-Curtis  | Piecewise linear   | 8         |
| Chebyshev        | Polynomial (DCT)   | 10        |
| Gauss-Patterson  | Nested Gaussian    | 6         |
| Maximum          | Piecewise linear   | 8         |
| NoBoundary       | Piecewise linear   | 8         |

## Routines

| Python binding  | Description                                    |
|-----------------|------------------------------------------------|
| `spvals`        | Compute hierarchical surpluses                 |
| `spinterp`      | Evaluate interpolant (and gradient) at points  |
| `spquad`        | Integrate interpolant over domain              |
| `spgrid`        | Return sparse grid node coordinates            |
| `spgetseq`      | Generate multi-index level sequences           |

## Quick example

```python
import spinterp

z = spinterp.spvals(lambda x, y, t: x**2 + y**2 - 2 * t, d=3)
f = spinterp.spinterp(z, 0.5, 0.2, 0.2)
q = spinterp.spquad(z)
```

## Installation

```bash
pip install spinterp
```

Requires Python 3.10+ and NumPy.

### Build from source

```bash
python bin/build.py install   # build and install via uv
python bin/build.py wheel     # produce a wheel in dist/
python bin/build.py clean     # remove all build artifacts
```

Requires gfortran, Meson, and f2py (`numpy`).

## Documentation

- [Theory](https://eggzec.github.io/spinterp/theory/) — hierarchical basis, Smolyak construction, algorithms
- [Quickstart](https://eggzec.github.io/spinterp/quickstart/) — runnable examples
- [API Reference](https://eggzec.github.io/spinterp/usage/) — function signatures and arguments

## License

GNU General Public License v3 (GPLv3) — see [LICENSE](LICENSE).
