# PDE Symmetry Tools

Symbolic detection of **permutation** and **parity** symmetries for systems of evolution equations (PDE/ODE), plus **group-based symmetrization** of expressions.

This repository is structured as:
- `README.md`
- `requirements.txt`
- `LICENSE`
- `examples/` with runnable demos and a `figures/` subfolder

## Features

- Permutation and parity transforms on **independent** and **dependent** variables.
- Symmetry checking up to equation **reordering** and **per-equation sign**.
- Brute-force search for symmetry groups (small systems).
- Symmetrization (orbit-averaging) of monomials/polynomials.

## Quickstart

```bash
python -m venv .venv
source .venv/bin/activate  # on Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

Run an example:
```bash
python examples/run_demo.py
```

## Minimal API

```python
from pde_symmetries import PDESymmetrySystem, Transform
import sympy as sp

x, y = sp.symbols('x y')
u, v = sp.Function('u'), sp.Function('v')

eqs = [
    sp.diff(u(x,y), x) + sp.diff(u(x,y), y),
    sp.diff(v(x,y), x) + sp.diff(v(x,y), y),
]

sys = PDESymmetrySystem([u,v], [x,y], eqs)

# Check a transform
T = Transform(dep_perm=(1,0), ind_perm=(1,0), dep_sign=(), ind_sign=())
print(sys.is_symmetry(T))

# Enumerate symmetries
group = sys.find_full_symmetries()
print("Group size:", len(group))
```



## Citation

If you use this code, please cite:
> Moataz M. Alghamdi, **Symbolic Detection of Permutation and Parity Symmetries of Evolution Equations**, KAUST MSc Thesis, 2017.
