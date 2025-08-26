# examples/run_demo.py
import sympy as sp
from pde_symmetries import PDESymmetrySystem, Transform

x, y = sp.symbols('x y')
u, v = sp.Function('u'), sp.Function('v')

eqs = [
    sp.diff(u(x,y), x) + sp.diff(u(x,y), y),
    sp.diff(v(x,y), x) + sp.diff(v(x,y), y),
]

sys = PDESymmetrySystem([u, v], [x, y], eqs)

print("== Dependent permutation symmetries ==")
for T in sys.find_dep_permutation_symmetries():
    print(T)

print("\n== Full symmetries (small system) ==")
G = sys.find_full_symmetries()
print("Group size:", len(G))

expr = u(x,y)**2
beta = sys.symmetrize(expr, G) if G else expr
print("Symmetrized u(x,y)^2 ->", sp.simplify(beta))

# Example transform check
T = Transform(dep_perm=(1,0), ind_perm=(1,0), dep_sign=(), ind_sign=())
print("Is swap(u,v) & swap(x,y) a symmetry? ->", sys.is_symmetry(T))
