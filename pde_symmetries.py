from __future__ import annotations
from dataclasses import dataclass
from itertools import permutations
from typing import List, Tuple, Sequence, Set, Optional, Dict
import sympy as sp

def _inverse_perm(p: Sequence[int]) -> List[int]:
    q = [0] * len(p); 
    for i, pi in enumerate(p): q[pi] = i
    return q

def _identity_perm(n: int) -> List[int]:
    return list(range(n))

def _as_set(idxs): 
    return set([] if idxs is None else idxs)

@dataclass(frozen=True)
class Transform:
    dep_perm: Tuple[int, ...]
    ind_perm: Tuple[int, ...]
    dep_sign: Tuple[int, ...] = ()
    ind_sign: Tuple[int, ...] = ()

class PDESymmetrySystem:
    def __init__(self, dep_funcs: Sequence[sp.Function], ind_vars: Sequence[sp.Symbol], equations: Sequence[sp.Expr]) -> None:
        self.dep_funcs = list(dep_funcs)
        self.ind_vars = list(ind_vars)
        self.equations = [sp.simplify(e) for e in equations]
        if not self.dep_funcs or not self.ind_vars:
            raise ValueError("dep_funcs and ind_vars must be non-empty.")
        for f in self.dep_funcs:
            if not isinstance(f, (sp.FunctionClass, sp.Function)):
                raise TypeError(f"Dependent entry {f} must be a sympy.Function.")
        for e in self.equations:
            if not isinstance(e, sp.Expr):
                raise TypeError("All equations must be sympy expressions.")
        self._dep_class_to_index = {f: i for i, f in enumerate(self.dep_funcs)}
        self._ind_index = {v: i for i, v in enumerate(self.ind_vars)}

    def transform_equations(self, T: Transform):
        dep_perm = list(T.dep_perm); ind_perm = list(T.ind_perm)
        dep_sign_set = _as_set(T.dep_sign); ind_sign_set = _as_set(T.ind_sign)
        if sorted(dep_perm) != list(range(len(self.dep_funcs))): raise ValueError("dep_perm invalid")
        if sorted(ind_perm) != list(range(len(self.ind_vars))):  raise ValueError("ind_perm invalid")
        inv_ind = _inverse_perm(ind_perm); m = len(self.ind_vars)
        ind_sign = [(-1 if i in ind_sign_set else 1) for i in range(m)]
        xprime = [ind_sign[i] * self.ind_vars[ind_perm[i]] for i in range(m)]
        bare_var_map = { self.ind_vars[j]: ind_sign[inv_ind[j]] * self.ind_vars[inv_ind[j]] for j in range(m)}

        def func_replacer(call: sp.Function) -> sp.Expr:
            fcls = call.func
            if fcls not in self._dep_class_to_index:
                new_args = [transform_expr(a) for a in call.args]
                return fcls(*new_args)
            i = self._dep_class_to_index[fcls]
            f_new = self.dep_funcs[dep_perm[i]]
            new_call = f_new(*xprime)
            if i in dep_sign_set: new_call = -new_call
            return new_call

        def transform_expr(e: sp.Expr) -> sp.Expr:
            if e.is_Number or e.is_NumberSymbol: return e
            if e.is_Symbol: return bare_var_map.get(e, e)
            if isinstance(e, sp.AppliedUndef): return func_replacer(e)
            if e.is_Derivative:
                base = transform_expr(e.expr); coeff = 1; new_vars = []
                for v in e.variables:
                    if v not in self._ind_index: new_vars.append(v); continue
                    j = self._ind_index[v]; k = inv_ind[j]
                    if ind_sign[k] == -1: coeff *= -1
                    new_vars.append(self.ind_vars[k])
                return coeff * sp.Derivative(base, *new_vars)
            if hasattr(e, 'args') and e.args:
                new_args = [transform_expr(a) for a in e.args]
                try: return e.func(*new_args)
                except Exception: return e.xreplace({a: na for a, na in zip(e.args, new_args)})
            return e

        return [sp.simplify(transform_expr(eq)) for eq in self.equations]

    def is_symmetry(self, T: Transform, simplify: bool=True) -> bool:
        transformed = self.transform_equations(T)
        original = self.equations
        if simplify:
            transformed = [sp.simplify(e) for e in transformed]
            original = [sp.simplify(e) for e in original]
        used = [False]*len(original)
        for te in transformed:
            matched = False
            for i, oe in enumerate(original):
                if used[i]: continue
                if sp.simplify(te - oe) == 0 or sp.simplify(te + oe) == 0:
                    used[i] = True; matched = True; break
            if not matched: return False
        return all(used)

    def find_dep_permutation_symmetries(self):
        k = len(self.dep_funcs); m = len(self.ind_vars)
        out = []
        for p in permutations(range(k)):
            T = Transform(dep_perm=tuple(p), ind_perm=tuple(_identity_perm(m)))
            if self.is_symmetry(T): out.append(T)
        return out

    def find_full_symmetries(self, search_dep_perm=True, search_ind_perm=True, allow_dep_parity=True, allow_ind_parity=True, max_results: Optional[int]=None):
        def all_parity(n: int):
            return [tuple(i for i in range(n) if (mask>>i)&1) for mask in range(1<<n)]
        k = len(self.dep_funcs); m = len(self.ind_vars)
        dep_perm_space = [tuple(_identity_perm(k))] if not search_dep_perm else list(permutations(range(k)))
        ind_perm_space = [tuple(_identity_perm(m))] if not search_ind_perm else list(permutations(range(m)))
        dep_parity_space = [()] if not allow_dep_parity else all_parity(k)
        ind_parity_space = [()] if not allow_ind_parity else all_parity(m)
        results = []
        for dp in dep_perm_space:
            for ip in ind_perm_space:
                for ds in dep_parity_space:
                    for isg in ind_parity_space:
                        T = Transform(dp, ip, ds, isg)
                        if self.is_symmetry(T):
                            results.append(T)
                            if max_results and len(results) >= max_results: return results
        return results

    def symmetrize(self, expr: sp.Expr, group):
        if not group: raise ValueError("Empty group")
        total = 0
        for T in group:
            # rebuild variable maps
            dep_perm = list(T.dep_perm); ind_perm = list(T.ind_perm)
            inv_ind = _inverse_perm(ind_perm); m = len(self.ind_vars)
            ind_sign = [(-1 if i in set(T.ind_sign) else 1) for i in range(m)]
            xprime = [ind_sign[i] * self.ind_vars[ind_perm[i]] for i in range(m)]
            bare_var_map = { self.ind_vars[j]: ind_sign[inv_ind[j]] * self.ind_vars[inv_ind[j]] for j in range(m)}

            def func_replacer(call: sp.Function) -> sp.Expr:
                fcls = call.func
                if fcls not in self._dep_class_to_index:
                    new_args = [transform_expr(a) for a in call.args]
                    return fcls(*new_args)
                i = self._dep_class_to_index[fcls]
                f_new = self.dep_funcs[dep_perm[i]]
                new_call = f_new(*xprime)
                if i in set(T.dep_sign): new_call = -new_call
                return new_call

            def transform_expr(e: sp.Expr) -> sp.Expr:
                if e.is_Number or e.is_NumberSymbol: return e
                if e.is_Symbol: return bare_var_map.get(e, e)
                if isinstance(e, sp.AppliedUndef): return func_replacer(e)
                if e.is_Derivative:
                    base = transform_expr(e.expr); coeff = 1; new_vars = []
                    for v in e.variables:
                        if v not in self._ind_index: new_vars.append(v); continue
                        j = self._ind_index[v]; k = inv_ind[j]
                        if ind_sign[k] == -1: coeff *= -1
                        new_vars.append(self.ind_vars[k])
                    return coeff * sp.Derivative(base, *new_vars)
                if hasattr(e, 'args') and e.args:
                    new_args = [transform_expr(a) for a in e.args]
                    try: return e.func(*new_args)
                    except Exception: return e.xreplace({a: na for a, na in zip(e.args, new_args)})
                return e

            total += transform_expr(expr)
        return sp.simplify(total/len(group))

if __name__ == "__main__":
    x, y = sp.symbols('x y')
    u, v = sp.Function('u'), sp.Function('v')
    eqs = [sp.diff(u(x,y), x) + sp.diff(u(x,y), y),
           sp.diff(v(x,y), x) + sp.diff(v(x,y), y)]
    sys = PDESymmetrySystem([u,v],[x,y],eqs)
    G = sys.find_full_symmetries()
    print("Found symmetries:", len(G))
    beta = sys.symmetrize(u(x,y)**2, G) if G else None
    if beta is not None: print("Symmetrized u^2:", beta)
