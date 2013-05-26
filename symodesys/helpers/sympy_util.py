#!/usr/bin/env python
# -*- coding: utf-8 -*-

# stdlib imports
from collections import defaultdict

# other imports
import sympy
import numpy as np
from sympy import C
from sympy.utilities.lambdify import implemented_function

from symodesys.helpers import cache

def get_new_symbs(expr, known_symbs):
    new_symbs = set()
    for atom in expr.atoms():
        if not atom in known_symbs and not atom.is_Number:
            new_symbs.add(atom)
    return new_symbs


def alt_ufuncify(args, expr, **kwargs):
    """
    This function mimics ufuncify in the autowrap module
    with the difference that all arguments are assumed to
    be arrays (of equal length)
    """
    from sympy.utilities.autowrap import autowrap
    y = C.IndexedBase(C.Dummy('y'))
    m = C.Dummy('m', integer=True)
    i = C.Dummy('i', integer=True)
    i = C.Idx(i, m)
    l = C.Lambda(args, expr)
    f = implemented_function('f', l)

    # first argument accepts an array
    internal_args = [C.IndexedBase(C.Dummy(a.name)) for a in args]
    return autowrap(C.Equality(y[i], f(*[a[i] for a in internal_args])), **kwargs)


@cache
def _array_subs_ufuncifier(expr, keys):
    """
    Special ufuncifies which supports sympy.Function and sympy.Derivative
    instances as arguments in expression

    returns ufunc-like callback and a tuple of keys which are
    used as function signature of ufunc

    This function is costly hence it is cached.
    """

    dummies = []
    def make_dummy(key, not_in):
        if isinstance(key, sympy.Function):
            candidate = sympy.Symbol(key.func.__name__, real=key.is_real)
        elif isinstance(key, sympy.Derivative):
            candidate = sympy.Symbol('d'+key.args[0].func.__name__+\
                                     ''.join('d'+x.name for x in key.args[1:]))
        elif isinstance(key, sympy.Symbol):
            candidate = key
        elif isinstance(key, str):
            candidate = sympy.Symbol(key)
        else:
            raise NotImplementedError
        # Now lets make sure the dummie symbol is unique
        if candidate not in dummies and not not_in.has(candidate):
            dummies.append(candidate)
            return candidate
        else:
            return make_dummy(candidate.name+'p', not_in)
    # First let's determine the used keys
    used_keys = []
    used_keys_symbs = []
    derivs = defaultdict(dict)
    for key in filter(expr.has, keys):
        used_keys.append(key)
        if isinstance(key, sympy.Symbol):
            used_keys_symbs.append(key)
        else:
            new_dummy = make_dummy(key, not_in=expr)
            used_keys_symbs.append(new_dummy)
            # Determine order of derivative
            if isinstance(key, sympy.Derivative):
                derivorder = len(key.args)-1
                derivs[derivorder][key] = new_dummy
            elif isinstance(key, sympy.Function):
                derivorder = 0
                derivs[derivorder][key] = new_dummy
            elif isinstance(key, sympy.Symbol):
                pass
            else:
                raise NotImplementedError
    for n in sorted(derivs.keys())[::-1]:
        # Substitutions of symbolified derivatives
        # need to be made for the highest degrees of
        # derivatives first
        for deriv, symb in derivs[n].items():
            expr = expr.subs({deriv: symb})

    # Now we are almost there..
    # But e.g. z(t) has included t into used_keys
    # we need to double check that t is truly used
    # now that z(t) has been turned into z
    truly_used_keys, truly_used_keys_symbs = [], []
    for k, ks in zip(used_keys, used_keys_symbs):
        if expr.has(ks):
            truly_used_keys.append(k)
            truly_used_keys_symbs.append(ks)

    ufunc = alt_ufuncify(tuple(truly_used_keys_symbs), expr)

    # autowrap is not called (from within alt_ufuncify) with args
    # specified, hence we need to sort alphabetically our used_keys
    # _based on_ alphabetic order of used_keys_symbs
    return ufunc, tuple(sorted(
        truly_used_keys,
        key=lambda x: truly_used_keys_symbs[truly_used_keys.index(x)]))

def array_subs(expr, array_dict):
    """
    E.g.
    >>> import numpy as np
    >>> from sympy import symbols
    >>> x, y = symbols('x y')
    >>> f = x+y**2
    >>> d = {x: np.arange(10), y: np.arange(10,20)}
    >>> array_subs(f, d)

    so far it looks like ufuncify, however, ufuncify fails for:

    >>> import numpy as np
    >>> from sympy import symbols
    >>> x, y = symbols('x y')
    >>> g = sympy.Function('g')(x)
    >>> f = x+y**2+g.diff(x)
    >>> d = {x: np.arange(10), y: np.arange(10,20), g.diff(x): np.arange(20,30)}
    >>> array_subs(f, d)

    which array_subs work around (to avoid limitations in ufuncify)
    by doing some variable substitions behind the scenes.
    """

    array0 = array_dict.values()[0]
    for array in array_dict.values():
        # No broadcasting supported
        assert array.shape == array0.shape
    from sympy.utilities.autowrap import ufuncify
    uf, used_keys = _array_subs_ufuncifier(expr, array_dict.keys())
    return uf(*(array_dict[k] for k in used_keys))

class SympyEvalr(object):
    """
    Specialized class to mimic the interface of IVP_Integrator
    class. It is used when an ODE has an analytic solution.

    Instances evaluates sympy expressions dependent on one variable and
    also substitutes parameters in expression using provided instance
    of FirstOrderODESystem
    The instances also evaluates the derivates in the independent variable
    up to requested order. (to facilitate interpolation)
    """

    default_dtype = np.float64

    def __init__(self, nderiv=0, dtype=None):
        """

        Arguments:
        - `exprs`: List of expressions to be evaluated
        - `indep_var_symb`: Sympy symbol of the independent variable
        - `params_by_symb`: Dictionary mapping sympy symbols of parameters to values
        - `nderiv`: set higher than 0 (default) in nderiv to also evaluate derivatives.
                   (output is sotred in self.Yout)
        """
        self.nderiv = nderiv
        if dtype == None: dtype = self.default_dtype
        self._dtype = dtype


    def configure(self, fo_odesys, param_vals):
        self._exprs = fo_odesys.solved_exprs
        self._indep_var_symb = fo_odesys.indepv
        self._params_by_symb = param_vals


    def eval_for_indep_array(self, arr, extra_params_by_symb):
        """
        Evaluate all expressions for values of indepedndent variable
        in array `arr` using self._params_symb and `extra_params_by_symb`
        for static substitution in sympy expressions in the list self._exprs
        """
        _Yout = np.empty((len(arr), len(self._exprs), self.nderiv+1), dtype=self._dtype)
        for expr_idx, expr in enumerate(self._exprs):
            for nderiv in range(self._nderiv):
                subs = self._params_by_symb
                subs.update(extra_params_by_symb)
                for i, t in enumerate(arr):
                    subs.update({self._indep_var_symb: t})
                    diff_expr = expr.diff(self._indep_var_symb, nderiv)
                    _Yout[i, expr_idx, nderiv] = diff_expr.subs(subs)
        self.Yout = _Yout
