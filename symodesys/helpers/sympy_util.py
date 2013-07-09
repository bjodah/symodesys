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

def reassign_const(expr, dest, known, source='C'):
    """
    e.g.
    >>> reassing_const(x*C1+C2, 'K', [x])
    x*K1+K2
    """
    def get_resymb(id_):
        tail = ''
        while sympy.Symbol(dest+id_+tail) in known:
            tail += 'A' # not pretty but should work..
        return sympy.Symbol(dest+id_+tail)

    new_symbs = get_new_symbs(expr, known)
    reassigned_symbs = set()
    new_not_reassigned = set()
    for symb in new_symbs:
        if source:
            # Only reassign matching source
            if symb.name.startswith(source):
                id_ = source.join(symb.name.split(source)[1:])
                resymb = get_resymb(id_)
                expr = expr.subs({symb: resymb})
                reassigned_symbs.add(resymb)
            else:
                # The new symb didn't match, store in
                # new_not_reassigned
                new_not_reassigned.add(symb)
        else:
            # All new symbs are to be renamed
            resymb = get_resymb(symb.name)
            expr.subs({symb: resymb})
            reassigned_symbs.append(resymb)
    return expr, reassigned_symbs, new_not_reassigned

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

    # all arguments accepts arrays
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




def _clean_args_from_piecewise(inexpr, undefined):
    """
    See docstring of `get_without_piecewise`
    """
    new_args = []
    for arg in inexpr.args:
        if isinstance(arg, sympy.Piecewise):
            # Found it!.. now get `True` condition
            for expr, cond in arg.args:
                if cond == True:
                    new_args.append(expr)
                else:
                    undefined.append(cond)
        else:
            if arg.has(sympy.Piecewise):
                # Recursive call!
                new_args.append(_clean_args_from_piecewise(
                    arg, undefined))
            else:
                new_args.append(arg)
    return inexpr.fromiter(new_args)

def get_without_piecewise(inexpr):
    """
    Returns expression without Piecewise elements
    and a list of conditions which were neglected
    (and hence are to be considered undefined from
    in the new expression)

    Arguments:
    - `inexpr`: the expression containing Piecewise terms

    Returns:
    A tuple of lenght 2, with new expression as first
    item and a list of neglected conditions as second.
    """

    undefined = []
    newexpr = _clean_args_from_piecewise(inexpr, undefined)
    return newexpr, undefined
