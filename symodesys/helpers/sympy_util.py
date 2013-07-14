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


def reassign_const(expr, dest, known, source='C'):
    """
    e.g.
    >>> reassing_const(x*C1+C2, 'K', [x], 'C')
    x*K1+K2
    """
    def get_resymb(id_):
        tail = ''
        while sympy.Symbol(dest+id_+tail) in known:
            tail += 'A' # not pretty but should work..
        return sympy.Symbol(dest+id_+tail)

    new_symbs = get_new_symbs(expr, known)
    reassigned_symbs = []
    new_not_reassigned = []
    for symb in new_symbs:
        if source:
            # Only reassign matching source
            if symb.name.startswith(source):
                id_ = source.join(symb.name.split(source)[1:])
                resymb = get_resymb(id_)
                expr = expr.subs({symb: resymb})
                reassigned_symbs.append(resymb)
            else:
                # The new symb didn't match, store in
                # new_not_reassigned
                new_not_reassigned.append(symb)
        else:
            # All new symbs are to be renamed
            resymb = get_resymb(symb.name)
            expr.subs({symb: resymb})
            reassigned_symbs.append(resymb)
    return expr, reassigned_symbs, new_not_reassigned


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
