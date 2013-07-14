#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import sympy

from symodesys.helpers import array_subs, get_without_piecewise, reassign_const


def test_get_without_piecewise():
    x, k = sympy.symbols('x k')
    f = sympy.Function('f')
    dfdx_expr = x*sympy.exp(-k*x)
    f_analytic = -(k*x+1)*sympy.exp(-k*x)/k**2
    sol = sympy.dsolve(f(x).diff(x)-dfdx_expr,f(x))
    without_piecewise, undefined = get_without_piecewise(sol.rhs)
    # k != 0 (default)
    assert (without_piecewise - f_analytic).simplify() == \
        sympy.Symbol('C1')
    # k == 0
    assert undefined[0].rhs == 0
    assert undefined[0].lhs.subs({k:0}) == 0

def test_reassign_const():
    x, C1, C2 = sympy.symbols('x C1 C2')
    new_expr, rea, not_rea = reassign_const(x*C1+C2, 'K', [x], 'C')
    assert not_rea == set([])
    K1, K2 = sympy.symbols('K1 K2')
    assert new_expr - x*K1 - K2 == 0
    assert rea == set([K1, K2])

if __name__ == '__main__':
    test_reassign_const()
    test_get_without_piecewise()
    test_array_subs___1()
    test_array_subs___2()
