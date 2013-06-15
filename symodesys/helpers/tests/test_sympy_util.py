#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import sympy

from symodesys.helpers import array_subs, get_without_piecewise, reassign_const

def test_array_subs___1():
    # Example that ufuncify would not be able to treat.
    x, y, z = sympy.symbols('x y z')
    g = sympy.Function('g')(x)
    f = x+y**2+g+g.diff(x)
    d = {x: np.arange(10), y: np.arange(10,20), g.diff(x): np.arange(20,30),
         g: np.arange(40,50), z: np.arange(10)}
    assert np.allclose(array_subs(f, d), np.array(
        [160, 184, 210, 238, 268, 300, 334, 370, 408, 448]))

def test_array_subs___2():
    # y'(t) = t**4; y(0)=1.0  ===> y = t**5/5 + 1.0
    # z(t) = log(y(t)) ===> y=exp(z(t)), dzdt = t**4*exp(-z(t))

    num_t = np.linspace(0.0,3.0,4)
    exact_y = num_t**5/5.0+1.0
    exact_dydt = num_t**4

    num_z = np.log(exact_y)
    num_dzdt = num_t**4*np.exp(-num_z)

    z = sympy.Function('z')
    y = sympy.Function('y')
    t = sympy.Symbol('t')

    y_in_z = sympy.exp(z(t))
    z_in_y = sympy.log(y(t))

    z_data = {t: num_t, z(t): num_z, z(t).diff(t): num_dzdt}

    num_y = array_subs(y_in_z, z_data)
    assert np.allclose(num_y, exact_y)

    dydt_in_z = y_in_z.diff(t)
    num_dydt = array_subs(dydt_in_z, z_data)
    assert np.allclose(num_dydt, exact_dydt)


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
    assert not_rea == []
    K1, K2 = sympy.symbols('K1 K2')
    assert new_expr - x*K1 - K2 == 0
    assert rea == [K1, K2]

if __name__ == '__main__':
    test_reassign_const()
    test_get_without_piecewise()
    test_array_subs___1()
    test_array_subs___2()
