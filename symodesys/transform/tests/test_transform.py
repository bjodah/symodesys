#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sympy
import numpy as np

from symodesys.transform import Transformer


def test_Transformer(tempdir=None):
    x, y = sympy.symbols('x y')
    t = Transformer([(x+1)**2,(x+1)**3-y],[x, y], tempdir=tempdir, save_temp=True)
    x_, y_ = np.linspace(0,10,10), np.linspace(10,20,10)
    out = t(x_, y_)
    assert np.allclose(out, np.vstack(((x_+1)**2, (x_+1)**3-y_)).transpose())

def test_Transformer__complex_argument_names(tempdir=None):
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
    dydt_in_z = y_in_z.diff(t)

    z_data = {t: num_t, z(t): num_z, z(t).diff(t): num_dzdt}

    exprs = [y_in_z, dydt_in_z]
    inp = [t, z(t), z(t).diff(t)]
    tfmr = Transformer(exprs, inp, tempdir=tempdir, save_temp=True)
    result = tfmr(*[z_data[k] for k in inp])

    num_y = result[:,0]
    assert np.allclose(num_y, exact_y)

    num_dydt = result[:,1]
    assert np.allclose(num_dydt, exact_dydt)


if __name__ == '__main__':
    test_Transformer__complex_argument_names('./tmp/')
    test_Transformer('./tmp/')
