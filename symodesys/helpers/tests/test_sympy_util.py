#!/usr/bin/env python
# -*- coding: utf-8 -*-

from symodesys.helpers import array_subs

def test_array_subs():
    import numpy as np
    import sympy
    x, y, z = sympy.symbols('x y z')
    g = sympy.Function('g')(x)
    f = x+y**2+g+g.diff(x)
    d = {x: np.arange(10), y: np.arange(10,20), g.diff(x): np.arange(20,30),
         g: np.arange(40,50), z: np.arange(10)}
    assert np.allclose(array_subs(f, d), np.array(
        [160, 184, 210, 238, 268, 300, 334, 370, 408, 448]))


if __name__ == '__main__':
    test_array_subs()
