#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

import os

import numpy as np

from sympy.abc import x, y, z

from symodesys.helpers.sympy_multi_ufuncify import multi_ufuncify, multi_autowrap


def test_multi_autowrap():
    from sympy.utilities.autowrap import autowrap
    exprs = [((x - y + z)**(13)).expand(),
             ((x / (1 + y + z))**(3)).expand()]
    binary_func = multi_autowrap(exprs)#, tempdir=os.path.join(os.path.dirname(__file__),'build/'))
    print(binary_func(1, 4, 2))
    #array([-1.0, 0.0029154518950437313])

def test_multi_ufucify():
    exprs = [((x - y + z)**(13)).expand(),
             ((x / (1 + y + z))**(3)).expand()]

    binary_func = multi_ufuncify([x,y,z], exprs, tempdir='build/ufunc/')
    x_ = np.linspace(1,2,3)
    y_ = np.linspace(4,5,3)
    z_ = np.linspace(2,3,3)
    print(binary_func(x_, y_, z_))


if __name__ == '__main__':
    test_multi_ufucify()
    test_multi_autowrap()
