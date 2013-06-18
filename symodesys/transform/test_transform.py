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


if __name__ == '__main__':
    test_Transformer('./tmp/')
