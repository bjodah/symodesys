#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import sympy
import numpy as np
import matplotlib.pyplot as plt

from symodesys.odesys import AnyOrderODESystem, FirstOrderODESystem
from symodesys.gui import get_chaco_viewer

def HarmonicOscillator():
    """
    Returns AnyOrderODESystem instance (mimics a class instantiation)
    of a second order system of differential equation(s).
    """
    x = sympy.Symbol('x', real = True)
    y = sympy.Function('y')(x)
    k = sympy.Symbol('k', real = True) # Spring constant

    harmonic_oscillator_eq = sympy.Eq(y.diff(x, 2), -k*y)
    return AnyOrderODESystem.from_list_of_eqs([harmonic_oscillator_eq])


def main(y0, t0, tend, params, N = 20):
    """
    Present a GUI with interactive control of initial values and
    parameters.
    """
    odesys = HarmonicOscillator()
    fo_odesys = odesys.reduce_to_sys_of_first_order(y0, 0.0)
    print(fo_odesys.all_depv)
    viewer = get_chaco_viewer(fo_odesys, y0, params, t0, tend, N)
    viewer.configure_traits()
    viewer.clean_up()


if __name__ == '__main__':
    main({'y': 1.0}, 0.0, 10.0, {'k': 1.0})
