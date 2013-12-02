#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import sympy
import numpy as np
import matplotlib.pyplot as plt

from symodesys.ivp import IVP
from symodesys.odesys import AnyOrderODESystem, FirstOrderODESystem
from harmonic_oscillator import HarmonicOscillator


"""
This example is work in progress...
"""

def main(init_y, init_dy, indep_var_lim, N = 0):
    """
    Solves a 1D harmonic oscillator analytically
    """

    odesys = HarmonicOscillator()
    y = odesys['y(x)']
    x = odesys['x']
    phi, A, omega = sympy.symbols('phi A omega')

    # Below it is obvious we want analytic reduction at very low level
    hypo_expr = A * sympy.sin(omega * x + phi)

    fo_odesys = odesys.reduce_to_sys_of_first_order()

    success = odesys.attempt_analytic_sol(
        y, hypo_expr, t0, y0, params, [A, omega, phi])

    assert success

    viewer = get_chaco_viewer(odesys, internal_y0, internal_params,
                              indep_var_lim[0], indep_var_lim[1], N)
    viewer.configure_traits()
    viewer.clean_up()


if __name__ == '__main__':
    main(1.0, 1.0, [0.0, 10.0], N = 0)
