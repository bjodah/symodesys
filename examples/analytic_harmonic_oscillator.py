#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import sympy
import numpy as np
import matplotlib.pyplot as plt

from symodesys.ivp import IVP
from symodesys.odesys import AnyOrderODESystem, FirstOrderODESystem
from harmonic_oscillator import HarmonicOscillator


def main(init_y, k_val, indep_var_lim, N = 0):
    """
    Solves a 1D harmonic oscillator analytically
    """

    odesys = HarmonicOscillator()
    y = odesys['y(x)']
    x = odesys['x']
    phi, A, omega = sympy.symbols('phi A omega')

    # Below it is obvious we want analytic reduction at very low level
    solved_odesys = odesys.attempt_analytic_sol(
        y, A * sympy.sin(omega * x + phi), [A, omega, phi])

    # TODO: FIX EVERYTHING BELOW



if __name__ == '__main__':
    main(1.0, 1.0, [0.0, 10.0], N = 0)
