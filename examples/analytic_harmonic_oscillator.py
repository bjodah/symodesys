#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import sympy
import numpy as np
import matplotlib.pyplot as plt

from symodesys.ivp import IVP
from symodesys.odesys import AnyOrderODESystem
from decay import Decay


def main(init_y, k_val, indep_var_lim, N = 0):
    """
    Integrate
    """
    x = sympy.Symbol('x', real = True)
    y = sympy.Function('y')(x)
    k = sympy.Symbol('k', real = True) # Spring constant

    harmonic_oscillator_eq = sympy.Eq(y.diff(x, 2), -k*y)

    odesys = AnyOrderODESystem.from_list_of_eqs([harmonic_oscillator_eq])

    fo_odesys = odesys.reduce_to_sys_of_first_order()

    param_vals_by_symb = {k: k_val}

    y0 = {y: init_y}
    t0, tend = indep_var_lim
    ivp = IVP(fo_odesys, y0, param_vals_by_symb, t0)

    # Attempt analytic reduction
    print ivp.recursive_analytic_reduction()

    ivp.integrate(tend, N = N)
    t, y = ivp.tout, ivp.yout

    ivp.plot(interpolate = True, show = False)
    plt.legend()
    plt.show()


if __name__ == '__main__':
    main(1.0, 1.0, [0.0, 10.0], N = 0)
