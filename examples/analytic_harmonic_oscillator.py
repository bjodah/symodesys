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
    # Below it is obvious we want analytic reduction at very low level
    solved_odesys = odesys.attempt_analytic_reduction()

    # TODO: FIX EVERYTHING BELOW
    ########################################
    odesys = odesys.reduce_to_sys_of_first_order()
    hlprs = odesys._1st_ordr_red_helper_fncs

    fo_odesys = FirstOrderODESystem(odesys) # FirstOrderODESystem is needed for IVP
    print fo_odesys.indep_var_symb
    param_vals_by_symb = {k: k_val}

    y0 = {k[2]: 0.0 for k in hlprs}
    y0.update({y: init_y})
    t0, tend = indep_var_lim

    print fo_odesys.dydt(t0, [y0[k] for k in fo_odesys.dep_var_func_symbs], [k_val])

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
