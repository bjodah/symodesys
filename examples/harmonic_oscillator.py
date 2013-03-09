#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import sympy
import numpy as np
import matplotlib.pyplot as plt

from symodesys.ivp import IVP
from symodesys.odesys import AnyOrderODESystem, FirstOrderODESystem

def HarmonicOscillator():
    """
    Returns AnyOrderODESystem instance (mimics a class instantiation)
    """
    x = sympy.Symbol('x', real = True)
    y = sympy.Function('y')(x)
    k = sympy.Symbol('k', real = True) # Spring constant

    harmonic_oscillator_eq = sympy.Eq(y.diff(x, 2), -k*y)
    return AnyOrderODESystem.from_list_of_eqs([harmonic_oscillator_eq])


def main(init_y, k_val, indep_var_lim, N = 0):
    """
    Integrate
    """

    odesys = HarmonicOscillator()

    odesys = odesys.reduce_to_sys_of_first_order()

    fo_odesys = FirstOrderODESystem(odesys) # FirstOrderODESystem is needed for IVP
    print fo_odesys.indep_var_symb
    param_vals_by_symb = {odesys['k']: k_val}

    # TODO:This part is cryptic, enhance!
    hlprs = odesys._1st_ordr_red_helper_fncs
    y0 = {k[2]: 0.0 for k in hlprs}
    # END TODO
    y0.update({odesys['y(x)']: init_y})
    t0, tend = indep_var_lim

    print fo_odesys.dydt(t0, [y0[k] for k in fo_odesys.all_depvar], [k_val])

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
