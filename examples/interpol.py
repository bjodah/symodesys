#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import sympy
import numpy as np
import matplotlib.pyplot as plt

from symodesys.odesys import SimpleFirstOrderODESystem
from symodesys.ivp import IVP
from symodesys.convenience import plot_numeric_vs_analytic

class X5(SimpleFirstOrderODESystem):

    depv_tokens = 'y',

    @property
    def expressions(self):
        return {self['y']: self.indepv ** 4}

    def analytic_y(self, indep_vals, y0, params, t0):
        return y0['y'] + indep_vals ** 5 / 5 - t0 ** 5

    analytic_sol = {'y': analytic_y}


if __name__ == '__main__':
    for i in range(3):
        ax=plt.subplot(310+i+1)
        ax=plot_numeric_vs_analytic(
            ODESys = X5,
            y0 = {'y': 1.0},
            params = {},
            t0 = 0,
            tend = 3.0,
            N = 2,
            nderiv=i,
            datapoints=True,
            ax=ax)
        if i > 0:
            box = ax.get_position()
            # Shrink height with 20%
            ax.set_position([box.x0, box.y0, box.width, box.height * 0.8])

        ax.legend(loc='upper left')
        plt.title('Using {} derivatives per point'.format(i))

    plt.show()
