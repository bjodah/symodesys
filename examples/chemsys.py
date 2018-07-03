#!/usr/bin/env python
# -*- coding: utf-8 -*-

from collections import defaultdict

import sys


from symodesys.gsl import GSL_IVP_Integrator as Integrator
from symodesys.odesys import SimpleFirstOrderODESystem
from symodesys.ivp import IVP


class ChemSys(SimpleFirstOrderODESystem):

    depv_tokens = 'A AcCN AcLB B C LB P Ti Ti2 TiA TiCN2 TiCNTiCN TiOAc2'.split()
    param_tokens = 'k1 k2f k2b k3 k4 k5f k5b k6f k6b k7 k8'.split()

    A_low, A_high, A_default = 0, 1, 1
    AcCN_low, AcCN_high, AcCN_default = 0, 2, 2
    LB_low, LB_high, LB_default = 0, 0.1, 0.1
    Ti2_low, Ti2_high, Ti2_default = 0, 0.05, 0.05
    AcLB_default = None
    B_default = None
    C_default = None
    P_default = None
    Ti_default = None
    TiA_default = None
    TiCN2_default = None
    TiCNTiCN_default = None
    TiOAc2_default = None

    k1_low, k1_high, k1_default = 50., 5000., 500.
    k2f_low, k2f_high, k2f_default = 1e5, 1e7, 1e6
    k2b_low, k2b_high, k2b_default = 2e4, 2e6, 2e5
    k3_low, k3_high, k3_default = 1., 100., 10.
    k4_low, k4_high, k4_default = 2e3, 2e5, 2e4
    k5f_low, k5f_high, k5f_default = 8., 800., 80.
    k5b_low, k5b_high, k5b_default = 0.8, 80., 8.
    k6f_low, k6f_high, k6f_default = 3e4, 3e6, 3e5
    k6b_low, k6b_high, k6b_default = 7e3, 7e5, 7e4
    k7_low, k7_high, k7_default = 3e2, 3e4, 3e3
    k8_low, k8_high, k8_default = 1e2, 1e4, 1e3

    @property
    def expressions(self):
        A, AcCN, AcLB, B, C, LB, P, Ti, Ti2, TiA, TiCN2, TiCNTiCN, TiOAc2 = [
            self[x] for x in self.depv_tokens]
        k1, k2f, k2b, k3, k4, k5f, k5b, k6f, k6b, k7, k8 = [
            self[x] for x in self.param_tokens]
        r1 = k1 * Ti2 * AcCN**2
        r2f = B * TiOAc2
        r2b = C * Ti2 * AcCN
        r3 = k3 * AcCN * LB
        r4 = k4 * AcLB * B
        r5f = k5f * Ti * A
        r5b = k5b * TiA
        r6f = k6f * Ti2
        r6b = k6b * Ti * Ti
        r7 = k7 * TiA * TiCN2
        r8 = k8 * TiCNTiCN * A
        return {
            A: -r8 - r5f + r5b,
            AcCN: 2*r2b - 2*r2f -2*r1 -r3,
            AcLB: r3 - r4,
            B: r8 + r7 - r4 - 2*r2f + 2*r2b,
            C: r2f - r2b,
            LB: r4 - r3,
            P: r4,
            Ti: 2*r6f - 2*r6b - r5f + r5b,
            Ti2: r6b - r6f + 2*r2f - 2*r2b - r1,
            TiA: r5f - r5b - r7,
            TiCN2: r1 - r7,
            TiCNTiCN: r4 - r8,
            TiOAc2: r1 - r2f + r2b,
        }


def main(y0, params, tend, t0 = 0.0, N = 50):
    """
    Example program integrating an IVP problem of van der Pol oscillator
    default is adaptive step size (N=0)
    """
    cs = ChemSys()
    ivp = IVP(cs, y0, params, t0, integrator=Integrator(
        save_temp=True, tempdir='build_chemsys/')
    )

    ivp.integrate(tend, N)
    ivp.plot(interpolate=False, datapoints=True, show=True)

y0 = {}
y0['A'] = 1.
y0['AcCN'] = 2.
y0['LB'] = 0.1
y0['Ti2'] = 0.05

for k in 'AcLB B C P Ti TiA TiCN2 TiCNTiCN TiOAc2'.split():
    y0[k] = 0.0

params = dict(k1=500, k2f=1e6, k2b=2e5, k3=10, k4=2e4, k5f=80,
              k5b=8, k6f=3e5, k6b=7e4, k7=3e3, k8=1e3)

if __name__ == '__main__':
    if len(sys.argv) > 1:
        mu = float(sys.argv[1])
    else:
        mu = 1.0

    main(y0=y0, params=params, tend=3.0)
