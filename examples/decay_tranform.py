#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sympy

from decay import Decay
from symodesys.ivp import IVP

def main(y0={'u':1.0}, params={'lambda_u':1.0}, t0=0.0, tend=5.0, N=20):
    """
    """
    d = Decay()
    ivp = IVP(d, y0, params, t0)
    u = d['u']
    lnu = ivp.mk_depv_symbol('lnu')
    ivp2 = ivp.use_internal_depv_trnsfm({lnu: sympy.log(u)}, {u: sympy.exp(lnu)})
    ivp2.integrate(tend, N, order=0)
    ivp2.plot(show=True)

if __name__ == '__main__':
    main()
