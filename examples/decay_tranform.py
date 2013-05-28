#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sympy
import numpy as np
from decay import Decay
from symodesys.ivp import IVP

def main(y0={'u':1.0}, params={'lambda_u':1.0}, t0=0.0, tend=5.0, N=20):
    """
    Solve decay with variable substitution
    """
    d = Decay()
    ivp = IVP(d, y0, params, t0)
    u = d['u']
    lnu = ivp.mk_depv_symbol('lnu')
    ivp2 = ivp.use_internal_depv_trnsfm({lnu: sympy.log(u)}, {u: sympy.exp(lnu)})
    ivp2.integrate(tend, N, nderiv=2)

    t, u = ivp2.indep_out(), ivp2.trajectories()[d['u']]

    expected_u = np.vstack((
        y0['u']*np.exp(-params['lambda_u']*t),
        y0['u']*np.exp(-params['lambda_u']*t)*-params['lambda_u'],
        y0['u']*np.exp(-params['lambda_u']*t)*params['lambda_u']**2,
        )).transpose()
    assert np.allclose(u, expected_u)

    ivp2.plot(show=True)

if __name__ == '__main__':
    main()
