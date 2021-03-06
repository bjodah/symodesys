from __future__ import print_function, division, absolute_import, unicode_literals

import numpy as np
import matplotlib.pyplot as plt

from symodesys.ivp import IVP
from symodesys.integrator import SciPy_IVP_Integrator

plotting_colors = 'k b r g m'.split()

def _get_default_integrator(Integrator=SciPy_IVP_Integrator,
                            abstol=1e-8, reltol=1e-8, nderiv=None):
    nderiv = nderiv or 1
    return Integrator(abstol=abstol, reltol=reltol, nderiv=nderiv)


def numeric_vs_analytic(ODESys, depv_init, params, indepv_init, indepv_end,
                        integrator=None, N=0, nderiv=None,
                        plot=True, subplot=True, show=True, acceptance_factor=0,
                        logger=None, interpolate=True, datapoints=False, **kwargs):
    """
    Used to plot and/or assert numeric solution in relation to analytic solution
    of odesystem
    In the case of plotting: in same (sub)plot, to also show absolute and relative error, see
    plot_numeric_error

    kwargs are passed to ivp.plot()
    """

    integrator = integrator or _get_default_integrator(nderiv=nderiv)
    odesys = ODESys()

    ls = ['-', '--']
    c = plotting_colors

    with IVP(odesys, depv_init, params,
             indepv_init, integrator, logger=logger) as ivp:
        ivp.integrate(indepv_end, N)

        # Analyse output
        if plot:
            plot_t = np.linspace(indepv_init, indepv_end, 50)
            if subplot:
                assert 'ax' not in kwargs
                ax = plt.subplot(311)
            else:
                ax = kwargs.pop('ax', None)
            ax = ivp.plot(interpolate=interpolate, show=False, ax=ax,
                          datapoints=datapoints,
                          nderiv=nderiv, ls=ls[:1], c=c, **kwargs)
        if acceptance_factor:
            abserr, relerr = {}, {}
        for i, depv in enumerate(odesys.all_depv):
            cb = odesys.analytic_sol[depv]
            analytic = cb(ivp.indepv_out(),
                          depv_init, params, indepv_init)
            numeric = ivp.trajectories(nderiv)[odesys[depv]][:,0]
            if plot:
                ax.plot(plot_t, cb(plot_t, depv_init,
                                   params, indepv_init),
                        label = 'Analytic {}'.format(depv),
                        c=c[i], ls=ls[1])
                if subplot:
                    plt.subplot(312)
                    plt.plot(ivp.indepv_out(),
                             (numeric-analytic)/integrator.abstol,
                             label='{}: abserr / abstol'.format(depv))
                    plt.legend()
                    plt.subplot(313)
                    plt.plot(ivp.indepv_out(),
                             (numeric-analytic)/analytic/integrator.reltol,
                               label='{}: relerr / reltol'.format(depv))
                    plt.legend()

            if acceptance_factor:
                abserr[depv] = numeric - analytic
                relerr[depv] = abserr[depv]/numeric

    if acceptance_factor:
        for depv, ae in abserr.items():
            if not np.all(ae/integrator.abstol < acceptance_factor):
                print(ae/integrator.abstol)
                raise RuntimeError("Too large errors in {}".format(depv))

        for depv, re in relerr.items():
            if not np.all(re/integrator.reltol < acceptance_factor):
                print(re/integrator.reltol)
                raise RuntimeError("Too large errors in {}".format(depv))

    if plot:
        if show: plt.show()
        return ax
    else:
        return
