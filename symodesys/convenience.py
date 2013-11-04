from __future__ import print_function, division, absolute_import, unicode_literals

import numpy as np
import matplotlib.pyplot as plt

from symodesys.ivp import IVP
from symodesys.integrator import SciPy_IVP_Integrator

def _get_default_integrator(Integrator=SciPy_IVP_Integrator,
                            abstol=1e-8, reltol=1e-8, nderiv=2):
    return Integrator(abstol=abstol, reltol=reltol, nderiv=nderiv)


def plot_numeric_error(ODESys, depv_init, params, indepv_init,
                       indepv_end, N=0, integrator=None):
    """
    Convenience function for instantiating ODESys class and assigning
    it to an associated IVP instance, run the integration and in the
    case of ODESys having 'analytic_sol', plot the absolute and
    relative errors made during the integration.
    """
    integrator = integrator or _get_default_integrator()
    odesys = ODESys()
    with IVP(odesys, depv_init, params,
             indepv_init, integrator) as ivp:
        ivp.integrate(indepv_end, N)
        t, y = ivp.indepv_out(), ivp.trajectories()

        ivp.plot(interpolate = True, show = False, ax=plt.subplot(311))

        for i, (k, cb) in enumerate(odesys.analytic_sol.items()):
            analytic = cb(odesys, t, depv_init, params, indepv_init)
            plt.subplot(312)
            plt.plot(t, (y[odesys[k]][:, 0] - analytic) \
                     / ivp.integrator.abstol,
                     label = k+': abserr / abstol')
            plt.legend()
            plt.subplot(313)
            plt.plot(t, (y[odesys[k]][:, 0] - analytic) \
                     / analytic / ivp.integrator.reltol,
                     label = k+': relerr / reltol')
            plt.legend()
    if show: plt.show()

def plot_numeric_vs_analytic(ODESys, depv_init, params, indepv_init, indepv_end,
                             integrator=None, N=0, nderiv=2,
                             show=True, **kwargs):
    """
    Used to plot numeric solution of odesystem and analytic solution
    in same (sub)plot, to also show absolute and relative error, see
    plot_numeric_error

    kwargs are passed to ivp.plot()
    """
    show = kwargs.pop('show', False)

    integrator = integrator or _get_default_integrator()
    odesys = ODESys()
    with IVP(odesys, depv_init, params,
             indepv_init, integrator) as ivp:
        ivp.integrate(indepv_end, N)

        # Anlyse output
        plot_t = np.linspace(indepv_init, indepv_end, 50)

        ax = ivp.plot(interpolate = True, show=False, **kwargs)
        for depv, cb in odesys.analytic_sol.items():
            ax.plot(plot_t, cb(odesys, plot_t, depv_init, params, indepv_init),
                    label = 'Analytic {}'.format(depv))
        plt.legend()
    if show: plt.show()
    return ax

def assert_acceptable_numeric_errors(ODESys, depv_init, params,
                                     indepv_init, indepv_end,
                                     integrator=None, N=0, nderiv=2,
                                     abstol=1e-8, reltol=1e-8,
                                     acceptance_factor=100):

    integrator = integrator or _get_default_integrator(
        abstol=abstol, reltol=reltol)
    odesys = ODESys()
    with IVP(odesys, depv_init, params,
             indepv_init, integrator) as ivp:
        ivp.integrate(indepv_end, N)

        abserr, relerr = {}, {}
        for i, (depv, cb) in enumerate(odesys.analytic_sol.items()):
            analytic = cb(odesys, ivp.indepv_out(),
                          depv_init, params, indepv_init)
            numeric = ivp.trajectories()[odesys[depv]][:,0]
            abserr[depv] = numeric - analytic
            relerr[depv] = abserr[depv]/numeric

    for depv, ae in abserr.items():
        if not np.all(ae/abstol < acceptance_factor):
            print(ae/abstol)
            raise RuntimeError("Too large errors in {}".format(depv))

    for depv, re in relerr.items():
        if not np.allclose(re/reltol < acceptance_factor):
            print(re/reltol)
            raise RuntimeError("Too large errors in {}".format(depv))
