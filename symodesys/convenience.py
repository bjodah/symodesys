import numpy as np
import matplotlib.pyplot as plt

from symodesys.ivp import IVP


def plot_numeric_error(ODESys, y0, params, t0, tend, integrator=None, N=0):
    """
    Convenience function for instantiating ODESys class and assigning
    it to an associated IVP instance, run the integration and in the case
    of ODESys having 'analytic_sol', plot the absolute and relative errors
    made during the integration.
    """
    odesys = ODESys()
    with IVP(odesys, y0, params, t0, integrator=integrator) as ivp:
        ivp.integrate(tend, N)
        t, y = ivp.indep_out(), ivp.trajectories()

        ivp.plot(interpolate = True, show = False, ax=plt.subplot(311))

        for i, (k, cb) in enumerate(odesys.analytic_sol.items()):
            analytic = cb(odesys, t, y0, params, t0)
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
    plt.show()

def plot_numeric_vs_analytic(ODESys, y0, params, t0, tend,
                             integrator=None, N=0, nderiv=2,
                             **kwargs):
    """
    Used to plot numeric solution of odesystem and analytic solution
    in same (sub)plot, to also show absolute and relative error, see
    plot_numeric_error
    """
    odesys = ODESys()

    # Solve
    with IVP(odesys, y0, params, t0, integrator=integrator) as ivp:
        ivp.integrate(tend, N=N, nderiv=nderiv)

        # Anlyse output
        plot_t = np.linspace(t0, tend, 50)

        ax = ivp.plot(interpolate = True, show = kwargs.pop('show', False), **kwargs)
        for depv, cb in odesys.analytic_sol.items():
            ax.plot(plot_t, cb(odesys, plot_t, y0, params, t0), label = 'Analytic {}'.format(depv))
    return ax
