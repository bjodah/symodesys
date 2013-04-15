import numpy as np
import matplotlib.pyplot as plt

from symodesys.ivp import IVP

def plot_numeric_vs_analytic(ODESys, y0, params, t0, tend, N = 0, ivp_kwargs=None):
    """
    Run compiled integrator (pointless for such a small ODE but
    useful for large systems)
    """
    ivp_kwargs = ivp_kwargs or {}
    odesys = ODESys()

    # Solve
    ivp = IVP(odesys, y0, params, t0, **ivp_kwargs)
    ivp.integrate(tend, N=N, order=2)

    # Anlyse output
    plot_t = np.linspace(t0, tend, 50)

    ivp.plot(interpolate = True, show = False)
    for depv, cb in odesys.analytic_sol.items():
        plt.plot(plot_t, cb(odesys, plot_t, y0, params), label = 'Analytic {}'.format(depv))
    plt.legend()
    plt.show()

