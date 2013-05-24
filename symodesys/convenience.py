import numpy as np
import matplotlib.pyplot as plt

from symodesys.ivp import IVP

def plot_numeric_vs_analytic(ODESys, y0, params, t0, tend, integrator=None, N=0):
    """
    Run compiled integrator (pointless for such a small ODE but
    useful for large systems)
    """
    odesys = ODESys()

    # Solve
    with IVP(odesys, y0, params, t0, integrator=integrator) as ivp:
        ivp.integrate(tend, N=N, nderiv=2)

        # Anlyse output
        plot_t = np.linspace(t0, tend, 50)

        ivp.plot(interpolate = True, show = False)
        for depv, cb in odesys.analytic_sol.items():
            plt.plot(plot_t, cb(odesys, plot_t, y0, params), label = 'Analytic {}'.format(depv))
        plt.legend()
        plt.show()
