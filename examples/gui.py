#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np

from enthought.chaco.api import ArrayPlotData, Plot
from enthought.enable.component_editor import ComponentEditor
from enthought.traits.api import HasTraits, Instance, Array, Property, Range, Float, Enum
from enthought.traits.ui.api import Item, View

from van_der_pol import VanDerPolOscillator
from symodesys.ivp import IVP

class ODESolViewer(HasTraits):

    plot = Instance(Plot)
    plotdata = Instance(ArrayPlotData, args=())

    t = Array
    u = Property(Array, depends_on = ['u0', 'v0', 'mu'])
    v = Property(Array, depends_on = ['u0', 'v0', 'mu'])
    u0 = Range(low = 0.0, high = 10.0, value = 1.0)
    v0 = Range(low = 0.0, high = 10.0, value = 1.0)
    mu = Range(low = 0.0, high = 10.0, value = 3.0)
    plot_type = Enum("line", "scatter")

    traits_view = View(
        Item('plot', editor=ComponentEditor(),
             show_label=False),
        Item(name = 'u0'),
        Item(name = 'v0'),
        Item(name = 'mu'),
        Item(name = 'plot_type'),
        width = 600,
        height = 800,
        resizable = True,
        title = "van der Pol oscillator",
        )

    def _plot_default(self):
        plot = Plot(self.plotdata)
        self.plotdata.set_data('t', self.t)
        self.plotdata.set_data('u', self.u)
        self.plotdata.set_data('v', self.v)
        plot.plot(("t", "u"), color = 'red', type_trait="plot_type", name = 'u')
        plot.plot(("y", "v"), color = 'blue', type_trait="plot_type", name = 'v')
        plot.legend.visible = True
        plot.title = "van der Pol oscillator"
        plot.x_axis.title = 't'
        return plot

    def _u_changed(self):
        self.plotdata.set_data('u', self._get_u())

    def _v_changed(self):
        self.plotdata.set_data('v', self._get_v())

    def _t_default(self):
        return self.t_default

    def _get_u(self):
        self.run_integration()
        return self.interpolated_yres['u']

    def _get_v(self):
        self.run_integration()
        return self.interpolated_yres['v']

    def __init__(self, ivp, indep_var_lim, N):
        super(ODESolViewer, self).__init__()
        self.t_default = np.linspace(indep_var_lim[0], indep_var_lim[1], 100)
        self.ivp = ivp
        self.N = N
        self.run_integration()

    def run_integration(self):
        self.ivp.integrate(self.t_default[-1], N = self.N)
        self.interpolated_yres = self.ivp.get_interpolated(self.t)

def get_gui(ODESys, indep_var_lim, depv0_by_token, param_vals, N=0):
    odesys = ODESys()
    params_by_symb = odesys.val_by_symb_from_by_token(param_vals)
    y0_by_symb = odesys.val_by_symb_from_by_token(depv0_by_token)

    t0, tend = indep_var_lim
    ivp = IVP(odesys, y0_by_symb, params_by_symb, t0)
    return ODESolViewer(ivp, indep_var_lim, N)


if __name__ == '__main__':
    viewer = get_gui(VanDerPolOscillator, (0.0, 10.0), {'u':1.0, 'v':1.0},
                     {'mu': 2.5}, N=0)
    viewer.configure_traits()
