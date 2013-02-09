#!/usr/bin/env python
# -*- coding: utf-8 -*-

from traits.api import HasTraits, Array, Range, Float, Enum, on_trait_change, Property
from traitsui.api import View, Item
from chaco.chaco_plot_editor import ChacoPlotItem
from numpy import arange

from van_der_pol import VanDerPolOscillator
from symodesys.integrator import SciPy_IVP_Integrator

class ODESolViewer(HasTraits):
    x = Array
    u = Property(Array, depends_on = ['u0', 'v0', 'mu'])
    v = Property(Array, depends_on = ['u0', 'v0', 'mu'])
    u0 = Range(low = -10.0, high = 10.0, value = 1.0)
    v0 = Range(low = -10.0, high = 10.0, value = 0.0)
    mu = Range(low = -10.0, high = 10.0, value = 3.0)
    plot_type = Enum("line", "scatter")

    traits_view = View(ChacoPlotItem("x", "u",
                                     type_trait = "plot_type",
                                     resizable = True,
                                     x_label = "X",
                                     y_label = "Y",
                                     x_bounds = ( 0, 10),
                                     x_auto = False,
                                     y_bounds = ( - 10, 10),
                                     y_auto = False,
                                     color = "blue",
                                     bgcolor = "white",
                                     border_visible = True,
                                     border_width = 1,
                                     title = 'Y vs. X',
                                     padding_bg_color = "lightgray"
                       ),
                       Item(name = 'u0'),
                       Item(name = 'v0'),
                       Item(name = 'mu'),
                       Item(name = 'plot_type'),
                       resizable = True,
                       buttons = ["OK"],
                       title = "van der Pol oscialltor",
                       width = 600, height = 800)


    def _x_default(self):
        return arange( 0, 10, 0.1)

    def _get_u(self):
        return integrate(self.u0, self.v0, self.mu)[:, 0]

    def _get_v(self):
        return integrate(self.u0, self.v0, self.mu)[:, 1]

def integrate(u0, v0, mu):
    vdpo = VanDerPolOscillator()
    intr = SciPy_IVP_Integrator(vdpo, [mu])

    int_kwargs = {'abstol': 1e-6,
                  'reltol': 1e-6}

    #N = 0 # adaptive stepsize controls output
    N = 101
    intr.integrate([u0, v0], 0.0, 10.0, N, **int_kwargs)
    return intr.yout

if __name__ == '__main__':
    viewer = ODESolViewer()
    viewer.configure_traits()
