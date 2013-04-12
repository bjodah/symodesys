#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

import numpy as np

from enthought.chaco.api import ArrayPlotData, Plot
from enthought.chaco.tools.api import PanTool, ZoomTool
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
        plot.plot(("t", "v"), color = 'blue', type_trait="plot_type", name = 'v')
        plot.legend.visible = True
        plot.title = "van der Pol oscillator"
        plot.x_axis.title = 't'

        # Add pan and zoom to the plot
        plot.tools.append(PanTool(plot, constrain_key="shift"))
        zoom = ZoomTool(plot)
        plot.overlays.append(zoom)
        return plot


    def _u_changed(self):
        self.plotdata.set_data('u', self._get_u())

    def _v_changed(self):
        self.plotdata.set_data('v', self._get_v())

    def _t_default(self):
        return self.t_default

    def _get_u(self):
        self.run_integration()
        return self.interpolated_yres[self.ivp.get_index_of_depv('u'),:]

    def _get_v(self):
        self.run_integration()
        return self.interpolated_yres[self.ivp.get_index_of_depv('v'),:]

    def __init__(self, ODESys, y0, params, t0, tend, N):
        super(ODESolViewer, self).__init__()
        for k, v in y0.items():
            if hasattr(self, k+'0'):
                setattr(self, k+'0', v)
            else:
                raise AttributeError('No init cond. {}'.format(k+'0'))
        for k, v in params.items():
            pass
        self.t_default = np.linspace(t0, tend, 500)
        self.ivp = IVP(ODESys(), y0, params, t0)
        self.N = N
        self.run_integration()

    @property
    def init_vals(self):
        return {'u': self.u0, 'v': self.v0}

    @property
    def param_vals(self):
        return {'mu': self.mu}

    def run_integration(self):
        print("Runnig integration..", end='')
        self.ivp.init_vals=self.init_vals
        self.ivp.param_vals=self.param_vals
        self.ivp.integrate(self.t_default[-1], N = self.N)
        self.interpolated_yres = self.ivp.get_interpolated(self.t)
        print("...DONE! {}".format(','.join([str(self.ivp.Yres()[-1,0,0]), str(self.ivp.tout[-1])])))
        print(self.ivp.Yres.cache)
        print([(id(k), id(v)) for k,v in self.ivp.Yres.cache.items()])


if __name__ == '__main__':
    ODESys=VanDerPolOscillator
    y0={'u':1.0, 'v':1.0}
    params={'mu': 2.5},
    t0=0.0
    tend=10.0
    N=0
    viewer = ODESolViewer(ODESys, y0, params, t0, tend, N)
    viewer.configure_traits()
