#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

import numpy as np

from chaco.api import ArrayPlotData, Plot
from chaco.tools.api import PanTool, ZoomTool
from enable.component_editor import ComponentEditor
from traits.api import HasTraits, Instance, Array, Property, Range, Float, Enum
from traitsui.api import Item, View

from van_der_pol import VanDerPolOscillator
from symodesys.ivp import IVP
from symodesys.gsl import GSL_IVP_Integrator as Binary_IVP_Integrator
#from symodesys.odepack import LSODES_IVP_Integrator as Binary_IVP_Integrator
from symodesys.helpers import cache

class ODESolViewer(HasTraits):

    plot = Instance(Plot)
    plotdata = Instance(ArrayPlotData, args=())

    indepv = Array
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
        self.plotdata.set_data('t', self.indepv)
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

    def _indepv_default(self):
        return self.indepv_default

    def _get_u(self):
        self.run_integration()
        return self.interpolated_yres[self.ivp.get_index_of_depv('u'),:]

    def _get_v(self):
        self.run_integration()
        return self.interpolated_yres[self.ivp.get_index_of_depv('v'),:]

    def __init__(self, ODESys, y0, params, t0, tend, N,
                 Integrator=Binary_IVP_Integrator):
        for k, v in y0.items():
            if hasattr(self, k+'0'):
                setattr(self, k+'0', v)
            else:
                raise AttributeError('No init cond. {}'.format(k+'0'))
        for k, v in params.items():
            if hasattr(self, k):
                setattr(self, k, v)
            else:
                raise AttributeError('No param {}'.format(k))
        self.indepv_default = np.linspace(t0, tend, 2048)
        self.ivp = IVP(ODESys(), y0, params, t0,
                       integrator=Integrator(nderiv=2, tempdir='tmp', save_temp=True))
        self.N = N
        super(ODESolViewer, self).__init__()
        self.run_integration()

    @property
    def depv_init(self):
        return {'u': self.u0, 'v': self.v0}

    @property
    def params(self):
        return {'mu': self.mu}

    def run_integration(self):
        self.interpolated_yres = self._integrate(
            self.depv_init, self.params, self.indepv_default[-1], self.N,
        )

    def _integrate(self, depv_init, params, tend, N):
        self.ivp.depv_init=depv_init
        self.ivp.params=params
        self.ivp.integrate(tend, N=N, step_type='bsimp')
        return self.ivp.get_interpolated(self.indepv)

if __name__ == '__main__':
    ODESys=VanDerPolOscillator
    y0={'u':1.0, 'v':1.0}
    params={'mu': 2.5}
    t0=0.0
    tend=10.0
    N=100
    viewer = ODESolViewer(ODESys, y0, params, t0, tend, N)
    viewer.configure_traits()
    viewer.ivp.clean()
