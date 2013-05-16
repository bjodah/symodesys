#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from symodesys.ivp import IVP
from symodesys.gsl import GSL_IVP_Integrator

def get_chaco_viewer(odesys, y0, params, t0, tend, N):
    """
    Usage:
        >>>viewer=get_chaco_viewer(MyOscillator, y0, params, t0, tend, N)
        >>>viewer.configure_traits()
    """
    from enthought.chaco.api import ArrayPlotData, Plot
    from enthought.chaco.tools.api import PanTool, ZoomTool
    from enthought.enable.component_editor import ComponentEditor
    from enthought.traits.api import HasTraits, Instance, Array, Property, Range, Float, Enum
    from enthought.traits.ui.api import Item, View

    class ODESolViewer(HasTraits):

        plot = Instance(Plot)
        plotdata = Instance(ArrayPlotData, args=())

        plot_type = Enum("line", "scatter")

        def _plot_default(self):
            plot = Plot(self.plotdata)
            self.plotdata.set_data(str(odesys.indepv), getattr(self, str(odesys.indepv)))
            for dv in odesys.all_depv:
                self.plotdata.set_data(dv.func.__name__, getattr(self, dv))
            for dv in odesys.all_depv:
                plot.plot((str(odesys.indepv), dv.func.__name__),
                          color = 'red', type_trait="plot_type", name = dv.func.__name__)
            plot.legend.visible = True
            plot.title = getattr(odesys, 'title', 'Solution')
            plot.x_axis.title = str(odesys.indepv)

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

        def __init__(self, odesys_, y0, params, t0, tend, N, Integrator=GSL_IVP_Integrator):
            for k, v in y0.items():
                if hasattr(self, 'init_'+k):
                    setattr(self, 'init_'+k, v)
                else:
                    raise AttributeError('No init cond. {}'.format('init_'+k))
            for k, v in params.items():
                if hasattr(self, k):
                    setattr(self, k, v)
                else:
                    raise AttributeError('No param {}'.format(k))
            self.t_default = np.linspace(t0, tend, 500)
            self.ivp = IVP(odesys_, y0, params, t0, integrator=GSL_IVP_Integrator())
            self.N = N
            super(ODESolViewer, self).__init__()
            self.run_integration()

        @property
        def init_vals(self):
            return {dv: getattr(self, 'init_'+str(dv.func.__name__)) for dv in odesys.all_depv}

        @property
        def param_vals(self):
            return {p: getattr(self, str(p)) for p in odesys.param_and_sol_symbs}

        def run_integration(self):
            # This may need caching
            self.interpolated_yres = self._integrate(self.init_vals, self.param_vals,
                                                     self.t_default[-1], self.N)

        def _integrate(self, init_vals, param_vals, tend, N):
            self.ivp.init_vals=init_vals
            self.ivp.param_vals=param_vals
            self.ivp.integrate(tend, N = N)
            print self.t, type(self.t)
            return self.ivp.get_interpolated(self.t)

    viewer = ODESolViewer
    all_init = ['init_' + depv.func.__name__ for depv in odesys.all_depv]
    all_param = odesys.param_and_sol_symbs

    items = [Item('plot', editor=ComponentEditor(),
                  show_label=False),
             Item(name='plot_type')]

    setattr(viewer, str(odesys.indepv), Array) # e.g. `t`
    setattr(viewer, '_'+str(odesys.indepv)+'_default', lambda self: self.t_default)

    for depv in odesys.all_depv:
        # TODO: enusre 'init_' + ... does not collide
        setattr(viewer, depv.func.__name__, Property(Array, depends_on=all_init+all_param))
        setattr(viewer, 'init_'+depv.func.__name__, Range(low=0.0, high=10.0, value=1.0))
        print 'setattr, init_{}'.format(depv.func.__name__) ###
        setattr(viewer, '_'+depv.func.__name__+'_changed',
                lambda self: self.plotdata.set_data(depv.func.__name__, getattr(self, '_get_'+depv.func.__name__)))
        def getter(self):
            self.run_integration()
            return self.interpolated_yres[self.ivp.get_index_of_depv(depv.func.__name__),:]
        setattr(viewer, '_get_'+depv.func.__name__, getter)
        items.append(Item(name=str('init_'+depv.func.__name__)))

    for param in odesys.param_and_sol_symbs:
        setattr(viewer, str(param), Range(low=0.0, high=10.0, value=1.0))
        items.append(Item(name=str(param)))

    view = View(*items,
            width = 600,
            height = 800,
            resizable = True,
            title = getattr(odesys,'title',"Solution"))
    setattr(viewer, 'traits_view', view)

    return viewer(odesys, y0, params, t0, tend, N)
