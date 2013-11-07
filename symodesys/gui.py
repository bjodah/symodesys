#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from symodesys.ivp import IVP
from symodesys.gsl import GSL_IVP_Integrator

COLR = ['red', 'green', 'blue', 'black']

def get_chaco_viewer(odesys, y0, params, t0, tend, N):
    """
    Dynamically builds a Traits enabled class which
    can produce interactive plotting through Chaco

    Usage:
        >>>viewer=get_chaco_viewer(MyOscillator, y0, params, t0, tend, N)
        >>>viewer.configure_traits()
    """
    from enthought.chaco.api import ArrayPlotData, Plot
    from enthought.chaco.tools.api import PanTool, ZoomTool
    from enthought.enable.component_editor import ComponentEditor
    from enthought.traits.api import HasTraits, Instance, Array, Property, Range, Float, Enum
    from enthought.traits.ui.api import Item, View

    class ODESolViewer(object):
        """
        ODESolViewer makes a chaco viewer instance.
        """
        # It uses some not so pretty meta programming
        # in the __new__ part in order to enable setting
        # attributes before HasTraits meta-class machinery
        # gets its hands on it. It works, but could be rewritten
        # to enhance brevity and extensibility.

        plot = Instance(Plot)
        plotdata = Instance(ArrayPlotData, args=())

        plot_type = Enum("line", "scatter")


        def _plot_default(self):
            plot = Plot(self.plotdata)
            self.plotdata.set_data(str(odesys.indepv), getattr(self, str(odesys.indepv)))
            for dv in odesys.all_depv:
                self.plotdata.set_data(dv.func.__name__, getattr(self, dv.func.__name__))
            for i, dv in enumerate(odesys.all_depv):
                plot.plot((str(odesys.indepv), dv.func.__name__),
                          color = COLR[i % len(COLR)], type_trait="plot_type", name = dv.func.__name__)
            plot.legend.visible = True
            plot.title = 'Solution'
            if hasattr(odesys, 'title'):
                if odesys.title: plot.title = odesys.title
            plot.x_axis.title = str(odesys.indepv)

            # Add pan and zoom to the plot
            plot.tools.append(PanTool(plot, constrain_key="shift"))
            zoom = ZoomTool(plot)
            plot.overlays.append(zoom)
            return plot

        def __new__(cls, *args, **kwargs):
            all_init = ['init_' + depv.func.__name__ for depv in odesys.all_depv]
            all_param = map(str,odesys.param_and_sol_symbs)

            items = [Item('plot', editor=ComponentEditor(),
                          show_label=False),
                     Item(name='plot_type')]

            setattr(cls, str(odesys.indepv), Array)
            setattr(cls, '_'+str(odesys.indepv)+'_default',
                    lambda self: getattr(self, str(odesys.indepv)+'_default'))

            for depv in odesys.all_depv:
                # TODO: enusre 'init_' + ... does not collide
                depv_str = depv.func.__name__
                setattr(cls, depv_str, Property(Array, depends_on=all_init+all_param))
                setattr(cls, 'init_'+depv_str, Range(low=0.0, high=10.0, value=1.0))

                def getter(self, depstr=depv_str):
                    # The keyword argument needs to be there not to be overwritten
                    # in the closure as the loop proceeds (python behaviour)
                    return self.interpolated_yres[self.ivp.get_index_of_depv(depstr),:]
                setattr(cls, '_get_'+depv_str, getter)

                def changed(self, depstr=depv_str):
                    self.plotdata.set_data(depstr, getattr(self, '_get_'+depstr)())
                setattr(cls, '_'+depv_str+'_changed', changed)

                def init_changed(self, depstr=depv_str):
                    self.run_integration()
                setattr(cls, '_init_'+depv_str+'_changed', init_changed)

                items.append(Item(name=str('init_'+depv_str)))

            for param in odesys.param_and_sol_symbs:
                setattr(cls, str(param), Range(low=0.0, high=10.0, value=1.0))
                items.append(Item(name=str(param)))
                def param_changed(self, depstr=str(param)):
                    self.run_integration()
                setattr(cls, '_'+str(param)+'_changed', param_changed)


            view = View(*items,
                        width = 600,
                        height = 800,
                        resizable = True)
                       #title = getattr(odesys, 'title', "Solution"))
            setattr(cls, 'traits_view', view)
            newcls = type(cls.__name__, (HasTraits,),
                          {k: v for k,v in cls.__dict__.items() if k != '__new__'})
            instance = newcls.__new__(newcls)
            instance.__dict__ = {}
            instance.__init__(*args, **kwargs)
            return instance

        def __init__(self, odesys_, y0, params, t0, tend, N, Integrator=GSL_IVP_Integrator):
            for k, v in y0.items():
                if hasattr(k, 'func'):
                    s = k.func.__name__
                else:
                    s = str(k)
                if hasattr(self, 'init_'+s):
                    setattr(self, 'init_'+s, v)
                else:
                    raise AttributeError('No init cond. {}'.format('init_'+s))
            for k, v in params.items():
                if hasattr(k, 'func'):
                    s = k.func.__name__
                else:
                    s = str(k)
                if hasattr(self, k):
                    setattr(self, k, v)
                else:
                    raise AttributeError('No param {}'.format(k))
            setattr(self, str(odesys.indepv) + '_default', np.linspace(t0, tend, 500))
            self.ivp = IVP(odesys_, y0, params, t0, integrator=GSL_IVP_Integrator())
            self.N = N
            self.old_args = []
            #super(ODESolViewer, self).__init__()
            self.run_integration()

        @property
        def depv_init(self):
            return {dv: getattr(self, 'init_'+str(dv.func.__name__)) for dv in odesys.all_depv}

        @property
        def params(self):
            return {p: getattr(self, str(p)) for p in odesys.param_and_sol_symbs}

        def run_integration(self):
            """
            Runs the actual numeric integration
            """

            # Only run if updated since last run
            args = []
            for depv in odesys.all_depv:
                args.append(getattr(self, 'init_'+depv.func.__name__))
            for param in odesys.param_and_sol_symbs:
                args.append(getattr(self, str(param)))
            if args == self.old_args: return
            self.old_args = args

            # Make a call to the heavy lifting.
            self.interpolated_yres = self._integrate(
                self.depv_init, self.params,
                getattr(self, str(odesys.indepv)+'_default')[-1], self.N)
            for depv in odesys.all_depv:
                depv_str = depv.func.__name__
                self.plotdata.set_data(depv_str, getattr(self, '_get_'+depv_str)())

        def _integrate(self, depv_init, params, tend, N):
            self.ivp.depv_init=depv_init
            self.ivp.params=params
            self.ivp.integrate(tend, N = N)
            return self.ivp.get_interpolated(getattr(self, str(odesys.indepv)))

        def clean_up(self):
            try:
                self.ivp.clean()
            except AttributeError:
                pass

    return ODESolViewer(odesys, y0, params, t0, tend, N)
