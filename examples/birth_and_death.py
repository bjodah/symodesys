#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

# stdlib imports
import sys
import logging
import time
from collections import OrderedDict


# External package imports
import sympy
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
from matplotlib.ticker import MaxNLocator
import argh

# Project internal imports
from symodesys.odesys import SimpleFirstOrderODESystem
from symodesys.ivp import IVP
from symodesys.odepack import LSODES_IVP_Integrator as Integrator

"""
This example illustrates the need of symodesys
to generate loop structures in code generation.
(and efficiently handle sparse systems)

TODO: use a banded solver here.
"""


class BirthDeathSystem(SimpleFirstOrderODESystem):

    # Following two lines are optional but useful for
    # automatic labeling when plotting:

    @classmethod
    def mk_of_dim(cls, dim):
        """
        Makes an instance of provided `dim`
        """
        cls.depv_tokens = ['y{}'.format(i) for i in range(dim)]
        cls.param_tokens   = ['l{}'.format(i) for i in range(dim)]
        cls.dim = dim
        return cls()

    @property
    def expressions(self):
        y = [self['y{}'.format(i)] for i in range(self.dim)]
        l = [self['l{}'.format(i)] for i in range(self.dim)]
        first =  (y[0], -l[0]*y[0])
        last  = (y[self.dim-1], l[self.dim-2]*y[self.dim-2])
        return OrderedDict([first] + [
            (y[i], y[i-1]*l[i-1]-y[i]*l[i]) for \
            i in range(1,self.dim-1)
        ] + [last])

def plot_population(ivp, depv_z_map):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    verts=[]
    nt = 1000
    t = np.linspace(ivp.indepv_out()[0],
                    ivp.indepv_out()[-1], nt)
    n = len(depv_z_map)
    maxval = 0

    # Get interpolated (and back-transformed) results
    tm = time.time()
    for depv, x in depv_z_map.items():
	interpol = ivp.get_interpolated(t, [depv])
	# Let polygons be aligned with bottom
        # (important: choose nt large!)
	interpol[0,0] = 0.0
	interpol[0,-1] = 0.0
        # Check fo new max value:
	maxval = max(maxval, np.max(interpol))
	stacked = zip(t,interpol.reshape(nt))
        verts.append(stacked)#.transpose())
    print("Time of interpolation and back transformation: {}".format(
        time.time()-tm))

    # Generate plot
    tm = time.time()
    poly = PolyCollection(verts, facecolors = [
        'rgbk'[i%4] for i in range(len(verts))])
    poly.set_alpha(0.7)
    ax.add_collection3d(poly, zs=[
        x+1 for x in depv_z_map.values()], zdir='y')

    ax.set_xlabel('t')
    ax.set_xlim3d(0, t[-1])
    ax.set_ylabel('y')
    ax.set_ylim3d(0.5, n+0.5)
    ax.set_zlabel('z')
    ax.set_zlim3d(0, maxval)

    ya = ax.get_yaxis()
    ya.set_major_locator(MaxNLocator(integer=True))

    print ("Time of plot generation: {}".format(time.time()-tm))
    plt.show()


def get_transformed_ivp(n, logger2intgrtr=False, logger=None):
    bd = BirthDeathSystem.mk_of_dim(n)
    print("Original system:")
    print('\n'.join([str(x) for x in bd.eqs]))
    y0 = {'y{}'.format(i): (n-i+1)*1.0 for i in range(n)}
    params = {'l{}'.format(i): (n-i+1)/n for i in range(n)}
    t0, tend, N = 0, 15.0, 500
    tm = time.time()
    integrator_logger = logger if logger2intgrtr else None
    ivp1 = IVP(bd, y0, params, t0,
              integrator=Integrator(
                  save_temp=True, tempdir='build_birth_and_death/',
                  logger=integrator_logger),
              logger=logger)
    print("Time to construct inital IVP = {} s".format(time.time()-tm))

    # Apply a logarithmic transform
    tm = time.time()
    trnsfm, inv_trnsfm = OrderedDict(), OrderedDict()
    for i in range(n):
        # Generate the transform (incl. inverse)
        lny = ivp1.fo_odesys.mk_depv('lny{}'.format(i))
        y = bd['y{}'.format(i)]
        trnsfm[lny] = sympy.log(y)
        inv_trnsfm[y] = sympy.exp(lny)
    ivp2 = ivp1.use_internal_depv_trnsfm(trnsfm, inv_trnsfm)
    print("Time of analytic forward transform = {} s".format(time.time() - tm))

    # Solve constants
    tm = time.time()
    ivp2.recursive_analytic_reduction(complexity=0)
    t_analytic = time.time() - tm
    print("Transformed:")
    print('\n'.join([str(x) for x in ivp2.fo_odesys.eqs]))

    # Run the integration
    tm = time.time()
    tm_int = ivp2.integrate(tend, N)
    print("Time of integration: {} s".format(tm_int))
    print("  + overhead (codegen, compilation): {} s".format(
        time.time()-tm))
    return ivp2


def main(n=5):
    """
    Plot the evolutoin of a birth death network
    """
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger(__file__)

    with get_transformed_ivp(n, logger=logger) as ivp:
        # ivp.plot(interpolate=True, datapoints=False, show=True)
        sorted_keys = sorted(
            ivp.all_depv, key=lambda x: int(x.func.__name__[1:]))
        plot_population(ivp, depv_z_map=OrderedDict(
            [(k, i) for i, k in enumerate(sorted_keys)]))


argh.dispatch_command(main)
