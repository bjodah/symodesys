#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

from van_der_pol import VanDerPolOscillator
from symodesys.gui import get_chaco_viewer

if __name__ == '__main__':
    ODESys = VanDerPolOscillator
    y0 = {'u':1.0, 'v':1.0}
    params = {'mu': 2.5}
    t0 = 0.0
    tend = 10.0
    N = 50
    viewer = get_chaco_viewer(ODESys(), y0, params, t0, tend, N)
    viewer.configure_traits()
    viewer.clean_up()
