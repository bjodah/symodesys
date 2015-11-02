#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

from chemsys import ChemSys, y0, params
from symodesys.gui import get_chaco_viewer

if __name__ == '__main__':
    ODESys = ChemSys
    t0 = 0.0
    tend = 3.0
    N = 50
    viewer = get_chaco_viewer(ODESys(), y0, params, t0, tend, N)
    viewer.configure_traits()
    viewer.clean_up()
