#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Demonstrates codeexport, compilation and
integration using either:
  Gnu Scientific Library (GSL)
  ODEPACK
  Sundials
"""

from __future__ import division

import os
import logging

import argh

from coupled_decay import coupled_decay_numeric_vs_analytic

def main(implementation='odepack'):
    impl = implementation.lower()
    if impl == 'odepack':
        from symodesys.odepack import LSODES_IVP_Integrator as IVP_Integrator
    elif impl == 'gsl':
        from symodesys.gsl import GSL_IVP_Integrator as IVP_Integrator
    elif impl == 'sundials':
        from symodesys.sundials import CVODE_IVP_Integrator as IVP_Integrator
    else:
        raise ValueError("Unkown implementation: {}".format(implementation))

    import logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    integrator_kwargs = {
        'tempdir': os.path.join(
            os.path.dirname(__file__),
            'build_codegen_'+impl
        ),
        'save_temp': True,
        'logger': logger,
    }
    coupled_decay_numeric_vs_analytic(
        IVP_Integrator,
        integrator_kwargs,
        logger=logger
    )
    fmtstr = 'To look at the code (and Py independent code) look in {}'
    print(fmtstr.format(integrator_kwargs['tempdir']))

argh.dispatch_command(main)
