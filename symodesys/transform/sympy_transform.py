#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Lightning fast binary callbacks using F90 for array expressions
"""

import os
from functools import partial

import sympy

from symodesys.codeexport import F90_Code
from symodesys.helpers.compilation import pyx2obj


class Transformer(F90_Code):
    """
    The the Transformer instance can export expressions
    (given in `exprs`) to C/Fortran and compile callbacks
    which takes input arrays/scalars (given in `inp`).
    The callback returns a two dimesional array with
    shape: (len(exprs), len(inp))
    """

    templates = ['transform_template.f90',
                 'transform_wrapper_template.pyx']

    _source_files = ['transform.f90',
                     'transform_wrapper.c']

    _obj_files = ['transform.o',
                  'transform_wrapper.o',]

    _so_file = 'transform_wrapper.so'

    def __init__(self, exprs, inp, **kwargs):
        self._basedir = os.path.dirname(__file__)
        self._exprs = exprs
        self._inp = inp
        self._obj_files = [x.replace('_template', '')[:-4]+'.o' for x
                           in self.templates]
        super(Transformer, self).__init__(**kwargs)
        self._binary = self.compile_and_import_binary()


    def variables(self):
        args = [str(v) for v in self._inp]
        cses, exprs_in_cse = sympy.cse(
            self._exprs, symbols=sympy.numbered_symbols('cse'))
        return {'ARGS_COMMA': ', '.join(args),
                'ARGS': args,
                'CSES': cses,
                'EXPRS_IN_CSE': [self.wcode(x) for x in exprs_in_cse],
                'N_EXPRS': len(self._exprs)}

    def _write_code(self):
        # first we need to export the template in this module
        super(Transformer, self)._write_code()


    def _compile_obj(self):
        pyxpath = [x for x in self._written_files if x.endswith('.pyx')][0]
        pyx2obj(pyxpath) # .pyx -> (.c) -> .o
        super(Transformer, self)._compile_obj(self._source_files[:-1])

    def __call__(self, *args):
        return self._binary.transform(*args)
