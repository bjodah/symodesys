#!/usr/bin/env python
# -*- coding: utf-8 -*-

from symodesys.codeexport import F90_Code

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

    _source_files = ['transform_template.f90',
                     'transform_wrapper.c']

    _obj_files = ['transform_template.o',
                  'transform_wrapper.o',]

    _so_file = 'transform_wrapper.so'

    def __init__(self, exprs, inp, **kwargs):
        self._exprs = exprs
        self._inp = inp
        self._obj_files = [x.replace('_template', '')[:-4]+'.o' for x
                           in self.templates]
        super(Transformer, self).__init__(**kwargs)
        self._binary = self.compile_and_import_binary()


    def _write_code(self):
        # first we need to export the template in this module
        super(Transformer, self)._write_code()


    def _compile_obj(self):
        pyxpath = filter(partial(str.endswith, suffix='.pyx'),
                         self._written_files)
        pyx2obj(pyxpath) # .pyx -> (.c) -> .o
        super(Transformer, self)._compile_obj()

    def __call__(self, *args):
        return self._binary.transform(*args)
