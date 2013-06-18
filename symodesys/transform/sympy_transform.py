#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Lightning fast binary callbacks using F90 for array expressions
"""

import os
from functools import partial
from collections import defaultdict

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

    _templates = ['transform_template.f90',
                 'transform_wrapper_template.pyx']

    _source_files = ['transform.f90',
                     'transform_wrapper.c']

    _obj_files = ['transform.o',
                  'transform_wrapper.o',]

    _so_file = 'transform_wrapper.so'

    def __init__(self, exprs, inp, **kwargs):
        self._cached_files = self._cached_files or []
        self._cached_files += [self._so_file]
        self._basedir = os.path.dirname(__file__)
        self._exprs = exprs
        self._inp = inp

        self._robust_exprs, self._robust_inp_symbs, self._robust_inp_dummies = \
            self.robustify(self._exprs, self._inp)

        self._obj_files = [x.replace('_template', '')[:-4]+'.o' for x
                           in self._templates]
        super(Transformer, self).__init__(**kwargs)
        self._binary = self.compile_and_import_binary()


    def robustify(self, exprs, inp):
        dummies = []
        def make_dummy(key, not_in_any):
            if isinstance(key, sympy.Function):
                candidate = sympy.Symbol(key.func.__name__, real=key.is_real)
            elif isinstance(key, sympy.Derivative):
                candidate = sympy.Symbol('d'+key.args[0].func.__name__+\
                                         ''.join('d'+x.name for x in key.args[1:]))
            elif isinstance(key, sympy.Symbol):
                candidate = key
            elif isinstance(key, str):
                candidate = sympy.Symbol(key)
            else:
                raise NotImplementedError
            # Now lets make sure the dummie symbol is unique
            if candidate not in dummies and not any(
                    [not_in.has(candidate) for not_in in not_in_any]):
                dummies.append(candidate)
                return candidate
            else:
                return make_dummy(candidate.name+'p', not_in)
        # First let's determine the used keys
        used_keys = []
        used_keys_dummies = []
        derivs = defaultdict(dict)
        any_expr_has = lambda x: any([expr.has(x) for expr in exprs])
        for key in filter(any_expr_has, inp):
            used_keys.append(key)
            if isinstance(key, sympy.Symbol):
                used_keys_dummies.append(key)
            else:
                new_dummy = make_dummy(key, not_in_any=exprs)
                used_keys_dummies.append(new_dummy)
                # Determine order of derivative
                if isinstance(key, sympy.Derivative):
                    derivorder = len(key.args)-1
                    derivs[derivorder][key] = new_dummy
                elif isinstance(key, sympy.Function):
                    derivorder = 0
                    derivs[derivorder][key] = new_dummy
                elif isinstance(key, sympy.Symbol):
                    pass
                else:
                    raise NotImplementedError
        for n in sorted(derivs.keys())[::-1]:
            # Substitutions of symbolified derivatives
            # need to be made for the highest degrees of
            # derivatives first
            for deriv, symb in derivs[n].items():
                new_exprs = []
                for expr in exprs:
                    new_exprs.append(expr.subs({deriv: symb}))
                exprs = new_exprs

        # Now we are almost there..
        # But e.g. z(t) has included t into used_keys
        # we need to double check that t is truly used
        # now that z(t) has been turned into z
        truly_used_keys, truly_used_keys_dummies = [], []
        for k, ks in zip(used_keys, used_keys_dummies):
            if any_expr_has(ks):
                truly_used_keys.append(k)
                truly_used_keys_dummies.append(ks)

        return exprs, truly_used_keys, truly_used_keys_dummies

        # return ufunc, tuple(sorted(
        #     truly_used_keys,
        #     key=lambda x: truly_used_keys_dummies[truly_used_keys.index(x)]))

    def __call__(self, *args):
        """
        Give arrays in as arguments in same order as `inp` at init.
        """
        dummies = dict(zip(self._robust_inp_symbs, self._robust_inp_dummies))
        idxs = [self._inp.index(symb) for symb in self._robust_inp_symbs]
        return self._binary.transform(*[args[i] for i in idxs])


    def variables(self):
        args = [str(v) for v in self._robust_inp_dummies]
        cses, exprs_in_cse = sympy.cse(
            self._robust_exprs, symbols=sympy.numbered_symbols('cse'))
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
