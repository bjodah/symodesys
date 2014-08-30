from __future__ import print_function, division

# For performance reasons it is preferable that the
# numeric integration is performed in a compiled language
# the codeexport module provide classes enabling FirstOrderODESystem
# instances to be used as blueprint for generating, compiling
# and importing a binary which performs the computations.

# Both C, C++ and Fortran is considered, but since
# the codegeneration uses templates, one can easily extend
# the functionality to other languages.

# stdlib imports
import tempfile
import shutil
import re
import os

from collections import OrderedDict
from functools import reduce, partial
from operator import add


# External imports
import sympy
import mako
import mako.template


# Intrapackage imports
from symodesys.integrator import IVP_Integrator
from pycodeexport.codeexport import Generic_Code, DummyGroup, ArrayifyGroup


class ODESys_Code(Generic_Code):
    """
    Wraps some sympy functionality of code generation from matrices
    returned by FirstOrderODESystem.dydt and FirstOrderODESystem.dydt_jac
    """

    tempdir_basename = "_symodesys_compile"
    compile_kwargs = {
        'options': ['warn','pic','fast']
    }

    def __init__(self, fo_odesys, **kwargs):
        self._fo_odesys = fo_odesys
        super(ODESys_Code, self).__init__(**kwargs)


    @property
    def NY(self):
        return len(self._fo_odesys.na_depv)


    @property
    def prog_param_symbs(self):
        """
        Bridge slight difference in meaning of parameter going from
        ODESystem to numerical integration routine (analytically
        solving part of ODESys introduces new 'parameters')
        """
        return self._fo_odesys.param_and_sol_symbs


    def variables(self):
        dummy_groups = (
            DummyGroup('depvdummies', self._fo_odesys.na_depv),
            DummyGroup('paramdummies', self._fo_odesys.param_and_sol_symbs),
        )

        arrayify_groups = (
            ArrayifyGroup('depvdummies', self.depv_tok, self.depv_offset),
            ArrayifyGroup('paramdummies', self.param_tok, self.param_offset)
        )

        code_func_cse, code_func_exprs = self.get_cse_code(
            self._fo_odesys.na_f.values(), 'csefunc', dummy_groups, arrayify_groups)

        y0 = ', '.join(['1.0'] * self.NY)
        y0_names = ', '.join(map(str, self._fo_odesys.na_depv))
        params = ', '.join(['1.0'] * len(self.prog_param_symbs))
        param_names = ', '.join(map(str, self.prog_param_symbs))

        indepv = self._fo_odesys.indepv

        # TODO: this is inefficient, implement linear scaling algo in
        # FirstOrderODESystem
        sparse_jac = OrderedDict(reduce(add, [
            [((i, j), expr) for j, expr in enumerate(row) if expr != 0]\
            for i, row in enumerate(self._fo_odesys.na_jac.tolist())
            ]))

        dfdt = self._fo_odesys.na_dfdt.values()

        code_jac_cse, code_jac_exprs = self.get_cse_code(
            sparse_jac.values() + dfdt, 'csejac', dummy_groups, arrayify_groups)

        code_pure_dfdt_cse, code_pure_dfdt_exprs = self.get_cse_code(
            dfdt, 'csedfdt', dummy_groups, arrayify_groups)

        code_dfdt_exprs = code_jac_exprs[len(sparse_jac):]
        code_jac_exprs = zip(sparse_jac.keys(),
                             code_jac_exprs[:len(sparse_jac)])


        # Populate ia, ja (sparse index specifiers using fortran indexing)
        # see documentation of LSODES in ODEPACK for definition
        # (Yale sparse matrix)
        # ja contains row indices of nonzero elements
        # ia contains index in ja where row i starts
        ia, ja = [1], []
        k = 1 # <--- index starts at 1 in fortran, not at 0 as in C/Python
        code_yale_jac_exprs = []
        code_yale_jac_cse = []
        for ci in range(self.NY):
            col_exprs = []
            cur_ja = []
            for ri in [sri for sri, sci in sparse_jac.keys() if sci == ci]:
                cur_ja.append(ri+1) # Fortran index
                k += 1
                col_exprs.append(sparse_jac[(ri,ci)])

            # Store indices in ja
            ja.extend(cur_ja)

            # Store indicies in ia
            ia.append(k)

            # Extract common subexpressions for this column
            cse_defs, cse_exprs = sympy.cse(
                col_exprs, symbols = sympy.numbered_symbols(
                    'csejaccol{}'.format(ci)))

            # Format code: expressions in cse terms
            code_exprs = zip(cur_ja, [
                self.as_arrayified_code(x, dummy_groups, arrayify_groups) for x \
                in cse_exprs
            ])
            code_yale_jac_exprs.append(code_exprs)

            # Format code: CSE definitions
            code_cse_defs=[]
            for var_name, var_expr in cse_defs:
                code_var_expr = self.as_arrayified_code(
                    var_expr, dummy_groups, arrayify_groups)
                code_cse_defs.append((var_name, code_var_expr))
            code_yale_jac_cse.append(code_cse_defs)
        ia = ia [:-1]

        return {'NY': self.NY,
                'NNZ': len(sparse_jac),
                'IA': ia,
                'JA': ja,
                'NPARAM': len(self.prog_param_symbs),
                'y0_comma_sep_str': y0,
                'y0_names': y0_names,
                'param_vals_comma_sep_str': params,
                'param_names': param_names,
                'cse_func': code_func_cse, 'f': code_func_exprs,
                'cse_jac': code_jac_cse, 'jac': code_jac_exprs,
                'yale_jac_exprs': code_yale_jac_exprs,
                'yale_jac_cse': code_yale_jac_cse,
                'dfdt': code_dfdt_exprs,
                'pure_dfdt_exprs': code_pure_dfdt_exprs,
                'pure_dfdt_cse': code_pure_dfdt_cse,
        }


class Binary_IVP_Integrator(IVP_Integrator):
    """
    Generic baseclass for writing integrators
    which generate code (set
    pycompilation.codeexport.Generic_Code derived
    class as CodeClass in derived class).
    """

    CodeClass = None
    _code = None # will hold instance (set_fo_odesys)

    def __init__(self, **kwargs):
        self.tempdir = kwargs.pop('tempdir', None)
        self.save_temp = kwargs.pop('save_temp', False)
        self.logger = kwargs.pop('logger', None)
        super(Binary_IVP_Integrator, self).__init__(**kwargs)


    def set_fo_odesys(self, fo_odesys):
        if fo_odesys != self._fo_odesys:
            # Don't recompile if we system is the same
            super(Binary_IVP_Integrator, self).set_fo_odesys(fo_odesys)
            self._binary_mod = None # <-- Clears cache
            self._code = self.CodeClass(
                fo_odesys = self._fo_odesys,
                tempdir = self.tempdir,
                save_temp = self.save_temp,
                logger=self.logger,
            )


    def clean(self):
        if self._code:
            self._code.clean()


    @property
    def binary_mod(self):
        """
        Returns compiled and imported module.
        Note: lazy caching is employed, set self._binary_mod equal to None to invalidate
        """
        if self._binary_mod == None:
            self._binary_mod = self._code.compile_and_import_binary()
        return self._binary_mod
