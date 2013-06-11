#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import shutil
import tempfile

# other imports
import sympy
import numpy as np
from sympy import C
from sympy.utilities.lambdify import implemented_function

from sympy.utilities.codegen import (
    get_code_generator, CodeGenArgumentListError, Routine
)

from sympy.utilities.autowrap import (
    CodeWrapper, F2PyCodeWrapper, CythonCodeWrapper,
    DummyWrapper,
)

class MultiCodeWrapper(CodeWrapper):
    """
    Extends CodeWrapper to support multiple output values
    """

    def _generate_code(self, main_routines, routines):
        routines.extend(main_routines)
        self.generator.write(
            routines, self.filename, True, self.include_header,
            self.include_empty)

    def wrap_code(self, routines, helpers=[]):
        workdir = self.filepath or tempfile.mkdtemp("_sympy_compile")
        if not os.access(workdir, os.F_OK):
            os.mkdir(workdir)
        oldwork = os.getcwd()
        os.chdir(workdir)
        try:
            sys.path.append(workdir)
            self._generate_code(routines, helpers)
            self._prepare_files(routines)
            self._process_files(routines) # <--- _process_files doesnt seem to use arg
            mod = __import__(self.module_name)
        finally:
            sys.path.remove(workdir)
            CodeWrapper._module_counter += 1
            os.chdir(oldwork)
            if not self.filepath:
                shutil.rmtree(workdir)

        return self._get_wrapped_function(mod, routunes)



class MultiF2PyCodeWrapper(MultiCodeWrapper, F2PyCodeWrapper):
    pass


class MultiCythonCodeWrapper(MultiCodeWrapper, CythonCodeWrapper):

    def _prepare_files(self, routines):
        pyxfilename = self.module_name + '.pyx'
        codefilename = "%s.%s" % (self.filename, self.generator.code_extension)

        # pyx
        with open(pyxfilename, 'w') as f:
            self.dump_pyx(routines, f, self.filename,
                self.include_header, self.include_empty)

        # setup.py
        ext_args = [repr(self.module_name), repr([pyxfilename, codefilename])]
        with open('setup.py', 'w') as f:
            print >> f, CythonCodeWrapper.setup_template % {
                'args': ", ".join(ext_args)}


def _get_code_wrapper_class(backend):
    wrappers = { 'F2PY': MultiF2PyCodeWrapper, 'CYTHON': MultiCythonCodeWrapper,
        'DUMMY': DummyWrapper}
    return wrappers[backend.upper()]

# multi_autowrap is essentially modified autowrap to support higher dimensional array
def multi_autowrap(
    exprs, language='C', backend='Cython', tempdir=None, args=None, flags=[],
        verbose=False, helpers=[]):
    """Generates python callable binaries based on the math expression.

    exprs
        The SymPy expressions that should be wrapped as a binary routine

    :Optional arguments:

    language
        The programming language to use, currently 'C' or 'F95'
    backend
        The wrapper backend to use, currently f2py or Cython
    tempdir
        Path to directory for temporary files.  If this argument is supplied,
        the generated code and the wrapper input files are left intact in the
        specified path.
    args
        Sequence of the formal parameters of the generated code, if ommited the
        function signature is determined by the code generator.
    flags
        Additional option flags that will be passed to the backend
    verbose
        If True, autowrap will not mute the command line backends.  This can be
        helpful for debugging.
    helpers
        Used to define auxillary expressions needed for the main expr.  If the
        main expression need to do call a specialized function it should be put
        in the ``helpers`` list.  Autowrap will then make sure that the compiled
        main expression can link to the helper routine.  Items should be tuples
        with (<funtion_name>, <sympy_expression>, <arguments>).  It is
        mandatory to supply an argument sequence to helper routines.

    >>> from sympy.abc import x, y, z
    >>> from sympy.utilities.autowrap import autowrap
    >>> exprs = [((x - y + z)**(13)).expand(),
                 ((x / (1 + y + z))**(3)).expand()]
    >>> binary_func = multi_autowrap(exprs)
    >>> binary_func(1, 4, 2)
    array([-1.0, 0.0029154518950437313])

    """

    code_generator = get_code_generator(language, "autowrap")
    CodeWrapperClass = _get_code_wrapper_class(backend)
    code_wrapper = CodeWrapperClass(code_generator, tempdir, flags, verbose)
    routines = []
    for i, expr in enumerate(exprs):
        try:
            routines.append(Routine('autofunc'+str(i), expr, args))
        except CodeGenArgumentListError, e:
            # if all missing arguments are for pure output, we simply attach them
            # at the end and try again, because the wrappers will silently convert
            # them to return values anyway.
            new_args = []
            for missing in e.missing_args:
                if not isinstance(missing, OutputArgument):
                    raise
                new_args.append(missing.name)
            routines.append(Routine('autofunc', expr, args + new_args))

    helps = []
    for name, expr, args in helpers:
        helps.append(Routine(name, expr, args))

    return code_wrapper.wrap_code(routines, helpers=helps)


def multi_ufuncify(args, exprs, **kwargs):
    """
    This function mimics ufuncify in the autowrap module
    with the difference that all arguments are assumed to
    be arrays (of equal length).

    Furthermore it wraps a number of exprs all of which
    depend on (a subset of) the args provided.

    It means that the returned array is 2 dimensional with
    the first dimension representing the index of the
    expressions provided and the second representing the
    index in the argument arrays. Hence the result has
    the dimension:
    (len(exprs),len(args[0]))

    **kwargs is passed onto autowrap
    """
    from sympy.utilities.autowrap import autowrap
    m = C.Dummy('m', integer=True) # denotes: len(args[0]) (resolved by autowrap)
    #n = C.Dummy('n', integer=True) # denotes: len(exprs)
    nexpr = len(exprs)
    y = [C.IndexedBase(C.Dummy('y'+str(k))) for k in range(nexpr)]
    i = C.Idx(C.Dummy('i', integer=True), nexpr)
    j = C.Idx(C.Dummy('j', integer=True), m)
    lmbs = [C.Lambda(args, expr) for expr in exprs]
    f = [implemented_function('f'+str(x), lmb) for x, lmb in enumerate(lmbs)]

    # all arguments accepts arrays
    internal_args = [C.IndexedBase(C.Dummy(a.name)) for a in args]
    return multi_autowrap([C.Equality(y[k][j], f[k](*[a[j] for a in internal_args])) for \
                     k in range(nexpr)], **kwargs)
