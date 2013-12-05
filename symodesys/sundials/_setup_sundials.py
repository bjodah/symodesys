import os
from pycompilation import pyx2obj, src2obj

USE_LAPACK = False

def main(dst, **kwargs):
    # Cythonize pyx file

    defmacros = ['SUNDIALS_DOUBLE_PRECISION']
    options = ['pic', 'warn', 'fast']
    if USE_LAPACK:
        defmacros += ['USE_LAPACK']
        options += ['lapack']

    return [
        pyx2obj('drivers_wrapper.pyx', dst,
                only_update=True, metadir=dst, **kwargs),
        src2obj(
            'drivers.c',
            objpath=dst,
            options=options,
            std='c99',
            libs=['m', 'sundials_cvode', 'sundials_nvecserial'],
            defmacros=defmacros,
            metadir=dst,
            only_update=True,
            **kwargs)
    ]
