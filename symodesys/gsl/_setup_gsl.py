import os
from pycompilation import pyx2obj, src2obj

def main(dst, **kwargs):
    return [
        pyx2obj('drivers_wrapper.pyx', dst, only_update=True,
                metadir=dst, **kwargs),
        src2obj(
            'drivers.c',
            objpath=dst,
            defmacros=['GSL_RANGE_CHECK_OFF', 'HAVE_INLINE'],
            options=['pic', 'warn', 'fast'],
            std='c99',
            metadir=dst,
            only_update=True,
            **kwargs
        )
    ]
