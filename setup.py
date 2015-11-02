#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from distutils.core import setup

pkg_name = 'symodesys'

version_ = '0.0.1'

classifiers = [
    "Development Status :: 2 - Pre-Alpha",
    "Intended Audience :: Developers",
    "License :: OSI Approved :: BSD License",
    "Operating System :: OS Independent",
    "Programming Language :: Python",
    "Programming Language :: C",
    "Programming Language :: Cython",
    "Programming Language :: Fortran",
    "Topic :: Software Development :: Code Generators",
    "Topic :: Software Development :: Compilers",
    "Topic :: Software Development :: Libraries :: Python Modules",
    "Topic :: Scientific/Engineering",
    "Topic :: Scientific/Engineering :: Mathematics",
]

if '--help'in sys.argv[1:] or sys.argv[1] in (
        '--help-commands', 'egg_info', 'clean', '--version'):
    cmdclass_ = {}
    ext_modules=ext_modules_,
else:
    from pycodeexport import pce_build_ext
    from symodesys.shared._setup_shared import get_shared_pce_ext
    from symodesys.odepack._setup_odepack import get_odepack_pce_ext
    from symodesys.gsl._setup_gsl import get_gsl_pce_ext
    from symodesys.sundials._setup_sundials import get_sundials_pce_ext

    ext_modules_ = [
        get_shared_pce_ext(pkg_name),
        #get_odepack_pce_ext(pkg_name),
        get_gsl_pce_ext(pkg_name),
        get_sundials_pce_ext(pkg_name),
    ]
    cmdclass_ = {'build_ext': pce_build_ext}


setup(
    name=pkg_name,
    version=version_,
    author='Björn Dahlgren',
    author_email='bjodah@DELETEMEgmail.com',
    description='Convenience functions for use with sympy.',
    license="BSD",
    url='https://github.com/bjodah/'+pkg_name,
    download_url='https://github.com/bjodah/'+pkg_name+'/archive/v'+version_+'.tar.gz',
    packages=['symodesys', 'symodesys.helpers'],
    cmdclass=cmdclass_,
    ext_modules=ext_modules_,
    classifiers=classifiers,
)
