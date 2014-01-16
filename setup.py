#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os

pkg_name = 'symodesys'

version_ = '0.0.1'

if '--help'in sys.argv[1:] or sys.argv[1] in (
        '--help-commands', 'egg_info', 'clean', '--version'):
    cmdclass_ = {}
    sub_pkgs = []
else:

    from pycompilation.util import make_dirs
    from symodesys.shared._setup_shared import prebuild as shared_prebuild
    from symodesys.odepack._setup_odepack import prebuild as odepack_prebuild
    from symodesys.gsl._setup_gsl import prebuild as gsl_prebuild
    from symodesys.sundials._setup_sundials import prebuild as sundials_prebuild

    sub_folders = ['shared', 'odepack', 'gsl', 'sundials']
    prebuilds = zip(sub_folders,
                [shared_prebuild, odepack_prebuild, gsl_prebuild, sundials_prebuild])
    sub_pkgs = ['symodesys.' + x for x in sub_folders]

    def run_prebuilds(build_lib):
        """
        Precompile some sources to object files
        and store in `prebuilt/` directories for
        speeding up meta-programming compilations.
        """
        import logging
        logging.basicConfig(level=logging.DEBUG)
        logger = logging.getLogger(__name__)

        for name, cb in prebuilds:
            destdir = os.path.join(build_lib, pkg_name, name)
            prebuilt_destdir = os.path.join(destdir, 'prebuilt')
            if not os.path.exists(prebuilt_destdir): make_dirs(prebuild_destdir)
            srcdir = os.path.join(os.path.dirname(__file__), pkg_name, name)
            cb(srcdir, destdir, logger=logger)

    from distutils.command.build_py import build_py as _build_py
    from distutils.core import setup

    class build_py(_build_py):
        """Specialized Python source builder."""
        def run(self):
            if not self.dry_run:
                run_prebuilds(self.build_lib)
            _build_py.run(self)

    cmdclass_ = {'build_py': build_py}

setup(
    name=pkg_name,
    version=version_,
    author='Björn Dahlgren',
    author_email='bjodah@DELETEMEgmail.com',
    description='Convenience functions for use with sympy.',
    license = "BSD",
    url='https://github.com/bjodah/'+pkg_name,
    download_url='https://github.com/bjodah/'+pkg_name+'/archive/v'+version_+'.tar.gz',
    packages=['symodesys', 'symodesys.helpers'] + sub_pkgs,
    # package_data={pkg_name: [
    #     'symodesys/'
    # ]}
    cmdclass = cmdclass_
)
