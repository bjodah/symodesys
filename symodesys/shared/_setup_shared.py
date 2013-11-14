# -*- coding: utf-8 -*-

import os
import shutil
from pycompilation import pyx2obj, CCompilerRunner, missing_or_other_newer

def main(cwd, logger):
    for f in ['symodesys_util.c']:
        fpath = os.path.join(cwd, f)
        dst = os.path.join(cwd, 'prebuilt/'+f.split('.')[0]+'.o')
        if missing_or_other_newer(dst, f):
            runner = CCompilerRunner(
                [fpath], dst, run_linker=False,
                cwd=cwd,
                options=['pic', 'warn', 'fast', 'c99', 'lapack'],
                flags=['-DUSE_LAPACK'],
                metadir='prebuilt/',
                logger=logger)
            if os.path.exists(dst):
                # make sure failed compilation kills the party..
                os.unlink(dst)
            runner.run()
        else:
            logger.info("{} newer than {}, did not recompile.".format(
                dst, fpath))
