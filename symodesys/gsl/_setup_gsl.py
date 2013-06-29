import os
from pycompilation import pyx2obj, CCompilerRunner, missing_or_other_newer

def main(cwd, logger):
    # Cythonize pyx file
    src = 'drivers_wrapper.pyx'
    dst = 'prebuilt/drivers_wrapper.o'
    pyx2obj(src, dst, cwd=cwd, logger=logger, only_update=True)

    f = 'drivers.c'
    fpath = os.path.join(cwd, f)
    dst = 'prebuilt/drivers.o'
    if missing_or_other_newer(dst, f):
        # Intel Fortran fails for opkda1.f, hence prefer `gnu`
        runner = CCompilerRunner(
            [fpath], dst, run_linker=False,
            cwd=cwd, options=['pic', 'warn', 'fast'],
            preferred_vendor='gnu', metadir='prebuilt/',
            logger=logger)
        runner.run()
    else:
        logger.info("{} newer than {}, did not recompile.".format(dst, fpath))