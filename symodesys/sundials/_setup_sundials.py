import os
from pycompilation import pyx2obj, CCompilerRunner, missing_or_other_newer

def main(cwd, logger):
    # Cythonize pyx file
    src = 'drivers_wrapper.pyx'
    dst = 'prebuilt/drivers_wrapper.o'
    pyx2obj(src, dst, cwd=cwd, logger=logger,
            only_update=True, metadir='prebuilt/')

    f = 'drivers.c'
    fpath = os.path.join(cwd, f)
    dst = os.path.join(cwd, 'prebuilt/drivers.o')
    if missing_or_other_newer(dst, f):
        runner = CCompilerRunner(
            [fpath], dst, run_linker=False,
            cwd=cwd,
            options=['pic', 'warn', 'fast', 'c99', 'lapack'],
            libs=['m', 'sundials_cvode', 'sundials_nvecserial'],
            defmacros=['SUNDIALS_DOUBLE_PRECISION'],
            metadir='prebuilt/',
            logger=logger)
        if os.path.exists(dst):
            # make sure failed compilation kills the party..
            os.unlink(dst)
        runner.run()
    else:
        logger.info("{} newer than {}, did not recompile.".format(
            dst, fpath))
