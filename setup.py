import os

import logging

from pycompilation.util import run_sub_setup

from symodesys.odepack._setup_odepack import main as odepack_main
from symodesys.gsl._setup_gsl import main as gsl_main


def main():
    """
    Precompile some sources to object files
    and store in `prebuilt/` directories for
    speeding up meta-programming compilations.
    """
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger(__name__)

    # ODEPACK
    cwd = os.path.join(os.path.abspath(
        os.path.dirname(__file__)),
                       'symodesys/odepack/')
    run_sub_setup(cwd, odepack_main, logger)

    # GSL
    cwd = os.path.join(os.path.abspath(
        os.path.dirname(__file__)),
                       'symodesys/gsl/')
    run_sub_setup(cwd, gsl_main, logger)

if __name__ == '__main__':
    main()
