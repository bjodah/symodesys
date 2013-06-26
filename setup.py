import os

import logging

from pycompilation.helpers import run_sub_setup

from symodesys.odepack._setup_odepack import main as odepack_main
from symodesys.transform._setup_transform import main as transform_main
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
    run_sub_setup('./symodesys/odepack/', odepack_main, logger)

    # GSL
    run_sub_setup('./symodesys/gsl/', gsl_main, logger)

if __name__ == '__main__':
    main()
