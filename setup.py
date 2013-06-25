import os

import logging

from symodesys.odepack._setup_odepack import main as odepack_main
from symodesys.transform._setup_transform import main as transform_main
from symodesys.gsl._setup_gsl import main as gsl_main

def _run_sub_setup(path, cb, logger):
    ori_dir = os.path.abspath(os.curdir)
    cwd = os.path.join(os.path.abspath(
        os.path.dirname(__file__)), path)
    os.chdir(cwd)
    cb(cwd, logger)
    os.chdir(ori_dir)


def main():
    """
    Precompile some sources to object files
    and store in `prebuilt/` directories for
    speeding up meta-programming compilations.
    """
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger(__name__)

    # ODEPACK
    _run_sub_setup('./symodesys/odepack/', odepack_main, logger)

    # Transform
    _run_sub_setup('./symodesys/transform/', transform_main, logger)

    # GSL
    _run_sub_setup('./symodesys/gsl/', gsl_main, logger)

if __name__ == '__main__':
    main()
