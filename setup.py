import os

import logging

import _setup_odepack
import _setup_transform
import _setup_gsl

def _run_sub_setup(path, module, logger):
    ori_dir = os.path.abspath(os.curdir)
    cwd = os.path.join(os.path.abspath(
        os.path.dirname(__file__)), path)
    os.chdir(cwd)
    module.main(cwd, logger)
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
    _run_sub_setup('./symodesys/odepack/', _setup_odepack, logger)

    # Transform
    _run_sub_setup('./symodesys/transform/', _setup_transform, logger)

    # GSL
    _run_sub_setup('./symodesys/gsl/', _setup_gsl, logger)

if __name__ == '__main__':
    main()
