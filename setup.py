import os

import logging

from pycompilation.util import run_sub_setup

from symodesys.shared._setup_shared import main as shared_main
from symodesys.odepack._setup_odepack import main as odepack_main
from symodesys.gsl._setup_gsl import main as gsl_main
from symodesys.sundials._setup_sundials import main as sundials_main


def main():
    """
    Precompile some sources to object files
    and store in `prebuilt/` directories for
    speeding up meta-programming compilations.
    """
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger(__name__)

    for name, cb in zip(['shared', 'odepack', 'gsl', 'sundials'],
                        [shared_main, odepack_main, gsl_main, sundials_main]):

        cwd = os.path.join(os.path.abspath(
            os.path.dirname(__file__)),
                           'symodesys/'+name+'/')
        run_sub_setup(cwd, cb, logger)


if __name__ == '__main__':
    main()
