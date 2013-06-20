import os

import logging

import _setup_odepack
import _setup_transform

# ODEPACK
opksrc = './symodesys/odepack/'

# Transform
tfmsrc = './symodesys/transform/'


def main():
    """
    Precompile some sources to object files
    and store in `prebuilt/` directories for
    speeding up meta-programming compilations.
    """
    ori_dir = os.path.abspath(os.curdir)

    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger(__name__)

    # ODEPACK
    cwd_odepack = os.path.join(os.path.abspath(
        os.path.dirname(__file__)), opksrc)
    os.chdir(cwd_odepack)
    _setup_odepack.main(cwd_odepack, logger)
    os.chdir(ori_dir)

    # Transform
    cwd_transform = os.path.join(os.path.abspath(
        os.path.dirname(__file__)), tfmsrc)
    os.chdir(cwd_transform)
    _setup_transform.main(cwd_transform, logger)
    os.chdir(ori_dir)

if __name__ == '__main__':
    main()
