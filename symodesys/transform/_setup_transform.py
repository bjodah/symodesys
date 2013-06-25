from pycompilation import pyx2obj

def main(cwd, logger):
    # Cythonize pyx file
    src = 'transform_wrapper.pyx'
    dst = 'prebuilt/transform_wrapper.o'
    pyx2obj(src, dst, cwd=cwd, logger=logger, only_update=True)
