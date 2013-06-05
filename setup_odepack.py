from __future__ import print_function, division

import os
import sys
import subprocess
import shutil

from symodesys.helpers import md5_of_file
from symodesys.helpers.compilation import FortranCompilerRunner, CCompilerRunner

optimize = True

opkfiles = ['opkdmain.f', 'opkda1.f', 'opkda2.f']
opksrc = './symodesys/odepack/'
opkweb = 'https://computation.llnl.gov/casc/odepack/software/'
md5sums = {
    'opkda1.f':    '00a675f71ab375376bb6108d24a33c0b',
    'opkda2.f':    'dd03a71ea1a5ac746169c0279aa4c551',
    'opkdmain.f':  '47d81cc73a1e82210f47a97c43daa8cf'
}

# Download sources ----------------------------------------
for f in opkfiles:
    if not os.path.exists(opksrc+f):
        import urllib2
        print('Downloading: {}'.format(opkweb+f))
        open(opksrc+f, 'wt').write(urllib2.urlopen(opkweb+f).read())
    fmd5 = md5_of_file(opksrc+f).hexdigest()
    if fmd5 != md5sums[f]:
        print("""Warning: MD5 sum of {} differs from that provided in setup.py.
        i.e. {} vs. {}""".format(opksrc+f, fmd5, md5sums[f]))


# (Pre)compile sources ----------------------------------------

# Distutils does not allow to use .o files in compilation
# (see http://bugs.python.org/issue5372)
# hence the compilation of ODEPACK is done once and for all and
# saved in prebuilt dir

cwd = os.path.join(os.path.abspath(
    os.path.dirname(__file__)), opksrc)
cur_dir = os.path.abspath(os.curdir)
os.chdir(cwd)
for f in opkfiles:
    # compiles to: prebuilt/opkd{a1,a2,main}.o
    name, ext = os.path.splitext(f)
    dst = os.path.join('prebuilt', name+'.o') # .ext -> .o
    if not os.path.exists(os.path.join(cwd, dst)):
        runner = FortranCompilerRunner([f], dst, run_linker=False, cwd=cwd,
            options=['pic', 'warn', 'fast'], verbose=True)
        out, err, exit_status = runner.run()
        if exit_status != 0:
            print(out)
            print(err)
        else:
            print('...Success!')
    else:
        print("Found {}, did not recompile.".format(dst))


# Cythonize pyx file
dst = 'pylsodes_bdf.c'
if not os.path.exists(os.path.join(cwd, dst)):
    from Cython.Compiler.Main import default_options, compile, CompilationOptions
    cy_sources = ['pylsodes_bdf.pyx']
    cy_options = CompilationOptions(default_options)
    compile(cy_sources, cy_options)
else:
    print("Found {}, did not re-cythonize.".format(dst))

# Compile cython generated .c file to object file
dst = 'prebuilt/pylsodes_bdf.o'
if not os.path.exists(os.path.join(cwd, dst)):
    from distutils.sysconfig import get_python_inc, get_config_vars
    import numpy
    includes = [get_python_inc(), numpy.get_include()]
    cc = " ".join(get_config_vars('CC', 'BASECFLAGS', 'OPT', 'CFLAGSFORSHARED'))
    compilern, flags = cc.split()[0], cc.split()[1:]
    print(includes)
    runner =CCompilerRunner(['pylsodes_bdf.c'], dst, flags, run_linker=False,
                            compiler=[compilern]*2, cwd=cwd,
                            inc_dirs=includes, verbose=True)
    out, err, exit_status = runner.run()
    if exit_status != 0:
        print(out)
        print(err)
    else:
        print('...Success!')
else:
    print("Found {}, did not recompile.".format(dst))
os.chdir(cur_dir)
