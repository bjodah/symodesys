from __future__ import print_function, division
import subprocess


def find_binary_of_command(candidates):
    for c in candidates:
        p = subprocess.Popen(['which', c], stdout=subprocess.PIPE)
        stdoutdata, stderrdata = p.communicate()
        p.wait()
        if p.returncode == 0:
            chosen_candidate = c
            binary_path = stdoutdata.split('\n')[0]
            return chosen_candidate, binary_path
    raise RuntimeError('No binary located for candidates: {}'.format(
        candidates))


def _uniquify(l):
    result = []
    for x in l:
        if not x in result:
            result.append(x)
    return result

class CompilerRunner(object):

    flag_dict = None # Lazy unified defaults for compilers

    def __init__(self, sources, out, flags=None, run_linker=True,
                 compiler=None, cwd=None, inc_dirs=None, libs=None,
                 lib_dirs=None,
                 options=None, verbose=False):

        self.sources = sources
        self.out = out
        self.flags = flags or []
        #self.run_linker = run_linker
        if compiler:
            self.compiler_name, self.compiler_binary = compiler
            assert self.compiler_name in self.flag_dict
        else:
            # Find a compiler
            self.compiler_name, self.compiler_binary = self.find_compiler()
        self.cwd = cwd
        self.inc_dirs = inc_dirs or []
        self.libs = libs or []
        self.lib_dirs = lib_dirs or []
        self.options = options or []
        self.verbose = verbose
        if run_linker:
            # both gcc and ifort have '-c' flag for disabling linker
            self.flags = filter(lambda x: x == '-c', self.flags)
        else:
            self.flags.append('-c')

        for inc_dir in self.inc_dirs:
            self.flags.append('-I'+inc_dir)

        for lib_dir in self.lib_dirs:
            self.lib_dirs.append('-L'+lib_dir)

        for opt in self.options:
            extra_flags = self.flag_dict[self.compiler_name][opt]
            self.flags.extend(extra_flags)


    def find_compiler(self):
        """
        Identify a suitable fortran compiler
        Currently only support *nix systems
        """
        return find_binary_of_command(self.flag_dict.keys())


    def run(self):
        self.flags = _uniquify(self.flags)

        # Append output flag and name to tail of flags
        self.flags.extend(['-o', self.out])

        cmd = [self.compiler_binary]+self.flags+self.sources+['-l'+x for x in self.libs]
        if self.verbose: print('Executing... : {}'.format(' '.join(cmd)))
        p = subprocess.Popen(cmd,
                             cwd=self.cwd,
                             #shell=True,
                             stdin= subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE
        )
        stdoutdata, stderrdata = p.communicate()
        return stdoutdata, stderrdata, p.returncode


class CCompilerRunner(CompilerRunner):

    flag_dict = {
        'gcc': {
            'pic': ('-fPIC',),
            'warn': ('-Wall', '-Wextra', '-Wimplicit-interface'),
            'fast': ('-O3', '-march=native', '-ffast-math', '-funroll-loops'),
        },
        'icc': {
            'pic': ('-fPIC',),
            'fast': ('-fast',),
        }
    }


class FortranCompilerRunner(CompilerRunner):

    flag_dict = {
        'gfortran': {
            'f90': ('-std=f2008',),
        },
        'ifort': {
            'warn': ('-warn', 'all',),
            'f90': ('-stand f95',),
        }
    }


    def __init__(self, *args, **kwargs):
        # gfortran takes a superset of gcc arguments
        for key, value in self.flag_dict.items():
            value.update(CCompilerRunner.flag_dict[{'gfortran': 'gcc', 'ifort': 'icc'}.get(key)])
        super(FortranCompilerRunner, self).__init__(*args, **kwargs)
