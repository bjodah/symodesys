from __future__ import print_function, division

import os
import subprocess
import pickle

from distutils.spawn import find_executable

# TODO: change print statements to logging statements.

def find_binary_of_command(candidates):
    """
    Currently only support *nix systems (invocation of which)
    """
    for c in candidates:
        binary_path = find_executable(c)
        if c:
            return c, binary_path
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
    metadata_filename = '.metadata_CompilerRunner'

    def __init__(self, sources, out, flags=None, run_linker=True,
                 compiler=None, cwd=None, inc_dirs=None, libs=None,
                 lib_dirs=None,
                 options=None, verbose=False, preferred_vendor=None,
                 metadir=None):
        """
        Arguments:
        - `preferred_vendor`: key of compiler_dict
        """

        self.sources = sources if hasattr(sources,'__iter__') else [sources]
        self.out = out
        self.flags = flags or []
        #self.run_linker = run_linker
        if compiler:
            self.compiler_name, self.compiler_binary = compiler
        else:
            # Find a compiler
            preferred_compiler_name = self.compiler_dict.get(preferred_vendor,None)
            self.compiler_name, self.compiler_binary = self.find_compiler(
                preferred_compiler_name, metadir or cwd)
            if self.compiler_binary == None:
                raise RuntimeError("No compiler found (searched: {})".format(
                    ', '.join(self.compiler_dict.values())))
        self.cwd = cwd
        self.inc_dirs = inc_dirs or []
        self.libs = libs or []
        self.lib_dirs = lib_dirs or []
        self.options = options or []
        self.verbose = verbose
        if run_linker:
            # both gcc and ifort have '-c' flag for disabling linker
            self.flags = filter(lambda x: x != '-c', self.flags)
        else:
            self.flags.append('-c')

        for inc_dir in self.inc_dirs:
            self.flags.append('-I'+inc_dir)

        for lib_dir in self.lib_dirs:
            self.lib_dirs.append('-L'+lib_dir)

        for opt in self.options:
            extra_flags = self.flag_dict[self.compiler_name][opt]
            self.flags.extend(extra_flags)


    @classmethod
    def find_compiler(cls, preferred_compiler_name=None, load_save_choice=None):
        """
        Identify a suitable C/fortran/other compiler

        When it is possible that the user (un)installs a compiler inbetween
        compilations of object files we want to catch that. This method
        allows compiler choice to be stored in a pickled metadata file.
        Provide load_save_choice a dirpath to make the class save choice
        there in a file with cls.metadata_filename as name.
        """
        if load_save_choice:
            name_path = cls.get_from_metadata_file(load_save_choice, 'compiler')
            if name_path != None: return name_path
        candidates = cls.flag_dict.keys()
        if preferred_compiler_name:
            if preferred_compiler_name in candidates:
                # Duplication doesn't matter
                candidates = [preferred_compiler_name] + candidates
        name_path = find_binary_of_command(candidates)
        if load_save_choice:
            cls.save_to_metadata_file(load_save_choice, 'compiler', name_path)
        return name_path


    @classmethod
    def _get_metadata_key(cls, kw):
        """ kw could be e.g. 'compiler' """
        return cls.__name__+'_'+kw

    @classmethod
    def get_from_metadata_file(cls, dirpath, key):
        """
        Get value of key in metadata file dict.
        """
        fullpath = os.path.join(dirpath, cls.metadata_filename)
        if os.path.exists(fullpath):
            d = pickle.load(open(fullpath,'r'))
            return d.get(cls._get_metadata_key(key), None)
        else:
            # Raise an exception instead?
            return None

    @classmethod
    def save_to_metadata_file(cls, dirpath, key, value):
        """
        Store `key: value` in metadata file dict.
        """
        fullpath = os.path.join(dirpath, cls.metadata_filename)
        if os.path.exists(fullpath):
            d = pickle.load(open(fullpath,'r'))
            d.update({key: value})
            pickle.dump(d, open(fullpath,'w'))
        else:
            print({key: value}, fullpath)
            pickle.dump({key: value}, open(fullpath,'w'))

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

    compiler_dict = {
        'gnu': 'gcc',
        'intel': 'icc',
    }

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

    compiler_dict = {
        'gnu': 'gfortran',
        'intel': 'ifort',
    }

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
