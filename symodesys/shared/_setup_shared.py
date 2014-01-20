# -*- coding: utf-8 -*-

from pycompilation.codeexport import C_Code, prebuild_Code

prebuild_sources = ['symodesys_util.c']

class SymodesysUtilCode(C_Code):
    compile_kwargs = {'std': 'c99'}
    source_files = prebuild_sources
    copy_files = 'main_ex_template.c symodesys_util.h symodesys_util.c'.split()

def prebuild(srcdir, destdir, build_temp, **kwargs):
    return prebuild_Code(
        srcdir, destdir, build_temp, SymodesysUtilCode,
        prebuild_sources, **kwargs)
