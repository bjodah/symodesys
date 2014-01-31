# -*- coding: utf-8 -*-

import os
from pycompilation.codeexport import make_CleverExtension_for_prebuilding_Code, C_Code

prebuild_sources = ['symodesys_util.c']

class SymodesysUtilCode(C_Code):
    compile_kwargs = {'std': 'c99'}
    source_files = prebuild_sources
    build_files = 'main_ex_template.c symodesys_util.h symodesys_util.c'.split()


def get_shared_clever_ext(basename):
    return make_CleverExtension_for_prebuilding_Code(
        basename+'.shared._FOO', SymodesysUtilCode,
        prebuild_sources,
        srcdir=os.path.join(basename, 'shared'),
        logger=True
    )
