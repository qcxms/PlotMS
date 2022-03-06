# This file is part of PlotMS.

os = host_machine.system()
fc = meson.get_compiler('fortran')
fc_id = fc.get_id()
static = get_option('default_library') == 'static' and get_option('fortran_link_args').contains('-static')

#> choose compiler options
if fc_id == 'gcc'
  add_project_arguments(
    '-ffree-line-length-none',
    '-fbacktrace',
    language: 'fortran',
  )
  if fc.version().version_compare('<8')
    error('GCC version 8 or higher is required.')
  endif
elif fc_id == 'intel'
  add_project_arguments(
    '-traceback',
    language: 'fortran',
  )
  add_project_arguments(
    '-DLINUX',
    language: 'c',
  )                                                                                                                                                                                                                                                                    
elif fc_id == 'pgi' or fc_id == 'nvidia_hpc'
  add_project_arguments(
    '-Mbackslash',
    '-Mallocatable=03',
    '-traceback',
    language: 'fortran',
  )
  error('This project does not support compilation with NVHPC yet')
endif
########################################################################


# Create the tool chain library as subproject
mctc_dep = dependency(
  'mctc-lib',
  version: '>=0.2',
  fallback: ['mctc-lib', 'mctc_dep'],
  default_options: ['default_library=static'],
  static: static,
)
exe_deps += mctc_dep
