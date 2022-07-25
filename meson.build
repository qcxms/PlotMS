# This build file is part of the PlotMS program 
# for plotting QCxMS files

project('plotms',
        'fortran',
        version: '6.2',
        license: 'LGPL-3.0-or-later',
        meson_version: '>=0.55',
        default_options: [
        'buildtype=debugoptimized',
        'default_library=static',
  ],
)
install = not (meson.is_subproject() and get_option('default_library') == 'static')


prog = []
srcs = []
exe_deps = []

# Collect sources
subdir('config')
subdir('src')

#plotms_exe = executable(
#executable('plotms',srcs)
#
# Create library target
plotms_lib = library(
  meson.project_name(),
  sources: srcs,
  dependencies: exe_deps,
)

# Export as dependency
plotms_inc = plotms_lib.private_dir_include()
plotms_dep = declare_dependency(
  link_with: plotms_lib,
  include_directories: plotms_inc,
  dependencies: exe_deps,
)


# Create executable target
plotms_exe = executable(
  meson.project_name(),
  sources: prog,
  dependencies: plotms_dep,
  install: install,
)


