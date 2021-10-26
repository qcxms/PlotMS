# This build file is part of the PlotMS program 
# for plotting QCxMS files

project('PlotMS',
        'fortran',
        version: '6.0',
        license: 'LGPL-3.0-or-later',
        meson_version: '>=0.55',
        default_options: [
        'buildtype=debugoptimized',
        'default_library=static',
  ],
)

srcs = ['PlotMS_v6.0.f90',
       'isotopes.f90',
       'boxmuller.f90',
       'bucket.f90',
       'readcommon.f90',
       'readl.f90',
       'accuracy.f90'
       ]

#plotms_exe = executable(
executable('plotms',srcs)


