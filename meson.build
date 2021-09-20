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
       'plotms_isotopes.f90',
       'qcxms_boxmuller.f90',
       'qcxms_readl.f90',
       'xtb_mctc_accuracy.f90'
       ]

#plotms_exe = executable(
executable('plotms',srcs)


