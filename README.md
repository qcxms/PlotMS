# PlotMS for QCxMS

**Plot** **M**ass **S**pectra (PlotMS) plotting program for the [QCxMS program](https://github.com/qcxms/QCxMS). 

**Installation**

Put the `plotms` [executable](https://github.com/qcxms/PlotMS/releases) and the `.mass_raw.agr` files into the `$HOME` folder. 

**Compiling**

Using **meson** as build system requires you to install a fairly new version like 0.57.2 or newer. 
To use the default backend of meson you have to install **ninja** version 1.10 or newer.

```bash
export FC=ifort CC=icc
meson setup build -Dfortran_link_args=-static
ninja -C build 
```
Copy the binary from the *build/plotms* file into a directory in your path, e.g. *~/bin/*.


**Running the script**

The files `qcxms.res` and `qcxms_cid.res` are produced then QCxMS is finished with the last run. Run `plotms` in your working directory. 
From version 6.0 and higher, **exact masses** are produced.

The following files are produced:
1) *mass.agr* - output in the XMGRACE format. The *m/z* and relative intensities can be found at the bottom of the file.
2) *result.jdx* - output in the JCAMP format.
3) *result.csv* - output in the CSV format. 
the *m/z* and the abundances values can be found. 

For the older QCEIMS program, use plotms version 4.2. Calculations done with QCEIMS can still be used with newer versions by renaming the `qceims.res` to `qcxms.res` or `qceims_cid.res` to `qcxms_cid.res` .  

**Documentation**

Find the documentation and other useful information for visualization at the [PlotMS docs](https://xtb-docs.readthedocs.io/en/latest/qcxms_doc/qcxms_plot.html).
