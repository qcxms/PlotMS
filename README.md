# PlotMS for QCxMS

**Plot** **M**ass **S**pectra (PlotMS) plotting program for the QCxMS program. 

**Installation**

Put the `plotms` [executable](https://github.com/qcxms/PlotMS/releases) into your `$HOME/bin/` folder. 
Put the *.mass_raw.agr* file into your *$HOME* folder. 

**Compiling**

Download the source code. Go into the created folder and run `make`.

**Running the script**

The files *qcxms.res* and *qcxms_cid.res* are produced then QCxMS is finished with the last run. Run `plotms` in your working directory. 
From version 6.0 and higher, **exact masses** are produced.

The following files are produced:
1) *mass.agr* - output in the XMGRACE format. The *m/z* and relative intensities can be found at the bottom of the file.
2) *result.jdx* - output in the JCAMP format.
3) *result.csv* - output in the CSV format. 
the *m/z* and the abundances values can be found. 


For the older QCEIMS program, use plotms version 4.2

**Documentation**

Find the documentation and other useful information for visualization at the [PlotMS docs](https://xtb-docs.readthedocs.io/en/latest/qcxms_doc/qcxms_plot.html).
