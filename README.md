# PlotMS for QCxMS

**Plot** **M**ass **S**pectra (PlotMS) plotting program for the QCxMS program. 

**Installation**

Put the `plotms` [executable](https://github.com/qcxms/PlotMS/releases) into your `$HOME/bin/` folder. 
Put the *.mass_raw.agr* file into your *$HOME* folder. 

**Compiling**

Download the source code. Go into the created folder and run `make`.

**Running the script**

The files *qcxms.res* and *qcxms_cid.res* are produced then QCxMS is finished with the last run. Run `plotms` in your working directory. This produces the *mass.agr* file in the XMGRACE format. At the bottom of this file, the *m/z* and the abundances values can be found. 

For the older QCEIMS program, use plotms version 4.2

**Documentation**

Find the documentation and other useful information for visualization at the [PlotMS docs](https://xtb-docs.readthedocs.io/en/latest/qcxms_doc/qcxms_plot.html).
