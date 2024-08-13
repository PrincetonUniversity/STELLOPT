Tutorial: Solovev NESCOIL Example
=============================================

This tutorial will walk the user through running the [NESCOIL](NESCOIL.nd) code for
generating a coil based on a Solovev equilibrium. Before beginning the user should run the [BNORM Solovev tutorial](BNORM Solovev Example.md) to generate a bnorm and
nescin file.

------------------------------------------------------------------------

**Run the NESCOIL code**

From the directory where the VMEC nescin and bnorm files are located 
execute the NESCOIL code.

```
>~/bin/xnescoil nescin.solovev_nfp5
```

**Examine the output.**

The NESCOIL code will generate a `nescout.solovev_nfp5` file containing
the output current potential of the nescoil code. The results can be plotted
using the pySTEL python routine `nescoil_util.py`.

![](images/nescout_solovev_nfp5.png)

![](images/nescout_3d_solovev_nfp5.png)

**Cut coils**

The python utility `nescoil_util.py` can be used to cut coils from the
NESCOIL current potential.

![](images/nescout_solovev_nfp5_coilcut.png)

![](images/nescout_coil_solovev_nfp5.png)
