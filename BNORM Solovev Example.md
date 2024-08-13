Tutorial: SOLOVEV BNORM Example
=============================================

This tutorial will walk the user through running the [BNORM](BNORM.nd) code for
generating a bnorm file based on a Solovev equilibrium. The wout file
can be found here: [wout_solovev_nfp5.nc](examples/wout_solovev_nfp5.nc)

------------------------------------------------------------------------

**Run the BNORM code**

From the directory where the VMEC wout file is located execute the BNORM
code with an offset boundary of 0.4 m.

```
>~/bin/xbnorm wout_solovev_nfp5.nc 0.4
```

**Examine the output.**

The BNORM code will generate a `bnorm.solovev_nfp5` file containing
the sine harmonics of the plasma B-field normal magnetic field component 
on the plasma surface. These harmonics are nomalized to the total polidal
current in the equilibrium. A `nescin.solovev_nfp5` NESCOIL input file
will also be generated.
