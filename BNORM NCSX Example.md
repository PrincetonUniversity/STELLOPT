Tutorial: NCSX BNORM Example
=============================================

This tutorial will walk the user through running the [BNORM](BNORM.nd) code for
generating a bnorm file based on NCSX. Before beginning the user should run the
[VMEC fixed boundary tutorial](VMEC Fixed Boundary Run.md) tutorial to 
generate an input and wout file.

------------------------------------------------------------------------

**Run the BNORM code**

From the directory where the VMEC wout file is located execute the BNORM
code with an offset boundary of 0.4 m.

```
>~/bin/xbnorm wout_ncsx_c09r00_fixed.nc 0.4
```

**Examine the output.**

The BNORM code will generate a `bnorm.ncsx_c09r00_fixed` file containing
the sine harmonics of the plasma B-field normal magnetic field component 
on the plasma surface. These harmonics are nomalized to the total polidal
current in the equilibrium. A `nescin.ncsx_c09r00_fixed` NESCOIL input file
will also be generated.
