NESCOIL
=======

The NESCOIL (NEumann Solver for fields produced by externals COILs)
([P. Merkel 1987 //Nucl. Fusion// \*\*27\*\* 867](http://dx.doi.org/10.1088/0029-5515/27/5/018))
code calculates a surface current on the exterior surface of two
toroidally closed surfaces such that the normal field on the interior
surface is minimized.

------------------------------------------------------------------------

### Theory

Here\'s the section to explain the theory behind the code.

------------------------------------------------------------------------

### Compilation

Here\'s the section to explain how to compile the code.

------------------------------------------------------------------------

### Input Data Format

The NESCOIL code takes an input file \'nescin.ext\' where ext is the
extension appended to each run. An input file is generated with each run
of the [BNORM](BNORM) code. There must also be a bnorm file with the
same extension in the directory from which the code is run. The input
file has the following format

    ------ Spatial dimensions ----
    nu, nv, nu1, nv1, npol, ntor, lasym_bn
     180 180 180 180 20 10 F

    ------ Fourier Dimensions ----
    mf, nf, md, nd (max in surf and bnorm files)
     10 10 20 20

    ------ Plasma information from VMEC ----
    np     iota_edge       phip_edge       curpol
     3 0.65460308521511945 -8.18670745572676883E-2 4.9781746242537697

    ------ Current Controls ----
    cut  cup  ibex(=1,use fixed background coils)
     0.E+0 1. 0

    ------ SVD controls -----
    mstrt, mstep, mkeep, mdspw, curwt, trgwt
     0 0 0 4 0.E+0 0.E+0

    ------ Output controls -----
    w_psurf w_csurf w_bnuv w_jsurf w_xerr w_svd
     1 1 1 1 1 1

    ------ Plasma Surface ----
    Number of fourier modes in table
     94
    Table of fourier coefficients
    m,n,crc,czs,cls,crs,czc,clc
          0     0  1.378200000000E+00  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
          0    -1 -4.145200000000E-03  7.163400000000E-03  3.337515916590E-02  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00

The following table explains each parameter as defined by the preceding
line:

\|\| Input Parameter Name \|\| Description \|\| \|\| nu \|\| Number of
poloidal grid points on current surface \|\| \|\| nv \|\| Number of
toroidal grid points on current surface \|\| \|\| nu1 \|\| Number of
poloidal grid points on inner surface \|\| \|\| nv1 \|\| Number of
toroidal grid points on inner surface \|\| \|\| npol \|\| Number of
segments of a modular or helical filament \|\| \|\| ntor \|\| Number of
filaments per period \|\| \|\| mf \|\| Number of poloidal modes for
current potential \|\| \|\| nf \|\| Number of toroidal modes for current
potential \|\| \|\| md \|\| Number of poloidal modes for surface shapes
(and B normal) \|\| \|\| nd \|\| Number of toroidal modes for surface
shapes (and B normal) \|\| \|\| np \|\| Number of field periods \|\|
\|\| iota\_edge \|\| Equilibrium rotational transform at edge \|\| \|\|
phip\_edge \|\| Equilibrium toroidal flux derivative at edge \|\| \|\|
curpol \|\| Equilibrium total poloidal current in Amps per field period
\|\| \|\| cut \|\| Net toroidal current +/-1 or 0 \|\| \|\| cup \|\| Net
poloidal current +/-1 or 0 \|\| \|\| ibex \|\| External field supplied
if IBEX = 1 \|\| \|\| mstrt \|\| Method + svdscan start if \> 1, \>=0
Berr, \<=0 Least square \|\| \|\| mstep \|\| Method + svdscan stepsize
\<=0 least square, =0 use F04ABE, no svd \|\| \|\| mkeep \|\| svd/scan
control 0 svdscan, else keep nkeep wgts \|\| \|\| mdspw \|\| 2 +
exponent of dsur multiplying bfn, ben \|\| \|\| curwt \|\| Weight for
surface current minimization (only in LSQ branch) \|\| \|\| trgwt \|\|
Not yet implemented \|\| \|\| w\_psurf \|\| Write plasma surface info
\|\| \|\| w\_csurf \|\| Write coil surface info \|\| \|\| w\_bnuv \|\|
Write Bnorm field info \|\| \|\| w\_jsurf \|\| Write J surface current
info \|\| \|\| w\_xerr \|\| Write X error (displacement) info \|\| \|\|
w\_svd \|\| Write SVD solution info \|\|

------------------------------------------------------------------------

### Execution

The code is executed from the command line by passing the name of the
\'nescin.ext\' file to the code

------------------------------------------------------------------------

### Output Data Format

A text file is output with the name \'nescout.ext\' where \'ext\' in the
same extension as the input file which generated the output. In the
file, the various input parameters are output in tables along with the
Fourier harmonics of the surface potential. Setting the w\_psurf,
w\_csurf, w\_bnuv, w\_jsurf, w\_xerr, and w\_svd flags can modify the
parameters which are output to this file. This file is essentially a
self-documenting text file.

------------------------------------------------------------------------

### Visualization

Explain how to visualize the data.

------------------------------------------------------------------------

### Tutorials

Put links to tutorial pages here.
