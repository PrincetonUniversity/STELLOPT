NESCOIL
=======

The NESCOIL (NEumann Solver for fields produced by externals COILs)
([P. Merkel 1987 //Nucl. Fusion// \*\*27\*\* 867](http://dx.doi.org/10.1088/0029-5515/27/5/018))
code calculates a surface current on the exterior surface of two
toroidally closed surfaces such that the normal field on the interior
surface is minimized.

------------------------------------------------------------------------

### Theory

We seek a current potential on a toroidal surface (D) which results in
the existance of an enclosed flux surface (R).

The field produced by current potential can be written as

$$ \vec{B}\left(\vec{x}\right) = \frac{\mu_0}{4\pi}\int\frac{\vec{j}\times\left(\vec{x}-\vec{x}'\right)}{\left|\vec{x}-\vec{x}'\right|^3}d^3x'$$

where we define the current density by a current potential $$\vec{j} = \hat{n}\times\nabla\Phi$$.
Which can be separated into secular and non-secular parts

$$\Phi\left(u,v\right) = \sum_{m=0}^M\sum_{n=-N}^N \Phi_{mn} sin\left(2\pi mu+2\pi nv\right) - \frac{I_{pol}}{N_p}v - I_{tor}u$$

Here $$I_{pol}$$ is the total poloidal current per field period of the equilbrium and
$$I_{tor}$$ is the total toroidal current. Note that the argument to
sine is not the same as that of [VMEC](VMEC). The plasma surface and 
current potential surface also use this formulation and not the [VMEC](VMEC)
one.

The two scalar parts of the potential represent magnetic fields arrising
from a toroidal field ($$\frac{I_{pol}}{N_p}v$$) and a vertical field
($$I_{tor}u$$). While these two quantities are free parameters, the
BNORM code normalizes to the poloidal current per field period, therefore
$$\frac{I_{pol}}{N_p}=1$$ is assumed. However, the user is free to modify
`CUP` and `CUT` to meet their needs.

------------------------------------------------------------------------

### Compilation

NESCOIL is distributed as part of the STELLOPT package of codes through
Git.

------------------------------------------------------------------------

### Input Data Format

The NSCOIL code takes an 'nescin.ext' file as input and optionally a
similarly named 'bnorm.ext' file (as produced by the [BNORM](BNORM) code).
If no 'bnorm.ext' file is found then it is assumed that the normal
plasma field is zero (vacuum condition). The 'nescin' file has the following
format

```
------ Spatial dimensions ----
nu, nv, nu1, nv1, npol, ntor, lasym_bn
         256         256         256         256          64          10 F

------ Fourier Dimensions ----
mf, nf, md, nd (max in surf and bnorm files)
          24          15          24          22

------ Plasma information from VMEC ----
np     iota_edge       phip_edge       curpol
           4   0.0000000000000000       0.32524904170259739       -16.089222854438592     

------ Current Controls ----
cut  cup  ibex(=1,use fixed background coils)
   0.0000000000000000        1.0000000000000000                0

------ SVD controls -----
mstrt, mstep, mkeep, mdspw, curwt, trgwt
           0           0           0           4   0.0000000000000000        0.0000000000000000     

------ Output controls -----
w_psurf w_csurf w_bnuv w_jsurf w_xerr w_svd
           0           0           0           0           0           0

------ Plasma Surface ---- 
Number of fourier modes in table
          10
Table of fourier coefficients
m,n,crc,czs,cls,crs,czc,clc
      0     0  4.000000000000E+00  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      1     0  1.000000000000E+00  1.000000000000E+00 -1.925732563734E-01  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      2     0  0.000000000000E+00  0.000000000000E+00 -5.090906726166E-02  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      3     0  0.000000000000E+00  0.000000000000E+00  6.458084248335E-03  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      4     0  0.000000000000E+00  0.000000000000E+00 -5.915484436642E-03  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      5     0  0.000000000000E+00  0.000000000000E+00  7.030074268844E-04  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      6     0  0.000000000000E+00  0.000000000000E+00 -1.024498213757E-03  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      7     0  0.000000000000E+00  0.000000000000E+00 -1.091616780958E-04  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      8     0  0.000000000000E+00  0.000000000000E+00  6.174851363263E-06  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      9     0  0.000000000000E+00  0.000000000000E+00 -1.438747872521E-06  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00

------ Current Surface: Coil-Plasma separation =   1.000000000000E+00 -----
Number of fourier modes in table
           2
Table of fourier coefficients
m,n,crc2,czs2,crs2,czc2
      0     0  4.000000000000E+00  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      1     0  2.000000000000E+00  2.000000000000E+00  0.000000000000E+00  0.000000000000E+00

```

The following table explains each parameter as defined by the preceding
line:

| Input Parameter Name | Description | 
|:------------- |:-------------:|:----- |
| nu | Number of poloidal grid points on current surface | 
| nv | Number of toroidal grid points on current surface | 
| nu1 | Number of poloidal grid points on plasma surface | 
| nv1 | Number of toroidal grid points on plasma surface | 
| npol | Number of poloidal Fourier modes in potential solution | 
| ntor | Number of toroidal Fourier modes in potential solution | 
| mf | Number of poloidal Fourier modes for plasma surface | 
| nf | Number of toroidal Fourier modes for plasma surface | 
| md | Number of poloidal Fourier modes for current potential surface (and B normal) | 
| nd | Number of toroidal Fourier modes for current potential surface (and B normal) | 
| np | Number of field periods |
| iota\_edge | Equilibrium rotational transform at edge | 
| phip\_edge | Equilibrium toroidal flux derivative at edge | 
| curpol | Equilibrium total poloidal current in Amps per field period | 
| cut | Normalized toroidal current | 
| cup | Normalized poloidal current | 
| ibex | Include external 1/R field ibex=1| 
| mstrt | Method + svdscan start if \> 1, \>=0 Berr, \<=0 Least square | 
| mstep | Method + svdscan stepsize \<=0 least square, =0 use F04ABE, no svd | 
| mkeep | svd/scan control 0 svdscan, else keep nkeep wgts | 
| mdspw | 2 + exponent of dsur multiplying bfn, ben | 
| curwt | Weight for surface current minimization (only in LSQ branch) | 
| trgwt | Not yet implemented | 
| | For the Output control values -2 means just option 2, +2 means 1 and 2 |
| w\_psurf | Write plasma surface info (1: R/Z, 2: X/Y/Z, 3: NX/NY/NZ, 4: dXdu/dYdu/dXdv/dYdv )| 
| w\_csurf | Write coil surface info (1: R/Z, 2: X/Y/Z, 3:NX/NY/NZ, 4: jx/jy/jz) | 
| w\_bnuv | Write Bnorm field info (2: BN_EXT) | 
| w\_jsurf | Write J surface current info (1: Potential, 2: Current) | 
| w\_xerr | Write X error (displacement) info | 
| w\_svd | Write SVD solution info |

The user should run the [BNORM](BNORM) code to generate a 'nescin' file
with offset winding surface.

------------------------------------------------------------------------

### Execution

To invoke NESCOIL 

```
> xnescoil nescin.example
```

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

To calculate the potential one simply need evaluate:

$$\Phi\left(u,v\right) = \sum_{m=0}^M\sum_{n=-N}^N \Phi_{mn} sin\left(2\pi mu+2\pi nv\right) - \frac{I_{pol}}{N_p}v - I_{tor}u$$

where $$u$$ and $$v$$ are defined on the unit circle from 0 to 1.
When cutting coils one simply should calculate contours of constant
potential. And use the $$u,v$$ paris to evaluate R,phi, and Z. Note
that $$v$$ is define over one field period. 

------------------------------------------------------------------------

### Tutorials

[NCSX NESCOIL Example](NESCOIL NCSX Example.md)
[Solovev NESCOIL Example](NESCOIL Solovev Example.md)

------------------------------------------------------------------------

### References
- [Merkel, P., Solution of stellarator boundary value problems with external currents, Nucl. Fusion 27 5 (1987) 867.](https://iopscience.iop.org/article/10.1088/0029-5515/27/5/018)
- [Landreman, M., Boozer, A.H., Efficient magnetic fields for supporting toroidal plasmas, Phys. Plasmas 23 3 (2016) 032506.](http://aip.scitation.org/doi/10.1063/1.4943201)
