CAS3D
=====

[toc](toc) The Code for the Analysis of the Stability of 3-D
equilibria\<ref\>[Schwab, C. \"Ideal magnetohydrodynamics: Global mode analysis of three-dimensional plasma configurations.\" Phys. Fluids B, 5(9): 3195, 1993.](@http://dx.doi.org/10.1063/1.860656)\</ref\>
(CAS3D) calculates the global stability properties of 3D equilibria.

------------------------------------------------------------------------

### Theory\[\[\#Theory\]\]

Here\'s the section to explain the theory behind the code.

------------------------------------------------------------------------

### Compilation\[\[\#Compilation\]\]

The CAS3D code is distributed as a suite of three codes (MC3D, PRERUN,
CAS3D). The codes build using
[AUTOCONF](@http://www.gnu.org/software/autoconf/autoconf.html) and the
[Make](@https://www.gnu.org/software/make/) build system. The
\'configure.ac\' file supplied with the code should be the only file one
needs to edit to get their system to run. One should check if their
machine type is already supported by running \'uname -n\' on their
system. One can then check the file to see if their build system is
supported. If not they need to add a section to the \'configure.ac\'
file which specifies the following:
[code format=\"autoconf\"](code format="autoconf") \[MACHINE\],
\[COMPILER=\"GNU\" FC=\"gfortran4\" ! Path to Fortran compiler
FC\_MPI=\"mpif90\" ! Path to MPI Fortran compiler AS\_IF(\[test
\"\$DEBUG\" = \"yes\" \], \[FCFLAGS=\" -C -g -I\\\$(OBJ\_DIR)
-J\\\$(OBJ\_DIR) -I\\\$(NETCDF\_DIR)/include \"\], ! Debug compile flags
\[FCFLAGS=\" -O2 -fdefault-real-8 -I\\\$(OBJ\_DIR) -J\\\$(OBJ\_DIR)
-I\\\$(NETCDF\_DIR)/include\"\]) ! Non-debug compile flags
FCFLAGS\_FIX=\$FCFLAGS ! Fixed format flags FCFLAGS=\$FCFLAGS\"
-ffree-form -ffree-line-length-none -ffixed-line-length-none\" ! Free
format flags FPPFLAGS\_FIX= FPPFLAGS= DPREF=-D ! Precompiler flag
(almost always -D) BITS= LD=\$FC ! Path to linker
LDFLAGS=\"-L\\\$(BLAS\_DIR)/lib -lblas -L\\\$(NAG\_DIR)/lib -lnag
-L\\\$(NETCDF\_DIR)/lib -lnetcdf -lnetcdff\" ! Linker flag, put
libraries here PLPLOT\_DIR=\"\\\$(PLPLOT\_DIR) \\\$(PLPLOT\_DIR)/lib\"
!PLPLOT directories (leave blank if PLPLOT not supported)
PLPLOT=\"-lplplotd \" !PLPLOT libraries (leave blank if PLPLOT not
supported) LIB\_DIR= LIB= IRAM=\"no\" ! yes: ARPACK included, otherwise
no \], [code](code) One the \'configure.ac\' file has been properly
setup, the user may issue the following commands:

[code](code) \> autoconf \> ./configure Making a RELEASE version
(default) Making a SERIAL version (default) Doing an ORDINARY make
(default) Using STELLARATOR symmetry (default) Using the KINETIC energy
as normalization (default) For EVEN-parity perturbation functions
(default) NOT using the PHASE-FACTOR transform (default)
sunfire10.pppl.gov sunfire configure: PLPLOT not set: preparing
executables without plotting configure: creating ./config.status
config.status: creating mercier/source/makefile config.status: creating
mercier/source/src/os.mk config.status: creating mc3d/source/makefile
config.status: creating mc3d/source/src/os.mk config.status: creating
cas3d/source/makefile config.status: creating cas3d/source/src/os.mk
config.status: creating ./makefile \> make [code](code) Alternatively
you may pass options to the configure command to control what type of
executable is built (type \'./configure -h\'):

[code](code) Optional Features: \--disable-option-checking ignore
unrecognized \--enable/\--with options \--disable-FEATURE do not include
FEATURE (same as \--enable-FEATURE=no) \--enable-FEATURE\[=ARG\] include
FEATURE \[ARG=yes\] \--enable-debug: Making a DEBUG version default:
release version \--enable-clean: CLEAN make (remakes existing .o files,
even if they are not out of date) default: ordinary make
\--enable-general: Using GENERAL toroidal symmetry (up-down asymmetric
tokamak) default: stellarator symmetry(up-down symm. tokamak)
\--enable-parallel: Making a PARALLEL version default: serial version
\--enable-odd: Using ODD-parity perturbation functions default:
even-parity \--enable-asym: Using GENERAL-parity perturbation functions
default: even-parity

Optional Packages: \--with-PACKAGE\[=ARG\] use PACKAGE \[ARG=yes\]
\--without-PACKAGE do not use PACKAGE (same as \--with-PACKAGE=no)
\--with-phasefactor: Using the PHASE-FACTOR transform default: no
phase-factor transform \--with-b1norm: Using the PERTURBED MAGNETIC
ENERGY as normalization default: kinetic energy \--with-c1boundary:
Using a BOUNDARY integral as normalization default: kinetic energy

Report bugs to the package provider. [code](code)

------------------------------------------------------------------------

### Input Data Format\[\[\#Input\]\]

To perform a stability analysis the data from a VMEC equilibrium must be
converted into straight field line coordinates by MC3D, become adjusted
by PRERUN, and then CAS3D may be run. The result is three separate data
input mechanisms.

#### MC3D

The MC3D code takes a fortran namelist (MC3D\_INDATA) as input, and a
file entitled \'fort.77\' which contains the name of the plotting files
if desired. For convenience this namelist can reside in the VMEC input
file. The namelist has the following format: [code](code) &mc3d\_indata
icode=2007 ! Type of input equilibrium (see table below)
woutfile=\"wout\_test.nc\" ! Name of VMEC output file icode=2007 ! Type
of wout file nsurf=19 ! Number of surfaces (for plotting) indexsurf = 5
10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 ! Array of
surfaces for plotting m0 =64 ! Maximum equilibrium poloidal mode number
(reset by wout file) n0 =32 ! Maximum poloidal mode number nu=258 !
Poloidal real space points nv=130 ! Toroidal real space points m0b = 64
! Fourier truncation of poloidal straight field line spectrum n0b = 32 !
Fourier truncation of toroidal straight field line spectrum conv=0.001 /

[code](code) The value of ICODE tells the code how to read the VMEC
\'wout\' file. The following table explains the possible values: \|\|
ICODE \|\| Meaning \|\| \|\| 91 \|\| VMEC form 1991 \|\| \|\| 98 \|\|
NEMEC Geiger (ASCII) \|\| \|\| 99 \|\| NEMEC99 (ASCII) \|\| \|\| 2000
\|\| VMEC2000 (ASCII) \|\| \|\| 2001 \|\| NEMEC2000 (ASCII fort.8) \|\|
\|\| 2002 \|\| NEMEC ERS (ASCII) \|\| \|\| 2003 \|\| VMEC2003 (ASCII)
\|\| \|\| 2005 \|\| TERPSICHORE (ASCII fort.18) \|\| \|\| 2007 \|\|
VMEC2000 (netCDF) \|\| \|\| 2009 \|\| VMEC2000 (ASCII) \|\|

#### PRERUN

The PRERUN code takes the output from MC3D (the magnetic coordinates)
and allows the user to adjust the radial grid while also calculating a
few new metric coefficients for CAS3D. The PRERUN code takes a name list
(PRERUN\_INDATA) and the \'for\_cas3d\_stp\_14.dat\' file as input. The
namelist has the format: [code](code) &PRERUN\_INDATA
PERFECT\_EQUILIBRIUM = .TRUE. NEW\_FORMAT = .TRUE. BACKTRANSFORM =
.TRUE. ILEFT = 0 IRIGHT = 128 NS = 128 NU = 258 NV = 130 NPER = 3 NTOPOL
= 3 binary\_input = \'\' / [code](code) CAS3D The CAS3D code takes the
\'for\_cas3d\_rzuv.dat\' from MC3D, the \'fort.12\' file from PRERUN,
and an \'input2\' file which controls execution. The input2 file has the
following format: [code](code) ncsx gamma phase-factor factor on
back-transform use of D3D-data use of HMSC-data plot of b1n gamma (M,N)
c2 back d3d hmsc for\_b1n / xis 1.6666 0 0 1. .true. .false. .false.
.false. matrix(yes=0) evp(1=no,3=shift) shift shiftmin shiftmax number 0
3 -7.0e-03 2.000e-02 0.00e-02 1 number density(axis) \[10\^20 m\^-3\]
partition of shifts: 1:linear 2:arsinh-partition 0.34 3 shape function:
0:s**(m/2)\*(1-s), 1:s**(m/2), 2:localized support,3:j 4: jumps 5: m=2
xis axis IRAND (IAA IBB IALP) \|\| icutl icutr (radial interval) 2 -7
3000 3 0 48 eigensolver (IVIT=inverse iteration or IRAM=implicitly
restarted Arnoldi) eigensolver (upper case letters only) IVIT

1.  svd-solution\|normal component of B1 on plasma boundary isvd
    value\_b1 (fraction of \|B\|(axis)) find mu tableau (only if
    gamma=0) 0 3.25e-05 .false. boundary condition: 0: free boundary, 1:
    fixed boundary, 2: free-without IFIX Fourier in integral equation
    M\_pot N\_pot POINT SPACE (field period) NU\_vac NV\_vac 0 25 15 75
    30 for the vacuum part: representation of conducting wall
    with\_wall: 0=no, 1=yes mnwall: number of harmonics (=lines to skip
    if with\_wall=0) 0 9 ellis008:3048 ls all\_profiles.dat
    for\_cas3d\_rzuv\_18.dat fort.77 physical\_profiles.dat
    xis\_harmonics.dat xis\_largest31-40.dat apar\_harmonics\_efield.dat
    fort.12 input2 plot\_fil\_0.meta xis\_largest01-10.dat
    b1\_boundary\_plasma\_even.dat fort.14 mercier.dat scan\_0.dat
    xis\_largest11-20.dat b1\_boundary\_vacuum\_even.dat fort.16
    phi\_harmonics\_efield.dat wkin\_0.dat xis\_largest21-30.dat
    ellis008:3049 cat input2 ncsx gamma phase-factor factor on
    back-transform use of D3D-data use of HMSC-data plot of b1n gamma
    (M,N) c2 back d3d hmsc for\_b1n / xis 1.6666 0 0 1. .true. .false.
    .false. .false. matrix(yes=0) evp(1=no,3=shift) shift shiftmin
    shiftmax number 0 3 -7.0e-03 2.000e-02 0.00e-02 1 number
    density(axis) \[10\^20 m\^-3\] partition of shifts: 1:linear
    2:arsinh-partition 0.34 3 shape function: 0:s**(m/2)\*(1-s),
    1:s**(m/2), 2:localized support,3:j 4: jumps 5: m=2 xis axis IRAND
    (IAA IBB IALP) \|\| icutl icutr (radial interval) 2 -7 3000 3 0 48
    eigensolver (IVIT=inverse iteration or IRAM=implicitly restarted
    Arnoldi) eigensolver (upper case letters only) IVIT
2.  svd-solution\|normal component of B1 on plasma boundary isvd
    value\_b1 (fraction of \|B\|(axis)) find mu tableau (only if
    gamma=0) 0 3.25e-05 .false. boundary condition: 0: free boundary, 1:
    fixed boundary, 2: free-without IFIX Fourier in integral equation
    M\_pot N\_pot POINT SPACE (field period) NU\_vac NV\_vac 0 25 15 75
    30 for the vacuum part: representation of conducting wall
    with\_wall: 0=no, 1=yes mnwall: number of harmonics (=lines to skip
    if with\_wall=0) 0 9 Maximum \# of inverse iterations \# of
    independent eigenfunctions to be saved ITNUM num\_eigen 1000 0
    printing of matrix elements: yes=1, no=0 IPRINT isprint 0 24
    relative distance of independent eigenvalues selection: plot only
    mdom relative\_distance MDOM coefficients 1.e-3 -1 s-interval for
    plots sa sb 0. 1. option for plots: 0=n0 1=yes form coef energy xis
    0 1 1 1 use grid perturbed or unperturbed field on\_v\_cut
    true=grid,false=shift true=unperturbed, false=perturbed
    true=on\_v\_cut, false=on\_phi\_cut .false. .false. .true. gourdon
    grid nf nfd nr nz r00 z00 dr0 dz0 itormin itormax ftol\_resi
    dist\_curve 240 240 100 128 5.5 0. 1.5 1.5 1 240 0.3 0.00001
    perturbation spectrum:isearch=0 (c2,islow=0!) =1 (sqrtg a) isearch 0
    for plot: factor on c3 plotting threshold(xis,eta,mu plots) 1. 0.01
    location of continua for selected partial modes number label
    curvature threshold 0 20. mass density and electron temperature
    profiles (=0: unity, \>0 number of terms in polynomial
    a\_0+a\_1\*s+a\_2\*s\*s+\... ) idens coeffs read\_temperature\_data
    (false: use p0=n0\*T0 .true.: insert polynomial data after rho
    polynomial) 0 .false. aspect ratio 11. c ntopol nper 3 3 c number of
    Fourier modes c xis eta mu nstop 100 100 100 0 c number of resonant
    modes c xis eta mu 0 0 0 c (m,n) index pairs c xis 0 20 -11 1 0

0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20

-11 0 0 0 0 0 0 0 0 0 0 0 82 83 84 85 86 87 88 89 90 91 -10 0 0 0 0 0 0
0 0 0 0 72 73 74 75 76 77 78 79 80 81 92 -9 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 -8 0 0 0 0 0 0 0 0 0 62 63 64 65 66 67 68 69 70 71 93 100
-7 0 0 0 0 0 0 0 52 53 54 55 56 57 58 59 60 61 94 95 0 0 -6 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -5 0 0 0 0 0 31 32 33 34 35 36 37 38 39 40
96 99 0 0 0 0 -4 0 0 0 21 22 23 24 25 26 27 28 29 30 97 0 0 0 0 0 0 0 -3
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 0 11 12 13 14 15 16 17 18
19 20 98 0 0 0 0 0 0 0 0 0 -1 51 1 2 3 4 5 6 7 8 9 10 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 41 42 43 44 45 46 47
48 49 50 0 0 0 0 0 0 0 0 0 0 c eta 0 20 -11 1 0

0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20

-11 0 0 0 0 0 0 0 0 0 0 0 82 83 84 85 86 87 88 89 90 91 -10 0 0 0 0 0 0
0 0 0 0 72 73 74 75 76 77 78 79 80 81 92 -9 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 -8 0 0 0 0 0 0 0 0 0 62 63 64 65 66 67 68 69 70 71 93 100
-7 0 0 0 0 0 0 0 52 53 54 55 56 57 58 59 60 61 94 95 0 0 -6 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -5 0 0 0 0 0 31 32 33 34 35 36 37 38 39 40
96 99 0 0 0 0 -4 0 0 0 21 22 23 24 25 26 27 28 29 30 97 0 0 0 0 0 0 0 -3
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 0 11 12 13 14 15 16 17 18
19 20 98 0 0 0 0 0 0 0 0 0 -1 51 1 2 3 4 5 6 7 8 9 10 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 41 42 43 44 45 46 47
48 49 50 0 0 0 0 0 0 0 0 0 0 c mu 0 20 -11 1 0

0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20

-11 0 0 0 0 0 0 0 0 0 0 0 82 83 84 85 86 87 88 89 90 91 -10 0 0 0 0 0 0
0 0 0 0 72 73 74 75 76 77 78 79 80 81 92 -9 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 -8 0 0 0 0 0 0 0 0 0 62 63 64 65 66 67 68 69 70 71 93 100
-7 0 0 0 0 0 0 0 52 53 54 55 56 57 58 59 60 61 94 95 0 0 -6 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -5 0 0 0 0 0 31 32 33 34 35 36 37 38 39 40
96 99 0 0 0 0 -4 0 0 0 21 22 23 24 25 26 27 28 29 30 97 0 0 0 0 0 0 0 -3
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 0 11 12 13 14 15 16 17 18
19 20 98 0 0 0 0 0 0 0 0 0 -1 51 1 2 3 4 5 6 7 8 9 10 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 41 42 43 44 45 46 47
48 49 50 0 0 0 0 0 0 0 0 0 0 m n RBC RBS ZBC ZBS 0 0 1.5587E+00
0.0000E+00 -7.4282E-03 0.0000E+00 1 0 9.1202E-01 1.2839E-02 1.2839E-02
1.6251E+00 2 0 1.2956E-01 7.5692E-03 4.7646E-03 -1.0589E-01 3 0
4.3674E-02 -2.5566E-03 -1.7452E-03 -1.8561E-02 4 0 -3.1006E-03
-5.4278E-04 -1.8013E-03 1.0895E-02 5 0 -1.6836E-04 1.2065E-03 6.6790E-04
4.1134E-03 6 0 3.1740E-03 -1.3182E-04 3.1708E-04 -3.0336E-03 7 0
-1.4565E-04 -4.4772E-04 -6.8330E-04 -3.7720E-04 8 0 -6.2512E-04
7.2924E-05 -1.6238E-04 1.6542E-03 c resonant b1n harmonics (on plasma
boundary) relative to average modB on axis c number of resonant
harmonics (up to 10) 0 c m n b1n\_res\_boundary / \<\|B0\|\>\_axis
\[default value: 0.0001\] c plot of the B1 normal component on the
resonant surface: true: phi E \[0, 1\] else \[-1/2, 1/2\] .false. c
resolution parameters for surface current m\_cur\_low m\_cur n\_cur
nu\_cur nv\_cur threshold\_svd 0 30 0 210 4 0. ! RHS delta1W
pressure\_change geometry\_change external error-field .false. .false.
.false. ! perturbed pressure: use polynomial (if false: pointwise)
.true. ! modified pressure profile (pressure in Pascal) ! ipres 1 !
data: (i, a\_i) OR (s\_i, p\_i) 0 0.

[code](code)

------------------------------------------------------------------------

### Execution\[\[\#Execution\]\]

Explain how to execute the code and what it produces.

------------------------------------------------------------------------

### Output Data Format\[\[\#Output\]\]

Explain how the output data is formatted.

------------------------------------------------------------------------

### Visualization\[\[\#Visualization\]\]

Explain how to visualize the data.

------------------------------------------------------------------------

### Tutorials\[\[\#Tutorials\]\]

Put links to tutorial pages here.
