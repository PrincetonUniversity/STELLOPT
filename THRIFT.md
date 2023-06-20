{% include head.html %}

THRIFT
=======

The THRIFT code is a current evolution code for three dimensional
equilibria. It calculates the temporal evolution of the current
profile given a variety of sources.

------------------------------------------------------------------------
Table of Contents
   * [Theory](#theory)
   * [Compilation](#compilation)
   * [Input Data Format](#input-data-format)
   * [Execution](#execution)
   * [Output Data Format](#output-data-format)
   * [Visualization](#visualization)
   * [Tutorials](#tutorials)

------------------------------------------------------------------------
### Theory

TBD

------------------------------------------------------------------------

### Compilation

THRIFT is distributed as part of the STELLOPT package of codes through
Git.

------------------------------------------------------------------------

### Input Data Format

The THRIFT code is controlled through command line inputs and an input
namelists which should be included in the input.ext file. The VMEC
and THRIFT_INPUT namelists are required at a minimum.  Additional
namelists for other codes (BOOTSJ, DIAGNO) should be included in this file.

```fortran
&INDATA
!----- Runtime Parameters -----
  DELT =  1.00000000000000E+00
  NITER = 20000
  NSTEP = 200
  TCON0 =  1.00000000000000E+00
  NS_ARRAY =  16  32  64  128
  FTOL_ARRAY =  1.0E-30  1.0E-30  1.0E-30  1.0E-12
  NITER_ARRAY = 1000 2000 4000 8000 8000 20000
!----- Grid Parameters -----
  LASYM = F
  NFP = 1
  MPOL = 6
  NTOR = 0
  PHIEDGE =  6.28
!----- Free Boundary Parameters -----
  LFREEB = F
!----- Pressure Parameters -----
  GAMMA =  0.00000000000000E+00
  BLOAT =  1.00000000000000E+00
  SPRES_PED =  1.00000000000000E+00
  PRES_SCALE = 1.00000000000000E+00
  PMASS_TYPE = 'power_series'
  AM =  1.00000000000000E+00 -1.00000000000000E+00  0.00000000000000E+00
        0.00000000000000E+00 -1.00000000000000E+00  1.00000000000000E+00
!----- Current/Iota Parameters -----
  CURTOR =  7.000000E+03
  NCURR = 1
  PIOTA_TYPE = 'power_series'
  AI =  0.60000000000000E+00 -0.20000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  PCURR_TYPE = 'power_series'
  AC =  1.00000000000000E+00 -1.00000000000000E+00  0.00000000000000E+00
        0.00000000000000E+00 -1.00000000000000E+00  1.00000000000000E+00
!----- Axis Parameters -----
  RAXIS =  1.00000000000000E+03
  ZAXIS =  0.00000000000000E+00
!----- Boundary Parameters -----
RBC(  0,  0) =  1.0000000000e+03  ZBS(  0,  0) =  0.0000000000e+00
RBC(  0,  1) =  1.0000000000e+00  ZBS(  0,  1) =  1.0000000000e+00
!----- Created by write_vmec 13-Jan-2012 15:12:17 -----
/
&THRIFT_INPUT
  !---------- GENERAL PARAMETERS ------------
  NPARALLEL_RUNS = 1
  BOOTSTRAP_TYPE = 'simple'
  JTOL = 0.01
  NPICARD = 10
  PICARD_FACTOR = 0.75
  LVERBJ = F
  !---------- GRID PARAMETERS ------------
  NRHO = 1024
  NTIMESTEPS = 201
  TSTART = 0
  TEND = 10
  ! BOOZER TRANSFORM
  !---------- BOOZER TRANSFORMATION ------------
  MBOZ = 32
  NBOZ = 16
/
&END
```

Additionaly, the code requries an HDF5 in which the temporal
evoultion of the plasma profiles is specified.

Various examples of how to run the code can be found in
BENCHMARKS subdirectory of the STELLOPT package.

------------------------------------------------------------------------

### Execution

The THRIFT code is controlled through a combination of command-line
inputs and an input namelist. The code requires at least two command
line inputs to be specified `-vmec ext` and `-prof prof.h5` which 
specify the extension of the input file and the HDF5 file containing
the profile information.

    xthrift -vmec <VMEC FILE> -prof <PROFILE FILE> -diagno -noverb -help

| Argument | Default | Description |
|:------------- |:-------------:|:----- |
| -vmec | NONE | VMEC input extension |
| -prof | NONE | Plasma profiles HDF5 file |
| -diagno | NONE | Compute magnetic diagnostic response |
| -noverb | NONE | Suppresses screen output |
| -help | NONE | Print help message. |

In it\'s simplest invokation the code requires
a VMEC input file.

    mpiexec -np 2 /Users/lazerson/bin/xthrift -vmec ATF -prof profiles_ATF_test.h5
    THRIFT Version  0.50
    -----  HDF5 Parameters  -----
       HDF5_version:   1.14 release: 00
    -----  MPI Parameters  -----
       MPI_version:   3.01
       Open MPI v4.1.4, package: Open MPI root@Lazersons-MacBook-Air.local Distribution ident: 4.1.4, repo rev: v4.1.4, May 26, 2022
       Nproc_total:         2
       Nproc_shared:        2
    -----  GIT Repository  -----
       Repository: git@github.com:PrincetonUniversity/STELLOPT.git
       Branch:     THRIFT
       Version:    v2.75-1600-ga77d8
       Built-on:   20.06.2023 08:57:51
       Hash:       a77d822100eaccd25c5419d5b11c16f8b0c45d60
    ----- THRIFT Input Parameters -----
       FILE:             input.ATF
       BOOTSTRAP MODEL:  simple
       NRHO:    1024
       NS:      1024
       NT:       201
       TSTART:   0.0000
       TEND:     2.0000
    ----- Reading Profile File -----
       FILE: profiles_ATF_test.h5
       RHO  = [    0.000,    1.000];  NRHO:  101
       T    = [    0.000,     2.00];    NT:   21
       Ne   = [    0.000,    0.490] E20 m^-3
       Te   = [    0.000,    1.700] keV
       Ni(1)= [    0.000,    0.117] E20 m^-3;  M:   1 amu;  Z:  1
       Ti(1)= [    0.000,    1.700] keV
       Ni(2)= [    0.000,    0.062] E20 m^-3;  M:  12 amu;  Z:  6
       Ti(2)= [    0.000,    1.700] keV
       P    = [    0.000,   18.217] kPa
    -----  MPI Params.   -----
       Parallel runs requested:            1
       Number of Processors:            2
       Shared memory groups:            1
       Processors per group:            2
       Workers per run:            2
       Parallel runs provided:            1
    
    NS =   16 NO. FOURIER MODES =   60 FTOLV =  1.000E-06 NITER =   4000
    PROCESSOR COUNT - RADIAL:    2
    
    ITER    FSQR      FSQZ      FSQL    RAX(v=0)    DELT       WMHD
    
      1  6.88E-01  6.91E-02  2.65E-01  2.000E+00  9.00E-01  3.2453E+01
    148  9.80E-07  8.46E-07  7.88E-08  2.064E+00  9.00E-01  2.9466E+01
    
    NS =   32 NO. FOURIER MODES =   60 FTOLV =  1.000E-08 NITER =   2000
    PROCESSOR COUNT - RADIAL:    2
    
    ITER    FSQR      FSQZ      FSQL    RAX(v=0)    DELT       WMHD
    
      1  2.34E-02  1.56E-02  7.60E-06  2.064E+00  9.00E-01  2.9463E+01
    250  1.12E-07  7.60E-08  6.23E-10  2.090E+00  9.00E-01  2.9462E+01
    330  8.46E-09  8.11E-09  6.58E-11  2.093E+00  9.00E-01  2.9462E+01
    
    NS =   64 NO. FOURIER MODES =   60 FTOLV =  1.000E-10 NITER =   4000
    PROCESSOR COUNT - RADIAL:    2
    
    ITER    FSQR      FSQZ      FSQL    RAX(v=0)    DELT       WMHD
    
      1  5.24E-02  3.02E-02  6.45E-07  2.093E+00  9.00E-01  2.9461E+01
    250  3.31E-07  1.07E-07  5.60E-11  2.094E+00  7.54E-01  2.9461E+01
    500  1.00E-09  3.71E-10  1.16E-12  2.097E+00  7.54E-01  2.9461E+01
    597  9.60E-11  3.58E-11  1.94E-13  2.097E+00  7.54E-01  2.9461E+01
    
    NS =  128 NO. FOURIER MODES =   60 FTOLV =  1.000E-12 NITER =  20000
    PROCESSOR COUNT - RADIAL:    2
    
    ITER    FSQR      FSQZ      FSQL    RAX(v=0)    DELT       WMHD
    
       1  5.45E-02  3.13E-02  1.32E-07  2.097E+00  9.00E-01  2.9461E+01
     250  3.48E-07  1.40E-07  5.11E-12  2.097E+00  5.83E-01  2.9461E+01
     500  8.62E-09  3.37E-09  3.62E-13  2.097E+00  5.83E-01  2.9461E+01
     750  4.86E-10  1.79E-10  3.38E-14  2.098E+00  5.83E-01  2.9461E+01
    1000  1.88E-11  5.55E-12  1.93E-15  2.098E+00  5.83E-01  2.9461E+01
    1224  9.63E-13  2.47E-13  1.46E-16  2.098E+00  5.83E-01  2.9461E+01
    
    EXECUTION TERMINATED NORMALLY
    
    FILE : ATF.001_001
    NUMBER OF JACOBIAN RESETS =    3
    
      TOTAL COMPUTATIONAL TIME (SEC)        15.78
      TIME TO INPUT/OUTPUT                   0.09
         READ IN DATA                        0.00
         WRITE OUT DATA TO WOUT              0.09
      TIME IN FUNCT3D                       15.53
         BCOVAR FIELDS                       1.86
         FOURIER TRANSFORM                   3.66
         INVERSE FOURIER TRANSFORM           3.03
         FORCES AND SYMMETRIZE               1.99
         RESIDUE                             2.83
         EQFORCE                             0.00
    
     T  NSUB  BETA        ITOR     IPLASMA       IBOOT MAX dJ/JOLD
    ===============================================================================
    0.000  1  0.00   0.000E+00   0.000E+00   0.000E+00   0.000E+0
    0.010  1  0.00  -9.982E+00   1.802E+01  -2.800E+01   2.544E+02
    0.010  2  0.00  -9.982E+00   1.802E+01  -2.800E+01   6.628E-04
    0.020  1  0.00  -1.830E+01   3.722E+01  -5.552E+01   5.256E+02
    0.020  2  0.00  -1.830E+01   3.722E+01  -5.552E+01   3.037E-03
    0.030  1  0.01  -2.579E+01   5.682E+01  -8.262E+01   8.007E+02
    0.030  2  0.01  -2.579E+01   5.682E+01  -8.262E+01   4.819E-04
    0.040  1  0.01  -3.282E+01   7.653E+01  -1.094E+02   1.077E+03
    0.040  2  0.01  -3.282E+01   7.653E+01  -1.094E+02   5.049E-04
    0.050  1  0.01  -3.957E+01   9.622E+01  -1.358E+02   1.355E+03
    0.050  2  0.01  -3.957E+01   9.622E+01  -1.358E+02   1.436E-03
    0.060  1  0.01  -4.614E+01   1.158E+02  -1.620E+02   1.633E+03
    0.060  2  0.01  -4.614E+01   1.158E+02  -1.620E+02   4.883E-04
    0.070  1  0.01  -5.262E+01   1.354E+02  -1.880E+02   1.914E+03
    0.070  2  0.01  -5.262E+01   1.354E+02  -1.880E+02   8.017E-04
    0.080  1  0.02  -5.908E+01   1.548E+02  -2.139E+02   2.196E+03
    0.080  2  0.02  -5.908E+01   1.548E+02  -2.139E+02   6.949E-04
    0.090  1  0.02  -6.559E+01   1.741E+02  -2.397E+02   2.482E+03
    0.090  2  0.02  -6.559E+01   1.741E+02  -2.397E+02   5.058E-04
    0.100  1  0.02  -7.220E+01   1.933E+02  -2.655E+02   2.772E+03
    0.100  2  0.02  -7.220E+01   1.933E+02  -2.655E+02   1.489E-03
    0.110  1  0.02  -7.894E+01   2.123E+02  -2.912E+02   3.065E+03
    0.110  2  0.02  -7.894E+01   2.123E+02  -2.912E+02   1.511E-03

The BENCHMARKS directory contains a set of tests for THRIFT.  The 
input files are located in the THRIFT_TEST subdirectory.  They can
all be invoked by calling make thrift_test from the BENCHMARKS 
directory. The comparrision scripts require Python.

------------------------------------------------------------------------

### Output Data Format

The data from each run is output into a HDF5 file where all datasets
located at the root level.  The following table helps to define the
variables (all values in mks units, angles in radians)

| Name | Type | Size | Description |
| :--- | :--- | :---: | :--- |
| VERSION | DOUBLE | 1 | THRIFT version information |
| lvmec | BOOLEAN | 1 | VMEC Equilibrium |
| leccd | BOOLEAN | 1 | Electron Cyclotron Current Drive included |
| lnbcd | BOOLEAN | 1 | Neutral Beam Current Drive included |
| lohmic | BOOLEAN | 1 | Ohmic current drive included |
| ntimesteps | INTEGER | 1 | Number of timesteps |
| nssize | INTEGER | 1 | Number of radial grid points (s) |
| nrho | INTEGER | 1 | Number of radial grid points (r/a) |
| npicard | INTEGER | 1 | Maximum number of Picard Iterations |
| jtol | DOUBLE | 1 | Current profile convergence factor % |
| picard_factor | DOUBLE | 1 | Picard iteration factor |
| THRIFT_RHO | DOUBLE | nrho | THRIFT radial grid (r/a) |
| THRIFT_S | DOUBLE | nssize | THRIFT radial grid (s) |
| THRIFT_T | DOUBLE | ntimesteps | THRIFT temporal grid (sec) |
| THRIFT_PHIEDGE | DOUBLE | ntimesteps | Enclosded toroidal flux (Wb) |
| THRIFT_J | DOUBLE | nssize,ntimesteps | THRIFT total current density A/m^2 |
| THRIFT_JBOOT | DOUBLE | nssize,ntimesteps | THRIFT bootstrap current density A/m^2 |
| THRIFT_JECCD | DOUBLE | nssize,ntimesteps | THRIFT ECCD current density A/m^2 |
| THRIFT_JNBCD | DOUBLE | nssize,ntimesteps | THRIFT NBCD current density A/m^2 |
| THRIFT_JOHMIC | DOUBLE | nssize,ntimesteps | THRIFT Ohmic current density A/m^2 |
| THRIFT_JPLASMA | DOUBLE | nssize,ntimesteps | THRIFT plasma response current density A/m^2 |
| THRIFT_JSOURCE | DOUBLE | nssize,ntimesteps | THRIFT total source current density A/m^2 |
| THRIFT_I | DOUBLE | nssize,ntimesteps | THRIFT total current A|
| THRIFT_IBOOT | DOUBLE | nssize,ntimesteps | THRIFT bootstrap current A |
| THRIFT_IECCD | DOUBLE | nssize,ntimesteps | THRIFT ECCD current A |
| THRIFT_INBCD | DOUBLE | nssize,ntimesteps | THRIFT NBCD current A |
| THRIFT_IOHMIC | DOUBLE | nssize,ntimesteps | THRIFT Ohmic current A |
| THRIFT_IPLASMA | DOUBLE | nssize,ntimesteps | THRIFT plasma response current A |
| THRIFT_ISOURCE | DOUBLE | nssize,ntimesteps | THRIFT total source current A |
| THRIFT_ETAPARA | DOUBLE | nssize,ntimesteps | Parallel electric resistivity |
| THRIFT_P | DOUBLE | nssize,ntimesteps | Pressure profile Pa |
| THRIFT_PPRIME | DOUBLE | nssize,ntimesteps | Pressure profile Pa derivative |

------------------------------------------------------------------------

### Visualization

The data from the HDF5 file may be easily plotted in many plotting
packages.

------------------------------------------------------------------------

### Tutorials

TBD

------------------------------------------------------------------------

### References

-   [van Ham, L. \"THRIFT THESIS\" TU/e (2023)]()
-   [Strand, P.I. and Houlberg, W.A. \"Magnetic flux evolution in highly shpaed plasmas\" Physics of Plasmas 8, 2782 (2001)](http://dx.doi.org/10.1063/1.1366618)
-   [Strand, P.I. and Houlberg, W.A. \"Solution techniques for magnetic flux evolution in toroidal plasmas.\" Fusion Techno. 39, 1091-1095](https://doi.org/10.13182/FST01-A11963389)



