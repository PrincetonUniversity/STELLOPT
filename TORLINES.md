TORLINES
========

The TORLINES code is a field line following code which follows magnetic
field lines on a grid concentric with a VMEC equilibrium. The code has
facilities to allow field lines to be followed or EMC3 grid data output.

------------------------------------------------------------------------

### Theory

![](images/torlines_ncsx_grid.jpg)

The TORLINES code follows field lines on a background grid concentric
with a given equilibrium ([VMEC](VMEC)). In the domain of the
equilibrium the magnetic field is supplied on a radial, poloidal and
toroidal grid. Outside the equilibrium domain a concentric grid is
constructed using similar routines to those found in the [BNORM](BNORM)
code. This allows the field lines to be followed on a more natural grid.
The magnetic field is represented in terms of radial, poloidal and
toroidal components.

The field line following is achieved using a common simplification
allowing the elimination of time in favor of toroidal angle. The ODE\'s
the code solves are [math](math) \\frac{\\partial \\rho}{\\partial
\\phi} = R\\frac{B\_\\rho}{B\_\\phi} [math](math) and [math](math)
\\frac{\\partial \\theta}{\\partial \\phi} =
R\\frac{B\_\\theta}{B\_\\phi} [math](math) The resulting ODE is solved
with a user determined step-size and accuracy. The available ODE
packages are:
[LSODE](@https://computation.llnl.gov/casc/odepack/odepack_home.html),
and
[NAG D02CJF](@http://www.nag.co.uk/numeric/fl/manual/pdf/D02/d02cjf.pdf)(requires
license). The calculation of field line trajectories is parallelized
over each field line. Thus each processor can follow each field line
independently to speed computation.

The code assembles fields from various sources. The vacuum component of
the field can be calculated directly from a coils file or an mgrid file.
In the later case the TORLINES grid must fall within the domain of the
mgrid file. The plasma field inside the equilibria domain is placed
directly on the background grid. The plasma response external to the
equilibria is calculated using a virtual casing principle. In order to
maintain accuracy near the surface of the plasma, the virtual casing
principle employs an adaptive integration scheme over the surface
current. This scheme is either handled by
[NAG D01EAF](http://www.nag.co.uk/numeric/fl/manual/pdf/D01/d01eaf.pdf)
(if available), or the
[DCUHRE](@http://www.math.wsu.edu/faculty/genz/software/software.html)
algorithm\<ref\>[J. Berntsen, T. O. Espelid and A. Genz, Trans. Math. Softw. 17 (1991), pp. 437-451.](@http://dx.doi.org/10.1145/210232.210233)\</ref\>
. The grids are treated in vector fashion allowing an fast scaleable
parallel execution over the vacuum region.

------------------------------------------------------------------------

### Compilation

TORLINES is a component of the STELLOPT suite of codes. Compilation of
the STELLOPT suite is discussed on the
[STELLOPT Compilation Page](vmecwiki/STELLOPT Compilation). To obtain
the code please contact the author Samuel A. Lazerson
(<lazerson@pppl.gov>).

------------------------------------------------------------------------

### Input Data Format

The TORLINES code is controlled through command line inputs and an input
namelist which should be placed in the input.ext file. While the entire
VMEC input name is not required the EXTCUR array is required. The name
lists should look like:

    #!fortran
    &INDATA
     ! VMEC input namelist (only need coil currents for usual runs)
     EXTCUR(1) = 10000.00
     EXTCUR(2) = 10000.00
     EXTCUR(3) = 12000.00
     EXTCUR(4) = 12000.00
     EXTCUR(5) = 6000.00
    /
    &TORLINES_INPUT
     K = 64      ! Number of radial gridpoints to use
     NU = 64   ! Number of poloidal gridpoints
     NV = 64   ! Number of toroidal gridpoints (per period)
     NPOINC = 60  ! How many times per field period to output a fieldline (for EMC3 stuff this is irrelevant)
     NU_VC = 128  ! Number of poloidal points for virtual casing (minimum of 4*(MPOL-1) )
     NV_VC = 128  ! Number of toroidal points for virtual casing (minimum of 4*NTOR)
     LVC_FIELD = T ! Use virtual casing (T: Virtual Casing, F: Volume Integration)
     VC_ADAPT_TOL = 1.0E-3 ! Tollerance for virtual casing
     FOLLOW_TOL = 1.0E-9 ! Field line following tollerance
     INT_TYPE = 'LSODE'  ! Integrator type (LSODE,NAG)
     BOUND_SEPARATION = 1.3 ! How much to scale up the boundary
     NZONE_EMC3 = 12 ! Number of toroidal zones (-emc3 option)
     NPOLO_EMC3 = 361 ! Number of poloidal grid points (-emc3 option)
     NTORO_EMC3 = 11 ! Number of toroidal grid points (-emc3 option)
     S0_EMC3 = 0.35 ! Inner plasma boundary for EMC3 (-emc3 option, in TORLINES coordinates)
     S1_EMC3 = 0.85 ! Outer plasma boundary for EMC3 (-emc3 option, in TORLINES coordinates)
    /
    &END

------------------------------------------------------------------------

### Execution

The TORLINES code is controlled through a combination of command-line
inputs and an input namelist. The input namelist must be in the
equilibrium input file. That file must also contain the VMEC INDATA
namelist (although only the EXTCUR array will be used). The TORLINES
code is run from the command line taking an equilibrium input file as a
necessary argument. This input file must have the INDATA (for EXTCUR)
and TORLINES\_INPUT namelists in it.

    XTORLINES -vmec <VMEC FILE> -coil <COIL FILE> -mgrid <MGRID FILE> -vessel <VESSEL FILE> -emc3 -auto -raw -noverb -help

\|\| Argument \|\| Default \|\| Description \|\| \|\| -vmec \|\| NONE
\|\| VMEC input extension \|\| \|\| -coil \|\| NONE \|\| Coils File \|\|
\|\| -mgrid \|\| NONE \|\| Makegrid style vacuum grid file \|\| \|\|
-vessel \|\| NONE \|\| First wall file \|\| \|\| -auto \|\| NONE \|\|
Auto calculate axis and edge maximum resolution \|\| \|\| -noverb \|\|
NONE \|\| Suppresses screen output \|\| \|\| -emc3 \|\| NONE \|\|
Outputs the B-Field on the cylindrical grid only. \|\| \|\| -raw \|\|
NONE \|\| Treats EXTCUR array as raw values (EXTCUR is a scale factor
applied to what\'s in the coils file). \|\| \|\| -help \|\| NONE \|\|
Print help message \|\| In it\'s simplest invokation the code requires a
VMEC input and wout file.

    >mpirun -np 64 ~/bin/xtorlines -vmec ncsx_c09r00_free -coil coils.c09r00 -auto
    TORLINES Version 1.00
     -----TORLINES File Parameters-----
              k:  70
             nu:  90   nv:  60
      bound_sep:    1.200
     -----VMEC File Parameters-----
        file: ncsx_c09r00_free
           m:  11   nu:  90
           n:   6   nv:  60
       mnmax:  137
         nfp:   3
          ns:  99
       Total Current: -178.653 [kA]
    ----- Vacuum Grid Info. -----
           Eq. Surface:     59
           Grid Points:  59400
    ----- COILS Information -----
       FILE: coils.c09r00
       Coil Periodicity:   3
       Current Systems:  10
       Current Type:      SCALED
       Num Coils  =    6  EXTCUR =  652.272 [kA]
       Num Coils  =    6  EXTCUR =  651.869 [kA]
       Num Coils  =    6  EXTCUR =  537.744 [kA]
       Num Coils  =    8  EXTCUR =    0.000 [A]
       Num Coils  =    8  EXTCUR =    0.000 [A]
       Num Coils  =    8  EXTCUR =   28.095 [kA]
       Num Coils  =   12  EXTCUR =  -54.805 [kA]
       Num Coils  =    4  EXTCUR =   30.123 [kA]
       Num Coils  =    2  EXTCUR =   94.241 [kA]
       Num Coils  =   18  EXTCUR =   45.514 [kA]
         Vacuum Field Calculation [ 99]%
    ----- Virtual Casing Information -----
       INTEGRAL TYPE: Surface Current
       MIN_GRID_DISTANCE =  6.0493E-02
       NORMAL_AREA =  2.4595E+01
       NR =    1;   NU =  128;  NV =  128;  NFP =   3
       NUVP =  49152
       ABS_TOL =  0.0000E+00;   REL_TOL =  5.0000E-03
       MIN_CLS =      0   (16777216)
         Plasma Field Calculation [ 99]%
    ---------- EXECUTION ----------
    ----     Rough GRID     ----
       rho   = [ 0.00000, 1.00000];
      ----- FOLLOWING FIELD LINES -----
          Method: NAG
           Lines:    201
           Steps:  180000   Delta-phi: 0.3491E-01
             Tol: 0.1000E-08  Type: M
         Fieldline Calculation [ 75]%
    ----     Edge GRID     ----
       rho   = [ 0.78824, 1.00000];
      ----- FOLLOWING FIELD LINES -----
          Method: NAG
           Lines:    201
           Steps:  180000   Delta-phi: 0.3491E-01
             Tol: 0.1000E-08  Type: M
         Fieldline Calculation [ 75]%
    ----     Final GRID     ----
       rho   = [ 0.00000, 1.00000];
      ----- FOLLOWING FIELD LINES -----
          Method: NAG
           Lines:    201
           Steps:  180000   Delta-phi: 0.3491E-01
             Tol: 0.1000E-08  Type: M
         Fieldline Calculation [ 83]%
    ---------- TORLINES DONE ----------

------------------------------------------------------------------------

### Output Data Format

Data is output in HDF5 format (where available). In the torlines\_ext.h5
file the background grid, magnetic field, and field line trajectories
can be found. Field line trajectories are output in terms of R,phi,Z
coordinates for ease of plotting in real space.

------------------------------------------------------------------------

### Visualization

The data can be visualized in many ways. The Poincare points are stored
as arrays of each trajectory and NPOINC cross sections per field period.

------------------------------------------------------------------------

### Tutorials

[NCSX TORLINES Example](NCSX TORLINES Example)