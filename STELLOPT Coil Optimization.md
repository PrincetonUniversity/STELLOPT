Tutorial: STELLOPT Coil Optimization Tutorial
==================================================

The [STELLOPT](STELLOPT) code can be used to optimize coil shapes using
a spline representation. The coils are defined by a set of cubic spline
knots which are defined in terms of rho, theta, and zeta. The rho value
is a physical distance from the equilibrium surface. The theta value is
the [VMEC](VMEC) poloidal angle. The zeta value is defined as running
from zero to 2 pi over a field period. These values are defined in the
OPIMUM namelist. [VMEC](VMEC) and [BNORM](BNORM) are run automatically 
and an option exists to produce [FIELDLINES](FIELDLINES) for each
minimum found.

------------------------------------------------------------------------

# Input namelists

When running a coil optimization run the code will run VMEC and BNORM
once in the beginning and reference those values for all subsequent
coil iterations. The only relevant parameter in the `&INDATA` namelist
is the `EXTCUR` array as that defines the coil currents. 

```fortran
&INDATA
  LFREEB = F ! For now free boundary coil optimization is not supported
  EXTCUR = 5*-1.0E6 ! This is used define the coil currents
  ! The rest of INDATA should be set to your target equilibrium.
/
&OPTIMUM
!-----------------------------------------------------------------------
!          OPTIMIZER RUN CONTROL PARAMETERS
!-----------------------------------------------------------------------
  NFUNC_MAX = 5
  EQUIL_TYPE = 'VMEC2000_ONEEQ'   ! So we run VMEC only once
  OPT_TYPE   = 'GADE'             ! You can run any optimizer
  FACTOR      =    0.5
  EPSFCN      =    0.3
  MODE        =  0002
  CR_STRATEGY =  0000
  LKEEP_MINS  = T
  NPOPULATION =  10
  NOPTIMIZERS = 5
!-----------------------------------------------------------------------
!          OPTIMIZED QUANTITIES
!-----------------------------------------------------------------------
  ! These options tell the code to vary all 3 knots
  LCOIL_KTS_OPT(1,1:8) = 8*T
  LCOIL_KTS_OPT(2,1:8) = 8*T
  LCOIL_KTS_OPT(3,1:8) = 8*T
  LCOIL_KTS_OPT(4,1:8) = 8*T
  LCOIL_KTS_OPT(5,1:8) = 8*T
!-----------------------------------------------------------------------
!          MULTI-FILAMENT OPTIONS (comment out for single filament)
!-----------------------------------------------------------------------
  ! This sets the number of windings in the binormal directions of the coil
  NW_COIL = 3
  NH_COIL = 3
  ! This sets the width of the coil windings
  WIDTH_COIL = 0.1
  HEIGHT_COIL = 0.1
!----------------------------------------------------------------------
!       Coil Spline Knots
!       RHO is distance from plasma surface in m.
!       THETA is poloidal angle in radians.
!       Zeta is field period toroidal angle in radians.
!       Note that theta will be closed by the code. Modular coils assumed
!----------------------------------------------------------------------
!----- COIL  1
  RHO_COIL_KTS(  1,:) =    3.000000000000E-001   3.000000000000E-001   3.000000000000E-001   3.000000000000E-001   3.000000000000E-001   3.000000000000E-001   3.000000000000E-001   3.000000000000E-001
  THETA_COIL_KTS(  1,:) =    0.000000000000E+000   7.853981600000E-001   1.570796330000E+000   2.356194490000E+000   3.141592650000E+000   3.926990820000E+000   4.712388980000E+000   5.497787140000E+000
  ZETA_COIL_KTS(  1,:) =    3.141592700000E-001   3.141592700000E-001   3.141592700000E-001   3.141592700000E-001   3.141592700000E-001   3.141592700000E-001   3.141592700000E-001   3.141592700000E-001
!----- COIL  2
  RHO_COIL_KTS(  2,:) =    3.000000000000E-001   3.000000000000E-001   3.000000000000E-001   3.000000000000E-001   3.000000000000E-001   3.000000000000E-001   3.000000000000E-001   3.000000000000E-001
  THETA_COIL_KTS(  2,:) =    0.000000000000E+000   7.853981600000E-001   1.570796330000E+000   2.356194490000E+000   3.141592650000E+000   3.926990820000E+000   4.712388980000E+000   5.497787140000E+000
  ZETA_COIL_KTS(  2,:) =    9.424778000000E-001   9.424778000000E-001   9.424778000000E-001   9.424778000000E-001   9.424778000000E-001   9.424778000000E-001   9.424778000000E-001   9.424778000000E-001
!----- COIL  3
  RHO_COIL_KTS(  3,:) =    3.000000000000E-001   3.000000000000E-001   3.000000000000E-001   3.000000000000E-001   3.000000000000E-001   3.000000000000E-001   3.000000000000E-001   3.000000000000E-001
  THETA_COIL_KTS(  3,:) =    0.000000000000E+000   7.853981600000E-001   1.570796330000E+000   2.356194490000E+000   3.141592650000E+000   3.926990820000E+000   4.712388980000E+000   5.497787140000E+000
  ZETA_COIL_KTS(  3,:) =    1.570796330000E+000   1.570796330000E+000   1.570796330000E+000   1.570796330000E+000   1.570796330000E+000   1.570796330000E+000   1.570796330000E+000   1.570796330000E+000
!----- COIL  4
  RHO_COIL_KTS(  4,:) =    3.000000000000E-001   3.000000000000E-001   3.000000000000E-001   3.000000000000E-001   3.000000000000E-001   3.000000000000E-001   3.000000000000E-001   3.000000000000E-001
  THETA_COIL_KTS(  4,:) =    0.000000000000E+000   7.853981600000E-001   1.570796330000E+000   2.356194490000E+000   3.141592650000E+000   3.926990820000E+000   4.712388980000E+000   5.497787140000E+000
  ZETA_COIL_KTS(  4,:) =    2.199114860000E+000   2.199114860000E+000   2.199114860000E+000   2.199114860000E+000   2.199114860000E+000   2.199114860000E+000   2.199114860000E+000   2.199114860000E+000
!----- COIL  5
  RHO_COIL_KTS(  5,:) =    3.000000000000E-001   3.000000000000E-001   3.000000000000E-001   3.000000000000E-001   3.000000000000E-001   3.000000000000E-001   3.000000000000E-001   3.000000000000E-001
  THETA_COIL_KTS(  5,:) =    0.000000000000E+000   7.853981600000E-001   1.570796330000E+000   2.356194490000E+000   3.141592650000E+000   3.926990820000E+000   4.712388980000E+000   5.497787140000E+000
  ZETA_COIL_KTS(  5,:) =    2.827433390000E+000   2.827433390000E+000   2.827433390000E+000   2.827433390000E+000   2.827433390000E+000   2.827433390000E+000   2.827433390000E+000   2.827433390000E+000
!----------------------------------------------------------------------
!          The following variable turn off variation of
!          RHO, THETA, or ZETA respectively if set to T.
!          Thus one can mimic COILOPT by setting LFIX_RHO_COIL=T.
!----------------------------------------------------------------------
  LFIX_RHO_COIL = F
  LFIX_THETA_COIL = F
  LFIX_ZETA_COIL = F
!----------------------------------------------------------------------
!          POINCARE PLOTS (after every minimium found)
!----------------------------------------------------------------------
  LPOINCARE = T
!----------------------------------------------------------------------
!          TARGET BNORMAL ON PLASMA SURFACE
!----------------------------------------------------------------------
  NU_BNORMAL = 128
  NV_BNORMAL = 64 ! Over field period
  TARGET_BNORMAL =    0.000000000000E+000
  SIGMA_BNORMAL =    1.000000000000E+002
!----------------------------------------------------------------------
!          TARGET COIL CURVATURE
!----------------------------------------------------------------------
  TARGET_COIL_CURVATURE =    0.000000000000E+000
  SIGMA_COIL_CURVATURE =    1.000000000000E+000
!----------------------------------------------------------------------
!          TARGET COIL TORSION
!----------------------------------------------------------------------
  TARGET_COIL_TORSION =    0.000000000000E+000
  SIGMA_COIL_TORSION =    1.000000000000E+000
!----------------------------------------------------------------------
!          TARGET COIL-COIL DISTANCE
!----------------------------------------------------------------------
  TARGET_COILCOIL_DISTANCE =    0.000000000000E+000
  SIGMA_COILCOIL_DISTANCE =    1.000000000000E+000
!----------------------------------------------------------------------
!          TARGET BNORMAL HARMONICS (n,m)
!----------------------------------------------------------------------
  TARGET_BNMNS(000,001) =    0.000000000000E+000  SIGMA_BNMNS(000,001) =    1.000000000000E+000
    TARGET_BNMNC(000,001) =    0.000000000000E+000    SIGMA_BNMNC(000,001) =    1.000000000000E+000
  TARGET_BNMNS(003,001) =    0.000000000000E+000  SIGMA_BNMNS(003,001) =    1.000000000000E+000
    TARGET_BNMNC(003,001) =    0.000000000000E+000    SIGMA_BNMNC(003,001) =    1.000000000000E+000
  TARGET_BNMNS(002,002) =    0.000000000000E+000  SIGMA_BNMNS(002,002) =    1.000000000000E+000
    TARGET_BNMNC(002,002) =    0.000000000000E+000    SIGMA_BNMNC(002,002) =    1.000000000000E+000
  TARGET_BNMNS(003,003) =    0.000000000000E+000  SIGMA_BNMNS(003,003) =    1.000000000000E+000
    TARGET_BNMNC(003,003) =    0.000000000000E+000    SIGMA_BNMNC(003,003) =    1.000000000000E+000
!----------------------------------------------------------------------
!          TARGET COIL B ON AXIS
!            This is the projection of the coil field onto the VMEC
!            magnetic axis. Only really useful for zero beta VMEC
!            runs. B.N/B = 1 if perfectly aligned.
!----------------------------------------------------------------------
  TARGET_COIL_BAXIS =    1.000000000000E+000
  SIGMA_COIL_BAXIS =    1.000000000000E+000
!----------------------------------------------------------------------
!          TARGET COIL LENGTH
!----------------------------------------------------------------------
  TARGET_COIL_LENGTH(1) =    1.000000000000E+001 SIGMA_COIL_LENGTH(1) =    1.000000000000E+000
  TARGET_COIL_LENGTH(2) =    1.000000000000E+001 SIGMA_COIL_LENGTH(2) =    1.000000000000E+000
  TARGET_COIL_LENGTH(3) =    1.000000000000E+001 SIGMA_COIL_LENGTH(3) =    1.000000000000E+000
  TARGET_COIL_LENGTH(4) =    1.000000000000E+001 SIGMA_COIL_LENGTH(4) =    1.000000000000E+000
  TARGET_COIL_LENGTH(5) =    1.000000000000E+001 SIGMA_COIL_LENGTH(5) =    1.000000000000E+000
!----------------------------------------------------------------------
!          TARGET COIL ENERGY
!----------------------------------------------------------------------
  TARGET_COIL_ENERGY(1) =    0.000000000000E+001 SIGMA_COIL_ENERGY(1) =    1.000000000000E+000
  TARGET_COIL_ENERGY(2) =    0.000000000000E+001 SIGMA_COIL_ENERGY(2) =    1.000000000000E+000
  TARGET_COIL_ENERGY(3) =    0.000000000000E+001 SIGMA_COIL_ENERGY(3) =    1.000000000000E+000
  TARGET_COIL_ENERGY(4) =    0.000000000000E+001 SIGMA_COIL_ENERGY(4) =    1.000000000000E+000
  TARGET_COIL_ENERGY(5) =    0.000000000000E+001 SIGMA_COIL_ENERGY(5) =    1.000000000000E+000
/ 
&FIELDLINES_INPUT
!-----------------------------------------------------------------------
!          FIELDLINES PARAMETERS (only used if LPOINCARE=T)
!-----------------------------------------------------------------------
 NR = 128
 NZ = 128
 NPHI = 120
 RMIN = 0.436
 RMAX = 2.436
 ZMIN = -1.0
 ZMAX = 1.0
 PHIMIN = 0.0
 PHIMAX = 2.09439510239
 MU = 0.0
 R_START   = 1.50  2.00
 Z_START   = 0.00  0.00
 PHI_START = 0.00  0.00
 PHI_END   = 2*6283.0
 NPOINC    = 6
 INT_TYPE  = 'LSODE'
 FOLLOW_TOL = 1.0E-9
 VC_ADAPT_TOL = 1.0E-7
/
&END
```

------------------------------------------------------------------------

# Execute the code

The code will execute like any other STELLOPT run.
```shell
> mpiexec --use-hw-threads xstelloptv2 input.myfilename
```
# Examine the output

When run in this mode various output files are generated.

### BNORM

These are the plasma bnormal harmonics as produced by 
the [BNORM](BNORM) code.

### BNORM_REAL

These text files contain information about the normal field on 
the plasma bounday. The first line is the total number of rows 
in the file. Each row contains the following values where U is 
a poloidal index and V is a toroidal index. The columns are:

 | Index | U | V | THETA | ZETA | PHI | R | Z | Nx | Ny | Nz | Bnormal (plasma) | Bnormal (coil) | Bnormal (total) |

Here N is the outward pointing normal vector from the surface.

### BNORM_HARM

These text files contain the Fourier decomposition of the normal 
magnetic field harmonics. The first value is the total number
of rows in the file. Only the total field is decomposed. The
colums are:

 | Index | M | N | BMN (cos) | BMN (sin) |

Note that for stellarator symmetry only the sine components
should be considered.

### BAXIS_REAL

These text files contain information about the magnetic field
along the VMEC magnetic axis. The first line is the total
number of rows in the file. The columns are:

 | V | ZETA | PHI | R | Z | Tx | Ty | Tz | Bx | By | Bz | B.N |

Here the T vector is the tangential vector to the VMEC
magnetic axis. The B vector is the coil field long the 
magnetic axis. And B.N is the projection of B along the axis.

### COIL_CURVATURE

These text file contain information about the coil curvature,
torsion, and derivatives. The first line is:

 | NCOILGROUPS | NW_COILS | NH_COILS | NS |

where NS is the number of points along each coil. The number 
of rows is then NCOIL_GROUPS*NW_COILS*NH_COILS*NS where
the columns are:

 | I | J | X | Y | Z | X' | Y' | Z' | X'' | Y'' | Z'' | X''' | Y''' | Z'''| CURVATURE | TORSION |

where I is the index over the coil groups and j is the index 
over the windings. The prime is with respect to the coil
trajectory. The curvature and torsion are local.

### COILS Files

These are the standard coils files as used by the 
[MAKEGRID](MAKEGRID) code.

