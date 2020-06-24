Tutorial: STELLOPT Energetic Particle Optimization
==================================================
![](images/BEAMS3D_orbits.png)

The [STELLOPT](STELLOPT) code can target energetic particle confinement
using the [BEAMS3D](BEAMS3D) neutral beam package as a collision less
gyro-particle integrator. When run in this manner the
[STELLOPT](STELLOPT) code initialises the particle starting points and
energies for [BEAMS3D](BEAMS3D) (skipping over the neutral beam
injection parts of the code). For equilibrium run in fixed boundary
mode, [STELLOPT](STELLOPT) constructs a tessellated wall at the
[VMEC](VMEC) boundary. All particles which hit this wall are considered
lost.

------------------------------------------------------------------------

1.  **Input namelists**
    In addition to the OPTIMUM name list
    the BEAMS3D\_INPUT name list must be included in the STELLOPT input
    file (also the INDATA namelist). The STELLOPT namelist should
    include the NPOPULATION parameter, this determines how to divide up
    the processors. For example if you has 256 processors available and
    NPOPULATION was set to 8, STELLOPT would execute with 8 processors
    doing the optimization while each of those 8 processors would have
    an addition 31 processors sitting around to help run BEAMS3D. Thus
    you could have up to 8 copies of BEAMS3D running with 32 processors
    each. Below is an example (truncated) set of namelists: >


```fortran
&INDATA
! Whatever initial equilibrium you want to start with
/
&OPTIMUM
!-----------------------------------------------------------------------
!          OPTIMIZER RUN CONTROL PARAMETERS
!-----------------------------------------------------------------------
  NFUNC_MAX = 1
  EQUIL_TYPE = 'VMEC2000'
  OPT_TYPE   = 'ONE_ITER'
  FTOL =  1.00000000000000E-06
  XTOL =  1.00000000000000E-30
  GTOL =  1.00000000000000E-30
  FACTOR =   100.0
  EPSFCN =   1.0E-05
  LKEEP_MINS = T
  NOPTIMIZERS = 1
!-----------------------------------------------------------------------
!          OPTIMIZED QUANTITIES
!-----------------------------------------------------------------------
  LCURTOR_OPT = T   CURTOR_MIN  = -1.0E6  CURTOR_MAX  = -1.0E4
!------------------------------------------------------------------------
!       Particle Transport
!       Mu = 0.5*m*v*v/B or (Thermal Energy)/B
!------------------------------------------------------------------------
  NU_ORBIT = 20
  NV_ORBIT = 20
  NP_ORBIT = 20
  VLL_ORBIT = 20*7.76E5  ! 50 [keV]
  VLL_ORBIT = 20*5.00E5
  MASS_ORBIT = 1.6726219E-27
  Z_ORBIT    = 1.0
  ! 5-> 100[ev]  num2str(((5:95/19:100).*1000)*1.60217733E-19./1.5,'%20.10E')
  MU_ORBIT  = 5.3405911000E-16    1.0681182200E-15    1.6021773300E-15    2.1362364400E-15    2.6702955500E-15
    3.2043546600E-15    3.7384137700E-15    4.2724728800E-15    4.8065319900E-15    5.3405911000E-15
    5.8746502100E-15    6.4087093200E-15    6.9427684300E-15    7.4768275400E-15    8.0108866500E-15
    8.5449457600E-15    9.0790048700E-15    9.6130639800E-15    1.0147123090E-14    1.0681182200E-14
  TARGET_ORBIT(10) = 0.0   SIGMA_ORBIT(10) = 1.0 
/
&BEAMS3D_INPUT
  NR = 128
  NZ = 128
  NPHI = 36
  RMIN =  4.36000000000000E-01
  RMAX =  2.43600000000000E+00
  ZMIN = -1.00000000000000E+00
  ZMAX =  1.00000000000000E+00
  PHIMIN =  0.00000000000000E+00
  PHIMAX =  2.09439510239000E+00
  NPOINC = 200
  T_END_IN = 2*0.003
  INT_TYPE = 'LSODE'
  FOLLOW_TOL =  1.00000000000000E-09
  VC_ADAPT_TOL =  1.00000000000000E-06
  R_START_IN(1) = 1.0 ! need this to turn off BEAM stuff
/
&END
```
In the above example the VMEC flux surface ns=10 is initialized at 400
    (20x20) unique poloidal and toroidal points. At each of these points
    20 particles are launched with mangetic moments (MU\_ORBIT) and
    parallel velocities (VLL\_ORBIT) as specified in associated arrays.
    In this way the user can set the pitch angles which are evaluated,
    the result being 8000 particles being launched from surface 10. If
    an additional surface had been specified (through TARGET_ORBIT and
    SIGMA\_ORBIT) the total number of particles followed would be 16,000
    (and so on). It should be noted that while the STELLOPT electron
    temperature and density will be read, it is not used as only
    collisionless particle orbits are currently followed. Alternatively,
    the user may set an array called VPERP\_ORBIT which then overrides
    MU\_ORBIT. The equilibirum modB, VPERP\_ORBIT and MASS\_ORBIT are
    then used to calculate MU\_ORBIT for the run.
2.  **Execute the code**
    The coupled BEAMS3D/STELLOPT codes
    will execute like any other STELLOPT run (note in this example
    we've used the SINGLE\_ITER optimization type and set
    NPOPULATION=1) >

```
> mpiexec -np 2 /Users/lazerson/bin/xstelloptv2 input.BEAMS3D
STELLOPT Version  2.70
  Equilibrium calculation provided by: 
  =================================================================================
  =========   Parallel Variational Moments Equilibrium Code (v 9.0)       =========
  =========                (S. Hirshman, J. Whitson)                      =========
  =========         http://vmecwiki.pppl.wikispaces.net/VMEC              =========
  =================================================================================
     
  Energetic Particle calculation provided by: 
  =================================================================================
  =========                      BEAMS3D (v 2.70)                         =========
  =========                  (M. McMillan, S. Lazerson)                   =========
  =========                       lazerson@pppl.gov                       =========
  =========          http://vmecwiki.pppl.wikispaces.net/BEAMS3D          =========
  =================================================================================
     
 -----  MPI Params.   -----
    MPI Version:  3.01
    Optimizers requested:            1
    Number of Processors:            2
    Shared memory groups:            1
    Processors per group:            2
   Workers per optimizer:            2
     Optimizers provided:            1
 -----  Optimization  -----
    =======VARS=======
     CURTOR:   Total Toroidal Current
    ======TARGETS=====
     Particle Orbits (BEAMS3D)
    ==================
    Number of Parameters:            1
    Number of Targets:               1
    !!!! EQUILIBRIUM RESTARTING NOT UTILIZED !!!!
     OPTIMIZER: SINGLE_ITERATION
     NFUNC_MAX:          100
 ---------------------------  EQUILIBRIUM CALCULATION  ------------------------

  NS =    9 NO. FOURIER MODES =   94 FTOLV =  1.000E-06 NITER =   1000
  PROCESSOR COUNT - RADIAL:    2
 INITIAL JACOBIAN CHANGED SIGN!
 TRYING TO IMPROVE INITIAL MAGNETIC AXIS GUESS
  ---- Improved AXIS Guess ----
      RAXIS_CC =    1.4888391525837754       0.11677945389900736       -2.7940380117869989E-003  -1.1018895172466669E-003  -2.6280314029721459E-005   1.8347864362562431E-003
      ZAXIS_CS =   -0.0000000000000000       -3.6771928287768728E-002   1.2315609423830821E-002   8.3228855092970868E-004  -8.0307918459706570E-003  -2.9798398326088293E-003
  -----------------------------


  ITER    FSQR      FSQZ      FSQL    RAX(v=0)    DELT       WMHD


    1  1.06E+00  1.63E-01  1.71E-01  1.604E+00  9.00E-01  3.7586E+00
  126  8.45E-07  2.61E-07  2.43E-07  1.582E+00  7.29E-01  3.6371E+00

  NS =   29 NO. FOURIER MODES =   94 FTOLV =  1.000E-08 NITER =   2000
  PROCESSOR COUNT - RADIAL:    2

  ITER    FSQR      FSQZ      FSQL    RAX(v=0)    DELT       WMHD

    1  4.11E-02  1.73E-02  3.13E-04  1.582E+00  9.00E-01  3.6367E+00
  200  7.22E-08  8.49E-09  6.01E-09  1.580E+00  7.29E-01  3.6365E+00
  286  9.46E-09  1.44E-09  1.38E-09  1.579E+00  7.29E-01  3.6365E+00

  NS =   49 NO. FOURIER MODES =   94 FTOLV =  1.000E-10 NITER =   4000
  PROCESSOR COUNT - RADIAL:    2

  ITER    FSQR      FSQZ      FSQL    RAX(v=0)    DELT       WMHD

    1  7.84E-03  4.35E-03  4.59E-06  1.579E+00  9.00E-01  3.6365E+00
  200  1.19E-08  2.81E-09  9.58E-10  1.579E+00  6.87E-01  3.6365E+00
  400  1.15E-09  2.06E-10  2.02E-10  1.578E+00  6.87E-01  3.6365E+00
  600  3.18E-10  5.68E-11  3.77E-11  1.578E+00  6.87E-01  3.6365E+00
  723  1.00E-10  1.79E-11  1.08E-11  1.578E+00  6.87E-01  3.6365E+00

  NS =   99 NO. FOURIER MODES =   94 FTOLV =  1.000E-12 NITER =  10000
  PROCESSOR COUNT - RADIAL:    2

  ITER    FSQR      FSQZ      FSQL    RAX(v=0)    DELT       WMHD

    1  2.21E-02  1.27E-02  1.79E-06  1.578E+00  9.00E-01  3.6365E+00
  200  2.46E-08  9.25E-09  1.18E-10  1.578E+00  5.57E-01  3.6365E+00
  400  9.75E-10  2.45E-10  1.71E-11  1.578E+00  5.57E-01  3.6365E+00
  600  1.24E-10  3.59E-11  8.24E-12  1.578E+00  5.57E-01  3.6365E+00
  800  3.70E-11  1.24E-11  3.43E-12  1.578E+00  5.57E-01  3.6365E+00
 1000  1.30E-11  4.65E-12  1.00E-12  1.578E+00  5.57E-01  3.6365E+00
 1200  4.45E-12  1.67E-12  2.64E-13  1.578E+00  5.57E-01  3.6365E+00
 1400  1.82E-12  6.66E-13  7.02E-14  1.578E+00  5.57E-01  3.6365E+00
 1522  9.97E-13  3.51E-13  2.99E-14  1.578E+00  5.57E-01  3.6365E+00

 EXECUTION TERMINATED NORMALLY

 FILE : BEAMS3D_opt0
 NUMBER OF JACOBIAN RESETS =    4

    TOTAL COMPUTATIONAL TIME (SEC)        16.89
    TIME TO INPUT/OUTPUT                   0.06
       READ IN DATA                        0.00
       WRITE OUT DATA TO WOUT              0.06
    TIME IN FUNCT3D                       16.65
       BCOVAR FIELDS                       1.97
       FOURIER TRANSFORM                   4.07
       INVERSE FOURIER TRANSFORM           3.44
       FORCES AND SYMMETRIZE               2.11
       RESIDUE                             2.69
       EQFORCE                             0.00
 -------------------------  PARAVMEC CALCULATION DONE  -----------------------
     ASPECT RATIO:    4.365
             BETA:    0.042  (total)
                      0.584  (poloidal)
                      0.046  (toroidal)
  TORIDAL CURRENT:   -0.174295786408E+06
     TORIDAL FLUX:    0.514
           VOLUME:    2.979
     MAJOR RADIUS:    1.422
     MINOR RADIUS:    0.326
       AXIS FIELD:    1.539
    STORED ENERGY:    0.192547218957E+06
 ---------------------------  ORBIT CALCULATION  -------------------------
 ----- Particle Initialization -----
   Z   =  1
   M   =   1.67E-27
   S   = [  0.09184,  0.09184];   NS:        1
   U   = [  0.00000,  5.96903];   NU:     20
   V   = [  0.00000,  2.98451];   NV:     20
   V_||= [  5.00E+05,  5.00E+05];  NP:     20
   Mu  = [  0.00000,  1.07E-14];  NP:     20
----- Profile Initialization -----
   Ne  = [  -0.00,   0.10] 10^19 [m^-3];  Nne:     50
   Te  = [  -0.00,   0.00] [keV];  Nte:     50
   Ti  = [  -0.00,   0.00] [keV];  Nti:     50
  Zeff = [  -1.00,   1.00];      NZeff:     50
BEAMS3D Version  2.70
----- Input Parameters -----
   R   = [  0.95222,  1.82145];  NR:    128
   PHI = [ 0.00000, 2.09440];  NPHI:   36
   Z   = [-0.66631, 0.66631];  NZ:    128
   # of Particles to Start:     8000
   MAGNETIC FIELD FROM PLASMA ONLY!
----- Profile Parameters -----
   Te   = [  0.00000,  0.00100] keV;  NTE:     50
   Ti   = [  0.00000,  0.00100] keV;  NTI:     50
   Ne   = [  0.00000,  0.01000] E20 m^-3;  NNE:     50
   Zeff = [  1.00000,  1.00000];  NZEFF:   50
   PLASMA_MASS =    1.00728 amu
   PLASMA_ZAVG =    1.00000 <Z>
   PLASMA_ZMEAN =    1.00000 [Z]
----- VMEC Information -----
   FILE: BEAMS3D_opt0
   R       = [  0.99173,  1.78194]
   Z       = [-0.62680, 0.62680]
   BETA    =   0.042;  I  =  -0.174 [MA]
   AMINOR  =   0.326 [m]
   PHIEDGE =   0.514 [Wb]
   VOLUME  =   2.979 [m^3]
 -----  Vessel Information  -----
   Wall Name :  HARMONICS
   Date      :  TODAY
   Faces     :   28800
                                    
----- Constructing Splines -----
   R   = [  0.95222,  1.82145];  NR:    128
   PHI = [ 0.00000, 2.09440];  NPHI:   36
   Z   = [-0.66631, 0.66631];  NZ:    128
   HERMITE FORM: 1
----- FOLLOWING PARTICLE TRAJECTORIES -----
      Method: LSODE
   Particles:      8000
       Steps:     29999   Delta-t:  0.1000E-06
      NPOINC:       200    dt_out:  0.1500E-04
         Tol:  0.1000E-08  Type: 10
     Trajectory Calculation [  0]%
```


2.  **Examine the output**
    The STELLOPT code will run as it
    usually does, however now BEAMS3D HDF5 output files will be
    produced. Given the large size of these files, only the first file
    HDF5 file produced will contain the full trajectory of particles
    followed (NPOINC locations along each trajectory). The additional
    files will just record the positions along the boundary where
    particles left the equilibrium. The 'stellopt' file will include
    the loss fraction from each surface.
