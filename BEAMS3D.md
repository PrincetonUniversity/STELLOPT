{% include head.html %}

BEAMS3D
=======

The BEAMS3D ([Matthew McMillan and Samuel A Lazerson 2014 *Plasma Phys.
Control. Fusion* **56** 095019](https://iopscience.iop.org/article/10.1088/0741-3335/56/9/095019/))
code is a guiding center particle code
capable of following both user defined ensembles of particles and
modeling neutral beam injection in parallel. The magnetic field is
represented by a three dimensional splines over a cylindrical grid. It
is currently interfaced to the MAKEGRID coils file, MAKEGRID output
file, and VMEC equilibria. Ionization and recombination models are
provided by ADAS.

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

![](images/depositing.png)

The BEAMS3D code follows the guiding center orbit equations on a
cylindrical grid
$$ \frac{d\vec{R}}{dt}=\frac{\hat{b}}{qB}\left(\mu\nabla B
+\frac{mv_{ll}^2}{2B}\left(\hat{b}\cdot\nabla\right)\vec{B}\right)+v_{ll}\hat{b} $$,
$$ \frac{dv_{ll}}{dt}=-\frac{\mu}{m}\hat{b}\cdot\left(\nabla
B\right) $$. These ODE\'s can be solved via a NAG routine,
LSODE, or Runge-Kutta algorithm. The magnetic field is splined over the
cylindrical grid (R,phi,Z). The initial position and velocity of the
particles can either be specified or modeled using a neutral beam model.
The neutral beam model relies on ADAS for ionization and recombination
physics. $$ \mu = \frac{mv_\perp^2}{2B} $$

------------------------------------------------------------------------

### Compilation

BEAMS3D is distributed as part of the STELLOPT package of codes through
Git.

------------------------------------------------------------------------

### Input Data Format

The BEAMS3D code is controlled through command line inputs and an input
namelist which should be placed in the input.ext file. While the entire
VMEC input name is not required some parts will be read. The name lists
should look like:

```fortran
&INDATA
 ! VMEC input namelist (only need coil currents for runs with vacuum fields)
 EXTCUR(1) = 10000.00
 EXTCUR(2) = 10000.00
 EXTCUR(3) = 12000.00
 EXTCUR(4) = 12000.00
 EXTCUR(5) = 6000.00
/
&BEAMS3D_INPUT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!            GRID PARAMETERS                                        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 NR = 201                          ! Number of radial gridpoints, overridden if using mgrid
 NPHI = 36                         ! Number of toroidal gridpoints, overridden if using mgrid
 NZ = 201                          ! Number of vertical gridpoints, overridden if using mgrid
 RMIN = 2.5                        ! Minimum extent of radial grid, overridden if using mgrid
 RMAX = 5.0                        ! Maximum extent of radial grid, overridden if using mgrid
 ZMIN = -1.5                       ! Minimum extent of vertical grid, overridden if using mgrid
 ZMAX = 1.5                        ! Maximum extent of radial grid, overridden if using mgrid
 PHIMIN = 0.0                      ! Minimum extent of toroidal grid, overridden if using mgrid
 PHIMAX = 1.2566370614             ! Maximum extent of toroidal grid, overridden if using mgrid
 VC_ADAPT_TOL = 1.0E-3             ! Virtual casing tolerance (for plasma field outside)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!            PLASMA PARAMETERS                                      !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 PLASMA_MASS  = 1.6726231E-27       ! Mean plasma mass [kg]
 PLASMA_ZAVG  = 1.0                 ! <Z>
 PLASMA_ZMEAN = 1.0                 ! [Z]
 TE_SCALE  = 1.0                    ! Electron Temperature Scaling factor
 TE_AUX_S  = 0.0 0.5 1.0            ! Electron Temperature Knots [0,1]
 TE_AUX_F  = 0.0 1.0 2.0            ! Electron Temperature [eV]
 NE_SCALE  = 1.0                    ! Electron Density Scaling factor
 NE_AUX_S  = 0.0 0.5 1.0            ! Electron Density Knots [0,1]
 NE_AUX_F  = 0.0 1.0 2.0            ! Electron Density [m^-3]
 TI_SCALE  = 1.0                    ! Ion Temperature Scaling factor
 TI_AUX_S  = 0.0 0.5 1.0            ! Ion Temperature Knots [0,1]
 TI_AUX_F  = 0.0 1.0 2.0            ! Ion Temperature [eV]
 POT_AUX_S  = 0.0 0.5 1.0           ! Electrostatic Potential Knots [0,1]
 POT_AUX_F  = 0.0 1.0 2.0           ! Electrostatic Potential [V] (Phi, not dPhi/dr)
 THERM_FACTOR = 1.5                 ! Factor at which to thermalize (Vtherm*THERM_FACTOR)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!            PLASMA PARAMETERS (MULTI-ION)                          !!
!!     Note: This overrides the above profile deffintion             !!
!!           of PLASMA_MASS, PLASMA_ZAVG, PLASMA_ZMEAN               !!
!!           UNITS OF M [kg]; UNITS OF NI [m^-3]                     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  NI_AUX_Z = 1 1 2 6 ! D T He C
  NI_AUX_M = 3.3435837724E-27 5.008267217094E-27 6.6464764063E-27 1.994423482E-26
  NI_AUX_S =     0.0     0.2     0.4     0.6     0.8    1.0
  NI_AUX_F(1,:) =  9.60E+18  8.00E+18  6.40E+18  4.80E+18  3.20E+18  1.60E+18
  NI_AUX_F(2,:) =  9.60E+18  8.00E+18  6.40E+18  4.80E+18  3.20E+18  1.60E+18
  NI_AUX_F(3,:) =  4.56E+18  3.80E+18  3.04E+18  2.28E+18  1.52E+18  7.60E+17
  NI_AUX_F(4,:) =  2.40E+17  2.00E+17  1.60E+17  1.20E+17  8.00E+16  4.00E+16
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!            DISTRIBUTION FUNCTION                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 NRHO_DIST   = 64                   ! Radial bins (0,1)
 NTHETA_DIST = 4                    ! Poloidal bins (0,2*pi)
 NZETA_DIST  = 4                    ! Toroidal bins (0,2*pi), 4*NFP is a good option
 NVPARA_DIST = 32                   ! Parallel velocity bins (-vmin,vmax)
 NVPERP_DIST = 64                   ! Perpendicular velocity bins (0,vmax)
 PARTVMAX    = 3.0E6                ! Maximum velocity in dist. (vmax)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!            PARTICLE INTEGRATION PARAMETERS                        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 INT_TYPE = 'LSODE'                 ! Particle trajectory integration method (NAG, RKH68, LSODE)
 FOLLOW_TOL = 1.0E-12               ! Trajectory following tolerance (NAG, LSODE)
 NPOINC = 100                       ! Number of trajector points to save per particle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!            PARTICLE INITIAL CONDITION (INDIVIDUAL)                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 R_START_IN =  3.6  3.7  3.8       ! Radial starting locations of particles [m]
 Z_START_IN =  0.0  0.0  0.0       ! Vertical starting locations of particles [m]
 PHI_START_IN =  0.0  0.0  0.0     ! Toroidal starting locations of particles (radians)
 VLL_START_IN =  1.0E6 1.0E6 1.0E6 ! Initial parallel velocity of particles [m/s]
 MU_START_IN =   3*1.0E-15         ! Particle magnetic moment [J/T] (0.5*m*v^2/B)
 CHARGE_IN   =   3*1.60217733E-19  ! Particle charge [C]
 MASS_IN     =   3*1.6726231E-27   ! Particle mass [kg]
 ZATOM_IN    =   3*1.0             ! Particle charge number
 T_END_IN    =   3*0.001           ! How long to follow particles [s]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!            PARTICLE INITIAL CONDITION (Neutral Beam)              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 NPARTICLES_START = 1024           ! Number of particles per beamline
 R_BEAMS(1,1:2) = 1.673 1.400      ! R chord deffinition (start,end) BEAM#1
 Z_BEAMS(1,1:2) = 0.000 0.000      ! Z chord deffinition (start,end) BEAM#1
 PHI_BEAMS(1,1:2) = 0.000 0.000    ! PHI chord deffinition (start,end) BEAM#1
 ASIZE_BEAMS(1)   = 0.15           ! Aperature Size [m] BEAM#1
 ADIST_BEAMS(1)   = 2.00           ! Aperature Distance [m] BEAM#1
 DIV_BEAMS(1)     = 0.1            ! Beam divergence [rad] BEAM#1
 E_BEAMS(1)       = 0.1            ! Beam Energy [J] BEAM#1
 MASS_BEAMS(1)    = 1.6726231E-27  ! Beam species mass [kg] BEAM#1
 CHARGE_BEAMS(1)  = 1.0            ! Beam species charge [C] BEAM#1
 ZATOM_BEAMS(1)   = 1.0            ! Beam species Z [norm] BEAM#1
 T_END_IN(1)      = 0.05           ! How long to follow particles [s]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!            PARTICLE INITIAL CONDITION (Fusion Reactions)          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUSION_SCALE     = 1.0            ! Scaleing Factor to apply to energy
/
&END
```

It is important to note that neutral beam lines are defined by two
points. The first index of the array is the beam number, the second
defines the two points. The first point (X,1) is the origin on the beam.
The second points (X,2) defines the beamline from that origin.

![ Diagram of neutral beam parameters showing how ASIZE\_BEAMS, DIV\_BEAMS, and ADIST\_BEAMS are used to initialise neutral particles.](images/beam_diagram.jpg )

------------------------------------------------------------------------

### Execution

The BEASM3D code is controlled through a combination of command-line
inputs and an input namelist. The input namelist must be in the
equilibrium input file. That file must also contain the VMEC INDATA
namelist (although only the EXTCUR array will be used). The BEASM3D code
is run from the command line taking an equilibrium input file as a
necessary argument. This input file must have the INDATA (for EXTCUR)
and BEASM3D\_INPUT namelists in it.

    BEASM3D -vmec <VMEC FILE> -coil <COIL FILE> -mgrid <MGRID FILE> -vessel <VESSEL FILE> -vac -full -noverb -help

| Argument | Default | Description |
|:------------- |:-------------:|:----- |
| -vmec | NONE | VMEC input extension |
| -coil | NONE | Coils File |
| -mgrid | NONE | Makegrid style vacuum grid file |
| -vessel | NONE | First wall file |
| -beamlet | NONE | Beamlet deffintion HDF5 file. |
| -restart | NONE | Restart run from particles in previous run (HDF5 file) |
| -vac | NONE | Only compute the vacuum field |
| -beam_simple | NONE | Assume monoenergetic beams (normally 1% variance around injection energy) |
| -collisions | NONE | Force use of slowing down/scattering operator. |
| -depo | NONE | Calculate deposition only |
| -field | NONE | Outputs the B-Field on the cylindrical grid only. |
| -ascot4 | NONE | Creates input HDF5 file for ASCOT4 (BBNBI, no particles) |
| -ascot5 | NONE | Creates input HDF5 file for ASCOT5. |
| -hitonly | NONE | Only save vessel strike points.|
| -plasma | NONE | Only compute fields inside the plasma domain (places wall at LCFS) |
| -raw | NONE | Treats EXTCUR array as raw values (EXTCUR is a scale factor applied to what\'s in the coils file). |
| -suzuki | NONE | Use Suzuki beam deposition model (default if no ADAS/PREACT). |
| -w7x | NONE | Use W7-X beam shape model. |
| -fusion | NONE | Use nuclear fusion thermal birth model. |
| -fusion_alpha | NONE | Use nuclear fusion thermal birth model (alphas only). |
| -noverb | NONE | Suppresses screen output |
| -help | NONE | Print help message. |

In it\'s simplest invokation the code requires
a VMEC input file.

    >~/bin/xbeams3d -vmec ncsx_c09r00_free -mgrid mgrid_c09r00.nc -vac
    BEAMS3D Version 1.00
    ----- Particle Initialization -----
       S   = [  0.24490,  0.24490];   NS:        1
       U   = [  0.00000,  5.96903];   NU:     20
       V   = [  0.00000,  5.96903];   NV:     20
       V_||= [*********,  2.00E+05];  NP:     10
       Mu  = [  0.00000,  2.00E-15];  NP:     10
    ----- Profile Initialization -----
       Ne  = [   0.00,   5.00] 10^19 [m^-3];  Nne:     50
       Te  = [   0.00,   4.58] [keV];  Nte:     50
       Ti  = [   0.00,   4.58] [keV];  Nti:     50
    ----- Input Parameters -----
       R   = [  0.95223,  1.82142];  NR:    201
       PHI = [ 0.00000, 2.09440];  NPHI:   60
       Z   = [-0.66634, 0.66634];  NZ:    201
       # of Particles to Start:   4000
    ----- Constructing Splines -----
       R   = [  0.95223,  1.82142];  NR:    201
       PHI = [ 0.00000, 2.09440];  NPHI:   60
       Z   = [-0.66634, 0.66634];  NZ:    201
       HERMITE FORM: 1
    ----- FOLLOWING PARTICLE TRAJECTORIES -----
          Method: NAG
      Particles: 4000
           Steps:     320   Delta-t: 0.3125E-06
          NPOINC:     320    dt_out: 0.1000E-06
             Tol: 0.1000E-08  Type: M
    ----- WRITING DATA TO FILE -----
       FILE: beams3d_ncsx_c09r00_free.h5
    ----- BEAMS3D DONE -----

The BENCHMARKS directory contains a set of tests for BEAMS3D based on an
axisymmetric VMEC tokamak equilibrium.  The input files are located
in the BEAMS3D_TEST subdirectory.  They can all be invoked by
calling make beams3d_test from the BENCHMARKS directory.  The
comparrision scripts require Python.

------------------------------------------------------------------------

### Output Data Format

The data from each run is output into a HDF5 file where all datasets
located at the root level.  The following table helps to define the
variables (all values in mks units, angles in radians)

| Name | Type | Size | Description |
| :--- | :--- | :---: | :--- |
| **Background Grid** |
| nr | INTEGER | 1 | Number of radial background gridpoints |
| nphi | INTEGER | 1 | Number of toroidal background gridpoints |
| nz | INTEGER | 1 | Number of vertical background gridpoints |
| nion | INTEGER | 1 | Number of background ion species |
| plasma_mass | DOUBLE | 1 | Average plasma background mass |
| plasma_zavg | DOUBLE | 1 | Average plasma \<Z\> |
| plasma_zavg | DOUBLE | 1 | Average plasma \[Z\] |
| raxis | DOUBLE | nr | R values of background grid |
| phiaxis | DOUBLE | nphi | Phi values of background grid |
| zaxis | DOUBLE | nz | Z values of background grid |
| BR_ARR | DOUBLE | nr,nphi,nz | Magnetic field (B_R) |
| BPHI_ARR | DOUBLE | nr,nphi,nz | Magnetic field (B_PHI) |
| BZ_ARR | DOUBLE | nr,nphi,nz | Magnetic field (B_Z) |
| S_ARR | DOUBLE | nr,nphi,nz | Normalized Toroidal Flux (s) |
| U_ARR | DOUBLE | nr,nphi,nz | Poloidal-like angle (u) |
| POT_ARR | DOUBLE | nr,nphi,nz | Electrostatic scalar potential |
| NE | DOUBLE | nr,nphi,nz | Electron number density |
| TE | DOUBLE | nr,nphi,nz | Electron Temperature eV |
| NI | DOUBLE | nion,nr,nphi,nz | Ion number density |
| TI | DOUBLE | nr,nphi,nz | Ion Temperature eV |
| ZEFF_ARR | DOUBLE | nr,nphi,nz | Zeff |
| **Marker Trajectory** |
| npoinc | INTEGER | 1 | Number of Timesteps Saved |
| nparticles | INTEGER | 1 | Number of markers Evolved |
| mass | DOUBLE | nparticles | Fast Ion Masses |
| charge | DOUBLE | nparticles | Fast Ion Charge |
| weight | DOUBLE | nparticles | Fast Weight (1/s) |
| Beam | INTEGER | nparticles | Population index |
| Zatom | INTEGER | nparticles | Particle Charge Number |
| end_state | INTEGER | nparticles | Particle End State (0: Orbiting; 1: Thermalized; 2: Wall Strike; 3: Shine-through; 4: Port-Load) |
| t_end | DOUBLE | nparticles | Particle Time of flight |
| R_lines | DOUBLE | npoinc+1,nparticles | R trajectory of markers. |
| PHI_lines | DOUBLE | npoinc+1,nparticles | Phi trajectory of markers. |
| Z_lines | DOUBLE | npoinc+1,nparticles | Z trajectory of markers. |
| vll_lines | DOUBLE | npoinc+1,nparticles | Parallel velocity trajectory of markers. |
| moment_lines | DOUBLE | npoinc+1,nparticles | Magnetic Moment trajectory of markers. |
| neut_lines | BOOLEAN | npoinc+1,nparticles | If true markers is a neutral at that point. |
| S_lines | DOUBLE | npoinc+1,nparticles | Normalized toroidal flux rajectory of markers. |
| U_lines | DOUBLE | npoinc+1,nparticles | Poloidal angle trajectory of markers. |
| B_lines | DOUBLE | npoinc+1,nparticles | mod(B) trajectory of markers. |
| R_lines | DOUBLE | npoinc+1,nparticles | R Trajectory of markers. |
| **Distribution Function** |
| nbeams | INTEGER | 1 | Number of fast ion populations |
| ns_prof1 | INTEGER | 1 | Number of radial distribution gridpoints |
| ns_prof2 | INTEGER | 1 | Number of poloidal distribution gridpoints |
| ns_prof3 | INTEGER | 1 | Number of toroidal distribution gridpoints |
| ns_prof4 | INTEGER | 1 | Number of parallel velocity distribution gridpoints |
| ns_prof5 | INTEGER | 1 | Number of perpendicular velocity distribution gridpoints |
| partvmax | DOUBLE | 1 | Maximum velocity of distribution function. |
| dist_prof | DOUBLE | nbeams,ns_prof1..5 | Distribution function. (part*m^-3*s^-3, no physical volume) |
| ndot_prof | DOUBLE | nbeams,ns_prof1 | Fast Ion Source (m^-3/s) |
| epower_prof | DOUBLE | nbeams,ns_prof1 | Electron Heating W/m^3 |
| ipower_prof | DOUBLE | nbeams,ns_prof1 | Ion Heating W/m^3 |
| j_prof | DOUBLE | nbeams,ns_prof1 | Fast Ion Current A/m^2 |
| dense_prof | DOUBLE | nbeams,ns_prof1 | Fast Ion Density m^-3 |
| Gfactor | DOUBLE | ns_prof1 | Correction factor (1-L31)/Zeff to obtain NBCD. |
| **NBI** |
| Energy | DOUBLE | nbeams | Fast Ion Population Initial Energy |
| V_neut | DOUBLE | 3,nparticles | Initial neutral velocity (Vx,Vy,Vz) |
| Shinethrough | DOUBLE | nbeams | Beam Shinethrough % |
| Shineport | DOUBLE | nbeams | Beam loss to port % |
| **Wall Model** |
| nvertex | INTEGER | 1 | Number of vertices in wall model |
| nface | INTEGER | 1 | Number of faces in wall model |
| wall_vertex | DOUBLE | nvertex,3 | Vertex locations (x,y,z) |
| wall_faces | INTEGER | nface,3 | Face Vertex Indicies (v1,v2,v3) |
| wall_strikes | INTEGER | nface | Number of marker strikes on a given wall model face. |
| wall_load | DOUBLE | nface | Heat flux to a given wall model face W/m^2 |
| wall_shine | DOUBLE | nface | Shinethrough heat flux to a given wall model face W/m^2 |

**Background Grids** These are the background quantities the code uses
for pushing the markers.  The `raxis`, `phiaxis`, and `zaxis` arrays
hold the locations of the gridpoints.

**Marker Trajectories** These arrays hold the state of every marker for
each timestep.  The timesteps are defined as the maximum time 
(`t_end` from the input) divided by `npoinc`.  For neutral beam
simulations the first time index is the initial position of the neutral,
the second index is either a collision with the wall or the point at
which the neutral ionizes, and the thrid index is the initial condition
for the gyrocenter.  For a slowing down run check the `end_state` array
if `t_end` was set long enough there should be no orbiting markers.
When a marker hits a wall structure we save two states.  The last index
is the position of the strike on the wall element.  The preceeding index
is the previous sub-timestep.  This allows us to reconstruct the ray
used for determining the wall strike.

**Distribution Functions** These arrays hold information from slowing 
down distribution information. These arrays are defined on a rho grid
where $$\rho=\sqrt{s}$$. The `ns_profX` variable indicate the number of
bins.  The extents are $$\rho=\left[0,1\right]$$, $$u=\left[0,2\pi\right]$$, $$\phi=\left[0,2\pi\right]$$,
$$v_{para}=\left[-partvmax,partvmax\right]$$, and $$v_{perp}=\left[0,partvmax\right]$$. Also
note that the 5D distribution fucntion is not normalized by volume,
the user should do this before attempting to interpolate to another
grid.

**Wall Model** The wall model is stored for each run.  The `wall_faces`
array contains the three indices into the `wall_vertex` array which
define each triangular face.

------------------------------------------------------------------------

### Visualization

The data from the HDF5 file may be easily plotted in many plotting
packages. The R , PHI, Z, and V\_PARALLEL coordinates of each particle
are save NPOINC times along the trajectory. To aid in plotting the
evolution of the distribution function the code outputs a textfile which
bins by VLL the particles at each NPOINC time step.

[Movie of Simulation](https://www.youtube.com/watch?v=8TdZouWdmNY)

------------------------------------------------------------------------

### Tutorials

[NCSX Neutral Beam Injection Example](NCSX Neutral Beam Injection Example.md)

[Benchmarking and Validation](BEAMS3D Validation and Benchmarking on HPC systems.md)

------------------------------------------------------------------------

### References

-   [McMillan, M. and Lazerson, S.A. \"BEAMS3D: Neutral beam injection model.\" Plasma Phys. and Control. Fusion 56, 095019 (2014)](https://doi.org/10.1088/0741-3335/56/9/095019)
-   [Lazerson, S.A. et al. \"Validation of the BEASM3D neutral beam deposition model on Wendelstein 7-X\" Nuclear Fusion 60, 706020 (2020)](https://doi.org/10.1088/1741-4326/ab8e61)
-   [Lazerson, S.A. et al. \"Modeling and measurement of energetic particle slowing down on Wendelstein 7-X\" Nuclear Fusion 61, 096006 (2021)](https://doi.org/10.1088/1741-4326/ac0771)

