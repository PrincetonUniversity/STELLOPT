{% include head.html %}

DKES
====

[toc](toc) The Drift Kinetic Equation Solver
(DKES)\<ref\>[W. I. van Rij and S. P. Hirshman, Variational Bounds for Transport Coefficients in Three-Dimensional Plasmas, Phys. Fluids B 1,563(1989)](http://dx.doi.org/10.1063/1.859116)\</ref\>
solves a set of 3D drift kinetic equations to obtain upper and lower
bounds for the diffusion coefficients of a prescribed toroidal plasma
equilibrium. The three dimensions are theta (poloidal angle), zeta
(toroidal angle), and alpha (pitch angle). Straight-line flux
coordinates are used to describe the equilibrium, which must satisfy the
stellarator symmetry conditions.

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

The DKES code is based on a variation principle for the linearized
drift-kinetic, Fokker-Planck equation. The formulation neglects the effect
of magnetic drifts on the orbits of the perturbed distribution function.
This neglect precludes the code from being able to study resonant
superbanana orbits, and systems with very low collisionality. Such
effects must also be considered when the radial electric field potential
satisfies $$ e\Phi/T < 1/A$$. The advantage of this formualtion is that it
reduces the problem from 5 dimension to 3. The conservation of symmetry
and conservative properties are also preserved allowing for a
variation approach to be taken. 

------------------------------------------------------------------------

### Compilation

DKES is a component of the STELLOPT suite of codes.

------------------------------------------------------------------------

### Input Data Format

The DKES code needs an 'input_dkes.ext' file which contains the `DKES_INDATA`
namelist such as

```fortran
&DKES_INDATA
 NZPERIOD=  3, ! Field Periods
 MPOL=17,      ! Number of poloidal modes in distribution
 NTOR= 26,     ! Number of toroidal modes in distribution
 MPOLB =  6,   ! Number of poloidal modes in equilibrium  
 NTORB =  8,   ! Number of toroidal modes in equilibrium
 LALPHA= 100,  ! Number of polynom. to use in pitch angle
 NRUN = 1,     ! Number of cmul/efield pairs to be calculated
 CMUL =   0.1000E-02,   ! Array of Collisionality
 EFIELD =   0.1000E+03, ! Array of Radial electric field values
 CHIP =  0.1403,        ! d(chi)/d(rho) Poloidal flux derivative
 PSIP =  0.2783,        ! d(psi)/d(rho) Toroidal flux derivative
 BTHETA =     4.6719062102E-03, ! I Toroidal current
 BZETA =     2.3228644337E+00,  ! J the Poloidal current
 IBBI = 1,     ! 1: BORBI is B; 2: BORBI is 1/B
 ! BORBI is the of nonzero Fourier coefficients for B (or 1/B)
 ! B(theta,zeta) = borbi(n,m)*cos[(m-1)*theta - nzperiod*nvalsb(n)*zeta]
 BORBI(0,0)=  0.15748E+01,
 BORBI(0,1)=  -.11885E+00,
 BORBI(0,2)=  -.19193E-01,
 BORBI(1,2)=  0.88849E-02,
 BORBI(0,3)=  0.48462E-02,
 BORBI(2,3)=  -.44171E-02,
 BORBI(-1,1)=  -.26255E-02,
 BORBI(1,1)=  -.23643E-02,
 BORBI(2,2)=  -.21185E-02,
 BORBI(3,3)=  0.20004E-02,
 BORBI(4,4)=  -.13090E-02,
 BORBI(1,0)=  -.11995E-02,
 BORBI(1,3)=  0.10712E-02,
 BORBI(1,4)=  -.99923E-03,
 BORBI(2,0)=  -.97714E-03,
 BORBI(4,3)=  -.95223E-03,
 BORBI(5,2)=  0.87245E-03,
 BORBI(3,4)=  0.86562E-03,
 BORBI(-2,1)=  -.86180E-03,
 BORBI(4,2)=  -.81482E-03,
 BORBI(3,2)=  0.77844E-03,
 BORBI(2,1)=  -.71725E-03,
 BORBI(-1,2)=  0.55397E-03,
 BORBI(3,0)=  0.51978E-03,
 BORBI(1,5)=  0.51622E-03,
 ! MMNN Matrix containing the distribution of toroidal mode numbers
!                 (expressed in units of nzperiod) in column 1, and the
!                 number of poloidal modes associated with each toroidal
!                 mode in column 2. the poloidal mode numbers associated
!                 with each toroidal mode are stored in the same row,
!                 starting in column 3. WARNING : this spectrum must
!                 encompass the spectrum of the input magnetic field.
 mmnn(1,1)=  0, 14,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
  13,
 mmnn(1,2)=  1, 15,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
  13, 14,
 mmnn(1,3)=  2, 17,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
  13, 14, 15, 16,
 mmnn(1,4)=  3, 16,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
  13, 14, 15,
 mmnn(1,5)=  4, 16,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
  13, 14, 15,
 mmnn(1,6)=  5, 16,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
  13, 14, 15,
 mmnn(1,7)=  6, 15,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
  13, 14,
 mmnn(1,8)=  7, 15,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
  13, 14,
 mmnn(1,9)=  8, 15,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
  13, 14,
 mmnn(1,10)=  9, 14,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
  13,
 mmnn(1,11)= 10, 13,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
 mmnn(1,12)= 11, 11,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
 mmnn(1,13)= 12, 10,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
 mmnn(1,14)= 13,  6,  4,  6,  7,  8,  9, 10,
 mmnn(1,15)= 14,  3,  6,  7,  8,
 mmnn(1,16)= 15,  1,  6,
 mmnn(1,17)= -1, 10,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
 mmnn(1,18)= -2,  8,  1,  2,  3,  4,  5,  6,  7,  8,
 mmnn(1,19)= -3,  8,  1,  2,  3,  4,  5,  6,  7,  8,
 mmnn(1,20)= -4,  6,  1,  2,  3,  4,  5,  6,
 mmnn(1,21)= -5,  5,  1,  2,  3,  4,  5,
 mmnn(1,22)= -6,  4,  1,  2,  3,  4,
 mmnn(1,23)= -7,  3,  1,  2,  3,
 mmnn(1,24)= -8,  3,  1,  2,  3,
 mmnn(1,25)= -9,  2,  1,  2,
 mmnn(1,26)=-10,  2,  1,  2,
/
```

This data can be generated from a 'boozmn' Boozer Transformation file.
The user can also provide DKES with such a file to ease implementation.

------------------------------------------------------------------------

### Execution

The simplest invokation of DKES is as follow
```
> xdkes input_dkes.example
```
where the file `input_dkes.example` contains the `DKES_INDATA` namelist
for that run.

Alternatively the user can invoke DKES as follows
```
> xdkes input.example 5 1.0E-3 1.0
```
where the files `wout_example.nc` (or `wout.example`) and `boozmn_example.nc`
are present. The next four numbers are the boozer radial surface,
collisionality (nu/v in m^-1), and radial electric field (E_s/v).
The user may provide a screen output flag after those values (T/F).
Additionally, the user may provide 3 additional command line arguments:
a filename modifer which is appended to the output files, the order of
the Legendre polynominals used in the code (100), and the number of iterations
for coupling (4). This invokation will automatically produce an `input_dkes.example`
file with the filename modifer appended if supplied.

------------------------------------------------------------------------

### Output Data Format

The code will produce three files 'dkesout', 'opt_dkes', and 'results'.
Additionally, if called using the alternative invokation an 'input_dkes'
file will be generated. The 'results' file contains a table of 19 values
with a row for each collisionality/efield pair. The 'opt_dkes' file
contains a simplified output used in the old version of STELLOPT. The
'dkesout' file contains a more detailed and nuanced version of the
output in 'results.'

------------------------------------------------------------------------

### Visualization

The 'results' file may be read by any plotting package capable of
importing space delimited data. Note the columns for `LXXm` and `LXXp`
are the variational bounds for the coefficients.

------------------------------------------------------------------------

### Tutorials

Put links to tutorial pages here.
