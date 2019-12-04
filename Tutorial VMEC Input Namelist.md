Tutorial: VMEC Input Namelist
============================================================================================================

This tutorial is designed to help the new user to VMEC understand the
VMEC input namelist and what options are available.

The VMEC code is controlled through a FORTRAN input namelist called
(INDATA)[VMEC Input Namelist (v8.47)]. A FORTRAN input namelist is simply a text file where variable
values may be set in a straightforward fashion. They may then be read
into a program to initialize various variables. It is important to note
that while variables in a given namelist may be omitted from a file the
file must not contain unknown variable names. The VMEC code looks for a
file (in the current directory) by the name of 'input.ext' where
'ext' is the extension passed to it as a command line option. This
file should contain an input namelist by the name of 'INDATA.'
Multiple input namelists may exist in a text file. Each list is denoted
by a '&' followed by the input namelist name. Each namelist is
terminated by the '/' (forward-slash) character. There are various
examples in the tutorial section.

The VMEC input namelist can be broken down into sections according to
what the different variables control. For the rest of this tutorial we
will examine each section giving detailed descriptions of how the choice
of variables affects the execution of the code. A full version of this
file can be found here:

The namelist declaration and Runtime control parameters The VMEC input
namelist is declared 'INDATA' and is the first line of the namelist.
Various parameters control the execution of the VMEC code. The DELT
parameter controls the evolution of the VMEC solution (in a sense the
stepsize between equilibria). The NITER parameter determines the maximum
numer of iterations for a given solution at a given radial resolution.
The solver will actually run for twice this value if it has not
converged, and print out a message to this effect. The NSTEP parameter
determines the number of iterations between outputs of the progress
towards convergence. The TCON0 parameter controls the constrained force
calculation. It's value is set to 1.0 for any value greater than 1.0.
The NS_ARRAY determines the number of radial grid points to use for
each equilibrium calculation. The FTOL_ARRAY parameter determines the
cutoff in the normalized force residuals for each radial grid. Once a
given value of FTOL_ARRAY is reached at a given NS_ARRAY radial mesh
the solver moves on to the next radial grid. The NITER_ARRAY determines
the number of iterations at a given radial resolution. If the value in
the corresponding FTOL_ARRAY has not been meet in NITER_ARRAY steps,
the code moves to the next radial grid. The LWOUTTXT logical variable
controls how the output is saved. If set to FALSE (default) the wout
file will be output in netCDF. If set TRUE then a text file will be
generated. If compiled without netCDF then code will output text wout
files.

```Fortran
&INDATA
DELT = 9.00000000000000E-01
NITER = 10000
NSTEP = 200
TCON0 = 1.00000000000000E+00
NS_ARRAY = 9 29 49 99
FTOL_ARRAY = 1.00000000000000E-08 1.00000000000000E-10 1.00000000000000E-12 1.00000E-14
NITER_ARRAY= 2000 3000 5000 10000
LWOUTTXT = F
```

The Grid Parameters The VMEC computational domain is controlled through
these parameters. The LASYM parmeter (T/F) determines if up/down
symmetry is to be violated. This value is defaulted to false
(stellarator symmetry). The NFP parameter controls the periodicity of
the simulation. This allows for a significant reduction in the number of
toroidal modes. The MPOL parameter controls the total number of poloidal
modes the simulation uses (m=0\...MPOL-1). The NTOR parameter controls
the total number of toroidal modes (n=-NTOR\...NTOR). Internally the NFP
variable and n are combined allowing for an efficient reduction in
computational effort. The PHIEDGE parameter controls the total enclosed
toroidal flux. In essence this controls the total volume of the plasma,
by scaling the choice of boundary coefficients to match the value found
in PHIEDGE. [code format=\"fortran\"](code format="fortran") LASYM = F
NFP = 10 MPOL = 8 NTOR = 6 PHIEDGE = 8.28000000000000E-01 [code](code)

The Free Boundary Parameters To preform a run in free boundary mode VMEC
must be supplied various parameters. The LFREEB parameter (T/F)
indicates if the code should be executed in free boundary mode. The
MGRID_FILE parameter indicates the location of the 'mgrid' file which
contains the vacuum magnetic field on a grid in R,Z and PHI. The NTHETA
parameter determines the number of points in theta to use to represent
the VMEC flux surfaces. This value defaults to 2*(MPOL)+6, and in the
strict sense is not just a free boundary parameters. The NZETA parameter
determines the number of gridpoints in zeta/phi to use. This value must
be equal to the number of phi planes in the 'mgrid' file, for free
boundary runs. The default value is 2*NTOR+4 unless NTOR=0 then it
defaults to 1. The EXTCUR parameter is an array specifying the current
in each current group found in the 'mgrid' file. The NVACSKIP
parameter determines the number of steps between update of the vacuum
solution.

```Fortran
LFREEB = T
MGRID_FILE = '/home/usr/coils/mgrid.lhd_ys4.msize'
NTHETA = 22
NZETA = 32
NVACSKIP = 6
EXTCUR
= 4.82137988281250E+03 3.14672241210938E+03 2.13679345703125E+03 -5.03201074218750E+03 -3.84028015136719E+02 3.77810131835938E+03
```

The Pressure Profile Parameters The pressure profile in VMEC is
specified as a function of radial flux space coordinates (s). Where for
the default case is the normalized toroidal flux (normalized to
PHIEDGE). If LRFP is set to TRUE, this coordinates becomes the
normalized poloidal flux. The AM parameter determines the polynomial
coefficients (0..10) used to calculate pressure profile: \$$ p=\sum_{n=0}^{10} am(n) * s^n . $$
The GAMMA parameter controls this profile by scaling it. In general the user should choose
it's value to be 0.0. The BLOAT parameter acts as a scaling factor to
s. Again it's value should be chosen to be 1.0. The SPRES_PED
parameter is that value of s at which the pressure is flat. This allows
modeling of pedestal pressure profiles. In general it should be set to
1.0 equating to the boundary in normalized coordinates. The PRES_SCALE
value is a scale factor applied to the profile allowing it to be scaled
up and down. Note the current version (8.47) provides the user with more
detailed ways of specifying the pressure profiles. The user is encourage
to examine the profile_functions.f file for details (see
[VMEC Advanced Profiles](VMEC Advanced Profiles)).

```Fortran
GAMMA = 0.00000000000000E+00
BLOAT = 1.00000000000000E+00
SPRES_PED = 1.00000000000000E+00
PRES_SCALE = 1.0
AM = 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00
```

The Toroidal Current / Rotational Transform Parameters The VMEC code
provides the user an option to specify either a rotational transform
radial profile or a toroidal current density profile. The NCURR
parameter determine which form of the profile to use (0: Rotational
Transform, 1: Toroidal Current Density). The AI parameter specifies the
polynomial coefficients (0..10) used to calculate the rotational
transform profile (NCURR=0) $$ \iota=\sum_{n=0}^{10} ai(n)
* s^n . $$ The AC_FORM parameter determines the form of the
current profiles used (NCURR=1). For AC_FORM=0 the toroidal current
profile is power series in s defined by the AC parameter
\$$ j=\sum_{n=0}^{10} ac(n) * s^n . $$

```Fortran
NCURR = 1
CURTOR = 0.00000000000000E+00
AI = 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00
AC = 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00
```

The Magnetic Axis Parameters The VMEC code needs an initial guess for
the magnetic axis. These values are specified as Fourier harmonics in
toroidal mode number (n=0..NTOR). The RAXIS parameter stores the cosine
then sine harmonics while the ZAXIS parameter stores the sine then
cosine harmonics.

```Fortran
RAXIS = 3.80000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00
ZAXIS = 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00
```

The Boundary Shape Parameters The VMEC boundary is specified in terms of
Fourier harmonics. The RBC parameter store the radial Fourier
coefficients while the ZBS parameter stores the vertical Fourier
coefficients. The RBS and ZBC parameters are provided for axisymmetric
runs.

```Fortran
RBC( 0,0) = 3.80
ZBS( 0,0) = 0.00
RBC( 0,1) = 0.30
ZBS( 0,1) = 0.60
RBC(-1,1) = -0.15
ZBS(-1,1) = 0.30
```
