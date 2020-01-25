VMEC Input Namelist (v6.9)
==========================

------------------------------------------------------------------------

The following table contains explanations for each variable in the input
namelist used by VMEC. Not all values need be specified in a file.
Examples can be found to guide the user in creating an input namelist
for their application. The file containing the input namelist is a text
file named input.run where 'run' is specified by the user and appended
to all output files for a given run. The 'Type' column indicates the
form of the data. Values denoted real(12,100) are reals defined by
selected_real_kind(12,100), (12 point percision, and 100 point range).
The 'Size' column indicates the size of the variable if it's an array
(denoted 1 if it's just a scalar value). | Variable | Type |
Size | Description | \n | mgrid_file | character | 100 |
Vacuum Green's Function Filename | \n | time_slice |
real(12,100) | 1 | Label for shot (used for output). | \n | nfp
| integer | 1 | Number of toroidal field periods | \n |
ncurr | integer | 1 | Switch for using (0) flux conservation or
(1) prescribed toroidal current | \n | nsin | integer | 1 |
Deprecated parameters. | \n | lasym | logical | 1 |
Non-stellarator symmetric configuration (T) | \n | niter | integer
| 1 | Total number of iterations | \n | nstep | integer |
1 | Output interval to screen and threed1 file | \n | nvacskip
| integer | 1 | Interval for full update of vacuum solutions
| \n | delt | real(12,100) | 1 | Blending parameter |
| ftol | real(12,100) | 1 | Array of value of residual(s) at
which each multigrid iteration ends | \n | gamma | real(12,100)
| 1 | Adiabatic index (compressional index) | \n | am |
real(12,100) | 0:10 | Mass or Pressure (gamma=0) expansion
coefficients | \n | ai | real(12,100) | 0:10 | Iota
expansion coefficients | \n | ac | real(12,100) | 0:10 |
Toroidal current density expansion coefficients | \n | aphi |
real(12,100) | 0:10 | Deprecated parameters. | \n | rbc |
real(12,100) | -61:61,0:60 | Boundary cosine coefficients for R
cos(m*theta-n*zeta) | \n | rbs | real(12,100) | -61:61,0:60
| Boundary sine coefficients for R | \n | zbc | real(12,100)
| -61:61,0:60 | Boundary cosine coefficients for Z | \n | zbs
| real(12,100) | -61:61,0:60 | Boundary sine coefficients for Z
| \n | raxis | real(12,100) | 0:61,2 | Radial position of
magnetic axis (R=raxis*cos(-n*zeta)) | \n | zaxis | real(12,100)
| 0:61,2 | Vertical position of magnetic axis | \n | spres_ped
| real(12,100) | 1 | Value (in s) beyond which pressure profile
is flat. (pedestal) | \n | mpol | integer | 1 | Poloidal
Mode Number (m) | \n | ntor | integer | 1 | Toroidal Mode
Number (n) | \n | ntheta | integer | 1 | Number of theta
grid points (\>=2*mpol+6) | \n | nzeta | integer | 1 |
Number of planes (in zeta) on which mgrid data has been calculated. |
| psa | real(12,100) | 100 | \n | | pfa | real(12,100
| 100 | \n | | isa | real(12,100) | 100 | \n | |
ifa | real(12,100) | 100 | \n | | ns_array | integer
| 100 | Number of radial grid points for each grid iteration |
| ftol_array | real(12,100) | 100 | Array of residual
values at which a given multigrid iteration ends. | \n | tcon0 |
real(12,100) | 1 | Weight factor for constrained force. (defaults
to 1.0 for values greater than 1.0) | \n | curtor | real(12,100)
| 1 | Edge value for toroidal current. | \n | sigma_current
| real(12,100) | 1 | Standard deviation in toroidal current
\[A\]. Standard deviations \<0 are interpreted as percent of respective
measurement. | \n | extcur | real(12,100) | 500 | Array of
currents in each external group for free boundary run. | \n | phiedge
| real(12,100) | 1 | Total enclosed toroidal flux \[Wb\]. |
| imatch_phiedge | integer | 1 | PHIEDGE switch: 0 use
pressure profile, (1) match value, 2 use LIMPOS data (in mgrid), 3 use
Ip (fixed boundary only) | \n | iopt_raxis | integer | 1 |
| \n | tensi | real(12,100) | 1 | Spline tension for iota
profile | \n | tensp | real(12,100) | 1 | Spline tension for
pressure profile. | \n | mseangle_offset | real(12,100) | 1
| Uniform experimental offset of MSE data (calibration offset) |
| imse | integer | 1 | Number of Motional Stark Effect (MSE)
data points | \n | isnodes | integer | 1 | Number of iota
spline points (computed internally unless specified explicitly) |
| rstark | real(12,100) | 100 | \n | | datastark |
real(12,100) | 100 | Pitch angle data from stark measurement |
| sigma_stark | real(12,100) | 100 | Standard deviation in
MSE data \[degrees\]. Standard deviations of \<0 are interpreted as a
percent of the respective measurement. | \n | itse | integer |
1 | Number of pressure profile data points. | \n | ipnodes |
integer | 1 | Number of pressure spline points (computed
internally unless specified explicitly) | \n | presfac |
real(12,100) | 1 | Number by which Thomson scattering data is
scaled to get actual pressure. | \n | pres_offset | real(12,100)
| 1 | Uniform arbitrary radial offset of pressure data. | \n |
rthom | real(12,100) | 100 | Radial coordinate Data for Thomson
scattering. (lpofr=.true. then in real space, lpofr=.false. then in flux
space) | \n | datathom | real(12,100) | 100 | Pressure data
from Thomson scattering. | \n | sigma_thom | real(12,100) |
100 | Standard deviation for pressure profile data \[Pa\]. Standard
deviations of \<0 are interpreted as a percent of the respective
measurement. | \n | phidiam | real(12,100) | 1 | Diamagnetic
toroidal flux \[Wb\]. | \n | sigma_delphid | real(12,100) | 1
| Standard deviation for pressure profile data \[Wb\]. Standard
deviations of \<0 are interpreted as a percent of the respective
measurement. | \n | tensi2 | real(12,100) | 1 | vbl spline
tension for iota | \n | fpolyi | real(12,100) | 1 | vbl
spline tension form factor (note if tensi!=tensi2 then tension(ith
point)=tensi+(tensi2-tensi)*(i/n-1)**fpolyi ) | \n | nflxs |
integer | 1 | Number of flux loop measurements used in matching.
| \n | indxflx | integer | 100 | Array giving index of flux
measurement in iconnect array. | \n | dsiobt | real(12,100) |
100 | Measured flux loop signals corresponding to the combination of
signals in iconnect array. | \n | sigma_flux | real(12,100) |
100 | Standard deviation for external poloidal flux data \[Wb\].
Standard deviations of \<0 are interpreted as a percent of the
respective measurement. | \n | nbfld | integer | 5 | Number
of selected external bfield measurements used in matching. | \n |
indxbfld | integer | 100,5 | Array giving index of bfield
measurement used in matching | \n | bbc | real(12,100) | 100,5
| Measured magnetic field at rbcoil(m,n) zbcoil(m,n) at the
orientation br*cos(abcoil) + bz*sin(abcoil) | \n | sigma_b |
real(12,100) | 100,5 | Standard deviation for external magnetic
field data \[T\]. Standard deviations of \<0 are interpreted as a
percent of the respective measurement. | \n | lpofr | logical |
1 | Switch for pressure data coordinates (.true. real space, .false.
flux space) | \n | lfreeb | logical | 1 | Switch for free
boundary run. | \n | lrecon | logical | 1 | Switch for
reconstruction run. | \n | lmca | logcial | 1 | \n | |
loldout | logical | 1 | Switch for old output format. | \n |
ldiagno | logical | 1 | Switch for production of the
'diagno_in.' file used for reconstruction with the DIAGNO routine.
| \n | ledge_dump | logical | 1 | Output edge values to
'FORT.NEDGE0'. | \n | lspectrum_dump | logical | 1 | \n |
| loptim | logical | 1 | \n | | bloat | real(12,100)
| 1 | Bloating factor for pressure and current profiles (bloats
the domain not the value). | \n | ac_form | integer | 1 |
Determines form of toroidal current (if NCURR=1). Expansion in toroidal
flux (1) otherwise include square root term. |**
