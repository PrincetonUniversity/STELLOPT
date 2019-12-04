VMEC Input Namelist (v8.47)
===========================

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
(denoted 1 if it's just a scalar value).

| Variable | Type | Size | Description |
|---|---| --- | --- | ---- |
 | nfp | integer | 1 | Number of toroidal field periods | 
 | ncurr | integer | 1 | Switch for using (0) flux conservation or (1) prescribed toroidal current |
 | nsin | integer | 1 | Deprecated parameters. |
 |niter | integer | 1 | Total number of iterations | 
 | nstep | integer | 1 | Output interval to screen and threed1 file | 
 | nvacskip | integer | 1 | Interval for full update of vacuum solutions | 
 | mpol | integer | 1 | Poloidal Mode Number (m) | 
 | ntor | integer | 1 | Toroidal Mode Number (n) | 
 | ntheta | integer | 1 | Number of theta grid points (\>=2*mpol+6) | 
 | nzeta | integer | 1 | Number of planes (in zeta) on which mgrid data has been calculated. | 
 | mfilter_fbdy | integer | 1 |  |
 | nfilter_fbdy | integer | 1 |  |
 | time_slice | real(12,100) | 1 | Time slice of equilibria (label) | 
 | curtor | real(12,100) | 1 | Total toroidal curent. \[A\] |
| delt | real(12,100) | 1 | Blending parameter | 
 | ftol | real(12,100) | 1 | Residual tolerance. | 
 | tcon0 | real(12,100) | 1 | Weight factor for constrained force. (defaults to 1.0 for values greater than 1.0) | 
 | gamma | real(12,100) | 1 | Adiabatic index (compressional index) | 
 | phiedge | real(12,100) | 1 | Total enclosed toroidal flux [Wb]. | 
 | spres_ped | real(12,100) | 1 | Value (in s) beyond which pressure profile is flat. (pedestal) | 
 | bloat | real(12,100) | 1 | Bloating factor for pressure and current profiles (bloats the domain not the value). | 
 | pres_scale | real(12,100) | 1 |  |
 | prec2d_threshold | real(12,100) | 1 |  |
 | am | real(12,100) | 0:20 | Mass or Pressure (gamma=0) expansion coefficients | 
 | ai | real(12,100) | 0:20 | Iota expansion coefficients | 
 | ac | real(12,100) | 0:20 | Toroidal current density expansion coefficients | 
 | aphi | real(12,100) | 0:20 | Deprecated parameters. | 
 | ns_array | integer | 100 | Number of radial grid points for each grid iteration | 
 | niter_array | integer | 100 | Maximum number of iterations for a given radial resolution | 
 | ftol_array | real(12,100) | 100 | Array of residual values at which a given multigrid iteration ends. | 
 | extcur | real(12,100) | 500 | Array of currents in each external group for free boundary run. | 
 | raxis | real(12,100) | 0:61 | see raxis_cc | 
 | zaxis | real(12,100) | 0:61 | see zaxis_cs | 
 | raxis_cc | real(12,100) | 0:61 | Radial Fourier Cosine Coefficients of magnetic Axis (R=raxis_cc*cos(-n*zeta)) | 
 | raxis_cs | real(12,100) | 0:61 | Radial Fourier Sine Coefficients of magnetic Axis (R=raxis_cs*sin(-n*zeta)) | 
 | zaxis_cc | real(12,100) | 0:61 | Vertical Fourier Cosine Coefficients of magnetic Axis (Z=zaxis_cc*cos(-n*zeta)) | 
 | zaxis_cs | real(12,100) | 0:61 | Vertical Fourier Sine Coefficients of magnetic Axis (Z=zaxis_cs*sin(-n*zeta)) | 
 | rbc | real(12,100) | -61:61,0:60 | Boundary cosine coefficients for R= cos(m*theta-n*zeta) | 
 | rbs | real(12,100) | -61:61,0:60 | Boundary sine coefficients for R | 
 | zbc | real(12,100) | -61:61,0:60 | Boundary cosine coefficients for Z | 
 | zbs | real(12,100) | -61:61,0:60 | Boundary sine coefficients for Z| 
 | lfreeb | logical | 1 | Switch for free boundary run.| 
 | lrecon | logical | 1 | Switch for reconstruction run.| 
 | lmac | logcial | 1 |  |
 | loldout | logical | 1 | Switch for old output format. This will produce a fort.8 file. | 
 | ledge_dump | logical | 1 | Output edge values to 'FORT.NEDGE0'. | 
 lspectrum_dump | logical | 1 | | 
 | lasym | logical | 1 | Non-stellarator symmetric configuration (T) | 
 | loptim | logical | 1 | Obsolete| 
 | lforbal | logical | 1 |  |
 | lrfp | logical | 1 | Switch for using poloidal flux as radial coordinate. |
 |l_spectrum_dump | 1 | Obsolete | |
 | mgrid_file | character | 100 | Vacuum Green's Function Filename | 
 | precon_type | character | 10 | Preconditioner type |
 | pcurr_type | character | 20 | Current Profile type | 
 | piota_type | character | 20 | Iota Profile type | 
 | pmass_type | character | 20 | Pressure (mass) profile type | 
 | phidiam | real(12,100) | 1 | Diamagnetic toroidal flux [Wb]. | 
 | psa | real(12,100) | 100 |  |
 | pfa | real(12,100 | 100 |  |
 | isa | real(12,100) | 100 |  |
 | ifa | real(12,100) | 100 |  |
 |sigma_current | real(12,100) | 1 | Standard deviation in toroidal current [A]. Standard deviations <0 are interpreted as percent of respective measurement. | 
 | imatch_phiedge | integer | 1 | PHIEDGE switch: 0 use pressure profile, (1) match value, (2) use LIMPOS data (in mgrid), (3) use Ip (fixed boundary only) |
 | iopt_raxis | integer | 1 |  |
 | tensi | real(12,100) | 1 | Spline tension for iota profile | 
 | tensp | real(12,100) | 1 | Spline tension for pressure profile. |
| mseangle_offset | real(12,100) | 1 | Uniform experimental offset of MSE data (calibration offset) | 
 | imse | integer | 1 | Number of Motional Stark Effect (MSE) data points | 
 | isnodes | integer | 1 | Number of iota spline points (computed internally unless specified explicitly) | 
 | rstark | real(12,100) | 100 |  |
 | datastark | real(12,100) | 100 | Pitch angle data from stark measurement | 
 | sigma_stark | real(12,100) | 100 | Standard deviation in MSE data [degrees]. Standard deviations of <0 are interpreted as a percent of the respective measurement. | 
 | itse | integer | 1 | Number of pressure profile data points. | 
 | ipnodes | integer | 1 | Number of pressure spline points (computed internally unless specified explicitly) | 
 | presfac | real(12,100) | 1 | Number by which Thomson scattering data is scaled to get actual pressure. | 
 | pres_offset | real(12,100) | 1 | Uniform arbitrary radial offset of pressure data. | 
 | rthom | real(12,100) | 100 | Radial coordinate Data for Thomson scattering. (lpofr=.true. then in real space, lpofr=.false. then in flux space) | 
 | datathom | real(12,100) | 100 | Pressure data from Thomson scattering. | 
 | sigma_thom | real(12,100) | 100 | Standard deviation for pressure profile data [Pa]. Standard deviations of <0 are interpreted as a percent of the respective measurement. | 
 | sigma_delphid | real(12,100) | 1 | Standard deviation for pressure profile data [Wb]. Standard deviations of <0 are interpreted as a percent of the respective measurement. |
| tensi2 | real(12,100) | 1 | vbl spline tension for iota | 
| fpolyi | real(12,100) | 1 | vbl spline tension form factor (note if tensi!=tensi2 then tension(ith point)=tensi+(tensi2-tensi)*(i/n-1)**fpolyi ) |
| nflxs | integer | 1 | Number of flux loop measurements used in matching.| 
 | indxflx | integer | 100 | Array giving index of flux measurement in iconnect array. | 
 | dsiobt | real(12,100) | 100 | Measured flux loop signals corresponding to the combination of signals in iconnect array. | 
 | sigma_flux | real(12,100) | 100 | Standard deviation for external poloidal flux data [Wb]. Standard deviations of <0 are interpreted as a percent of the respective measurement. | 
 | nbfld | integer | 5 | Number of selected external bfield measurements used in matching. | 
 | indxbfld | integer | 100,5 | Array giving index of bfield measurement used in matching | 
 | bbc | real(12,100) | 100,5 | Measured magnetic field at rbcoil(m,n) zbcoil(m,n) at the orientation br*cos(abcoil) + bz*sin(abcoil) | 
 | sigma_b | real(12,100) | 100,5 | Standard deviation for external magnetic field data [T]. Standard deviations of <0 are interpreted as a percent of the respective measurement. | 
 | lpofr | logical | 1 | Switch for pressure data coordinates (.true. real space, .false. flux space) |
 
