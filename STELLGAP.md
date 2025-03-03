STELLGAP
======

Stellgap calculates the shear Alfvén gap structure for 3D configurations (stellarators, RFPs, 3D tokamaks)

These codes are used to calculate shear Alfven continua for 3D configurations, both with and without sound wave coupling effects. The associated paper is [D. A. Spong, R. Sanchez, A. Weller, "Shear Alfvén continua in stellarators," Phys. Plasmas 10 (2003) 3217–3224.](https://doi.org/10.1063/1.1590316)

------------------------------------------------------------------------

### Theory

The Alfvén continuum equation for 3D stellarator equilibira in low
plasma beta, incompressible limite is written:

$$\mu_0\rho\omega^2\frac{|\nabla\psi|^2}{B^2}E_\psi+\vec{B}\cdot\nabla\left[\frac{|\nabla\psi|^2}{B^2}\left(\vec{B}\cdot\nabla\right)E_\psi\right]=0$$

The STELLGAP code reformulates this equation in terms of a eigenvalue
equation which is sovled using the Lapack routine [DGGEV](https://netlib.org/lapack/explore-html-3.6.1/d9/d8e/group__double_g_eeigen_gab3b93851a33e592f5705fdcbc876c186.html) routine.

The following table provides a list of Alfvén couplings which are of
interest to 3D equilibria.

| Abbreviation | Name | $\delta_m$ | $\delta_n$ |
| :---: | :---: | :---: | :---: |
| GAE | Global Alfvén eigenmode | 0 | 0 |
| TAE | Toroidal Alfvén eigenmode | ±1 | 0 |
| EAE | Elliptical Alfvén eigenmode | ±2 | 0 |
| NAE | Noncircular Alfvén eigenmode | >2 | 0 |
| MAE | Mirror Alfvén eigenmode | 0 | ±1, ±2, ... |
| HAE | Helical Alfvén eigenmode | >0 | ±1, ±2, ... | 


------------------------------------------------------------------------

### Compilation

The STELLGAP code can be found in the ORNL repository [Github:STELLGAP](https://github.com/ORNL-Fusion/Stellgap). Once you download the code there
are shell scripts for building the code, otherwise you can use the
following makefile:

```makefile

PRECOMP=gcc -traditional-cpp -cpp 
FC=mpif90 -g -fbacktrace -fexternal-blas
BUILD_TYPE = -DSERIAL

LIBSTELL_DIR = $(STELLOPT_PATH)LIBSTELL/Release
LIBSTELL_LIB = ~/bin/libstell.a

FFLAGS= -I${LIBSTELL_DIR} $(shell nf-config --fflags)
LDFLAGS=${LIBSTELL_LIB} $(shell nf-config --flibs) -framework Accelerate 


OBJ=fitpack.o Fourier_lib_convolve.o 

%.o: %.f
	$(PRECOMP) -E -P $(BUILD_TYPE) $^ > temp.f
	$(FC) -c -o $@ temp.f $(FFLAGS)

xmetric: metric_element_create_ver8.46.o
	$(FC) -o $@ $^ $(LDFLAGS)

xstgap: stellgap_ver5.o $(OBJ) 
	$(FC) -o $@ $^ $(LDFLAGS)

xstgap_snd_ver6: stellgap_soundwave_lagrng_ver6.o $(OBJ) 
	$(FC) -o $@ $^ $(LDFLAGS)

xstgap_snd_ver7: stellgap_soundwave_lagrng_ver7.o $(OBJ) 
	$(FC) -o $@ $^ $(LDFLAGS)

clean:
	@rm -rf *.o xmetric xstgap xstgap_snd_ver6 xstgap_snd_ver7 temp.f
```

Note that the parallel version can be built by changing the line
`BUILD_TYPE` to `-DPARALLEL`. It is important that whatever compiler was
used to build LIBSTELL be used to build STELLGAP.


------------------------------------------------------------------------

### Input Data Format

The STELLGAP code requires that the VMEC equilibrium be transformed into
Boozer Coordinates by the [BOOZ_XFORM] code. The output of that code is
then transformed into input files for STELLGAP. There are two files the
user must provide (`fourier.dat` and `plasma.dat`).

The `plasma.dat` code contains a Fortran input namelist which defines the
density profile for the computation.

```fortran
&PLASMA_INPUT
! ION MASS in AMU
 ion_to_proton_mass = 2.
! Core Ion Density in m^3
ion_density_0 = 2.5e+20
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Ion Profiles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ion profile type iota
! ion_profile = 0 ! [iota(rho)/iota(0)]**2
! Ion profile type polynomial
! ion_profile = 1
! nion = 1.0 -2.0 1.0 ! Similar to AM/AI/AC array
! Ion profile type constant
! ion_profile = 2
! Ion profile type two power
 ion_profile = 3 ! [1 - aion*rho**bion]**cion
 aion = 0.99
 bion = 0.50
 cion = 1.00
!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Unused AE3D variables
!!!!!!!!!!!!!!!!!!!!!!!!!!
 jdqz_data = T
 egnout_form = 'asci'
/
```

The `fourier.dat` file is a text file containing the grid and spectrum
information it has the following format (all quantities are integers).

```
NFP  NTHETA  NZETA  MODE_FAMILY
NMODES
N ML MU
N ML MU
N ML MU
...
```

 * `NFP`: number of field periods
 * `NTHETA`: Number of poloidal grid points (see surface_area_elements)
 * `NZETA`: Number of toroidal grid points (see surface_area_elements)
 * `MODE_FAMILY`: Toroidal mode number about which the eigenfunctions are built (not utilized by code)
 * `NMODES`: Number of toroidal modes used (should be multiple of NFP)
 * `N`: Toroidal mode number
 * `ML`: Starting Poloidal mode number
 * `MU`: Ending Poloidal mode number

The choice of toroidal mode numbers to investigate should follow the
mode coupleing paradigm. Specically `N=n0+NFP*k` and `N=n0-NPF*k`. For
example, and `n0=1` modes for an NFP=5 device would be -11,-9,-6,-4,4,6,9,11 
where both postive and negative `n0` are accounted for.

------------------------------------------------------------------------

### Execution

To execute the STELLGAP code you must start from a [VMEC] equilibrium,
transform the equilibrium to Boozer coordiantes using [BOOZ_XFORM],
process the Boozer coordinates with the `xmetric` routine, generate the
`fourier.dat` and `profiles.dat` files, and finally run a version of
STELLGAP (`xstgap`, `xstgap_snd_ver6`, `xstgap_snd_ver7`).

1. Run the [VMEC] code to produce a `wout` file.
2. Generate the `in_booz` file for your VMEC for all radial gridpoints in your equilibrium.
3. Run the `xmetric` utility to generate input files. This routine takes the boozmn file suffix as input on the command line.
4. Generate the `fourier.dat` and `profiles.dat` files for your run.
5. Run STELLGAP (`xstgap`, `xstgap_snd_ver6`, `xstgap_snd_ver7`). STELLGAP takes two values command line input (`irad` and `ir_fine_scl`)
 * `irad`: should be set to the number of flux surface in original file - 2.
 * `ir_fine_scl`: Number of surfaces in the output should be greater than `irad`. 

------------------------------------------------------------------------

### Output Data Format

The following files are produced by the code.

 * XMETRIC
   * ae_metric.dat : First line is `(nsd-2), izeta, itheta, nznt`, then followed by profiles and finally various metric elements.
   * cart_coords : X,Y,Z data.
   * full_torus_coords : X,Y,Z over entire torus.
   * metrics : Metric elements
   * surf_area_elements : First line is `nfp, izeta, itheta`, then zeta, theta, and surface area.
   * tae_data_boozer : Each block starts with `ks,iotac,phipc,jtorc,jpolc` followed by metric information.
   * torus_slice_coords : X,Y,Z data.
 * XSTELGAP
   * alfven_spec : rho, alfr, alfi, beta, m, n
   * coef_array : rm, rn, f1, f3a, f3b, f3c
   * data_post : iopt, mn, if_fine_scale, isym_opt
   * ion_profile : rho, n_ion, iota, va
   * modes : Table of modes


------------------------------------------------------------------------

### Visualization

STELLGAP provides two routines `post_process.f` and `post_process_snd.f`
to aide in plotting values from `alfven_spec`. These routine process the
`alfven_spec` file defining lam_r = alfr/beta and lam_i = alfi/beta. Then
the complex number omega2 = (lam_r,lam_i) is defined. Finally, omega_r
is computed as the real part of the square root of omega2. 

------------------------------------------------------------------------

### Tutorials

TBD

------------------------------------------------------------------------

### References

 - [D. A. Spong, R. Sanchez, A. Weller, "Shear Alfvén continua in stellarators," Phys. Plasmas 10 (2003) 3217–3224.](https://doi.org/10.1063/1.1590316)
 - [A.G. Ghiozzi et al., "Modeling of frequency-sweeping Alfvén modes in the TJ-II stellarator," Nucl. Fusion 64 036005 (2004)](https://doi.org/10.1088/1741-4326/ad1c93)
