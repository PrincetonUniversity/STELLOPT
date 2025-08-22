Preparing FIDASIM simulations with BEAMS3D
=============================================

This tutorial will walk the user through running the BEAMS3D code and generate the necessary output files for [FIDASIM](https://d3denergetic.github.io/FIDASIM/index.html).

------------------------------------------------------------------------
Table of Contents
   * [Distribution](#dist)
   * [Input](#namelist)
   * [Output](#output)

------------------------------------------------------------------------

### Distribution

BEAMS3D can calculate the steady-state slowing down distribution, which can be used by FIDASIM to generate synthetic FIDA spectra which are directly comparable to experimental data and thus serve as a valuable tool for validation and physics studies. 
To generate the FIDASIM files, the code has to be run with either the `-fidasim` or `-fidasim_cyl` flag.  The following explanation is based on the [BEAMS3D validation paper](https://doi.org/10.1088/1741-4326/adeda2).

The native distribution function in BEAMS3D is collected as a 6-dimensional array, where the dimensions are: number of beams and beam energies, radial bins (units of toroidal flux), poloidal bins (flux aligned), toroidal bins, parallel velocity bins, and perpendicular velocity bins. The resulting distribution is normalized to the real $$[\mathrm{m}^{-3}]$$ and phase-space volumes $$[(\mathrm{m^3/s^3})^{-1}]$$, thus giving units of $$[\mathrm{s^3/m^6}]$$. 
FIDASIM requires the distribution to be given on a cylindrical real-space grid (as opposed to the flux-aligned bins of the BEAMS3D distribution) as well as energy and pitch phase-space coordinates, explicitly in units $$[\mathrm{1/cm^3/keV}]$$. As such, a conversion both in real and phase space is necessary. The phase space conversion is achieved when running with the `-fidasim` flag by

  $$f_{\mathrm{FIDASIM}}(E,p) = 2\pi \frac{e \cdot v}{10^3 \cdot m} \cdot f_{\mathrm{B3D}}(v_{\parallel},v_{\perp}),$$

where $$e$$ is the elementary charge in [C], $$v$$ is the total velocity in [m/s] and $$m$$ is the mass of the fast ion species in [kg]. This has been implemented in BEAMS3D using nearest-neighbor binning to avoid a memory-intensive interpolation. 

Alternatively, the distribution function for FIDASIM can also be collected directly alongside the regular distribution (using the `-fidasim_cyl` flag), resulting in increased memory usage during simulations but significantly increased precision of the distribution, while not significantly affecting overall run time. There are many grid cells in the cylindrical real-space grid outside the plasma, which have to be allocated in memory but are largely free of any fast ions. This memory inefficiency is the main reason why BEAMS3D usually collects the distribution on a flux-aligned grid. However, for the application in the FIDASIM context the increased accuracy by skipping the renormalization and the independence from conversion to flux coordinates makes this approach attractive. This is especially true for cases where continuous flux coordinates are difficult or impossible to construct, as in cases where magnetic islands are included in the equilibrium.

### Input

The additional input namelist in the BEAMS3D namelist file has to contain all relevant information for generating the FIDASIM files. Detailed information on the inputs can be found in the [FIDASIM documentation](https://d3denergetic.github.io/FIDASIM/page/03_technical/01_prefida_inputs.html).

```fortran
&fidasim_inputs_b3d
    ab = 1.0
    adist = 586.108695237747, 667.7809604002308
    ai = 1.0
    alpha = -1.6792548743640188
    aoffy = 19.48, -18.2
    aoffz = 9.102, 0.186
    ashape = 1, 1
    awidy = 30.0, 30.0
    awidz = 40.0, 40.0
    axis_nbi = -0.107, -0.99, -0.085
    axis_spec(:,1) = 0.63102, 0.43877, -0.63976
    axis_spec(:,2) = -0.99542, -0.079458, -0.053215
    axis_spec(:,3) = 0.58196, 0.43525, -0.68694
    axis_spec(:,4) = -0.9939, -0.095986, -0.054339
    beta = 0.08527719127138363
    calc_bes = 1
    calc_birth = 0
    calc_brems = 1
    calc_cold = 1
    calc_dcx = 1
    calc_fida = 1
    calc_fida_wght = 1
    calc_halo = 1
    calc_neutron = 0
    calc_npa = 0
    calc_npa_wght = 0
    calc_pfida = 1
    calc_pnpa = 0
    comment = 'BEAMS3D Distribution function'
    current_fractions = 0.5, 0.3, 0.2
    device = 'Test device'
    divy = 0.0099, 0.010, 0.012
    divz = 0.0099, 0.010, 0.012
    einj = 55.0
    emax_wght = 65.0
    focy = 650.0
    focz = 850.0
    gamma = 0.0
    id = 'Channel 1', 'Channel 2', 'Channel 3', 'Channel 4'
    impurity_charge = 5
    lambdamax = 669.0
    lambdamax_wght = 665.0
    lambdamin = 645.0
    lambdamin_wght = 647.0
    lens(:,1) = -52., 541., 106.
    lens(:,2) = 210., 610., 37.
    lens(:,3) = -52., 541., 105.
    lens(:,4) = 210., 610., 37.
    n_birth = 10000
    n_dcx = 500000
    n_fida = 5000000
    n_halo = 50000
    n_nbi = 50000
    n_npa = 5000000
    n_pfida = 50000000
    n_pnpa = 50000000
    name = 'W7X K21:7'
    naperture = 2
    nbi_data_source = 'generated'
    nchan = 48
    ne_wght = 10
    nlambda = 1024
    nlambda_wght = 10
    np_wght = 10
    nphi_wght = 8
    nx = 90
    ny = 90
    nz = 60
    origin = 33., 669., 34.
    pinj = 1.7
    radius = 595., 596., 594., 593.
    result_dir = '' !for current directory
    runid = 'fidasim_runid'   !! runID
    shape = 1
    shot = 999999
    sigma_pi = 0.9, 0.9, 0.9, 0.9
    spec_data_source = 'matlab'
    spot_size = 0.0, 0.0, 0.0, 0.0
    src = 104., 1321., 90.
    system = 'FIDA system'
    tables_file = 'atomic_tables.h5'
    time = 1.0
    widy = 11.0
    widz = 25.0
    xmax = 180.0
    xmin = 0.0
    ymax = 180.0
    ymin = -80.0
    zmax = 60.0
    zmin = -60.0
/
```
------------------------------------------------------------------------

### Output

The resulting output files are

- fidasim_runid_distribution.h5
- fidasim_runid_geometry.h5
- fidasim_runid_equilibrium.h5
- fidasim_runid_inputs.dat

These can be directly used to run FIDASIM as outlined in its [documentation](https://d3denergetic.github.io/FIDASIM/page/01_getting_started/03_running.html).

The outputs can be analyzed using the functions included in [matlabVMEC](https://github.com/kudav/matlabVMEC/tree/Kulla/FIDASIM), as well as the python interface included in FIDASIM.