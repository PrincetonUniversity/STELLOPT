Tutorial: Using STELLOPT to minimize turbulent transport
========================================================

The STELLOPT code can optimize equilibria for reduced turbulent
transport using either analytic proxy functions
\<ref\>[H. Mynick, P. Xanthopoulos, B. Faber, M. Lucia, M. Rorvig, and J. Talmadge \"Turbulent optimization of toroidal configurations.\" Plasma Phys. Control. Fusion 56 (2014)](http://dx.doi.org/10.1088/0741-3335/56/9/094001)\</ref\>
or linear growth rates as calculated by the [GENE](http://genecode.org/)
code. This page describes the various proxy functions available to the
user and how to compile STELLOPT with the [GENE](http://genecode.org/)
code. It should also be noted that this feature of STELLOPT replicates
the functionality of the
[GIST](http://www2.ipp.mpg.de/~pax/GIST/gist.html) code for interfacing
the VMEC equilibria to GENE.

------------------------------------------------------------------------

Turbulent transport proxy functions
-----------------------------------

The STELLOPT code can optimize equilibria for reduced turbulent
transport using various proxy functions. These functions do rely on a
direct calculation of turbulence, instead they return quantities which
are related to the growth rates. This results in much simpler (faster)
calculations. The idea being that if the quantities which influence
growth rates are minimized, so will the growth rates. The following
table lists the available proxies \|\| Name \|\| Description \|\| \|\|
prox1 \|\| Ion Temperature Gradient (ITG) proxy \|\| \|\| prox1b \|\|
Ion Temperature Gradient (ITG) proxy \|\| \|\| prox1c \|\| Ion
Temperature Gradient (ITG) proxy \|\| \|\| prox1d \|\| Ion Temperature
Gradient (ITG) proxy \|\| \|\| prox1e \|\| Ion Temperature Gradient
(ITG) proxy \|\| \|\| prox1f \|\| Ion Temperature Gradient (ITG) proxy
\|\| \|\| prox1g \|\| Ion Temperature Gradient (ITG) proxy \|\| \|\|
tem\_overlap \|\| Trapped Electron Mode (TEM) proxy (feedback on
curvature) \|\| \|\| tem\_bounce \|\| Trapped Electron Mode (TEM) proxy
(feedback on curvature and magnetic well) \|\| \|\| tem\_bounce\_tau
\|\| Trapped Electron Mode (TEM) proxy (feedback on curv. and mag well.
weighted) \|\| \|\| gene \|\| Runs the GENE code in linear mode in
serial mode (see next section) \|\| \|\| gene\_parallel \|\| Runs the
GENE code in linear mode in parallel mode (see next section) \|\| At
each radial location you wish to evaluate the turbulent transport proxy
you\'ll need to set a target (TARGET\_TXPORT), sigma (SIGMA\_TXPORT),
and normalized toroidal flux value (S\_TXPORT). Additional name lists
for the GIST code will also need to be included in the input file an
abridged version of the input file is shown below

    &INDATA
    ......
    /
    &OPTIMUM
    ......
      TXPORT_PROXY = 'PROX1D'
      TXPORT_NZ = 128 ! Must match nz0 in parameters file (for GENE)
      LTXPORT_GLOBAL = F
      NALPHA = 1 ! Number of flux tubes
      ALPHA_START = -0.6283 ! Start of flux tube
      ALPHA_END   =  0.6283 ! END of flux tube
      TARGET_TXPORT(01) = 0.0E+00  SIGMA_TXPORT(01) = 1.0E-03  S_TXPORT(01) = 0.25
      TARGET_TXPORT(02) = 0.0E+00  SIGMA_TXPORT(02) = 1.0E-03  S_TXPORT(02) = 0.50
      TARGET_TXPORT(03) = 0.0E+00  SIGMA_TXPORT(03) = 1.0E-03  S_TXPORT(03) = 0.66
    ......
    /

Note that the STELLOPT version of GIST only supports PEST coordinates.
While the full version of the GIST input namelists can be read in,
STELLOPT overrides many of the namelist parameters (for example S0 in
the &SETUP namelist will be over-ridden by the values in S\_TXPORT).

------------------------------------------------------------------------

Compiling STELLOPT with GENE
----------------------------

Please note that the STELLOPT setup script must be modified to properly
link [GENE](http://www2.ipp.mpg.de/~fsj/gene/) into STELLOPT (currently
supported on the PPPL cluster and NERSC-hopper). The following steps
outline what must be done to compile the STELLOPT code with GENE
support.

1.  Obtain and compile the GENE code.
2.  Set an environment variable titles GENE\_PATH which points to the
    \<GENE\>/bin/\<obj\> directory where the .o and .mod files are.
3.  Run the STELLOPT setup script. Please note that you must have TXPORT
    support in your version. If the code properly compiles you should
    see a message indicating the TXPORT was detected along with GENE
    support. This happens before the codes begin compilation.

------------------------------------------------------------------------

Running STELLOPT with GENE
--------------------------

To optimize and equilibrium with the GENE code you will need to set the
TXPORT\_PROXY namelist parameter to \'gene\' or \'gene\_parallel.\' You
will also need a parameters file in your run directory (GENE reads this
file). If you choose to run the GENE code in parallel mode you will need
to specify the NPOPULATION namelist parameter. GENE must be run with
multiples of 2. For example, say you request a run with 40 processors
(ex. mpirun -np 40) then your value for NPOPULATION should be 5.
STELLOPT will then optimize as if 5 processors were requested and use
the additional processors to run GENE (so that each GENE run uses 8
processors, 5x8 = 40). \_\_Please note: n\_procs\_z (&parallelization),
nz0 (&box), npopulation (&optimum), and nz0 (&SETUP) must all be
consistent with the total number of processors (nz0 must be evenly
divisible by n\_procs\_z).\_\_

Some clusters may require that the environment variable
MP\_DEBUG\_NOTIMEOUT=yes to avoid a segmentation fault as some
processors may have to wait at a barrier statement for a few minutes.

Here is an example parameters file:

    &parallelization
    ! Defines number of MPI processes handling each dimension
    ! must evenly divide the resolution number (value)
    n_procs_s =   1  ! Species parallelization (n_spec)
    n_procs_v =   1  ! Vll parallelization (nv0)
    n_procs_w =   1  ! Mu (magnetic moment) parallelization (nw0)
    n_procs_y =   1  ! Poloidal wavenumber parallelization (nky0)
    n_procs_z =   1  ! Field line coordinate parallelization (nz0)
    n_procs_x =   1
    /
    &box
    ! Simulation domain definition
    !    Linear runs will reduce nx0 by 1 if an even number is input set nx0 = 2^n-1 to avoid error message.
    n_spec = 1  ! Number of species
    nx0 =  15  ! Number gridpoints in radial direction (power of two)
    nky0 = 1   ! Number of Fourier modes in poloidal direction
    nz0 =  128 ! Number of grid points along field line (must be even)
    nv0 =  32  ! Number of grid points in parallel direction (must be even)
    nw0 =   8  ! Number of grid points in Mu direction (must be even)
    lx =        125.628  ! Extension of the box in x direction (gyroradii)
    kymin =     0.9 ! Extension of the box in polodial direction (inverse gyroradii)
    lv =        3.0 ! Extension of the parallel velocity direction (thermal velocity)
    lw =        9.0 ! Extension of the Mu direction (T/Bref)
    adapt_lx = .t.  ! Set's the optimal lx =1/(ky*shat) (ignored for nonlinear)
    /
    &in_out
    ! Control of input and output
    diagdir = './'         ! Diagnotic dir (set to ./ to trick GENE into running)
    write_checkpoint = .t. ! Write checkpoint files at end of run
    istep_field =   100    ! Steps between entries in field.dat
    istep_mom =     100    ! Steps between entries in mom.dat
    istep_nrg =      10    ! Steps between entries in nrg.dat
    istep_vsp =     5000   ! Steps between entries in vsp.dat
    istep_schpt =  10000   ! Steps between entries in s_checkpoint
    /
    &general
    ! Control of how the algorithm runs
    comp_type = 'EV'       ! Eigenvalue (EV-linear), Initial Value (IV), or Neoclassical (NC)
    nonlinear = .f.        ! Linear vs. Non-linear runs
    ! Initial Value Runs
    ntimesteps = 50000     ! Maximum number of timesteps
    timelim  = 500        ! Limit computation to this number of seconds
    simtimeline = 490      ! Control stop of code in seconds
    omega_prec =  1.0E-4   ! Desired precision of growth rate and frequency
    overflow_limit = 1.0E30 ! Maximum number of NRG moments
    underflow_limit = 1E-8 ! Minimum value of NRG moments
    dt_max = 0.0385        ! Maximum Timestep
    calc_dt   = .t.        ! Automatic determination of dt_max
    timescheme = 'RK4'     ! RK3, RK4, RK4M, IE1p, IE1s Timeschemes
    coll_split_scheme = 'RKCa' ! EE1, RKC2, RKC3, RKC4, RKCa, none Collision Operatore time scheme
    ! Eigenvalue Runs
    ev_prec = 1.0E-2       ! Desired Precision
    which_ev = 'jd'        ! Desired solver
    n_ev = 1               ! Number of Eignevalues to be computed
    ev_max_it = 500        ! Maximum number of EV iterations
    ev_shift  = (10.0,0.0) ! Shift in the complex plane
    ! Neoclassical
    include_f0_contr = .F. ! Calculated neoclassical contribution
    ! Initialization
    init_cond = 'alm'      ! Initialization of modes (all modes)
    ! Species independent physical parameters
    beta = 0.0             ! Plasma beta
    collision_op = 'none'  ! Collision Operator
    coll_cons_model = 'default' ! Collision back-reaction term
    coll = 1.0E-10         ! Normalized collision frequency
    Zeff = 1.0             ! Z-effective
    debye2 = 0.0           ! Squared debye wavelength normalized to rho_ref
    ! Hyper Diffusion Settings
    ! hyp_XXX : Hyper diffusion coefficient in XXX direction (X,Y,Z,V)
    ! hyp_XXX_order: Exponent of the hyper diffusion operator in a given direction
    hyp_z = 0.25
    hyp_z_order = 4
    hyp_v = 0.2
    hyp_v_order = 4
    ! Other
    bpar = .f.             ! Calculation of parallel magnetic fluctuations (beta >0)
    delzonal = .f.         ! Suppress zonal flows (phi)
    delazonal = .f.        ! Suppress zonal flows (A)

    courant = 0.4
    pressure_off=.f.
    /
    &geometry
    ! Magnetic Geometry Namelist
    magn_geometry = 'gist'  ! Type of data
    ! 'slab','s_alpha','circular'
    q0      = 1.4           ! Central q
    shat    = 0.6565        ! Magnetic Shear
    !amhd    = 0.0           ! Alpha parameter
    !major_r = 1.6           ! Major radius
    !minor_r = 0.3           ! Minor radius
    !trpeps =  0.18          ! Inverse aspect ratio at flux tube
    ! 'miller'
    !kappa = 1.0             ! Elongation
    !delta = 0.0             ! Triangularity
    !zeta = 0.0              ! Squareness
    !s_kappa = 0.0           ! Shear
    !s_delta = 0.0
    !s_zeta  = 0.0
    !drR     = 0.0           ! Major radius shift
    ! 'slab'
    !shat = 0.6565           ! Magnetic shear
    !parscale = 1.0          ! Parallel scale length
    ! 'tracer', 'tracer-efit', 'gist'
    geomdir = '.'           ! Directory of geometry file (STELLOPT controlled)
    geomfile = '.'          ! Name of geometry file (STELLOPT controlled)
    ! 'chease'
    !geomfile = '.'          ! Name of geometry file
    !x_def = 'arho_t'        ! Definition of x variable
    ! Optional
    !dpdx_pm = 0              ! Amplitude of negative pressure gradient
    dpdx_term = 'full_drift'  ! Treatment of pressure gradient term
    /
    ! The code will read each SPECIES namelist Sequentially based on the value n_spec
    ! name = 'string' ! A arbitrary string which defines a species
    ! passive = .f.   ! Treat as tracer species with no influence on fields
    ! omn = 0.0       ! Normalized density gradient
    ! omt = 4.0       ! Normalized temperature gradient
    ! mass = 1.0      ! Mass normalized to mass scale (mref)
    ! charge = 1      ! Charge number of present species (-1 for electrons)
    ! temp = 1.0      ! Temperature normalized to (Tref)
    ! dens = 1.0      ! Density normalied to electron density
    ! prof_type = 0   ! Profile type (0: use above values, otherwise they're calculated)
    &species
    name =   'ions'
    omn =    0.0
    omt =    4.0
    mass =    1.0
    charge =    1
    temp =    1.0
    dens =    1.0
    /
    &species
    name =   'electrons'
    omn =    2.22
    omt =    6.92
    mass =   0.0025
    charge =   -1
    temp =    1.0
    dens =    1.0
    /
