MODULE AJAX_MOD
!-------------------------------------------------------------------------------
!AJAX-(Take your pick)
!    -Algebraic Jacobians for Advanced eXperiments
!    -The son of Telamon of Salamis and a warrior of great stature and prowess
!     who fought against Troy
!    -The son of Ileus of Locris and a warrior of small stature and arrogant
!     character who fought against Troy
!    -Brand name of a cleaning agent
!
!AJAX_MOD is an F90 module of routines that serve as an interface between a
!  transport code and MHD equilibrium solutions.
!
!References:
!
!  W.A.Houlberg, F90 free format 8/2004
!
!Contains PUBLIC routines:
!
!  Input:
!
!    AJAX_LOAD_RZBDY   -loads approx to 2D MHD equilibria from boundary values
!    AJAX_LOAD_RZLAM   -loads 2D/3D MHD equilibria in inverse coordinate form
!    AJAX_LOAD_MAGFLUX -loads magnetic flux data
!
!  Coordinate conversions:
!
!    AJAX_CAR2CYL      -converts from Cartesian to cylindrical coordinates
!    AJAX_CYL2CAR      -converts from cylindrical to Cartesian coordinates
!    AJAX_CYL2FLX      -converts from cylindrical to flux coordinates
!    AJAX_FLX2CYL      -converts from flux to cylindrical coordinates
!
!  Local magnetics:
!
!    AJAX_B            -gets components of B in various forms
!    AJAX_LAMBDA       -gets stream function and its derivatives
!
!  Flux surface quantities:
!
!    AJAX_FLUXAV_B     -gets flux surface quantities that depend on B
!    AJAX_FLUXAV_G     -gets flux surface quantities that depend on geometry
!    AJAX_I            -gets enclosed toroidal and external poloidal currents
!                       and maximum and minimum values
!    AJAX_MAGFLUX      -gets poloidal and toroidal magnetic fluxes
!    AJAX_SHAPE        -gets shift, elongation, triangularity, Rmax, Rmin, etc
!
!  Miscellaneous:
!
!    AJAX_GLOBALS      -gets global characteristics of the data
!    AJAX_MINMAX_RZ    -gets maximum or minimum of R or Z
!
!Comments:
!
!  This module is designed for use in a transport code that calculates MHD
!    equilibria (flux surface geometry) at infrequent intervals and updates the
!    the magnetic flux (e.g., safety factor or rotational transform) more
!    frequently (e.g., every time step) in a Grad-Hogan approach:
!    1) After a new equilibrium (geometry and magnetic flux) call:
!       AJAX_LOAD_RZLAM or AJAX_LOAD_RZBDY
!       AJAX_LOAD_MAGFLUX
!    2) After each additional timestep (magnetic flux) call:
!       AJAX_LOAD_MAGFLUX
!
!  The flux surface averaging routines are split to provide efficient
!    evaluation:
!    1) After a new equilibrium (geometry dependent) call:
!       AJAX_FLUXAV_G
!       AJAX_SHAPE
!    2) After each timestep (geometry and magnetic flux dependent) call:
!       AJAX_FLUXAV_B
!
!  There is extensive use of optional arguments to:
!    1) Provide flexibility of input variables
!    2) Give the user control over the fineness of the internal grids
!
!  The modernization of the code structure into an F90 module takes advantage of
!    some of the more attractive features of F90:
!    -use of KIND for precision declarations
!    -optional arguments for I/O
!    -generic names for all intrinsic functions
!    -compilation using either free or fixed form
!    -no common blocks or other deprecated Fortran features
!    -dynamic and automatic alocation of variables
!    -array syntax for vector operations
!-------------------------------------------------------------------------------
USE SPEC_KIND_MOD
USE LINEAR1_MOD
USE SPLINE1_MOD
IMPLICIT NONE

!-------------------------------------------------------------------------------
!Private procedures
!-------------------------------------------------------------------------------
PRIVATE :: &
  AJAX_INIT,           & !initializes or resets private data and its allocations
                         !  called from AJAX_LOAD_RZLAM
  AJAX_INIT_FLUXAV_B,  & !initializes magnetic field arrays
                         !  called from AJAX_FLUXAV_B
  AJAX_INIT_FLUXAV_G,  & !initializes geometry arrays
                         !  called from AJAX_FLUXAV_B
                         !  called from AJAX_FLUXAV_G
  AJAX_INIT_LAMBDA,    & !calculates the magnetic stream function from R,Z data
                         !  version is only valid for axisymmetric plasmas
                         !  called from AJAX_LOAD_RZLAM
  AJAX_LOAD_LAMBDA       !loads the magnetic stream function
                         !  called from AJAX_LOAD_RZLAM

!-------------------------------------------------------------------------------
!Private data
!-------------------------------------------------------------------------------
!Logical switches
LOGICAL, PRIVATE, SAVE :: &     
  l_fluxavb_3d,        & !option for whether b-dependent integrals set up [-]
  l_fluxavg_3d,        & !option for whether geom-dependent integrals set up [-]
  l_mfilter_3d           !option for mode filtering [logical]

!Poloidal and toroidal mode expansions
INTEGER, PRIVATE, SAVE :: &     
  krz_3d,              & !no. of modes in RZ expansion [-]
  klam_3d,             & !no. of modes in lambda expansion [-]
  km0n0_3d,            & !index of (m,n)=(0,0) mode [-]
  km1n0_3d,            & !index of (m,n)=(1,0) mode [-]
  nper_3d                !no. of toroidal periods (3-D) [-]

INTEGER, PRIVATE, SAVE, ALLOCATABLE :: &   
  m_3d(:),             & !poloidal modes [-]
  n_3d(:),             & !toroidal modes [-]
  mabs_3d(:)             !|m| in normalization, but restricted by l_mfilter [-]

!Radial, poloidal and toroidal grids
INTEGER, PRIVATE, SAVE :: & 
  nrho_3d,             & !no. of radial nodes [-]
  ntheta_3d,           & !no. of poloidal nodes [-]
  nzeta_3d               !no. of toroidal nodes [-]

REAL(KIND=rspec), PRIVATE, SAVE :: &
  dtheta_3d,           & !poloidal step size for integrals [rad]
  dzeta_3d               !toroidal step size for integrals [rad]

REAL(KIND=rspec), PRIVATE, SAVE, ALLOCATABLE :: &  
  rho_3d(:),           & !radial grid prop to sqrt(tor flux) [-]
  wtheta_3d(:),        & !pol weighting factor for integrals [-]
  wzeta_3d(:)            !tor weighting factor for integrals [-]

!Toroidal flux and magnetic field directions
REAL(KIND=rspec), PRIVATE, SAVE :: &    
  phitot_3d,           & !total toroidal flux [Wb]
  signbp_3d,           & !sign of the pol magnetic field [-]
  signbt_3d              !sign of the toroidal magnetic field [-]

!Splined quantities in radial coordinate:
REAL(KIND=rspec), PRIVATE, SAVE, ALLOCATABLE :: &  
  iotabar_3d(:,:),     & !rotational transform/2pi [-]
  r_3d(:,:,:),         & !expansion coefficients for R [m]
  z_3d(:,:,:),         & !expansion coefficients for Z [m]
  lam_3d(:,:,:)          !expansion coefficients for lambda [-]
  
!Non-stellarator symmetric Splined quantities in radial coordinate:     !SAL 9/9/11
REAL(KIND=rspec), PRIVATE, SAVE, ALLOCATABLE :: &  
  rs_3d(:,:,:),         & !expansion coefficients for R [m] (sin)
  zc_3d(:,:,:),         & !expansion coefficients for Z [m] (cos)
  lamc_3d(:,:,:)          !expansion coefficients for lambda [-] (cos)

!Arrays for flux surface averaging
REAL(KIND=rspec), PRIVATE, SAVE, ALLOCATABLE :: &  
  az_3d(:),            & !temp tor array [arb]
  gt_3d(:),            & !temp pol array [arb]
  vp_3d(:),            & !V'(rho) [m**3/rho]
  rcyl_3d(:,:,:),      & !R coords on internal grid [m]
  zcyl_3d(:,:,:),      & !Z coords on internal grid [m]
  gsqrt_3d(:,:,:),     & !3D Jacobian on internal grid [m**3/rho]
  eltheta_3d(:,:,:),   & !d(lambda)/d(theta) on internal grid [/rad]
  elzeta_3d(:,:,:),    & !d(lambda)/d(zeta) on internal grid [/rad]
  b_3d(:,:,:),         & !total B on internal grid [T]
  btheta_3d(:,:,:),    & !contravariant poloidal B on internal grid [T/m]
  gcyl_3d(:,:,:,:)       !R,Z derivatives
                         !=(R_rho,R_theta,R_zeta,Z_rho,Z_theta,Z_zeta)
                         ! [m/rho,m,      m,     m/rho,m,      m     ]

!Major and minor radius quantities
REAL(KIND=rspec), PRIVATE, SAVE :: &    
  r000_3d,             & !coefficient of (m,n)=(0,0) mode for R at axis [m]
  rhomin_3d,           & !inner radial boundary of R,Z domain [rho]
  rhomax_3d,           & !outer radial boundary of R,Z domain [rho]
  rhores_3d              !resolution factor for rho [rho]

!Physical and conversion constants
REAL(KIND=rspec), PRIVATE, PARAMETER :: &   
  z_mu0=1.2566e-06,    & !permeability of free space [H/m]
  z_pi=3.141592654       !pi [-]

REAL(KIND=rspec), PRIVATE, SAVE :: &    
  z_large,             & !largest real number [-]
  z_precision,         & !machine precision [-]
  z_small                !smallest  real number [-]

!-------------------------------------------------------------------------------
! Public procedures
!-------------------------------------------------------------------------------
CONTAINS

SUBROUTINE AJAX_LOAD_RZBDY(r0,a0,s0,e0,e1,d1, &
                           iflag,message, & 
                           NRHO_AJAX,NTHETA_AJAX,NKLAM_AJAX,RHOMAX_AJAX)
!-------------------------------------------------------------------------------
!AJAX_LOAD_RZBDY converts geometric boundary and axial constraints to values of
!  the R,Z coordinates then calls AJAX_LOAD_RZLAM to fill in the profiles
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &    
  r0,                  & !major radius of geometric center [m]
  a0,                  & !minor radius in midplane [m]
  s0,                  & !axis shift normalized to a0 [-]
  e0,                  & !axis elongation normalized to a0 [-]
  e1,                  & !edge elongation normalized to a0 [-]
  d1                     !edge triangularity normalized to a0 [-]

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error

!Declaration of optional input variables
INTEGER, INTENT(IN), OPTIONAL :: &    
  NRHO_AJAX,           & !no. of radial nodes in internal data [-]
                         !=21 default
  NTHETA_AJAX,         & !no. of poloidal nodes in internal data [odd]
                         !=21 default
  NKLAM_AJAX             !no. of lambda modes in nternal data [-]
                         !=5 default

REAL(KIND=rspec), INTENT(IN), OPTIONAL :: &   
  RHOMAX_AJAX            !value of internal radial grid at R,Z boundary [rho]
                         !=1.0 default

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &      
  i,k,nk_rz,nk_lam,nr_rz,nthetal

INTEGER, ALLOCATABLE :: &     
  m(:),n(:)

REAL(KIND=rspec) :: &     
  c,dm1,em1,r2,rhomaxl

REAL(KIND=rspec), ALLOCATABLE :: &    
  rho_rz(:),r(:,:),z(:,:)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Null output
iflag=0
message=''

!Check optional input
IF(PRESENT(RHOMAX_AJAX)) THEN

  rhomaxl=RHOMAX_AJAX

ELSE

  rhomaxl=1

ENDIF

IF(PRESENT(NRHO_AJAX)) THEN

  nr_rz=NRHO_AJAX

ELSE

  nr_rz=31

ENDIF

IF(PRESENT(NKLAM_AJAX)) THEN

  nk_lam=NKLAM_AJAX

ELSE

  nk_lam=8

ENDIF

IF(PRESENT(NTHETA_AJAX)) THEN

  nthetal=NTHETA_AJAX

ELSE

  nthetal=4*nk_lam+1

ENDIF

!Allocate and set poloidal mode information
ALLOCATE(m(nk_lam), &     
         n(nk_lam))

  nk_rz=3
  m(:)=(/ (k-1,k=1,nk_lam) /)
  n(:)=0

!-------------------------------------------------------------------------------
!Convert from geometric to moments quantities     
!-------------------------------------------------------------------------------
!Moments triangularity ~geometric triangularity/4
dm1=d1/4
c=0

!Iterate for triangularity value, convergence is fast and robust
DO i=1,10 !Over iteration

  c=4*dm1/(SQRT(1.0+32*dm1**2)+1.0)
  dm1=d1/(4.0-6*c**2)

ENDDO !Over iteration

!Elongation
em1=e1/(SQRT(1.0-c**2)*(1.0+2*dm1*c))

!-------------------------------------------------------------------------------
!Fill in radial values of R,Z expansion coefficients
!-------------------------------------------------------------------------------
!Allocate radial grid and R,Z arrays
ALLOCATE(rho_rz(nr_rz), &     
         r(nr_rz,nk_rz), &    
         z(nr_rz,nk_rz))

  rho_rz(:)=0
  r(:,:)=0
  z(:,:)=0

!Set radial grid
rho_rz(:)=(/ (REAL(i-1,rspec)/REAL(nr_rz-1,rspec),i=1,nr_rz) /)

!Radial variation
DO i=1,nr_rz !Over radial nodes

  r2=rho_rz(i)**2
  r(i,1)=r0+a0*(s0*(1.0-r2)-dm1*r2)
  r(i,2)=a0*rho_rz(i)
  z(i,2)=-r(i,2)*(e0*(1.0-r2)+em1*r2)
  r(i,3)=a0*r2*dm1
  z(i,3)=r(i,3)*em1
 
ENDDO !Over radial nodes

!Load private data
CALL AJAX_LOAD_RZLAM(nr_rz,nk_rz,rho_rz,m,n,r,z, &
                     iflag,message, &    
                     NRHO_AJAX=nr_rz, &   
                     NTHETA_AJAX=nthetal, &   
                     RHOMAX_AJAX=rhomaxl, &   
                     NK_LAM=nk_lam)

!Check messages
IF(iflag /= 0) message='AJAX_LOAD_RZBDY/'//message

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

IF(ALLOCATED(m)) THEN

  !Deallocate mode information
  DEALLOCATE(m, &      
             n)

ENDIF

IF(ALLOCATED(rho_rz)) THEN

  !Deallocate grid and R,Z arrays
  DEALLOCATE(rho_rz, &     
             r, &      
             z)

ENDIF

END SUBROUTINE AJAX_LOAD_RZBDY

SUBROUTINE AJAX_LOAD_RZLAM(nr_rz,nk_rz,rho_rz,m,n,r,z, &
                           iflag,message, &    
                           K_GRID,L_MFILTER_AJAX,NRHO_AJAX,NTHETA_AJAX, & 
                           NZETA_AJAX,RHOMAX_AJAX,NR_LAM,NK_LAM,RHO_LAM,LAM, &
                           R_S,Z_C,LAM_C)
!-------------------------------------------------------------------------------
!AJAX_LOAD_RZLAM loads 2D/3D MHD equilibria in inverse coordinate form
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!
!Comments:
!  Knowledge of the form of the input radial grid (see K_GRID) is necessary for
!    accurate mapping to the internal grid
!  The internal radial grid, rho_3d, is proportional to the square root of the
!    toroidal flux
!  The maximum value of the R,Z radial grid defines rho/rhomax=1
!  The internal radial grid will be scaled to rhomax if specified, otherwise it
!    will default to the domain [0,1]
!  The input poloidal and toroidal mode numbers are the same for the R,Z and
!    lambda expansions and dimensioned to the greater of the two; e.g., if there
!    are more poloidal and toroidal modes for the input lambda they should be
!    appended to the end of the sequence of values for the R,Z expansions
!  Mode filtering may be tried to improve calculations near the axis and the
!    outer boundary
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &     
  nr_rz,               & !no. of radial nodes in the input R,Z [-]
  nk_rz,               & !no. of poloidal & toroidal modes in the input R,Z [-]
  m(:),                & !poloidal mode numbers [-]
  n(:)                   !toroidal mode numbers [-]
     
REAL(KIND=rspec), INTENT(IN) :: &                                      
  rho_rz(:),           & !radial nodes in the input R,Z [arb]
  r(:,:),              & !expansion coeffs for R [m] (cos)
  z(:,:)                 !expansion coeffs for Z [m] (sin)

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error
  
!Declaration of optional input variables
LOGICAL, INTENT(IN), OPTIONAL :: &    
  L_MFILTER_AJAX         !option for mode filtering [logical]

INTEGER, INTENT(IN), OPTIONAL :: &    
  K_GRID,              & !option designating type of input radial grid [-]
                         !=0 proportional to sqrt(toroidal flux)
                         !=1 proportional to toroidal flux
                         !=else not allowed
  NRHO_AJAX,           & !no. of radial nodes in internal data [-]
                         !=21 default
  NTHETA_AJAX,         & !no. of poloidal nodes in internal data [odd]
                         !=21 default
  NZETA_AJAX,          & !no. of toroidal nodes per period in internal data [odd]
                         !=11 default for non-axisymmetric plasma
                         !=0 default for axisymmetric plasma
  NR_LAM,              & !no. of radial nodes in the input lambda [-]
  NK_LAM                 !no. of poloidal & toroidal modes in the input lambda [-]

REAL(KIND=rspec), INTENT(IN), OPTIONAL :: &   
  RHOMAX_AJAX,         & !value of internal radial grid at R,Z boundary [rho]
                         !=1.0 default
  RHO_LAM(:),          & !radial nodes in the input lambda [arb]
  LAM(:,:)               !expansion coeffs for lambda [-] (sin)

REAL(KIND=rspec), INTENT(IN), OPTIONAL :: & !SAL 9/9/11
  R_S(:,:),            & !expansion coeffs for R [m] (sin)
  Z_C(:,:),            & !expansion coeffs for Z [m] (cos)
  LAM_C(:,:)             !expansion coeffs for lambda [-] (cos)

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &      
  i,j,k,kmax,nset1,k_grid_l,k_vopt(1:3)=(/1,0,0/),k_bc1=3,k_bcn=0
  
REAL(KIND=rspec), ALLOCATABLE :: &    
  rho(:),rmn(:,:),zmn(:,:),fspl(:,:),values(:,:)
  
REAL(KIND=rspec), ALLOCATABLE :: &               !SAL 9/9/11
  rmns(:,:),zmnc(:,:)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Null output
iflag=0
message=''

!Set number of angular modes
krz_3d=nk_rz

IF(PRESENT(NK_LAM)) THEN

  klam_3d=NK_LAM

ELSE

  klam_3d=MAX(krz_3d,8)

ENDIF

kmax=MAX(krz_3d,klam_3d)

!Set the number of toroidal field periods, < 100 periods
i=100

DO k=1,kmax !Over modes

  !Find lowest non-zero toroidal mode number
  IF(n(k) /= 0) i=MIN(i,ABS(n(k)))

ENDDO !Over modes

IF(i == 100) THEN

  !Found no peiodicity in data, set to axisymmetric
   nper_3d=0

ELSE

  !Non-axisymmetric 
  nper_3d=i

  !Check periodicity
  DO k=1,kmax !Over modes

    IF(n(k) /= 0) THEN

      j=MOD(n(k),nper_3d)

      IF(j /= 0) THEN

        !Some mode does not fit the periodicity
        iflag=1
        message='AJAX_LOAD_RZLAM/ERROR(1):periodicity not found'
        GOTO 9999

      ENDIF

    ENDIF

  ENDDO !Over modes

ENDIF

!-------------------------------------------------------------------------------
!Check optional input for internal grid initialization
!-------------------------------------------------------------------------------
!Set the number of internal radial nodes
IF(PRESENT(NRHO_AJAX)) THEN

  !Use input value
  nrho_3d=NRHO_AJAX

ELSE

  !Use default value
  nrho_3d=31

ENDIF

!Set the number of internal poloidal nodes
IF(PRESENT(NTHETA_AJAX)) THEN

  !Use input value
  ntheta_3d=NTHETA_AJAX

  !Make sure the number of modes is odd for 4th order Simpson integration
  IF(MOD(ntheta_3d,2) == 0) THEN

    ntheta_3d=ntheta_3d+1
    iflag=-1
    message='AJAX_LOAD_RZLAM/WARNING(2):poloidal nodes reset'

  ENDIF

ELSE

  IF(nper_3d == 0) THEN
  
    !Axisymmetric
    !4k+1 typically yields <0.1% error in B_zeta = RB_t ~ const
    ntheta_3d=4*klam_3d+1

  ELSE

    !Non-axisymmetric
    !Error estimate yet to be quantified, typically |m| <= 8
    ntheta_3d=33

  ENDIF

ENDIF

!Set the number of internal toroidal nodes
IF(PRESENT(NZETA_AJAX)) THEN

  !Use input value
  nzeta_3d=NZETA_AJAX

  IF(nper_3d == 0 .AND. nzeta_3d /= 1) THEN

    nzeta_3d=1
    iflag=-1
    message='AJAX_LOAD_RZLAM/WARNING(3):toroidal nodes reset to 1'

  ENDIF

  !Make sure the number of modes is odd for 4th order Simpson integration
  IF(MOD(nzeta_3d,2) == 0) THEN

    nzeta_3d=nzeta_3d+1
    iflag=-1
    message='AJAX_LOAD_RZLAM/WARNING(4):toroidal nodes reset'

  ENDIF

ELSEIF(nper_3d == 0) THEN

  !Use axisymmetric value
  nzeta_3d=1

ELSE

  !Use non-axisymmetric default value
  !Error estimate yet to be quantified, typically |n/nper| <= 5
  nzeta_3d=21

ENDIF

!Set the option for mode filtering
IF(PRESENT(L_MFILTER_AJAX)) THEN

  !Use input value
  l_mfilter_3d=L_MFILTER_AJAX

ELSE

  !Use default value
  l_mfilter_3d=.TRUE.

ENDIF

!-------------------------------------------------------------------------------
!Set general information and initialize arrays
!-------------------------------------------------------------------------------
!Scale the outer radial boundary of R,Z to rhomax
IF(PRESENT(RHOMAX_AJAX)) THEN

  !Use input value
  rhomax_3d=RHOMAX_AJAX

ELSE

  !Use default value
  rhomax_3d=1

ENDIF

!Set the mode numbers of the (m,n)=(0,0) and (m,n)=(1,0) modes
km0n0_3d=0
km1n0_3d=0

LOOP_K: DO k=1,kmax

  IF(ABS(m(k)) == 0 .AND. n(k) == 0) km0n0_3d=k
  IF(ABS(m(k)) == 1 .AND. n(k) == 0) km1n0_3d=k
  IF(km0n0_3d /= 0 .AND. km1n0_3d /= 0) EXIT LOOP_K

ENDDO LOOP_K !Over modes

!Call initialization routine
CALL AJAX_INIT

!-------------------------------------------------------------------------------
!Poloidal and toroidal mode data
!-------------------------------------------------------------------------------
!Initialization
m_3d(:)=0
n_3d(:)=0

!Copy data
m_3d(1:kmax)=m(1:kmax)
n_3d(1:kmax)=n(1:kmax)

!Set |m| with mode filtering for rho**|m| normalization
DO k=1,kmax !Over modes

  mabs_3d(k)=ABS(m_3d(k))

  !Check for filtering
  IF(l_mfilter_3d) THEN

    IF(mabs_3d(k) > 3) THEN

      !Filtering for |m| > 3
      IF(MOD(mabs_3d(k),2) == 0) THEN

        !Even mode, remove rho**2
        mabs_3d(k)=2

      ELSE

        !Odd mode, remove rho**3
        mabs_3d(k)=3

      ENDIF

    ENDIF

  ENDIF

ENDDO !Over modes

!-------------------------------------------------------------------------------
!Load R and Z expansion data
!-------------------------------------------------------------------------------
!Initialization
r_3d(:,:,:)=0
z_3d(:,:,:)=0

!Make copy of radial grid and expansion coefficients
!If the user does not provide an axial value for R,Z data, need to fill in
IF(rho_rz(1) > rhores_3d) THEN

  !Allocate radial grid and R,Z arrays, add radial node at axis
  nset1=nr_rz+1
  ALLOCATE(rho(nset1), &     
           rmn(nset1,nk_rz), &    
           zmn(nset1,nk_rz))

    rho(:)=0
    rmn(:,:)=0
    zmn(:,:)=0
    rho(2:nset1)=rho_rz(1:nr_rz)/rho_rz(nr_rz)
    rmn(2:nset1,1:nk_rz)=r(1:nr_rz,1:nk_rz)
    zmn(2:nset1,1:nk_rz)=z(1:nr_rz,1:nk_rz)
           
  ! SAL 9/9/11
  IF (PRESENT(R_S)) THEN
     ALLOCATE(rmns(nset1,nk_rz), &
              zmnc(nset1,nk_rz))
    rmns(2:nset1,1:nk_rz)=R_S(1:nr_rz,1:nk_rz)
    zmnc(2:nset1,1:nk_rz)=Z_C(1:nr_rz,1:nk_rz)
  ENDIF

ELSE

  !Allocate radial grid and R,Z arrays, use input grid
  nset1=nr_rz
  ALLOCATE(rho(nset1), &     
           rmn(nset1,nk_rz), &    
           zmn(nset1,nk_rz))

    rho(:)=0
    rmn(:,:)=0
    zmn(:,:)=0
    rho(1:nset1)=rho_rz(1:nset1)/rho_rz(nset1)
    rmn(1:nset1,1:nk_rz)=r(1:nset1,1:nk_rz)
    zmn(1:nset1,1:nk_rz)=z(1:nset1,1:nk_rz)
           
  ! SAL 9/9/11
  IF (PRESENT(R_S)) THEN
     ALLOCATE(rmns(nset1,nk_rz), &
              zmnc(nset1,nk_rz))
    rmns(1:nset1,1:nk_rz)=R_S(1:nset1,1:nk_rz)
    zmnc(1:nset1,1:nk_rz)=Z_C(1:nset1,1:nk_rz)
  ENDIF

ENDIF

!Convert rho to sqrt(toroidal flux) if necesssary and scale
k_grid_l=0
IF(PRESENT(K_GRID)) k_grid_l=K_GRID

IF(k_grid_l == 1) THEN

  !~toroidal flux
  rho(:)=SQRT(rho(:))*rhomax_3d

ELSEIF(k_grid_l == 0) THEN

  !~sqrt(toroidal flux)
  rho(:)=rho(:)*rhomax_3d

ELSE

  !Unallowed choice of input grid
  iflag=1
  message='AJAX_LOAD_RZLAM/ERROR(5):unallowed choice of K_GRID'
  GOTO 9999

ENDIF

!Set the inner radial boundary of R,Z for checking axial extrapolation
rhomin_3d=rho(2)

!Allocate spline arrays
ALLOCATE(fspl(4,nset1), &     
         values(3,nrho_3d))

  fspl(:,:)=0
  values(:,:)=0

!Normalize the expansion coefficients to rho**m
DO k=1,krz_3d !Over modes

  DO i=2,nset1 !Over radial nodes

    rmn(i,k)=rmn(i,k)/rho(i)**mabs_3d(k)
    zmn(i,k)=zmn(i,k)/rho(i)**mabs_3d(k)
    
    IF (PRESENT(R_S)) THEN  ! 9/9/11 SAL
       rmns(i,k)=rmns(i,k)/rho(i)**mabs_3d(k)
       zmnc(i,k)=zmnc(i,k)/rho(i)**mabs_3d(k)
    END IF

  ENDDO !Over radial nodes

!Parabolic extrapolation to axis (user values ignored)
  rmn(1,k)=(rmn(2,k)*rho(3)**2-rmn(3,k)*rho(2)**2)/(rho(3)**2-rho(2)**2)
  zmn(1,k)=(zmn(2,k)*rho(3)**2-zmn(3,k)*rho(2)**2)/(rho(3)**2-rho(2)**2)
  
  IF (PRESENT(R_S)) THEN  ! 9/9/11 SAL
    rmns(1,k)=(rmns(2,k)*rho(3)**2-rmns(3,k)*rho(2)**2)/(rho(3)**2-rho(2)**2)
    zmnc(1,k)=(zmnc(2,k)*rho(3)**2-zmnc(3,k)*rho(2)**2)/(rho(3)**2-rho(2)**2)
  END IF

!Map R_mn to internal radial grid
  fspl(1,1:nset1)=rmn(1:nset1,k)
  iflag=0
  message=''
  CALL SPLINE1_INTERP(k_vopt,nset1,rho,fspl,nrho_3d,rho_3d,values, &
                      iflag,message, &    
                      K_BC1=k_bc1, &    
                      K_BCN=k_bcn)

  !Check messages
  IF(iflag > 0) THEN

    message='AJAX_LOAD_RZLAM(6)/'//message
    GOTO 9999

  ENDIF

  !Respline the R coeffs for internal storage
  r_3d(1,1:nrho_3d,k)=values(1,1:nrho_3d)
  CALL SPLINE1_FIT(nrho_3d,rho_3d,r_3d(:,:,k), &  
                   K_BC1=k_bc1, &    
                   K_BCN=k_bcn)

  IF (PRESENT(R_S)) THEN ! 9/9/11 SAL
  !Map Rs_mn to internal radial grid
    fspl(1,1:nset1)=rmns(1:nset1,k)
    iflag=0
    message=''
    CALL SPLINE1_INTERP(k_vopt,nset1,rho,fspl,nrho_3d,rho_3d,values, &
                      iflag,message, &    
                      K_BC1=k_bc1, &    
                      K_BCN=k_bcn)

    !Check messages
    IF(iflag > 0) THEN

      message='AJAX_LOAD_RZLAM(6)/'//message
      GOTO 9999

    ENDIF

    !Respline the R coeffs for internal storage
    rs_3d(1,1:nrho_3d,k)=values(1,1:nrho_3d)
    CALL SPLINE1_FIT(nrho_3d,rho_3d,rs_3d(:,:,k), &  
                   K_BC1=k_bc1, &    
                   K_BCN=k_bcn)
  END IF
   
!Map Z_mn to internal radial grid
  fspl(1,1:nset1)=zmn(1:nset1,k)
  iflag=0
  message=''
  CALL SPLINE1_INTERP(k_vopt,nset1,rho,fspl,nrho_3d,rho_3d, & 
                      values,iflag,message, &   
                      K_BC1=k_bc1, &    
                      K_BCN=k_bcn)

  !Check messages
  IF(iflag > 0) THEN

    message='AJAX_LOAD_RZLAM(7)/'//message
    GOTO 9999

  ENDIF

  !Respline the Z coeffs for internal storage
  z_3d(1,1:nrho_3d,k)=values(1,1:nrho_3d)        
  CALL SPLINE1_FIT(nrho_3d,rho_3d,z_3d(:,:,k), &  
                   K_BC1=k_bc1, &    
                   K_BCN=k_bcn)
                   
  IF (PRESENT(R_S)) THEN ! 9/9/11 SAL
    !Map Zc_mn to internal radial grid
    fspl(1,1:nset1)=zmnc(1:nset1,k)
    iflag=0
    message=''
    CALL SPLINE1_INTERP(k_vopt,nset1,rho,fspl,nrho_3d,rho_3d,values, &
                      iflag,message, &    
                      K_BC1=k_bc1, &    
                      K_BCN=k_bcn)

    !Check messages
    IF(iflag > 0) THEN

      message='AJAX_LOAD_RZLAM(6)/'//message
      GOTO 9999

    ENDIF

    !Respline the R coeffs for internal storage
    zc_3d(1,1:nrho_3d,k)=values(1,1:nrho_3d)
    CALL SPLINE1_FIT(nrho_3d,rho_3d,zc_3d(:,:,k), &  
                   K_BC1=k_bc1, &    
                   K_BCN=k_bcn)
  END IF

ENDDO !Over modes

!Set the length scale factors
r000_3d=r_3d(1,1,km0n0_3d)

!-------------------------------------------------------------------------------
!Magnetic stream function
!-------------------------------------------------------------------------------
IF(PRESENT(NR_LAM) .AND. &     
   PRESENT(NK_LAM) .AND. &     
   PRESENT(RHO_LAM) .AND. &     
   PRESENT(LAM)) THEN

  !Load lambda from input
  iflag=0
  message=''
  
!  SAL 09/09/11
!  CALL AJAX_LOAD_LAMBDA(k_grid_l,NR_LAM,NK_LAM,rho_rz(nr_rz), & 
!                        RHO_LAM,LAM,iflag,message)
  IF(PRESENT(LAM_C)) THEN
    CALL AJAX_LOAD_LAMBDA(k_grid_l,NR_LAM,NK_LAM,rho_rz(nr_rz), & 
                        RHO_LAM,LAM,iflag,message,LAM_C2=LAM_C)
  ELSE
    CALL AJAX_LOAD_LAMBDA(k_grid_l,NR_LAM,NK_LAM,rho_rz(nr_rz), & 
                        RHO_LAM,LAM,iflag,message)
  ENDIF
      

  !Check messages
  IF(iflag /= 0) message='AJAX_LOAD_RZLAM(8)/'//message

ELSEIF(nper_3d == 0) THEN

  !Calculate lambdas for axisymmetric plasma
  !Set up 3d grid for flux surface integrals
  iflag=0
  message=''
  CALL AJAX_INIT_FLUXAV_G(iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_LOAD_RZLAM(9)/'//message
    IF(iflag > 0) GOTO 9999

  ENDIF


  !Calculate lambda
  CALL AJAX_INIT_LAMBDA

ELSE

  !Need to specify stream function for non-axisymmetric plasmas
  iflag=1
  message='AJAX_LOAD_RZLAM/ERROR(10):need lambdas'

ENDIF

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

IF(ALLOCATED(rho)) THEN

  !Deallocate radial grid and R,Z arrays
  DEALLOCATE(rho, &      
             rmn, &      
             zmn)

ENDIF

IF (ALLOCATED(rmns)) THEN    ! 9/9/11 SAL
   DEALLOCATE(rmns, &
              zmnc)
ENDIF

IF(ALLOCATED(fspl)) THEN

  !Deallocate spline arrays
  DEALLOCATE(fspl, &      
             values)

ENDIF

END SUBROUTINE AJAX_LOAD_RZLAM

SUBROUTINE AJAX_LOAD_MAGFLUX(phitot,k_pflx,nr_pflx,rho_pflx,pflx, &
                             iflag,message)
!-------------------------------------------------------------------------------
!AJAX_LOAD_MAGFLUX loads magnetic flux data
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &     
  k_pflx,              & !option for represention of poloidal flux [-]
                         !=1 iotabar (rotational transform)
                         !=2 d(Psi)/d(rho)
                         !=else q (safety factor)
  nr_pflx                !no. of radial nodes in input flux [-]

REAL(KIND=rspec), INTENT(IN) :: &    
  rho_pflx(:),         & !radial nodes in input flux [rho]
  pflx(:),             & !poloidal magnetic flux function
                         !=q [-]
                         !=iotabar [-]
                         !=d(Psi)/d(rho) [Wb/rho]
  phitot                 !total toroidal flux [Wb]

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER, PARAMETER :: &     
  k_vopt(1:3)=(/1,0,0/),k_bc1=3,k_bcn=0 

INTEGER :: &      
  nset1
 
REAL(KIND=rspec) :: &     
  fspl(1:4,1:nr_pflx),values(1:3,1:nrho_3d)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Null output
iflag=0
message=''

!-------------------------------------------------------------------------------
!Iotabar allocation
!-------------------------------------------------------------------------------
IF(ALLOCATED(iotabar_3d)) THEN

  !Radial dimension
  nset1=SIZE(iotabar_3d,2)

  !If storage requirements have changed, reallocate
  IF(nset1 /= nrho_3d) THEN

    !Reallocate iotabar
    DEALLOCATE(iotabar_3d)
    ALLOCATE(iotabar_3d(4,nrho_3d))

      iotabar_3d(:,:)=0

  ENDIF

ELSE

  !Allocate iotabar
  ALLOCATE(iotabar_3d(4,nrho_3d))

    iotabar_3d(:,:)=0

ENDIF

!-------------------------------------------------------------------------------
!Total toroidal flux in private data
!-------------------------------------------------------------------------------
phitot_3d=phitot

!-------------------------------------------------------------------------------
!Set iotabar in private data
!-------------------------------------------------------------------------------
!Initialize
fspl(:,:)=0
values(:,:)=0

IF(k_pflx == 1) THEN

  !iotabar input
  fspl(1,1:nr_pflx)=pflx(1:nr_pflx)

ELSEIF(k_pflx == 2) THEN

  !d(Psi)/d(rho) input
  fspl(1,1:nr_pflx)=pflx(1:nr_pflx)/rho_pflx(1:nr_pflx)
  fspl(1,1:nr_pflx)=rhomax_3d**2/phitot_3d/2*fspl(1,1:nr_pflx)

ELSE

  !q input
  fspl(1,1:nr_pflx)=1/pflx(1:nr_pflx)

ENDIF

!Interpolate iotabar onto internal grid
iflag=0
message=''
CALL SPLINE1_INTERP(k_vopt,nr_pflx,rho_pflx,fspl,nrho_3d,rho_3d, &
                    values,iflag,message, &   
                    K_BC1=k_bc1, &    
                    K_BCN=k_bcn)

!Check messages
IF(iflag > 0) THEN

  message='AJAX_LOAD_MAGFLUX(2)/'//message
  GOTO 9999

ENDIF

!Get spline coefficients on internal grid
iotabar_3d(:,:)=0
iotabar_3d(1,1:nrho_3d)=values(1,1:nrho_3d)

!If value at origin not supplied, overwrite with approximation
IF(rho_pflx(1) > rhores_3d) THEN

  iotabar_3d(1,1)=iotabar_3d(2,1)-rho_3d(2)*(iotabar_3d(3,1)-iotabar_3d(2,1)) &  
                                  /(rho_3d(3)-rho_3d(2))

ENDIF

CALL SPLINE1_FIT(nrho_3d,rho_3d,iotabar_3d, &   
                 K_BC1=k_bc1, &    
                 K_BCN=k_bcn)

!-------------------------------------------------------------------------------
!Set signs of poloidal and toroidal magnetic fields
!-------------------------------------------------------------------------------
signbt_3d=SIGN(1.0_rspec,phitot_3d)
signbp_3d=SIGN(1.0_rspec,iotabar_3d(1,nrho_3d))*signbt_3d

!-------------------------------------------------------------------------------
!Set logical switches for flux surface averaging
!-------------------------------------------------------------------------------
l_fluxavb_3d=.FALSE.

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE AJAX_LOAD_MAGFLUX

SUBROUTINE AJAX_CAR2CYL(r_car, &
                        r_cyl)
!-------------------------------------------------------------------------------
!AJAX_CAR2CYL converts from Cartesian to cylindrical coordinates
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &    
  r_car(:)               !Cartesian coordinates (x,y,z) [m,m,m]

!Declaration of output variables
REAL(KIND=rspec), INTENT(OUT) :: &    
  r_cyl(:)               !cylindrical coordinates (R,phi,Z) [m,rad,m]

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Null output
r_cyl(:)=0

!-------------------------------------------------------------------------------
!Cartesian to cylindrical conversion
!-------------------------------------------------------------------------------
r_cyl(1)=SQRT(r_car(1)**2+r_car(2)**2) !R
r_cyl(2)=ATAN2(r_car(2),r_car(1))      !phi
r_cyl(3)=r_car(3)                      !Z

!Ensure 0 <= phi <= 2*pi
IF(r_cyl(2) < 0.0) r_cyl(2)=r_cyl(2)+2*z_pi

END SUBROUTINE AJAX_CAR2CYL

SUBROUTINE AJAX_CYL2CAR(r_cyl, &
                        r_car)
!-------------------------------------------------------------------------------
!AJAX_CYL2CAR converts from cylindrical to Cartesian coordinates
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &    
  r_cyl(:)               !cylindrical coordinates (R,phi,Z) [m,rad,m]

!Declaration of output variables
REAL(KIND=rspec), INTENT(OUT) :: &    
  r_car(:)               !Cartesian coordinates (x,y,z) [m,m,m]

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Null output
r_car(:)=0

!-------------------------------------------------------------------------------
!Cylindrical to Cartesian conversion
!-------------------------------------------------------------------------------
r_car(1)=r_cyl(1)*COS(r_cyl(2))   !x
r_car(2)=r_cyl(1)*SIN(r_cyl(2))   !y
r_car(3)=r_cyl(3)                 !z

END SUBROUTINE AJAX_CYL2CAR

SUBROUTINE AJAX_CYL2FLX(r_cyl, &
                        r_flx, &
                        iflag,message, &  
                        G_CYL,GSQRT,TAU)
!-------------------------------------------------------------------------------
!AJAX_CYL2FLX converts from cylindrical to flux coordinates
!
!References:
!  S.E.Attenberger, W.A.Houlberg, S.P.Hirshman J Comp Phys 72 (1987) 435
!  W.A.Houlberg, F90 free format 8/2004
!
!Comments:
!  The basic method is a Newton iteration in 2 dimensions
!  Input values of r_flx are used as an initial guess
!  Convergence is typically one order of magnitude reduction in the
!    error in R and Z per iteration - should rarely exceed 5 iterations
!  The toroidal angle in flux coordinates is the same as in cylindrical
!    coordinates giving a right-handed system with positive Jacobian
!  Point 0 is the previous best guess and point 1 is a trial point
!  If point 1 is further from (R,Z) than point 0, a new trial point is
!    generated by halving the step and using interpolated derivatives
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &    
  r_cyl(:)               !cylindrical coordinates (R,phi,Z) [m,rad,m]

!Declaration of input/output variables
REAL(KIND=rspec), INTENT(INOUT) :: &    
  r_flx(:)               !flux coordinates (rho,theta,zeta) [rho,rad,rad]
      
!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error

!Declaration of optional output variables
REAL(KIND=rspec), INTENT(OUT), OPTIONAL :: &                           
  G_CYL(:),            & !R,Z derivatives
                         !=(R_rho,R_theta,R_zeta,Z_rho,Z_theta,Z_zeta)
                         ! [m/rho,m,      m,     m/rho,m,      m     ]
  GSQRT,               & !3D Jacobian [m**3/rho]
  TAU                    !2D Jacobian in phi=zeta=constant plane [m**2/rho]
     
!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &      
  it,itmax,nh
 
REAL(KIND=rspec) :: &     
  dr,dr0,dr1,dz,dz0,dz1,drho,dtheta,err0,err1,gsqrt1,tau0,tau1, &                                 
  taut,tol

REAL(KIND=rspec) :: &     
  r_flx0(1:3),g_cyl0(1:6),r_flx1(1:3),r_cyl1(1:3),g_cyl1(1:6), & 
  g_cylt(1:6)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Null output
iflag=0
message=''

!Coordinates and metrics for iteration
r_flx0(:)=0
g_cyl0(:)=0
r_flx1(:)=0
r_cyl1(:)=0
g_cyl1(:)=0
g_cylt(:)=0

!Tolerance for convergence
tol=0.1*rhores_3d/rhomax_3d

!Maximum iterations
itmax=20
   
!Grid halving index
nh=1

!Set flux coordinate toroidal angle equal to cylindrical angle
r_flx(3)=r_cyl(2)

!Set starting points and parameters
tau0=0
err0=0.1*z_large
r_flx0(1:3)=r_flx(1:3)
r_flx1(1:3)=r_flx(1:3)

!Move rho away from axis if necessary
IF(r_flx1(1) < rhores_3d) r_flx1(1)=rhores_3d

!-------------------------------------------------------------------------------
!Iterate to find flux coordinates
!-------------------------------------------------------------------------------
DO it=1,itmax !Over iteration

  !Get cylindrical coordinates at point 1
  iflag=0
  message=''
  CALL AJAX_FLX2CYL(r_flx1,r_cyl1,iflag,message, &  
                    G_CYL=g_cyl1, &    
                    GSQRT=gsqrt1, &    
                    TAU=tau1)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_CYL2FLX(1)/'//message
    IF(iflag > 0) GOTO 9999

  ENDIF

  dr1=(r_cyl(1)-r_cyl1(1))
  dz1=(r_cyl(3)-r_cyl1(3))
  err1=(dr1**2+dz1**2)/r000_3d**2
  !Check convergence
  IF((ABS(dr1) <= tol*r000_3d .AND. ABS(dz1) <= tol*r000_3d) & 
    .OR. (ABS(err1-err0) < 5.0e-6)) THEN
!  IF((ABS(dr1) <= tol*r000_3d .AND. ABS(dz1) <= tol*r000_3d) & 
!    .OR. (ABS(err1-err0) < 1.0d-60)) THEN

    !Converged, but first check if solution is in the R,Z domain
    IF(r_flx1(1) > rhomax_3d) THEN

      iflag=-1
      message='AJAX_CYL2FLX(2)/WARNING:outside R,Z domain'

    ENDIF

    !Solution is at point 1, copy to output arrays and exit
    r_flx(1:3)=r_flx1(1:3)
    IF(r_flx(2) < 0.0) r_flx(2)=r_flx(2)+2*z_pi
    r_flx(2)=MOD(r_flx(2),2*z_pi)
    IF(PRESENT(GSQRT)) GSQRT=gsqrt1
    IF(PRESENT(G_CYL)) G_CYL(:)=g_cyl1(:)
    IF(PRESENT(TAU)) TAU=tau1
    GOTO 9999

  ENDIF

  !Need improved estimate for rho and theta and get new Jacobian
  IF(r_flx1(1) < rhores_3d) r_flx1(1)=rhores_3d

  IF(ABS(tau1) < 10*z_small) THEN

    !No solution exists for the Newton method
    !Assume the solution is at the origin where tau=0
    iflag=-1
    message='AJAX_CYL2FLX(3)/WARNING:solution at origin'
    r_flx(1)=0
    r_flx(2)=0

  ENDIF

  !Check consistency of sign of Jacobian
  IF(SIGN(1.0_rspec,tau1) /= SIGN(1.0_rspec,tau0) .AND. tau0 /= 0.0) THEN

    !Bad data
    iflag=1
    message='AJAX_CYL2FLX(4)/ERROR:bad Jacobian'
    GOTO 9999

  ENDIF

  !Check whether convergence is improving
  IF(err1 < err0) THEN

    !Converging, double step if possible and set point 0 = point 1
    IF(nh > 1) nh=nh/2
    r_flx0(:)=r_flx1(:)
    g_cyl0(:)=g_cyl1(:)
    err0=err1
    tau0=tau1
    dr0=dr1
    dz0=dz1

  ELSE

    !Not converging - halve step
    nh=nh*2

  ENDIF

  !Calculate new rho and theta steps
  !Use empirical combination of last two steps to avoid oscillation
  !0 and 1 may coincide here
  g_cylt(:)=0.75*g_cyl0(:)+0.25*g_cyl1(:)
  taut=0.75*tau0+0.25*tau1
  dr=dr0/nh
  dz=dz0/nh

  !Project to desired point using Jacobian information
  !drho=(R_theta*dZ-Z_theta*dR)/tau
  drho=(g_cylt(2)*dz-g_cylt(5)*dr)/taut
  !dtheta=(Z_rho*dR-R_rho*dZ)/tau
  dtheta=(g_cylt(4)*dr-g_cylt(1)*dz)/taut
  IF(ABS(dtheta) > z_pi/4) dtheta=SIGN(z_pi/4,dtheta)

  !Set flux coordinates for new point 1
  r_flx1(1)=r_flx0(1)+drho
  r_flx1(2)=r_flx0(2)+dtheta

  IF(r_flx1(1) < 0.0) THEN

    !Stepped past minor axis, flip flux coordinates
    r_flx1(1)=-r_flx1(1)
    r_flx1(2)=r_flx1(2)+z_pi-2*dtheta

  ENDIF

ENDDO !Over iteration

!Exceeded maximum iterations
iflag=1
message='AJAX_CYL2FLX(5)/ERROR:max iterations exceeded'

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE AJAX_CYL2FLX

SUBROUTINE AJAX_FLX2CYL(r_flx, &
                        r_cyl,iflag,message, &  
                        G_CYL,GSQRT,TAU)
!-------------------------------------------------------------------------------
!AJAX_FLX2CYL converts from flux to cylindrical coordinates
!
!References:
!  S.E.Attenberger, W.A.Houlberg, S.P.Hirshman J Comp Phys 72 (1987) 435
!  W.A.Houlberg, F90 free format 8/2004
!
!Comments:
!  The representation is extended beyond the rho_3d grid by linear
!    extrapolation of the m=1, n=0 term in rho, with all other terms
!    held fixed at the edge value
!  This permits unique representation of all space for flux surfaces
!    that are everywhere convex (does not include bean-shapes)
!  rho**m is factored out of the Fourier coefficients prior to spline
!    fitting for increased accuracy near the origin
!  Mode filtering is used to improve calculations near the axis and the outer
!    boundary
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &    
  r_flx(:)               !flux coordinates (rho,theta,zeta) [rho,rad,rad]

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error

REAL(KIND=rspec), INTENT(OUT) :: &    
  r_cyl(:)               !cylindrical coordinates (R,phi,Z) [m,rad,m]

!Declaration of optional output variables
REAL(KIND=rspec), INTENT(OUT), OPTIONAL :: &   
  G_CYL(:),            & !R,Z derivatives
                         !=(R_rho,R_theta,R_zeta,Z_rho,Z_theta,Z_zeta)
                         ! [m/rho,m,      m,     m/rho,m,      m     ]
  GSQRT,               & !3D Jacobian [m**3/rho]
  TAU                    !2D Jacobian in phi=zeta=constant plane [m**2/rho]
    
!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &      
  k

INTEGER, SAVE :: &      
  i=1

INTEGER, PARAMETER :: &     
  k_vopt(1:3)=(/1,1,0/)

REAL(KIND=rspec) :: &     
  ct,st,rho,rhom,drhom,rmnx,drmnx,zmnx,dzmnx,g_cylt(1:6),value(1:3)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Null output
iflag=0
message=''
r_cyl(:)=0

!Null local values
g_cylt(:)=0
value(:)=0

!phi = zeta
r_cyl(2)=r_flx(3)

!Limit to R,Z domain
rho=MIN(r_flx(1),rhomax_3d)

!Axial resolution
rho=MAX(rho,rhores_3d)

!-------------------------------------------------------------------------------
!Evaluate R,Z and derivatives
!-------------------------------------------------------------------------------
!Loop over modes
DO k=1,krz_3d !Over modes

  !Set sine and cosine values
  ct=COS(m_3d(k)*r_flx(2)-n_3d(k)*r_flx(3))
  st=SIN(m_3d(k)*r_flx(2)-n_3d(k)*r_flx(3))

  !Calculate rho**m and its derivative
  IF(mabs_3d(k) == 0) THEN

    !rho**0
    rhom=1
    drhom=0

  ELSEIF(r_flx(1) < rhomax_3d+rhores_3d) THEN

    !Inside last closed surface, rho**|m|
    rhom=rho**mabs_3d(k)
    drhom=mabs_3d(k)*rho**(mabs_3d(k)-1)

  ELSE

    !Outside last closed surface, rhomax**|m|
    rhom=rhomax_3d**mabs_3d(k)
    drhom=0

  ENDIF

  !R_mn and dR_mn/drho from spline fits
  CALL SPLINE1_EVAL(k_vopt,nrho_3d,rho,rho_3d, &  
                     r_3d(1:4,1:nrho_3d,k),i,value)
  rmnx=value(1)
  drmnx=value(2)

  !Z_mn and dZ_mn/drho from spline fits
  CALL SPLINE1_EVAL(k_vopt,nrho_3d,rho,rho_3d,z_3d(1:4,1:nrho_3d,k),i,value)
  zmnx=value(1)
  dzmnx=value(2)

  !R = sum_mn [ R_mn * cos(m*theta-n*zeta) * rho**m ]
  r_cyl(1)=r_cyl(1)+rmnx*ct*rhom

  !Z = sum_mn [ Z_mn * sin(m*theta-n*zeta) * rho**m ]
  r_cyl(3)=r_cyl(3)+zmnx*st*rhom

  !dR/dtheta = sum_mn [ -m * R_mn * sin(m*theta-n*zeta) * rho**m ]
  g_cylt(2)=g_cylt(2)-m_3d(k)*rmnx*st*rhom

  !dR/dzeta = sum_mn [ n * R_mn * sin(m*theta-n*zeta) * rho**m ]
  g_cylt(3)=g_cylt(3)+n_3d(k)*rmnx*st*rhom

  !dZ/dtheta = sum_mn [ m * Z_mn * cos(m*theta-n*zeta) * rho**m ]
  g_cylt(5)=g_cylt(5)+m_3d(k)*zmnx*ct*rhom

  !dZ/dtheta = sum_mn [ -n * Z_mn * cos(m*theta-n*zeta) * rho**m ]
  g_cylt(6)=g_cylt(6)-n_3d(k)*zmnx*ct*rhom

  !Radial derivatives inside R,Z domain
  IF(r_flx(1) <= rhomax_3d+rhores_3d) THEN

    !dR/drho = sum_mn [ rho**m * dR_mn/drho + R_mn *d(rho**m)/drho ]
    !                   * cos(m*theta-n*zeta)
    g_cylt(1)=g_cylt(1)+(drmnx*rhom+rmnx*drhom)*ct

    !dZ/drho = sum_mn [ rho**m * dZ_mn/drho + Z_mn *d(rho**m)/drho ]
    !                   * sin(m*theta-n*zeta)
    g_cylt(4)=g_cylt(4)+(dzmnx*rhom+zmnx*drhom)*st

  ENDIF
  
  IF (ALLOCATED(rs_3d) .AND. ALLOCATED(zc_3d)) THEN !SAL 09/09/11
    !R_mn and dR_mn/drho from spline fits
    CALL SPLINE1_EVAL(k_vopt,nrho_3d,rho,rho_3d, &  
                     rs_3d(1:4,1:nrho_3d,k),i,value)
    rmnx=value(1)
    drmnx=value(2)                !SAL 09/09/11
    !Z_mn and dZ_mn/drho from spline fits
    CALL SPLINE1_EVAL(k_vopt,nrho_3d,rho,rho_3d,zc_3d(1:4,1:nrho_3d,k),i,value)
    zmnx=value(1)
    dzmnx=value(2)
    !R = sum_mn [ Rs_mn * sin(m*theta-n*zeta) * rho**m ]
    r_cyl(1)=r_cyl(1)+rmnx*st*rhom
    !Z = sum_mn [ Zc_mn * cos(m*theta-n*zeta) * rho**m ]
    r_cyl(3)=r_cyl(3)+zmnx*ct*rhom
    !dR/dtheta = sum_mn [ m * Rs_mn * cos(m*theta-n*zeta) * rho**m ]
    g_cylt(2)=g_cylt(2)+m_3d(k)*rmnx*ct*rhom
    !dR/dzeta = sum_mn [ -n * Rs_mn * cos(m*theta-n*zeta) * rho**m ]
    g_cylt(3)=g_cylt(3)-n_3d(k)*rmnx*ct*rhom
    !dZ/dtheta = sum_mn [ -m * Zc_mn * sin(m*theta-n*zeta) * rho**m ]
    g_cylt(5)=g_cylt(5)-m_3d(k)*zmnx*st*rhom
    !dZ/dtheta = sum_mn [ n * Zc_mn * sin(m*theta-n*zeta) * rho**m ]
    g_cylt(6)=g_cylt(6)+n_3d(k)*zmnx*st*rhom

    !Radial derivatives inside R,Z domain
    IF(r_flx(1) <= rhomax_3d+rhores_3d) THEN

      !dR/drho = sum_mn [ rho**m * dRs_mn/drho + Rs_mn *d(rho**m)/drho ]
      !                   * sin(m*theta-n*zeta)
      g_cylt(1)=g_cylt(1)+(drmnx*rhom+rmnx*drhom)*st

      !dZ/drho = sum_mn [ rho**m * dZc_mn/drho + Zc_mn *d(rho**m)/drho ]
      !                   * cos(m*theta-n*zeta)
      g_cylt(4)=g_cylt(4)+(dzmnx*rhom+zmnx*drhom)*ct

    ENDIF
  END IF

ENDDO !Over modes

!Radial derivatives outside R,Z domain
IF(r_flx(1) > rhomax_3d+rhores_3d) THEN

  !This point is off the grid toward the wall
  ct=COS(m_3d(km1n0_3d)*r_flx(2))
  st=SIN(m_3d(km1n0_3d)*r_flx(2))
  r_cyl(1)=r_cyl(1)+(r_flx(1)-rhomax_3d)*r_3d(1,nrho_3d,km1n0_3d)*ct
  r_cyl(3)=r_cyl(3)+(r_flx(1)-rhomax_3d)*z_3d(1,nrho_3d,km1n0_3d)*st
  g_cylt(1)=r_3d(1,nrho_3d,km1n0_3d)*ct
  g_cylt(2)=g_cylt(2)-(r_flx(1)-rhomax_3d)*r_3d(1,nrho_3d,km1n0_3d)*st & 
                      *m_3d(km1n0_3d)
  g_cylt(4)=z_3d(1,nrho_3d,km1n0_3d)*st
  g_cylt(5)=g_cylt(5)+(r_flx(1)-rhomax_3d)*z_3d(1,nrho_3d,km1n0_3d)*ct & 
                      *m_3d(km1n0_3d)
                      
  IF (ALLOCATED(rs_3d) .AND. ALLOCATED(zc_3d)) THEN !SAL 09/09/11
    r_cyl(1)=r_cyl(1)+(r_flx(1)-rhomax_3d)*rs_3d(1,nrho_3d,km1n0_3d)*st
    r_cyl(3)=r_cyl(3)+(r_flx(1)-rhomax_3d)*zc_3d(1,nrho_3d,km1n0_3d)*ct
    g_cylt(1)=g_cylt(1)+rs_3d(1,nrho_3d,km1n0_3d)*st
    g_cylt(2)=g_cylt(2)+(r_flx(1)-rhomax_3d)*rs_3d(1,nrho_3d,km1n0_3d)*ct & 
                      *m_3d(km1n0_3d)
    g_cylt(4)=g_cylt(4)+zc_3d(1,nrho_3d,km1n0_3d)*ct
    g_cylt(5)=g_cylt(5)-(r_flx(1)-rhomax_3d)*zc_3d(1,nrho_3d,km1n0_3d)*st & 
                      *m_3d(km1n0_3d)
  END IF

ENDIF

!-------------------------------------------------------------------------------
!Optional output
!-------------------------------------------------------------------------------
!R,Z derivatives
IF(PRESENT(G_CYL)) G_CYL(:)=g_cylt(:)

!2D Jacobian = dR/dtheta * dZ/drho - dR/drho * dZ/dtheta
IF(PRESENT(TAU)) TAU=g_cylt(2)*g_cylt(4)-g_cylt(1)*g_cylt(5)

!3D Jacobian = R * tau
IF(PRESENT(GSQRT)) GSQRT=r_cyl(1)*(g_cylt(2)*g_cylt(4)-g_cylt(1)*g_cylt(5))

END SUBROUTINE AJAX_FLX2CYL

SUBROUTINE AJAX_B(r_flx,r_cyl,g_cyl, &
                  iflag,message, &  
                  B_CON,B_CO,B_CYL,B_CAR,B_POL,B_TOR,B_MOD)
!-------------------------------------------------------------------------------
!AJAX_B gets components of B in various forms
!
!References:
!  S.E.Attenberger, W.A.Houlberg, S.P.Hirshman, J Comp Phys 72 (1987) 435
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &    
  r_flx(:),            & !flux coordinates (rho,theta,zeta) [rho,rad,rad]
  r_cyl(:),            & !cylindrical coordinates (R,phi,Z) [m,rad,m]
  g_cyl(:)               !R,Z derivatives
                         !=(R_rho,R_theta,R_zeta,Z_rho,Z_theta,Z_zeta)
                         ! [m/rho,m,      m,     m/rho,m,      m     ]

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error

!Declaration of optional output variables
REAL(KIND=rspec), INTENT(OUT), OPTIONAL :: &   
  B_CON(:),            & !contravariant (B^rho,B^theta,B^zeta) [T*rho/m,T/m,T/m]
  B_CO(:),             & !covariant (B_rho,B_theta,B_zeta) [T*m/rho,T*m,T*m]
  B_CYL(:),            & !cylindrical (B_R,B_phi,B_Z) [T,T,T]
  B_CAR(:),            & !Cartesian (B_x,B_y,B_z) [T,T,T]
  B_POL,               & !poloidal field [T]
  B_TOR,               & !toroidal field [T]
  B_MOD                  !|B| [T]

!-------------------------------------------------------------------------------
!Declaration of local variables
REAL(KIND=rspec) :: &     
  gsqrt,iotabar(1),phiprm(1),lam_theta,lam_zeta,b_cont(1:3), & 
  b_cylt(1:3)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Null output
iflag=0
message=''

!-------------------------------------------------------------------------------
!Set primary variables
!-------------------------------------------------------------------------------
!3D Jacobian = R * (dR/dtheta * dZ/drho - dR/drho * dZ/dtheta)
gsqrt=r_cyl(1)*(g_cyl(2)*g_cyl(4)-g_cyl(1)*g_cyl(5))

!Rotational transform and derivative of toroidal flux
CALL AJAX_MAGFLUX(1,r_flx, &
                  iflag,message, &   
                  IOTABAR_R=iotabar, &    
                  PHIPRM_R=phiprm)
                  
!Magnetic stream function derivatives
CALL AJAX_LAMBDA(r_flx, &
                 lam_theta,lam_zeta)


!Check messages
IF(iflag /= 0) THEN

  message='AJAX_B/'//message
  IF (iflag > 0 ) GOTO 9999

ENDIF

!-------------------------------------------------------------------------------
!Contravariant components
!-------------------------------------------------------------------------------
!B^rho
b_cont(1)=0

!B^theta
!b_cont(2)=phiprm(1)*(iotabar(1)-lam_zeta)/(2*z_pi*gsqrt)
b_cont(2)=phiprm(1)*(iotabar(1)-lam_zeta)/(2*z_pi*gsqrt)

!B^zeta
b_cont(3)=phiprm(1)*(1.0_rspec+lam_theta)/(2*z_pi*gsqrt)

PRINT *,'r_cyl(1:3)',r_cyl(1),r_cyl(2),r_cyl(3)
PRINT *,'r_flx(1:3)',r_flx(1),r_flx(2),r_flx(3)
PRINT *,'phiprm(1)',phiprm(1)
PRINT *,'phiprm(1)',phiprm(1),'iotabar',iotabar(1),&
        'lam_theta',lam_theta,'lam_zeta',lam_zeta,&
        'gsqrt',gsqrt
PRINT *,'b_cont(1:3)',b_cont(1),b_cont(2),b_cont(3)
PRINT *,'g_cyl(1:6)',g_cyl(1),g_cyl(2),g_cyl(3),&
                     g_cyl(4),g_cyl(5),g_cyl(6)

IF(PRESENT(B_CON)) B_CON(1:3)=b_cont(1:3)

!-------------------------------------------------------------------------------
!Cylindrical components
!-------------------------------------------------------------------------------
IF(PRESENT(B_CYL) .OR. &     
   PRESENT(B_CO) .OR. &     
   PRESENT(B_CAR) .OR. &     
   PRESENT(B_POL) .OR. &     
   PRESENT(B_TOR) .OR. &     
   PRESENT(B_MOD)) THEN

  !B_R
  b_cylt(1)=g_cyl(2)*b_cont(2)+g_cyl(3)*b_cont(3)

  !B_phi
  b_cylt(2)=r_cyl(1)*b_cont(3)

  !B_Z
  b_cylt(3)=g_cyl(5)*b_cont(2)+g_cyl(6)*b_cont(3)

  PRINT *,'b_cylt(1:3)',b_cylt(1),b_cylt(2),b_cylt(3)

  IF(PRESENT(B_CYL)) B_CYL(1:3)=b_cylt(1:3)

  !-------------------------------------------------------------------------------
  !Covariant components
  !-------------------------------------------------------------------------------
  IF(PRESENT(B_CO)) THEN

    !B_rho
    B_CO(1)=b_cylt(1)*g_cyl(1)+b_cylt(3)*g_cyl(4)

    !B_theta
    B_CO(2)=b_cylt(1)*g_cyl(2)+b_cylt(3)*g_cyl(5)

    !B_zeta
    B_CO(3)=b_cylt(1)*g_cyl(3)+b_cylt(2)*r_cyl(1)+b_cylt(3)*g_cyl(6)

  ENDIF

  !-------------------------------------------------------------------------------
  !Cartesian components
  !-------------------------------------------------------------------------------
  IF(PRESENT(B_CAR)) THEN

    !B_x
    B_CAR(1)=b_cylt(1)*COS(r_cyl(2))-b_cylt(2)*SIN(r_cyl(2))

    !B_y
    B_CAR(2)=b_cylt(1)*SIN(r_cyl(2))+b_cylt(2)*COS(r_cyl(2))

    !B_z
    B_CAR(3)=b_cylt(3)

  ENDIF

  !-------------------------------------------------------------------------------
  !Poloidal field
  !-------------------------------------------------------------------------------
  IF(PRESENT(B_POL)) THEN

    !B_pol
    B_POL=signbp_3d*SQRT(b_cylt(1)**2+b_cylt(3)**2)

  ENDIF

  !-------------------------------------------------------------------------------
  !Toroidal field
  !-------------------------------------------------------------------------------
  IF(PRESENT(B_TOR)) THEN

    !B_tor
    B_TOR=b_cylt(2)

  ENDIF

  !-------------------------------------------------------------------------------
  !|B|
  !-------------------------------------------------------------------------------
  IF(PRESENT(B_MOD)) THEN

    !B_mod
    B_MOD=SQRT(b_cylt(1)**2+b_cylt(2)**2+b_cylt(3)**2)

  ENDIF

ENDIF

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE AJAX_B

SUBROUTINE AJAX_LAMBDA(r_flx, &
                       lam_theta,lam_zeta)
!-------------------------------------------------------------------------------
!AJAX_LAMBDA gets the stream function and its derivatives
!
!References:
!  S.E.Attenberger, W.A.Houlberg, S.P.Hirshman, J Comp Phys 72 (1987) 435
!  W.A.Houlberg, F90 free format 8/2004
!
!Comments:
!  Mode filtering is used to improve calculations near the axis and the outer
!    boundary
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &    
  r_flx(:)               !flux coordinates (rho,theta,zeta) [rho,rad,rad]

!Declaration of output variables
REAL(KIND=rspec), INTENT(OUT) :: &    
  lam_theta,           & !d(lambda)/d(theta) [/rad]
  lam_zeta               !d(lambda)/d(zeta) [/rad]
     
!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER, PARAMETER :: &     
  k_vopt(1:3)=(/1,0,0/)

INTEGER, SAVE :: &      
  i=1

INTEGER  :: &      
  k

REAL(KIND=rspec) :: &     
  value(1:3),cosk,lam,rho
  
  
REAL(KIND=rspec) :: &     ! SAL 09/09/11
  sink

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Null output
lam_theta=0
lam_zeta=0

!Spline values
value(:)=0

!Limit to R,Z domain
rho=MIN(r_flx(1),rhomax_3d)

!Axial resolution
rho=MAX(rho,rhores_3d)

!-------------------------------------------------------------------------------
!Calculate derivatives of the magnetic stream function
!-------------------------------------------------------------------------------
!Loop over the modes     
DO k=1,klam_3d !Over modes

  cosk=COS(m_3d(k)*r_flx(2)-n_3d(k)*r_flx(3))
  CALL SPLINE1_EVAL(k_vopt,nrho_3d,rho,rho_3d,lam_3d(1:4,1:nrho_3d,k),i,value)
  !Reintroduce normalization
  lam=value(1)*rho**mabs_3d(k)
  lam_theta=lam_theta+lam*m_3d(k)*cosk
  lam_zeta=lam_zeta-lam*n_3d(k)*cosk
  
  IF (ALLOCATED(lamc_3d)) THEN  !SAL 09/09/11
    sink=SIN(m_3d(k)*r_flx(2)-n_3d(k)*r_flx(3))
    CALL SPLINE1_EVAL(k_vopt,nrho_3d,rho,rho_3d,lamc_3d(1:4,1:nrho_3d,k),i,value)
    !Reintroduce normalization
    lam=value(1)*rho**mabs_3d(k)
    lam_theta=lam_theta-lam*m_3d(k)*sink
    lam_zeta=lam_zeta+lam*n_3d(k)*sink
  ENDIF

ENDDO !Over modes

END SUBROUTINE AJAX_LAMBDA

SUBROUTINE AJAX_FLUXAV_B(nrho_r,rho_r, &
                         iflag,message, &  
                         B2_R,BM2_R,FM_R,FTRAP_R,GR2BM2_R,GRTH_R,SUS11_R, &
                         SUS12_R,SUS21_R,SUS22_R)
!-------------------------------------------------------------------------------
!AJAX_FLUXAV_B gets flux surface quantities that depend on B
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &     
  nrho_r                 !no. of radial nodes [-]

REAL(KIND=rspec), INTENT(IN) :: &    
  rho_r(:)               !radial nodes [rho]

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error

!Declaration of optional output variables
REAL(KIND=rspec), INTENT(OUT), OPTIONAL :: &   
  B2_R(:),             & !<B**2> [T**2]
  BM2_R(:),            & !<1/B**2> [/T**2]
  FM_R(:,:),           & !poloidal moments of <[n.grad(B)]**2>/<B**2> [-]
  FTRAP_R(:),          & !trapped particle fraction [-]
  GR2BM2_R(:),         & !<grad(rho)**2/B**2> [rho**2/m**2/T**2]
  GRTH_R(:),           & !n.grad(Theta) [/m]
  SUS11_R(:),          & !suceptance matrix element 11 (pol-pol) [-]
  SUS12_R(:),          & !suceptance matrix element 12 (pol-tor) [-]
  SUS21_R(:),          & !suceptance matrix element 21 (tor-pol) [-]
  SUS22_R(:)             !suceptance matrix element 22 (tor-tor) [-]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &      
  i,j,k,m,nr

REAL(KIND=rspec) :: &     
  dbdt,h

REAL(KIND=rspec) :: &     
  b(1:nrho_3d),b2(1:nrho_3d),bmax(1:nrho_3d), &   
  captheta(1:nrho_3d,1:ntheta_3d), &    
  gam(1:nrho_3d),v1(1:nrho_3d),v2(1:nrho_3d)

REAL(KIND=rspec) :: &     
  v_r(1:nrho_r)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Null output
iflag=0
message=''

!-------------------------------------------------------------------------------
!Check whether flux surface averaging arrays have been set
!-------------------------------------------------------------------------------
IF(.NOT. l_fluxavg_3d) THEN

  iflag=0
  message=''
  CALL AJAX_INIT_FLUXAV_G(iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_FLUXAV_B(1)/'//message
    GOTO 9999

  ENDIF

ENDIF

IF(.NOT. l_fluxavb_3d) CALL AJAX_INIT_FLUXAV_B

!-------------------------------------------------------------------------------
!Check whether values outside R,Z domain are requested (nr is last point inside)
!-------------------------------------------------------------------------------
LOOP_I1: DO i=nrho_r,1,-1

  nr=i
  IF(rho_r(nr) < rhomax_3d+rhores_3d) EXIT LOOP_I1

ENDDO LOOP_I1

!-------------------------------------------------------------------------------
!<B**2> on internal grid for b2_r, ftrap_r, and fm_r
!-------------------------------------------------------------------------------
IF(PRESENT(B2_R) .OR. &     
   PRESENT(FM_R) .OR. &     
   PRESENT(FTRAP_R)) THEN

  !Initialization
  b2(:)=0
  gt_3d(:)=0
  az_3d(:)=0

  !Calculate on internal grid 
  DO i=2,nrho_3d !Over radial nodes

    DO k=1,nzeta_3d !Over toroidal nodes

      gt_3d(1:ntheta_3d)=gsqrt_3d(i,1:ntheta_3d,k)*b_3d(i,1:ntheta_3d,k)**2
      az_3d(k)=SUM(wtheta_3d*gt_3d) !Over poloidal nodes

    ENDDO !Over toroidal nodes

    b2(i)=SUM(wzeta_3d*az_3d)/vp_3d(i) !Over toroidal nodes

  ENDDO !Over radial nodes

  !Extrapolate to axis
  b2(1)=b2(2)-rho_3d(2)*(b2(3)-b2(2))/(rho_3d(3)-rho_3d(2))

ENDIF

!-------------------------------------------------------------------------------
!Theta and gam=n.grad(Theta) for Fm and n.grad(Theta) in axisymmetric plasmas
!-------------------------------------------------------------------------------
IF((PRESENT(FM_R) .OR. &     
    PRESENT(GRTH_R)) .AND. &     
    nper_3d == 0) THEN

  !Initialization
  captheta(:,:)=0
  gam(:)=0

  !Calculate on internal grid 
  DO i=2,nrho_3d !Over radial nodes

    DO j=2,ntheta_3d !Over poloidal nodes

      captheta(i,j)=captheta(i,j-1)+dtheta_3d/2 &
                    *(b_3d(i,j-1,1)/btheta_3d(i,j-1,1) &  
                     +b_3d(i,j,1)/btheta_3d(i,j,1))

    ENDDO !Over poloidal nodes

    gam(i)=2*z_pi/captheta(i,ntheta_3d)
    captheta(i,1:ntheta_3d)=gam(i)*captheta(i,1:ntheta_3d)

  ENDDO !Over radial nodes

  !Extrapolate to axis
  gam(1)=gam(2)-rho_3d(2)*(gam(3)-gam(2))/(rho_3d(3)-rho_3d(2))

ENDIF

!-------------------------------------------------------------------------------
!<B**2>
!-------------------------------------------------------------------------------
IF(PRESENT(B2_R)) THEN

  !Initialization
  B2_R(1:nrho_r)=0

  !Interpolate to user grid inside R,Z domain
  iflag=0
  message=''
  CALL LINEAR1_INTERP(nrho_3d,rho_3d,b2,nr,rho_r,b2_r,iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_FLUXAV_B(2)/'//message
    IF(iflag > 0) GOTO 9999

  ENDIF

  IF(nr < nrho_r) THEN

    !Extrapolate to user grid outside R,Z domain
    B2_R(nr+1:nrho_r)=B2_R(nr)

  ENDIF

ENDIF

!-------------------------------------------------------------------------------
!<1/B**2>
!-------------------------------------------------------------------------------
IF(PRESENT(BM2_R)) THEN

  !Initialization
  BM2_R(1:nrho_r)=0
  gt_3d(:)=0
  az_3d(:)=0
  v1(:)=0

  !Calculate on internal grid 
  DO i=2,nrho_3d !Over radial nodes

    DO k=1,nzeta_3d !Over toroidal nodes

      gt_3d(1:ntheta_3d)=gsqrt_3d(i,1:ntheta_3d,k)/b_3d(i,1:ntheta_3d,k)**2
      az_3d(k)=SUM(wtheta_3d*gt_3d) !Over poloidal nodes

    ENDDO !Over toroidal nodes

    v1(i)=SUM(wzeta_3d*az_3d)/vp_3d(i) !Over toroidal nodes

  ENDDO !Over radial nodes

  !Extrapolate to axis
  v1(1)=v1(2)-rho_3d(2)*(v1(3)-v1(2))/(rho_3d(3)-rho_3d(2))

  !Interpolate to user grid inside R,Z domain
  iflag=0
  message=''
  CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,BM2_R,iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_FLUXAV_B(3)/'//message
    IF(iflag > 0) GOTO 9999

  ENDIF

  IF(nr < nrho_r) THEN

  !Extrapolate to user grid outside R,Z domain
  BM2_R(nr+1:nrho_r)=BM2_R(nr)

  ENDIF

ENDIF

!-------------------------------------------------------------------------------
!<grad(rho)**2/B**2>
!-------------------------------------------------------------------------------
IF(PRESENT(GR2BM2_R)) THEN

  !Initialization
  GR2BM2_R(1:nrho_r)=0
  gt_3d(:)=0
  az_3d(:)=0
  v1(:)=0

  DO i=2,nrho_3d !Over radial nodes

    DO k=1,nzeta_3d !Over toroidal nodes

      gt_3d(1:ntheta_3d)=((gcyl_3d(3,i,1:ntheta_3d,k) &  
                          *gcyl_3d(5,i,1:ntheta_3d,k) &  
                          -gcyl_3d(2,i,1:ntheta_3d,k) &  
                          *gcyl_3d(6,i,1:ntheta_3d,k))**2 & 
                          +rcyl_3d(i,1:ntheta_3d,k)**2 & 
                         *(gcyl_3d(2,i,1:ntheta_3d,k)**2 & 
                          +gcyl_3d(5,i,1:ntheta_3d,k)**2)) & 
                         /gsqrt_3d(i,1:ntheta_3d,k)/b_3d(i,1:ntheta_3d,k)**2
      az_3d(k)=SUM(wtheta_3d*gt_3d) !Over poloidal nodes

    ENDDO !Over toroidal nodes

    v1(i)=SUM(wzeta_3d*az_3d)/vp_3d(i) !Over toroidal nodes

  ENDDO !Over radial nodes

  !Extrapolate to axis
  v1(1)=v1(2)-rho_3d(2)*(v1(3)-v1(2))/(rho_3d(3)-rho_3d(2))

  !Interpolate to user grid inside R,Z domain
  iflag=0
  message=''
  CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,GR2BM2_r,iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_FLUXAV_B(4)/'//message
    IF(iflag > 0) GOTO 9999

  ENDIF

  IF(nr < nrho_r) THEN

  !Extrapolate to user grid outside R,Z domain
  GR2BM2_R(nr+1:nrho_r)=GR2BM2_R(nr)

  ENDIF

ENDIF

!-------------------------------------------------------------------------------
!f_trap
!-------------------------------------------------------------------------------
IF(PRESENT(FTRAP_R)) THEN

  !Initialization
  FTRAP_R(1:nrho_r)=0
  gt_3d(:)=0
  az_3d(:)=0
  b(:)=0
  bmax(:)=0
  v1(:)=0

  !<|B|> and B_max
  DO i=2,nrho_3d !Over radial nodes

    DO k=1,nzeta_3d !Over toroidal nodes

      DO j=1,ntheta_3d !Over poloidal nodes

        gt_3d(j)=b_3d(i,j,k)*gsqrt_3d(i,j,k)
        IF(b_3d(i,j,k) > bmax(i)) bmax(i)=b_3d(i,j,k)

      ENDDO !Over poloidal nodes

      az_3d(k)=SUM(wtheta_3d*gt_3d) !Over poloidal nodes

    ENDDO !Over toroidal nodes

    b(i)=SUM(wzeta_3d*az_3d)/vp_3d(i) !Over toroidal nodes

  ENDDO !Over radial nodes

  !Extrapolate to axis
  b(1)=b(2)-rho_3d(2)*(b(3)-b(2))/(rho_3d(3)-rho_3d(2))
  bmax(1)=bmax(2)-rho_3d(2)*(bmax(3)-bmax(2))/(rho_3d(3)-rho_3d(2))

  !Upper and lower trapped fractions bound solution
  DO i=2,nrho_3d !Over radial nodes

    !Flux surface average for upper estimate of trapped fraction
    DO k=1,nzeta_3d !Over toroidal nodes

      DO j=1,ntheta_3d !Over poloidal nodes

        h=b_3d(i,j,k)/bmax(i)
        gt_3d(j)=(1.0-SQRT(1.0-h)*(1.0+h/2))/h**2*gsqrt_3d(i,j,k)

      ENDDO !Over poloidal nodes

      az_3d(k)=SUM(wtheta_3d*gt_3d) !Over poloidal nodes

    ENDDO !Over toroidal nodes

    v1(i)=SUM(wzeta_3d*az_3d)/vp_3d(i) !Over toroidal nodes

    !Trapped fraction normalized to sqrt(rho) = 0.25*f_trap_low + 0.75*f_trap_upper
    h=b(i)/bmax(i)
    v1(i)=((1.0-b2(i)/bmax(i)**2*v1(i))/4+0.75*(1.0-b2(i)/b(i)**2 &  
          *(1.0-SQRT(1.0-h)*(1.0+h/2))))/SQRT(rho_3d(i))

  ENDDO !Over radial nodes


  !Extrapolate to axis using constant normalized trapped fraction 
  !Find index of first point away from axis in R,Z domain
  LOOP_I2: DO i=1,nrho_r !Over radial nodes

    j=i
    IF(rho_3d(j) > rhomin_3d) EXIT LOOP_I2

  ENDDO LOOP_I2 !Over radial nodes

  IF(j > 1) v1(1:j-1)=v1(j)

  !Interpolate to user grid inside R,Z domain
  iflag=0
  message=''
  CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,FTRAP_R,iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_FLUXAV_B(5)/'//message
    IF(iflag > 0) GOTO 9999

  ENDIF

  IF(nr < nrho_r) THEN

    !Extrapolate to user grid outside R,Z domain
    FTRAP_R(nr+1:nrho_r)=FTRAP_R(nr)

  ENDIF

  !Reintroduce dominant radial dependence
  FTRAP_R(1:nrho_r)=FTRAP_R(1:nrho_r)*SQRT(rho_r(1:nrho_r))

ENDIF

!-------------------------------------------------------------------------------
!F_m
!-------------------------------------------------------------------------------
IF(PRESENT(FM_R)) THEN

  !Initialization
  FM_R(:,:)=0
  gt_3d(:)=0
  v1(:)=0
  v2(:)=0
  v_r(:)=0

  IF(nper_3d /= 0) THEN

    !Non-axisymmetric plasma
    !q(rho) for fm_1 interpolation
    DO i=2,nrho_3d !Over radial nodes

      v1(i)=1/ABS(iotabar_3d(1,i))

    ENDDO !Over radial nodes

    !Extrapolate to axis
    v1(1)=v1(2)-rho_3d(2)*(v1(3)-v1(2))/(rho_3d(3)-rho_3d(2))

    !Interpolate to user grid inside R,Z domain
    iflag=0
    message=''
    CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,v_r,iflag,message)

    !Check messages
    IF(iflag /= 0) THEN

      message='AJAX_FLUXAV_B(6)/'//message
      IF(iflag > 0) GOTO 9999

    ENDIF

    IF(nr < nrho_r) THEN

      !Extrapolate to user grid outside R,Z domain
      v_r(nr+1:nrho_r)=v_r(nr)

    ENDIF

    !Reintroduce dominant radial dependence
    !If first node is not at the axis, include in calculation
    j=2
    IF(rho_r(1) > rhores_3d) j=1

    DO i=j,nrho_r !Over radial nodes

      !fm_r(2)=eps**1.5
      FM_R(2,i)=(rho_r(i)/r000_3d)**1.5

      !fm_r(1)=R_0*q/eps**1.5
      FM_R(1,i)=r000_3d*v_r(i)/fm_r(2,i)

    ENDDO !Over radial nodes

  ELSE

    !Axisymmetric plasma
    DO m=1,SIZE(FM_R,1) !Over poloidal moments

      DO i=2,nrho_3d !Over radial nodes

        !sin(Theta) terms
        DO j=1,ntheta_3d !Over poloidal nodes

          IF(j == 1 .OR. j == ntheta_3d) THEN

            !d(B)/d(theta) at end points noting periodicity
            dbdt=(b_3d(i,2,1)-b_3d(i,ntheta_3d-1,1))/2/dtheta_3d

          ELSE

            dbdt=(b_3d(i,j+1,1)-b_3d(i,j-1,1))/2/dtheta_3d

          ENDIF

          !sin(m*Theta)*B^theta/B*(dB/dtheta)
          gt_3d(j)=gsqrt_3d(i,j,1)*btheta_3d(i,j,1)/b_3d(i,j,1) & 
                   *SIN(m*captheta(i,j))*dbdt

        ENDDO !Over poloidal nodes

        v1(i)=2*z_pi*SUM(wtheta_3d*gt_3d)/vp_3d(i) !Over poloidal nodes

        DO j=1,ntheta_3d !Over poloidal nodes

          !sin(m*Theta)*B^theta*(dB/dtheta)
          gt_3d(j)=gt_3d(j)*b_3d(i,j,1)*gam(i)

        ENDDO !Over poloidal nodes

        v2(i)=v1(i)*2*z_pi*SUM(wtheta_3d*gt_3d)/vp_3d(i) !Over poloidal nodes

        !cos(Theta) terms
        DO j=1,ntheta_3d !Over poloidal nodes

          IF(j == 1 .OR. j == ntheta_3d) THEN

            !d(B)/d(theta) at end points noting periodicity
            dbdt=(b_3d(i,2,1)-b_3d(i,ntheta_3d-1,1))/(2*dtheta_3d)

          ELSE

            dbdt=(b_3d(i,j+1,1)-b_3d(i,j-1,1))/2/dtheta_3d

          ENDIF

          !cos(m*Theta)*B^theta/B*(dB/dtheta)
          gt_3d(j)=gsqrt_3d(i,j,1)*btheta_3d(i,j,1)/b_3d(i,j,1) & 
                   *COS(m*captheta(i,j))*dbdt

        ENDDO !Over poloidal nodes

        v1(i)=2*z_pi*SUM(wtheta_3d*gt_3d) &  
              /vp_3d(i) !Over poloidal nodes

        !cos(m*Theta)*B^theta*(dB/dtheta)
        gt_3d(:)=gt_3d(:)*b_3d(i,:,1)*gam(i)
        v2(i)=v2(i)+v1(i)*2*z_pi*SUM(wtheta_3d*gt_3d)/vp_3d(i) !Over poloidal nodes

        !<B^theta>
        gt_3d(:)=gsqrt_3d(i,:,1)*btheta_3d(i,:,1)
        v1(i)=2*z_pi*SUM(wtheta_3d*gt_3d)/vp_3d(i) !Over poloidal nodes

        !Multiply by 2/<B^theta>/<B**2> and normalize to rho**1.5
        v2(i)=v2(i)*2/v1(i)/b2(i)/rho_3d(i)**1.5

      ENDDO !Over radial nodes

      !Extrapolate to axis
      v2(1)=v2(2)-rho_3d(2)*(v2(3)-v2(2))/(rho_3d(3)-rho_3d(2))

      !Interpolate to user grid inside R,Z domain
      iflag=0
      message=''
      CALL LINEAR1_INTERP(nrho_3d,rho_3d,v2,nr,rho_r,v_r,iflag,message)

      !Check messages
      IF(iflag /= 0) THEN

        message='AJAX_FLUXAV_B(7)/'//message
        IF(iflag > 0) GOTO 9999

      ENDIF

      IF(nr < nrho_r) THEN

        !Extrapolate to user grid outside R,Z domain
        v_r(nr+1:nrho_r)=v_r(nr)

      ENDIF

      !Reintroduce dominant radial dependence
      FM_R(m,1:nrho_r)=v_r(1:nrho_r)*rho_r(1:nrho_r)**1.5

    ENDDO !Over poloidal moments

  ENDIF

ENDIF

!-------------------------------------------------------------------------------
!n.grad(Theta)
!-------------------------------------------------------------------------------
IF(PRESENT(GRTH_R)) THEN

  !Initialization
  GRTH_R(1:nrho_r)=0

  !Interpolate to user grid inside R,Z domain
  IF(nper_3d == 0) THEN

    iflag=0
    message=''
    CALL LINEAR1_INTERP(nrho_3d,rho_3d,gam,nr,rho_r,GRTH_R,iflag,message)

    !Check messages
    IF(iflag /= 0) THEN

      message='AJAX_FLUXAV_B(8)/'//message
      IF(iflag > 0) GOTO 9999

    ENDIF

  ENDIF

  IF(nr < nrho_r) THEN

    !Extrapolate to user grid outside R,Z domain
    GRTH_R(nr+1:nrho_r)=GRTH_R(nr)

  ENDIF

ENDIF

!-------------------------------------------------------------------------------
!Susceptance matrix element S11
!-------------------------------------------------------------------------------
IF(PRESENT(SUS11_R)) THEN

  !Initialization
  SUS11_R(1:nrho_r)=0
  gt_3d(:)=0
  az_3d(:)=0
  v1(:)=0

  !Calculate on internal grid
  DO i=2,nrho_3d !Over radial nodes

    DO k=1,nzeta_3d !Over toroidal nodes

      gt_3d(1:ntheta_3d)=( gcyl_3d(2,i,1:ntheta_3d,k)**2 & 
                          +gcyl_3d(5,i,1:ntheta_3d,k)**2 ) & 
                         /gsqrt_3d(i,1:ntheta_3d,k)
      az_3d(k)=SUM(wtheta_3d*gt_3d) !Over poloidal nodes

    ENDDO !Over toroidal nodes

    !Remove dominant radial dependence
    v1(i)=SUM(wzeta_3d*az_3d)/rho_3d(i) !Over toroidal nodes

  ENDDO !Over radial nodes

  !Extrapolate to axis
  v1(1)=v1(2)-rho_3d(2)*(v1(3)-v1(2))/(rho_3d(3)-rho_3d(2))

  !Interpolate to user grid inside R,Z domain
  iflag=0
  message=''
  CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,SUS11_R,iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_FLUXAV_B(9)/'//message
    IF(iflag > 0) GOTO 9999

  ENDIF

  IF(nr < nrho_r) THEN

    !Extrapolate to user grid outside R,Z domain
    SUS11_R(nrho_r+1:nr)=SUS11_R(nr)

  ENDIF

  !Reintroduce dominant radial dependence
  SUS11_R(1:nrho_r)=SUS11_R(1:nrho_r)*rho_r(1:nrho_r)/4/z_pi**2

ENDIF

!-------------------------------------------------------------------------------
!Susceptance matrix element S12
!-------------------------------------------------------------------------------
IF(PRESENT(SUS12_R)) THEN

  !Initialization
  SUS12_R(1:nrho_r)=0
  gt_3d(:)=0
  az_3d(:)=0
  v1(:)=0

  !Calculate on internal grid
  DO i=2,nrho_3d !Over radial nodes

    DO k=1,nzeta_3d !Over toroidal nodes

      gt_3d(1:ntheta_3d)=( ( gcyl_3d(2,i,1:ntheta_3d,k) & 
                            *gcyl_3d(3,i,1:ntheta_3d,k) &  
                            +gcyl_3d(5,i,1:ntheta_3d,k) &  
                            *gcyl_3d(6,i,1:ntheta_3d,k)) & 
                           *(1.0+eltheta_3d(i,1:ntheta_3d,k)) & 
                          -( gcyl_3d(2,i,1:ntheta_3d,k)**2 & 
                            +gcyl_3d(5,i,1:ntheta_3d,k)**2 ) & 
                           *elzeta_3d(i,1:ntheta_3d,k) ) & 
                         /gsqrt_3d(i,1:ntheta_3d,k)
      az_3d(k)=SUM(wtheta_3d*gt_3d) !Over poloidal nodes

    ENDDO !Over toroidal nodes

    !Remove dominant radial dependence
    v1(i)=SUM(wzeta_3d*az_3d)/rho_3d(i) !Over toroidal nodes

  ENDDO !Over radial nodes

  !Extrapolate to axis
  v1(1)=v1(2)-rho_3d(2)*(v1(3)-v1(2))/(rho_3d(3)-rho_3d(2))

  !Interpolate to user grid inside R,Z domain
  iflag=0
  message=''
  CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,SUS12_R,iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_FLUXAV_B(10)/'//message
    IF(iflag > 0) GOTO 9999

  ENDIF

  IF(nr < nrho_r) THEN

    !Extrapolate to user grid outside R,Z domain
    sus12_r(nr+1:nrho_r)=SUS12_R(nr)

  ENDIF

  !Reintroduce dominant radial dependence
  SUS12_R(1:nrho_r)=SUS12_R(1:nrho_r)*rho_r(1:nrho_r)/4/z_pi**2

ENDIF

!-------------------------------------------------------------------------------
!Susceptance matrix element S21
!-------------------------------------------------------------------------------
IF(PRESENT(SUS21_R)) THEN

  !Initialization
  SUS21_R(1:nrho_r)=0
  gt_3d(:)=0
  az_3d(:)=0
  v1(:)=0

  !Calculate on internal grid
  DO i=2,nrho_3d !Over radial nodes

    DO k=1,nzeta_3d !Over toroidal nodes

      gt_3d(1:ntheta_3d)=( gcyl_3d(2,i,1:ntheta_3d,k) &  
                          *gcyl_3d(3,i,1:ntheta_3d,k) &  
                          +gcyl_3d(5,i,1:ntheta_3d,k) &  
                          *gcyl_3d(6,i,1:ntheta_3d,k) ) &  
                         /gsqrt_3d(i,1:ntheta_3d,k)
      az_3d(k)=SUM(wtheta_3d*gt_3d) !Over poloidal nodes

    ENDDO !Over toroidal nodes

    !Remove dominant radial dependence
    v1(i)=SUM(wzeta_3d*az_3d)/rho_3d(i) !Over toroidal nodes

  ENDDO !Over radial nodes

  !Extrapolate to axis
  v1(1)=v1(2)-rho_3d(2)*(v1(3)-v1(2))/(rho_3d(3)-rho_3d(2))

  !Interpolate to user grid inside R,Z domain
  iflag=0
  message=''
  CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,SUS21_R,iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_FLUXAV_B(11)/'//message
    IF(iflag > 0) GOTO 9999

  ENDIF

  IF(nr < nrho_r) THEN

    !Extrapolate to user grid outside R,Z domain
    SUS21_R(nr+1:nrho_r)=SUS21_R(nr)

  ENDIF

  !Reintroduce dominant radial dependence
  SUS21_R(1:nrho_r)=SUS21_R(1:nrho_r)*rho_r(1:nrho_r)/4/z_pi**2

ENDIF

!-------------------------------------------------------------------------------
!Susceptance matrix element S22
!-------------------------------------------------------------------------------
IF(PRESENT(SUS22_R)) THEN

  !Initialization
  SUS22_R(1:nrho_r)=0
  gt_3d(:)=0
  az_3d(:)=0
  v1(:)=0

  !Calculate on internal grid
  DO i=2,nrho_3d !Over radial nodes

    DO k=1,nzeta_3d !Over toroidal nodes

      gt_3d(1:ntheta_3d)=( ( rcyl_3d(i,1:ntheta_3d,k)**2 & 
                            +gcyl_3d(3,i,1:ntheta_3d,k)**2 & 
                            +gcyl_3d(6,i,1:ntheta_3d,k)**2 ) & 
                           *(1.0+eltheta_3d(i,1:ntheta_3d,k)) &
                           -( gcyl_3d(2,i,1:ntheta_3d,k) & 
                             *gcyl_3d(3,i,1:ntheta_3d,k) & 
                             +gcyl_3d(5,i,1:ntheta_3d,k) & 
                             *gcyl_3d(6,i,1:ntheta_3d,k) ) & 
                            *elzeta_3d(i,1:ntheta_3d,k) ) & 
                         /gsqrt_3d(i,1:ntheta_3d,k)
      az_3d(k)=SUM(wtheta_3d*gt_3d) !Over poloidal nodes

    ENDDO !Over toroidal nodes

    !Remove dominant radial dependence
    v1(i)=SUM(wzeta_3d*az_3d)*rho_3d(i) !Over toroidal nodes

  ENDDO !Over radial nodes

  !Extrapolate to axis
  v1(1)=v1(2)-rho_3d(2)*(v1(3)-v1(2))/(rho_3d(3)-rho_3d(2))

  !Interpolate to user grid inside R,Z domain
  iflag=0
  message=''
  CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,SUS22_R,iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_FLUXAV_B(12)/'//message
    IF(iflag > 0) GOTO 9999

  ENDIF

  IF(nr < nrho_r) THEN

    !Extrapolate to user grid outside R,Z domain
    SUS22_R(nr+1:nrho_r)=SUS22_R(nr)

  ENDIF

  !Reintroduce dominant radial dependence
  j=2
  !If first node is not at the axis, include in calculation
  IF(rho_r(1) > rhores_3d) j=1

  SUS22_R(j:nrho_r)=SUS22_R(j:nrho_r)/rho_r(j:nrho_r)/4/z_pi**2

ENDIF

IF(j == 2) SUS22_R(1)=SUS22_R(2)

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE AJAX_FLUXAV_B

SUBROUTINE AJAX_FLUXAV_G(nrho_r,rho_r, &
                         iflag,message, &  
                         AREA_R,DVOL_R,GRHO1_R,GRHO2_R,GRHO2RM2_R,RM2_R,VOL_R, &
                         VP_R)
!-------------------------------------------------------------------------------
!AJAX_FLUXAV_G gets flux surface quantities that depend on geometry
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &     
  nrho_r                 !no. of radial nodes [-]

REAL(KIND=rspec), INTENT(IN) :: &    
  rho_r(:)               !radial nodes [rho]

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error

!Declaration of optional output variables
REAL(KIND=rspec), INTENT(OUT), OPTIONAL :: &   
  AREA_R(:),           & !surface area [m**2]
  DVOL_R(:),           & !cell volume between rho_r(i) and rho_r(i+1) [m**3]
  GRHO1_R(:),          & !<|grad(rho)|> [rho/m]
  GRHO2_R(:),          & !<grad(rho)**2> [rho**2/m**4]
  GRHO2RM2_R(:),       & !<grad(rho)**2/R**2>) [rho**2/m**2]
  RM2_R(:),            & !<1/R**2> [/m**2]
  VOL_R(:),            & !enclosed volume [-]
  VP_R(:)                !dV/drho [m**3/rho]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &      
  i,k,nr

REAL(KIND=rspec) :: &     
  gr(1:nrho_r),v1(1:nrho_3d),vp(1:nrho_r)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Null output
iflag=0
message=''

!Check whether flux surface averaging arrays have been set
IF(.NOT. l_fluxavg_3d) THEN

  CALL AJAX_INIT_FLUXAV_G(iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_FLUXAV_G(1)/'//message
    GOTO 9999

  ENDIF

ENDIF

!Check whether values outside R,Z domain are requested (nr is last point inside)
LOOP_I: DO i=nrho_r,1,-1 !Over nodes

  nr=i
  IF(rho_r(nr) < rhomax_3d+rhores_3d) EXIT LOOP_I

ENDDO LOOP_I !Over nodes

!-------------------------------------------------------------------------------
!d(V)/d(rho) on user grid for area_r, dvol_r, vol_r and vp_r
!-------------------------------------------------------------------------------
IF(PRESENT(AREA_R) .OR. &     
   PRESENT(DVOL_R) .OR. &     
   PRESENT(VOL_R) .OR. &     
   PRESENT(VP_R)) THEN

  !Initialization
  vp(:)=0
  v1(:)=0

  !Remove dominant radial dependence
  v1(2:nrho_3d)=vp_3d(2:nrho_3d)/rho_3d(2:nrho_3d)

  !Extrapolate to axis
  v1(1)=v1(2)-rho_3d(2)*(v1(3)-v1(2))/(rho_3d(3)-rho_3d(2))

  !Interpolate to user grid inside R,Z domain
  iflag=0
  message=''
  CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,vp,iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_FLUXAV_G(2)/'//message
    IF(iflag > 0) GOTO 9999

  ENDIF

  IF(nr < nrho_r) THEN

  !Extrapolate to user grid outside R,Z domain
  vp(nr+1:nrho_r)=vp(nr)

  ENDIF

  !Reintroduce dominant radial dependence
  vp(1:nrho_r)=vp(1:nrho_r)*rho_r(1:nrho_r)

ENDIF

!-------------------------------------------------------------------------------
!<|grad(rho)|> on user grid for area_r and grho1_r
!-------------------------------------------------------------------------------
IF(PRESENT(AREA_R) .OR. &     
   PRESENT(GRHO1_R)) THEN

  !Initialization
  gr(:)=0
  gt_3d(:)=0
  az_3d(:)=0
  v1(:)=0

  !Calculate on internal grid
  DO i=2,nrho_3d !Over radial nodes

    DO k=1,nzeta_3d !Over toroidal nodes

      gt_3d(1:ntheta_3d)=SQRT(( ( gcyl_3d(3,i,1:ntheta_3d,k) & 
                                 *gcyl_3d(5,i,1:ntheta_3d,k) & 
                                 -gcyl_3d(2,i,1:ntheta_3d,k) & 
                                 *gcyl_3d(6,i,1:ntheta_3d,k))**2 & 
                                +rcyl_3d(i,1:ntheta_3d,k)**2 & 
                                *(gcyl_3d(2,i,1:ntheta_3d,k)**2 & 
                                 +gcyl_3d(5,i,1:ntheta_3d,k)**2) ))
      az_3d(k)=SUM(wtheta_3d*gt_3d) !Over poloidal nodes

    ENDDO !Over toroidal nodes

    v1(i)=SUM(wzeta_3d*az_3d)/vp_3d(i) !Over toroidal nodes

  ENDDO !Over radial nodes

  v1(1)=v1(2)

  !Interpolate to user grid inside R,Z domain
  iflag=0
  message=''
  CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,gr,iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_FLUXAV_G(3)/'//message
    IF(iflag > 0) GOTO 9999

  ENDIF

  IF(nr < nrho_r) THEN

    !Extrapolate to user grid outside R,Z domain
    gr(nr+1:nrho_r)=gr(nr)

  ENDIF

ENDIF

!-------------------------------------------------------------------------------
!Area
!-------------------------------------------------------------------------------
IF(PRESENT(AREA_R)) THEN

  !Initialization
  AREA_R(:)=0

  !Use temporary arrays
  AREA_R(1:nrho_r)=vp(1:nrho_r)*gr(1:nrho_r)

ENDIF

!-------------------------------------------------------------------------------
!delta(Volume)
!-------------------------------------------------------------------------------
IF(PRESENT(DVOL_R)) THEN

  !Initialization
  DVOL_R(:)=0

  !Calculate from integral rho*drho*(V'/rho), where (V'/rho) is a weak function
  !Allow for first node to be away from the axis
  IF(rho_r(1) > rhores_3d) THEN

    !First node is off axis
    k=1

  ELSE

    !First node is on axis
    k=2
    DVOL_R(1)=vp(2)*rho_r(2)/2

  ENDIF

  DO i=k,nrho_r-1 !Over radial nodes

    DVOL_R(i)=(vp(i)/rho_r(i)+vp(i+1)/rho_r(i+1))/2 &
              *(rho_r(i+1)**2-rho_r(i)**2)/2

  ENDDO !Over radial nodes

  !Set ghost node value equal to last node
  DVOL_R(nrho_r)=DVOL_R(nrho_r-1)

ENDIF

!-------------------------------------------------------------------------------
!<|grad(rho)|>
!-------------------------------------------------------------------------------
IF(PRESENT(GRHO1_R)) THEN

  !Initialization
  GRHO1_R(:)=0

  !Copy from temporary array
  GRHO1_R(1:nrho_r)=gr(1:nrho_r)

ENDIF

!-------------------------------------------------------------------------------
!<grad(rho)**2>
!-------------------------------------------------------------------------------
IF(PRESENT(GRHO2_R)) THEN

  !Initialization
  GRHO2_R(:)=0
  gt_3d(:)=0
  az_3d(:)=0
  v1(:)=0

  !Calculate on internal grid
  DO i=2,nrho_3d !Over radial nodes

    DO k=1,nzeta_3d !Over toroidal nodes

      gt_3d(1:ntheta_3d)=( (gcyl_3d(3,i,1:ntheta_3d,k) &  
                           *gcyl_3d(5,i,1:ntheta_3d,k) &  
                           -gcyl_3d(2,i,1:ntheta_3d,k) &  
                           *gcyl_3d(6,i,1:ntheta_3d,k))**2 & 
                           +rcyl_3d(i,1:ntheta_3d,k)**2 &  
                           *(gcyl_3d(2,i,1:ntheta_3d,k)**2 & 
                           +gcyl_3d(5,i,1:ntheta_3d,k)**2) ) & 
                         /gsqrt_3d(i,1:ntheta_3d,k)
      az_3d(k)=SUM(wtheta_3d*gt_3d) !Over poloidal nodes

    ENDDO !Over toroidal nodes

    v1(i)=SUM(wzeta_3d*az_3d)/vp_3d(i) !Over toroidal nodes

  ENDDO !Over radial nodes

  !Extrapolate to axis
  v1(1)=v1(2)

  !Interpolate to user grid inside R,Z domain
  iflag=0
  message=''
  CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,GRHO2_R,iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_FLUXAV_G(4)/'//message
    IF(iflag > 0) GOTO 9999

  ENDIF

  IF(nr < nrho_r) THEN

    !Extrapolate to user grid outside R,Z domain
    GRHO2_R(nr+1:nrho_r)=GRHO2_R(nr)

  ENDIF

ENDIF

!-------------------------------------------------------------------------------
!<grad(rho)**2/R**2>
!-------------------------------------------------------------------------------
IF(PRESENT(GRHO2RM2_R)) THEN

  !Initialization
  GRHO2RM2_r(:)=0
  gt_3d(:)=0
  az_3d(:)=0
  v1(:)=0

  !Calculate on internal grid
  DO i=2,nrho_3d !Over radial nodes

    DO k=1,nzeta_3d !Over toroidal nodes

      gt_3d(1:ntheta_3d)=( (gcyl_3d(3,i,1:ntheta_3d,k) &  
                           *gcyl_3d(5,i,1:ntheta_3d,k) &  
                           -gcyl_3d(2,i,1:ntheta_3d,k) &  
                           *gcyl_3d(6,i,1:ntheta_3d,k))**2 & 
                           /rcyl_3d(i,1:ntheta_3d,k)**2 &  
                           +(gcyl_3d(2,i,1:ntheta_3d,k)**2 & 
                           +gcyl_3d(5,i,1:ntheta_3d,k)**2) ) & 
                         /gsqrt_3d(i,1:ntheta_3d,k)
      az_3d(k)=SUM(wtheta_3d*gt_3d) !Over poloidal nodes

    ENDDO !Over toroidal nodes

    v1(i)=SUM(wzeta_3d*az_3d)/vp_3d(i) !Over toroidal nodes

  ENDDO !Over radial nodes

  !Extrapolate to axis
  v1(1)=v1(2)

  !Interpolate to user grid inside R,Z domain
  iflag=0
  message=''
  CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,GRHO2RM2_R,iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_FLUXAV_G(5)/'//message
    IF(iflag > 0) GOTO 9999

  ENDIF

  IF(nr < nrho_r) THEN

    !Extrapolate to user grid outside R,Z domain
    GRHO2RM2_R(nr+1:nrho_r)=GRHO2RM2_R(nr)

  ENDIF

ENDIF

!-------------------------------------------------------------------------------
!<1/R**2>
!-------------------------------------------------------------------------------
IF(PRESENT(RM2_R)) THEN

  !Initialization
  RM2_R(:)=0
  gt_3d(:)=0
  az_3d(:)=0
  v1(:)=0

  !Calculate on internal grid
  DO i=2,nrho_3d !Over radial nodes

    DO k=1,nzeta_3d !Over toroidal nodes

      gt_3d(1:ntheta_3d)=gsqrt_3d(i,1:ntheta_3d,k) &  
                         /rcyl_3d(i,1:ntheta_3d,k)**2
      az_3d(k)=SUM(wtheta_3d*gt_3d) !Over poloidal nodes

    ENDDO !Over toroidal nodes

    v1(i)=SUM(wzeta_3d*az_3d)/vp_3d(i) !Over toroidal nodes

  ENDDO !Over radial nodes

  !Extrapolate to axis
  v1(1)=v1(2)

  !Interpolate to user grid inside R,Z domain
  iflag=0
  message=''
  CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,RM2_R,iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_FLUXAV_G(6)/'//message
    IF(iflag > 0) GOTO 9999

  ENDIF

  IF(nr < nrho_r) THEN

    !Extrapolate to user grid outside R,Z domain
    RM2_R(nr+1:nrho_r)=RM2_R(nr)

  ENDIF

ENDIF

!-------------------------------------------------------------------------------
!Volume
!-------------------------------------------------------------------------------
IF(PRESENT(VOL_R)) THEN

  !Initialization
  VOL_R(:)=0

  !Sum up dvol elements
  !Allow for first node to be away from the axis
  IF(rho_r(1) > rhores_3d) THEN

    !First node is off axis
    k=1
    vol_r(1)=vp(1)*rho_r(1)/2

  ELSE

    !First node is on axis
    k=2
    VOL_R(1)=0
    VOL_R(2)=vp(2)*rho_r(2)/2

  ENDIF

  DO i=k+1,nrho_r !Over radial nodes

    VOL_R(i)=VOL_R(i-1)+(vp(i-1)/rho_r(i-1)+vp(i)/rho_r(i))/2 & 
                        *(rho_r(i)**2-rho_r(i-1)**2)/2

  ENDDO !Over radial nodes

ENDIF

!-------------------------------------------------------------------------------
!d(Volume)/d(rho)
!-------------------------------------------------------------------------------
IF(PRESENT(VP_R)) THEN

  !Initialization
  VP_R(:)=0

  !Copy from temporary array
  VP_R(1:nrho_r)=vp(1:nrho_r)

ENDIF

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE AJAX_FLUXAV_G

SUBROUTINE AJAX_I(nrho_r,rho_r, &
                  iflag,message, &   
                  CUR_I_R,CUR_F_R,CUR_IMN_R,CUR_IMX_R,CUR_FMN_R,CUR_FMX_R)
!-------------------------------------------------------------------------------
!AJAX_I gets the enclosed toroidal and external poloidal currents for a set of
!  flux surfaces from the covariant components of B, as well as maximum and
!  minimum values for use as a diagnostic on the accuracy of the stream function
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &     
  nrho_r                 !no. of radial nodes [-]

REAL(KIND=rspec), INTENT(IN) :: &    
  rho_r(:)               !radial nodes [rho]

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error

!Declaration of optional output variables
REAL(KIND=rspec), INTENT(OUT), OPTIONAL :: &   
  CUR_I_R(:),          & !enclosed toroidal current [A]
  CUR_F_R(:),          & !external poloidal current [A]
  CUR_IMN_R(:),        & !minimum toroidal current on AJAX toroidal grid [A]
  CUR_IMX_R(:),        & !maximum toroidal current on AJAX toroidal grid [A]
  CUR_FMN_R(:),        & !minimum poloidal current on AJAX poloidal grid [A]
  CUR_FMX_R(:)           !maximum poloidal current on AJAX poloidal grid [A]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &      
  i,j,k,imin

REAL(KIND=rspec) :: &     
  cmax,cmin,r_flx(1:3),r_cyl(1:3),g_cyl(1:6)

REAL(KIND=rspec), ALLOCATABLE :: &    
  b_co(:,:,:),gt(:),gz(:)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Null output
iflag=0
message=''

!Allocate covariant B field array
ALLOCATE(b_co(3,ntheta_3d,nzeta_3d))

  b_co(:,:,:)=0

IF(PRESENT(CUR_I_R) .OR. &     
   PRESENT(CUR_IMN_R) .OR. &     
   PRESENT(CUR_IMX_R)) THEN

  ALLOCATE(gz(nzeta_3d))

    gz(:)=0

ENDIF

IF(PRESENT(CUR_F_R) .OR. &     
   PRESENT(CUR_FMN_R) .OR. &     
   PRESENT(CUR_FMX_R)) THEN

  ALLOCATE(gt(ntheta_3d))

    gt(:)=0

ENDIF

IF(PRESENT(CUR_I_R)) CUR_I_R(1:nrho_r)=0
IF(PRESENT(CUR_IMN_R)) CUR_IMN_R(1:nrho_r)=0
IF(PRESENT(CUR_IMX_R)) CUR_IMX_R(1:nrho_r)=0
IF(PRESENT(CUR_F_R)) CUR_F_R(1:nrho_r)=0
IF(PRESENT(CUR_FMN_R)) CUR_FMN_R(1:nrho_r)=0
IF(PRESENT(CUR_FMX_R)) CUR_FMX_R(1:nrho_r)=0

imin=1
IF(rho_r(1) <= rhomin_3d) imin=2

!-------------------------------------------------------------------------------
!Perform integrals of covariant B components on each surface
!-------------------------------------------------------------------------------
DO i=imin,nrho_r !Over radial nodes

  r_flx(1)=rho_r(i)

  DO k=1,nzeta_3d !Over toroidal nodes

    r_flx(3)=(k-1)*dzeta_3d

    DO j=1,ntheta_3d !Over poloidal nodes

      r_flx(2)=(j-1)*dtheta_3d

      !Get cylindrical coordinates and metrics
      iflag=0
      message=''
      CALL AJAX_FLX2CYL(r_flx,r_cyl,iflag,message,G_CYL=g_cyl)

      !Check messages
      IF(iflag /= 0) THEN

        message='AJAX_I(1)/'//message
        IF(iflag > 0) GOTO 9999

      ENDIF

      !Get covariant B components
      iflag=0
      message=''
      CALL AJAX_B(r_flx,r_cyl,g_cyl,iflag,message,B_CO=b_co(:,j,k))

      !Check messages
      IF(iflag /= 0) THEN

        message='AJAX_I(2)/'//message
        IF(iflag > 0) GOTO 9999

      ENDIF

    ENDDO !Over poloidal nodes

  ENDDO !Over toroidal nodes

  !Toroidal current
  IF(PRESENT(CUR_I_R) .OR. &     
     PRESENT(CUR_IMN_R) .OR. &    
     PRESENT(CUR_IMX_R)) THEN

    cmax=-z_large
    cmin=z_large

    DO k=1,nzeta_3d !Over toroidal nodes

      gz(k)=SUM(wtheta_3d*b_co(2,1:ntheta_3d,k))/z_mu0
      IF(gz(k) < cmin) cmin=gz(k)
      IF(gz(k) > cmax) cmax=gz(k)

      IF(PRESENT(CUR_I_R)) CUR_I_R(i)=SUM(wzeta_3d*gz)/2/z_pi
      IF(PRESENT(CUR_IMN_R)) CUR_IMN_R(i)=cmin
      IF(PRESENT(CUR_IMX_R)) CUR_IMX_R(i)=cmax

    ENDDO !Over toroidal nodes

  ENDIF

  !Poloidal current
  IF(PRESENT(CUR_F_R) .OR. &     
     PRESENT(CUR_FMN_R) .OR. &    
     PRESENT(CUR_FMX_R)) THEN

    cmax=-z_large
    cmin=z_large

    DO j=1,ntheta_3d !Over poloidal nodes

      gt(j)=SUM(wzeta_3d*b_co(3,j,1:nzeta_3d))/z_mu0
      IF(gt(j) < cmin) cmin=gt(j)
      IF(gt(j) > cmax) cmax=gt(j)

      IF(PRESENT(CUR_F_R)) CUR_F_R(i)=SUM(wtheta_3d*gt)/2/z_pi
      IF(PRESENT(CUR_FMN_R)) CUR_FMN_R(i)=cmin
      IF(PRESENT(CUR_FMX_R)) CUR_FMX_R(i)=cmax

    ENDDO !Over toroidal nodes

  ENDIF

ENDDO !Over radial nodes

!Check whether axial values are requested
IF(imin == 2) THEN

  !Use value of second node
  IF(PRESENT(CUR_I_R)) CUR_I_R(1)=0
  IF(PRESENT(CUR_IMN_R)) CUR_IMN_R(1)=0
  IF(PRESENT(CUR_IMX_R)) CUR_IMX_R(1)=0
  IF(PRESENT(CUR_F_R)) CUR_F_R(1)=CUR_F_R(2)
  IF(PRESENT(CUR_FMN_R)) CUR_FMN_R(1)=CUR_FMN_R(2)
  IF(PRESENT(CUR_FMX_R)) CUR_FMX_R(1)=CUR_FMX_R(2)

ENDIF

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

!Deallocate arrays
DEALLOCATE(b_co)
IF(ALLOCATED(gt)) DEALLOCATE(gt)
IF(ALLOCATED(gz)) DEALLOCATE(gz)

END SUBROUTINE AJAX_I

SUBROUTINE AJAX_MAGFLUX(nrho_r,rho_r, &
                        iflag,message, &  
                        IOTABAR_R,Q_R,PHIPRM_R,PSIPRM_R,PHI_R,PSI_R)
!-------------------------------------------------------------------------------
!AJAX_MAGFLUX gets poloidal and toroidal magnetic fluxes
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &     
  nrho_r                 !no. of radial nodes [-]

REAL(KIND=rspec), INTENT(IN) :: &    
  rho_r(:)               !radial nodes [rho]

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error

!Declaration of output variables
REAL(KIND=rspec), INTENT(OUT), OPTIONAL :: &   
  IOTABAR_R(:),        & !rotational transform = d(Psi)/d(Phi) [-]
  Q_R(:),              & !safety factor = d(Phi)/d(Psi) [-]
  PHIPRM_R(:),         & !d(Phi)/d(rho) [Wb/rho]
  PSIPRM_R(:),         & !d(Psi)/d(rho) [Wb/rho]
  PHI_R(:),            & !toroidal flux [Wb]
  PSI_R(:)               !poloidal flux [Wb]

!-------------------------------------------------------------------------------
!Declaration of local variables             
INTEGER :: &      
  i,j,nr

INTEGER, PARAMETER :: &     
  k_vopt(1:3)=(/1,0,0/)

REAL(KIND=rspec) :: &     
  iotabar_t(1:nrho_r),phiprm_t(1:nrho_r),value(1:3), &  
  values(1:nrho_r)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Null output
iflag=0
message=''

!Internal values
phiprm_t(:)=0
iotabar_t(:)=0

!Check whether values outside R,Z domain are requested (nr is last point inside)
LOOP_I: DO i=nrho_r,1,-1 !Over nodes

  nr=i
  IF(rho_r(nr) < rhomax_3d+rhores_3d) EXIT LOOP_I

ENDDO LOOP_I !Over nodes

!-------------------------------------------------------------------------------
!Get temporary base quantities for phiprm and iotabar
!-------------------------------------------------------------------------------
phiprm_t(1:nrho_r)=2*rho_r(1:nrho_r)*phitot_3d/rhomax_3d**2
IF(rho_r(1) < rhores_3d) phiprm_t(1)=2*rhores_3d*phitot_3d/rhomax_3d**2

  !Inside MHD solution domain
  j=1

DO i=1,nr !Over radial nodes

  CALL SPLINE1_EVAL(k_vopt,nrho_3d,rho_r(i),rho_3d,iotabar_3d,j,value)
  iotabar_t(i)=value(1)

ENDDO !Over radial nodes

IF(nr < nrho_r) THEN

  !Extrapolate outside with constant slope in iotabar
  iotabar_t(nr+1:nrho_r)=iotabar_3d(1,nrho_3d) &  
                         +(rho_r(nr+1:nrho_r)-rhomax_3d) & 
                         *(iotabar_3d(1,nrho_3d) &  
                         -iotabar_3d(1,nrho_3d-1)) &  
                         /(rhomax_3d-rho_3d(nrho_3d-1))

ENDIF

!-------------------------------------------------------------------------------
!Return the appropriate quantity
!-------------------------------------------------------------------------------
!Rotational transform
IF(PRESENT(IOTABAR_R)) THEN

  IOTABAR_R(:)=0
  IOTABAR_R(1:nrho_r)=iotabar_t(1:nrho_r)

ENDIF

!Safety factor
IF(PRESENT(Q_R)) THEN

  Q_R(:)=0
  Q_R(1:nrho_r)=1/iotabar_t(1:nrho_r)

ENDIF

!d(Phi)/d(rho)
IF(PRESENT(PHIPRM_R)) THEN

  PHIPRM_R(:)=0
  PHIPRM_R(1:nrho_r)=phiprm_t(1:nrho_r)

ENDIF

!d(Psi)/d(rho)
IF(PRESENT(PSIPRM_R)) THEN

  PSIPRM_R(:)=0
  PSIPRM_R(1:nrho_r)=iotabar_t(1:nrho_r)*phiprm_t(1:nrho_r)

ENDIF

!Total toroidal flux
IF(PRESENT(PHI_R)) THEN

  PHI_R(:)=0
  PHI_R(1:nrho_r)=phitot_3d*(rho_r(1:nrho_r)/rhomax_3d)**2

ENDIF

!Total poloidal flux
IF(PRESENT(PSI_R)) THEN

  PSI_R(:)=0
  iflag=0
  message=''
  CALL SPLINE1_INTEG(1,nrho_3d,rho_3d,iotabar_3d,nr,rho_r,values,iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_MAGFLUX/'//message
    IF(iflag > 0) GOTO 9999

  ENDIF

  PSI_R(1:nr)=values(1:nr)*2*phitot_3d/rhomax_3d**2

  IF(nr < nrho_r) THEN

    DO i=nr+1,nrho_r !Over radial nodes

      PSI_R(i)=PSI_R(i-1)+phitot_3d/rhomax_3d**2 &  
               *(rho_r(i-1)*iotabar_t(i-1)+rho_r(i)*iotabar_t(i)) &
               *(rho_r(i)-rho_r(i-1))

    ENDDO !Over radial nodes

  ENDIF

ENDIF

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE AJAX_MAGFLUX

SUBROUTINE AJAX_SHAPE(nrho_r,rho_r, &
                      iflag,message, &  
                      ZETA, &
                      SHIFT_R,ELONG_R,TRIANG_R,RMAX_R,RMIN_R,ZMAX_R,ZMIN_R, &
                      RBOT_R,RTOP_R)
!-------------------------------------------------------------------------------
!AJAX_SHAPE gets shift, elongation, and triangularity
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &     
  nrho_r                 !no. of radial nodes [-]

REAL(KIND=rspec), INTENT(IN) :: &    
  rho_r(:)               !radial nodes [rho]

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error

!Declaration of optional input variables
REAL(KIND=rspec), INTENT(IN), OPTIONAL :: &   
  ZETA                   !toroidal plane for the shape [-]
!                        !=0 by default

!Declaration of optional output variables
REAL(KIND=rspec), INTENT(OUT), OPTIONAL :: &   
  SHIFT_R(:),          & !shift [-]
  ELONG_R(:),          & !elongation [-]
  TRIANG_R(:),         & !triangularity [-]
  RMAX_R(:),           & !maximum R of plasma in zeta plane [m]
  RMIN_R(:),           & !minimum R of plasma in zeta plane [m]
  ZMAX_R(:),           & !maximum Z of plasma in zeta plane [m]
  ZMIN_R(:),           & !minimum Z of plasma in zeta plane [m]
  RBOT_R(:),           & !R at ZMIN_R in zeta plane [m]
  RTOP_R(:)              !R at ZMAX_R in zeta plane [m]

!-------------------------------------------------------------------------------
!Declaration of local variables
LOGICAL :: &      
  l_rminmax

INTEGER :: &      
  i,nr

REAL(KIND=rspec) :: &     
  a0,rg0,r_cyl(1:3),r_flx(1:3),v1(1:nrho_3d),rbot(1:nrho_3d), & 
  rtop(1:nrho_3d),rmax(1:nrho_3d),rmin(1:nrho_3d),zmax(1:nrho_3d), &
  zmin(1:nrho_3d)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Null output
iflag=0
message=''

!Null local values
r_flx(:)=0
r_cyl(:)=0
rbot(:)=0
rtop(:)=0
rmax(:)=0
rmin(:)=0
zmax(:)=0
zmin(:)=0
v1(:)=0

!Set the toroidal angle for the calculations
IF(PRESENT(ZETA)) r_flx(3)=ZETA

!Check whether values outside R,Z domain are requested (nr is last point inside)
LOOP_I: DO i=nrho_r,1,-1 !Over nodes

  nr=i
  IF(rho_r(nr) < rhomax_3d+rhores_3d) EXIT LOOP_I

ENDDO LOOP_I !Over nodes

!-------------------------------------------------------------------------------
!Half diameter and center of plasma boundary
!-------------------------------------------------------------------------------
!Find coordinates where dR/dtheta=0 on inside for R_min
l_rminmax=.TRUE.
r_flx(1)=rhomax_3d
r_flx(2)=z_pi
iflag=0
message=''
CALL AJAX_MINMAX_RZ(l_rminmax,r_flx,r_cyl,iflag,message)

!Check messages
IF(iflag /= 0) THEN

  message='AJAX_SHAPE(1)/'//message
  IF(iflag > 0) GOTO 9999

ENDIF

rmin(nrho_3d)=r_cyl(1)

!Find coordinates where dR/dtheta=0 on outside for R_max
l_rminmax=.TRUE.
r_flx(1)=rhomax_3d
r_flx(2)=0
iflag=0
message=''
CALL AJAX_MINMAX_RZ(l_rminmax,r_flx,r_cyl,iflag,message)

!Check messages
IF(iflag /= 0) THEN

  message='AJAX_SHAPE(2)/'//message
  IF(iflag > 0) GOTO 9999

ENDIF

rmax(nrho_3d)=r_cyl(1)

!Major radius of center of plasma boundary
rg0=(rmin(nrho_3d)+rmax(nrho_3d))/2

!Minor radius = half diameter
a0=(rmax(nrho_3d)-rmin(nrho_3d))/2

!-------------------------------------------------------------------------------
!Find inside, outside (dR/theta=0) and top, bottom (dZ/dtheta=0) of each surface
!-------------------------------------------------------------------------------
DO i=2,nrho_3d !Over radial nodes

  !dR/dtheta=0 on inside
  l_rminmax=.TRUE.
  r_flx(1)=rho_3d(i)
  r_flx(2)=z_pi
  iflag=0
  message=''
  CALL AJAX_MINMAX_RZ(l_rminmax,r_flx,r_cyl,iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_SHAPE(3)/'//message
    IF(iflag > 0) GOTO 9999

  ENDIF

  rmin(i)=r_cyl(1)

  !dR/dtheta=0 on outside
  l_rminmax=.TRUE.
  r_flx(1)=rho_3d(i)
  r_flx(2)=0
  iflag=0
  message=''
  CALL AJAX_MINMAX_RZ(l_rminmax,r_flx,r_cyl,iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_SHAPE(4)/'//message
    IF(iflag > 0) GOTO 9999

  ENDIF

  rmax(i)=r_cyl(1)

  !dZ/dtheta=0 on bottom
  l_rminmax=.FALSE.
  r_flx(1)=rho_3d(i)
  r_flx(2)=z_pi/2
  iflag=0
  message=''
  CALL AJAX_MINMAX_RZ(l_rminmax,r_flx,r_cyl,iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_SHAPE(5)/'//message
    IF(iflag > 0) GOTO 9999

  ENDIF

  rbot(i)=r_cyl(1)
  zmin(i)=r_cyl(3)

  !dZ/dtheta=0 on top
  l_rminmax=.FALSE.
  r_flx(1)=rho_3d(i)
  r_flx(2)=-z_pi/2
  iflag=0
  message=''
  CALL AJAX_MINMAX_RZ(l_rminmax,r_flx,r_cyl,iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_SHAPE(6)/'//message
    IF(iflag > 0) GOTO 9999

  ENDIF

  rtop(i)=r_cyl(1)
  zmax(i)=r_cyl(3)

ENDDO !Over radial nodes

!-------------------------------------------------------------------------------
!Shift
!-------------------------------------------------------------------------------
IF(PRESENT(SHIFT_R)) THEN

  !Initialization
  SHIFT_R(1:nrho_r)=0

  !Shift relative to center of outer surface   
  v1(2:nrho_3d)=((rmax(2:nrho_3d)+rmin(2:nrho_3d))/2-rg0)/a0

  !Set axial value
  v1(1)=v1(2)-rho_3d(2)*(v1(3)-v1(2))/(rho_3d(3)-rho_3d(2))

  !Interpolate to external grid points inside R,Z domain
  iflag=0
  message=''
  CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,SHIFT_R,iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_SHAPE(7)/'//message
    GOTO 9999

  ENDIF

  IF(nr < nrho_r) THEN

    !Extrapolate to external grid points outside R,Z domain
    SHIFT_R(nr+1:nrho_r)=SHIFT_R(nr)

  ENDIF

ENDIF

!-------------------------------------------------------------------------------
!Elongation
!-------------------------------------------------------------------------------
IF(PRESENT(ELONG_R)) THEN

  !Initialization
  ELONG_R(1:nrho_r)=0

  !Elongation is total height over total width
  v1(2:nrho_3d)=(zmax(2:nrho_3d)-zmin(2:nrho_3d)) &  
                /(rmax(2:nrho_3d)-rmin(2:nrho_3d))

  !Set axial value
  v1(1)=v1(2)-rho_3d(2)*(v1(3)-v1(2))/(rho_3d(3)-rho_3d(2))

  !Interpolate to external grid points inside R,Z domain
  iflag=0
  message=''
  CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,ELONG_R,iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_SHAPE(8)/'//message
    GOTO 9999

  ENDIF

  IF(nr < nrho_r) THEN

    !Extrapolate to external grid points outside R,Z domain
    ELONG_R(nr+1:nrho_r)=ELONG_R(nr)

  ENDIF

ENDIF

!-------------------------------------------------------------------------------
!Triangularity
!-------------------------------------------------------------------------------
IF(PRESENT(TRIANG_R)) THEN

  !Initialization
  TRIANG_R(1:nrho_r)=0

  !Triangularity - average of upper and lower
  v1(2:nrho_3d)=((rmax(2:nrho_3d)+rmin(2:nrho_3d)) &  
                -(rbot(2:nrho_3d)+rtop(2:nrho_3d))) &  
                /(rmax(2:nrho_3d)-rmin(2:nrho_3d))

  !Set axial value
  v1(1)=0

  !Interpolate to external grid points inside R,Z domain
  iflag=0
  message=''
  CALL LINEAR1_INTERP(nrho_3d,rho_3d,v1,nr,rho_r,TRIANG_R,iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_SHAPE(9)/'//message
    GOTO 9999

  ENDIF

  IF(nr < nrho_r) THEN

    !Extrapolate to external grid points outside R,Z domain
    TRIANG_R(nr+1:nrho_r)=TRIANG_R(nr)

  ENDIF

ENDIF

!-------------------------------------------------------------------------------
!Maximum R
!-------------------------------------------------------------------------------
IF(PRESENT(RMAX_R)) THEN

  !Initialization
  RMAX_R(1:nrho_r)=0

  !To cover odd-shaped plasmas, look for maxima near pi/4 and 7pi/4
  DO i=2,nrho_3d !Over radial nodes

    l_rminmax=.TRUE.
    r_flx(1)=rho_3d(i)
    r_flx(2)=z_pi/4
    iflag=0
    message=''
    CALL AJAX_MINMAX_RZ(l_rminmax,r_flx,r_cyl,iflag,message)

    !Check messages
    IF(iflag /= 0) THEN

      message='AJAX_SHAPE(10)/'//message
      GOTO 9999

    ENDIF

    IF(r_cyl(1) > rmax(i)) rmax(i)=r_cyl(1)
    r_flx(2)=7*z_pi/4
    iflag=0
    message=''
    CALL AJAX_MINMAX_RZ(l_rminmax,r_flx,r_cyl,iflag,message)

    !Check messages
    IF(iflag /= 0) THEN

      message='AJAX_SHAPE(11)/'//message
      GOTO 9999

    ENDIF

    IF(r_cyl(1) > rmax(i)) rmax(i)=r_cyl(1)

  ENDDO !Over radial nodes

  !Set axial value
  rmax(1)=rmax(2)-rho_3d(2)*(rmax(3)-rmax(2))/(rho_3d(3)-rho_3d(2))

  !Interpolate to external grid points inside R,Z domain
  iflag=0
  message=''
  CALL LINEAR1_INTERP(nrho_3d,rho_3d,rmax,nr,rho_r,RMAX_R,iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_SHAPE(12)/'//message
    GOTO 9999

  ENDIF

ENDIF

!-------------------------------------------------------------------------------
!Minimum R
!-------------------------------------------------------------------------------
IF(PRESENT(RMIN_R)) THEN

  !Initialization
  RMIN_R(1:nrho_r)=0

  !To cover odd-shaped plasmas, look for minima near 3pi/4 and 5pi/4
  DO i=2,nrho_3d !Over radial nodes

    l_rminmax=.TRUE.
    r_flx(1)=rho_3d(i)
    r_flx(2)=3*z_pi/4
    iflag=0
    message=''
    CALL AJAX_MINMAX_RZ(l_rminmax,r_flx,r_cyl,iflag,message)

    !Check messages
    IF(iflag /= 0) THEN

      message='AJAX_SHAPE(13)/'//message
      GOTO 9999

    ENDIF

    IF(r_cyl(1) < rmin(i)) rmin(i)=r_cyl(1)
    r_flx(2)=5*z_pi/4
    iflag=0
    message=''
    CALL AJAX_MINMAX_RZ(l_rminmax,r_flx,r_cyl,iflag,message)

    !Check messages
    IF(iflag /= 0) THEN

      message='AJAX_SHAPE(14)/'//message
      GOTO 9999

    ENDIF

    IF(r_cyl(1) < rmin(i)) rmin(i)=r_cyl(1)

  ENDDO !Over radial nodes

  !Set axial value
  rmin(1)=rmin(2)-rho_3d(2)*(rmin(3)-rmin(2))/(rho_3d(3)-rho_3d(2))

  !Interpolate to external grid points inside R,Z domain
  iflag=0
  message=''
  CALL LINEAR1_INTERP(nrho_3d,rho_3d,rmin,nr,rho_r,RMIN_R,iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_SHAPE(15)/'//message
    GOTO 9999

  ENDIF

ENDIF

!-------------------------------------------------------------------------------
!Maximum Z
!-------------------------------------------------------------------------------
IF(PRESENT(ZMAX_R)) THEN

  !Initialization
  ZMAX_R(1:nrho_r)=0

  !Set axial value
  zmax(1)=zmax(2)-rho_3d(2)*(zmax(3)-zmax(2))/(rho_3d(3)-rho_3d(2))

  !Interpolate to external grid points inside R,Z domain
  iflag=0
  message=''
  CALL LINEAR1_INTERP(nrho_3d,rho_3d,zmax,nr,rho_r,ZMAX_R,iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_SHAPE(16)/'//message
    GOTO 9999

  ENDIF

ENDIF

!-------------------------------------------------------------------------------
!Minimum Z
!-------------------------------------------------------------------------------
IF(PRESENT(ZMIN_R)) THEN

  !Initialization
  ZMIN_R(1:nrho_r)=0

  !Set axial value
  zmin(1)=zmin(2)-rho_3d(2)*(zmin(3)-zmin(2))/(rho_3d(3)-rho_3d(2))

  !Interpolate to external grid points inside R,Z domain
  iflag=0
  message=''
  CALL LINEAR1_INTERP(nrho_3d,rho_3d,zmin,nr,rho_r,ZMIN_R,iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_SHAPE(17)/'//message
    GOTO 9999

  ENDIF

ENDIF

!-------------------------------------------------------------------------------
!R at maximum Z
!-------------------------------------------------------------------------------
IF(PRESENT(RTOP_R)) THEN

  !Initialization
  RTOP_R(1:nrho_r)=0

  !Set axial value
  rtop(1)=rtop(2)-rho_3d(2)*(rtop(3)-rtop(2))/(rho_3d(3)-rho_3d(2))

  !Interpolate to external grid points inside R,Z domain
  iflag=0
  message=''
  CALL LINEAR1_INTERP(nrho_3d,rho_3d,rtop,nr,rho_r,RTOP_R,iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_SHAPE(18)/'//message
    GOTO 9999

  ENDIF

ENDIF

!-------------------------------------------------------------------------------
!R at minimum Z
!-------------------------------------------------------------------------------
IF(PRESENT(RBOT_R)) THEN

  !Initialization
  RBOT_R(1:nrho_r)=0

  !Set axial value
  rbot(1)=rbot(2)-rho_3d(2)*(rbot(3)-rbot(2))/(rho_3d(3)-rho_3d(2))

  !Interpolate to external grid points inside R,Z domain
  iflag=0
  message=''
  CALL LINEAR1_INTERP(nrho_3d,rho_3d,rbot,nr,rho_r,RBOT_R,iflag,message)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_SHAPE(19)/'//message
    GOTO 9999

  ENDIF

ENDIF

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE AJAX_SHAPE

SUBROUTINE AJAX_GLOBALS(iflag,message, &   
                        R000,PHITOT,RHOMAX,RHOMIN,SIGNBP,SIGNBT,NPER,NRHO, &
                        NTHETA,NZETA)
!-------------------------------------------------------------------------------
!AJAX_GLOBALS gets global characteristics of the data
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------
      
!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error

!Declaration of optional output variables
INTEGER, INTENT(OUT), OPTIONAL :: &    
  NPER,                & !number of field periods, =0 for axisymmetry [-]
  NRHO,                & !number of radial nodes in internal data [-]
  NTHETA,              & !number of poloidal nodes in internal data [-]
  NZETA                  !number of toroidal nodes in internal data [-]

REAL(KIND=rspec), INTENT(OUT), OPTIONAL :: &   
  R000,                & !coefficient of (m,n)=(0,0) mode for R at axis [m]
  PHITOT,              & !total toroidal flux [Wb]
  RHOMIN,              & !inner radial boundary of R,Z domain [rho]
  RHOMAX,              & !outer radial boundary of R,Z domain [rho]
  SIGNBP,              & !sign of the poloidal magnetic field [-]
                         !positive is down on outside
  SIGNBT                 !sign of the toroidal magnetic field [-]
                         !positive is counterclockwise from top view

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Null output
iflag=0
message=''

!-------------------------------------------------------------------------------
!Return values
!-------------------------------------------------------------------------------
IF(PRESENT(R000)) R000=r000_3d
IF(PRESENT(PHITOT)) PHITOT=phitot_3d
IF(PRESENT(RHOMAX)) RHOMAX=rhomax_3d
IF(PRESENT(RHOMIN)) RHOMIN=rhomin_3d
IF(PRESENT(SIGNBP)) SIGNBP=signbp_3d
IF(PRESENT(SIGNBT)) SIGNBT=signbt_3d
IF(PRESENT(NPER)) NPER=nper_3d
IF(PRESENT(NRHO)) NRHO=nrho_3d
IF(PRESENT(NTHETA)) NTHETA=ntheta_3d
IF(PRESENT(NZETA)) NZETA=nzeta_3d

END SUBROUTINE AJAX_GLOBALS

SUBROUTINE AJAX_MINMAX_RZ(l_rminmax, &
                          r_flx, &
                          r_cyl,iflag,message)
!-------------------------------------------------------------------------------
!AJAX_MINMAX_RZ gets maximum or minimum of R or Z
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
LOGICAL, INTENT(IN) :: &     
  l_rminmax              !switch to select search for R or Z [logical]
                         !=.TRUE. find maximum or minimum R
                         !=.FALSE. find maximum or minimum Z

!Declaration of input/output variables
REAL(KIND=rspec), INTENT(INOUT) :: &    
  r_flx(:)               !flux coordinates (rho,theta,zeta) [rho,rad,rad]

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error

REAL(KIND=rspec), INTENT(OUT) :: &    
  r_cyl(:)               !cylindrical coordinates (R,phi,Z) [m,rad,m]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &      
  i,ig

REAL(KIND=rspec) :: &     
  g_cyl0(1:6),g_cyl1(1:6),r_flx0(1:3),r_flx1(1:3)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Null output
iflag=0
message=''

!Null local values
r_flx0(:)=0
r_flx1(:)=0
g_cyl0(:)=0
g_cyl1(:)=0

!Index for derivative
IF(l_rminmax) THEN

  !d(R)/d(theta)
  ig=2

ELSE

  !d(Z)/d(theta)
  ig=5

ENDIF

!-------------------------------------------------------------------------------
!Find extrema
!-------------------------------------------------------------------------------
r_flx0(1:3)=r_flx(1:3)
r_flx0(2)=r_flx0(2)-0.05*z_pi
iflag=0
message=''
CALL AJAX_FLX2CYL(r_flx0,r_cyl,iflag,message, &
                  G_CYL=g_cyl0)

!Check messages
IF(iflag /= 0) THEN

  message='AJAX_MINMAX_RZ(1)/'//message
  IF(iflag > 0) GOTO 9999

ENDIF

r_flx1(1:3)=r_flx(1:3)
r_flx1(2)=r_flx1(2)+0.05*z_pi
iflag=0
message=''
CALL AJAX_FLX2CYL(r_flx1,r_cyl,iflag,message, &   
                  G_CYL=g_cyl1)

!Check messages
IF(iflag /= 0) THEN

  message='AJAX_MINMAX_RZ(2)/'//message
  IF(iflag > 0) GOTO 9999

ENDIF

!Normally 1 or 2 iterations are required
LOOP_I: DO i=1,10 !Over iteration

  r_flx(2)=r_flx0(2)-g_cyl0(ig)*(r_flx1(2)-r_flx0(2))/(g_cyl1(ig)-g_cyl0(ig))
  r_flx0(1:3)=r_flx1(1:3)
  g_cyl0(:)=g_cyl1(:)
  r_flx1(1:3)=r_flx(1:3)
  iflag=0
  message=''
  CALL AJAX_FLX2CYL(r_flx,r_cyl,iflag,message, &  
                    G_CYL=g_cyl1)

  !Check messages
  IF(iflag /= 0) THEN

    message='AJAX_MINMAX_RZ(3)/'//message
    IF(iflag > 0) GOTO 9999

  ENDIF

  IF(ABS(g_cyl1(ig)) < 1.0e-5*r000_3d/z_pi) EXIT LOOP_I

ENDDO LOOP_I !Over iteration

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE AJAX_MINMAX_RZ

SUBROUTINE AJAX_INIT
!-------------------------------------------------------------------------------
!AJAX_INIT initializes or resets private data for the radial, poloidal and
!  toroidal grids, and allocates storage for the R, Z, and stream function
!  expansions
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of local variables
INTEGER :: &      
  i,nset1,nset2,kmax

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Machine constants
z_large=HUGE(1.0_rspec)
z_precision=EPSILON(1.0_rspec)
z_small=TINY(1.0_rspec)

!Resolution of rho
rhores_3d=1.0e-3*rhomax_3d

!-------------------------------------------------------------------------------
!Radial grid and R,Z
!-------------------------------------------------------------------------------
IF(ALLOCATED(rho_3d)) THEN

  !Check dimensions
  nset1=SIZE(r_3d,2)
  nset2=SIZE(r_3d,3)

  !If storage requirements have changed, reallocate
  IF(nset1 /= nrho_3d .OR. &     
     nset2 /= krz_3d) THEN

    !Reallocate radial grid and R,Z
    DEALLOCATE(rho_3d, &     
               r_3d, &     
               z_3d)
    ALLOCATE(rho_3d(nrho_3d), &    
             r_3d(4,nrho_3d,krz_3d), &   
             z_3d(4,nrho_3d,krz_3d))

      rho_3d(:)=0
      r_3d(:,:,:)=0
      z_3d(:,:,:)=0
      
    !Reallocate radial grid and Rs,Zc 09/09/11 SAL
    DEALLOCATE(rs_3d, &     
               zc_3d)
    ALLOCATE(rs_3d(4,nrho_3d,krz_3d), &   
             zc_3d(4,nrho_3d,krz_3d))

      rs_3d(:,:,:)=0
      zc_3d(:,:,:)=0

  ENDIF

ELSE

  !Allocate radial grid and R,Z
  ALLOCATE(rho_3d(nrho_3d), &    
           r_3d(4,nrho_3d,krz_3d), &   
           z_3d(4,nrho_3d,krz_3d))

    rho_3d(:)=0
    r_3d(:,:,:)=0
    z_3d(:,:,:)=0
                                                                        
  !Allocate radial grid and Rs,Zc  09/09/11 SAL                                     
  ALLOCATE(rs_3d(4,nrho_3d,krz_3d), &   
           zc_3d(4,nrho_3d,krz_3d))

    rs_3d(:,:,:)=0
    zc_3d(:,:,:)=0

ENDIF

!Set radial grid
rho_3d(:)=REAL((/ (i-1,i=1,nrho_3d) /),rspec)*rhomax_3d/(nrho_3d-1)

!-------------------------------------------------------------------------------
!Lambda
!-------------------------------------------------------------------------------
IF(ALLOCATED(lam_3d)) THEN

  !Check dimensions
  nset1=SIZE(lam_3d,2)
  nset2=SIZE(lam_3d,3)

  !If storage requirements have changed, reallocate
  IF(nset1 /= nrho_3d .OR. &     
     nset2 /= klam_3d) THEN

    !Realloate lambda
    DEALLOCATE(lam_3d)
    ALLOCATE(lam_3d(4,nrho_3d,klam_3d))

      lam_3d(:,:,:)=0

    !Realloate lambda 09/09/11 SAL
    DEALLOCATE(lamc_3d)
    ALLOCATE(lamc_3d(4,nrho_3d,klam_3d))

      lamc_3d(:,:,:)=0

  ENDIF

ELSE

  !Allocate lambda
  ALLOCATE(lam_3d(4,nrho_3d,klam_3d))

    lam_3d(:,:,:)=0
  
  !Allocate lambda 09/09/11 SAL
  ALLOCATE(lamc_3d(4,nrho_3d,klam_3d))

    lamc_3d(:,:,:)=0

ENDIF

!-------------------------------------------------------------------------------
!Toroidal and poloidal modes
!-------------------------------------------------------------------------------
kmax=MAX(krz_3d,klam_3d)

IF(ALLOCATED(m_3d)) THEN

  !Check dimension
  nset1=SIZE(m_3d)

  !If storage requirements have changed, reallocate
  IF(nset1 /= kmax) THEN

    !Reallocate toroidal and poloidal modes
    DEALLOCATE(m_3d, &     
               n_3d, &     
               mabs_3d)
    ALLOCATE(m_3d(kmax), &     
             n_3d(kmax), &     
             mabs_3d(kmax))

      m_3d(:)=0
      n_3d(:)=0
      mabs_3d(:)=0

  ENDIF

ELSE

  !Allocate toroidal and poloidal modes
  ALLOCATE(m_3d(kmax), &     
           n_3d(kmax), &     
           mabs_3d(kmax))

    m_3d(:)=0
    n_3d(:)=0
    mabs_3d(:)=0

ENDIF

!-------------------------------------------------------------------------------
!Poloidal grid spacing and weighting
!-------------------------------------------------------------------------------
dtheta_3d=2*z_pi/(ntheta_3d-1)

IF(ALLOCATED(wtheta_3d)) THEN

  !Check dimension
  nset1=SIZE(wtheta_3d)

  IF(nset1 /= ntheta_3d) THEN

    !Reallocate poloidal weighting
    DEALLOCATE(wtheta_3d)
    ALLOCATE(wtheta_3d(ntheta_3d))

      wtheta_3d(:)=0    

    !Set poloidal weights
    wtheta_3d(1)=1
    wtheta_3d(ntheta_3d)=1

    DO i=2,ntheta_3d-1 !Over poloidal nodes

      IF(MOD(i,2) == 0) THEN

        !Even node
        wtheta_3d(i)=4

      ELSE

        !Odd mode
        wtheta_3d(i)=2

      ENDIF

    ENDDO !Over poloidal nodes

    wtheta_3d(:)=wtheta_3d(:)*2*z_pi/3/(ntheta_3d-1)

  ENDIF

ELSE

  !Allocate poloidal weighting
  ALLOCATE(wtheta_3d(ntheta_3d))

    wtheta_3d(:)=0

  !Set poloidal weights
  wtheta_3d(1)=1
  wtheta_3d(ntheta_3d)=1

  DO i=2,ntheta_3d-1 !Over poloidal nodes

    IF(MOD(i,2) == 0) THEN

      !Even node
      wtheta_3d(i)=4

    ELSE

      !Odd mode
      wtheta_3d(i)=2

    ENDIF

  ENDDO !Over poloidal nodes

  wtheta_3d(:)=wtheta_3d(:)*2*z_pi/3/(ntheta_3d-1)

ENDIF

!-------------------------------------------------------------------------------
!Poloidal grid spacing and weighting
!-------------------------------------------------------------------------------
IF(nzeta_3d == 1) THEN

  dzeta_3d=2*z_pi

ELSE

  dzeta_3d=2*z_pi/(nzeta_3d-1)/nper_3d

ENDIF

!Check allocation of toroidal weighting
IF(ALLOCATED(wzeta_3d)) THEN

  !Check dimension
  nset1=SIZE(wzeta_3d)

  IF(nset1 /= nzeta_3d) THEN

    !Reallocate toroidal weighting
    DEALLOCATE(wzeta_3d)
    ALLOCATE(wzeta_3d(nzeta_3d))

      wzeta_3d(:)=0

    !Set toroidal weights
    wzeta_3d(1)=dzeta_3d

    IF(nzeta_3d > 1) THEN

      !Non-axisymmetric plasma
      wzeta_3d(1)=1
      wzeta_3d(nzeta_3d)=1

      DO i=2,nzeta_3d-1 !Over toroidal nodes

        IF(MOD(i,2) == 0) THEN

          !Even node
          wzeta_3d(i)=4

        ELSE

          !Odd mode
          wzeta_3d(i)=2

        ENDIF

      ENDDO !Over toroidal nodes

      wzeta_3d(:)=wzeta_3d(:)*nper_3d*dzeta_3d/3

    ENDIF

  ENDIF

ELSE

  !Allocate toroidal weighting
  ALLOCATE(wzeta_3d(nzeta_3d))

    wzeta_3d(:)=0

  !Set toroidal weights
  wzeta_3d(1)=dzeta_3d

  IF(nzeta_3d > 1) THEN

    !Non-axisymmetric plasma
    wzeta_3d(1)=1
    wzeta_3d(nzeta_3d)=1

    DO i=2,nzeta_3d-1 !Over toroidal nodes

      IF(MOD(i,2) == 0) THEN

        !Even node
        wzeta_3d(i)=4

      ELSE

        !Odd mode
        wzeta_3d(i)=2

      ENDIF

    ENDDO !Over toroidal nodes

    wzeta_3d(:)=wzeta_3d(:)*nper_3d*dzeta_3d/3

  ENDIF

ENDIF

!-------------------------------------------------------------------------------
!Set logical switches for flux surface averaging
!-------------------------------------------------------------------------------
l_fluxavb_3d=.FALSE.
l_fluxavg_3d=.FALSE.

END SUBROUTINE AJAX_INIT

SUBROUTINE AJAX_INIT_FLUXAV_B
!-------------------------------------------------------------------------------
!AJAX_INIT_FLUXAV_B initializes magnetic field arrays
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of local variables
INTEGER :: &      
  i,j,k,nset1,nset2,nset3

REAL(KIND=rspec) :: &     
  br,bz,bzeta,phiprm,r_flx(1:3)

!-------------------------------------------------------------------------------
!Load lambda derivatives for flux surface averaging
!-------------------------------------------------------------------------------
!Check allocation
IF(ALLOCATED(eltheta_3d)) THEN

  !Check dimensions
  nset1=SIZE(eltheta_3d,1)
  nset2=SIZE(eltheta_3d,2)
  nset3=SIZE(eltheta_3d,3)

  !If storage requirements have changed, reallocate
  IF(nset1 /= nrho_3d .OR. &    
     nset2 /= ntheta_3d .OR. &    
     nset3 /= nzeta_3d) THEN

    !Reallocate lambda derivatives
    DEALLOCATE(eltheta_3d, &    
               elzeta_3d)
    ALLOCATE(eltheta_3d(nrho_3d,ntheta_3d,nzeta_3d), & 
             elzeta_3d(nrho_3d,ntheta_3d,nzeta_3d))

  ENDIF

ELSE

  !Allocate lambda derivatives
  ALLOCATE(eltheta_3d(nrho_3d,ntheta_3d,nzeta_3d), & 
           elzeta_3d(nrho_3d,ntheta_3d,nzeta_3d))

ENDIF

!Initialization
r_flx(:)=0
eltheta_3d(:,:,:)=0
elzeta_3d(:,:,:)=0

!Fill arrays
DO i=2,nrho_3d !Over radial nodes

  r_flx(1)=rho_3d(i)

  DO j=1,ntheta_3d !Over poloidal nodes

    r_flx(2)=(j-1)*dtheta_3d

    DO k=1,nzeta_3d !Over toroidal nodes

      r_flx(3)=(k-1)*dzeta_3d

      CALL AJAX_LAMBDA(r_flx, &
                       eltheta_3d(i,j,k),elzeta_3d(i,j,k))

    ENDDO !Over toroidal nodes

  ENDDO !Over poloidal nodes

ENDDO !Over radial nodes

!-------------------------------------------------------------------------------
!Load B data for flux surface averaging
!-------------------------------------------------------------------------------
!Check allocation      
IF(ALLOCATED(b_3d)) THEN

  !Check dimensions
  nset1=SIZE(b_3d,1)
  nset2=SIZE(b_3d,2)
  nset3=SIZE(b_3d,3)

  !If storage requirements have changed, reallocate
  IF(nset1 /= nrho_3d .OR. &     
     nset2 /= ntheta_3d .OR. &    
     nset3 /= nzeta_3d) THEN

    !Reallocate B data
    DEALLOCATE(b_3d, &     
               btheta_3d)
    ALLOCATE(b_3d(nrho_3d,ntheta_3d,nzeta_3d), &  
             btheta_3d(nrho_3d,ntheta_3d,nzeta_3d))
  ENDIF

ELSE

  !Allocate B data
  ALLOCATE(b_3d(nrho_3d,ntheta_3d,nzeta_3d), &  
           btheta_3d(nrho_3d,ntheta_3d,nzeta_3d))

ENDIF

!Initialization
btheta_3d(:,:,:)=0
b_3d(:,:,:)=0

!Evaluate 3D arrays
DO i=2,nrho_3d !Over radial nodes

  phiprm=2*rho_3d(i)*phitot_3d/rhomax_3d**2

  DO j=1,ntheta_3d !Over poloidal nodes

    DO k=1,nzeta_3d !Over toroidal nodes

      btheta_3d(i,j,k)=(iotabar_3d(1,i)-elzeta_3d(i,j,k))
      bzeta=(1.0+eltheta_3d(i,j,k))
      br=gcyl_3d(2,i,j,k)*btheta_3d(i,j,k)+gcyl_3d(3,i,j,k)*bzeta
      bz=gcyl_3d(5,i,j,k)*btheta_3d(i,j,k)+gcyl_3d(6,i,j,k)*bzeta
      b_3d(i,j,k)=SQRT(br**2+bz**2+(rcyl_3d(i,j,k)*bzeta)**2)/gsqrt_3d(i,j,k)
      btheta_3d(i,j,k)=btheta_3d(i,j,k)/gsqrt_3d(i,j,k)

    ENDDO !Over toroidal nodes

  ENDDO !Over poloidal nodes

  btheta_3d(i,:,:)=phiprm*btheta_3d(i,:,:)
  b_3d(i,:,:)=ABS(phiprm)*b_3d(i,:,:)

ENDDO !Over radial nodes

!Poloidal values at origin are linear extrapolation of averages at 2 and 3
DO k=1,nzeta_3d

  btheta_3d(1,1:ntheta_3d,k)=(SUM(btheta_3d(2,1:ntheta_3d-1,k))*rho_3d(3) &   
                             -SUM(btheta_3d(3,1:ntheta_3d-1,k)) & 
                              *rho_3d(2))/(ntheta_3d-1)/(rho_3d(3)-rho_3d(2))
  b_3d(1,1:ntheta_3d,k)=(SUM(b_3d(2,1:ntheta_3d-1,k))*rho_3d(3) &    
                        -SUM(b_3d(3,1:ntheta_3d-1,k)) &  
                         *rho_3d(2))/(ntheta_3d-1)/(rho_3d(3)-rho_3d(2))

ENDDO

!Divide by 2*pi
btheta_3d(:,:,:)=btheta_3d(:,:,:)/2/z_pi
b_3d(:,:,:)=b_3d(:,:,:)/2/z_pi

!-------------------------------------------------------------------------------
!Set initialization flag
!-------------------------------------------------------------------------------
l_fluxavb_3d=.TRUE.

END SUBROUTINE AJAX_INIT_FLUXAV_B

SUBROUTINE AJAX_INIT_FLUXAV_G(iflag,message)
!-------------------------------------------------------------------------------
!AJAX_INIT_FLUXAV_G initializes geometry arrays
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &      
  i,j,k,nset1,nset2,nset3

REAL(KIND=rspec) :: &     
  r_cyl(1:3),r_flx(1:3)

!-------------------------------------------------------------------------------
!Null output
!-------------------------------------------------------------------------------
!Error flag and message
iflag=0
message=''

!-------------------------------------------------------------------------------
!Load R,Z and other geometric data for flux surface averaging
!-------------------------------------------------------------------------------
!Check allocation
IF(ALLOCATED(rcyl_3d)) THEN

  !Check dimensions
        nset1=SIZE(rcyl_3d,1)
        nset2=SIZE(rcyl_3d,2)
        nset3=SIZE(rcyl_3d,3)

  !If storage requirements have changed, reallocate
  IF(nset1 /= nrho_3d .OR. &     
     nset2 /= ntheta_3d .OR. &    
     nset3 /= nzeta_3d) THEN

    !Reallocate R,Z and other geometric data
    DEALLOCATE(rcyl_3d, &     
               zcyl_3d, &     
               gsqrt_3d, &     
               gcyl_3d, &     
               vp_3d, &     
               gt_3d, &     
               az_3d)
    ALLOCATE(rcyl_3d(nrho_3d,ntheta_3d,nzeta_3d), &
             zcyl_3d(nrho_3d,ntheta_3d,nzeta_3d), &
             gsqrt_3d(nrho_3d,ntheta_3d,nzeta_3d), &
             gcyl_3d(6,nrho_3d,ntheta_3d,nzeta_3d), &
             vp_3d(nrho_3d), &
             gt_3d(ntheta_3d), &
             az_3d(nzeta_3d))

  ENDIF

ELSE

  !Allocate R,Z and other geometric data
  ALLOCATE(rcyl_3d(nrho_3d,ntheta_3d,nzeta_3d), &
           zcyl_3d(nrho_3d,ntheta_3d,nzeta_3d), &
           gsqrt_3d(nrho_3d,ntheta_3d,nzeta_3d), &
           gcyl_3d(6,nrho_3d,ntheta_3d,nzeta_3d), &
           vp_3d(nrho_3d), &
           gt_3d(ntheta_3d), &
           az_3d(nzeta_3d))

ENDIF

!Initialization
rcyl_3d(:,:,:)=0
zcyl_3d(:,:,:)=0
gsqrt_3d(:,:,:)=0
gcyl_3d(:,:,:,:)=0
vp_3d(:)=0
gt_3d(:)=0
az_3d(:)=0
r_flx(:)=0
r_cyl(:)=0

!Evaluate 3D cylindrical coordinates, metrics, and Jacobian
DO i=2,nrho_3d !Over radial nodes

  DO j=1,ntheta_3d !Over poloidal nodes

    DO k=1,nzeta_3d !Over toroidal nodes

      r_flx(1)=rho_3d(i)
      r_flx(2)=(j-1)*dtheta_3d
      r_flx(3)=(k-1)*dzeta_3d
      iflag=0
      message=''
      CALL AJAX_FLX2CYL(r_flx,r_cyl,iflag,message, &  
                        G_CYL=gcyl_3d(1:6,i,j,k), &  
                        GSQRT=gsqrt_3d(i,j,k))

      !Check messages
      IF(iflag /= 0) THEN

        message='AJAX_INIT_FLUXAV_G/'//message
        IF(iflag > 0) GOTO 9999

      ENDIF

      rcyl_3d(i,j,k)=r_cyl(1)
      zcyl_3d(i,j,k)=r_cyl(3)

    ENDDO !Over toroidal nodes

  ENDDO !Over poloidal nodes

ENDDO !Over radial nodes

!Calculate V' on internal grid 
DO i=2,nrho_3d !Over radial nodes

  DO k=1,nzeta_3d !Over toroidal nodes

    gt_3d(1:ntheta_3d)=gsqrt_3d(i,1:ntheta_3d,k)
    az_3d(k)=SUM(wtheta_3d*gt_3d) !Over poloidal nodes

  ENDDO !Over toroidal nodes

  vp_3d(i)=SUM(wzeta_3d*az_3d) !Over toroidal nodes

ENDDO !Over radial nodes

!-------------------------------------------------------------------------------
!Set logical switches for flux surface averaging
!-------------------------------------------------------------------------------
l_fluxavg_3d=.TRUE.

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE AJAX_INIT_FLUXAV_G

SUBROUTINE AJAX_INIT_LAMBDA
!-------------------------------------------------------------------------------
!AJAX_INIT_LAMBDA calculates the magnetic stream function expansion coefficients
!  for an axisymmetric plasma
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of local variables
INTEGER :: &      
  i,j,k,k_bc1=3,k_bcn=0

REAL(KIND=rspec) :: &     
  theta,g1(1:ntheta_3d),g2(1:ntheta_3d)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
lam_3d(:,:,:)=0
g1(:)=0
g2(:)=0

!-------------------------------------------------------------------------------
!Calculate lambda values
!-------------------------------------------------------------------------------
!Loop over harmonics
DO k=1,klam_3d !Over modes

  IF(k /= km0n0_3d) THEN

    !Loop over radial grid
    DO i=2,nrho_3d !Over radial nodes

      !Loop over theta grid
      DO j=1,ntheta_3d !Over poloidal nodes

              g1(j)=gsqrt_3d(i,j,1)/rcyl_3d(i,j,1)**2
              theta=(j-1)*dtheta_3d
              g2(j)=g1(j)*COS(m_3d(k)*theta)

      ENDDO !Over poloidal nodes

      !Theta averages to find lambda coefficients
      lam_3d(1,i,k)=2*SUM(wtheta_3d*g2)/SUM(wtheta_3d*g1)/m_3d(k) !Over poloidal nodes

    ENDDO !Over radial nodes

    lam_3d(1,2:nrho_3d,k)=lam_3d(1,2:nrho_3d,k)/rho_3d(2:nrho_3d)**mabs_3d(k)

    !Make a parabolic fit to axis
    lam_3d(1,1,k)=(lam_3d(1,2,k)*rho_3d(3)**2-lam_3d(1,3,k)*rho_3d(2)**2) &  
                  /(rho_3d(3)**2-rho_3d(2)**2)
 
    !Spline the Lmn coefficients for internal storage (not-a-knot edge BC)
    CALL SPLINE1_FIT(nrho_3d,rho_3d,lam_3d(:,:,k), &  
                     K_BC1=k_bc1, &    
                     K_BCN=k_bcn)

  ENDIF

ENDDO !Over modes

END SUBROUTINE AJAX_INIT_LAMBDA

SUBROUTINE AJAX_LOAD_LAMBDA(k_grid,nr_lam,nk_lam,rhorz,rho_lam,lam, & 
                            iflag,message,LAM_C2)
!-------------------------------------------------------------------------------
!AJAX_LOAD_LAMBDA loads the magnetic stream function
!
!References:
!  W.A.Houlberg, F90 free format 8/2004
!
!Comments:
!  The radial grid is assumed to be the same as for the R,Z data, so the scale
!    factor, rhorz, and grid type, k_grid are needed input
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &     
  k_grid,              & !option designating type of input radial grid [-]
                         !=0 proportional to sqrt(toroidal flux)
                         !=1 proportional to toroidal flux
                         !=else not allowed
  nr_lam,              & !no. of radial nodes in the input lambda [-]
  nk_lam                 !no. of poloidal & toroidal modes in the input lambda [-]

REAL(KIND=rspec), INTENT(IN) :: &    
  rhorz,               & !max rho of the R,Z data [arb]
  rho_lam(:),          & !radial nodes in the input lambda [arb]
  lam(:,:)               !expansion coeffs for lambda [-] (sin)
  
  
REAL(KIND=rspec), INTENT(IN), OPTIONAL :: &  ! SAL 09/09/11
  LAM_C2(:,:)             !expansion coeffs for lambda [-] (cos)

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag                  !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error

!-------------------------------------------------------------------------------
!Declaration of local variables
LOGICAL :: &      
  l_edge

INTEGER :: &      
  j,k,nset1,nset2,k_vopt(1:3)=(/1,0,0/),k_bc1=3,k_bcn=0

REAL(KIND=rspec), ALLOCATABLE :: &    
  rho(:),lmn(:,:),fspl(:,:),values(:,:)
  
REAL(KIND=rspec), ALLOCATABLE :: &    ! 09/09/11 SAL
  lmnc(:,:)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Null output
iflag=0
message=''

!Stream function
lam_3d(:,:,:)=0

!Stream function 09/09/11 SAL
lamc_3d(:,:,:)=0

!-------------------------------------------------------------------------------
!Load lambda expansion data
!-------------------------------------------------------------------------------
!Copy temporary radial grid and expansion coefficients
nset2=nk_lam

!Check whether inner boundary is at axis
IF(rho_lam(1) > 10*z_precision*rhorz) THEN

  !Radial node needs to be added at axis, check whether boundary node is present
  IF(rho_lam(nr_lam) < (1.0+10*z_precision)*rhorz) THEN

    !Add radial nodes at center and edge
    nset1=nr_lam+2
    ALLOCATE(rho(nset1), &     
             lmn(nset1,nset2))

      rho(:)=0
      lmn(:,:)=0

    rho(1)=0
    rho(2:nset1-1)=rho_lam(1:nr_lam)/rhorz
    rho(nset1)=1
    lmn(2:nset1-1,1:nset2)=lam(1:nr_lam,1:nset2)
    l_edge=.TRUE.
    
    !SAL 09/09/11
    IF (PRESENT(LAM_C2)) THEN
      ALLOCATE(lmnc(nset1,nset2))
        lmnc(:,:)=0
      lmnc(2:nset1-1,1:nset2)=LAM_C2(1:nr_lam,1:nset2)
    ENDIF

  ELSE

    !Add radial node at center
    nset1=nr_lam+1
    ALLOCATE(rho(nset1), &     
             lmn(nset1,nset2))

      rho(:)=0
      lmn(:,:)=0

    rho(1)=0
    rho(2:nset1)=rho_lam(1:nr_lam)/rhorz
    lmn(2:nset1,1:nset2)=lam(1:nr_lam,1:nset2)
    l_edge=.FALSE.
    
    !SAL 09/09/11
    IF (PRESENT(LAM_C2)) THEN
      ALLOCATE(lmnc(nset1,nset2))
        lmnc(:,:)=0
      lmnc(2:nset1-1,1:nset2)=LAM_C2(1:nr_lam,1:nset2)
    ENDIF

  ENDIF

ELSE

  !Radial node is present at axis
  !Check whether boundary node is present
  IF(rho_lam(nr_lam) < (1.0+10*z_precision)*rhorz) THEN

    !Add radial node at edge
    nset1=nr_lam+1
    ALLOCATE(rho(nset1), &     
             lmn(nset1,nset2))

      rho(:)=0
      lmn(:,:)=0

    rho(1:nset1-1)=rho_lam(1:nset1-1)/rhorz
    rho(nset1)=1
    lmn(1:nset1-1,1:nset2)=lam(1:nset1-1,1:nset2)
    l_edge=.TRUE.

    !SAL 09/09/11
    IF (PRESENT(LAM_C2)) THEN
      ALLOCATE(lmnc(nset1,nset2))
        lmnc(:,:)=0
      lmnc(1:nset1-1,1:nset2)=LAM_C2(1:nset1-1,1:nset2)
    ENDIF

  ELSE

    !Use input grid
    nset1=nr_lam
    ALLOCATE(rho(nset1), &     
             lmn(nset1,nset2))

      rho(:)=0
      lmn(:,:)=0

    rho(1:nset1)=rho_lam(1:nset1)/rhorz
    lmn(1:nset1,1:nset2)=lam(1:nset1,1:nset2)
    l_edge=.FALSE.

    !SAL 09/09/11
    IF (PRESENT(LAM_C2)) THEN
      ALLOCATE(lmnc(nset1,nset2))
        lmnc(:,:)=0
      lmnc(1:nset1,1:nset2)=LAM_C2(1:nset1,1:nset2)
    ENDIF

  ENDIF

ENDIF

!Convert rho to sqrt(toroidal flux) if necesssary and scale
IF(k_grid == 1) THEN

  !~toroidal flux
  rho(:)=SQRT(rho(:))*rhomax_3d

ELSE

  !~sqrt(toroidal flux)
  rho(:)=rho(:)*rhomax_3d

ENDIF

!Allocate temporary work arrays
ALLOCATE(fspl(4,nset1), &     
         values(3,nrho_3d))

  fspl(:,:)=0
  values(:,:)=0

!Normalize the expansion coeffs to rho**m
DO k=1,klam_3d !Over modes

  IF(l_edge) THEN

    !Extrapolate to edge
    DO j=2,nset1-1 !Over radial nodes

      lmn(j,k)=lmn(j,k)/rho(j)**mabs_3d(k)

    ENDDO !Over radial nodes

    lmn(nset1,k)=(lmn(nset1-1,k)*(rho(nset1)-rho(nset1-2)) & 
                 -lmn(nset1-2,k)*(rho(nset1)-rho(nset1-1))) & 
                  /(rho(nset1-1)-rho(nset1-2))
            
    !Extrapolate to edge SAL 09/09/11
    IF (PRESENT(LAM_C2)) THEN
      DO j=2,nset1-1 !Over radial nodes

        lmnc(j,k)=lmnc(j,k)/rho(j)**mabs_3d(k)

      ENDDO !Over radial nodes

      lmnc(nset1,k)=(lmnc(nset1-1,k)*(rho(nset1)-rho(nset1-2)) & 
                 -lmnc(nset1-2,k)*(rho(nset1)-rho(nset1-1))) & 
                  /(rho(nset1-1)-rho(nset1-2))
    ENDIF
    
  ELSE

    DO j=2,nset1 !Over radial nodes

      lmn(j,k)=lmn(j,k)/rho(j)**mabs_3d(k)

    ENDDO !Over radial nodes

    IF (PRESENT(LAM_C2)) THEN       !SAL 09/09/11
      DO j=2,nset1 !Over radial nodes
        lmnc(j,k)=lmnc(j,k)/rho(j)**mabs_3d(k)
      ENDDO !Over radial nodes
    ENDIF

  ENDIF

  !Make a parabolic fit to axis
  lmn(1,k)=(lmn(2,k)*rho(3)**2-lmn(3,k)*rho(2)**2)/(rho(3)**2-rho(2)**2)
  
  !Map lambda_mn to internal radial grid
  fspl(1,1:nset1)=lmn(1:nset1,k)
  iflag=0
  message=''
  CALL SPLINE1_INTERP(k_vopt,nset1,rho,fspl,nrho_3d,rho_3d,values, &
                      iflag,message, &    
                      K_BC1=k_bc1, &    
                      K_BCN=k_bcn)

  !Check messages
  IF(iflag > 0) THEN

    message='AJAX_LOAD_LAM/'//message
    GOTO 9999

  ENDIF

  !Respline the Lmn coefficients for internal storage (not-a-knot now)
  lam_3d(1,1:nrho_3d,k)=values(1,1:nrho_3d)
  CALL SPLINE1_FIT(nrho_3d,rho_3d,lam_3d(:,:,k), &  
                   K_BC1=k_bc1, &    
                   K_BCN=k_bcn)
  
  IF (PRESENT(LAM_C2)) THEN       !SAL 09/09/11
    !Make a parabolic fit to axis
    lmnc(1,k)=(lmnc(2,k)*rho(3)**2-lmnc(3,k)*rho(2)**2)&
             /(rho(3)**2-rho(2)**2)
    !Map lambda_mn to internal radial grid
    fspl(1,1:nset1)=lmnc(1:nset1,k)
    iflag=0
    message=''
    CALL SPLINE1_INTERP(k_vopt,nset1,rho,fspl,nrho_3d,rho_3d,values, &
                      iflag,message, &    
                      K_BC1=k_bc1, &    
                      K_BCN=k_bcn)

    !Check messages
    IF(iflag > 0) THEN
      message='AJAX_LOAD_LAM/'//message
      GOTO 9999
    ENDIF

    !Respline the Lmn coefficients for internal storage (not-a-knot now)
    lamc_3d(1,1:nrho_3d,k)=values(1,1:nrho_3d)
    CALL SPLINE1_FIT(nrho_3d,rho_3d,lamc_3d(:,:,k), &  
                   K_BC1=k_bc1, &    
                   K_BCN=k_bcn)
  END IF

ENDDO !Over modes

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

IF(ALLOCATED(rho)) THEN

  !Deallocate radial grid and lambda arrays
  DEALLOCATE(rho, &      
             lmn)

ENDIF

IF (ALLOCATED(lmnc)) THEN         !SAL 09/09/11
   DEALLOCATE(lmnc)
ENDIF

IF(ALLOCATED(fspl)) THEN

  !Deallocate spline arrays
  DEALLOCATE(fspl, &      
             values)

ENDIF

END SUBROUTINE AJAX_LOAD_LAMBDA

END MODULE AJAX_MOD
