MODULE TRACK_MOD
!-------------------------------------------------------------------------------
!TRACK-TRACKs line segments through a torus
!
!TRACK_MOD is an F90 module of routines to obtain plasma information along a
!  segmented path
!
!References:
!
!  Attenberger, Houlberg, Hirshman, J Comp Phys 72 (1987) 435
!  W.A.Houlberg, P.I.Strand, D.McCune 4/2002
!  W.A.Houlberg, F90 free format 8/2004
!
!Contains PUBLIC routines:
!
!  TRACK               -intersections of a segmented path with flux surfaces
!  TRACK_D             -coordinates of point a distance d along segment
!  TRACK_G             -Cartesian gradients along segment
!
!Comments:
!
!  Version 2.0 is a major upgrade of a previous package that has been in wide
!    use for about 15 years for RF, pellet, diagnostic and other tokamak and
!    stellarator applications.
!
!  The revisions include a faster and more robust algorithm for the
!    intersections.
!
!  This module is designed to map segmented straight-line trajectories in real
!    space to flux coordinates in order to find intersections with a
!    user-specified set of flux surfaces (TRACK).
!
!  Two auxiliary routines (TRACK_D and TRACK_G) are also accessible to the user
!    to find other information about a segment.
!
!  TRACK is meant to be used in transport analysis codes for mapping RF ray,
!    NBI, pellet, diagnostic or other trajectories to plasma coordinates for
!    plasma-neutral, plasma-wave, or other interactions that require information
!    from different coordinate systems.
!
!  The general tasks the TRACK modue performs are inherently 3D, so no
!    axisymmetric assumptions are made - that possibe restriction is left to the
!    MHD equilibrium interface module that actually performs the coordinate
!    conversions.
!
!  Two MHD equilibrium interface options are available:
!    AJAX              -2D or 3D inverse coordinate representations implemented
!                       by direct calls to AJAX routines
!    XPLASMA           -2D inverse coordinate or Psi(R,Z) representations
!                       implemented by calls to a module that maps AJAX calls to
!                       XPLASMA calls (AJAX_XPLASMA_MOD)
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

!Uncomment the appropriate MHD equilibrium interface, comment out the other
!For XPLASMA use mapping from XPLASMA to AJAX calls
!USE AJAX_XPLASMA_MOD

!For AJAX use direct calls
USE AJAX_MOD

IMPLICIT NONE

REAL(kind=rspec), PARAMETER :: one=1, zero=0

CONTAINS

SUBROUTINE TRACK(n_rho,rho,n_seg,r_seg, &
                 n_int,irho_int,s_int,iflag,message, &
                 K_SEG,IZONE_INT,RFLX_INT,RCYL_INT,RCAR_INT,SDOTB_INT, &
                 SDOTPHI_INT)
!-------------------------------------------------------------------------------
!TRACK finds the intersections of a segmented path with a set of flux surfaces
!
!References:
!  Attenberger, Houlberg, Hirshman, J Comp Phys 72 (1987) 435
!  W.A.Houlberg, P.I.Strand, D.McCune 4/2002
!  W.A.Houlberg, F90 free format 8/2004
!
!Comments:
!  It is assumed that rho(i) is increasing with i so that the i'th volume
!    element is bounded by rho(i) and rho(i+1)
!  The region outside rho(n_rho) is volume element n_rho
!-------------------------------------------------------------------------------

!Declaration of input variables
INTEGER, INTENT(IN) :: &
  n_rho,               & !number of radial nodes [-]
  n_seg                  !number of points defining the segments = segments+1 [-]

REAL(KIND=rspec), INTENT(IN) :: & 
  rho(:),              & !radial nodes [-]
  r_seg(:,:)             !coordinates defining the segments (see K_SEG)

!Declaration of optional input variables
INTEGER, INTENT(IN), OPTIONAL :: &
  K_SEG                  !option for coordinates specifying the segments [-]
                         !=1 flux coordinates (rho,theta,zeta)
                         !=2 Cartesian coordinates (x,y,z)
                         !=else default cylindrical coordinates (R,phi,Z)

!Declaration of output variables
CHARACTER(len=*), INTENT(OUT) :: &
  message                !warning or error message [character]

INTEGER, INTENT(OUT) :: &
  iflag,               & !error and warning flag [-]
                         !=-1 warning
                         !=0 none
                         !=1 error
  n_int,               & !number of nodes + intersections [-]
  irho_int(:)            !radial nodes of intersections [-]
                         !=0 if point is a node of the ray

REAL(KIND=rspec), INTENT(OUT) :: &
  s_int(:)               !length along path to intersections [m]

!Declaration of optional output variables
INTEGER, INTENT(OUT), OPTIONAL :: &
  IZONE_INT(:)           !zone being entered [-]

REAL(KIND=rspec), INTENT(OUT), OPTIONAL :: &
  RFLX_INT(:,:),       & !flux coordinates [rho,rad,rad]
  RCYL_INT(:,:),       & !cylindrical coordinates [m,rad,m]
  RCAR_INT(:,:),       & !Cartesian coordinates [m,m,m]
  SDOTB_INT(:),        & !cos of angle between path and B [-]
  SDOTPHI_INT(:)         !cos of angle between path and phi [-]

!-------------------------------------------------------------------------------
!Declaration of local variables
INTEGER :: &
  i,ilo,ihit,k,mhalf

REAL(KIND=rspec) :: &
  d0,d1,dseg,drds0,drds1,drds1f,drds_min,drho_min,ds,ds_max,ds_min,r000, &
  rhomax,sseg,tau,tol,tol1

REAL(KIND=rspec) :: &
  r_car0(1:3),r_cyl0(1:3),r_flx0(1:3),g_cyl0(1:6), &
  r_car1(1:3),r_cyl1(1:3),r_flx1(1:3),g_cyl1(1:6), &
  r_cars(1:3),g_cars(1:3),r_cylc(1:3,1:n_seg),b_car(1:3)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!General tolerence for min step sizes, resolution, convergence
tol=1.0e-4
tol1=1.0+tol

!Get rho and characterstic major radius
iflag=0
message=''
CALL AJAX_GLOBALS(iflag,message, & 
                  RHOMAX=rhomax, &
                  R000=r000)

!Check messages
IF(iflag /= 0) THEN

  message='TRACK(1)/'//message
  IF(iflag > 0) GOTO 9999

ENDIF

!Set limits on step sizes relative to scale of torus
ds_min=tol*r000
ds_max=0.02*r000
drho_min=tol*rhomax
drds_min=drho_min/ds_max

!Initialize local variables
drds0=0
drds1=0

!Initialize output arrays
irho_int(:)=0
s_int(:)=0
IF(PRESENT(IZONE_INT)) IZONE_INT(:)=0
IF(PRESENT(RFLX_INT)) RFLX_INT(1:3,:)=0
IF(PRESENT(RCYL_INT)) RCYL_INT(1:3,:)=0
IF(PRESENT(RCAR_INT)) RCAR_INT(1:3,:)=0
IF(PRESENT(SDOTB_INT)) SDOTB_INT(:)=0
IF(PRESENT(SDOTPHI_INT)) SDOTPHI_INT(:)=0

!Set option for coordinates defining segments
IF(PRESENT(K_SEG)) THEN

  k=K_SEG

ELSE

  k=0

ENDIF

!Set cylindrical coordinates defining segments
IF(k == 1) THEN

  !Flux to cylindrical conversion
  DO i=1,n_seg !Over chord points

    iflag=0
    message=''
    CALL AJAX_FLX2CYL(r_seg(1:3,i), &
                      r_cylc(1:3,i),iflag,message)

    !Check messages
    IF(iflag /= 0) THEN

      message='TRACK(2)/'//message
      IF(iflag > 0) GOTO 9999

    ENDIF

  ENDDO !Over chord points

ELSEIF(k == 2) THEN

  !Cartesian to cylindrical conversion
  DO i=1,n_seg !Over chord points

    CALL AJAX_CAR2CYL(r_seg(1:3,i),r_cylc(1:3,i))

  ENDDO !Over chord points

ELSE

  !Default cylindrical
  r_cylc(1:3,1:n_seg)=r_seg(1:3,1:n_seg)

ENDIF

!-------------------------------------------------------------------------------
!Set starting point of first segment
!-------------------------------------------------------------------------------
!Distance along entire path
sseg=0

!Cylindrical coordinates
r_cyl0(:)=r_cylc(:,1)

!Cartesian coordinates
CALL AJAX_CYL2CAR(r_cyl0,r_car0)

!Flux coordinates
r_flx0(1)=rhomax
r_flx0(2)=ATAN2(-r_cyl0(3),r_cyl0(1)-r000)
r_flx0(3)=r_cyl0(2)
iflag=0
message=''
CALL AJAX_CYL2FLX(r_cyl0, &
                  r_flx0,iflag,message, &
                  G_CYL=g_cyl0)

!Check messages
IF(iflag > 0) THEN

  !Assume starting point is outside plasma with bad Jacobian and proceed
  r_flx0(1)=1.2*rhomax
  r_flx0(2)=ATAN2(-r_cyl0(3),r_cyl0(1)-r000)
  r_flx0(3)=r_cyl0(2)
  g_cyl0(:)=(/one,zero,zero,zero,-1.2*rhomax,zero/)
  iflag=0
  message=''

ENDIF

!Find nearest flux surface between starting point and axis (may be on surface)
IF(r_flx0(1) > rho(n_rho)+drho_min) THEN

  !Outside largest surface of interest to user
  ilo=n_rho
  irho_int(1)=0

ELSE

  !Find radial starting radial zone or intersection
  LOOP_I: DO i=n_rho,1,-1 !Over radial nodes

    ilo=i

    IF(ABS(r_flx0(1)-rho(ilo)) <= drho_min) THEN

      !On a surface
      irho_int(1)=ilo
      EXIT LOOP_I

    ELSEIF(r_flx0(1) >= rho(ilo)) THEN

      !ilo is the lower surface
      irho_int(1)=0
      EXIT LOOP_I

    ENDIF

  ENDDO LOOP_I !Over radial nodes

ENDIF

!Output starting point information
n_int=1
s_int(n_int)=0
IF(PRESENT(IZONE_INT)) IZONE_INT(n_int)=ilo
IF(PRESENT(RFLX_INT)) RFLX_INT(1:3,n_int)=r_flx0(1:3)
IF(PRESENT(RCYL_INT)) RCYL_INT(1:3,n_int)=r_cyl0(1:3)
IF(PRESENT(RCAR_INT)) RCAR_INT(1:3,n_int)=r_car0(1:3)

CALL TRACK_G(r_cylc(1,1),r_cylc(1,2),dseg,g_cars)

IF(PRESENT(SDOTB_INT)) THEN

  iflag=0
  message=''
  b_car(:)=0
  CALL AJAX_B(r_flx0,r_cyl0,g_cyl0, &
              iflag,message, &
              B_CAR=b_car)

  !Check messages
  IF(iflag > 0) THEN

    message='TRACK(3)/ERROR:failed to get B components'
    GOTO 9999

  ENDIF

  IF(SUM(b_car(1:3)) /= 0.0) THEN

    SDOTB_INT(n_int)=SUM(g_cars(1:3)*b_car(1:3)) &
                     /SQRT(b_car(1)**2+b_car(2)**2+b_car(3)**2)
     
  ENDIF

ENDIF

IF(PRESENT(SDOTPHI_INT)) THEN

  SDOTPHI_INT(n_int)=-g_cars(1)*SIN(r_cyl0(2))+g_cars(2)*COS(r_cyl0(2))

ENDIF

!Initialize flux coordinates at 1
r_flx1(:)=r_flx0(:)

!-------------------------------------------------------------------------------
!Loop over segments
!-------------------------------------------------------------------------------
DO i=1,n_seg-1 !Over segments

  !Set Cartesian gradients of segment
  CALL TRACK_G(r_cylc(1,i),r_cylc(1,i+1),dseg,g_cars)

  !Cartesian coordinates of segment start
  CALL AJAX_CYL2CAR(r_cylc(:,i),r_cars)

  !Distance from start of segment
  d0=0

  !Cylindrical coordinates and other information
  iflag=0
  message=''
  CALL TRACK_D(d0,r_cars,g_cars, &
               r_flx0, &
               r_car0,r_cyl0,drds0,g_cyl0,tau,iflag,message)

  !Make sure gradient exceeds minimum
  IF(ABS(drds0) <= drds_min) drds0=SIGN(one,drds0)*tol1*drds_min

  !Check messages
  IF(iflag > 0 .OR. &
     tau < 0.0) THEN

    !Assume d0 is outside plasma with bad Jacobian and proceed
    r_flx0(:)=(/1.1*rhomax,zero,zero/)
    drds0=0
    iflag=0
    message=''

  ENDIF

  !Set step halving counter, and initialize position of point 1
  mhalf=0
  d1=0

!-------------------------------------------------------------------------------
!Step along segment: 0 is beginning of each step, 1 is end of step
!-------------------------------------------------------------------------------
  LOOP_K: DO k=1,1000 !Over steps in segment

    !Check limiting number of steps
    IF(k == 1000) THEN

      !Too many steps
      iflag=1
      message='TRACK(4)/ERROR:too many steps'
      GOTO 9999

    ENDIF

    !Check conditions for termination
    IF(mhalf > 5) THEN

      !Too many step halvings
      iflag=1
      message='TRACK(5)/ERROR:too many step halvings'
      GOTO 9999

    ENDIF

    IF(d1 > dseg-ds_min) EXIT LOOP_K

    !Set step size
    !Set ds from drho/ds at 0
    IF(ABS(drds0) < drds_min) THEN

      !Weak slope, limit step size to maximum
      ds=ds_max

    ELSEIF(drds0 < 0.0) THEN

      !Heading inward
      IF(ilo == 0) THEN

        !Special case where innermost surface is not axis
        ds=ds_max

      ELSE

        !Try to step past ilo
        ds=1.2*ABS((r_flx0(1)-rho(ilo))/drds0)

      ENDIF

    ELSE

      IF(ilo < n_rho) THEN

        !Heading outward, try to step past ilo+1
        ds=1.2*ABS((rho(ilo+1)-r_flx0(1))/drds0)

      ELSE

        !Outside plasma, limit step to maximum
        ds=ds_max

      ENDIF

    ENDIF

    !Check restrictions on step size
    ds=MIN(ds,dseg-d0,ds_max)
    ds=ds/REAL(2**mhalf,rspec)
    ds=MAX(ds,ds_min)

    !Set position in segment at end of step
    d1=d0+ds

    !Get parameters at 1
    r_flx1(:)=r_flx0(:)
    iflag=0
    message=''
    CALL TRACK_D(d1,r_cars,g_cars, &
                 r_flx1, &
                 r_car1,r_cyl1,drds1,g_cyl1,tau,iflag,message)

    !Make sure gradient exceeds minimum
    IF(ABS(drds1) <= drds_min) drds1=SIGN(one,drds1)*tol1*drds_min

    !Check messages
    IF(iflag > 0 .OR. &
       tau < 0.0) THEN

      !Assume d1 is outside plasma with bad Jacobian and proceed
      r_flx1(:)=(/1.1*rhomax,zero,r_cyl1(2)/)
      drds1=0
      iflag=0
      message=''

    ENDIF

!-------------------------------------------------------------------------------
!Checks for intersections and tangencies depend on drho/ds at endpoints
!-------------------------------------------------------------------------------
    IF(ABS(drds0) < drds_min) THEN

!-------------------------------------------------------------------------------
!Case I: Invalid soln at 0

      IF(ABS(drds1) < drds_min) THEN

!---------
!Case I.A: Invalid solns at 0 and 1
        !Proceed with next step
        ihit=0
        mhalf=0

      ELSE

        IF(r_flx1(1) > rho(n_rho)) THEN

!---------
!Case I.B&C.1: Outside target surface
          !Proceed with next step
          ihit=0
          mhalf=0

        ELSE

!---------
!Case I.B&C.2: Crossed target surface
          !Reduce step size to get good data point outside plasma
          ihit=0
          mhalf=1

        ENDIF

      ENDIF

    ELSEIF(drds0 > 0.0) THEN

!-------------------------------------------------------------------------------
!Case II: Valid outward soln at 0

      IF(ABS(drds1) < drds_min) THEN

!---------
!Case II.A: Valid outward soln at 0, invalid soln at 1

        IF(ilo == n_rho) THEN

!Case II.A.1: Outside plasma
          !Proceed with next step
          ihit=0
          mhalf=0

        ELSE

!Case II.A.2: Inside plasma
          !Reduce step to find intersection
          ihit=0
          mhalf=mhalf+1

        ENDIF

      ELSEIF(drds1 > 0.0) THEN

!---------
!Case II.B: Valid outward soln at 0, valid outward soln at 1
        IF(ilo == n_rho) THEN

!Case II.B.1: Outside plasma
          !Proceed with next step
          ihit=0
          mhalf=0

        ELSEIF(ABS(r_flx1(1)-rho(ilo+1)) < drho_min) THEN

!Case II.B.2: Hit surface ilo+1
          !Record hit and proceed with next step
          ihit=ilo+1
          mhalf=0
          ilo=ilo+1
          d1=d1+MAX(ds_min,(rho(ihit)-r_flx1(1))/drds1)
          ds=d1-d0

        ELSEIF(ilo < n_rho-1 .AND. &
               r_flx1(1) > rho(ilo+2)) THEN

!Case II.B.3: Crossed more than one surface
          !Halve step
          ihit=0
          mhalf=mhalf+1

        ELSEIF(r_flx1(1) > rho(ilo+1)) THEN

!Case II.B.4: Crossed surface ilo+1
          !Get intersection, record hit, and proceed with next step
          ihit=ilo+1
          mhalf=0
          ilo=ilo+1
          d1=d0+MAX(ds_min,(ds*(rho(ihit)-r_flx0(1)) &
             /(r_flx1(1)-r_flx0(1))))
          ds=d1-d0

        ELSE
!Case II.B.5: Crossed no surfaces
          !Proceed with next step
          ihit=0
          mhalf=0

        ENDIF

      ELSE

!---------
!Case II.C: Valid outward soln at 0, valid inward soln at 1, passed tangency
        !Get position of tangency
        d1=MAX(d0+ds_min,(d1*drds0-d0*drds1)/(drds0-drds1))
        ds=d1-d0

        !Get parameters at tangency
        r_flx1(:)=r_flx0(:)
        iflag=0
        message=''
        CALL TRACK_D(d1,r_cars,g_cars, &
                     r_flx1, &
                     r_car1,r_cyl1,drds1,g_cyl1,tau,iflag,message)

        !Check messages
        IF(iflag > 0 .OR. &
           tau < 0.0) THEN

          message='TRACK(6)/'//message
          GOTO 9999

        ENDIF

        !Force gradient to minimum beyond tangency
        drds1=-drds_min

        !Check for intersection at tangency
        IF(ilo == n_rho) THEN

!Case II.C.1: Tangency outside outermost surface
          !Proceed with next step
          ihit=0
          mhalf=0

        ELSEIF(ABS(r_flx1(1)-rho(ilo+1)) < drho_min) THEN

!Case II.C.2: Tangent to surface ilo+1 (degenerate case)
          !Record hit and proceed with next step
          ihit=ilo+1
          mhalf=0

        ELSEIF(r_flx1(1) > rho(ilo+1)) THEN

!Case II.C.3: Crossed surface ilo+1 before tangency
          !Get intersection and proceed with next step
          ihit=ilo+1
          mhalf=0
          ilo=ilo+1
          d1=d0+MAX(ds_min,(ds*(rho(ihit)-r_flx0(1)) &
             /(r_flx1(1)-r_flx0(1))))
          ds=d1-d0

          !Unforce gradient and proceed outward
          drds1=2*drds_min

        ELSE

!Case II.C.4: Crossed no surface
          !Proceed with next step
          ihit=0
          mhalf=0

        ENDIF

      ENDIF

    ELSE

!-------------------------------------------------------------------------------
!Case III: Valid inward soln at 0

      IF(ABS(drds1) < drds_min) THEN

!---------
!Case III.A: Valid inward soln at 0, invalid soln at 1
        IF(ilo == n_rho) THEN

!Case III.A.1: Outside plasma
          !Proceed with next step
          ihit=0
          mhalf=0

        ELSE

!Case III.A.2: Inside plasma
          !Failure
          iflag=1
          message='TRACK(7)/ERROR:conversion failure in plasma'
          GOTO 9999

        ENDIF

      ELSEIF(drds1 > 0) THEN

!---------
!Case III.B: Valid inward soln at 0, valid outward soln at 1, passed tangency
        !Get position of tangency
        d1=MAX(d0+ds_min,(d1*drds0-d0*drds1)/(drds0-drds1))
        ds=d1-d0

        !Get parameters at tangency
        r_flx1(:)=r_flx0(:)
        iflag=0
        message=''
        CALL TRACK_D(d1,r_cars,g_cars, &
                     r_flx1, &
                     r_car1,r_cyl1,drds1,g_cyl1,tau,iflag,message)

        !Force gradient to minimum beyond tangency
         drds1=drds_min

        !Check messages
        IF(iflag > 0 .OR. &
           tau < 0.0) THEN

          message='TRACK(8)/'//message
          GOTO 9999

        ENDIF

        !Check for intersection at tangency
        IF(ilo == 0) THEN

!Case III.B.1: Tangency inside innermost surface
          !Proceed with next step
          ihit=0
          mhalf=0

        ELSEIF(rho(ilo) < drho_min) THEN

!Case III.B.2: Don't look for tangency near axis
          !Proceed with next step
          ihit=0
          mhalf=0

        ELSEIF(ABS(r_flx1(1)-rho(ilo)) < drho_min) THEN

!Case III.B.3: Tangent to surface ilo (degenerate case)
          !Record hit and proceed with next step
          ihit=ilo
          mhalf=0

        ELSEIF(r_flx1(1) < rho(ilo)) THEN

!Case III.B.4: Crossed surface ilo before tangency
          !Get and record intersection and proceed with next step
          ihit=ilo
          mhalf=0
          ilo=ilo-1
          d1=d0+MAX(ds_min,(ds*(rho(ihit)-r_flx0(1)) &
             /(r_flx1(1)-r_flx0(1))))
          ds=d1-d0

          !Unforce gradient and proceed inward
          drds1=-2*drds_min

        ELSE

!Case III.B.5: Crossed no surface
          !Proceed with next step
          ihit=0
          mhalf=0

        ENDIF

      ELSE

!---------
!Case III.C: Valid inward soln at 0, valid inward soln at 1
        IF(ilo == 0) THEN

!Case III.C.1: Hasn't passed tangency to re-hit lowest surface
          !Proceed with next step
          ihit=0
          mhalf=0

        ELSEIF(rho(ilo) < drho_min) THEN

!Case III.C.2: Don't look for intersection near axis
          !Proceed with next step
          ihit=0
          mhalf=0

        ELSEIF(ABS(r_flx1(1)-rho(ilo)) < drho_min) THEN

!Case III.C.3: Hit surface ilo
          !Record hit and proceed with next step
          ihit=ilo
          mhalf=0
          ilo=ilo-1
          d1=d1+MAX(ds_min,(rho(ihit)-r_flx1(1))/drds1)
          ds=d1-d0

        ELSEIF(ilo > 1 .AND. &
               r_flx1(1) < rho(ilo-1)) THEN

!Case III.C.4: Crossed more than one surface
          !Halve step
          ihit=0
          mhalf=mhalf+1

        ELSEIF(r_flx1(1) < rho(ilo)) THEN

!Case III.C.5: Crossed surface ilo
          !Get intersection, record hit, and proceed with next step
          ihit=ilo
          mhalf=0
          ilo=ilo-1
          d1=d0+MAX(ds_min,(ds*(rho(ihit)-r_flx0(1)) &
             /(r_flx1(1)-r_flx0(1))))
          ds=d1-d0

        ELSE

!Case III.C.6: Crossed no surface
          !Proceed with next step
          ihit=0
          mhalf=0

        ENDIF

      ENDIF

    ENDIF

!-------------------------------------------------------------------------------
!Step completed
!-------------------------------------------------------------------------------

    IF(ihit /= 0) THEN

      !Save in case of forced end gradients
      drds1f=drds1

      !Record hit
      iflag=0
      message=''
      CALL TRACK_D(d1,r_cars,g_cars, &
                   r_flx1, &
                   r_car1,r_cyl1,drds1,g_cyl1,tau,iflag,message)

      IF(ABS(drds1) <= tol1*drds_min) drds1=drds1f

      !Check messages
      IF(iflag > 0 .OR. &
         tau < 0.0) THEN

        message='TRACK(9)/'//message
        GOTO 9999

      ENDIF

      n_int=n_int+1
      r_flx1(1)=rho(ihit)
      irho_int(n_int)=ihit
      s_int(n_int)=sseg+d1
      IF(PRESENT(IZONE_INT) .AND. &
         drds0 > 0.0) IZONE_INT(n_int)=ihit
      IF(PRESENT(IZONE_INT) .AND. &
         drds0 < 0.0) IZONE_INT(n_int)=ihit-1
      IF(PRESENT(RFLX_INT)) RFLX_INT(1:3,n_int)=r_flx1(1:3)
      IF(PRESENT(RCYL_INT)) RCYL_INT(1:3,n_int)=r_cyl1(1:3)
      IF(PRESENT(RCAR_INT)) RCAR_INT(1:3,n_int)=r_car1(1:3)

      IF(PRESENT(SDOTB_INT)) THEN

        iflag=0
        message=''
        b_car(:)=0
        CALL AJAX_B(r_flx1,r_cyl1,g_cyl1, &
                    iflag,message, &
                    B_CAR=b_car)

        !Check messages
        IF(iflag > 0) THEN

          message='TRACK(10)/ERROR:failed to get B components'
          GOTO 9999

        ENDIF

        SDOTB_INT(n_int)=SUM(g_cars(1:3)*b_car(1:3)) &
                         /SQRT(b_car(1)**2+b_car(2)**2+b_car(3)**2)

      ENDIF

      IF(PRESENT(SDOTPHI_INT)) THEN

        SDOTPHI_INT(n_int)=-g_cars(1)*SIN(r_cyl1(2))+g_cars(2) &
                           *COS(r_cyl1(2))

      ENDIF

      ihit=0

    ENDIF

    IF(mhalf == 0) THEN

      !Advance point 0 to point 1
      r_flx0(1:3)=r_flx1(1:3)
      r_cyl0(1:3)=r_cyl1(1:3)
      r_car0(1:3)=r_car1(1:3)
      g_cyl0(1:6)=g_cyl1(1:6)
      drds0=drds1
      d0=d1

    ENDIF

  ENDDO LOOP_K !Over steps in segment
   
!-------------------------------------------------------------------------------
!Segment completed
!-------------------------------------------------------------------------------
  sseg=sseg+dseg
  n_int=n_int+1
  irho_int(n_int)=0
  s_int(n_int)=sseg
  r_cyl1(1:3)=r_cylc(1:3,i+1)

  iflag=0
  message=''
  CALL AJAX_CYL2FLX(r_cyl1, &
                    r_flx1,iflag,message, &
                    G_CYL=g_cyl1)

  !Check messages
  IF(iflag > 0) THEN

    !Assume d1 is outside plasma with bad Jacobian and proceed
    r_flx1(:)=(/1.1*rhomax,zero,zero/)
    drds1=0
    iflag=-1
    message='TRACK(11)/WARNING:end point assumed outside plasma'

  ENDIF

  CALL AJAX_CYL2CAR(r_cyl1,r_car1)

  IF(PRESENT(IZONE_INT)) IZONE_INT(n_int)=ilo
  IF(PRESENT(RFLX_INT)) RFLX_INT(1:3,n_int)=r_flx1(1:3)
  IF(PRESENT(RCYL_INT)) RCYL_INT(1:3,n_int)=r_cyl1(1:3)
  IF(PRESENT(RCAR_INT)) RCAR_INT(1:3,n_int)=r_car1(1:3)

  IF(PRESENT(SDOTB_INT)) THEN

    iflag=0
    message=''
    CALL AJAX_B(r_flx1,r_cyl1,g_cyl1, &
                iflag,message, &
                B_CAR=b_car)

    !Check messages
    IF(iflag > 0) THEN

      message='TRACK(12)/ERROR:failed to get B components'
      GOTO 9999

    ENDIF

    SDOTB_INT(n_int)=SUM(g_cars(1:3)*b_car(1:3)) &
                     /SQRT(b_car(1)**2+b_car(2)**2+b_car(3)**2)

  ENDIF

  IF(PRESENT(SDOTPHI_INT)) THEN

    SDOTPHI_INT(n_int)=-g_cars(1)*SIN(r_cylc(2,i+1)) &
                       +g_cars(2)*COS(r_cylc(2,i+1))

  ENDIF

ENDDO !Over segments

!Make sure intersection distances are monotonic
DO k=1,n_int-1

   IF(s_int(k+1) <= s_int(k)) THEN

      iflag=1
      message='TRACK(12)/ERROR:intercept lengths out of order.'
      GOTO 9999

   ENDIF

ENDDO

!-------------------------------------------------------------------------------
!Cleanup and exit
!-------------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE TRACK

SUBROUTINE TRACK_D(d,r_cars,g_cars, &
                   r_flx, &
                   r_car,r_cyl,drhods,g_cyl,tau,iflag,message)
!-------------------------------------------------------------------------------
!TRACK_D finds the Cartesian, cylindrical and flux coordinates of a point
!  located a distance d along a line segment specified in Cartesian coordinates,
!  and returns other relevant information for that point
!
!References:
!  Attenberger, Houlberg, Hirshman J Comp Phys 72 (1987) 435
!  W.A.Houlberg, P.I.Strand 11/2001
!  W.A.Houlberg, F90 free format 8/2004
!
!Comments:
!  For efficiency, r_flx should be a good initial guess on entry
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  d,                   & !distance along segment [m]
  r_cars(:),           & !Cartesian coordinates of start of segment [m]
  g_cars(:)              !unit length vector along segment [-]

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
  r_car(:),            & !Cartesian coordinates (x,y,z) [m,m,m]
  r_cyl(:),            & !cylindrical coordinates (R,phi,Z) [m,rad,m]
  g_cyl(:),            & !R,Z derivatives
                         !=(R_rho,R_theta,R_zeta,Z_rho,Z_theta,Z_zeta)
                         ! [m/rho,m,      m,     m/rho,m,      m     ]
  drhods,              & !drho/ds at position d along the chord [rho/m]
  tau                    !2-D Jacobian in phi=zeta=constant plane [m**2/rho]

!-------------------------------------------------------------------------------
!Declaration of local variables
REAL(KIND=rspec) :: &
  rhop,rhor,rhoz

!-------------------------------------------------------------------------------
!Get coordinates
!-------------------------------------------------------------------------------
!Cartesian coordinates of end point
r_car(1:3)=r_cars(1:3)+g_cars(1:3)*d

!Cylindrical coordinates
CALL AJAX_CAR2CYL(r_car, &
                  r_cyl)

!Flux coordinates
iflag=0
message=''
CALL AJAX_CYL2FLX(r_cyl, &
                  r_flx,iflag,message, &
                  G_CYL=g_cyl, &
                  TAU=tau)

!Check messages
IF(iflag >0) THEN

  message='TRACK_D/'//message

ENDIF

!drho/ds
rhor=-g_cyl(5)/tau
rhoz=g_cyl(2)/tau
rhop=(g_cyl(3)*g_cyl(5)-g_cyl(2)*g_cyl(6))/tau
drhods=(rhor*r_car(1)-rhop*r_car(2)/r_cyl(1))*g_cars(1)/r_cyl(1) &
       +(rhor*r_car(2)+rhop*r_car(1)/r_cyl(1))*g_cars(2)/r_cyl(1)+rhoz*g_cars(3)

END SUBROUTINE TRACK_D

SUBROUTINE TRACK_G(r_cyl0,r_cyl1, &
                   d,g_car)
!-------------------------------------------------------------------------------
!TRACK_G evaluates the Cartesian gradients along a segment defined by two sets
!  of cylindrical coordinates
!
!References:
!  Attenberger, Houlberg, Hirshman, J Comp Phys 72 (1987) 435
!  W.A.Houlberg 8/2001
!  W.A.Houlberg, F90 free format 8/2004
!-------------------------------------------------------------------------------

!Declaration of input variables
REAL(KIND=rspec), INTENT(IN) :: &
  r_cyl0(3),           & !cylindrical coordinates at start of segment [m,rad,m]
  r_cyl1(3)              !cylindrical coordinates at end of segment [m,rad,m]

!Declaration of output variables
REAL(KIND=rspec), INTENT(OUT) :: &
  d,                   & !length of segment [m]
  g_car(3)               !unit length vector along segment [-]

!-------------------------------------------------------------------------------
!Declaration of local variables
REAL(KIND=rspec) :: &
  r_car0(3),r_car1(3)

!-------------------------------------------------------------------------------
!Get gradients
!-------------------------------------------------------------------------------
!Cartesian coordinates for endpoints
CALL AJAX_CYL2CAR(r_cyl0,r_car0)
CALL AJAX_CYL2CAR(r_cyl1,r_car1)

!Increments in Cartesian coordinates
g_car(:)=r_car1(:)-r_car0(:)

!Segment length
d=SQRT((r_car1(1)-r_car0(1))**2+(r_car1(2)-r_car0(2))**2 &
       +(r_car1(3)-r_car0(3))**2)

!Normalized gradients
g_car(:)=g_car(:)/d

END SUBROUTINE TRACK_G

END MODULE TRACK_MOD
