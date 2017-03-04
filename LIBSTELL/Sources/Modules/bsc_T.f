!*******************************************************************************
!  File bsc_mod.f
!  Contains module bsc
!  Under CVS control at logjam.gat.com

!*******************************************************************************
!  MODULE bsc
!    (Biot-Savart Coils)
! SECTION I.    VARIABLE DECLARATIONS
! SECTION II.   DERIVED TYPE (STRUCTURE) DECLARATIONS
! SECTION III.  INTERFACE BLOCKS
! SECTION IV.   CONSTRUCTION SUBROUTINES
! SECTION V.    DESTRUCTION SUBROUTINES
! SECTION VI.   ASSIGNMENT SUBROUTINES
! SECTION VII.  COIL COLLECTION MANIPULATION SUBROUTINES
! SECTION VIII. ROTATE AND SHIFT SUBROUTINES
! SECTION IX.   VECTOR POTENTIAL SUBROUTINES
! SECTION X.    MAGNETIC FIELD SUBROUTINES
! SECTION XI.   MAGNETIC FLUX INTEGRAL SUBROUTINES
! SECTION XII.  AUXILIARY SUBROUTINES
! SECTION XIII. DEBUGGING SUBROUTINES
! SECTION XIV.  SPECIAL FUNCTIONS
! SECTION XV.   DUPLICATE CODING FOR TESTING
! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
      MODULE bsc_T

      IMPLICIT NONE
!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Type declarations - lengths of reals, integers, and complexes.
!-------------------------------------------------------------------------------
      INTEGER, PARAMETER :: rprec = SELECTED_REAL_KIND(12,100)
      INTEGER, PARAMETER :: iprec = SELECTED_INT_KIND(8)
      INTEGER, PARAMETER :: cprec = KIND((1.0_rprec,1.0_rprec))

!-------------------------------------------------------------------------------
!  Frequently used mathematical constants, lots of extra precision.
!-------------------------------------------------------------------------------
      REAL(rprec), PARAMETER :: pi=3.14159265358979323846264338328_rprec
      REAL(rprec), PARAMETER :: twopi=6.28318530717958647692528677_rprec
      REAL(rprec), PARAMETER :: one = 1.0_rprec, zero = 0.0_rprec

!-------------------------------------------------------------------------------
!  Physical Constants
!-------------------------------------------------------------------------------
      REAL(rprec), PARAMETER :: bsc_k2_def = 1.0e-7_rprec
      REAL(rprec), PARAMETER :: bsc_k2inv_def = 1.0e+7_rprec

!-------------------------------------------------------------------------------
!  Computational Constants
!-------------------------------------------------------------------------------
      REAL(rprec), PARAMETER :: bsc_mach_eps = EPSILON(one)
      
!-------------------------------------------------------------------------------
!  Make type declarations and constants Private, so there are no conflicts.
!-------------------------------------------------------------------------------
      PRIVATE rprec, iprec, cprec, pi, twopi, one, zero, bsc_k2_def,           &
     &   bsc_k2inv_def, bsc_mach_eps

!-------------------------------------------------------------------------------
!  Tuning Parameters (Private)
!-------------------------------------------------------------------------------
!    These are parameters that are pretty well determined, so they
!    can be made private. They are not available outside the bsc module.
!    For testing purposes, comment out the PRIVATE statement.
      REAL(rprec) :: bsc_emcut = 0.01_rprec

!  bsc_emcut    cut off value of m,for A and B for circular coils
!    switches between power series and elliptic integrals
!  put it here, so that can change in test code

!      PRIVATE bsc_emcut
      
!-------------------------------------------------------------------------------
!  Tuning Parameters (Not Private)
!-------------------------------------------------------------------------------
!    These are parameters that are I want to be able to change
!    from outside the module, for testing purposes, or for tuning
!    the algorithms.

! (Nothing here right now)

!*******************************************************************************
! SECTION II. DERIVED TYPE (STRUCTURE) DECLARATIONS
!   Type to describe all coils:
!     bsc_coil  
!   Type of coil specified by  % c_type.
!   Allowable values of c_type:
!     fil_loop  filamentary loop, consisting of straight line segments
!     fil_circ  filamentary circular loop, uses complete elliptic integrals
!     floop     Same as fil_loop. Use is deprecated.
!     fcirc     Same as fil_circ. Use is deprecated.
!     fil_rogo  Filamentary (partial) Rogowski coil
!               (other c_types not yet implemented)
!
!   Coil Collection
!     bsc_coilcoll  collection of bsc_coil's
!
!   Structure for Rotate and Shift Information
!     bsc_rs
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Declare type bsc_coil
!    Common to all c_types                                                     
!       c_type   character, type of coil
!       s_name   character, short name of coil
!       l_name   character, long name of coil
!       current  current in the coil
!       eps_sq   pedestal to avoid numerical singularity if observation point lies on coil
!       raux     real auxiliary variable. Carried by bsc, but not used within bsc.
!    Used for c_type = 'fil_circ'
!       rcirc    radius of circle
!       xcent    array(3), Cartesian center of circle
!       enhat    array(3), unit vector, normal to plane of circle
!    Used for c_type = 'fil_loop'
!       xnod     array(3,-) of cartesian node positions
!       dxnod    array(3,-) of xnod differences
!       ehnod    array(3,-) of ehat vectors (normalized dxnod)
!       lsqnod   array(-) of square lengths of straight line segments (dxnod)
!       lnod     array(-) of lengths of straight line segments (dxnod)
!    Used for c_type = 'fil_rogo'
!       ave_n_area  average value of number of turns per unit length times
!                   cross-sectional area of turns
!-------------------------------------------------------------------------------
      TYPE bsc_coil
         CHARACTER (len=8) :: c_type
         CHARACTER (len=30) :: s_name                                 
         CHARACTER (len=80) :: l_name
         REAL(rprec) :: eps_sq
         REAL(rprec) :: current
         REAL(rprec) :: raux
         REAL(rprec) :: rcirc
         REAL(rprec) :: ave_n_area
         REAL(rprec), DIMENSION(3) :: xcent, enhat
         REAL(rprec), DIMENSION(:,:), POINTER :: xnod => null()
         REAL(rprec), DIMENSION(:,:), POINTER :: dxnod => null()
         REAL(rprec), DIMENSION(:,:), POINTER :: ehnod => null()
         REAL(rprec), DIMENSION(:),   POINTER :: lsqnod => null()
         REAL(rprec), DIMENSION(:),   POINTER :: lnod => null()
      END TYPE bsc_coil
!-------------------------------------------------------------------------------
!  Declare type bsc_coilcoll                                                    
!       s_name   character, short name of coil collection
!       l_name   character, long name of coil collection
!       ncoil    number of coils
!       coils    array of coils
!-------------------------------------------------------------------------------
      TYPE bsc_coilcoll
         CHARACTER (len=30) :: s_name
         CHARACTER (len=80) :: l_name
         INTEGER(iprec) :: ncoil
         TYPE (bsc_coil), DIMENSION(:), POINTER :: coils => null()    
      END TYPE bsc_coilcoll
!-------------------------------------------------------------------------------
!  Declare type bsc_rs
!       rot_matrix  real, array(3,3), rotation matrix for the specific coil
!       c_of_rot    real, array(3), center of rotation for the specified coil 
!       shift       real, array(3), shift vector for the specified coil
!-------------------------------------------------------------------------------
      TYPE bsc_rs
         REAL(rprec), DIMENSION(3,3) :: rot_matrix
         REAL(rprec), DIMENSION(3)   :: c_of_rot
         REAL(rprec), DIMENSION(3)   :: shift
      END TYPE bsc_rs
!*******************************************************************************
! SECTION III. INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Assignment for bsc_coil
!-------------------------------------------------------------------------------
      INTERFACE ASSIGNMENT (=)
         MODULE PROCEDURE bsc_coil_to_coil, bsc_coil_a_to_coil_a
      END INTERFACE

!-------------------------------------------------------------------------------
!  Generic construct
!-------------------------------------------------------------------------------
      INTERFACE bsc_construct
         MODULE PROCEDURE bsc_construct_coil, bsc_construct_coilcoll,          &
     &                    bsc_construct_rs
      END INTERFACE

!-------------------------------------------------------------------------------
!  Generic destroy
!-------------------------------------------------------------------------------
      INTERFACE bsc_destroy
         MODULE PROCEDURE bsc_destroy_coil, bsc_destroy_coilcoll,              &
     &                    bsc_destroy_coil_a
      END INTERFACE

!-------------------------------------------------------------------------------
!  Generic bsc_rot_shift (Rotation and Shift)
!-------------------------------------------------------------------------------
      INTERFACE bsc_rot_shift
         MODULE PROCEDURE bsc_rot_shift_pt, bsc_rot_shift_pts,                 &
     &                    bsc_rot_shift_coil,                                  & 
     &                    bsc_rot_shift_coil_a, bsc_rot_shift_coilcoll
      END INTERFACE

!-------------------------------------------------------------------------------
!  Generic bsc_a (Vector Potential)
!-------------------------------------------------------------------------------
      INTERFACE bsc_a
         MODULE PROCEDURE bsc_a_coil, bsc_a_coil_a, bsc_a_coilcoll
      END INTERFACE

!-------------------------------------------------------------------------------
!  Generic bsc_b (Magnetic Field)
!-------------------------------------------------------------------------------
      INTERFACE bsc_b
         MODULE PROCEDURE bsc_b_coil, bsc_b_coil_a, bsc_b_coilcoll
      END INTERFACE

!-------------------------------------------------------------------------------
!  Generic bsc_fluxba (Magnetic Flux of a through coil b)
!-------------------------------------------------------------------------------
      INTERFACE bsc_fluxba
         MODULE PROCEDURE bsc_fluxba_coil, bsc_fluxba_coil_a,                  &
     &   bsc_fluxba_coilcoll
      END INTERFACE
               
!-------------------------------------------------------------------------------
!  Interface block for testing goes here. See SECTION XII.
!-------------------------------------------------------------------------------

      CONTAINS
!*******************************************************************************
! SECTION IV. CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a coil
!
!  For c_type = 'fil_loop' (filamentary loops)
!    Note that not all the structure components are needed as arguments.
!      dxnod, ehnod, lsqnod, and lnod are computed during construction.
!
!    Note that the number of nodes in the coil is determined by the SIZE
!      of xnod, as passed into bsc_construct_coil.
!
!    Note that the coil as input is assumed to be NOT CLOSED.
!    JDH To Do: add coding to CHECK this. Want to avoid closing a loop
!       that is already closed. Better yet, why not prune out very short
!       segments, say < 10^-8 times total length.
!
!  For c_type = 'fil_circ' (filamentary circular loops)
!     Unit length for enhat is enforced on construction
!
!  For c_type = 'fil_rogo' (filamentary Rogowski coils)
!     Much of the coding is the same as for fil_loop coils.
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_construct_coil(this,c_type,s_name,l_name,current, 
     &   xnod,rcirc,xcent,enhat,raux,anturns,xsarea)

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (bsc_coil), INTENT (inout) :: this
      CHARACTER (len=*), INTENT(in)   :: c_type
      CHARACTER (len=*), INTENT(in)   :: s_name
      CHARACTER (len=*), INTENT(in)   :: l_name
      REAL(rprec), INTENT(in)  :: current
      REAL(rprec), DIMENSION(:,:), INTENT(in), OPTIONAL :: xnod
      REAL(rprec), INTENT(in), OPTIONAL :: rcirc
      REAL(rprec), DIMENSION(3), INTENT(in), OPTIONAL :: xcent
      REAL(rprec), DIMENSION(3), INTENT(in), OPTIONAL :: enhat
      REAL(rprec), INTENT(in), OPTIONAL :: raux

!  Declare Arguments that aren't bsc_coil components 
      REAL(rprec), INTENT(in), OPTIONAL :: anturns
      REAL(rprec), INTENT(in), OPTIONAL :: xsarea
!  anturns    Total number of turns in Rogowski coil. 
!              NB. Declared as REAL, not INTEGER
!  xsarea    Cross-sectional area of turns in Rogowski coil

!  Declare local variables
      INTEGER(iprec) :: n_xnod_1, n_xnod_2, i, n, nm1
      REAL(rprec) :: enlength

!  Local Variables added 2010-07-06 JDH      
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: xnod_temp
      REAL(rprec), DIMENSION(3) :: vec_temp
      INTEGER(iprec) :: itemp
      REAL(rprec) :: lsqnod_temp

!  Start of executable code

!      WRITE(*,*) ' Executing bsc_construct_coil'

!  First, destroy the coil
      CALL bsc_destroy(this)

!  Scalar assignments, common to all c_types
      this % s_name = s_name
      this % l_name = l_name
      this % current = current
      IF (PRESENT(raux)) THEN
         this % raux = raux
      ELSE
         this % raux = zero
      END IF

!  Different coding, depending on c_type
      SELECT CASE (c_type)
      CASE ('fil_loop','floop','fil_rogo')
         this % c_type = c_type

!  Check for presence of necessary arguments
         IF (.not. PRESENT(xnod)) THEN
            WRITE(*,*) 'FATAL: bsc_construct_coil '
            WRITE(*,*) 'argument xnod not present for c_type =', c_type
            STOP
         END IF

!  Check lengths of xnod
         n_xnod_1 = SIZE(xnod,1)
         n_xnod_2 = SIZE(xnod,2)
         IF (n_xnod_1 .ne. 3) THEN
            STOP ' FATAL:sub. bsc_construct: n_xnod_1 .ne. 3'
         ELSE IF (n_xnod_2 .lt. 2) THEN
            STOP ' FATAL:sub. bsc_construct: n_xnod_2 < 2'
         ENDIF

!  JDH 2010-07-06
!  First, Eliminate all zero-length segments.
         ALLOCATE(xnod_temp(3,n_xnod_2 + 1))
         
         itemp = 1
         xnod_temp(1:3,1) = xnod(1:3,1)
         DO i = 2, n_xnod_2
            vec_temp(1:3) = xnod_temp(1:3,itemp) - xnod(1:3,i)
            lsqnod_temp = DOT_PRODUCT(vec_temp,vec_temp)
            IF (lsqnod_temp .eq. zero) THEN
               CYCLE
            ELSE
               itemp = itemp + 1
               xnod_temp(1:3,itemp) = xnod(1:3,i)
            ENDIF
         END DO

!  Check that wrap-around segment is not zero length
!  (This effectively unwraps the loop, if it came in wrapped)
!  JDH 2012-01-23. Make sure this does NOT apply to Rogowskis
!    (because the addition of the wrapping segment only happens for
!    fil_loops)
         IF ((c_type .eq.'fil_loop') .or. (c_type .eq.'floop')) THEN
            vec_temp(1:3) = xnod_temp(1:3,itemp) - xnod(1:3,1)
            lsqnod_temp = DOT_PRODUCT(vec_temp,vec_temp)
            IF (lsqnod_temp .eq. zero) THEN
               itemp = itemp - 1
            ENDIF
         ENDIF
         
         IF (itemp .eq. 1) THEN
            STOP ' FATAL:sub. bsc_construct: itemp .eq. 1'
         ENDIF
         
!  Close (wrap) fil_loops
         IF ((c_type .eq.'fil_loop') .or. (c_type .eq.'floop')) THEN
            IF (itemp .eq. 2) THEN
               WRITE (6,*) ' WARNING: bsc_construct: one straight ',           &
     &                  ' filament for coil: ', TRIM(s_name)
               WRITE (6,*) ' (this is a straight coil filament)'
            ELSE  !  Here - need to close (wrap) the coil
               itemp = itemp + 1
               xnod_temp(1:3,itemp) = xnod(1:3,1)
            ENDIF
         ENDIF
         
         n = itemp
         nm1 = n - 1
!  JDH 2010-07-06

!  Allocate space for pointers. No need to deallocate space
!  as this was done in bsc_destroy (called at
!  start of bsc_construct_coil)
         ALLOCATE(this % xnod(3,n))    
         ALLOCATE(this % dxnod(3,nm1))           
         ALLOCATE(this % ehnod(3,nm1))           
         ALLOCATE(this % lsqnod(nm1))            
         ALLOCATE(this % lnod(nm1))              

!  Copy xnod to this % xnod.
!  Modified 2010-07-06 JDH
!         this % xnod(1:3,1:nm1) = xnod(1:3,1:nm1)
         this % xnod(1:3,1:n) = xnod_temp(1:3,1:n)
         DEALLOCATE(xnod_temp)
        
!  Calculations for the other arrays (not included as arguments)
!  Compute dxnod
         this % dxnod(1:3,1:nm1) = this % xnod(1:3,2:n)                        &
     &                           - this % xnod(1:3,1:nm1)

!  Compute lsqnod = dxnod*dxnod
         this % lsqnod(1:nm1) = this % dxnod(1,1:nm1)**2 +                     &
     &      this % dxnod(2,1:nm1)**2 + this % dxnod(3,1:nm1)**2          
         IF (ANY(this % lsqnod(1:nm1) .eq. zero))                              &
     &      STOP 'FATAL: bsc_construct_coil : lsqnod must be nonzero'
         this % eps_sq = bsc_mach_eps *                                        &
     &                   MINVAL(this % lsqnod(1:nm1))
!  JDH 11-21-03. Change EPSILON(-) to bsc_mach_eps
!  JDH 11-21-03. Not sure just how this % eps_sq will get used.

!  Compute lnod
         this % lnod(1:nm1) = SQRT(this % lsqnod(1:nm1))              

!  Compute ehnod
         DO i = 1,3
            this % ehnod(i,1:nm1) = this % dxnod(i,1:nm1) /                    &                     
     &         this % lnod(1:nm1)                                     
         END DO
      
      CASE ('fil_circ','fcirc')
         this % c_type = c_type

!  Check for presence of necessary arguments
         IF (.not. PRESENT(rcirc)) THEN
            WRITE(6,*) 'FATAL: bsc_construct_coil '
            WRITE(6,*) 'arg rcirc not present for c_type =', c_type
            STOP
         ELSE IF (.not. PRESENT(xcent)) THEN
            WRITE(6,*) 'FATAL: bsc_construct_coil '
            WRITE(6,*) 'arg xcent not present for c_type =', c_type
            STOP
         ELSE IF (.not. PRESENT(enhat)) THEN
            WRITE(6,*) 'FATAL: bsc_construct_coil '
            WRITE(6,*) 'arg enhat not present for c_type =', c_type
            STOP
         END IF

!  Copy arguments
         this % rcirc = rcirc
         this % xcent(1:3) = xcent(1:3)

!  Normalize the enhat unit vector
!  ZZZ Think about 1.d-40 number here
         enlength = SQRT(DOT_PRODUCT(enhat,enhat))
         IF (enlength .le. 1.e-40_rprec) THEN
            this % enhat(1) = 0            ! JDH Integers get converted correctly
            this % enhat(2) = 0
            this % enhat(3) = 1
            WRITE(*,*) 'WARN: bsc_contruct: enhat set to (0,0,1)'
         ELSE
            this % enhat = enhat / enlength
         END IF ! normalize the enhat vector

      CASE DEFAULT
         WRITE(*,*) 'FATAL: bsc_contruct: unrecognized c_type = ',c_type
         STOP
      END SELECT ! Different coding depending on c_type
      
!  More stuff, for Rogowskis
      IF (this % c_type .eq. 'fil_rogo') THEN
         IF (.not. PRESENT(anturns)) THEN
            WRITE(6,*) 'FATAL: bsc_construct_coil '
            WRITE(6,*) 'arg anturns not present for c_type =', c_type
            STOP
         ELSE IF (.not. PRESENT(xsarea)) THEN
            WRITE(6,*) 'FATAL: bsc_construct_coil '
            WRITE(6,*) 'arg xsarea not present for c_type =', c_type
            STOP
         END IF
         
         this % ave_n_area = xsarea * anturns / sum(this % lnod(1:nm1))
      END IF

      END SUBROUTINE bsc_construct_coil

!-------------------------------------------------------------------------------
!  Construct a coil collection: bsc_coilcoll
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_construct_coilcoll(this,s_name,l_name,ncoil_init)
      IMPLICIT NONE

!  Required Arguments 
      TYPE (bsc_coilcoll), INTENT (inout) :: this
      CHARACTER (len=*), INTENT(in)  :: s_name
      CHARACTER (len=*), INTENT(in)  :: l_name

!  Optional Arguments
!    ncoil_init allows user to override default initial size of this % coils array
      INTEGER(iprec), INTENT(in), OPTIONAL :: ncoil_init

!  Declare local variables
      INTEGER(iprec), PARAMETER :: ncoil_init_def = 10
      INTEGER(iprec) :: ncoil_init_use

!  Start of executable code
!  Check coil pointer. If it is associated, then destroy new_coilcoll
!  and start over

      IF (ASSOCIATED(this % coils)) CALL bsc_destroy(this)
      
!  Initial assignments and allocation
      this % s_name = s_name
      this % l_name = l_name
      this % ncoil = 0

      IF (PRESENT(ncoil_init)) THEN
         ncoil_init_use = MAX(2,ncoil_init)
      ELSE
         ncoil_init_use = ncoil_init_def
      ENDIF

      ALLOCATE(this % coils(ncoil_init_use))
      
      END SUBROUTINE bsc_construct_coilcoll

!-------------------------------------------------------------------------------
!  Construct a coil rotation and shift type bsc_rs  
!
!  Required Arguments
!  this            :   bsc_rs type to create. 
!                      on exit, contains the rotation matrix, the center of
!                      rotation, and the shift vector
!  theta           :   real(rprec), theta angle in spherical coordinates to 
!                      indicate the direction of the rotation axis vector                 
!  phi             :   real(rprec), phi angle in spherical coordinates to 
!                      indicate the direction of the rotation axis vector
!  rot_ang         :   real(rprec), angle specifying rigid-body rotation with 
!                      respect to the axis of rotation (left-hand rule).
!
!  Optional Arguments 
!  c_of_rot        :   real(rprec), array (size 3) center of rigid-body center 
!                      of mass(com) shifts   
!  shift           :   real(rprec), array (size 3) shift vector for the 
!                      translation of the coil   
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_construct_rs(this,theta,phi,rot_ang,                      &
     &                                 c_of_rot,shift)
      IMPLICIT NONE

!  Argument Declaration
!  Required Arguments
      TYPE (bsc_rs), INTENT(inout) :: this
      REAL(rprec), INTENT(in) :: theta
      REAL(rprec), INTENT(in) :: phi
      REAL(rprec), INTENT(in) :: rot_ang

!  Optional Arguments
      REAL(rprec), DIMENSION(3), OPTIONAL, INTENT(in) :: c_of_rot
      REAL(rprec), DIMENSION(3), OPTIONAL, INTENT(in) :: shift
      
!  Local Variable Declaration
      REAL(rprec) :: omega(3)
      REAL(rprec) :: cosrot, sinrot, onemcos
  
!  Start of executable code
!*******************************************************************************
!     Apply rotation about an axis of rotation formula 
!     (Ref: "Classical Mechanics", Goldstein,  (first ed. P-162)
!                                              (second ed. pp 164-65)
!
!     x(rot) = [x(in) dot OMEGA]OMEGA  
!            + cos(rot_ang)[X(in) - (X(in) dot OMEGA)OMEGA] 
!            + sin(rot_ang) X(in) cross OMEGA
!
!            == R * X(in)
!
!     where OMEGA = (sin(theta)cos(phi) xhat + sin(theta)*sin(phi) yhat + 
!                    cos(theta) zhat
!     is the unit rotation axis vector. Note that only the component of X(in) 
!     NORMAL to OMEGA is rotated.
!
!     NOTE: R(inv) = R(transpose) = R(-rot_ang)
!     NOTE: Goldstein's convention is a "left-handed" one. If you point your _left_
!        thumb along the OMMEGA vector, the fingers of your left hand indicate
!        the direction of rotation.
!
!*******************************************************************************

      omega(1) = SIN(theta) * COS(phi)
      omega(2) = SIN(theta) * SIN(phi)
      omega(3) = COS(theta)

      cosrot = COS(rot_ang);  sinrot = SIN(rot_ang)
      onemcos = 1 - cosrot

      this % rot_matrix(1,1) = cosrot + onemcos * omega(1) ** 2
      this % rot_matrix(1,2) = sinrot * omega(3) +                             &
     &                  onemcos * omega(1) * omega(2)
      this % rot_matrix(1,3) = -sinrot * omega(2) +                            &
     &                  onemcos * omega(1) * omega(3)

      this % rot_matrix(2,1) = -sinrot * omega(3) +                            &
     &                  onemcos * omega(1) * omega(2)
      this % rot_matrix(2,2) = cosrot + onemcos * omega(2) ** 2
      this % rot_matrix(2,3) = sinrot * omega(1) +                             &
     &                  onemcos * omega(2) * omega(3)

      this % rot_matrix(3,1) = sinrot * omega(2) +                             & 
     &                  onemcos * omega(1) * omega(3)
      this % rot_matrix(3,2) = -sinrot * omega(1) +                            &
     &                  onemcos * omega(2) * omega(3)
      this % rot_matrix(3,3) = cosrot + onemcos * omega(3) ** 2


!  Assignment of the c_of_rot vector, if present gives the present value of the 
!  c_of_rot or gives the default value zero

      IF (PRESENT(c_of_rot)) THEN 
         this % c_of_rot = c_of_rot
      ELSE
         this % c_of_rot = zero
      END IF 
                       
!  Assignment of the shift vector, if not present gives the default value zero
      IF (PRESENT(shift)) THEN 
         this % shift = shift
      ELSE
         this % shift = zero        
      END IF 

      END SUBROUTINE bsc_construct_rs

!*******************************************************************************
! SECTION V. DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Destroy a coil
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_destroy_coil(this)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (bsc_coil), INTENT(inout) :: this

!  Start of executable code

!  Get rid of all components
      this % c_type = ''
      this % s_name = ''
      this % l_name = ''
      this % current = zero
      this % eps_sq = zero
      this % raux = zero
      this % rcirc = zero
      this % xcent = zero
      this % enhat = zero
      this % ave_n_area = zero
      
      IF (ASSOCIATED(this % xnod)) DEALLOCATE(this % xnod)      
      IF (ASSOCIATED(this % dxnod)) DEALLOCATE(this % dxnod)      
      IF (ASSOCIATED(this % ehnod)) DEALLOCATE(this % ehnod)      
      IF (ASSOCIATED(this % lsqnod)) DEALLOCATE(this % lsqnod)
      IF (ASSOCIATED(this % lnod)) DEALLOCATE(this % lnod)

      END SUBROUTINE bsc_destroy_coil

!-------------------------------------------------------------------------------
!  Destroy an array of coils
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_destroy_coil_a(this)

!  Declare Arguments 
      TYPE (bsc_coil), DIMENSION(:), INTENT(inout) :: this
      INTEGER :: n, nsize
      
      nsize = SIZE(this)
      DO n = 1, nsize
         CALL bsc_destroy_coil(this(n))
      END DO
      
      END SUBROUTINE bsc_destroy_coil_a
      
!-------------------------------------------------------------------------------
!  Destroy a coilcoll
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_destroy_coilcoll(this)

!  Declare Arguments 
      TYPE (bsc_coilcoll), INTENT(inout) :: this

!  Declare local variables
      INTEGER(iprec) :: ncoild

!  Start of executable code

      this % ncoil = 0
      this % s_name = ''
      this % l_name = ''

!  Get rid of coils. Destroy them one by one, to avoid memory leaks.
      IF (ASSOCIATED(this % coils)) THEN
         ncoild = SIZE(this % coils)
         CALL bsc_destroy(this % coils(1:ncoild))         
         DEALLOCATE(this % coils)
      END IF

      END SUBROUTINE bsc_destroy_coilcoll

!*******************************************************************************
! SECTION VI. ASSIGNMENT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Assignment for coils
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_coil_to_coil(left,right)
      IMPLICIT NONE
      
!  Declare Arguments 
      TYPE (bsc_coil), INTENT (out) :: left
      TYPE (bsc_coil), INTENT (in) :: right
      
!  Declare temporary variables
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: temp2
      REAL(rprec), DIMENSION(:), ALLOCATABLE:: temp1

!  Declare local variables
      INTEGER(iprec) :: n, nm1
         
!  Start of executable code

!      WRITE(*,*) ' Executing bsc_coil_to_coil'
      
!  Non-pointer variables.
!  Copy them over, for all types of coils
      left % c_type = right % c_type
      left % s_name = right % s_name
      left % l_name = right % l_name
      left % current = right % current
      left % eps_sq = right % eps_sq
      left % raux = right % raux
      left % rcirc = right % rcirc
      left % ave_n_area = right % ave_n_area
      left % xcent = right % xcent 
      left % enhat = right % enhat

!  Pointers components of bsc_coil.
!  Only bother with this coding if coil is a fil_loop or fil_rogo
      IF ((right % c_type .eq. 'fil_loop') .or.                                &
     &   (right % c_type .eq. 'floop') .or.                                    &
     &   (right % c_type .eq. 'fil_rogo')) THEN
       
!  Find the SIZE of the pointer arrays.
         n   = SIZE(right % xnod,2)
         nm1 = n - 1

!  Allocate space for the left components, and copy stuff from the right.

!  If left and right are identical, then naive coding could
!  accidentally deallocate right before it can get copied.
!  To avoid this, I first copy the right components into temporary space
!  and then deal with the left component deallocation and allocation.

!  Two dimensional pointers
         ALLOCATE(temp2(3,n))

         temp2(1:3,1:n) =  right % xnod(1:3,1:n) 
         IF (ASSOCIATED(left % xnod)) DEALLOCATE(left % xnod)
         ALLOCATE(left % xnod(3,n))
         left % xnod(1:3,1:n) = temp2(1:3,1:n)

         temp2(1:3,1:nm1) =  right % dxnod(1:3,1:nm1)                
         IF (ASSOCIATED(left % dxnod)) DEALLOCATE(left % dxnod)
         ALLOCATE(left % dxnod(3,nm1))                               
         left % dxnod(1:3,1:nm1) = temp2(1:3,1:nm1)                  

         temp2(1:3,1:nm1) =  right % ehnod(1:3,1:nm1)                
         IF (ASSOCIATED(left % ehnod)) DEALLOCATE(left % ehnod)
         ALLOCATE(left % ehnod(3,nm1))                               
         left % ehnod(1:3,1:nm1) = temp2(1:3,1:nm1)                  

         DEALLOCATE(temp2)
      
!  One dimensional pointers
         ALLOCATE(temp1(n))

         temp1(1:nm1) =  right % lsqnod(1:nm1)                          
         IF (ASSOCIATED(left % lsqnod)) DEALLOCATE(left % lsqnod)
         ALLOCATE(left % lsqnod(nm1))                                   
         left % lsqnod(1:nm1) = temp1(1:nm1)                            

         temp1(1:nm1) =  right % lnod(1:nm1) 
         IF (ASSOCIATED(left % lnod)) DEALLOCATE(left % lnod)
         ALLOCATE(left % lnod(nm1))                                     
         left % lnod(1:nm1) = temp1(1:nm1)                              

         DEALLOCATE(temp1)

      END IF ! c_type .eq. fil_loop or fil_rogo
         
      END SUBROUTINE bsc_coil_to_coil

!-------------------------------------------------------------------------------
!  Assignment for arrays of type bsc_coil
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_coil_a_to_coil_a(left,right)
      
!  Declare Arguments 
      TYPE (bsc_coil), DIMENSION(:), INTENT (out) :: left
      TYPE (bsc_coil), DIMENSION(:), INTENT (in) :: right
      
!  Declare temporary variables
      INTEGER(iprec) :: nleft, nright, i
         
!  Start of executable code

!      WRITE(*,*) ' Executing bsc_coil_a_to_coil_a'

      nleft = SIZE(left)
      nright = SIZE(right)
      IF (nleft .ne. nright) THEN
         WRITE(*,*) 'FATAL in bsc_coil_a_to_coil_a. nleft .ne. nright'
         STOP
      END IF
      
!  Assignment, one by one
      DO i = 1,nleft
         left(i) = right(i)
      END DO

      END SUBROUTINE bsc_coil_a_to_coil_a


!  ZZZ What about scalar = array(1)?
!  ZZZ What about broadcast, array(i:j) = scalar?

!*******************************************************************************
! SECTION VII. COIL COLLECTION MANIPULATION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Append a coil to a coil collection
!
!    Note that this has some clunky coding, (within the IF test)
!      so that the number of coils in a coil collection can be
!      arbitrarily large.
!  ZZZ Think about ways to do this better.
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_append(this, newcoil)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (bsc_coilcoll), INTENT(inout) :: this
      TYPE (bsc_coil), INTENT(in) :: newcoil

!  Declare local variables
      INTEGER(iprec) :: icoil, ncoild
      INTEGER(iprec), PARAMETER :: nincr = 10

      TYPE (bsc_coil), DIMENSION(:), ALLOCATABLE :: coils_temp

!  Start of executable code

!  Check to see status of coilcoll      
      IF (.NOT. ASSOCIATED(this % coils)) THEN
         CALL bsc_construct(this,'id from bsc_append', '')
      END IF
      ncoild = SIZE(this % coils)

!  Check to see if need to increment size of coils array
      IF (this % ncoil + 1 .gt. ncoild) THEN
      
!  Make some temporary space, and copy coils to it
         ALLOCATE(coils_temp(ncoild))

!  JDH 03.27.03. SPH had coding in here for different compilers
!  Problem was that one compiler did not like array syntax for
!         coils_temp(1:ncoild) = this % coils(1:ncoild)
!   "ON IBM RISC, NEED V7.1.1. OR LATER FOR ASSIGNMENT OF ARRAY TO WORK CORRECTLY"
!  I just eliminated array syntax coding. Do loop works, and avoids conditional
!  compilation.
         DO icoil = 1, ncoild
            coils_temp(icoil) = this % coils(icoil)
         END DO

!  Destroy the existing coils and free up the space
         CALL bsc_destroy(this % coils(1:ncoild))         
         DEALLOCATE(this % coils)
         
!  Increase the size of the coils array, and copy old stuff there
         ALLOCATE(this % coils(ncoild + nincr))
         this % coils(1:ncoild) = coils_temp(1:ncoild)

!  Get rid of the temporary coils
         CALL bsc_destroy(coils_temp(1:ncoild))
         DEALLOCATE(coils_temp)

      END IF
!  End of Check to see if need to increment SIZE of coils
      
!  Copy newcoil onto end of ncoils
      this % ncoil = this % ncoil + 1
      this % coils(this % ncoil) = newcoil

      END SUBROUTINE bsc_append
!*******************************************************************************
! SECTION VIII. ROTATION AND SHIFT SUBROUTINES
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Generic Subroutine: bsc_rot_shift
!    First argument is the object to be rotated and shifted
!    Second argument is a bsc_rs, which carries the rotation and shift information
!
!  Different versions depending on first argument:
!     bsc_rot_shift_pts       First argument is a two dimensional array of points
!     bsc_rot_shift_pt        First argument is a single point (Dimension(3))
!     bsc_rot_shift_coil      First argument is a bsc_coil 
!     bsc_rot_shift_coil_a    First argument is an array of bsc_coil
!     bsc_rot_shift_coilcoll  First argument is a bsc_coilcoll
!
!  Particular subroutines, called by __coil, depending on c_type. 
!     bsc_rot_shift_coil_fil_loop
!     bsc_rot_shift_coil_fil_circ
!     c_type = 'fil_rogo': also calls bsc_rot_shift_coil_fil_loop
!-------------------------------------------------------------------------------
!  Rotation and Shift for an array of points
!    Optional third argument xyz_dim is an integer variable that specifies
!    which index of the 2 d array is the coordinate index. Legal values are 1 and 2.
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_rot_shift_pts(this,my_rs,xyz_dim)
      IMPLICIT NONE

!  Argument Declaration
!  Required Arguments
      REAL(rprec), DIMENSION(:,:), INTENT(inout) :: this
      TYPE (bsc_rs), INTENT(in) :: my_rs
      
!  Optional Arguments
      INTEGER(iprec), OPTIONAL, INTENT(in) :: xyz_dim

!  Local Variable Declaration
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: this_temp 
      REAL(rprec), DIMENSION(3) :: shift_2
      INTEGER(iprec) :: xyz_dim_use
      INTEGER(iprec) :: is1, js2, i, j, nm1
            
!  Start of executable code

      is1 = SIZE(this, 1)
      js2 = SIZE(this, 2)

!  Allocates the temporary array of rotated vectors with the same size as this
      ALLOCATE(this_temp(is1,js2))
      
!  Logic to determinate the value of xyz_dim_use 

      IF (PRESENT(xyz_dim)) THEN
         SELECT CASE (xyz_dim)
         CASE (1)
            IF (is1 == 3) THEN
               xyz_dim_use = 1
            ELSE
               WRITE(*,*) 'ERROR: bsc_rot_shift_pts: xyz_dim invalid'         
               STOP
            END IF   
         
         CASE (2)
            IF (js2 == 3) THEN
               xyz_dim_use = 2
            ELSE
               WRITE(*,*) 'ERROR: bsc_rot_shift_pts: xyz_dim invalid'         
               STOP
            END IF   
         
         CASE DEFAULT 
            WRITE(*,*) 'WARNING: bsc_rot_shift_pts: xyz_dim has not a',        &
     &                  'valid value'         
         END SELECT
      
      ELSE
         IF (is1 == 3) THEN
            xyz_dim_use = 1
         ELSE IF (js2 == 3) THEN
            xyz_dim_use = 2
         ELSE
            WRITE(*,*) 'ERROR: bsc_rot_shift_pts: points have no ',            &
     &                  'dimension of length 3'         
            STOP
         END IF      
      
      END IF

!  Computes the shifts and the rotations
      shift_2(1:3) = my_rs % c_of_rot(1:3) + my_rs % shift(1:3)

      SELECT CASE (xyz_dim_use)
      CASE (1)
         this(1:3,1:js2) = this(1:3,1:js2) -                                   &
     &                     SPREAD(my_rs % c_of_rot,2,js2)
         this_temp = MATMUL(my_rs % rot_matrix, this)
         this(1:3,1:js2) = this_temp(1:3,1:js2) + SPREAD(shift_2,2,js2)

      CASE (2)
         this(1:is1,1:3) = this(1:is1,1:3) -                                   &
     &                     SPREAD(my_rs % c_of_rot,1,is1)     
         this_temp = MATMUL(this,TRANSPOSE(my_rs % rot_matrix))
         this(1:is1,1:3) = this_temp(1:is1,1:3) + SPREAD(shift_2,1,is1)

      CASE DEFAULT
         WRITE(*,*) 'FATAL ERROR: bsc_rot_shift_pts: xyz_dim_use is',          &
     &                 ' not a valid value (1 or 2)'         

      END SELECT
         
!  Deallocates the temporary array of rotated vectors
      DEALLOCATE(this_temp)

      END SUBROUTINE bsc_rot_shift_pts

!-------------------------------------------------------------------------------
!  Rotation and Shift for single point
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_rot_shift_pt(this,my_rs)
      IMPLICIT NONE

!  Argument Declaration
!  Required Arguments
      REAL(rprec), DIMENSION(3), INTENT(inout) :: this
      TYPE (bsc_rs), INTENT(in) :: my_rs
      
!  Local Variable Declaration
      REAL(rprec), DIMENSION(3) :: this_temp 
            
!  Start of executable code
      this = this - my_rs % c_of_rot
      this_temp = MATMUL(my_rs % rot_matrix, this)
      this = this_temp + my_rs % c_of_rot + my_rs % shift 

      END SUBROUTINE bsc_rot_shift_pt

!-------------------------------------------------------------------------------
!  Rotation and Shift for single coil
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_rot_shift_coil(this,my_rs)
      IMPLICIT NONE

!  Argument Declaration
!  Required Arguments
      TYPE (bsc_coil), INTENT(inout) :: this
      TYPE (bsc_rs), INTENT(in) :: my_rs

!  Local Variable Declaration
            
!  Start of executable code

      SELECT CASE (this % c_type)
      CASE ('fil_loop','floop','fil_rogo')
         CALL bsc_rot_shift_coil_fil_loop(this,my_rs)

      CASE ('fil_circ','fcirc')
         CALL bsc_rot_shift_coil_fil_circ(this,my_rs)

      CASE DEFAULT 
         WRITE(*,*) 'FATAL: bsc_rot_shift_coil: c_type unrecognized:',         &
     &     this % c_type
         STOP
      END SELECT

      END SUBROUTINE bsc_rot_shift_coil

!-------------------------------------------------------------------------------
!  Rotation and Shift for single filamentary loop
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_rot_shift_coil_fil_loop(this,my_rs)
!  Rotation and shift for single filamentary loop
!  Should only be called from bsc_rot_shift_coil

!  Argument Declaration
!  Required Arguments
      TYPE (bsc_coil), INTENT(inout) :: this
      TYPE (bsc_rs), INTENT(in) :: my_rs

!  Local Variable Declaration
      INTEGER(iprec) :: nwire, nm1, i
            
!  Start of executable code

      nwire = SIZE(this % xnod,2)
      nm1 = MAX(1, nwire-1)

!  Rotates the points that form the coil

      CALL bsc_rot_shift_pts(this % xnod,my_rs,xyz_dim = 1_iprec)

!  Recompute dxnod, ehat vectors for rotated coils
!  lnod and lsqnod should be invariant (check this?)

      this % dxnod(1:3,1:nm1) = this % xnod(1:3,2:nwire)                       &
     &                        - this % xnod(1:3,1:nm1)

      DO i = 1,3
         this % ehnod(i,1:nm1) = this % dxnod(i,1:nm1) /                       &                     
     &         this % lnod(1:nm1)                                     
      END DO
      
      END SUBROUTINE bsc_rot_shift_coil_fil_loop

!-------------------------------------------------------------------------------
!  Rotation and shift for single filamentary circle
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_rot_shift_coil_fil_circ(this,my_rs)
      IMPLICIT NONE
!  Rotation and shift for single filamentary circle
!  Should only be called from bsc_rot_shift_coil

!  Argument Declaration
!  Required Arguments
      TYPE (bsc_coil), INTENT(inout) :: this
      TYPE (bsc_rs), INTENT(in) :: my_rs

!  Local Variable Declaration
      REAL(rprec), DIMENSION(3,2) :: rot_vectors

!  Start of executable code

!  Rotates the circular coil and its normal vector with respect to the
!  center of rotation

      rot_vectors(1:3,1) = this % xcent(1:3)
      rot_vectors(1:3,2) = this % enhat(1:3) + this % xcent(1:3)

      CALL bsc_rot_shift_pts(rot_vectors,my_rs,xyz_dim = 1_iprec)
      
      this % xcent(1:3) = rot_vectors(1:3,1)
      this % enhat(1:3) = rot_vectors(1:3,2) - this % xcent(1:3)

      END SUBROUTINE bsc_rot_shift_coil_fil_circ

!-------------------------------------------------------------------------------
!  Rotation and Shift for array of coils
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_rot_shift_coil_a(this,my_rs)
      IMPLICIT NONE

!  Argument Declaration
!  Required Arguments
      TYPE (bsc_coil), DIMENSION(:), INTENT(inout) :: this
      TYPE (bsc_rs), INTENT(in) :: my_rs
      
!  Local Variable Declaration
      INTEGER(iprec) :: n, nsize
            
!  Start of executable code
      nsize = SIZE(this)

      DO n = 1,nsize
         CALL bsc_rot_shift_coil(this(n),my_rs)
      END DO

      END SUBROUTINE bsc_rot_shift_coil_a

!-------------------------------------------------------------------------------
!  Rotation and Shift for a coil collection
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_rot_shift_coilcoll(this,my_rs)

!  Argument Declaration
!  Required Arguments
      TYPE (bsc_coilcoll), INTENT(inout) :: this
      TYPE (bsc_rs), INTENT(in) :: my_rs
      
!  Start of executable code
      CALL bsc_rot_shift_coil_a(this % coils(1:this % ncoil),my_rs)

      END SUBROUTINE bsc_rot_shift_coilcoll

!*******************************************************************************
! SECTION IX. VECTOR POTENTIAL SUBROUTINES
!*******************************************************************************
!  Generic Subroutine: bsc_a
!  Different versions depending on first argument:
!     bsc_a_coil       First argument a bsc_coil. All calls come through here.
!     bsc_a_coil_a     First argument an array of bsc_coil
!     bsc_a_coilcoll   First argument a bsc_coilcoll
!  Particular to a c_type. Called from bsc_a_coil
!     bsc_a_coil_fil_loop
!     bsc_a_coil_fil_circ
!  NB c_type = 'fil_rogo' calls B field: bsc_b_coil_fil_loop
!-------------------------------------------------------------------------------
!  Vector Potential field for single coil
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_a_coil(this,x,a,bsc_k2)
      IMPLICIT NONE

!  Argument Declaration
!  Required Arguments
      TYPE (bsc_coil), INTENT(in) :: this
      REAL(rprec), DIMENSION(3), INTENT(in) :: x
      REAL(rprec), DIMENSION(3), INTENT(out) :: a
      
!  Optional Arguments
      REAL(rprec), OPTIONAL, INTENT(in) :: bsc_k2
      
!  Local Variable Declaration
      REAL(rprec) :: bsc_k2_use
            
!  Start of executable code
      SELECT CASE (this % c_type)
      CASE ('fil_loop','floop')
         CALL bsc_a_coil_fil_loop(this,x,a)

      CASE ('fil_circ','fcirc')
         CALL bsc_a_coil_fil_circ(this,x,a)
      
      CASE ('fil_rogo')
!  Rogowski. Compute the B field from the nodes.
         CALL bsc_b_coil_fil_loop(this,x,a)
!  Scale by the correct factor.
         a(1:3) = this % ave_n_area * a(1:3)

      CASE DEFAULT 
         WRITE(*,*) 'FATAL: bsc_a_coil: c_type unrecognized:',                 &
     &     this % c_type
         STOP
      END SELECT

! Scale with current and k2 constant
      IF (PRESENT(bsc_k2)) THEN
         bsc_k2_use = bsc_k2
      ELSE
         bsc_k2_use = bsc_k2_def
      ENDIF
      
      a(1:3) = this % current * bsc_k2_use * a(1:3)
      
      END SUBROUTINE bsc_a_coil

!-------------------------------------------------------------------------------
!  Vector Potential field for single filamentary loop
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_a_coil_fil_loop(this,x,a)
!  Vector Potential field for single filamentary loop
!  Should only be called from bsc_a_coil

!  Argument Declaration
!  Required Arguments
      TYPE (bsc_coil), INTENT(IN) :: this
      REAL(rprec), DIMENSION(3), INTENT(in) :: x
      REAL(rprec), DIMENSION(3), INTENT(out) :: a
      
!  Local Variable Declaration
      INTEGER(iprec) :: n, nm1, i
      REAL(rprec), DIMENSION(SIZE(this % xnod,2))   :: capR
      REAL(rprec), DIMENSION(SIZE(this % xnod,2)-1) :: lnfactor
            
!  Start of executable code
      n = SIZE(this % xnod,2)
      nm1 = n - 1

!  Form array of relative vector lengths
      capR = SQRT((x(1) - this % xnod(1,1:n))**2                               &
     &     +      (x(2) - this % xnod(2,1:n))**2                               &
     &     +      (x(3) - this % xnod(3,1:n))**2)

!  Form lnfactor
      lnfactor(1:nm1) = this % lnod(1:nm1)/(capR(1:nm1) + capR(2:n))
      CALL log_eps (lnfactor)

!      lnfactor(1:nm1) = log((capR(1:nm1) + capR(2:n)
!     &   + this % lnod(1:nm1)) / (capR(1:nm1) + capR(2:n)
!     &   - this % lnod(1:nm1)))

! Sum for A field
      DO i = 1,3
         a(i) = DOT_PRODUCT(this % ehnod(i,1:nm1),lnfactor(1:nm1))
      END DO
      
      END SUBROUTINE bsc_a_coil_fil_loop

!-------------------------------------------------------------------------------
!  Vector Potential field for single filamentary circle
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_a_coil_fil_circ(this,x,a)
      IMPLICIT NONE
!  Vector Potential field for single filamentary circle
!  Should only be called from bsc_a_coil

!  Argument Declaration
!  Required Arguments
      TYPE (bsc_coil), INTENT(in) :: this
      REAL(rprec), DIMENSION(3), INTENT(in) :: x
      REAL(rprec), DIMENSION(3), INTENT(out) :: a
      
!  Local Variable Declaration
      REAL(rprec), PARAMETER :: two_third = 2._rprec / 3
      REAL(rprec), PARAMETER :: pio16 = pi / 16
      REAL(rprec), DIMENSION(3) :: rprime, rhoprime, rhophat, phiphat
      REAL(rprec) :: zprime, fsq, em, emone
      REAL(rprec) :: rhopmsq, rhopmag, f, rf, rd, brackg, aphi, radcc

      REAL(rprec), PARAMETER :: c_gb1 = one
      REAL(rprec), PARAMETER :: c_gb2 = 3._rprec / 4
      REAL(rprec), PARAMETER :: c_gb3 = 75._rprec / 128
      REAL(rprec), PARAMETER :: c_gb4 = 245._rprec / 512
      REAL(rprec), PARAMETER :: c_gb5 = 6615._rprec / 16384
      REAL(rprec), PARAMETER :: c_gb6 = 22869._rprec / 65536
      REAL(rprec), PARAMETER :: c_gb7 = 1288287._rprec / 4194304

            
!  Start of executable code
!  Transform to local (primed) coordinates
      rprime = x(1:3) - this % xcent(1:3)
      zprime = DOT_PRODUCT(rprime,this % enhat)
      rhoprime = rprime(1:3) - zprime * this % enhat(1:3)
      rhopmsq = DOT_PRODUCT(rhoprime,rhoprime)
!  ZZZ Think about numbers 1.e-30 and 1.e-15
      IF (rhopmsq .lt. 1.e-30_rprec) THEN
         rhoprime(1) = 1.e-15_rprec
         rhoprime(2) = zero
         rhoprime(3) = zero
         rhopmsq = DOT_PRODUCT(rhoprime,rhoprime)
      END IF
      rhopmag = SQRT(rhopmsq)
      rhophat(1:3) = rhoprime(1:3) / rhopmag
      radcc = this % rcirc

!  various factors. 
      fsq = one / ((rhopmag + radcc) ** 2 + zprime ** 2)
      f = SQRT(fsq)
      em = 4 * rhopmag * radcc * fsq
      emone = one - em
      
      IF (em .gt. bsc_emcut) THEN           ! large m, use elliptic integrals
!  All the calculation of the complete elliptic integrals is localized
!  in bsc_cei
         CALL bsc_cei(emone,rf,rd)
         brackg = two_third * rd - rf
      ELSE                                  ! small m, use power series in m
         brackg = pio16 * em * (c_gb1 + em * (c_gb2 + em *                     &
     &      (c_gb3  + em * (c_gb4 + em * (c_gb5 + em * (c_gb6 + em *           &
     &      c_gb7))))))
      END IF
      aphi = 4 * radcc * f * brackg

!  Now convert to global cylindrical coordinates.
!  First, find the phi_prime_hat by taking the cross
!  product of z_prime_hat and rho_prime_hat
      phiphat(1) = this % enhat(2) * rhophat(3) -                              &
     &   this % enhat(3) * rhophat(2)
      phiphat(2) = this % enhat(3) * rhophat(1) -                              &
     &   this % enhat(1) * rhophat(3)
      phiphat(3) = this % enhat(1) * rhophat(2) -                              &
     &   this % enhat(2) * rhophat(1)
      a(1:3) = phiphat(1:3) * aphi

      END SUBROUTINE bsc_a_coil_fil_circ

!-------------------------------------------------------------------------------
!  Vector Potential field for array of coils
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_a_coil_a(this,x,a,bsc_k2)
      IMPLICIT NONE

!  Argument Declaration
!  Required Arguments
      TYPE (bsc_coil), DIMENSION(:), INTENT(IN) :: this
      REAL(rprec), DIMENSION(3), INTENT(IN) :: x
      REAL(rprec), DIMENSION(3), INTENT(OUT) :: a
      
!  Optional Arguments
      REAL(rprec), OPTIONAL, INTENT(in) :: bsc_k2
      
!  Local Variable Declaration
      INTEGER(iprec) :: i, n
      REAL(rprec), DIMENSION(3,SIZE(this)) :: aarray
            
!  Start of executable code

!      WRITE(*,*) 'Executing bsc_a_coil_a'

!  calls to bsc_a_coil
      n = SIZE(this)
      DO i = 1,n
         CALL bsc_a_coil(this(i),x,aarray(1:3,i))
      END DO

! Sum for A field
      a(1:3) = SUM(aarray(1:3,1:n),2)
      
!  Rescale if bsc_k2 present
      IF (PRESENT(bsc_k2)) THEN
         a(1:3) = a(1:3) * bsc_k2 * bsc_k2inv_def
      ENDIF
      
      RETURN
      END SUBROUTINE bsc_a_coil_a

!-------------------------------------------------------------------------------
!  Vector Potential for a coil collection
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_a_coilcoll(this,x,a,bsc_k2)

!  Argument Declaration
!  Required Arguments
      TYPE (bsc_coilcoll), INTENT(in) :: this
      REAL(rprec), DIMENSION(3), INTENT(in) :: x
      REAL(rprec), DIMENSION(3), INTENT(out) :: a
      
!  Optional Arguments
      REAL(rprec), OPTIONAL, INTENT(in) :: bsc_k2
      
!  Local Variable Declaration
      INTEGER(iprec) :: n
            
!  Start of executable code
!  Calculate a from coils in coilcoll
      n = this % ncoil
      IF (n .gt. 0) then
         CALL bsc_a_coil_a(this % coils(1:n),x,a(1:3))
      ELSE
         a(1:3) = zero
      END IF

!  Rescale if bsc_k2 present
      IF (PRESENT(bsc_k2)) THEN
         a(1:3) = a(1:3) * bsc_k2 * bsc_k2inv_def
      ENDIF
      
      END SUBROUTINE bsc_a_coilcoll

!*******************************************************************************
! SECTION X. MAGNETIC FIELD SUBROUTINES
!*******************************************************************************
!  Generic Subroutine: bsc_b
!  Different versions depending on first argument:
!     bsc_b_coil       First argument a bsc_coil. All calls come through here.
!     bsc_b_coil_a     First argument an array of bsc_coil
!     bsc_b_coilcoll   First argument a bsc_coilcoll
!  Particular to a c_type. Called from bsc_b_coil
!     bsc_b_coil_fil_loop
!     bsc_b_coil_fil_circ
!  NB: c_type = 'fil_rogo' doesn't compute magnetic field.
!-------------------------------------------------------------------------------
!  Magnetic field for single coil
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_b_coil (this,x,b,bsc_k2)
      IMPLICIT NONE

!  Argument Declaration
!  Required Arguments
      TYPE (bsc_coil), INTENT(in) :: this
      REAL(rprec), DIMENSION(3), INTENT(in) :: x
      REAL(rprec), DIMENSION(3), INTENT(out) :: b
      
!  Optional Arguments
      REAL(rprec), OPTIONAL, INTENT(in) :: bsc_k2
      
!  Local Variable Declaration
      REAL(rprec) :: bsc_k2_use


!  Start of executable code
      SELECT CASE (this % c_type)
      CASE ('fil_loop','floop')
         CALL bsc_b_coil_fil_loop(this,x,b)

      CASE ('fil_circ','fcirc')
         CALL bsc_b_coil_fil_circ(this,x,b)
      
      CASE ('fil_rogo')
!  Rogowski. Not yet implemented
         WRITE(*,*) 'WARNING: bsc_b_coil: NOT YET IMPLEMENTED',                &
     &     this % c_type

      CASE DEFAULT 
         WRITE(*,*) 'FATAL: bsc_b_coil: c_type unrecognized:',                 &
     &     this % c_type
         STOP
      END SELECT

! Scale with current and k2 constant
      IF (PRESENT(bsc_k2)) THEN
         bsc_k2_use = bsc_k2
      ELSE
         bsc_k2_use = bsc_k2_def
      ENDIF
      
      b(1:3) = this % current * bsc_k2_use * b(1:3)
      
      END SUBROUTINE bsc_b_coil

!-------------------------------------------------------------------------------
!  Magnetic field for single filamentary loop
!-------------------------------------------------------------------------------
!2 4 6 8(1)2 4 6 8(2)2 4 6 8(3)2 4 6 8(4)2 4 6 8(5)2 4 6 8(6)2 4 6 8(7)2 4 6 8(8
      SUBROUTINE bsc_b_coil_fil_loop (this,x,b)
      IMPLICIT NONE
!  Should only be called from bsc_b_coil
!  Can also be called from bsc_a_coil with c_type = 'fil_rogo'

!  Argument Declaration
!  Required Arguments
      TYPE (bsc_coil), INTENT(in) :: this
      REAL(rprec), DIMENSION(3), INTENT(in) :: x
      REAL(rprec), DIMENSION(3), INTENT(out) :: b
      
!  Local Variable Declaration
      INTEGER(iprec) :: nm1, n, i, j, k
      REAL(rprec), DIMENSION(SIZE(this % xnod,2))     :: capR
      REAL(rprec), DIMENSION(3,SIZE(this % xnod,2))   :: capRv  
      REAL(rprec), DIMENSION(SIZE(this % xnod,2)-1)   :: Rfactor, R1p2
      REAL(rprec), DIMENSION(3,SIZE(this % xnod,2)-1) :: crossv 
            
!  Start of executable code
      n = SIZE(this % xnod,2)
      nm1 = n - 1

!  Form array of vectors relative to observation point x(i)
      DO i = 1,3                                                              
         capRv(i,1:n) = x(i) - this % xnod(i,1:n)                             
      END DO                                                                  
      
!  Form array of relative vector lengths
!JDH Quick Fix 2007-05-24
      capR(1:n) = SQRT(MAX(this % eps_sq,capRv(1,1:n) * capRv(1,1:n) +         &
     &                 capRv(2,1:n) * capRv(2,1:n) +                           &
     &                 capRv(3,1:n) * capRv(3,1:n)))                                          
!      capR(1:n) = SQRT(capRv(1,1:n) * capRv(1,1:n) +                           &
!     &                 capRv(2,1:n) * capRv(2,1:n) +                           &
!     &                 capRv(3,1:n) * capRv(3,1:n))                                          

!  Form Cross Product
      DO i = 1, 3
         j = mod(i,3_iprec) + 1
         k = mod(j,3_iprec) + 1
         crossv(i,1:nm1) = this % dxnod(j,1:nm1) * capRv(k,1:nm1)              &
     &                   - this % dxnod(k,1:nm1) * capRv(j,1:nm1)                                
      END DO

      R1p2(1:nm1) = capR(1:nm1) + capR(2:n)
      Rfactor(1:nm1) = 2 * R1p2(1:nm1) / (capR(1:nm1) * capR(2:n) *            &
     &      MAX(R1p2(1:nm1)*R1p2(1:nm1) - this % lsqnod(1:nm1),                &
     &          this % eps_sq))

! Sum for B field
      DO i = 1,3                                                              
         b(i) = DOT_PRODUCT(crossv(i,1:nm1),Rfactor(1:nm1))
      END DO
            
      END SUBROUTINE bsc_b_coil_fil_loop

!-------------------------------------------------------------------------------
! Magnetic Field for single filamentary circle
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_b_coil_fil_circ(this,x,b)
      IMPLICIT NONE
!  Magnetic field for single filamentary circle
!  Should only be CALLed from bsc_b_coil

!  Argument Declaration
!  Required Arguments
      TYPE (bsc_coil), INTENT(in) :: this
      REAL(rprec), DIMENSION(3), INTENT(in) :: x
      REAL(rprec), DIMENSION(3), INTENT(out) :: b
      
!  Local Variable Declaration
      REAL(rprec), PARAMETER :: third = one / 3
      REAL(rprec), PARAMETER :: sixth = one / 6
      REAL(rprec), PARAMETER :: two_third = 2._rprec / 3

      REAL(rprec), DIMENSION(3) :: rprime, rhoprime, rhophat, phiphat
      REAL(rprec) :: zprime, fsq, fcube, em, emone, cfcube, geofac1
      REAL(rprec) :: rhopmsq, rhopmag, f, rf, rd, brackg, brackh,              &
     &   brfac1, brho, bz, aphi, radcc, add_on, fac1, fac2
     
!  Coefficients for power series in m
      REAL(rprec) :: c_bz0, c_bz1, c_bz2, c_bz3, c_bz4, c_bz5, c_bz6
      REAL(rprec) :: cb_bz0, cb_bz1, cb_bz2, cb_bz3, cb_bz4, cb_bz5,           &
     &   cb_bz6
      REAL(rprec), PARAMETER :: c_brho1 = one
      REAL(rprec), PARAMETER :: c_brho2 = 5._rprec / 4
      REAL(rprec), PARAMETER :: c_brho3 = 175._rprec / 128
      REAL(rprec), PARAMETER :: c_brho4 = 735._rprec / 512
      REAL(rprec), PARAMETER :: c_brho5 = 24255._rprec / 16384
      REAL(rprec), PARAMETER :: c_brho6 = 99099._rprec / 65536
      REAL(rprec), PARAMETER :: c_brho7 = 6441435._rprec / 4194304
      REAL(rprec), PARAMETER :: ca_bz1 =  3._rprec / 4
      REAL(rprec), PARAMETER :: ca_bz2 =  75._rprec / 128
      REAL(rprec), PARAMETER :: ca_bz3 =  245._rprec / 512
      REAL(rprec), PARAMETER :: ca_bz4 =  6615._rprec / 16384
      REAL(rprec), PARAMETER :: ca_bz5 =  22869._rprec / 65536
      REAL(rprec), PARAMETER :: ca_bz6 =  1288287._rprec / 4194304
      REAL(rprec), PARAMETER :: pio4 =  pi / 4
      REAL(rprec), PARAMETER :: pi3o16 =  3._rprec * pi / 16
            
!  Start of executable code
!  Transform to local (primed) coordinates
      rprime = x(1:3) - this % xcent(1:3)
      zprime = DOT_PRODUCT(rprime,this % enhat)
      rhoprime = rprime(1:3) - zprime * this % enhat(1:3)
      rhopmsq = DOT_PRODUCT(rhoprime,rhoprime)
!  ZZZ Think about numbers 1.e-30 and 1.e-15
      IF (rhopmsq .lt. 1.e-30_rprec) THEN
         rhoprime(1) = 1.e-15_rprec
         rhoprime(2) = zero
         rhoprime(3) = zero
         rhopmsq = DOT_PRODUCT(rhoprime,rhoprime)
      END IF
      rhopmag = SQRT(rhopmsq)
      rhophat(1:3) = rhoprime(1:3) / rhopmag
      radcc = this % rcirc

!  various factors. 
      fsq = one / ((rhopmag + radcc) ** 2 + zprime ** 2)
      f = SQRT(fsq)
      em = 4 * rhopmag * radcc * fsq
      emone = one - em
      cfcube = 4 * radcc * f * fsq 
!  Current gets multiplied in bsc_b_coil
      geofac1 = (radcc ** 2 + zprime ** 2) / rhopmag
      
      IF (em .gt. bsc_emcut) THEN                ! large m, use elliptic integrals
!  All the calculation of the complete elliptic integrals is localized
!  in bsc_cei
         CALL bsc_cei(emone,rf,rd)
         brackg = two_third * rd - rf
         brackh = sixth * (-(1 + 3 * emone) * rd +                             &
     &      3 * (1 + emone) * rf) / emone
         brfac1 = brackg + 2 * brackh
         brho = cfcube * zprime * brfac1
         bz = cfcube * (brackg * (geofac1 + radcc)                             &
     &      + brackh * (geofac1 - rhopmag))
      ELSE                                          ! small m, use power series
         fac2 = pi3o16 * zprime * cfcube
         brho = fac2 * em * (c_brho1 + em * (c_brho2 +                         &
     &      em * (c_brho3 + em * (c_brho4 + em * (c_brho5 + em *               &
     &      (c_brho6 + em * c_brho7))))))
         fac1 = cfcube * radcc * fsq * pio4
         add_on = zprime ** 2 + (radcc + rhopmag) * (radcc - rhopmag)
         cb_bz0 = (radcc + rhopmag) * rhopmag + 2 * add_on
         cb_bz1 = cb_bz0 + add_on
         cb_bz2 = cb_bz1 + add_on
         cb_bz3 = cb_bz2 + add_on
         cb_bz4 = cb_bz3 + add_on
         cb_bz5 = cb_bz4 + add_on
         cb_bz6 = cb_bz5 + add_on
         c_bz0 = cb_bz0
         c_bz1 = ca_bz1 * cb_bz1
         c_bz2 = ca_bz2 * cb_bz2
         c_bz3 = ca_bz3 * cb_bz3
         c_bz4 = ca_bz4 * cb_bz4
         c_bz5 = ca_bz5 * cb_bz5
         c_bz6 = ca_bz6 * cb_bz6
         bz = fac1 * (c_bz0 + em * (c_bz1 + em * (c_bz2 +                      &
     &      em * (c_bz3 + em * (c_bz4 + em * (c_bz5 + em * c_bz6))))))
      END IF

!  Now convert to global cylindrical coordinates.
      b(1:3) = brho * rhophat(1:3) + bz * this % enhat(1:3)

      END SUBROUTINE bsc_b_coil_fil_circ

!-------------------------------------------------------------------------------
!  Magnetic field for array of coils
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_b_coil_a(this,x,b,bsc_k2)
      IMPLICIT NONE

!  Argument Declaration
!  Required Arguments
      TYPE (bsc_coil), DIMENSION(:), INTENT(IN) :: this
      REAL(rprec), DIMENSION(3), INTENT(IN) :: x
      REAL(rprec), DIMENSION(3), INTENT(OUT) :: b
      
!  Optional Arguments
      REAL(rprec), OPTIONAL, INTENT(in) :: bsc_k2
      
!  Local Variable Declaration
      INTEGER(iprec) :: i, n
      REAL(rprec), DIMENSION(3,SIZE(this)) :: barray
            
!  Start of executable code
!  Calls to bsc_b_coil
      n = SIZE(this)
      DO i = 1,n
         CALL bsc_b_coil(this(i),x,barray(1:3,i))
      END DO

! Sum for B field
      b = SUM(barray(1:3,1:n),2)
      
!  Rescale if bsc_k2 present
      IF (PRESENT(bsc_k2)) THEN
         b(1:3) = b(1:3) * bsc_k2 * bsc_k2inv_def
      ENDIF
      
      END SUBROUTINE bsc_b_coil_a

!-------------------------------------------------------------------------------
!  Magnetic field for a coil collection
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_b_coilcoll(this,x,b,bsc_k2)
      IMPLICIT NONE

!  Argument Declaration
!  Required Arguments
      TYPE (bsc_coilcoll), INTENT(in) :: this
      REAL(rprec), DIMENSION(3), INTENT(in) :: x
      REAL(rprec), DIMENSION(3), INTENT(out) :: b
      
!  Optional Arguments
      REAL(rprec), OPTIONAL, INTENT(in) :: bsc_k2
      
!  Local Variable Declaration
      INTEGER(iprec) :: n
            
!  Start of executable code
!  Calculate b from coils
      n = this % ncoil
      IF (n .gt. 0) then
         CALL bsc_b_coil_a(this % coils(1:n),x,b(1:3))
      ELSE
         b(1:3) = zero
      END IF

!  Rescale if bsc_k2 present
      IF (PRESENT(bsc_k2)) THEN
         b(1:3) = b(1:3) * bsc_k2 * bsc_k2inv_def
      ENDIF
      
      END SUBROUTINE bsc_b_coilcoll

!*******************************************************************************
! SECTION XI.  MAGNETIC FLUX INTEGRAL SUBROUTINES
!*******************************************************************************
!  No mutual induction subroutine yet - use magnetic flux instead
!  Mutual inductance can be computed from the flux, by dividing by the current
!  in coil_a, and multiplying by the appropriate number of turns for coil_a and
!  coil_b
!  Note: bsc_fluxba is the generic subroutine. There are implementations for
!  the first argument as a bsc_coil, as an array of bsc_coils, and as a bsc_coilcoll
!  The three subroutine should look very similar, since the difference is only
!  in the vector potential.
!  All of the work related to coil_b is put into subroutine bsc_flux_pos, so that 
!  it does not have to be duplicated in each of the three subroutines.
!  2012-01-17 JDH. Add len_integrate argument

!-------------------------------------------------------------------------------
!  Magnetic Flux for a single coil
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_fluxba_coil(coil_a,coil_b,len_integerate,flux,            &
     &   bsc_k2)
      IMPLICIT NONE

!  This subroutine calculates the magnetic flux due to coil_a through coil_b.
!  The flux is computed using the vector potential due to coil_a, with a line
!  integral around coil_b. 

!  Argument Declaration

!  Required Arguments
      TYPE (bsc_coil), INTENT(in) :: coil_a
      TYPE (bsc_coil), INTENT(in) :: coil_b
      REAL(rprec), INTENT(IN)     :: len_integerate
      REAL(rprec), INTENT(out)    :: flux
      
!  Optional Arguments
      REAL(rprec), OPTIONAL, INTENT(in) :: bsc_k2
      
!  Local Variable Declaration
      REAL(rprec), DIMENSION(:,:), POINTER :: positions => null(),             &
     &  tangents => null(), avecs => null()
      INTEGER(iprec) :: i, npoints
            
!  Start of executable code

!  Get array of positions at which to evaluate the vector potential
!  Subroutine also allocates space to the pointers: positions, tangents, and
!  avecs.
      CALL bsc_flux_pos(coil_b,len_integerate,positions,tangents,              &
     &   avecs,npoints)

!  Calculate vector potentials
      SELECT CASE (coil_b % c_type)
      CASE DEFAULT
         DO i = 1,npoints
            CALL bsc_a(coil_a,positions(1:3,i),avecs(1:3,i))
         END DO
      CASE ('fil_rogo')
         DO i = 1,npoints
            CALL bsc_b(coil_a,positions(1:3,i),avecs(1:3,i))
         END DO
      END SELECT

!  Do summations
      CALL bsc_flux_sum(coil_b,positions,tangents,avecs,npoints,flux)

!  Rescale if bsc_k2 present
      IF (PRESENT(bsc_k2)) THEN
         flux = flux * bsc_k2 * bsc_k2inv_def
      ENDIF

!  Deallocate space. 
      DEALLOCATE(avecs,positions,tangents)

!  That's all
      END SUBROUTINE bsc_fluxba_coil

!-------------------------------------------------------------------------------
!  Magnetic Flux for an array of coils
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_fluxba_coil_a(coil_a,coil_b,len_integerate,flux,          &
     &   bsc_k2)
      IMPLICIT NONE

!  This subroutine calculates the magnetic flux due to coil_a through coil_b.
!  The flux is computed using the vector potential due to coil_a, with a line
!  integral around coil_b. 

!  Argument Declaration

!  Required Arguments
      TYPE (bsc_coil), DIMENSION(:), INTENT(IN) :: coil_a
      TYPE (bsc_coil), INTENT(in) :: coil_b
      REAL(rprec), INTENT(IN)     :: len_integerate
      REAL(rprec), INTENT(out) :: flux
      
!  Optional Arguments
      REAL(rprec), OPTIONAL, INTENT(in) :: bsc_k2
      
!  Local Variable Declaration
      REAL(rprec), DIMENSION(:,:), POINTER :: positions => null(),             &
     &  tangents => null(), avecs => null()
      INTEGER(iprec) :: i, npoints
            
!  Start of executable code

!  Get array of positions at which to evaluate the vector potential
!  Subroutine also allocates space to the pointers: positions, tangents, and
!  avecs.
      CALL bsc_flux_pos(coil_b,len_integerate,positions,tangents,              &
     &   avecs,npoints)

!  Calculate vector potentials
      SELECT CASE (coil_b % c_type)
      CASE DEFAULT
         DO i = 1,npoints
            CALL bsc_a(coil_a,positions(1:3,i),avecs(1:3,i))
         END DO
      CASE ('fil_rogo')
         DO i = 1,npoints
            CALL bsc_b(coil_a,positions(1:3,i),avecs(1:3,i))
         END DO
      END SELECT

!  Do summations
      CALL bsc_flux_sum(coil_b,positions,tangents,avecs,npoints,flux)

!  Rescale if bsc_k2 present
      IF (PRESENT(bsc_k2)) THEN
         flux = flux * bsc_k2 * bsc_k2inv_def
      ENDIF

!  Deallocate space. 
      DEALLOCATE(avecs,positions,tangents)

!  That's all
      END SUBROUTINE bsc_fluxba_coil_a
 
!-------------------------------------------------------------------------------
!  Magnetic Flux for a coil collection
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_fluxba_coilcoll(coil_a,coil_b,len_integerate,flux,        &
     &   bsc_k2)
      IMPLICIT NONE

!  This subroutine calculates the magnetic flux due to coil_a through coil_b.
!  The flux is computed using the vector potential due to coil_a, with a line
!  integral around coil_b. 

!  Argument Declaration

!  Required Arguments
      TYPE (bsc_coilcoll), INTENT(IN) :: coil_a
      TYPE (bsc_coil), INTENT(in) :: coil_b
      REAL(rprec), INTENT(IN)     :: len_integerate
      REAL(rprec), INTENT(out) :: flux
      
!  Optional Arguments
      REAL(rprec), OPTIONAL, INTENT(in) :: bsc_k2
      
!  Local Variable Declaration
      REAL(rprec), DIMENSION(:,:), POINTER :: positions => null(),             &
     &  tangents => null(), avecs => null()
      INTEGER(iprec) :: i, npoints
            
!  Start of executable code

!  Get array of positions at which to evaluate the vector potential
!  Subroutine also allocates space to the pointers: positions, tangents, and
!  avecs.
      CALL bsc_flux_pos(coil_b,len_integerate,positions,tangents,              &
     &   avecs,npoints)

!  Calculate vector potentials
      SELECT CASE (coil_b % c_type)
      CASE DEFAULT
         DO i = 1,npoints
            CALL bsc_a(coil_a,positions(1:3,i),avecs(1:3,i))
         END DO
      CASE ('fil_rogo')
         DO i = 1,npoints
            CALL bsc_b(coil_a,positions(1:3,i),avecs(1:3,i))
         END DO
      END SELECT

!  Do summations
      CALL bsc_flux_sum(coil_b,positions,tangents,avecs,npoints,flux)

!  Rescale if bsc_k2 present
      IF (PRESENT(bsc_k2)) THEN
         flux = flux * bsc_k2 * bsc_k2inv_def
      ENDIF

!  Deallocate space. 
      DEALLOCATE(avecs,positions,tangents)

!  That's all
      END SUBROUTINE bsc_fluxba_coilcoll
        
        
!-------------------------------------------------------------------------------
!  Subroutine to find positions and tangents for second coil
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_flux_pos(coil_b,len_integrate,positions,                  &
     &   tangents,avecs,npoints)
      IMPLICIT NONE
!  This subroutine will return positions and tangents to the coil coil_b.
!  It is called by the bsc_flux_ subroutines.
!  For now, the positions chosen lie ON the coil. In future, this should be
!  changed, so that the evaluations points are nearby (about a wire radius away)
!  2012-01-17 JDH. Add len_integrate argument.

!     Argument Declaration
!  coil_b         Coil to determine integration positions for
!  len_integrate  Integration length. If zero or negative, use one point
!                 per segment.
!  positions      Array of positions for integration
!  tangents       Array of tangent vectors

!  Required Arguments
      TYPE (bsc_coil), INTENT(in) :: coil_b
      REAL(rprec), INTENT(in) :: len_integrate
      REAL(rprec), DIMENSION(:,:), POINTER :: positions, tangents,             &
     &   avecs
      INTEGER(iprec), INTENT(out) ::  npoints
      
!  Local Variable Declaration
      INTEGER(iprec) :: i, j, iseg, npoints_this_segment, n_denom
      INTEGER(iprec) :: nseg
      REAL(rprec), DIMENSION(3) :: xhatlocal, yhatlocal, zhatlocal,            &
     &  rhohat, phihat
      REAL(rprec) :: dphi, phi, cphi, sphi, frac_denom
      INTEGER(iprec), PARAMETER :: npcirc = 60   ! zzz - Picked out of a hat

!  Start of executable code
      
!  Define length of arrays
      SELECT CASE (coil_b % c_type)
      CASE ('fil_loop','floop','fil_rogo')
         IF (len_integrate .le. zero) THEN
            npoints = SIZE(coil_b % lnod)
         ELSE
            npoints = 0
            DO iseg = 1,SIZE(coil_b % lnod)
               npoints_this_segment = coil_b % lnod(iseg) /                    &
     &             len_integrate + 1
               npoints = npoints + npoints_this_segment
            END DO
         ENDIF
      CASE ('fil_circ','fcirc')
         IF (len_integrate .le. zero) THEN
            npoints = npcirc
         ELSE
            npoints = twopi * coil_b % rcirc / len_integrate + 1
            npoints = MAX(npcirc,npoints)
         ENDIF
      END SELECT

!  Allocate space
      IF (ASSOCIATED(positions)) DEALLOCATE(positions)
      IF (ASSOCIATED(tangents)) DEALLOCATE(tangents)
      IF (ASSOCIATED(avecs)) DEALLOCATE(avecs)
      ALLOCATE(positions(3,npoints),tangents(3,npoints),                       &
     &   avecs(3,npoints))

!  For fil_loops, points spaced throughout the segment
!  First and last points in a segment weighted at half the weight of all other
!  segments.
!    points per segment       position fraction              weights
!          1                        1/2                        1
!          2                    1/4, 3/4                      1/2, 1/2
!          3                    1/8, 1/2, 7/8               1/4, 1/2, 1/4
!          4                    1/12, 4/12, 8/12, 11/12     1/6, 1/3, 1/3, 1/6
      SELECT CASE (coil_b % c_type)
      CASE ('fil_loop','floop','fil_rogo')
         i = 1
         nseg = SIZE(coil_b % lnod)
         DO iseg = 1,nseg
            IF (len_integrate .le. zero) THEN
               npoints_this_segment = 1
            ELSE
               npoints_this_segment = coil_b % lnod(iseg) /                    &
     &             len_integrate + 1
            ENDIF
            IF (npoints_this_segment .eq. 1) THEN
               positions(1:3,i) = coil_b % xnod(1:3,iseg) +                    &
     &            0.5 * coil_b % dxnod(1:3,iseg)
               tangents(1:3,i) = coil_b % dxnod(1:3,iseg)
               i = i + 1
            ELSEIF (npoints_this_segment .eq. 2) THEN
               positions(1:3,i) = coil_b % xnod(1:3,iseg) +                    &
     &             0.25 * coil_b % dxnod(1:3,iseg)
               tangents(1:3,i) = 0.5 * coil_b % dxnod(1:3,iseg)
               i = i + 1
               positions(1:3,i) = coil_b % xnod(1:3,iseg) +                    &
     &            0.75 * coil_b % dxnod(1:3,iseg)
               tangents(1:3,i) = 0.5 * coil_b % dxnod(1:3,iseg)
               i = i + 1
            ELSE                            ! Here, 3 or more points per segment
               n_denom = npoints_this_segment - 1
               frac_denom = one / n_denom
               positions(1:3,i) = coil_b % xnod(1:3,iseg) +                    &
     &            0.25 * frac_denom * coil_b % dxnod(1:3,iseg)
               tangents(1:3,i) = 0.5 * frac_denom *                            &
     &            coil_b % dxnod(1:3,iseg)
               i = i + 1               
               DO j = 2,npoints_this_segment - 1
                  positions(1:3,i) = coil_b % xnod(1:3,iseg) +                 &
     &               frac_denom * (j-1) * coil_b % dxnod(1:3,iseg)   
                  tangents(1:3,i) = frac_denom *                               &
     &                coil_b % dxnod(1:3,iseg)
                  i = i + 1
               END DO
               positions(1:3,i) = coil_b % xnod(1:3,iseg) +                    &
     &            (1 - 0.25 * frac_denom) * coil_b % dxnod(1:3,iseg)
               tangents(1:3,i) = 0.5 * frac_denom *                            &
     &            coil_b % dxnod(1:3,iseg)
               i = i + 1
            ENDIF
         END DO
         IF (i - 1 .ne. npoints) THEN
            STOP 'ERROR 1 in bsc_flux_pos'
         ENDIF
      CASE ('fil_circ','fcirc')
!  For circles, npcirc points on the circle
         CALL bsc_triplet(coil_b % enhat,xhatlocal,yhatlocal,zhatlocal)
         dphi = twopi / npoints
         DO i = 1,npoints
            phi = i * dphi
            cphi = cos(phi)
            sphi = sin(phi)
            rhohat = cphi * xhatlocal + sphi * yhatlocal
            phihat = - sphi * xhatlocal + cphi * yhatlocal
            positions(1:3,i) = coil_b % xcent(1:3) + coil_b % rcirc *          &
     &         rhohat
            tangents(1:3,i) = coil_b % rcirc * dphi * phihat
         END DO
      END SELECT
      
      END SUBROUTINE bsc_flux_pos      

!-------------------------------------------------------------------------------
!  Subroutine to find do summations for flux
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_flux_sum(coil_b,positions,tangents,avecs,npoints,         &
     &   flux)
      IMPLICIT NONE
!  This subroutine will compute some sums
!  It is called by the bsc_flux_ subroutines

!  Argument Declaration

!  Required Arguments
      TYPE (bsc_coil), INTENT(in) :: coil_b
      REAL(rprec), DIMENSION(:,:), POINTER :: positions, tangents,             &
     &   avecs
      INTEGER(iprec), INTENT(in) ::  npoints
      REAL(rprec), INTENT(out) :: flux
      
!  Local Variable Declaration
      INTEGER(iprec) :: i

!  Start of executable code
      
!  Sum dot products
      flux = zero
      DO i = 1,npoints
         flux = flux + DOT_PRODUCT(avecs(1:3,i),tangents(1:3,i))
      END DO

!  Extra Factor for coil_b a Rogowski coil
      IF (coil_b % c_type .eq. 'fil_rogo') THEN
         flux = flux * coil_b % ave_n_area
      END IF
      
      END SUBROUTINE bsc_flux_sum      

!*******************************************************************************
! SECTION XII.  AUXILIARY SUBROUTINES
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Subroutine to compute complete elliptic integrals needed for circular loops
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_cei(xarg, rf, rd)
      IMPLICIT NONE
!  
!  Should only be called from bsc_a_coil_fil_circ subroutine
!  or from the bsc_b_coil_fil_circ subroutine

!  This subroutine calculates Rf(0,xarg,1) and Rd(0,xarg,1)
!  where Rf(x,y,z) and Rd(x,y,z) are the Carlson forms of the
!  elliptic integral functions. The algorithm is based on those in section
!  6.11 of Numerical Recipes, 2nd edition. The iteration loops have been
!  combined. The first two iterations have been done explicitly, since
!  the Rf and Rd arguments x and z are 0 and 1, respectively.

!  Argument Declaration
!  Required Arguments
      REAL(rprec), INTENT(IN) :: xarg
      REAL(rprec), INTENT(out) :: rf, rd
      
!  Local Variable Declaration
      REAL(rprec) :: xt, yt, zt, sum, x
      REAL(rprec) :: alamb, fac, sqrtx, sqrty, sqrtz, avef, recavef,           &
     &     delxf, delyf, delzf, e2f, e3f, aved, recaved, delxd,                &
     &     delyd, delzd, ead, ebd, ecd, edd, eed
      INTEGER(iprec) :: iter
      INTEGER(iprec), PARAMETER :: niter = 5

!   parameters for the algorithm
      REAL(rprec), PARAMETER :: third = one / 3._rprec
      REAL(rprec), PARAMETER :: c1f = one / 24._rprec
      REAL(rprec), PARAMETER :: c2f = 0.1_rprec
      REAL(rprec), PARAMETER :: c3f = 3._rprec / 44._rprec
      REAL(rprec), PARAMETER :: c4f = one / 14._rprec
      REAL(rprec), PARAMETER :: sixth = one / 6._rprec
      REAL(rprec), PARAMETER :: twelfth = one / 12._rprec
      REAL(rprec), PARAMETER :: c1d = 3._rprec / 14._rprec
      REAL(rprec), PARAMETER :: c2d = one / 6._rprec
      REAL(rprec), PARAMETER :: c3d = 9._rprec / 22._rprec
      REAL(rprec), PARAMETER :: c4d = 3._rprec / 26._rprec
      REAL(rprec), PARAMETER :: c5d = .25_rprec * c3d
      REAL(rprec), PARAMETER :: c6d = 1.5_rprec * c4d
            
!  Start of executable code
!  First, make sure that 0 < x <= 1
      x = min(max(xarg,1.e-12_rprec),one)

!  Do first iteration explicitly, since xt, yt, and zt are known
!  (Saves two square roots)
      alamb = SQRT(x)
      fac = .25_rprec
      xt = fac * alamb
      yt = fac * (x + alamb)
      zt = fac * (one + alamb)
      sum = one / (one + alamb)

!  Do other iterations
      do iter = 2,niter
         sqrtx = SQRT(xt)
         sqrty = SQRT(yt)
         sqrtz = SQRT(zt)
         alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz
         sum = sum + fac / (sqrtz * (zt + alamb))
         xt = .25_rprec * (xt + alamb)
         yt = .25_rprec * (yt + alamb)
         zt = .25_rprec * (zt + alamb)
         fac = .25_rprec * fac
      end do

!  Done with iterations. Now get stuff for power series evaluation
!  of Rf and Rd.

      aved = .2_rprec * (xt + yt + 3._rprec * zt)
      recaved = one / aved
      delxd = one - xt * recaved
      delyd = one - yt * recaved
      delzd = one - zt * recaved
      ead = delxd * delyd
      ebd = delzd * delzd
      ecd = ead - ebd
      edd = ead - 6._rprec * ebd
      eed = edd + ecd + ecd
      rd = 3._rprec * sum + fac * (one + edd * (- c1d + c5d * edd -            &
     &   c6d * delzd * eed) + delzd * (c2d * eed + delzd *                     &
     &   (- c3d * ecd + delzd * c4d * ead))) * recaved *                       &
     &   SQRT(recaved)

      avef = third * (xt + yt + zt)
      recavef = one / avef
      delxf = one - xt * recavef
      delyf = one - yt * recavef
      delzf = one - zt * recavef
      e2f = delxf * delyf - delzf ** 2
      e3f = delxf * delyf * delzf
      rf = (one + (c1f * e2f - c2f - c3f * e3f) *                              &
     &   e2f + c4f * e3f) * SQRT(recavef)

      END SUBROUTINE bsc_cei

!-------------------------------------------------------------------------------
!  Subroutine to compute local unit vectors, given an initial z direction
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_triplet(zlocal,xhatlocal,yhatlocal,zhatlocal)
      IMPLICIT NONE
!  Subroutine to compute a right-handed set of orthogonal unit vectors (x,y,z)hatlocal
!  where zhatlocal is parallel to the input argument zlocal.

!  Argument Declaration

!  Required Arguments
      REAL(rprec), DIMENSION(3), INTENT(in) :: zlocal
      REAL(rprec), DIMENSION(3), INTENT(out) :: xhatlocal, yhatlocal,          &
     &   zhatlocal
      
!  Local Variable Declaration
      REAL(rprec), DIMENSION(3) :: zloc, xloc
      REAL(rprec) :: fac, lsq

!  Start of executable code
!  Normalize zlocal. Copy first, so zlocal can be INTENT(in).
      zloc = zlocal
      lsq = DOT_PRODUCT(zloc,zloc)
      IF (lsq .le. 0.) THEN
         lsq = 1._rprec
         zloc(3) = 1._rprec
      END IF
      fac = 1._rprec / SQRT(lsq)
      zhatlocal = fac * zloc
      
!  Find a direction that is not parallel to zhatlocal
      IF (abs(zhatlocal(1)) .le. 0.8_rprec) THEN
         xloc = (/ 1._rprec, 0._rprec, 0._rprec /)
      ELSE 
         xloc = (/ 0._rprec, 1._rprec, 0._rprec /)
      END IF

!  Subtract off piece of xloc that is parallel to zhatlocal
      fac = DOT_PRODUCT (xloc,zhatlocal)
      xloc = xloc - fac * zhatlocal

!  Normalize xloc, call it xhatlocal 
      lsq = DOT_PRODUCT(xloc,xloc)     
      fac = 1._rprec / SQRT(lsq)
      xhatlocal = fac * xloc

!  Cross product to determine yhatlocal
      yhatlocal(1) = zhatlocal(2) * xhatlocal(3) -                             &
     &   zhatlocal(3) * xhatlocal(2)
      yhatlocal(2) = zhatlocal(3) * xhatlocal(1) -                             &
     &   zhatlocal(1) * xhatlocal(3)
      yhatlocal(3) = zhatlocal(1) * xhatlocal(2) -                             &
     &   zhatlocal(2) * xhatlocal(1)

      END SUBROUTINE bsc_triplet
      
!-------------------------------------------------------------------------------
!  Subroutine to compute the mean position of a coil  
!    this            :   bsc_coil to compute its mean position 
!    mean_r          :   real(rprec), array (size 3) output specifying the
!                        mean position of the coil
!
!  Difference between mean_xnod and mean_r is in the weighting.
!    mean_xnod is just a simple average of the node positions
!    mean_r is length weighted.
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_mean_r(this,mean_r)
      IMPLICIT NONE

!  Argument Declaration
!  Required Arguments
      TYPE (bsc_coil), INTENT(in) :: this 
      REAL(rprec), DIMENSION(3), INTENT(out) :: mean_r
      
!  Local Variable Declaration
      INTEGER :: nwire, nm1, iwire, i
      REAL(rprec) :: coil_length

!  Start of executable code

!  Computes the mean position depending on the type

      SELECT CASE (this % c_type)
      CASE ('fil_circ','fcirc')
         mean_r(1:3) = this % xcent(1:3)   
      CASE ('fil_loop','floop','fil_rogo')
         nwire = SIZE(this % xnod,2)
         nm1 = MAX(1, nwire-1)
         coil_length =  SUM(this % lnod(1:nm1))
         DO i = 1,3
            mean_r(i) = DOT_PRODUCT(this % lnod(1:nm1),                        &
     &         this % xnod(i,1:nm1) + 0.5 * this % dxnod(i,1:nm1)) /           & 
     &         coil_length
         END DO
      CASE DEFAULT 
         WRITE(*,*) 'FATAL: bsc_mean_r:                                        & 
     &   c_type unrecognized:', this % c_type
         STOP
      END SELECT

      END SUBROUTINE bsc_mean_r
      
!-------------------------------------------------------------------------------
!  Subroutine to compute the mean xnod of a coil  
!    this            :   bsc_coil to compute its mean xnod 
!    mean_xnod       :   real(rprec), array (size 3) output specifying the
!                        mean position of the coil. Average of Node positions.
!
!  Difference between mean_xnod and mean_r is in the weighting.
!    mean_xnod is just a simple average of the node positions
!    mean_r is length weighted.
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_mean_xnod(this,mean_xnod)
      IMPLICIT NONE

!  Argument Declaration
!  Required Arguments
      TYPE (bsc_coil), INTENT(in) :: this 
      REAL(rprec), DIMENSION(3), INTENT(out) :: mean_xnod
      
!  Local Variable Declaration
      INTEGER :: nwire, nm1, iwire, i

!  Start of executable code

!  Computes the mean position depending on the type

      SELECT CASE (this % c_type)
      CASE ('fil_circ','fcirc')
         mean_xnod(1:3) = this % xcent(1:3)   
      CASE ('fil_loop','floop','fil_rogo')
         nwire = SIZE(this % xnod,2)
         nm1 = MAX(1, nwire-1)
         DO i = 1,3
            mean_xnod(i) = SUM(this % xnod(i,1:nm1)) / nm1
         END DO
      CASE DEFAULT 
         WRITE(*,*) 'FATAL: bsc_mean_r:                                        & 
     &   c_type unrecognized:', this % c_type
         STOP
      END SELECT

      END SUBROUTINE bsc_mean_xnod
 
          
!*******************************************************************************
! SECTION XIII.  DEBUGGING SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Print out the contents of a coil
!-------------------------------------------------------------------------------

      SUBROUTINE bsc_spill_coil(this,identifier)
      IMPLICIT NONE
!  Subroutine to print out the contents of a coil structure

      TYPE (bsc_coil), INTENT (in) :: this
      CHARACTER (len=*) :: identifier
      INTEGER(iprec) :: i, n

!  start of executable code
      WRITE(*,*)
      WRITE(*,*)    'sub spill_coil CALL with -- ',identifier,' --'

! Components common to all c_types
      WRITE(*,*)    '  c_type = ', this % c_type
      WRITE(*,*)    '  s_name = ', this % s_name
      WRITE(*,*)    '  l_name = ', this % l_name
      WRITE(*,*)    '  eps_sq = ', this % eps_sq
      WRITE(*,*)    '  current = ', this % current
      WRITE(*,*)    '  raux = ', this % raux
      IF (this % c_type .eq. 'fil_rogo') THEN
         WRITE(*,*)    '  ave_n_area = ', this % ave_n_area
      END IF

      SELECT CASE (this % c_type)
      CASE ('fil_loop','floop','fil_rogo')
         IF (ASSOCIATED(this % xnod)) THEN
            n = SIZE(this % xnod,2)
            WRITE(*,*)    '  i xnod(1,i), xnod(2,i), xnod(3,i)'
            DO i = 1,n
               WRITE(*,*) i, this % xnod(1:3,i)
            END DO
!  Comment out Printout of dxnod, lsqnod, lnod and ehnod - redundant info
!            WRITE(*,*)    '  i dxnod(1-3,i), lsqnod(i)'
!            DO i = 1,n-1
!               WRITE(*,*) i, this % dxnod(1:3,i), this % lsqnod(i)
!            END DO
!            WRITE(*,*)    '  i ehnod(1-3,i), lnod(i)'
!            DO i = 1,n-1
!               WRITE(*,*) i, this % ehnod(1:3,i), this % lnod(i)
!            END DO
         ELSE
            WRITE(*,*)   '  xnod is NOT ASSOCIATED'
         END IF ! this % xnod ASSOCIATED
      CASE ('fil_circ','fcirc')
         WRITE(*,*)    '  rcirc = ', this % rcirc
         WRITE(*,*)    '  xcent = ', this % xcent(1:3)
         WRITE(*,*)    '  enhat = ', this % enhat(1:3)
      END SELECT

      END SUBROUTINE bsc_spill_coil

!-------------------------------------------------------------------------------
!  Print out lots of stuff about a coil collection
!-------------------------------------------------------------------------------
      SUBROUTINE bsc_spill_coilcoll(this,identifier)
      IMPLICIT NONE
!  Subroutine to print out the contents of a coil structure

      TYPE (bsc_coilcoll), INTENT (in) :: this
      CHARACTER (len=*) :: identifier
      INTEGER(iprec) :: i, n, ncoild
      

!  start of executable code
      WRITE(*,*)
      WRITE(*,*)    'sub spill_coilcoll CALL with -- ',identifier,' --'
      WRITE(*,*)    '  s_name = ', this % s_name
      WRITE(*,*)    '  l_name = ', this % l_name
      WRITE(*,*)    '  ncoil = ', this % ncoil
      ncoild = SIZE(this % coils)
      WRITE(*,*)    '  ncoild = ', ncoild
      
!  Note that loop below ru ns up through dimension of this % coils array
!  So that expect blank s_names for ncoil < i <= ncoild
      DO i = 1,ncoild
         WRITE(*,*) i,' this coil, s_name = ',this % coils(i) % s_name
         IF (ASSOCIATED(this % coils(i) % xnod)) THEN
            n = SIZE(this % coils(i) % xnod,2)
            WRITE(*,*)    ' xnod associated. SIZE(2) = ',n
         ELSE
            WRITE(*,*)   '  xnod is NOT ASSOCIATED'
         END IF
      END DO

!  JDH 12.03.02
!  Here's an alternate loop over the coils in the coilcoll:
!      DO i = 1,this % ncoil
!         CALL bsc_spill_coil(this % coils(i),'In bsc_spill_coilcoll')
!      END DO
      
      WRITE(*,*) 'End of spill_coilcoll'
      
      END SUBROUTINE bsc_spill_coilcoll


!*******************************************************************************
! SECTION XIV.  SPECIAL FUNCTION(S)
!*******************************************************************************
      SUBROUTINE log_eps (x1)
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(RPREC), INTENT(INOUT) :: x1(:)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: c3 = one/3,                                    &
     &   c5 = one/5, c7 = one/7, c9 = one/9, c11 = one/11
      REAL(rprec), PARAMETER :: threshold = 0.1_rprec
             !(MAX ALLOWED: 1/3, LOG(1+x), x < 1)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i
      REAL(rprec) :: x2
!-----------------------------------------------
!
!     FAST COMPUTATION OF LOG[(1+x)/(1-x)] USING TAYLOR OR PADE APPROXIMATION
!     FOR X < THRESHOLD
!
!     ON ENTRY,  x1 = E < 1  [ = LNOD/(RI + RF) ]
!     ON RETURN, x1 = LOG((1+E)/(1-E))
!
! JDH Notes, 11-21-2003:
! 1) Simplified coding, 11/2003. Previously, used for testing.
! 2) threshold value set from testing. Could vary with compiler, etc.
! 3) Other options tested (and coding eliminated) 
!       - Pade approximation instead of Taylor Series
!       - Where loop instead of do loop  
!
!-----------------------------------------------

      DO i = 1, SIZE(x1)
         IF (x1(i) <= threshold) THEN
            x2 = x1(i)*x1(i)

!     TAYLOR SERIES 
!     (TEST CASE: 11.2 s ON PC, MAX ERR = 1.305E-10 FOR THRESHOLD = 0.2,
!                                         1.6E-14   FOR THRESHOLD = 0.1)

             x1(i) = 2 * x1(i) * (one + x2 * (c3 + x2 * (c5 + x2 * (c7         &
     &                                + x2 * (c9 + x2 * c11)))))
         ELSE

!     THE MAX BELOW IS NECESSARY TO AVOID REAL SINGULARITY WHEN 
!     OBSERVATION POINT IS ON A COIL SEGMENT
!  bsc_mach_eps is machine epsilon parameter

            x1(i) = LOG((1 + x1(i))/MAX(1 - x1(i) ,bsc_mach_eps))

         END IF
      END DO

      RETURN

      END SUBROUTINE log_eps

!*******************************************************************************
! SECTION XV.  DUPLICATE CODING FOR TESTING
!*******************************************************************************
!  File 6.3: Duplicate magnetic field subroutine, and alter coding
!  so that can do speed tests.
!  The Interface block is above, in SECTION III
!  File 6.4 - eliminated.
!*******************************************************************************
! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!  
!  JDH 10.03.02
!     Adding power series in m to magnetic field and vector potential
!     subroutines, for circular coils.
!     Subroutines bsc_a_coil_fil_circ and bsc_b_coil_fil_circ
!  JDH 10.11.02
!     Fixed some bugs, extended power series in m to higher order. Adjusted
!     bsc_emcut appropriately.
!  JDH 10.23.02
!     A minor bug with fcirc fixed (in flux calculation), and a few comments changed.
!  JDH 11.19.2002
!     Added % raux component to bsc_coil definition. Used to carry auxiliary information
!     Not used within bsc.
!  JDH 11.20.2002
!     Changed file name from bsc_mod.f to bsc.f
!  JDH 12.03.02
!     1) Steve Hirshman made some minor changes, eliminating unused variables, and adding
!        lots of IMPLICIT NONE statements.
!     2) Updated bsc_spill_coil and bsc_spill_coilcoll
!  JDH 03.27.03
!    1) Added c_type='fil_rogo'
!    2) Changed many IF - ELSEIF structures to SELECT CASE structures.
!  JDH 05.15.03
!    1) Fixed bug with assignment of fil_rogo. (Forgot to copy '% ave_n_area')
!    2) In bsc_spill_coil, now print out ave_n_area.
!  JDH 05.16.03
!     Cleaned up intial variable declaration a bit
!     Revised the flux subroutines, so that they are shorter. More gets done in 
!     bsc_flux_pos and the new subroutine bsc_flux_sum
!  JDH 10.06.03
!     Fixed a bug in subroutine bsc_triplet. Problem was in the choice of a direction
!     for the initial xhatlocal.
!  JDH 11.22.03
!     Simplified coding in log_eps. Little changes throughout to make both free and fixed 
!     format compatable. Added bsc_mach_eps parameter, and replaced EPSILON calls.
!  JDH 07-14-03
!     Added Coil Rotations.
!     Original Coding by Steve Hirshman and Dennis Strickler.
!     Function calls and structure declarations changed from their original.
!     Parameterization of rotation is unchanged.
!     Also note the subroutines bsc_mean_r and bsc_mean_xnod, used to compute the mean
!     postion of a coil. These are in the Auxiliary Subroutines section.
!     Jorge Munoz did most of coding.
!     I made some further changes.
!
!  JDH 2007-06-13
!    Lines 1481-1487, in subroutine bsc_b_coil_fil_loop. Did a quick fix to 
!    eliminate a division by zero. (Change made 2007-05-04).

!  >>>  TO DO NEXT:
!    1) Fix up flux_ba, more points, with tuning parameters

!  JDH 2010-07-06
!    Changed coding for wrapping of coils. Need to eliminate zero length segments
!    without crashing, for W7X diagnostics.
!    2010-07-07. Passed test code from J. Hebert.
           
      END MODULE bsc_T
