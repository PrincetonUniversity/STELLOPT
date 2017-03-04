!*******************************************************************************
!  File recon_param_T.f
!  Contains module recon_param_T
!  Defines derived-types: recon_param, recon_cnstrnts
!   recon_param
!    This deals with the parameters that get changed during the course of
!    reconstruction. They are related to the variable parameters of the model.
!    Now (8/2006) the relation is simple - a 1-1 correpondence between
!    recon_param's and eq_param_var components. This may change in the future.
!   recon_cnstrnts
!     (New 2007-10-03)
!     This deals with constraints (relations) amongst the model parameters
!
!*******************************************************************************
!  MODULE recon_param_T
!    (Reconstruction Parameter Type Definition, for the V3FIT code)
! SECTION I.    VARIABLE DECLARATIONS
! SECTION II.   DERIVED-TYPE DECLARATIONS
! SECTION III.  INTERFACE BLOCKS
! SECTION IV.   CONSTRUCTION SUBROUTINES
! SECTION V.    DESTRUCTION SUBROUTINES
! SECTION VI.   CHANGE VALUE SUBROUTINES
! SECTION VII.  ASSIGNMENT SUBROUTINES
! SECTION VIII. OUTPUT SUBROUTINES

! SECTION XII.  AUXILIARY SUBROUTINES
! SECTION XIII. DEBUGGING SUBROUTINES
! SECTION XV.   DUPLICATE CODING FOR TESTING
! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
      MODULE recon_param_T

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Type declarations - lengths of reals, integers, and complexes.
!  Frequently used mathematical constants, lots of extra precision.
!-------------------------------------------------------------------------------
      USE stel_kinds, only : rprec, cprec
      USE stel_constants, only : pi, twopi, one, zero
     
!-------------------------------------------------------------------------------
!  Use Statements for other structures, V3 Utilities
!-------------------------------------------------------------------------------
      USE v3_utilities

!-------------------------------------------------------------------------------
!  Implicit None comes after USE statements, before other declarations
!-------------------------------------------------------------------------------
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Make type declarations and constants Private, so there are no conflicts.
!-------------------------------------------------------------------------------
      PRIVATE rprec, cprec, pi, twopi, one, zero

!-------------------------------------------------------------------------------
!  Lengths of Character Variables
!-------------------------------------------------------------------------------
      INTEGER, PARAMETER, PRIVATE :: type_len=20
      INTEGER, PARAMETER, PRIVATE :: sn_len=30
      INTEGER, PARAMETER, PRIVATE :: ln_len=80
      INTEGER, PARAMETER, PRIVATE :: units_len=30

!*******************************************************************************
! SECTION II. DERIVED-TYPE DECLARATIONS
! 1)   Reconstruction Parameter:
!         recon_param  
!     Type of reconstruction parameter specified by  % p_type.
!     Allowable values of p_type:
!       ac          - VMEC2000 current profile specification
!       ac_aux_s    _ VMEC2000 current profile specification, aux_s
!       ac_aux_f    _ VMEC2000 current profile specification, aux_f
!       ai          - VMEC2000 iota profile specification
!       ai_aux_s    _ VMEC2000 iota profile specification, aux_s
!       ai_aux_f    _ VMEC2000 iota profile specification, aux_f
!       am          - VMEC2000 pressure profile expansion coefficient
!       am_aux_s    _ VMEC2000 pressure profile specification, aux_s
!       am_aux_f    _ VMEC2000 pressure profile specification, aux_f
!       bloat       _ VMEC2000 - used to expand computation region
!       extcur      - VMEC2000 external coil current
!       curtor      - VMEC2000 torodial current
!       phiedge     - VMEC2000 toroidal flux
!       pres_scale  - VMEC2000 pressure scaling factor
!       density_max - max plasma density (assumed to be at magnetic axis)
!       density_tau - exponent on radial flux for density distribution of form:
!                   - ne = ne_max*(1-s^tau)^kappa +ne_ambient
! 2)   Reconstruction Constraint
!         recon_cnstrnts
!     Type of constraint is specified by % c_type
!     Allowable values of c_type:
!       ac_edge     - VMEC200 current profile (I') at edge specification
!       am_edge     - VMEC2000 pressure profile at edge specification
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Declare type recon_param
!       p_type          character, type of reconstruction parameter
!       index           integer, index to an array, if needed
!       value           value of the reconstruction parameter
!       vrnc            variance, or scaling value
!       fdf             finite difference factor - for determining the
!                       step size for finite difference partial derivatives
!                       (Added 2008-01-21)
!-------------------------------------------------------------------------------
      TYPE recon_param
         CHARACTER (len=type_len)         :: p_type
         INTEGER                          :: index
         REAL(rprec)                      :: value
         REAL(rprec)                      :: vrnc
         REAL(rprec)                      :: fdf
      END TYPE recon_param
!-------------------------------------------------------------------------------
!  Declare type recon_cnstrnts
!       c_type          character, type of reconstruction constraint
!       index           integer, index to an array, if needed
!       value           value of the reconstruction constraint
! Coding Notes
!   I want to put in some arrays here, so all the constraints are in one place
!   I don't want to bother with allocatable arrays, so I will just dimension
!   them of length 10, which should be plenty long enough.
!-------------------------------------------------------------------------------
      TYPE recon_cnstrnts
         INTEGER                                       :: n_cnstrnts
         CHARACTER(len=type_len), DIMENSION(10)        :: c_type
         INTEGER, DIMENSION(10)                        :: index
         REAL(rprec), DIMENSION(10)                    :: value
      END TYPE recon_cnstrnts

!*******************************************************************************
! SECTION III. INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Change Value Subroutine, with either scalar or vector
!-------------------------------------------------------------------------------
      INTERFACE recon_param_change_value
         MODULE PROCEDURE recon_param_change_value_s,                          & 
     &      recon_param_change_value_a
      END INTERFACE

      CONTAINS
!*******************************************************************************
! SECTION IV. CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a recon_param
!
!  TO DO - make index and vrnc optional arguments
!-------------------------------------------------------------------------------
      SUBROUTINE recon_param_construct(this,p_type,value,index,vrnc,           &
     &   fdf)

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (recon_param), INTENT(inout)          :: this
      CHARACTER (len=*), INTENT(in)              :: p_type
      REAL(rprec), INTENT(in)                    :: value
      INTEGER, INTENT(in)                        :: index
      REAL(rprec), INTENT(in)                    :: vrnc
      REAL(rprec), INTENT(in), OPTIONAL          :: fdf

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'recon_param_construct: '

!  Start of executable code

!  Assignments 
      this % p_type = TRIM(ADJUSTL(p_type))
      this % index = index
      this % value = value
      this % vrnc = vrnc
      IF(PRESENT(fdf)) THEN
         this % fdf = fdf
      ELSE
         this % fdf = one
      ENDIF

!  Different coding, depending on p_type
      SELECT CASE (TRIM(ADJUSTL(p_type)))
      CASE ('ac')
      CASE ('ac_aux_s')
      CASE ('ac_aux_f')
      CASE ('ai')
      CASE ('ai_aux_s')
      CASE ('ai_aux_f')
      CASE ('am')
      CASE ('am_aux_s')
      CASE ('am_aux_f')
      CASE ('bloat')
      CASE ('extcur')
      CASE ('curtor')
      CASE ('phiedge')
      CASE ('pres_scale')
      CASE ('density_max')
      CASE ('density_tau')

      CASE DEFAULT
         CALL err_fatal(sub_name // 'unrecognized p_type: ',                   &
     &      char=p_type)
      END SELECT ! Different coding depending on p_type
      
      END SUBROUTINE recon_param_construct
!-------------------------------------------------------------------------------
!  Construct a recon_cnstrnts
!
!  For c_type = 'ac_edge' (VMEC200 current profile (I') at edge specification)
!             = 'am_edge' (VMEC2000 pressure profile at edge specification)
!  Construct - make zero length. _append to get constraints into the type
!-------------------------------------------------------------------------------
      SUBROUTINE recon_cnstrnts_construct(this)

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (recon_cnstrnts), INTENT(inout)        :: this

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'recon_cnstrnt_construct: '

!  Start of executable code

!  Assignments
      this % n_cnstrnts = 0
      this % c_type = ''
      this % index = 0
      this % value = zero
      
      END SUBROUTINE recon_cnstrnts_construct
      
!*******************************************************************************
! SECTION IV.B APPEND SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Append to a recon_cnstrnts
!-------------------------------------------------------------------------------
!
!  For c_type = 'ac_edge' (VMEC200 current profile (I') at edge specification)
!             = 'am_edge' (VMEC2000 pressure profile at edge specification)
!-------------------------------------------------------------------------------
      SUBROUTINE recon_cnstrnts_append(this,c_type,value,index)

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (recon_cnstrnts), INTENT(inout)       :: this
      CHARACTER (len=*), INTENT(in)              :: c_type
      REAL(rprec), INTENT(in)                    :: value
      INTEGER, INTENT(in)                        :: index

!  Declare local variables
      INTEGER              :: n
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'recon_cnstrnts_append: '

!  Start of executable code

!  Check if too many
      this % n_cnstrnts = this % n_cnstrnts + 1
      n = this % n_cnstrnts
      IF (n .gt. SIZE(this % c_type)) THEN
         CALL err_warn('Too many recon-cnstrnts' // sub_name,int=n)
      ENDIF

!  Different coding, depending on c_type
      SELECT CASE (TRIM(ADJUSTL(c_type)))
      CASE ('ac_edge')
      CASE ('am_edge')

      CASE DEFAULT
         CALL err_fatal(sub_name // 'unrecognized c_type: ',                   &
     &      char=c_type)
      END SELECT ! Different coding depending on c_type

!  Actual Assignments    
      this % c_type(n) = TRIM(ADJUSTL(c_type))
      this % index = index
      this % value = value
      
      END SUBROUTINE recon_cnstrnts_append
      
!*******************************************************************************
! SECTION V. DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Destroy a recon_param
!-------------------------------------------------------------------------------
      SUBROUTINE recon_param_destroy(this)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (recon_param), INTENT(inout) :: this

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'recon_param_destroy: '

!  Start of executable code

!  Get rid of all components
      this % p_type = ''
      this % index = 0
      this % value = zero
      this % vrnc = one
      this % fdf = one

!  Different coding, depending on p_type
      SELECT CASE (TRIM(ADJUSTL(this % p_type)))
      CASE ('ac')
      CASE ('ac_aux_s')
      CASE ('ac_aux_f')
      CASE ('ai')
      CASE ('ai_aux_s')
      CASE ('ai_aux_f')
      CASE ('am')
      CASE ('am_aux_s')
      CASE ('am_aux_f')
      CASE ('bloat')
      CASE ('extcur')
      CASE ('curtor')
      CASE ('phiedge')
      CASE ('pres_scale')
      CASE ('density_max')
      CASE ('density_tau')

      CASE DEFAULT
         CALL err_fatal(sub_name // 'unrecognized p_type: ',                   &
     &      char=this % p_type)
      END SELECT ! Different coding depending on p_type

      END SUBROUTINE recon_param_destroy
!-------------------------------------------------------------------------------
!  Destroy a recon_cnstrnts
!-------------------------------------------------------------------------------
      SUBROUTINE recon_cnstrnts_destroy(this)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (recon_cnstrnts), INTENT(inout) :: this

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'recon_cnstrnt_destroy: '

!  Start of executable code

!  Get rid of all components
      this % n_cnstrnts = 0
      this % c_type = ''
      this % index = 0
      this % value = zero

      END SUBROUTINE recon_cnstrnts_destroy

!*******************************************************************************
! SECTION VI. CHANGE VALUE SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Change the value of a recon_param
!-------------------------------------------------------------------------------
      SUBROUTINE recon_param_change_value_s(this,new_value,delta_value)
!  Either a new value, or a change to the value, can be supplied.

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (recon_param), INTENT(inout)          :: this
      REAL(rprec), INTENT(in), OPTIONAL          :: new_value
      REAL(rprec), INTENT(in), OPTIONAL          :: delta_value

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'recon_param_change_value_s: '

!  Start of executable code

      IF(PRESENT(new_value) .AND. PRESENT(delta_value)) THEN
         CALL err_warn(sub_name // 'Both new_ and delta_ present')
      ENDIF
      
      IF (.NOT.(PRESENT(new_value) .OR. PRESENT(delta_value))) THEN
         CALL err_warn(sub_name // 'Neither new_ nor delta_ present')
      ENDIF
      
      IF (PRESENT(new_value)) THEN
         this % value = new_value
      ENDIF
      
      IF (PRESENT(delta_value)) THEN
         this % value = this % value + delta_value
      ENDIF

      END SUBROUTINE recon_param_change_value_s
      
!-------------------------------------------------------------------------------
!  Change the value of an array of recon_param
!-------------------------------------------------------------------------------      
      SUBROUTINE recon_param_change_value_a(this,new_value,delta_value)
!  Subroutine to change the values of an array of recon_param.
!  Either new values, or changes to the values, can be supplied.

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (recon_param), DIMENSION(:), INTENT(inout)   :: this
      REAL(rprec), DIMENSION(:), INTENT(in), OPTIONAL   :: new_value
      REAL(rprec), DIMENSION(:), INTENT(in), OPTIONAL   :: delta_value

!  Declare local variables
      INTEGER :: n1,n2, n
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'recon_param_change_value_a: '

!  Start of executable code
      n1 = SIZE(this,1)

      IF(PRESENT(new_value) .AND. PRESENT(delta_value)) THEN
         CALL err_warn(sub_name // 'Both new_ and delta_ present')
      ENDIF
      
      IF (.NOT.(PRESENT(new_value) .OR. PRESENT(delta_value))) THEN
         CALL err_warn(sub_name // 'Neither new_ nor delta_ present')
      ENDIF
      
      IF (PRESENT(new_value)) THEN
         n2 = SIZE(new_value,1)
         CALL assert_eq(n1,n2,sub_name,err_class='Warning')
         n = MIN(n1,n2)
         this(1:n) % value = new_value(1:n)
      ENDIF
      
      IF (PRESENT(delta_value)) THEN
         n2 = SIZE(delta_value,1)
         CALL assert_eq(n1,n2,sub_name,err_class='Warning')
         n = MIN(n1,n2)
         this(1:n) % value = this(1:n) % value + delta_value(1:n)
      ENDIF
      
      END SUBROUTINE recon_param_change_value_a
!-------------------------------------------------------------------------------
!  No need to change the value of a recon_cnstrnts
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  Increment a value, for finite differencing
!-------------------------------------------------------------------------------      
      SUBROUTINE recon_param_fd_increment(this,delta_rp)
!  Subroutine to change the value of a recon_param, for finite difference

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (recon_param), INTENT(inout)   :: this
      REAL(rprec), INTENT(out)            :: delta_rp

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'recon_param_fd_increment: '

!  Start of executable code
      delta_rp = this % fdf * this % vrnc
      this % value = this % value + delta_rp
      
      END SUBROUTINE recon_param_fd_increment
!-------------------------------------------------------------------------------
!  Get a finite difference factor from a recon_param
!-------------------------------------------------------------------------------      
      SUBROUTINE recon_param_fdf_get(this,fdf)
!  Subroutine to get the finite difference factor

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (recon_param), INTENT(inout)   :: this
      REAL(rprec), INTENT(out)            :: fdf

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'recon_param_fdf_get: '

!  Start of executable code
      fdf = this %  fdf 
      
      END SUBROUTINE recon_param_fdf_get
!-------------------------------------------------------------------------------
!  Put a finite difference factor to a recon_param
!-------------------------------------------------------------------------------      
      SUBROUTINE recon_param_fdf_put(this,fdf)
!  Subroutine to get the finite difference factor

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (recon_param), INTENT(inout)   :: this
      REAL(rprec), INTENT(in)            :: fdf

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'recon_param_fdf_put: '

!  Start of executable code
      this % fdf = fdf 
      
      END SUBROUTINE recon_param_fdf_put

!*******************************************************************************
! SECTION VII. ASSIGNMENT SUBROUTINES
!*******************************************************************************
!  recon_param is pretty simple, so the generic assign should work just fine.
!  recon_cnstrnts is pretty simple, so the generic assign should work just fine.

!*******************************************************************************
! SECTION VIII.  OUTPUT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Write out the contents of an array of recon_param
!-------------------------------------------------------------------------------

!  JDH 09-08-06. Need to make generic, called with scalar uses array of length 1.

      SUBROUTINE recon_param_write(this,identifier,unit,verbose)
      IMPLICIT NONE

!  Declare Arguments
      TYPE (recon_param), DIMENSION(:), INTENT (in) :: this
      CHARACTER (len=*), INTENT(in), OPTIONAL :: identifier
      INTEGER, INTENT(in), OPTIONAL :: unit
      INTEGER, INTENT(in), OPTIONAL :: verbose
!  identifier   character variable, also written out
!  unit         I/O unit number to write to
!  verbose      integer, to specify verbosity level of write

!  Declare local variables and constants
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'recon_param_write: '
      INTEGER :: iv_default = 1
      INTEGER :: iv
      INTEGER :: iou_default = 6
      INTEGER :: iou
      INTEGER :: i
      CHARACTER (len=60) :: id

!  Declare Format array
      CHARACTER(len=*), PARAMETER, DIMENSION(4) :: fmt1 = (/                   &
     & '(" start recon_param array write, called with id = ",a)    ',          &
     & '(" i     %p_type   %index    %value          %vrnc fdf")   ',          &
     & '(i4,2x,a10,2x,i4,2x,es12.5,2x,es12.5,2x,es12.5)            ',          &
     & '(" end recon_param array write, called with id = ",a)      '           &
     &  /) 

!  start of executable code
!  Check for arguments present
      IF (PRESENT(identifier)) THEN
         id = identifier
      ELSE
         id = ' '
      END IF

      IF (PRESENT(unit)) THEN
         iou = unit
      ELSE
         iou = iou_default
      END IF

      IF (PRESENT(verbose)) THEN
         iv = verbose
      ELSE
         iv = iv_default
      END IF
      
!  Select Case of Verbosity Level
      SELECT CASE(iv)
      CASE( :0)  ! VERY Terse
         WRITE(iou,fmt1(2))
         DO i = 1,size(this)
            WRITE(iou,fmt1(3)) i, this(i)%p_type, this(i)%index,               &
     &         this(i)%value, this(i)%vrnc, this(i)%fdf
         END DO
      
      CASE(1:)    ! Default, more verbose (add id at beginning and end)
         WRITE(iou,fmt1(1)) id
         WRITE(iou,fmt1(2))
         DO i = 1,size(this)
            WRITE(iou,fmt1(3)) i, this(i)%p_type, this(i)%index,               &
     &         this(i)%value, this(i)%vrnc, this(i)%fdf
         END DO
         WRITE(iou,fmt1(4)) id
      
      END SELECT

      END SUBROUTINE recon_param_write

!-------------------------------------------------------------------------------
!  Write out old, new, difference for array of recon_params
!-------------------------------------------------------------------------------

      SUBROUTINE recon_param_write_ond(old,new,identifier,unit,verbose)
      IMPLICIT NONE

!  Declare Arguments
      TYPE (recon_param), DIMENSION(:), INTENT (in) :: old, new
      CHARACTER (len=*), INTENT(in), OPTIONAL :: identifier
      INTEGER, INTENT(in), OPTIONAL :: unit
      INTEGER, INTENT(in), OPTIONAL :: verbose
!  identifier   character variable, also written out
!  unit         I/O unit number to write to
!  verbose      integer, to specify verbosity level of write

!  Declare local variables and constants
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'recon_param_write_ond: '
      INTEGER :: iv_default = 1
      INTEGER :: iv
      INTEGER :: iou_default = 6
      INTEGER :: iou
      INTEGER :: i, n_old, n_new, n
      CHARACTER (len=60) :: id

!  Declare Format array
!  (Formats did not fit on one line - use FORMAT statements)

!  start of executable code
!  Check for arguments present
      IF (PRESENT(identifier)) THEN
         id = identifier
      ELSE
         id = ' '
      END IF

      IF (PRESENT(unit)) THEN
         iou = unit
      ELSE
         iou = iou_default
      END IF

      IF (PRESENT(verbose)) THEN
         iv = verbose
      ELSE
         iv = iv_default
      END IF
      
!  Select Case of Verbosity Level
      SELECT CASE(iv)
      CASE( :-10)  ! VERY Terse
      
      CASE(-9:)    ! Default
         n_old = SIZE(old)
         n_new = SIZE(new)
         CALL assert_eq(n_old,n_new,sub_name // 
     &     ' old-new length discrepancy','Warn')
         n = MIN(n_old,n_new)
         WRITE(iou,1000)
         WRITE(iou,1100) (i,old(i) % p_type,old(i) % index,                    &
     &      old(i) % value, new(i) % value,                                    &
     &      new(i) % value - old(i) % value,i=1,n)
         WRITE(iou,1200)
         WRITE(iou,1100) (i,old(i) % p_type,old(i) % index,                    &
     &      old(i) % value / old(i) % vrnc,                                    &
     &      new(i) % value / new(i) % vrnc,                                    &
     &      (new(i) % value - old(i) % value) /                                &
     &      new(i) % vrnc,i=1,n)
1000  FORMAT(/'Reconstruction Parameters, old, new, difference '/              &
     &   'irp',t7,'type',t16,'index',t26,'old',t39,                            &
     &   'new',t52,'diff')
1100  FORMAT(i3,2x,a10,2x,i3,2x,es11.4,2x,es11.4,2x,es11.4)
1200  FORMAT(/'Normalized recon_params, old, new, difference'/                 &
     &   'irp',t7,'type',t16,'index',t26,'old_n',t39,                          &
     &   'new_n',t52,'diff_n')
      
      END SELECT

      END SUBROUTINE recon_param_write_ond
!-------------------------------------------------------------------------------
!  Write out the contents of a recon_cnstrnts
!-------------------------------------------------------------------------------

      SUBROUTINE recon_cnstrnts_write(this,identifier,unit,verbose)
      IMPLICIT NONE

!  Declare Arguments
      TYPE (recon_cnstrnts), INTENT (in) :: this
      CHARACTER (len=*), INTENT(in), OPTIONAL :: identifier
      INTEGER, INTENT(in), OPTIONAL :: unit
      INTEGER, INTENT(in), OPTIONAL :: verbose
!  identifier   character variable, also written out
!  unit         I/O unit number to write to
!  verbose      integer, to specify verbosity level of write

!  Declare local variables and constants
      INTEGER :: iv_default = 1
      INTEGER :: iv
      INTEGER :: iou_default = 6
      INTEGER :: iou
      INTEGER :: i
      CHARACTER (len=60) :: id

!  Declare Format array
      CHARACTER(len=*), PARAMETER, DIMENSION(4) :: fmt1 = (/                   &
     & '(" start recon_cnstrnts array write, called with id = ",a) ',          &
     & '(" i     %c_type   %index    %value               ")       ',          &
     & '(i4,2x,a10,2x,i4,2x,es12.5)                                ',          &
     & '(" end recon_cnstrnts array write, called with id = ",a)   '           &
     &  /) 

!  start of executable code
!  Check for arguments present
      IF (PRESENT(identifier)) THEN
         id = identifier
      ELSE
         id = ' '
      END IF

      IF (PRESENT(unit)) THEN
         iou = unit
      ELSE
         iou = iou_default
      END IF

      IF (PRESENT(verbose)) THEN
         iv = verbose
      ELSE
         iv = iv_default
      END IF
      
!  Select Case of Verbosity Level
      SELECT CASE(iv)
      CASE( :0)  ! VERY Terse
         WRITE(iou,fmt1(2))
         DO i = 1,this % n_cnstrnts
            WRITE(iou,fmt1(3)) i, this%c_type(i), this%index(i),               &
     &         this%value(i)
         END DO
      
      CASE(1:)    ! Default, more verbose (add id at beginning and end)
         WRITE(iou,fmt1(1)) id
         WRITE(iou,fmt1(2))
         DO i = 1,this % n_cnstrnts
            WRITE(iou,fmt1(3)) i, this%c_type(i), this%index(i),               &
     &         this%value(i)
         END DO
         WRITE(iou,fmt1(4)) id
      
      END SELECT

      END SUBROUTINE recon_cnstrnts_write

!*******************************************************************************
! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!  
!  JDH 08-16-2006
!     Modifying signal_T.f to get parameter_T.f
!
!  JDH 08-22-2006
!     Changed module name to recon_param_T
!
!  JDH 09-08-2006
!      Added recon_param_write, 
!
!  JDH 11-29-2006
!      Added p_types 'am' and 'pres_scale'
!
!  JDH 12-28-2006
!      Added change value subroutines
!
!  JDH 12-29-2006
!      Added p_types 'phiedge' and 'extcur'
!
!  JMS 7-25-2007
!      Added p_types 'density_max' and 'density_tau'
!      also increased type_len from 10 to 20
!
!  JDH 2007-10-03
!      Added type recon_cnstrnt, construct, destroy, write
!
!  JDH 2007-10-06
!      Changed to type recon_cnstrnts, put all constraints into one type.
!
!  JDH 2008-01-21
!      SPH eliminated iprec ~January 2008. Finish elimination.
!      Add finite difference factor component to the recon_param type
!      Subroutines recon_param_fd_increment, recon_param_fdf_get, recon_param_fdf_put
!
!  JDH 2010-12-06
!      Added more p_type s

      END MODULE recon_param_T
