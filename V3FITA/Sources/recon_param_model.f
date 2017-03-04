!*******************************************************************************
!  File recon_param_model.f
!  Contains module recon_param_model

!  Module is part of the V3FIT stellarator equilibrium reconstruction code.
!  It deals with the interaction of reconstruction parameters and models

!*******************************************************************************
!  MODULE recon_param_model
!    (recon_param - model Calculations)
! SECTION I.     VARIABLE DECLARATIONS
! SECTION II.    INTERFACE BLOCKS
! SECTION III.   GET AND PUT SUBROUTINES
! SECTION X.     COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
      MODULE recon_param_model

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Type declarations - lengths of reals, integers, and complexes.
!  Frequently used mathematical constants, lots of extra precision.
!-------------------------------------------------------------------------------
      USE stel_kinds, only : rprec, iprec, cprec
      USE stel_constants, only : pi, twopi, one, zero
      
!-------------------------------------------------------------------------------
!  Reconstruction Parameter Derived Types
!-------------------------------------------------------------------------------
      USE recon_param_T
      
!-------------------------------------------------------------------------------
!  Model Derived Types
!-------------------------------------------------------------------------------
      USE model_T
      
!-------------------------------------------------------------------------------
!  V3 Utilities
!-------------------------------------------------------------------------------
      USE v3_utilities

!-------------------------------------------------------------------------------
!  Equilibrium Interface
!-------------------------------------------------------------------------------
      USE eq_interface

      IMPLICIT NONE
      
!*******************************************************************************
! SECTION II.  INTERFACE BLOCKS
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Put parameters into the model, scalar and array of recon_param
!-------------------------------------------------------------------------------
      INTERFACE recon_param_model_put
         MODULE PROCEDURE recon_param_model_put_s,                             & 
     &      recon_param_model_put_a
      END INTERFACE
      CONTAINS
          
!*******************************************************************************
! SECTION III.  GET AND PUT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get values of an array of recon_param from a model
!-------------------------------------------------------------------------------
      SUBROUTINE recon_param_model_get(rparam,a_model)

!-------------------------------------------------------------------------------
!  Declarations
!-------------------------------------------------------------------------------
!  Declare Arguments 
      TYPE (recon_param), DIMENSION(:), INTENT(inout) :: rparam
      TYPE (model), TARGET, INTENT(in) :: a_model

!  Declare local variables
      TYPE (eq_state), POINTER :: state => null()
      INTEGER                  :: i, iac, iai, iam
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'recon_param_model_get: '

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
!
!  Point state to the model equilibrium state
      state => a_model % eqstate

!-------------------------------------------------------------------------------
!  Loop over rparam array
!-------------------------------------------------------------------------------

      DO i = 1, size(rparam)
         SELECT CASE (TRIM(ADJUSTL(rparam(i) % p_type)))
         CASE ('ac')
            iac = rparam(i) % index
            CALL assert_eq(0,MIN(0,iac),sub_name // 'iac range min1')
            CALL assert_eq(20,MAX(20,iac),sub_name // 'iac range max2')
            rparam(i) % value = state % varp % ac(iac)

         CASE ('ac_aux_s')
            iac = rparam(i) % index
            CALL assert_eq(1,MIN(1,iac),sub_name // 'iac range min3')
            CALL assert_eq(101,MAX(101,iac),sub_name //                        &
     &         'iac range max4')
            rparam(i) % value = state % varp % ac_aux_s(iac)

         CASE ('ac_aux_f')
            iac = rparam(i) % index
            CALL assert_eq(1,MIN(1,iac),sub_name // 'iac range min5')
            CALL assert_eq(101,MAX(101,iac),sub_name //                        &
     &         'iac range max6')
            rparam(i) % value = state % varp % ac_aux_f(iac)

         CASE ('ai')
            iai = rparam(i) % index
            CALL assert_eq(0,MIN(0,iai),sub_name // 'iai range min7')
            CALL assert_eq(20,MAX(20,iai),sub_name // 'iai range max8')
            rparam(i) % value = state % varp % ai(iai)

         CASE ('ai_aux_s')
            iai = rparam(i) % index
            CALL assert_eq(1,MIN(1,iai),sub_name // 'iai range min9')
            CALL assert_eq(101,MAX(101,iai),sub_name //                        &
     &         'iai range max10')
            rparam(i) % value = state % varp % ai_aux_s(iai)

         CASE ('ai_aux_f')
            iai = rparam(i) % index
            CALL assert_eq(1,MIN(1,iai),sub_name // 'iai range min11')
            CALL assert_eq(101,MAX(101,iai),sub_name //                        &
     &         'iai range max12')
            rparam(i) % value = state % varp % ai_aux_f(iai)

         CASE ('am')
            iam = rparam(i) % index
            CALL assert_eq(0,MIN(0,iam),sub_name // 'iam range min13')
            CALL assert_eq(20,MAX(20,iam),sub_name // 'iam rang max14')
            rparam(i) % value = state % varp % am(iam)

         CASE ('am_aux_s')
            iam = rparam(i) % index
            CALL assert_eq(1,MIN(1,iam),sub_name // 'iam range min15')
            CALL assert_eq(101,MAX(101,iam),sub_name //                        &
     &         'iam range max16')
            rparam(i) % value = state % varp % am_aux_s(iam)

         CASE ('am_aux_f')
            iam = rparam(i) % index
            CALL assert_eq(1,MIN(1,iam),sub_name // 'iam range min17')
            CALL assert_eq(101,MAX(101,iam),sub_name //                        &
     &         'iam range max18')
            rparam(i) % value = state % varp % am_aux_f(iam)

         CASE ('bloat')
            rparam(i) % value = state % varp % bloat

          CASE ('curtor')
            rparam(i) % value = state % varp % curtor
                       
         CASE ('pres_scale')
            rparam(i) % value = state % varp % pres_scale
            
         CASE ('phiedge')
            rparam(i) % value = state % varp % phiedge

         CASE ('extcur')
            iam = rparam(i) % index
            CALL assert_eq(0,MIN(0,iam),sub_name // 'iam range min')
            CALL assert_eq(100,MAX(100,iam),sub_name //'iam range max')
            rparam(i) % value = state % varp % extcur(iam)
            
         CASE ('density_max')
            rparam(i) % value = a_model % density % ne_max

         CASE ('density_tau')
            rparam(i) % value = a_model % density % tau

         CASE DEFAULT
            CALL err_fatal(sub_name // 'unrecognized p_type: ',                &
     &         char=rparam(i) % p_type)
         END SELECT ! Different coding depending on p_type
      END DO

      RETURN
      END SUBROUTINE recon_param_model_get

!-------------------------------------------------------------------------------
!  Put values of an array of recon_param into a model
!-------------------------------------------------------------------------------
      SUBROUTINE recon_param_model_put_a(rparam_a,a_model,rcons,cvmi)

!-------------------------------------------------------------------------------
!  Declarations
!-------------------------------------------------------------------------------
!  rparam_a    array of recon_param
!  a_model     a model
!  rcons       a recon_cnstrnts - reconstruction constraints
!  cvmi        optional character argument. IF .eq. 'VMI', then parameters are
!                propagated to the VMEC internal state.

!  Declare Arguments 
      TYPE (recon_param), DIMENSION(:), INTENT(in) :: rparam_a
      TYPE (model), TARGET, INTENT (inout) :: a_model
      TYPE (recon_cnstrnts), INTENT(in) :: rcons
      CHARACTER(len=*), INTENT(in), OPTIONAL :: cvmi

!  Declare local variables
      TYPE (eq_state), POINTER :: state => null()
      INTEGER                  :: i, iac, iai, iam
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'recon_param_model_put: '

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
!  Point state to the model equilibrium state
      state => a_model % eqstate

!-------------------------------------------------------------------------------
!  Loop over rparam_a array
!-------------------------------------------------------------------------------

      DO i = 1, size(rparam_a)
         SELECT CASE (TRIM(ADJUSTL(rparam_a(i) % p_type)))

         CASE ('ac')
            iac = rparam_a(i) % index
            CALL assert_eq(0,MIN(0,iac),sub_name // 'iac range min1')
            CALL assert_eq(20,MAX(20,iac),sub_name // 'iac range max2')
            state % varp % ac(iac) = rparam_a(i) % value

         CASE ('ac_aux_s')
            iac = rparam_a(i) % index
            CALL assert_eq(1,MIN(1,iac),sub_name // 'iac range min3')
            CALL assert_eq(101,MAX(101,iac),sub_name //                        &
     &         'iac range max4')
            state % varp % ac_aux_s(iac) = rparam_a(i) % value

         CASE ('ac_aux_f')
            iac = rparam_a(i) % index
            CALL assert_eq(1,MIN(1,iac),sub_name // 'iac range min5')
            CALL assert_eq(101,MAX(101,iac),sub_name //                        &
     &         'iac range max6')
            state % varp % ac_aux_f(iac) = rparam_a(i) % value

         CASE ('ai')
            iai = rparam_a(i) % index
            CALL assert_eq(0,MIN(0,iai),sub_name // 'iai range min7')
            CALL assert_eq(20,MAX(20,iai),sub_name // 'iai range max8')
            state % varp % ai(iai) = rparam_a(i) % value

         CASE ('ai_aux_s')
            iai = rparam_a(i) % index
            CALL assert_eq(1,MIN(1,iai),sub_name // 'iai range min9')
            CALL assert_eq(101,MAX(101,iai),sub_name //                        &
     &         'iai range max10')
            state % varp % ai_aux_s(iai) = rparam_a(i) % value

         CASE ('ai_aux_f')
            iai = rparam_a(i) % index
            CALL assert_eq(1,MIN(1,iai),sub_name // 'iai range min11')
            CALL assert_eq(101,MAX(101,iai),sub_name //                        &
     &         'iai range max12')
            state % varp % ai_aux_f(iai) = rparam_a(i) % value

         CASE ('am')
            iam = rparam_a(i) % index
            CALL assert_eq(0,MIN(0,iam),sub_name // 'iam range min13')
            CALL assert_eq(20,MAX(20,iam),sub_name // 'iam rang max14')
            state % varp % am(iam) = rparam_a(i) % value

         CASE ('am_aux_s')
            iam = rparam_a(i) % index
            CALL assert_eq(1,MIN(1,iam),sub_name // 'iam range min15')
            CALL assert_eq(101,MAX(101,iam),sub_name //                        &
     &         'iam range max16')
            state % varp % am_aux_s(iam) = rparam_a(i) % value

         CASE ('am_aux_f')
            iam = rparam_a(i) % index
            CALL assert_eq(1,MIN(1,iam),sub_name // 'iam range min17')
            CALL assert_eq(101,MAX(101,iam),sub_name //                        &
     &         'iam range max18')
            state % varp % am_aux_f(iam) = rparam_a(i) % value

         CASE ('bloat')
            state % varp % bloat = rparam_a(i) % value

         CASE ('curtor')
            state % varp % curtor = rparam_a(i) % value
            
         CASE ('pres_scale')
            state % varp % pres_scale = rparam_a(i) % value
            
         CASE ('phiedge')
            state % varp % phiedge = rparam_a(i) % value
            
         CASE ('extcur')
            iam = rparam_a(i) % index
            CALL assert_eq(0,MIN(0,iam),sub_name // 'iam range min')
            CALL assert_eq(100,MAX(100,iam),sub_name //'iam range max')
            state % varp % extcur(iam) = rparam_a(i) % value

         CASE ('density_max')
            a_model % density % ne_max = rparam_a(i) % value

         CASE ('density_tau')
            a_model % density % tau = rparam_a(i) % value

         CASE DEFAULT
            CALL err_fatal(sub_name // 'unrecognized p_type: ',                &
     &         char=rparam_a(i) % p_type)
         END SELECT ! Different coding depending on p_type
         
      END DO
      
!-------------------------------------------------------------------------------
!  Impose Constraints
!-------------------------------------------------------------------------------

!  JDH 2010-12-06 To DO - Check for correct pxxx_type before do constraints.

      DO i = 1, rcons % n_cnstrnts
!  Different coding, depending on c_type
         SELECT CASE (TRIM(ADJUSTL(rcons % c_type(i))))
         CASE ('ac_edge')
            state % varp % ac(rcons % index(i)) = rcons % value(i) -           &
     &         sum(state % varp % ac(0:rcons % index(i)-1))
            state % varp % ac(rcons % index(i)+1:) = zero
         CASE ('am_edge')
            state % varp % am(rcons % index(i)) = rcons % value(i) -           &
     &         sum(state % varp % am(0:rcons % index(i)-1))
            state % varp % am(rcons % index(i)+1:) = zero

         CASE DEFAULT
            CALL err_fatal(sub_name // 'unrecognized c_type: ',                &
     &         char=rcons % c_type(i), int=i)
         END SELECT ! Different coding depending on c_type
      END DO
         
!  Check last argument, to see if need to propagate parameters to the
!  VMEC internal state
      IF (PRESENT(cvmi)) THEN
         IF (TRIM(ADJUSTL(cvmi)) .eq. 'VMI') THEN
            CALL eq_change_vmi_varp(state)
         ELSE
             WRITE(*,*) sub_name, ' cvmi = ', cvmi
         ENDIF
      ELSE
         WRITE(*,*) sub_name, ' Argument cvmi not present'
      ENDIF

      RETURN
      END SUBROUTINE recon_param_model_put_a

!-------------------------------------------------------------------------------
!  Put values of a scalar recon_param into a model
!-------------------------------------------------------------------------------
      SUBROUTINE recon_param_model_put_s(rparam_s,a_model,rcons,cvmi)

!-------------------------------------------------------------------------------
!  Declarations
!-------------------------------------------------------------------------------
!  rparam_s    scalar recon_param
!  a_model     a model
!  rcons       a recon_cnstrnts - reconstruction constraints
!  cvmi        optional character argument. IF .eq. 'VMI', then parameters are
!                propagated to the VMEC internal state.

!  Declare Arguments 
      TYPE (recon_param), INTENT(in) :: rparam_s
      TYPE (model), TARGET, INTENT (inout) :: a_model
      TYPE (recon_cnstrnts), INTENT(in) :: rcons
      CHARACTER(len=*), INTENT(in), OPTIONAL :: cvmi

!  Declare local variables
      TYPE (recon_param), DIMENSION(1) :: rparam_a
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'recon_param_model_put_s: '

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
!  Copy the scalar into an array of length 1, and call the array subroutine
      rparam_a(1) = rparam_s
      CALL recon_param_model_put_a(rparam_a,a_model,rcons,cvmi)

      RETURN
      END SUBROUTINE recon_param_model_put_s

!*******************************************************************************
! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!  
!  JDH 09-08-2006
!     First version of module
!
!  JDH 11-28-2006
!     Added am array and pres_scale
!
!  JDH 11-29-06
!     index bug with am fixed.
!
!  JDH 12-29-06
!     Added phiedge and extcur
!
!  JMS 07-25-07
!     Added density_max and density_tau
!
!  JDH 2007-10-06
!     Added argument for reconstruction constraints to recon_param_model_put
!
!  JDH 2008-01-19 - Added => null() to pointer declarations
!
!  JDH 2008-08-08 - Changed recon_param_model_put to recon_param_model_put_a
!      Added recon_param_model_put_s. Scalars and arrays.
!
!  JDH 2010-12-06
!      Added more p_type s.
           
      END MODULE recon_param_model
