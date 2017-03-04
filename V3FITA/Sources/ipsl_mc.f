!*******************************************************************************
!  File ipsl_mc.f
!  Contains module ipsl_mc

!  Module is part of the V3FIT stellarator equilibrium reconstruction code.
!  The postfix _mc in the name indicates Model Compute
!  The module contains the subroutine ipsl_mc_model_compute
!  Both an ipsl_desc and a model are needed for this computation.

!*******************************************************************************
!  MODULE ipsl_mc
!    (MDDC - Model Compute Calculations)
! SECTION I.     VARIABLE DECLARATIONS
! SECTION II.    INTERFACE BLOCKS
! SECTION III.   MAIN COMPUTATIONAL SUBROUTINE
!*******************************************************************************
      MODULE ipsl_mc

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Type declarations, constants, utilities
!-------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
      USE v3_utilities

!-------------------------------------------------------------------------------
!  MDDC Derived Types
!-------------------------------------------------------------------------------
      USE ipsl_T
      
!-------------------------------------------------------------------------------
!  Model Derived Types
!-------------------------------------------------------------------------------
      USE model_T

!-------------------------------------------------------------------------------
!  IPSL specific modules
!-------------------------------------------------------------------------------
!      USE ip_beamline
      USE ipsl_T
      USE ip_global
      USE faraday_mod, ONLY: ipsl_integral

      IMPLICIT NONE
      
!*******************************************************************************
! SECTION II.  INTERFACE BLOCKS
!*******************************************************************************

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      CONTAINS
          
!*******************************************************************************
! SECTION III.  MAIN COMPUTATIONAL SUBROUTINE
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Compute an MDDC signal
!
!    Information comes from the ipsl_desc and the model
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  Actual computation of the model signal is in this subroutine
!    s_type = diagnostic
!      d_type = ipsl
!        Magnetic Diagnostic-Dot Coil
!        signal_model_compute_ipsl 
!           % data(1) - total flux
!           % data(2) - flux from the plasma currents
!           % data(3) - flux from all the external coils
!           % data(4) - data(3 + nextcur) - flux from individual external coils
!
!-------------------------------------------------------------------------------

      SUBROUTINE ipsl_mc_model_compute(a_ipsl,a_model,arr_data,                &
     &   arr_sigma)
!-------------------------------------------------------------------------
!      
! FUNCTION: computes signals for interferometry/polarimetry straight
!           line" signals  (s_type = ipsl)
!
!           % data(1) -  Faraday Rotation angle      (ipsl_type = 'polarim')  OR
!                     -  Interferometry phase angle  (ipsl_type = 'interf')
!
!! created 6/19/07 by J. Shields
!-------------------------------------------------------------------------

!..........passed Arguments................................! 
      TYPE (ipsl_desc), INTENT (inout), TARGET :: a_ipsl
      TYPE (model), INTENT (inout), TARGET :: a_model    ! equilbrium info
      REAL(rprec), POINTER, DIMENSION(:) :: arr_data, arr_sigma
!  Note - convention is that arr_data and arr_sigma
!    1) are deallocated in this subroutine, and then
!    2) are allocated in this subroutine

!.........Local Variables...............................!
      TYPE (eq_state), POINTER :: state => null()
      LOGICAL :: lfirst
      INTEGER :: n_data
      INTEGER :: i, k, mc, index
      REAL(rprec) :: line_integral_out
      INTEGER :: istat1, istat2
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'ipsl_mc_model_compute: '

!      write(*,*) 'Now entering ', sub_name

!  Allocate data and sigma
      n_data = 1
      IF (ASSOCIATED(arr_data)) THEN
         DEALLOCATE(arr_data,arr_sigma, stat=istat1)
         CALL assert_eq(0,istat1,sub_name // 'arr_data, sigma dealloc')
      ENDIF
      ALLOCATE(arr_data(n_data),arr_sigma(n_data),stat=istat1)
      CALL assert_eq(0,istat1,sub_name // 'arr_data, sigma alloc')

      call IPSL_INTEGRAL(a_ipsl, line_integral_out)
      arr_data(1) = line_integral_out
      arr_sigma(1) = zero

      write(*,*) sub_name, 'output line integral = ', arr_data(1)

      RETURN
      END SUBROUTINE ipsl_mc_model_compute

!===================================================================================
      SUBROUTINE ipsl_model_init()
!  2009-03-03 JDH. Could this subroutine be eliminated?

!--------------------------------------------------------------------------------
!
! Function:  Tests VMEC/V3FIT interfacing routines for B transformations
!
! 
! 
!
!
!
!----------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
!      USE v3f_global   ! module containing V3FIT global data
      IMPLICIT NONE

!!..............dummy variables.............................................!
!      INTEGER, INTENT(in) :: i_beam     ! index of interferometer beam being integrated over

!............local variables.................................................!
      INTEGER      :: i, j, kp_junk
      REAL(rprec)               ::   blah
      CHARACTER(len=*), PARAMETER :: subname = 'V3F_GLOBAL_TEST: '




!      write(*,*) '===================================================='    
!      write(*,*) ' '
      write(*,*) 'Hi.  Now in ', subname

!     call VMEC_B_INIT


      RETURN
      END SUBROUTINE ipsl_model_init

!*******************************************************************************
! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!
!  JDH 2009-03-03. First version of module. Based on coding from previous
!    module signal_model
!
      END MODULE ipsl_mc
