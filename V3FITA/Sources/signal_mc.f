!*******************************************************************************
!  File signal_mc.f
!  Contains module signal_mc

!  Module is part of the V3FIT stellarator equilibrium reconstruction code.
!  The postfix _mc in the name indicates Model Compute
!  The module contains the subroutine signal_mc_model_compute
!  Both an signal_desc and a model are needed for this computation.

!*******************************************************************************
!  MODULE signal_mc
!    (Signal - Model Calculations)
! SECTION I.     VARIABLE DECLARATIONS
! SECTION II.    INTERFACE BLOCKS
! SECTION III.   MAIN COMPUTATIONAL SUBROUTINE
!*******************************************************************************
      MODULE signal_mc

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
!  Signal Derived Types
!-------------------------------------------------------------------------------
      USE signal_T
      
!-------------------------------------------------------------------------------
!  Model Derived Types
!-------------------------------------------------------------------------------
      USE model_T

!-------------------------------------------------------------------------------
!  Extra Modules, for IPSL stuff
!-------------------------------------------------------------------------------
      USE read_wout_mod, ONLY: read_wout_file      
      USE b_transform_vmec, ONLY: VMEC_B_INIT

!-------------------------------------------------------------------------------
!  Other _mc modules
!-------------------------------------------------------------------------------
      USE mddc_mc
      USE ipsl_mc
      USE sxrc_mc           !GJH 2010-01-20
      USE edge_limit_mc
      
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
!  Compute a model signal
!
!    Information comes from the signal and the equilibrium
!-------------------------------------------------------------------------------
      SUBROUTINE signal_mc_model_compute(sdata,sdesc,a_model)

!-------------------------------------------------------------------------------
!  Declarations
!-------------------------------------------------------------------------------
!  Declare Arguments 
      TYPE (signal_data), INTENT (inout) :: sdata
      TYPE (signal_desc), TARGET, INTENT (inout) :: sdesc
      TYPE (model), INTENT (inout) :: a_model

!  Declare local variables
      INTEGER :: ierr
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'signal_mc_model_compute: '

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Things common to all s_type
!-------------------------------------------------------------------------------

!  Deallocate data and sigma
      IF (ASSOCIATED(sdata % data)) THEN
         DEALLOCATE(sdata % data,sdata % sigma)
      ENDIF

!  Set signal_desc pointer
      sdata % desc => sdesc

!  Set sd_type
      sdata % sd_type = 'model'

!-------------------------------------------------------------------------------
!  Different coding, depending on s_type
!-------------------------------------------------------------------------------
      SELECT CASE (TRIM(ADJUSTL(sdesc % s_type)))
      CASE ('diagnostic')
         SELECT CASE (TRIM(ADJUSTL(sdesc % diag % d_type)))
         CASE ('mddc')
             CALL mddc_mc_model_compute(sdesc % diag % mddc,a_model,           &
     &               sdata % data,sdata % sigma)
         CASE ('ipsl')
             CALL read_wout_file('wout_.nc', ierr)
             if (ierr .ne. 0 ) write(*,*) sub_name // 'read_wout error'
             call VMEC_B_INIT(1)
             CALL ipsl_mc_model_compute(sdesc % diag % ipsl,a_model,           &
     &               sdata % data,sdata % sigma)
         CASE ('sxrc')
             CALL sxrc_mc_model_compute(sdesc % diag % sxrc,a_model,           &
     &               sdata % data, sdata % sigma)
         CASE DEFAULT
            CALL err_fatal(sub_name // 'unrecognized d_type: ',                &
     &         sdesc % diag % d_type)
         END SELECT ! Select Case on d_type
      
      CASE ('geometric')
         SELECT CASE (TRIM(ADJUSTL(sdesc % geom % g_type)))
         CASE ('edge_limit')
             CALL edge_limit_mc_model_compute(sdesc % geom % el_desc,          &
     &               a_model,sdata % data,sdata % sigma)
         CASE DEFAULT
            CALL err_fatal(sub_name // 'unrecognized g_type: ',                &
     &         sdesc % geom % g_type)
         END SELECT ! Select Case on g_type
      
      CASE DEFAULT
         CALL err_fatal(sub_name // 'unrecognized s_type: ',                   &
     &      char=sdesc % s_type)
      END SELECT ! Different coding depending on d_type

      RETURN
      END SUBROUTINE signal_mc_model_compute


!*******************************************************************************
! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!  
!  JDH 11-28-2004
!     First version of module
!
!  JDH 12-14-2004
!     Moved from signal_eq to signal_model
!
!  JDH 2007-06-12
!     Changed from signal_mrf to mddc_mrf located in diagnostic_desc
!
!  JDH 2007-10-05
!     Slight change to error trapping in compute_mddc
!
!  JDH 2008-01-19 - Added => null() to pointer declarations
!
!  JDH 2009-01-30. Added subroutines signal_model_compute_edge_limit and
!     signal_model_aux_rz_one.
!
!  JDH 2009-03-03. Revised structure, removed stuff, changed name
!    from signal_model to signal_mc
!
!  JDH 2010-06-11. Added GJH coding for sxrc.
!

      END MODULE signal_mc
