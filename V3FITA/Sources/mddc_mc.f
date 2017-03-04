!*******************************************************************************
!  File mddc_mc.f
!  Contains module mddc_mc

!  Module is part of the V3FIT stellarator equilibrium reconstruction code.
!  The postfix _mc in the name indicates Model Compute
!  The module contains the subroutine mddc_mc_model_compute
!  Both an mddc_desc and a model are needed for this computation.

!*******************************************************************************
!  MODULE mddc_mc
!    (MDDC - Model Compute Calculations)
! SECTION I.     VARIABLE DECLARATIONS
! SECTION II.    INTERFACE BLOCKS
! SECTION III.   MAIN COMPUTATIONAL SUBROUTINE
!*******************************************************************************
      MODULE mddc_mc

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
      USE mddc_T
      
!-------------------------------------------------------------------------------
!  Model Derived Types
!-------------------------------------------------------------------------------
      USE model_T

!-------------------------------------------------------------------------------
!  Module for bivariate interpolation, in LIBSTELL/Sources/Modules
!-------------------------------------------------------------------------------
      USE bivariate

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
!    Information comes from the mddc_desc and the model
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  Actual computation of the model signal is in this subroutine
!    s_type = diagnostic
!      d_type = mddc
!        Magnetic Diagnostic-Dot Coil
!        signal_model_compute_mddc 
!           % data(1) - total flux
!           % data(2) - flux from the plasma currents
!           % data(3) - flux from all the external coils
!           % data(4) - data(3 + nextcur) - flux from individual external coils
!-------------------------------------------------------------------------------

      SUBROUTINE mddc_mc_model_compute(a_mddc,a_model,arr_data,                &
     &   arr_sigma)
      
!-------------------------------------------------------------------------------
!  Declarations
!-------------------------------------------------------------------------------
!  Declare Arguments 
      TYPE (mddc_desc), INTENT (inout), TARGET :: a_mddc
      TYPE (model), INTENT (inout), TARGET :: a_model
      REAL(rprec), POINTER, DIMENSION(:) :: arr_data, arr_sigma
!  Note - convention is that arr_data and arr_sigma
!    1) are deallocated in this subroutine, and then
!    2) are allocated in this subroutine

!  Declare local variables
      TYPE (eq_state), POINTER :: state => null()
      TYPE (mddc_mrf), POINTER :: mrf_local => null()
      LOGICAL :: lfirst
      INTEGER:: ib
!--   rgrid             : R on cylindrical grid in m                 --
!--   zgrid             : Z on cylindrical grid in m                 --
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: rgrid
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: zgrid

      REAL(rprec), DIMENSION(:,:), POINTER :: pl_response => null()
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: pl_respsuv
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: sumsuv
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: sum1k

      INTEGER :: nextcur, nextcur_2, nextcur_3
      INTEGER :: n_data
      INTEGER :: i, k, mc, index
      INTEGER :: ir, jz, kp, kp_store, ns, ju, kv
      REAL(rprec) :: fperiod, delr, delz, delp, dels, delu, delv,              &
     &    deluv, delsuv
      REAL(rprec) :: sumtot
      INTEGER :: istat1, istat2
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'mddc_mc_model_compute: '

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------

!  Point state to the model equilibrium state
      state => a_model % eqstate

!  ASSUME that mddc_mrf is defined. and is a TARGET
!  Point the local mrf to it.
      mrf_local => a_mddc % mrf

!  Check that aux2 is defined. (DEFINE of aux2 will define aux1 if needed)
!  (Should kv argument be kp or kpstore?)
      IF (.not. state % l_def_aux2) THEN
         CALL eq_aux_2_define(state,mrf_local % kp)  !! NEED TO FIX
      END IF

!  Check nextcur. Only need to worry about nextcur for free-boundary equilibria
      nextcur = state % fixp % nextcur
      nextcur_2 = size(mrf_local % rdiag_coilg_1,1)
      nextcur_3 = mrf_local % n_field_cg
      IF (state % fixp % lfreeb) THEN
         CALL assert_eq(nextcur,nextcur_2,nextcur_3,                           & 
     &      sub_name // 'nextcur discrepancy','Warn' )
      ENDIF
      nextcur = Min(nextcur,nextcur_2,nextcur_3)
      
!  Check for lstell_sym - lasym problem
!  (Only a problem when both are T, diagnostics assuming a symmetric equilibrium
!  but the equilibrium is asymmetric. Not a problem when the equilibrium is 
!  symmetric [lasym = F] yet the diagnostic does not assume a symmetric 
!  equilibrium [lstell_sym = F])
      IF (mrf_local % lstell_sym .and. state % fixp % lasym) THEN
         CALL err_fatal(sub_name // ' lstell_sym is T and lasym is T')
      ENDIF

!  Allocate data and sigma
      n_data = nextcur + 3
      IF (ASSOCIATED(arr_data)) THEN
         DEALLOCATE(arr_data,arr_sigma, stat=istat1)
         CALL assert_eq(0,istat1,sub_name // 'arr_data, sigma dealloc')
      ENDIF
      ALLOCATE(arr_data(n_data),arr_sigma(n_data),stat=istat1)
      CALL assert_eq(0,istat1,sub_name // 'arr_data, sigma alloc')

!----------------------------------------------------------------------
!  Flux due to Plasma
!----------------------------------------------------------------------
!  Much of this coding is copied from v3post.
!  First, set some of v3post's logical variables
!    Reset the grid
      lfirst = .true.
!  Do a volume, not surface integration
!  Comment out lsurf = .true. coding 
!  When try to add back lsurf = .true. coding, there may be more necessary than
!  just the commented out pieces.
      ib = 1

!----------------------------------------------------------------------
!-- set up (R,P,Z) cylindrical grids                                 --
!----------------------------------------------------------------------
!  Integers
      ir = mrf_local % ir
      jz = mrf_local % jz
      kp = mrf_local % kp
      kp_store = mrf_local % kp_store
      ns = size(state % aux2 % rsuv,1) ! Might this not be better set to currusuv?
      ju = size(state % aux2 % rsuv,2)
      kv = size(state % aux2 % rsuv,3)

      CALL assert_eq(size(state % aux2 % rsuv,3),kv,sub_name //                &
     &   'kv size discrepancy')

!  Reals
      fperiod = twopi / mrf_local % n_field_periods
      delr = (mrf_local % rmax -mrf_local %  rmin) /                           &
     &   (ir - 1)
      delz = (mrf_local % zmax - mrf_local % zmin) /                           &
     &   (jz - 1)
      delp = fperiod / kp
      dels = one / (ns - 1)
      delu = twopi / ju
      delv = fperiod / kv

!  Allocation
      ALLOCATE (rgrid(ir),stat = istat1)
      ALLOCATE (zgrid(jz),stat = istat2)
      CALL assert_eq(0,istat1,istat2,sub_name // 'r,z grid alloc')
      ALLOCATE(pl_respsuv(ns,ju,3),stat=istat1)
      ALLOCATE (sumsuv(ns, ju),stat=istat2)
      CALL assert_eq(0,istat1,istat2,sub_name // 'pl_, sum suv alloc')
      ALLOCATE (sum1k(kp_store),stat=istat1)
      CALL assert_eq(0,istat1,sub_name // 'sum1k alloc')

!  Fill arrays
      DO i = 1, ir
        rgrid(i) = mrf_local % rmin + (i - 1) * delr
      END DO
      DO i = 1, mrf_local % jz
        zgrid(i) = mrf_local % zmin + (i - 1) * delz
      END DO

!----------------------------------------------------------------------
!-- Loop over toroidal planes                                        --
!----------------------------------------------------------------------
      sumtot = zero
      TOR_PLANE: DO k = 1, kp_store

         DO mc = 1, 3
            SELECT CASE (mc)
               CASE (1) 
                  pl_response => mrf_local % a_r(:,:,k)
               CASE (2) 
                  pl_response => mrf_local % a_f(:,:,k)
               CASE (3) 
                  pl_response => mrf_local % a_z(:,:,k)
            END SELECT
!----------------------------------------------------------------------
!-- bivariate (4 pt) interpolation response tables                  ---
!----------------------------------------------------------------------
            IF (lfirst .and. mc.eq.1) THEN
               CALL setbivariate (state % aux2 % rsuv(:,:,k),                  &
     &                            state % aux2 % zsuv(:,:,k),                  &
     &                            rgrid, zgrid, k, kv, ib)
            END IF
!            IF (.not.(lsurf)) THEN
               CALL bivariate4pt (pl_response, pl_respsuv(1,1,mc), k)
!            ELSE
!               CALL bivariate4pt (pl_response, pl_respns(1), k)
!               pl_respsuv(ns,:,mc) = pl_respns(:)
!            END IF
         END DO ! mc
!----------------------------------------------------------------------
!-- form plasma response integral (R=1, PHI=2, Z=3)                  --
!----------------------------------------------------------------------
!         IF (lsurf) THEN
!            sumsuv(ns,:) = pl_respsuv(ns,:,1)                                  &
!     &                    *(bsubuns(:,k)*rvsuv(ns,:,k)                         & 
!     &                    - bsubvns(:,k)*rusuv(ns,:,k))                        & 
!     &                    + pl_respsuv(ns,:,2)                                 &
!     &                    * bsubuns(:,k)*rsuv(ns,:,k)                          &
!     &                    + pl_respsuv(ns,:,3)                                 &
!     &                    *(bsubuns(:,k)*zvsuv(ns,:,k)                         & 
!     &                    - bsubvns(:,k)*zusuv(ns,:,k))&
!            sum1k(k) = SUM(sumsuv(ns,:))
!          ELSE
             sumsuv(1:ns,1:ju) = state % aux2 % currusuv(:,:,k)                &
     &          *(state % aux2 % rusuv(:,:,k)*pl_respsuv(:,:,1)                & 
     &          + state % aux2 % zusuv(:,:,k)*pl_respsuv(:,:,3))               & 
     &          + state % aux2 % currvsuv(:,:,k)                               & 
     &          *(state % aux2 % rvsuv(:,:,k)*pl_respsuv(:,:,1)                & 
     &          + state % aux2 % zvsuv(:,:,k)*pl_respsuv(:,:,3)                & 
     &          + state % aux2 % rsuv (:,:,k)*pl_respsuv(:,:,2))
             sum1k(k) =  SUM(sumsuv(2:ns-1,1:ju)) +                            &
     &          SUM(sumsuv(1,1:ju) + sumsuv(ns,1:ju))/2
!          END IF
      END DO TOR_PLANE
        
      sumtot = SUM(sum1k)
      
      IF (mrf_local % lstell_sym) THEN
         sumtot = sumtot - (sum1k(1) + sum1k(kp_store)) / 2 
      END IF
      
      deluv = delu*delv
      delsuv = dels*deluv

!      IF (lsurf) THEN
!         sumtot = -sumtot*deluv/mu0
!      ELSE
         sumtot = -sumtot*delsuv
!      END IF

      arr_data(2) = sumtot

!----------------------------------------------------------------------
!  Flux due to Coils - only need for free-boundary case
!----------------------------------------------------------------------
!  mgrid_mode       variable to specify mode for "Scaled' or "Raw"
!       In "S" mode, B-fields are "unit" current responses, so EXTCUR
!       has units of A. The responses are obtained from SUM(inductance * EXTCUR)
!       In "R" mode, B-fields are the true fields corresponding to
!       the currents in the coils-dot file, so EXTCUR ~ unity dimensionless
!       multiplier. Thus, the responses are SUM(rdiag_coil * EXTCUR)

      IF (state % fixp % lfreeb) THEN
         SELECT CASE (state % aux1 % mgrid_mode)
         CASE ('S','s', 'N', 'n')  ! Scaled. Per Ed Lazarus Request, (V3POST 
                                ! comment of 07.15.04), also treat 'N' this way
                                ! Uses "inductance"
            DO i = 1,nextcur
              index = 3 + i
              arr_data(index) = mrf_local % rdiag_coilg_1(i) *                 & 
     &            state % varp % extcur(i) / mrf_local % extcur_mg(i)
            END DO
      
         CASE DEFAULT ! Includes Raw ('R','r') 
            DO i = 1,nextcur
              index = 3 + i
              arr_data(index) = mrf_local % rdiag_coilg_1(i) *                 & 
     &            state % varp % extcur(i)
            END DO
         END SELECT
      ENDIF
      
!----------------------------------------------------------------------
!  Total Flux, Coil Flux, Sigmas
!----------------------------------------------------------------------
      
      IF (state % fixp % lfreeb) THEN
         arr_data(3) = SUM(arr_data(4:n_data))
      ELSE
         arr_data(3) = zero
      ENDIF
      arr_data(1) = arr_data(2) + arr_data(3)
      arr_sigma(1:n_data) = zero

!----------------------------------------------------------------------
!  DEALLOCATE space
!----------------------------------------------------------------------
      DEALLOCATE(rgrid,stat = istat1)
      CALL assert_eq(0,istat1,sub_name // 'dealloc rgrid')

      DEALLOCATE(zgrid,stat = istat1)
      CALL assert_eq(0,istat1,sub_name // 'dealloc zgrid')

      DEALLOCATE(pl_respsuv,stat = istat1)
      CALL assert_eq(0,istat1,sub_name // 'dealloc pl_respsuv')

      DEALLOCATE(sumsuv,stat = istat1)
      CALL assert_eq(0,istat1,sub_name // 'dealloc sumsuv')
      
      DEALLOCATE(sum1k,stat = istat1)
      CALL assert_eq(0,istat1,sub_name // 'dealloc sum1k')
      
      RETURN
      END SUBROUTINE mddc_mc_model_compute

!*******************************************************************************
! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!
!  JDH 2010-11-15. Added logic for VMEC fixed boundary case (lfreeb = .false.) 
!
!  JDH 2009-03-03. First version of module. Based on coding from previous
!    module signal_model
!
      END MODULE mddc_mc
