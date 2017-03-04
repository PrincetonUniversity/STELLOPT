!     SPH011408: Replace all INTEGER(iprec) with INTEGER
!*******************************************************************************
!  File eq_interface.f
!  Contains module eq_interface

!  Module is part of the V3FIT stellarator equilibrium reconstruction code.
!  It deals with the interface with the equilibrium code.
!  Right now (6/2004), the VMEC code is the only equilibrium code.

!*******************************************************************************
!  MODULE eq_mod
!    (Equilibrium)
! SECTION I.      VARIABLE DECLARATIONS
! SECTION II.     INTERFACE BLOCKS
! SECTION III.    EQUILIBRIUM INITIALIZATION SUBROUTINES
! SECTION IV.     EQUILIBRIUM STEP ROUTINES
! SECTION V.      EQUILIBRIUM CHANGE SUBROUTINES
! SECTION VI.     EQUILIBRIUM AUXILIARY CALCULATION ROUTINES
! SECTION VII.    EQUILIBRIUM WRITE ROUTINES
!*******************************************************************************
      MODULE eq_interface

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Type declarations - lengths of reals, integers, and complexes.
!-------------------------------------------------------------------------------
      USE stel_kinds

!-------------------------------------------------------------------------------
!  Frequently used mathematical constants, lots of extra precision.
!-------------------------------------------------------------------------------
      USE stel_constants
      
!-------------------------------------------------------------------------------
!  Equilibrium Derived Types
!-------------------------------------------------------------------------------
      USE eq_T
      
!-------------------------------------------------------------------------------
!  V3 Utilities
!-------------------------------------------------------------------------------
      USE v3_utilities
      
!-------------------------------------------------------------------------------
!  V3 Global
!-------------------------------------------------------------------------------
      USE v3f_global, ONLY: l_zero_xcdot
!  JDH 2008-08-04
!     The USE v3f_global is a KLUDGE. This is an UGLY way to get control 
!     information to the eq_interface. Use for testing, and try to GET RID OF IT

      IMPLICIT NONE
      
!*******************************************************************************
! SECTION II.  INTERFACE BLOCKS
!*******************************************************************************
!SPH  Interface to external runvmec routine
      INTERFACE
         SUBROUTINE runvmec(ictrl_array, input_file, lscreen,                  &
     &                      reset_file_name)
         INTEGER, INTENT(inout), TARGET :: ictrl_array(5)
         LOGICAL, INTENT(in)            :: lscreen
         CHARACTER(LEN=*), INTENT(in)   :: input_file
         CHARACTER(LEN=*), OPTIONAL     :: reset_file_name
         END SUBROUTINE runvmec

      END INTERFACE

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      CONTAINS
          
!*******************************************************************************
! SECTION III.  EQUILIBRIUM INITIALIZATION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Initialize the equilibrium solver - from a file
!
!    Reads information in from the file 
!    Initializes the equilibrium solver
!    Constructs the equilibrium structures
!
!    For now, specific to the VMEC code.
!-------------------------------------------------------------------------------
      SUBROUTINE eq_init_file(filename,state)

!-------------------------------------------------------------------------------
!  USE statements
!-------------------------------------------------------------------------------

! VMEC input variables
      USE vmec_input, ONLY: lasym, lfreeb, lforbal, lrfp,                      &
     &   mgrid_file, nfp,                                                      &
     &   ncurr, nvacskip, mpol, ntor, ntheta, nzeta, ns_array,                 &
     &   pcurr_type, piota_type, pmass_type,                                   &
     &   niter, nstep, delt, ftol_array, gamma, tcon0,                         &
     &   ac, ac_aux_s, ac_aux_f, ai, ai_aux_s, ai_aux_f,                       &
     &   am, am_aux_s, am_aux_f, bloat, rbc, zbs, zbc, rbs,                    &     
     &   extcur, curtor, phiedge, pres_scale, l_v3fit

! VMEC dimensions
!  10-21-04 Check with Steve Hirshman, are there any gotchas with this
      USE vmec_dim, ONLY: mns

! VMEC flags for communicating with runvmec
      USE vmec_params, ONLY: restart_flag, readin_flag, timestep_flag,         &
     &   more_iter_flag, norm_term_flag

! VMEC fsq - for monitoring convergence
      USE vmec_main, ONLY: fsqr, fsqz, fsql

! VMEC xstuff - state arrays
      USE xstuff, ONLY: xc_vmec_1 => xc, xcdot_vmec_1 => xcdot

! VMEC nextcur - number of external currents
      USE mgrid_mod, ONLY: nextcur

!  VMEC_HISTORY
      USE vmec_history

!-------------------------------------------------------------------------------
!  Declarations
!-------------------------------------------------------------------------------
!  Declare Arguments 
      CHARACTER(len=*), INTENT(inout) :: filename
      TYPE (eq_state), INTENT (inout) :: state

!  Declare local variables
      TYPE (eq_param_fix) :: fixp
      TYPE (eq_param_var) :: varp
      CHARACTER(len=*), PARAMETER :: code = "VMEC2000"
!  Should get version from the VMEC code itself, rather than hard-wire it here
      CHARACTER(len=*), PARAMETER :: version = "8.45"
      CHARACTER(len=30) :: s_id = " "
      CHARACTER(len=80) :: l_id = " "
      INTEGER :: ns_index
      REAL(rprec) :: fsq_max

      INTEGER :: ictrl_array(5), ier_flag
      LOGICAL, PARAMETER :: lscreen     = .false.
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'eq_init_file: '

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------

!  Initialize VMEC: read in the file filename to load module vmec_input variables
!  AND initialize xc_vmec state array. Flags = reset(1) + initialize(2) + 
!  timestep(4) = 7.

      ictrl_array = 0
      ictrl_array(1) = restart_flag + readin_flag + timestep_flag
      ictrl_array(3) = 1           ! Run one step to initialize xc and xcdot

!  JDH 11-30-04. For now, limit VMEC to use only the first value in the 
!  ns_array. (i.e., no multigrid, just a single grid: ns_index = 1)
      ns_index = 1
      ictrl_array(4) = ns_index        

!  Set the l_v3fit logical
      l_v3fit = .true.

!  Set vmec_history printing off, so that VMEC doesn't control
!  Set vmec history integers - indicates where called from
      CALL vmec_history_print_flag_off
      CALL eq_history_set(i1 = 0,i2 = 0) 

      CALL runvmec(ictrl_array, filename, lscreen)
      ier_flag = ictrl_array(2)
         
!  Check for abnormal terminations from runvmec
      IF ((ier_flag .ne. norm_term_flag) .and.                                 &
     &    (ier_flag .ne. more_iter_flag)) THEN
             CALL err_fatal(sub_name // 'Bad ier_flag', int = ier_flag)
      ENDIF

!      CALL assert_eq(0,ier_flag,sub_name // 'Reading vmec input data')
!      IF (ier_flag .ne. 0) STOP 'Error reading in vmec input data!'

      fsq_max= MAX(fsqr,fsqz,fsql)
      WRITE(*,*) 'In eq_init_file, fsq_max = ', fsq_max

!  Construct the structures
      CALL eq_param_fix_construct(fixp,                                        &
     &   lasym, lfreeb, lforbal, lrfp, mgrid_file, nfp, ncurr,                 &
     &   nvacskip, mpol, ntor, ntheta, nzeta, mns, ns_array, nextcur,          &
     &   pcurr_type, piota_type, pmass_type,                                   &
     &   niter, nstep, delt, ftol_array, gamma, tcon0)
     
      CALL eq_param_var_construct(varp,ns_index,                               &
     &   ac, ac_aux_s, ac_aux_f, ai, ai_aux_s, ai_aux_f,                       &
     &   am, am_aux_s, am_aux_f, bloat, rbc, zbs, zbc, rbs,                    &     
     &   extcur, curtor, phiedge, pres_scale)
     
      CALL eq_state_construct(state,code,version,s_id,l_id,                    &
     &   fixp,varp,xc_vmec_1,xcdot_vmec_1,' ',fsq_max)

      RETURN
      END SUBROUTINE eq_init_file

!-------------------------------------------------------------------------------
!  Initialize the equilibrium solver - from equilibrium structures
!
!    Reads information in from the structures
!    Initializes the equilibrium solver
!
!    For now, specific to the VMEC code.
!-------------------------------------------------------------------------------
      SUBROUTINE eq_init_structure(state)

! Note to SH. I think you'll want some USE statements here, to get access to the
! VMEC variables

!  Declare Arguments 
      TYPE (eq_state), INTENT (inout) :: state

!  Declare local variables

!  Start of executable code

!  Initialize VMEC, using the information in fixp, varp, and state

      RETURN
      END SUBROUTINE eq_init_structure

!*******************************************************************************
! SECTION IV.  EQUILIBRIUM STEP SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Step the equilibrium solver
!    Takes as input the number of iterations numiter and the eq_state state
!    The state that is iterated is the  * internal * state of the equilibrium solver
!    The eq_state argument state contains the state on output.
!    The equilibrium solver MUST have been initialized already.
!
!    For now, specific to the VMEC code.
!-------------------------------------------------------------------------------

      SUBROUTINE eq_step(numiter,state,lconverged,iter_es,                     &
     &   ier_flag_vmec,vmec_flag_name,delta_iter2)

!  JDH 2009-09-10. Changed err_fatal to err_warn for various VMEC flags
!  JDH 2009-07-16. Revised and cleaned up coding

!-------------------------------------------------------------------------------
!  USE statements
!-------------------------------------------------------------------------------

! VMEC xstuff - state arrays
      USE xstuff, ONLY: xc_vmec_2 => xc, xcdot_vmec_2 => xcdot

! VMEC flags for communicating with runvmec
      USE vmec_params

! VMEC fsq - for monitoring convergence
!      iter2 - vmec internal iteration counter
!      SPH013108: add iter1 here
      USE vmec_main, ONLY: fsqr, fsqz, fsql, iter1, iter2

!  VMEC history - iteration counter
!      USE vmec_history, ONLY: iter_ha
!     ^  JDH 2010-08-03. Use iter2 from vmec_main for iter_es argument

!  IOU for runlog file
      USE v3f_global, ONLY: iou_runlog

!-------------------------------------------------------------------------------
!  Declarations
!-------------------------------------------------------------------------------
!  numiter         number of iterations for vmec to perform
!  state           Type eq_state, the equilibrium state to iterate
!  lconverged      Logical, whether or not equilibrium has converged
!  iter_es         integer, equilibrium solver internal iteration counter.
!                  (07-24-06 - use iter_ha from vmec_history)
!                  (2010-08-03 - use iter2 from vmec_main)
!  ier_flag_vmec   integer ier_flag returned from runvmec
!  vmec_flag_name  name of the ier_flag from VMEC
!  delta_iter2     Change in VMEC iter2 counter, due to runvmec call
!  
!  Declare Arguments 
      INTEGER, INTENT(inout) :: numiter
      TYPE (eq_state), INTENT (inout) :: state
      LOGICAL, INTENT(out), OPTIONAL :: lconverged
      INTEGER, INTENT(out), OPTIONAL :: iter_es
      INTEGER, INTENT(out), OPTIONAL :: ier_flag_vmec
      CHARACTER	(len=*), INTENT(inout), OPTIONAL :: vmec_flag_name
      INTEGER, INTENT(out), OPTIONAL :: delta_iter2     

!  Declare local variables
      INTEGER :: ictrl_array(5), ier_flag, ier1, iter2_before
      LOGICAL, PARAMETER :: lscreen     = .false.
      LOGICAL :: lconverged_local
      CHARACTER	(len=80) :: vmec_flag_name_local
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'eq_step: '

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  Check consistency of state variables and VMI variables
!-------------------------------------------------------------------------------
      CALL assert(ASSOCIATED(state % xc),ASSOCIATED(state % xcdot),            &
     &   sub_name // 'state variables not associated')
      CALL assert_eq(SIZE(state % xc),SIZE(xc_vmec_2),                         &
     &   sub_name // 'xc and state % xc not the same size')
      CALL assert_eq(SIZE(state % xcdot),SIZE(xcdot_vmec_2),                   &
     &   sub_name // 'xcdot and state % xcdot not the same size')

!  Internal VMEC state is unchanged before call to runvmec.
 
!-------------------------------------------------------------------------------
!  Setup call to runvmec
!-------------------------------------------------------------------------------
!  Use ns_index from the state
      ictrl_array = 0
      ictrl_array(1) = timestep_flag + output_flag + reset_jacdt_flag
      ictrl_array(3) = numiter
      ictrl_array(4) = state % varp % ns_index
      iter2_before = iter2
      iter1        = iter2
      
      IF (l_zero_xcdot) CALL eq_change_vmi_zero_xcdot
!  l_zero_xcdot gets here via a USE v3f_global
!  For testing 2008-08-04

      CALL runvmec(ictrl_array, ' ', lscreen)

!-------------------------------------------------------------------------------
!  Set error flag. (Flags declared in vmec_params.f)
!-------------------------------------------------------------------------------

      ier_flag = ictrl_array(2)
      SELECT CASE (ier_flag)
      CASE (norm_term_flag)
         vmec_flag_name_local = 'norm_term_flag'       ! = 0
      CASE (bad_jacobian_flag)
         vmec_flag_name_local = 'bad_jacobian_flag'    ! = 1
      CASE (more_iter_flag)
         vmec_flag_name_local = 'more_iter_flag'       ! = 2
      CASE (jac75_flag)
         vmec_flag_name_local = 'jac75_flag'           ! = 4
      CASE (input_error_flag)
         vmec_flag_name_local = 'input_error_flag'     ! = 5
      CASE (phiedge_error_flag)
         vmec_flag_name_local = 'phiedge_error_flag'   ! = 7
      CASE (ns_error_flag)
         vmec_flag_name_local = 'ns_error_flag'        ! = 8
      CASE (misc_error_flag)
         vmec_flag_name_local = 'misc_error_flag'      ! = 9
      CASE (successful_term_flag)
         vmec_flag_name_local = 'successful_term_flag' ! = 11
      CASE DEFAULT
         vmec_flag_name_local = 'default'
      END SELECT

      IF (PRESENT(vmec_flag_name)) vmec_flag_name =                            &
     &      vmec_flag_name_local

      IF (PRESENT(ier_flag_vmec)) ier_flag_vmec = ier_flag
         
!  Check for abnormal terminations from runvmec. Only warn
      IF ((ier_flag .ne. norm_term_flag) .and.                                 &
     &       (ier_flag .ne. successful_term_flag) .and.                        &
     &       (ier_flag .ne. more_iter_flag)) THEN
          CALL err_warn(sub_name // 'Bad ier_flag ' //                         &
     &            vmec_flag_name_local, int = ier_flag)
      ENDIF
      
!-------------------------------------------------------------------------------
!  Copy the VMEC variables into the state variables
!-------------------------------------------------------------------------------
!  First, need to check to see if sizes have changed
!  (They would have changed if ns_index was increased)
      IF (SIZE(state % xc) .ne. SIZE(xc_vmec_2)) THEN
         DEALLOCATE(state % xc,STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'dealloc xc')
         ALLOCATE(state % xc(1:SIZE(xc_vmec_2)),STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'alloc xc')
      ENDIF
      IF (SIZE(state % xcdot) .ne. SIZE(xcdot_vmec_2)) THEN
         DEALLOCATE(state % xcdot,STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'dealloc xcdot')
         ALLOCATE(state % xcdot(1:SIZE(xcdot_vmec_2)),STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'alloc xcdot')
      ENDIF

      state % xc = xc_vmec_2
      state % xcdot = xcdot_vmec_2
      state % fsq_max = MAX(fsqr,fsqz,fsql)
!  Need to get correct wout filename.
      state % wout_filename = 'wout_.nc'

!-------------------------------------------------------------------------------
!  Logic to see if have converged
!-------------------------------------------------------------------------------
      lconverged_local = .false.
      IF (state%fsq_max .le.                                                   &
     &   state%fixp%ftol_array(state%varp%ns_index)) THEN
         lconverged_local = .true.
      ENDIF
      IF (PRESENT(lconverged)) lconverged = lconverged_local
      
!-------------------------------------------------------------------------------
!  Diagnostic Printout
!-------------------------------------------------------------------------------
      WRITE(*,1000) iter2, iter2 - iter2_before,                               &
     &   vmec_flag_name_local, ier_flag, lconverged_local
      WRITE(iou_runlog,1000) iter2, iter2 - iter2_before,                      &
     &   vmec_flag_name_local, ier_flag, lconverged_local
1000  FORMAT('   eq_step:', i7,1x,i7,1x,a20,1x,i3,1x,l1)

!  Print out more information when have not converged
      IF (.not. lconverged_local) THEN
         WRITE(*,*) ' lconverged_local = ', lconverged_local,                  &
     &      ' in eq_step'
         WRITE(*,*) ' state%fsq_max = ',state%fsq_max
         WRITE(*,*) ' state%fixp%ftol_array(state%varp%ns_index) = ',          &
     &      state%fixp%ftol_array(state%varp%ns_index)
         WRITE(*,*) ' vmec_flag_name_local = ', vmec_flag_name_local
         WRITE(*,*) ' ier_flag = ', ier_flag
         WRITE(iou_runlog,*) ' lconverged_local = ', lconverged_local,         &
     &      ' in eq_step'
         WRITE(iou_runlog,*) ' state%fsq_max = ',state%fsq_max
         WRITE(iou_runlog,*)                                                   &
     &      ' state%fixp%ftol_array(state%varp%ns_index) = ',                  &
     &      state%fixp%ftol_array(state%varp%ns_index)
         WRITE(iou_runlog,*) ' vmec_flag_name_local = ',                       &
     &      vmec_flag_name_local
         WRITE(iou_runlog,*) ' ier_flag = ', ier_flag
      ENDIF

!-------------------------------------------------------------------------------
!  Clean up loose ends
!-------------------------------------------------------------------------------
         
      IF (PRESENT(delta_iter2)) delta_iter2 = iter2 - iter2_before

!  Logicals for auxiliarys defined in eq_state set to false
      state % l_def_aux1 = .false.
      state % l_def_aux2 = .false.

!  Cumulative Internal iteration counter
!      iter_es = iter_ha JDH 2010-08-03
      iter_es = iter2
      
      RETURN
      
      END SUBROUTINE eq_step
!*******************************************************************************
! SECTION V.  EQUILIBRIUM CHANGE SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Change the internal state of the equilibrium solver
!
!    For now, specific to the VMEC code.
!-------------------------------------------------------------------------------

      SUBROUTINE eq_change_vmi_varp(state)

!-------------------------------------------------------------------------------
!  USE statements
!-------------------------------------------------------------------------------

! VMEC xstuff - state arrays
      USE xstuff, ONLY: xc_vmec_3 => xc, xcdot_vmec_3 => xcdot,                &
     &   xstore_vmec_3 => xstore

! VMEC input variables. mgrid_file, nfp to reread for extcur change
      USE vmec_input, ONLY:                                                    &
     &   ac, ac_aux_s, ac_aux_f, ai, ai_aux_s, ai_aux_f,                       &
     &   am, am_aux_s, am_aux_f, bloat, rbc, zbs, zbc, rbs,                    &     
     &   extcur, curtor, phiedge, pres_scale,                                  &
     &   mgrid_file, nfp

! VMEC  - icurv array, for checking. currv - for curtor
!   irzloff - for call to profil3d
! VMEC ivac, fsqr, fsqz for ensuring that extcur changes are absorbed
      USE vmec_main, ONLY: icurv, currv, irzloff, ivac, fsqr, fsqz

!  mu0 travels through vparams, and vmec_main to get to readin
      USE stel_constants, ONLY: mu0

!  nv is in vacmod0, used as argument to read_mgrid
      USE vacmod0, ONLY: nv

!  subroutine read_mgrid is in the mgrid_mod module
!  Need mgrid_path_old to get around a test in read_mgrid
      USE mgrid_mod, ONLY: read_mgrid, mgrid_path_old

!-------------------------------------------------------------------------------
!  Declarations
!-------------------------------------------------------------------------------
!  state       Type eq_state. Change the internal state of vmec to match the varp
!  
!  Declare Arguments 
      TYPE (eq_state), INTENT (inout) :: state

!  Declare local variables
!  lextcur_diff  logical to determine if extcur has changed.
!  lscreen_arg   argument for read_mgrid
!  ier_flag_arg  argument for read_mgrid
      LOGICAL, PARAMETER :: l_reset_xc     = .false.
      LOGICAL            :: lextcur_diff, lscreen_arg
      INTEGER            :: ier_flag_arg
      REAL(rprec)        :: t1, t2, t3
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'eq_change_vmi_varp: '

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Copy the variable parameters into the VMEC variables
!  Sizes were checked on construction, so don't bother here
!-------------------------------------------------------------------------------
 
!  JDH 2010-12-06 Add more _var variables to those copied to VMI

      ac = state % varp % ac
      ac_aux_s = state % varp % ac_aux_s
      ac_aux_f = state % varp % ac_aux_f
      ai = state % varp % ai
      ai_aux_s = state % varp % ai_aux_s
      ai_aux_f = state % varp % ai_aux_f
      am = state % varp % am
      am_aux_s = state % varp % am_aux_s
      am_aux_f = state % varp % am_aux_f
      bloat = state % varp % bloat
      curtor = state % varp % curtor
      pres_scale = state % varp % pres_scale

!  Extra logic and coding for change of phiedge
!  JDH 2010-11-24. Added when VMEC lrfp coding added
!  SPH email --- The reason for the change is that in the new code, instead
!  of phip*(1+lambda) we have phip + lambda, and lambda needs to be scaled
!  by this factor to keep the initial lambda force small.---
!  JDH 2010-11-29. Added coding to als0 rescale xstore
      IF (phiedge .ne. state % varp % phiedge) THEN
         IF (phiedge .ne. 0._rprec) THEN
            t1 = state % varp % phiedge / phiedge
            xc_vmec_3(1+2*irzloff:3*irzloff) = t1 *                            &
     &         xc_vmec_3(1+2*irzloff:3*irzloff)
            xstore_vmec_3(1+2*irzloff:3*irzloff) = t1 *                        &
     &         xstore_vmec_3(1+2*irzloff:3*irzloff)
         END IF
         phiedge = state % varp % phiedge
      ENDIF

!  See if extcur has changed, set lextcur_diff
      t1 = dot_product(extcur - state % varp % extcur,                         &
     &   extcur - state % varp % extcur)
      t2 = dot_product(extcur,extcur)
      t3 = dot_product(state % varp % extcur,state % varp % extcur)
      IF (t1 .eq. 0) THEN
         lextcur_diff = .FALSE.
      ELSE
!  JDH 2007-12-21 : Change test from 1.e-9 to 1.e-20
         lextcur_diff = (t1 / (t2 + t3)) .ge. 1.e-20
      ENDIF
      
      IF (lextcur_diff) THEN
         extcur = state % varp % extcur
      ENDIF

!-------------------------------------------------------------------------------
!  Make sure vmec "absorbs" these changes.
!-------------------------------------------------------------------------------
!  Execute relevant statements from readin
!  (For curtor, gets put into currv)
      currv = mu0*curtor              !Convert to Internal units

!  More 'absorption' in profil1d
      CALL profil1d (xc_vmec_3, xcdot_vmec_3, l_reset_xc)

!  Change to phiedge needs profil3d called.
!  Mimic call from initialize_radial
      CALL profil3d(xc_vmec_3(1), xc_vmec_3(1+irzloff), l_reset_xc)

!  Have to reread the mgrid file if the external currents have been changed
!  Have to change mgrid_path_old (module variable in mgrid_mod) to avoid a trap
!  in read_mgrid, that detects if the mgrid file has been read already.
      IF (lextcur_diff) THEN
         mgrid_path_old = ' '
         lscreen_arg = .FALSE.
         ier_flag_arg = 0
         CALL read_mgrid(mgrid_file,extcur,nv, nfp,lscreen_arg,                &
     &      ier_flag_arg)
         CALL assert_eq(0,ier_flag_arg,sub_name // 'bad read_mgrid')
      ENDIF

!  Coding added 2007-12-21 JDH
!  Also have to ensure that changed current is used in the first VMEC iteration
!  See VMEC coding in subroutine funct3d.f
!  Want to ensure that coding inside IF statements at lines
!  154 - IF (lfreeb .and. iter2.gt.1 .and. iequi.eq.0) THEN
!    and
!  157  - IF (ivac .ge. 0) THEN
!     and
!  166  - IF (ivac .le. 2) ivacskip = 0
!  gets executed.
!  JDH 2008-01-10
!  Try to prevent execution of 
!  155  - IF ((fsqr + fsqz) .le. 1.e-3_dp) ivac = ivac + 1
!  JDH 2008-01-14 - Did not work without fsqr=fsqz=1.
      IF (lextcur_diff) THEN
         ivac = 1
         fsqr = one
         fsqz = one
       ENDIF

      RETURN
      
      END SUBROUTINE eq_change_vmi_varp
      
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      SUBROUTINE eq_change_vmi_zero_xcdot
!  Subroutine to set the VMEC xcdot array to zero

!-------------------------------------------------------------------------------
!  USE statements
!-------------------------------------------------------------------------------

! VMEC xstuff - state arrays
      USE xstuff, ONLY: xcdot_vmec_4 => xcdot

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------

      xcdot_vmec_4 = zero
      
      END SUBROUTINE eq_change_vmi_zero_xcdot

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      SUBROUTINE eq_change_vmi_cp_xc(state)
!  Subroutine to copy the xc array from an eq_state to the VMEC internal state

!-------------------------------------------------------------------------------
!  USE statements
!-------------------------------------------------------------------------------

! VMEC xstuff - state arrays
      USE xstuff, ONLY: xc_vmec_5 => xc

!-------------------------------------------------------------------------------
!  Declarations
!-------------------------------------------------------------------------------
!  state       Type eq_state.
!  
!  Declare Arguments 
      TYPE (eq_state), INTENT (in) :: state

!  Declare local variables
!  n_xc_vmec      size of the VMEC internal xc array
!  n_xc_state     size of the state xc array
      INTEGER            :: n_xc_vmec, n_xc_state
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'eq_change_cp_xc: '

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------

! Get the array sizes
      n_xc_vmec = SIZE(xc_vmec_5)
      n_xc_state = SIZE(state % xc)
      
      IF (n_xc_vmec .eq. n_xc_state) THEN
         xc_vmec_5 = state % xc
      ELSE
         WRITE(*,*) 'WARNING from ', sub_name
         WRITE(*,*) ' xc arrays of different sizes'
         WRITE(*,*) ' n_xc_vmec, n_xc_state = ', n_xc_vmec, n_xc_state
         WRITE(*,*) 'END WARNING'
      ENDIF
            
      END SUBROUTINE eq_change_vmi_cp_xc

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      SUBROUTINE eq_change_vmi_precon2d_off
!  Subroutine to turn off 2d preconditioning in the VMEC internal state

!-------------------------------------------------------------------------------
!  USE statements
!-------------------------------------------------------------------------------

! VMEC xstuff 
       USE precon2d, ONLY: ictrl_prec2d

!-------------------------------------------------------------------------------
!  Declarations
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------

! Get the array sizes
      ictrl_prec2d = 0
            
      END SUBROUTINE eq_change_vmi_precon2d_off

!*******************************************************************************
! SECTION VI.  EQUILIBRIUM AUXILIARY CALCULATION SUBROUTINES
!*******************************************************************************


!*******************************************************************************
! SECTION VI.  EQUILIBRIUM WRITE SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Print out history information
!    Module vmec_history, in VMEC2000/Sources/TimeStep store some FTOL 
!    information after each iteration.
!    This subroutine calls the vmec_hisotry subroutine that prints the
!    information out.
!
!    Specific to the VMEC code.
!-------------------------------------------------------------------------------
      SUBROUTINE eq_history_print

      USE vmec_history

!  Start of executable code

!  write out the vmec_history
      CALL vmec_history_print_flag_on
      CALL vmec_history_print

      RETURN
      END SUBROUTINE eq_history_print

!-------------------------------------------------------------------------------
!  Set integers for vmec_history
!    Use i1 for reconstruction iteration number
!    Use i2 to specify which reconstruction parameter, when calculating jacobian
!
!    Specific to the VMEC code.
!-------------------------------------------------------------------------------
      SUBROUTINE eq_history_set(i1,i2)

      USE vmec_history

!  Declare Arguments 
      INTEGER, OPTIONAL :: i1, i2

!  Start of executable code

      IF (PRESENT(i1)) THEN
         IF (PRESENT(i2)) THEN
            CALL vmec_history_set(i1,i2)  !  both i1 and i2
         ELSE
            CALL vmec_history_set(i1)     !  i1 only
         ENDIF
      ELSE
         IF (PRESENT(i2)) THEN
            CALL vmec_history_set(i2=i2)  !  i2 only
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE eq_history_set

!-------------------------------------------------------------------------------
!  Get integers for vmec_history
!
!    Specific to the VMEC code.
!-------------------------------------------------------------------------------
      SUBROUTINE eq_history_get(i1,i2)

      USE vmec_history

!  Declare Arguments 
      INTEGER :: i1, i2

!  Start of executable code
      CALL vmec_history_get(i1,i2)

      RETURN
      END SUBROUTINE eq_history_get

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      SUBROUTINE eq_interface_nli_write(input_file_end,istat)
!  Write out a VMEC namelist-input file
      USE vmec_input
      USE mgrid_mod, ONLY:  curlabel, nextcur
      USE safe_open_mod
      USE write_array_generic
      USE vacfield_mod, ONLY: write_invac
      USE date_and_computer
      USE vmec_params, ONLY: version_
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(out) :: istat
      CHARACTER(LEN=*) :: input_file_end
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(1) :: ins
      INTEGER :: iftol, m, n, i, k, ns0, ns1, iunit,                           &
     &    nextcur_max, iout=23
      CHARACTER(LEN=10) :: date0, time0, zone0
      CHARACTER(LEN=50) :: dateloc, Version
      INTEGER :: imon
      CHARACTER(LEN=80), PARAMETER ::                                          &
     &   banner =                                                              &
     &'! THIS IS WRITTEN BY V3FIT USING VMEC2000, VERSION '
!-----------------------------------------------
      CALL GetComputerInfo
      CALL DATE_AND_TIME(date0,time0,zone0)
      READ(date0(5:6),'(i2)')imon
      WRITE(dateloc,100)months(imon),date0(7:8),date0(1:4),                    &
     &  time0(1:2),time0(3:4),time0(5:6)
 100  FORMAT('! DATE = ',a3,' ',a2,',',a4,' ',' TIME = ',2(a2,':'),a2)

      iunit = iout+2
      CALL safe_open(iunit, istat, input_file_end,                             & 
     &	'replace', 'formatted')
      IF (istat .ne. 0) THEN
         istat = -3
         RETURN
      END IF

      ins = MAXLOC(ns_array)
      ns0 = ns_array(ins(1))

      iftol =0 
      DO WHILE(ftol_array(iftol+1).ne.zero .and. iftol.lt.100)
         iftol = iftol + 1
      END DO

      IF (nextcur .eq. 0) THEN
         DO nextcur = SIZE(extcur), 1, -1
            IF (extcur(nextcur) .ne. zero) EXIT
         END DO
         nextcur = MAX(0, nextcur)
      END IF

!
!     Start writing out formatted VMEC data with &INDATA heading
!     so that it can be read in under the indata NAMELIST
!
      WRITE (iunit, '(a7)') '&INDATA'
      WRITE (iunit, '(2x,3a)') "MGRID_FILE = '", TRIM(mgrid_file), "'"
      WRITE (iunit, 1007) 'LFREEB = ', lfreeb
      WRITE (iunit, 1007) 'LASYM = ', lasym
      WRITE (iunit, 1007) 'LFORBAL = ', lforbal
      WRITE (iunit, '(2x,a,1p,e10.2)') 'DELT = ', delt
      WRITE (iunit, '(2x,a,1p,e10.2)') 'TCON0 = ', tcon0
      WRITE (iunit,'(2x,3a)') "PRECON_TYPE = '", TRIM(precon_type),"'" 
      WRITE (iunit,'(2x,a,1p,e14.6)') "PREC2D_THRESHOLD = ",                   &
     &                                prec2d_threshold
      WRITE (iunit, '(2x,a,i4)') 'NFP = ', nfp
      WRITE (iunit, '(2x,a,i4)') 'NCURR = ', ncurr
      WRITE (iunit, '(2(2x,a,i4))') 'MPOL = ', mpol, 'NTOR = ', ntor
      IF (ntheta.gt.0) WRITE (iunit, '(2x,a,i4)') 'NTHETA = ',ntheta
      IF (nzeta.gt.0)  WRITE (iunit, '(2x,a,i4)') 'NZETA  = ',nzeta
      IF (mfilter_fbdy.gt.0)                                                   &
     &    WRITE (iunit, '(2x,a,i4)') 'mfilter_fbdy  = ',mfilter_fbdy
      IF (nfilter_fbdy.gt.0)                                                   & 
     &    WRITE (iunit, '(2x,a,i4)') 'nfilter_fbdy  = ',nfilter_fbdy

      WRITE (iunit, '(2x,a)') 'NS_ARRAY = '
      WRITE (iunit, 980) (ns_array(i),i=1,ins(1))
!  JDH 2010-07-20 i -> i6, single line below, to make gfortran happy
      WRITE (iunit, '(2x,a,i6)') 'NITER = ', niter
      WRITE (iunit, '(2x,a,i4)') 'NSTEP = ', nstep
      WRITE (iunit, '(2x,a,i4)') 'NVACSKIP = ', nvacskip
      WRITE (iunit, 994) gamma
      WRITE (iunit, '(2x,a)') 'FTOL_ARRAY ='
      WRITE (iunit, 995, err=2000) (ftol_array(i),i=1,iftol )
      WRITE (iunit, 996) phiedge
      WRITE (iunit, 997) pres_scale
      WRITE (iunit, 998) curtor
      WRITE (iunit, '(2x,a,1p,e14.6)') 'BLOAT = ', bloat
      WRITE (iunit, '(2x,3a)') "PMASS_TYPE = '", pmass_type, "'"      
      WRITE (iunit, '(2x,a,/,2x,3(1pe22.14))')'AM = ', am
      WRITE (iunit, '(2x,3a)') "PCURR_TYPE = '", pcurr_type, "'"      
      WRITE (iunit, '(2x,a,/,2x,3(1pe22.14))')'AC = ', ac
      WRITE (iunit, '(2x,a,/,2x,3(1pe22.14))')'AI = ', ai
      IF (nextcur .gt. 0) then
        DO i=1,nextcur 
          WRITE(iunit,'(a,i2.2,a,1p,e22.14,2a)')                               &
     &		"  EXTCUR(",i,") = ", extcur(i),"	!",TRIM(curlabel(i))
!     1    CALL write_array(iunit, 'EXTCUR', extcur, nextcur)
        ENDDO
      ENDIF

      CALL write_rbzb (iunit, istat)
      IF (istat .ne. 0) GOTO 2000

  980 FORMAT(2x,16i5)
  981 FORMAT(2x,a,10(i5))
  990 FORMAT(2x,a,6(1p,e22.14))
  991 FORMAT(1x,a,5(1p,e122.14))
  993 FORMAT(2(2x,a,1p,e22.14))
  994 FORMAT(2x,'GAMMA = ',1p,e14.6)
  995 FORMAT(2x,1p,4e14.6)
  996 FORMAT(2x,'PHIEDGE = ',1p,e22.14)
  997 FORMAT(2x,'PRES_SCALE = ',1p,e22.14)
  998 FORMAT(2x,'CURTOR = ',1p,e22.14)
 1000 FORMAT(5(1p,e22.14))
 1007 FORMAT(2x,a,L1)
 1008 FORMAT(2x,a,L1,2x,a,i2)
 1009 FORMAT("/",a)
      WRITE (iunit, 1009) 
      CALL GetComputerInfo
      Version = TRIM(ADJUSTL(version_))
      WRITE(iunit,'(a,x,a,/,a,x,a,/,a,/,a,x,a)')                               &
     &  TRIM(banner),TRIM(Version),                                            &
     &     '! COMPUTER: ', TRIM(computer),                                     &
     &  TRIM(dateloc)
      CLOSE (iunit)
      RETURN
     
!
!     Handle errors here
!
 2000 istat = -5

      END  SUBROUTINE eq_interface_nli_write

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      SUBROUTINE write_rbzb(iunit, istat)
      USE vmec_input
      USE write_array_generic
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: iunit, istat, m, n
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      CHARACTER(LEN=*), DIMENSION(7), PARAMETER :: outfmt =                    &
     &   (/ "(2(a6,i1,a1,i1,a4,1p,e22.14,3x))",                                &
     &      "(2(a6,i1,a1,i2,a4,1p,e21.14,3x))",                                &
     &      "(2(a6,i2,a1,i1,a4,1p,e21.14,3x))",                                &
     &      "(2(a6,i2,a1,i2,a4,1p,e21.14,3x))",                                &
     &      "(2(a6,i3,a1,i1,a4,1p,e21.14,3x))",                                &
     &      "(2(a6,i3,a1,i2,a4,1p,e21.14,3x))",                                &
     &      "(2(a6,i ,a1,i ,a4,1p,e21.14,3x))" /)
      CHARACTER(len=LEN(outfmt(1))) :: outcfmt
      CALL write_array(iunit, 'RAXIS_CC', raxis_cc, ntor+1, low_index=0)
      CALL write_array(iunit, 'ZAXIS_CS', zaxis_cs, ntor+1, low_index=0)
      CALL write_array(iunit, 'RAXIS_CS', raxis_cs, ntor+1, low_index=0)
      CALL write_array(iunit, 'ZAXIS_CC', zaxis_cc, ntor+1, low_index=0)
      DO m = 0, mpol-1
         DO n = -ntor, ntor
            IF ((rbc(n,m).ne.zero) .or. (zbs(n,m).ne.zero)) THEN
!
!     handle formatting for up to 2 digit n by 3 digit m.
!     while this is probably overkill, we at least have to handle up
!     thru 2 digit n by 2 digit m to handle known cases.
!
               outcfmt = outfmt(1)
               IF( m > 9) outcfmt = outfmt(2)
               IF(( n>-10 .and. n<0 ) .or. (n > 9 .and. n < 100)) THEN
                  outcfmt = outfmt(3)
                  IF( m > 9) outcfmt = outfmt(4)
               ELSE IF( n>-100 .and. n< -9 ) THEN
                  outcfmt = outfmt(5)
                  IF( m > 9) outcfmt = outfmt(6)
               ENDIF

               IF( n>= 100 .or. n <= -100 .or. m >=100) THEN
                  outcfmt = outfmt(7)
               ENDIF
             IF(lasym) then
               WRITE (iunit, outcfmt, err=2000)                                &
     &            '  RBC(', n, ',', m, ') = ', rbc(n,m),                       &
     &            '  ZBS(', n, ',', m, ') = ', zbs(n,m),                       &
     &            '  RBS(', n, ',', m, ') = ', rbs(n,m),                       &
     &            '  ZBC(', n, ',', m, ') = ', zbc(n,m)
             ELSE
               WRITE (iunit, outcfmt, err=2000)                                &
     &            '  RBC(', n, ',', m, ') = ', rbc(n,m),                       &
     &            '  ZBS(', n, ',', m, ') = ', zbs(n,m)
             ENDIF
            ENDIF
         END DO
      END DO

 100  FORMAT(a,(1p,4e22.14))

      istat = 0
      RETURN

 2000 istat = -5

      END SUBROUTINE write_rbzb
!*******************************************************************************
! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!  
!  JDH 09-21-2004
!     First version of module
!
!  JDH 11-02-2004
!     First version of eq_step
!
!  JDH 11-26-2004
!     Modified for revised eq_T module
!
!  JDH 11-28-2004
!     Added nextcur access (in module mgrid_mod). Set state logical variables 
!     in eq_step.
!
!  JDH 12-06-2004
!    ns_index issues, added specification of ictrl_array(4) to eq_step, Change
!    size of xc and xcdot if needed.
!
!  JDH 06-28-2006
!    Added iter_es argument to eq_step, equilibrium solver iteration counter
!
!  JDH 07-18-06
!    Added subroutine eq_history_print
!
!  JDH 07-20-06
!    Modified eq_step, so that don't copy varp items into vmec at start.
!    Add subroutine eq_change_vmi_varp
!
!  JDH 07-24-06
!    Modified eq_step a bit more: only step, don't copy from state to vmec
!    internal variables. Changes to eq_change_vmi_varp: change ac also.
!
!  JDH 11-28-2006
!    Changed eq_change_vmi_varp. Changes pres_scale and am (I hope)
!
!  JDH 12-27-2006
!    Changed reset_jdt_flag to reset_jacdt_flag, consistency with latest runvmec.f
!
!  JDH 12-29-2006
!    Added coding to eq_change_vmi_varp for phiedge and extcur
!
!  JDH 2007-12-21
!    Changed test value for extcur change detection (lextcur_diff)
!    from 1.e-9 to 1.e-20, in subroutine eq_change_vmi_varp
!
!  JDH 2008-01-14 (started 2007-12-21)
!    Additional changes so that extcur change "takes" in VMEC
!
!  JDH 2008-08-04
!    Added subroutines eq_change_vmi_zero_xcdot and eq_change_vmi_cp_xc(state)
!
!  JDH 2009-03-24
!    Changes from Joan Jimenez via Steve Hirshman. Peculiarity of compiler on
!    origin2000. Complaints about repeated name aliases in USE statements in
!    different module subroutines.
!
!  JDH 2009-05-15
!    Set the l_v3fit logical in eq_init_file. Add pcurr_type and pmass_type
!    variables to eq_interface_nli_write.
!
!  JDH 2009-08-27
!    Significant changes to eq_step. Actual coding changes written 2009-07-16
!    Revised and cleaned up code.
!
!  JDH 2009-11-16
!    Slight changes to coding
!
!  JDH 2010-11-24. 
!    Change to eq_change_vmi_varp, extra logic with VMEC lrfp coding. when
!    phiedge changes.
!
!  JDH 2010-12-06. 
!    Changed argument sequence for calls to eq_param_fix_construct and 
!    eq_param_var_construct, in eq_init_file. Corresponding changes in
!    eq_T.
!    Added more _var variables to be copied in eq_change_vmi_varp
          
      END MODULE eq_interface
