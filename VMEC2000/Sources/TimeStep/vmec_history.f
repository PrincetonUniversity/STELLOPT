!*******************************************************************************
!  File vmec_history.f
!  Contains module vmec_history

!*******************************************************************************
!  MODULE vmec_history
!    (history of a vmec run)
! SECTION I.    Variable Declarations
! SECTION II.   Subroutines
! SECTION III.  Comments - version history
!*******************************************************************************

      MODULE vmec_history

!*******************************************************************************
! SECTION I. Variable Declarations
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Type declarations - lengths of reals, integers, and complexes.
!  Frequently used mathematical constants, lots of extra precision.
!  Make type declarations and constants Private, so there are no conflicts.
!-------------------------------------------------------------------------------

      USE stel_kinds, only : rprec
      USE safe_open_mod
      USE xstuff

      IMPLICIT NONE

      PRIVATE rprec

!-------------------------------------------------------------------------------
!  Module Variables - Scalars
!-------------------------------------------------------------------------------
!  vmh_dim                Length of History Arrays
!  vmh_index              index to the history arrays, and call counter
!  vmh_save_i1   		  Integer stored value. Value set with subroutine
!                             vmec_history_set. Used by V3FIT for reconstruction
!                             iteration number
!  vmh_save_i2   		  Integer stored value. Value set with subroutine
!                             vmec_history_set. Used by V3FIT for jacobian 
!                             calculation loop (reconstruction parameter number)
!  vmh_print_flag         Logical to control printing
!    Print Flag Usage - Now (2010-08-12)
!      1) Initialized to .False.
!      2) vmec_history_print not called from VMEC
!      3) Flag changed to .TRUE. in V3FIT
!      4) vmec_history_print called at end of V3FIT run
!    Print Flag Usage - Future - If want VMEC run alone to also call vmec_history_print
!      1) Initialized to .TRUE.
!      2) vmec_history_print from VMEC
!      3) When run from V3FIT
!           a) have V3FIT change flag to .False. at start
!                 Will need to fix up eq_interface
!           b) change flag to .TRUE. at end, and call vmec_history_print
!  vmh_time_zero         Real variable to store the initial result of the second0 call

      INTEGER, PARAMETER :: vmh_dim = 100000
      INTEGER            :: vmh_index = 0
      INTEGER            :: vmh_save_i1 = - 1
      INTEGER            :: vmh_save_i2 = - 1
!      LOGICAL            :: vmh_print_flag = .TRUE.           
      LOGICAL            :: vmh_print_flag = .FALSE.
      REAL(rprec)        :: vmh_time_zero = 0
      PRIVATE vmh_dim, vmh_index, vmh_save_i1, vmh_save_i2,                    &
     &   vmh_print_flag

!-------------------------------------------------------------------------------
!  Module Variables - Integer arrays
!-------------------------------------------------------------------------------
!  vmh_iterc           VMEC's iterc 
!  vmh_iter2m1         VMEC's iter2 - iter1 
!  vmh_ns              VMEC's ns 
!  vmh_nvacskip        VMEC's nvacskip
!  vmh_ivac            VMEC's ivac
!  vmh_ictrl_prec2d    VMEC's ictrl_prec2d
!  vmh_i1              V3FIT, reconstruction iteration number
!  vmh_i2              V3FIT, jacobian calculation, reconstruction parameter #

      INTEGER, DIMENSION(vmh_dim) :: vmh_iterc, vmh_iter2m1, vmh_ns,           &
     &    vmh_nvacskip, vmh_ivac, vmh_ictrl_prec2d, vmh_i1, vmh_i2

!-------------------------------------------------------------------------------
!  Module Variables - Real arrays
!-------------------------------------------------------------------------------
! vmh_fsq-     		Convergence diagnostics
! vmh_time_step		VMEC delt0 values
! vmh_time          SYSTEM_TIME, via LIBSTELL's second0 function

      REAL(rprec), DIMENSION(vmh_dim) :: vmh_time_step, vmh_fsqr,              &
     &   vmh_fsqz, vmh_fsql, vmh_fedge, vmh_time

!*******************************************************************************
! SECTION II. Subroutines
!*******************************************************************************
      CONTAINS

!*******************************************************************************
!*******************************************************************************
!  vmec_history_store
!     subroutine to store local vmec and variables into the history arrays.
!     Should be called ONLY from subroutine eqsolve, right after iterc
!     is incremented.
!-------------------------------------------------------------------------------

      SUBROUTINE vmec_history_store(time_step)

!-------------------------------------------------------------------------------
!  USE statements
!-------------------------------------------------------------------------------

      USE vmec_main, ONLY: iter1, iter2, iterc, fsqr, fsqz, fsql,              &
     &   fedge, ivac
      USE vmec_dim, ONLY: ns
      USE vmec_input, ONLY: nvacskip
      USE precon2d, ONLY: ictrl_prec2d

!-------------------------------------------------------------------------------
!  ARGUMENT declaration
!-------------------------------------------------------------------------------
      REAL(rprec), INTENT(in) :: time_step
!  delt0 is a local variable in runvmec
!    runvmec: CALL eqsolve(.,delt0,...)
!       dummy argument of eqsolve is called delt0
!    eqsolve: CALL vmec_history_store(delt0)
!    eqsolve: CALL evolve(delt0,...)
!       dummy argument of evolve is time_step

!-------------------------------------------------------------------------------
!  Local Variables
!-------------------------------------------------------------------------------
      REAL(rprec) :: time_now

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      IF (vmh_index .eq. 0) THEN
         CALL second0(vmh_time_zero)
      ENDIF
      vmh_index = vmh_index + 1
      IF (vmh_index .le. vmh_dim) THEN
         vmh_iterc(vmh_index) = iterc
         vmh_iter2m1(vmh_index) = iter2 - iter1
         vmh_ns(vmh_index) = ns
         vmh_nvacskip(vmh_index) = nvacskip
         vmh_ivac(vmh_index) = ivac
         vmh_ictrl_prec2d(vmh_index) = ictrl_prec2d
         vmh_i1(vmh_index) = vmh_save_i1
         vmh_i2(vmh_index) = vmh_save_i2
         vmh_time_step(vmh_index) = time_step
         vmh_fsqr(vmh_index) = fsqr
         vmh_fsqz(vmh_index) = fsqz
         vmh_fsql(vmh_index) = fsql
         vmh_fedge(vmh_index) = fedge
         CALL second0(time_now)
         vmh_time(vmh_index) = time_now - vmh_time_zero
      ENDIF
      RETURN
      END SUBROUTINE vmec_history_store

!*******************************************************************************
!*******************************************************************************
!  vmec_history_print
!     subroutine to print out the accumulated history
!-------------------------------------------------------------------------------

      SUBROUTINE vmec_history_print

!-------------------------------------------------------------------------------
!  USE statements
!-------------------------------------------------------------------------------
      USE vmec_input, ONLY: input_extension

!-------------------------------------------------------------------------------
!  Local Variables
!-------------------------------------------------------------------------------
      INTEGER  :: vmh_iou = 73
      INTEGER  :: istat, i
      CHARACTER(LEN=120) :: vmh_history_file_name
      CHARACTER(LEN=80) :: vmh_format2 = 
     &   '(3(i5,1x),i4,1x,i3,1x,i5,1x,3(i3,1x),7(2x,es9.2))'
      CHARACTER(LEN=150) :: vmh_header 

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
      IF (.NOT. vmh_print_flag) RETURN

      vmh_history_file_name = TRIM('vmec_history.' // input_extension)
      CALL safe_open(vmh_iou,istat,TRIM(vmh_history_file_name),                &
     &   'replace','formatted',delim_in='none',record_in=150)
      IF (istat .ne. 0) THEN
         WRITE(*,*) 'In subroutine vmec_history_print: Error from'
         WRITE(*,*) 'call to safe_open. istat = ', istat
         STOP ' (source file vmec_history.f)'
      ENDIF
      
      WRITE(vmh_iou,*) 'History arrays are dimensioned ',vmh_dim
      WRITE(vmh_iou,*) 'Subroutine vmec_history_store was called ',            &
     &   vmh_index, ' times'
      WRITE(vmh_iou,*) 

      vmh_header = '    i iterc  2m1   ns nvac ivac ictrl_ i1 i2' //           &
     &  '    time_step  fsqr      fsqz       fsql      max(fsq)' //            &
     &  '    fedge      sys-time'
      WRITE(vmh_iou,*) TRIM(vmh_header)
      WRITE(vmh_iou,*) '                      skip      prec2d'

      DO i = 1,MIN(vmh_index,vmh_dim)
         WRITE(vmh_iou,vmh_format2)                                            &
     &      i, vmh_iterc(i), vmh_iter2m1(i), vmh_ns(i),                        &
     &      vmh_nvacskip(i), vmh_ivac(i),                                      &
     &      vmh_ictrl_prec2d(i), vmh_i1(i), vmh_i2(i),                         &
     &      vmh_time_step(i), vmh_fsqr(i), vmh_fsqz(i), vmh_fsql(i),           &
     &      MAX(vmh_fsqr(i),vmh_fsqz(i),vmh_fsql(i)),                          &
     &      vmh_fedge(i), vmh_time(i)
      END DO

      RETURN
      END SUBROUTINE vmec_history_print

!*******************************************************************************
!*******************************************************************************
!  vmec_history_set
!     Subroutine to set values of vmh_save_i1 and/or vmh_save_i2
!-------------------------------------------------------------------------------

      SUBROUTINE vmec_history_set(i1,i2)

!  Declare Arguments 
      INTEGER, OPTIONAL :: i1, i2

!  Start of executable code

      IF (PRESENT(i1)) vmh_save_i1 = i1
      IF (PRESENT(i2)) vmh_save_i2 = i2

      RETURN
      END SUBROUTINE vmec_history_set

!*******************************************************************************
!*******************************************************************************
!  vmec_history_get
!     Subroutine to get values of vmh_save_i1 and/or vmh_save_i2
!-------------------------------------------------------------------------------

      SUBROUTINE vmec_history_get(i1,i2)

!  Declare Arguments 
      INTEGER :: i1, i2

!  Start of executable code

      i1 = vmh_save_i1
      i2 = vmh_save_i2

      RETURN
      END SUBROUTINE vmec_history_get

!*******************************************************************************
!*******************************************************************************
!  vmec_history_print_flag_off
!     Subroutine to turn on the print flag.
!-------------------------------------------------------------------------------

      SUBROUTINE vmec_history_print_flag_off
      vmh_print_flag = .FALSE.
      RETURN
      END SUBROUTINE vmec_history_print_flag_off

!*******************************************************************************
!*******************************************************************************
!  vmec_history_print_flag_on
!     Subroutine to turn off the print flag.
!-------------------------------------------------------------------------------

      SUBROUTINE vmec_history_print_flag_on
      vmh_print_flag = .TRUE.
      RETURN
      END SUBROUTINE vmec_history_print_flag_on

      END MODULE vmec_history

!*******************************************************************************
! SECTION III. Comments - version history
!*******************************************************************************
!
!  JDH 07-12-2006. First version. 
!  Module to store history information about a vmec run
!
!  07-24-2006 JDH
!    Changed so that iter_ha is index of last assigned value in history arrays.
!    So, iter_ha increment happens BEFORE assignments in vmec_history_store.
!
!  2010-08-03 JDH
!    Significant Revisions, added more variables to store. Still needs V3F interface
!     Use variable prefix vmh_
!  2010-08-05 JDH  Add mechanism to turn off and on the printing
!  2010-08-10 JDH  Add time_step
!  2010-08-12 JDH  Add vmec_history_get, clean up headings, comments.
!  2011-02-08 JDH  Added vmh_fsq and vmh_fedge
!  2011-02-15 JDH  Removed vmh_fsq
!  2011-02-18 JDH  Added vmh_time - system time via LIBSTELL's second0
!  2011-02-19 JDH  Fixed up issues with vmh_time
!  2012-06-20 JDH  Changed iter2 to iterc for cumulative counter



