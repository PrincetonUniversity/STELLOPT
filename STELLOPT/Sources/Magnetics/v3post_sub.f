!*******************************************************************************
!  CODE V3POST
!
! SECTION I.   MAIN PROGRAM
! SECTION II.  PRIMARY SUBROUTINES (CALLED FROM MAIN)
! SECTION III.    AUXILIARY SUBROUTINES
! SECTION IV.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************

!*******************************************************************************
! SECTION I. MAIN PROGRAM
!*******************************************************************************
!CALL v3post_sub( TRIM(v3post_in), TRIM(extension), ierr)
      SUBROUTINE v3post_sub(v3post_in, extension, ierr)
      USE stel_kinds
      USE v3post_rfun
      IMPLICIT NONE
!----------------------------------------------------------------------
! D U M M Y Arguments Declarations
!----------------------------------------------------------------------
      CHARACTER (len=*) :: v3post_in, extension
      INTEGER :: ierr
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec) :: t1, t2
      REAL(rprec) :: ton, toff
      CHARACTER (len=120) :: plasfcn=' ', coilfcn=' '
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
      CALL second0(t1)
      ierr=0
!----------------------------------------------------------------------
!-- Open input file (v3post.ext), and read plasma equilibrium and v3post namelist--
!----------------------------------------------------------------------
      CALL get_input_v3p (TRIM(v3post_in), TRIM(extension),
     1   plasfcn, coilfcn)    !  reads the namelist
!----------------------------------------------------------------------
!-- Compute magnetic signals due to external coils                   --
!----------------------------------------------------------------------
! Note that second0 is a timing routine in LIBSTELL
      CALL second0(ton)
      IF (freeb) CALL get_coil_v3p (coilfcn)
      IF (.not.freeb) CALL no_coil (plasfcn)
      CALL second0(toff)
!----------------------------------------------------------------------
!-- Compute plasma contribution                                      --
!----------------------------------------------------------------------
      CALL get_plasma_v3p (plasfcn, ierr)
!----------------------------------------------------------------------
!-- Outputs                                                          --
!----------------------------------------------------------------------
      IF (ierr .eq. 0) CALL write_out_v3p
      CALL second0(t2)
!      PRINT 120, t2-t1, n_diagn_c

 100  FORMAT(a,1pe12.3)
 120  FORMAT(
     1 /,' TIME IN V3POST CODE:',1pe12.2,' SEC FOR ',i8.4,' SENSORS')
!      CALL flush(6)
      END SUBROUTINE v3post_sub
