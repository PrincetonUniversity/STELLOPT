      SUBROUTINE get_input_v3p (input_id, ideq_cmdline,
     1   plasfcn, coilfcn)
      USE stel_kinds
      USE stel_constants
      USE v3post_rfun
!DEC$ IF DEFINED (MPI_OPT)
      USE read_response
!DEC$ ELSE
      USE read_response_nompi
!DEC$ ENDIF
      USE safe_open_mod
      USE read_wout_mod, ONLY: read_wout_file, ns, presf
      IMPLICIT NONE
!----------------------------------------------------------------------
! D U M M Y Arguments Declarations
!----------------------------------------------------------------------
      CHARACTER (len=*) :: plasfcn, coilfcn
      CHARACTER (len=*)  ::  input_id, ideq_cmdline
!----------------------------------------------------------------------
! L O C A L Variable Declarations
!----------------------------------------------------------------------
      INTEGER ::  istat, ierr, nin=100
      CHARACTER (len=100)  :: listin
      REAL(rprec) :: sump
      NAMELIST/v3p_in/listin, idrfun, eqtype, ideqfile, myid, freeb,
     &                lsurf

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
!----------------------------------------------------------------------
      CALL safe_open(nin, istat, TRIM(input_id), 'old', 'formatted')
      IF (istat .ne. 0) STOP 'safe_open iostat != 0 in v3post'
!----------------------------------------------------------------------
!-- read inputs from namelist                                        --
!----------------------------------------------------------------------
      READ (nin, nml=v3p_in,iostat=istat)
      CLOSE(nin)
      ideqfile = TRIM(ideq_cmdline)
      IF (TRIM(myid) .eq. 'none') myid=' '
      plasfcn = TRIM(listin)
      coilfcn = "crfun_"//TRIM(idrfun)//".nc"
      IF (myid(1:1) .ne. ' ') coilfcn=
     .  "crfun"//TRIM(myid)//"_"//TRIM(idrfun)//".nc"
!
!----------------------------------------------------------------------
!-- read VMEC equilibrium data from wout                                  --
!----------------------------------------------------------------------
        eq_file = 'wout_' // TRIM(ideqfile)
        CALL  read_wout_file(eq_file, ierr)
        IF (ierr .ne. 0) THEN
           WRITE (6,*) 'Error reading wout file: ierr = ', ierr
           STOP
        ENDIF
!
!Added by EALazarus (2005): Accounts for negative pressure regions arising during optimization
!        ppenalty = 0
!        sump = SUM(ABS(presf))
!        IF(sump .ne. ppenalty)
!     &       ppenalty = -ns*(SUM(presf)-SUM(ABS(presf)))/sump

      END SUBROUTINE get_input_v3p
