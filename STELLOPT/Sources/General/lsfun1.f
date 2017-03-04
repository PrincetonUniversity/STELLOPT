      SUBROUTINE lsfun1(nopt, nvar, xc_opt, fvec, iflag, niter)
      USE optim
      USE safe_open_mod
      USE mpi_params
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      INCLUDE 'mpif.h'                                                          ! MPI
!DEC$ ENDIF
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER nopt, nvar, iflag, niter
      REAL(rprec), DIMENSION(nvar), INTENT(in) :: xc_opt
      REAL(rprec), DIMENSION(nopt), INTENT(out) :: fvec
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: iflag_clean = -100, iflag_singletask = -1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: istat, ncnt, ikey
      CHARACTER(len=200) :: temp_input
      LOGICAL :: lscreen, lvmec_succeed
C-----------------------------------------------
!     ON ENTRY:
!        NITER:           NUMBER OF LEVENBERG ITERATIONS (=0 INITIALLY)
!        FVEC:            CONTAINS NOTHING THAT SHOULD BE USED
!        XC_OPT:          CONTAINS PERTURBED X VALUES
!        IFLAG = -1,      THIS IS A SINGLE-TASKED CALL, PRINT TO SCREEN
!        IFLAG 1 to NOPT  MULTI-TASKED CALL, DO NOT PRINT TO SCREEN
!        IFLAG = -100,    SPECIAL CLEAN UP CODE CALLED
!
!     ON EXIT:
!        FVEC:         CONTAINS F(X)
!        XC_OPT:       UNCHANGED
!        IFLAG = 0     FINISHED WITHOUT A PROBLEM DETECTED
!        IFLAG < 0     ERROR SENT TO CALLER
!
!     SUGGESTIONS FOR IMPROVING CONVERGENCE
!     1. USE A REASONABLY SMALL VALUE FOR EPSFCN (3.E-4)
!     2. TRY 'LOOSENING' THE TOLERANCE IN FTOL ARRAY, i.e., RAISE
!        FTOL TO 1-5*E-9.
!     3. TRY 'TIGHTENING' THE TOLERANCE, i.e., FTOL <= 5.E-10
!     4. TRY PERTURBING R(m=1,n=0) and/or Z(m=1,n=0) A LITTLE
!
      IF (iflag <= iflag_clean) THEN
         CALL clean_up (nvar,.false., iflag)    !PPPL
!         CALL clean_up (nvar, iflag)    ! ORNL
         RETURN
      ENDIF
      CALL FLUSH(6)
C      PRINT *,'MYID: ',myid,' IN FCN'
      fvec(:nopt) = 0

!************************************************************************
!     COMPUTE UNIQUE EXTENSION FOR PARALLELIZED OPERATION
!************************************************************************
      istat = iflag
      IF (iflag .eq. iflag_singletask) istat = 0
      WRITE (temp_input,'(i5)') istat
      opt_ext = TRIM(seq_ext)//'_opt'//TRIM(ADJUSTL(temp_input))

      input_file = 'input.' // TRIM(opt_ext)

      IF (iflag .ge. 0) THEN
         ncnt = niter + iflag
         lscreen = .false.
      ELSE
         ncnt = niter
         lscreen = .true.
      END IF

!************************************************************************
!     RUN VMEC TO COMPUTE TARGET CRITERIA.
!************************************************************************
C      PRINT *,'MYID: ',myid,' IN FCN LOAD_PARAMS'
      CALL load_params (xc_opt, nvar, iflag, lscreen)
      lvmec_succeed = (ierr_vmec.eq.0 .and. iflag.eq.0)
C      PRINT *,'MYID: ',myid,' IN FCN LOAD_PARAMS (DONE)'

!DEC$ IF DEFINED (MPI_OPT)
!************************************************************************
!     Set up DYNAMIC worker communication group for processors that had a 
!     successful VMEC run (SUBGROUP of MPI_COMM_WORKERS)
!************************************************************************
C      CALL MPI_COMM_RANK(MPI_COMM_WORKERS, worker_id, ierr_mpi)
C      IF (ierr_mpi .ne. 0) STOP 'MPI error in lsfun1!'
C      IF (.not.lvmec_succeed) THEN
C         ikey = MPI_UNDEFINED
C      ELSE
C         ikey = WORKER_SPLIT_KEY+1
C      END IF
C      CALL MPI_COMM_SPLIT(MPI_COMM_WORKERS, ikey, worker_id, 
C     1                    MPI_COMM_WORKERS_OK, ierr_mpi)
C      IF (ierr_mpi .ne. 0) STOP 'MPI_COMM_SPLIT error in lsfun1'
!DEC$ ENDIF

!************************************************************************
!     LOAD USER-DEFINED TARGET FUNCTIONALS
!     EVEN IF VMEC FAILED (IERR_VMEC != 0)
!************************************************************************
      IF (lvmec_succeed) THEN
         iunit_opt_local = iunit_opt + 1
         CALL safe_open(iunit_opt_local, istat, 'output.' //
     1                  TRIM(opt_ext), 'replace', 'formatted')
C         PRINT *,'MYID: ',myid,' AFTER SAFEOPEN istat=',istat
C         IF (istat .ne. 0) THEN
C            iflag = -4
C            RETURN
C         END IF

!************************************************************************
!     CHECK FAILURE OF VMEC CONVERGENCE
!************************************************************************
C      PRINT *,'MYID: ',myid,' IN FCN LOAD_TARGETS'
         CALL load_target (xc_opt, fvec, opt_ext, nopt, ncnt,
     1                     iflag, lscreen)
         CLOSE (iunit_opt_local, iostat=istat)
C      PRINT *,'MYID: ',myid,' IN FCN LOAD_TARGETS (DONE)'
!DEC$ IF DEFINED (MPI_OPT)
C         CALL MPI_COMM_FREE(MPI_COMM_WORKERS_OK, ierr_mpi)
!DEC$ ENDIF
      ENDIF

!************************************************************************
!     IFLAG MAY HAVE FAILED IN LOAD_TARGET CALL, MUST CHECK IFLAG AGAIN
!
!     IF VMEC DID NOT CONVERGE, THEN DISCOURAGE THIS SEARCH DIRECTION 
!     BY IMPOSING A LARGE VALUE OF CHI-SQ.
!     IF LSCREEN (FIRST CALL), THEN STOP: VMEC COULD NOT GET STARTED
!************************************************************************
      IF (ierr_vmec.ne.0 .or. iflag.ne.0) THEN
         IF( nopt_alg == 0) THEN
            fvec(:nopt) = 10*SQRT(chisq_min/nopt)
         ELSE
            fvec(:nopt) = 10*SQRT(bigno/nopt)
         ENDIF
         IF (ierr_vmec.ne.0 .and. .not.lscreen) iflag = 0
      END IF
      RETURN

      END SUBROUTINE lsfun1
