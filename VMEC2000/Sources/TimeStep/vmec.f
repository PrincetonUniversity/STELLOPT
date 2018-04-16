      PROGRAM vmec
      USE vmec_input
      USE vmec_seq
      USE safe_open_mod
!      USE precon2d, ONLY: ScratchFile
      USE vparams, ONLY: nlog, nlog0, nthreed
      USE vmec_params, ONLY: more_iter_flag,
     1                       bad_jacobian_flag,
     2    restart_flag, readin_flag, timestep_flag,
     3    output_flag, cleanup_flag,
     4    norm_term_flag, successful_term_flag ! J Geiger: for more iterations and full 3D1-output
      USE parallel_include_module
      USE parallel_vmec_module, ONLY: MyEnvVariables, 
     1    InitializeParallel, FinalizeParallel
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: nseq0 = 12
      CHARACTER(LEN=*), PARAMETER :: increase_niter = 
     1   "Try increasing NITER",
     2    bad_jacobian = "The jacobian was non-definite!",
     3    full_3d1output_request = "Full threed1-output request!"
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: numargs, ierr_vmec, index_end,
     1   iopen, isnml, iread, iseq, index_seq,
     2   index_dat, iunit, ncount, nsteps, i
      INTEGER :: ictrl(5)
      CHARACTER(LEN=120) :: input_file, seq_ext, reset_file_name, arg
      CHARACTER(LEN=120) :: log_file
      CHARACTER(LEN=120), DIMENSION(10) :: command_arg
      LOGICAL :: lscreen
      INTEGER :: RVC_COMM

#if defined(SKS)
      REAL(dp) :: ton, toff
      REAL(dp) :: totalton, totaltoff 
#endif
C-----------------------------------------------
!***
!                              D   I   S   C   L   A   I   M   E   R
!
!       You are using a beta version of the PROGRAM VMEC, which is currently
!       under development by S. P. Hirshman at the Fusion Energy Division,
!       Oak Ridge National Laboratory.  Please report any problems or comments
!       to him.  As a beta version, this program is subject to change
!       and improvement without notice.
!
!       1. CODE SYNOPSIS
!
!       THIS PROGRAM - VMEC (Variational Moments Equilibrium Code)  -
!       SOLVES THREE-DIMENSIONAL MHD EQUILIBRIUM EQUATIONS USING
!       FOURIER SPECTRAL (MOMENTS) METHODS. A CYLINDRICAL COORDINATE
!       REPRESENTATION IS USED (R-Z COORDINATES). THE POLOIDAL
!       ANGLE VARIABLE IS RENORMALIZED THROUGH THE STREAM FUNCTION
!       LAMBDA, WHICH IS SELF-CONSISTENTLY DETERMINED AND DIFFERENCED
!       VARIATIONALLY ON THE HALF-RADIAL MESH. THE POLOIDAL ANGLE IS
!       DETERMINED BY MINIMIZING <M> = m**2 S(m) , WHERE S(m) =
!       Rm**2 + Zm**2 . AN EVEN-ODD DECOMPOSITION IN THE POLOIDAL MODE
!       NO. OF R,Z, AND LAMDA IS USED TO IMPROVE RADIAL RESOLUTION.
!       A FREE-BOUNDARY OPTION IS AVAILABLE (FOR lfreeb=T), WITH A
!       USER-SUPPLIED DATA-FILE "MGRID" NEEDED TO COMPUTE THE PLASMA
!       VACUUM FIELD COMPONENTS BR, BPHI, BZ (see SUBROUTINE BECOIL)
!
!       THE MAGNETIC FIELD IS REPRESENTED INTERNALLY AS FOLLOWS:
!
!       B(s,u,v) = grad(phiT) X ( grad(u) + grad(lambda) ) +
!
!                  iota(s) * grad(v) X grad(phiT)
!
!       WHERE phiT is the toroidal flux (called phi in code) and
!       u,v are the poloidal, toroidal angles, respectively.
!
!       2. ADDITIONAL CODES REQUIRED
!       For the fixed boundary calculation, the user must provide the Fourier
!       coefficients for the plasma boundary (the last surface outside of which
!       the pressure gradient vanishes). For ALL but the simplest geometry, the
!       SCRUNCH code (available from R. Wieland), based on the DESCUR curve-fitting
!       code, can be used to produce the optimized VMEC Fourier representation for
!       an arbritrary closed boundary (it need not be a 'star-like' DOmain, nor
!       need it possess vertical, or 'stellarator', symmetry).
!
!       For the free boundary calculation, the MAKEGRID code (available upon
!       request) is needed to create a binary Green''s FUNCTION table for the
!       vacuum magnetic field(s) and, IF data analysis is to be done, flux and
!       field loops as well. The user provides a SUBROUTINE (BFIELD) which can be
!       called at an arbitrary spatial location and which should RETURN the three
!       cylindrical components of the vacuum field at that point. (Similary,
!       locations of diagnostic flux loops, Rogowski coils, etc. are required IF
!       equilibrium reconstruction is to be done.)
!
!       Plotting is handled by a stand-alone package, PROUT.NCARG (written by
!       R. M. Wieland). It uses NCAR-graphics calls and reads the primary VMEC output
!       file, WOUT.EXT, WHERE 'EXT' is the command-line extension of the INPUT file.
!
!
!       3. UNIX SCRIPT SETUP PARAMETERS
!       The VMEC source code (vmec.lsqh) is actually a UNIX script file which uses
!       the C-precompiler to produce both the machine-specific Fortran source and a
!       make-file specific to ANY one of the following platforms:
!
!       IBM-RISC6000, CRAY, ALPHA (DEC-STATION), HP-UX WORKSTATION,
!       WINDOWS-NT, DEC-VMS
!
!       Additional platforms are easy to add to the existing script as required.
!
!
!       4. FORTRAN PARAMETER STATEMENTS set by user
!       In the Fortran-90 version of VMEC these PARAMETER statements have
!       been replaced by dynamic memory allocation. So the user should set the
!       run-time parameters ns (through ns_array), mpol, ntor in the NAMELIST INDATA.
!
!
!       Added features since last edition (see vmec_params for revision history list)
!       1. Implemented preconditioning algorithm for R,Z
!       2. The physical (unpreconditioned) residuals are used
!          to determine the level of convergence
!       3. The original (MOMCON) scaling of lambda is used, i.e.,
!          Bsupu = phip*(iota - lamda[sub]v)/SQRT(g). This is needed to
!          maintain consistency with the time-stepper for arbitrary PHIP.
!
!       WRITTEN BY S. P. HIRSHMAN (8/28/85 - REVISED 3/1/86) BASED ON
!       1. S. P. Hirshman and J. C. Whitson, Phys. Fluids 26, 3553 (1983).
!       2. S. P. Hirshman and H. K. Meier, Phys. Fluids 28, 1387 (1985).
!       3. S. P. Hirshman and D. K. Lee, Comp. Phys. Comm. 39, 161 (1986).
!

!     Local variables
!
!     ictrl:   array(5) of control variables for running "runvmec" routine
!              see "runvmec" for a description
!

!
!     Read in command-line arguments to get input file or sequence file,
!     screen display information, and restart information
!
      INTERFACE
         SUBROUTINE runvmec(ictrl_array, input_file0, 
     1                      lscreen, RVC_COMM, reset_file_name)
         IMPLICIT NONE
         INTEGER, INTENT(inout), TARGET :: ictrl_array(5)
         LOGICAL, INTENT(in) :: lscreen
         CHARACTER(LEN=*), INTENT(in) :: input_file0
         INTEGER, INTENT(in), OPTIONAL :: RVC_COMM
         CHARACTER(LEN=*), OPTIONAL :: reset_file_name
         END SUBROUTINE runvmec
      END INTERFACE

#if defined(SKS)
      CALL MyEnvVariables
      CALL InitializeParallel
      CALL MPI_COMM_DUP(MPI_COMM_WORLD,RVC_COMM,MPI_ERR)
      CALL second0(totalton)
      ton = totalton
#endif
      CALL getcarg(1, command_arg(1), numargs)
      DO iseq = 2, numargs
         CALL getcarg(iseq, command_arg(iseq), numargs)
      END DO
#if defined(SKS)
      CALL second0(toff)
      get_args_time = get_args_time + (toff -ton)
#endif
      lscreen = .false.
      IF(grank.EQ.0) lscreen = .true.
      reset_file_name = " "

      IF (numargs .lt. 1) THEN
         STOP 'Invalid command line'
      ELSE IF (command_arg(1).eq.'-h' .or. command_arg(1).eq.'/h') THEN
         PRINT *,
     1   ' ENTER INPUT FILE NAME OR INPUT-FILE SUFFIX ON COMMAND LINE'
         PRINT *
         PRINT *,' For example: '
         PRINT *,'    xvmec input.tftr OR xvmec tftr ',
     1           'OR xvmec ../input.tftr'
         PRINT *
         PRINT *,' Sequence files, containing a list of input files',
     1           ' are also allowed. For example: '
         PRINT *,'    xvmec input.tftr_runs'
         PRINT *
         PRINT *,' Here, input.tftr_runs contains a &VSEQ namelist',
     1           ' entry'
         PRINT *
         PRINT *,' Additional (optional) command arguments are',
     1           ' allowed:'
         PRINT *
         PRINT *,'  xvmec <filename> [noscreen] [reset=reset_wout_file]'
         PRINT *
         PRINT *,' noscreen: supresses all output to screen ',
     1           ' (default, or "screen", displays output)'
         PRINT *,' name of reset wout file (defaults to none)'

         STOP
      ELSE
         DO iseq = 2, MIN(numargs,10)
            arg = command_arg(iseq)
            IF (TRIM(arg).eq.'noscreen' .or. TRIM(arg).eq.'NOSCREEN')
     1         lscreen = .false.
            index_end = INDEX(arg, "reset=")
            index_seq = MAX(INDEX(arg, "RESET="), index_end)
            IF (index_seq .gt. 0) reset_file_name = arg(index_seq+6:)
         END DO
      END IF

!
!     Determine type of file opened (sequential or input-data)
!     ARG1 (char var)
!          By DEFAULT, ARG1 obtained from the command
!          line is parsed as follows to determine the input data file(s):
!               a. Attempt to OPEN file ARG1 (full path + file name).
!                  Look for the VSEQ NAMELIST to obtain nseq, nseq_select, and
!                  extension array. If they exist and nseq>0, VMEC will run
!                  sequentially using input determined from the array EXTENSION[i]
!                  or input.EXTENSION[i]
!               b. If the command argument is not a sequence NAMELIST, THEN the data file
!                  ARG1 or input.ARG1 is READ directly, with NSEQ=1.
!
      arg = command_arg(1)
      index_dat = INDEX(arg,'.')
      index_end = LEN_TRIM(arg)
      IF (index_dat .gt. 0) THEN
         seq_ext  = arg(index_dat+1:index_end)
         input_file = TRIM(arg)
      ELSE
         seq_ext = TRIM(arg)
         input_file = 'input.'//TRIM(seq_ext)
      END IF

      nseq = 1
      nseq_select(1) = 1
      extension(1) = input_file
!
!     READ IN NAMELIST VSEQ TO GET ARRAY
!     OF INPUT FILE EXTENSIONS AND INDEXING ARRAY, NSEQ_SELECT
!
      nlog = nlog0
      iunit = nseq0
      DO iseq = 1, 2
         IF (iseq .EQ. 1) THEN
           arg = input_file
         ELSE
           arg = seq_ext
         END IF
#if defined(SKS)
         CALL second0(ton)
#endif
         CALL safe_open(iunit, iopen, TRIM(arg), 'old', 'formatted')
#if defined(SKS)
         CALL second0(toff)
         safe_open_time = safe_open_time + (toff - ton)
#endif
         IF (iopen .eq. 0) THEN
           DO ncount = 1, nseqmax
              nseq_select(ncount) = ncount
           END DO
#if defined(SKS)
           CALL second0(ton)
#endif
           CALL read_namelist (iunit, isnml, 'vseq')
#if defined(SKS)
           CALL second0(toff)
           read_namelist_time = read_namelist_time + (toff - ton)
#endif
!
!       OPEN FILE FOR STORING SEQUENTIAL RUN HISTORY
!
           IF (isnml .eq. 0) THEN
              IF (nseq .gt. nseqmax) STOP 'NSEQ>NSEQMAX'
              log_file = 'log.'//seq_ext
#if defined(SKS)
              CALL second0(ton)
#endif
              CALL safe_open(nlog, iread, log_file, 'replace',
     1           'formatted')
#if defined(SKS)
              CALL second0(toff)
              safe_open_time = safe_open_time + (toff - ton)
#endif
              IF (iread .NE. 0) THEN
                 PRINT *, log_file,
     1           ' LOG FILE IS INACCESSIBLE: IOSTAT= ',iread
                 STOP 3
              ELSE
                 EXIT        !!Break out of loop
              END IF
           ENDIF
        ENDIF
        CLOSE (iunit)
      END DO

!
!     CALL EQUILIBRIUM SOLVER
!
!     nseq_select:      If sequence file (VSEQ NAMELIST given with nseq >0)
!                       array giving indices into EXTENSION array prescribing
!                       the order in which the input files are run by VMEC
!     nseq:             number of sequential VMEC runs to make
!
!
!     CALL VMEC WITH POSSIBLE SEQUENCE EXTENSION (SEQ_EXT)
!     AND ARRAY OF INPUT FILE EXTENSIONS (EXTENSION)
!
      ictrl = 0

!     GOTO 200  !ENABLE THIS STATEMENT TO TEST REVERSE-COMMUNICATION STOP/RESTART CODING

      SEQ: DO iseq = 1, nseq
         index_seq = nseq_select(iseq)
         ictrl(1) = restart_flag+readin_flag+timestep_flag
     1            + output_flag+cleanup_flag                !Sets all flags
         ictrl(2) = 0
!         ictrl(3) = 100
!         ictrl(4) = 2
         ictrl(5) = iseq-1
         ncount = 0
         IF (iseq .GT. 1) THEN
            reset_file_name =
#ifdef NETCDF
     1       'wout_' // TRIM(extension(index_seq-1)) // ".nc"
#else
     1       'wout.' // TRIM(extension(index_seq-1))
#endif
         END IF
!    
!     SET UP A "REVERSE-COMMUNICATION" LOOP FOR RUNNING VMEC
!

 100     CONTINUE

         RVCCALLNUM=1
         CALL runvmec (ictrl, extension(index_seq), 
     1                 lscreen, RVC_COMM, reset_file_name)

         ierr_vmec = ictrl(2)

         SELECT CASE (ierr_vmec)
         CASE (more_iter_flag)                                !Need a few more iterations to converge
            IF (grank .EQ. 0) THEN
               IF(lscreen ) WRITE (6, '(1x,a)') increase_niter
               WRITE (nthreed, '(1x,a)') increase_niter
               WRITE (nthreed, '(1x,a)') "PARVMEC aborting..."
               CALL FLUSH(nthreed)
#if defined(SKS)
               CALL MPI_Abort(MPI_COMM_WORLD, MPI_ERR)
#endif
            END IF
! J Geiger: if lmoreiter and lfull3d1out are false
!           the o-lines (original) are the only
!           ones to be executed.
            IF (lmoreiter) THEN                                 ! J Geiger: --start--
              DO i=2,max_main_iterations                        ! Changes to run
                ictrl(1) = timestep_flag                        ! some more iterations if requested
                ictrl(3) = niter                                ! - this is the number of iterations
                RVCCALLNUM=2
                CALL runvmec (ictrl, extension(1), lscreen,
     1       RVC_COMM, reset_file_name)    ! - the second iteration run with ictrl(3) iterations
                IF (ictrl(2).EQ.more_iter_flag) THEN
                  IF (grank .EQ. 0) THEN
                    WRITE (nthreed, '(1x,a)') increase_niter
                    IF(lscreen) WRITE (6, '(1x,a)') 
     1                            increase_niter
                  ENDIF
                ENDIF
              ENDDO
              ictrl(1) = output_flag+cleanup_flag               ! - Output, cleanup
              IF(ictrl(2).ne.successful_term_flag)
     &                ictrl(2)=successful_term_flag             ! - force success flag to get full threed1-output!
              ictrl(3) = 0                                      ! - this is the number of iterations
              RVCCALLNUM=3
              CALL runvmec (ictrl, extension(1), lscreen,       ! - final call.
     1        RVC_COMM, reset_file_name)
            ELSE                                                ! else-branch contains original code. 
#if defined(SKS)		                                               
              CALL MPI_Barrier(RVC_COMM,MPI_ERR)
#endif
              IF(grank.EQ.0) 
     1           PRINT *,"vmec.f:There is a bug in this branch."
              STOP 

              ictrl(1) = output_flag+cleanup_flag               !Output, cleanup  ! o-lines
              ictrl(2) = 0                                      ! o-lines
              IF(lfull3d1out) THEN
                ictrl(2)=successful_term_flag
                IF (grank .EQ. 0) THEN
                  WRITE(6,'(1x,a)') full_3d1output_request
                  WRITE(nthreed,'(1x,a)') full_3d1output_request
                END IF
              ENDIF
              !STOP 111
              RVCTRIGGER=.TRUE.
              RVCCALLNUM=4
              CALL runvmec (ictrl, extension(1), lscreen,        ! o-lines
     1       RVC_COMM, reset_file_name)
              RVCTRIGGER=.FALSE.
            ENDIF                                                ! J Geiger: -- end --

         CASE (bad_jacobian_flag)                             !Bad jacobian even after axis reset and ns->3
           IF(grank.EQ.0) THEN
             IF (lscreen) WRITE (6, '(/,1x,a)') bad_jacobian
             WRITE (nthreed, '(/,1x,a)') bad_jacobian
           END IF
         CASE DEFAULT
         END SELECT
      END DO SEQ

      GOTO 300

!REVERSE-COMMUNICATION TEST LOOP: SHOULD BE OFF FOR NORMAL VMEC2000 RUN
!ONLY SHOWS HOW TO IMPLEMENT EXTERNAL STOP/START OF VMEC
 200  CONTINUE

      nsteps = 50
      ictrl(1) = restart_flag+readin_flag          !Initialization only
      ictrl(3) = nsteps
      ictrl(4) = 1                                 !Go to fine grid directly (assumes it is grid index=1)
      RVCCALLNUM=5
      CALL runvmec (ictrl, extension(1), lscreen,
     1       RVC_COMM, reset_file_name)

      ictrl(1) = timestep_flag  + output_flag      !Set timestep flag (output_flag, too, if wout needed)
      ictrl(1) = timestep_flag                     !Set timestep flag (output_flag, too, if wout needed)

      DO iopen = 1, 3                              !Scan through grids (3, for this example)
!      DO iopen = 1,1                              !Scan through grids
         ictrl(4) = iopen          
         DO ncount = 1, MAX(1,niter/nsteps)
           RVCCALLNUM=6
            CALL runvmec (ictrl, extension(1), lscreen,
     1       RVC_COMM, reset_file_name)
            PRINT *,' BREAK HERE'
            ierr_vmec = ictrl(2)
            IF (ierr_vmec .ne. more_iter_flag) EXIT
         END DO
      END DO

      ictrl(1) = output_flag+cleanup_flag      !Output, cleanup
      RVCCALLNUM=7
      CALL runvmec (ictrl, extension(1), lscreen,
     1       RVC_COMM, reset_file_name)

 300  CONTINUE

      CLOSE (nlog)

!      IF (ScratchFile .ne. "") THEN
!         OPEN(unit=20, file=ScratchFile, iostat=iopen)
!         IF (iopen .eq. 0) CLOSE (20, status='delete')
!      END IF

#if defined(SKS)
      CALL second0(totaltoff)
      total_time = total_time + (totaltoff-totalton)
      toff = totaltoff
      !IF (.NOT.LV3FITCALL.AND.lactive) CALL PrintTimes
      IF (.NOT.LV3FITCALL.AND.lactive) CALL WriteTimes('timings.txt')
      CALL FinalizeParallel
#endif

      END PROGRAM vmec
