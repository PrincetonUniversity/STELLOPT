      PROGRAM boozer_xform
      USE booz_params
      USE safe_open_mod
      USE mpi_params, ONLY: master, myid, numprocs, ierr_mpi,
     1                      MPI_COMM_BOOZER, MPI_CALC_MYRANGE
      USE mpi_inc
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: istat, jrad, iread, numargs
      INTEGER :: mystart, myend, jsurf, ncount
      REAL(rprec) :: t1, t2
      CHARACTER(LEN=50) :: arg1, arg2
      CHARACTER(LEN=120) :: extension
C-----------------------------------------------
!
!     driver: reads from command line file the wout file extension and surfaces
!     (half-radial) at which the boozer coordinates are required
!     writes the boozer coordinates to a file, boozmn.extension
!
!     call this as follows:
!
!          xbooz_xform input.boz [T or F]
!
!     WHERE input.boz CONTAINS the mboz, nboz, wout file extension and the jrad values (as a
!     blank-delimited list, not necessarily on a single line):
!
!     mboz   nboz
!     FILE_EXTENSION (does NOT have to include the .nc or .txt extension for netcdf,text file)
!     1  3   5   10  12
!
!     The OPTIONAL command line argument, (T) or (F), Allows the user to turn off screen
!     output IF set to F.
!
!     CALL xbooz_xform -h brings up a help screen
!
!     MPI parallelization: the flux surfaces selected by lsurf_boz are
!     distributed across ranks in MPI_COMM_BOOZER.  Each rank computes
!     its subset of surfaces and writes into the corresponding columns
!     of the packed output arrays (rmncb, zmnsb, ...).  The remaining
!     columns are left at zero (from allocate_boozer), so a single
!     MPI_REDUCE with MPI_SUM gathers the full result to the master
!     rank, which writes boozmn.
!

!
!     INITIALIZE MPI
!
      myid = master
      numprocs = 1
#if defined(MPI_OPT)
      CALL MPI_INIT(ierr_mpi)
      CALL MPI_COMM_DUP(MPI_COMM_WORLD, MPI_COMM_BOOZER, ierr_mpi)
      CALL MPI_COMM_RANK(MPI_COMM_BOOZER, myid, ierr_mpi)
      CALL MPI_COMM_SIZE(MPI_COMM_BOOZER, numprocs, ierr_mpi)
#endif

      lscreen = .true.                  !!Default, write to screen
      IF (myid .ne. master) lscreen = .false.

!
!     Read command line argument to get input file name
!
      CALL getcarg(1, arg1, numargs)
      IF (numargs .gt. 1) CALL getcarg(2, arg2, numargs)

      IF (numargs .lt. 1) THEN
         IF (myid .eq. master) THEN
            PRINT *,'Invalid command line in calling xbooz_xform'
            PRINT *,'Type xbooz_xform -h to get more information'
         END IF
         GOTO 2000
      ELSE IF (arg1 .eq. '-h' .or. arg1 .eq. '/h') THEN
         IF (myid .eq. master) THEN
            PRINT *,' ENTER INPUT FILE NAME ON COMMAND LINE'
            PRINT *,' For example: xbooz_xform in_booz.ext'
            PRINT *
            PRINT *,' WHERE in_booz.ext is the input file'
            PRINT *
            PRINT *,' Optional command line argument'
            PRINT *,' xbooz_xform <infile> (T or F)'
            PRINT *
            PRINT *,' where F suppresses output to the screen'
         END IF
         GOTO 2000
      ELSE IF (numargs .gt. 1) THEN
         IF (arg2(1:1).eq.'f' .or. arg2(1:1).eq.'F') lscreen = .false.
      ENDIF

!
!     Every rank opens the (read-only) input and wout files
!     independently; this avoids broadcasting the many arrays filled
!     by read_wout_booz.
!
      iread = unit_booz-1
      CALL safe_open (iread, istat, TRIM(arg1), 'old', 'formatted')
      IF (istat .ne. 0) THEN
         IF (myid .eq. master)
     1      PRINT *,'Error opening input file in booz_xform'
         GOTO 2000
      END IF

      READ (iread, *, iostat = istat) mboz, nboz
      READ (iread, *, iostat = istat) extension
      IF (istat .ne. 0) THEN
         IF (myid .eq. master)
     1      PRINT *,'Error reading input file in booz_xform'
         CLOSE(unit=iread)
         GOTO 2000
      END IF

!
!     READ IN PARAMETERS, DATA FROM WOUT FILE
!
      CALL read_wout_booz(extension, iread, istat)
      IF (istat .ne. 0) THEN
         IF (myid .eq. master)
     1      PRINT *,' ierr_vmec !=0 in booz_xform read_wout_booz'
         CLOSE(unit=iread)
         GOTO 1010
      END IF

      CLOSE (unit=iread)
      CALL second0(t1)

!
!     ONE-TIME BOOZER GRID / TRANSFORM SETUP (executed on every rank;
!     populates persistent arrays used by every call to boozer_coords)
!
      CALL boozer_setup

!
!     DISTRIBUTE PACKED SURFACE INDICES (1..jsize) ACROSS RANKS.
!     MPI_CALC_MYRANGE returns mystart>myend on ranks that do not
!     receive any work, which makes the DO loop a no-op there.
!
      mystart = 1
      myend   = jsize
#if defined(MPI_OPT)
      CALL MPI_CALC_MYRANGE(MPI_COMM_BOOZER, 1, jsize, mystart, myend)
#endif

!
!     COMPUTE BOOZER TRANSFORM ON THIS RANK'S SUBSET OF SURFACES
!
      DO jsurf = mystart, myend
         jrad = jlist(jsurf)
         IF (jrad .ge. 1 .and. jrad .le. ns) THEN
            IF (lsurf_boz(jrad)) CALL boozer_coords(jrad, jsurf)
         END IF
      END DO

#if defined(MPI_OPT)
!
!     GATHER PACKED OUTPUT ARRAYS TO MASTER VIA SUM-REDUCE.  Columns
!     not filled by a given rank remain zero (see allocate_boozer),
!     so summing over ranks reproduces the full serial result.
!
      ncount = mnboz*jsize
      IF (myid .eq. master) THEN
         CALL MPI_REDUCE(MPI_IN_PLACE, bmncb, ncount,
     1        MPI_DOUBLE_PRECISION, MPI_SUM, master,
     2        MPI_COMM_BOOZER, ierr_mpi)
         CALL MPI_REDUCE(MPI_IN_PLACE, rmncb, ncount,
     1        MPI_DOUBLE_PRECISION, MPI_SUM, master,
     2        MPI_COMM_BOOZER, ierr_mpi)
         CALL MPI_REDUCE(MPI_IN_PLACE, zmnsb, ncount,
     1        MPI_DOUBLE_PRECISION, MPI_SUM, master,
     2        MPI_COMM_BOOZER, ierr_mpi)
         CALL MPI_REDUCE(MPI_IN_PLACE, pmnsb, ncount,
     1        MPI_DOUBLE_PRECISION, MPI_SUM, master,
     2        MPI_COMM_BOOZER, ierr_mpi)
         CALL MPI_REDUCE(MPI_IN_PLACE, gmncb, ncount,
     1        MPI_DOUBLE_PRECISION, MPI_SUM, master,
     2        MPI_COMM_BOOZER, ierr_mpi)
      ELSE
         CALL MPI_REDUCE(bmncb, bmncb, ncount,
     1        MPI_DOUBLE_PRECISION, MPI_SUM, master,
     2        MPI_COMM_BOOZER, ierr_mpi)
         CALL MPI_REDUCE(rmncb, rmncb, ncount,
     1        MPI_DOUBLE_PRECISION, MPI_SUM, master,
     2        MPI_COMM_BOOZER, ierr_mpi)
         CALL MPI_REDUCE(zmnsb, zmnsb, ncount,
     1        MPI_DOUBLE_PRECISION, MPI_SUM, master,
     2        MPI_COMM_BOOZER, ierr_mpi)
         CALL MPI_REDUCE(pmnsb, pmnsb, ncount,
     1        MPI_DOUBLE_PRECISION, MPI_SUM, master,
     2        MPI_COMM_BOOZER, ierr_mpi)
         CALL MPI_REDUCE(gmncb, gmncb, ncount,
     1        MPI_DOUBLE_PRECISION, MPI_SUM, master,
     2        MPI_COMM_BOOZER, ierr_mpi)
      END IF

      IF (lasym_b) THEN
         IF (myid .eq. master) THEN
            CALL MPI_REDUCE(MPI_IN_PLACE, bmnsb, ncount,
     1           MPI_DOUBLE_PRECISION, MPI_SUM, master,
     2           MPI_COMM_BOOZER, ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE, rmnsb, ncount,
     1           MPI_DOUBLE_PRECISION, MPI_SUM, master,
     2           MPI_COMM_BOOZER, ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE, zmncb, ncount,
     1           MPI_DOUBLE_PRECISION, MPI_SUM, master,
     2           MPI_COMM_BOOZER, ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE, pmncb, ncount,
     1           MPI_DOUBLE_PRECISION, MPI_SUM, master,
     2           MPI_COMM_BOOZER, ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE, gmnsb, ncount,
     1           MPI_DOUBLE_PRECISION, MPI_SUM, master,
     2           MPI_COMM_BOOZER, ierr_mpi)
         ELSE
            CALL MPI_REDUCE(bmnsb, bmnsb, ncount,
     1           MPI_DOUBLE_PRECISION, MPI_SUM, master,
     2           MPI_COMM_BOOZER, ierr_mpi)
            CALL MPI_REDUCE(rmnsb, rmnsb, ncount,
     1           MPI_DOUBLE_PRECISION, MPI_SUM, master,
     2           MPI_COMM_BOOZER, ierr_mpi)
            CALL MPI_REDUCE(zmncb, zmncb, ncount,
     1           MPI_DOUBLE_PRECISION, MPI_SUM, master,
     2           MPI_COMM_BOOZER, ierr_mpi)
            CALL MPI_REDUCE(pmncb, pmncb, ncount,
     1           MPI_DOUBLE_PRECISION, MPI_SUM, master,
     2           MPI_COMM_BOOZER, ierr_mpi)
            CALL MPI_REDUCE(gmnsb, gmnsb, ncount,
     1           MPI_DOUBLE_PRECISION, MPI_SUM, master,
     2           MPI_COMM_BOOZER, ierr_mpi)
         END IF
      END IF
#endif

!
!     WRITE OUT CONVERTED RESULTS (master only)
!
      IF (myid .eq. master) CALL write_boozmn(extension)

 1010 CONTINUE
!
!     FREE MEMORY : USER MUST CALL BOOZER_COORDS AT END WITH LDEALLOC = TRUE
!     OTHERWISE MEMORY ALLOCATED WILL NOT BE AVAILABLE
!

      CALL free_mem_boozer

      CALL second0(t2)

      IF (myid .eq. master) PRINT 120, t2-t1
 120  FORMAT(/,' TIME IN BOOZER TRANSFORM CODE:',1pe12.2,' SEC')

 2000 CONTINUE
#if defined(MPI_OPT)
      CALL MPI_FINALIZE(ierr_mpi)
#endif

      END PROGRAM boozer_xform
