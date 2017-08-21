      MODULE lmpar_mod
      USE stel_kinds
      INTEGER :: nscan, ldfjac
      INTEGER, DIMENSION(:), POINTER :: ipvt
      REAL(rprec) :: pnorm, fnorm1, delta, par, spread_ratio
      REAL(rprec), DIMENSION(:), POINTER :: wa2p, wa3p, wa4p
      REAL(rprec), DIMENSION(:), POINTER :: diag, qtf
      REAL(rprec), DIMENSION(:,:), POINTER :: fjac
      LOGICAL :: lfirst_lm

      CONTAINS

      SUBROUTINE levmarq_param_mp(x, wa1, wa2, wa3, wa4,
     1     nfev, m, n, iflag, fcn, lev_step_range,fnorm_min,
     2     xvmin,xvmax)
!      SUBROUTINE levmarq_param_mp(x, wa1, wa2, wa3, wa4,
!     1     nfev, m, n, iflag, fcn)
      USE fdjac_mod, ONLY: flag_cleanup, flag_cleanup_lev, jac_index,
     1                     n_red
      USE mpi_params
      USE safe_open_mod
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      INCLUDE 'mpif.h'                                       !mpi stuff
!DEC$ ENDIF
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: n, m
      INTEGER :: iflag, nfev, lev_step_range                            !PPPL
!      INTEGER :: iflag, nfev
      REAL(rprec), INTENT(in) :: fnorm_min
      REAL(rprec), INTENT(in) :: x(n)
      REAL(rprec) :: wa1(n), wa2(n), wa3(n), wa4(m)
      REAL(rprec), INTENT(in), OPTIONAL :: xvmin(n)
      REAL(rprec), INTENT(in), OPTIONAL :: xvmax(n)
      EXTERNAL fcn
!DEC$ IF DEFINED (MPI_OPT)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1
C      REAL(rprec), DIMENSION(11), PARAMETER :: factors =
C     1  (/ 1.0_dp, 0.5_dp, 0.25_dp, 0.128_dp, 0.75_dp,
C     2      1.25_dp, 1.5_dp, 0.9_dp, 1.1_dp, 1.75_dp, 2.1_dp /)
!      real(rprec), dimension(11), parameter :: factors =
!     1  (/ 1.0_dp, 0.5_dp, 0.25_dp, 0.128_dp, 2.1_dp, 0.75_dp,
!     2      1.25_dp, 1.5_dp, 0.9_dp, 1.1_dp, 1.75_dp /)                 !PPPL
      real(rprec), dimension(11), parameter :: factors =
     1  (/ 0.128_dp, 0.250_dp, 0.500_dp, 0.750_dp, 0.900_dp,
     2     1.000_dp, 0.110_dp, 1.250_dp, 1.500_dp, 1.750_dp,
     3     2.000_dp /)                                                 !SAL
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: iproc, iproc_min, nfact, num_lev, istat, ierr, j,
     1           iflag_min(1), i, iunit, nfev_output
      INTEGER, ALLOCATABLE, DIMENSION(:) :: iflag_array
      REAL(rprec) :: scale_factor
!      REAL(rprec), DIMENSION(:), ALLOCATABLE :: fnorm_array
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: fnorm_array,
     1                                          delta_array,
     2                                          par_array,
     3                                          diag_red !PPPL
      CHARACTER(LEN=1) :: ext, low_mark
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      REAL(rprec), EXTERNAL :: enorm
C-----------------------------------------------
!     Perform Numprocs different function calls (in parallel) to find the minimum norm (chi-sq).
!     MPI calls are used to determine which processor has the minimum norm and then the
!     information is sent to all processors using MPI Broadcasts (modifications made by D. A. Spong 8/27/2000).
!
!     Uses processors 0,...,numprocs-1 (which belong to the user-defined MPI_COMM_WORKERS communicator)

      nfact = SIZE(factors)
      num_lev = numprocs
      iproc = myid + 1
      ALLOCATE (iflag_array(numprocs), fnorm_array(numprocs),
     1          delta_array(numprocs), par_array(numprocs),
     2          diag_red(n_red),stat=istat)   !PPPL
      IF (istat .ne. 0) THEN
         IF (myid == master) THEN
            WRITE(6,*) 'ALLOCATION ERROR in levmarq_param_mp'
            WRITE(6,*) 'STAT = ',istat
            CALL FLUSH(6)
         END IF
         CALL MPI_BARRIER(MPI_COMM_STEL,ierr)
         CALL mpi_stel_abort(0)
      END IF
            
      
!
!       Do an exponential spread the first call (lfirst_lm=.true.) to see where we are
!     
      IF (lfirst_lm .and. num_lev > 2) THEN
         scale_factor = exp((iproc-1)*log(spread_ratio)/num_lev)          !PPPL
         lfirst_lm = .false.
      ELSE IF (num_lev > 2*nfact) THEN
         scale_factor = (iproc*MAXVAL(factors))/num_lev
      ELSE IF (iproc .lt. nfact) THEN
         scale_factor = factors(iproc)
      ELSE
         scale_factor = ((iproc-nfact)*MINVAL(factors))/
     1                   (num_lev-nfact)
      END IF

      delta = scale_factor*delta
      
      DO i = 1, n_red
         j = jac_index(i)
         diag_red(i) = diag(j)
      END DO

      CALL lmpar (n_red, fjac, ldfjac, ipvt, diag_red, qtf,
     1            delta, par, wa1, wa2, wa3, wa4)
     
!
!     store the direction p and x + p. calculate the norm of p.
!

      IF (par .eq. zero) wa1 = wa1*scale_factor
      wa1 = -wa1
      wa2 = x
      wa3 = 0
      DO i = 1, n_red
         j = jac_index(i)
         wa2(j) = x(j) + wa1(i)
         wa3(j) = diag(j)*wa1(i)
      END DO
      ! Begin MJL additions
      print *,"i=",i
      print *,"j=",j
      print *,"wa2=",wa2
      print *,"wa3=",wa3
      print *,"About to try present(xvmin)"
      IF (PRESENT(xvmin)) THEN
         print *,"AAA Here comes size(xvmin):"
         print *,size(xvmin)
         print *,"Here comes shape(xvmin):"
         print *,shape(xvmin)
         print *,"Here comes xvmin:"
         print *,xvmin
         WHERE(wa2 < xvmin) wa3 = 0
         print *,"BBB"
         WHERE(wa2 < xvmin) wa2 = xvmin
         print *,"CCC"
      END IF
      print *,"Done with present(xvmin)"
      IF (PRESENT(xvmax)) THEN
         WHERE(wa2 > xvmax) wa3 = 0
         WHERE(wa2 > xvmax) wa2 = xvmax
      END IF
      print *,"Done with present(xvmax)"
      pnorm = enorm(n,wa3)
      
!
!     evaluate the function at x + p and calculate its norm.
!     Only do for 0 < myid < n processors (the MPI_COMM_WORKERS group),
      iflag = iproc
      CALL MPI_BARRIER(MPI_COMM_STEL,ierr)
      IF (ierr .ne. 0) CALL mpi_stel_abort(ierr)   
      CALL fcn (m, n, wa2, wa4, iflag, nfev)
      CALL MPI_BARRIER(MPI_COMM_STEL,ierr)
      IF (ierr .ne. 0) CALL mpi_stel_abort(ierr)      
      fnorm1 = enorm(m,wa4)
     

!
!     Create the xvec.dat file
!
      j=0; iunit = 27; istat = 0
      DO j = 0, numprocs-1
         IF (myid == j) THEN
            CALL safe_open(iunit,istat,'xvec.dat','unknown',
     1                     'formatted', ACCESS_IN='APPEND')
            WRITE(iunit,'(2(2X,I5.5))') n,nfev+j+1 ! +1 necessary to because j is PID not iteration number
            WRITE(iunit,'(10ES22.12E3)') wa2(1:n)
            WRITE(iunit,'(ES22.12E3)') fnorm1
            CLOSE(iunit)
         END IF
         CALL MPI_BARRIER(MPI_COMM_STEL,ierr)
         IF (ierr .ne. 0) CALL mpi_stel_abort(ierr)
      END DO

!
!     Gather iflag information to all processors and check for iflag < 0
!
      CALL MPI_ALLGATHER(iflag, 1, MPI_INTEGER, iflag_array, 1,
     1     MPI_INTEGER, MPI_COMM_STEL, ierr)
      IF (ierr .ne. 0) CALL mpi_stel_abort(ierr)

      iflag = minval(iflag_array)
      if (iflag .lt. 0) return

!
!     Find processor with minimum fnorm1 value
!
      CALL MPI_ALLGATHER(fnorm1, 1, MPI_REAL8, fnorm_array, 1,
     1     MPI_REAL8, MPI_COMM_STEL, ierr)
      IF (ierr .ne. 0) CALL mpi_stel_abort(ierr)
      iflag_array(1:1) = MINLOC(fnorm_array)
      iproc_min = iflag_array(1) - 1
      CALL MPI_ALLGATHER(delta, 1, MPI_REAL8, delta_array, 1,
     1     MPI_REAL8, MPI_COMM_STEL, ierr)                             !PPPL
      IF (ierr .ne. 0) CALL mpi_stel_abort(ierr)             !PPPL
      CALL MPI_ALLGATHER(par, 1, MPI_REAL8, par_array, 1,
     1     MPI_REAL8, MPI_COMM_STEL, ierr)                             !PPPL
      IF (ierr .ne. 0) CALL mpi_stel_abort(ierr)             !PPPL

      !PPPL WAY
      IF( myid .eq. master) THEN
         DO j=1, num_lev
            ext = ' '
            IF (j == 0) ext = '*'
            low_mark = ' '
            IF (j == iproc_min+1) low_mark = '*'

            WRITE(6, '(2x,i6,8x,i3,4x,2(3x,es12.4,a),(3x,es12.4))')
     1         j+nfev, j, fnorm_array(j)**2, low_mark,
     2         par_array(j), ext, delta_array(j)
         END DO
         WRITE(6, '(a)') '  '

         CALL flush(6)
      END IF
      
      fnorm1 = fnorm_array(iproc_min+1)                                 !PPPL
      !delta  = delta_array(iproc_min+1)                                 !PPPL
      !par    = par_array(iproc_min+1)                                   !PPPL
                                  !PPPL
      
      
!
!     Calc lev_step_range (PPPL)
!
      lev_step_range = 0                                                
      iflag_array(1:1) = maxloc(delta_array)                             
      if( iproc_min == iflag_array(1) - 1 ) then
         lev_step_range = 1
      else if( iproc_min == iflag_array(numprocs) - 1 ) then
         lev_step_range = -1
      endif
      
      ! Only update delta and par if actual minimum found
      IF (fnorm1**2 < 1.0E12) 
     1      delta  = delta_array(iproc_min+1)                                 !PPPL
      IF (fnorm1**2 < 1.0E12)
     1      par    = par_array(iproc_min+1) 
      IF (fnorm1**2 >= 1.0E12)  lev_step_range = 0     

 100  CONTINUE
 
      DEALLOCATE (iflag_array, fnorm_array, delta_array, par_array,
     1            diag_red)     !PPPL

!
!     Broadcast all relevant scalars and arrays from the
!     processor with minimum fnorm1 to the other processors,
!     overwriting their data. Note: diag, ipvt are same already on
!     ALL processors. wa3 is overwritten...
!
      CALL MPI_BCAST(pnorm,1,MPI_REAL8,iproc_min,
     1     MPI_COMM_STEL,ierr)
      IF (ierr .ne. 0) GOTO 3000
      CALL MPI_BCAST(x,n,MPI_REAL8,iproc_min,
     1     MPI_COMM_STEL,ierr)
      IF (ierr .ne. 0) GOTO 3000
      CALL MPI_BCAST(wa1,n,MPI_REAL8,iproc_min,
     1     MPI_COMM_STEL,ierr)
      IF (ierr .ne. 0) GOTO 3000
      CALL MPI_BCAST(wa2,n,MPI_REAL8,iproc_min,
     1     MPI_COMM_STEL,ierr)
      IF (ierr .ne. 0) GOTO 3000
      CALL MPI_BCAST(wa4,m,MPI_REAL8,iproc_min,
     1     MPI_COMM_STEL,ierr)
      IF (ierr .ne. 0) GOTO 3000

!
!     CLEANUP AFTER LEVENBERG-MARQUARDT LOOP AS NEEDED (WA4 IS NOT CHANGED)
!
   
      nfev_output = nfev
      IF (myid == iproc_min .and. fnorm_min > fnorm1) THEN
         iflag = flag_cleanup_lev        ! SAL 6/26/12
         nfev_output  = nfev+iproc_min+1
      ELSE
         iflag = -3  ! This is a dummy so everyone else does nothing.
      END IF
      CALL fcn (m, n, wa2, wa4, iflag, nfev_output)             !Contains Bcast Barrier
      CALL MPI_BARRIER(MPI_COMM_STEL,ierr)
      IF (ierr .ne. 0) CALL mpi_stel_abort(ierr)

      nfev = nfev + num_lev


      RETURN

 3000 CONTINUE

      WRITE (6, *) 'MPI_BCAST error in LEVMARQ_PARAM_MP, ierr = ', ierr
      CALL mpi_stel_abort(ierr)
!DEC$ ENDIF
      END SUBROUTINE levmarq_param_mp

      SUBROUTINE levmarq_param(x, wa1, wa2, wa3, wa4,
     1     time, nfev, m1, n1, iflag, fcn)
      USE fdjac_mod, mp => m, np => n
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: nfev, n1, m1, iflag
      REAL(rprec) :: time
      REAL(rprec), TARGET :: x(n1), wa1(n1), wa2(n1), wa3(n1), wa4(m1)
      EXTERNAL fcn
!DEC$ IF .NOT.DEFINED (MPI_OPT)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: j, istat, iread, ic1, ic2, irate, count_max,
     1     jmin
      REAL(rprec), DIMENSION(num_lm_params) ::
     1      fnorm_min, pnorm_min, delta_min, par_min
      CHARACTER(LEN=1) :: ext, low_mark
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL multiprocess, lmpar_parallel
C-----------------------------------------------
!
!     Initialize module variables for use in lmpar_parallel
!
      xp => x
      wap => wa1
      wa2p => wa2
      wa3p => wa3
      wa4p => wa4
      np = n1
      mp = m1
      ncnt = nfev

      CALL SYSTEM_CLOCK(ic1, irate)

      CALL multiprocess(num_lm_params, max_processors,
     1    lmpar_parallel, fcn)

      CALL SYSTEM_CLOCK(ic2, irate, count_max)
      IF (ic2 .lt. ic1) ic2 = ic2 + count_max

      nfev = nfev + num_lm_params

!
!     Read in optimal wa1, wa2, wa4, par, delta, pnorm, fnorm1 value from file
!
      DO j = 1, num_lm_params

        READ (j+1000, iostat=iread) istat, iflag, pnorm_min(j),
     1        fnorm_min(j), par_min(j), delta_min(j)
        IF (iread .ne. 0) THEN
           WRITE (6, *) 'Error reading from file fort.', j+1000,
     1       ' in levmarq_param', ' IOSTAT = ', iread
           iflag = -15
        ELSE IF (j .ne. istat) THEN
           WRITE (6, *)
     1        'Incorrect value READ in for INDEX j in levmarq_param'
           iflag = -15
        END IF

        IF (iflag .ne. 0) RETURN

        IF (j .eq. 1) fnorm1 = fnorm_min(j)
        IF (fnorm_min(j) .le. fnorm1) THEN
           jmin = j
           fnorm1 = fnorm_min(jmin)
           pnorm  = pnorm_min(jmin)
           par    = par_min(jmin)
           delta  = delta_min(jmin)
!DEC$ IF DEFINED (CRAY)
           DO k = 1, n1
              READ (j+1000) wa1(k), wa2(k)
              DO istat = 1, n1
                 READ (j+1000) fjac(k, istat)
              END DO
           END DO
           DO k = 1, m1
              READ (j+1000) wa4(k)
           END DO
!DEC$ ELSE
           READ (j+1000) wa1, wa2, wa4, fjac(1:n1, 1:n1)
!DEC$ ENDIF

        END IF

        CLOSE (j+1000, status='delete')                        !!Needed to run correctly in multi-tasking...

      END DO

      DO j = 1, num_lm_params
         ext = ' '
         low_mark = ' '
         IF (j .eq. 1) ext = '*'
         IF (j .eq. jmin) low_mark = '*'
         WRITE (6, '(2x,i6,4x,2(3x,1es12.4,a),3x,1es12.4)') j+ncnt,
     1         fnorm_min(j)**2, low_mark, par_min(j), ext, delta_min(j)
      END DO

!
!     Do any special cleanup now for IFLAG = flag_cleanup. WA4 LEFT UNCHANGED
!
      iflag = flag_cleanup
      CALL fcn(m1, n1, x, wa4, iflag, ncnt)

      time = time + REAL(ic2 - ic1)/REAL(irate)                !!Time in multi-process CALL
!DEC$ ENDIF
      END SUBROUTINE levmarq_param

      END MODULE lmpar_mod
