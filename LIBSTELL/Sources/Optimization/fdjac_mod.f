      MODULE fdjac_mod
      USE stel_kinds

      INTEGER, PARAMETER :: flag_singletask = -1, flag_cleanup = -100
      INTEGER, PARAMETER :: flag_cleanup_jac = -100
      INTEGER, PARAMETER :: flag_cleanup_lev = -101
      INTEGER, PARAMETER :: flag_cleanup_lbfgsb = -102

      INTEGER :: m, n, ncnt, max_processors, num_lm_params
      INTEGER :: ix_min, jac_count, n_red                                !PPPL
      REAL(rprec), DIMENSION(:), POINTER :: xp, wap
      LOGICAL, DIMENSION(:), ALLOCATABLE :: flip
      integer, dimension(:), allocatable :: jac_order                    !PPPL
      INTEGER, DIMENSION(:), ALLOCATABLE :: jac_err, jac_index
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: h_order
      REAL(rprec) :: eps

      CONTAINS
 
      SUBROUTINE fdjac2_mp(fcn, m, n, x, fvec, fjac, ldfjac,
     1           iflag, ncnt, epsfcn, fnorm_min, x_min, wa)
      USE stel_kinds
      USE mpi_params
      USE mpi_inc
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: m, n, ldfjac, ncnt
      INTEGER :: iflag
      REAL(rprec), INTENT(in) :: epsfcn
      REAL(rprec) :: x(n)
      REAL(rprec), INTENT(in) :: fvec(m)
      REAL(rprec) :: fjac(ldfjac,n)
      REAL(rprec) :: fnorm_min, x_min(n), wa(m)
      EXTERNAL fcn
!DEC$ IF DEFINED (MPI_OPT)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, j, iproc_min, num_split, joff
      INTEGER :: iflag_array(1)
      REAL(rprec) :: eps, epsmch, dpmpar, temp, enorm, cur_norm
      REAL(rprec), PARAMETER :: one = 1, zero = 0
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: buffer, h,
     1     fnorm_array                                       !mpi stuff
      REAL(rprec), ALLOCATABLE :: fjac_overflow(:,:), fnorm_overflow(:)
      INTEGER :: istat, ierr_flag, ikey
      LOGICAL :: lflip
      LOGICAL, ALLOCATABLE     :: flip_overflow(:)
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL dpmpar, enorm
C-----------------------------------------------
c
c     SUBROUTINE fdjac2
c
c     this SUBROUTINE computes a forward-difference approximation
c     to the m by n jacobian matrix ASSOCIATED with a specified
c     problem of m functions in n variables.
c
c     the SUBROUTINE statement is
c
c       SUBROUTINE fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)
c
c     WHERE
c
c       fcn is the name of the user-supplied SUBROUTINE which
c         calculates the functions. fcn must be declared
c         in an EXTERNAL statement in the user calling
c         program, and should be written as follows.
c
c         SUBROUTINE fcn(m,n,x,fvec,iflag)
c         INTEGER m,n,iflag
c         REAL(rprec) x(n),fvec(m)
c         ----------
c         calculate the functions at x and
c         RETURN this vector in fvec.
c         ----------
c         RETURN
c         END
c
c         the value of iflag should not be changed by fcn unless
c         the user wants to terminate execution of fdjac2.
c         in this CASE set iflag to a negative INTEGER.
c
c       m is a positive INTEGER input variable set to the number
c         of functions.
c
c       n is a positive INTEGER input variable set to the number
c         of variables. n must not exceed m.
c
c       x is an input array of length n.
c
c       fvec is an input array of length m which must contain the
c         functions evaluated at x.
c
c       fjac is an output m by n array which contains the
c         approximation to the jacobian matrix evaluated at x.
c
c       ldfjac is a positive INTEGER input variable not less than m
c         which specifies the leading DIMENSION of the array fjac.
c
c       iflag is an INTEGER variable which can be used to terminate
c         the execution of fdjac2. see description of fcn.
c
c       epsfcn is an input variable used in determining a suitable
c         step length for the forward-difference approximation. this
c         approximation assumes that the relative errors in the
c         functions are of the order of epsfcn. If epsfcn is less
c         than the machine precision, it is assumed that the relative
c         errors in the functions are of the order of the machine
c         precision.
c
c       wa is a work array of length m.
c
c     subprograms called
c
c       user-supplied ...... fcn
c
c       MINpack-supplied ... dpmpar
c
c       fortran-supplied ... abs,max,sqrt
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      ALLOCATE (fnorm_array(n), h(n), buffer(m), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error in fdjac2_mp'
c
c     epsmch is the machine precision.
c
      epsmch = dpmpar(1)
c
      eps = SQRT(MAX(epsfcn,epsmch))
      cur_norm = enorm(m,fvec)

      h(1:n) = eps*ABS(x(1:n))
      WHERE (h .eq. zero) h=eps
      WHERE (flip) h = -h
      ierr_flag = 0


!     NOTE: if numprocs > n, num_split = 0 and we go to the overflow 
!           ("left-over") loop for num_split=mod(n,numprocs) == n
      num_split = n/numprocs
      MPI_COMM_WORKERS = MPI_COMM_STEL
      joff = 1

      DO i = 1, num_split
         j    = myid + joff
         temp = x(j)
         x(j) = temp + h(j)

         iflag = j
         CALL fcn(m, n, x, wa, iflag, ncnt)

         x(j) = temp
         buffer(:m) = (wa(:m) - fvec(:m))/h(j)

!
!        STORE FNORM OF PERTURBED STATE (X + H)
!
         temp = enorm(m,wa)
         IF (temp > cur_norm) lflip = .not. flip(j)

!        WRITE ITERATION, PROCESSOR, CHISQ TO SCREEN
         WRITE (6, '(2x,i6,8x,i3,7x,1es12.4)') ncnt+j,
     1         myid+1, temp**2
         CALL flush(6)
      
         CALL MPI_ALLGATHER(buffer, m, MPI_REAL8, 
     1        fjac(1,joff), m, MPI_REAL8, MPI_COMM_STEL, ierr_mpi)
         CALL MPI_ALLGATHER(temp, 1, MPI_REAL8, 
     1        fnorm_array(joff), 1, MPI_REAL8, MPI_COMM_STEL, ierr_mpi)
         CALL MPI_ALLGATHER(lflip, 1, MPI_LOGICAL, 
     1        flip(joff), 1, MPI_LOGICAL, MPI_COMM_STEL, ierr_mpi)

         joff = joff + numprocs
 
         IF (ierr_mpi .ne. 0) GO TO 100
      END DO

!
!     account for left-over variables OR special case when numprocs>n: 
!     must redefine the MPIC_COMM_WORDERS group to consist of MOD(n,numprocs) only
!
      num_split = MOD(n, numprocs)          !overflow...
      IF (num_split .ne. 0) THEN
      
         ALLOCATE(fjac_overflow(m,numprocs), fnorm_overflow(numprocs),
     1            flip_overflow(numprocs), stat = istat)
         IF (istat .ne. 0) THEN
            PRINT *,' istat = ',istat,' in fdjac_mod, myid: ', myid
            STOP
         END IF
         j    = myid + joff
         IF (j .le. n) THEN
            temp = x(j)
            x(j) = temp + h(j)
            ikey = WORKER_SPLIT_KEY+1
         ELSE
            temp = 0
            ikey = MPI_UNDEFINED
         END IF

!     Set new WORKERS communicator to include first MOD(n,numprocs) only
         CALL MPI_COMM_SPLIT(MPI_COMM_STEL, ikey, worker_id, 
     1                       MPI_COMM_WORKERS, ierr_mpi)
         IF (ierr_mpi .ne. 0) 
     1       STOP 'MPI_COMM_SPLIT error in fdjac_mp_queue'


         IF (j .le. n) THEN            !(MPI_COMM_WORKERS ONLY ACCESS HERE)
            iflag = j
            CALL fcn(m, n, x, wa, iflag, ncnt)
            x(j) = temp
            buffer(:m) = (wa(:m) - fvec(:m))/h(j)
            temp = enorm(m,wa)
            IF (temp > cur_norm) lflip = .not. flip(j)

            WRITE (6, '(2x,i6,8x,i3,7x,1es12.4)') ncnt+j,
     1             myid+1, temp**2
      
            CALL MPI_COMM_FREE(MPI_COMM_WORKERS, ierr_mpi)
         ELSE                          !(ALL OTHER PROCESSORS GO HERE)
            buffer = 0; lflip = .true.
         END IF

         CALL MPI_ALLGATHER(buffer, m, MPI_REAL8, 
     1        fjac_overflow, m, MPI_REAL8, MPI_COMM_STEL, ierr_mpi)
         CALL MPI_ALLGATHER(temp, 1, MPI_REAL8, 
     1        fnorm_overflow, 1, MPI_REAL8, MPI_COMM_STEL, ierr_mpi)
         CALL MPI_ALLGATHER(lflip, 1, MPI_LOGICAL, 
     1        flip_overflow, 1, MPI_LOGICAL, MPI_COMM_STEL, ierr_mpi)

         fjac(:,joff:n)      = fjac_overflow(:,1:num_split)
         fnorm_array(joff:n) = fnorm_overflow(1:num_split)
         flip(joff:n)        = flip_overflow(1:num_split)

         DEALLOCATE(fjac_overflow, fnorm_overflow, flip_overflow)
 
      END IF
      
!
!     Find processor with minimum fnorm_min value and broadcast wa (=fvec_min), x_min, fnorm_min
!
      fnorm_min   = MINVAL(fnorm_array)
      iflag_array = MINLOC(fnorm_array)
      iproc_min = iflag_array(1)
      IF (iproc_min .le. 0 .or. iproc_min .gt. n) THEN
         PRINT *,'IPROC_MIN=',iproc_min,' out of range in fdjac2_mp'
         STOP
      END IF
      wa(:) = fjac(:, iproc_min)*h(iproc_min) + fvec(:)
      x_min(:) = x(:)
      x_min(iproc_min) = x(iproc_min) + h(iproc_min)

      DEALLOCATE (h, fnorm_array, buffer)
!
!     Do any special cleanup now for iflag = flag_cleanup
!     Barrier appears in fcn call
!
      ikey = flag_cleanup
      MPI_COMM_WORKERS = MPI_COMM_STEL
      CALL fcn(m, n, x, wa, ikey, ncnt)

      RETURN

 100  CONTINUE
      WRITE (6, *) ' MPI_BCAST error in FDJAC2_MP: IERR=', ierr_mpi,
     1             ' PROCESSOR: ',myid

!DEC$ ENDIF
      END SUBROUTINE fdjac2_mp


      SUBROUTINE fdjac2_mp_queue(fcn, fnorm, m, n, x, fvec, fjac,
     1           ldfjac, iflag, ncnt, epsfcn, fnorm_min, x_min, wa,
     2           fnorm_array,clean_flag)      !PPPL (adds fnorm)
      USE stel_kinds
      USE mpi_params
      USE safe_open_mod, ONLY: safe_open
      USE mpi_inc
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: m, n, ldfjac
      INTEGER, INTENT(inout) :: ncnt
      INTEGER :: iflag
      REAL(rprec), INTENT(in) :: epsfcn, fnorm
      REAL(rprec) :: x(n)
      REAL(rprec), INTENT(in) :: fvec(m)
      REAL(rprec) :: fjac(ldfjac,n)
      REAL(rprec) :: fnorm_min, x_min(n), wa(m)
      REAL(rprec) :: fnorm_array(n)
      logical     :: clean_flag
      EXTERNAL fcn
!DEC$ IF DEFINED (MPI_OPT)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      logical, dimension(n) :: lmask
      INTEGER :: i, j, iunit,iproc_min, ix_temp
      integer, dimension(1) :: isort
      integer, dimension(numprocs) :: iflag_array
      REAL(rprec) :: eps, epsmch, dpmpar, temp, enorm, temp_norm
      REAL(rprec), PARAMETER :: one = 1, zero = 0
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: h
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: x_global,fvec_array
      CHARACTER(16) ::  temp_string
      CHARACTER(256) ::  jac_file
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL dpmpar, enorm
C-----------------------------------------------
c
c     SUBROUTINE fdjac2
c
c     this SUBROUTINE computes a forward-difference approximation
c     to the m by n jacobian matrix ASSOCIATED with a specified
c     problem of m functions in n variables.
c
c     the SUBROUTINE statement is
c
c       SUBROUTINE fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)
c
c     WHERE
c
c       fcn is the name of the user-supplied SUBROUTINE which
c         calculates the functions. fcn must be declared
c         in an EXTERNAL statement in the user calling
c         program, and should be written as follows.
c
c         SUBROUTINE fcn(m,n,x,fvec,iflag)
c         INTEGER m,n,iflag
c         REAL(rprec) x(n),fvec(m)
c         ----------
c         calculate the functions at x and
c         RETURN this vector in fvec.
c         ----------
c         RETURN
c         END
c
c         the value of iflag should not be changed by fcn unless
c         the user wants to terminate execution of fdjac2.
c         in this case set iflag to a negative integer.
c
c       m is a positive INTEGER input variable set to the number
c         of functions.
c
c       n is a positive INTEGER input variable set to the number
c         of variables. n must not exceed m.
c
c       x is an input array of length n.
c
c       fvec is an input array of length m which must contain the
c         functions evaluated at x.
c
c       fjac is an output m by n array which contains the
c         approximation to the jacobian matrix evaluated at x.
c
c       ldfjac is a positive INTEGER input variable not less than m
c         which specifies the leading DIMENSION of the array fjac.
c
c       iflag is an INTEGER variable which can be used to terminate
c         the execution of fdjac2. see description of fcn.
c
c       epsfcn is an input variable used in determining a suitable
c         step length for the forward-difference approximation. this
c         approximation assumes that the relative errors in the
c         functions are of the order of epsfcn. If epsfcn is less
c         than the machine precision, it is assumed that the relative
c         errors in the functions are of the order of the machine
c         precision.
c
c       wa is a work array of length m.
c
c     subprograms called
c
c       user-supplied ...... fcn
c
c       MINpack-supplied ... dpmpar
c
c       fortran-supplied ... abs,max,sqrt
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
c*************** Updated by SAL****************************************
      ALLOCATE (x_global(n,n),h(n),fvec_array(m,n))

      epsmch = dpmpar(1)
      eps = SQRT(MAX(epsfcn,epsmch))
      if( epsfcn < 0 ) eps = - eps

      h = eps*ABS(x)
      WHERE (h .le. eps/100.) h = eps/100.
      WHERE (h .eq. zero) h=eps
      WHERE (flip) h = -h
      h_order = h

      DO i = 1, n
         x_global(:,i) = x(:)
         x_global(i,i) = x(i) + h(i)
      END DO

!
!     Evaluate finite difference jacobian
!
      CALL eval_x_queued(fcn,m,n,n,x_global,fvec_array,ncnt,
     1                   MPI_COMM_STEL)
      CALL MPI_BCAST(fvec_array,m*n,MPI_DOUBLE_PRECISION,
     1               master,MPI_COMM_STEL,ierr_mpi)

!
!     Calculate Jacobian
!
      fnorm_array = SQRT(SUM(fvec_array*fvec_array,DIM=1))
      DO i = 1, n
         fjac(:,i) = (fvec_array(:,i) - fvec(:))/h(i)
         temp_norm = fnorm_array(i)*fnorm_array(i)
         IF ((temp_norm >= 1.0E12).or.(temp_norm/=temp_norm)) THEN
            jac_err(i) = 0
            fjac(:m,i) = 0.0
            flip(i) = .not. flip(i)
         ELSE
            jac_err(i) = 1
            IF (temp_norm > fnorm) flip(i) = .not. flip(i)
            END IF
      END DO

!
!     Check to make sure there is a search direction (PPPL)
!
      IF (ALL(fnorm_array >= 1.0E12)) iflag = 327
      IF (iflag == 327) RETURN

!
!     Output Function Evaluations, Jacobian, and reorder jacobian
!
      IF (myid .eq. master) THEN
         ! ADDED by S. Lazerson (dump feval then Jacobian to file)
         IF (ncnt == 0) THEN
            WRITE(temp_string,'(i6.6)') ncnt
         ELSE
            WRITE(temp_string,'(i6.6)') ncnt+1
         END IF
         jac_file = 'fevals.' // TRIM(temp_string)
         iunit = 27; j=0;
         CALL safe_open(iunit,j,TRIM(jac_file),'new','formatted')
         WRITE(iunit,'(2X,i6,2X,i6)') m,n
         WRITE(iunit,'(1p,4e22.14)') (fvec(i), i=1,m)
         DO i = 1, m
            WRITE(iunit,'(1p,4e22.14)') (fvec_array(i,j), j=1,n)
         END DO
         CLOSE(iunit)
         jac_file = ''
         jac_file = 'jacobian.' // TRIM(temp_string)
         CALL safe_open(iunit,j,TRIM(jac_file),'new','formatted')
         WRITE(iunit,'(2X,i6,2X,i6)') m,n
         DO i = 1, m
            WRITE(iunit,'(1p,4e22.14)') (fjac(i,j), j=1,n)
         END DO
         CLOSE(iunit)
      END IF

!
!     Reorder jacobian
!
      DO i = 1, n
         jac_index(i) = i
      END DO
      i = 1
      DO WHILE(i <= COUNT(jac_err > 0))
         IF (jac_err(i) == 0) THEN
            DO j = i, n-1
               fjac(1:m,j) = fjac(1:m,j+1)
               jac_index(j) = jac_index(j+1)
               fnorm_array(j) = fnorm_array(j+1)
               jac_err(j) = jac_err(j+1)
               h_order(j) = h_order(j+1)
            END DO
            fjac(1:m,n) = 0.0
            fnorm_array(n) = 1.0E30
            jac_err(n) = 0
            jac_index(n) = 0
            h_order(n)   = 0.0
         ELSE
            i = i + 1
         END IF
      END DO
      n_red = COUNT(jac_err > 0)

!
!     Find minimum
!
      jac_order = 0
      ix_temp = 1
      temp = 0
      lmask = .true.
      DO WHILE (ix_temp <= n_red)
         isort = MINLOC(fnorm_array, MASK=lmask)
         temp = fnorm_array(isort(1))
         jac_order(ix_temp) = isort(1)
         IF(isort(1) <= 0 .or. isort(1) > n_red) THEN
            EXIT
         ELSE IF(fnorm_array(isort(1)) > fnorm) THEN
            EXIT
         ELSE
            lmask(isort(1)) = .false.
            ix_temp = ix_temp + 1
         END IF
      END DO
      ix_min = jac_order(1)
      IF(ix_temp <= n) jac_order(ix_temp:) = 0
      jac_count = ix_temp - 1
      IF (ix_min .le. 0 .or. ix_min .gt. n) THEN
         PRINT *,' IX_MIN = ',ix_min,' out of range'
         STOP
      END IF
      j = jac_index(ix_min)
      fnorm_min = fnorm_array(ix_min)
      wa(:) = fjac(:, ix_min)*h(j) + fvec(:)                     ! Note wa ~ fvec_min
      x_min(:) = x(:)
      x_min(j) = x(j) + h(j)

!
!     Cleanup
!
      i = flag_cleanup_jac
      if (clean_flag) CALL fcn(m, n, x, wa, i, ncnt)
      DEALLOCATE(x_global,h,fvec_array)
      
      RETURN
c**********************************************************************
!DEC$ ENDIF
      END SUBROUTINE fdjac2_mp_queue


      SUBROUTINE fdjac2(fcn, m1, n1, x, fvec, fjac, ldfjac, iflag,
     1    ncnt1, epsfcn, wa, time, fnorm_min, x_min, fvec_min)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: m1, n1, ldfjac, ncnt1
      INTEGER, TARGET :: iflag
      REAL(rprec) ::  epsfcn, time
      REAL(rprec), TARGET :: x(n), wa(m)
      REAL(rprec), INTENT(in) :: fvec(m)
      REAL(rprec), INTENT(out) :: fjac(ldfjac,n)
      REAL(rprec), INTENT(out) :: fnorm_min, x_min(n), fvec_min(m)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: j, istat, iread, ic1, ic2, irate, count_max
      REAL(rprec) :: eps, epsmch, h, dpmpar, temp, cur_norm
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL fcn, dpmpar, multiprocess, fdjac_parallel
      REAL(rprec), EXTERNAL :: enorm
C-----------------------------------------------
      fjac = 0
      fnorm_min = 0
      x_min = 0
      fvec_min = 0
!DEC$ IF .NOT.DEFINED (MPI_OPT)
c
c     SUBROUTINE fdjac2
c
c     this SUBROUTINE computes a forward-difference approximation
c     to the m by n jacobian matrix ASSOCIATED with a specified
c     problem of m functions in n variables.
c
c     Here
c
c       fcn is the name of the user-supplied SUBROUTINE which
c         calculates the functions. fcn must be declared
c         in an external statement in the user calling
c         program (see LMDIF1 for documentation), and should be written as follows:
c
c         SUBROUTINE fcn(m,n,x,fvec,iflag,ncnt)
c         INTEGER m,n,iflag
c         REAL(rprec) x(n),fvec(m)
c         ----------
c         calculate the functions at x and
c         RETURN this vector in fvec.
c         ----------
c         RETURN
c         END
c
c       fjac is an output m by n array which CONTAINS the
c         approximation to the jacobian matrix evaluated at x.
c
c       ldfjac is a positive INTEGER input variable not less than m
c         which specifies the leading DIMENSION of the array fjac.
c
c       iflag is an INTEGER variable which can be used to terminate
c         the execution of fdjac2. see description of fcn.
c
c       epsfcn is an input variable used in determining a suitable
c         step length for the forward-difference approximation. this
c         approximation assumes that the relative errors in the
c         functions are of the order of epsfcn. IF epsfcn is less
c         than the machine precision, it is assumed that the relative
c         errors in the functions are of the order of the machine
c         precision.
c
c       wa is a work array of length m.
c
c     subprograms called
c
c       user-supplied ...... fcn
c
c       MINpack-supplied ... dpmpar
c
c       fortran-supplied ... ABS,max,sqrt
c
c     argonne national laboratory. MINpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********

c
c     epsmch is the machine precision.
c
      fjac = 0
      fnorm_min = 0
      x_min = 0
      fvec_min = 0
      epsmch = dpmpar(1)
c
      eps = SQRT(MAX(epsfcn,epsmch))
!
!     Load fdjac_mod values. pointers will automatically update targets...
!     Prepare for multi-processing...
!
      m = m1
      n = n1
      ncnt = ncnt1
      xp => x
      wap => wa

!     Find min chisq = fnorm**2 state for this jacobian evaluation
!     (Do NOT retain from previous evaluations, or could get into a non-converging loop...)

      fnorm_min = HUGE(fnorm_min)

      CALL SYSTEM_CLOCK(ic1, irate)
      CALL multiprocess(n, max_processors, fdjac_parallel, fcn)
      CALL SYSTEM_CLOCK(ic2, irate, count_max)
      IF (ic2 .lt. ic1) ic2 = ic2 + count_max

c     DO j = 1, n
c         temp = x(j)
c         h = eps*ABS(temp)
c         IF (h .eq. zero) h = eps
c         x(j) = temp + h
c         CALL fcn (m, n, x, wa, iflag, ncnt)
c         IF (iflag .lt. 0) EXIT
c         x(j) = temp
c         fjac(:m,j) = (wa - fvec)/h
c      END DO


      cur_norm = enorm(m,fvec)      ! where are we now?

      DO j = 1, n

        READ (j+1000, iostat=iread) istat, iflag, h, temp
        IF (iread .ne. 0) THEN
           WRITE (6, *) 'Error reading from file fort.', j+1000,
     1     ' in fdjac2: IOSTAT = ', iread
           iflag = -14
        ELSE IF (j .ne. istat) THEN
           WRITE (6, *) 'Wrong value for INDEX j READ in fdjac2'
           iflag = -14
        END IF

        IF (iflag .ne. 0) EXIT

!DEC$ IF DEFINED (CRAY)
        DO k = 1, m
           READ (j+1000) wa(k)
        END DO
!DEC$ ELSE
        READ (j+1000) wa
!DEC$ ENDIF
        fjac(:m,j) = (wa - fvec)/h

        IF( temp > cur_norm) flip(j) = .not. flip(j)  ! flip for next time

        IF( temp < fnorm_min) THEN
           fnorm_min = temp
           fvec_min = wa
!DEC$ IF DEFINED (CRAY)
           DO k = 1, n
              READ (j+1000) x_min(k)
           END DO
!DEC$ ELSE
           READ (j+1000) x_min
!DEC$ ENDIF
        END IF

        CLOSE (j+1000, status='delete')                        !!Needed to run correctly in multi-tasking...

      END DO

!
!     Do any special cleanup now for iflag = flag_cleanup
!
      iflag = flag_cleanup
      CALL fcn(m, n, x, wa, iflag, ncnt)

      time = time + REAL(ic2 - ic1)/REAL(irate)                !!Time in multi-process CALL
!DEC$ ENDIF
      END SUBROUTINE fdjac2


      END MODULE fdjac_mod

