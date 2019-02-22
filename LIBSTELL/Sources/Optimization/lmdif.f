      SUBROUTINE lmdif(fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev,
     1   epsfcn, diag, mode, factor, nprint, info, nfev, fjac,
     2   ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4, xvmin, xvmax)
      USE stel_kinds
      USE lmpar_mod, fjac_mod=>fjac, ldfjac_mod=>ldfjac,
     1   ipvt_mod=>ipvt, qtf_mod=>qtf, diag_mod=>diag
!DEC$ IF DEFINED (MPI_OPT)
      USE fdjac_mod, ONLY: flip,flag_singletask,flag_cleanup,
     1                     fdjac2_mp_queue, jac_order, jac_count,
     2                     ix_min, h_order, flag_cleanup_lev,
     3                     jac_err, jac_index, n_red !PPPL   
!DEC$ ELSE
      USE fdjac_mod, ONLY: max_processors, flip, flag_singletask, 
     1                     flag_cleanup, fdjac2, jac_order, h_order,
     2                     jac_err, jac_index, flag_cleanup_lev,
     3                     ix_min, jac_count, n_red
!DEC$ ENDIF
      USE mpi_params
      USE safe_open_mod
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      INCLUDE 'mpif.h'                                       !mpi stuff
!DEC$ ENDIF
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: m, n, maxfev, mode, nprint, info, ldfjac
      INTEGER, INTENT(inout)  :: nfev
      REAL(rprec), INTENT(in) ::  ftol, xtol, gtol, factor
      REAL(rprec), INTENT(inout) :: epsfcn
      REAL(rprec), DIMENSION(n) :: x, wa1, wa2, wa3
      REAL(rprec), DIMENSION(m) :: fvec, wa4
      INTEGER, DIMENSION(n), TARGET :: ipvt
      REAL(rprec), DIMENSION(n), TARGET :: diag, qtf
      REAL(rprec), DIMENSION(ldfjac,n), TARGET :: fjac
      REAL(rprec), DIMENSION(n), INTENT(in), OPTIONAL :: xvmin
      REAL(rprec), DIMENSION(n), INTENT(in), OPTIONAL :: xvmax
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1,
     1   p1=0.1_dp, p5=0.5_dp, p25=0.25_dp, p75=0.75_dp, p0001=1.e-4_dp
      CHARACTER(LEN=130), DIMENSION(0:13) :: info_array 
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: iflag, iter, j, l, istat, ikey, lev_step_range,iunit         !PPPL
!      INTEGER :: iflag, iter, j, l, istat, ikey
      INTEGER :: cycle_count, subcycle,i,num_jacs
      REAL(rprec) :: actred, dirder, epsmch, fnorm,
     1   gnorm, prered, ratio, sum0, temp,
     2   temp1, temp2, xnorm, delta_old, actred_lev, par_old,fnorm_old  !PPPL
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: fjac_save
!DEC$ IF .NOT.DEFINED (MPI_OPT)
      REAL(rprec) :: wall_time, wall_time_lev
!DEC$ ENDIF
      REAL(rprec) :: fnorm_min, epsfcn0, epsfcn_temp
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: x_min, fvec_min
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: fnorm_array
      logical :: step_improved, first_jacobian, lsmall_step, lredo_diag !PPPL
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL fcn
      REAL(rprec), EXTERNAL :: dpmpar, enorm
C-----------------------------------------------
c
c     SUBROUTINE lmdif
c
c     the purpose of lmdif is to minimize the sum of the squares of
c     m nonlinear functions in n variables by a modification of
c     the levenberg-marquardt algorithm. the user must provide a
c     subroutine which calculates the functions. the jacobian is
c     then calculated by a forward-difference approximation.
c
c     the subroutine statement is
c
c       subroutine lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,
c                        diag,mode,factor,nprint,info,nfev,fjac,
c                        ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
c
c     where
c
c       fcn is the name of the user-supplied subroutine which
c         calculates the functions. fcn must be declared
c         in an external statement in the user calling. see lmdif1 for
c         documentation and should be written as follows.
c
c         subroutine fcn(m, n, x, fvec, iflag, ncnt)
c         integer m,n,iflag
c         real(rprec) x(n),fvec(m)
c         ----------
c         calculate the functions at x and
c         return this vector in fvec.
c         ----------
c         return
c         end
c
c
c       ftol is a nonnegative input variable. termination
c         occurs when both the actual and predicted relative
c         reductions in the sum of squares are at most ftol.
c         therefore, ftol measures the relative error desired
c         in the sum of squares.
c
c       xtol is a nonnegative input variable. termination
c         occurs when the relative error between two consecutive
c         iterates is at most xtol. therefore, xtol measures the
c         relative error desired in the approximate solution.
c
c       gtol is a nonnegative input variable. termination
c         occurs when the cosine of the angle between fvec and
c         any column of the jacobian is at most gtol in absolute
c         value. therefore, gtol measures the orthogonality
c         desired between the function vector and the columns
c         of the jacobian.
c
c       maxfev is a positive integer input variable. termination
c         occurs when the number of calls to fcn is at least
c         maxfev by the end of an iteration.
c
c       epsfcn is an input variable used in determining a suitable
c         step length for the forward-difference approximation. this
c         approximation assumes that the relative errors in the
c         functions are of the order of epsfcn. if epsfcn is less
c         than the machine precision, it is assumed that the relative
c         errors in the functions are of the order of the machine
c         precision.
c
c       diag is an array of length n. if mode = 1 (see
c         below), diag is internally set. if mode = 2, diag
c         must contain positive entries that serve as
c         multiplicative scale factors for the variables.
c
c       mode is an integer input variable. if mode = 1, the
c         variables will be scaled internally. if mode = 2,
c         the scaling is specified by the input diag. other
c         values of mode are equivalent to mode = 1.
c
c       factor is a positive input variable used in determining the
c         initial step bound. this bound is set to the product of
c         factor and the euclidean norm of diag*x if nonzero, or else
c         to factor itself. in most cases factor should lie in the
c         interval (.1,100.). 100. is a generally recommended value.
c
c       nprint is an integer input variable that enables controlled
c         printing of iterates if it is positive. in this case,
c         fcn is called with iflag = 0 at the beginning of the first
c         iteration and every nprint iterations thereafter and
c         immediately prior to return, with x and fvec available
c         for printing. if nprint is not positive, no special calls
c         of fcn with iflag = 0 are made.
c
c       info is an integer output variable. if the user has
c         terminated execution, info is set to the (negative)
c         value of iflag. see description of fcn. otherwise,
c         info is set
c       nfev is an integer output variable set to the number of
c         calls to fcn.
c
c       fjac is an output m by n array. the upper n by n submatrix
c         of fjac contains an upper triangular matrix r with
c         diagonal elements of nonincreasing magnitude such that
c
c                t     t           t
c               p *(jac *jac)*p = r *r,
c
c         where p is a permutation matrix and jac is the final
c         calculated jacobian. column j of p is column ipvt(j)
c         (see below) of the identity matrix. the lower trapezoidal
c         part of fjac contains information generated during
c         the computation of r.
c
c       ldfjac is a positive integer input variable not less than m
c         which specifies the leading dimension of the array fjac.
c
c       ipvt is an integer output array of length n. ipvt
c         defines a permutation matrix p such that jac*p = q*r,
c         where jac is the final calculated jacobian, q is
c         orthogonal (not stored), and r is upper triangular
c         with diagonal elements of nonincreasing magnitude.
c         column j of p is column ipvt(j) of the identity matrix.
c
c       qtf is an output array of length n which contains
c         the first n elements of the vector (q transpose)*fvec.
c
c       wa1, wa2, and wa3 are work arrays of length n.
c
c       wa4 is a work array of length m.
c
c     subprograms called
c
c       user-supplied ...... fcn
c
c       minpack-supplied ... dpmpar,enorm,fdjac2,lmpar,qrfac
c
c       fortran-supplied ... ABS,max,min,sqrt,mod
c
c     argonne national laboratory. MINpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********

!DEC$ IF DEFINED (MPI_OPT)
c     Get mpi parameters
      CALL MPI_COMM_RANK (MPI_COMM_STEL, myid, ierr_mpi)       !mpi stuff
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
      IF (ierr_mpi .ne. 0) STOP 'MPI_COMM_RANK error in LMDIF'
      CALL MPI_COMM_SIZE (MPI_COMM_STEL, numprocs, ierr_mpi)   !mpi stuff
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
      IF (ierr_mpi .ne. 0) STOP 'MPI_COMM_SIZE error in LMDIF'
!DEC$ ENDIF

      !info = 0;      iflag = 0;      nfev = 0;      cycle_count = 0
      info = 0;      iflag = 0;      cycle_count = 0

      info_array(9) = 
     1"improper input parameters (number constraints MUST " //
     1"be greater than number variables (>0)."

      info_array(1) = 
     1"algorithm estimates that the relative error in the " //
     1"sum of squares is at most ftol."

      info_array(2) = 
     2"algorithm estimates that the relative error between" //
     2" x and the solution is at most xtol."

      info_array(3) = 
     3"algorithm estimates that the relative error in the " //
     3"sum of squares and between x and the solution is at" //
     3" most ftol and xtol."

      info_array(4) = 
     4"the cosine of the angle between fvec and any column" //
     4" of the jacobian is at most gtol in absolute value."

      info_array(5) = 
     5"number of calls to fcn has reached or exceeded maxfev."

      info_array(6) = 
     6"ftol is too small. no further reduction in the sum " //
     6"of squares is possible."

      info_array(7) = 
     7"xtol is too small. no further improvement in the " // 
     7"approximate solution x is possible."

      info_array(8) = 
     8"gtol is too small. fvec is orthogonal to the columns" //
     8" of the jacobian to machine precision."

      info_array(0) = 
     9"levenberg-marquardt optimizer terminated properly." 

      info_array(10) = 
     A"Levenberg-Marquardt terminated prematurely due to " //
     A"function evaluation error." 

      info_array(11) = 
     A"Levenberg-Marquardt terminated prematurely: " //
     A"must request more than one (1) processor" 

      info_array(12) = 
     A"Levenberg-Marquardt terminated prematurely: " //
     A"because of bad Jacobian" 

      info_array(13) = 
     1"improper input parameters " //
     1"X_MIN < X_MAX"


!DEC$ IF DEFINED (MPI_OPT)
      IF (numprocs > n) THEN
         IF (myid .eq. master) THEN
            WRITE (6, *)'Warning: more processors have been requested',
     1      ' than the maximum (nvar) required = ',n
         END IF
      ELSE IF (numprocs <= 1) THEN
         info = 11
         GOTO 400
      END IF
!DEC$ ENDIF

!
!     check the input parameters for errors.
!
      ALLOCATE (x_min(n), fvec_min(m), flip(n),jac_order(n),
     1          fnorm_array(n),h_order(n),jac_err(n),jac_index(n),
     2          stat=istat)             
      IF (istat .ne. 0) STOP 'Allocation error in lmdif'

!
!     epsmch is the machine precision.
!     flip is control for direction flipping in fdjac2!
      epsmch = dpmpar(1)
      flip = .false.
      lsmall_step = .false.  !SAL for EPSMCH step
      lredo_diag = .false.    !SAL recompute DIAG scaling
      jac_order = 0                                                     !PPPL
!      info = 0; iflag = 0; nfev = 0; cycle_count = 0;                   !PPPL
      info = 0;      iflag = 0;      cycle_count = 0
      delta = 0                                                         !PPPL
!DEC$ IF .NOT.DEFINED (MPI_OPT)
      myid = 0;      wall_time = 0;      wall_time_lev = 0
!DEC$ ENDIF
!
!     ASSIGN VALUES TO MODULE VARIABLES (FACILITATES PASSING TO SUBROUTINES)
!
      ldfjac_mod = ldfjac
      ipvt_mod => ipvt
      fjac_mod => fjac
      diag_mod => diag
      qtf_mod => qtf
!
!     check the input parameters for errors.
!
      IF (n.le.0 .or. m.lt.n .or. ldfjac.lt.m .or. ftol.lt.zero
     1   .or. xtol.lt. zero .or. gtol.lt.zero .or. maxfev.le.0
     2   .or. factor.le.zero) THEN
         info = 9
         GOTO 400
      END IF

      IF (mode .eq. 2) THEN
         DO j = 1, n
            IF (diag(j) .le. zero) GOTO 300
         END DO
      END IF

      IF (PRESENT(xvmin) .AND. PRESENT(xvmax)) THEN
         DO j = 1, n
            IF (xvmin(j) .ge. xvmax(j)) info = 13
         END DO
         IF (info > 0) GOTO 400
      END IF

!     Set up workers communicator (only master processor here) for initial run
!DEC$ IF DEFINED (MPI_OPT)
      IF (myid .ne. master) THEN
         ikey = MPI_UNDEFINED
      ELSE
         ikey = WORKER_SPLIT_KEY+1
      END IF
      CALL MPI_COMM_SPLIT(MPI_COMM_STEL, ikey, worker_id, 
     1                    MPI_COMM_WORKERS, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF

!     evaluate the function at the starting point and calculate its norm.
      IF (myid .eq. master) THEN
         iflag = flag_singletask
         if (nfev .ne. 0) iflag = 0
         CALL fcn (m, n, x, fvec, iflag, nfev)
         IF (iflag .ne. 0) THEN
            WRITE(6,*) "FIRST RUN FAILS!  IMPROVE INPUT!"
            STOP "ERRROR!"
         END IF
         fnorm = enorm(m,fvec)
         iunit = 12; istat = 0
         CALL safe_open(iunit,istat,'xvec.dat','unknown','formatted',
     1                  ACCESS_IN='APPEND')
         WRITE(iunit,'(2(2X,I5.5))') n,0
         WRITE(iunit,'(10ES22.12E3)') x(1:n)
         WRITE(iunit,'(ES22.12E3)') fnorm
         CLOSE(iunit)
!DEC$ IF DEFINED (MPI_OPT)
         CALL MPI_COMM_FREE(MPI_COMM_WORKERS, ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)      
!DEC$ ENDIF
      END IF
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BCAST(x,n, MPI_REAL8, master, 
     1               MPI_COMM_STEL, ierr_mpi)
      !CALL MPI_BARRIER(MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF
      iflag = FLAG_CLEANUP
      IF (myid == master) iflag = flag_cleanup_lev
      call fcn (m, n, x, fvec, iflag, nfev)

!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BCAST(iflag,1, MPI_INTEGER, master, 
     1               MPI_COMM_STEL, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
      IF (iflag .ge. 0) CALL
     1     MPI_BCAST(fvec, m, MPI_REAL8, master, 
     2               MPI_COMM_STEL, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF

      IF (iflag .lt. 0) GOTO 300
      fnorm = enorm(m,fvec)
      IF (nfev.ge.maxfev .or. maxfev.eq.1) info = 5
      IF (info .ne. 0) GOTO 300
!
!     initialize levenberg-marquardt parameter (par) and iteration counter.
!
      par = 0
      iter = 1
      lfirst_lm = .true.
!DEC$ IF DEFINED (MPI_OPT)
      IF (myid .eq. master) WRITE (6, 1000) numprocs
 1000 FORMAT (/,' Beginning Levenberg-Marquardt Iterations',/,
     1        ' Number of Processors: ',i4,//,
!     1        ' Number of Processors: ',i4,' (1 controller proc)',//,
     2        70('='),/,2x,'Iteration',3x,'Processor',7x,'Chi-Sq',7x,
     3       'LM Parameter',6x,'Delta Tol'/,70('='))
!DEC$ ELSE
      WRITE (6, 1000) max_processors
 1000 FORMAT (/,' Beginning Levenberg-Marquardt Iterations',/,
     1        ' Number processors requested: ', i4,//,
     1        59('='),/,2x,'Iteration',8x,'Chi-Sq',7x,
     2       'LM Parameter',6x,'Delta Tol',/,59('='))
!DEC$ ENDIF
      IF (myid == master)  WRITE(6, '(2x,i6,8x,i3,7x,1es12.4)')
     1                           0, myid, fnorm*fnorm
      CALL FLUSH(6)

!
!     beginning of the outer loop.
!
      first_jacobian = .true.                                            !PPPL
      outerloop: DO WHILE (nfev .lt. maxfev)                         
         delta_old = delta                                               !PPPL
         par_old = par                                                   !PPPL
         fnorm_old = fnorm
!DEC$ IF DEFINED (MPI_OPT)
         CALL MPI_BARRIER(MPI_COMM_STEL,ierr_mpi)                       !PPPL -SAL
         IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF
         
!
!        calculate the jacobian matrix.
!
!DEC$ IF DEFINED (MPI_OPT)
         CALL fdjac2_mp_queue(fcn, fnorm, m, n, x, fvec, fjac, ldfjac,
     1           iflag, nfev, epsfcn, fnorm_min, x_min, fvec_min,
     2           fnorm_array,.true.)        !PPPL
     
!DEC$ ELSE
         iflag = 2
         CALL fdjac2(fcn, m, n, x, fvec, fjac, ldfjac, iflag, nfev,
     1               epsfcn, wa4, wall_time, fnorm_min, x_min, fvec_min)
!DEC$ ENDIF
         nfev = nfev + n
         
         ! The Jacobian evaluation was all failed evaluations
         IF (iflag == 327) THEN
            epsfcn = epsfcn*0.1
            iflag = 0
            IF (myid == master) WRITE(6,*) 
     1             '  Improving Jacobian: Small stepsize', epsfcn
            CYCLE outerloop
         END IF
         IF (iflag .lt. 0) EXIT outerloop

!
!        compute the qr factorization of the jacobian.
!
!         CALL qrfac(m, n, fjac, ldfjac, .true., ipvt,
!     1              n, wa1, wa2, wa3)
         ipvt = 0
         wa1  = 0
         wa2  = 0
         wa3  = 0
         CALL qrfac(m, n_red, fjac, ldfjac, .true., ipvt,
     1              n_red, wa1, wa2, wa3)

!
!        on the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.
!        also on the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.
!
         IF (iter .eq. 1 .or. lredo_diag) THEN
            IF (mode .ne. 2) THEN
               diag = 1.0
               DO i = 1, n_red
                  j = jac_index(i)
                  diag(j) = wa2(i)
                  IF (wa2(i) == zero) diag(j) = one
               END DO
               !WHERE (wa2 .eq. zero) diag = one
               !IF (myid == master) THEN
               !   WRITE(6,*) '  ----- DIAG -----'
               !   DO i =1, n
               !      WRITE(6,*) i,diag(i)
               !   END DO
               !   WRITE(6,*) '  ----------------'
               !END IF
            END IF
            wa3 = diag*x
            xnorm = enorm(n,wa3)
            delta = factor*xnorm
            IF (delta .eq. zero) delta = factor
            lredo_diag = .false.
         END IF

!
!        form (q TRANSPOSE)*fvec and store the first n components in qtf.
!
        
         wa4 = fvec
         DO j = 1, n_red
            IF (fjac(j,j) .ne. zero) THEN
               sum0 = SUM(fjac(j:m,j)*wa4(j:m))
               temp = -sum0/fjac(j,j)
               wa4(j:m) = wa4(j:m) + fjac(j:m,j)*temp
            END IF
            fjac(j,j) = wa1(j)
            qtf(j) = wa4(j)
         END DO

!
!        compute the norm of the scaled gradient.
!
         gnorm = zero
         IF (fnorm .ne. zero) THEN
            DO j = 1, n_red
               l = ipvt(j)
               IF (wa2(l) .ne. zero) THEN
                  sum0 = SUM(fjac(1:j,j)*(qtf(1:j)/fnorm))
                  gnorm = MAX(gnorm,ABS(sum0/wa2(l)))
               END IF
            END DO
         END IF

!
!        test for convergence of the gradient norm.
!
         IF (gnorm .le. gtol) info = 4
         IF (info .ne. 0) EXIT outerloop

!
!        rescale if necessary.
!
         IF (mode .ne. 2) THEN
            DO i = 1, n_red
               j = jac_index(i)
               diag(j) = MAX(diag(j),wa2(i))
            END DO
         END IF
         !IF (mode .ne. 2) diag = MAX(diag,wa2)

!
!        If maxfev == 2 then we calculate the Jacobian and stop the code
!        10/15/12 - SAL (PPPL)
!
         IF (maxfev.eq.2) info = 5
         IF (info .ne. 0) GOTO 300

!
!        set up for inner loop (levmarqloop) to determine x update.
!
         subcycle = 0
         ratio = 0

         ALLOCATE (fjac_save(n,n), stat=istat)
         IF (istat .ne. 0) STOP 'Fjac_save allocation error'
         fjac_save(:n,:n) = fjac(:n,:n)

!
!        Levenburg Inner Loop
!                                     
         levmarqloop: DO WHILE (ratio.lt.p0001 .and. subcycle.lt.2)       

           subcycle = subcycle + 1
           fjac(:n,:n) = fjac_save(:n,:n)
           spread_ratio = abs(epsfcn/delta/10)                           !PPPL

!
!        Determine the levenberg-marquardt PARAMETER.
!
!DEC$ IF DEFINED (MPI_OPT)
           IF (PRESENT(xvmin) .and. PRESENT(xvmax)) THEN
              CALL levmarq_param_mp (x, wa1, wa2, wa3, wa4,
     1                               nfev, m, n, iflag, fcn, 
     2                               lev_step_range,fnorm_min,
     3                               xvmin,xvmax)   !PPPL
           ELSE
              CALL levmarq_param_mp (x, wa1, wa2, wa3, wa4,
     1                               nfev, m, n, iflag, fcn, 
     2                               lev_step_range,fnorm_min)   !PPPL
           END IF
!DEC$ ELSE
           CALL levmarq_param(x, wa1, wa2, wa3, wa4,
     1          wall_time_lev, nfev, m, n, iflag, fcn)
!DEC$ ENDIF
           IF (iflag .lt. 0) EXIT

!
!        compute the scaled actual reduction.
!
           actred = -one
           IF (p1*fnorm1 .lt. fnorm) actred = one - (fnorm1/fnorm)**2

!
!        compute the scaled predicted reduction (prered) and
!        the scaled directional derivative (dirder).
!
           DO j = 1, n_red
              wa3(j) = zero
              l = ipvt(j)
              i = jac_index(l)
              temp = wa1(l)
              wa3(1:j) = wa3(1:j) + fjac(1:j,j)*temp
           END DO
           temp1 = enorm(n_red,wa3)/fnorm
           temp2 = (SQRT(par)*pnorm)/fnorm
           prered = temp1**2 + temp2**2/p5
           dirder = -(temp1**2 + temp2**2)

!
!        compute the ratio of the actual to the predicted reduction.
!
           actred_lev = actred                                           !PPPL
           ratio = zero
           IF (prered .ne. zero) ratio = actred/prered

!
!        test if Jacobian calculation gave best improvement
!        (added by MZ and SH) (always false for non-mpi runs)
!           
!DEC$ IF DEFINED (MPI_OPT)
           step_improved = fnorm_min < MIN(fnorm,fnorm1)
!DEC$ ENDIF
           
           ! Jacobian gives best improvment so take orthagonal step
           IF (step_improved) THEN
              IF (myid == master)
     1           jac_count = MIN(count(jac_order(:) .ne. 0),
     2                           INT(SQRT(REAL(numprocs))))
!DEC$ IF DEFINED (MPI_OPT)
              call MPI_BCAST(jac_count, 1, MPI_INTEGER, master,
     1                       MPI_COMM_STEL, ierr_mpi)
              IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
              IF (ierr_mpi .ne. 0) GOTO 3000
!DEC$ ENDIF
              
              IF(myid .ne. master) jac_order = 0
              
!DEC$ IF DEFINED (MPI_OPT)
              CALL MPI_BCAST(jac_order, jac_count, MPI_INTEGER, master,
     1                       MPI_COMM_STEL, ierr_mpi)
              IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
              IF (ierr_mpi .ne. 0) GOTO 3000
              
              CALL MPI_BARRIER(MPI_COMM_STEL,ierr_mpi)                  !PPPL -SAL
              IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF
              
              CALL stepopt_mp(fcn, wa2, wa4, m, n, x_min, fvec_min,
     1                         fnorm_min, iflag, nfev, epsfcn)
              
              IF (myid .eq. master) THEN
                 WRITE(6,'(a,i6,a)')
     1           ' Using minimum from Jacobian/step improvement (',
     2           ix_min, ')'
              END IF
               
              wa2 = x_min
              wa4 = fvec_min
              fnorm1 = fnorm_min
              IF ( cycle_count > 1) THEN
                 delta = MAX(delta, delta_old)
                 par = MIN(par, par_old)
              END IF

              actred = 1 - (fnorm1/fnorm)**2
              ratio = actred/prered
              lredo_diag = .true.
              mode = 1
           ! If there was no direction of improvement found try
           ! recomputing jacobian with signs of displacements flipped
           ELSE IF (first_jacobian .and. actred < 0) THEN
              first_jacobian = .false.
              deallocate(fjac_save)
              lredo_diag = .true.
              mode = 1
              IF (myid == master) WRITE(6,*) 
     1             '  Improving Jacobian: flipping direction'
              cycle outerloop
           ! Try a very small step size
           ELSE IF (.not. lsmall_step .and. actred < 0) THEN
              epsfcn = epsfcn / 100.
              lsmall_step = .true.
              first_jacobian = .true.
              deallocate(fjac_save)
              lredo_diag = .true.
              mode = 1
              IF (myid == master) WRITE(6,*) 
     1             '  Improving Jacobian: Small stepsize'
              cycle outerloop
           ! In the multiprocessor case, if we have enough samples,
           ! defer to the best delta parameter found.  Because of the
           ! asymmetry in the 'factors' (in levmarg_param), we need
           ! to perturb ('goose') delta when it is at the top end of
           ! the range.
           ! Not sure if this works or not
!           ELSE IF (numprocs > 8) THEN
!              IF (actred < 0) THEN
!                 delta = 0.25 * MIN(delta, delta_old)
!                 par   = 4 * par
!              ELSE IF (lev_step_range == 1) THEN
!                 delta = 2 * delta
!                 par   = 0.5 * par
!              END IF
           ! Update the step bound, if Levenberg Step gave the best
           ! improvment.
           ELSE 
              IF (ratio .le. p25) THEN
                 IF (actred .ge. zero) THEN
                    temp = p5
                 ELSE
                    temp = p5*dirder/(dirder + p5*actred)
                 END IF
                 IF (p1*fnorm1.ge.fnorm .or. temp.lt.p1) temp = p1
                 delta = temp*MIN(delta,pnorm/p1)
!                 delta = max( temp*min(delta,pnorm/p1), p1*delta)       !PPPL
                 par = par/temp
              ELSE IF (par.eq.zero .or. ratio.ge.p75) THEN
                 delta = pnorm/p5
                 par = p5*par
              END IF
           END IF

!
!        test for successful iteration.
!        update x, fvec, and their norms.
!
           if( n_red > numprocs .and. subcycle < 2 .and.
     1         ratio .lt. p0001 ) cycle levmarqloop                     !PPPL
     
           IF (ratio .ge. p0001) THEN
              x = wa2
              ! New way to calculated xnorm
              wa2 = 0.0
              DO i = 1, n_red
                 j = jac_index(i)
                 wa2(j) = diag(j)*x(j)
              END DO
              !wa2 = diag*x
              fvec = wa4
              xnorm = enorm(n,wa2)
              fnorm = fnorm1
              iter = iter + 1
              IF (myid .eq. master) WRITE(6,'(/,3(2x,a,1es10.3)/)')
     1          'new minimum = ', fnorm**2,'lm-par = ', par,
     2          'delta-tol = ', delta
           END IF

           cycle_count = cycle_count + 1

!
!        tests for convergence.
!
           IF (ABS(actred).le.ftol .and. prered.le.ftol
     1        .and. p5*ratio.le.one) info = 1
!
!        next test made more stringent to avoid premature declaration of
!        completion observed on some problems.
!
           IF (delta.le.xtol*xnorm/10 .and. cycle_count>2
     1         .and. .not. step_improved ) info = 2                     !PPPL
           IF (ABS(actred).le.ftol .and. prered.le.ftol
     1        .and. p5*ratio.le.one .and. info.eq.2) info = 3
           IF (info .ne. 0) EXIT levmarqloop

!
!        tests for termination and stringent tolerances.
!
           IF (nfev .ge. maxfev) info = 5
           IF (ABS(actred).le.epsmch .and. prered.le.epsmch
     1        .and. p5*ratio.le.one) info = 6
           IF (delta .le. epsmch*xnorm) info = 7
           IF (gnorm .le. epsmch) info = 8
           IF (info .ne. 0) EXIT levmarqloop
!
!        END of the inner loop. repeat if iteration unsuccessful (ratio < p0001)
!
         END DO levmarqloop

         DEALLOCATE (fjac_save)
         IF (info.ne.0 .or. iflag.ne.0) EXIT outerloop
         
         first_jacobian = .true.                                         !PPPL

      END DO outerloop
c
c     termination, either normal or user imposed.
c
 300  CONTINUE

      IF (info.eq.0 .and. iflag.ne.0) info = 10

 400  CONTINUE

      IF (myid .eq. master) THEN                                         ! MPI
         IF (info.ge.LBOUND(info_array,1) .and. 
     1       info.le.UBOUND(info_array,1)) THEN
            WRITE (6, '(2(/,1x,a))')
     1  'Levenberg-Marquardt optimizer status: ',TRIM(info_array(info))
         ELSE
            WRITE (6, '(a,i5)')' LM info out of bounds: ', info
         END IF
      ENDIF                                                              ! MPI

      IF (ALLOCATED(x_min)) DEALLOCATE(x_min)
      IF (ALLOCATED(fvec_min)) DEALLOCATE(fvec_min)
      IF (ALLOCATED(flip)) DEALLOCATE(flip)
      IF (ALLOCATED(jac_order)) DEALLOCATE(jac_order)
      IF (ALLOCATED(fnorm_array)) DEALLOCATE(fnorm_array)
      IF (ALLOCATED(h_order)) DEALLOCATE(h_order)

      IF (iflag .lt. 0) info = iflag

      IF (nfev.le.1 .and. info.ne.11) THEN
         nfev = 1
         iflag = flag_cleanup                        !!Clean-up last time through for master
         CALL fcn (m, n, x, fvec, iflag, nfev)
      END IF

!DEC$ IF DEFINED (MPI_OPT)
      RETURN

 3000 CONTINUE
      WRITE (6, *) 'MPI_BCAST error in LMDIF, ierr = ', ierr_mpi
!DEC$ ELSE
      WRITE(*, '(2(/,a, f10.2))')
     1     ' Total wall clock time in jacobian multi-process call  = ',
     2     wall_time,
     3     ' Total wall clock time in lev param multi-process call = ',
     4     wall_time_lev
!DEC$ ENDIF
      END SUBROUTINE lmdif
