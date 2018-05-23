
SUBROUTINE lmanalytic(fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, &
        diag, mode, factor, nprint, info, nfev, fjac, &
        ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4, xvmin, xvmax)

      USE mpi_params
      USE safe_open_mod

			USE stel_kinds
			USE vparams, ONLY: sfincs_nmax, sfincs_mmax, ndatafmax
      USE lmpar_mod, ONLY: fjac_mod=>fjac, ldfjac_mod=>ldfjac, &
        ipvt_mod=>ipvt, qtf_mod=>qtf, diag_mod=>diag, delta, fnorm1, &
				lfirst_lm, par, pnorm, spread_ratio, levmarq_param_mp
!DEC$ IF DEFINED (MPI_OPT)
      USE fdjac_mod, ONLY: flip,flag_singletask,flag_cleanup, &
                          fdjac2_mp_queue, jac_order, jac_count, &
                           ix_min, h_order, flag_cleanup_lev, &
                         jac_err, jac_index, n_red !PPPL
!DEC$ ELSE
      USE fdjac_mod, ONLY: max_processors, flip, flag_singletask, &
                         flag_cleanup, fdjac2, jac_order, h_order, &
                          jac_err, jac_index, flag_cleanup_lev, &
												 ix_min, jac_count, n_red
!DEC$ ENDIF

			USE jacfcn_mod, ONLY: jac_analytic, fjac_curr, jacfcn

      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      INCLUDE 'mpif.h'                                       !mpi stuff
!DEC$ ENDIF
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
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
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1, &
        p1=0.1_dp, p5=0.5_dp, p25=0.25_dp, p75=0.75_dp, p0001=1.e-4_dp
      CHARACTER(LEN=130), DIMENSION(0:12) :: info_array
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: iflag, iter, j, l, istat, ikey, lev_step_range,iunit         !PPPL
!      INTEGER :: iflag, iter, j, l, istat, ikey
      INTEGER :: cycle_count, subcycle,i,num_jacs
      REAL(rprec) :: actred, dirder, epsmch, fnorm, &
        gnorm, prered, ratio, sum0, temp, &
        temp1, temp2, xnorm, delta_old, actred_lev, par_old,fnorm_old  !PPPL
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: fjac_save
!DEC$ IF .NOT.DEFINED (MPI_OPT)
      REAL(rprec) :: wall_time, wall_time_lev
!DEC$ ENDIF
      REAL(rprec) :: fnorm_min, epsfcn0, epsfcn_temp
      REAL(rprec) :: x_min(n), fvec_min(m)
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: fnorm_array
      logical :: step_improved, first_jacobian, lsmall_step, lredo_diag !PPPL
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      EXTERNAL fcn
      REAL(rprec), EXTERNAL :: dpmpar, enorm
!-----------------------------------------------

!     SUBROUTINE lmanalytic

!     the purpose of lmanalytic is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of
!     the levenberg-marquardt algorithm. the user must provide a
!     subroutine which calculates the functions. the jacobian is
!     computed using analytic derivatives
!
!     the subroutine statement is
!
!       subroutine lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,
!                        diag,mode,factor,nprint,info,nfev,fjac,
!                        ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling. see lmdif1 for
!				 documentation and should be written as follows.
!
!         subroutine fcn(m, n, x, fvec, iflag, ncnt)
!         integer m,n,iflag
!         real(rprec) x(n),fvec(m)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ----------
!         return
!         end
!
!
!       ftol is a nonnegative input variable. termination
!         occurs when both the actual and predicted relative
!         reductions in the sum of squares are at most ftol.
!         therefore, ftol measures the relative error desired
!         in the sum of squares.
!
!       xtol is a nonnegative input variable. termination
!         occurs when the relative error between two consecutive
!         iterates is at most xtol. therefore, xtol measures the
!         relative error desired in the approximate solution.
!
!       gtol is a nonnegative input variable. termination
!        occurs when the cosine of the angle between fvec and
!         any column of the jacobian is at most gtol in absolute
!         value. therefore, gtol measures the orthogonality
!         desired between the function vector and the columns
!         of the jacobian.
!
!       maxfev is a positive integer input variable. termination
!         occurs when the number of calls to fcn is at least
!         maxfev by the end of an iteration.
!
!				jacfcn is the name of a user-supplied subroutine which
!				calculates the jacobian matrix.
!
!       diag is an array of length n. if mode = 1 (see
!         below), diag is internally set. if mode = 2, diag
!         must contain positive entries that serve as
!         multiplicative scale factors for the variables.
!
!       mode is an integer input variable. if mode = 1, the
!         variables will be scaled internally. if mode = 2,
!         the scaling is specified by the input diag. other
!         values of mode are equivalent to mode = 1.
!
!       factor is a positive input variable used in determining the
!         initial step bound. this bound is set to the product of
!         factor and the euclidean norm of diag*x if nonzero, or else
!         to factor itself. in most cases factor should lie in the
!         interval (.1,100.). 100. is a generally recommended value.
!
!       nprint is an integer input variable that enables controlled
!         printing of iterates if it is positive. in this case,
!         fcn is called with iflag = 0 at the beginning of the first
!         iteration and every nprint iterations thereafter and
!         immediately prior to return, with x and fvec available
!         for printing. if nprint is not positive, no special calls
!         of fcn with iflag = 0 are made.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set
!       nfev is an integer output variable set to the number of
!         calls to fcn.
!
!       fjac is an output m by n array. the upper n by n submatrix
!         of fjac contains an upper triangular matrix r with
!         diagonal elements of nonincreasing magnitude such that
!
!                t     t           t
!               p *(jac *jac)*p = r *r,
!
!         where p is a permutation matrix and jac is the final
!         calculated jacobian. column j of p is column ipvt(j)
!         (see below) of the identity matrix. the lower trapezoidal
!         part of fjac contains information generated during
!         the computation of r.
!
!       ldfjac is a positive integer input variable not less than m
!         which specifies the leading dimension of the array fjac.
!
!       ipvt is an integer output array of length n. ipvt
!         defines a permutation matrix p such that jac*p = q*r,
!         where jac is the final calculated jacobian, q is
!         orthogonal (not stored), and r is upper triangular
!         with diagonal elements of nonincreasing magnitude.
!         column j of p is column ipvt(j) of the identity matrix.
!
!       qtf is an output array of length n which contains
!         the first n elements of the vector (q transpose)*fvec.
!
!       wa1, wa2, and wa3 are work arrays of length n.
!
!       wa4 is a work array of length m.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... dpmpar,enorm,fdjac2,lmpar,qrfac
!
!       fortran-supplied ... ABS,max,min,sqrt,mod
!
!     argonne national laboratory. MINpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********

!DEC$ IF DEFINED (MPI_OPT)
!c     Get mpi parameters
      CALL MPI_COMM_RANK (MPI_COMM_STEL, myid, ierr_mpi)       !mpi stuff
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
      IF (ierr_mpi .ne. 0) STOP 'MPI_COMM_RANK error in LMDIF'
      CALL MPI_COMM_SIZE (MPI_COMM_STEL, numprocs, ierr_mpi)   !mpi stuff
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
      IF (ierr_mpi .ne. 0) STOP 'MPI_COMM_SIZE error in LMDIF'
!DEC$ ENDIF

      !info = 0;      iflag = 0;      nfev = 0;      cycle_count = 0
      info = 0;      iflag = 0;      cycle_count = 0

      info_array(9) = &
     "improper input parameters (number constraints MUST  &
     &be greater than number variables (>0)."

      info_array(1) = &
     "algorithm estimates that the relative error in the &
     &sum of squares is at most ftol."

      info_array(2) = &
     "algorithm estimates that the relative error between &
     & x and the solution is at most xtol."

      info_array(3) = &
     "algorithm estimates that the relative error in the &
     &sum of squares and between x and the solution is at &
     &most ftol and xtol."

      info_array(4) = &
     "the cosine of the angle between fvec and any column &
     &of the jacobian is at most gtol in absolute value."

      info_array(5) = &
     "number of calls to fcn has reached or exceeded maxfev."

      info_array(6) = &
     "ftol is too small. no further reduction in the sum &
     &of squares is possible."

      info_array(7) = &
     "xtol is too small. no further improvement in the &
     &approximate solution x is possible."

      info_array(8) = &
     "gtol is too small. fvec is orthogonal to the columns &
     &of the jacobian to machine precision."

      info_array(0) = &
     "levenberg-marquardt optimizer terminated properly."

      info_array(10) = &
     "Levenberg-Marquardt terminated prematurely due to &
     &function evaluation error."

      info_array(11) = &
     "Levenberg-Marquardt terminated prematurely: &
     &must request more than one (1) processor"

      info_array(12) = &
     "Levenberg-Marquardt terminated prematurely: &
     &because of bad Jacobian"


!DEC$ IF DEFINED (MPI_OPT)
      IF (numprocs > n) THEN
         IF (myid .eq. master) THEN
            WRITE (6, *)'Warning: more processors have been requested&
           & than the maximum (nvar) required = ',n
         END IF
      ELSE IF (numprocs <= 1) THEN
         info = 11
         GOTO 400
      END IF
!DEC$ ENDIF

!     check the input parameters for errors.

      ALLOCATE (flip(n),jac_order(n), fnorm_array(n),h_order(n),jac_err(n),jac_index(n),fjac_curr(m,n), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error in lmdif'

			jac_analytic = .TRUE.

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

!     ASSIGN VALUES TO MODULE VARIABLES (FACILITATES PASSING TO SUBROUTINES)

      ldfjac_mod = ldfjac
      ipvt_mod => ipvt
      fjac_mod => fjac
      diag_mod => diag
      qtf_mod => qtf

!     check the input parameters for errors.

      IF (n.le.0 .or. m.lt.n .or. ldfjac.lt.m .or. ftol.lt.zero &
        .or. xtol.lt. zero .or. gtol.lt.zero .or. maxfev.le.0 &
        .or. factor.le.zero) THEN
         info = 9
         GOTO 400
      END IF

      IF (mode .eq. 2) THEN
         DO j = 1, n
            IF (diag(j) .le. zero) GOTO 300
         END DO
      END IF

!     Set up workers communicator (only master processor here) for initial run
!DEC$ IF DEFINED (MPI_OPT)
      IF (myid .ne. master) THEN
         ikey = MPI_UNDEFINED
      ELSE
         ikey = WORKER_SPLIT_KEY+1
      END IF
      CALL MPI_COMM_SPLIT(MPI_COMM_STEL, ikey, worker_id, MPI_COMM_WORKERS, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF

!     evaluate the function at the starting point and calculate its norm.
      IF (myid .eq. master) THEN
         iflag = flag_singletask
         if (nfev .ne. 0) iflag = 0
         !CALL fcn (m, n, x, fvec, iflag, nfev, fjac)
						CALL fcn (m, n, x, fvec, iflag, nfev)
				 IF (iflag .ne. 0) THEN
            WRITE(6,*) "FIRST RUN FAILS!  IMPROVE INPUT!"
            STOP "ERRROR!"
         END IF
         fnorm = enorm(m,fvec)
         iunit = 12; istat = 0
         CALL safe_open(iunit,istat,'xvec.dat','unknown','formatted', ACCESS_IN='APPEND')
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
      CALL MPI_BCAST(x,n, MPI_REAL8, master, &
                    MPI_COMM_STEL, ierr_mpi)
			!CALL MPI_BARRIER(MPI_COMM_STEL,ierr_mpi)
			IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
			CALL MPI_BCAST(fjac,m*n,MPI_DOUBLE_PRECISION,master,MPI_COMM_STEL,ierr_mpi)
			IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF
      iflag = FLAG_CLEANUP
      IF (myid == master) iflag = flag_cleanup_lev
      !call fcn (m, n, x, fvec, iflag, nfev, fjac)
			call fcn (m, n, x, fvec, iflag, nfev)

!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BCAST(iflag,1, MPI_INTEGER, master, &
                    MPI_COMM_STEL, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
      IF (iflag .ge. 0) CALL &
          MPI_BCAST(fvec, m, MPI_REAL8, master, &
                    MPI_COMM_STEL, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
			! Jacobian computed in call to fcn
			IF (iflag .ge. 0) CALL &
          MPI_BCAST(fjac, m*n, MPI_REAL8, master, &
                    MPI_COMM_STEL, ierr_mpi)
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
 1000 FORMAT (/,' Beginning Levenberg-Marquardt Iterations',/, &
             ' Number of Processors: ',i4,//, &
!     1        ' Number of Processors: ',i4,' (1 controller proc)',//,
             70('='),/,2x,'Iteration',3x,'Processor',7x,'Chi-Sq',7x, &
            'LM Parameter',6x,'Delta Tol'/,70('='))
!DEC$ ELSE
      WRITE (6, 1000) max_processors
 1000 FORMAT (/,' Beginning Levenberg-Marquardt Iterations',/, &
             ' Number processors requested: ', i4,//, &
             59('='),/,2x,'Iteration',8x,'Chi-Sq',7x, &
            'LM Parameter',6x,'Delta Tol',/,59('='))
!DEC$ ENDIF
      IF (myid == master)  WRITE(6, '(2x,i6,8x,i3,7x,1es12.4)') &
                                0, myid, fnorm*fnorm
      CALL FLUSH(6)

!     beginning of the outer loop.

      first_jacobian = .true.                                            !PPPL
      outerloop: DO WHILE (nfev .lt. maxfev)
         delta_old = delta                                               !PPPL
         par_old = par                                                   !PPPL
         fnorm_old = fnorm
!DEC$ IF DEFINED (MPI_OPT)
         CALL MPI_BARRIER(MPI_COMM_STEL,ierr_mpi)                       !PPPL -SAL
         IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF

!!!!        calculate the jacobian matrix.
!				IF (myid .eq. master) THEN
				CALL jacfcn(m,n,x,fvec,fjac,x_min,fnorm,nfev,fvec_min,fnorm_min,epsfcn)
!				END IF

!!DEC$ IF DEFINED (MPI_OPT)
!				CALL MPI_BCAST(fjac,m*n,MPI_DOUBLE_PRECISION, &
!											master,MPI_COMM_STEL,ierr_mpi)
!				CALL MPI_BCAST(x_min,n,MPI_DOUBLE_PRECISION, &
!						master,MPI_COMM_STEL,ierr_mpi)
!				CALL MPI_BCAST(fvec_min,m,MPI_DOUBLE_PRECISION, &
!						master,MPI_COMM_STEL,ierr_mpi)
!				CALL MPI_BCAST(jac_index,n,MPI_DOUBLE_PRECISION, &
!						master,MPI_COMM_STEL,ierr_mpi)
!				CALL MPI_BCAST(jac_err,n,MPI_INTEGER, &
!						master,MPI_COMM_STEL,ierr_mpi)
!				CALL MPI_BCAST(n_red,1,MPI_INTEGER, &
!						master,MPI_COMM_STEL,ierr_mpi)
!				CALL MPI_BCAST(jac_order,n,MPI_INTEGER, &
!						master,MPI_COMM_STEL,ierr_mpi)
!				CALL MPI_BCAST(ix_min,1,MPI_INTEGER, &
!						master,MPI_COMM_STEL,ierr_mpi)
!				CALL MPI_BCAST(fnorm,1,MPI_DOUBLE_PRECISION, &
!						master,MPI_COMM_STEL,ierr_mpi)
!				CALL MPI_BCAST(h_order,n,MPI_DOUBLE_PRECISION, &
!						master,MPI_COMM_STEL,ierr_mpi)
!				CALL MPI_BCAST(fnorm_min,1,MPI_DOUBLE_PRECISION, &
!						master,MPI_COMM_STEL,ierr_mpi)
!!DEC$ ENDIF

				 IF (iflag .lt. 0) EXIT outerloop

				 ipvt = 0
				 wa1  = 0
				 wa2  = 0
				 wa3  = 0

				 CALL qrfac(m, n_red, fjac, ldfjac, .true., ipvt, &
									 n_red, wa1, wa2, wa3)

!        on the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.
!        also on the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.

				 IF (iter .eq. 1 .or. lredo_diag) THEN
            IF (mode .ne. 2) THEN
               diag = 1.0
               DO i = 1, n_red
                  j = jac_index(i)
                  diag(j) = wa2(i)
                  IF (wa2(i) == zero) diag(j) = one
               END DO
            END IF
            wa3 = diag*x
            xnorm = enorm(n,wa3)
            delta = factor*xnorm
            IF (delta .eq. zero) delta = factor
            lredo_diag = .false.
         END IF

!        form (q TRANSPOSE)*fvec and store the first n components in qtf.
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

!        compute the norm of the scaled gradient.
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

!        test for convergence of the gradient norm.
         IF (gnorm .le. gtol) info = 4
         IF (info .ne. 0) EXIT outerloop

!        rescale if necessary.
         IF (mode .ne. 2) THEN
            DO i = 1, n_red
               j = jac_index(i)
               diag(j) = MAX(diag(j),wa2(i))
            END DO
         END IF
         !IF (mode .ne. 2) diag = MAX(diag,wa2)

!        If maxfev == 1 (for analytic gradients) then we calculate the Jacobian and stop the code
!        10/15/12 - SAL (PPPL)

         IF (maxfev.eq.1) info = 5
         IF (info .ne. 0) GOTO 300

!        set up for inner loop (levmarqloop) to determine x update.
         subcycle = 0
         ratio = 0

         ALLOCATE (fjac_save(n,n), stat=istat)
         IF (istat .ne. 0) STOP 'Fjac_save allocation error'
         fjac_save(:n,:n) = fjac(:n,:n)

!        Levenburg Inner Loop
         levmarqloop: DO WHILE (ratio.lt.p0001 .and. subcycle.lt.2)

           subcycle = subcycle + 1
           fjac(:n,:n) = fjac_save(:n,:n)
           spread_ratio = abs(epsfcn/delta/10)                           !PPPL


!        Determine the levenberg-marquardt PARAMETER.
!DEC$ IF DEFINED (MPI_OPT)
           IF (PRESENT(xvmin) .and. PRESENT(xvmax)) THEN
              ! Begin MJL
              print *,"Calling levmarq_param_mp WITH xvmin,xvmax."
              print *,"Here comes xvmin:"
              print *,xvmin
              ! End MJL
							CALL levmarq_param_mp (x, wa1, wa2, wa3, wa4, &
                                   nfev, m, n, iflag, fcn, &
																	 lev_step_range,fnorm_min, &
																	 xvmin,xvmax)
           ELSE
              print *,"Calling levmarq_param_mp WITHOUT xvmin,xvmax."
              CALL levmarq_param_mp (x, wa1, wa2, wa3, wa4, &
                                    nfev, m, n, iflag, fcn, &
                                    lev_step_range,fnorm_min)   !PPPL
           END IF
!DEC$ ELSE
					 CALL levmarq_param(x, wa1, wa2, wa3, wa4, &
								wall_time_lev, nfev, m, n, iflag, fcn)
!DEC$ ENDIF
           IF (iflag .lt. 0) EXIT


!        compute the scaled actual reduction.

           actred = -one
           IF (p1*fnorm1 .lt. fnorm) actred = one - (fnorm1/fnorm)**2

!        compute the scaled predicted reduction (prered) and
!        the scaled directional derivative (dirder).

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

!        compute the ratio of the actual to the predicted reduction.

           actred_lev = actred                                           !PPPL
           ratio = zero
           IF (prered .ne. zero) ratio = actred/prered

!        test if Jacobian calculation gave best improvement
!        (added by MZ and SH) (always false for non-mpi runs)


!DEC$ IF DEFINED (MPI_OPT)
           step_improved = fnorm_min < MIN(fnorm,fnorm1)
!DEC$ ENDIF

           ! Jacobian gives best improvment so take orthagonal step
           IF (step_improved) THEN
              IF (myid == master) &
                jac_count = MIN(count(jac_order(:) .ne. 0), &
                                INT(SQRT(REAL(numprocs))))
!DEC$ IF DEFINED (MPI_OPT)
              call MPI_BCAST(jac_count, 1, MPI_INTEGER, master, &
                            MPI_COMM_STEL, ierr_mpi)
              IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
              IF (ierr_mpi .ne. 0) GOTO 3000
!DEC$ ENDIF

              IF(myid .ne. master) jac_order = 0

!DEC$ IF DEFINED (MPI_OPT)
              CALL MPI_BCAST(jac_order, jac_count, MPI_INTEGER, master, &
                            MPI_COMM_STEL, ierr_mpi)
              IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
              IF (ierr_mpi .ne. 0) GOTO 3000

              CALL MPI_BARRIER(MPI_COMM_STEL,ierr_mpi)                  !PPPL -SAL
              IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF

              CALL stepopt_mp(fcn, wa2, wa4, m, n, x_min, fvec_min, &
															 fnorm_min, iflag, nfev, epsfcn)

              IF (myid .eq. master) THEN
                 WRITE(6,'(a,i6,a)') &
                ' Using minimum from Jacobian/step improvement (', &
                ix_min, ')'
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
           ELSE
              IF (ratio .le. p25) THEN
                 IF (actred .ge. zero) THEN
                    temp = p5
                 ELSE
                    temp = p5*dirder/(dirder + p5*actred)
                 END IF
                 IF (p1*fnorm1.ge.fnorm .or. temp.lt.p1) temp = p1
                 delta = temp*MIN(delta,pnorm/p1)
                 par = par/temp
              ELSE IF (par.eq.zero .or. ratio.ge.p75) THEN
                 delta = pnorm/p5
                 par = p5*par
              END IF
           END IF

!        test for successful iteration.
!        update x, fvec, and their norms.

           if( n_red > numprocs .and. subcycle < 2 .and. &
              ratio .lt. p0001 ) cycle levmarqloop                     !PPPL

           IF (ratio .ge. p0001) THEN
              x = wa2
              ! New way to calculated xnorm
              wa2 = 0.0
              DO i = 1, n_red
                 j = jac_index(i)
                 wa2(j) = diag(j)*x(j)
              END DO
              fvec = wa4
              xnorm = enorm(n,wa2)
              fnorm = fnorm1
              iter = iter + 1
              IF (myid .eq. master) WRITE(6,'(/,3(2x,a,1es10.3)/)') &
               'new minimum = ', fnorm**2,'lm-par = ', par, &
               'delta-tol = ', delta
           END IF

           cycle_count = cycle_count + 1

!        tests for convergence.

           IF (ABS(actred).le.ftol .and. prered.le.ftol &
             .and. p5*ratio.le.one) info = 1

!        next test made more stringent to avoid premature declaration of
!        completion observed on some problems.

           IF (delta.le.xtol*xnorm/10 .and. cycle_count>2 &
              .and. .not. step_improved ) info = 2                     !PPPL
           IF (ABS(actred).le.ftol .and. prered.le.ftol &
             .and. p5*ratio.le.one .and. info.eq.2) info = 3
           IF (info .ne. 0) EXIT levmarqloop

!
!        tests for termination and stringent tolerances.
!
           IF (nfev .ge. maxfev) info = 5
           IF (ABS(actred).le.epsmch .and. prered.le.epsmch &
             .and. p5*ratio.le.one) info = 6
           IF (delta .le. epsmch*xnorm) info = 7
           IF (gnorm .le. epsmch) info = 8
           IF (info .ne. 0) EXIT levmarqloop
!
!        END of the inner loop. repeat if iteration unsuccessful (ratio < p0001)
!
         END DO levmarqloop

         DEALLOCATE (fjac_save)
         IF (info.ne.0 .or. iflag.ne.0) EXIT outerloop

         first_jacobian = .true.  !PPPL

      END DO outerloop

!     termination, either normal or user imposed.

 300  CONTINUE

      IF (info.eq.0 .and. iflag.ne.0) info = 10

 400  CONTINUE

      IF (myid .eq. master) THEN                                         ! MPI
         IF (info.ge.LBOUND(info_array,1) .and. &
            info.le.UBOUND(info_array,1)) THEN
            WRITE (6, '(2(/,1x,a))') &
       'Levenberg-Marquardt optimizer status: ',TRIM(info_array(info))
         ELSE
            WRITE (6, '(a,i5)')' LM info out of bounds: ', info
         END IF
      ENDIF                                                              ! MPI

!      IF (ALLOCATED(x_min)) DEALLOCATE(x_min)
!      IF (ALLOCATED(fvec_min)) DEALLOCATE(fvec_min)
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
      WRITE(*, '(2(/,a, f10.2))') &
          ' Total wall clock time in jacobian multi-process call  = ', &
          wall_time, &
          ' Total wall clock time in lev param multi-process call = ', &
          wall_time_lev
!DEC$ ENDIF
END SUBROUTINE lmanalytic
