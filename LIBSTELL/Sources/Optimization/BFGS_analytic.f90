SUBROUTINE BFGS_analytic(fcn,m,n,x,ftol,gtol,maxfev,nfev,alpha_backtrack,&
      c_armijo,rho_backtrack,beta_hess,alpha_min)

  USE mpi_params
  USE safe_open_mod
  USE stel_kinds
  USE fdjac_mod, ONLY: FLAG_CLEANUP, FLAG_CLEANUP_LEV, flag_singletask
  ! flag_singletask = -1, flag_cleanup = -100, flag_cleanup_lev = -101
  USE jacfcn_mod, ONLY: fjac_curr, write_jacobian, jac_analytic

!-----------------------------------------------------------------------
!     Subroutine:    BFGS_analytic
!     Authors:       E. J. Paul (ejpaul@umd.edu)
!     Date:          05/26/2018
!     Description:   This subroutine performs a BFGS quasi-newton
!                    optimization of user-supplied fcn beginning at
!                    x. This is based off of Nocedal & Wright Algorithm 6.1.
!
!    Input parameters
!    fcn             User-specified subroutine to compute objective function
!    m               Number of target functions.
!    n               Number of variables.
!    x               Vector of variables.
!    ftol            Tolerance on norm of function value.
!    gtol            Tolerance on norm of gradient.
!    maxfev          Maximum number of function evals.
!    alpha_backtrack Initial step size for line search. Should typically be 1.
!    c_armijo        Parameter used to test significant decrease condition.
!    rho_backtrack   Factor by which step size is shrunk in backtracking line
!                    search.
!    beta_hess       Scaling used for initial Hessian approximation. This
!                    sets the initial stepsize and is used to scale the
!                    stepsize if a reasonable descent direction is not
!                    obtained.
!		 alpha_min	     Minimum allowed step size for line search.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!DEC$ IF DEFINED (MPI_OPT)
INCLUDE 'mpif.h'
!DEC$ ENDIF

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  INTEGER :: m, n, maxfev, nfev
  REAL(rprec), INTENT(in) ::  ftol, gtol
  REAL(rprec), DIMENSION(n) :: x
  REAL(rprec), INTENT(IN) :: alpha_backtrack, c_armijo, &
      rho_backtrack, beta_hess, alpha_min

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  REAL(rprec) :: fnorm, f_new, fvec(m)
  INTEGER :: iflag, istat, iunit, ikey, nvar, iter, nvar1, nvar2, exit_code
  REAL(rprec) :: hess(n,n), f_curr, grad_curr(n), y_curr(m), s_curr(n)
  REAL(rprec) :: eye(n,n), rho_curr, p_curr(n), x_new(n)
  REAL(rprec) :: s_s_outer(n,n), s_y_outer(n,n), mat1(n,n), mat2(n,n)
  REAL(rprec) :: grad_new(n)

!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
  EXTERNAL fcn
  REAL(rprec), EXTERNAL :: dpmpar, enorm

  ALLOCATE(fjac_curr(m,n))
  jac_analytic = .TRUE.

!DEC$ IF DEFINED (MPI_OPT)
   ! Get mpi parameters
  CALL MPI_COMM_RANK (MPI_COMM_STEL, myid, ierr_mpi)
  IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
  IF (ierr_mpi .ne. 0) STOP 'MPI_COMM_RANK error in LMDIF'
  CALL MPI_COMM_SIZE (MPI_COMM_STEL, numprocs, ierr_mpi)
  IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
  IF (ierr_mpi .ne. 0) STOP 'MPI_COMM_SIZE error in LMDIF'
!DEC$ ENDIF

  IF (numprocs>1) THEN
    print *,"WARNING! BFGS_analytic is a serial optimizer. The additional NOPTIMIZERS will be wasted."
  END IF

  ! Check the input parameters for errors.

  IF ((n<0) .or. (m<0) .or. (ftol<0) .or. (gtol<0) &
    .or. (maxfev<0) .or. (alpha_backtrack<0) .or. (c_armijo>1) .or. &
    (c_armijo<0) .or. (rho_backtrack<0) .or. (rho_backtrack>1) &
    .or. (beta_hess<0)) THEN
     STOP "Error! BFGS_analytic called with improper arguments."
  END IF

  ! Evaluate function at intial point
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
     CALL safe_open(iunit,istat,'xvec.dat','unknown','formatted', ACCESS_IN='APPEND')
     WRITE(iunit,'(2(2X,I5.5))') n,0
     WRITE(iunit,'(10ES22.12E3)') x(1:n)
     WRITE(iunit,'(ES22.12E3)') fnorm
     CLOSE(iunit)

    ! Write jacobian
    CALL write_jacobian(nfev,m,n)

    ! Cleanup after fcn call
    iflag = FLAG_CLEANUP_LEV
    CALL fcn (m, n, x, fvec, iflag, nfev)
    nfev = nfev+1

    ! Compute value of objective function
    f_curr = fnorm**2
    write(*,"(A,E22.14)") "Initial function value: ", f_curr

    ! Compute gradient of objective function
    ! f_curr = sum(fvec^2)
    ! fvec = (vals(1:m)-targets(1:m))/abs(sigmas(1:m))
    ! fjac = d(fvec)/d(vars)
    ! d(f_curr)/d(vars) = sum(2*fvec*fjac)
    grad_curr = 2*matmul(fvec,fjac_curr)

    ! Initial Hessian set to multiple of identity
    eye = 0
    DO nvar = 1,n
      eye(nvar,nvar) = 1
    END DO
    hess = beta_hess*eye

    iter = 0
    ! Main BFGS loop - See algorithm 6.1 Nocedal & Wright
		exit_code = -1
    do while ((enorm(n,grad_curr)>gtol) .and. (f_curr > ftol) .and. (nfev < maxfev))
			if (enorm(n,grad_curr)<gtol) then
				exit_code = 1
				exit
			end if
			if (f_curr<ftol) then
				exit_code = 2
				exit
			end if
      ! Compute search direction
      p_curr = -matmul(hess,grad_curr)

      ! Line search to compute x_new, grad_new, f_new
      CALL backtrack(m,n,nfev,fcn,p_curr,f_curr,x, &
        grad_curr,f_new,x_new,grad_new,alpha_backtrack,c_armijo,rho_backtrack,beta_hess,alpha_min)

      ! Form quantities needed for Hessian update
      s_curr = x_new - x
      y_curr = grad_new - grad_curr
      rho_curr = 1/dot_product(y_curr,s_curr)

      ! Form outer products needed for Hessian update
      s_y_outer = 0
      s_s_outer = 0
      do nvar1 = 1,n
        do nvar2 = 1,n
          s_y_outer(nvar1,nvar2) = s_curr(nvar1)*y_curr(nvar2)
          s_s_outer(nvar1,nvar2) = s_curr(nvar1)*s_curr(nvar2)
        end do
      end do

      ! Rescale if this is the initial Hessian (see Nocedal & Wright 6.20)
      if (iter==0) then
        hess = (dot_product(y_curr,s_curr)/dot_product(y_curr,y_curr))*eye
      end if

      ! Update approximate Hessian (See Nocedal & Wright 6.17)
      mat1 = eye - rho_curr*s_y_outer
      mat2 = eye - rho_curr*transpose(s_y_outer)
      hess = matmul(mat1,matmul(hess,mat2)) + rho_curr*s_s_outer

      ! Update x, f_curr, and grad_curr
      x = x_new
      f_curr = f_new
      grad_curr = grad_new
      iter = iter + 1
    end do
		if (exit_code==1) then
			write(*,"(A)") "BFGS_analytic terminated due to gtol."
		elseif (exit_code==2) then
			write(*,"(A)") "BFGS_analaytic terminated due to ftol."
		elseif (exit_code==-1) then
			write(*,"(A)") "BFGS_analytic terminated due to nfev."
		else
			write(*,"(A)") "BFGS_analytic terminated for an unknown reason."
		end if

    write(*,"(A,I2,A)") "BFGS_analytic terminated in ", iter, " iterations."
    write(*,"(A,E22.14)") "Function value: ", f_curr
    write(*,"(A,E22.14)") "New gradient norm: ", enorm(n,grad_curr)
    write(*,"(A,I2,A)") "Function evaluations: ", nfev
  END IF

END SUBROUTINE BFGS_analytic
