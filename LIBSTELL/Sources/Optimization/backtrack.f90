
  SUBROUTINE backtrack(m,n,nfev,fcn,p_curr,f_curr,x,grad_curr, &
             f_new,x_new,grad_new,alpha_backtrack,c_armijo,&
             rho_backtrack,beta_hess)
  USE stel_kinds
  USE fdjac_mod, ONLY: FLAG_CLEANUP, FLAG_CLEANUP_LEV, flag_singletask
  USE mpi_params
  USE safe_open_mod
  USE jacfcn_mod, ONLY: write_jacobian, fjac_curr

!-----------------------------------------------------------------------
!     Subroutine:    backtrack
!     Authors:       E. J. Paul (ejpaul@umd.edu)
!     Date:          05/26/2018
!     Description:   This subroutine performs a backtracking line search
!                    given a descent direction p_curr and a function fcn.
!                    This was originally written to be called by
!                    BFGS_analytic, returning a new point satisfying the
!                    Armijo condition. This was based off of Nocedal
!                    and Wright Algorithm 3.1.
!
!     Input parameters
!     m               Number of target functions.
!     n               Number of variables.
!     fcn             User-specified subroutine to compute objective function.
!     p_curr          Descent direction chosen by BFGS.
!     f_curr          Current function value.
!     x               Vector of variables.
!     grad_curr       Vector of objective function gradients.
!     alpha_backtrack Initial step size for line search. Should typically be 1.
!     c_armijo        Parameter used to test significant decrease condition.
!     rho_backtrack   Factor by which step size is shrunk.
!     beta_hess       Scaling of steepest descent direction. This is only
!                     used if p_curr is not a descent direction.
!
!     Output parameters
!     nfev            Number of function evaluations.
!     f_new           Function value at optimal point in line search.
!     x_new           Vector of variables at optimal point in line search.
!     grad_new        Objective function gradient at optimal point in line search.
!-----------------------------------------------------------------------

  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  INTEGER, INTENT(IN) :: m, n
  INTEGER, INTENT(INOUT) :: nfev
  REAL(rprec), INTENT(IN) :: x(n),grad_curr(n), f_curr
  REAL(rprec), INTENT(INOUT) :: p_curr(n)
  REAL(rprec), INTENT(OUT) :: x_new(n), grad_new(n), f_new
  REAL(rprec), INTENT(IN) :: alpha_backtrack, c_armijo, rho_backtrack, beta_hess

  !-----------------------------------------------
  !   E x t e r n a l   F u n c t i o n s
  !-----------------------------------------------
  EXTERNAL fcn
  REAL(rprec), EXTERNAL :: enorm

  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  REAL(rprec) :: fvec_new(m), grad_dot_p, fnorm_new, alpha
  INTEGER :: iflag, niter, nvar

  ! Set initial step length
  alpha = alpha_backtrack

  grad_dot_p = dot_product(grad_curr,p_curr)
  print *,"grad_dot_p: ", grad_dot_p
  if (grad_dot_p >= 0) then
    write(*,"(A)") "Warning. Backtrack was not called with a descent direction. A steepest descent step will be taken (scaled by beta_hess)"
    p_curr = -grad_curr*beta_hess
  end if

  x_new = x + alpha*p_curr
  iflag = myid
  call fcn(m,n,x_new,fvec_new,iflag,nfev)

  ! Save to xvec.dat
  fnorm_new = enorm(m,fvec_new)
  iunit = 12; istat = 0
  CALL safe_open(iunit,istat,'xvec.dat','unknown','formatted', ACCESS_IN='APPEND')
  WRITE(iunit,'(2(2X,I5.5))') n,nfev
  WRITE(iunit,'(10ES22.12E3)') x_new(1:n)
  WRITE(iunit,'(ES22.12E3)') fnorm_new
  CLOSE(iunit)

  ! Write jacobian
  CALL write_jacobian(nfev,m,n)

  ! Cleanup after fcn call
  iflag = FLAG_CLEANUP_LEV
  CALL fcn (m,n,x_new,fvec_new,iflag,nfev)
  nfev = nfev+1

  f_new = (enorm(m,fvec_new))**2
  write(*,"(A,E22.14)") "New function value: ",f_new

  niter = 0
  ! Test Armijo condition
  print *,"f_curr + c_armijo*alpha*grad_dot_p: ", f_curr + c_armijo*alpha*grad_dot_p
  print *,"f_new: ", f_new
  do while (f_new > (f_curr + c_armijo*alpha*grad_dot_p))
    alpha = rho_backtrack*alpha

    x_new = x + alpha*p_curr
    call fcn(m,n,x_new,fvec_new,iflag,nfev)

    ! Save to xvec.dat
    fnorm_new = enorm(m,fvec_new)
    iunit = 12; istat = 0
    CALL safe_open(iunit,istat,'xvec.dat','unknown','formatted', ACCESS_IN='APPEND')
    WRITE(iunit,'(2(2X,I5.5))') n,nfev
    WRITE(iunit,'(10ES22.12E3)') x_new(1:n)
    WRITE(iunit,'(ES22.12E3)') fnorm_new
    CLOSE(iunit)

    ! Write jacobian
    CALL write_jacobian(nfev,m,n)

    ! Cleanup after fcn call
    iflag =  FLAG_CLEANUP_LEV
    CALL fcn (m,n,x_new,fvec_new,iflag,nfev)
    nfev = nfev+1
    print *,"niter: ", niter
    print *,"alpha: ", alpha
    print *,"f_curr + c_armijo*alpha*grad_dot_p: ", f_curr + c_armijo*alpha*grad_dot_p
    print *,"f_new: ", f_new
    f_new = (enorm(m,fvec_new))**2
    niter = niter + 1
  end do
  write(*,"(A,I2,A)") "Backtracking line search completed in ",niter," iterations."
  write(*,"(A,E22.14)") "New function value: ",f_new
  ! Compute new gradient
  grad_new = 2*matmul(fvec_new,fjac_curr)
  write(*,"(A,E22.14)") "New gradient norm: ", enorm(n,grad_new)

  END SUBROUTINE backtrack
