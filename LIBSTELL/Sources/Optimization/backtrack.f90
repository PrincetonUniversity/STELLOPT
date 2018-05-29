
  SUBROUTINE backtrack(m,n,nfev,fcn,p_curr, &
            f_curr,x,grad_curr,f_new,x_new,grad_new,alpha_init,c,rho)
  USE stel_kinds
  USE jacfcn_mod, ONLY: fjac_curr
  USE fdjac_mod, ONLY: FLAG_CLEANUP, FLAG_CLEANUP_LEV, flag_singletask
  USE mpi_params
  USE safe_open_mod
  USE jacfcn_mod, ONLY: write_jacobian

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
!-----------------------------------------------------------------------

  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  INTEGER, INTENT(IN) :: m, n
  INTEGER, INTENT(INOUT) :: nfev
  REAL(rprec), INTENT(IN) :: x(n),grad_curr(n), f_curr
  REAL(rprec), INTENT(OUT) :: x_new(n), grad_new(n)
  REAL(rprec), INTENT(OUT) :: f_new, p_curr(n)
  REAL(rprec), INTENT(IN) :: alpha_init, c, rho

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
  alpha = alpha_init

  grad_dot_p = sum(grad_curr*p_curr)

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

  f_new = (sum(fvec_new**2))
  write(*,"(A,E22.14)") "New function value: ",f_new

  niter = 0
  ! Test Armijo condition
  do while (f_new > (f_curr + c*alpha*grad_dot_p))
    alpha = rho*alpha

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

    f_new = (sum(fvec_new**2))
    niter = niter + 1
  end do
  write(*,"(A,I2,A)") "Backtracking line search completed in ",niter," iterations."
  write(*,"(A,E22.14)") "New function value: ",f_new
  write(*,"(A,E22.14)") "New minimum: ", x_new
  ! Compute new gradient
  DO nvar = 1,n
    grad_new(nvar) = 2*sum(fvec_new*fjac_curr(:,nvar))
  END DO
  write(*,"(A,E22.14)") "New gradient norm: ", enorm(n,grad_new)

  END SUBROUTINE backtrack
