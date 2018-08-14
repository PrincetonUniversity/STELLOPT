SUBROUTINE BFGS_backtrack(m,n,nfev,fcn,p_curr,f_curr,x,grad_curr, &
             f_new,x_new,grad_new,alpha_backtrack,c_armijo, &
             rho_backtrack,beta_hessian,alpha_min, dx_init, &
             x_min, fvec_min, fnorm_array)

  USE BFGS_params
  IMPLICIT NONE
!-----------------------------------------------------------------------
!     Subroutine:    BFGS_backtrack
!     Authors:       E. J. Paul (ejpaul@umd.edu)
!                    J.C. Schmitt (jcschmitt@auburn.edu)
!     Date:          05/26/2018, 08/02/2018
!     Description:   This subroutine performs a backtracking line search
!                    given a descent direction p_curr and a function fcn.
!                    This is to be called by
!                    BFGS_FD, returning a new point satisfying the
!                    Armijo condition. This is based off of Nocedal
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
!     alpha_min       Minimum allowed step size for line search.
!     c_armijo        Parameter used to test significant decrease condition.
!     rho_backtrack   Factor by which step size is shrunk.
!     beta_hessian       Scaling of steepest descent direction. This is only
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
  REAL(rprec), INTENT(IN) :: x(n),grad_curr(n), f_curr, dx_init
  REAL(rprec), INTENT(INOUT) :: p_curr(n)
  REAL(rprec), INTENT(OUT) :: x_new(n), grad_new(n), f_new
  REAL(rprec), INTENT(IN) :: alpha_backtrack, c_armijo, rho_backtrack, &
                             beta_hessian
  REAL(rprec), INTENT(IN) :: alpha_min
  REAL(rprec), INTENT(IN) :: x_min(n), fvec_min(m), fnorm_array(n)

  !-----------------------------------------------
  !   E x t e r n a l   F u n c t i o n s
  !-----------------------------------------------
  EXTERNAL fcn
  REAL(rprec), EXTERNAL :: enorm

  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  REAL(rprec) :: fvec_new(m), grad_dot_p, fnorm_new, alpha
  INTEGER :: iflag, niter, nvar, istat, iunit
  REAL(rprec) :: fjac_curr(m,n)  ! the current jacobian
  REAL(rprec) :: fnorm_min !, epsfcn0, epsfcn_temp
 
 
  ! Set initial step length
  alpha = alpha_backtrack

  grad_dot_p = dot_product(grad_curr,p_curr)

  if (DEBUG_BFGS .and. (myid .eq. master)) then
    print *,"<---Backtrack Step 0a: grad_dot_p: ", grad_dot_p
  end if

  niter = 0

  if (grad_dot_p >= 0) then
    if (myid .eq. master) then
      write(*,"(A)") "Warning. Backtrack was not called with a descent direction. A steepest descent step will be taken (scaled by beta_hessian)"
    end if
    p_curr = -grad_curr*beta_hessian
  end if

  x_new = x + alpha*p_curr
  iflag = myid
  call fcn(m, n, x_new, fvec_new, iflag, nfev)
  nfev = nfev+1

  ! Save to xvec.dat
  fnorm_new = enorm(m,fvec_new)
!  if (myid .eq. master) then
!    iunit = 12; istat = 0
!    CALL safe_open(iunit,istat,'xvec.dat','unknown','formatted', ACCESS_IN='APPEND')
!    WRITE(iunit,'(2(2X,I5.5))') n,nfev
!    WRITE(iunit,'(10ES22.12E3)') x_new(1:n)
!    WRITE(iunit,'(ES22.12E3)') fnorm_new
!    CLOSE(iunit)
!  end if

  if (myid .eq. master) then
    write (6, 56) 'Entering Armijo Backtracking Search', &
                  'Previous chi^2 = ', f_curr, 'Target chi^2 <=  ', &
                   f_curr+c_armijo*alpha*grad_dot_p, &
                  'nfev', 'Iteration', 'Alpha', 'Chi-Sq'
56  FORMAT (/, A, /, A, 1es12.4, /, A, 1ES12.4, /, 50('='), /, &
            4x, A, 3x, A, 6X, A, 8X, A, /, 50('='))
  end if

  if (myid .eq. master) then
    write(6, '(2X,I6,2X,I6,4X,1ES12.4,3X,1ES12.4)') &
        nfev, niter, alpha, fnorm_new*fnorm_new
  end if

  ! Cleanup after fcn call
  !IF (myid .eq. master)  iflag = FLAG_CLEANUP_BFGS
  !nfev = nfev+1
  !CALL fcn (m,n,x_new,fvec_new,iflag,nfev)
  ! nfev = nfev+1

  f_new = (enorm(m,fvec_new))**2
  
  if (DEBUG_BFGS .and. (myid .eq. master)) then
    print *, '<---------------------------->'
    print *, "<---Backtrack Step 0b"
    print *, "<-x", x
    print *, "<-alpha", alpha
    print *, "<-p_curr", p_curr
    print *, "<-x_new", x_new
    print *, "<-fvec_new: ",fvec_new
    print *, "<-f_curr: ",f_curr
    print *, "<-c_armijo: ",c_armijo
    print *, "<-grad_dot_p: ",grad_dot_p
    print *, "<-f_new: ",f_new
    print *, '<---------------------------->'
  end if

  ! Test Armijo condition
  if (DEBUG_BFGS .and. (myid .eq. master)) then
    print *, '<---------------------------->>>>'
    print *, "<---Backtrack Step 0c"
    print *, "<--f_curr: ", f_curr
    print *, "<--c_armijo: ", c_armijo
    print *, "<--alpha: ", alpha
    print *, "<--grad_dot_p: ", grad_dot_p
    print *, "<--f_curr + c_armijo*alpha*grad_dot_p: ", f_curr + c_armijo*alpha*grad_dot_p
    !print *, "<--fvec_new: ",fvec_new
    print *, "<--f_new: ", f_new
    print *, "<--Desire f_new<f_curr+c_armijo*alpha*grad_dot_p"
    print *, '<---------------------------->>>>'
    print *, '<========Entering alpha search'
  end if
  ! Initialize the exit status to 'normal'
  bt_exit_flag = BT_EXIT_NORMAL

  do while (f_new > (f_curr + c_armijo*alpha*grad_dot_p))
    niter = niter + 1
    alpha = rho_backtrack*alpha
    if (alpha < alpha_min) then
      bt_exit_flag = BT_EXIT_STEPTOOSMALL
      exit
    end if

    x_new = x + alpha*p_curr
    call fcn(m,n,x_new,fvec_new,iflag,nfev)
    nfev = nfev + 1

    if (myid .eq. master) then
      write(6, '(2X,I6,2X,I6,4X,1ES12.4,3X,1ES12.4)') &
            nfev, niter, alpha, fnorm_new*fnorm_new
    end if

  if (DEBUG_BFGS .and. (myid .eq. master)) then
      print *, '<---------------------------->>>>>>>>'
      print *, "<---Backtrack Step 1"
      print *, '<----Tried alpha=', alpha
      print *, '<----with c_armijo=', c_armijo
      print *, '<----and grad_dot_p=', grad_dot_p
      print *, '<----x_new=', x_new
      print *, "<----fvec_new: ",fvec_new
      print *, "<----f_new: ", f_new
      print *, "<----f_curr: ", f_curr
      print *, "<----f_curr + c_armijo*alpha*grad_dot_p: ", f_curr + c_armijo*alpha*grad_dot_p
      print *, '<---------------------------->>>>>>>>'
    end if
    bt_exit_flag = BT_EXIT_NORMAL


  if (DEBUG_BFGS .and. (myid .eq. master)) then
      print *, '<---------------------------->>>>>>>>>>>>>>>'
      print *, "<---Backtrack 2"
      print *, "<---niter: ", niter
      print *, "<---alpha: ", alpha
      print *, "<---f_curr + c_armijo*alpha*grad_dot_p: ", f_curr + c_armijo*alpha*grad_dot_p
      print *, '<---x_new=', x_new
      print *, "<---fvec_new: ",fvec_new
      print *, "<---f_new: ", f_new
      print *, '<---------------------------->>>>>>>>>>>>>>>'
    end if
  end do

  select case (bt_exit_flag)
    case (BT_EXIT_NORMAL)
      if (myid .eq. master) then
        write(*,"(A)") "<----Linesearch terminated: Armijo condition satisfied."
        ! Save to xvec.dat
        fnorm_new = enorm(m,fvec_new)
        nfev = nfev + niter
        if (myid .eq. master) then
          iunit = 12; istat = 0
          CALL safe_open(iunit,istat,'xvec.dat','unknown','formatted', ACCESS_IN='APPEND')
          WRITE(iunit,'(2(2X,I5.5))') n,nfev
          WRITE(iunit,'(10ES22.12E3)') x_new(1:n)
          WRITE(iunit,'(ES22.12E3)') fnorm_new
          CLOSE(iunit)
        end if


        ! Cleanup after fcn call
        iflag =  FLAG_CLEANUP_BFGS
        !nfev = nfev+1
        CALL fcn (m,n,x_new,fvec_new,iflag,nfev)
        f_new = (enorm(m,fvec_new))**2
      end if

  if (DEBUG_BFGS .and. (myid .eq. master)) then
    print *, '<========Exiting alpha search'
    write(*,"(A,I6,A)") "Backtracking line search completed in ",niter," iterations."
    write(*,"(A,E22.14)") "New function value: ",f_new
    print *, '<----Calculating Finite Differences'
    print *, '<----x_new=', x_new
    print *, '<----dx_init=', dx_init
  end if

    IF (myid .eq. master) then
       WRITE (6, 500) 
500    FORMAT (/,' Performing BFGS-Finite Difference Iterations', /, &
                40('='), / ,2x, &
               'Iteration', 3x, 'Processor', 7x, 'Chi-Sq', 7x, &
               /, 40('='))
        WRITE(6, '(2x,i6,8x,i3,7x,1es12.4)'), 0, myid, f_new
    END IF

  ! Compute new gradient
  CALL fdjac2_mp_queue(fcn, fnorm_new, m, n, x_new, fvec_new, fjac_curr, &
                       m, iflag, nfev, dx_init, fnorm_min, x_min, &
                       fvec_min, fnorm_array, .true.)
  !nfev = nfev + n
  grad_new = 2*matmul(fvec_new, fjac_curr)

  ! Check the flag to see if 'flipping' is permitted
  if (enable_flip .eqv. .false.) flip = .false. 

!  if ((myid .eq. master)) then
!      print *, "<----Finite Differences Calculated"
!  end if
  if (DEBUG_BFGS .and. (myid .eq. master)) then
      print *, '<---------------------------------'
      print *, "<----Finite Differences Calculated"
      print *, '<----x_new=', x_new
      print *, "<----fvec_new: ", fvec_new
      print *, "<----f_new: ", f_new
      print *, "<----fjac_curr: ", fjac_curr
      print *, "<----grad_new: ", grad_new
      write(*,"(A,E22.14)") "New gradient norm: ", enorm(n,grad_new)
  end if
    case (BT_EXIT_STEPTOOSMALL)
      if (myid .eq. master) then
        write(*,"(A)") "<----Linesearch terminated: Stepsize too small."
      end if
  end select

  END SUBROUTINE BFGS_backtrack
