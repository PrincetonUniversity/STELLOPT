SUBROUTINE BFGS_FD(fcn, m, n, x, ftol, gtol, maxfev, nfev, &
                   alpha_backtrack, c_armijo, rho_backtrack, &
                   beta_hessian, alpha_min, dx_init)

  USE bfgs_params

!-----------------------------------------------------------------------
!     Subroutine:    BFGS_FD
!     Authors:       E. J. Paul (ejpaul@umd.edu)
!                    J.C. Schmitt (jcschmitt@auburn.edu)
!     Date:          05/26/2018, 08/02/2018
!     Description:   This subroutine performs a BFGS quasi-newton
!                    optimization of user-supplied fcn beginning at
!                    x. This is based off of Nocedal & Wright Algorithm 6.1.
!                    Initially written to handle analytic derivatives. (EJP)
!                    Modified to handle finite differences. (JCS)
!
!    Input parameters
!    fcn             User-specified subroutine to compute objective function
!    m               Number of target functions.
!    n               Number of variables.
!    x               Vector of variables.
!    ftol            Tolerance on norm of function value.
!    gtol            Tolerance on norm of gradient.
!    maxfev          Maximum number of function evals.
!    alpha_backtrack Initial step size for line search.
!                    Should typically be 1.
!    c_armijo        Parameter used to test significant decrease condition.
!    rho_backtrack   Factor by which step size is shrunk in
!                    backtracking line search.
!    beta_hessian    Scaling used for initial Hessian approximation. This
!                    sets the initial stepsize and is used to scale the
!                    stepsize if a reasonable descent direction is not
!                    obtained.
!    alpha_min	     Minimum allowed step size for line search.
!    dx_init         The inital step 'dx' for the change in variables 
!                    used for the inital Jacobian calculation
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
      rho_backtrack, beta_hessian, alpha_min, dx_init

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  REAL(rprec) :: fnorm, f_new, fvec(m)
  INTEGER :: iflag, istat, iunit, ikey, nvar, iter, nvar1, nvar2
  INTEGER :: ii, jj
  ! nvar: used to loop over the 'n' variables
  ! iter: counts the number of iterations through the BFGS loop
  REAL(rprec) :: hessian(n,n), f_curr, grad_curr(n), y_curr(m), s_curr(n)
  REAL(rprec) :: invHessian(n,n) ! The inverse Hessian (inverse of 'hessian')
  REAL(rprec) :: eye(n,n), rho_curr, p_curr(n), x_curr(n), x_new(n)
  ! eye: Identigy matrix
  REAL(rprec) :: fjac_curr(m,n)  ! the 'current' jacobian
  REAL(rprec) :: fjac_new(m,n)  ! the 'new' jacobian
  REAL(rprec) :: s_s_outer(n,n), s_y_outer(n,n), mat1(n,n), mat2(n,n)
  REAL(rprec) :: grad_new(n)
  REAL(rprec) :: fnorm_min !, epsfcn0, epsfcn_temp
  REAL(rprec), ALLOCATABLE, DIMENSION(:) :: x_min, fvec_min
  REAL(rprec), ALLOCATABLE, DIMENSION(:) :: fnorm_array
  CHARACTER(16) :: string1
  CHARACTER(256) :: string2

!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
  EXTERNAL fcn
  REAL(rprec), EXTERNAL :: dpmpar, enorm

!DEC$ IF DEFINED (MPI_OPT)
    ! Get mpi parameters - check for mpi errors after each call 
    ! to help identify bugs, errors and other problems

    ! Assign the process rank to variable 'myid'
    CALL MPI_COMM_RANK (MPI_COMM_STEL, myid, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
    IF (ierr_mpi .ne. 0) STOP 'MPI_COMM_RANK error in BFGS_FD'

    ! Assign the size of the group to 'numprocs'
    CALL MPI_COMM_SIZE (MPI_COMM_STEL, numprocs, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
    IF (ierr_mpi .ne. 0) STOP 'MPI_COMM_SIZE error in BFGS_FD'
!DEC$ ELSE
    ! If not using mpi, assign a meaningful value to 'numprocs'
    numprocs = 1
    ! 'myid' is set to the default value of 'master' in mpi_params
!DEC$ ENDIF

!DEC$ IF DEFINED (MPI_OPT)
    IF ((numprocs > n) .and. (myid .eq. master)) THEN
      WRITE (6, '(2A,2X,I5)'), &
             'K====Warning: more processors have been requested ', &
             'than the maximum (nvar) required = ', n
    END IF
!DEC$ ENDIF

  ! Check the input parameters for errors.
  IF ( (n<0) .or. (m<0) .or. (ftol<0) .or. (gtol<0) &
    .or. (maxfev<0) .or. (alpha_backtrack<0) .or. (c_armijo>1) .or. &
    (c_armijo<0) .or. (rho_backtrack<0) .or. (rho_backtrack>1) &
    .or. (beta_hessian<0) ) THEN
     STOP "K====Error! BFGS_FD called with improper input arguments."
  END IF

! Alllocate memory for fdjac2_mp_queue
  ALLOCATE (x_min(n), fvec_min(m), flip(n), jac_order(n), &
            fnorm_array(n), h_order(n), jac_err(n), jac_index(n), &
            stat=istat)
  IF (istat .ne. 0) STOP 'K====ALLOCATION ERROR IN BFGS_FD'

!  Initialize the newly allocated variables (at least the ones that need it)

  flip = .false. 
  jac_order = 0

!     Set up workers communicator (only master performs this task)
!     for the initial run only
!DEC$ IF DEFINED (MPI_OPT)
    IF (myid .ne. master) THEN
       ikey = MPI_UNDEFINED
    ELSE
       ikey = WORKER_INITIAL_RUN_SPLIT_KEY
    END IF
    CALL MPI_COMM_SPLIT(MPI_COMM_STEL, ikey, worker_id, &
                        MPI_COMM_WORKERS, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF

  ! Evaluate function at intial point and calculate its norm.
  IF (myid .eq. master) THEN

    ! Set iflag to flag_singletask to tell fcn what to do
    iflag = FLAG_SINGLETASK

    ! JCS: How would nfev ever NOT be equal to zero?
    IF (nfev .ne. 0) then
      STOP '<----NFEV IS NOT ZERO!!!'
      iflag = 0
    END IF
    CALL fcn (m, n, x, fvec, iflag, nfev)

    IF (iflag .ne. 0) THEN
      WRITE(6,*) "<----FIRST RUN FAILS!  IMPROVE INPUT!"
      STOP "K====ERROR in BFGS_FD!"
    END IF

    ! Calculate the Euclidean norm here
    fnorm = enorm(m,fvec)

    ! Write useful information to 'xvec.dat'
    iunit = 12; istat = 0
    CALL safe_open(iunit,istat,'xvec.dat','unknown','formatted', &
                   ACCESS_IN='APPEND')
    ! Number of variables, 'n', followed by a '0' (iteration is 0)
    WRITE(iunit,'(2(2X,I5.5))') n,0
    ! The variables, 'x'
    WRITE(iunit,'(10ES22.12E3)') x(1:n)
    ! The norm, 'fnorm'
    WRITE(iunit,'(ES22.12E3)') fnorm
    CLOSE(iunit)

!DEC$ IF DEFINED (MPI_OPT)
      ! Mark the communicator for deallocation
      CALL MPI_COMM_FREE(MPI_COMM_WORKERS, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)      
!DEC$ ENDIF
  END IF ! End of initial evaluation

!DEC$ IF DEFINED (MPI_OPT)

    CALL MPI_BCAST(x,n, MPI_REAL8, master, &
                   MPI_COMM_STEL, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF

  ! Cleanup after fcn call
  IF (myid .eq. master) iflag = FLAG_CLEANUP_BFGS

  if (myid .eq. master)  print *, "<----In BFGS_FD. Master will do first bfgs cleanup"
  ! The master will do the bfgs cleanup.
  ! The workers will do 'something else'?.
  !if (myid .eq. master)  write (6,'(A30,i5,a6,i5)'), '<----Before cleanup: myid=', myid, 'iflag=', iflag
  CALL fcn (m, n, x, fvec, iflag, nfev)
  !if (myid .eq. master)  write (6,'(A30,i5,a6,i5)'), '<----After cleanup: myid=', myid, 'iflag=', iflag
  if (myid .eq. master)  print *, "<----In BFGS_FD. Master just did bfgs cleanup"

  ! Increment nfev by 1 (why? JCS)
  nfev = nfev + 1

!DEC$ IF DEFINED (MPI_OPT)
    ! The master will broadcast its value of iflag to all other workers
    CALL MPI_BCAST(iflag, 1, MPI_INTEGER, master, &
                   MPI_COMM_STEL, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
    ! If the value of iflag >= 0, then master will broadcast its 'fvec'
    ! array to all other workers (length of 'm')
    IF (iflag .ge. 0) CALL &
        MPI_BCAST(fvec, m, MPI_REAL8, master,  &
                  MPI_COMM_STEL, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF


  ! Now that all of the processes have 'fvec', make sure they
  ! also have 'fnorm'
  fnorm = enorm(m, fvec)

  ! Compute value of objective function
   f_curr = fnorm**2
  ! if not using fnorm squrared,  just using fnorm
  ! f_curr = fnorm
  
  if (DEBUG_BFGS .and. (myid .eq. master)) then
    write(*,"(A,1ES12.4E2)") "K=Initial function value: ", f_curr
  end if

    IF (myid .eq. master) then
       WRITE (6, 500) numprocs
500    FORMAT (/,' Beginning BFGS-Finite Difference Iterations', /, &
               ' Number of Processors: ', i6, //, 40('='), / ,2x, &
               'Iteration', 3x, 'Processor', 7x, 'Chi-Sq', 7x, &
               /, 40('='))
        WRITE(6, '(2x,i6,8x,i3,7x,1es12.4)'), 0, myid, fnorm*fnorm
    END IF
 
!DEC$ IF DEFINED (MPI_OPT)
  ! Sync all of the processors at this point - otherwise
  ! some of the processors won't have the data yet...??? (this may not be
  ! neccessary)
      CALL MPI_BARRIER(MPI_COMM_STEL, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF

  CALL fdjac2_mp_queue(fcn, fnorm, m, n, x, fvec, fjac_curr, &
                       m, iflag, nfev, dx_init, fnorm_min, x_min, &
                       fvec_min, fnorm_array, .true.)
  nfev = nfev + n

  ! Check the flag to see if 'flipping' is permitted
  if (enable_flip .eqv. .false.) flip = .false. 
 
  ! Compute gradient of objective function
  ! f_curr = sum(fvec^2)
  ! fvec = (vals(1:m)-targets(1:m))/abs(sigmas(1:m))
  ! fjac = d(fvec)/d(vars)
  ! d(f_curr)/d(vars) = sum(2*fvec*fjac)
  grad_curr = 2*matmul(fvec, fjac_curr)
  !grad_curr = sum(fjac_curr, 1)

  ! to do: print out the grad to a file
  ! to do: calculate the regular hessian, and print it out to a file


  ! Initial Hessian set to multiple of identity
  eye = 0
  DO nvar = 1,n
    eye(nvar,nvar) = 1
  END DO
  hessian = (1.0/beta_hessian)*eye
  
  ! Initial inverse Hessian set to multiple of identity
  invHessian = beta_hessian*eye

  iter = 0
  if (DEBUG_BFGS .and. (myid .eq. master)) then
    print *, 'K====Entering BFGS main loop'
    print *, 'K====iter = ', iter
    print *, 'K====x = ', x
    print *, 'K====f_curr = ', f_curr
    print *, 'K====fvec = ', fvec
    print *, 'K====fjac_curr = ', fjac_curr
    print *, 'K====grad_curr = ', grad_curr
    print *, 'K====norm(grad_curr) = ', enorm(n,grad_curr)
    print *, 'K====Inverse Hessian=', invHessian
    print *, 'K=============================================================='
  end if

  ! Main BFGS loop - See algorithm 6.1 Nocedal & Wright
  do while ((enorm(n,grad_curr)>gtol) .and. ((f_curr*f_curr) > ftol) &
            .and. (nfev < maxfev))
    ! Compute search direction, Eq. 6.18 (or 3.2)
    p_curr = -matmul(invHessian, grad_curr)

    if (DEBUG_BFGS .and. (myid .eq. master)) then
      print *, '<================================================='
      print *, '<----In BFGS main loop'
      print *, '<----iter = ', iter
      print *, '<----x = ', x
      print *, '<----f_curr = ', f_curr
      print *, '<----fvec = ', fvec
      print *, '<----fjac_curr = ', fjac_curr
      print *, '<----grad_curr = ', grad_curr
      print *, '<----norm(grad_curr) = ', enorm(n,grad_curr)
      print *, '<----Inverse Hessian=', invHessian
      print *, '<----p_curr=', p_curr
      print *, '<================================================='
    end if
    
    ! Line search to compute x_new, grad_new, f_new
    CALL BFGS_backtrack(m, n, nfev, fcn, p_curr, f_curr, x, grad_curr, &
                        f_new, x_new, grad_new, alpha_backtrack, &
                        c_armijo, rho_backtrack, beta_hessian, &
                        alpha_min, dx_init, x_min, fvec_min, fnorm_array)

    if (BT_EXIT_FLAG .ne. BT_EXIT_NORMAL) EXIT

    if (DEBUG_BFGS .and. (myid .eq. master)) then
      print *, '<======================================='
      print *, '<----Backtrack exit flag: ', BT_EXIT_FLAG
      print *, '<----nfev=', nfev
      print *, '<======================================='
    end if

    ! Form quantities needed for update to the inverse Hessian 
    ! Eqns 6.5
    s_curr = x_new - x
    y_curr = grad_new - grad_curr
    ! Eqn 6.14
    rho_curr = 1/dot_product(y_curr,s_curr)
    if (DEBUG_BFGS .and. (myid .eq. master)) then
      print *, '<======================================='
      print *, '<----Preparing to update Inverse Hessian'
      print *, '<----x_old=', x
      print *, '<----x_new=', x_new
      print *, '<----s_curr = ', s_curr
      print *, '<----y_curr =', y_curr
      print *, '<----rho_curr = ', rho_curr
      print *, '<======================================='
    end if

    ! Form outer products needed for inverse Hessian update
    s_y_outer = 0
    s_s_outer = 0
    do nvar1 = 1,n
      do nvar2 = 1,n
        s_y_outer(nvar1,nvar2) = s_curr(nvar1)*y_curr(nvar2)
        s_s_outer(nvar1,nvar2) = s_curr(nvar1)*s_curr(nvar2)
      end do
    end do

    ! Rescale if this is the initial inverse Hessian
    !               (see Nocedal & Wright 6.20)
    if (iter==0) then
      !OLD:hessian = (dot_product(y_curr,s_curr) / &
      !           dot_product(y_curr,y_curr)) * eye
      invHessian = (dot_product(y_curr,s_curr) / &
                    dot_product(y_curr,y_curr)) * eye

      ! write the inverse Hessian
      if (myid .eq. master) then
        write(string1,'(A)') '-1'
        string2 = 'invHess.' // TRIM(string1)
        iunit = 28; ii=0
        CALL safe_open(iunit, ii, TRIM(string2), 'new', 'formatted')
        WRITE(iunit,'(2X,i6,2x,i6)'), m, n
        DO jj = 1,m
          WRITE(iunit,'(1p,4e22.14)'), (invHessian(jj,ii), ii=1,n)
        END DO
        CLOSE(iunit)
      end if

    end if

    ! Update approximate inverse Hessian (See Nocedal & Wright 6.17)
    mat1 = eye - rho_curr*s_y_outer
    mat2 = eye - rho_curr*transpose(s_y_outer)
    if (DEBUG_BFGS .and. (myid .eq. master)) then
      print *, '<==========================---------============='
      print *, '<----Updating Inverse Hessian'
      print *, '<----iter = ', iter
      print *, '<----s_y_outer = ', s_y_outer
      print *, '<----s_s_outer = ', s_s_outer
      print *, '<----mat1 = ', mat1
      print *, '<----mat2 = ', mat2
      print *, '<----invHessian = ', invHessian
      print *, '<==========================---------============='
    end if

    invHessian = matmul(mat1,matmul(invHessian,mat2)) + &
                 rho_curr*s_s_outer

    if (DEBUG_BFGS .and. (myid .eq. master)) then
      print *, '<======================================='
      print *, '<----Printing out inverse Hessian'
      print *, '<----invHessian = ', invHessian
      print *, '<======================================='
    end if

    ! write the inverse Hessian
    if (myid .eq. master) then
      write(string1,'(i6.6)') nfev
      string2 = 'invHess.' // TRIM(string1)
      iunit = 28; ii=0
      CALL safe_open(iunit, ii, TRIM(string2), 'new', 'formatted')
      WRITE(iunit,'(2X,i6,2x,i6)'), m,n
      DO jj = 1,m
        WRITE(iunit,'(1p,4e22.14)'), (invHessian(jj,ii), ii=1,n)
      END DO
      CLOSE(iunit)
    end if

    ! Update x, f_curr, and grad_curr
    x = x_new
    f_curr = f_new
    grad_curr = grad_new
    iter = iter + 1
    nfev = nfev+n
  end do

99    IF (ALLOCATED(x_min)) DEALLOCATE(x_min)
      IF (ALLOCATED(fvec_min)) DEALLOCATE(fvec_min)
      IF (ALLOCATED(flip)) DEALLOCATE(flip)
      IF (ALLOCATED(jac_order)) DEALLOCATE(jac_order)
      IF (ALLOCATED(fnorm_array)) DEALLOCATE(fnorm_array)
      IF (ALLOCATED(h_order)) DEALLOCATE(h_order)
  if (myid .eq. master) then
    write(*,"(A,I2,A)") "<----BFGS_FD terminated after ", iter, " iterations."
    write(*,"(A,10E22.14)") "<----New x value: ", x
    write(*,"(A,10E22.14)") "<----Function value: ", f_curr
    write(*,"(A,10E22.14)") "<----Function graidant: ", grad_curr
    write(*,"(A,10E22.14)") "<----New gradient norm: ", enorm(n,grad_curr)
    write(*,"(A,I6)") "<----Function evaluations: ", nfev
  end if
END SUBROUTINE BFGS_FD
