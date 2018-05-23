
MODULE jacfcn_mod

			USE stel_kinds
			USE fdjac_mod, ONLY: jac_index, jac_err, n_red, jac_order, ix_min, h_order, flip
			USE safe_open_mod, ONLY: safe_open
			USE	vparams, ONLY: sfincs_nmax, sfincs_mmax, ndatafmax

			! Value of Jacobian obtained from the most recent call to fcn
			REAL(rprec), DIMENSION(:,:), allocatable :: fjac_curr
			! Are analytic gradients avaiable? (maybe should have same dimension as
			! fjac_curr eventually
			LOGICAL :: jac_analytic = .FALSE.

!-----------------------------------------------------------------------
!     Subroutine:    stellopt_jacfcn
!     Authors:       E. J. Paul (ejpaul@umd.edu)
!     Date:          05/17/2018
!     Description:   This subroutine calculates the jacobian matrix for the
!										 function which is minimized by STELLOPT. This is
! 									 intended to be used with the lmanalytic subroutine
!										 in the case that the jacobian matrix is not computed
!										 purely with forward differencing. When this subroutine
!										 is called, fjac should contain the jacobian matrix
! 									 which is reordered here.
!-----------------------------------------------------------------------
		CONTAINS

      SUBROUTINE jacfcn(m, n, x, fvec, fjac, x_min, fnorm, ncnt, wa, fnorm_min, epsfcn)

			EXTERNAL dpmpar

!-----------------------------------------------------------------------
!   I've used the same naming convention as fdjac_mod
!
!			Input Variables
!        m       is a positive INTEGER input variable set to the number
!         			 of functions.
!        n       is a positive INTEGER input variable set to the number
!         			 of variables. n must not exceed m.
!        x       is an input array of length n.
!        fvec    is an input array of length m which must contain the
!         			 functions evaluated at x.
!				 fjac		 jacobian of the function.
! 			 ncnt		 number of function evaluations
!				 fnorm	 norm of fcn at point x
!				 epsfcn  is an input variable used in determining a suitable
!         			 step length for the forward-difference approximation. this
!         			 approximation assumes that the relative errors in the
!         			 functions are of the order of epsfcn. If epsfcn is less
!         			 than the machine precision, it is assumed that the relative
!         			 errors in the functions are of the order of the machine
!         			 precision.
!
!			Output Variables
!				 x_min	 finite diff direction with minimal fcnß
!			   wa		   working array of length m.
!				 fnorm_min norm of fcn evaluated at x_min
!----------------------------------------------------------------------

			! Input
      INTEGER, INTENT(in)      ::  m, n
      INTEGER, INTENT(inout)   :: ncnt
      REAL(rprec), INTENT(inout), DIMENSION(n)  :: x
			REAL(rprec), INTENT(inout), DIMENSION(m,n) :: fjac
			REAL(rprec), INTENT(in), DIMENSION(m) :: fvec
			REAL(rprec), INTENT(in) :: fnorm
			REAL(rprec), INTENT(out) :: wa(m), x_min(n)
			REAL(rprec), INTENT(out) :: fnorm_min
			REAL(rprec) ::  epsfcn

		 ! Local variables
			INTEGER :: iunit, j, i, ix_temp
			CHARACTER(16) ::  temp_string
			CHARACTER(256) ::  jac_file
			REAL(rprec) :: temp
			REAL(rprec), DIMENSION(n) :: jacnorm_array
			logical, dimension(n) :: lmask
			integer, dimension(1) :: isort
			REAL(rprec), DIMENSION(n) :: h
			REAL(rprec) :: eps, epsmch, temp_norm
			REAL(rprec), DIMENSION(m,n) :: fvec_array
			REAL(rprec), DIMENSION(n) :: fnorm_array
			INTEGER :: ierr_mpi

			! Use jacobian from current iteration
			fjac = fjac_curr

			! Instead of using fnorm_array, we could look at norm of colomns of Jacobian
			! for sorting
			!			jacnorm_array = SQRT(SUM(fjac*fjac,DIM=1))

			! h needs to be computed for fvec_array
			epsmch = dpmpar(1)
      eps = SQRT(MAX(epsfcn,epsmch))
      h(1:n) = eps*ABS(x(1:n))
      WHERE (h .eq. zero) h=eps
			WHERE (flip) h = -h

			! fvec_array can be computed from fjac
			DO nvar = 1,n
				DO mtarget = 1,m
					fvec_array(mtarget,nvar) = fjac(mtarget,nvar)*h(nvar) + fvec(mtarget)
				END DO
			END DO
			fnorm_array = SQRT(SUM(fvec_array*fvec_array,DIM=1))

			! The following is modified from fdjac_mod.f
		  ! Check fjac for errors
			DO i = 1, n
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
!     Output Function Evaluations, Jacobian, and reorder jacobian
!
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

!			! If no errors in computing jac
!			jac_err = 1
!			n_red = n

      DO i = 1, n
         jac_index(i) = i
      END DO
		! The following should not modify Jacobian if free of errors
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
!     Find minimum and reorder
!
      jac_order = 0
      ix_temp = 1
      temp = 0
      lmask = .true.
      DO WHILE (ix_temp <= n_red)
         isort = MINLOC(fnorm_array, MASK=lmask) ! find element w/ min. value
!				 isort = MINLOC(jacnorm_array, MASK=lmask)
         temp = jacnorm_array(isort(1)) ! min value
         jac_order(ix_temp) = isort(1)
         IF(isort(1) <= 0 .or. isort(1) > n_red) THEN
            EXIT
				 ! fnorm_array - norm of forward diff step function evaluation
				 ! fnorm - norm of function evaluation
				 ! If smallest element of fnorm_array is larger than fnorm
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

			END SUBROUTINE jacfcn

END MODULE jacfcn_mod
