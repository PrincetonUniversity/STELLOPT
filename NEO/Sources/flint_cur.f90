
SUBROUTINE flint_cur()
! **********************************************************************
! Modules
! **********************************************************************
  USE neo_precision
  USE neo_units
  USE neo_exchange
  USE neo_output
  USE partpa_cur
  USE neo_rk4dcur
! **********************************************************************
! Local definitions
! **********************************************************************
  IMPLICIT NONE
! **********************************************************************
  INTEGER                                          ::  i,n,j1,nfl
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE, TARGET ::  y
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE         ::  y_s
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE         ::  t

  REAL(kind=dp)                                    ::  theta0, phi0
  REAL(kind=dp)                                    ::  theta, phi
  REAL(kind=dp)                                    ::  hphi
  REAL(kind=dp)                                    ::  tmin, tmax, ht

  REAL(kind=dp),               POINTER             ::  p_theta, p_bm2
  REAL(kind=dp),               POINTER             ::  p_bm2gv
  REAL(kind=dp),               POINTER             ::  p_lamps
  REAL(kind=dp),               POINTER             ::  p_lamb1n, p_lamb1d
  REAL(kind=dp),               POINTER             ::  p_lamb2n, p_lamb2d
  REAL(kind=dp), DIMENSION(:), POINTER             ::  p_l, p_k1, p_k
! **********************************************************************
! Allocation
! **********************************************************************
  ndim_cur = npq_cur + 3 * npart_cur
  ALLOCATE( t(npart_cur) )
  ALLOCATE( y_part(npart_cur) )
  ALLOCATE( yfac(npart_cur) )
  ALLOCATE( sqyfac(npart_cur) )
  ALLOCATE( y(ndim_cur) )
  ALLOCATE( y_s(ndim_cur) )
  ALLOCATE( k_fac1(npart_cur) )
  ALLOCATE( k_fac2(npart_cur) )
! **********************************************************************
! Pointer
! **********************************************************************
  p_theta  => y(1)
  p_bm2    => y(2)
  p_bm2gv  => y(3)
  p_lamps  => y(4)
  p_lamb1n => y(5)
  p_lamb1d => y(6)
  p_lamb2n => y(7)
  p_lamb2d => y(8)
  p_l      => y(npq_cur+1:npq_cur+npart_cur)
  p_k1     => y(npq_cur+npart_cur+1:npq_cur+2*npart_cur)
  p_k      => y(npq_cur+2*npart_cur+1:npq_cur+3*npart_cur)
! **********************************************************************
! Initial settings
! **********************************************************************
  ierr   = 0
  hphi   = twopi/(nstep_per*nper)
  bmod0  = bmref
!
  theta0 = theta_bmax
  phi0   = phi_bmax
! **********************************************************************
! Definition of t-values for particles
! **********************************************************************
  tmin = 1.0e-3_dp
  tmax = 1
  ht = (tmax-tmin)/(npart_cur-1)
  DO i=1,npart_cur
     t(i)=tmin+ht*(i-1)
  ENDDO
  y_part = ONE - t**alpha_cur
!!$  PRINT *, 'hphi, ht ',hphi, ht
!!$  PRINT *, 'bmod0 ', bmod0
!!$  PRINT *, 'ndim_cur ', ndim_cur
!!$  PRINT *, 'theta0 ', theta0
!!$  PRINT *, 'phi0 ', phi0
!!$  PAUSE
!!$  PRINT *, 't ',t
!!$  PAUSE
!!$  PRINT *, 'y_part ',y_part
!!$  PAUSE
! **********************************************************************
  IF (write_cur_inte == 1) THEN
     OPEN(unit=w_u8,file='current.dat',status='replace',form='formatted')
  ENDIF
! **********************************************************************
  phi     = phi0
  y       = 0
  p_theta = theta0
! **********************************************************************
! Integration steps (Summation)
! **********************************************************************
! ATTENTION JUST FOR TEST
  hit_rat = 0
  nintfp  = 10000
! ATTENTION JUST FOR TEST
  IF (hit_rat .EQ. 0) THEN
     DO n=1,nintfp
        DO j1=1,nstep_per
           CALL RK4D_CUR(phi,y,hphi)
        END DO
        IF (write_cur_inte == 1) THEN
           CALL calccur(n,t,ht,y)
        ENDIF
     END DO
     IF (write_cur_inte == 0) THEN
        CALL calccur(n,t,ht,y)
     END IF
  ELSE
    y_s = 0
    DO nfl = 0, nfl_rat
       phi     = phi0
       y       = 0
       theta   = theta0 + nfl * delta_theta_rat
       p_theta = theta
       DO n=1,nfp_rat
          DO j1=1,nstep_per
             CALL RK4D_CUR(phi,y,hphi)
          END DO
          IF (write_cur_inte == 1) THEN
             CALL calccur(nfl*nfp_rat+n,t,ht,y)
          ENDIF
       END DO
       y_s = y_s + y
    END DO
    y = y_s
    CALL calccur(nfl*nfp_rat+n+1,t,ht,y)
  END IF
! **********************************************************************
  lambda_b   = lambda_b   / avnabpsi
  lambda_b1  = lambda_b1  / avnabpsi
  lambda_b2  = lambda_b2  / avnabpsi
  lambda_ps1 = lambda_ps1 / avnabpsi
  lambda_ps2 = lambda_ps2 / avnabpsi
! **********************************************************************
  IF (write_cur_inte == 1) THEN
     CLOSE(unit=w_u8)
  ENDIF
! **********************************************************************
! Deallocation
! **********************************************************************
  DEALLOCATE( t )
  DEALLOCATE( y_part )
  DEALLOCATE( yfac )
  DEALLOCATE( sqyfac )
  DEALLOCATE( y )
  DEALLOCATE( y_s )
  DEALLOCATE( k_fac1 )
  DEALLOCATE( k_fac2 )
!
  RETURN
END SUBROUTINE flint_cur
