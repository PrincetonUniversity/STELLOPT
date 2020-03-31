
SUBROUTINE calccur(n,t,ht,y)
! **********************************************************************
! Modules
! **********************************************************************
  USE neo_precision
  USE neo_exchange
  USE neo_output
  USE neo_units
  USE partpa_cur
! **********************************************************************
! Local definitions
! **********************************************************************
  IMPLICIT NONE
! **********************************************************************
  REAL(kind=dp), DIMENSION(ndim_cur),  INTENT(in), TARGET ::  y
  REAL(kind=dp), DIMENSION(npart_cur), INTENT(in)         ::  t
  REAL(kind=dp),                       INTENT(in)         ::  ht
  INTEGER                                                 ::  i,n
  REAL(kind=dp),               POINTER                    ::  p_theta, p_bm2
  REAL(kind=dp),               POINTER                    ::  p_bm2gv
  REAL(kind=dp),               POINTER                    ::  p_lamps
  REAL(kind=dp),               POINTER                    ::  p_lamb1n, p_lamb1d
  REAL(kind=dp),               POINTER                    ::  p_lamb2n, p_lamb2d
  REAL(kind=dp), DIMENSION(:), POINTER                    ::  p_l, p_k1, p_k
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
  lambda_b = 0
! **********************************************************************
  avnabpsi = p_bm2gv / p_bm2
  DO i = 1, npart_cur
     lambda_b = lambda_b + t(i)**(alpha_cur-1) * y_part(i)*y_part(i) *   &
                p_k(i) / p_l(i)
  END DO
  lambda_b   = - 3.0_dp / 8.0_dp * lambda_b * alpha_cur * ht
!
  lambda_ps1 = 2.0_dp * p_lamb1n / p_lamb1d
  lambda_ps2 = 2.0_dp * p_lamb2n / p_lamb2d
!
  lambda_b1  = lambda_b + lambda_ps1
  lambda_b2  = lambda_b + lambda_ps2
! **********************************************************************
  IF (write_cur_inte .EQ. 1) THEN
     WRITE(w_u8,'(1(1x,i8),23(1x,e17.10))')                            &
                 n,                                                    & ! 1
                 avnabpsi,lambda_b,                                    & ! 2-3
                 lambda_ps1,lambda_ps2,                                & ! 4-5
                 lambda_b1,lambda_b2,                                  & ! 6-7
                 p_lamps,p_lamb1n,p_lamb1d,p_lamb2n,p_lamb2d,          & ! 8-12
                 p_k(1)/p_l(1), p_l(1),p_k(1),p_k1(1),                 & ! 13-16
                 p_k(50)/p_l(50), p_l(50),p_k(50),p_k1(50),            & ! 17-20
                 p_k(100)/p_l(100), p_l(100),p_k(100),p_k1(100)          ! 21-24
  END IF
END SUBROUTINE calccur
