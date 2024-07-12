
SUBROUTINE rhs_cur1(phi,y,dery)
! Right Hand Side of Differential Equation
!
!   y(1)                     - $\theta$
!   y(2)                     - $\int \rd \phi / B^2$
!   y(3)                     - $\int \rd \phi |\nabla \psi| / B^2$
!   y(4)                     - $\int \rd \phi K_G / B^3$
!   y(4+1)-y(4+npart)        - $\int \rd \phi \frac{\rd I_{fj}}{rd \phi}$
!   y(4+npart)-y(4+2*npart)  - $\int \rd \phi \frac{\rd H_{fj}}{rd \phi}$
!
  USE neo_precision
  USE partpa_cur
  USE neo_exchange
!
  IMPLICIT NONE
!
  REAL(kind=dp),               INTENT(in)          ::  phi
  REAL(kind=dp)                                    ::  theta,bmod,gval,qval
  REAL(kind=dp)                                    ::  bmodm2,bmodm3
  REAL(kind=dp)                                    ::  geodcu,pardeb,bra
  REAL(kind=dp)                                    ::  curfac
  REAL(kind=dp), DIMENSION(:)      , TARGET, INTENT(in)  ::  y
  REAL(kind=dp), DIMENSION(SIZE(y)), TARGET, INTENT(out) ::  dery
  REAL(kind=dp),               POINTER             ::  pd_iota, pd_bm2
  REAL(kind=dp),               POINTER             ::  pd_bm2gv
  REAL(kind=dp),               POINTER             ::  pd_lamps
  REAL(kind=dp),               POINTER             ::  pd_lamb1n, pd_lamb1d
  REAL(kind=dp),               POINTER             ::  pd_lamb2n, pd_lamb2d
  REAL(kind=dp), DIMENSION(:), POINTER             ::  pd_l, pd_k1, pd_k
  REAL(kind=dp), DIMENSION(:), POINTER             ::  p_l, p_k1, p_k
  INTEGER      , PARAMETER                         ::  k_meth = 2
!
  pd_iota   => dery(1)
  pd_bm2    => dery(2)
  pd_bm2gv  => dery(3)
  pd_lamps  => dery(4)
  pd_lamb1n => dery(5)
  pd_lamb1d => dery(6)
  pd_lamb2n => dery(7)
  pd_lamb2d => dery(8)
  pd_l      => dery(npq_cur+1:npq_cur+npart_cur)
  pd_k1     => dery(npq_cur+npart_cur+1:npq_cur+2*npart_cur)
  pd_k      => dery(npq_cur+2*npart_cur+1:npq_cur+3*npart_cur)

  p_l       => y(npq_cur+1:npq_cur+npart_cur)
  p_k1      => y(npq_cur+npart_cur+1:npq_cur+2*npart_cur)
  p_k       => y(npq_cur+2*npart_cur+1:npq_cur+3*npart_cur)

  theta=y(1)
! Evaluation of Splines
  CALL neo_eval(theta,phi,bmod,gval,geodcu,pardeb,qval)
!
  bra = bmod / bmod0
  bmodm2 = ONE / bra**2
  bmodm3 = bmodm2 / bra
  curfac = geodcu * fac / bmod0

!!$  PRINT *, 'phi ',phi
!!$  PRINT *, 'theta ',theta
!!$  PRINT *, 'bra ',bra
!!$  PRINT *, 'bmodm2 ',bmodm2
!!$  PRINT *, 'bmodm3 ',bmodm3
!!$  PRINT *, 'curfac ',curfac
!!$  PRINT *, 'geodcu ',geodcu
!!$  PRINT *, 'fac ',fac
!!$  PRINT *, 'bmod0 ',bmod0
!!$  PRINT *, 'bmod ',bmod
!!$  PRINT *, 'gvql ',gval
!!$  PRINT *, 'qvql ',qval

  yfac   = ONE - y_part*bra
  sqyfac = SQRT(yfac)
!
  pd_iota   = iota(psi_ind)
  pd_bm2    = bmodm2
  pd_bm2gv  = bmodm2 * gval
  pd_lamps  = curfac * bmodm3
  pd_lamb1n = qval * bmodm2 * pd_lamps
  pd_lamb1d = qval * bmodm2
  pd_lamb2n = pd_lamps
  pd_lamb2d = 1

  yfac      = ONE - y_part*bra
  sqyfac    = SQRT(yfac)
  pd_l      = bmodm2 * sqyfac
  IF (k_meth .EQ. 1) THEN
     pd_k1     = curfac / yfac / sqyfac
     pd_k      = pd_l * p_k1
  ELSE
     k_fac1 = -(2.0_dp/bra + y_part / 2.0_dp / yfac  ) * pardeb / bmod0
     k_fac2 = bmodm2 / yfac * curfac
     pd_k1  = k_fac1 * p_k1 + k_fac2
     pd_k   = p_k1
  END IF
!
  RETURN
END SUBROUTINE rhs_cur1
