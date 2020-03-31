
  SUBROUTINE rhs_bo1(phi,y,dery)
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
    USE partpa_bo
    USE neo_exchange
!
    IMPLICIT NONE
!
    REAL(kind=dp),               INTENT(in)          ::  phi
    REAL(kind=dp)                                    ::  theta,bmod,gval,qval
    REAL(kind=dp)                                    ::  bmodm2,bmodm3
    REAL(kind=dp)                                    ::  geodcu,pardeb,bra,subsq,sq,sqeta
    REAL(kind=dp), DIMENSION(:),         INTENT(in)  ::  y
    REAL(kind=dp), DIMENSION(SIZE(y)), TARGET, INTENT(out) ::  dery
    INTEGER                                          ::  ipass,i
    REAL(kind=dp),               POINTER             ::  p_iota, p_bm2
    REAL(kind=dp),               POINTER             ::  p_bm2gv, p_bm3ge
    REAL(kind=dp), DIMENSION(:), POINTER             ::  p_i, p_h
!
    p_iota  => dery(1)
    p_bm2   => dery(2)
    p_bm2gv => dery(3)
    p_bm3ge => dery(4)
    p_i     => dery(npq+1:npq+npart)
    p_h     => dery(npq+npart+1:npq+2*npart)

    theta=y(1)
! Evaluation of Splines
    CALL neo_eval(theta,phi,bmod,gval,geodcu,pardeb,qval)
!
    bmodm2 = ONE / bmod**2
    bmodm3 = bmodm2 / bmod
    bra = bmod / bmod0
!
    IF(pardeb*pard0.LE.ZERO .AND. pardeb.GT.ZERO) THEN
       ipass=1
    ELSE
       ipass=0
    ENDIF
    IF(ipmax.EQ.0.AND.pardeb*pard0.LE.ZERO.AND.pardeb.LT.ZERO) ipmax=1
    pard0=pardeb
!
    p_iota  = iota(psi_ind)
    p_bm2   = bmodm2
    p_bm2gv = bmodm2 * gval
    p_bm3ge = geodcu * bmodm3
!
    DO i=1,npart
       subsq = ONE - bra/eta(i)
       sqeta = sqrt(eta(i))
       IF(subsq.GT.ZERO) THEN
          isw(i) = 1
          icount(i) = icount(i) + 1
          ipa(i) = ipa(i)+ipass
          sq = sqrt(subsq)*bmodm2
          p_i(i) = sq
          p_h(i) = sq*(4.0_dp/bra-ONE/eta(i))*geodcu/sqeta
       ELSE
          sq = 0
          IF(isw(i).EQ.1) THEN
             isw(i) = 2
          ELSEIF(isw(i).EQ.2) THEN
             CONTINUE
          ELSE
             isw(i) = 0
          ENDIF
          p_i(i) = 0
          p_h(i) = 0
       ENDIF
    ENDDO
!
    RETURN
  END SUBROUTINE rhs_bo1
