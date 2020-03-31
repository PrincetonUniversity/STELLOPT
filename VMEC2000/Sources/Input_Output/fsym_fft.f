      SUBROUTINE fext_fft (bout, bs_s, bs_a)
      USE vmec_main
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(nzeta,ntheta3), INTENT(out) :: bout
      REAL(rprec), DIMENSION(nzeta,ntheta2), INTENT(in) ::  bs_s, bs_a
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ir, i, kz, kzr
C-----------------------------------------------
!     
!     Extends bs from ntheta2 interval to full ntheta3 interval in angle u
!     bs_s ~ cos(mu-nv) (symmetric);   bs_a ~ sin(mu-nv) (anti-symmetric)
!     ntheta2 = pi
!
      bout(:,1:ntheta2) = bs_s(:,1:ntheta2) + bs_a(:,1:ntheta2)
      DO i = 1+ntheta2, ntheta3
         ir = ntheta1+2-i                     !-theta
         DO kz= 1, nzeta
!            kzr = ireflect(kz*ns)/ns          !-zeta
            kzr = nzeta+2-kz
            IF (kz .eq. 1) kzr = 1
            bout(kz,i) = bs_s(kzr,ir) - bs_a(kzr,ir)
         END DO
      END DO

      END SUBROUTINE fext_fft


      SUBROUTINE fsym_fft (bs, bu, bv, bs_s, bu_s, bv_s,
     1    bs_a, bu_a, bv_a)
      USE vmec_main
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(nzeta,ntheta3), INTENT(in)      :: bs
      REAL(rprec), DIMENSION(nzeta,ntheta3,0:1), INTENT(in)  :: bu, bv
      REAL(rprec), DIMENSION(nzeta,ntheta2,0:1), INTENT(out) ::
     1   bu_s, bv_s, bu_a, bv_a
      REAL(rprec), DIMENSION(nzeta,ntheta2) :: bs_s, bs_a
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ir, i, kz, kzr
C-----------------------------------------------
!       CONTRACTS bs,bu,bv FROM FULL nu INTERVAL TO HALF-U INTERVAL
!       SO COS,SIN INTEGRALS CAN BE PERFORMED ON HALF-U INTERVAL
!
!       bs_s(v,u) = .5*( bs(v,u) - bs(-v,-u) )     ! * SIN(mu - nv)
!       bs_a(v,u) = .5*( bs(v,u) + bs(-v,-u) )     ! * COS(mu - nv)
!
!       bu, bv have opposite parity       
!
      DO i = 1, ntheta2
         ir = ntheta1+2-i                     !-theta
         IF (i == 1) ir = 1
         DO kz = 1, nzeta
!           kzr = ireflect(ns*kz)/ns          !-zeta
            kzr = nzeta+2-kz
            IF (kz .eq. 1) kzr = 1
            bs_a(kz,i)      = cp5*(bs(kz,i)+bs(kzr,ir))
            bs_s(kz,i)      = cp5*(bs(kz,i)-bs(kzr,ir))
            bu_a(kz,i,:) = cp5*(bu(kz,i,:)-bu(kzr,ir,:))
            bu_s(kz,i,:) = cp5*(bu(kz,i,:)+bu(kzr,ir,:))
            bv_a(kz,i,:) = cp5*(bv(kz,i,:)-bv(kzr,ir,:))
            bv_s(kz,i,:) = cp5*(bv(kz,i,:)+bv(kzr,ir,:))
            END DO
         END DO

      END SUBROUTINE fsym_fft
