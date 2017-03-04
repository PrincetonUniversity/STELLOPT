      SUBROUTINE flux_init(phipog)
      USE vmec_main
      USE vmec_params, ONLY: signgs
      USE vforces, ONLY : r12=>armn_o, gsqrt=>azmn_o, orsq=>blmn_o
      USE vsvd
      USE realspace
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(*) :: phipog
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: p5 = 0.5_dp, c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: l, js
C-----------------------------------------------

!
!     COMPUTE OBSERVATION POINT - INVARIANT FUNCTIONS OF RADIUS
!     CURRENT = PHIP * INTEGRAL(0,2pi)[ guu / gsqrt]   (on full mesh)
!     RM2     = < R**(-2) >
!     VRM2    = V` * RM2    (V` = 2pi * VP)
!     ORSQ    = SQRT(G) * R**(-2)   (on half mesh)
!
!     MUST HAVE GONE THROUGH NEWPROFILE DETERMINATION OF IOTA AT
!     LEAST ONCE, OTHERWISE IOTAS IS UNDEFINED!
!
      IF (iresidue .le. 0) RETURN
      current(1) = zero
      presint(1) = one

      DO l = 2,nrzt-1
        orsq(l) = p5*( phipog(l) + phipog(l+1) ) *
     1  (ru0(l)*ru0(l) + zu0(l)*zu0(l))
      ENDDO
      DO l = ns,nrzt,ns
        orsq(l) = ( c1p5*phipog(l) - p5*phipog(l-1) ) *
     1  (ru0(l)*ru0(l) + zu0(l)*zu0(l))
      ENDDO

      DO js = 2, ns
         current(js) = twopi*SUM(orsq(js:nrzt:ns)*wint(js:nrzt:ns))
         presint(js) = one
      END DO

      DO l = 2, nrzt
         orsq(l) = gsqrt(l)/r12(l)**2
      END DO
      DO js = 2, ns
         vrm2(js) = twopi*SUM(orsq(js:nrzt:ns)*wint(js:nrzt:ns))
         rm2(js) = vrm2(js)/(twopi*signgs*vp(js))
         ovrm2(js) = one/vrm2(js)
         ochip(js) = one/(phip(js)*iotas(js))
         presph(js) = presf(js) - presf(js - 1)
      END DO

      END SUBROUTINE flux_init
