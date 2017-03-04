      SUBROUTINE setup_booz(ntorsum, ns, mnmax, ohs, 
     1   xmb, xnb, sfull, scl, mboz, nboz, mnboz, nu2_b, 
     2   nu_boz, nv_boz, nfp, lasym)
      USE stel_kinds
      USE booz_persistent, ONLY: xm
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: ns, mnmax, mboz, nboz, nu2_b, 
     1                       nu_boz, nv_boz, nfp
      INTEGER, INTENT(inout) :: mnboz
      REAL(rprec) :: ohs
      INTEGER, INTENT(out), DIMENSION(0:1)   :: ntorsum
      REAL(rprec), INTENT(out),  DIMENSION(mnboz) :: xmb, xnb, scl
      REAL(rprec), INTENT(out), DIMENSION(ns) :: sfull
      LOGICAL, INTENT(in) :: lasym
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mnboz0, n2, m, n1, n, i
      REAL(rprec) :: fac, hs
C-----------------------------------------------
!
!     SETUP BOOZER M,N ARRAYS
!
      mnboz0 = 0
      n2 = nboz
      DO m = 0, mboz-1
         n1 = -nboz
         IF (m .eq. 0) n1 = 0
         DO n = n1, n2
            mnboz0 = mnboz0 + 1
            IF (mnboz0 .gt. mnboz) THEN
               STOP 'mnboz exceeds limit in booz xform'
            END IF
            xnb(mnboz0) = n*nfp
            xmb(mnboz0) = m
         END DO
      END DO

      IF (mnboz0 .ne. mnboz) mnboz = mnboz0

!     SCALE FACTOR FOR NORMALIZATION OF FOURIER TRANSFORMS

      IF (lasym) THEN
         fac = 2.0_dp/(nu_boz*nv_boz)
      ELSE
         fac = 2.0_dp/((nu2_b-1)*nv_boz)
      END IF

      scl(:mnboz) = fac

      WHERE (NINT(xnb(:mnboz)).eq.0 .and. NINT(xmb(:mnboz)).eq.0)
     1   scl(:mnboz) = fac/2

      ntorsum(0) = 0
      ntorsum(1) = 0

      DO i=1,mnmax
        IF (NINT(xm(i)) .eq. 0) ntorsum(0) = ntorsum(0)+1
        IF (NINT(xm(i)) .le. 1) ntorsum(1) = ntorsum(1)+1
      END DO

      ohs = (ns-1)
      hs = one/ohs
      DO i=2,ns
        sfull(i) = SQRT(hs*(i-1))
      END DO

      END SUBROUTINE setup_booz
