      SUBROUTINE surface(rc, rs, zs, zc, xm, xn, mnmax)
      USE vacmod
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER mnmax
      REAL(rprec), DIMENSION(mnmax) :: rc, rs, zs, zc, xm, xn
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
      INTEGER :: i, mn, m, n, n1
      REAL(rprec), ALLOCATABLE, DIMENSION(:) ::
     1   ruu, ruv, rvv, zuu, zuv, zvv
      REAL(rprec) :: cosmn1, sinmn1
!-----------------------------------------------
!
!       THIS ROUTINE COMPUTES THE SURFACE VALUES OF R,Z AND DERIVATIVES
!
!
!       Compute R & Z (and their derivatives) on surface
!
!       R = SUM [RC(m,n)*COS(mu - nv) + RS(m,n)*SIN(mu - nv)]
!       Z = SUM [ZS(m,n)*SIN(mu - nv) + ZC(m,n)*COS(mu - nv)]
!
!       NOTE: u, v here are actual angles (0, 2pi), NOT the normalized
!             variables used in PKM paper
!
      ALLOCATE (ruu(nuv2), ruv(nuv2), rvv(nuv2), zuu(nuv2), zuv(nuv2),
     1          zvv(nuv2), stat = i)
      IF (i .NE. 0) STOP 'Allocation error in SURFACE'

      r1b = 0;   rub = 0;   rvb = 0;  ruu = 0; ruv = 0; rvv = 0
      z1b = 0;   zub = 0;   zvb = 0;  zuu = 0; zuv = 0; zvv = 0

      DO mn = 1, mnmax
         m = NINT(xm(mn))
         n = NINT(xn(mn)/(nfper))
         n1 = ABS(n)
         DO i = 1, nuv2
            cosmn1 = cosu1(i,m)*cosv1(i,n1) + csign(n)*sinu1(i,m)*
     1               sinv1(i,n1)
            sinmn1 = sinu1(i,m)*cosv1(i,n1) - csign(n)*cosu1(i,m)*
     1               sinv1(i,n1)
            r1b(i) = r1b(i) + rc(mn) * cosmn1
            rub(i) = rub(i) - xm(mn) * rc(mn) * sinmn1
            rvb(i) = rvb(i) + xn(mn) * rc(mn) * sinmn1
            z1b(i) = z1b(i) + zs(mn) * sinmn1
            zub(i) = zub(i) + xm(mn) * zs(mn) * cosmn1
            zvb(i) = zvb(i) - xn(mn) * zs(mn) * cosmn1
            ruu(i) = ruu(i) - xm(mn)*xm(mn)*rc(mn) * cosmn1
            ruv(i) = ruv(i) + xm(mn)*xn(mn)*rc(mn) * cosmn1
            rvv(i) = rvv(i) - xn(mn)*xn(mn)*rc(mn) * cosmn1
            zuu(i) = zuu(i) - xm(mn)*xm(mn)*zs(mn) * sinmn1
            zuv(i) = zuv(i) + xm(mn)*xn(mn)*zs(mn) * sinmn1
            zvv(i) = zvv(i) - xn(mn)*xn(mn)*zs(mn) * sinmn1
         IF (lasym) THEN
            r1b(i) = r1b(i) + rs(mn) * sinmn1
            rub(i) = rub(i) + xm(mn) * rs(mn) * cosmn1
            rvb(i) = rvb(i) - xn(mn) * rs(mn) * cosmn1
            z1b(i) = z1b(i) + zc(mn) * cosmn1
            zub(i) = zub(i) - xm(mn) * zc(mn) * sinmn1
            zvb(i) = zvb(i) + xn(mn) * zc(mn) * sinmn1
            ruu(i) = ruu(i) - xm(mn)*xm(mn)*rs(mn) * sinmn1
            ruv(i) = ruv(i) + xm(mn)*xn(mn)*rs(mn) * sinmn1
            rvv(i) = rvv(i) - xn(mn)*xn(mn)*rs(mn) * sinmn1
            zuu(i) = zuu(i) - xm(mn)*xm(mn)*zc(mn) * cosmn1
            zuv(i) = zuv(i) + xm(mn)*xn(mn)*zc(mn) * cosmn1
            zvv(i) = zvv(i) - xn(mn)*xn(mn)*zc(mn) * cosmn1
         END IF
      END DO
      END DO

!
!     COMPUTE METRIC COEFFICIENTS GIJ_B AND SURFACE NORMAL COMPONENTS
!     [SNR, SNV, SNZ] = NP*[Xu cross Xv]
!
!     NOTE: These should be multiplied by -signgs to point OUTWARD from vacuum INTO plasma 
!           for either handed-ness of the coordinate system
!
!           Eq. 2.4 in PKM has wrong sign for a left-handed coordinate system
!
!     NOTE: guv = .5*np guv_b; gvv = np*np* gvv_b, where GUV, GVV are the
!           REAL metric elements. CAP(A), etc. defined in Eq. (2.13) of PKM paper
!
!           AUU == NP*CAP(A) = .5*Xuu dot [Xu cross Xv] * NP
!
!           AUV == 2*NP*CAP(B) =  Xuv dot [Xu cross Xv] * NP
!
!           AVV == NP*CAP(C) = .5*Xvv dot [Xu cross Xv] * NP
!
      DO i = 1,nuv2
        guu_b(i) = rub(i)*rub(i) + zub(i)*zub(i)
        guv_b(i) = (rub(i)*rvb(i)+ zub(i)*zvb(i))*onp*2
        gvv_b(i) = (rvb(i)*rvb(i)+ zvb(i)*zvb(i)+(r1b(i)*r1b(i)))*onp2
        rzb2(i) = r1b(i)*r1b(i) + z1b(i)*z1b(i)
        snr(i) = signgs*r1b(i)*zub(i)
        snv(i) = signgs*(rub(i)*zvb(i) - rvb(i)*zub(i))
        snz(i) =-signgs*r1b(i)*rub(i)
        drv(i) = -(r1b(i)*snr(i) + z1b(i)*snz(i))
        auu(i) = p5*(snr(i)*ruu(i) + snz(i)*zuu(i))
        auv(i) = (snr(i)*ruv(i) + snv(i)*rub(i) + snz(i)*zuv(i))*onp
        avv(i) = (snv(i)*rvb(i) + p5*(snr(i)*(rvv(i) - r1b(i))
     1                          +     snz(i)* zvv(i)))*onp2
      END DO
      IF (.NOT.lasym) THEN
         DO i = 1 + nv, nuv2 - nv
            rzb2(imirr(i)) = rzb2(i)
            r1b(imirr(i))  = r1b(i)
            z1b(imirr(i))  =-z1b(i)
         END DO
      END IF

      DO i = 1,nuv
        rcosuv(i) = r1b(i)*cosuv(i)
        rsinuv(i) = r1b(i)*sinuv(i)
      END DO

      DEALLOCATE (ruu, ruv, rvv, zuu, zuv, zvv, stat=i)

      END SUBROUTINE surface
