      SUBROUTINE surface(rc, rs, zs, zc, xm, xn, mnmax)
      USE vacmod
      USE parallel_include_module
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER mnmax
      REAL(dp), DIMENSION(mnmax) :: rc, rs, zs, zc, xm, xn
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, mn, m, n, n1
      REAL(dp), ALLOCATABLE, DIMENSION(:) ::
     1   ruu, ruv, rvv, zuu, zuv, zvv, cosmn1, sinmn1
      REAL(dp) :: tsurfon, tsurfoff
C-----------------------------------------------
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
      CALL second0(tsurfon)

      ALLOCATE (ruu(nuv3), ruv(nuv3), rvv(nuv3), zuu(nuv3), zuv(nuv3),
     1          zvv(nuv3), cosmn1(nuv3), sinmn1(nuv3), stat = i)
      IF (i .ne. 0) STOP 'Allocation error in SURFACE'

      r1b = 0; z1b = 0
      DO i = nuv3min, nuv3max  
         zub(i) = 0;   zvb(i) = 0;  zuu(i) = 0; zuv(i) = 0; zvv(i) = 0
         rub(i) = 0;   rvb(i) = 0;  ruu(i) = 0; ruv(i) = 0; rvv(i) = 0
      END DO

      DO mn = 1, mnmax
         m = NINT(xm(mn))
         n = NINT(xn(mn)/(nfper))
         n1 = ABS(n)
         cosmn1(:) = cosu1(:,m)*cosv1(:,n1) + csign(n)*sinu1(:,m)*
     1               sinv1(:,n1)
         sinmn1(:) = sinu1(:,m)*cosv1(:,n1) - csign(n)*cosu1(:,m)*
     1               sinv1(:,n1)
         DO i = 1, nuv3
            r1b(i) = r1b(i) + rc(mn) * cosmn1(i)
            z1b(i) = z1b(i) + zs(mn) * sinmn1(i)
         END DO
         DO i = nuv3min, nuv3max
            rub(i) = rub(i) - xm(mn) * rc(mn) * sinmn1(i)
            rvb(i) = rvb(i) + xn(mn) * rc(mn) * sinmn1(i)
            zub(i) = zub(i) + xm(mn) * zs(mn) * cosmn1(i)
            zvb(i) = zvb(i) - xn(mn) * zs(mn) * cosmn1(i)
            ruu(i) = ruu(i) - xm(mn)*xm(mn)*rc(mn) * cosmn1(i)
            ruv(i) = ruv(i) + xm(mn)*xn(mn)*rc(mn) * cosmn1(i)
            rvv(i) = rvv(i) - xn(mn)*xn(mn)*rc(mn) * cosmn1(i)
            zuu(i) = zuu(i) - xm(mn)*xm(mn)*zs(mn) * sinmn1(i)
            zuv(i) = zuv(i) + xm(mn)*xn(mn)*zs(mn) * sinmn1(i)
            zvv(i) = zvv(i) - xn(mn)*xn(mn)*zs(mn) * sinmn1(i)
         END DO

         IF (.NOT.lasym) CYCLE

         DO i = 1, nuv3
            r1b(i) = r1b(i) + rs(mn) * sinmn1(i)
            z1b(i) = z1b(i) + zc(mn) * cosmn1(i)
         END DO
         DO i = nuv3min, nuv3max
            rub(i) = rub(i) + xm(mn) * rs(mn) * cosmn1(i)
            rvb(i) = rvb(i) - xn(mn) * rs(mn) * cosmn1(i)
            zub(i) = zub(i) - xm(mn) * zc(mn) * sinmn1(i)
            zvb(i) = zvb(i) + xn(mn) * zc(mn) * sinmn1(i)
            ruu(i) = ruu(i) - xm(mn)*xm(mn)*rs(mn) * sinmn1(i)
            ruv(i) = ruv(i) + xm(mn)*xn(mn)*rs(mn) * sinmn1(i)
            rvv(i) = rvv(i) - xn(mn)*xn(mn)*rs(mn) * sinmn1(i)
            zuu(i) = zuu(i) - xm(mn)*xm(mn)*zc(mn) * cosmn1(i)
            zuv(i) = zuv(i) + xm(mn)*xn(mn)*zc(mn) * cosmn1(i)
            zvv(i) = zvv(i) - xn(mn)*xn(mn)*zc(mn) * cosmn1(i)
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
      DO i = nuv3min, nuv3max
        guu_b(i) = rub(i)*rub(i) + zub(i)*zub(i)
        guv_b(i) = (rub(i)*rvb(i)+ zub(i)*zvb(i))*onp*2
        gvv_b(i) = (rvb(i)*rvb(i)+ zvb(i)*zvb(i)+(r1b(i)*r1b(i)))*onp2
        snr(i) = signgs*r1b(i)*zub(i)
        snv(i) = signgs*(rub(i)*zvb(i) - rvb(i)*zub(i))
        snz(i) =-signgs*r1b(i)*rub(i)
        drv(i) = -(r1b(i)*snr(i) + z1b(i)*snz(i))
        auu(i) = p5*(snr(i)*ruu(i) + snz(i)*zuu(i))
        auv(i) = (snr(i)*ruv(i) + snv(i)*rub(i) + snz(i)*zuv(i))*onp
        avv(i) = (snv(i)*rvb(i) + p5*(snr(i)*(rvv(i) - r1b(i))
     1                          +     snz(i)* zvv(i)))*onp2
      END DO

      DO i = 1, nuv3
         rzb2(i) = r1b(i)*r1b(i) + z1b(i)*z1b(i)
      END DO
      IF (.NOT.lasym) THEN
         DO i = 1 + nv, nuv3 - nv
            rzb2(imirr(i)) = rzb2(i)
            r1b(imirr(i))  = r1b(i)
            z1b(imirr(i))  =-z1b(i)
         END DO
      END IF

      DO i = 1, nuv
        rcosuv(i) = r1b(i)*cosuv(i)
        rsinuv(i) = r1b(i)*sinuv(i)
      END DO

      DEALLOCATE (ruu, ruv, rvv, zuu, zuv, zvv, cosmn1, sinmn1, stat=i)

      CALL second0(tsurfoff)
      surface_time = surface_time + (tsurfoff-tsurfon)

      END SUBROUTINE surface

