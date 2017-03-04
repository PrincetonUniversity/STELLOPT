      SUBROUTINE vcoords_rz(rmnc, zmns, lmns, rmns, zmnc, lmnc, xm, xn,
     1   ntorsum, ns, jrad, mnmax, r, z, lt, lz, lam, sfull, 
     2   nparity, nznt, nfp, lasym)
      USE stel_kinds
      USE booz_persistent, ONLY: cosm_b, sinm_b, cosn_b, sinn_b
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: jrad, ns, mnmax, nparity, nznt, nfp
      INTEGER, DIMENSION(0:1) :: ntorsum
      REAL(rprec), DIMENSION(mnmax,ns) :: rmnc, zmns, lmns
      REAL(rprec), DIMENSION(mnmax,ns) :: rmns, zmnc, lmnc
      REAL(rprec), DIMENSION(mnmax), INTENT(in) :: xm, xn
      REAL(rprec), DIMENSION(nznt), INTENT(out) :: r, z
      REAL(rprec), DIMENSION(nznt), INTENT(out) :: lam, lt, lz
      REAL(rprec), DIMENSION(ns), INTENT(in) :: sfull
      LOGICAL, INTENT(in)        :: lasym
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: js, js1, mn, m, n
      REAL(rprec) :: t1, t2, rc, zs, sgn, rs, zc
      REAL(rprec), DIMENSION(nznt) :: tsin, tcos
C-----------------------------------------------
      js = jrad
      js1= js-1
      IF (js .le. 1) STOP 'js must be > 1!'

      r  = zero
      z  = zero

!
!       Compute Reven, Rodd and Zeven, Zodd in Real Space
!       on full radial grid at js and js1 (will average onto half-mesh later)
!       Lambda is on half grid
!       (even, nparity = 0; odd, nparity = 1)
!
      IF (nparity .eq. 0) THEN
         t1 = one
         t2 = one
         lt  = zero
         lz  = zero
         lam = zero
      ELSE IF (js .gt. 2) THEN
         t1 = one/sfull(js)
         t2 = one/sfull(js1)
      ELSE
         t1 = one/sfull(2)
         t2 = one
         rmnc(1+ntorsum(0):ntorsum(1),1) = 2*rmnc(1+ntorsum(0):
     1      ntorsum(1),2)/sfull(2) - rmnc(1+ntorsum(0):ntorsum(1),3)/
     2      sfull(3)
         zmns(1+ntorsum(0):ntorsum(1),1) = 2*zmns(1+ntorsum(0):
     1      ntorsum(1),2)/sfull(2) - zmns(1+ntorsum(0):ntorsum(1),3)/
     2      sfull(3)
         IF (lasym) THEN
            rmns(1+ntorsum(0):ntorsum(1),1) = 2*rmns(1+ntorsum(0):
     1      ntorsum(1),2)/sfull(2) - rmns(1+ntorsum(0):ntorsum(1),3)/
     2      sfull(3)
            zmnc(1+ntorsum(0):ntorsum(1),1) = 2*zmnc(1+ntorsum(0):
     1      ntorsum(1),2)/sfull(2) - zmnc(1+ntorsum(0):ntorsum(1),3)/
     2      sfull(3)
         ENDIF
      ENDIF

      t1 = t1/2
      t2 = t2/2

      DO mn = 1, mnmax
         m = NINT(xm(mn))
         IF (MOD(m,2) .ne. nparity) CYCLE

         n = NINT(ABS(xn(mn)/nfp))
         sgn = SIGN(one,xn(mn))

         tcos = cosm_b(:,m)*cosn_b(:,n)
     1        + sinm_b(:,m)*sinn_b(:,n)*sgn
         tsin = sinm_b(:,m)*cosn_b(:,n)
     1        - cosm_b(:,m)*sinn_b(:,n)*sgn
         rc = t1*rmnc(mn,js)+t2*rmnc(mn,js1)
         zs = t1*zmns(mn,js)+t2*zmns(mn,js1)
         r   = r   + tcos*rc
         z   = z   + tsin*zs
         lt  = lt  + tcos*lmns(mn,js)*xm(mn)
         lz  = lz  - tcos*lmns(mn,js)*xn(mn)
         lam = lam + tsin*lmns(mn,js)
         IF (lasym) THEN
            rs = t1*rmns(mn,js) + t2*rmns(mn,js)
            zc = t1*zmnc(mn,js) + t2*zmnc(mn,js)
            r   = r   + tsin*rs
            z   = z   + tcos*zc
            lt  = lt  - tsin*lmnc(mn,js)*xm(mn)
            lz  = lz  + tsin*lmnc(mn,js)*xn(mn)
            lam = lam + tcos*lmnc(mn,js)
         END IF
      END DO

      END SUBROUTINE vcoords_rz

      SUBROUTINE vcoords_w(bmnc, bmns, pmns, pmnc, xm, xn, jrad, 
     1    mnmax, bmod, wt, wz, w, nznt, nfp, lasym)
      USE stel_kinds
      USE booz_persistent, ONLY: cosm_w => cosm_nyq, sinm_w => sinm_nyq,
     1                           cosn_w => cosn_nyq, sinn_w => sinn_nyq
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: jrad, mnmax, nznt, nfp
      REAL(rprec), DIMENSION(mnmax), INTENT(in) :: xm, xn, pmns, bmnc
      REAL(rprec), DIMENSION(mnmax), INTENT(in) :: pmnc, bmns
      REAL(rprec), DIMENSION(nznt), INTENT(out) :: w, wt, wz, bmod
      LOGICAL, INTENT(in) :: lasym
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mn, m, n
      REAL(rprec) :: sgn
      REAL(rprec), DIMENSION(nznt) :: tsin, tcos
C-----------------------------------------------
      IF (jrad .le. 1) STOP 'jrad must be > 1!'

!
!     Compute w and derivatives (p transformation of right-side in Eq.(10)) 
!     and |B| on half radial grid in REAL space
!
      w    = zero
      wt   = zero
      wz   = zero
      bmod = zero

      DO mn = 1, mnmax
         m = NINT(xm(mn))
         n = NINT(ABS(xn(mn)))/nfp
         sgn = SIGN(one,xn(mn))

         tcos = cosm_w(:,m)*cosn_w(:,n)
     1        + sinm_w(:,m)*sinn_w(:,n)*sgn
         tsin = sinm_w(:,m)*cosn_w(:,n)
     1        - cosm_w(:,m)*sinn_w(:,n)*sgn
         w   = w   + tsin*pmns(mn)
         wt  = wt  + tcos*pmns(mn)*xm(mn)
         wz  = wz  - tcos*pmns(mn)*xn(mn)
         bmod= bmod+ tcos*bmnc(mn)
         
         IF (.not.lasym) CYCLE
         
         w   = w   + tcos*pmnc(mn)
         wt  = wt  - tsin*pmnc(mn)*xm(mn)
         wz  = wz  + tsin*pmnc(mn)*xn(mn)
         bmod= bmod+ tsin*bmns(mn)

      END DO

      END SUBROUTINE vcoords_w


      SUBROUTINE vcoords_rzb(rmnc, zmns, rmns, zmnc, xmb, xnb,
     1                       cosm_boz, sinm_boz, cosn_boz, sinn_boz,
     2            mboz, nboz,  mnmax, jsurf, ns, r, z, nznt, nfp, lasym)
      USE stel_kinds
      USE booz_persistent, ONLY: thgrd, ztgrd
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: jsurf, ns, mnmax, nznt, nfp, mboz, nboz
      REAL(rprec), DIMENSION(mnmax,ns), INTENT(in) :: rmnc, zmns    !BOOZER COORDINATES
      REAL(rprec), DIMENSION(mnmax,ns), INTENT(in) :: rmns, zmnc
      REAL(rprec), DIMENSION(mnmax), INTENT(in) :: xmb, xnb
      REAL(rprec), DIMENSION(nznt), INTENT(out) :: r, z
      REAL(rprec), DIMENSION(nznt,0:mboz) :: cosm_boz, sinm_boz
      REAL(rprec), DIMENSION(nznt,0:nboz) :: cosn_boz, sinn_boz
      LOGICAL, INTENT(in)        :: lasym
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: js, mn, m, n
      REAL(rprec) :: rc, zs, sgn, rs, zc
      REAL(rprec), DIMENSION(nznt) :: tsin, tcos
C-----------------------------------------------
      js = jsurf

      r = 0
      z = 0

!
!     theta-boz = thgrd:  now uniform angle mesh in BOOZER space
!     zeta-boz  = ztgrd
!
!     IF user wants to check xform in VMEC space, change
!     thgrd -> thgrd + uboz,  ztgrd -> ztgrd in call to trigfunc
!
      CALL trigfunc (thgrd, ztgrd, cosm_boz, sinm_boz, cosn_boz, 
     1               sinn_boz, mboz, nboz, nznt)

!
!       Compute Reven, Rodd and Zeven, Zodd in Real BOOZER Space
!       on half radial grid (rmncb, zmnsb, etc ARE ALREADY on half mesh
!

      DO mn = 1, mnmax
         m = NINT(xmb(mn))
         n = NINT(ABS(xnb(mn)/nfp))
         sgn = SIGN(one,xnb(mn))

         tcos = cosm_boz(:,m)*cosn_boz(:,n)
     1        + sinm_boz(:,m)*sinn_boz(:,n)*sgn
         tsin = sinm_boz(:,m)*cosn_boz(:,n)
     1        - cosm_boz(:,m)*sinn_boz(:,n)*sgn
         rc = rmnc(mn,js)
         zs = zmns(mn,js)
         r  = r + tcos*rc
         z  = z + tsin*zs
         IF (lasym) THEN
            rs = rmns(mn,js)
            zc = zmnc(mn,js)
            r = r + tsin(:nznt)*rs
            z = z + tcos(:nznt)*zc
         END IF
      END DO

      END SUBROUTINE vcoords_rzb
