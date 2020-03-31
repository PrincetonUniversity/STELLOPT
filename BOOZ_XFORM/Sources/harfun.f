      SUBROUTINE harfun(jacfac, hiota, gpsi, ipsi, js, nznt, xlt,
     1   xlz, xl, wt, wz, w, uboz, vboz, xjac)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: js, nznt
      REAL(rprec) :: jacfac
      REAL(rprec), DIMENSION(*) :: hiota, gpsi, ipsi
      REAL(rprec), DIMENSION(nznt) ::
     1   xl, xlt, xlz, w, wt, wz, uboz, vboz, xjac
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: one = 1, zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), ALLOCATABLE, DIMENSION(:) ::
     1   bsupu, bsupv, psubu, psubv
      REAL(rprec) :: dem, dem2, gpsi1, hiota1, ipsi1
	INTEGER     :: istat
C-----------------------------------------------
!
!     HERE, W IS THE PART OF THE TRANSFORMATION FUNCTION P
!     WHICH DEPENDS ON THE SURFACE-AVERAGED COVARIANT B COMPONENTS [RIGHT-SIDE IN Eq(10)],
!     COMPUTED IN transpmn AND vcoord_w.
!
!     THE FULL TRANSFORMATION FUNCTION IS THEN [Eq(10 in paper], with D=gpsi+iota*Ipsi:
!
!     P = (w - Ipsi*Lambda) / D
!
!     ALSO [Eq.(3) in paper]:
!
!     uboz = lambda + iota*P   (non-secular piece of boozer theta in VMEC coordinates)
!
!          = (gpsi*lambda + iota*w)/D
!
!     vboz = P                 (non-secular piece of boozer phi in VMEC coordinates)
!
!     FINALLY, XJAC IS THE JACOBIAN BETWEEN BOOZER, VMEC COORDINATES
!
!
!     NOTE THAT LAMBDA = xl, d(lambda)/du = xlt, d(lambda)/dv = xlz
!
      ALLOCATE(bsupu(nznt), bsupv(nznt), psubu(nznt), psubv(nznt),
     1         stat=istat)
	IF (istat .ne. 0) STOP 'Allocation error in harfun!'

      jacfac = gpsi(js) + hiota(js)*ipsi(js)
      IF (jacfac .eq. zero) 
     1   STOP 'Boozer coordinate XFORM failed, jacfac = 0!'
      dem = one/jacfac
      gpsi1 = gpsi(js)*dem
      hiota1 = hiota(js)*dem
      ipsi1 = ipsi(js)*dem

      vboz = dem*w - ipsi1*xl        !TOTAL p in Eq.(10)
      uboz = xl + hiota(js)*vboz
      psubv = dem*wz - ipsi1*xlz
      psubu = dem*wt - ipsi1*xlt
      bsupv = 1 + xlt
      bsupu = hiota(js) - xlz

!
!     Eq. (12) 
!
      xjac = bsupv*(1+psubv) + bsupu*psubu

      dem = MINVAL(xjac)
      dem2= MAXVAL(xjac)
!     SAL 07/06/16
!      IF (dem*dem2 .le. zero) PRINT *,
!     1   ' Jacobian xjac changed sign in harfun in xbooz_xform'

      DEALLOCATE(bsupu, bsupv, psubu, psubv, stat=istat)

      END SUBROUTINE harfun

      SUBROUTINE modbooz(bmnc, bmns, bmod, xmb, xnb,
     1   u_b, v_b, cosmm, sinmm, cosnn, sinnn,
     2   mnmax, mboz, nboz, nfp, lasym)
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: mnmax, mboz, nboz, nfp
      REAL(rprec), DIMENSION(mnmax), INTENT(in) ::
     1   bmnc, bmns
      REAL(rprec), DIMENSION(mnmax), INTENT(in) :: xmb, xnb
      REAL(rprec), DIMENSION(4), INTENT(in)  :: u_b, v_b
      REAL(rprec), DIMENSION(4), INTENT(out) :: bmod
      REAL(rprec), DIMENSION(0:mboz) :: cosmm, sinmm
      REAL(rprec), DIMENSION(0:nboz) :: cosnn, sinnn
      LOGICAL, INTENT(in) :: lasym
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: one = 1.0_dp, zero = 0.0_dp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: mn, m, n, angles
      REAL(rprec) :: cost, sint
      REAL(rprec) :: sgn
!-----------------------------------------------
      bmod = zero

      ANGLE_LOOP: DO angles=1,4

      cosmm(0) = one
      sinmm(0) = zero
      cosmm(1) = COS(u_b(angles))
      sinmm(1) = SIN(u_b(angles))

      cosnn(0) = one
      sinnn(0) = zero
      IF (nboz .ge. 1) THEN
      cosnn(1) = COS(v_b(angles)*nfp)
      sinnn(1) = SIN(v_b(angles)*nfp)
      END IF

      DO m = 2,mboz
         cosmm(m) = cosmm(m-1)*cosmm(1)
     1            - sinmm(m-1)*sinmm(1)
         sinmm(m) = sinmm(m-1)*cosmm(1)
     1            + cosmm(m-1)*sinmm(1)
      END DO

      DO n = 2,nboz
         cosnn(n) = cosnn(n-1)*cosnn(1)
     1            - sinnn(n-1)*sinnn(1)
         sinnn(n) = sinnn(n-1)*cosnn(1)
     1            + cosnn(n-1)*sinnn(1)
      END DO

      DO mn=1,mnmax
        m = NINT(xmb(mn))
        n = NINT(ABS(xnb(mn)))/nfp
        sgn = SIGN(one,xnb(mn))
        cost = cosmm(m)*cosnn(n)
     1       + sinmm(m)*sinnn(n)*sgn
        bmod(angles) = bmod(angles) + bmnc(mn)*cost
        IF (lasym) THEN
        sint = sinmm(m)*cosnn(n)
     1       - cosmm(m)*sinnn(n)*sgn
        bmod(angles) = bmod(angles) + bmns(mn)*sint
        END IF
      END DO

      END DO ANGLE_LOOP

      END SUBROUTINE modbooz
