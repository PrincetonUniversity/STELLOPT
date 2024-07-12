      SUBROUTINE boozer(thgrd, ztgrd, bmod, rad, zee, xmb, xnb,
     1   bmncb, rmncb, zmnsb, pmnsb, gmncb, 
     2   bmnsb, rmnsb, zmncb, pmncb, gmnsb,
     3   scl, uboz, vboz, xjac, 
     4   cosmm, sinmm, cosnn, sinnn, mnmax, nznt, mboz, nboz, nfp, 
     5   nu2, nv, jacfac)
C...MODIFIED 6/98 by A. WARE to speed up by factor of 8
C
      USE stel_kinds
      USE booz_params, ONLY: lasym_b
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: mnmax, nznt, mboz, nboz, nfp, nu2, nv
      REAL(rprec), DIMENSION(nznt), INTENT(in) ::
     1   thgrd, ztgrd, bmod, rad, zee, uboz, vboz, xjac
      REAL(rprec), DIMENSION(mnmax), INTENT(in) ::
     1   xmb, xnb, scl
      REAL(rprec), DIMENSION(nznt,0:mboz), INTENT(out) ::
     1   cosmm, sinmm
      REAL(rprec), DIMENSION(nznt,0:nboz), INTENT(out) ::
     1   cosnn, sinnn
      REAL(rprec), DIMENSION(mnmax), INTENT(out) ::
     1   bmncb, rmncb, zmnsb, pmnsb, gmncb
      REAL(rprec), DIMENSION(mnmax), INTENT(out) ::
     1   bmnsb, rmnsb, zmncb, pmncb, gmnsb
      REAL(rprec), INTENT(in) :: jacfac
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: one=1, zero=0, p5=0.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mn, m, n, imax, i
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: cost, sint, uang, 
     1                                          vang, bbjac
      REAL(rprec) :: sgn
C-----------------------------------------------
!
!     theta-boz = thgrd + uboz              
!     zeta-boz  = ztgrd + vboz
!
      ALLOCATE (cost(nznt), sint(nznt), uang(nznt), 
     1          vang(nznt), bbjac(nznt), stat=i)
      IF (i .ne. 0) STOP 'Allocation error in boozer!'

      uang = thgrd+uboz
      vang = ztgrd+vboz
      CALL trigfunc (uang, vang, cosmm, sinmm, cosnn, sinnn, 
     1               mboz, nboz, nznt)

      IF (.not.lasym_b) THEN
!     ONLY INTEGRATE IN U HALF WAY AROUND (FOR LASYM=F)
         i = nv*(nu2-1)+1                               !u=pi interval: i:imax
         imax = i-1+nv
         DO m = 0,mboz
            cosmm(1:nv,m)    = p5*cosmm(1:nv,m)         !u=0
            cosmm(i:imax,m)  = p5*cosmm(i:imax,m)       !u=pi
            sinmm(1:nv,m)    = p5*sinmm(1:nv,m)         !should be zeroes
            sinmm(i:imax,m)  = p5*sinmm(i:imax,m)       !should be zeroes
         END DO
      END IF

!     jacobian from VMEC to Boozer coords, with SPECIAL
!     radial variable s = (toroidal flux)/twopi (phip = 1)
!     cost = cos(mu-nv);  sint = sin(mu-nv)
      bbjac = jacfac/(bmod*bmod)

      DO mn = 1,mnmax
        m = NINT(xmb(mn))
        n = NINT(ABS(xnb(mn)))/nfp
        sgn = SIGN(one,xnb(mn))
        cost = (cosmm(:,m)*cosnn(:,n)
     1       +  sinmm(:,m)*sinnn(:,n)*sgn)*xjac
        sint = (sinmm(:,m)*cosnn(:,n)
     1       -  cosmm(:,m)*sinnn(:,n)*sgn)*xjac

        bmncb(mn) = DOT_PRODUCT(bmod,cost)
        rmncb(mn) = DOT_PRODUCT(rad, cost)
        zmnsb(mn) = DOT_PRODUCT(zee, sint)
        pmnsb(mn) =-DOT_PRODUCT(vboz,sint)
        gmncb(mn) = DOT_PRODUCT(bbjac, cost)

        IF (.not.lasym_b) CYCLE

        bmnsb(mn) = DOT_PRODUCT(bmod,sint)
        rmnsb(mn) = DOT_PRODUCT(rad ,sint)
        zmncb(mn) = DOT_PRODUCT(zee, cost)
        pmncb(mn) =-DOT_PRODUCT(vboz,cost)
        gmnsb(mn) = DOT_PRODUCT(bbjac, sint)

      END DO

      DEALLOCATE (cost, sint, uang, vang, bbjac, stat=i)

      bmncb = scl*bmncb
      rmncb = scl*rmncb
      zmnsb = scl*zmnsb
      pmnsb = scl*pmnsb
      gmncb = scl*gmncb

      IF (.not.lasym_b) RETURN

      bmnsb = scl*bmnsb
      rmnsb = scl*rmnsb
      zmncb = scl*zmncb
      pmncb = scl*pmncb
      gmnsb = scl*gmnsb

      END SUBROUTINE boozer
