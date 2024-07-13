      SUBROUTINE coil_curve (nu, uc, us, nv, vc, vs, la, lb, lc,
     1   ls, nf, ns, mnum, nnum, rmn, zmn, nfil, s, alf, bet,
     2   deln, delt, u, v, r, c, x, y, z, ierr)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE stel_constants
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y  A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(out) :: ierr
      INTEGER, INTENT(in)  :: nu, nv, nf, ns, nfil
      REAL(rprec), INTENT(in)  :: s, alf, bet, deln, delt
      REAL(rprec), INTENT(out) :: r
      REAL(rprec), INTENT(out) :: u, v
      REAL(rprec), DIMENSION(0:*), INTENT(inout) :: uc, us
      REAL(rprec), DIMENSION(0:*), INTENT(inout) :: vc, vs
      INTEGER, DIMENSION(*), INTENT(in) :: nnum, mnum
      REAL(rprec), DIMENSION(*), INTENT(in) :: rmn, zmn
      REAL(rprec), DIMENSION(2,3,5), INTENT(out) :: x, y, z
      REAL(rprec), DIMENSION(2), INTENT(out) :: c
      LOGICAL, INTENT(in) :: la, lb, lc, ls
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, k, n, np, mk, nk, ifail
      REAL(rprec), DIMENSION(2) :: x0, x1, x2
      REAL(rprec), DIMENSION(2) :: y0, y1, y2
      REAL(rprec), DIMENSION(2) :: tcx, tcy, tcz
      REAL(rprec), DIMENSION(2) :: tsx, tsy, tsz
      REAL(rprec), DIMENSION(2) :: xu, xv, yu, yv
      REAL(rprec), DIMENSION(2) :: nx, ny, nz, nr, nphi, dn, d1
      REAL(rprec), DIMENSION(2) :: ph0,  ph1,  ph2, phv
      REAL(rprec), DIMENSION(4) :: fu, fv
      REAL(rprec) :: s0, s1, dph
      REAL(rprec) :: ck, sk
      REAL(rprec) :: r0, r1, r2
      REAL(rprec) :: z0, z1, z2
      REAL(rprec) :: u0,  u1,  u2
      REAL(rprec) :: v0,  v1,  v2
      REAL(rprec) :: ru, rv, zu, zv
      REAL(rprec) :: sinph0, cosph0
!-----------------------------------------------

      dph = twopi/nf
      s0 = twopi*s
      s1 = twopi

!     compute toroidal variable v and derivatives wrt s
      IF (ls) THEN

         ! b-spline series representation
         IF (lc) THEN
            v0 = 0
         ELSE
            v0 = vc(0)
         END IF
         v1 = 0
         v2 = 0
         ! coefficients
         IF (nv .lt. 7) THEN
            ierr = 1
         ELSE
            ! breakpoint boundary conditions
            DO k = 1, 4
               vs(k) = zero
               vs(nv + k) = one
            END DO
            ! control point boundary conditions
!           IF (lc) THEN
!              vc(nv) = vc(1)
!              vc(nv-1) = 2.0_dp*vc(1) - vc(2)
!              vc(nv-2) = 2.0_dp*vc(1) - vc(3)
!           END IF
            CALL speval (nv, vs, vc, s, fv, ifail)
            v0 = v0 + fv(1)
            v1 = fv(2)
            v2 = fv(3)
         END IF

      ELSE

         ! fourier series representation
         v0 = 0
         v1 = 0
         v2 = 0
         DO k = 0,nv
            ck = COS(k*s0)
            sk = SIN(k*s0)
            v0 = v0 +  vc(k)*ck + vs(k)*sk
            v1 = v1 + (vs(k)*ck - vc(k)*sk)*k
            v2 = v2 + (vs(k)*sk + vc(k)*ck)*k**2
         END DO
         v1 = s1*v1
         v2 = -s1**2*v2

      END IF           ! IF (ls)

      IF (lb .and. (.not. lc)) THEN
!     lb = T => helical or wavy_pf
         v0 = v0 + s
         v1 = v1 + 1
      END IF

!     compute poloidal variable u and derivatives wrt s
      IF (ls) THEN

         ! b-spline series representation
         u0 = 0
         u1 = 0
         u2 = 0
         ! coefficients
         IF (nu .lt. 7) THEN
            ierr = 1
         ELSE
            ! breakpoint boundary conditions
            DO k = 1, 4
               us(k) = zero
               us(nu + k) = one
            END DO
            ! control point boundary conditions
!           IF (lc) THEN
!              uc(nu) = one + uc(1)
!              uc(nu-1) = uc(nu) - uc(2) + uc(1)
!              uc(nu-2) = uc(nu) - uc(3) + uc(1)
!           END IF
            CALL speval (nu, us, uc, s, fu, ifail)
            u0 = fu(1)
            u1 = fu(2)
            u2 = fu(3)
         END IF

      ELSE

         ! fourier series representation
         u0 = 0
         u1 = 0
         u2 = 0
         DO k = 0,nu
            ck = COS(k*s0)
            sk = SIN(k*s0)
            u0 = u0 +  uc(k)*ck + us(k)*sk
            u1 = u1 + (us(k)*ck - uc(k)*sk)*k
            u2 = u2 + (us(k)*sk + uc(k)*ck)*k**2
         END DO

         u1 = s1*u1
         u2 = -s1**2*u2

      END IF           ! IF (ls)

      IF (la .and. (.not. lc)) THEN
!        la = T => modular
         u0 = u0 + s
         u1 = u1 + 1
      END IF

!     compute R, Z and derivatives wrt s

      r0 = 0
      r1 = 0
      r2 = 0
      z0 = 0
      z1 = 0
      z2 = 0
      ru = 0
      rv = 0
      zu = 0
      zv = 0
      DO k = 1, ns
         mk = mnum(k)
         nk = nnum(k)
         ck = COS(twopi*(mk*u0 + nk*v0))    !!Half time used by these 2 lines on PC
         sk = SIN(twopi*(mk*u0 + nk*v0))
!     R ...
         r0 = r0 + rmn(k)*ck
         r1 = r1 - rmn(k)*(twopi*(mk*u1 + nk*v1))*sk
         r2 = r2 - rmn(k)*((twopi*(mk*u2 + nk*v2))*sk
     1      + (twopi*(mk*u1 + nk*v1))**2*ck)
         ru = ru - rmn(k)*mk*sk*twopi
         rv = rv - rmn(k)*nk*sk*twopi
!     Z ...
         z0 = z0 + zmn(k)*sk
         z1 = z1 + zmn(k)*(twopi*(mk*u1 + nk*v1))*ck
         z2 = z2 + zmn(k)*((twopi*(mk*u2 + nk*v2))*ck
     1      - (twopi*(mk*u1 + nk*v1))**2*sk)
         zu = zu + zmn(k)*mk*ck*twopi
         zv = zv + zmn(k)*nk*ck*twopi
      END DO

!     derivatives of toroidal angle wrt s
      ph1(1) =  twopi*v1/nf
      ph1(2) = -twopi*v1/nf
      ph2(1) =  twopi*v2/nf
      ph2(2) = -twopi*v2/nf
      phv(1) =  twopi/nf
      phv(2) = -twopi/nf

!     field period invariance
      DO np = 1,nf

!        toroidal angle in field period np
         ph0(1) =  twopi*v0/nf + (np - 1)*dph
         ph0(2) = -twopi*v0/nf + (np - 1)*dph

!        stellarator symmetry
         DO i = 1,2
!        X and derivatives wrt s
         sinph0 = SIN(ph0(i))
         cosph0 = COS(ph0(i))
         x0(i) = r0*cosph0
         x1(i) = r1*cosph0 - r0*ph1(i)*sinph0
         x2(i) = r2*cosph0 - r1*ph1(i)*sinph0
     1         - (r0*ph1(i)**2*cosph0
     2         + (r0*ph2(i) + r1*ph1(i))*sinph0)

!        X derivatives wrt u, v

         xu(i) = ru*cosph0
         xv(i) = rv*cosph0 - r0*sinph0*phv(i)

!        Y and derivatives wrt s

         y0(i) = r0*sinph0
         y1(i) = r1*sinph0 + r0*ph1(i)*cosph0
         y2(i) = r2*sinph0 + r1*ph1(i)*cosph0
     1         + (-r0*ph1(i)**2*sinph0
     2         + (r0*ph2(i) + r1*ph1(i))*cosph0)

!        Y derivatives wrt u, v

         yu(i) = ru*sinph0
         yv(i) = rv*sinph0 + r0*cosph0*phv(i)

!        compute |dr/ds| = dl/ds

         d1(i) = SQRT (x1(i)*x1(i) + y1(i)*y1(i) + z1*z1)

!        END loop for stellarator symmetry
         END DO

!        Tangent vector (to curve) tcx=dx/ds, tcy=dy/ds, tcz=dz/ds

         tcx(1) =  x1(1)/d1(1)
         tcy(1) =  y1(1)/d1(1)
         tcz(1) =  z1/d1(1)
!        stellarator symmetry
         tcx(2) =  x1(2)/d1(2)
         tcy(2) =  y1(2)/d1(2)
         tcz(2) = -z1/d1(2)

!        components of surface normal vector (cylindrical) nr, nphi, nz

         nr(1) =  r0*zu
         nphi(1) =  ru*zv - rv*zu
         nz(1) = -r0*ru
         dn(1) =  SQRT (nr(1)*nr(1) + nphi(1)*nphi(1) + nz(1)*nz(1))
         nr(1) =  nr(1)/dn(1)
         nphi(1) = nphi(1)/dn(1)
         nz(1) =  nz(1)/dn(1)
!        stellarator symmetry
         nr(2) = -r0*zu
         nphi(2) = -ru*zv + rv*zu
         nz(2) = -r0*ru
         dn(2) =  SQRT (nr(2)*nr(2) + nphi(2)*nphi(2) + nz(2)*nz(2))
         nr(2) =  nr(2)/dn(2)
         nphi(2) = nphi(2)/dn(2)
         nz(2) =  nz(2)/dn(2)

!        components of surface normal vector (cartesian) nx, ny, nz

         nx(1) =  yu(1)*zv - zu*yv(1)
         ny(1) =  zu*xv(1) - zv*xu(1)
         nz(1) =  xu(1)*yv(1) - xv(1)*yu(1)
         dn(1) =  SQRT (nx(1)*nx(1) + ny(1)*ny(1) + nz(1)*nz(1))
         nx(1) =  nx(1)/dn(1)
         ny(1) =  ny(1)/dn(1)
         nz(1) =  nz(1)/dn(1)
!        stellarator symmetry
         nx(2) = -yu(2)*zv + zu*yv(2)
         ny(2) = -zu*xv(2) + zv*xu(2)
         nz(2) =  xu(2)*yv(2) - xv(2)*yu(2)
         dn(2) =  SQRT (nx(2)*nx(2) + ny(2)*ny(2) + nz(2)*nz(2))
         nx(2) =  nx(2)/dn(2)
         ny(2) =  ny(2)/dn(2)
         nz(2) =  nz(2)/dn(2)

!        TANgent vector (to surface) tsx, tsy, tsz

         tsx(1) = ny(1)*tcz(1) - tcy(1)*nz(1)
         tsy(1) = nz(1)*tcx(1) - tcz(1)*nx(1)
         tsz(1) = nx(1)*tcy(1) - tcx(1)*ny(1)
         dn(1) = SQRT (tsx(1)*tsx(1) + tsy(1)*tsy(1) + tsz(1)*tsz(1))
         tsx(1) = tsx(1)/dn(1)
         tsy(1) = tsy(1)/dn(1)
         tsz(1) = tsz(1)/dn(1)
!        stellarator symmetry
         tsx(2) = ny(2)*tcz(2) - tcy(2)*nz(2)
         tsy(2) = nz(2)*tcx(2) - tcz(2)*nx(2)
         tsz(2) = nx(2)*tcy(2) - tcx(2)*ny(2)
         dn(2) = SQRT (tsx(2)*tsx(2) + tsy(2)*tsy(2) + tsz(2)*tsz(2))
         tsx(2) = tsx(2)/dn(2)
         tsy(2) = tsy(2)/dn(2)
         tsz(2) = tsz(2)/dn(2)

!        RETURN curvature c = |dr/ds X d2r/ds2|/|dr/ds|^3

         c(1) = SQRT ((y1(1)*z2 - y2(1)*z1)**2
     1        + (z1*x2(1) - z2*x1(1))**2
     2        + (x1(1)*y2(1) - x2(1)*y1(1))**2)/d1(1)**3
!        stellarator symmetry
         c(2) = SQRT ((-y1(2)*z2 + y2(2)*z1)**2
     1        + (-z1*x2(2) + z2*x1(2))**2
     2        + (x1(2)*y2(2) - x2(2)*y1(2))**2)/d1(2)**3

!        RETURN multifilament representation for x, y, z
!        filament 1
         n = 1
         x(1,np,n) =  x0(1)
         y(1,np,n) =  y0(1)
         z(1,np,n) =  z0
!        stellarator symmetry
         x(2,np,n) =  x0(2)
         y(2,np,n) =  y0(2)
         z(2,np,n) = -z0
         IF (nfil .lt. 4) THEN
!           Three filament model
!           filament 2
            n = 2
            x(1,np,n) =  x0(1) + deln*(alf*nx(1) + bet*tsx(1))
            y(1,np,n) =  y0(1) + deln*(alf*ny(1) + bet*tsy(1))
            z(1,np,n) =  z0    + deln*(alf*nz(1) + bet*tsz(1))
!           stellarator symmetry
            x(2,np,n) =  x0(2) + deln*(alf*nx(2) + bet*tsx(2))
            y(2,np,n) =  y0(2) + deln*(alf*ny(2) + bet*tsy(2))
            z(2,np,n) = -z0    + deln*(alf*nz(2) + bet*tsz(2))
!           filament 3
            n = 3
            x(1,np,n) =  x0(1) - deln*(alf*nx(1) + bet*tsx(1))
            y(1,np,n) =  y0(1) - deln*(alf*ny(1) + bet*tsy(1))
            z(1,np,n) =  z0    - deln*(alf*nz(1) + bet*tsz(1))
!           stellarator symmetry
            x(2,np,n) =  x0(2) - deln*(alf*nx(2) + bet*tsx(2))
            y(2,np,n) =  y0(2) - deln*(alf*ny(2) + bet*tsy(2))
            z(2,np,n) = -z0    - deln*(alf*nz(2) + bet*tsz(2))
         ELSE
!           Four filament model (no current in central filament 1)
!           filament 2
            n = 2
            x(1,np,n) =  x0(1) + deln*nx(1) + delt*tsx(1)
            y(1,np,n) =  y0(1) + deln*ny(1) + delt*tsy(1)
            z(1,np,n) =  z0    + deln*nz(1) + delt*tsz(1)
!           stellarator symmetry
            x(2,np,n) =  x0(2) + deln*nx(2) + delt*tsx(2)
            y(2,np,n) =  y0(2) + deln*ny(2) + delt*tsy(2)
            z(2,np,n) = -z0    + deln*nz(2) + delt*tsz(2)
!           filament 3
            n = 3
            x(1,np,n) =  x0(1) - deln*nx(1) + delt*tsx(1)
            y(1,np,n) =  y0(1) - deln*ny(1) + delt*tsy(1)
            z(1,np,n) =  z0    - deln*nz(1) + delt*tsz(1)
!           stellarator symmetry
            x(2,np,n) =  x0(2) - deln*nx(2) + delt*tsx(2)
            y(2,np,n) =  y0(2) - deln*ny(2) + delt*tsy(2)
            z(2,np,n) = -z0    - deln*nz(2) + delt*tsz(2)
!           filament 4
            n = 4
            x(1,np,n) =  x0(1) - deln*nx(1) - delt*tsx(1)
            y(1,np,n) =  y0(1) - deln*ny(1) - delt*tsy(1)
            z(1,np,n) =  z0    - deln*nz(1) - delt*tsz(1)
!           stellarator symmetry
            x(2,np,n) =  x0(2) - deln*nx(2) - delt*tsx(2)
            y(2,np,n) =  y0(2) - deln*ny(2) - delt*tsy(2)
            z(2,np,n) = -z0    - deln*nz(2) - delt*tsz(2)
!           filament 5
            n = 5
            x(1,np,n) =  x0(1) + deln*nx(1) - delt*tsx(1)
            y(1,np,n) =  y0(1) + deln*ny(1) - delt*tsy(1)
            z(1,np,n) =  z0    + deln*nz(1) - delt*tsz(1)
!           stellarator symmetry
            x(2,np,n) =  x0(2) + deln*nx(2) - delt*tsx(2)
            y(2,np,n) =  y0(2) + deln*ny(2) - delt*tsy(2)
            z(2,np,n) = -z0    + deln*nz(2) - delt*tsz(2)
         END IF      ! number of filaments (nfil)

!     END loop for field period invariance
      END DO

!     RETURN u, v, r
      u = u0
      v = v0
      r = r0

      END SUBROUTINE coil_curve
