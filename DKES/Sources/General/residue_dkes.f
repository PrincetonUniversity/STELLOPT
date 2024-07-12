      SUBROUTINE residue_dkes(a, bm1, bp1, cm2, cp2, csave, fz1, fz3,
     1   srces, rsd1, rsd3, g11, g33, g31, g13, crs1, crs3)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE Vnamecl2
      USE dkes_input, ONLY: idisk, lalpha
      USE dkes_realspace, ONLY: diagl, diagle, srces0, mpnt
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), INTENT(out) ::
     1     rsd1, rsd3, g11, g33, g31, g13, crs1, crs3
      REAL(rprec), DIMENSION(mpnt,0:lalpha), INTENT(in) :: fz1, fz3
      REAL(rprec), DIMENSION(mpnt,mpnt) :: a, bm1, bp1, cm2, cp2, csave
      REAL(rprec), INTENT(in), DIMENSION(mpnt,4,2,2) :: srces
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: k, kkount, mnp, j, ll
      REAL(rprec) :: fnrm1, fnrm3, fac, fnx1, fnx3, gnx1
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: src1, src3, xrc1, xrc3
!-----------------------------------------------------------------------
!
!  This SUBROUTINE calculates the solution residuals.
!  See Eq.(24) and following in W. I. van Rij, et. al.
!
!  cm2 * f(l-2) + bm1 * f(l-1) + a * f(l) + bp1 * f(l+1) + cp2 * f(l+2) = source(l)
!
!  WHERE:   bp1(l) = bm1(l+1) (transpose);     cp2(l) = cm1(l+2) (transpose)
!
!
!     READ in f for tri-diagonal code
!
!     READ (33, iostat=k) fz1(1:mpnt,1:), fz3(1:mpnt,1:)
!     READ (33, iostat=k) fz1(1:mpnt,0), fz3(1:mpnt,0)
!     IF (k .ne. 0) STOP 'Unable to READ in data from tri-diag unit33'

      ALLOCATE (src1(mpnt), src3(mpnt), xrc1(mpnt), xrc3(mpnt), stat=k)
      IF (k .ne. 0) STOP 'allocation error in DKES residue!'

      fnrm1 = 0;  fnrm3 = 0;  rsd1 = 0;  rsd3 = 0
      g11 = 0;    g33 = 0;    g31 = 0;   g13 = 0

      k = idisk*lam6 + 1
      ll = lap1 - k
      kkount = lam3 - k

!  l >= 2

      CALL blox (ll, a, bm1, cm2, src1, src3, srces0)
      IF (idisk .eq. 0) THEN                                !!ll = lalpha block row
         fnrm1 = fnrm1 + SUM(src1*src1)
         fnrm3 = fnrm3 + SUM(src3*src3)
         DO mnp = 1, mpnt
            src1(:) = src1(:) - a(:,mnp)*fz1(mnp,ll) - bm1(:,mnp)
     1         *fz1(mnp,ll-1) - cm2(:,mnp)*fz1(mnp,ll-2)
            src3(:) = src3(:) - a(:,mnp)*fz3(mnp,ll) - bm1(:,mnp)
     1         *fz3(mnp,ll-1) - cm2(:,mnp)*fz3(mnp,ll-2)
         END DO
         rsd1 = rsd1 + SUM(src1*src1)
         rsd3 = rsd3 + SUM(src3*src3)
         g11 = g11 - SUM(src1*fz1(:,ll))
         g33 = g33 - SUM(src3*fz3(:,ll))
         g31 = g31 - SUM(src1*fz3(:,ll))
         g13 = g13 - SUM(src3*fz1(:,ll))
      ENDIF

      bp1 = bm1
      cp2 = cm2

      k = k + 1
      ll = ll - 1
      CALL blox (ll, a, bm1, cm2, src1, src3, srces0)
      IF (idisk .eq. 0) THEN
         fnrm1 = fnrm1 + SUM(src1*src1)
         fnrm3 = fnrm3 + SUM(src3*src3)
         DO mnp = 1, mpnt
            src1(:) = src1(:) - bp1(mnp,:)*fz1(mnp,ll+1) -
     1         a(:,mnp)*fz1(mnp,ll) - bm1(:,mnp)*fz1(mnp,ll-1) -
     2         cm2(:,mnp)*fz1(mnp,ll-2)
            src3(:) = src3(:) - bp1(mnp,:)*fz3(mnp,ll+1) -
     1         a(:,mnp)*fz3(mnp,ll) - bm1(:,mnp)*fz3(mnp,ll-1) -
     2         cm2(:,mnp)*fz3(mnp,ll-2)
         END DO
         rsd1 = rsd1 + SUM(src1*src1)
         rsd3 = rsd3 + SUM(src3*src3)
         g11 = g11 - SUM(src1*fz1(:,ll))
         g33 = g33 - SUM(src3*fz3(:,ll))
         g31 = g31 - SUM(src1*fz3(:,ll))
         g13 = g13 - SUM(src3*fz1(:,ll))
      ENDIF

      bp1 = bm1
      csave = cm2

      DO j = 1, kkount
         k = k + 1
         ll = ll -1
         CALL blox (ll, a, bm1, cm2, src1, src3, srces0)
         fnrm1 = fnrm1 + SUM(src1*src1)
         fnrm3 = fnrm3 + SUM(src3*src3)
         DO mnp = 1, mpnt
            src1(:) = src1(:) - cp2(mnp,:)*fz1(mnp,ll+2) -
     1         bp1(mnp,:)*fz1(mnp,ll+1) - a(:,mnp)*fz1(mnp,ll) -
     2         bm1(:,mnp)*fz1(mnp,ll-1) - cm2(:,mnp)*fz1(mnp,ll-2)
            src3(:) = src3(:) - cp2(mnp,:)*fz3(mnp,ll+2) -
     1         bp1(mnp,:)*fz3(mnp,ll+1) - a(:,mnp)*fz3(mnp,ll) -
     2         bm1(:,mnp)*fz3(mnp,ll-1) - cm2(:,mnp)*fz3(mnp,ll-2)
         END DO
         rsd1 = rsd1 + SUM(src1*src1)
         rsd3 = rsd3 + SUM(src3*src3)
         g11 = g11 - SUM(src1*fz1(:,ll))
         g33 = g33 - SUM(src3*fz3(:,ll))
         g31 = g31 - SUM(src1*fz3(:,ll))
         g13 = g13 - SUM(src3*fz1(:,ll))
         bp1 = bm1
         cp2 = csave
         csave = cm2
      END DO

!  l = 1

      k = k + 1
      ll = ll - 1

      CALL blox (ll, a, bm1, cm2, src1, src3, srces0)
      fnrm1 = fnrm1 + SUM(src1*src1)
      fnrm3 = fnrm3 + SUM(src3*src3)

      DO mnp = 1, mpnt
         src1(:) = src1(:) - diagl(:,mnp,iswpm)*fz1(mnp,0) -
     1      cp2(mnp,:)*fz1(mnp,ll+2) -
     2      bp1(mnp,:)*fz1(mnp,ll+1) - a(:,mnp)*fz1(mnp,ll) -
     3      bm1(:,mnp)*fz1(mnp,ll-1)
         src3(:) = src3(:) - diagl(:,mnp,iswpm)*fz3(mnp,0) -
     1      cp2(mnp,:)*fz3(mnp,ll+2) -
     2      bp1(mnp,:)*fz3(mnp,ll+1) - a(:,mnp)*fz3(mnp,ll) -
     3      bm1(:,mnp)*fz3(mnp,ll-1)
      END DO
      rsd1 = rsd1 + SUM(src1*src1)
      rsd3 = rsd3 + SUM(src3*src3)
      g11 = g11 - SUM(src1*fz1(:,ll))
      g33 = g33 - SUM(src3*fz3(:,ll))
      g31 = g31 - SUM(src1*fz3(:,ll))
      g13 = g13 - SUM(src3*fz1(:,ll))
      bp1 = bm1
      cp2 = csave

!  l = 0

      k = k + 1
      ll = ll - 1
      CALL blox (ll, a, bm1, cm2, src1, src3, srces0)
      fnrm1 = fnrm1 + SUM(src1*src1)
      fnrm3 = fnrm3 + SUM(src3*src3)
      fac = iswpm - 1
      fnx1 = 0
      fnx3 = 0

      DO mnp = 1, mpnt
         src1(:) = src1(:) - diagle(:,mnp,iswpm)*fz1(mnp,0)
     1           - cp2(mnp,:)*fz1(mnp,ll+2) - bp1(mnp,:)
     2           * fz1(mnp,ll+1) - a(:,mnp)*fz1(mnp,ll)
         src3(:) = src3(:) - diagle(:,mnp,iswpm)*fz3(mnp,0)
     1           - cp2(mnp,:)*fz3(mnp,ll+2) - bp1(mnp,:)
     2           *fz3(mnp,ll+1) - a(:,mnp)*fz3(mnp,ll)
      END DO


!  Particle conservation. Here, The diagle, diagl coefficients (for V(l=0)) are
!  -transpose of the ones used in the vf(l=0) contributions above.

      xrc1(:) = MATMUL(TRANSPOSE(diagle(:,:,iswpm)),fz1(:,ll))
     1        + MATMUL(TRANSPOSE(diagl(:,:,iswpm)), fz1(:,ll+1))
     2        + fac*srces(:,1,1,1)
      xrc3(:) = MATMUL(TRANSPOSE(diagle(:,:,iswpm)),fz3(:,ll))
     1        + MATMUL(TRANSPOSE(diagl(:,:,iswpm)), fz3(:,ll+1))
      fnx1 = SUM((ABS(MATMUL(TRANSPOSE(diagle(:,:,1)),fz1(:,ll)))
     1        + ABS(MATMUL(TRANSPOSE(diagl(:,:,1)),fz1(:,ll+1))))**2)
      fnx3 = SUM((ABS(MATMUL(TRANSPOSE(diagle(:,:,1)),fz3(:,ll)))
     1        + ABS(MATMUL(TRANSPOSE(diagl(:,:,1)),fz3(:,ll+1))))**2)

      crs1 = SUM(xrc1*xrc1)
      crs3 = SUM(xrc3*xrc3)
      gnx1 = fac*SUM(srces(:,1,1,1)*srces(:,1,1,1))

      IF (gnx1 .ne. zero) fnx1 = gnx1
      IF (fnx1 .eq. zero) fnx1 = one
      IF (fnx3 .eq. zero) fnx3 = one
      IF (fnrm1 .eq. zero) fnrm1 = EPSILON(fnrm1)
      IF (fnrm3 .eq. zero) fnrm3 = EPSILON(fnrm3)
      crs1 = SQRT(crs1/fnx1)
      crs3 = SQRT(crs3/fnx3)
      IF (efield1.eq.zero .and. iswpm.eq.2) crs3 = zero
      rsd1 = SQRT((rsd1 + SUM(src1*src1))/fnrm1)
      rsd3 = SQRT((rsd3 + SUM(src3*src3))/fnrm3)

!     ll = 1 here (SPH060910: replaced fz1,3(1,ll) with fz1,3(:,ll))
      g11 = (g11 - SUM(src1*fz1(:,1)) - SUM(xrc1*fz1(:,0)))/vp
      g33 =  g33 - SUM(src3*fz3(:,1)) - SUM(xrc3*fz3(:,0))/vp
      g31 = (g31 - SUM(src1*fz3(:,1)) - SUM(xrc1*fz3(:,0)))/vp
      g13 = (g13 - SUM(src3*fz1(:,1)) - SUM(xrc3*fz1(:,0)))/vp

      DEALLOCATE (src1, src3, xrc1, xrc3, stat=k)

      END SUBROUTINE residue_dkes
