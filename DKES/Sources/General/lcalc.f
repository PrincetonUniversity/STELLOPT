      SUBROUTINE lcalc
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE Vnamecl2
      USE dkes_realspace, ONLY: diagl
      USE dkes_input, ONLY: lalpha
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: l1 = 2
      REAL(rprec), PARAMETER :: half = 0.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: l, ll
      REAL(rprec) :: ap, am
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: cinv
C-----------------------------------------------
c  This SUBROUTINE calculates l (degree of pitch-angle Legendre
c  polynomial) arrays that are stored once and reused each cmul,efield run.
      lap1 = lalpha + 1
      lam1 = lalpha - 1
      lam3 = lalpha - 3
      lam6 = lalpha - 6

      ALLOCATE (omgl(0:lap1), cols(lalpha), al1(lalpha), al2(lalpha),
     1   al3(lalpha), al4(lalpha), bl1(lalpha), bl2(lalpha),
     2   bl3(lalpha), bl4(lalpha), cl1(lalpha), cl2(lalpha),
     3   cl3(lalpha), cl4(lalpha), cols0(lap1), al01(lalpha),
     4   al02(lalpha), al03(lalpha), al04(lalpha), bl01(lalpha),
     5   bl02(lalpha), bl03(lalpha), bl04(lalpha), cl01(lalpha),
     6   cl02(lalpha), cl03(lalpha), cl04(lalpha), cinv(0:lap1),
     9   stat=l)
      IF (l .ne. 0) STOP 'allocation error in LCALC'

      omgl(0) = 0;   omgl(1) = 0;   cols0(1) = 0
      cinv(0) = 0;   cinv(1) = 0
      DO l = 2, lap1
         ll = l - 1
         omgl(l)  = half*ll/SQRT(4*ll*ll - one)
         cols0(l) = half*ll*l                                            !.5*l*(l+1)
         cinv(l)  = one/cols0(l)                                         !1/nu(l), l>0
      END DO

      DO l = 1, lalpha
         ll = l - 1
         ap = cinv(l+1)*omgl(l+1)**2
         am = cinv(l-1)*omgl(l)**2
!
!        coefficients of BMAT1 - BMAT6 comprising A(l) [multiplies f(l)]
!
         al01(l) = 4*(ap + am)
         al02(l) = ap*ll**2 + am*l**2
         al03(l) = cinv(l)
         al04(l) = 2*(am*l - ap*ll)
         bl01(l) = 2*cinv(l-1)*omgl(l)
         bl02(l) = 2*cinv(l)*omgl(l)
         bl03(l) = cinv(l-1)*l*omgl(l)
         bl04(l) =-cinv(l)*(ll - 1)*omgl(l)
         cl01(l) = cinv(l-1)*omgl(l-1)*omgl(l)                           !w(l) * w(l-1) /-nu(l-1)
!
!        coefficients of BMAT1 - BMAT6 comprising C-(l) [multiplies f(l-2)]
!
         cl02(l) =-(ll - 2)*l*cl01(l)
         cl03(l) = 2*l*cl01(l)
         cl04(l) =-2*(ll - 2)*cl01(l)

         cl01(l) = 4*cl01(l)

      END DO

      DEALLOCATE (cinv)

c   Matrix elements for particle conservation. These are the coefficients
c   of Vf(l=0) (l=0,1 contributions), and are related by -TRANSPOSE to V(l=0) f
c   RecALL that -efield (Bsubv d/du - Bsubu d/dv) is the electric drift term

      diagl(:,:,1)  =-2*omgl(l1) * diagl(:,:,1)
      diagl(:,:,2)  = 2*omgl(l1) * diagl(:,:,2)

      END SUBROUTINE lcalc
