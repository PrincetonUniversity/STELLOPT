
      subroutine caprsh2(irho)
c--
c  Evaluate CAPR, CAPS, and H2.  Here, m is the poloidal mode number
c  and n is the toroidal mode number/periods.
c--LAB--changed calculation of H2 to use proper 1/b**2 Jacobian for
c  flux surface average
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      use parambs
      use vmec0
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer irho
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nh, m, i, j
      real(rprec) :: h2top_th, h2top_phi, den, qn, e2
      real(rprec) ::  b2, h2top
      real(rprec) :: h2top_tmp, den_tmp, sin_eps
C-----------------------------------------------------------------------

c  the following code was added to evaluate H2 LAB
      h2top = zero
      den   = zero
      do i=1, nthetah
         do j=1, nzetah
            h2top_th = zero
            h2top_phi = zero
            den_tmp = zero
            do nh = 0, nbuse            !  evaluate i,j terms in num and denom of H2
               qn = qsafety(irho)*periods*nh*zetasign
               do m = -mbuse, mbuse
                  sin_eps = amnfit(irho,m,nh)*
     1            (sinmi(m,i)*cosnj(nh,j)+cosmi(m,i)*sinnj(nh,j))
                  h2top_th = h2top_th - m*sin_eps
                  h2top_phi = h2top_phi - qn*sin_eps
                  den_tmp = den_tmp - (m + qn)*sin_eps
               enddo
            enddo
            b2 = bfield(i,j)**2
            h2top = h2top + (h2top_th**2 - h2top_phi**2)/b2
            den  = den + den_tmp**2/b2
         enddo
      enddo

      ! Modified by SAL so code won't stop
      !if (den .eq. zero) stop 'den = 0 in caprsh2'
      IF (den .eq. zero) THEN
         IF (lscreen) THEN
            WRITE(6,'(A)')      'WARNING(BOOTSJ): den = 0 in caprsh2'
            WRITE(6,'(A,I3.3)') '                 irho: ',irho
         END IF
         h2(irho) = 0.0
      ELSE
         h2(irho) = h2top/den
      END IF
      capr(irho) = (one - h2(irho))/(2*qsafety(irho))
      caps(irho) = (one + h2(irho))/2

      end subroutine caprsh2
