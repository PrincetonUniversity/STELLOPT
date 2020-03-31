
      subroutine grad(gradbs1, gradbs2, gradbs3, gradbs4, irho)
c
c     calculate gradient factors, gradbs1 and gradbs2
c     gradbs1 - due to electron density gradient
c     gradbs2 - due to ion density gradient
c     gradbs3 - due to electron temperature gradient
c     gradbs4 - due to ion temperature gradient
c
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      use parambs
      use vmec0
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: irho
      real(rprec) :: gradbs1, gradbs2, gradbs3, gradbs4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec) :: gradbs, de, de1, di, di1, p
C-----------------------------------------------

      gradbs =-2.5_dp*gpsi(irho)*qsafety(irho)*betar(irho)*sign_jacobian
      de = dense(irho)                               !electron density
      de1 = densrho1
      di = dense(irho)/zeff1                         !ion density
      di1 = densrho1/zeff1
      p = de*tempe1(irho) + di*tempi1(irho) + 1.E-36_dp     !plasma pressure
      gradbs1 = gradbs*de1*tempe1(irho)/p  !due to electron density gradient
      gradbs2 = gradbs*di1*tempi1(irho)/p  !due to ion density gradient
c                                     !due to electron temperature gradi
      gradbs3 = alphae*gradbs*temperho1*de/p
c                                       !due to ion temperature gradient
      gradbs4 = alphai*gradbs*tempirho1*di/p

      end subroutine grad
