!-----------------------------------------------------------------------
!     Subroutine:    vmec_ohkawa
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          07/31/2012
!     Description:   This subroutine calculates the Ohkawa correction
!                    factor using the Hirshman-Sigmar formulation.
!                    	K.C. Shaing, B.A. Carreras, N. Dominguez, 
!                           V.E. Lynch, J.S. Tolliver
!                           "Bootstrap current control in stellarators", 
!                           Phys. Fluids B1, 1663 (1989). 
!                           aip.scitation.org/doi/10.1063/1.858945
!                    The Total correction is
!                    j_nbcd = j_f * [1 - Z_f/Zeff*(1-l31)]
!                    we return (1-l31)/Zeff
!-----------------------------------------------------------------------
      SUBROUTINE vmec_ohkawa(sflx,zeff,Gfac)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE read_wout_mod, ONLY: ns, mnmax_nyq, bmnc, gmnc, lasym, &
                               bmns, gmns, xm_nyq, xn_nyq, &
                               lwout_opened, nfp
!-----------------------------------------------------------------------
!     Arguments
!          sflx			Normalized Toroidal Flux
!          zeff         Effective plasma charge
!          Gfac         Correction Factor (1-l31)/zeff
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: sflx, zeff
      REAL(rprec), INTENT(out) :: Gfac
!-----------------------------------------------------------------------
!     PARAMETERS
!          nu			Number of poloidal gridpoints
!          nv           Number of toroidal gridpoints
!          nla          Number of lambda gridpoints
!          pi2          2*pi
!-----------------------------------------------------------------------
      INTEGER, PARAMETER :: nu  = 64
      INTEGER, PARAMETER :: nv  = 64
      INTEGER, PARAMETER :: nla = 64
      REAL(rprec), PARAMETER      :: pi2 = 6.283185482025146D+00
!-----------------------------------------------------------------------
!     Variables
!          ns1			Radial index of first grid point
!          ns2          Radial index of second grid point
!          mn           Linear harmonic index
!          m/n          Poloidal/Toroidal Harmonic
!          zt           Linear theta/zeta index
!          hs           Radial grid spacing
!          zt           Linear array of theta/zeta
!          wlo          Interpolation Weights (whi)
!          wlo_odd      Interpolation Weights, odd (whi_odd)
!          bmn          Fourier Harmonics of |B| (cos)
!          gmn          Fourier Haromincs of Jacobian sqrt(g) (cos)
!          bmn2         Fourier Harmonics of |B| (sin)
!          gmn2         Fourier Haromincs of Jacobian sqrt(g) (sin)
!          u/v          Poloidal/Toroidal var [0,1]
!          kernel       Kernel (m*theta-n*zeta)
!          b            Linear array of |B|
!          g            Linear array of Jacobian sqrt(g)
!          lam          Lambda Array (mu*B_max/B<1, passing)
!          sumg         Sum over Jacobian
!          b2avg        <B^2>
!          bmax         B_max
!          avgbobm2     <B^2>/(B_max*B_max)
!          avgbpov      <sqrt(1-lam*B/B_max)>
!          ft           Trapped Particle fraction
!          fp           Passing particle fraction
!          x            ft/fp
!          z2           Zeff*Zeff
!          d            Demoninator in l31 equation
!          a            Numerator in l31 equation
!-----------------------------------------------------------------------
      INTEGER :: ns1, ns2, zt, mn, m, n
      REAL(rprec) :: hs, wlo, whi, wlo_odd, whi_odd, sumg, b2avg, &
                     bmax, fp ,ft, x, z2, d, a, avgbobm2, kernel, u, v
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: bmn, gmn, b, g, lam, &
                                                avgbpov, bmn2, gmn2
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: tz
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Intializations
      Gfac = 0
      ns1 = 0; ns2 = 0;
      z2 = zeff*zeff

      ! Checks
      IF (.not.lwout_opened) &
         STOP "ERROR: READ_WOUT not loaded!"
      IF (sflx>1 .or. sflx<0) RETURN

      ! Interpolate or extrapolate
      hs  = 1.0 / (ns - 1)
      ns1  = INT(1.5 + sflx * (ns - 1))
      ns2  = ns1 + 1
      wlo = (hs * (ns2 - 1.5) - sflx) * (ns - 1)
      whi = 1 - wlo
      IF (ns1 == ns) THEN
         ns2 = ns1-1
         wlo = 1+whi
         whi = -whi
      ELSE IF (ns1 == 1) THEN
         ns1 = 2
      END IF
      whi_odd = whi * SQRT(sflx / (hs * (ns2 - 1.5)))
      IF (ns1 /= 1) THEN
         wlo_odd = wlo * SQRT( sflx / (hs * (ns1 - 1.5)))
      ELSE
         wlo_odd = 0
         whi_odd = SQRT(sflx / (hs * (ns2 - 1.5)))
      END IF

      ! Extrapolate Harmonics
      ALLOCATE(bmn(mnmax_nyq),gmn(mnmax_nyq))
      gmn = 0
      bmn = 0
      WHERE (MOD(NINT(xm_nyq(:)),2) == 0)
         bmn = wlo * bmnc(:,ns1) + whi * bmnc(:,ns2)
         gmn = wlo * gmnc(:,ns1) + whi * gmnc(:,ns2)
      ELSEWHERE
         bmn = wlo_odd * bmnc(:,ns1) + whi_odd * bmnc(:,ns2)
         gmn = wlo_odd * gmnc(:,ns1) + whi_odd * gmnc(:,ns2)
      END WHERE

      IF (lasym) THEN
         ALLOCATE(bmn2(mnmax_nyq),gmn2(mnmax_nyq))
         gmn = 0; bmn = 0
         WHERE (MOD(NINT(xm_nyq(:)),2) == 0)
            bmn2 = wlo * bmns(:,ns1) + whi * bmns(:,ns2)
            gmn2 = wlo * gmns(:,ns1) + whi * gmns(:,ns2)
         ELSEWHERE
            bmn2 = wlo_odd * bmns(:,ns1) + whi_odd * bmns(:,ns2)
            gmn2 = wlo_odd * gmns(:,ns1) + whi_odd * gmns(:,ns2)
         END WHERE
      END IF

      ! Setup theta/zeta grid
      ALLOCATE(tz(nu*nv,2))
      tz = 0
      DO zt = 1, nu*nv
         u = MOD(zt-1,nu)/REAL(nu)
         v = FLOOR(zt/REAL(nu+1))/REAL(nv)
         tz(zt,1) = u*pi2
         tz(zt,2) = v*pi2/nfp
      END DO

      ! Transform quantities (linear OK)
      ALLOCATE(b(nu*nv),g(nu*nv))
      b = 0; g = 0;
      DO zt = 1, nu*nv
         DO mn = 1, mnmax_nyq
            m = xm_nyq(mn)
            n = xn_nyq(mn)
            kernel = cos(m*tz(zt,1)-n*tz(zt,2))
            b(zt)  = b(zt) + bmn(mn)*kernel
            g(zt)  = g(zt) + gmn(mn)*kernel
         END DO
      END DO
      DEALLOCATE(bmn,gmn)

      IF (ALLOCATED(bmn2)) THEN
         DO zt = 1, nu*nv
            DO mn = 1, mnmax_nyq
               m = xm_nyq(mn)
               n = xn_nyq(mn)
               kernel = sin(m*tz(zt,1)-n*tz(zt,2))
               b(zt)  = b(zt) + bmn2(mn)*kernel
               g(zt)  = g(zt) + gmn2(mn)*kernel
            END DO
         END DO
         DEALLOCATE(bmn2,gmn2)
      END IF
      DEALLOCATE(tz)


      ! Adjust
      b = ABS(b)
      g = -g ! Jacobain negative in VMEC

      ! Flux surface quantities
      sumg = SUM(g)
      b2avg = SUM(b*b*g)/sumg
      bmax = MAXVAL(b)
      avgbobm2 = b2avg/(bmax*bmax)

      ! Normalize B
      b = b / bmax
      WHERE (b > 1) b = 1;

      ! Setup Integral
      ALLOCATE(lam(nla),avgbpov(nla))
      DO mn = 1, nla
         lam(mn) = REAL(mn-1)/REAL(nla-1)
         avgbpov(mn) = sum(sqrt(abs(1 - lam(mn) * b)) * g)
      END DO

      ! Calculate trapped and passing fractions
      fp = 0.75 * avgbobm2 * sum(lam/avgbpov)*sumg/(nla-1)
      fp = MAX(MIN(fp,1.0),0.0)
      IF (fp == 0) RETURN
      x = (1-fp)/fp

      ! Calculate G~l31
      d = 1.414*zeff +z2 &
         + x * (0.754 + 2.657*zeff + 2 * z2) &
         + x  * x * (0.348 + 1.243*zeff +     z2)
      a = 0.754 + 2.21 * zeff + z2 + x * (0.348 + 1.243 * zeff + z2)
      Gfac = (1 - x*a/d)/zeff

      !STOP

      ! Cleanup
      DEALLOCATE(b,g)
      DEALLOCATE(lam,avgbpov)

      RETURN

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE vmec_ohkawa