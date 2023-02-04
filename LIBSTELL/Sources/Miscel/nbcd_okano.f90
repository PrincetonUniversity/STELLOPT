!-----------------------------------------------------------------------
!     Subroutine:    nbcd_okano
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          02/03/2023
!     Description:   This subroutine calculates the neutral beam
!                    current drive using the Okano model.  The details
!                    of which can be found in:
!                    K. Okano "Neoclassical formula for netural beam
!                       current drive." 1990 Nucl. Fusion 30 423
!                    https://dx.doi.org/10.1088/0029-5515/30/3/004
!                    It returns j/p [A*m/W]
!                               j [A/m^2]
!                               p [W/m^3]
!-----------------------------------------------------------------------
      SUBROUTINE nbcd_okano(ne,Te,Ti,mass,Zeff,InvAspect,G,Abeam,Vbeam,Zbeam,alpha,theta,jopb)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     Arguments
!          ne			Electron Density [m^-3]
!          Te           Electron Temperature [eV]
!          Ti           Ion Temperature [eV]
!          mass         Plasma mass
!          Zeff         Sum(ni*Zi*Zi)/ne
!          InvAspect    Inverse Aspect Ratio r/R
!          G            Trapped Electron Correction Factor
!          Abeam        Beam mass number
!          Vbeam        Beam velocity [m/s]
!          Zbeam        Beam Charge Number
!          alpha        Beam pitch vll/vtotal
!          theta        Beam poloidal angle
!          jopb         j/pb Current density over beam Power [A*m/W]
!-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: ne, Te, Ti, mass, Zeff, InvAspect, G, Abeam, Vbeam, Zbeam, alpha, theta
      DOUBLE PRECISION, INTENT(out) :: jopb
!-----------------------------------------------------------------------
!     PARAMETERS
!        pi2            2*pi
!        c              speed of light [m/s]
!        electron_mass  Mass of an electron [kg]
!        sqrt_pi        pi^(1/2)
!        bm             Beta_m parameters (eq 8b from paper)
!-----------------------------------------------------------------------
      DOUBLE PRECISION, PARAMETER :: pi2 = 6.283185482025146D+00
      DOUBLE PRECISION, PARAMETER :: c   = 299792458D+00
      DOUBLE PRECISION, PARAMETER :: electron_mass = 9.10938356D-31
      DOUBLE PRECISION, PARAMETER :: e_charge      = 1.60217662E-19
      DOUBLE PRECISION, PARAMETER :: sqrt_pi       = 1.7724538509
      DOUBLE PRECISION, PARAMETER, DIMENSION(9) :: bm = &
        (/-0.02338, 8.973,-29.37, &
          36.15, 152.2,-566.6, &
          773.0,-479.7, 112.3 /)
!-----------------------------------------------------------------------
!     Variables
!          m            Summation index over beta_m
!          x            Normalized Velocity (vbeam/vcritical)
!          y            y parameter (eq 4)
!          v_c          Critical Velocity Wesson formulation
!          eps0         Pitch angle parameter (eq 5)
!          t            Neoclassical parameter
!          C_sigma      Sigma coefficient
!          sigma        Neoclassical parameter
!          C_b          B coefficient
!          b            Neoclassical parameter
!          Fnc          Neoclassical factor 1-b*(r/R)**sigma
!          Fed          Energy diffusion factor
!          Fcx          Charge Exchange Factor (not implemented)
!          J0           Current Function
!-----------------------------------------------------------------------
      INTEGER          :: m
      DOUBLE PRECISION :: x, y, v_c, eps0, t, C_sigma, sigma, C_b, b, &
                          Fnc, Fed, Fcx, J0
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      jopb = 0
      v_c = SQRT(Te)*SQRT(2*e_charge/mass) * &
            (0.75*sqrt_pi*sqrt(mass/electron_mass))**(1.0/3.0) ! Wesson pg 226 5.4.9
      x = Vbeam/v_c
      y = 4 * Zeff / (5 * Abeam)
      ! Fast ion pitch angle
      eps0 = SQRT(1-(1-alpha*alpha)*(1-InvAspect)/(1-cos(theta)))
      ! Neoclassical correction factor
      t = (1-eps0)/(1.273 - 0.240*x + 0.00667*x*x)
      IF (y < 0.8) THEN
         C_sigma = (0.0709/y - 0.107 + 0.0238*y)*x + 0.294/y + 0.538 + 0.119*y
      ELSE IF (y <= 1.2) THEN
         C_sigma = 0.365/y + 0.430 + 0.143*y
      ELSE
      END IF
      sigma = C_sigma*(0.8-0.05*x+0.01*x*x)
      IF (y < 0.8) THEN
         C_b = (0.148*(8.42-x)*(0.8-y)*((0.94+x)*(1.2-y)+0.541)*(1/eps0 - 1)+1)**sigma
      ELSE IF (y <= 1.2) THEN
         C_b = ((28-x+3*x*x)/23.4*(eps0-0.6)*EXP(-(7-10*eps0)**2)-(1-eps0)*(10-x+x*x)/10)*(y-0.8) + 1
      END IF
      DO m = 1, 9
         b = b + bm(m)*t**(m-1)
      END DO
      b = (b*b + 0.775*x**(-0.85) + 0.125*x)*C_b
      Fnc = MAX(1 - b * InvAspect**sigma,0.0)
      ! Energy Diffusion
      Fed = 1+(-0.002207/x**6 + 0.01545/x**4 + 0.03727/x**2)*(1+0.072*(y-0.4)*EXP(-0.55*x**0.95))*(1+4000*(Ti/Te - 1)*EXP(-9*x**0.3))
      ! Charge Exchange
      Fcx = 1
      ! Current Function
      J0 = x*x/(4 + 3*y + x*x*(x+1.39+1.61*y**0.7))
      ! j/Pb
      jopb = 1.58*Te*1D-3*eps0*Fnc*(1-Zbeam*(1-G)/Zeff)*Fed*Fcx*J0/(ne*1D-20)

      RETURN

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE nbcd_okano


