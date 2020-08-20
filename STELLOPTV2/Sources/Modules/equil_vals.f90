!-----------------------------------------------------------------------
!     Module:        equil_vals
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/24/2012
!     Description:   This module contains the various values which
!                    define an equilibria.
!-----------------------------------------------------------------------
      MODULE equil_vals
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE EZspline_obj
!-----------------------------------------------------------------------
!     Module Variables
!            aspect    Equilibrium Aspect Ratio
!              beta    Equilibrium Beta
!            curtor    Net Toroidal Current
!           phiedge    Total Enclosed toroidal flux (edge value)
!            volume    Equilibrium Volume
!                wp    Equilibrium Stored Energy
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL     ::  lasym
      INTEGER     ::  bcs0(2) = (/ 0, 0/)
      INTEGER     ::  bcs1(2) = (/-1,-1/)
      INTEGER     ::  nrad, nfp
      REAL(rprec) ::  aspect, betat, curtor, phiedge, volume, wp, drho,&
                      rbtor, r0, z0, iota_res_tgt, betap, beta, Rmajor, &
                      Aminor, mach0, kx_gene, kink_omega, Baxis
      REAL(rprec),ALLOCATABLE :: rho(:), shat(:), extcur(:), eff_ripple(:), &
                                 orbit_lost_frac(:), radto_ece(:,:), radtx_ece(:,:)
      REAL(rprec),ALLOCATABLE :: balloon_grate(:,:,:)
      REAL(rprec),ALLOCATABLE :: txport_q(:,:,:), txport_q_all(:,:,:,:)
      TYPE(EZspline1_r8) :: pres_spl, iota_spl, phi_spl, ip_spl, V_spl, &
                            te_spl, ne_spl, ti_spl, th_spl, ah_spl, &
                            jdotb_spl, zeff_spl, jcurv_spl, omega_spl,&
                            nustar_spl, emis_xics_spl
      REAL(rprec), ALLOCATABLE :: wp_kink(:), wk_kink(:), omega_kink(:),&
                     growth_kink(:)

      END MODULE equil_vals
