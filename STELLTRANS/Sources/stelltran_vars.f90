!-----------------------------------------------------------------------
!     Module:        stelltran_vars
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          08/29/2015
!     Description:   This module contains various internal variables.
!-----------------------------------------------------------------------
      MODULE stelltran_vars
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
!-----------------------------------------------------------------------
!     Module Variables
!          prof_length  Length of internal profiles.
!          ec           Electron charge (abs).
!          Qe           Heat flux density of electrons
!          Xe           Heat flux Coefficient of electrons
!          Qi           Heat flux density of ions
!          Xi           Heat flux coefficient of ions
!          Ge           Electron particle flux density
!          De           Electron particle flux coefficient
!          Gi           Ion particle flux density
!          Di           Ion particle flux coefficient
!          te           Electron temperature profile
!          ti           Ion temperature profile
!          ne           Electron density profile
!          ni           Ion density profile
!          Vp           V' profile
!          rho_coord    Rho coordinate corresponding to the normalized magnetic flux coordinate
!          gradrho      Gradient of rho = r/a coordinate
!          gradrhosq    Squared gradient of rho = r/a coordinate
!          Er           Ambipolar electric field profile
!          zeff         Z-effective profile
!          johm         Ohmic current profile
!          jboot        Bootstrap current profile
!          jrf          RF Current profile
!          pe_rf        Electron power from RF
!          pi_col       Ion power from collisions with electrons
!          S_pe          Source of electrons
!          S_pi         Source of ions
!          jbeam        Beam driven current profile
!----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, PARAMETER :: prof_length = 64
      REAL(rprec), PARAMETER :: ec  = 1.6021766208D-19
      REAL(rprec), PARAMETER :: me  = 9.10938356D-31
      REAL(rprec), PARAMETER :: e0  = 8.854187817D-12
      
      INTEGER :: mboz, nboz
      REAL(rprec), DIMENSION(prof_length,2) :: te, ti, ne, ni, zeff, Qe, S_pe, S_pi, Ge, Gi, Qi, Er
      REAL(rprec), DIMENSION(prof_length,2) :: De, Di, Xe, Xi, Vp, gradrhosq, Ptot
      REAL(rprec), DIMENSION(prof_length,2) :: johm,jboot,jrf,jbeam,jecrh
      REAL(rprec), DIMENSION(prof_length,2) :: dPdV_erf,dPdV_irf,dPdV_ohm, dPdV_beam, pi_col
      REAL(rprec), DIMENSION(prof_length,2) :: te_old, ti_old, ne_old, zeff_old
      REAL(rprec), DIMENSION(prof_length,2) :: johm_old,jboot_old,jrf_old,jbeam_old
      
      END MODULE stelltran_vars
