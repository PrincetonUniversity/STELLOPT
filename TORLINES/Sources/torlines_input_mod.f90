!-----------------------------------------------------------------------
!     Module:        pies_input_mod
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/7/2011
!     Description:   This module contains the PIES input namelist and
!                    subroutine which initializes and reads the
!                    PIES input namelist.
!-----------------------------------------------------------------------
      MODULE torlines_input_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds
      USE torlines_runtime
      USE torlines_background
      USE torlines_realspace, ONLY: nu, nv, k
      USE torlines_fieldlines, ONLY: nlines
      USE virtual_casing_mod, ONLY: nu_vc, nv_vc
      
!-----------------------------------------------------------------------
!     Module Variables
!         
!-----------------------------------------------------------------------
      IMPLICIT NONE
      
!-----------------------------------------------------------------------
!     Input Namelists
!         &fieldlines_input
!            k              Number of radial gridpoints
!            nu             Number of poloidal gridpoints
!            nv             Number of toroidal gridpoints
!            r_start        Radial starting locations for fieldlines [m]
!            phi_start      Toroidal starting locations for fieldlines [radians]
!            phi_end        Toroidal ending locations for fieldlines [radians]
!            z_start        Vertical starting locations for fieldlines [m]
!            npoinc         Number of points (per field period) in which data is saved
!            dphi           Fieldlines following stepsize [radians]
!            r_hc           Initial guess for homoclinic tangle [m]
!            z_hc           Initial guess for homoclinic tangle [m]
!            phi_hc         Initial guess for homoclinic tangle [radians]
!            num_hcp        Number of points for homoclinic tangle
!            delta_hcp      Length of initial line for homoclinic tangle
!            follow_tol     Tollerance for fieldline following (LSODE and NAG)
!            vc_adapt_tol   Tollerance for adaptive integration using Virtual casing
!                           (note set to negative value to use non-adaptive integration)
!            int_type       Field line integration method
!                           'NAG','LSODE','RKH68'
!
!            NOTE:  Some grid parameters may be overriden (such as
!                   phimin and phimax) to properly represent a given
!                   field period.
!-----------------------------------------------------------------------
      NAMELIST /torlines_input/   k, nu, nv, bound_separation, nu_vc, nv_vc,&
                                  r_start, phi_start, phi_end,z_start,&
                                  npoinc, dphi, follow_tol,&
                                  vc_adapt_tol, int_type, lvc_field, &
                                  r_hc, phi_hc, z_hc, num_hcp, delta_hc, &
                                  NZONE_EMC3, NPOLO_EMC3, NTORO_EMC3, &
                                  S0_EMC3, S1_EMC3
!-----------------------------------------------------------------------
!     Subroutines
!         read_pies_input:   Reads pies_input namelist
!-----------------------------------------------------------------------
      CONTAINS
      
      SUBROUTINE read_torlines_input(iunit, istat)
      IMPLICIT NONE
      INTEGER :: iunit, istat
      ! Initializations
      istat = 0
      k = 50
      nu = -1
      nv = -1
      nu_vc = -1
      nv_vc = -1
      r_start   = -1
      z_start   = -1
      phi_start = -1
      phi_end   = -1
      r_hc      = -1
      z_hc      = -1
      phi_hc    = -1
      follow_tol = 1e-7
      lfreeb = .true.
      lvc_field = .true.
      bound_separation = 1.0
      int_type = 'NAG'
      phi_end = -1
      npoinc = 36
      NZONE_EMC3 = 12
      NPOLO_EMC3 = 361
      NTORO_EMC3 = 11
      S0_EMC3    = 0.5 
      S1_EMC3    = 0.8
      ! Read namelist
      READ(iunit,NML=torlines_input,IOSTAT=istat)
      IF (istat /= 0) CALL handle_err(NAMELIST_READ_ERR,'torlines_input in: input.'//TRIM(id_string),istat)
      ! Adjust EMC3 VARS
      NZONET   = NZONE_EMC3
      SRF_POLO = NPOLO_EMC3
      SRF_TORO = NTORO_EMC3
      s_inner  = S0_EMC3
      s_outer  = S1_EMC3
      END SUBROUTINE read_torlines_input
      
      END MODULE torlines_input_mod
