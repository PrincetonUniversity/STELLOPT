!-----------------------------------------------------------------------
!     Module:        pies_input_mod
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/7/2011
!     Description:   This module contains the PIES input namelist and
!                    subroutine which initializes and reads the
!                    PIES input namelist.
!-----------------------------------------------------------------------
      MODULE pies_input_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds
      USE pies_runtime
      USE pies_realspace, ONLY: nu, nv
      USE pies_background, ONLY: m, n, k
      USE pies_fieldlines, ONLY: follow_tol, magaxis_abs, magaxis_eta,&
                                 dkmin,nu_spline,nv_spline
      USE virtual_casing_mod, ONLY: nu_vc, nv_vc
      
!-----------------------------------------------------------------------
!     Module Variables
!         
!-----------------------------------------------------------------------
      IMPLICIT NONE
      
!-----------------------------------------------------------------------
!     Input Namelists
!         pies_input
!-----------------------------------------------------------------------
      NAMELIST /pies_input/ nu,nv,m,n,k,lfreeb,extsurfs,lideal,&
                            magaxis_abs, magaxis_eta, follow_tol,&
                            dkmin,nu_spline,nv_spline, nu_vc, nv_vc
      
!-----------------------------------------------------------------------
!     Subroutines
!         read_pies_input:   Reads pies_input namelist
!-----------------------------------------------------------------------
      CONTAINS
      
      SUBROUTINE read_pies_input(iunit, istat)
      INTEGER :: iunit, istat
      ! Initializations
      istat = 0
      iter = 0
      k = -1
      m = -1
      n = -1
      nu = -1
      nv = -1
      nu_vc = -1
      nv_vc = -1
      follow_tol = 1e-7
      magaxis_abs = 1e-7
      dkmin = 0.004
      magaxis_eta = 0.0
      lfreeb = .true.
      extsurfs = 3
      lideal = .true.
      nu_spline = 90
      nv_spline = 90
      ! Read namelist
      READ(iunit,NML=pies_input,IOSTAT=istat)
      IF (istat /= 0) CALL handle_err(NAMELIST_READ_ERR,'pies_input in: input.'//TRIM(id_string),istat)
      END SUBROUTINE read_pies_input
      
      END MODULE pies_input_mod
