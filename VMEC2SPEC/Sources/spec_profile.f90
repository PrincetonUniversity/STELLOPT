!-----------------------------------------------------------------------
!     Module:        spec_profile
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          04/09/2012
!     Description:   This module contains the various profile quantities
!                    which define the equilibrium
!-----------------------------------------------------------------------
      MODULE spec_profile
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds
      USE EZspline_obj
      USE EZspline
!-----------------------------------------------------------------------
!     Module Variables
!          torflux_edge Toridal flux at the edge of the plasma
!          press        Pressure
!          torflux      Toroidal Flux
!          iprime       I'=dI/dPsi
!          p_spl        Pressure Spline Object
!          i_spl        I' Spline Object
!----------------------------------------------------------------------
      IMPLICIT NONE
      
      REAL(rprec) :: torflux_edge
      REAL(rprec) :: curtor
      REAL(rprec) :: pscale
      REAL(rprec), ALLOCATABLE :: press(:)
      REAL(rprec), ALLOCATABLE :: adiab(:)
      REAL(rprec), ALLOCATABLE :: currf(:)
      REAL(rprec), ALLOCATABLE :: iotaf(:)
      REAL(rprec), ALLOCATABLE :: torflux_pies(:)
      REAL(rprec), ALLOCATABLE :: torflux_norm(:)
      REAL(rprec), ALLOCATABLE :: iprime(:)
      REAL(rprec), ALLOCATABLE :: p_cubspl(:,:)
      REAL(rprec), ALLOCATABLE :: ip_cubspl(:,:)
      TYPE(EZspline1_r8) :: p_spl
      TYPE(EZspline1_r8) :: ip_spl
      TYPE(EZspline1_r8) :: iota_spl
      
      END MODULE spec_profile
