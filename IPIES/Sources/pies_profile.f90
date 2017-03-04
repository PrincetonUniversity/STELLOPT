!-----------------------------------------------------------------------
!     Module:        pies_profile
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/2/2011
!     Description:   This module contains the various profile quantities
!                    which define the equilibrium
!-----------------------------------------------------------------------
      MODULE pies_profile
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
!         
!----------------------------------------------------------------------
      IMPLICIT NONE
      
      REAL(rprec) :: torflux_edge
      REAL(rprec) :: curtor
      REAL(rprec), ALLOCATABLE :: press(:)
      REAL(rprec), ALLOCATABLE :: torflux(:)
      REAL(rprec), ALLOCATABLE :: torflux_norm(:)
      REAL(rprec), ALLOCATABLE :: iprime(:)
      REAL(rprec), ALLOCATABLE :: iota(:)
      TYPE(EZspline1_r8) :: p_spl
      TYPE(EZspline1_r8) :: ip_spl
      
      END MODULE pies_profile
