!-----------------------------------------------------------------------
!     Module:        pies_realspace
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/2/2011
!     Description:   This module contains the realspace arrays.
!-----------------------------------------------------------------------
      MODULE pies_realspace
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds
!-----------------------------------------------------------------------
!     Module Variables
!          nu         Number of real space poloidal points
!          nv         Number of real space toroidal points
!          nuv        nu*nv
!          nuvp       nu*nv*np
!          nv_vmec    Number of toroidal points in vacuum field file
!          xu         Poloidal Points Number Array [0,1]
!          xv         Toroidal Points Number Array [0,1]
!          rreal      Real Space R
!          zreal      Real Space Z
!          bsreal     Real Space B^s
!          bureal     Real Space B^u
!          bvreal     Real Space B^v
!          rs         Real Space dR/ds
!          ru         Real Space dR/du
!          rv         Real Space dR/dv
!          zs         Real Space dZ/ds
!          zu         Real Space dZ/du
!          zv         Real Space dZ/dv
!         
!----------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER                  :: nu, nv, nuv, nuvp, nv_vmec
      REAL(rprec), ALLOCATABLE :: xu(:), xv(:)
      REAL(rprec), ALLOCATABLE :: rreal(:,:,:), zreal(:,:,:)
      REAL(rprec), ALLOCATABLE :: bsreal(:,:,:), bureal(:,:,:)
      REAL(rprec), ALLOCATABLE :: bvreal(:,:,:)
      REAL(rprec), ALLOCATABLE :: rs(:,:,:),zs(:,:,:)
      REAL(rprec), ALLOCATABLE :: ru(:,:,:),rv(:,:,:)
      REAL(rprec), ALLOCATABLE :: zu(:,:,:),zv(:,:,:)
      REAL(rprec), ALLOCATABLE :: g11(:,:,:),g12(:,:,:),g13(:,:,:)
      REAL(rprec), ALLOCATABLE :: g22(:,:,:),g23(:,:,:),g33(:,:,:)
      REAL(rprec), ALLOCATABLE :: detg(:,:,:)
      REAL(rprec), ALLOCATABLE :: jsreal(:,:,:), jureal(:,:,:)
      REAL(rprec), ALLOCATABLE :: jvreal(:,:,:)
      
      END MODULE pies_realspace
