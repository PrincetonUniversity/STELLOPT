!-----------------------------------------------------------------------
!     Module:        stelltran_data
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          08/21/2015
!     Description:   This module contains various parameters related to
!                    the code parameters.
!-----------------------------------------------------------------------
      MODULE stelltran_data
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
!-----------------------------------------------------------------------
!     Module Variables
!          ntimesteps   Number of timesteps in database file
!         
!          te_R/Z/PHI   Location of Electron Temperature (R/Z/PHI)
!          te_f         Electron Temperature
!          ti_R/Z/PHI   Location of Ion Temperature (R/Z/PHI)
!          ti_f         Ion Temperature
!          ne_R/Z/PHI   Location of Electron Density (R/Z/PHI)
!          ne_f         Electron Density
!          zeff_R/Z/PHI Location of Z-effective (R/Z/PHI)
!          zeff_f       Z-effective
!----------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER :: ntimesteps, nne, nte, nti, nzeff, nprof, necrh
      REAL(rprec) :: dt
      REAL(rprec), ALLOCATABLE :: Vloop(:)
      REAL(rprec), ALLOCATABLE :: Ip(:)
      REAL(rprec), ALLOCATABLE :: te_s(:), ti_s(:), ne_s(:), zeff_s(:)
      REAL(rprec), ALLOCATABLE :: te_f(:,:),ti_f(:,:),ne_f(:,:),zeff_f(:,:)
      REAL(rprec), ALLOCATABLE :: P_ecrh(:,:)
      REAL(rprec), ALLOCATABLE :: nu_star(:,:)
      REAL(rprec), ALLOCATABLE :: johm_sav(:,:),jboot_sav(:,:),jecrh_sav(:,:),jbeam_sav(:,:)
      REAL(rprec), ALLOCATABLE :: dPdVohm_e(:,:),dPdVohm_i(:,:)
      REAL(rprec), ALLOCATABLE :: dPdVecrh_e(:,:),dPdVicrh_i(:,:)
      REAL(rprec), ALLOCATABLE :: dPdVbeam_e(:,:),dPdVbeam_i(:,:)
      
      END MODULE stelltran_data
