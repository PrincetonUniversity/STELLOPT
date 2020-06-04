!-----------------------------------------------------------------------
!     Module:        beams3d_lines
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/21/2012
!     Description:   This module contains the FIELDLINES field line
!                    variables.
!-----------------------------------------------------------------------
      MODULE beams3d_lines
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      
!-----------------------------------------------------------------------
!     Module Variables
!          myline    Dummy index
!          nparticles    Number of Particles
!          nsteps    Number of integration steps along fieldline
!          R_lines   Radial locations along fieldline [m] (npoinc per field period)
!          Z_lines   Vertical locations along field line [m]
!          PHI_lines Toroidal locations along field line [radians]
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL  ::  ltherm
      INTEGER,PARAMETER  :: ns_prof1 = 128 !rho-Grid
      INTEGER,PARAMETER  :: ns_prof2 = 16 !u-grid
      INTEGER,PARAMETER  :: ns_prof3 = 16 !v-grid
      INTEGER,PARAMETER  :: ns_prof4 = 64 !vll-grid
      INTEGER,PARAMETER  :: ns_prof5 = 32 !vperp-grid
      INTEGER  :: nparticles, nsteps, myline, mybeam, mytdex, myend, mystart_save,myend_save
      REAL(rprec) :: xlast,ylast,zlast ! for storing position
      REAL(rprec) :: moment, mycharge, myZ, mymass, myv_neut(3), &
                     B_temp(4), rand_prob, cum_prob, tau, next_t, &
                     dt_out, partvmax, fact_crit, fact_pa, fact_vsound, partpmax
      LOGICAL, ALLOCATABLE     :: neut_lines(:,:)
      INTEGER, ALLOCATABLE     :: end_state(:)
      REAL(rprec), ALLOCATABLE :: shine_through(:), shine_port(:)
      REAL(rprec), ALLOCATABLE :: ndot_prof(:,:),epower_prof(:,:), &
                                  ipower_prof(:,:),j_prof(:,:), dense_prof(:,:), dist2d_prof(:,:,:)
!      REAL(rprec), ALLOCATABLE :: dist_prof(:,:,:,:,:,:)
      REAL(rprec), ALLOCATABLE :: R_lines(:,:),Z_lines(:,:),PHI_lines(:,:),vll_lines(:,:),moment_lines(:,:),&
                                  S_lines(:,:),U_lines(:,:),B_lines(:,:)

      END MODULE beams3d_lines
