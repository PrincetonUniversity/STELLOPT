!===============================================================================
      MODULE ip_global
!-----------------------------------------------------------------------------
! module to define globally available variables for the Faraday Rotation module
!
! integrand_toggle : set integrand to Faraday Rotation (1) or Phase Shift (2)
! ibeam : index to specify the particular interferometer/polarimeter beamline
! intpol : TYPE to store all relevant info about all beams
!
!-------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
      USE v3_utilities       ! assert_eq is here
      USE plasma_edge
      USE ip_beamline
      USE signal_T
      USE diagnostic_T
      IMPLICIT NONE

      INTEGER, SAVE                             :: integrand_toggle
      INTEGER, SAVE                             :: ibeam
      INTEGER, SAVE                             :: n_d_desc_g
      INTEGER, SAVE                             :: n_s_desc_g
      TYPE(int_pol_coll), SAVE                       :: ipcoll_g
      TYPE(int_pol), SAVE                            :: intpol
      TYPE(ip_beam), SAVE                            :: ipbeam_g
      TYPE(pv_edge), ALLOCATABLE, DIMENSION(:), SAVE  :: myedge
      TYPE(diagnostic_desc), DIMENSION(:), ALLOCATABLE, TARGET ::              &
     &    d_desc_g
      TYPE(signal_desc), DIMENSION(:), ALLOCATABLE, TARGET ::                  &
     &    s_desc_g

      END MODULE ip_global
