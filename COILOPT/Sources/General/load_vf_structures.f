      SUBROUTINE load_vf_structures
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vf_coils
      USE coils
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i
!-----------------------------------------------

!     load the coil parameters with values from variables

      DO i=1, num_vf
         vertical(i)%current = cc_vf(i)
      END DO

      END SUBROUTINE load_vf_structures
