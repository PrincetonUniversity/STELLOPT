      SUBROUTINE allocate_wires
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE biotsavart, ONLY: single_coil
      USE Vcoilpts
      USE Vwire, ONLY: nwire, nwire1
      USE tf_coils
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: status
!-----------------------------------------------
      nwire = nwdim
      nwire1 = nwire + 1
      nwdim1 = nwdim + 1
      
      ALLOCATE (tfc_x(ncdim, nwdim1), tfc_y(ncdim, nwdim1), 
     2          tfc_z(ncdim, nwdim1), single_coil, stat=status)

      IF (status .ne. 0) 
     1   STOP 'Allocation error in COILOPT allocate_wires'

    
      END SUBROUTINE allocate_wires
