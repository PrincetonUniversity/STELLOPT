      SUBROUTINE allocate_vf_coils
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vf_coils
      USE coils
      USE mpi_params                                         !mpi stuff
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: status
!-----------------------------------------------

!     ALLOCATE vf coil arrays

      ALLOCATE (x_vf(nwdim1,3,ncdim), y_vf(nwdim1,3,ncdim), 
     1          z_vf(nwdim1,3,ncdim), vertical(num_vf), stat=status)
      IF( myid .eq. master ) THEN                            !mpi stuff
        IF(status /= 0) STOP "Cannot ALLOCATE vertical"
      END IF     !IF( myid .eq. master )                     !mpi stuff

      END SUBROUTINE allocate_vf_coils
