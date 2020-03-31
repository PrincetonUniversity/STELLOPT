      SUBROUTINE read_bcoils
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE bcoils_mod
      USE Vcoilpts
      USE safe_open_mod
      USE biotsavart
      USE mpi_params                                         !mpi stuff
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ierr, icoils=21
      INTEGER :: i, n, mbw
      LOGICAL :: lcoil_open
!-----------------------------------------------

      CALL cleanup_biotsavart
      INQUIRE(FILE=TRIM(bcoil_file),NUMBER=icoils,
     1        OPENED=lcoil_open)
      IF (lcoil_open) CLOSE(icoils)
      CALL parse_coils_file(TRIM(bcoil_file))
      mbcoils = SIZE(coil_group) !SAL
      IF (myid .eq. master) THEN
         WRITE(6,'(A)')   '----- COILS Information -----'
         WRITE(6,'(A,A)') '   FILE: ',TRIM(bcoil_file)
         WRITE(6,'(A,I3)')'   Coil Periodicity: ',nfp_bs
         WRITE(6,'(A,I3)')'   Current Systems: ',mbcoils
      END IF
!      CALL safe_open (icoils, ierr, TRIM(bcoil_file), 'old','formatted')
!      IF (ierr .ne. 0) THEN
!         IF (myid .EQ. master) PRINT *, 
!     1   'Error opening bcoil_file: ierr = ', ierr,
!     2   ' Setting lbcoil=F'
!!         STOP
!         lbcoil=.FALSE.
!         RETURN
!      END IF

!      ALLOCATE (bcoil_x(ncdim, nwdim1), bcoil_y(ncdim, nwdim1), 
!     1          bcoil_z(ncdim, nwdim1), stat=ierr)
!      IF (ierr .ne. 0) STOP 'Allocation error in read_bcoils'

!     READ coil filament coordinates

!      READ (icoils, *) mbcoils
!      IF (mbcoils .gt. ncdim) THEN
!         IF( myid .eq. master ) THEN                         !mpi stuff
!           PRINT *, 'Error: no. of bg coils > ncdim ', mbcoils
!         END IF     !IF( myid .eq. master )                  !mpi stuff
!         STOP
!      END IF

!      DO n=1,mbcoils
!         READ (icoils, *) mbwires(n)
!         mbw = mbwires(n) + 1
!         IF (mbw .gt. nwdim1) THEN
!            IF( myid .eq. master ) THEN                         !mpi stuff
!               PRINT *, 'Error: no. of bg wires > nwdim ', mbw
!            END IF     !IF( myid .eq. master )                  !mpi stuff
!            STOP
!         END IF
!         DO i=1,mbw
!            READ (icoils, *) bcoil_x(n,i), bcoil_y(n,i), bcoil_z(n,i)
!         END DO
!      END DO

!      CLOSE (icoils)

      END SUBROUTINE read_bcoils
