      SUBROUTINE read_bnorm_coefs
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE boundary
      USE bnorm_mod
      USE mpi_params                                         !mpi stuff
      USE safe_open_mod
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l  P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: imnmax = 1000
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ierr, ibnorm = 22
      INTEGER :: i, im, in
      REAL (rprec) :: bn
!-----------------------------------------------
      IF (myid .eq. master)
     1  WRITE(6,'(3X,A,A)') 'FILE: ',TRIM(bnorm_file)
      CALL safe_open(ibnorm, ierr, bnorm_file, 'old', 'formatted')
      IF (ierr .ne. 0) THEN
         IF( myid .eq. master ) THEN
           PRINT *, 'Error opening ' // TRIM(bnorm_file)
     1                               // ': ierr =', ierr
           PRINT *, 'Run xnorm wout.filename to generate bnorm file'
         END IF
         STOP
      END IF

!     READ and load coefficients into boundary

      mnbn_max = 0
      DO i=1,imnmax
          READ (ibnorm, *, END=120, iostat=ierr) im, in, bn
          mnbn_max = i
          xbn_m(i) = im
          xbn_n(i) = in*nfp                                  ! note nfp
          bn_coef(i) = bn
      END DO
  120 CONTINUE

      CLOSE(ibnorm)
      IF (myid .eq. master)
     1  WRITE(6,'(3X,A,I6)') '# of coefs: ',mnbn_max
!      IF( myid .eq. master )
!     1   PRINT *, 'number of bnorm coefficients READ = ', mnbn_max

      END SUBROUTINE read_bnorm_coefs
