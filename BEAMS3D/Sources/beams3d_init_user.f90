!-----------------------------------------------------------------------
!     Module:        beams3d_init_user
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          05/21/2025
!     Description:   This subroutine allows one to create a user
!                    defined scenario for doing example plots.
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_init_user
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE beams3d_runtime
      USE beams3d_grid, ONLY: raxis,phiaxis,zaxis, nr, nphi, nz, &
                                 rmin, rmax, zmin, zmax, phimin, &
                                 phimax, B_R, B_Z, B_PHI, S_ARR
      USE mpi_params  
      USE mpi_inc      
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, PARAMETER :: BYTE_8 = SELECTED_INT_KIND (8)
#if defined(MPI_OPT)
      INTEGER(KIND=BYTE_8),ALLOCATABLE :: mnum(:), moffsets(:)
      INTEGER :: numprocs_local, mylocalid, mylocalmaster
      INTEGER :: MPI_COMM_LOCAL
#endif
      INTEGER(KIND=BYTE_8) :: chunk
      INTEGER :: ier, iunit, s, i, j, mystart, myend, k, ik, ig, coil_dex
      REAL(rprec)  :: br, bphi, bz, current, current_first, &
                      br_temp, bphi_temp, bz_temp
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------

      ! Divide up Work
#if defined(MPI_OPT)
      CALL MPI_COMM_DUP( MPI_COMM_SHARMEM, MPI_COMM_LOCAL, ierr_mpi)
      CALL MPI_COMM_RANK( MPI_COMM_LOCAL, mylocalid, ierr_mpi )              ! MPI
      CALL MPI_COMM_SIZE( MPI_COMM_LOCAL, numprocs_local, ierr_mpi )          ! MPI
#endif
      mylocalmaster = master

      IF (lverb) THEN
         WRITE(6,'(A)')   '----- USER Defined Scenario -----'
         CALL FLUSH(6)
      END IF

      ! Break up the Work
      CALL MPI_CALC_MYRANGE(MPI_COMM_LOCAL, 1, nr*nphi*nz, mystart, myend)
      
      ! Here is where you can set the mangetic field.
      DO s = mystart, myend
         i = MOD(s-1,nr)+1
         j = MOD(s-1,nr*nphi)
         j = FLOOR(REAL(j) / REAL(nr))+1
         k = CEILING(REAL(s) / REAL(nr*nphi))
         ! Pure vertical field
         B_R(i,j,k)   = 0.0
         B_PHI(i,j,k) = 0.0
         B_Z(i,j,k)   = 1.0+10.0*ABS(zaxis(k)/zaxis(nz))**2.0
         S_ARR(i,j,k) = 0.0
         IF (lverb .and. (MOD(s,nr) == 0)) THEN
            CALL backspace_out(6,6)
            WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*s)/(myend-mystart+1)),']%'
            CALL FLUSH(6)
         END IF
      END DO
      
      ! Clean up the progress bar
      IF (lverb) THEN
         CALL backspace_out(6,36)
         WRITE(6,'(38X)',ADVANCE='no')
         CALL backspace_out(6,36)
         WRITE(6,*)
         CALL FLUSH(6)
      END IF    

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
      CALL MPI_COMM_FREE(MPI_COMM_LOCAL,ierr_mpi)
      CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'beams3d_init_user',ierr_mpi)
#endif
      
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE beams3d_init_user
