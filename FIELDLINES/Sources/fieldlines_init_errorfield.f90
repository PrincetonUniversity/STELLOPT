!-----------------------------------------------------------------------
!     Module:        fieldlines_init_errorfield
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          05/14/2021
!     Description:   This subroutine adds and errorfield to the existing
!                    background field.
!-----------------------------------------------------------------------
      SUBROUTINE fieldlines_init_errorfield
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE fieldlines_runtime
      USE fieldlines_grid, ONLY: raxis,phiaxis,zaxis, nr, nphi, nz, &
                                 rmin, rmax, zmin, zmax, phimin, &
                                 phimax, B_R, B_Z, B_PHI
      USE mpi_params
      USE mpi_inc            
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, PARAMETER :: BYTE_8 = SELECTED_INT_KIND (8)
      INTEGER(KIND=BYTE_8),ALLOCATABLE :: mnum(:), moffsets(:)
      INTEGER :: numprocs_local, mylocalid, mylocalmaster
      INTEGER :: MPI_COMM_LOCAL
      INTEGER(KIND=BYTE_8) :: chunk
      INTEGER :: ier, iunit, s, i, j, mystart, myend, k, ig
      REAL(rprec)  :: bx_err, by_err
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------

      ! Divide up Work
      numprocs_local = 1; mylocalid = master
#if defined(MPI_OPT)
      CALL MPI_COMM_DUP( MPI_COMM_SHARMEM, MPI_COMM_LOCAL, ierr_mpi)
      CALL MPI_COMM_RANK( MPI_COMM_LOCAL, mylocalid, ierr_mpi )              ! MPI
      CALL MPI_COMM_SIZE( MPI_COMM_LOCAL, numprocs_local, ierr_mpi )          ! MPI
#endif
      mylocalmaster = master

      IF (lverb) THEN
         WRITE(6,'(A)')   '----- ERROR FIELD Information -----'
         DO ig = 1, 20
            IF (errorfield_amp(ig) .eq. 0) CYCLE
            WRITE(6,'(A,I2.2,A,E11.4,A,F5.3,A)') 'n = ',ig,'; B = ',errorfield_amp(ig),&
                          ' T; phase = ',errorfield_phase(ig), 'rad'
         END DO
         CALL FLUSH(6)
      END IF

      ! Break up the Work
      CALL MPI_CALC_MYRANGE(MPI_COMM_LOCAL,1, nr*nphi*nz, mystart, myend)

      ! Loop
      DO s = mystart, myend
         i = MOD(s-1,nr)+1
         j = MOD(s-1,nr*nphi)
         j = FLOOR(REAL(j) / REAL(nr))+1
         k = CEILING(REAL(s) / REAL(nr*nphi))
         bx_err = 0
         by_err = 0
         DO ig = 1, 20
            bx_err = errorfield_amp(ig) * COS(ig*phiaxis(j)+errorfield_phase(ig))
            by_err = errorfield_amp(ig) * SIN(ig*phiaxis(j)+errorfield_phase(ig))
            B_R(i,j,k) = B_R(i,j,k) + bx_err*cos(phiaxis(j)) + by_err*sin(phiaxis(j))
            B_PHI(i,j,k) = B_PHI(i,j,k) + by_err*cos(phiaxis(j)) - bx_err*sin(phiaxis(j))
         END DO
      END DO 

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init_coil1',ierr_mpi)
      CALL MPI_COMM_FREE(MPI_COMM_LOCAL,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'fieldlines_init_coil: MPI_COMM_LOCAL',ierr_mpi)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init_coil',ierr_mpi)
#endif
      
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE fieldlines_init_errorfield
