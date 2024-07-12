!-----------------------------------------------------------------------
!     Subroutine:    fieldlines_init_nescoil
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          04/28/2016
!     Description:   This subroutine reads the NESCOIL output file
!                    and calculates the field from that file.
!-----------------------------------------------------------------------
      SUBROUTINE fieldlines_init_nescoil
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE safe_open_mod
      USE fieldlines_runtime
      USE fieldlines_grid, ONLY: raxis_g => raxis, phiaxis, &
                                 zaxis_g => zaxis, nr, nphi, nz, &
                                 rmin, rmax, zmin, zmax, phimin, &
                                 phimax, B_R, B_Z, B_PHI,&
                                 BR_spl, BZ_spl
      USE fieldlines_lines, ONLY: nlines
      USE read_nescoil_mod
      USE mpi_params
      USE mpi_inc
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: mylocalid, mylocalmaster
      INTEGER :: MPI_COMM_LOCAL
      INTEGER :: ier, iunit, s, i, j, k, mystart, myend
      REAL(rprec)  :: x, y, z, bx, by, bz
      INTEGER, PARAMETER :: nu_local = 128
      INTEGER, PARAMETER :: nv_local = 90
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------

      ! Divide up Work
      mylocalid = myworkid
#if defined(MPI_OPT)
      CALL MPI_COMM_DUP( MPI_COMM_SHARMEM, MPI_COMM_LOCAL, ierr_mpi)
      CALL MPI_COMM_RANK( MPI_COMM_LOCAL, mylocalid, ierr_mpi )              ! MPI
#endif
      mylocalmaster = master

      ! Read the NESCOIL FILE
      IF (mylocalid == mylocalmaster) THEN
         CALL read_nescout('nescout.'//TRIM(nescoil_string),ier)
         IF (ier /= 0) STOP "ERROR reading nescout file."
      END IF
      CALL nescoil_bfield_init(nu_local,nv_local,MPI_COMM_LOCAL)
      
      IF (lverb) THEN
         CALL nescoil_info(6)
         WRITE(6,'(5X,A,I3.3,A)',ADVANCE='no') 'Vacuum Field Calculation [',0,']%'
         CALL FLUSH(6)
      END IF
      
      ! Break up the Work
      CALL MPI_CALC_MYRANGE(MPI_COMM_LOCAL,1, nr*nphi*nz, mystart, myend)

      IF (lafield_only) THEN
            ! This is not supported
      ELSE
         DO s = mystart, myend
            i = MOD(s-1,nr)+1
            j = MOD(s-1,nr*nphi)
            j = FLOOR(REAL(j) / REAL(nr))+1
            k = CEILING(REAL(s) / REAL(nr*nphi))
            x = raxis_g(i)*COS(phiaxis(j))
            y = raxis_g(i)*SIN(phiaxis(j))
            z = zaxis_g(k)
            bx = 0; by = 0; bz = 0;
            CALL nescoil_bfield(x,y,z,bx,by,bz)
            B_R(i,j,k)   = bx * COS(phiaxis(j)) + by *SIN(phiaxis(j))
            B_PHI(i,j,k) = by * COS(phiaxis(j)) - bx * SIN(phiaxis(j))
            B_Z(i,j,k)   = bz
            IF (lverb .and. (MOD(s,nr) == 0)) THEN
               CALL backspace_out(6,6)
               WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*s)/(myend-mystart+1)),']%'
               CALL FLUSH(6)
            END IF
         END DO
      END IF
      
      IF (lverb) THEN
         CALL backspace_out(6,36)
         CALL FLUSH(6)
         WRITE(6,'(36X)',ADVANCE='no')
         CALL FLUSH(6)
         CALL backspace_out(6,36)
         CALL FLUSH(6)
      END IF

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init_nescoil1',ierr_mpi)
#endif

      ! Deallocate
      CALL read_nescout_deallocate(MPI_COMM_LOCAL)

      ! Fix any B_PHI==0 points
      IF (mylocalid == mylocalmaster) WHERE(B_PHI == 0) B_PHI = 1.0

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init_nescoil2',ierr_mpi)
      CALL MPI_COMM_FREE(MPI_COMM_LOCAL,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'fieldlines_init_nescoil3: MPI_COMM_LOCAL',ierr_mpi)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init_nescoil4',ierr_mpi)
#endif
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE fieldlines_init_nescoil
