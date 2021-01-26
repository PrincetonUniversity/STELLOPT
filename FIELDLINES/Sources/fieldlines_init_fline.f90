!-----------------------------------------------------------------------
!     Module:        fieldlines_init_fline
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          10/30/2020
!     Description:   This subrotine loads a previous fieldlines run and
!                    generates a set of fieldline starting points
!                    from a given field line trace.
!-----------------------------------------------------------------------
      SUBROUTINE fieldlines_init_fline
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE fieldlines_lines, ONLY: nlines
      USE fieldlines_runtime
      USE ez_hdf5
      USE mpi_params
      USE mpi_inc
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, i
      INTEGER :: npoinc_old, nlines_old
      REAL(rprec) :: dphi_old
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: R_old, Z_old, PHI_old

!      INTEGER, PARAMETER :: line_select = 64
      INTEGER, PARAMETER :: nnew_lines = 2**17
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      IF (lverb) THEN
         WRITE(6,'(A)')  '----- Generating Starting Points -----'
         WRITE(6,'(A)')  '   FILE: '//TRIM(restart_string)
         WRITE(6,'(A,I6)') '   NEW_LINES: ',nnew_lines*2
         WRITE(6,'(A,I3)') '      LINE #: ',line_select
      END IF
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
#endif
      ! Setup PHI_END before going on
      i = MAXLOC(ABS(PHI_end),DIM=1)
      PHI_end(1:nnew_lines) = PHI_end(i)
      ! Initialize fieldine data
      r_start   = -1
      z_start   = -1
      phi_start =  0
      phi_end(nnew_lines+1:MAXLINES)   =  0
#if defined(LHDF5)
      IF (myworkid == master) THEN
         CALL open_hdf5(TRIM(restart_string),fid,ier,LCREATE=.false.)
         CALL read_scalar_hdf5(fid,'nsteps',ier,INTVAR=npoinc_old)
         CALL read_scalar_hdf5(fid,'nlines',ier,INTVAR=nlines_old)
         IF (lverb) THEN
            WRITE(6,'(A,I3)') '   OLD_LINES: ',nlines_old
            WRITE(6,'(A,I6)') '      NSTEPS: ',npoinc_old
         END IF
         ALLOCATE(R_old(nlines_old,npoinc_old+1), &
                  Z_old(nlines_old,npoinc_old+1), &
                  PHI_old(nlines_old,npoinc_old+1))
         CALL read_var_hdf5(fid,'R_lines',nlines_old,npoinc_old+1,ier,DBLVAR=R_old)
         CALL read_var_hdf5(fid,'Z_lines',nlines_old,npoinc_old+1,ier,DBLVAR=Z_old)
         CALL read_var_hdf5(fid,'PHI_lines',nlines_old,npoinc_old+1,ier,DBLVAR=PHI_old)
         CALL close_hdf5(fid,ier)
         npoinc_old = npoinc_old - 1
         ! Setup new PHI
         dphi_old = PHI_old(line_select,2) - PHI_old(line_select,1)
         CALL RANDOM_NUMBER(PHI_start(1:nnew_lines))
         PHI_start = PHI_start * PHI_old(line_select,npoinc_old)
         CALL spline_it(npoinc_old, &
                        PHI_old(line_select,1:npoinc_old), &
                        R_old(line_select,1:npoinc_old), &
                        nnew_lines, &
                        PHI_start(1:nnew_lines), &
                        R_start(1:nnew_lines), 0)
         CALL spline_it(npoinc_old, &
                        PHI_old(line_select,1:npoinc_old), &
                        Z_old(line_select,1:npoinc_old), &
                        nnew_lines, &
                        PHI_start(1:nnew_lines), &
                        Z_start(1:nnew_lines), 0)
         DEALLOCATE(R_old, PHI_old, Z_old)
         ! Now setup reverse fieldline trace
         i = nnew_lines+1
         PHI_start(i:2*nnew_lines) = PHI_start(1:nnew_lines)
         R_start(i:2*nnew_lines) = R_start(1:nnew_lines)
         Z_start(i:2*nnew_lines) = Z_start(1:nnew_lines)
         PHI_end(i:2*nnew_lines) = -PHI_end(1:nnew_lines)
         nlines = 2*nnew_lines
      END IF
#endif
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_BCAST(nlines,           1, MPI_INTEGER, master, MPI_COMM_FIELDLINES, ierr_mpi)
      CALL MPI_BCAST(r_start,   MAXLINES,   MPI_REAL8, master, MPI_COMM_FIELDLINES, ierr_mpi)
      CALL MPI_BCAST(z_start,   MAXLINES,   MPI_REAL8, master, MPI_COMM_FIELDLINES, ierr_mpi)
      CALL MPI_BCAST(phi_start, MAXLINES,   MPI_REAL8, master, MPI_COMM_FIELDLINES, ierr_mpi)
      CALL MPI_BCAST(phi_end,   MAXLINES,   MPI_REAL8, master, MPI_COMM_FIELDLINES, ierr_mpi)
#endif
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE fieldlines_init_fline
