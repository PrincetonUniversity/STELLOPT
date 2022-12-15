!-----------------------------------------------------------------------
!     Module:        fieldlines_init_restart
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/27/2012
!     Description:   This subroutine loads the fields from a restart file.
!                    This is still under development.
!-----------------------------------------------------------------------
      SUBROUTINE fieldlines_init_restart
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE fieldlines_runtime
      USE fieldlines_grid
      USE fieldlines_lines, ONLY: nlines
      USE ez_hdf5
      USE mpi_params
      USE mpi_inc
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, i
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      IF (lverb) THEN
         WRITE(6,'(A)')  '----- Reading Restart File -----'
         WRITE(6,'(A)')  '   FILE: '//TRIM(restart_string)
      END IF
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
#endif
#if defined(LHDF5)
      IF (myworkid == master) THEN
      CALL open_hdf5(TRIM(restart_string),fid,ier)
      PRINT *,'ier',ier
      CALL read_scalar_hdf5(fid,'nr',ier,INTVAR=nr)
      CALL read_scalar_hdf5(fid,'nphi',ier,INTVAR=nphi)
      CALL read_scalar_hdf5(fid,'nz',ier,INTVAR=nz)
      ALLOCATE(raxis(nr),phiaxis(nphi),zaxis(nz))
      CALL read_var_hdf5(fid,'raxis',nr,ier,DBLVAR=raxis)
      CALL read_var_hdf5(fid,'phiaxis',nphi,ier,DBLVAR=phiaxis)
      CALL read_var_hdf5(fid,'zaxis',nz,ier,DBLVAR=zaxis)
      ALLOCATE(B_R(nr,nphi,nz),B_PHI(nr,nphi,nz),B_Z(nr,nphi,nz))
      CALL read_var_hdf5(fid,'B_R',nr,nphi,nz,ier,DBLVAR=B_R)
      CALL read_var_hdf5(fid,'B_PHI',nr,nphi,nz,ier,DBLVAR=B_PHI)
      CALL read_var_hdf5(fid,'B_Z',nr,nphi,nz,ier,DBLVAR=B_Z)
      CALL close_hdf5(fid,ier)
      END IF
#endif
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (myworkid /= master) THEN
           ! DEALLOCATE(raxis, phiaxis, zaxis, &
           ! B_R, B_PHI, B_Z)
            ALLOCATE(raxis(nr), phiaxis(nphi), zaxis(nz), &
                       B_R(nr,nphi,nz), B_PHI(nr,nphi,nz), B_Z(nr,nphi,nz))
      END IF
         CALL MPI_BCAST(nr,1,MPI_INTEGER, master, MPI_COMM_FIELDLINES,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_restart',ierr_mpi)
         CALL MPI_BCAST(nphi,1,MPI_INTEGER, master, MPI_COMM_FIELDLINES,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_restart',ierr_mpi)
         CALL MPI_BCAST(nz,1,MPI_INTEGER, master, MPI_COMM_FIELDLINES,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_restart',ierr_mpi)
         CALL MPI_BCAST(raxis,nr,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_restart',ierr_mpi)
         CALL MPI_BCAST(phiaxis,nphi,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_restart',ierr_mpi)
         CALL MPI_BCAST(zaxis,nz,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_restart',ierr_mpi)
         CALL MPI_BCAST(B_R,nr*nphi*nz,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_restart',ierr_mpi)
         CALL MPI_BCAST(B_PHI,nr*nphi*nz,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_restart',ierr_mpi)
         CALL MPI_BCAST(B_Z,nr*nphi*nz,MPI_REAL8, master, MPI_COMM_FIELDLINES,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_restart',ierr_mpi)
#endif
      rmin = MINVAL(raxis)
      rmax = MAXVAL(raxis)
      phimin = MINVAL(phiaxis)
      phimax = MAXVAL(phiaxis)
      zmin   = MINVAL(zaxis)
      zmax   = MAXVAL(zaxis)
      IF (lverb) THEN
         WRITE(6,'(A,F8.5,A,F8.5,A,I4)') '   R   = [',rmin,',',rmax,'];  NR:   ',nr
         WRITE(6,'(A,F8.5,A,F8.5,A,I4)') '   PHI = [',phimin,',',phimax,'];  NPHI: ',nphi
         WRITE(6,'(A,F8.5,A,F8.5,A,I4)') '   Z   = [',zmin,',',zmax,'];  NZ:   ',nz
         WRITE(6,'(A,I4)')               '   # of Fieldlines: ',nlines
         CALL FLUSH(6)
      END IF
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE fieldlines_init_restart
