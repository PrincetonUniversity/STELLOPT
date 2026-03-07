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
   USE mpi_sharmem
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
      CALL close_hdf5(fid,ier)
   end if
#endif
#if defined(MPI_OPT)
   CALL MPI_BCAST(nr,   1, MPI_INTEGER, master, MPI_COMM_FIELDLINES, ierr_mpi)
   CALL MPI_BCAST(nphi, 1, MPI_INTEGER, master, MPI_COMM_FIELDLINES, ierr_mpi)
   CALL MPI_BCAST(nz,   1, MPI_INTEGER, master, MPI_COMM_FIELDLINES, ierr_mpi)
   CALL mpialloc(raxis,   nr,   myid_sharmem, master, MPI_COMM_SHARMEM, win_raxis)
   CALL mpialloc(phiaxis, nphi, myid_sharmem, master, MPI_COMM_SHARMEM, win_phiaxis)
   CALL mpialloc(zaxis,   nz,   myid_sharmem, master, MPI_COMM_SHARMEM, win_zaxis)

   CALL mpialloc(B_R,   nr, nphi, nz, myid_sharmem, master, MPI_COMM_SHARMEM, win_B_R)
   CALL mpialloc(B_PHI, nr, nphi, nz, myid_sharmem, master, MPI_COMM_SHARMEM, win_B_PHI)
   CALL mpialloc(B_Z,   nr, nphi, nz, myid_sharmem, master, MPI_COMM_SHARMEM, win_B_Z)
#endif
#if defined(LHDF5)
   IF (myworkid == master) THEN
      CALL open_hdf5(TRIM(restart_string),fid,ier)
      CALL read_var_hdf5(fid,'raxis',nr,ier,DBLVAR=raxis)
      CALL read_var_hdf5(fid,'phiaxis',nphi,ier,DBLVAR=phiaxis)
      CALL read_var_hdf5(fid,'zaxis',nz,ier,DBLVAR=zaxis)
      CALL read_var_hdf5(fid,'B_R',nr,nphi,nz,ier,DBLVAR=B_R)
      CALL read_var_hdf5(fid,'B_PHI',nr,nphi,nz,ier,DBLVAR=B_PHI)
      CALL read_var_hdf5(fid,'B_Z',nr,nphi,nz,ier,DBLVAR=B_Z)
      CALL close_hdf5(fid,ier)
   END IF
#endif

#if defined(MPI_OPT)
   CALL MPI_WIN_SYNC(win_raxis, ierr_mpi)
   CALL MPI_WIN_SYNC(win_phiaxis, ierr_mpi)
   CALL MPI_WIN_SYNC(win_zaxis, ierr_mpi)
   CALL MPI_WIN_SYNC(win_B_R, ierr_mpi)
   CALL MPI_WIN_SYNC(win_B_PHI, ierr_mpi)
   CALL MPI_WIN_SYNC(win_B_Z, ierr_mpi)

   CALL MPI_BARRIER(MPI_COMM_SHARMEM, ierr_mpi)
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
end subroutine fieldlines_init_restart
