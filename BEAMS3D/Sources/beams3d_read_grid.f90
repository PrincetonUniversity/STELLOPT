!-----------------------------------------------------------------------
!     Module:        beams3d_read_grid
!     Authors:       D. Kulla (david.kulla@ipp.mpg.de)
!     Date:          30/05/2023
!     Description:   This subroutine reads a BEAMS3D grid into memory.
!-----------------------------------------------------------------------
SUBROUTINE beams3d_read_grid(file_ext)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
   USE stel_kinds, ONLY: rprec
#if defined(LHDF5)
   USE ez_hdf5
#endif
   USE beams3d_lines
   USE beams3d_grid
   USE beams3d_runtime
   USE wall_mod
   USE safe_open_mod, ONLY: safe_open
   USE wall_mod, ONLY: nface,nvertex,face,vertex,ihit_array
   USE mpi_sharmem
   USE mpi_params
   USE mpi_inc
!-----------------------------------------------------------------------
!     Input Variables
!          file_ext     Extension of file (beams3d_ext.h5)
!-----------------------------------------------------------------------
   IMPLICIT NONE
   CHARACTER(LEN=*), INTENT(in)           :: file_ext
!-----------------------------------------------------------------------
!     Local Variables
!          ier          Error Flag
!          iunit        File ID
!-----------------------------------------------------------------------
   INTEGER :: ier, iunit
   REAL(rprec) :: ver_temp, t_end_restart
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
   IF (lverb) THEN
      WRITE(6,'(A)')  '----- READING GRID FROM FILE -----'
   END IF
#if defined(LHDF5)

   IF (myid_sharmem == master) THEN
      IF (lverb) WRITE(6,'(A)')  '   FILE: '//'beams3d_'//TRIM(file_ext)//'.h5'
      CALL open_hdf5('beams3d_'//TRIM(file_ext)//'.h5',fid,ier,LCREATE=.false.)
      IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,'beams3d_'//TRIM(file_ext)//'.h5',ier)
      CALL read_scalar_hdf5(fid,'nr',ier,INTVAR=nr)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'nr',ier)
      CALL read_scalar_hdf5(fid,'nphi',ier,INTVAR=nphi)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'nphi',ier)
      CALL read_scalar_hdf5(fid,'nz',ier,INTVAR=nz)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'nz',ier)
   END IF
   ! Grid
   IF (ASSOCIATED(raxis))   DEALLOCATE(raxis)
   IF (ASSOCIATED(zaxis))   DEALLOCATE(zaxis)
   IF (ASSOCIATED(phiaxis)) DEALLOCATE(phiaxis)
   IF (ASSOCIATED(B_R))     DEALLOCATE(B_R)
   IF (ASSOCIATED(B_Z))     DEALLOCATE(B_Z)
   IF (ASSOCIATED(B_PHI))   DEALLOCATE(B_PHI)
   IF (ASSOCIATED(S_ARR))   DEALLOCATE(S_ARR)
   IF (ASSOCIATED(U_ARR))   DEALLOCATE(U_ARR)
   IF (ASSOCIATED(POT_ARR)) DEALLOCATE(POT_ARR)
   IF (ASSOCIATED(TE))     DEALLOCATE(TE)
   IF (ASSOCIATED(NE))     DEALLOCATE(NE)
   IF (ASSOCIATED(TI))     DEALLOCATE(TI)
   IF (ASSOCIATED(NI))     DEALLOCATE(NI)
   IF (ASSOCIATED(ZEFF_ARR))     DEALLOCATE(ZEFF_ARR)
   CALL MPI_BARRIER(MPI_COMM_SHARMEM, ier)
   CALL mpialloc(raxis, nr, myid_sharmem, 0, MPI_COMM_SHARMEM, win_raxis)
   CALL mpialloc(phiaxis, nphi, myid_sharmem, 0, MPI_COMM_SHARMEM, win_phiaxis)
   CALL mpialloc(zaxis, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_zaxis)
   CALL mpialloc(B_R, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_B_R)
   CALL mpialloc(B_PHI, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_B_PHI)
   CALL mpialloc(B_Z, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_B_Z)
   CALL mpialloc(TE, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_TE)
   CALL mpialloc(NE, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_NE)
   CALL mpialloc(TI, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_TI)
   CALL mpialloc(NI, NION, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_NI)
   CALL mpialloc(ZEFF_ARR, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_ZEFF_ARR)
   CALL mpialloc(POT_ARR, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_POT_ARR)
   CALL mpialloc(S_ARR, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_S_ARR)
   CALL mpialloc(U_ARR, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_U_ARR)
   IF (myid_sharmem == master) THEN
      ! ALLOCATE(raxis(nr))
      ! ALLOCATE(phiaxis(nphi))
      ! ALLOCATE(zaxis(nz))
      ! ALLOCATE(B_R(nr,nphi,nz))
      ! ALLOCATE(B_PHI(nr,nphi,nz))
      ! ALLOCATE(B_Z(nr,nphi,nz))
      ! ALLOCATE(MODB(nr,nphi,nz))
      ! ALLOCATE(POT_ARR(nr,nphi,nz))
      ! ALLOCATE(S_ARR(nr,nphi,nz))
      ! ALLOCATE(U_ARR(nr,nphi,nz))
      CALL read_var_hdf5(fid,'raxis',nr,ier,DBLVAR=raxis)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'raxis',ier)
      CALL read_var_hdf5(fid,'zaxis',nz,ier,DBLVAR=zaxis)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'zaxis',ier)
      CALL read_var_hdf5(fid,'phiaxis',nphi,ier,DBLVAR=phiaxis)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'phiaxis',ier)
      CALL read_var_hdf5(fid,'B_R',nr,nphi,nz,ier,DBLVAR=B_R)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'B_R',ier)
      CALL read_var_hdf5(fid,'B_PHI',nr,nphi,nz,ier,DBLVAR=B_PHI)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'B_PHI',ier)
      CALL read_var_hdf5(fid,'B_Z',nr,nphi,nz,ier,DBLVAR=B_Z)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'B_Z',ier)
      CALL read_var_hdf5(fid,'S_ARR',nr,nphi,nz,ier,DBLVAR=S_ARR)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'S_ARR',ier)
      CALL read_var_hdf5(fid,'U_ARR',nr,nphi,nz,ier,DBLVAR=U_ARR)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'U_ARR',ier)
      CALL read_var_hdf5(fid,'POT_ARR',nr,nphi,nz,ier,DBLVAR=POT_ARR)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'POT_ARR',ier)

      ! ALLOCATE(TE(nr,nphi,nz))
      ! ALLOCATE(NE(nr,nphi,nz))
      ! ALLOCATE(TI(nr,nphi,nz))
      ! ALLOCATE(NI(NION,nr,nphi,nz))
      ! ALLOCATE(ZEFF_ARR(nr,nphi,nz))
      rmax = MAXVAL(raxis)
      rmin = MINVAL(raxis)
      phimax = MAXVAL(phiaxis)
      phimin = MINVAL(phiaxis)
      zmax = MAXVAL(zaxis)
      zmin = MINVAL(zaxis)
      CALL read_var_hdf5(fid,'TE',nr,nphi,nz,ier,DBLVAR=TE)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'TE',ier)
      CALL read_var_hdf5(fid,'NE',nr,nphi,nz,ier,DBLVAR=NE)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'NE',ier)
      WRITE (27, '(EN12.3)') MAXVAL(NE)
      CALL read_var_hdf5(fid,'TI',nr,nphi,nz,ier,DBLVAR=TI)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'TI',ier)
      CALL read_var_hdf5(fid,'NI',nion,nr,nphi,nz,ier,DBLVAR=NI)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'NI',ier)
      CALL read_var_hdf5(fid,'ZEFF_ARR',nr,nphi,nz,ier,DBLVAR=ZEFF_ARR)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'ZEFF_ARR',ier)
      ! CALL read_var_hdf5(fid,'dist_rhoaxis',ns_prof1,ier,DBLVAR=dist_rhoaxis)
      ! IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'dist_rhoaxis',ier)
      ! CALL read_var_hdf5(fid,'dist_uaxis',ns_prof2,ier,DBLVAR=dist_uaxis)
      ! IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'dist_uaxis',ier)
      ! CALL read_var_hdf5(fid,'dist_paxis',ns_prof2,ier,DBLVAR=dist_paxis)
      ! IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'dist_paxis',ier)
      ! CALL read_var_hdf5(fid,'dist_Vaxis',ns_prof4,DBLVAR=dist_Vaxis)
      ! IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'dist_Vaxis',ier)
      ! CALL read_var_hdf5(fid,'dist_Waxis',ns_prof5,DBLVAR=dist_Waxis)
      ! IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'dist_Waxis',ier)

      ! ! Try to read the faces
      ! CALL read_scalar_hdf5(fid,'nvertex',ier,INTVAR=nvertex)
      ! IF (ier == 0) THEN
      !    ALLOCATE(vertex(nvertex,3))
      !    CALL read_var_hdf5(fid,'wall_vertex',nvertex,3,ier,DBLVAR=vertex)
      !    IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'wall_vertex',ier)
      ! END IF
      ! ier = 0
      ! CALL read_scalar_hdf5(fid,'nface',ier,INTVAR=nface)
      ! IF (ier == 0) THEN
      !    ALLOCATE(face(nface,3),ihit_array(nface))
      !    CALL read_var_hdf5(fid,'wall_faces',nface,3,ier,INTVAR=face)
      !    IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'wall_faces',ier)
      !    CALL read_var_hdf5(fid,'wall_strikes',nface,ier,INTVAR=ihit_array)
      !    !IF (ier /= 0) DEALLOCATE(ihit_array)
      !    !lwall_loaded=.true.
      ! END IF
      ! ier = 0
      ! Close the file
      CALL close_hdf5(fid,ier)
      IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,'beams3d_'//TRIM(file_ext)//'.h5',ier)
   END IF

#else
   ! To be done
   IF (lverb) WRITE(6,*) 'ERROR: Reading from non-HDF5 not implemented!'
#endif

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
END SUBROUTINE beams3d_read_grid
