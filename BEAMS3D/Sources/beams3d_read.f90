!-----------------------------------------------------------------------
!     Module:        beams3d_read
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          01/09/2014
!     Description:   This subroutine reads a BEAMS3D run into memory.
!-----------------------------------------------------------------------
SUBROUTINE beams3d_read(file_ext)
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
      WRITE(6,'(A)')  '----- READING DATA FROM FILE -----'
   END IF
#if defined(LHDF5)
   IF (lverb) WRITE(6,'(A)')  '   FILE: '//'beams3d_'//TRIM(file_ext)//'.h5'
   CALL open_hdf5('beams3d_'//TRIM(file_ext)//'.h5',fid,ier,LCREATE=.false.)
   IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,'beams3d_'//TRIM(file_ext)//'.h5',ier)
   ! Runtime
   CALL read_scalar_hdf5(fid,'VERSION',ier,DBLVAR=ver_temp)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'VERSION',ier)
   CALL read_scalar_hdf5(fid,'lvmec',ier,BOOVAR=lvmec)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'lvmec',ier)
   CALL read_scalar_hdf5(fid,'lpies',ier,BOOVAR=lpies)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'lpies',ier)
   CALL read_scalar_hdf5(fid,'lspec',ier,BOOVAR=lspec)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'lspec',ier)
   CALL read_scalar_hdf5(fid,'leqdsk',ier,BOOVAR=leqdsk)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'leqdsk',ier)
   CALL read_scalar_hdf5(fid,'lcoil',ier,BOOVAR=lcoil)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'lcoil',ier)
   CALL read_scalar_hdf5(fid,'lmgrid',ier,BOOVAR=lmgrid)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'lmgrid',ier)
   CALL read_scalar_hdf5(fid,'lvessel',ier,BOOVAR=lvessel)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'lvessel',ier)
   CALL read_scalar_hdf5(fid,'lvac',ier,BOOVAR=lvac)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'lvac',ier)
   CALL read_scalar_hdf5(fid,'lbeam_simple',ier,BOOVAR=lbeam_simple)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'lbeam_simple',ier)
   CALL read_scalar_hdf5(fid,'ldepo',ier,BOOVAR=ldepo)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'ldepo',ier)
   CALL read_scalar_hdf5(fid,'lbeam',ier,BOOVAR=lbeam)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'lbeam',ier)
   CALL read_scalar_hdf5(fid,'lcollision',ier,BOOVAR=lcollision)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'lcollision',ier)
   CALL read_scalar_hdf5(fid,'lascot',ier,BOOVAR=lascot)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'lascot',ier)
   CALL read_scalar_hdf5(fid,'lfusion',ier,BOOVAR=lfusion)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'lfusion',ier)
   CALL read_scalar_hdf5(fid,'lhitonly',ier,BOOVAR=lhitonly)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'lhitonly',ier)
   ! Trajectories
   CALL read_scalar_hdf5(fid,'nparticles',ier,INTVAR=nparticles)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'nparticles',ier)
   CALL read_scalar_hdf5(fid,'nbeams',ier,INTVAR=nbeams)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'nbeams',ier)
   CALL read_scalar_hdf5(fid,'nsteps',ier,INTVAR=nsteps)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'nsteps',ier)
   CALL read_scalar_hdf5(fid,'npoinc',ier,INTVAR=npoinc)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'npoinc',ier)
   CALL read_scalar_hdf5(fid,'nbeams',ier,INTVAR=nbeams)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'nbeams',ier)
   CALL read_scalar_hdf5(fid,'ns_prof1',ier,INTVAR=ns_prof1)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'ns_prof1',ier)
   CALL read_scalar_hdf5(fid,'ns_prof2',ier,INTVAR=ns_prof2)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'ns_prof2',ier)
   CALL read_scalar_hdf5(fid,'ns_prof3',ier,INTVAR=ns_prof3)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'ns_prof3',ier)
   CALL read_scalar_hdf5(fid,'ns_prof4',ier,INTVAR=ns_prof4)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'ns_prof4',ier)
   CALL read_scalar_hdf5(fid,'ns_prof5',ier,INTVAR=ns_prof5)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'ns_prof5',ier)
   CALL read_scalar_hdf5(fid,'t_end_in',ier,DBLVAR=t_end_restart)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'t_end',ier)
   t_end_in = t_end_restart !Broadcast to t_end_in
   IF (ALLOCATED(t_end)) DEALLOCATE(t_end)
   IF (ALLOCATED(mass)) DEALLOCATE(mass)
   IF (ALLOCATED(charge)) DEALLOCATE(charge)
   IF (ALLOCATED(Zatom)) DEALLOCATE(Zatom)
   IF (ALLOCATED(Weight)) DEALLOCATE(Weight)
   IF (ALLOCATED(Beam)) DEALLOCATE(Beam)
   IF (ALLOCATED(end_state)) DEALLOCATE(end_state)
   ALLOCATE(t_end(nparticles),mass(nparticles),charge(nparticles), &
      Zatom(nparticles),end_state(nparticles), Weight(nparticles), &
      Beam(nparticles))
   CALL read_var_hdf5(fid,'t_end',nparticles,ier,DBLVAR=t_end)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'t_end',ier)
   CALL read_var_hdf5(fid,'mass',nparticles,ier,DBLVAR=mass)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'mass',ier)
   CALL read_var_hdf5(fid,'charge',nparticles,ier,DBLVAR=charge)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'charge',ier)
   CALL read_var_hdf5(fid,'Zatom',nparticles,ier,DBLVAR=Zatom)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'Zatom',ier)
   CALL read_var_hdf5(fid,'end_state',nparticles,ier,INTVAR=end_state)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'end_state',ier)
   IF (ALLOCATED(R_lines)) DEALLOCATE(R_lines)
   IF (ALLOCATED(Z_lines)) DEALLOCATE(Z_lines)
   IF (ALLOCATED(PHI_lines)) DEALLOCATE(PHI_lines)
   IF (ALLOCATED(vll_lines)) DEALLOCATE(vll_lines)
   IF (ALLOCATED(neut_lines)) DEALLOCATE(neut_lines)
   IF (ALLOCATED(moment_lines)) DEALLOCATE(moment_lines)
   IF (ALLOCATED(v_neut)) DEALLOCATE(v_neut)
   IF (ALLOCATED(S_lines)) DEALLOCATE(S_lines)
   IF (ALLOCATED(U_lines)) DEALLOCATE(U_lines)
   IF (ALLOCATED(B_lines)) DEALLOCATE(B_lines)
   ALLOCATE(R_lines(0:npoinc,nparticles),Z_lines(0:npoinc,nparticles),PHI_lines(0:npoinc,nparticles),&
      vll_lines(0:npoinc,nparticles),neut_lines(0:npoinc,nparticles),moment_lines(0:npoinc,nparticles),&
      v_neut(3,nparticles))
   ALLOCATE(S_lines(0:npoinc,nparticles),U_lines(0:npoinc,nparticles),B_lines(0:npoinc,nparticles))
   CALL read_var_hdf5(fid,'R_lines',npoinc+1,nparticles,ier,DBLVAR=R_lines)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'R_lines',ier)
   CALL read_var_hdf5(fid,'Z_lines',npoinc+1,nparticles,ier,DBLVAR=Z_lines)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'Z_lines',ier)
   CALL read_var_hdf5(fid,'PHI_lines',npoinc+1,nparticles,ier,DBLVAR=PHI_lines)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'PHI_lines',ier)
   CALL read_var_hdf5(fid,'vll_lines',npoinc+1,nparticles,ier,DBLVAR=vll_lines)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'vll_lines',ier)
   CALL read_var_hdf5(fid,'moment_lines',npoinc+1,nparticles,ier,DBLVAR=moment_lines)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'moment_lines',ier)
   CALL read_var_hdf5(fid,'neut_lines',npoinc+1,nparticles,ier,BOOVAR=neut_lines)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'neut_lines',ier)
   CALL read_var_hdf5(fid,'S_lines',npoinc+1,nparticles,ier,DBLVAR=S_lines)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'S_lines',ier)
   CALL read_var_hdf5(fid,'U_lines',npoinc+1,nparticles,ier,DBLVAR=U_lines)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'U_lines',ier)
   CALL read_var_hdf5(fid,'B_lines',npoinc+1,nparticles,ier,DBLVAR=B_lines)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'B_lines',ier)
   CALL read_var_hdf5(fid,'V_NEUT',3,nparticles,ier,DBLVAR=v_neut)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'v_neut',ier)
   ! Particle parameters
   CALL read_var_hdf5(fid,'Weight',nparticles,ier,DBLVAR=weight)
   IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'weight',ier)
   IF (lbeam) THEN
      IF (ALLOCATED(beam)) DEALLOCATE(beam)
      IF (ALLOCATED(v_neut)) DEALLOCATE(v_neut)
      IF (ALLOCATED(shine_through)) DEALLOCATE(shine_through)
      ALLOCATE(shine_through(nbeams))
      ALLOCATE(beam(nparticles),v_neut(3,nparticles))
      CALL read_var_hdf5(fid,'Shinethrough',nbeams,ier,DBLVAR=shine_through)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'shine_through',ier)
      CALL read_var_hdf5(fid,'Beam',nparticles,ier,INTVAR=beam)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'beam',ier)
      CALL read_var_hdf5(fid,'Energy',nbeams,ier,DBLVAR=e_beams)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'e_beams',ier)
      CALL read_var_hdf5(fid,'V_NEUT',3,nparticles,ier,DBLVAR=v_neut)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'V_NEUT',ier)
      IF (.not. ldepo) THEN
!            Need to discuss how this is used
!            IF (ALLOCATED(ndot_prof)) DEALLOCATE(ndot_prof)
!            IF (ALLOCATED(epower_prof)) DEALLOCATE(epower_prof)
!            IF (ALLOCATED(ipower_prof)) DEALLOCATE(ipower_prof)
!            IF (ALLOCATED(j_prof)) DEALLOCATE(j_prof)
!            !IF (ALLOCATED(dist_prof)) DEALLOCATE(dist_prof)
!            IF (ALLOCATED(dist2d_prof)) DEALLOCATE(dist2d_prof)
!            ALLOCATE(ndot_prof(nbeams,ns_prof1),epower_prof(nbeams,ns_prof1),&
!               ipower_prof(nbeams,ns_prof1),j_prof(nbeams,ns_prof1))
!            !ALLOCATE(dist_prof(nbeams,ns_prof1,ns_prof2,ns_prof3,ns_prof4,ns_prof5))
!            ALLOCATE(dist2d_prof(nbeams,ns_prof4,ns_prof5))
!            CALL read_var_hdf5(fid, 'ndot_prof',   nbeams, ns_prof1, ier, DBLVAR=ndot_prof)
!            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'ndot_prof',ier)
!            CALL read_var_hdf5(fid, 'epower_prof', nbeams, ns_prof1, ier, DBLVAR=epower_prof)
!            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'epower_prof',ier)
!            CALL read_var_hdf5(fid, 'ipower_prof', nbeams, ns_prof1, ier, DBLVAR=ipower_prof)
!            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'ipower_prof',ier)
!            CALL read_var_hdf5(fid, 'j_prof',      nbeams, ns_prof1, ier, DBLVAR=j_prof)
!            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'j_prof',ier)
!            !CALL read_var_hdf5(fid, 'dist_prof',   nbeams, ns_prof1, ns_prof2,ns_prof3,ns_prof4,ns_prof5,ier,DBLVAR=dist_prof)
!            !IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'dist_prof',ier)
!            CALL read_var_hdf5(fid, 'dist2d_prof',   nbeams, ns_prof4, ns_prof5, ier, DBLVAR=dist2d_prof)
!            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'dist2d_prof',ier)
      END IF
   END IF
   ! Close the file
   CALL close_hdf5(fid,ier)
   IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,'beams3d_'//TRIM(file_ext)//'.h5',ier)
#else
   ! To be done
   IF (lverb) WRITE(6,*) 'ERROR: Reading from non-HDF5 not implemented!'
#endif

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
END SUBROUTINE beams3d_read
