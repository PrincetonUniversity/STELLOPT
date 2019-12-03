!-----------------------------------------------------------------------
!     Module:        beams3d_read
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          01/09/2014
!     Description:   This subroutine reads a BEAMS3D run into memory.
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_read
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
!DEC$ IF DEFINED (LHDF5)
      USE ez_hdf5
!DEC$ ENDIF  
      USE beams3d_lines
      USE beams3d_grid
      USE beams3d_runtime
      USE safe_open_mod, ONLY: safe_open
      USE wall_mod, ONLY: nface,nvertex,face,vertex,ihit_array
      USE mpi_sharmem
!-----------------------------------------------------------------------
!     Local Variables
!          ier          Error Flag
!          iunit        File ID
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, iunit
      REAL(rprec) :: ver_temp
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      IF (lverb) THEN
         WRITE(6,'(A)')  '----- READING DATA FROM FILE -----'
      END IF
!DEC$ IF DEFINED (LHDF5)
      IF (lverb) WRITE(6,'(A)')  '   FILE: '//'beams3d_'//TRIM(id_string)//'.h5'
      CALL open_hdf5('beams3d_'//TRIM(id_string)//'.h5',fid,ier,LCREATE=.false.)
      IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,'beams3d_'//TRIM(id_string)//'.h5',ier)
      ! Runtime
      CALL read_scalar_hdf5(fid,'VERSION',ier,DBLVAR=ver_temp)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'VERSION',ier)
      CALL read_scalar_hdf5(fid,'lvmec',ier,BOOVAR=lvmec)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'lvmec',ier)
      CALL read_scalar_hdf5(fid,'lpies',ier,BOOVAR=lpies)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'lpies',ier)
      CALL read_scalar_hdf5(fid,'lspec',ier,BOOVAR=lspec)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'lspec',ier)
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
      CALL read_scalar_hdf5(fid,'lflux',ier,BOOVAR=lflux)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'lflux',ier)
      CALL read_scalar_hdf5(fid,'ldepo',ier,BOOVAR=ldepo)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'ldepo',ier)
      CALL read_scalar_hdf5(fid,'lbeam',ier,BOOVAR=lbeam)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'lbeam',ier)
      ! Trajectories
      CALL read_scalar_hdf5(fid,'nparticles',ier,INTVAR=nparticles)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'nparticles',ier)
      CALL read_scalar_hdf5(fid,'nbeams',ier,INTVAR=nbeams)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'nbeams',ier)
      CALL read_scalar_hdf5(fid,'nsteps',ier,INTVAR=nsteps)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'nsteps',ier)
      CALL read_scalar_hdf5(fid,'npoinc',ier,INTVAR=npoinc)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'npoinc',ier)
      IF (ALLOCATED(t_end)) DEALLOCATE(t_end)
      IF (ALLOCATED(mass)) DEALLOCATE(mass)
      IF (ALLOCATED(charge)) DEALLOCATE(charge)
      IF (ALLOCATED(Zatom)) DEALLOCATE(Zatom)
      ALLOCATE(t_end(nparticles),mass(nparticles),charge(nparticles),Zatom(nparticles))
      CALL read_var_hdf5(fid,'t_end',nparticles,ier,DBLVAR=t_end)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'t_end',ier)
      CALL read_var_hdf5(fid,'mass',nparticles,ier,DBLVAR=mass)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'mass',ier)
      CALL read_var_hdf5(fid,'charge',nparticles,ier,DBLVAR=charge)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'charge',ier)
      CALL read_var_hdf5(fid,'Zatom',nparticles,ier,DBLVAR=Zatom)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'Zatom',ier)
      IF (ALLOCATED(R_lines)) DEALLOCATE(R_lines)
      IF (ALLOCATED(Z_lines)) DEALLOCATE(Z_lines)
      IF (ALLOCATED(PHI_lines)) DEALLOCATE(PHI_lines)
      IF (ALLOCATED(vll_lines)) DEALLOCATE(vll_lines)
      IF (ALLOCATED(neut_lines)) DEALLOCATE(neut_lines)
      IF (ALLOCATED(moment_lines)) DEALLOCATE(moment_lines)
      IF (ALLOCATED(S_lines)) DEALLOCATE(S_lines)
      IF (ALLOCATED(U_lines)) DEALLOCATE(U_lines)
      ALLOCATE(R_lines(0:npoinc,nparticles),Z_lines(0:npoinc,nparticles),PHI_lines(0:npoinc,nparticles),&
            vll_lines(0:npoinc,nparticles),neut_lines(0:npoinc,nparticles),moment_lines(0:npoinc,nparticles))
      ALLOCATE(S_lines(0:npoinc,nparticles),U_lines(0:npoinc,nparticles),B_lines(0:npoinc,nparticles))
      ALLOCATE(PE_lines(0:npoinc,nparticles),PI_lines(0:npoinc,nparticles))
      CALL read_var_hdf5(fid,'R_lines',npoinc+1,nparticles,ier,DBLVAR=R_lines)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'R_lines',ier)
      CALL read_var_hdf5(fid,'Z_lines',npoinc+1,nparticles,ier,DBLVAR=Z_lines)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'Z_lines',ier)
      CALL read_var_hdf5(fid,'PHI_lines',npoinc+1,nparticles,ier,DBLVAR=PHI_lines)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'PHI_lines',ier)
      CALL read_var_hdf5(fid,'vll_lines',npoinc+1,nparticles,ier,DBLVAR=vll_lines)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'vll_lines',ier)
      CALL read_var_hdf5(fid,'neut_lines',npoinc+1,nparticles,ier,BOOVAR=neut_lines)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'neut_lines',ier)
      CALL read_var_hdf5(fid,'moment_lines',npoinc+1,nparticles,ier,DBLVAR=moment_lines)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'moment_lines',ier)
      CALL read_var_hdf5(fid,'S_lines',npoinc+1,nparticles,ier,DBLVAR=S_lines)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'S_lines',ier)
      CALL read_var_hdf5(fid,'U_lines',npoinc+1,nparticles,ier,DBLVAR=U_lines)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'U_lines',ier)
      CALL read_var_hdf5(fid,'B_lines',npoinc+1,nparticles,ier,DBLVAR=B_lines)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'B_lines',ier)
      CALL read_var_hdf5(fid,'PE_lines',npoinc+1,nparticles,ier,DBLVAR=PE_lines)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'PE_lines',ier)
      CALL read_var_hdf5(fid,'PI_lines',npoinc+1,nparticles,ier,DBLVAR=PI_lines)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'PI_lines',ier)
      IF (lflux) THEN
         ! DIAGNOSTICS
         IF (ALLOCATED(shine_through)) DEALLOCATE(shine_through)
         ALLOCATE(shine_through(nbeams))
         CALL read_var_hdf5(fid,'Shinethrough',nbeams,ier,DBLVAR=shine_through)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'shine_through',ier)
         IF (lbeam .and. .not.ldepo) THEN
            IF (ALLOCATED(ndot_prof)) DEALLOCATE(ndot_prof)
            IF (ALLOCATED(epower_prof)) DEALLOCATE(epower_prof)
            IF (ALLOCATED(ipower_prof)) DEALLOCATE(ipower_prof)
            IF (ALLOCATED(j_prof)) DEALLOCATE(j_prof)
            CALL read_scalar_hdf5(fid,'ns_prof',ier,INTVAR=ns_prof)
            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'ns_prof',ier)
            CALL read_var_hdf5(fid,'ndot_prof',nbeams,ns_prof,ier,DBLVAR=ndot_prof)
            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'ndot_prof',ier)
            CALL read_var_hdf5(fid,'epower_prof',nbeams,ns_prof,ier,DBLVAR=epower_prof)
            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'epower_prof',ier)
            CALL read_var_hdf5(fid,'ipower_prof',nbeams,ns_prof,ier,DBLVAR=ipower_prof)
            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'ipower_prof',ier)
            CALL read_var_hdf5(fid,'j_prof',nbeams,ns_prof,ier,DBLVAR=j_prof)
            IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'j_prof',ier)
         END IF
      END IF
      IF (lbeam) THEN
         IF (ALLOCATED(weight)) DEALLOCATE(weight)
         IF (ALLOCATED(beam)) DEALLOCATE(beam)
         IF (ALLOCATED(v_neut)) DEALLOCATE(v_neut)
         ALLOCATE(weight(nparticles_start,nbeams),beam(nparticles),v_neut(3,nparticles))
         CALL read_var_hdf5(fid,'Weight',nparticles_start,nbeams,ier,DBLVAR=weight)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'weight',ier)
         CALL read_var_hdf5(fid,'Beam',nparticles,ier,INTVAR=beam)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'beam',ier)
         CALL read_var_hdf5(fid,'Energy',nbeams,ier,DBLVAR=e_beams)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'e_beams',ier)
         CALL read_var_hdf5(fid,'V_NEUT',3,nparticles,ier,DBLVAR=v_neut)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'e_beams',ier)
      END IF
      ! Grid
      CALL read_scalar_hdf5(fid,'nr',ier,INTVAR=nr)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'nr',ier)
      CALL read_scalar_hdf5(fid,'nphi',ier,INTVAR=nphi)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'nphi',ier)
      CALL read_scalar_hdf5(fid,'nz',ier,INTVAR=nz)
      IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'nz',ier)
      IF (ASSOCIATED(raxis))   DEALLOCATE(raxis)
      IF (ASSOCIATED(zaxis))   DEALLOCATE(zaxis)
      IF (ASSOCIATED(phiaxis)) DEALLOCATE(phiaxis)
      IF (ASSOCIATED(B_R))     DEALLOCATE(B_R)
      IF (ASSOCIATED(B_Z))     DEALLOCATE(B_Z)
      IF (ASSOCIATED(B_PHI))   DEALLOCATE(B_PHI)
      IF (ASSOCIATED(S_ARR))   DEALLOCATE(S_ARR)
      IF (ASSOCIATED(U_ARR))   DEALLOCATE(U_ARR)
      IF (ASSOCIATED(POT_ARR)) DEALLOCATE(POT_ARR)
      ALLOCATE(raxis(nr))
      ALLOCATE(phiaxis(nphi))
      ALLOCATE(zaxis(nz))
      ALLOCATE(B_R(nr,nphi,nz))
      ALLOCATE(B_PHI(nr,nphi,nz))
      ALLOCATE(B_Z(nr,nphi,nz))
      ALLOCATE(MODB(nr,nphi,nz))
      ALLOCATE(POT_ARR(nr,nphi,nz))
      ALLOCATE(S_ARR(nr,nphi,nz))
      ALLOCATE(U_ARR(nr,nphi,nz))
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
      ! Try to read the faces
      CALL read_scalar_hdf5(fid,'nvertex',ier,INTVAR=nvertex)
      IF (ier == 0) THEN
         ALLOCATE(vertex(nvertex,3))
         CALL read_var_hdf5(fid,'wall_vertex',nvertex,3,ier,DBLVAR=vertex)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'wall_vertex',ier)
      END IF
      ier = 0
      CALL read_scalar_hdf5(fid,'nface',ier,INTVAR=nface)
      IF (ier == 0) THEN
         ALLOCATE(face(nface,3),ihit_array(nface))
         CALL read_var_hdf5(fid,'wall_faces',nface,3,ier,INTVAR=face)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'wall_faces',ier)
         CALL read_var_hdf5(fid,'wall_strikes',nface,ier,INTVAR=ihit_array)
         !IF (ier /= 0) DEALLOCATE(ihit_array)
      END IF
      ier = 0

      ! Close the file
      CALL close_hdf5(fid,ier)
      IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,'beams3d_'//TRIM(id_string)//'.h5',ier)




!DEC$ ELSE

!DEC$ ENDIF  

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE beams3d_read
