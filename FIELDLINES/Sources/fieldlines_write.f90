!-----------------------------------------------------------------------
!     Module:        fieldlines_write
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/22/2012
!     Description:   This subroutine outputs the fieldline data to an
!                    HDF5 file or binary file.
!-----------------------------------------------------------------------
      SUBROUTINE fieldlines_write
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE ez_hdf5
      USE fieldlines_lines
      USE fieldlines_grid, ONLY: nr, nphi, nz, B_R, B_PHI, B_Z, raxis, &
                                 zaxis, phiaxis, n_qshep, x_qshep, &
                                 y_qshep, z_qshep
      USE fieldlines_runtime, ONLY: id_string, npoinc, lverb, lvmec, &
                                    lpies, lspec, lcoil, lmgrid, lmu, &
                                    lvessel, lvac, laxis_i, handle_err,&
                                    HDF5_OPEN_ERR,HDF5_WRITE_ERR,&
                                    HDF5_CLOSE_ERR, FIELDLINES_VERSION,&
                                    ladvanced, lbfield_only, lreverse,&
                                    lafield_only, lemc3, lmodb
      USE safe_open_mod, ONLY: safe_open
      USE wall_mod, ONLY: nface,nvertex,face,vertex,ihit_array
!-----------------------------------------------------------------------
!     Local Variables
!          ier          Error Flag
!          iunit        File ID
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, iunit
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      IF (lverb) THEN
         WRITE(6,'(A)')  '----- WRITING DATA TO FILE -----'
      END IF
!DEC$ IF DEFINED (LHDF5)
      WRITE(6,'(A)')  '   FILE: '//'fieldlines_'//TRIM(id_string)//'.h5'
      CALL open_hdf5('fieldlines_'//TRIM(id_string)//'.h5',fid,ier,LCREATE=.true.)
      IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,'fieldlines_'//TRIM(id_string)//'.h5',ier)
      ! Runtime
      CALL write_scalar_hdf5(fid,'VERSION',ier,DBLVAR=FIELDLINES_VERSION,ATT='Version Number',ATT_NAME='description')
      IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'VERSION',ier)
      CALL write_scalar_hdf5(fid,'lvmec',ier,BOOVAR=lvmec,ATT='VMEC input',ATT_NAME='description')
      IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lvmec',ier)
      CALL write_scalar_hdf5(fid,'lpies',ier,BOOVAR=lpies,ATT='PIES input',ATT_NAME='description')
      IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lpies',ier)
      CALL write_scalar_hdf5(fid,'lspec',ier,BOOVAR=lspec,ATT='SPEC input',ATT_NAME='description')
      IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lspec',ier)
      CALL write_scalar_hdf5(fid,'lcoil',ier,BOOVAR=lcoil,ATT='Coil input',ATT_NAME='description')
      IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lcoil',ier)
      CALL write_scalar_hdf5(fid,'lmgrid',ier,BOOVAR=lmgrid,ATT='MGRID input',ATT_NAME='description')
      IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lmgrid',ier)
      CALL write_scalar_hdf5(fid,'lmu',ier,BOOVAR=lmu,ATT='Diffusion Logical',ATT_NAME='description')
      IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lmu',ier)
      CALL write_scalar_hdf5(fid,'lvessel',ier,BOOVAR=lvessel,ATT='Vessel input',ATT_NAME='description')
      IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lvessel',ier)
      CALL write_scalar_hdf5(fid,'lvac',ier,BOOVAR=lvac,ATT='Vacuum calc',ATT_NAME='description')
      IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lvac',ier)
      CALL write_scalar_hdf5(fid,'laxis_i',ier,BOOVAR=laxis_i,ATT='Axis calc',ATT_NAME='description')
      IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'laxis_i',ier)
      CALL write_scalar_hdf5(fid,'ladvanced',ier,BOOVAR=ladvanced,ATT='Advanced Grid Flag',ATT_NAME='description')
      IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'ladvanced',ier)
      CALL write_scalar_hdf5(fid,'lreverse',ier,BOOVAR=lreverse,ATT='VMEC input',ATT_NAME='description')
      IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lreverse',ier)
      ! Fieldlines
      IF (.not. lbfield_only .and. .not. lafield_only .and. .not. lemc3) THEN
         CALL write_scalar_hdf5(fid,'nlines',ier,INTVAR=nlines,ATT='Number of Fieldlines',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nlines',ier)
         CALL write_scalar_hdf5(fid,'nsteps',ier,INTVAR=nsteps+1,ATT='Number of Steps Along Fieldline',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nsteps',ier)
         CALL write_scalar_hdf5(fid,'npoinc',ier,INTVAR=npoinc,ATT='Number of steps per field period',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'npoinc',ier)
         CALL write_var_hdf5(fid,'R_lines',nlines,nsteps+1,ier,DBLVAR=R_lines,ATT='Cylindrical R of Fieldline [m]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'R_lines',ier)
         CALL write_var_hdf5(fid,'Z_lines',nlines,nsteps+1,ier,DBLVAR=Z_lines,ATT='Cylindrical Z of Fieldline [m]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'Z_lines',ier)
         CALL write_var_hdf5(fid,'PHI_lines',nlines,nsteps+1,ier,DBLVAR=PHI_lines,ATT='Cylindrical Phi of Fieldline [rad]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'PHI_lines',ier)
         IF (lmodb) THEN
            CALL write_var_hdf5(fid,'B_lines',nlines,nsteps+1,ier,DBLVAR=B_lines,ATT='|B| Along fieldline [T]',ATT_NAME='description')
            IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'B_lines',ier)
         END IF
      END IF
      ! Wall Data
      IF (ALLOCATED(vertex)) THEN
         CALL write_var_hdf5(fid,'wall_vertex',nvertex,3,ier,DBLVAR=vertex,ATT='Wall Verticies (x,y,z) [m]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'wall_vertex',ier)
         DEALLOCATE(vertex)
      END IF
      IF (ALLOCATED(face)) THEN
         CALL write_var_hdf5(fid,'wall_faces',nface,3,ier,INTVAR=face,ATT='Wall Faces',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'wall_faces',ier)
         DEALLOCATE(face)
      END IF
      IF (ALLOCATED(ihit_array)) THEN
         CALL write_var_hdf5(fid,'wall_strikes',nface,ier,INTVAR=ihit_array,&
                                   ATT='Wall Strikes',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'wall_strikes',ier)
         DEALLOCATE(ihit_array)
      END IF
      ! Homocline
      IF (ALLOCATED(Rhc_lines)) THEN
         nlines = SIZE(Rhc_lines,1)
         nsteps = SIZE(Rhc_lines,2)
         CALL write_var_hdf5(fid,'Rhc_lines',nlines,nsteps,ier,DBLVAR=Rhc_lines,ATT='Cylindrical R of Homocline Fieldline [m]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'R_lines',ier)
         CALL write_var_hdf5(fid,'Zhc_lines',nlines,nsteps,ier,DBLVAR=Zhc_lines,ATT='Cylindrical Z of Homocline Fieldline [m]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'Z_lines',ier)
         
      END IF
      ! Here we output the grid
      IF (ladvanced) THEN
         CALL write_scalar_hdf5(fid,'n_grid',ier,INTVAR=n_qshep,ATT='Number of Gridpoints',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'n_grid',ier)
         CALL write_var_hdf5(fid,'r_grid',n_qshep,ier,DBLVAR=x_qshep,ATT='Radial Grid [m]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'r_grid',ier)
         CALL write_var_hdf5(fid,'phi_grid',n_qshep,ier,DBLVAR=y_qshep,ATT='Toroidal Grid [rad]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'phi_grid',ier)
         CALL write_var_hdf5(fid,'z_grid',n_qshep,ier,DBLVAR=z_qshep,ATT='Vertical Grid [m]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'z_grid',ier)
      ELSEIF (lafield_only) THEN
         CALL write_scalar_hdf5(fid,'nr',ier,INTVAR=nr,ATT='Number of Radial Gridpoints',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nr',ier)
         CALL write_scalar_hdf5(fid,'nphi',ier,INTVAR=nphi,ATT='Number of Toroidal Gridpoints',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nphi',ier)
         CALL write_scalar_hdf5(fid,'nz',ier,INTVAR=nz,ATT='Number of Vertical Gridpoints',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nz',ier)
         CALL write_var_hdf5(fid,'raxis',nr,ier,DBLVAR=raxis,ATT='Radial Axis [m]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'raxis',ier)
         CALL write_var_hdf5(fid,'phiaxis',nphi,ier,DBLVAR=phiaxis,ATT='Toroidal Axis [rad]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'phiaxis',ier)
         CALL write_var_hdf5(fid,'zaxis',nz,ier,DBLVAR=zaxis,ATT='Vertical Axis [rad]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'zaxis',ier)
         CALL write_var_hdf5(fid,'A_R',nr,nphi,nz,ier,DBLVAR=B_R,ATT='Radial Fieldline Eq. (AR)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'A_R',ier)
         CALL write_var_hdf5(fid,'A_Z',nr,nphi,nz,ier,DBLVAR=B_Z,ATT='Vertical Fieldline Eq. (AZ)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'A_Z',ier)
         CALL write_var_hdf5(fid,'A_PHI',nr,nphi,nz,ier,DBLVAR=B_PHI,ATT='Toroidal Vector Potential (APHI)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'A_PHI',ier)
      ELSEIF (lbfield_only) THEN
         CALL write_scalar_hdf5(fid,'nr',ier,INTVAR=nr,ATT='Number of Radial Gridpoints',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nr',ier)
         CALL write_scalar_hdf5(fid,'nphi',ier,INTVAR=nphi,ATT='Number of Toroidal Gridpoints',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nphi',ier)
         CALL write_scalar_hdf5(fid,'nz',ier,INTVAR=nz,ATT='Number of Vertical Gridpoints',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nz',ier)
         CALL write_var_hdf5(fid,'raxis',nr,ier,DBLVAR=raxis,ATT='Radial Axis [m]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'raxis',ier)
         CALL write_var_hdf5(fid,'phiaxis',nphi,ier,DBLVAR=phiaxis,ATT='Toroidal Axis [rad]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'phiaxis',ier)
         CALL write_var_hdf5(fid,'zaxis',nz,ier,DBLVAR=zaxis,ATT='Vertical Axis [rad]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'zaxis',ier)
         CALL write_var_hdf5(fid,'B_R',nr,nphi,nz,ier,DBLVAR=B_R,ATT='Radial Fieldline Eq. (BR)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'B_R',ier)
         CALL write_var_hdf5(fid,'B_Z',nr,nphi,nz,ier,DBLVAR=B_Z,ATT='Vertical Fieldline Eq. (BZ)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'B_Z',ier)
         CALL write_var_hdf5(fid,'B_PHI',nr,nphi,nz,ier,DBLVAR=B_PHI,ATT='Toroidal Field (BPHI)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'B_PHI',ier)
      ELSE
         CALL write_scalar_hdf5(fid,'nr',ier,INTVAR=nr,ATT='Number of Radial Gridpoints',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nr',ier)
         CALL write_scalar_hdf5(fid,'nphi',ier,INTVAR=nphi,ATT='Number of Toroidal Gridpoints',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nphi',ier)
         CALL write_scalar_hdf5(fid,'nz',ier,INTVAR=nz,ATT='Number of Vertical Gridpoints',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nz',ier)
         CALL write_var_hdf5(fid,'raxis',nr,ier,DBLVAR=raxis,ATT='Radial Axis [m]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'raxis',ier)
         CALL write_var_hdf5(fid,'phiaxis',nphi,ier,DBLVAR=phiaxis,ATT='Toroidal Axis [rad]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'phiaxis',ier)
         CALL write_var_hdf5(fid,'zaxis',nz,ier,DBLVAR=zaxis,ATT='Vertical Axis [rad]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'zaxis',ier)
         CALL write_var_hdf5(fid,'B_R',nr,nphi,nz,ier,DBLVAR=B_R,ATT='Radial Fieldline Eq. (BR/BPHI)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'B_R',ier)
         CALL write_var_hdf5(fid,'B_Z',nr,nphi,nz,ier,DBLVAR=B_Z,ATT='Vertical Fieldline Eq. (BZ/BPHI)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'B_Z',ier)
         CALL write_var_hdf5(fid,'B_PHI',nr,nphi,nz,ier,DBLVAR=B_PHI,ATT='Toroidal Field (BPHI)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'B_PHI',ier)
      END IF
      CALL close_hdf5(fid,ier)
      IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,'fieldlines_'//TRIM(id_string)//'.h5',ier)
!DEC$ ELSE
      iunit = 100
      WRITE(6,'(A)')  '   FILE: '//'fieldlines_'//TRIM(id_string)//'.bin'
      CALL safe_open(iunit,ier,'fieldlines_'//TRIM(id_string)//'.bin','replace','unformatted')
      WRITE(iunit) FIELDLINES_VERSION
      WRITE(iunit) lvmec,lpies,lspec,lcoil,lmgrid,lmu,lvessel,lvac,laxis_i,ladvanced,lreverse
      IF (.not. lbfield_only .and. .not. lafield_only .and. .not. lemc3) THEN
         WRITE(iunit) nlines,nsteps,npoinc
         WRITE(iunit) R_lines
         WRITE(iunit) Z_lines
         WRITE(iunit) PHI_lines
      END IF
      IF (ALLOCATED(Rhc_lines)) THEN
         nlines = SIZE(Rhc_lines,1)
         nsteps = SIZE(Rhc_lines,2)
         WRITE(iunit) nlines,nsteps
         WRITE(iunit) Rhc_lines
         WRITE(iunit) Zhc_lines
      END IF
      WRITE(iunit) nr,nphi,nz
      WRITE(iunit) raxis
      WRITE(iunit) phiaxis
      WRITE(iunit) zaxis
      WRITE(iunit) B_R
      WRITE(iunit) B_Z
      WRITE(iunit) B_PHI
      CLOSE(iunit)
!DEC$ ENDIF  

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE fieldlines_write
