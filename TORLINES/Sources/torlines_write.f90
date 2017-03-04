!-----------------------------------------------------------------------
!     Subroutine:    torlines_write
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          06/9/2015
!     Description:   This subroutine writes out the fieldline data
!-----------------------------------------------------------------------
      SUBROUTINE torlines_write(write_type)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE torlines_fieldlines
      USE torlines_background
      USE torlines_realspace
      USE torlines_runtime
      USE ez_hdf5
      USE mpi_params
!-----------------------------------------------------------------------
!     Input Variables
!          write_type  Type of write to preform
!-----------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN):: write_type
!-----------------------------------------------------------------------
!     Local Variables
!          ierr        Error flag
!-----------------------------------------------------------------------
!DEC$ IF DEFINED (MPI_OPT)
      INCLUDE 'mpif.h'
!DEC$ ENDIF  
      INTEGER :: ier
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ier = 0
      IF (myid == master) THEN
         SELECT CASE (TRIM(write_type))
            CASE('GRID_INIT')
               CALL open_hdf5('torlines_'//TRIM(id_string)//'.h5',fid,ier,LCREATE=.true.)
               IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,'torlines_'//TRIM(id_string)//'.h5',ier)
               CALL write_scalar_hdf5(fid,'VERSION',ier,DBLVAR=TORLINES_VERSION,ATT='Version Number',ATT_NAME='description')
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
               CALL write_scalar_hdf5(fid,'lvac',ier,BOOVAR=lvac,ATT='Vacuum calc',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lvac',ier)
               CALL write_scalar_hdf5(fid,'k',ier,INTVAR=k,ATT='Number of Rho Gridpoints',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'k',ier)
               CALL write_scalar_hdf5(fid,'nu',ier,INTVAR=nu,ATT='Number of Poloidal Gridpoints',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nu',ier)
               CALL write_scalar_hdf5(fid,'nv',ier,INTVAR=nv,ATT='Number of Toroidal Gridpoints',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nv',ier)
               CALL write_scalar_hdf5(fid,'nfp',ier,INTVAR=nfp,ATT='Number of Field Periods',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nfp',ier)
               CALL write_var_hdf5(fid,'rreal',k,nu,nv,ier,DBLVAR=rreal,ATT='R Gridpoints',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'rreal',ier)
               CALL write_var_hdf5(fid,'zreal',k,nu,nv,ier,DBLVAR=zreal,ATT='Z Gridpoints',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'zreal',ier)
            CASE('FIELD')
               CALL open_hdf5('torlines_'//TRIM(id_string)//'.h5',fid,ier,LCREATE=.false.)
               IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,'torlines_'//TRIM(id_string)//'.h5',ier)
               CALL write_var_hdf5(fid,'bsreal',k,nu,nv,ier,DBLVAR=bsreal,ATT='B^S Real',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'bsreal',ier)
               CALL write_var_hdf5(fid,'bureal',k,nu,nv,ier,DBLVAR=bureal,ATT='B^U Real',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'bureal',ier)
               CALL write_var_hdf5(fid,'bvreal',k,nu,nv,ier,DBLVAR=bvreal,ATT='B^V Real',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'bvreal',ier)
               CALL write_var_hdf5(fid,'breal',k,nu,nv,ier,DBLVAR=breal,ATT='|B| Real',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'breal',ier)
            CASE('FIELDLINES')
               CALL open_hdf5('torlines_'//TRIM(id_string)//'.h5',fid,ier,LCREATE=.false.)
               IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,'torlines_'//TRIM(id_string)//'.h5',ier)
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
               CALL write_var_hdf5(fid,'B_lines',nlines,nsteps+1,ier,DBLVAR=B_lines,ATT='|B| of Fieldline [T]',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'PHI_lines',ier)
               CALL write_var_hdf5(fid,'U_lines',nlines,nsteps+1,ier,DBLVAR=U_lines,ATT='Poloidal-like angle of Fieldlines [rad]',ATT_NAME='description')
               IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'PHI_lines',ier)
         END SELECT
         CALL close_hdf5(fid,ier)
         IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,'torlines_'//TRIM(id_string)//'.h5',ier)
      END IF
!DEC$ IF DEFINED (MPI_OPT)
      ! Make sure we wait till we're done
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'torlines_write',ierr_mpi)
!DEC$ ENDIF  
   
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE torlines_write
