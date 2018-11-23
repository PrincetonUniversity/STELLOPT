!-----------------------------------------------------------------------
!     Module:        beams3d_write_parhdf5
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/03/2017
!     Description:   Write trajectory info using parallel HDF5. Note,
!                    you must be on a parallel file system for this to
!                    work.
!-----------------------------------------------------------------------
      MODULE beams3d_write_par
      CONTAINS
      SUBROUTINE beams3d_write_parhdf5(n1,n2,m1,m2,mystart,myend,var_name,INTVAR,FLTVAR,DBLVAR)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE mpi_params
      USE hdf5
      USE beams3d_runtime, ONLY: id_string, nprocs_beams, handle_err, &
                                 MPI_BARRIER_ERR
!-----------------------------------------------------------------------
!     Input Variables
!          var_name     Name of variable
!          n1,n2        Extent of first dimesion of global array
!          m1,m2        Extent of second dimension of global array
!          mystart      Starting point of data (2nd dim)
!          myend        Ending point of data (2nd dim)
!          INTVAR       Array (INTEGER)
!          FLTVAR       Array (REAL)
!          DBLVAR       Array (DOUBLE PRCISION)
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER                                :: n1,n2,m1,m2,mystart,myend
      CHARACTER(LEN=*), INTENT(in)           :: var_name
      INTEGER, INTENT(in),  OPTIONAL         :: INTVAR(n1:n2,mystart:myend)
      REAL, INTENT(in),  OPTIONAL            :: FLTVAR(n1:n2,mystart:myend)
      DOUBLE PRECISION, INTENT(in), OPTIONAL :: DBLVAR(n1:n2,mystart:myend)
!-----------------------------------------------------------------------
!     Local Variables
!          ier          Error Flag
!          info         ?
!          plist_id     Property list
!          file_id      HDF5 file ID
!          rank         Data space rank
!          dimsf        Data space dimensions
!          fspace_id    Data space identifier
!          dset_id      Dataset identifier
!          mspace_id    Memory identifier
!          dataspace    Data identifier
!-----------------------------------------------------------------------
      INCLUDE 'mpif.h'
      LOGICAL :: livar, lfvar, ldvar
      INTEGER :: ier, info, rank, i
      INTEGER(HID_T) :: file_id, fspace_id, dset_id, mspace_id, &
                        fapl_id, dcpl_id, dxpl_id, driver_id
      INTEGER(HSIZE_T), ALLOCATABLE :: dimsf(:), counts(:), chunk_dims(:),&
                                       offset(:)

!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Handle types
      livar = .false.; lfvar=.false.; ldvar=.false.
      IF (PRESENT(INTVAR)) livar=.true.
      IF (PRESENT(FLTVAR)) lfvar=.true.
      IF (PRESENT(DBLVAR)) ldvar=.true.

      ! Setup Helper Arrays
      rank =2;
      ALLOCATE(dimsf(rank),chunk_dims(rank),counts(rank),offset(rank))
      dimsf(1) = n2-n1+1
      dimsf(2) = m2-m1+1
      chunk_dims(1) = n2-n1+1
      chunk_dims(2) = myend-mystart+1
      counts(1) = 1
      counts(2) = 1
      offset(1) = 0
      offset(2) = mystart-1

!DEC$ IF DEFINED (HDF5_PAR)
      ier = 0
      info = MPI_INFO_NULL
      ! Initialize
      CALL h5open_f(ier)
      ! Setup File access
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ier)
      CALL h5pset_fapl_mpio_f(fapl_id, MPI_COMM_BEAMS, info, ier)
      CALL h5pget_driver_f(fapl_id, driver_id, ier)
      ! Open file
      CALL h5fopen_f('beams3d_'//TRIM(id_string)//'.h5', H5F_ACC_RDWR_F, file_id, ier, access_prp = fapl_id)
!!!!!!! Begin writing

      WRITE(6,*) myworkid,dimsf,chunk_dims,counts,offset; CALL FLUSH(6)

      ! Create Spaces
      CALL h5screate_simple_f(rank, dimsf, fspace_id, ier)
      CALL h5screate_simple_f(rank, dimsf, mspace_id, ier)
!      CALL h5screate_simple_f(rank, chunk_dims, mspace_id, ier)


      ! Enable Chunking
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ier)
      CALL h5pset_chunk_f(dcpl_id, rank, chunk_dims, ier)

      ! Create Dataset
      IF (livar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_INTEGER, fspace_id, dset_id, ier, dcpl_id)
      IF (lfvar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_DOUBLE, fspace_id, dset_id, ier, dcpl_id)
      IF (ldvar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_DOUBLE, fspace_id, dset_id, ier, dcpl_id)


      ! Select Hyperslab in memory
      CALL h5sselect_hyperslab_f(mspace_id, H5S_SELECT_SET_F, offset, chunk_dims, ier)

      ! Select Hyperslab in File
      CALL h5sselect_hyperslab_f(fspace_id, H5S_SELECT_SET_F, offset, chunk_dims, ier)

      ! Create Properties
      CALL h5pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ier)
      CALL h5pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F, ier)


      ! Write dataset
      IF (livar) CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, INTVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id, xfer_prp = dxpl_id)
      IF (lfvar) CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, FLTVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id, xfer_prp = dxpl_id)
      IF (ldvar) CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, DBLVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id, xfer_prp = dxpl_id)

      ! Close Property list
      CALL h5pclose_f(fapl_id, ier)
      CALL h5pclose_f(dcpl_id, ier)
      CALL h5pclose_f(dxpl_id, ier)
      CALL h5sclose_f(mspace_id, ier)
      CALL h5sclose_f(fspace_id, ier)
      CALL h5dclose_f(dset_id, ier)

      ! Deallocate Helpers
      DEALLOCATE(dimsf,chunk_dims,counts,offset)

!!!!!!!CLOSE FILE
      ! Close the file
      CALL h5fclose_f(file_id, ier)
      ! Close the fortran interface
      CALL h5close_f(ier)

!DEC$ ELSE

      DO i = 0, nprocs_beams-1
         IF (myworkid == i) THEN
            ! Open the fotran interface
            CALL h5open_f(ier)
            !PRINT *,'h5open ',ier

            ! Setup File access
            CALL h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ier)
            !PRINT *,'h5pcreate_f ',ier

            ! Open file
            CALL h5fopen_f('beams3d_'//TRIM(id_string)//'.h5', H5F_ACC_RDWR_F, file_id, ier, access_prp = H5P_DEFAULT_F)
            !PRINT *,'h5fopen_f ',ier

            ! Open or create the dataset and get/create dataspace identifer
            IF  (myworkid == master) THEN
               CALL h5screate_simple_f(rank, dimsf, fspace_id, ier)
            !PRINT *,'h5screate_simple_f ',ier
               IF (livar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_INTEGER, fspace_id, dset_id, ier)
               IF (lfvar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_DOUBLE, fspace_id, dset_id, ier)
               IF (ldvar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_DOUBLE, fspace_id, dset_id, ier)
            !PRINT *,'h5dcreate_f ',ier
            ELSE
               CALL h5dopen_f(file_id, TRIM(var_name), dset_id, ier)
            !PRINT *,'h5dopen_f ',ier
               CALL h5dget_space_f(dset_id, fspace_id, ier)
            !PRINT *,'h5dget_space_f ',ier
            END IF

            ! Select Hyperslab
            CALL h5sselect_hyperslab_f(fspace_id, H5S_SELECT_SET_F,offset, chunk_dims, ier) 
            !PRINT *,'h5sselect_hyperslab_f ',ier
 
            ! Get Memory Space
            CALL h5screate_simple_f(rank, chunk_dims, mspace_id, ier)
            !PRINT *,'h5screate_simple_f ',ier
            
            ! Enable Chunking
            CALL h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ier)
            !PRINT *,'h5pcreate_f ',ier
            CALL h5pset_chunk_f(dcpl_id, rank, chunk_dims, ier)
            !PRINT *,'h5pset_chunk_f ',ier

            ! Write dataset
            IF (livar) CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, INTVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id)
            IF (lfvar) CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, FLTVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id)
            IF (ldvar) CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, DBLVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id)
            !PRINT *,'h5dwrite_f ',ier
            
            ! Close down
            CALL h5sclose_f(fspace_id, ier)
            CALL h5sclose_f(mspace_id, ier)
            CALL h5dclose_f(dset_id, ier)
      
            ! Close the file
            CALL h5fclose_f(file_id, ier)
      
            ! Close the fortran interface
            CALL h5close_f(ier)
         END IF
         CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
         IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'beams3d_init_coil',ierr_mpi)
         !IF (i == 1) STOP
      END DO

!DEC$ ENDIF

      ! Deallocate Helpers
      DEALLOCATE(dimsf,chunk_dims,counts,offset)
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE beams3d_write_parhdf5
      END MODULE beams3d_write_par
