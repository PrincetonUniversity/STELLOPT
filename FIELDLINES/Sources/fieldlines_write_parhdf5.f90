!-----------------------------------------------------------------------
!     Module:        fieldlines_write_par
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/11/2017
!     Description:   Write trajectory info using parallel HDF5
!-----------------------------------------------------------------------
      MODULE fieldlines_write_par
      CONTAINS
      SUBROUTINE fieldlines_write1d_parhdf5(n1,n2,mystart,myend,var_name,INTVAR,FLTVAR,DBLVAR)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE mpi_params
      USE hdf5
      USE fieldlines_runtime, ONLY: id_string, handle_err, &
                                 MPI_BARRIER_ERR, nprocs_fieldlines
      USE mpi_inc
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
      INTEGER                                :: n1,n2,mystart,myend
      CHARACTER(LEN=*), INTENT(in)           :: var_name
      INTEGER, INTENT(in),  OPTIONAL         :: INTVAR(mystart:myend)
      REAL, INTENT(in),  OPTIONAL            :: FLTVAR(mystart:myend)
      DOUBLE PRECISION, INTENT(in), OPTIONAL :: DBLVAR(mystart:myend)
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
      LOGICAL :: livar, lfvar, ldvar
      INTEGER :: ier, info, rank, i
      INTEGER(HID_T) :: file_id, fspace_id, dset_id, mspace_id, &
                        fapl_id, dcpl_id, dxpl_id, driver_id
      INTEGER(HSIZE_T), ALLOCATABLE :: dimsf(:), counts(:), chunk_dims(:),&
                                       offset(:), block(:), stride(:)

!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Handle types
      livar = .false.; lfvar=.false.; ldvar=.false.
      IF (PRESENT(INTVAR)) livar=.true.
      IF (PRESENT(FLTVAR)) lfvar=.true.
      IF (PRESENT(DBLVAR)) ldvar=.true.

      ! Setup Helper Arrays
      rank =1;
      ALLOCATE(dimsf(rank),chunk_dims(rank),counts(rank),offset(rank))
      dimsf(1) = n2-n1+1
      chunk_dims(1) = myend-mystart+1
      counts(1) = 1
      offset(1) = mystart-1

!DEC$ IF DEFINED (HDF5_PAR)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,chunk_dims(1),1,MPI_INTEGER,MPI_MAX,MPI_COMM_FIELDLINES,ier)
      ! Setup Helper Arrays
      ALLOCATE(stride(rank),block(rank))
      stride(1) = 1
      block(1)  = chunk_dims(1)

      ier = 0
      info = MPI_INFO_NULL

      ! Opent the fortran interface
      CALL h5open_f(ier)

      ! Setup the File access
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ier)

      ! Setup the parallel communicator
      CALL h5pset_fapl_mpio_f(fapl_id, MPI_COMM_FIELDLINES, info, ier)

      ! Open file
      CALL h5fopen_f('fieldlines_'//TRIM(id_string)//'.h5', H5F_ACC_RDWR_F, file_id, ier, access_prp = fapl_id)
      CALL h5pclose_f(fapl_id,ier)

      ! Create File Space
      CALL h5screate_simple_f(rank, dimsf, fspace_id, ier)

      ! Enable Chunking
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ier)
      CALL h5pset_chunk_f(dcpl_id, rank, chunk_dims, ier)

      ! Create Dataset
      IF (livar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_INTEGER, fspace_id, dset_id, ier, dcpl_id)
      IF (lfvar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_DOUBLE, fspace_id, dset_id, ier, dcpl_id)
      IF (ldvar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_DOUBLE, fspace_id, dset_id, ier, dcpl_id)
      CALL h5pclose_f(dcpl_id, ier)


      ! Create Memory Space
      chunk_dims(1) = myend-mystart+1
      CALL h5screate_simple_f(rank, chunk_dims, mspace_id, ier)

      ! Select Hyperslab in File
      !CALL h5sselect_hyperslab_f(fspace_id, H5S_SELECT_SET_F, offset, chunk_dims, ier)
      CALL h5sselect_hyperslab_f(fspace_id, H5S_SELECT_SET_F, offset, counts, ier,stride,block)

      ! Create Properties
      CALL h5pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ier)
      CALL h5pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F, ier)

      ! Write dataset
      IF (livar) CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, INTVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id, xfer_prp = dxpl_id)
      IF (lfvar) CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, FLTVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id, xfer_prp = dxpl_id)
      IF (ldvar) CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, DBLVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id, xfer_prp = dxpl_id)
      CALL h5pclose_f(dxpl_id, ier)


      ! Close Property list
      CALL h5sclose_f(fspace_id, ier)
      CALL h5sclose_f(mspace_id, ier)
      CALL h5dclose_f(dset_id, ier)

      ! Close the file
      CALL h5fclose_f(file_id, ier)

      ! Close the fortran interface
      CALL h5close_f(ier)

      DEALLOCATE(stride,block)

!DEC$ ELSE

      DO i = 0, nprocs_fieldlines-1
         IF (myworkid == i) THEN
            ! Open the fotran interface
            CALL h5open_f(ier)

            ! Setup File access
            CALL h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ier)

            ! Open file
            CALL h5fopen_f('fieldlines_'//TRIM(id_string)//'.h5', H5F_ACC_RDWR_F, file_id, ier, access_prp = H5P_DEFAULT_F)

            ! Open or create the dataset and get/create dataspace identifer
            IF  (myworkid == master) THEN
               CALL h5screate_simple_f(rank, dimsf, fspace_id, ier)
               IF (livar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_INTEGER, fspace_id, dset_id, ier)
               IF (lfvar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_DOUBLE, fspace_id, dset_id, ier)
               IF (ldvar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_DOUBLE, fspace_id, dset_id, ier)
            ELSE
               CALL h5dopen_f(file_id, TRIM(var_name), dset_id, ier)
               CALL h5dget_space_f(dset_id, fspace_id, ier)
            END IF

            ! Select Hyperslab
            CALL h5sselect_hyperslab_f(fspace_id, H5S_SELECT_SET_F,offset, chunk_dims, ier) 
 
            ! Get Memory Space
            CALL h5screate_simple_f(rank, chunk_dims, mspace_id, ier)
            
            ! Enable Chunking
            CALL h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ier)
            CALL h5pset_chunk_f(dcpl_id, rank, chunk_dims, ier)

            ! Write dataset
            IF (livar) CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, INTVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id)
            IF (lfvar) CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, FLTVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id)
            IF (ldvar) CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, DBLVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id)
            
            ! Close down
            CALL h5sclose_f(fspace_id, ier)
            CALL h5sclose_f(mspace_id, ier)
            CALL h5dclose_f(dset_id, ier)
      
            ! Close the file
            CALL h5fclose_f(file_id, ier)
      
            ! Close the fortran interface
            CALL h5close_f(ier)
         END IF
         CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
         IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_write_parhdf5',ierr_mpi)
      END DO

!DEC$ ENDIF

      ! Deallocate Helpers
      DEALLOCATE(dimsf,chunk_dims,counts,offset)
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE fieldlines_write1d_parhdf5

      SUBROUTINE fieldlines_write2d_parhdf5(n1,n2,m1,m2,mystart,myend,var_name,INTVAR,FLTVAR,DBLVAR)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE mpi_params
      USE hdf5
      USE fieldlines_runtime, ONLY: id_string, handle_err, &
                                 MPI_BARRIER_ERR, nprocs_fieldlines
      USE mpi_inc
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
      INTEGER, INTENT(in),  OPTIONAL         :: INTVAR(mystart:myend,m1:m2)
      REAL, INTENT(in),  OPTIONAL            :: FLTVAR(mystart:myend,m1:m2)
      DOUBLE PRECISION, INTENT(in), OPTIONAL :: DBLVAR(mystart:myend,m1:m2)
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
      LOGICAL :: livar, lfvar, ldvar
      INTEGER :: ier, info, rank, i
      INTEGER(HID_T) :: file_id, fspace_id, dset_id, mspace_id, &
                        fapl_id, dcpl_id, dxpl_id, driver_id
      INTEGER(HSIZE_T), ALLOCATABLE :: dimsf(:), counts(:), chunk_dims(:),&
                                       offset(:), block(:), stride(:)

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
      chunk_dims(1) = myend-mystart+1
      chunk_dims(2) = m2-m1+1
      counts(1) = 1
      counts(2) = 1
      offset(1) = mystart-1
      offset(2) = 0

!DEC$ IF DEFINED (HDF5_PAR)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,chunk_dims(1),1,MPI_INTEGER,MPI_MAX,MPI_COMM_FIELDLINES,ier)
      ! Setup Helper Arrays
      ALLOCATE(stride(rank),block(rank))
      stride(1) = 1
      stride(2) = 1
      block(1)  = chunk_dims(1)
      block(2)  = chunk_dims(2)

      ier = 0
      info = MPI_INFO_NULL

      ! Opent the fortran interface
      CALL h5open_f(ier)

      ! Setup the File access
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ier)

      ! Setup the parallel communicator
      CALL h5pset_fapl_mpio_f(fapl_id, MPI_COMM_FIELDLINES, info, ier)

      ! Open file
      CALL h5fopen_f('fieldlines_'//TRIM(id_string)//'.h5', H5F_ACC_RDWR_F, file_id, ier, access_prp = fapl_id)
      CALL h5pclose_f(fapl_id,ier)

      ! Create File Space
      CALL h5screate_simple_f(rank, dimsf, fspace_id, ier)

      ! Enable Chunking
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ier)
      CALL h5pset_chunk_f(dcpl_id, rank, chunk_dims, ier)

      ! Create Dataset
      IF (livar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_INTEGER, fspace_id, dset_id, ier, dcpl_id)
      IF (lfvar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_DOUBLE, fspace_id, dset_id, ier, dcpl_id)
      IF (ldvar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_DOUBLE, fspace_id, dset_id, ier, dcpl_id)

      ! Close the file space
      CALL h5sclose_f(fspace_id, ier)

      ! Create Memory Space
      chunk_dims(1) = myend-mystart+1
      block(1) = chunk_dims(1)
      CALL h5screate_simple_f(rank, chunk_dims, mspace_id, ier)

      ! Select Hyperslab in File
      CALL h5dget_space_f(dset_id, fspace_id, ier)
      CALL h5sselect_hyperslab_f(fspace_id, H5S_SELECT_SET_F, offset, counts, ier,stride,block)

      ! Create Properties
      CALL h5pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ier)
      CALL h5pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F, ier)

      ! Write dataset
      IF (livar) CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, INTVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id, xfer_prp = dxpl_id)
      IF (lfvar) CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, FLTVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id, xfer_prp = dxpl_id)
      IF (ldvar) CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, DBLVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id, xfer_prp = dxpl_id)
      CALL h5pclose_f(dxpl_id, ier)


      ! Close Property list
      CALL h5pclose_f(fapl_id, ier)
      CALL h5pclose_f(dcpl_id, ier)
      CALL h5pclose_f(dxpl_id, ier)
      CALL h5sclose_f(mspace_id, ier)
      CALL h5sclose_f(fspace_id, ier)
      CALL h5dclose_f(dset_id, ier)

      ! Close the file
      CALL h5fclose_f(file_id, ier)

      ! Close the fortran interface
      CALL h5close_f(ier)

      DEALLOCATE(stride,block)

!DEC$ ELSE

      DO i = 0, nprocs_fieldlines-1
         IF (myworkid == i) THEN
            ! Open the fotran interface
            CALL h5open_f(ier)

            ! Setup File access
            CALL h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ier)

            ! Open file
            CALL h5fopen_f('fieldlines_'//TRIM(id_string)//'.h5', H5F_ACC_RDWR_F, file_id, ier, access_prp = H5P_DEFAULT_F)

            ! Open or create the dataset and get/create dataspace identifer
            IF  (myworkid == master) THEN
               CALL h5screate_simple_f(rank, dimsf, fspace_id, ier)
               IF (livar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_INTEGER, fspace_id, dset_id, ier)
               IF (lfvar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_DOUBLE, fspace_id, dset_id, ier)
               IF (ldvar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_DOUBLE, fspace_id, dset_id, ier)
            ELSE
               CALL h5dopen_f(file_id, TRIM(var_name), dset_id, ier)
               CALL h5dget_space_f(dset_id, fspace_id, ier)
            END IF

            ! Select Hyperslab
            CALL h5sselect_hyperslab_f(fspace_id, H5S_SELECT_SET_F,offset, chunk_dims, ier) 
 
            ! Get Memory Space
            CALL h5screate_simple_f(rank, chunk_dims, mspace_id, ier)
            
            ! Enable Chunking
            CALL h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ier)
            CALL h5pset_chunk_f(dcpl_id, rank, chunk_dims, ier)

            ! Write dataset
            IF (livar) CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, INTVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id)
            IF (lfvar) CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, FLTVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id)
            IF (ldvar) CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, DBLVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id)
            
            ! Close down
            CALL h5sclose_f(fspace_id, ier)
            CALL h5sclose_f(mspace_id, ier)
            CALL h5dclose_f(dset_id, ier)
      
            ! Close the file
            CALL h5fclose_f(file_id, ier)
      
            ! Close the fortran interface
            CALL h5close_f(ier)
         END IF
         CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
         IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_write_parhdf5',ierr_mpi)
      END DO

!DEC$ ENDIF

      ! Deallocate Helpers
      DEALLOCATE(dimsf,chunk_dims,counts,offset)
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE fieldlines_write2d_parhdf5
      END MODULE fieldlines_write_par
