!-----------------------------------------------------------------------
!     Module:        beams3d_write_parhdf5
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/03/2017
!     Description:   Write trajectory info using parallel HDF5. Note,
!                    you must be on a parallel file system for this to
!                    work.
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_write_parhdf5(n1,n2,m1,m2,mystart,myend,Farr,var_name)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE mpi_params
      USE hdf5
      USE beams3d_runtime, ONLY: id_string
!-----------------------------------------------------------------------
!     Input Variables
!          mystart      Starting point of data
!          myend        Ending point of data
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: n1,n2,m1,m2,mystart,myend
      DOUBLE PRECISION :: Farr(n1:n2,mystart:myend)
      CHARACTER(LEN=*), INTENT(in)  :: var_name
!-----------------------------------------------------------------------
!     Local Variables
!          ier          Error Flag
!          info         ?
!          plist_id     Property list
!          file_id      HDF5 file ID
!          rank         Data space rank
!          dimsf        Data space dimensions
!          filespace    Data space identifier
!          dset_id      Dataset identifier
!          memspace     Memory identifier
!          dataspace    Data identifier
!-----------------------------------------------------------------------
      INCLUDE 'mpif.h'
      INTEGER :: ier, info, rank
      INTEGER(HID_T) :: plist_id, file_id, filespace, dset_id, memspace, &
                        dataspace
      INTEGER(HSIZE_T), ALLOCATABLE :: dimsf(:), counts(:)
      INTEGER(HSSIZE_T), ALLOCATABLE :: offset(:)

!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ier = 0
      info = MPI_INFO_NULL
      ! Initialize
      CALL h5open_f(ier)
      ! Setup File access
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ier)
      CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_BEAMS, info, ier)
      ! Open file
      CALL h5fopen_f('beams3d_'//TRIM(id_string)//'.h5', H5F_ACC_RDWR_F, file_id, ier, access_prp = plist_id)
      ! Close property list      CALL h5pclose_f(plist_id, ier)
!!!!!!! Begin writing

      ! Create Data Space
      rank =2;
      ALLOCATE(dimsf(rank))
      dimsf(1) = n2-n1+1
      dimsf(2) = m2-m1+1
      CALL h5screate_simple_f(rank, dimsf, filespace, ier)
      ! Create the Dataset      CALL h5dcreate_f(file_id, var_name, H5T_NATIVE_DOUBLE, filespace, dset_id, ier)
      ! Close the filespace
      CALL h5sclose_f(filespace, ier)

      ! Create the hyperslab in memory
      ALLOCATE(counts(rank))
      ALLOCATE(offset(rank))
      counts(1) = SIZE(Farr,1)
      counts(2) = SIZE(Farr,2)
      offset(1) = 0
      offset(2) = mystart-1
      CALL h5screate_simple_f(rank, counts, memspace, ier)
      ! Select the hyperslab      CALL h5dget_space_f(dset_id, dataspace, ier)      CALL h5sselect_hyperslab_f (dataspace, H5S_SELECT_SET_F, offset, counts, ier)

      ! Create a property list for the dataset      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ier)       CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, ier)

      ! Write dataset
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, Farr, dimsf, ier, file_space_id = dataspace, mem_space_id = memspace, xfer_prp = plist_id)

      ! Close Property list      CALL h5pclose_f(plist_id, error)
      ! Close memory space      CALL h5sclose_f(memspace, ier)
      ! Close dataspace      CALL h5sclose_f(dataspace, ier)
      ! Close Dataset      CALL h5dclose_f(dset_id, error)
      ! DEALLOCATE the helpers
      DEALLOCATE(dimsf,counts,offset)

!!!!!!!CLOSE FILE
      ! Close the file      CALL h5fclose_f(file_id, ier)
      ! Close the fortran interface
      CALL h5close_f(ier)
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE beams3d_write_parhdf5
