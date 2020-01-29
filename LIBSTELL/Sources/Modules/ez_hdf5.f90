!-----------------------------------------------------------------------
!     Module:        ez_hdf5
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          12/5/2012
!     Description:   This module helps simplfy the interface to hdf5
!                    files.  The HDF5 module must be loaded for your
!                    compiler. The following compiler commands should be
!                    added to your makefiles to properly compile this
!                    file:
!
!                    -L$(HDF5_HOME)/lib -I$(HDF5_HOME)/include 
!                    -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 
!                    -lpthread -lz -lm
!
!                    At the end of this file there is a commented
!                    section which contains an example showing how to
!                    call these routines.  
!-----------------------------------------------------------------------
      MODULE ez_hdf5
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
#ifdef LHDF5
   
      USE hdf5
!-----------------------------------------------------------------------
!     Module Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER(HID_T) :: fid
      INTEGER(HID_T) :: plistid
!-----------------------------------------------------------------------
!     Module Subroutines
!          open_hdf5         Opens and hdf5 file for input/output
!          write_scalr_hdf5  Outputs a scalar to the HDF5 file.
!          write_var_hdf5    Outputs a vector/array to the HDF5 file.
!          close_hdf5        Closes the HDF5 file
!-----------------------------------------------------------------------
      INTERFACE write_var_hdf5
         MODULE PROCEDURE write_scalar_hdf5, write_vector_hdf5,&
                          write_arr2d_hdf5, write_arr3d_hdf5
!         MODULE PROCEDURE write_vector_hdf5,&
!                          write_arr2d_hdf5, write_arr3d_hdf5
      END INTERFACE
      
      INTERFACE read_var_hdf5
         ! GCC complains
         MODULE PROCEDURE read_scalar_hdf5, read_vector_hdf5,&
                          read_arr2d_hdf5, read_arr3d_hdf5
!         MODULE PROCEDURE read_vector_hdf5,&
!                          read_arr2d_hdf5, read_arr3d_hdf5
      END INTERFACE
      
      
      CONTAINS
      
      !-----------------------------------------------------------------
      SUBROUTINE open_hdf5(filename,file_id,ier, lcreate,comm,info)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER(HID_T), INTENT(out)  :: file_id
      INTEGER, INTENT(out)         :: ier
      LOGICAL, OPTIONAL, INTENT(in) :: lcreate
      INTEGER, OPTIONAL, INTENT(in) :: comm
      INTEGER, OPTIONAL, INTENT(in) :: info
      INTEGER             ::  strlen, dex
      CHARACTER(LEN=1024) ::  file_temp
      file_id   = -1
      file_temp = TRIM(filename)
      strlen    = LEN(TRIM(file_temp))
      dex=INDEX(file_temp,'.h5',BACK=.TRUE.)
      IF (dex /= strlen-2) file_temp = TRIM(file_temp) // '.h5'
      CALL h5open_f(ier)
      IF (ier /=0) RETURN
#ifdef HDF5_PAR
      IF (PRESENT(comm)) THEN
         CALL h5pcreate_f(H5P_FILE_ACCESS_F, plistid, ier)
         CALL h5pset_fapl_mpio_f(plistid, comm, info, ier)
         IF (ier /=0) RETURN
      ELSE
         plistid = H5P_DEFAULT_F
      END IF
#else
      plistid = H5P_DEFAULT_F
#endif
      IF (PRESENT(lcreate)) THEN
         IF (lcreate) THEN
            CALL h5fcreate_f(TRIM(file_temp),H5F_ACC_TRUNC_F, file_id,ier, access_prp = plistid)
         ELSE
            CALL h5fopen_f(file_temp, H5F_ACC_RDWR_F, file_id, ier, access_prp = plistid)
         END IF
      ELSE
         CALL h5fopen_f(file_temp, H5F_ACC_RDWR_F, file_id, ier, access_prp = plistid)
      END IF
#ifdef HDF5_PAR
      IF (PRESENT(comm)) CALL h5pclose_f(plistid,ier)
#endif
      END SUBROUTINE open_hdf5
      !-----------------------------------------------------------------
         
      !-----------------------------------------------------------------
      SUBROUTINE close_hdf5(file_id,ier)
      IMPLICIT NONE
      INTEGER(HID_T), INTENT(inout):: file_id
      INTEGER, INTENT(out)         :: ier
      ier = 0
#ifdef HDF5_PAR
      IF (plistid /= H5P_DEFAULT_F) CALL h5pclose_f(plistid, ier)
#endif
      IF (ier /=0) RETURN
      CALL h5fclose_f(fid,ier)
      IF (ier /=0) RETURN
      CALL h5close_f(ier)
      IF (ier /=0) RETURN
      END SUBROUTINE close_hdf5
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      SUBROUTINE write_att_hdf5(loc_id,ATT_NAME,ATT,ierr)
      IMPLICIT NONE
      INTEGER(HID_T), INTENT(in)    :: loc_id
      CHARACTER(LEN=*), INTENT(in)  :: ATT
      CHARACTER(LEN=*), INTENT(in)  :: ATT_NAME
      INTEGER, INTENT(out)          :: ierr
      INTEGER        :: drank = 1
      INTEGER        :: arank = 1
      INTEGER(HID_T) :: dset_id
      INTEGER(HID_T) :: attr_id  
      INTEGER(HID_T) :: dspace_id
      INTEGER(HID_T) :: aspace_id
      INTEGER(HID_T) :: atype_id
      INTEGER(SIZE_T) :: attrlen
      INTEGER(HSIZE_T), DIMENSION(1) :: ddims = (/1/)
      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/)
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims


      ierr = 0
      attrlen=LEN(TRIM(ATT))
      data_dims(1)=1
      !CALL h5screate_simple_f(drank, ddims, aspace_id, ierr)
      CALL h5screate_f(H5S_SCALAR_F,aspace_id,ierr)
      IF (ierr /=0) RETURN
      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, ierr)
      IF (ierr /=0) RETURN
      CALL h5tset_size_f(atype_id,attrlen,ierr)
      IF (ierr /=0) RETURN
      CALL h5acreate_f(loc_id,TRIM(ATT_NAME), atype_id, aspace_id, attr_id,ierr)
      IF (ierr /=0) RETURN
      CALL h5awrite_f(attr_id,atype_id,TRIM(ATT),data_dims,ierr)
      IF (ierr /=0) RETURN
      CALL h5aclose_f(attr_id,ierr)
      IF (ierr /=0) RETURN
      CALL h5tclose_f(atype_id, ierr)
      IF (ierr /=0) RETURN
      CALL h5sclose_f(aspace_id,ierr)
      IF (ierr /=0) RETURN

      RETURN
      END SUBROUTINE write_att_hdf5   
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      SUBROUTINE write_att_int_hdf5(loc_id,ATT_NAME,ATT,ierr)
      IMPLICIT NONE
      INTEGER(HID_T), INTENT(in)    :: loc_id
      INTEGER, INTENT(in)    :: ATT
      CHARACTER(LEN=*), INTENT(in)  :: ATT_NAME
      INTEGER, INTENT(out)          :: ierr
      INTEGER        :: drank = 1
      INTEGER        :: arank = 1
      INTEGER(HID_T) :: dset_id
      INTEGER(HID_T) :: attr_id  
      INTEGER(HID_T) :: dspace_id
      INTEGER(HID_T) :: aspace_id
      INTEGER(HID_T) :: atype_id
      INTEGER(SIZE_T) :: attrlen
      INTEGER(HSIZE_T), DIMENSION(1) :: ddims = (/1/)
      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/)
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims


      ierr = 0
      data_dims(1)=1
      !CALL h5screate_simple_f(drank, ddims, aspace_id, ierr)
      CALL h5screate_f(H5S_SCALAR_F,aspace_id,ierr)
      IF (ierr /=0) RETURN
      CALL h5acreate_f(loc_id,TRIM(ATT_NAME), H5T_NATIVE_INTEGER, aspace_id, attr_id,ierr)
      IF (ierr /=0) RETURN
      CALL h5awrite_f(attr_id,H5T_NATIVE_INTEGER,ATT,data_dims,ierr)
      IF (ierr /=0) RETURN
      CALL h5aclose_f(attr_id,ierr)
      IF (ierr /=0) RETURN
      CALL h5sclose_f(aspace_id,ierr)
      IF (ierr /=0) RETURN

      RETURN
      END SUBROUTINE write_att_int_hdf5   
      !-----------------------------------------------------------------
      
      !-----------------------------------------------------------------
      SUBROUTINE write_scalar_hdf5(file_id,var,ierr,BOOVAR,INTVAR,FLTVAR,DBLVAR,ATT,ATT_NAME)
      IMPLICIT NONE
      INTEGER(HID_T), INTENT(in)    :: file_id
      CHARACTER(LEN=*), INTENT(in)  :: var
      INTEGER, INTENT(out)          :: ierr
      LOGICAL, INTENT(in), OPTIONAL :: BOOVAR
      INTEGER, INTENT(in), OPTIONAL :: INTVAR
      REAL, INTENT(in), OPTIONAL    :: FLTVAR   
      DOUBLE PRECISION, INTENT(in), OPTIONAL    :: DBLVAR  
      CHARACTER(LEN=*), INTENT(in), OPTIONAL :: ATT
      CHARACTER(LEN=*), INTENT(in), OPTIONAL :: ATT_NAME
!      INTEGER, INTENT(in), OPTIONAL    :: C1
!      INTEGER, INTENT(in), OPTIONAL    :: C2
      INTEGER        :: boo_temp
      INTEGER        :: drank = 1
      INTEGER        :: arank = 1
      INTEGER(HID_T) :: dset_id
      INTEGER(HID_T) :: attr_id  
      INTEGER(HID_T) :: dspace_id
      INTEGER(HID_T) :: aspace_id
      INTEGER(HID_T) :: atype_id
      INTEGER(SIZE_T) :: attrlen
      INTEGER(HSIZE_T), DIMENSION(1) :: ddims = (/1/)
      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/)
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
      CHARACTER(LEN=1024) ::  att_temp
      
      ierr = 0
      CALL h5screate_simple_f(drank,ddims,dspace_id,ierr)
      IF (ierr /=0) RETURN
      IF (PRESENT(INTVAR)) THEN
         CALL h5dcreate_f(file_id,TRIM(var),H5T_NATIVE_INTEGER,dspace_id,dset_id,ierr)
!         IF (plistid /= H5P_DEFAULT_F) CALL h5pset_dxpl_mpio_f(plistid, H5FD_MPIO_COLLECTIVE_F, ierr)
         CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, INTVAR, ddims, ierr, xfer_prp = plistid)
      ELSE IF (PRESENT(FLTVAR)) THEN
         CALL h5dcreate_f(file_id,TRIM(var),H5T_NATIVE_DOUBLE,dspace_id,dset_id,ierr)
!         IF (plistid /= H5P_DEFAULT_F) CALL h5pset_dxpl_mpio_f(plistid, H5FD_MPIO_COLLECTIVE_F, ierr)
         CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, FLTVAR, ddims, ierr, xfer_prp = plistid)
      ELSE IF (PRESENT(DBLVAR)) THEN
         CALL h5dcreate_f(file_id,TRIM(var),H5T_NATIVE_DOUBLE,dspace_id,dset_id,ierr)
!         IF (plistid /= H5P_DEFAULT_F) CALL h5pset_dxpl_mpio_f(plistid, H5FD_MPIO_COLLECTIVE_F, ierr)
         CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, DBLVAR, ddims, ierr, xfer_prp = plistid)
      ELSE IF (PRESENT(BOOVAR)) THEN
         boo_temp = 0
         IF (BOOVAR) boo_temp = 1
         CALL h5dcreate_f(file_id,TRIM(var),H5T_NATIVE_INTEGER,dspace_id,dset_id,ierr)
!         IF (plistid /= H5P_DEFAULT_F) CALL h5pset_dxpl_mpio_f(plistid, H5FD_MPIO_COLLECTIVE_F, ierr)
         CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, boo_temp, ddims, ierr, xfer_prp = plistid)
      ELSE
         ierr=-2
      END IF
      IF (PRESENT(ATT)) THEN
         IF (PRESENT(ATT_NAME)) THEN
            att_temp=TRIM(ATT_NAME)
         ELSE
            att_temp=TRIM('_att')
         END IF
         data_dims(1) = 1
         CALL h5screate_simple_f(arank,adims,aspace_id,ierr)
         IF (ierr /=0) RETURN
         CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id,ierr)
         IF (ierr /=0) RETURN
         attrlen=LEN(TRIM(ATT))
         CALL h5tset_size_f(atype_id,attrlen,ierr)
         IF (ierr /=0) RETURN
         CALL h5acreate_f(dset_id,TRIM(var) // '_' // TRIM(att_temp), atype_id, aspace_id, attr_id,ierr)
         IF (ierr /=0) RETURN
         data_dims(1)=1
         CALL h5awrite_f(attr_id,atype_id,TRIM(ATT),data_dims,ierr)
         IF (ierr /=0) RETURN
         CALL h5aclose_f(attr_id,ierr)
         IF (ierr /=0) RETURN
         CALL h5sclose_f(aspace_id,ierr)
         IF (ierr /=0) RETURN
      END IF
      IF (ierr /=0) RETURN
      CALL h5dclose_f(dset_id,ierr)
      IF (ierr /=0) RETURN
      CALL h5sclose_f(dspace_id,ierr)
      IF (ierr /=0) RETURN
      END SUBROUTINE write_scalar_hdf5   
      !-----------------------------------------------------------------
      
      !-----------------------------------------------------------------
      SUBROUTINE write_vector_hdf5(file_id,var,n,ier,BOOVAR,INTVAR,FLTVAR,DBLVAR,ATT,ATT_NAME)
      IMPLICIT NONE
      INTEGER(HID_T), INTENT(in)    :: file_id
      CHARACTER(LEN=*), INTENT(in)  :: var
      INTEGER, INTENT(in)           :: n
      INTEGER, INTENT(out)          :: ier
      INTEGER, INTENT(in), OPTIONAL :: INTVAR(n)
      LOGICAL, INTENT(in), OPTIONAL :: BOOVAR(n)
      REAL, INTENT(in), OPTIONAL    :: FLTVAR(n)
      DOUBLE PRECISION, INTENT(in), OPTIONAL    :: DBLVAR(n)
      CHARACTER(LEN=*), INTENT(in), OPTIONAL :: ATT
      CHARACTER(LEN=*), INTENT(in), OPTIONAL :: ATT_NAME
!      INTEGER, INTENT(in), OPTIONAL    :: C1
!      INTEGER, INTENT(in), OPTIONAL    :: C2
      INTEGER        :: drank = 1
      INTEGER        :: arank = 1
      INTEGER        :: boo_temp(n)
      INTEGER(HID_T) :: dset_id
      INTEGER(HID_T) :: attr_id  
      INTEGER(HID_T) :: dspace_id
      INTEGER(HID_T) :: mspace_id
      INTEGER(HID_T) :: aspace_id
      INTEGER(HID_T) :: atype_id
      INTEGER(SIZE_T) :: attrlen
      INTEGER(HSIZE_T), DIMENSION(1) :: ddims
      INTEGER(HSIZE_T), DIMENSION(1) :: cdims
      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/)
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
      INTEGER(HSIZE_T), DIMENSION(1) :: counter
      INTEGER(HSIZE_T), DIMENSION(1) :: offset
      INTEGER(HSIZE_T), DIMENSION(1) :: stride
      INTEGER(HSIZE_T), DIMENSION(1) :: block
      CHARACTER(LEN=1024) ::  att_temp
      
      ier = 0
      ddims(1)=n
      !IF (PRESENT(C1)) THEN
      !   IF (PRESENT(INTVAR)) CALL h5dcreate_f(file_id, TRIM(var), H5T_NATIVE_INTEGER, dspace_id, dset_id, ier, plistid)
      !   IF (PRESENT(FLTVAR)) CALL h5dcreate_f(file_id, TRIM(var), H5T_NATIVE_DOUBLE, dspace_id, dset_id, ier, plistid)
      !   IF (PRESENT(DBLVAR)) CALL h5dcreate_f(file_id, TRIM(var), H5T_NATIVE_DOUBLE, dspace_id, dset_id, ier, plistid)
      !   IF (PRESENT(BOOVAR)) CALL h5dcreate_f(file_id, TRIM(var), H5T_NATIVE_INTEGER, dspace_id, dset_id, ier, plistid)
      !   CALL h5sclose_f(dspace_id,ier)
      !   IF (ier /=0) RETURN
      !   cdims(1) = C2-C1+1 ! Length of my segment
      !   offset(1) = C1-1   ! Where I stick it
      !   CALL h5screate_simple_f(drank,cdims,mspace_id,ier)
      !   CALL h5dget_space_f(dset_id, dspace_id, ier)
      !   CALL h5sselect_hyperslab_f (dspace_id, H5S_SELECT_SET_F, offset, cdims, ier)
      !   CALL h5pcreate_f(H5P_DATASET_XFER_F, plistid, ier)
      !   CALL h5pset_dxpl_mpio_f(plistid, H5FD_MPIO_COLLECTIVE_F, ierr)
      !   CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, INTVAR, ddims, ier,&
      !                   file_space_id = dspace_id, mem_space_id = mspace_id, xfer_prp = plist_id)
      !   CALL h5sclose_f(mspace_id, error)
      !END IF
      CALL h5screate_simple_f(drank,ddims,dspace_id,ier)
      IF (ier /=0) RETURN
      IF (PRESENT(INTVAR)) THEN
         CALL h5dcreate_f(file_id,TRIM(var),H5T_NATIVE_INTEGER,dspace_id,dset_id,ier, plistid)
         CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, INTVAR, ddims, ier)
      ELSE IF (PRESENT(FLTVAR)) THEN
         CALL h5dcreate_f(file_id,TRIM(var),H5T_NATIVE_DOUBLE,dspace_id,dset_id,ier, plistid)
         CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, FLTVAR, ddims, ier)
      ELSE IF (PRESENT(DBLVAR)) THEN
         CALL h5dcreate_f(file_id,TRIM(var),H5T_NATIVE_DOUBLE,dspace_id,dset_id,ier, plistid)
         CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, DBLVAR, ddims, ier)
      ELSE IF (PRESENT(BOOVAR)) THEN
         boo_temp = 0
         WHERE(BOOVAR) boo_temp = 1
         CALL h5dcreate_f(file_id,TRIM(var),H5T_NATIVE_INTEGER,dspace_id,dset_id,ier)
         CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, boo_temp, ddims, ier)
      ELSE
         ier=-2
      END IF
      IF (PRESENT(ATT)) THEN
         IF (PRESENT(ATT_NAME)) THEN
            att_temp=TRIM(ATT_NAME)
         ELSE
            att_temp=TRIM('_att')
         END IF
         data_dims(1) = 1
         CALL h5screate_simple_f(arank,adims,aspace_id,ier)
         IF (ier /=0) RETURN
         CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id,ier)
         IF (ier /=0) RETURN
         attrlen=LEN(TRIM(ATT))
         CALL h5tset_size_f(atype_id,attrlen,ier)
         IF (ier /=0) RETURN
         CALL h5acreate_f(dset_id,TRIM(var) // '_' // TRIM(att_temp), atype_id, aspace_id, attr_id,ier)
         IF (ier /=0) RETURN
         data_dims(1)=1
         CALL h5awrite_f(attr_id,atype_id,TRIM(ATT),data_dims,ier)
         IF (ier /=0) RETURN
         CALL h5aclose_f(attr_id,ier)
         IF (ier /=0) RETURN
         CALL h5sclose_f(aspace_id,ier)
         IF (ier /=0) RETURN
      END IF
      IF (ier /=0) RETURN
      CALL h5dclose_f(dset_id,ier)
      IF (ier /=0) RETURN
      CALL h5sclose_f(dspace_id,ier)
      IF (ier /=0) RETURN
      
      END SUBROUTINE write_vector_hdf5   
      !-----------------------------------------------------------------
      
      !-----------------------------------------------------------------
      SUBROUTINE write_arr2d_hdf5(file_id,var,n1,n2,ier,BOOVAR,INTVAR,FLTVAR,DBLVAR,ATT,ATT_NAME)
      IMPLICIT NONE
      INTEGER(HID_T), INTENT(in)    :: file_id
      CHARACTER(LEN=*), INTENT(in)  :: var
      INTEGER, INTENT(in)           :: n1
      INTEGER, INTENT(in)           :: n2
      INTEGER, INTENT(out)          :: ier
      INTEGER, INTENT(in), OPTIONAL :: INTVAR(n1,n2)
      LOGICAL, INTENT(in), OPTIONAL :: BOOVAR(n1,n2)
      REAL, INTENT(in), OPTIONAL    :: FLTVAR(n1,n2)
      DOUBLE PRECISION, INTENT(in), OPTIONAL  :: DBLVAR(n1,n2)
      CHARACTER(LEN=*), INTENT(in), OPTIONAL :: ATT
      CHARACTER(LEN=*), INTENT(in), OPTIONAL :: ATT_NAME
      INTEGER        :: drank = 2
      INTEGER        :: arank = 1
      INTEGER        :: boo_temp(n1,n2)
      INTEGER(HID_T) :: dset_id
      INTEGER(HID_T) :: attr_id  
      INTEGER(HID_T) :: dspace_id
      INTEGER(HID_T) :: aspace_id
      INTEGER(HID_T) :: atype_id
      INTEGER(SIZE_T) :: attrlen
      INTEGER(HSIZE_T), DIMENSION(2) :: ddims
      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/)
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
      CHARACTER(LEN=1024) ::  att_temp
      
      ier = 0
      ddims(1)=n1
      ddims(2)=n2
      CALL h5screate_simple_f(drank,ddims,dspace_id,ier)
      IF (ier /=0) RETURN
      IF (PRESENT(INTVAR)) THEN
         CALL h5dcreate_f(file_id,TRIM(var),H5T_NATIVE_INTEGER,dspace_id,dset_id,ier)
         CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, INTVAR, ddims, ier)
      ELSE IF (PRESENT(FLTVAR)) THEN
         CALL h5dcreate_f(file_id,TRIM(var),H5T_NATIVE_DOUBLE,dspace_id,dset_id,ier)
         CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, FLTVAR, ddims, ier)
      ELSE IF (PRESENT(DBLVAR)) THEN
         CALL h5dcreate_f(file_id,TRIM(var),H5T_NATIVE_DOUBLE,dspace_id,dset_id,ier)
         CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, DBLVAR, ddims, ier)
      ELSE IF (PRESENT(BOOVAR)) THEN
         boo_temp = 0
         WHERE(BOOVAR) boo_temp = 1
         CALL h5dcreate_f(file_id,TRIM(var),H5T_NATIVE_INTEGER,dspace_id,dset_id,ier)
         CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, boo_temp, ddims, ier)
      ELSE
         ier=-2
      END IF
      IF (PRESENT(ATT)) THEN
         IF (PRESENT(ATT_NAME)) THEN
            att_temp=TRIM(ATT_NAME)
         ELSE
            att_temp=TRIM('_att')
         END IF
         data_dims(1) = 1
         CALL h5screate_simple_f(arank,adims,aspace_id,ier)
         IF (ier /=0) RETURN
         CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id,ier)
         IF (ier /=0) RETURN
         attrlen=LEN(TRIM(ATT))
         CALL h5tset_size_f(atype_id,attrlen,ier)
         IF (ier /=0) RETURN
         CALL h5acreate_f(dset_id,TRIM(var) // '_' // TRIM(att_temp), atype_id, aspace_id, attr_id,ier)
         IF (ier /=0) RETURN
         data_dims(1)=1
         CALL h5awrite_f(attr_id,atype_id,TRIM(ATT),data_dims,ier)
         IF (ier /=0) RETURN
         CALL h5aclose_f(attr_id,ier)
         IF (ier /=0) RETURN
         CALL h5sclose_f(aspace_id,ier)
         IF (ier /=0) RETURN
      END IF
      IF (ier /=0) RETURN
      CALL h5dclose_f(dset_id,ier)
      IF (ier /=0) RETURN
      CALL h5sclose_f(dspace_id,ier)
      IF (ier /=0) RETURN
      
      END SUBROUTINE write_arr2d_hdf5   
      !-----------------------------------------------------------------
      
      !-----------------------------------------------------------------
      SUBROUTINE write_arr3d_hdf5(file_id,var,n1,n2,n3,ier,BOOVAR,INTVAR,FLTVAR,DBLVAR,ATT,ATT_NAME)
      IMPLICIT NONE
      INTEGER(HID_T), INTENT(in)    :: file_id
      CHARACTER(LEN=*), INTENT(in)  :: var
      INTEGER, INTENT(in)           :: n1
      INTEGER, INTENT(in)           :: n2
      INTEGER, INTENT(in)           :: n3
      INTEGER, INTENT(out)          :: ier
      INTEGER, INTENT(in), OPTIONAL :: INTVAR(n1,n2,n3)
      LOGICAL, INTENT(in), OPTIONAL :: BOOVAR(n1,n2,n3)
      REAL, INTENT(in), OPTIONAL    :: FLTVAR(n1,n2,n3)
      DOUBLE PRECISION, INTENT(in), OPTIONAL  :: DBLVAR(n1,n2,n3)
      CHARACTER(LEN=*), INTENT(in), OPTIONAL :: ATT
      CHARACTER(LEN=*), INTENT(in), OPTIONAL :: ATT_NAME
      INTEGER        :: drank = 3
      INTEGER        :: arank = 1
      INTEGER        :: boo_temp(n1,n2,n3)
      INTEGER(HID_T) :: dset_id
      INTEGER(HID_T) :: attr_id  
      INTEGER(HID_T) :: dspace_id
      INTEGER(HID_T) :: aspace_id
      INTEGER(HID_T) :: atype_id
      INTEGER(SIZE_T) :: attrlen
      INTEGER(HSIZE_T), DIMENSION(3) :: ddims
      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/)
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
      CHARACTER(LEN=1024) ::  att_temp
      
      ier = 0
      ddims(1)=n1
      ddims(2)=n2
      ddims(3)=n3
      CALL h5screate_simple_f(drank,ddims,dspace_id,ier)
      IF (ier /=0) RETURN
      IF (PRESENT(INTVAR)) THEN
         CALL h5dcreate_f(file_id,TRIM(var),H5T_NATIVE_INTEGER,dspace_id,dset_id,ier)
         CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, INTVAR, ddims, ier)
      ELSE IF (PRESENT(FLTVAR)) THEN
         CALL h5dcreate_f(file_id,TRIM(var),H5T_NATIVE_DOUBLE,dspace_id,dset_id,ier)
         CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, FLTVAR, ddims, ier)
      ELSE IF (PRESENT(DBLVAR)) THEN
         CALL h5dcreate_f(file_id,TRIM(var),H5T_NATIVE_DOUBLE,dspace_id,dset_id,ier)
         CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, DBLVAR, ddims, ier)
      ELSE IF (PRESENT(BOOVAR)) THEN
         boo_temp = 0
         WHERE(BOOVAR) boo_temp = 1
         CALL h5dcreate_f(file_id,TRIM(var),H5T_NATIVE_INTEGER,dspace_id,dset_id,ier)
         CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, boo_temp, ddims, ier)
      ELSE
         ier=-2
      END IF
      IF (PRESENT(ATT)) THEN
         IF (PRESENT(ATT_NAME)) THEN
            att_temp=TRIM(ATT_NAME)
         ELSE
            att_temp=TRIM('_att')
         END IF
         data_dims(1) = 1
         CALL h5screate_simple_f(arank,adims,aspace_id,ier)
         IF (ier /=0) RETURN
         CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id,ier)
         IF (ier /=0) RETURN
         attrlen=LEN(TRIM(ATT))
         CALL h5tset_size_f(atype_id,attrlen,ier)
         IF (ier /=0) RETURN
         CALL h5acreate_f(dset_id,TRIM(var) // '_' // TRIM(att_temp), atype_id, aspace_id, attr_id,ier)
         IF (ier /=0) RETURN
         data_dims(1)=1
         CALL h5awrite_f(attr_id,atype_id,TRIM(ATT),data_dims,ier)
         IF (ier /=0) RETURN
         CALL h5aclose_f(attr_id,ier)
         IF (ier /=0) RETURN
         CALL h5sclose_f(aspace_id,ier)
         IF (ier /=0) RETURN
      END IF
      IF (ier /=0) RETURN
      CALL h5dclose_f(dset_id,ier)
      IF (ier /=0) RETURN
      CALL h5sclose_f(dspace_id,ier)
      IF (ier /=0) RETURN
      
      END SUBROUTINE write_arr3d_hdf5   
      !-----------------------------------------------------------------
      
      !-----------------------------------------------------------------
      SUBROUTINE read_scalar_hdf5(file_id,var,ierr,BOOVAR,INTVAR,FLTVAR,DBLVAR,ATT,ATT_NAME)
      IMPLICIT NONE
      INTEGER(HID_T), INTENT(in)    :: file_id
      CHARACTER(LEN=*), INTENT(in)  :: var
      INTEGER, INTENT(out)          :: ierr
      LOGICAL, INTENT(out), OPTIONAL :: BOOVAR
      INTEGER, INTENT(out), OPTIONAL :: INTVAR
      REAL, INTENT(out), OPTIONAL    :: FLTVAR   
      DOUBLE PRECISION, INTENT(out), OPTIONAL    :: DBLVAR  
      CHARACTER(LEN=*), INTENT(out), OPTIONAL :: ATT
      CHARACTER(LEN=*), INTENT(out), OPTIONAL :: ATT_NAME
      INTEGER        :: boo_temp
      INTEGER        :: drank = 1
      INTEGER        :: arank = 1
      INTEGER(HID_T) :: dset_id
      INTEGER(HID_T) :: attr_id  
      INTEGER(HID_T) :: dspace_id
      INTEGER(HID_T) :: aspace_id
      INTEGER(HID_T) :: atype_id
      INTEGER(SIZE_T) :: attrlen
      INTEGER(HSIZE_T), DIMENSION(1) :: ddims = (/1/)
      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/)
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
      CHARACTER(LEN=1024) ::  att_temp
      
      ierr = 1
      CALL h5dopen_f(file_id,TRIM(var), dset_id,ierr)
      IF (ierr /=0) RETURN
      IF (PRESENT(INTVAR)) THEN
         INTVAR = 0
         CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, INTVAR, ddims, ierr)
      ELSE IF (PRESENT(FLTVAR)) THEN
         FLTVAR = 0
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, FLTVAR, ddims, ierr)
      ELSE IF (PRESENT(DBLVAR)) THEN
         DBLVAR = 0
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, DBLVAR, ddims, ierr)
      ELSE IF (PRESENT(BOOVAR)) THEN
         boo_temp = 0
         CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, boo_temp, ddims, ierr)
         BOOVAR = .FALSE.
         IF (boo_temp == 1) BOOVAR = .TRUE.
      ELSE
         ierr=-2
      END IF
      CALL h5dclose_f(dset_id,ierr)
      IF (ierr /=0) RETURN
      END SUBROUTINE read_scalar_hdf5   
      !-----------------------------------------------------------------
      
      !-----------------------------------------------------------------
      SUBROUTINE read_vector_hdf5(file_id,var,n,ier,BOOVAR,INTVAR,FLTVAR,DBLVAR,ATT,ATT_NAME)
      IMPLICIT NONE
      INTEGER(HID_T), INTENT(in)    :: file_id
      CHARACTER(LEN=*), INTENT(in)  :: var
      INTEGER, INTENT(in)           :: n
      INTEGER, INTENT(out)          :: ier
      INTEGER, INTENT(out), OPTIONAL :: INTVAR(n)
      LOGICAL, INTENT(out), OPTIONAL :: BOOVAR(n)
      REAL, INTENT(out), OPTIONAL    :: FLTVAR(n)
      DOUBLE PRECISION, INTENT(out), OPTIONAL    :: DBLVAR(n)
      CHARACTER(LEN=*), INTENT(out), OPTIONAL :: ATT
      CHARACTER(LEN=*), INTENT(out), OPTIONAL :: ATT_NAME
      INTEGER        :: drank = 1
      INTEGER        :: arank = 1
      INTEGER        :: boo_temp(n)
      INTEGER(HID_T) :: dset_id
      INTEGER(HID_T) :: attr_id  
      INTEGER(HID_T) :: dspace_id
      INTEGER(HID_T) :: aspace_id
      INTEGER(HID_T) :: atype_id
      INTEGER(SIZE_T) :: attrlen
      INTEGER(HSIZE_T), DIMENSION(1) :: ddims
      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/)
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
      CHARACTER(LEN=1024) ::  att_temp
      
      ier = 0
      ddims(1)=n
      !ATT=''
      !ATT_NAME=''
      CALL h5dopen_f(file_id,TRIM(var), dset_id,ier)
      IF (ier /=0) RETURN
      IF (PRESENT(INTVAR)) THEN
         INTVAR = 0
         CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, INTVAR, ddims, ier)
      ELSE IF (PRESENT(FLTVAR)) THEN
         FLTVAR = 0
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, FLTVAR, ddims, ier)
      ELSE IF (PRESENT(DBLVAR)) THEN
         DBLVAR = 0
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, DBLVAR, ddims, ier)
      ELSE IF (PRESENT(BOOVAR)) THEN
         boo_temp = 0
         CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, boo_temp, ddims, ier)
         BOOVAR = .FALSE.
         WHERE(boo_temp == 1) BOOVAR = .TRUE.
      ELSE
         ier=-2
      END IF
      CALL h5dclose_f(dset_id,ier)
      IF (ier /=0) RETURN
      !CALL h5sclose_f(dspace_id,ier)
      !IF (ier /=0) RETURN
      
      END SUBROUTINE read_vector_hdf5   
      !-----------------------------------------------------------------
      
      !-----------------------------------------------------------------
      SUBROUTINE read_arr2d_hdf5(file_id,var,n1,n2,ier,BOOVAR,INTVAR,FLTVAR,DBLVAR,ATT,ATT_NAME)
      IMPLICIT NONE
      INTEGER(HID_T), INTENT(in)    :: file_id
      CHARACTER(LEN=*), INTENT(in)  :: var
      INTEGER, INTENT(in)           :: n1
      INTEGER, INTENT(in)           :: n2
      INTEGER, INTENT(out)          :: ier
      INTEGER, INTENT(out), OPTIONAL :: INTVAR(n1,n2)
      LOGICAL, INTENT(out), OPTIONAL :: BOOVAR(n1,n2)
      REAL, INTENT(out), OPTIONAL    :: FLTVAR(n1,n2)
      DOUBLE PRECISION, INTENT(out), OPTIONAL  :: DBLVAR(n1,n2)
      CHARACTER(LEN=*), INTENT(out), OPTIONAL :: ATT
      CHARACTER(LEN=*), INTENT(out), OPTIONAL :: ATT_NAME
      INTEGER        :: drank = 2
      INTEGER        :: arank = 1
      INTEGER        :: boo_temp(n1,n2)
      INTEGER(HID_T) :: dset_id
      INTEGER(HID_T) :: attr_id  
      INTEGER(HID_T) :: dspace_id
      INTEGER(HID_T) :: aspace_id
      INTEGER(HID_T) :: atype_id
      INTEGER(SIZE_T) :: attrlen
      INTEGER(HSIZE_T), DIMENSION(2) :: ddims
      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/)
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
      CHARACTER(LEN=1024) ::  att_temp
      
      ier = 0
      ddims(1)=n1
      ddims(2)=n2
      !ATT=''
      !ATT_NAME=''
      CALL h5dopen_f(file_id,TRIM(var), dset_id,ier)
      IF (ier /=0) RETURN
      IF (PRESENT(INTVAR)) THEN
         INTVAR = 0
         CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, INTVAR, ddims, ier)
      ELSE IF (PRESENT(FLTVAR)) THEN
         FLTVAR = 0
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, FLTVAR, ddims, ier)
      ELSE IF (PRESENT(DBLVAR)) THEN
         DBLVAR = 0
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, DBLVAR, ddims, ier)
      ELSE IF (PRESENT(BOOVAR)) THEN
         boo_temp = 0
         CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, boo_temp, ddims, ier)
         BOOVAR = .FALSE.
         WHERE(boo_temp == 1) BOOVAR = .TRUE.
      ELSE
         ier=-2
      END IF
      IF (ier /=0) RETURN
      CALL h5dclose_f(dset_id,ier)
      IF (ier /=0) RETURN
      !CALL h5sclose_f(dspace_id,ier)
      !IF (ier /=0) RETURN
      
      END SUBROUTINE read_arr2d_hdf5   
      !-----------------------------------------------------------------
      
      !-----------------------------------------------------------------
      SUBROUTINE read_arr3d_hdf5(file_id,var,n1,n2,n3,ier,BOOVAR,INTVAR,FLTVAR,DBLVAR,ATT,ATT_NAME)
      IMPLICIT NONE
      INTEGER(HID_T), INTENT(in)    :: file_id
      CHARACTER(LEN=*), INTENT(in)  :: var
      INTEGER, INTENT(in)           :: n1
      INTEGER, INTENT(in)           :: n2
      INTEGER, INTENT(in)           :: n3
      INTEGER, INTENT(out)          :: ier
      INTEGER, INTENT(out), OPTIONAL :: INTVAR(n1,n2,n3)
      LOGICAL, INTENT(out), OPTIONAL :: BOOVAR(n1,n2,n3)
      REAL, INTENT(out), OPTIONAL    :: FLTVAR(n1,n2,n3)
      DOUBLE PRECISION, INTENT(out), OPTIONAL  :: DBLVAR(n1,n2,n3)
      CHARACTER(LEN=*), INTENT(out), OPTIONAL :: ATT
      CHARACTER(LEN=*), INTENT(out), OPTIONAL :: ATT_NAME
      INTEGER        :: drank = 3
      INTEGER        :: arank = 1
      INTEGER        :: boo_temp(n1,n2,n3)
      INTEGER(HID_T) :: dset_id
      INTEGER(HID_T) :: attr_id  
      INTEGER(HID_T) :: dspace_id
      INTEGER(HID_T) :: aspace_id
      INTEGER(HID_T) :: atype_id
      INTEGER(SIZE_T) :: attrlen
      INTEGER(HSIZE_T), DIMENSION(3) :: ddims
      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/)
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
      CHARACTER(LEN=1024) ::  att_temp
      
      ier = 0
      ddims(1)=n1
      ddims(2)=n2
      ddims(3)=n3
      !ATT=''
      !ATT_NAME=''
      CALL h5dopen_f(file_id,TRIM(var), dset_id,ier)
      IF (ier /=0) RETURN
      IF (PRESENT(INTVAR)) THEN
         INTVAR = 0
         CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, INTVAR, ddims, ier)
      ELSE IF (PRESENT(FLTVAR)) THEN
         FLTVAR = 0
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, FLTVAR, ddims, ier)
      ELSE IF (PRESENT(DBLVAR)) THEN
         DBLVAR = 0
         CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, DBLVAR, ddims, ier)
      ELSE IF (PRESENT(BOOVAR)) THEN
         boo_temp = 0
         CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, boo_temp, ddims, ier)
         BOOVAR = .FALSE.
         WHERE(boo_temp == 1) BOOVAR = .TRUE.
      ELSE
         ier=-2
      END IF
      IF (ier /=0) RETURN
      CALL h5dclose_f(dset_id,ier)
      IF (ier /=0) RETURN
      !CALL h5sclose_f(dspace_id,ier)
      !IF (ier /=0) RETURN
      
      END SUBROUTINE read_arr3d_hdf5   
      !-----------------------------------------------------------------
         
#endif
!-----------------------------------------------------------------------
!     End Module
!-----------------------------------------------------------------------
      END MODULE ez_hdf5
!
!
!
!-----  B E G I N  S A M P L E  C O D E  -------------------------------
!      SUBROUTINE WRITE_HDF5_TEST
!      USE ez_hdf5
!      IMPLICIT NONE
!      INTEGER :: ier
!      REAL    :: test_arr(4,6)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
!      test_arr(:,:) = 1.56
!      !-----  OPEN HDF5 FILE
!      PRINT *,'WRITING TO FILE: ','test.h5'
!      CALL open_hdf5('test.h5',fid,ier)
!      !-----  WRITE SCALARS
!      CALL write_var_hdf5(fid,'int_var',ier,INTVAR=314,ATT='An Integer Variable',ATT_NAME='Description')
!      CALL write_var_hdf5(fid,'flt_var',ier,FLTVAR=3.14159,ATT='Approximation of Pi')
!      CALL write_var_hdf5(fid,'boo_var',ier,BOOVAR=.TRUE.,ATT='A Boolean Variable')
!      !-----  WRITE VECTORS
!      CALL write_var_hdf5(fid,'int_arr',5,ier,INTVAR=/0,1,2,3,4/)
!      !-----  WRITE ARRAYS
!      CALL write_var_hdf5(fid,'rmnc',4,6,ier,FLTVAR=test_arr,ATT='4x6 Test Array')
!      !-----  CLOSE HDF5 FILE
!      CALL close_hdf5(fid,ier)
!      END SUBROUTINE WRITE_HDF5_TEST
!-----------------------------------------------------------------------
