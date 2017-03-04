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
      USE HDF5
!-----------------------------------------------------------------------
!     Module Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER(HID_T) :: fid
!-----------------------------------------------------------------------
!     Module Subroutines
!          open_hdf    Opens and hdf5 file for input/output
!          write_scalar_hdf5    Outputs a scalar to the HDF5 file.
!          write_vector_hdf5    Outputs a vector to the HDF5 file.
!          write_arr2d_hdf5     Outputs a 2D array to the HDF5 file.
!-----------------------------------------------------------------------
      CONTAINS
      
      !-----------------------------------------------------------------
      SUBROUTINE open_hdf5(filename,file_id,ier)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER(HID_T), INTENT(out)  :: file_id
      INTEGER, INTENT(out)         :: ier
      INTEGER             ::  strlen, dex
      CHARACTER(LEN=1024) ::  file_temp
      file_id   = -1
      file_temp = TRIM(filename)
      strlen    = LEN(TRIM(file_temp))
      dex=INDEX(file_temp,'.h5',BACK=.TRUE.)
      IF (dex /= strlen-2) file_temp = TRIM(file_temp) // '.h5'
      CALL h5open_f(ier)
      IF (ier /=0) RETURN
      CALL h5fcreate_f(TRIM(file_temp),H5F_ACC_TRUNC_F, file_id,ier)
      END SUBROUTINE open_hdf5
      !-----------------------------------------------------------------
         
      !-----------------------------------------------------------------
      SUBROUTINE close_hdf5(file_id,ier)
      IMPLICIT NONE
      INTEGER(HID_T), INTENT(out)  :: file_id
      INTEGER, INTENT(out)         :: ier
      CALL h5fclose_f(fid,ier)
      IF (ier /=0) RETURN
      CALL h5close_f(ier)
      IF (ier /=0) RETURN
      END SUBROUTINE close_hdf5
      !-----------------------------------------------------------------
      
      !-----------------------------------------------------------------
      SUBROUTINE write_scalar_hdf5(file_id,var,ier,BOOVAR,INTVAR,FLTVAR,DBLVAR,ATT,ATT_NAME)
      IMPLICIT NONE
      INTEGER(HID_T), INTENT(in)    :: file_id
      CHARACTER(LEN=*), INTENT(in)  :: var
      INTEGER, INTENT(out)          :: ier
      LOGICAL, INTENT(in), OPTIONAL :: BOOVAR
      INTEGER, INTENT(in), OPTIONAL :: INTVAR
      REAL, INTENT(in), OPTIONAL    :: FLTVAR   
      DOUBLE PRECISION, INTENT(in), OPTIONAL    :: DBLVAR  
      CHARACTER(LEN=*), INTENT(in), OPTIONAL :: ATT
      CHARACTER(LEN=*), INTENT(in), OPTIONAL :: ATT_NAME
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
         IF (BOOVAR) boo_temp = 1
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
      
      ddims(1)=n
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
!      CALL write_scalar_hdf5(fid,'int_var',ier,INTVAR=314,ATT='An Integer Variable',ATT_NAME='Description')
!      CALL write_scalar_hdf5(fid,'flt_var',ier,FLTVAR=3.14159,ATT='Approximation of Pi')
!      CALL write_scalar_hdf5(fid,'boo_var',ier,BOOVAR=.TRUE.,ATT='A Boolean Variable')
!      !-----  WRITE VECTORS
!      CALL write_vector_hdf5(fid,'int_arr',5,ier,INTVAR=/0,1,2,3,4/)
!      !-----  WRITE ARRAYS
!      CALL write_arr2d_hdf5(fid,'rmnc',4,6,ier,FLTVAR=test_arr,ATT='4x6 Test Array')
!      !-----  CLOSE HDF5 FILE
!      CALL close_hdf5(fid,ier)
!      END SUBROUTINE WRITE_HDF5_TEST
!-----------------------------------------------------------------------
