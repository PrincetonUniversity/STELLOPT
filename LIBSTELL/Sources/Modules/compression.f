!*******************************************************************************
!>  @file compression.f
!>  @brief Contains module @ref compression.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref compression_class. This class
!>  contains the code and buffers to hold compressed and uncompressed data. 2D
!>  matrix data is compressed using a singular value decomposition. See
!>  del-Castillo-Negrete et. al. doi:10.1016/j.jcp.2006.07.022.
!*******************************************************************************

      MODULE compression
      USE stel_kinds
      USE mpi_inc

      IMPLICIT NONE

!*******************************************************************************
!  reconstruction module parameters
!*******************************************************************************
!>  Maximum compressed size before uncompressed buffers are stored.
      REAL (rprec), PARAMETER :: compression_max_percentage = 100.0

!>  Postfix for the data buffer.
      CHARACTER (len=*), PARAMETER ::                                          &
     &   compression_data_buffer_post = "data_buffer"
!>  Postfix for the data buffer.
      CHARACTER (len=*), PARAMETER ::                                          &
     &   compression_u_buffer_post = "u_buffer"
!>  Postfix for the data buffer.
      CHARACTER (len=*), PARAMETER ::                                          &
     &   compression_wvt_buffer_post = "wvt_buffer"
!>  Postfix for the data buffer.
      CHARACTER (len=*), PARAMETER ::                                          &
     &   compression_data_dim_post = "data_dim"

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) compression base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class containing buffers for compressed and uncompressed data.
!-------------------------------------------------------------------------------
      TYPE compression_class
!  Uncompressed buffer
         REAL (rprec), DIMENSION(:,:), POINTER :: data_buffer => null()

!  Singular Value U buffer
         REAL (rprec), DIMENSION(:,:), POINTER :: u_buffer => null()
!  Singular Value V transpose buffer
         REAL (rprec), DIMENSION(:,:), POINTER :: wvt_buffer => null()

!  Dimensions of the data buffer. This is used in the special case that the
!  orginal data buffer is all zeros. Otherwise it remains unallocated.
         INTEGER, DIMENSION(:), POINTER        :: data_dim => null()
      END TYPE

!-------------------------------------------------------------------------------
!>  Pointer to a compression object. Used for creating arrays of compression
!>  pointers. This is needed because fortran does not allow arrays of pointers
!>  directly.
!-------------------------------------------------------------------------------
      TYPE compression_pointer
         TYPE (compression_class), POINTER :: p => null()
      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
      INTERFACE compression_construct
         MODULE PROCEDURE compression_construct_new,                           &
     &                    compression_construct_netcdf
      END INTERFACE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref compression_class object.
!>
!>  Constructs compression object and compresses the data buffer if compression
!>  meets the specified criteria.
!>
!>  @param[in] data_in     Data buffer to compress. This gets copied to the
!>                         internal storage buffer.
!>  @param[in] svn_cut_off Cutoff value for determining the number of singular
!>                         values to use for compression.
!>  @returns A pointer to a constructed @ref compression_class object.
!-------------------------------------------------------------------------------
      FUNCTION compression_construct_new(data_in, svd_cut_off)
      USE v3_utilities

      IMPLICIT NONE

!  Declare Arguments
      TYPE (compression_class), POINTER :: compression_construct_new
      REAL (rprec), DIMENSION(:,:), INTENT(in) :: data_in
      REAL (rprec), INTENT(in)                 :: svd_cut_off

!  local arguments
      INTEGER                                  :: m, n, i
      REAL (rprec), DIMENSION(:), ALLOCATABLE  :: w_svd
      REAL (rprec), DIMENSION(:), ALLOCATABLE  :: svd_work
      INTEGER                                  :: work_size
      INTEGER                                  :: status
      REAL (rprec)                             :: e_zero
      REAL (rprec)                             :: e_r

!  Start of executable code
      ALLOCATE(compression_construct_new)

!  Special case if all the data is zero, there is no need to retain any singular
!  values. This will be detected in the uncompress by having the data_dim array
!  allocated.
      IF (ALL(data_in .eq. 0.0)) THEN
         ALLOCATE(compression_construct_new%data_dim(2))
         compression_construct_new%data_dim = SHAPE(data_in)
         RETURN
      END IF

      m = SIZE(data_in, 1)
      n = SIZE(data_in, 2)
      ALLOCATE(compression_construct_new%data_buffer(m,n))
      compression_construct_new%data_buffer = data_in

      IF (svd_cut_off .eq. 0.0) THEN
         RETURN
      END iF

      ALLOCATE(compression_construct_new%wvt_buffer(n,n))
      ALLOCATE(compression_construct_new%u_buffer(m,m))
      ALLOCATE(svd_work(1))
      ALLOCATE(w_svd(MIN(m,n)))

      compression_construct_new%wvt_buffer = 0.0
      compression_construct_new%u_buffer = 0.0
      w_svd = 0.0
      svd_work = 0.0

!  Find the optimal work size.
      CALL dgesvd('All', 'All', m, n,                                          &
     &            compression_construct_new%data_buffer, m,                    &
     &            w_svd, compression_construct_new%u_buffer, m,                &
     &            compression_construct_new%wvt_buffer, n, svd_work, -1,       &
     &            status)
      work_size = INT(svd_work(1))
      DEALLOCATE(svd_work)
      ALLOCATE(svd_work(work_size))
      svd_work = 0.0

!  Factor the matrix to M = U * W * V^T
      CALL dgesvd('All', 'All', m, n,                                          &
     &            compression_construct_new%data_buffer, m,                    &
     &            w_svd, compression_construct_new%u_buffer, m,                &
     &            compression_construct_new%wvt_buffer, n, svd_work,           &
     &            work_size, status)
      CALL assert_eq(0, status, 'dgesvd problem when compressing ' //          &
     &                          'buffer')

      DEALLOCATE(svd_work)

!  Determine the number of singular values to keep. If the number of singular
!  values is greater than half, do not compress the buffer. First start by
!  setting any negative values to zero.
      WHERE (w_svd .lt. 0.0)
         w_svd = 0.0
      END WHERE

!  Find the minimum number of singular values to keep based on
!  1 - n(r) = svd_cut_off where n(r) is equation (7) in del-Castillo-Negrete
!  et. al. doi:10.1016/j.jcp.2006.07.022
      e_zero = DOT_PRODUCT(w_svd, w_svd)

      e_r = 0.0
!  We only need to check up to one more than half the number of singular values.
!  Anymore than that will use more data storage. Store the number of singular
!  values retained in the work_size.
      DO i = 1, SIZE(w_svd)/2 + 1
         e_r = e_r + w_svd(i)*w_svd(i)
         work_size = i
         IF ((1.0 - e_r/e_zero) .le. svd_cut_off) THEN
            EXIT
         END IF
      END DO

!  Check the compression percentage. The compression is defined as the
!  percentage of the orginal data size. The compression percentage is defined as
!
!    compressed_size/uncompressed_size*100.0                                 (1)
!
!  where
!
!    compressed_size = r + r*m + r*n                                         (2)
!
!  and
!
!    uncompressed_size = n*m                                                 (3)
!
!  n and m are the dimensions of the orginal matrix and r is the number of
!  singular values retained.
!
!  However you can gain even more storage space. The compressed data is '
!  uncompressed reconstructing the matrix from U.W.V^T. By precomputing W.V^T
!  and storing that buffer instead instead of the W and V^T buffers separately,
!  the an array the size of the number of signular values is eliminiated. This
!  also has the added benifit of speed when uncompressing the data. The new
!  compressed size becomes
!
!    compressed_size = r*m + r*n                                             (4)
!
!  Barrow the e_ variables to compute th compression
!  percentage. Use e_r for the compressed size and e_zero of the uncompressed
!  size.

      e_r = work_size*m + work_size*n
      e_zero = n*m

      IF (100.0*e_r/e_zero .ge. compression_max_percentage) THEN
!  Use the uncompressed buffers.
         compression_construct_new%data_buffer = data_in

         DEALLOCATE(compression_construct_new%u_buffer)
         compression_construct_new%u_buffer => null()

         DEALLOCATE(compression_construct_new%wvt_buffer)
         compression_construct_new%wvt_buffer => null()
      ELSE
!  Use the compressed buffers. Barrow the data_buffer as temp buffer for
!  resizing the storage buffers.
         DEALLOCATE(compression_construct_new%data_buffer)

!  Store the U buffer. The reduced U buffer now only needs m x r storage space.
         ALLOCATE(compression_construct_new%data_buffer(m,work_size))
         compression_construct_new%data_buffer =                               &
     &      compression_construct_new%u_buffer(:,1:work_size)
         DEALLOCATE(compression_construct_new%u_buffer)
         ALLOCATE(compression_construct_new%u_buffer(m,work_size))
         compression_construct_new%u_buffer =                                  &
     &      compression_construct_new%data_buffer
         DEALLOCATE(compression_construct_new%data_buffer)

!  Store the VT buffer. The reduced VT buffer now only needs r x n storage
!  space.
         ALLOCATE(compression_construct_new%data_buffer(work_size,n))
         compression_construct_new%data_buffer =                               &
     &      compression_construct_new%wvt_buffer(1:work_size,:)
         DEALLOCATE(compression_construct_new%wvt_buffer)
         ALLOCATE(compression_construct_new%wvt_buffer(work_size,n))

!  Multiply W.VT and store in the WVT buffer.
         DO i = 1, n
            compression_construct_new%wvt_buffer(:,i) =                        &
     &         w_svd(1:work_size) *                                            &
     &         compression_construct_new%data_buffer(:,i)
         END DO

         DEALLOCATE(compression_construct_new%data_buffer)
         compression_construct_new%data_buffer => null()
      END IF

      DEALLOCATE(w_svd)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref compression_class object by reading from an netcdf.
!>  file.
!>
!>  Constructs data from a precompressed netcdf file.
!>
!>  @param[in] ncid A netcdf file id.
!>  @param[in] name A netcdf variable name id.
!>  @param[in] name Flag discribing how the data was compressed.
!>  @returns A pointer to a constructed @ref compression_class object.
!-------------------------------------------------------------------------------
      FUNCTION compression_construct_netcdf(ncid, name)
      USE ezcdf

      IMPLICIT NONE

!  Declare Arguments
      TYPE (compression_class), POINTER :: compression_construct_netcdf
      INTEGER, INTENT(in)               :: ncid
      CHARACTER (len=*), INTENT(in)     :: name

!  local variables
      INTEGER, DIMENSION(2)             :: dim_lengths
      INTEGER                           :: error

!  Start of executable code
      ALLOCATE(compression_construct_netcdf)

      dim_lengths = -1

      CALL cdf_inquire(ncid, TRIM(name) // compression_data_buffer_post,       &
     &                 dim_lengths)
      IF (ALL(dim_lengths .gt. 0)) THEN
         ALLOCATE(                                                             &
     &      compression_construct_netcdf%data_buffer(dim_lengths(1),           &
     &                                               dim_lengths(2)))
         CALL cdf_read(ncid,                                                   &
     &                 TRIM(name) // compression_data_buffer_post,             &
     &                 compression_construct_netcdf%data_buffer)
         RETURN
      END IF

      CALL cdf_inquire(ncid, TRIM(name) // compression_u_buffer_post,          &
     &                 dim_lengths)
      IF (ALL(dim_lengths .gt. 0)) THEN
         ALLOCATE(                                                             &
     &      compression_construct_netcdf%u_buffer(dim_lengths(1),              &
     &                                            dim_lengths(2)))
         CALL cdf_read(ncid, TRIM(name) // compression_u_buffer_post,          &
     &                 compression_construct_netcdf%u_buffer)

         CALL cdf_inquire(ncid,                                                &
     &                    TRIM(name) // compression_wvt_buffer_post,           &
     &                    dim_lengths)
         ALLOCATE(                                                             &
     &      compression_construct_netcdf%wvt_buffer(dim_lengths(1),            &
     &                                              dim_lengths(2)))
         CALL cdf_read(ncid, TRIM(name) // compression_wvt_buffer_post,        &
     &                 compression_construct_netcdf%wvt_buffer)
         RETURN
      END IF

      CALL cdf_inquire(ncid, TRIM(name) // compression_data_dim_post,          &
     &                 dim_lengths(1:1))
      ALLOCATE( compression_construct_netcdf%data_dim(dim_lengths(1)))
      CALL cdf_read(ncid, TRIM(name) // compression_data_dim_post,             &
     &              compression_construct_netcdf%data_dim)

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref compression_class object.
!>
!>  Deallocates memory and uninitializes a @ref compression_class object.
!>
!>  @param[inout] this A @ref compression_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE compression_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (compression_class), POINTER :: this

!  Start of executable code
      IF (ASSOCIATED(this%data_buffer)) THEN
         DEALLOCATE(this%data_buffer)
         this%data_buffer => null()
      END IF

      IF (ASSOCIATED(this%u_buffer)) THEN
         DEALLOCATE(this%u_buffer)
         this%u_buffer => null()
      END IF

      IF (ASSOCIATED(this%wvt_buffer)) THEN
         DEALLOCATE(this%wvt_buffer)
         this%wvt_buffer => null()
      END IF

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Get the ith dimension length.
!>
!>  Retrives the uncompressed dimension of the first data buffer index. This
!>  can be used without needing to decompress the data first.
!>
!>  @param[inout] this A @ref compression_class instance.
!>  @returns The size of the ith dimension of the uncompressed data buffer.
!-------------------------------------------------------------------------------
      FUNCTION compression_get_dimension1(this)

      IMPLICIT NONE

!  Declare Aruments
      INTEGER                           :: compression_get_dimension1
      TYPE (compression_class), POINTER :: this

!  Start of executable code
      IF (ASSOCIATED(this%data_buffer)) THEN
         compression_get_dimension1 = SIZE(this%data_buffer, 1)
      ELSE IF (ASSOCIATED(this%u_buffer)) THEN
         compression_get_dimension1 = SIZE(this%u_buffer, 1)
      ELSE IF (ASSOCIATED(this%data_dim)) THEN
         compression_get_dimension1 = this%data_dim(1)
      ELSE
         compression_get_dimension1 = 0
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get the jth dimension length.
!>
!>  Retrives the uncompressed dimension of the second data buffer index. This
!>  can be used without needing to decompress the data first.
!>
!>  @param[inout] this A @ref compression_class instance.
!>  @returns The size of the jth dimension of the uncompressed data buffer.
!-------------------------------------------------------------------------------
      FUNCTION compression_get_dimension2(this)

      IMPLICIT NONE

!  Declare Aruments
      INTEGER                           :: compression_get_dimension2
      TYPE (compression_class), POINTER :: this

!  Start of executable code
      IF (ASSOCIATED(this%data_buffer)) THEN
         compression_get_dimension2 = SIZE(this%data_buffer, 2)
      ELSE IF (ASSOCIATED(this%u_buffer)) THEN
         compression_get_dimension2 = SIZE(this%wvt_buffer, 2)
      ELSE IF (ASSOCIATED(this%data_dim)) THEN
         compression_get_dimension2 = this%data_dim(2)
      ELSE
         compression_get_dimension2 = 0
      END IF

      END FUNCTION

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Decompress the data.
!>
!>  Checks if the U and WVT arrays are allocated. If they are not the
!>  data_buffer was left uncompressed. If the data dimensions is also
!>  unallocated, then the orginal data buffer was all zeros otherwise data is
!>  uncompressed by matrix a multiplication.
!>
!>  @param[inout] this A @ref compression_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE compression_decompress(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (compression_class), POINTER :: this

!  Start of executable code
      IF (ASSOCIATED(this%data_buffer)) THEN
!  Data is stored uncompressed. No need to do anything.
         RETURN
      END IF

      IF (ASSOCIATED(this%data_dim)) THEN
!  Data is all zeros. Allocate the data_buffer and set to zero.
         ALLOCATE(this%data_buffer(this%data_dim(1),this%data_dim(2)))
         this%data_buffer = 0.0
         RETURN
      END IF

!  If control has reached this point, the data is compressed. Uncompress it by
!  U.WVT.
      ALLOCATE(this%data_buffer(SIZE(this%u_buffer, 1),                        &
     &                          SIZE(this%wvt_buffer, 2)))

      this%data_buffer = MATMUL(this%u_buffer, this%wvt_buffer)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Cleanup decompressed data.
!>
!>  If the data dimensions arrays or the U and WVT is allocated, then delete the
!>  data buffer. Otherwise the data is stored uncompressed.
!>
!>  @param[inout] this A @ref compression_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE compression_cleanup(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (compression_class), POINTER :: this

!  Start of executable code
      IF (ASSOCIATED(this%data_dim) .or. ASSOCIATED(this%u_buffer)) THEN
         DEALLOCATE(this%data_buffer)
         this%data_buffer => null()
      END IF

      END SUBROUTINE

!*******************************************************************************
!  NETCDF SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Define variables for the compressed data.
!>
!>  Defines the netcdf variables for writting compressed data to a netcdf file.
!>
!>  @param[in] this A @ref compression_class instance.
!>  @param[in] ncid Netcdf file id.
!>  @param[in] name Base name of variables to compress.
!>  @returns The flag indicating which buffer was stored.
!-------------------------------------------------------------------------------
      SUBROUTINE compression_define(this, ncid, name)
      USE ezcdf

      IMPLICIT NONE

!  Declare Arguments
      TYPE (compression_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: ncid
      CHARACTER (len=*), INTENT(in)        :: name

!  Start of executable code
      IF (ASSOCIATED(this%data_buffer)) THEN
         CALL cdf_define(ncid,                                                 &
     &                   TRIM(name) // compression_data_buffer_post,           &
     &                   this%data_buffer)
      ELSE IF (ASSOCIATED(this%u_buffer)) THEN
         CALL cdf_define(ncid, TRIM(name) // compression_u_buffer_post,        &
     &                   this%u_buffer)
         CALL cdf_define(ncid,                                                 &
     &                   TRIM(name) // compression_wvt_buffer_post,            &
     &                   this%wvt_buffer)
      ELSE
         CALL cdf_define(ncid, TRIM(name) // compression_data_dim_post,        &
     &                   this%data_dim)
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Cleanup decompressed data.
!>
!>  If the data dimensions arrays or the U and WVT is allocated, then delete the
!>  data buffer. Otherwise the data is stored uncompressed.
!>
!>  @param[inout] this A @ref compression_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE compression_write(this, ncid, name)
      USE ezcdf

      IMPLICIT NONE

!  Declare Arguments
      TYPE (compression_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: ncid
      CHARACTER (len=*), INTENT(in)        :: name

!  Start of executable code
      IF (ASSOCIATED(this%data_buffer)) THEN

         CALL cdf_write(ncid,                                                  &
     &                  TRIM(name) // compression_data_buffer_post,            &
     &                  this%data_buffer)

      ELSE IF (ASSOCIATED(this%u_buffer)) THEN

         CALL cdf_write(ncid, TRIM(name) // compression_u_buffer_post,         &
     &                  this%u_buffer)
         CALL cdf_write(ncid, TRIM(name) // compression_wvt_buffer_post,       &
     &                  this%wvt_buffer)

      ELSE

         CALL cdf_write(ncid, TRIM(name) // compression_data_dim_post,         &
     &                  this%data_dim)

      END IF

      END SUBROUTINE

      END MODULE
