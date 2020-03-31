!*******************************************************************************
!>  @file file_opts.f
!>   @brief Contains module @ref file_opts
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Contains cross platform routines for manipulating files on the file system.
!>  Defines a functions to move, copy and delete a file that is cross platform.
!>  While, system contains a crossplatform way to call a command line
!>  commands, the command used will vary by platform.
!*******************************************************************************
      MODULE file_opts
      USE SYSTEM_MOD

      IMPLICIT NONE

!*******************************************************************************
!  file_opts module parameters
!*******************************************************************************
!>  Length of file paths.
      INTEGER, PARAMETER :: path_length = 300

      CONTAINS

!-------------------------------------------------------------------------------
!>  @brief Moves the source file to the destination.
!>
!>  Moves a file by calling the appropriate command for the platform. If the
!>  destination is an existing file, that file will be overwritten.
!>
!>  @param[in]  file_source The file to move.
!>  @param[in]  file_dest   The destination of the moved file.
!>  @param[out] error       The result of the move.
!-------------------------------------------------------------------------------
      SUBROUTINE move_file(file_source, file_dest, error)
#if defined(__INTEL_COMPILER)
      USE IFLPORT
#endif

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=*), INTENT(in) :: file_source
      CHARACTER (len=*), INTENT(in) :: file_dest
      INTEGER, INTENT(out)          :: error

!  local parameters
!  To avoid any possible interactive command lines, make sure the command
!  overwrites the file if it exists.
#if defined(WIN32)
      CHARACTER (len=*), PARAMETER :: cmd = 'move /y '
#else
      CHARACTER (len=*), PARAMETER :: cmd = 'mv -f '
#endif

!  Start of executable code

#if defined(__GFORTRAN__) || defined(__INTEL_COMPILER)
      error = RENAME(TRIM(file_source), TRIM(file_dest))
#else
      CALL system(cmd // TRIM(file_source) // ' ' // file_dest,           &
     &            error)
#endif
      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Copies the source file to the destination.
!>
!>  Copies a file by calling the appropriate command for the platform. If the
!>  destination is an existing file, that file will be overwritten.
!>
!>  @param[in]  file_source The file to copy.
!>  @param[in]  file_dest   The destination of the copied file.
!>  @param[out] error       The result of the copy.
!-------------------------------------------------------------------------------
      SUBROUTINE copy_file(file_source, file_dest, error)
      USE safe_open_mod
#if defined(FAST_COPY)
      USE, INTRINSIC :: ISO_C_BINDING
#endif

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=*), INTENT(in)                :: file_source
      CHARACTER (len=*), INTENT(in)                :: file_dest
      INTEGER, INTENT(out)                         :: error

!  local variables
      LOGICAL                                      :: l_exists
      LOGICAL                                      :: l_opened
      INTEGER                                      :: file_size
      INTEGER                                      :: max_block_size
      INTEGER                                      :: io_source
      INTEGER                                      :: io_dest
      INTEGER                                      :: n_rec
      INTEGER                                      :: i_rec
      CHARACTER (len=1), ALLOCATABLE, DIMENSION(:) :: buffer

!  local parameters
      INTEGER, PARAMETER            :: block_size = 1024

#if defined(FAST_COPY)
      INTERFACE
         FUNCTION copy_file_c(src, dest) BIND(C)
         USE, INTRINSIC :: ISO_C_BINDING

         IMPLICIT NONE

         INTEGER, PARAMETER :: length = 300
         INTEGER (c_int)    :: copy_file_c
         CHARACTER (kind=c_char,len=1), DIMENSION(length), INTENT(in) ::       &
     &      src
         CHARACTER (kind=c_char,len=1), DIMENSION(length), INTENT(in) ::       &
     &      dest
         END FUNCTION
      END INTERFACE
#endif

!  Start of executable code
#if defined(FAST_COPY)
      error = copy_file_c(TRIM(file_source) // C_NULL_CHAR,                    &
     &                    TRIM(file_dest) // C_NULL_CHAR)
#else

!  Check source file.
      INQUIRE (file=file_source, exist=l_exists, opened=l_opened)
      CALL getfilesize(file_source, file_size)

      IF (.not.l_exists) THEN
         WRITE (*,1001) TRIM(file_source)
         error = -1
         RETURN
      ELSE IF (l_opened) THEN
         WRITE (*,1000) TRIM(file_source)
         error = -1
         RETURN
      END IF

!  Check destination file.
      INQUIRE (file=file_dest, exist=l_exists, opened=l_opened)
      IF (l_opened) THEN
         WRITE (*,1000) TRIM(file_dest)
         error = -1
         RETURN
      ELSE
         CALL delete_file(file_dest, error)
      END IF

      io_source = 10
      io_dest = 11

      max_block_size = MIN(block_size, file_size)

      CALL safe_open(io_source, error, file_source, 'old',                     &
     &               'unformatted', record_in=max_block_size,                  &
     &               access_in='direct')
      IF (error .ne. 0) WRITE (*,1002) TRIM(file_source)

      CALL safe_open(io_dest, error, file_dest, 'new', 'unformatted',          &
     &               record_in=max_block_size, access_in='direct')
      IF (error .ne. 0) WRITE (*,1003) TRIM(file_dest)

      ALLOCATE (buffer(max_block_size), stat=error)
      IF (error .ne. 0) WRITE (*,1004)

!  Write file in block sized chunks.
      n_rec = MAX(1, file_size/max_block_size)
      DO i_rec = 1, n_rec
         READ (io_source, rec=i_rec) buffer
         WRITE (io_dest, rec=i_rec) buffer
      END DO

      CLOSE (io_source)
      CLOSE (io_dest)

!  The block size may not bet an even multiple of the file size. Copy the
!  remaining part in 1 byte chunks.
      IF (file_size - (i_rec - 1)*max_block_size .gt. 0) THEN
         OPEN (unit=io_source, file=file_source, status='old',                 &
     &         access='direct', recl=1)
         OPEN (unit=io_dest, file=file_dest, status='old',                     &
     &         access='direct', recl=1)

         DO i_rec = i_rec, file_size
            READ (io_source, rec=i_rec) buffer(1)
            WRITE (io_dest, rec=i_rec) buffer(1)
         END DO

         CLOSE (io_source)
         CLOSE (io_dest)
      END IF

      DEALLOCATE (buffer)

1000  FORMAT('ERROR: File ',a,' is already opened.')
1001  FORMAT('ERROR: File ',a,' does not exist.')
1002  FORMAT('ERROR: Fail to open ',a,'.')
1003  FORMAT('ERROR: Fail to create ',a,'.')
1004  FORMAT('ERROR: Fail to allocate copy buffer.')
#endif
      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Deletes a file
!>
!>  Deletes a file by calling the appropriate command for the platform. Ignores
!>  errors if the file doesn't exisit.
!>
!>  @param[in]  file_source The file to delete.
!>  @param[out] error       The result of the copy.
!-------------------------------------------------------------------------------
      SUBROUTINE delete_file(file_source, error)
#if defined(__INTEL_COMPILER)
      USE IFLPORT
#endif
      USE safe_open_mod

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=*), INTENT(in) :: file_source
      INTEGER, INTENT(out)          :: error

!  local variables
      INTEGER                       :: io_unit

!  Start of executable code
#if defined(__GFORTRAN__) || defined(__INTEL_COMPILER)
      error = UNLINK(file_source)
#else
      io_unit = 10
      CALL safe_open(io_unit, error, file_source, 'old', 'formatted')
      CLOSE (io_unit, iostat=error, status='delete')
#endif

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Makes a directory
!>
!>  Creates new directory by calling the appropriate command for the plaform.
!>
!>  @param[in]  file_source The file to delete.
!>  @param[out] error       The result of the copy.
!-------------------------------------------------------------------------------
      SUBROUTINE create_directory(directory_source, error)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=*), INTENT(in) :: directory_source
      INTEGER, INTENT(out)          :: error

!  local parameters
#if defined(WIN32)
      CHARACTER (len=*), PARAMETER :: cmd = 'mkdir '
#else
      CHARACTER (len=*), PARAMETER :: cmd = 'mkdir -p '
#endif

!  Start of executable code

      CALL system(cmd // TRIM(directory_source), error)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Deletes a directory
!>
!>  Deletes directory by calling the appropriate command for the plaform. This
!>  deletes all files inside the directory as well.
!>
!>  @param[in]  file_source The file to delete.
!>  @param[out] error       The result of the copy.
!-------------------------------------------------------------------------------
      SUBROUTINE delete_directory(directory_source, error)
      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=*), INTENT(in) :: directory_source
      INTEGER, INTENT(out)          :: error

!  local parameters
#if defined(WIN32)
      CHARACTER (len=*), PARAMETER :: cmd = 'rmdir /Q /S '
#else
      CHARACTER (len=*), PARAMETER :: cmd = 'rm -rf '
#endif

!  Start of executable code

      CALL system(cmd // TRIM(directory_source), error)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Queries if the a path is absoulte.
!>
!>  Determines if the root of a path exists or if the path is relative. If the
!>  path is absoulte return true. If the path is relative return false. An empty
!>  path is treated as absoulte.
!>
!>  @param[in] path Path to test.
!>  @returns True if the path is absolute.
!-------------------------------------------------------------------------------
      FUNCTION is_absolute_path(path)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: is_absolute_path
      CHARACTER (len=*), INTENT(in) :: path

!  Start of executable code
!  Treat empty paths as absolute.
      IF (path .eq. '') THEN
         is_absolute_path = .true.
         RETURN
      END IF

#if defined(WIN32)
!  On windows systems, root paths start with a variey of cases. Taken from
!
!  http://msdn.microsoft.com/en-us/library/windows/desktop/aa365247(v=vs.85).aspx#paths
!
!    * A UNC name of any format, which always start with two backslash characters
!      ("\\"). For more information, see the next section.
!    * A disk designator with a backslash, for example "C:\" or "d:\".
!    * A single backslash, for example, "\directory" or "\file.txt". This is
!      also referred to as an absolute path.
!
      is_absolute_path = (path(1:2)         .eq. '\\') .or.                    &
     &                   (INDEX(path, ':\') .eq. 0)    .or.                    &
     &                   (path(1:1)         .eq. '\')
#else
!  On unix systems, an absolute path starts with either / or ~
      is_absolute_path = (path(1:1) .eq. '/') .or.                             &
     &                   (path(1:1) .eq. '~')
#endif

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Returns the directory of a file.
!>
!>  Searches a path for last path separator and returns the directory path.
!>
!>  @param[in] path Path to file.
!>  @returns The path of the file.
!-------------------------------------------------------------------------------
      FUNCTION get_path_of_file(file)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=path_length)   :: get_path_of_file
      CHARACTER (len=*), INTENT(in) :: file

!  Start of executable code
#if defined(WIN32)
      get_path_of_file =                                                       &
     &   file(1:(INDEX(TRIM(file), '\', back = .true.) - 1))
#else
      get_path_of_file =                                                       &
     &   file(1:(INDEX(TRIM(file), '/', back = .true.) - 1))
#endif

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Builds a path.
!>
!>  Builds a path by concatenating a string onto a path.
!>
!>  @param[in] path Existing path.
!>  @param[in] name String to add to path.
!>  @returns The new path.
!-------------------------------------------------------------------------------
      FUNCTION build_path(path, name)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=path_length)   :: build_path
      CHARACTER (len=*), INTENT(in) :: path
      CHARACTER (len=*), INTENT(in) :: name

!  Start of executable code
      IF (TRIM(path) .eq. '') THEN
         build_path = name
         RETURN
      END IF

#if defined(WIN32)
      build_path = TRIM(path) // '\' // TRIM(name)
#else
      build_path = TRIM(path) // '/' // TRIM(name)
#endif

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Change working directory
!>
!>  Changes the current working directory of the running process. CHDIR is an
!>  extension on many compilers. Check for compatability.
!>
!>  @param[in] path Path of the working directory to change to.
!-------------------------------------------------------------------------------
      SUBROUTINE change_directory(path, error)
#if defined(__INTEL_COMPILER)
      USE IFLPORT
#endif

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=*), INTENT(in) :: path
      INTEGER, INTENT(out)          :: error

!  local parameters
      error = CHDIR(TRIM(path))

      END SUBROUTINE

      END MODULE
