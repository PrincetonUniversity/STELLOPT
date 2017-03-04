      MODULE read_v3post_mod
!
!     Use READ_V3POST_MOD to include variables dynamically allocated
!     in this module
!     Call DEALLOCATE_READ_V3POST to free this memory when it is no longer needed
!
      USE stel_kinds

      IMPLICIT NONE
!DEC$ IF DEFINED (NETCDF)
C-----------------------------------------------
C   L O C A L   P A R A M E T E R S
C-----------------------------------------------
! Variable names (vn_...) : put eventually into library, used by read_wout too...
      CHARACTER (LEN=*), PARAMETER :: 
     1  vn_cal = 'signal_diag_cal',
     2  vn_cext = 'signal_diag_cext',
     3  vn_plasma = 'signal_diag_plasma',
     4  vn_sname = 'signal_diag_sname'
      CHARACTER(LEN=*), PARAMETER, DIMENSION(1) ::
     1             d1dim = (/'num_diagno'/)
      CHARACTER(LEN=*), DIMENSION(2), PARAMETER :: 
     1             d2dim = (/'str_len   ','num_diagno'/)
!DEC$ ENDIF
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  :: num_diagno
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: signal_diag_cext,
     1                                          signal_diag_plasma
      CHARACTER (LEN=30), ALLOCATABLE, DIMENSION(:) :: signal_sname
!-----------------------------------------------

      CONTAINS

      SUBROUTINE read_v3post_file (file_or_extension, filename, ierr)
      USE safe_open_mod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(out) :: ierr
      CHARACTER(LEN=*), INTENT(in)  :: file_or_extension
      CHARACTER(LEN=*), INTENT(out) :: filename
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, PARAMETER :: iunit_init = 10
      INTEGER :: iunit
      LOGICAL :: isnc
C-----------------------------------------------
!
!     THIS SUBROUTINE READS THE WOUT FILE CREATED BY THE VMEC CODE
!     AND STORES THE DATA IN THE READ_WOUT MODULE
!
!     FIRST, CHECK IF THIS IS A FULLY-QUALIFIED PATH NAME
!     MAKE SURE wout IS NOT EMBEDDED IN THE NAME (PERVERSE USER...)
!
      filename = 'v3post'
      CALL parse_extension(filename, file_or_extension, isnc)
      IF (isnc) THEN
!DEC$ IF DEFINED (NETCDF)
         CALL read_v3post_nc (filename, ierr)
!DEC$ ELSE
         PRINT *, "NETCDF wout file can not be opened on this platform"
         ierr = -100
!DEC$ ENDIF
      ELSE
         iunit = iunit_init
         CALL safe_open (iunit, ierr, filename, 'old', 'formatted')
         IF (ierr .eq. 0) CALL read_v3post_text(iunit, ierr)
         CLOSE(unit=iunit)
      END IF

      END SUBROUTINE read_v3post_file


      SUBROUTINE read_v3post_text(iunit, ierr)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: iunit, ierr
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      STOP 'read_v3post TEXT file read not implemented!'
 
      END SUBROUTINE read_v3post_text


!DEC$ IF DEFINED (NETCDF)
      SUBROUTINE read_v3post_nc(filename, ierr)
      USE ezcdf
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(out) :: ierr
      CHARACTER(LEN=*), INTENT(in) :: filename
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ncdf, istat
      INTEGER, DIMENSION(2)   :: dimlens
C-----------------------------------------------
! Open cdf File
      call cdf_open(ncdf,filename,'r', ierr)
      IF (ierr .ne. 0) THEN
         PRINT *,' Error opening v3post .nc file'
         RETURN
      END IF

! Be sure all arrays are deallocated
      CALL read_v3post_deallocate

      CALL cdf_inquire(ncdf, vn_plasma, dimlens)
      num_diagno = dimlens(1)

      ALLOCATE (signal_diag_cext(num_diagno), 
     1          signal_diag_plasma(num_diagno), 
     1	signal_sname(num_diagno), stat=istat)

      IF (istat .ne. 0) STOP 'Allocation error in read_v3post'
      signal_diag_cext = 0 ; signal_diag_plasma = 0
      CALL cdf_read(ncdf, vn_cext, signal_diag_cext)
      CALL cdf_read(ncdf, vn_plasma, signal_diag_plasma)
      CALL cdf_read(ncdf, vn_sname, signal_sname)
      CALL cdf_close(ncdf)

      END SUBROUTINE read_v3post_nc
!DEC$ ENDIF

      SUBROUTINE read_v3post_deallocate
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: istat
!-----------------------------------------------
      istat = 0

      IF (ALLOCATED(signal_diag_cext)) DEALLOCATE (signal_diag_cext,
     1 signal_sname, signal_diag_plasma, stat=istat)

      END SUBROUTINE read_v3post_deallocate

      END MODULE read_v3post_mod
