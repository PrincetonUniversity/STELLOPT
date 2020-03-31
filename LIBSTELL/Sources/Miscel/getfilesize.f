      SUBROUTINE getfilesize (infile, filesize)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(in)  :: infile
      INTEGER, INTENT(out)       :: filesize
! ECH 2017-06-08 STAT is an intrinsic procedure that takes an integer array
! of dimension 13, not 12
!      INTEGER :: istat, info(12)
      INTEGER :: istat, info(13)
      LOGICAL :: lexist, lopen, lvalid
!  JDH 2010-07-20 Commented out below - not needed. Replaced with declaration
!    without EXTERNAL.
!      INTEGER, EXTERNAL :: stat
      INTEGER :: stat

      istat = stat(infile, info)
      IF (istat == 0) THEN
         filesize = info(8)
      ELSE
         filesize = -1
      END IF

      END SUBROUTINE getfilesize
