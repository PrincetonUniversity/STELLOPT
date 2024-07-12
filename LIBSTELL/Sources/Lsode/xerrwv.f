      SUBROUTINE xerrwv(msg, nmes, nerr, level, ni, i1, i2, nr, r1, r2)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: nmes, nerr, level, ni, i1, i2, nr
      REAL(rprec) r1, r2
      CHARACTER(LEN=*) :: msg
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
!      INTEGER :: lun, ncpw, nch, nwds
      INTEGER, PARAMETER ::  lunit = 6
c-----------------------------------------------------------------------
c Subroutines xerrwv, xsetf, and xsetun, as given here, constitute
c a simplified version of the slatec error handling package.
c written by a. c. hindmarsh at llnl.  version of march 30, 1987.
c the following is valid for the pdp-11, or vax with 2-byte integers.
c-----------------------------------------------------------------------
c WRITE the message. ---------------------------------------------------
      WRITE (lunit, '(a)') TRIM(msg)
      IF (ni == 1) WRITE (lunit, 20) i1
   20 FORMAT(6x,'in above message,  i1 =',i10)
      IF (ni == 2) WRITE (lunit, 30) i1, i2
   30 FORMAT(6x,'in above message,  i1 =',i10,3x,'i2 =',i10)
      IF (nr == 1) WRITE (lunit, 40) r1
   40 FORMAT(6x,'in above message,  r1 =',d14.6)
      IF (nr == 2) WRITE (lunit, 50) r1, r2
   50 FORMAT(6x,'in above,  r1 =',d21.13,3x,'r2 =',d14.6)

c abort the run if level = 2. ------------------------------------------
      IF (level /= 2) RETURN
      STOP

      END SUBROUTINE xerrwv
