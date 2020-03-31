      SUBROUTINE rddisk(iunit, a, now, incnow, irec, ierr)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER iunit, now, incnow, irec, ierr
      REAL(rprec), DIMENSION(now) :: a
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ig, irech, il
C-----------------------------------------------
c
c     Author:    H. Maassberg       Sept. 1988
c
c     auxilary routine for routines BLKTRD / BLK5D (DKES code)
c     for disk i/o on CRAY at Garching
c
c-----------------------------------------------------------------------
c
c     READ NOW words of vector A from disk (fortran IUNIT) with
c     direct access in record IREC
c
c-----------------------------------------------------------------------
      ig = 0
      irech = irec - 1

      DO WHILE (ig < now)
         il = ig + 1
         irech = irech + 1
         ig = MIN(now,ig + incnow)
         READ (iunit, rec=irech, err=10) a(il:ig)
      END DO
      ierr = 0
      RETURN
   10 CONTINUE
      ierr = 1
      PRINT *, ' error detected in disk i/o  (routine RDDISK)'

      END SUBROUTINE rddisk
