      SUBROUTINE ga_restart(i,istart,kount,filename, myid)
c#######################################################################
c
c  this subroutine writes restart information to the ga.restart file.
c
      USE ga_mod
      USE mpi_params, ONLY: master
      USE safe_open_mod
      IMPLICIT NONE
      INTEGER :: kount, i, j, l, istart, istat, myid
      CHARACTER(LEN=100) :: filename
      SAVE

      IF (myid .ne. master) RETURN

      kount=kount+1
      IF(i.eq.maxgen+istart-1 .or. kount.eq.kountmx) THEN
         iunit_ga_restart = 25
         CALL safe_open(iunit_ga_restart, istat,
     1                'ga_restart.'//TRIM(filename),
     2               'unknown', 'formatted')
         REWIND iunit_ga_restart
         WRITE(iunit_ga_restart,*) i+1,npopsiz
         DO 80 j=1,npopsiz
            WRITE(iunit_ga_restart,*) j,(iparent(l,j),l=1,nchrome)
c        IF(nchrome .le. 60) THEN
c           WRITE(iunit_ga_restart,1500) j,(iparent(l,j),l=1,nchrome)
c        ELSE
c           WRITE(iunit_ga_restart,1500) j,(iparent(l,j),l=1,60)
c           WRITE(iunit_ga_restart,1501) (iparent(l,j),l=61,nchrome)
c        END IF
 80      CONTINUE
         CLOSE (iunit_ga_restart)
         kount=0
      ENDif
c
 1500 FORMAT(i5,3x,60i2)
 1501 FORMAT(5x,3x,60i2)

      END SUBROUTINE ga_restart
