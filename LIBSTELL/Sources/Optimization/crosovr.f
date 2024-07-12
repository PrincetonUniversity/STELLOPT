      SUBROUTINE crosovr(ncross,j,mate1,mate2)
c#######################################################################
c
c  Subroutine for crossover between the randomly selected pair.
      USE stel_kinds
      USE ga_mod
      IMPLICIT NONE
      INTEGER :: ncross, j, mate1, mate2
      INTEGER :: n, icross
      REAL(rprec) :: rand
      SAVE
c
      IF (iunifrm.eq.0) THEN
c  Single-point crossover at a random chromosome point.
         CALL ran3(1,rand)
         IF(rand.gt.pcross) GOTO 69
         ncross=ncross+1
         CALL ran3(1,rand)
         icross=2+INT((nchrome-1)*rand)
         DO 50 n=icross,nchrome
            ichild(n,j)=iparent(n,mate2)
            IF(nchild.eq.2) ichild(n,j+1)=iparent(n,mate1)
 50      CONTINUE
      ELSE
c  Perform uniform crossover between the randomly selected pair.
         DO 60 n=1,nchrome
            CALL ran3(1,rand)
            IF(rand.le.pcross) THEN
               ncross=ncross+1
               ichild(n,j)=iparent(n,mate2)
               IF(nchild.eq.2) ichild(n,j+1)=iparent(n,mate1)
            END IF
 60      CONTINUE
      END IF
 69   CONTINUE

      END SUBROUTINE crosovr
