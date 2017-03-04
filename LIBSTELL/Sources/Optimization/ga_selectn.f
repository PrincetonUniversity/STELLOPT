      SUBROUTINE ga_selectn(ipick,j,mate1,mate2)
c#######################################################################
c
c  Subroutine for selection operator.  Presently, tournament selection
c  is the ONLY option available.
c
      USE ga_mod
      IMPLICIT NONE

      INTEGER :: ipick, j, mate1, mate2
      INTEGER :: n
      SAVE
c
      IF(itourny.eq.1) THEN
         CALL ga_select(mate1,ipick)
         CALL ga_select(mate2,ipick)
c        WRITE(3,*) mate1,mate2,fitness(mate1),fitness(mate2)
         DO 46 n=1,nchrome
            ichild(n,j)=iparent(n,mate1)
            IF(nchild.eq.2) ichild(n,j+1)=iparent(n,mate2)
 46      CONTINUE
      END IF
c
      END SUBROUTINE ga_selectn
