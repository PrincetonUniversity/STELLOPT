      SUBROUTINE ga_shuffle(ipick)
c#######################################################################
c
c  This routine shuffles the parent array and its corresponding fitness
c
      USE ga_mod
      IMPLICIT NONE
      INTEGER :: ipick, j, n, itemp, iother
      REAL(rprec) :: rand, temp
      SAVE
c
      ipick=1
      DO 10 j=1,npopsiz-1
         CALL ran3(1,rand)
         iother=j+1+INT((npopsiz-j)*rand)
         DO 20 n=1,nchrome
            itemp=iparent(n,iother)
            iparent(n,iother)=iparent(n,j)
            iparent(n,j)=itemp
 20      CONTINUE
         temp=fitness(iother)
         fitness(iother)=fitness(j)
         fitness(j)=temp
 10   CONTINUE
c
      END SUBROUTINE ga_shuffle
