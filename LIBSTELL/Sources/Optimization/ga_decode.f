      SUBROUTINE ga_decode(i,array,iarray)
c#######################################################################
c
c  This routine decodes a binary string to a REAL number.
c
      USE ga_mod
      IMPLICIT NONE
      INTEGER :: i, j, m, l, k, iparam
      REAL(rprec), DIMENSION(nparmax,indmax) :: array
      INTEGER, DIMENSION(nchrmax,indmax) :: iarray
      SAVE

      l=1
      DO 10 k=1,nparam
         iparam=0
         m=l
         DO 20 j=m,m+ig2(k)-1
            l=l+1
            iparam=iparam+iarray(j,i)*(2**(m+ig2(k)-1-j))
 20      CONTINUE
         array(k,i)=g0(k)+g1(k)*iparam
 10   CONTINUE

      END SUBROUTINE ga_decode
