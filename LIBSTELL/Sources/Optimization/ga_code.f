      SUBROUTINE ga_code(j,k,array,iarray)
c#######################################################################
c
c
c  This routine codes a PARAMETER into a binary string.
c
      USE ga_mod
      IMPLICIT NONE
      INTEGER :: j, k, i, istart, iparam, m
      REAL(rprec), DIMENSION(nparmax,indmax) :: array
      INTEGER, DIMENSION(nchrmax,indmax) :: iarray

      SAVE
c
c  First, establish the beginning location of the PARAMETER string of
c  interest.
      istart=1
      DO 10 i=1,k-1
         istart=istart+ig2(i)
 10   CONTINUE
c
c  Find the equivalent coded PARAMETER value, and back out the binary
c  string by factors of two.
      m=ig2(k)-1
      IF (g1(k).eq.0.0_dp) RETURN
      iparam=NINT((array(k,j)-g0(k))/g1(k))
      DO 20 i=istart,istart+ig2(k)-1
         iarray(i,j)=0
         IF ((iparam+1).gt.(2**m)) THEN
            iarray(i,j)=1
            iparam=iparam-2**m
         END IF
         m=m-1
 20   CONTINUE
c     WRITE(3,*)array(k,j),iparam,(iarray(i,j),i=istart,istart+ig2(k)-1)
c
      END SUBROUTINE ga_code
