      SUBROUTINE ga_newgen(npossum,ig2sum,myid)
c#######################################################################
c
c  Write child array back into parent array for new generation.  Check
c  to see IF the best parent was replicated; IF not, and IF ielite=1,
c  THEN reproduce the best parent into a random slot.
c
      USE ga_mod
      USE mpi_params, ONLY: master
      IMPLICIT NONE
      INTEGER :: npossum, ig2sum, kelite, jelite
      INTEGER :: irand, j, n, myid
      REAL(rprec) :: rand
      SAVE
c
      IF (npossum.lt.ig2sum) CALL ga_possibl(child,ichild,myid)
      kelite=0
      DO 94 j=1,npopsiz
         jelite=0
         DO 95 n=1,nchrome
            iparent(n,j)=ichild(n,j)
            IF (iparent(n,j).eq.ibest(n)) jelite=jelite+1
            IF (jelite.eq.nchrome) kelite=1
 95      CONTINUE
 94   CONTINUE
      IF (ielite.ne.0 .and. kelite.eq.0) THEN
         CALL ran3(1,rand)
         irand=1+INT(npopsiz*rand)
         iparent(1:nchrome,irand)=ibest(1:nchrome)
         IF (myid .eq. master) WRITE(iunit_ga_out,1260) irand
      END IF
c
 1260 FORMAT('  Elitist Reproduction on Individual ',i4)
c
      END SUBROUTINE ga_newgen
