      SUBROUTINE ga_micro(i,npossum,ig2sum, myid)
c#######################################################################
c
c  Micro-GA implementation SUBROUTINE
c
      USE ga_mod
      USE mpi_params, ONLY: master
      IMPLICIT NONE
      INTEGER :: i, npossum, ig2sum, myid
      INTEGER :: icount, j, n
      REAL(rprec) :: diffrac, rand
      SAVE
c
c
c  First, check for convergence of micro population.
c  If converged, start a new generation with best individual and fill
c  the remainder of the population with new randomly generated parents.
c
c  Count number of different bits from best member in micro-population
      icount=0
      DO 81 j=1,npopsiz
         DO 82 n=1,nchrome
            IF(iparent(n,j).ne.ibest(n)) icount=icount+1
 82      CONTINUE
 81   CONTINUE
c
c  If icount less than 5% of number of bits, THEN consider population
c  to be converged.  Restart with best individual and random others.
      diffrac=REAL(icount,rprec)/((npopsiz-1)*nchrome)
      IF (diffrac.lt.0.05_dp) THEN
         DO 87 n=1,nchrome
            iparent(n,1)=ibest(n)
 87      CONTINUE
         DO 88 j=2,npopsiz
            DO 89 n=1,nchrome
               CALL ran3(1,rand)
               iparent(n,j)=1
               IF(rand.lt.0.5_dp) iparent(n,j)=0
 89         CONTINUE
 88      CONTINUE
         IF (npossum.lt.ig2sum) CALL ga_possibl(parent,iparent,myid)
         IF (myid .eq. master) THEN
           WRITE(6,1375) i
            WRITE(iunit_ga_out,1375) i
         END IF
      END IF

 1375 FORMAT(//'%%%%%%%  Restart micro-population at generation',
     +       i5,'  %%%%%%%')
c
      END SUBROUTINE ga_micro
