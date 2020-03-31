      SUBROUTINE ga_possibl(array, iarray, myid)
c#######################################################################
c
c  This SUBROUTINE determines whether or not ALL parameters are within
c  the specified range of possibility.  If not, the PARAMETER is
c  randomly reassigned within the range.  This SUBROUTINE is ONLY
c  necessary when the number of possibilities per PARAMETER is not
c  optimized to be 2**n, i.e. IF nposSUM < ig2sum.
c
      USE ga_mod
      USE mpi_params, ONLY: master
      IMPLICIT NONE

      REAL(rprec), DIMENSION(nparmax,indmax) :: array
      INTEGER, DIMENSION(nchrmax,indmax) :: iarray
      INTEGER :: i, j, n2ig2j, irand, myid
      REAL(rprec) :: rand

      SAVE
c
      DO 10 i=1,npopsiz
         CALL ga_decode(i,array,iarray)
         DO 20 j=1,nparam
            n2ig2j=ig2(j)
            IF(nposibl(j).ne.n2ig2j .and. array(j,i).gt.par_max(j)) THEN
               CALL ran3(1,rand)
               irand=INT((2**nposibl(j))*rand)
               array(j,i)=g0(j)+irand*g1(j)
               CALL ga_code(i,j,array,iarray)
               IF (nowrite.eq.0 .and. myid.eq.master) THEN
                  WRITE(6,1000) i,j
                  WRITE(iunit_ga_out,1000) i,j
               END IF
            END IF
 20      CONTINUE
 10   CONTINUE

 1000 FORMAT('*** Parameter adjustment to individual     ',i4,
     1       ', PARAMETER  ',i3,' ***')

      END SUBROUTINE ga_possibl
