      SUBROUTINE ga_mutate (myid)
c#######################################################################
c
      USE ga_mod
      USE mpi_params, ONLY: master
      IMPLICIT NONE
      INTEGER :: nmutate, ncreep, j, k, myid
      REAL(rprec) :: rand, creep
      SAVE
c
c  This SUBROUTINE performs mutations on the children generation.
c  Perform random jump mutation IF a random number is less than pmutate.
c  Perform random creep mutation IF a different random number is less
c  than pcreep.
      nmutate=0
      ncreep=0
      DO 70 j=1,npopsiz
         DO 75 k=1,nchrome
c  Jump mutation
            CALL ran3(1,rand)
            IF (rand.le.pmutate) THEN
               nmutate=nmutate+1
               IF(ichild(k,j).eq.0) THEN
                  ichild(k,j)=1
               ELSE
                  ichild(k,j)=0
               END IF
               IF (nowrite.eq.0 .and. myid.eq.master) THEN
                  WRITE(6,1300) j,k
                  WRITE(iunit_ga_out,1300) j,k
               END IF
            END IF
 75      CONTINUE
c  Creep mutation (one discrete position away).
         IF (icreep.ne.0) THEN
            DO 76 k=1,nparam
               CALL ran3(1,rand)
               IF(rand.le.pcreep) THEN
                  CALL ga_decode(j,child,ichild)
                  ncreep=ncreep+1
                  creep=1
                  CALL ran3(1,rand)
                  IF (rand.lt.0.5_dp) creep=-1
                  child(k,j)=child(k,j)+g1(k)*creep
                  IF (child(k,j).gt.par_max(k)) THEN
                     child(k,j)=par_max(k)-1.0d0*g1(k)
                  ELSEIF (child(k,j).lt.par_min(k)) THEN
                     child(k,j)=par_min(k)+g1(k)
                  END IF
                  CALL ga_code(j,k,child,ichild)
                  IF (nowrite.eq.0 .and. myid.eq.master) THEN
                     WRITE(6,1350) j,k
                     WRITE(iunit_ga_out,1350) j,k
                  END IF
               END IF
 76         CONTINUE
         END IF
 70   CONTINUE
      IF (myid .eq. master) THEN
         WRITE(6,1250) nmutate,ncreep
         WRITE(iunit_ga_out,1250) nmutate,ncreep
      END IF

 1250 FORMAT(/'  Number of Jump Mutations  =',i5/
     +        '  Number of Creep Mutations =',i5)
 1300 FORMAT('*** Jump mutation performed on individual  ',i4,
     +       ', chromosome ',i3,' ***')
 1350 FORMAT('*** Creep mutation performed on individual ',i4,
     +       ', PARAMETER  ',i3,' ***')
c
      END SUBROUTINE ga_mutate
