      SUBROUTINE ga_niche(myid)
c#######################################################################
c
c  Implement "niching" through Goldberg''s multidimensional phenotypic
c  sharing scheme with a triangular sharing FUNCTION.  To find the
c  multidimensional distance from the best individual, normalize ALL
c  PARAMETER differences.
c
      USE ga_mod
      USE mpi_params, ONLY: master
      IMPLICIT NONE
      REAL(rprec) :: alpha, del, del2, sigshar, SUMshar, share
      INTEGER :: nniche, jj, ii, j, k, myid
      SAVE
c
c   Variable definitions:
c
c  alpha   = power law exponent for sharing function; typically = 1.0
c  del     = normalized multidimensional distance between ii and ALL
c            other members of the population
c            (equals the square root of del2)
c  del2    = sum of the squares of the normalized multidimensional
c            distance between member ii and all other members of
c            the population
c  nniche  = number of niched parameters
c  sigshar = normalized distance to be compared with del; in some sense,
c            1/sigshar can be viewed as the number of regions over which
c            the sharing function should focus, e.g. with sigshar=0.1,
c            the sharing function will try to clump in ten distinct
c            regions of the phase space.  a value of sigshar on the
c            order of 0.1 seems to work best.
c  share   = sharing function between individual ii and j
c  sumshar = sum of the sharing functions for individual ii
c
      alpha=1
      sigshar=0.1_dp
      nniche=0
      DO 33 jj=1,nparam
         nniche=nniche+nichflg(jj)
 33   CONTINUE
      IF (nniche.eq.0) THEN
         IF (myid .eq. master) THEN
            WRITE(6,1900)
            WRITE(iunit_ga_out,1900)
            CLOSE(iunit_ga_out)
         END IF
         STOP
      END IF
      DO 34 ii=1,npopsiz
         SUMshar=0
         DO 35 j=1,npopsiz
            del2=0
            DO 36 k=1,nparam
               IF (nichflg(k).ne.0) THEN
                  del2=del2+((parent(k,j)-parent(k,ii))/pardel(k))**2
               END IF
 36         CONTINUE
            del=SQRT(del2)/nniche
            IF (del.lt.sigshar) THEN
c               share=1.0-((del/sigshar)**alpha)
               share=1-(del/sigshar)
            ELSE
               share=0
            END IF
            SUMshar=sumshar+share/npopsiz
 35      CONTINUE
         IF (sumshar.ne.0.0_dp) fitness(ii)=fitness(ii)/sumshar
 34   CONTINUE

 1900 FORMAT(1x,'ERROR: iniche=1 and ALL values in nichflg array = 0'/
     1       1x,'       Do you want to niche or not?')

      END SUBROUTINE ga_niche
