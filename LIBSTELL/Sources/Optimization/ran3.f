      SUBROUTINE ran3(idum,rand)
c#######################################################################
c
c  Returns a uniform random deviate between 0.0 and 1.0.  Set idum to
c  ANY negative value to initialize or reinitialize the sequence.
c  This FUNCTION is taken from W.H. Press, "Numerical Recipes" p. 199.
c
      USE stel_kinds
      IMPLICIT NONE
      INTEGER :: idum, IFf, inext, inextp, i, ii, k
      REAL(rprec) :: rand, ma, mj, mk
      SAVE
      REAL(rprec), PARAMETER :: mbig=4000000, mseed=1618033,
     1           mz=0, fac=1._dp/mbig
c     PARAMETER (mbig=1000000000,mseed=161803398,mz=0,fac=1./mbig)
c
c  According to Knuth, ANY large mbig, and ANY smaller (but still large)
c  mseed can be substituted for the above values.
      DIMENSION ma(55)
      data IFf /0/
      IF (idum.lt.0 .or. IFf.eq.0) THEN
         IFf=1
         mj=mseed-ABS(idum)
         mj=MOD(mj,mbig)
         ma(55)=mj
         mk=1
         DO 11 i=1,54
            ii=MOD(21*i,55)
            ma(ii)=mk
            mk=mj-mk
            IF(mk.lt.mz) mk=mk+mbig
            mj=ma(ii)
 11      CONTINUE
         DO 13 k=1,4
            DO 12 i=1,55
               ma(i)=ma(i)-ma(1+MOD(i+30,55))
               IF(ma(i).lt.mz) ma(i)=ma(i)+mbig
 12         CONTINUE
 13      CONTINUE
         inext=0
         inextp=31
         idum=1
      ENDif
      inext=inext+1
      IF(inext.eq.56) inext=1
      inextp=inextp+1
      IF(inextp.eq.56) inextp=1
      mj=ma(inext)-ma(inextp)
      IF(mj.lt.mz) mj=mj+mbig
      ma(inext)=mj
      rand=mj*fac

      END SUBROUTINE ran3
