      SUBROUTINE fft99a(a, work, trigs, inc, jump, n, lot)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER inc, jump, n, lot
      REAL(rprec), DIMENSION(*) :: a, work, trigs
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER::nh,nx,ink,ia,ib,ja,jb,l,iabase,ibbase,jabase,jbbase,k
      REAL(rprec) :: c, s
C-----------------------------------------------
c
c     SUBROUTINE fft99a - preprocessing step for fft99, isign=+1
c     (spectral to gridpoint transform)
c
      nh=n/2
      nx=n+1
      ink=inc+inc
c
c     a(0) and a(n/2)
      ia=1
      ib=n*inc+1
      ja=1
      jb=2
cdir$ ivdep
      DO 10 l=1,lot
      work(ja)=a(ia)+a(ib)
      work(jb)=a(ia)-a(ib)
      ia=ia+jump
      ib=ib+jump
      ja=ja+nx
      jb=jb+nx
   10 CONTINUE
c
c     remaining wavenumbers
      iabase=2*inc+1
      ibbase=(n-2)*inc+1
      jabase=3
      jbbase=n-1
c
      DO 30 k=3,nh,2
      ia=iabase
      ib=ibbase
      ja=jabase
      jb=jbbase
      c=trigs(n+k)
      s=trigs(n+k+1)
cdir$ ivdep
      DO 20 l=1,lot
      work(ja)=(a(ia)+a(ib))-
     *    (s*(a(ia)-a(ib))+c*(a(ia+inc)+a(ib+inc)))
      work(jb)=(a(ia)+a(ib))+
     *    (s*(a(ia)-a(ib))+c*(a(ia+inc)+a(ib+inc)))
      work(ja+1)=(c*(a(ia)-a(ib))-s*(a(ia+inc)+a(ib+inc)))+
     *    (a(ia+inc)-a(ib+inc))
      work(jb+1)=(c*(a(ia)-a(ib))-s*(a(ia+inc)+a(ib+inc)))-
     *    (a(ia+inc)-a(ib+inc))
      ia=ia+jump
      ib=ib+jump
      ja=ja+nx
      jb=jb+nx
   20 CONTINUE
      iabase=iabase+ink
      ibbase=ibbase-ink
      jabase=jabase+2
      jbbase=jbbase-2
   30 CONTINUE
c
      IF (iabase.ne.ibbase) GOTO 50
c     wavenumber n/4 (IF it exists)
      ia=iabase
      ja=jabase
cdir$ ivdep
      DO 40 l=1,lot
      work(ja)=2.0*a(ia)
      work(ja+1)=-2.0*a(ia+inc)
      ia=ia+jump
      ja=ja+nx
   40 CONTINUE
c
   50 CONTINUE

      END SUBROUTINE fft99a
