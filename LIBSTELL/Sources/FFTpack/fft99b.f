      SUBROUTINE fft99b(work, a, trigs, inc, jump, n, lot)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER inc, jump, n, lot
      REAL(rprec), DIMENSION(*) :: work, a, trigs
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: one = 1, zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER::nh,nx,ink,ia,ib,ja,jb,l,iabase,ibbase,jabase,jbbase,k
      REAL(rprec) :: scale, c, s
C-----------------------------------------------
c
c     SUBROUTINE fft99b - postprocessing step for fft99, isign=-1
c     (gridpoint to spectral transform)
c
      nh=n/2
      nx=n+1
      ink=inc+inc
c
c     a(0) and a(n/2)
      scale = one/n
      ia=1
      ib=2
      ja=1
      jb=n*inc+1
cdir$ ivdep
      DO 10 l=1,lot
      a(ja)=scale*(work(ia)+work(ib))
      a(jb)=scale*(work(ia)-work(ib))
      a(ja+inc) = zero
      a(jb+inc) = zero
      ia=ia+nx
      ib=ib+nx
      ja=ja+jump
      jb=jb+jump
   10 CONTINUE
c
c     remaining wavenumbers
      scale=0.5_dp*scale
      iabase=3
      ibbase=n-1
      jabase=2*inc+1
      jbbase=(n-2)*inc+1
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
      a(ja)=scale*((work(ia)+work(ib))
     *   +(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
      a(jb)=scale*((work(ia)+work(ib))
     *   -(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
      a(ja+inc)=scale*((c*(work(ia)-work(ib))-s*(work(ia+1)+work(ib+1)))
     *    +(work(ib+1)-work(ia+1)))
      a(jb+inc)=scale*((c*(work(ia)-work(ib))-s*(work(ia+1)+work(ib+1)))
     *    -(work(ib+1)-work(ia+1)))
      ia=ia+nx
      ib=ib+nx
      ja=ja+jump
      jb=jb+jump
   20 CONTINUE
      iabase=iabase+2
      ibbase=ibbase-2
      jabase=jabase+ink
      jbbase=jbbase-ink
   30 CONTINUE
c
      IF (iabase.ne.ibbase) GOTO 50
c     wavenumber n/4 (IF it exists)
      ia=iabase
      ja=jabase
      scale=2*scale
cdir$ ivdep
      DO 40 l=1,lot
      a(ja)=scale*work(ia)
      a(ja+inc)=-scale*work(ia+1)
      ia=ia+nx
      ja=ja+jump
   40 CONTINUE
c
   50 CONTINUE

      END SUBROUTINE fft99b
