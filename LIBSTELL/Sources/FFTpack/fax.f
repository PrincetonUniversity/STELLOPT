      SUBROUTINE fax(ifax, n, mode)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER n, mode
      INTEGER, DIMENSION(13) :: ifax
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: nn, k, l, inc, nfax, ii, istop, i, item
C-----------------------------------------------
      nn=n
      IF (IABS(mode).eq.1) GOTO 10
      IF (IABS(mode).eq.8) GOTO 10
      nn=n/2
      IF ((nn+nn).eq.n) GOTO 10
      ifax(1)=-99
      RETURN
   10 k=1
c     test for factors of 4
   20 IF (MOD(nn,4).ne.0) GOTO 30
      k=k+1
      ifax(k)=4
      nn=nn/4
      IF (nn.eq.1) GOTO 80
      GOTO 20
c     test for extra factor of 2
   30 IF (MOD(nn,2).ne.0) GOTO 40
      k=k+1
      ifax(k)=2
      nn=nn/2
      IF (nn.eq.1) GOTO 80
c     test for factors of 3
   40 IF (MOD(nn,3).ne.0) GOTO 50
      k=k+1
      ifax(k)=3
      nn=nn/3
      IF (nn.eq.1) GOTO 80
      GOTO 40
c     now find remaining factors
   50 l=5
      inc=2
c     inc alternately takes on values 2 and 4
   60 IF (MOD(nn,l).ne.0) GOTO 70
      k=k+1
      ifax(k)=l
      nn=nn/l
      IF (nn.eq.1) GOTO 80
      GOTO 60
   70 l=l+inc
      inc=6-inc
      GOTO 60
   80 ifax(1)=k-1
c     ifax(1) CONTAINS number of factors
      nfax=ifax(1)
c     sort factors into ascending order
      IF (nfax.eq.1) GOTO 110
      DO 100 ii=2,nfax
      istop=nfax+2-ii
      DO 90 i=2,istop
      IF (ifax(i+1).ge.ifax(i)) GOTO 90
      item=ifax(i)
      ifax(i)=ifax(i+1)
      ifax(i+1)=item
   90 CONTINUE
  100 CONTINUE
  110 CONTINUE

      END SUBROUTINE fax
