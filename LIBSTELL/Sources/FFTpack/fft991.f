      SUBROUTINE fft991(a, work, trigs, ifax, inc, jump, n, lot, isign)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER inc, jump, n, lot, isign
      INTEGER, DIMENSION(13) :: ifax
      REAL(rprec), DIMENSION(*) :: a, work, trigs
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER::nfax,nx,nh,ink,igo,ibase,jbase,l,i,j,m,ia,la,k,ib
C-----------------------------------------------
c
c     SUBROUTINE 'fft991' - multiple REAL/half-complex periodic
c     fast fourier transform
c
c     same as fft99 except that ordering of data corresponds to
c     that in mrfft2
c
c     PROCEDURE used to convert to half-length complex transform
c     is given by cooley, lewis and welch (j. sound vib., vol. 12
c     (1970), 315-337)
c
c     a is the array containing input and output data
c     work is an area of SIZE (n+1)*lot
c     trigs is a previously prepared list of trig FUNCTION values
c     ifax is a previously prepared list of factors of n/2
c     inc is the increment within each data "vector"
c         (e.g. inc=1 for consecutively stored data)
c     jump is the increment between the start of each data vector
c     n is the length of the data vectors
c     lot is the number of data vectors
c     isign = +1 for transform from spectral to gridpoint
c           = -1 for transform from gridpoint to spectral
c
c     ordering of coefficients:
c         a(0),b(0),a(1),b(1),a(2),b(2),...,a(n/2),b(n/2)
c         WHERE b(0)=b(n/2)=0; (n+2) locations required
c
c     ordering of data:
c         x(0),x(1),x(2),...,x(n-1)
c
c     vectorization is achieved on cray by DOing the transforms in
c     parallel
c
c     *** n.b. n is assumed to be an even number
c
c     definition of transforms:
c     -------------------------
c
c     isign=+1: x(j)=SUM(k=0,...,n-1)(c(k)*EXP(2*i*j*k*pi/n))
c         WHERE c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)
c
c     isign=-1: a(k)=(1/n)*SUM(j=0,...,n-1)(x(j)*COS(2*j*k*pi/n))
c               b(k)=-(1/n)*SUM(j=0,...,n-1)(x(j)*SIN(2*j*k*pi/n))
c
      nfax=ifax(1)
      nx=n+1
      nh=n/2
      ink=inc+inc
      IF (isign.eq.+1) GOTO 30
c
c     IF necessary, transfer data to work area
      igo=50
      IF (MOD(nfax,2).eq.1) GOTO 40
      ibase=1
      jbase=1
      DO 20 l=1,lot
      i=ibase
      j=jbase
cdir$ ivdep
      DO 10 m=1,n
      work(j)=a(i)
      i=i+inc
      j=j+1
   10 CONTINUE
      ibase=ibase+jump
      jbase=jbase+nx
   20 CONTINUE
c
      igo=60
      GOTO 40
c
c     preprocessing (isign=+1)
c     ------------------------
c
   30 CONTINUE
      CALL fft99a(a,work,trigs,inc,jump,n,lot)
      igo=60
c
c     complex transform
c     -----------------
c
   40 CONTINUE
      ia=1
      la=1
      DO 80 k=1,nfax
      IF (igo.eq.60) GOTO 60
   50 CONTINUE
      CALL vpassm(a(ia),a(ia+inc),work(1),work(2),trigs,
     *   ink,2,jump,nx,lot,nh,ifax(k+1),la)
      igo=60
      GOTO 70
   60 CONTINUE
      CALL vpassm(work(1),work(2),a(ia),a(ia+inc),trigs,
     *    2,ink,nx,jump,lot,nh,ifax(k+1),la)
      igo=50
   70 CONTINUE
      la=la*ifax(k+1)
   80 CONTINUE
c
      IF (isign.eq.-1) GOTO 130
c
c     IF necessary, transfer data from work area
      IF (MOD(nfax,2).eq.1) GOTO 110
      ibase=1
      jbase=1
      DO 100 l=1,lot
      i=ibase
      j=jbase
cdir$ ivdep
      DO 90 m=1,n
      a(j)=work(i)
      i=i+1
      j=j+inc
   90 CONTINUE
      ibase=ibase+nx
      jbase=jbase+jump
  100 CONTINUE
c
c     fill in zeros at END
  110 CONTINUE
      ib=n*inc+1
cdir$ ivdep
      DO 120 l=1,lot
      a(ib)=0.0
      a(ib+inc)=0.0
      ib=ib+jump
  120 CONTINUE
      GOTO 140
c
c     postprocessing (isign=-1):
c     --------------------------
c
  130 CONTINUE
      CALL fft99b(work,a,trigs,inc,jump,n,lot)
c
  140 CONTINUE

      END SUBROUTINE fft991
