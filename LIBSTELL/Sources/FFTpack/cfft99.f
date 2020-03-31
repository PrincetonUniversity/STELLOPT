      SUBROUTINE cfft99(a, work, trigs, ifax, inc, jump, n, lot, isign)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER inc, jump, n, lot, isign
      INTEGER, DIMENSION(*) :: ifax
      REAL(rprec) :: a(*), work(*), trigs(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: nn, ink, jum, nfax, jnk, jst, ibase, ilast, nh, l, i1,
     1   i2, m, jbase, i, j, igo, la, k
      REAL(rprec) :: hREAL, himag
C-----------------------------------------------
c
c     SUBROUTINE 'cfft99' - multiple fast complex fourier transform
c
c     a is the array containing input and output data
c     work is an area of SIZE n*lot
c     trigs is a previously prepared list of trig FUNCTION values
c     ifax is a previously prepared list of factors of n
c     inc is the increment within each data "vector"
c         (e.g. inc=1 for consecutively stored data)
c     jump is the increment between the start of each data vector
c     n is the length of the data vectors
c     lot is the number of data vectors
c     isign = +1 for transform from spectral to gridpoint
c           = -1 for transform from gridpoint to spectral
c
c
c     vectorization is achieved on cray by DOing the transforms in
c     parallel.
c
c
      nn = n+n
      ink=inc+inc
      jum = jump+jump
      nfax=ifax(1)
      jnk = 2
      jst = 2
      IF (isign.ge.0) GOTO 30
      jnk = -2
      jst = nn-2
      IF (MOD(nfax,2).eq.1) GOTO 40
      ibase = 1
      ilast = (n-1)*ink
      nh = n/2
      DO 20 l=1,lot
      i1 = ibase+ink
      i2 = ibase+ilast
cdir$ ivdep
      DO 10 m=1,nh
c     swap REAL and imaginary portions
      hREAL = a(i1)
      himag = a(i1+1)
      a(i1) = a(i2)
      a(i1+1) = a(i2+1)
      a(i2) = hreal
      a(i2+1) = himag
      i1 = i1+ink
      i2 = i2-ink
   10 CONTINUE
      ibase = ibase+jum
   20 CONTINUE
      GOTO 100
c
   30 CONTINUE
      IF (MOD(nfax,2).eq.0) GOTO 100
c
   40 CONTINUE
c
c     during the transform process, nfax steps are taken, and the
c     results are stored alternately in work and in a.  IF nfax is
c     odd, the input data are first moved to work so that the final
c     result (after nfax steps) is stored in array a.
c
      ibase=1
      jbase=1
      DO 60 l=1,lot
c     move REAL and imaginary portions of element zero
      work(jbase) = a(ibase)
      work(jbase+1) = a(ibase+1)
      i=ibase+ink
      j=jbase+jst
cdir$ ivdep
      DO 50 m=2,n
c     move REAL and imaginary portions of other elements (possibly in
c     reverse order, depending on jst and jnk)
      work(j) = a(i)
      work(j+1) = a(i+1)
      i=i+ink
      j=j+jnk
   50 CONTINUE
      ibase=ibase+jum
      jbase=jbase+nn
   60 CONTINUE
c
  100 CONTINUE
c
c     perform the transform passes, one pass for each factor.  during
c     each pass the data are moved from a to work or from work to a.
c
c     for nfax even, the first pass moves from a to work
      igo = 110
c     for nfax odd, the first pass moves from work to a
      IF (MOD(nfax,2).eq.1) igo = 120
      la=1
      DO 140 k=1,nfax
      IF (igo.eq.120) GOTO 120
  110 CONTINUE
      CALL vpassm(a(1),a(2),work(1),work(2),trigs,
     *   ink,2,jum,nn,lot,n,ifax(k+1),la)
      igo=120
      GOTO 130
  120 CONTINUE
      CALL vpassm(work(1),work(2),a(1),a(2),trigs,
     *    2,ink,nn,jum,lot,n,ifax(k+1),la)
      igo=110
  130 CONTINUE
      la=la*ifax(k+1)
  140 CONTINUE
c
c     at this point the final transform result is stored in a.
c
      END SUBROUTINE cfft99
