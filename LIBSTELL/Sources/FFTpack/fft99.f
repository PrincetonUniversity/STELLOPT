      SUBROUTINE fft99(a, work, trigs, ifax, inc, jump, n, lot, isign)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER inc, jump, n, lot, isign
      INTEGER :: ifax(13)
      REAL(rprec) :: a(lot*(n+2)), work(lot*(n+1)), trigs(3*n/2+1)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER::nfax,nx,nh,ink,igo,ibase,jbase,l,i,j,m,ia,la,k,ib
C-----------------------------------------------
c
c purpose      performs multiple fast fourier transforms.  this package
c              will perform a number of simultaneous REAL/half-complex
c              periodic fourier transforms or corresponding inverse
c              transforms, i.e.  given a set of REAL data vectors, the
c              package returns a set of "half-complex" fourier
c              coefficient vectors, or vice versa.  the length of the
c              transforms must be an even number greater than 4 that has
c              no other factors except possibly powers of 2, 3, and 5.
c              this is an ALL fortran version of the craylib package
c              that is mostly written in cal.
c
c              the package fft99f CONTAINS several user-level routines:
c
c            SUBROUTINE fftfax
c                an initialization routine that must be called once
c                before a sequence of calls to the fft routines
c                (provided that n is not changed).
c
c            SUBROUTINEs fft99 and fft991
c                two fft routines that RETURN slightly different
c                arrangements of the data in gridpoint space.
c
c
c access       this fortran version may be accessed with
c
c                   *fortran,p=xlib,sn=fft99f
c
c              to access the cray object code, CALLing the user ENTRY
c              points from a cray PROGRAM is sufficient.  the source
c              fortran and cal code for the craylib version may be
c              accessed using
c
c                   fetch p=craylib,sn=fft99
c
c usage        let n be of the form 2**p * 3**q * 5**r, WHERE p .ge. 1,
c              q .ge. 0, and r .ge. 0.  THEN a typical sequence of
c              calls to transform a given set of REAL vectors of length
c              n to a set of "half-complex" fourier coefficient vectors
c              of length n is
c
c                   DIMENSION ifax(13),trigs(3*n/2+1),a(m*(n+2)),
c                  +          work(m*(n+1))
c
c                   CALL fftfax (n, ifax, trigs)
c                   CALL fft99 (a,work,trigs,ifax,inc,jump,n,m,isign)
c
c              see the individual WRITE-ups for fftfax, fft99, and
c              fft991 below, for a detailed description of the
c              arguments.
c
c history      the package was written by clive temperton at ecmwf in
c              november, 1978.  it was modified, documented, and tested
c              for ncar by russ rew in september, 1980.
c
c-----------------------------------------------------------------------
c
c SUBROUTINE fftfax (n,ifax,trigs)
c
c purpose      a set-up routine for fft99 and fft991.  it need ONLY be
c              called once before a sequence of calls to the fft
c              routines (provided that n is not changed).
c
c argument     ifax(13),trigs(3*n/2+1)
c dimensions
c
c arguments
c
c on input     n
c               an even number greater than 4 that has no prime factor
c               greater than 5.  n is the length of the transforms (see
c               the documentation for fft99 and fft991 for the
c               definitions of the transforms).
c
c              ifax
c               an INTEGER array.  the number of elements actually used
c               will depend on the factorization of n.  dimensioning
c               ifax for 13 suffices for ALL n less than a million.
c
c              trigs
c               a floating point array of DIMENSION 3*n/2 IF n/2 is
c               even, or 3*n/2+1 IF n/2 is odd.
c
c on output    ifax
c               CONTAINS the factorization of n/2.  ifax(1) is the
c               number of factors, and the factors themselves are stored
c               in ifax(2),ifax(3),...  IF fftfax is called with n odd,
c               or IF n has ANY prime factors greater than 5, ifax(1)
c               is set to -99.
c
c              trigs
c               an array of trignomentric FUNCTION values subsequently
c               used by the fft routines.
c
c-----------------------------------------------------------------------
c
c SUBROUTINE fft991 (a,work,trigs,ifax,inc,jump,n,m,isign)
c                       and
c SUBROUTINE fft99 (a,work,trigs,ifax,inc,jump,n,m,isign)
c
c purpose      perform a number of simultaneous REAL/half-complex
c              periodic fourier transforms or corresponding inverse
c              transforms, using ordinary spatial order of gridpoint
c              values (fft991) or explicit cyclic continuity in the
c              gridpoint values (fft99).  given a set
c              of REAL data vectors, the package returns a set of
c              "half-complex" fourier coefficient vectors, or vice
c              versa.  the length of the transforms must be an even
c              number that has no other factors except possibly powers
c              of 2, 3, and 5.  these version of fft991 and fft99 are
c              optimized for USE on the cray-1.
c
c argument     a(m*(n+2)), work(m*(n+1)), trigs(3*n/2+1), ifax(13)
c dimensions
c
c arguments
c
c on input     a
c               an array of length m*(n+2) containing the input data
c               or coefficient vectors.  this array is overwritten by
c               the results.
c
c              work
c               a work array of DIMENSION m*(n+1)
c
c              trigs
c               an array set up by fftfax, which must be called first.
c
c              ifax
c               an array set up by fftfax, which must be called first.
c
c              inc
c               the increment (in words) between successive elements of
c               each data or coefficient vector (e.g.  inc=1 for
c               consecutively stored data).
c
c              jump
c               the increment (in words) between the first elements of
c               successive data or coefficient vectors.  on the cray-1,
c               try to arrange data so that jump is not a multiple of 8
c               (to avoid memory bank conflicts).  for clarification of
c               inc and jump, see the examples below.
c
c              n
c               the length of each transform (see definition of
c               transforms, below).
c
c              m
c               the number of transforms to be done simultaneously.
c
c              isign
c               = +1 for a transform from fourier coefficients to
c                    gridpoint values.
c               = -1 for a transform from gridpoint values to fourier
c                    coefficients.
c
c on output    a
c               IF isign = +1, and m coefficient vectors are supplied
c               each containing the sequence:
c
c               a(0),b(0),a(1),b(1),...,a(n/2),b(n/2)  (n+2 values)
c
c               THEN the result consists of m data vectors each
c               containing the corresponding n+2 gridpoint values:
c
c               for fft991, x(0), x(1), x(2),...,x(n-1),0,0.
c               for fft99, x(n-1),x(0),x(1),x(2),...,x(n-1),x(0).
c                   (explicit cyclic continuity)
c
c               when isign = +1, the transform is defined by:
c                 x(j)=SUM(k=0,...,n-1)(c(k)*EXP(2*i*j*k*pi/n))
c                 WHERE c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)
c                 and i=SQRT (-1)
c
c               IF isign = -1, and m data vectors are supplied each
c               containing a sequence of gridpoint values x(j) as
c               defined above, THEN the result consists of m vectors
c               each containing the corresponding fourier cofficients
c               a(k), b(k), 0 .le. k .le n/2.
c
c               when isign = -1, the inverse transform is defined by:
c                 c(k)=(1/n)*SUM(j=0,...,n-1)(x(j)*EXP(-2*i*j*k*pi/n))
c                 WHERE c(k)=a(k)+i*b(k) and i=SQRT(-1)
c
c               a CALL with isign=+1 followed by a CALL with isign=-1
c               (or vice versa) returns the original data.
c
c               note: the fact that the gridpoint values x(j) are REAL
c               implies that b(0)=b(n/2)=0.  for a CALL with isign=+1,
c               it is not actually necessary to supply these zeros.
c
c examples      given 19 data vectors each of length 64 (+2 for explicit
c               cyclic continuity), compute the corresponding vectors of
c               fourier coefficients.  the data may, for example, be
c               arranged like this:
c
c first data   a(1)=    . . .                a(66)=             a(70)
c vector       x(63) x(0) x(1) x(2) ... x(63) x(0)  (4 empty locations)
c
c second data  a(71)=   . . .                                  a(140)
c vector       x(63) x(0) x(1) x(2) ... x(63) x(0)  (4 empty locations)
c
c               and so on.  here inc=1, jump=70, n=64, m=19, isign=-1,
c               and fft99 should be used (because of the explicit cyclic
c               continuity).
c
c               alternatively the data may be arranged like this:
c
c                first         second                          last
c                data          data                            data
c                vector        vector                          vector
c
c                 a(1)=         a(2)=                           a(19)=
c
c                 x(63)         x(63)       . . .               x(63)
c        a(20)=   x(0)          x(0)        . . .               x(0)
c        a(39)=   x(1)          x(1)        . . .               x(1)
c                  .             .                               .
c                  .             .                               .
c                  .             .                               .
c
c               in which CASE we have inc=19, jump=1, and the remaining
c               parameters are the same as before.  in either CASE, each
c               coefficient vector overwrites the corresponding input
c               data vector.
c
c-----------------------------------------------------------------------
c
c     SUBROUTINE 'fft99' - multiple fast REAL periodic transform
c     corresponding to old scalar routine fft9
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
c         x(n-1),x(0),x(1),x(2),...,x(n),x(0)
c         i.e. explicit cyclic continuity; (n+2) locations required
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
      ibase=inc+1
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
      ia=inc+1
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
      jbase=ia
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
c     fill in cyclic boundary points
  110 CONTINUE
      ia=1
      ib=n*inc+1
cdir$ ivdep
      DO 120 l=1,lot
      a(ia)=a(ib)
      a(ib+inc)=a(ia+inc)
      ia=ia+jump
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

      END SUBROUTINE fft99
