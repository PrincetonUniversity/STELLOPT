      SUBROUTINE vpassm(a, b, c, d, trigs, inc1, inc2, inc3, inc4, lot,
     1   n, ifac, la)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: SIN36=0.587785252292473_dp
      REAL(rprec), PARAMETER :: COS36=0.809016994374947_dp
      REAL(rprec), PARAMETER :: SIN72=0.951056516295154_dp
      REAL(rprec), PARAMETER :: COS72=0.309016994374947_dp
      REAL(rprec), PARAMETER :: SIN60=0.866025403784437_dp
      REAL(rprec), PARAMETER :: p5 = 0.5_dp
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER inc1, inc2, inc3, inc4, lot, n, ifac, la
      REAL(rprec), DIMENSION(*) :: a, b, c, d
      REAL(rprec), DIMENSION(*) :: trigs
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: m, iink, jink, jump, ibase, jbase, igo, ia, ja, ib, jb
     1   , l, i, j, ijk, la1, k, kb, ic, jc, kc, id, jd, kd, ie, je, ke
      REAL(rprec) :: c1,s1,c2,s2,c3,s3,c4,s4
C-----------------------------------------------
c
c     SUBROUTINE 'vpassm' - multiple version of 'vpassa'
c     performs one pass through data
c     as part of multiple complex (inverse) fft routine
c     a is first real input vector
c     b is first imaginary input vector
c     c is first real output vector
c     d is first imaginary output vector
c     trigs is precalculated table of sines & cosines
c     inc1 is addressing increment for a and b
c     inc2 is addressing increment for c and d
c     inc3 is addressing increment between a's & b's
c     inc4 is addressing increment between c's & d's
c     lot is the number of vectors
c     n is length of vectors
c     ifac is current factor of n
c     la is product of previous factors
c
c
      m=n/ifac
      iink=m*inc1
      jink=la*inc2
      jump=(ifac-1)*jink
      ibase=0
      jbase=0
      igo=ifac-1
      IF (igo .gt. 4) RETURN
      GOTO (10,50,90,130),igo
c
c     coding for factor 2
c
   10 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      DO 20 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      DO 15 ijk=1,lot
      c(ja+j)=a(ia+i)+a(ib+i)
      d(ja+j)=b(ia+i)+b(ib+i)
      c(jb+j)=a(ia+i)-a(ib+i)
      d(jb+j)=b(ia+i)-b(ib+i)
      i=i+inc3
      j=j+inc4
   15 CONTINUE
      ibase=ibase+inc1
      jbase=jbase+inc2
   20 CONTINUE
      IF (la.eq.m) RETURN
      la1=la+1
      jbase=jbase+jump
      DO 40 k=la1,m,la
      kb=k+k-2
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      DO 30 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      DO 25 ijk=1,lot
      c(ja+j)=a(ia+i)+a(ib+i)
      d(ja+j)=b(ia+i)+b(ib+i)
      c(jb+j)=c1*(a(ia+i)-a(ib+i))-s1*(b(ia+i)-b(ib+i))
      d(jb+j)=s1*(a(ia+i)-a(ib+i))+c1*(b(ia+i)-b(ib+i))
      i=i+inc3
      j=j+inc4
   25 CONTINUE
      ibase=ibase+inc1
      jbase=jbase+inc2
   30 CONTINUE
      jbase=jbase+jump
   40 CONTINUE
      RETURN
c
c     coding for factor 3
c
   50 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      DO 60 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      DO 55 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
      c(jb+j)=(a(ia+i)-p5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))
      c(jc+j)=(a(ia+i)-p5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))
      d(jb+j)=(b(ia+i)-p5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i)))
      d(jc+j)=(b(ia+i)-p5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i)))
      i=i+inc3
      j=j+inc4
   55 CONTINUE
      ibase=ibase+inc1
      jbase=jbase+inc2
   60 CONTINUE
      IF (la.eq.m) RETURN
      la1=la+1
      jbase=jbase+jump
      DO 80 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      DO 70 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      DO 65 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
      c(jb+j)=
     *    c1*((a(ia+i)-p5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i))))
     *   -s1*((b(ia+i)-p5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
      d(jb+j)=
     *    s1*((a(ia+i)-p5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i))))
     *   +c1*((b(ia+i)-p5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
      c(jc+j)=
     *    c2*((a(ia+i)-p5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i))))
     *   -s2*((b(ia+i)-p5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
      d(jc+j)=
     *    s2*((a(ia+i)-p5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i))))
     *   +c2*((b(ia+i)-p5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
      i=i+inc3
      j=j+inc4
   65 CONTINUE
      ibase=ibase+inc1
      jbase=jbase+inc2
   70 CONTINUE
      jbase=jbase+jump
   80 CONTINUE
      RETURN
c
c     coding for factor 4
c
   90 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      id=ic+iink
      jd=jc+jink
      DO 100 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      DO 95 ijk=1,lot
      c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
      c(jc+j)=(a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))
      d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
      d(jc+j)=(b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i))
      c(jb+j)=(a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))
      c(jd+j)=(a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))
      d(jb+j)=(b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i))
      d(jd+j)=(b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i))
      i=i+inc3
      j=j+inc4
   95 CONTINUE
      ibase=ibase+inc1
      jbase=jbase+inc2
  100 CONTINUE
      IF (la.eq.m) RETURN
      la1=la+1
      jbase=jbase+jump
      DO 120 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      kd=kc+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      c3=trigs(kd+1)
      s3=trigs(kd+2)
      DO 110 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      DO 105 ijk=1,lot
      c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
      d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
      c(jc+j)=
     *    c2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i)))
     *   -s2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
      d(jc+j)=
     *    s2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i)))
     *   +c2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
      c(jb+j)=
     *    c1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i)))
     *   -s1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
      d(jb+j)=
     *    s1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i)))
     *   +c1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
      c(jd+j)=
     *    c3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i)))
     *   -s3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
      d(jd+j)=
     *    s3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i)))
     *   +c3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
      i=i+inc3
      j=j+inc4
  105 CONTINUE
      ibase=ibase+inc1
      jbase=jbase+inc2
  110 CONTINUE
      jbase=jbase+jump
  120 CONTINUE
      RETURN
c
c     coding for factor 5
c
  130 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      id=ic+iink
      jd=jc+jink
      ie=id+iink
      je=jd+jink
      DO 140 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      DO 135 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
      c(jb+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *  -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
      c(je+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *  +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
      d(jb+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *  +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
      d(je+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *  -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
      c(jc+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *  -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
      c(jd+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *  +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
      d(jc+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *  +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
      d(jd+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *  -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
      i=i+inc3
      j=j+inc4
  135 CONTINUE
      ibase=ibase+inc1
      jbase=jbase+inc2
  140 CONTINUE
      IF (la.eq.m) RETURN
      la1=la+1
      jbase=jbase+jump
      DO 160 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      kd=kc+kb
      ke=kd+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      c3=trigs(kd+1)
      s3=trigs(kd+2)
      c4=trigs(ke+1)
      s4=trigs(ke+2)
      DO 150 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      DO 145 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
      c(jb+j)=
     *    c1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *      -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))
     *   -s1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *      +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      d(jb+j)=
     *    s1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *      -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))
     *   +c1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *      +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      c(je+j)=
     *    c4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *      +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))
     *   -s4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *      -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      d(je+j)=
     *    s4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *      +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))
     *   +c4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *      -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      c(jc+j)=
     *    c2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *      -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))
     *   -s2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *      +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      d(jc+j)=
     *    s2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *      -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))
     *   +c2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *      +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      c(jd+j)=
     *    c3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *      +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))
     *   -s3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *      -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      d(jd+j)=
     *    s3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *      +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))
     *   +c3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *      -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      i=i+inc3
      j=j+inc4
  145 CONTINUE
      ibase=ibase+inc1
      jbase=jbase+inc2
  150 CONTINUE
      jbase=jbase+jump
  160 CONTINUE

      END SUBROUTINE vpassm
