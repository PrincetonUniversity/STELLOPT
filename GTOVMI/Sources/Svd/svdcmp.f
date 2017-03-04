      SUBROUTINE svdcmp(a,m,n,mp,np,w,v)
      USE precision
      IMPLICIT NONE
      INTEGER :: m,mp,n,np
      REAL(rprec) :: a(mp,np),v(np,np),w(np)
      INTEGER, PARAMETER :: nmax=252
cu    uses pythag
      INTEGER i,its,j,jj,k,l,nm
      REAL(rprec) anorm,c,f,g,h,s,scaly,x,y,z,rv1(nmax),pythag
      g=0
      scaly=0
      anorm=0
      DO 25 i=1,n
        l=i+1
        rv1(i)=scaly*g
        g=0
        s=0
        scaly=0
        IF(i.le.m)then
          DO 11 k=i,m
            scaly=scaly+ABS(a(k,i))
11        CONTINUE
          IF(scaly.ne.0.0)then
            DO 12 k=i,m
              a(k,i)=a(k,i)/scaly
              s=s+a(k,i)*a(k,i)
12          CONTINUE
            f=a(i,i)
            g=-SIGN(SQRT(s),f)
            h=f*g-s
            a(i,i)=f-g
            DO 15 j=l,n
              s=0
              DO 13 k=i,m
                s=s+a(k,i)*a(k,j)
13            CONTINUE
              f=s/h
              DO 14 k=i,m
                a(k,j)=a(k,j)+f*a(k,i)
14            CONTINUE
15          CONTINUE
            DO 16 k=i,m
              a(k,i)=scaly*a(k,i)
16          CONTINUE
          ENDIF
        ENDIF
        w(i)=scaly *g
        g=0
        s=0
        scaly=0
        IF((i.le.m).and.(i.ne.n))then
          DO 17 k=l,n
            scaly=scaly+ABS(a(i,k))
17        CONTINUE
          IF(scaly.ne.0.0)then
            DO 18 k=l,n
              a(i,k)=a(i,k)/scaly
              s=s+a(i,k)*a(i,k)
18          CONTINUE
            f=a(i,l)
            g=-SIGN(SQRT(s),f)
            h=f*g-s
            a(i,l)=f-g
            DO 19 k=l,n
              rv1(k)=a(i,k)/h
19          CONTINUE
            DO 23 j=l,m
              s=0
              DO 21 k=l,n
                s=s+a(j,k)*a(i,k)
21            CONTINUE
              DO 22 k=l,n
                a(j,k)=a(j,k)+s*rv1(k)
22            CONTINUE
23          CONTINUE
            DO 24 k=l,n
              a(i,k)=scaly*a(i,k)
24          CONTINUE
          ENDIF
        ENDIF
        anorm=max(anorm,(ABS(w(i))+ABS(rv1(i))))
25    CONTINUE
      DO 32 i=n,1,-1
        IF(i.lt.n)then
          IF(g.ne.0.0)then
            DO 26 j=l,n
              v(j,i)=(a(i,j)/a(i,l))/g
26          CONTINUE
            DO 29 j=l,n
              s=0
              DO 27 k=l,n
                s=s+a(i,k)*v(k,j)
27            CONTINUE
              DO 28 k=l,n
                v(k,j)=v(k,j)+s*v(k,i)
28            CONTINUE
29          CONTINUE
          ENDIF
          DO 31 j=l,n
            v(i,j)=0
            v(j,i)=0
31        CONTINUE
        ENDIF
        v(i,i)=1.0
        g=rv1(i)
        l=i
32    CONTINUE
      DO 39 i=min(m,n),1,-1
        l=i+1
        g=w(i)
        DO 33 j=l,n
          a(i,j)=0
33      CONTINUE
        IF(g.ne.0.0)then
          g=1.0/g
          DO 36 j=l,n
            s=0
            DO 34 k=l,m
              s=s+a(k,i)*a(k,j)
34          CONTINUE
            f=(s/a(i,i))*g
            DO 35 k=i,m
              a(k,j)=a(k,j)+f*a(k,i)
35          CONTINUE
36        CONTINUE
          DO 37 j=i,m
            a(j,i)=a(j,i)*g
37        CONTINUE
        ELSE
          DO 38 j= i,m
            a(j,i)=0
38        CONTINUE
        ENDIF
        a(i,i)=a(i,i)+1.0
39    CONTINUE
      DO 49 k=n,1,-1
        DO 48 its=1,30
          DO 41 l=k,1,-1
            nm=l-1
            IF((ABS(rv1(l))+anorm).eq.anorm)  GOTO 2
            IF((ABS(w(nm))+anorm).eq.anorm)  GOTO 1
41        CONTINUE
1         c=0
          s=1.0
          DO 43 i=l,k
            f=s*rv1(i)
            rv1(i)=c*rv1(i)
            IF((ABS(f)+anorm).eq.anorm) GOTO 2
            g=w(i)
            h=pythag(f,g)
            w(i)=h
            h=1.0/h
            c= (g*h)
            s=-(f*h)
            DO 42 j=1,m
              y=a(j,nm)
              z=a(j,i)
              a(j,nm)=(y*c)+(z*s)
              a(j,i)=-(y*s)+(z*c)
42          CONTINUE
43        CONTINUE
2         z=w(k)
          IF(l.eq.k)then
            IF(z.lt.0.0)then
              w(k)=-z
              DO 44 j=1,n
                v(j,k)=-v(j,k)
44            CONTINUE
            ENDIF
            GOTO 3
          ENDIF
          IF(its.eq.30) PAUSE 'no convergence in svdcmp'
          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
          g=pythag(f,1.0_dbl)
          f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x
          c=1.0
          s=1.0
          DO 47 j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            z=pythag(f,h)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            DO 45 jj=1,n
              x=v(jj,j)
              z=v(jj,i)
              v(jj,j)= (x*c)+(z*s)
              v(jj,i)=-(x*s)+(z*c)
45          CONTINUE
            z=pythag(f,h)
            w(j)=z
            IF(z.ne.0.0)then
              z=1.0/z
              c=f*z
              s=h*z
            ENDIF
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            DO 46 jj=1,m
              y=a(jj,j)
              z=a(jj,i)
              a(jj,j)= (y*c)+(z*s)
              a(jj,i)=-(y*s)+(z*c)
46          CONTINUE
47        CONTINUE
          rv1(l)=0
          rv1(k)=f
          w(k)=x
48      CONTINUE
3       CONTINUE
49    CONTINUE
      RETURN
      END SUBROUTINE svdcmp
      FUNCTION pythag(a,b)
      USE precision
      IMPLICIT NONE
      REAL(rprec) a,b,pythag
      REAL(rprec) absa,absb
      absa=ABS(a)
      absb=ABS(b)
      IF(absa.gt.absb)then
        pythag=absa*SQRT(1.+(absb/absa)**2)
      ELSE
        IF(absb.eq.0.)then
          pythag=0.
        ELSE
          pythag=absb*SQRT(1.+(absa/absb)**2)
        ENDIF
      ENDIF
      RETURN
      END FUNCTION pythag

      SUBROUTINE svbksb(u,w,v,m,n,mp,np,b,x)
      USE precision
      IMPLICIT NONE
      INTEGER m,mp,n,np
      REAL(rprec) b(mp),u(mp,np),v(np,np),w(np),x(np)
      INTEGER, PARAMETER :: nmax=252
      INTEGER i,j,jj
      REAL s,tmp(nmax)
      DO 12 j=1,n
        s=0.
        IF(w(j).ne.0.)then
          DO 11 i=1,m
            s=s+u(i,j)*b(i)
11        CONTINUE
          s=s/w(j)
        ENDIF
        tmp(j)=s
12    CONTINUE
      DO 14 j=1,n
        s=0.
        DO 13 jj=1,n
          s=s+v(j,jj)*tmp(jj)
13      CONTINUE
        x(j)=s
14    CONTINUE
      RETURN
      END SUBROUTINE svbksb



