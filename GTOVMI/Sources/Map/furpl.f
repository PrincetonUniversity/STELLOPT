      SUBROUTINE furplm(p,grsq,x,y,v,xp,yp,gsq,np,nx,ny,ntx,nw,nh,is,js)
c ***********************************************************************
c given an x-y grid with a variable p defined on the grid this routine
c finds the level p curves of value v.  that is, it finds np points
c (xp(k),yp(k)) such that p(xp(k),yp(k)) is CLOSE to v.
c id is an INTEGER between 1 and 12 or 15 which identifies the option
c USEd to find the last point.
c IF in=1 the curve is moving up and IF in=0 the curve is moving DOwn.
c IF il=1 the curve is moving left and IF il=0 the curve is moving right.
c IF the last poINT was a grid poINT id=15 and igp is an INTEGER between
c 1 and 12 which identifies the grid poINT option USEd to find the last
c point.
c IF ig=1 the grid poINT was entered from the left and IF ig=0 the grid
c poINT was entered from the right.
c iq is an indicator USEd in that part of the PROGRAM which searches the
c grid squares surrounding a grid point.
c IF ntx is less than zero search for the first poINT starts from the
c top of the grid and goes DOwn.
c IF np eq -1 the search starts in the center of the grid
c IF ntx gt 0 the search starts at the bottom of the grid
c IF you would like further details see joanne helton.
c ***********************************************************************
      USE precision
      IMPLICIT NONE
      INTEGER :: js, is, ntx, ny, nx, np, nw, nh, nt, nx2, ny2
      INTEGER :: jst, jdec, kprt, j, i, k, in, iq, id, i1, jy,
     .	ix, igp, l, iyy, ixx, ig, il
      REAL(rprec) ::  xp(*),yp(*)
      REAL(rprec) :: grsq(nw,nh),gsq(*)
      REAL(rprec) ::  p(nw,*),x(*),y(*)
      REAL(rprec) :: a, b, c, f, yq, v, dx, dy, dp, d, e, ff, ws,
     .	xl, yl, x1, y1
      a=0.
      b=0.
      c=0.
      f=0.
      yq=0.
      nt=ntx
      jst=js
      jdec=1
      IF(np.ne.-1) go to 1999
      jst=ny/2+1
      jdec=1
      nt=iABS(nt)
 1999 IF(nt.gt.0) go to 2004
      jst=ny-1
      jdec=-1
      nt=-nt
 2004 CONTINUE
      nx2=nx/2
      ny2=ny/2
      dx=x(2)-x(1)
      dy=y(2)-y(1)
      dp=ABS(p(nx2,ny2))
      d=10.**(-1*8)*dp
      e=10.**(-1*4)*dx*3.0
      ff=10.**(-1*6)*dx
      kprt=0
2001  CONTINUE
      j=jst
      i=is
      k=0
      in=0
      iq=0
c ---------------------
c we find a first point
c ---------------------
   13 IF((v-p(i,j))*(p(i+1,j)-v).gt.0.) go to 30
      IF(ABS(p(i,j)-v).ge.d) go to 32
      in=1
      ig=1
      go to 15
   32 i=i+1
      IF(i.eq.nx) go to 25
      go to 13
   25 j=j+jdec
      IF(jdec.lt.0) go to 2005
      IF(j.gt.ny) go to 10
      i=1
      go to 13
 2005 IF(j.lt.1) go to 10
      i=1
      go to 13
   30 IF(i.gt.(nx-2)) go to 74
      IF(i.eq.1) go to 88
      xp(1)=x(i)
      CALL fit(2,x(i-1),x(i),x(i+1),x(i+2),p(i-1,j),p(i,j),p(i+1,j),
     1         p(i+2,j),xp(1),v,yq)
      ws=xp(1)
      CALL cubic(x(i-1),x(i),x(i+1),x(i+2),ws,a,b,c,f)
      gsq(1)=a*grsq(i-1,j)+b*grsq(i,j)+c*grsq(i+1,j)+f*grsq(i+2,j)
      go to 75
   88 xp(1)=x(2)
      CALL fit(2,x(1),x(2),x(3),x(4),p(1,j),p(2,j),p(3,j),p(4,j),
     1         xp(1),v,yq)
      ws=xp(1)
      CALL cubic(x(1),x(2),x(3),x(4),ws,a,b,c,f)
      gsq(1)=a*grsq(1,j)+b*grsq(2,j)+c*grsq(3,j)+f*grsq(4,j)
      go to 75
   74 xp(1)=x(nx-2)
      CALL fit(2,x(nx-3),x(nx-2),x(nx-1),x(nx),p(nx-3,j),p(nx-2,j),
     1         p(nx-1,j),p(nx,j),xp(1),v,yq)
      ws=xp(1)
      CALL cubic(x(nx-3),x(nx-2),x(nx-1),x(nx),ws,a,b,c,f)
      gsq(1)=a*grsq(nx-3,j)+b*grsq(nx-2,j)+c*grsq(nx-1,j)+f*grsq(nx,j)
   75 yp(1)=y(j)
      in=1
      IF(xp(1).le.x(i)) xp(1)=x(i)+ff
      IF(xp(1).ge.x(i+1)) xp(1)=x(i+1)-ff
      id=20
      il=1
      IF(jdec.lt.0) il=0
      k=1
      IF(kprt.eq.1) WRITE(66,102)k,xp(k),yp(k),p(i,j),i,j,ix,jy,ig,in,id
     x      ,iq,il
  102 FORMAT(i5,3e13.6,10i5)
      go to 5
c ---------------------------------------------------------------
c in the block below we determine WHERE the next poINT is located
c ---------------------------------------------------------------
c ---------------------------------------------------------------------
c the curve is entering the grid square from the right or from the left
c IF il equals 1 the curve is moving left
c ---------------------------------------------------------------------
    4 IF(il.eq.1) go to 21
c --------------------------------------------------
c the following options are for a curve moving right
c --------------------------------------------------
      IF(i+1.gt.nx) go to 10
c --------------------------
c first we check grid points
c --------------------------
      IF(ABS(v-p(i+1,j+1)).ge.d) go to 45
      ig=1
      in=1
      i=i+1
      j=j+1
      igp=11
      go to 15
   45 IF(ABS(v-p(i+1,j)).ge.d) go to 33
      ig=1
      in=0
      i=i+1
      igp=12
      go to 15
c --------------------------------
c now we check between grid points
c --------------------------------
   33 IF((v-p(i,j))*(p(i+1,j)-v).le.0.) go to 50
      ix=i
      jy=j
      id=1
      go to 1
   50 IF((v-p(i,j+1))*(p(i+1,j+1)-v).le.0) go to 51
      ix=i
      jy=j+1
      id=2
      go to 1
   51 IF((v-p(i+1,j))*(p(i+1,j+1)-v).le.0.) go to 52
      ix=i+1
      jy=j
      id=3
      go to 7
   52 CONTINUE
c -------------------------------------------------
c the following options are for a curve moving left
c -------------------------------------------------
      IF(i-1.lt.0) go to 10
c --------------------------
c first we check grid points
c --------------------------
   21 IF(ABS(v-p(i-1,j)).ge.d) go to 46
      in=0
      ig=0
      i=i-1
      igp=1
      go to 15
   46 IF(ABS(v-p(i-1,j+1)).ge.d) go to 34
      ig=0
      in=1
      i=i-1
      j=j+1
      igp=2
      go to 15
c --------------------------------
c now we check between grid points
c --------------------------------
   34 IF((v-p(i,j+1))*(p(i-1,j+1)-v).le.0.) go to 53
      ix=i-1
      jy=j+1
      id=4
      go to 1
   53 IF((v-p(i-1,j+1))*(p(i-1,j)-v).le.0) go to 54
      ix=i-1
      jy=j
      id=5
      go to 7
   54 IF((v-p(i,j))*(p(i-1,j)-v).le.0.) go to 55
      ix=i-1
      jy=j
      id=6
      go to 1
   55 CONTINUE
c ---------------------------------------------------------------------
c the curve is entering the grid square from the top or from the bottom
c or from a grid point
c IF ig equals 1 the grid poINT was entered from the left
c IF in equals 1 the curve is moving up
c ---------------------------------------------------------------------
    5 IF(id.eq.15.and.il.eq.1) i=i-1
   72 iq=0
   62 IF(in.eq.1) go to 12
c -------------------------------------------------
c the following options are for a curve moving DOwn
c -------------------------------------------------
      IF(j-1.lt.1) go to 10
   11 IF(ABS(v-p(i+1,j-1)).ge.d) go to 42
      ig=1
      in=0
      i=i+1
      j=j-1
      igp=3
      go to 15
   42 IF(ABS(v-p(i,j-1)).ge.d) go to 17
      ig=0
      in=0
      j=j-1
      igp=5
      go to 15
c --------------------------------
c now we check between grid points
c --------------------------------
   17 IF((v-p(i+1,j-1))*(p(i+1,j)-v).le.0.) go to 56
      IF(ABS(v-p(i+1,j)).le.d) go to 56
      IF(id.eq.15.and.igp.eq.7) go to 56
      ix=i+1
      jy=j-1
      id=7
      go to 7
   56 IF((v-p(i,j-1))*(p(i+1,j-1)-v).le.0.) go to 57
      ix=i
      jy=j-1
      id=8
      go to 1
   57 IF((v-p(i,j))*(p(i,j-1)-v).le.0.) go to 44
      IF(ABS(v-p(i,j)).le.d) go to 44
      ix=i
      jy=j-1
      id=9
      go to 7
   44 IF(iq.ne.1) go to 48
      in=1
      IF(il.eq.1) i=i+1
      IF(il.eq.0) i=i-1
      go to 72
   48 iq=1
      in=1
      go to 62
c -----------------------------------------------
c the following options are for a curve moving up
c -----------------------------------------------
   12 IF(j+1.gt.ny) go to 10
      IF(i.eq.nx) go to 10
c --------------------------
c first we check grid points
c --------------------------
      IF(ABS(v-p(i+1,j+1)).ge.d) go to 40
      ig=1
      in=1
      i=i+1
      j=j+1
      igp=7
      go to 15
   40 IF(ABS(v-p(i,j+1)).ge.d) go to 61
      ig=0
      in=1
      j=j+1
      igp=9
      go to 15
c --------------------------------
c now we check between grid points
c --------------------------------
   61 IF((v-p(i,j+1))*(p(i+1,j+1)-v).le.0) go to 59
      ix=i
      jy=j+1
      id=10
      go to 1
   59 IF((v-p(i+1,j))*(p(i+1,j+1)-v).le.0.) go to 60
      IF(ABS(v-p(i+1,j)).le.d) go to 60
      ix=i+1
      jy=j
      id=11
      go to 7
   60 IF((v-p(i,j+1))*(p(i,j)-v).le.0.) go to 41
      IF(ABS(v-p(i,j)).le.d) go to 41
      ix=i
      jy=j
      id=12
      go to 7
   41 IF(iq.ne.1) go to 49
      in=0
      IF(il.eq.1) i=i+1
      IF(il.eq.0) i=i-1
      go to 72
   49 iq=1
      in=0
      go to 62
c ---------------------------------------------------------------------
c in the block below the x and y values of the poINT found are computed
c ---------------------------------------------------------------------
c and stored and control is RETURNed to the appropriate poINT in the
c block above
    1 k=k+1
      IF(k.eq.nt) WRITE(66,106)
  106 FORMAT(/,"***nt points found in furplm***",/)
      IF(k.eq.nt) go to 2000
      IF(ix.eq.1) go to 2
      IF(ix.gt.(nx-2)) go to 3
      xp(k)=x(ix)
      CALL fit(2,x(ix-1),x(ix),x(ix+1),x(ix+2),p(ix-1,jy),p(ix,jy),
     1         p(ix+1,jy),p(ix+2,jy),xp(k),v,yq)
      ws=xp(k)
      CALL cubic(x(ix-1),x(ix),x(ix+1),x(ix+2),ws,a,b,c,f)
 2009 CONTINUE
      gsq(k)=a*grsq(ix-1,jy)+b*grsq(ix,jy)+c*grsq(ix+1,jy)
     1+f*grsq(ix+2,jy)
 2010 CONTINUE
      go to 6
    2 xp(k)=x(2)
      CALL fit(2,x(1),x(2),x(3),x(4),p(1,jy),p(2,jy),p(3,jy),p(4,jy),
     1         xp(k),v,yq)
      ws=xp(k)
      CALL cubic(x(1),x(2),x(3),x(4),ws,a,b,c,f)
      gsq(k)=a*grsq(1,jy)+b*grsq(2,jy)+c*grsq(3,jy)+f*grsq(4,jy)
      go to 6
    3 xp(k)=x(nx-2)
      CALL fit(2,x(nx-3),x(nx-2),x(nx-1),x(nx),p(nx-3,jy),p(nx-2,jy),
     1         p(nx-1,jy),p(nx,jy),xp(k),v,yq)
      ws=xp(k)
      CALL cubic(x(nx-3),x(nx-2),x(nx-1),x(nx),ws,a,b,c,f)
      gsq(k)=a*grsq(nx-3,jy)+b*grsq(nx-2,jy)+c*grsq(nx-1,jy)
     1+f*grsq(nx,jy)
    6 yp(k)=y(jy)
      IF((id.eq.4).or.(id.eq.6)) i=i-1
      IF(xp(k).le.x(ix)) xp(k)=x(ix)+ff
      IF(xp(k).ge.x(ix+1)) xp(k)=x(ix+1)-ff
      IF(k.le.3) go to 20
      IF((ABS(xp(k)-xp(1)).lt.e.and.ABS(yp(k)-yp(1)).lt.e)
     1.and.k.ne.1) go to 10
   20 IF(jy.eq.1.and.(((id.eq.1).or.(id.eq.6)).or.(id.eq.8))) go to 10
      IF(jy.eq.ny.and.(((id.eq.2).or.(id.eq.4)).or.(id.eq.10))) go to 10
      IF(((id.eq.2).or.(id.eq.4)).or.(id.eq.10)) j=j+1
      IF(id.eq.8) j=j-1
      in=1
      IF(yp(k).lt.yp(k-1)) in=0
      il=1
      IF(xp(k).gt.xp(k-1)) il=0
      IF(k > 2) THEN ! avoid reverence to element (0) that does not exist
      IF((ABS(xp(k)-xp(k-2)).lt.e.and.ABS(yp(k)-yp(k-2)).lt.e).and.k.ge.
     x      3) WRITE(66,107)
  107 FORMAT(/,"***poINT same as poINT before last***",/)
      IF((ABS(xp(k)-xp(k-2)).lt.e.and.ABS(yp(k)-yp(k-2)).lt.e)
     1.and.k.ge.3) go to 2000
      ENDIF
      IF(kprt.eq.1) WRITE(66,102)k,xp(k),yp(k),p(i,j),i,j,ix,jy,ig,in,id
     x      ,iq,il
      go to 5
    7 k=k+1
      IF(k.eq.nt) WRITE(66,106)
      IF(k.eq.nt) go to 2000
      xp(k)=x(ix)
      IF(jy.eq.1) go to 8
      IF(jy.gt.(ny-2)) go to 9
      yp(k)=y(jy)
      CALL fit(2,y(jy-1),y(jy),y(jy+1),y(jy+2),p(ix,jy-1),p(ix,jy),
     1         p(ix,jy+1),p(ix,jy+2),yp(k),v,yq)
      ws=yp(k)
      CALL cubic(y(jy-1),y(jy),y(jy+1),y(jy+2),ws,a,b,c,f)
      gsq(k)=a*grsq(ix,jy-1)+b*grsq(ix,jy)+c*grsq(ix,jy+1)
     1+f*grsq(ix,jy+2)
      go to 14
    8 yp(k)=y(2)
      CALL fit(2,y(1),y(2),y(3),y(4),p(ix,1),p(ix,2),p(ix,3),p(ix,4),
     1         yp(k),v,yq)
      ws=yp(k)
      CALL cubic(y(1),y(2),y(3),y(4),ws,a,b,c,f)
      gsq(k)=a*grsq(ix,1)+b*grsq(ix,2)+c*grsq(ix,3)+f*grsq(ix,4)
      go to 14
    9 yp(k)=y(ny-2)
      CALL fit(2,y(ny-3),y(ny-2),y(ny-1),y(ny),p(ix,ny-3),p(ix,ny-2),
     1         p(ix,ny-1),p(ix,ny),yp(k),v,yq)
      ws=yp(k)
      CALL cubic(y(ny-3),y(ny-2),y(ny-1),y(ny),ws,a,b,c,f)
      gsq(k)=a*grsq(ix,ny-3)+b*grsq(ix,ny-2)+c*grsq(ix,ny-1)
     1+f*grsq(ix,ny)
   14 IF(yp(k).le.y(jy)) yp(k)=y(jy)+ff
      IF(yp(k).ge.y(jy+1)) yp(k)=y(jy+1)-ff
      IF(k.le.3) go to 19
      IF((ABS(xp(k)-xp(1)).lt.e.and.ABS(yp(k)-yp(1)).lt.e)
     1.and.k.ne.1) go to 10
   19 IF(ix.eq.1.and.(((id.eq.5).or.(id.eq.9)).or.(id.eq.12))) go to 10
      IF(ix.eq.nx.and.(((id.eq.3).or.(id.eq.7)).or.(id.eq.11))) go to 10
      IF(id.eq.5) i=i-1
      IF((id.eq.7).or.(id.eq.9)) j=j-1
      in=1
      IF(yp(k).lt.yp(k-1)) in=0
      il=1
      IF(xp(k).gt.xp(k-1)) il=0
      IF(((id.eq.3).or.(id.eq.7.and.il.eq.0)).or.(id.eq.11.and.il.eq.0))
     1i=i+1
      IF(kprt.eq.1) WRITE(66,102)k,xp(k),yp(k),p(i,j),i,j,ix,jy,ig,in,id
     x      ,iq,il
      IF(k > 2) THEN ! avoid reverence to element (0) that does not exist
      IF((ABS(xp(k)-xp(k-2)).lt.e.and.ABS(yp(k)-yp(k-2)).lt.e).and.k.ge.
     x      3) WRITE(66,107)! avoid reverence to element (0) that does not exist
      IF((ABS(xp(k)-xp(k-2)).lt.e.and.ABS(yp(k)-yp(k-2)).lt.e)
     1.and.k.ge.3) go to 2000
      ENDIF
      go to 4
   15 k=k+1
      IF(k.eq.nt) WRITE(66,106)
      IF(k.eq.nt) go to 2000
      xp(k)=x(i)
      yp(k)=y(j)
      IF(k > 2) THEN ! avoid reverence to element (0) that does not exist
      IF((ABS(xp(k)-xp(k-2)).lt.e.and.ABS(yp(k)-yp(k-2)).lt.e).and.k.ge.
     x      3) WRITE(66,107)
      IF((ABS(xp(k)-xp(k-2)).lt.e.and.ABS(yp(k)-yp(k-2)).lt.e)
     1.and.k.ge.3) go to 2000
      ENDIF
      gsq(k)=grsq(i,j)
      id=15
      IF(k.le.3) go to 18
      IF((ABS(xp(k)-xp(1)).lt.e.and.ABS(yp(k)-yp(1)).lt.e)
     1.and.k.ne.1) go to 10
   18 IF(((i.eq.1.or.i.eq.nx).and.(j.eq.1.or.j.eq.ny)).and.k.ne.1)
     1go to 10
      in=1
      IF(yp(k).lt.yp(k-1)) in=0
      il=1
      IF(xp(k).gt.xp(k-1)) il=0
      IF(kprt.eq.1) WRITE(66,102)k,xp(k),yp(k),p(i,j),i,j,ix,jy,ig,in,id
     x      ,iq,il,igp
      go to 5
   10 CONTINUE
      np=k
      RETURN
 2000 CONTINUE
c test IF looping on separatrix
      DO 2020 i=1,k
      x1=xp(i)
      y1=yp(i)
      DO 2008 j=i+1,k
      IF(ABS(xp(j)-x1).gt.e) go to 2008
      IF(ABS(yp(j)-y1).le.e) go to 2030
 2008 CONTINUE
 2020 CONTINUE
      go to 2100
 2030 CONTINUE
      k=j-i+1
      DO  j=1,k
         l=i+j-1
         gsq(j)=gsq(l)
         xp(j)=xp(l)
         yp(j)=yp(l)
      END DO
      go to 10
 2100 CONTINUE
      WRITE(6,*)'error in furpl--fur1'
      STOP
      WRITE(66,2002) v,(ixx,x(ixx),ixx=1,nx),(iyy,y(iyy),iyy=1,ny)
2002  FORMAT("  furplm diagnostic center",/,
     .     "  value =",e14.6,/,"   x and y grids follow",/,
     .(i5,e14.6,i5,e14.6,i5,e14.6,i5,e14.6))
      WRITE(66,2003)
2003  FORMAT(//,"  k  ","       xp          yp         p(i,j)   ",
     ."  i    j    ix   jy   ig   in   id   iq   il")
      kprt=1
      go to 2001
      END SUBROUTINE furplm
      SUBROUTINE cubic(p1,p2,p3,p4,p,a,b,c,d)
      USE precision
      IMPLICIT NONE
      REAL(rprec) :: p1, p2, p3, p4, p, a, b, c, d, pp1, pp2, pp3, pp4
      REAL(rprec) :: zero, p1p2, p2p3, p1p3, p1p4, p2p4, p3p4
      INTEGER :: k12, k23, k34
      data  zero/1.e-6_dbl/
cray  lcm (p1),(p2),(p3),(p4)
      p1p2=p1-p2
      k12=SIGN(1.5_dbl,p1p2)
      IF(ABS(p1p2).lt.zero) go to 1
      p2p3=p2-p3
      k23=SIGN(1.5_dbl,p2p3)
      IF(ABS(p2p3).lt.zero) go to 2
      p1p3=p1-p3
      IF(ABS(p1p3).lt.zero) go to 1
      p1p4=p1-p4
      IF (ABS(p1p4).lt.zero) go to 3
      p2p4=p2-p4
      IF (ABS(p2p4).lt.zero) go to 4
      p3p4=p3-p4
      k34=SIGN(1.5_dbl,p3p4)
      IF(ABS(p3p4).lt.zero) go to 4
      IF (iABS(k12+k23+k34).eq.3)  go to 5
      IF (iABS(k23+k34).eq.2) go to 1
      IF (iABS(k12+k23).eq.2) go to 4
5     pp1=p-p1
      pp2=p-p2
      pp3=p-p3
      pp4=p-p4
      a=pp2*pp3*pp4/(p1p2*p1p3*p1p4)
      b=-pp1*pp3*pp4/(p1p2*p2p3*p2p4)
      c=pp1*pp2*pp4/(p1p3*p2p3*p3p4)
      d=-pp1*pp2*pp3/(p1p4*p2p4*p3p4)
      RETURN
1     a=0.
      CALL parab(p2,p3,p4,p,b,c,d)
      RETURN
2     a=0.
      b=1.
      c=0.
      d=0.
      RETURN
3     IF (iABS(k12+k23).eq.2) go to 4
      go to 1
4     d=0.
      CALL parab(p1,p2,p3,p,a,b,c)
      RETURN
      END SUBROUTINE cubic
      SUBROUTINE parab(p1,p2,p3,p,a,b,c)
      USE precision
      IMPLICIT NONE
      REAL(rprec) :: p1, p2, p3, p, a, b, c, zero, yp, y, xp, x,
     .	 p1p2, p2p3, p1p3, p1p4, p2p4, p3p4, pp1, pp2, pp3, pp4
      INTEGER :: k12, k23, k34
      data  zero/1.e-6/
cray  lcm (p1),(p2),(p3)
      p1p2=p1-p2
      p2p3=p2-p3
      p1p3=p1-p3
      pp1=p-p1
      pp2=p-p2
      pp3=p-p3
      k12=SIGN(1.5_dbl,p1p2)
      IF(ABS(p1p2).lt.zero) go to 1
      k23=SIGN(1.5_dbl,p2p3)
      IF(ABS(p2p3).lt.zero) go to 2
      IF(ABS(p1p3).lt.zero) go to 1
      IF (iABS(k12+k23).ne.2) go to 1
      a=pp2*pp3/(p1p2*p1p3)
      b=-pp1*pp3/(p1p2*p2p3)
      c=pp1*pp2/(p1p3*p2p3)
      RETURN
1     a=0.
      b=pp3/p2p3
      c=-pp2/p2p3
      RETURN
2     a=pp3/p1p3
      b=0.
      c=-pp1/p1p3
      RETURN
      END SUBROUTINE parab
      SUBROUTINE fit(k,x1,x2,x3,x4,y1,y2,y3,y4,x,y,yp)
c --------------------------------------------
c this routine is taken from an nyu PROGRAM
c set k=1 to find y,yp and k=2 to find x,yp
c --------------------------------------------
      USE precision
      IMPLICIT NONE
      REAL(rprec) :: x1, x2, x3, x4, y1, y2, y3, y4, x, y, yp, f,
     .	c1, c2, c3, c4, d1, d2, d3, d4, d12, d13, d14, d23, d24, d34,
     . crit, xa, xb, ya, yb, dydx
      INTEGER :: k, iturn, i
      iturn=0
      c1=y1/((x1-x2)*(x1-x3)*(x1-x4))
      c2=y2/((x2-x1)*(x2-x3)*(x2-x4))
      c3=y3/((x3-x1)*(x3-x2)*(x3-x4))
      c4=y4/((x4-x1)*(x4-x2)*(x4-x3))
      IF(k.eq.2) go to 2
   1  d1=x-x1
      d2=x-x2
      d3=x-x3
      d4=x-x4
      d12=d1*d2
      d13=d1*d3
      d14=d1*d4
      d23=d2*d3
      d24=d2*d4
      d34=d3*d4
      f=(c1*d23+c2*d13+c3*d12)*d4+c4*d12*d3
      yp=c1*(d23+d24+d34)+c2*(d13+d14+d34)
     +  +c3*(d12+d14+d24)+c4*(d12+d13+d23)
      IF(k.eq.2) go to 3
      y=f
      RETURN
    2 IF(y.ge.min(y1,y2,y3,y4).and.y.le.max(y1,y2,y3,y4)) GOTO 4
   21 CONTINUE
      WRITE(6,11)x1,x2,x3,x4,y1,y2,y3,y4,y
      x=0.
      WRITE(6,500)
  500 FORMAT("error in sub fit at label 21")
      WRITE(6,*)'error in furpl--fit1'
      STOP
      RETURN
  11  FORMAT(" fit",9e14.5)
   4  crit=(ABS(y1)+ABS(y2)+ABS(y3)+ABS(y4))*1.e-05
      i=0
      xa=x1
      xb=x2
      ya=y1
      yb=y2
      IF((x-xa)*(x-xb).lt.0.) go to 10
      xa=x2
      xb=x3
      ya=y2
      yb=y3
      IF((x-xa)*(x-xb).lt.0.) go to 10
      xa=x3
      xb=x4
      ya=y3
      yb=y4
      IF((x-xa)*(x-xb).lt.0.) go to 10
  12  xa=x2
      ya=y2
      xb=x3
      yb=y3
      x=(xa+xb)/2.
      IF((y-ya)*(y-yb).lt.0.) go to 1
      xa=x1
      ya=y1
      xb=x2
      yb=y2
      x=(xa+xb)/2.
      IF((y-ya)*(y-yb).lt.0.) go to 1
      xa=x3
      ya=y3
      xb=x4
      yb=y4
      x=(xa+xb)/2.
      IF((y-ya)*(y-yb).lt.0.) go to 1
      IF(y.ne.y1) go to 13
      x=x1
      go to 1
   13 IF(y.ne.y2) go to 14
      x=x2
      go to 1
   14 IF(y.ne.y3) go to 15
      x=x3
      go to 1
   15 IF(y.ne.y4) go to 16
      x=x4
      go to 1
   16 CONTINUE
      WRITE(6,11)x,y
      go to 21
3     IF(ABS(f-y).lt.crit) iturn=1
      IF(i.eq.1) go to 7
      dydx=(yb-ya)/(xb-xa)
      IF(ABS(yp-dydx).lt..2*ABS(yp)) go to 7
      IF((f-y)*(ya-y).lt.0.) go to 5
      xa=x
      ya=f
      go to 6
   5  xb=x
      yb=f
   6  x=(xa+xb)/2.
      i=1
      go to 1
    7 IF((f-y)*(ya-y).lt.0.) go to 8
      xa=x
      ya=f
      go to 9
   8  xb=x
      yb=f
   9  dydx=(yb-ya)/(xb-xa)
      IF(ABS(yp-dydx).lt..2*ABS(yp))dydx=yp
      x=x-(f-y)/dydx
      IF(iturn.eq.1) RETURN
      i=0
      go to 1
   10 IF((y-ya)*(y-yb).lt.0.) go to 1
      go to 12
      END SUBROUTINE fit
