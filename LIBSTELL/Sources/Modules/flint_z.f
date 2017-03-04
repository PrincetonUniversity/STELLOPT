C******************** START FILE FLINT_Z.FOR ; GROUP SIGS2 ********************
C==============================================================================
C
C  FLINT_Z
C   INTERPOLATE ON A 2D ARRAY
C
       SUBROUTINE FLINT_Z(X,Y,F,N,MTX,XMIN,XMAX,XLR,NX,YMIN,YMAX,YLR,
     > NY,IWARN)
C
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER n,nx,ny,it,iwarn,ix,iy,ixp1,iyp1
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 f,mtx,x,y,xmin,xmax,xlr,ymin,ymax,ylr,zx,zy,zzx,zlogx,
     > zzy,zlogy,zf00,zf01,zf10,zf11
C============
      DIMENSION MTX(NX,NY)
      DIMENSION X(N)
      DIMENSION Y(N)
      DIMENSION F(N)
c
c Interpolates the table MTX(NX,NY) to find the value at (x,y).
c If XLR>0, then the X grid is equally spaced on a logarithmic scale:
c
c   X(i) = XMIN * (exp(XLR)) ** ( (i-1)/(NX-1) )
c
c   So X(i) ranges from XMIN to XMIN*exp(XLR), and XLR=log(xmax/xmin).
c
c For XLR<=0, the X grid is equally spaced on a linear scale from
c XMIN to XMAX.
c
c For YLR>0, Y is logarithmically spaced from YMIN to YMIN*exp(YLR)
c and for YLR<=0, Y is linearly spaced from YMIN to YMAX.
C
C Modified 7/26/00 (PCR) - vectorized.
C
      integer idebug
      common/flint_debug/ idebug
C
C
      iwarn=0
C
      do 70 it=1,n
C
         ZX=X(it)
         ZY=Y(it)
C
         if(zx .lt. xmin .or. zx/xmin .le. 0.0) then
            zzx=0.0D0
         else
 
            IF(XLR.GT.0.0D0) THEN
               ZLOGX=log(ZX/XMIN)
               ZZX=1.0D0+ZLOGX/XLR * (NX-1)
            ELSE
               ZZX=1.0D0+(ZX-XMIN)/(XMAX-XMIN) * (NX-1)
            ENDIF
 
         endif
C
         if(zy .lt. ymin .or. zy/ymin .le. 0.0) then
            zzy=0.0D0
         else
 
            IF(YLR.GT.0.0D0) THEN
               ZLOGY=log(ZY/YMIN)
               ZZY=1.0D0+ZLOGY/YLR * (NY-1)
            ELSE
               ZZY=1.0D0+(ZY-YMIN)/(YMAX-YMIN) * (NY-1)
            ENDIF
 
         endif
C
         IX=ZZX
         IY=ZZY
C  CHECK FOR POINT OUT OF BOUNDS OF ARRAY
         IF(IX.GE.1) GO TO 10
C
         IF(ZX.LT.(XMIN - 5.0D-7)) THEN
            IWARN=1
         ENDIF
         IX=1
         ZZX=1.0D0
C
 10      CONTINUE
         IF(IY.GE.1) GO TO 20
C
         IF(ZY.LT.(YMIN - 5.0D-7)) THEN
            IWARN=1
         ENDIF
         IY=1
         ZZY=1.0D0
C
 20      CONTINUE
         IF(IX.LT.NX) GO TO 30
C
         IF(ZX.GT.(XMAX + 5.0D-7)) THEN
            IWARN=1
         ENDIF
         IX=NX-1
         ZZX=NX
C
 30      CONTINUE
         IF(IY.LT.NY) GO TO 40
C
         IF(ZY.GT.(YMAX + 5.0D-7)) THEN
            IWARN=1
         ENDIF
         IY=NY-1
         ZZY=NY
C
C---------------------------
C  OK
C
 40      CONTINUE
         IXP1=IX+1
         IYP1=IY+1
         ZZX=ZZX-IX
         ZZY=ZZY-IY
         ZF00=(1.D0-ZZX)*(1.D0-ZZY)
         ZF01=(1.D0-ZZX)*ZZY
         ZF10=ZZX*(1.D0-ZZY)
         ZF11=ZZX*ZZY
C
C  INTERPOLATE
C
         F(it)=ZF00*MTX(IX,IY)+ZF01*MTX(IX,IYP1)+
     >      ZF10*MTX(IXP1,IY)+ZF11*MTX(IXP1,IYP1)
c
         if(idebug.eq.99) then
            write(6,*) 'flin1_z:  zx,zy=',zx,zy,' ix,iy=',ix,iy
            write(6,*) '  mtx(ix,iy:iyp1) = ',mtx(ix,iy:iyp1)
            write(6,*) '  mtx(ixp1,iy:iyp1) = ',mtx(ixp1,iy:iyp1)
            write(6,*) '  it,F(it)=',it,F(it)
         endif
C
 70   continue
C
      idebug=0  ! clear debug flag
C
      END
C
C******************** END FILE FLINT.FOR ; GROUP SIGS2 ************************
! 26Jul2000 fgtok -s r8_precision.sub "r8con.csh conversion"
