C******************** START FILE FLIN1_Z.FOR ; GROUP SIGS2 ********************
C==============================================================================
C
C  FLIN1_Z
C   INTERPOLATE ON A 1D ARRAY
C
      SUBROUTINE FLIN1_Z(X,F,N,VTR,XMIN,XMAX,XLR,NX,IWARN)
c
c Interpolates the table VTR(NX) to find the value at x.
c If XLR>0, then the X grid is equally spaced on a logarithmic scale:
c
c   X(i) = XMIN * (exp(XLR)) ** ( (i-1)/(NX-1) )
c
c   So X(i) ranges from XMIN to XMIN*exp(XLR), and XLR=log(xmax/xmin).
c
c For XLR<=0, the X grid is equally spaced on a linear scale from
c XMIN to XMAX.
C
C Modified 7/26/00 (PCR) - vectorized.
C
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER n,nx,it,iwarn,ix,ixp1
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 f,vtr,x,xmin,xmax,xlr,zx,zy,zzx,zlogx
C============
      DIMENSION VTR(NX)
      DIMENSION X(N)
      DIMENSION F(N)
C
C  INTERPOLATE TO X  LOGARITHMIC SPACING IN X
C     LINEAR IF XLR .LE. 0.0
C
      integer idebug
      common/flint_debug/ idebug
C
      iwarn=0
C
      do 70 it=1,n
C
         ZX=X(it)
C
         if(zx .lt. xmin) then
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
         IX=ZZX
C
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
C
         IF(IX.LT.NX) GO TO 30
C
         IF(ZX.GT.(XMAX + 5.0D-7)) THEN
            IWARN=1
         ENDIF
         IX=NX-1
         ZZX=NX
C
C---------------------------
C  OK
C
 30      CONTINUE
C
         IXP1=IX+1
         ZZX=ZZX-IX
C
C  INTERPOLATE
C
         F(it)=(1.0D0-ZZX)*VTR(IX)+ZZX*VTR(IXP1)
c
         if(idebug.eq.99) then
            write(6,*) 'flin1_z:  zx=',zx,' ix=',ix
            write(6,*) '  vtr(ix:ixp1) = ',vtr(ix:ixp1)
            write(6,*) '  it,F(it)=',it,F(it)
         endif
C
 70   continue
C
      idebug=0  ! clear debug flag
C
      END
C
C******************** END FILE FLIN1.FOR ; GROUP SIGS2 ************************
! 26Jul2000 fgtok -s r8_precision.sub "r8con.csh conversion"
