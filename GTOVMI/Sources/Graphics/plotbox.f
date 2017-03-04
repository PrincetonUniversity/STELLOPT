        SUBROUTINE graf0 (y,n,lx,ly,lt,runlbl)
         IMPLICIT NONE
         INTEGER n, i
         REAL, DIMENSION(n), INTENT(in) ::  y
         REAL, DIMENSION(n) :: x
         REAL :: ymin, ymax, xmin, xmax, siz
         CHARACTER*(*) lx,ly,lt,runlbl
         CALL pgsave
         CALL pgbbuf
         CALL pgsci(1)
         xmin=1; xmax=n
         DO i=1,n 
          x(i)=i
         ENDDO
         ymin=MINVAL(y);ymax=MAXVAL(y)
         IF (ymax  .EQ.  ymin) THEN
          IF (ymax  .GT.  0.) ymin=.99*ymax
          IF (ymax  .LT.  0.) ymax=.99*ymin
          IF (ymax  .EQ.  0.) ymax=1.e-3;ymin=-1.e-3
         ENDIF
         CALL pgenv(MINVAL(x),MAXVAL(x),ymin,ymax,0,0)
         CALL pglab(TRIM(lx),TRIM(ly),TRIM(lt))
         CALL pgqch(siz)
         CALL pgsch(siz-.5)
         CALL pgmtxt('B',4.19,0.5,0.5,TRIM(runlbl))
         CALL pgsch(siz)
         CALL pgsci(2)
         CALL pgslw(3)
         CALL pgpt(n,x,y,2)
         CALL pgebuf
         CALL pgunsa
        END SUBROUTINE graf0
        SUBROUTINE graf1 (x,y,n,lx,ly,lt,runlbl)
         IMPLICIT NONE
         INTEGER n
         REAL, DIMENSION(n), INTENT(in) :: x, y
         REAL :: ymin, ymax, xmin, xmax, siz
         CHARACTER*(*) lx,ly,lt,runlbl
         CALL pgsave
         CALL pgbbuf
         CALL pgsci(1)
         xmin=MINVAL(x);xmax=MAXVAL(x)
         ymin=MINVAL(y);ymax=MAXVAL(y)
         IF (ymax  .EQ.  ymin) THEN
          IF (ymax  .GT.  0.) ymin=.99*ymax
          IF (ymax  .LT.  0.) ymax=.99*ymin
          IF (ymax  .EQ.  0.) ymax=1.e-3;ymin=-1.e-3
         ENDIF
         CALL pgenv(MINVAL(x),MAXVAL(x),ymin,ymax,0,0)
         CALL pglab(TRIM(lx),TRIM(ly),TRIM(lt))
         CALL pgqch(siz)
         CALL pgsch(siz-.5)
         CALL pgmtxt('B',4.19,0.5,0.5,TRIM(runlbl))
         CALL pgsch(siz)
         CALL pgsci(2)
         CALL pgslw(3)
         CALL pgline(n,x,y)
         CALL pgebuf
         CALL pgunsa
        END SUBROUTINE graf1
        SUBROUTINE graf1x (x,y,n,lx,ly,lt,runlbl)
         IMPLICIT NONE
         INTEGER n
         REAL, DIMENSION(n), INTENT(in) :: x, y
         REAL :: ymin, ymax, xmin, xmax, siz
         CHARACTER*(*) lx,ly,lt,runlbl
         CALL pgsave
         CALL pgbbuf
         CALL pgsci(1)
         xmin=MINVAL(x);xmax=MAXVAL(x)
         ymin=MINVAL(y);ymax=MAXVAL(y)
         IF (ymax  .EQ.  ymin) THEN
          IF (ymax  .GT.  0.) ymin=.99*ymax
          IF (ymax  .LT.  0.) ymax=.99*ymin
          IF (ymax  .EQ.  0.) ymax=1.e-3;ymin=-1.e-3
         ENDIF
         CALL pgenv(MINVAL(x),MAXVAL(x),ymin,ymax,0,1)
         CALL pglab(TRIM(lx),TRIM(ly),TRIM(lt))
         CALL pgqch(siz)
         CALL pgsch(siz-.5)
         CALL pgmtxt('B',4.19,0.5,0.5,TRIM(runlbl))
         CALL pgsch(siz)
         CALL pgsci(2)
         CALL pgslw(3)
         CALL pgline(n,x,y)
         CALL pgebuf
         CALL pgunsa
        END SUBROUTINE graf1x

        SUBROUTINE graf1pt (x,y,n,lx,ly,lt,runlbl)
         IMPLICIT NONE
         INTEGER n
         REAL, DIMENSION(n), INTENT(in)  :: x, y
         REAL :: ymin, ymax, siz
         CHARACTER*(*) lx,ly,lt,runlbl
         CALL pgsave
         CALL pgsci(1)
         ymin=MINVAL(y);ymax=MAXVAL(y)
         IF (ymax  .EQ.  ymin) THEN
          IF (ymax  .GT.  0.) ymin=.99*ymax
          IF (ymax  .LT.  0.) ymax=.99*ymin
          IF (ymax  .EQ.  0.) ymax=1.e-3;ymin=-1.e-3
         ENDIF
         CALL pgpage
         CALL pgvstd
         CALL pgswin(MINVAL(x),MAXVAL(x),ymin,ymax)
         CALL pgbox('BCNT',0.,0,'BCNTP1',0.,0)
         CALL pglab(TRIM(lx),TRIM(ly),TRIM(lt))
         CALL pgqch(siz)
         CALL pgsch(siz-.5)
         CALL pgmtxt('B',4.19,0.5,0.5,TRIM(runlbl))
         CALL pgsch(siz)
         CALL pgsci(2)
         CALL pgpt(n,x,y,2)
         CALL pgunsa
        END SUBROUTINE graf1pt
        SUBROUTINE graf2 (x,y,yfit,n,lx,ly,lt,runlbl)
!	two lines on a SINgle scale
         IMPLICIT NONE
         INTEGER i, n
         REAL :: ymin, ymax, siz
         REAL, DIMENSION(n) :: x,y,yfit
         CHARACTER*(*) lx,ly,lt,runlbl
         CALL pgsave
         CALL pgbbuf
         CALL pgsci(1)
         ymin=MIN(MINVAL(y),MINVAL(yfit))
         ymax=MAX(MAXVAL(y),MAXVAL(yfit))
         IF (ymax  .EQ.  ymin) THEN
          IF (ymax  .GT.  0.) ymin=.99*ymax
          IF (ymax  .LT.  0.) ymax=.99*ymin
          IF (ymax  .EQ.  0.) ymax=1.e-3;ymin=-1.e-3
         ENDIF
         CALL pgenv(MINVAL(x),MAXVAL(x),ymin,ymax,0,0)
         CALL pglab(TRIM(lx),TRIM(ly),TRIM(lt))
         CALL pgqch(siz)
         CALL pgsch(siz-.5)
         CALL pgmtxt('B',4.19,0.5,0.5,TRIM(runlbl))
         CALL pgsch(siz)
         CALL pgsci(3)
         CALL pgpt(n,x,y,-2)
         CALL pgslw(3)
         CALL pgsls(1)
         CALL pgsci(2)
         CALL pgline(n,x,yfit)
         CALL pgsci(3)
         CALL pgpt(n,x,y,2)
         CALL pgebuf
         CALL pgunsa
        END SUBROUTINE graf2

        SUBROUTINE graf2pt (x1,x2,y1,y2,n,lx,ly1,ly2,lt,runlbl)
         INTEGER n
         REAL, DIMENSION(n), INTENT(in)  :: x1, x2, y1, y2
         REAL :: ymin, ymax, siz
         CHARACTER*(*) lx,ly1,ly2,lt,runlbl
         CALL pgsave
         CALL pgslw(1)
         CALL pgsci(1)
         xmin=MINVAL(x1);xmax=MAXVAL(x1)
         xmin=MIN(xmin,MINVAL(x2))
         xmax=MAX(xmax,MAXVAL(x2))
         ymin=MINVAL(y1);ymax=MAXVAL(y1)
         ymin=MIN(ymin,MINVAL(y2))
         ymax=MAX(ymax,MAXVAL(y2))
         IF (ymax  .EQ.  ymin) THEN
          IF (ymax  .GT.  0.) ymin=.99*ymax
          IF (ymax  .LT.  0.) ymax=.99*ymin
          IF (ymax  .EQ.  0.) ymax=1.e-3;ymin=-1.e-3
         ENDIF
         CALL pgpage
         CALL pgvstd
         CALL pgswin(xmin,xmax,ymin,ymax)
         CALL pgbox('BCNT',0.,0,'BCNTP',0.,0)
         CALL pglab(TRIM(lx),TRIM(ly1//"   "//ly2),TRIM(lt))
         CALL pgqch(siz)
         CALL pgmtxt('T',.919,0.5,0.5,TRIM(runlbl))
         CALL pgsch(siz)
         CALL pgsci(4)
         CALL pgslw(1)
!         CALL pgsch(siz-.5)
         CALL pgpt(n,x1,y1,5)
         CALL pgsch(siz)
!        CALL pgmtxt('L',2.65,0.,0.,TRIM(ly1))
         CALL pgsci(2)
         CALL pgslw(1)
         CALL pgpt(n,x2,y2,4)
         CALL pgsch(siz)
!        CALL pgmtxt('L',2.65,1.0,1.0,TRIM(ly2))
         CALL pgunsa
        END SUBROUTINE graf2pt

        SUBROUTINE graf3pt (x,y1,y2,y3,n,lx,ly1,ly2,ly3,lt,runlbl)
         INTEGER n
         REAL, DIMENSION(n), INTENT(in)  :: x, y1, y2, y3
         REAL :: ymin, ymax, siz
         CHARACTER*(*) lx,ly1,ly2,ly3,lt,runlbl
         CALL pgsave
         CALL pgsci(1)
         ymin=MINVAL(y1);ymax=MAXVAL(y1)
         ymin=MIN(ymin,MINVAL(y2))
         ymax=MAX(ymax,MAXVAL(y2))
         ymin=MIN(ymin,MINVAL(y3))
         ymax=MAX(ymax,MAXVAL(y3))
         IF (ymax  .EQ.  ymin) THEN
          IF (ymax  .GT.  0.) ymin=.99*ymax
          IF (ymax  .LT.  0.) ymax=.99*ymin
          IF (ymax  .EQ.  0.) ymax=1.e-3;ymin=-1.e-3
         ENDIF
         CALL pgenv(MINVAL(x),MAXVAL(x),ymin,ymax,0,0)
         CALL pglab(TRIM(lx),' ',TRIM(lt))
         CALL pgqch(siz)
         CALL pgsch(siz-.5)
         CALL pgmtxt('B',4.19,0.5,0.5,TRIM(runlbl))
         CALL pgsch(siz)
         CALL pgsci(1)
         CALL pgpt(n,x,y1,2)
         CALL pgmtxt('L',2.65,0.,0.,TRIM(ly1))
         CALL pgsci(2)
         CALL pgpt(n,x,y2,3)
         CALL pgmtxt('L',2.65,0.5,0.5,TRIM(ly2))
         CALL pgsci(4)
         CALL pgpt(n,x,y3,4)
         CALL pgmtxt('L',2.65,1.,1.,TRIM(ly3))
         CALL pgunsa
        END SUBROUTINE graf3pt

        SUBROUTINE graf2x (x,y,yfit,n,lx,lyl,lyr,lt,runlbl)
!        two lines in a box with y axis left and right
         IMPLICIT NONE
         INTEGER i, n
         REAL :: ymin, ymax, dx, siz
         REAL, DIMENSION(n) :: x,y,yfit
         CHARACTER*(*) lx,lyl,lyr,lt,runlbl
         CHARACTER*8 xopt,yopt
         dx=x(3)-x(1)
         IF (n  .LT.  30)dx=x(2)-x(1)
         CALL pgsave
         CALL pgbbuf
         CALL pgsci(1)
         CALL pgpage
         CALL pgvstd
         CALL pgswin(MINVAL(x),MAXVAL(x),MINVAL(y),MAXVAL(y))
         xopt='BCNST';yopt='B'
         CALL pgbox(TRIM(xopt),0.,0,TRIM(yopt),0.,0)
         CALL pglab(TRIM(lx),' ',TRIM(lt))
         CALL pgqch(siz)
         CALL pgsch(siz-.5)
         CALL pgmtxt('B',4.19,0.5,0.5,TRIM(runlbl))
         CALL pgsch(siz)
         CALL pgsci(2)
         yopt='BNST'
         CALL pgbox(' ',0.,0,TRIM(yopt),0.,0)
         CALL pgmtxt('L',2.25,0.5,0.5,TRIM(lyl))
         CALL pgsls(1)
         CALL pgslw(3)
         CALL pgline(n,x,y)
         xopt=' ';yopt='CMST'
         CALL pgsci(4)
         CALL pgswin(MINVAL(x),MAXVAL(x),MINVAL(yfit),MAXVAL(yfit))
         CALL pgbox(' ',0.,0,TRIM(yopt),0.,0)
         CALL pgmtxt('R',2.25,0.5,0.5,TRIM(lyr))
         CALL pgsls(4)
         CALL pgslw(8)
         CALL pgline(n,x,yfit)
         CALL pgebuf
         CALL pgunsa
        END SUBROUTINE graf2x

