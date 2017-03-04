      SUBROUTINE packk(gp,xp,yp,mp,delta)
c-----------------------------------------------------------------------
c the INTerpolation and extrapolation routines are inaccurate IF some of
c points found by furplm are too CLOSE together.  this routine chooses
c subset of those points which are well spaced and redefines the arrays
c gr and drp to apply ONLY to these points.
c last changes made july 2 by fjh
c-----------------------------------------------------------------------
      USE precision
      IMPLICIT NONE
      REAL(rprec) :: xp(*),yp(*),gp(*)
      REAL(rprec) :: xs, ys, d, dd, delta
      INTEGER :: mp1, mpi, i, mp, i1, j
      mp1=mp-1
      mpi=mp1
      dd=delta*delta
      DO 1 i=1,mp1
      IF (i.gt.mpi) GOTO 3
6     xs=(xp(i+1)-xp(i))**2
      ys=(yp(i+1)-yp(i))**2
      d=xs+ys
      IF (d.ge.dd) GOTO 1
      IF (i.eq.mpi) GOTO 4
      i1=i+1
      DO  j=i1,mpi
         gp(j)=gp(j+1)
         xp(j)=xp(j+1)
         yp(j)=yp(j+1)
      END DO
      mpi=mpi-1
      GOTO 6
1     CONTINUE
3     mp=mpi+1
      GOTO 7
4     xp(mpi)=xp(mpi+1)
      yp(mpi)=yp(mpi+1)
      gp(mpi)=gp(mpi+1)
      mp=mpi
    7 DO 10 i=1,mp
      gp(i)=ABS(gp(i))
   10 CONTINUE
      END SUBROUTINE packk
