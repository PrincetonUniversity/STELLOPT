
!-------------------------------------------------------------
      subroutine bexter(xp, yp, zp, bx, by, bz)

!................................................................
c     Purpose:
c     Calculate B field from background coils at xp,yp,zp
c................................................................

      USE Vmeshes
      USE NumParams
      USE Vprecal1, ONLY: np
      USE Vvacuum3, ONLY: cup
      implicit none
c................................................................
C   D u m m y   A r g u m e n t s
c................................................................
      real(rprec) ::  xp, yp, zp, bx, by, bz
c................................................................
C   L o c a l   V a r i a b l e s
c................................................................
      real(rprec) :: bxy
c................................................................
c
c     FOR NOW, JUST SIMULATE EXTERNAL 1/R FIELD

      bxy =  cup / (xp*xp + yp*yp)
      bx  =  yp * bxy
      by  = -xp * bxy
      bz  =  zero

      end subroutine bexter
