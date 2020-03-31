!-----------------------------------------------------------------------
!     Function:      dflux_midpoint
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/02/2012
!     Description:   This subroutine calculates the flux due to a line
!                    segment using midpoint formula.
!-----------------------------------------------------------------------
      REAL FUNCTION dflux_midpoint(x0,y0,z0,x1,y1,z1,cg)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE virtual_casing_mod
      USE biotsavart
      USE bsc_T
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      INTEGER, intent(in) :: CG
      real(rprec)   , intent(in) :: x0,y0,z0,x1,y1,z1
      INTEGER     :: ier
      real(rprec) :: dx,dy,dz,xp,yp,zp,ax,ay,az
      REAL(rprec) :: xvec(3),avec(3)
!-----------------------------------------------------------------------
!     Begin Function
!-----------------------------------------------------------------------
      dx = x1-x0
      dy = y1-y0
      dz = z1-z0
      xp = .5*(x0+x1)
      yp = .5*(y0+y1)
      zp = .5*(z0+z1)
      IF (cg > 0) THEN
         xvec(1)=xp
         xvec(2)=yp
         xvec(3)=zp
         CALL bsc_a(coil_group(cg),xvec,avec)
         dflux_midpoint = avec(1)*dx+avec(2)*dy+avec(3)*dz
      ELSE
         ier = 1
         !CALL vecpot_virtual_casing_adapt(xp,yp,zp,ax,ay,az,ier)
         !CALL vecpot_volint_adapt(xp,yp,zp,ax,ay,az,ier)
         CALL vecpot_vc(xp,yp,zp,ax,ay,az,ier)
         dflux_midpoint = ax*dx+ay*dy+az*dz
      END IF
      RETURN
!-----------------------------------------------------------------------
!     End Function
!-----------------------------------------------------------------------    
      END FUNCTION dflux_midpoint
