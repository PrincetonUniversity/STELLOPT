!-----------------------------------------------------------------------
!     Function:      db_midpoint
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/02/2012
!     Description:   This subroutine calculates the flux due to a line
!                    segment using midpoint formula.
!-----------------------------------------------------------------------
      REAL FUNCTION db_midpoint(x0,y0,z0,x1,y1,z1,cg,ldb)
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
      LOGICAL, intent(in), optional :: ldb
      LOGICAL     :: ldb_local
      INTEGER     :: ier
      real(rprec) :: dx,dy,dz,xp,yp,zp,bx,by,bz
      REAL(rprec) :: xvec(3),bvec(3)
!-----------------------------------------------------------------------
!     Begin Function
!-----------------------------------------------------------------------
      dx = x1-x0
      dy = y1-y0
      dz = z1-z0
      xp = .5*(x0+x1)
      yp = .5*(y0+y1)
      zp = .5*(z0+z1)
      dl = 1.0
      ldb_local=.false.
      IF (PRESENT(ldb)) ldb_local = ldb
      IF (cg > 0) THEN
         xvec(1)=xp
         xvec(2)=yp
         xvec(3)=zp
         CALL bsc_b(coil_group(cg),xvec,bvec)
         dflux_midpoint = bvec(1)*dx+bvec(2)*dy+bvec(3)*dz
      ELSE
         ier = 1
         CALL bfield_vc(xp,yp,zp,bx,by,bz,ier)
         !CALL bfield_virtual_casing_adapt(xp,yp,zp,bx,by,bz,ier)
         !CALL bfield_volint_adapt(xp,yp,zp,bx,by,bz,ier)
         db_midpoint = bx*dx+by*dy+bz*dz
      END IF
      IF (ldb_local) db_midpoint = db_midpoint/SQRT(dx*dx+dy*dy+dz*dz) 
      RETURN
!-----------------------------------------------------------------------
!     End Function
!-----------------------------------------------------------------------    
      END FUNCTION db_midpoint
