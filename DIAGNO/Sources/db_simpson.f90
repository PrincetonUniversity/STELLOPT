!-----------------------------------------------------------------------
!     Function:      db_simpson
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/02/2012
!     Description:   This subroutine calculates the field due to a line
!                    segment using simpson formula.
!-----------------------------------------------------------------------
      REAL FUNCTION db_simpson(x0,y0,z0,x1,y1,z1,cg,ldb)
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
      IMPLICIT NONE
      INTEGER, intent(in) :: cg
      REAL(rprec)   , intent(in) :: x0,y0,z0,x1,y1,z1
      LOGICAL, intent(in), optional :: ldb
      LOGICAL     :: ldb_local
      REAL(rprec) :: delf,delt,dx,dy,dz,xp,yp,zp,bx,by,bz,dl
      REAL(rprec) :: xvec(3),bvec(3)
      INTEGER     :: i,ier
      integer, parameter :: nop=3
      real(rprec), dimension(nop), parameter :: ci=(/1./6.,2./3.,1./6./)
!-----------------------------------------------------------------------
!     Begin Function
!-----------------------------------------------------------------------

      delf = 0.0
      delt = 1./float(nop-1)
      dx = x1-x0
      dy = y1-y0
      dz = z1-z0
      dl = 1.0
      ldb_local=.false.
      IF (cg > 0) THEN
         DO i = 1, nop
            xvec(1) = x0+(i-1)*dx*delt
            xvec(2) = y0+(i-1)*dy*delt
            xvec(3) = z0+(i-1)*dz*delt
            CALL bsc_b(coil_group(cg),xvec,bvec)
            delf = delf+ ci(i)*(bvec(1)*dx+bvec(2)*dy+bvec(3)*dz)
         END DO
      ELSE
         DO i = 1, nop
            ier = 1; bx = 0.0; by = 0.0; bz = 0.0
            xp = x0+(i-1)*dx*delt
            yp = y0+(i-1)*dy*delt
            zp = z0+(i-1)*dz*delt
            CALL bfield_vc(xp,yp,zp,bx,by,bz,ier)
            !CALL bfield_virtual_casing_adapt(xp,yp,zp,bx,by,bz,ier)
            !CALL bfield_volint_adapt(xp,yp,zp,bx,by,bz,ier)
            delf = delf+ ci(i)*(bx*dx+by*dy+bz*dz)
         END DO
      END IF
      IF (ldb_local) THEN
         dl = 0.0
         DO i = 1, nop
            xp = x0+(i-1)*dx*delt
            yp = y0+(i-1)*dy*delt
            zp = z0+(i-1)*dz*delt
            dl = dl+ ci(i)*SQRT(dx*dx+dy*dy+dz*dz)
         END DO
      END IF
      db_simpson = delf/dl
      RETURN
!-----------------------------------------------------------------------
!     End Function
!-----------------------------------------------------------------------    
      END FUNCTION db_simpson
