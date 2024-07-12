!-----------------------------------------------------------------------
!     Function:      dflux_simpson
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/02/2012
!     Description:   This subroutine calculates the flux due to a line
!                    segment using simpson formula.
!-----------------------------------------------------------------------
      REAL FUNCTION dflux_simpson(x0,y0,z0,x1,y1,z1,cg)
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
      REAL(rprec) :: delf,delt,dx,dy,dz,xp,yp,zp,ax,ay,az
      REAL(rprec) :: xvec(3),avec(3)
      INTEGER     :: i,ier
      integer, parameter :: nop=3
      real(rprec), dimension(nop), parameter :: ci=(/1./6.,2./3.,1./6./)
!-----------------------------------------------------------------------
!     Begin Function
!-----------------------------------------------------------------------

      delf = 0
      delt = REAL(1)/REAL(nop-1)
      dx = x1-x0
      dy = y1-y0
      dz = z1-z0
      IF (cg > 0) THEN
         DO i = 1, nop
            xvec(1) = x0+(i-1)*dx*delt
            xvec(2) = y0+(i-1)*dy*delt
            xvec(3) = z0+(i-1)*dz*delt
            CALL bsc_a(coil_group(cg),xvec,avec)
            delf = delf+ ci(i)*(avec(1)*dx+avec(2)*dy+avec(3)*dz)
         END DO
      ELSE
         DO i = 1, nop
            ier = 1; ax = 0; ay = 0; az = 0
            xp = x0+(i-1)*dx*delt
            yp = y0+(i-1)*dy*delt
            zp = z0+(i-1)*dz*delt
            CALL vecpot_vc(xp,yp,zp,ax,ay,az,ier)
            delf = delf+ ci(i)*(ax*dx+ay*dy+az*dz)
         END DO
      END IF
      dflux_simpson = delf
      RETURN
!-----------------------------------------------------------------------
!     End Function
!-----------------------------------------------------------------------    
      END FUNCTION dflux_simpson
