!-----------------------------------------------------------------------
!     SUBROUTINE:     LIM_DIST
!
!     PURPOSE:        This subroutine calculates the minimum distance
!                     between a point (P) and a line segment defined by 
!                     two additional points (A and B).  If the cord of
!                     minimum distance falls between A and B then the
!                     minimum distance is returned, otherwise the
!                     value 1.0E+30 is returned.
!
!     INPUTS:         rp        Radial position of point P
!                     zp        Vertical position of point P
!                     r_lim(1)  Radial position of point A
!                     r_lim(2)  Radial position of point B
!                     z_lim(1)  Vertical position of point A
!                     z_lim(2)  Vertical position of point B
!
!     OUTPUTS:        d         Miminum distance between P and chord AB
!
!     LIBRARIES:      lib_opt.a - kind_spec
!
!     WRITTEN BY:     Neil Pomphrey (pomphrey@pppl.gov)
!
!     DATE:           01/11/11
!-----------------------------------------------------------------------
      subroutine lim_dist( rp, zp, r_lim, z_lim, d)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      use stel_kinds
!-----------------------------------------------------------------------
!     Input Arguments (see above)
!-----------------------------------------------------------------------
      implicit none
      real(rprec) :: rp, zp, d
      real(rprec), dimension(2) :: r_lim, z_lim
!-----------------------------------------------------------------------
!     Local Variables
!          dlr     Radial distance between A and B
!          dlz     Vertical distance between A and B
!          dl      Square of the length of segment AB (dlr^2+dlz^2)
!          u       Projection of AP onto AB (AP dot AB)/|AB|
!          dr      Radial point where the minimum distance lies on AB
!          dz      Vertical point where the mimimum distance lies on AB
!          d1      Square of distance from minimum point to A
!          d2      Square of distance from minimum point to B
!-----------------------------------------------------------------------
      real(rprec) :: dlr, dlz, dl, u, dr, dz, d1, d2
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------

      dlr = r_lim(2) - r_lim(1)
      dlz = z_lim(2) - z_lim(1)
      dl = (dlr**2 + dlz**2)

      u = ((rp - r_lim(1))*dlr + (zp - z_lim(1))*dlz)/dl

      dr = r_lim(1) + u*dlr
      dz = z_lim(1) + u*dlz

      d = (rp-dr)**2 + (zp-dz)**2

      d1 = (dr-r_lim(1))**2 + (dz-z_lim(1))**2
      d2 = (dr-r_lim(2))**2 + (dz-z_lim(2))**2
!     Check to see that the mimum point lies between A and B
      if( d1 > dl .or. d2 > dl) then
         d = 1.e30
      else
!         d = sign(sqrt(d), dlr*(zp-dz) - dlz*(rp-dr))
          d = sqrt(d)
      endif
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      end subroutine lim_dist
