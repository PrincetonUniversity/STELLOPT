

      function Map_v(v)
      use stel_kinds
      USE torlines_background, only: u0, v0, vb_ws
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) :: v, Map_v
C-----------------------------------------------

!
!     Finds phi (=vb_ws/Np) at a displaced winding surface corresponding to a
!     value of phi (=v/Np) on the plasma surface. It then computes the difference
!     compared with original phi (=v0/Np) on the plasma surface
!
      
      call normal_vector(u0, v)
c
      Map_v = vb_ws - v0                               !v0 is ORIGINAL Np*(Real Phi) on plasma surface

      end function Map_v
