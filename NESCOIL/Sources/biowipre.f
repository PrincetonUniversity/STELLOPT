
! ----------------------------------------------------------------------
      subroutine biowipre
! ----------------------------------------------------------------------
c                                                             05/01/89
c     purpose:
c
c
c ----------------------------------------------------------------------
      use Vmeshes
      USE Vsurfac13
      USE Vculine1
      USE Vbiopre, ONLY: vx, vy, vz, dx, dy, dz
      implicit none
C-----------------------------------------------
      dx(:nw-1) = (xw(2:nw)-xw(:nw-1))*curre(:nw-1)
      dy(:nw-1) = (yw(2:nw)-yw(:nw-1))*curre(:nw-1)
      dz(:nw-1) = (zw(2:nw)-zw(:nw-1))*curre(:nw-1)
      vx(:nw-1) = yw(:nw-1)*dz(:nw-1) - zw(:nw-1)*dy(:nw-1)
      vy(:nw-1) = zw(:nw-1)*dx(:nw-1) - xw(:nw-1)*dz(:nw-1)
      vz(:nw-1) = xw(:nw-1)*dy(:nw-1) - yw(:nw-1)*dx(:nw-1)
c ----------------------------------------------------------------------

      end subroutine biowipre
