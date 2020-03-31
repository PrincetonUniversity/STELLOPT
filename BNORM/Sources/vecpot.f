
      subroutine vecpot (ax, ay, az, i, nlo, nhi)
c ----------------------------------------------------------------------
      use stel_kinds
      use bnvariables, only: indu, indv, guu, guv, gvv, djx, djy, djz,
     1    x, y, z, tanu, tanv
      implicit none
      integer, intent(in) :: i, nlo, nhi
      integer :: ip, istat
      real(rprec), intent(inout) :: ax, ay, az
      real(rprec), target, allocatable :: dx(:), dy(:), dz(:)
      real(rprec), pointer :: sq(:), sqs(:), du(:), dv(:)
      real(rprec) :: sqsum
c ----------------------------------------------------------------------
      ip = nhi - nlo + 1
      if (ip .le. 0) return

      allocate (dx(nlo:nhi), dy(nlo:nhi), dz(nlo:nhi), stat=istat)
      if (istat .ne. 0) stop 'Allocation error in BNORM routine vecpot'

      sq => dz
      du => dx
      dv => dy
      sqs => dy

      dx    = x(i) - x(nlo:nhi)
      dy    = y(i) - y(nlo:nhi)
      dz    = z(i) - z(nlo:nhi)
      sq    = 1._dp/sqrt(dx*dx + dy*dy + dz*dz)
      du    = tanu(indu(nlo:nhi) - indu(i))
      dv    = tanv(indv(nlo:nhi) - indv(i))
      sqs   = 1._dp/sqrt(guu(i)*du*du + 2._dp*guv(i)*du*dv
     1                 + gvv(i)*dv*dv)

      sqsum = sum(sqs)

      ax = ax + sum(djx(nlo:nhi)*sq) - djx(i)*sqsum
      ay = ay + sum(djy(nlo:nhi)*sq) - djy(i)*sqsum
      az = az + sum(djz(nlo:nhi)*sq) - djz(i)*sqsum

      deallocate (dx, dy, dz)

      end subroutine vecpot
