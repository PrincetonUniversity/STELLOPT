
! ----------------------------------------------------------------------
      subroutine matrix(a, ip)
! ----------------------------------------------------------------------
c                                                             01.01.89
c     purpose:
c
c
C-----------------------------------------------
      use Vmeshes
      use NumParams
      USE Vprecal1
      USE Vprecal4
      USE Vsurface1
      USE Vsurface7, ONLY: snx, sny, snz
      USE Vsurface8, ONLY: xcur, ycur, zcur
      USE Vsurface9
      USE Vsurfac14, ONLY: bex, bey, bez
      USE LoopCtrl
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ip
      real(rprec), dimension(3*nuv) :: a
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: in, in1, in2, k, i
      real(rprec) :: sux, suy, suz
      real(rprec), allocatable, dimension(:) ::
     1     dax, day, daz, dbx, dby, dbz, dx, dy, dz, sq2, sq3, pn
c ----------------------------------------------------------------------

      i=0
      allocate (dax(nuv), day(nuv), daz(nuv), dbx(nuv),
     1          dby(nuv), dbz(nuv), dx(nuv), dy(nuv), dz(nuv),
     2          pn(nuv), sq2(nuv), sq3(nuv), stat = i)
      if (i .ne. 0) stop 'Allocation error in nescoil matrix'
      dax(:) = zero
      day(:) = zero
      daz(:) = zero
      dbx(:) = zero
      dby(:) = zero
      dbz(:) = zero
      dx(:) = zero
      dy(:) = zero
      dz(:) = zero
      pn(:) = zero
      sq2(:) = zero
      sq3(:) = zero

      in  = 0
      in1 = in + nuv
      in2 = in1 + nuv
      dz = z  - z1(ip)
      do k = 0, np - 1
         dx = x  - x1(ip)*cok(k) + y1(ip)*sik(k)
         dy = y  - y1(ip)*cok(k) - x1(ip)*sik(k)
         sq2 = one/(dx*dx + dy*dy + dz*dz)
         sq3 = sq2*sqrt(sq2)
         pn = three*sq2*(dx*snx + dy*sny + dz*snz )
         dax  = (snx  - pn*dx)*sq3
         day  = (sny  - pn*dy)*sq3
         daz  = (snz  - pn*dz)*sq3

         a(in +1:in1) = a(in +1:in1) + cok(k)*dax  + sik(k)*day
         a(in1+1:in2) = a(in1+1:in2) - sik(k)*dax  + cok(k)*day
         a(in2+1:in2+nuv) = a(in2+1:in2+nuv) + daz

         dbx  = (ycur*dz - zcur*dy)*sq3
         sux = sum(dbx(:nuv))*fnuv

         dby  = (zcur*dx - xcur*dz)*sq3
         suy = sum(dby(:nuv))*fnuv

         dbz  = (xcur*dy - ycur*dx)*sq3
         suz = sum(dbz(:nuv))*fnuv

         bex(ip) = bex(ip) + cok(k)*sux + sik(k)*suy
         bey(ip) = bey(ip) - sik(k)*sux + cok(k)*suy
         bez(ip) = bez(ip) + suz
      end do

c........Deallocate matrix arrays
      deallocate (dax, day, daz, dbx, dby, dbz, pn,
     1            sq2, sq3, dx, dy, dz)

      end subroutine matrix
