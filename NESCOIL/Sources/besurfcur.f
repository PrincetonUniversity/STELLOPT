
! ----------------------------------------------------------------------
      subroutine besurfcur(spx, spy, spz)
! ----------------------------------------------------------------------
c                                                             01/01/89
c     purpose:
c
c
c ----------------------------------------------------------------------
      use Vmeshes
      use NumParams
      USE Vprecal1
      USE Vprecal4
      USE Vsurface1, ONLY: x, y, z
      USE Vsurface7, ONLY: snx, sny, snz
      USE Vsurface8, ONLY: xcur, ycur, zcur
      USE Vdiagno3, ONLY: pote
      USE Vvacuum3
      USE Vcen, ONLY: xp, yp, zp, bmod, bmod1, bx, by, bz
      USE LoopCtrl
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(nuv) :: spx, spy, spz
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ico = 0, i, k
      real(rprec), allocatable, dimension(:) ::
     1     dbx, dby, dbz, dx, dy, dz, sq2, sq3, dn
      real(rprec) :: bxk, byk, bxadd, byadd, bzadd
c ----------------------------------------------------------------------
c
      i=0
      if (.not.allocated(dbx))
     1 allocate (dbx(nuv), dby(nuv), dbz(nuv), dn(nuv), dx(nuv),
     1   dy(nuv), dz(nuv), sq2(nuv), sq3(nuv), stat=i)
      if (i .ne. 0) stop 'allocation error in NESCOIL besurfcur'
      dbx(:) = zero
      dby(:) = zero
      dbz(:) = zero
      dn(:) = zero
      dx(:) = zero
      dy(:) = zero
      dz(:) = zero
      sq2(:) = zero
      sq3(:) = zero

      if(ico .eq. 0) then
         ico    = 1
         spx = pote*snx
         spy = pote*sny
         spz = pote*snz
      endif

      bx     = zero
      by     = zero
      bz     = zero

      do k = 0,np-1
         dx  = xp*cok(k) + yp*sik(k) - x
         dy  = yp*cok(k) - xp*sik(k) - y
         dz  = zp  - z
         sq2 = one/(dx*dx + dy*dy + dz*dz)
         sq3 = sq2*sqrt(sq2)
         dn  = three*sq2*(spx*dx + spy*dy + spz*dz)
         dbx = (spx - dn*dx + ycur*dz - zcur*dy)*sq3
         dby = (spy - dn*dy + zcur*dx - xcur*dz)*sq3
         dbz = (spz - dn*dz + xcur*dy - ycur*dx)*sq3
c
         bxk    = sum(dbx(1:nuv))*fnuv
         byk    = sum(dby(1:nuv))*fnuv
         bx     = bx + (bxk*cok(k) - byk*sik(k))
         by     = by + (byk*cok(k) + bxk*sik(k))
         bz     = bz + sum(dbz(1:nuv))*fnuv
      end do
C
      bx = -bx
      by = -by
      bz = -bz
C
c     TURNED THIS OFF, ADD BEN(i) IN CALLING ROUTINE
      if (ibex .ne. 0) then
         call bexter (xp, yp, zp, bxadd, byadd, bzadd)
         bx = bx + bxadd
         by = by + byadd
         bz = bz + bzadd
      endif
c
      bmod   = sqrt(bx*bx+by*by+bz*bz)
      bmod1  = one/bmod

      if (iloop .le. 0)
     1 deallocate (dbx, dby, dbz, dx, dy, dz, sq2, sq3, dn)

      end subroutine besurfcur
