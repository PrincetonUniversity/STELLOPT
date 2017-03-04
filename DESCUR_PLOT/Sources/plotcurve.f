      program plotcurve
      use vplot
      integer :: istat
      real, dimension(:), allocatable :: rin, zin
      real, dimension(mnd) :: rbc, zbs, rbs, zbc

      open(unit=10,file = 'plotout',status = 'old', iostat=istat)
      if (istat .ne. 0) stop 'XDES_PLOT could not open plotout file'

      read(10,1990)mpol,ntheta,nphi,mpol1,nphi2,nfp,mpnt

      allocate (rin(ntheta*nphi), zin(ntheta*nphi), stat=istat)
      if (istat .ne. 0) stop 'XDES_PLOT allocation error'

      do i = 1,ntheta*nphi
         read(10,2000) rin(i),zin(i)
      enddo

 1990   format(7i10)
 2000   format(1p2e12.4)


      call plotter(rin,zin,rbc,zbs,rbs,zbc)

      deallocate (rin, zin)

      call system("ctrans -d X11 gmeta")

      stop

      end program plotcurve
