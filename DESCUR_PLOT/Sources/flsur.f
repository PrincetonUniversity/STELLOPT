      subroutine flsur(rin,zin,r,z,rmax,rmin,zmax,zmin,kphi,ncurve1)
      use vplot
      real :: rin(ntheta,*), zin(ntheta,*), r(ncurve1), z(ncurve1)
      real :: x(ntheta), y(ntheta)
      character*40 :: mtitle

      pi = 4*atan(1.0)
      xav = rmin+rmax
      yav = zmin+zmax
      grsize = 5.
      dx = rmax-rmin
      dy = zmax-zmin
      step = max(dx,dy)/(.9*grsize)
      del = step*grsize
      xor = .5*(xav-del)
      yor = .5*(yav-del)
      xmax = .5*(xav+del)
      ymax = .5*(yav+del)
      deg = 360.*(kphi-1.)/real(nfp*nphi)
      write (mtitle,10)deg
 10   format('PHI = ',f6.2,' deg$')
      
!DEC$ IF .NOT.DEFINED (WIN32)
      call agsetr('X/MINIMUM.',xor)
      call agsetr('X/MAXIMUM.',xmax)
      call agsetr('X/NICE.',0.0)
      call agsetr('Y/MINIMUM.',yor)
      call agsetr('Y/MAXIMUM.',ymax)
      call agsetr('Y/NICE.',0.0)
      call agseti('FRAME.',2)
      call anotat('R$','Z$',0,0,0,0)
        call ezxy (r,z,ncurve1,mtitle)
        do 50 j = 1,ntheta
        x(j) = rin(j,kphi)
 50     y(j) = zin(j,kphi)
      call gsmksc(1.3)
        call points(x,y,ntheta,-4,0)
        call frame
!DEC$ ENDIF
        end subroutine flsur
