      subroutine plotter(rin,zin,rbc,zbs,rbs,zbc)
      use vplot
      integer, parameter :: ncurve=99, ncurve1=ncurve+1
      integer :: lmodes, l
      real :: bound
      real, dimension(0:mpol-1, -nphi2:nphi2) :: rbc, zbs, rbs, zbc
      real, dimension(ncurve1) :: r, z
      real, dimension(ntheta*nphi) :: rin, zin
      real :: xc(4*mnd)

!DEC$ IF .NOT.DEFINED (WIN32)
        xc = 0

        call opngks
        print *,' PLOTTING SURFACES..... '

        bound = 1.e-10
        l = 0
        m1 = 0;    n1 = 0
        do 10 m = 0,mpol-1
        do 10 n = -nphi2,nphi2
        read(10,2000)rbc(m,n),zbs(m,n),rbs(m,n),zbc(m,n)
        if( abs(rbc(m,n)).lt.bound .and. abs(zbs(m,n)).lt.bound .and.
     1      abs(rbs(m,n)).lt.bound .and. abs(zbc(m,n)).lt.bound )
     2      goto 10
        l = l+1
        m1(l) = m
        n1(l) = n*nfp
        xc(l) = rbc(m,n)
        xc(l+mpnt) = rbs(m,n)
        xc(l+2*mpnt) = zbc(m,n)
        xc(l+3*mpnt) = zbs(m,n)
 10     continue
 2000   format(1p4e12.4)
        close(unit = 10)
        lmodes = l
c
c       Compute "Bubble-Plot" of Mode-Amplitudes
c
        call bubble(mpnt,lmodes,m1,n1,xc,'R$')
        call bubble(mpnt,lmodes,m1,n1,xc(1+2*mpnt),'Z$')
c
c       Compute flux surfaces
        twopi = 8*atan(1.0)
        pit = twopi/ncurve
        do 20 kphi = 1,nphi
        call totz(xc,xc(1+mpnt),xc(1+2*mpnt),xc(1+3*mpnt),
     >  pit,kphi,r,z,ncurve,lmodes)
        if(kphi.ne.1) GOTO 30
        rmax = r(1)
        rmin = r(1)
        zmax = z(1)
        zmin = z(1)
 30     continue
        do 40 kt = 1,ncurve
        rmax = max(rmax,r(kt))
        rmin = min(rmin,r(kt))
        zmax = max(zmax,z(kt))
        zmin = min(zmin,z(kt))
 40     continue
 20     continue
        do 50 kphi = 1,nphi
        call totz(xc,xc(1+mpnt),xc(1+2*mpnt),xc(1+3*mpnt),
     >  pit,kphi,r,z,ncurve,lmodes)
        call flsur(rin,zin,r,z,rmax,rmin,zmax,zmin,kphi,ncurve1)
 50     continue
        call clsgks
cq      call system('ctrans -d ps.mono gmeta | lpr -Plr50')

!DEC$ ENDIF
        end subroutine plotter
