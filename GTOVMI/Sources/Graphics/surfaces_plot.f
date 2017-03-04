      subroutine surfaces_plot(g1,g2,device)
      use precision
      use Vname0
      use gfile             
      use mapout
      implicit none
      type(geqdsk) :: g1  
      type(rszs) :: g2
      character*(*) :: device
      real*4, dimension(:), allocatable :: sqphin,thetap,xp,yp
      real, pointer :: xpnt(:,:),ypnt(:,:)
      real(rprec) :: pi
      real*4 xmin, xmax, zmin, zmax
      integer :: id0, id1, id, pgopen
      integer :: jpsi,  itht, j, ij, m
      character(len=60) :: mlabel
      interface
       subroutine evals_boundary(g2,xp,yp,icount)
        use Vname0
        use mapout
        real(rprec), dimension(0:mu-1,0:0)  :: rbc,zbs,rbs,zbc 
        real*4, dimension(:), allocatable  :: rb, zb
        real*4, dimension(*) :: xp, yp
        type(rszs) :: g2
       end subroutine evals_boundary
      end interface
      pi = 4 * atan(1._dbl)
      call pgask(.true.)
      zmin=-165.        !       allow for D3D coils
      zmax=165.
      xmin=80.
      xmax=270.
!      call  pgpap (11., (zmax-zmin)/(xmax-xmin)/2)     ! 
!      if(INDEX(TRIM(device),'p')/=0)
     	call  pgpap (8., (zmax-zmin)/(xmax-xmin)/2)
      call pgsubp(2,1)
      call  pgsch (2.0)                        ! set height of text
      call pgsci(1.)
      call  pgenv (xmin, xmax, zmin, zmax,1,0)      ! define the subplot area
      call pglab('R(cm)','Z(cm)','mapcode')
      call  pgscf (1  )                             ! select simple typeface
      call pgslw(5)
      m=g1%limitr;allocate(xp(m),yp(m))
      xp(1:m)=100*g1%xlim(1:m);yp(1:m)=100*g1%ylim(1:m)
      call pgline(m,xp,yp)
      if( allocated(xp) )deallocate(xp,yp)
      allocate(xp(nu),yp(nu))
      xp=0;yp=0
      call evals_boundary(g2,xp,yp,m)
      xp=xp*100;yp=yp*100
      call pgpt(m,xp,yp,5)
      if( allocated(xp) )deallocate(xp,yp)
      call pltcol
      call pgslw(2)
      call pgsci(3)
      m=g1%nbdry;allocate(xp(m),yp(m))
      xp(1:m)=100*g1%rbdry(1:m);yp(1:m)=100*g1%zbdry(1:m)
      call pgpt(m,xp,yp,3)
      if( allocated(xp) )deallocate(xp,yp)
      call pgslw(2)
      itht=g2%nthet
      jpsi=g2%npsi
      allocate(xp(g2%nthet),yp(g2%nthet))
      call pgsci(4)
      do j=1,jpsi,32
         xp(1:itht)=100*g2%rs(j,1:itht)
         yp(1:itht)=100*g2%zs(j,1:itht)
         call pgline(itht,xp,yp)
      enddo
      if( allocated(xp) )deallocate(xp,yp)
      allocate(xp(g2%npsi),yp(g2%npsi))
      call pgsci(5)
      ij=itht/16
      do j=1,itht,ij
         xp(1:jpsi)=100*g2%rs(1:jpsi,j)
         yp(1:jpsi)=100*g2%zs(1:jpsi,j)
         call pgline(jpsi,xp,yp)
      enddo
      call pgsci(2)
      j=1
      xp(1:jpsi)=100*g2%rs(1:jpsi,j)
      yp(1:jpsi)=100*g2%zs(1:jpsi,j)
      call pgline(jpsi,xp,yp)
      j=itht/2
      xp(1:jpsi)=100*g2%rs(1:jpsi,j)
      yp(1:jpsi)=100*g2%zs(1:jpsi,j)
      call pgline(jpsi,xp,yp)
      call pgsci(1)
c second plot
      call  pgenv (xmin, xmax, zmin, zmax,1,0)      ! define the subplot area
      call pglab('R(cm)','Z(cm)','DESCUR v mapped eqdsk')
      write(mlabel,fmt='("% boundary flux=",f8.4)')
     .		100.*g2%fraction_bndry
      CALL pgmtxt('B',2.25,0.5,0.5,TRIM(mlabel))
      write(mlabel,fmt='("mpol=",i3.3)')mu
      CALL pgmtxt('T',0.325,0.5,0.5,TRIM(mlabel))
      call  pgscf (1  )                             ! select simple typeface
      call  pgsch (2.0)                        ! set height of text
      call pgslw(5)
      if( allocated(xp) )deallocate(xp,yp)
      m=g1%limitr;allocate(xp(m),yp(m))
      xp(1:m)=100*g1%xlim(1:m);yp(1:m)=100*g1%ylim(1:m)
      call pgline(m,xp,yp)
      if( allocated(xp) )deallocate(xp,yp)
      call pltcol
      call pgslw(2)
      if( allocated(xp) )deallocate(xp,yp)
      allocate(xp(g2%nthet),yp(g2%nthet))
      call pgsci(4)
       j=g2%npsi
         xp(1:itht)=100*g2%rs(j,1:itht)
         yp(1:itht)=100*g2%zs(j,1:itht)
      call pgpt(itht,xp,yp,5)
      if( allocated(xp) )deallocate(xp,yp)
      call pgsci(3)
      if( allocated(xp) )deallocate(xp,yp)
      allocate(xp(nu),yp(nu))
      xp=0;yp=0
      call evals_boundary(g2,xp,yp,m)
      xp=xp*100;yp=yp*100
      call pgslw(5)
      call pgsci(2)
      call pgline(m,xp,yp)
      end subroutine surfaces_plot
      subroutine pltcol
      parameter (ncoil=18)
      dimension rc(ncoil), zc(ncoil), wc(ncoil),
     .hc(ncoil),xx(5),yy(5)
      dimension ac(ncoil),ac2(ncoil)
c
      data rc/
     &.8608,.8614,.8628,.8611,1.0041,2.6124,
     &2.3733,1.2518,1.6890,.8608,.8607,.8611,
     &.8630,1.0025,2.6124,2.3834,1.2524,1.6889/
      data zc/
     &.16830,.50810,.84910,1.1899,1.5169,0.4376,
     &1.1171,1.6019,1.5874,-.17370,-.51350,-.85430,
     &-1.1957,-1.5169,-0.4376,-1.1171,-1.6027,-1.5780/
      data wc/
     &.0508,.0508,.0508,.0508,.13920,0.17320,
     &0.1880,.23490,.16940,.0508,.0508,.0508,
     &.0508,.13920,0.17320,0.1880,.23490,.16940/
      data hc/
     &.32106,.32106,.32106,.32106,.11940,0.1946,
     &0.16920,.08510,.13310,.32106,.32106,.32106,
     &.32106,.11940,0.1946,0.16920,.08510,.13310/
      data ac/
     &0.,0.,0.,0.,45.0,0.,
     &0.,0.,0.,0.,0.,0.,
     &0.,-45.0,0.,0.,0.,0./
      data ac2/
     &0.,0.,0.,0.,0.,92.40,
     &108.06,0.,0.,0.,0.,0.,
     &0.,0.,-92.40,-108.06,0.,0./
      do i=1,ncoil
         x= rc(i)
         y= zc(i)
         dx=wc(i)/2.
         dy=hc(i)/2.
c         sn1=sind(ac(i))
         sn1=sin(pi*ac(i)/180)
c         cos1=cosd(ac(i))
         cos1=cos(pi*ac(i)/180)
         tac=sn1/cos1
c         sn2=sind(ac2(i))
         sn2=sin(pi*ac2(i)/180)
c         cos2=cosd(ac2(i))
         cos2=cos(pi*ac2(i)/180)
         tac2=sn2/cos2
         if(ac2(i).eq.0) go to 40
         xx(1)=x-dx-dy/tac2
         xx(3)=x+dx+dy/tac2
         xx(5)=xx(1)
         xx(2)=xx(3)-wc(i)
         xx(4)=xx(1)+wc(i)
c
         yy(1)=y-dy
         yy(2)=y+dy
         yy(3)=y+dy
         yy(4)=y-dy
         yy(5)=y-dy
         go to 50
   40    continue
c
         dac=0.5*wc(i)*tac
         xx(1)=x-dx
         xx(2)=x-dx
         xx(3)=x+dx
         xx(4)=x+dx
         xx(5)=x-dx
c
         yy(1)=y-dy-dac
         yy(2)=y+dy-dac
         yy(3)=y+dy+dac
         yy(4)=y-dy+dac
         yy(5)=y-dy-dac
   50    continue
         do is=1,5
           xx(is)=xx(is)*100
           yy(is)=yy(is)*100
         enddo
         call pgsci(1)
         call pgslw(3)
         call pgline(5,xx,yy)
         call pgsci(2)
         call pgpoly(5,xx,yy)
      enddo
      end subroutine pltcol
       subroutine arc_integral(rv,zv,jpsi,itht,arcl)
       real arcl(jpsi-1)
       real, dimension(jpsi,itht) :: rv,zv
       real, dimension(:), allocatable :: dl
       if(.not.allocated(dl))allocate(dl(itht-1))
       pi=4 * atan(1.0)
       itht=itht
       do k=1,jpsi
        dl=sqrt( (rv(k,2:itht)-rv(k,1:itht-1))**2 + 
     &	      (zv(k,2:itht)-zv(k,1:itht-1))**2 	)
        arcl(k)=Sum(dl)
       enddo
       deallocate(dl)
       end subroutine arc_integral
