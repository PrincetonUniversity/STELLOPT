      subroutine rdeqdsk(g1)
      USE precision
      USE gfile
      USE mapg_mod
c-----------------------------------------------------------------------
c     changed dpsi to dpsiv: rlm 7/3/96
c-----------------------------------------------------------------------
      IMPLICIT NONE
      TYPE(geqdsk) :: g1
      CHARACTER*6 ntitle(5)
      CHARACTER*6 dat
      INTEGER :: ipestg, limitr, mbdry
      REAL(rprec) :: xdim,zdim,rc,redge,zmid
      REAL(rprec) :: btor
      REAL(rprec) :: totcur,psimx(2),xax(2)
      REAL(rprec) :: zax(2),psisep,xsep,zsep
      REAL(rprec) :: dpsiv
      INTEGER :: i,j,isave
      neqdsk=118
c
      open (neqdsk,file=filename,status='old', iostat=isave)
      if(isave /= 0) then
        print*,"iostat=",isave," opening file=",trim(filename),
     .		" on unit=",neqdsk
        STOP 'EQDSK'
      endif
c     read eqdsk file
      write(6,'("begin read eqdsk")')      
      READ(neqdsk,200,err=211)(ntitle(i),i=1,5),dat,ipestg,nx,nz
      GOTO 212
 211  print*,' READ ERROR :ipestg,nx,nz=',ipestg,nx,nz
      STOP 'eqdsk 1st line' 
 212  continue
 200  format(6a8,3i4)
      call eqdsk_allocate(nx,nz)
      read(neqdsk,300)xdim,zdim,rc,redge,zmid
      read(neqdsk,300)xaxis,zaxis,psiaxis,psilim,btor
      read(neqdsk,300)totcur,psimx(1),psimx(2),xax(1),xax(2)
      read(neqdsk,300)zax(1),zax(2),psisep,xsep,zsep
      read(neqdsk,300)(sf(i),i=1,nx)
      read(neqdsk,300)(sp(i),i=1,nx)
      read(neqdsk,300)(sffp(i),i=1,nx)
      read(neqdsk,300)(spp(i),i=1,nx)
      read(neqdsk,300)((psixz(i,j),i=1,nx),j=1,nz)
 300  format(5e16.9)
      kvtor=0
      nmass=0
      read(neqdsk,300,end=500) (qpsi(i),i=1,nx)
      read(neqdsk,'(2i5)') nbndry,nlim
      limitr=nlim; mbdry=nbndry
      allocate(xbndry(nbndry),zbndry(nbndry))
      read(neqdsk,300) (xbndry(i),zbndry(i),i=1,nbndry)
      allocate(xlim(nlim),zlim(nlim))
      read(neqdsk,300) (xlim(i),zlim(i),i=1,nlim)
      write(6,'("finished eqdsk")')
      if(rotate.ne.0) then
         read(neqdsk,'(i5,e16.9,i5)',end=500,err=500) kvtor,rvtor,nmass
         if(kvtor.gt.0) then
            read(neqdsk,300) (pressw(i),i=1,nx)
            read(neqdsk,300) (pwprim(i),i=1,nx)
c     guard against negative or zero pressw
            isave=1
            do i=1,nx
               if(pressw(i).le.0.) then
                  isave=i
                  go to 75
               endif
            end do
 75         continue
            if(isave.ne.1) then
               do i=isave,nx
                  pressw(i)=pressw(isave-1)+(i-isave)/(nx-isave)*
     $                 (1.-pressw(isave-1))
                  pwprim(i)=(pressw(isave-1)-1.)/
     $                 (psilim-psiaxis)*(nx-isave)/(nx-1)
               end do
            endif
         endif
         if(nmass.gt.0) then
            read(neqdsk,300) (rho0(i),i=1,nx)
            dpsiv=(psilim-psiaxis)/(nx-1.)
            do i=2,nx-1
               rho0p(i)=(rho0(i+1)-rho0(i-1))/(2.*dpsiv)
            end do
            rho0p(1)=(-3.*rho0(1)+4.*rho0(2)-rho0(3))/(2.*dpsiv)
            rho0p(nx)=(3.*rho0(nx)-4.*rho0(nx-1)+rho0(nx-2))/(2.*dpsiv)
         endif
      endif
c     generate x,z gridd
 500  continue
      dx=xdim/(nx-1.)
      dz=zdim/(nz-1.)
      do i=1,nx
         xgrid(i)=redge+(i-1.)*dx
      end do
      do i=1,nz
         zgrid(i)=-0.5*zdim+(i-1.)*dz
      end do
c     if sf is negative change it's sign
c     I don't know what the purpose of a negative B field is.
      if(sf(nx).lt.0.) then
         do i=1,nx
            sf(i)=-sf(i)
         end do
      endif
      close(neqdsk)
      call getbpsq(psixz,nxd,nzd,xgrid,dx,dz,nx,nz,bpsq)
!	load g
      write(6,'("begin loading geqdsk")')
      call gfile_allocate(nx,nz,mbdry,limitr,g1)
      read(ntitle(4)(1:6),fmt='(i6)')i
      g1%shot=i
      read(ntitle(5)(1:4),fmt='(i4)')j
      g1%time=j
      g1%source=trim(filename)
      g1%xdim=xdim
      g1%zdim=zdim
      g1%mw=nx
      g1%mh=nz
      g1%rzero=rc
      g1%rgrid1=redge
      g1%zmid=zmid
      g1%rmaxis=xaxis
      g1%zmaxis=zaxis
      g1%ssimag=psiaxis
      g1%ssibry=psilim
      g1%bcentr=btor
      g1%cpasma=totcur
      g1%fpol=sf
      g1%pres=sp
      g1%ffprim=sffp
      g1%pprime=spp
      g1%psirz=psixz(1:nx,1:nz)
      g1%qpsi=qpsi
      g1%nbdry=mbdry
      g1%limitr=limitr
      g1%rbdry(1:mbdry)=xbndry(1:mbdry)
      g1%zbdry(1:mbdry)=zbndry(1:mbdry)
      g1%xlim(1:limitr)=xlim(1:limitr)
      g1%ylim(1:limitr)=zlim(1:limitr)
      write(6,'("end loading geqdsk")')
      end subroutine rdeqdsk
      subroutine  eqdsk_allocate(mw,mh)
      USE mapg_mod
      INTEGER,  INTENT(IN) :: mw,mh
      IF(ALLOCATED(psixz))STOP 'allocation error in reqdsk'
      ALLOCATE(psixz(mw,mh),bpsq(mw,mh),STAT=istat)
        if(istat.ne.0)print*,"STAT=",istat
      ALLOCATE(sp(mw),spp(mw),sf(mw),sffp(mw),qpsi(mw)
     .	,STAT=istat)
        if(istat.ne.0)print*,"STAT=",istat
      ALLOCATE(xgrid(mw),zgrid(mw))
      ALLOCATE(pressw(mw),pwprim(mw),rho0(mw),rho0p(mw),STAT=istat)
        if(istat.ne.0)print*,"STAT=",istat
!      parameter(nxd=129,nzd=129,nxzd=6*(nxd+nzd),nh2=2*nzd,
!     $     nwk=2*nxd*nzd+nh2)
      nxd=mw; nzd=mh;nxzd=6*(nxd+nzd);nh2=2*nzd;nwk=2*nxd*nzd+nh2
! allocate splining arrays for eqdsk
      ALLOCATE(csplpsi(4,nxd,nzd),STAT=istat)
        if(istat.ne.0)print*,"STAT=",istat
! initialize
      pressw=0;pwprim=0;rho0=0;rho0p=0;csplpsi=0;xgrid=0;zgrid=0
      psixz=0;bpsq=0;sp=0;spp=0;sf=0;sffp=0;qpsi=0
      end subroutine  eqdsk_allocate

