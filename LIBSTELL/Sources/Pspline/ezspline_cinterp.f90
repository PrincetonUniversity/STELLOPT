subroutine EZspline_cinterp(spline_o, k, p1, p2, p3, f, ier)
  ! list of coordinate triplets
  ! this version of interp that has been inlined for better performance on the
  ! crays (pletzer@pppl.gov) Thu Jun 15 08:52:46 PDT 2000
  ! NOTE:
  ! This routine does not require the -dp switch on Crays since
  ! all constant have been made kind=r8.
 
  ! DMC Apr 2007 -- only the tricubic SPLINE evaluation has been inlined;
  ! if a different method interpolation is used, the standard f77-style
  ! pspline subroutine is called.
 
  use EZspline_obj
  implicit none
  type(EZspline3_r8) spline_o
  integer, intent(in) :: k
  real(ezspline_r8), intent(in) :: p1(k), p2(k), p3(k)  ! location arrays
  real(ezspline_r8), intent(out):: f(k)  ! interpolant array
  integer, intent(out) :: ier
 
  integer :: ifail
  integer, parameter :: ict(10)=(/1,0,0,0,0,0,0,0,0,0/)
  integer:: iwarn=0
 
      integer ix(k)                  ! zone indices {j}
      REAL(EZSPLINE_R8) dxn(k)       ! normalized displacements w/in zones
      REAL(EZSPLINE_R8) hx(k)        ! h(j) vector
      REAL(EZSPLINE_R8) hxi(k)       ! 1/h(j) vector
!
      integer iy(k)                  ! zone indices {j}
      REAL(EZSPLINE_R8) dyn(k)       ! normalized displacements w/in zones
      REAL(EZSPLINE_R8) hy(k)        ! h(j) vector
      REAL(EZSPLINE_R8) hyi(k)       ! 1/h(j) vector
!
      integer iz(k)                  ! zone indices {j}
      REAL(EZSPLINE_R8) dzn(k)       ! normalized displacements w/in zones
      REAL(EZSPLINE_R8) hz(k)        ! h(j) vector
      REAL(EZSPLINE_R8) hzi(k)       ! 1/h(j) vector
      integer v
!
      REAL(EZSPLINE_R8) , parameter ::  sixth=0.166666666666666667_ezspline_r8
      REAL(EZSPLINE_R8) , parameter ::  one=1.0_ezspline_r8
      REAL(EZSPLINE_R8) , parameter ::  three=3.0_ezspline_r8
 
      integer j1, j2 , j3
      REAL(EZSPLINE_R8) z36th,z216th,xp,xpi,xp2,xpi2,cx,cxi,hx2,cxd,cxdi,yp,       &
     & ypi,yp2,ypi2,cy,cyi,hy2,cyd,cydi,zp,zpi,zp2,zpi2,cz,czi,hz2,     &
     & czd,czdi,somme
 
  ier = 0
  ifail=0
  if(spline_o%isReady /= 1) then
     ier = 94
     return
  endif
 
  if (spline_o%isLinear == 1) then

     call r8vecpc3(ict, k, p1, p2, p3, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn,ifail)

  else if (spline_o%isHybrid == 1) then
  
     call r8vecintrp3d(ict, k, p1, p2, p3, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%hspline, &
          & spline_o%fspl(1,1,1,1), size(spline_o%fspl,1), &
          & size(spline_o%fspl,2), size(spline_o%fspl,3), &
          & size(spline_o%fspl,4), iwarn, ifail)

  else if (spline_o%isHermite ==0) then
     !
!!$     call r8vectricub(ict, k, p1, p2, p3, k, f, &
!!$          & spline_o%n1,spline_o%x1pkg(1,1), &
!!$          & spline_o%n2,spline_o%x2pkg(1,1), &
!!$          & spline_o%n3,spline_o%x3pkg(1,1), &
!!$          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
!!$          & iwarn, ifail)
 
      call r8xlookup(k,p1,spline_o%n1,spline_o%x1pkg(1,1),2,ix,dxn,hx,hxi,ifail)
      call r8xlookup(k,p2,spline_o%n2,spline_o%x2pkg(1,1),2,iy,dyn,hy,hyi,ifail)
      call r8xlookup(k,p3,spline_o%n3,spline_o%x3pkg(1,1),2,iz,dzn,hz,hzi,ifail)
 
!!$      call r8fvtricub(ict,k,k,f,ix,iy,iz,dxn,dyn,dzn,           &
!!$     &   hx,hxi,hy,hyi,hz,hzi,spline_o%fspl(1,1,1,1), &
!!$     &   spline_o%n1, spline_o%n2, spline_o%n3)
 
 
      z36th=sixth*sixth
      z216th=sixth*sixth*sixth
!
!  prepare useful parameters...
!
      do v=1,k
         j1=ix(v)
         j2=iy(v)
         j3=iz(v)
!
!   ...in x direction
!
         xp=dxn(v)
         xpi=one-xp
         xp2=xp*xp
         xpi2=xpi*xpi
!
!!$         if((ict(1).eq.1).or.(ict(3).eq.1).or.(ict(4).eq.1).or.         &
!!$     &      (ict(6).eq.1).or.(ict(7).eq.1).or.(ict(10).eq.1)) then
            cx=xp*(xp2-ONE)
            cxi=xpi*(xpi2-ONE)
            hx2=hx(v)*hx(v)
!!$         endif
!!$         if((ict(2).eq.1).or.(ict(8).eq.1).or.(ict(9).eq.1)) then
!!$            cxd=THREE*xp2-ONE
!!$            cxdi=-THREE*xpi2+ONE
!!$         endif
!
!   ...and in y direction
!
         yp=dyn(v)
         ypi=ONE-yp
         yp2=yp*yp
         ypi2=ypi*ypi
!
!!$         if((ict(1).eq.1).or.(ict(2).eq.1).or.(ict(4).eq.1).or.         &
!!$     &      (ict(5).eq.1).or.(ict(7).eq.1).or.(ict(9).eq.1)) then
            cy=yp*(yp2-ONE)
            cyi=ypi*(ypi2-ONE)
            hy2=hy(v)*hy(v)
!!$         endif
!!$         if((ict(3).eq.1).or.(ict(8).eq.1).or.(ict(10).eq.1)) then
!!$            cyd=THREE*yp2-ONE
!!$            cydi=-THREE*ypi2+ONE
!!$         endif
!
!   ...and in z direction
!
         zp=dzn(v)
         zpi=ONE-zp
         zp2=zp*zp
         zpi2=zpi*zpi
!
!!$         if((ict(1).eq.1).or.(ict(2).eq.1).or.(ict(3).eq.1).or.         &
!!$     &      (ict(5).eq.1).or.(ict(6).eq.1).or.(ict(8).eq.1)) then
            cz=zp*(zp2-ONE)
            czi=zpi*(zpi2-ONE)
            hz2=hz(v)*hz(v)
!!$         endif
!!$         if((ict(4).eq.1).or.(ict(9).eq.1).or.(ict(10).eq.1)) then
!!$            czd=THREE*zp2-ONE
!!$            czdi=-THREE*zpi2+ONE
!!$         endif
!
!  get desired values:
!
         if(ict(1).eq.1) then
!
!  function value:
!
            somme=(                                                       &
     &         zpi*(                                                    &
     &           xpi*(ypi*spline_o%fspl(1,j1,j2,j3)  +yp*spline_o%fspl(1,j1,j2+1,j3))+            &
     &            xp*(ypi*spline_o%fspl(1,j1+1,j2,j3)+yp*spline_o%fspl(1,j1+1,j2+1,j3)))          &
     &        +zp*(                                                     &
     &           xpi*(ypi*spline_o%fspl(1,j1,j2,j3+1)  +yp*spline_o%fspl(1,j1,j2+1,j3+1))+        &
     &            xp*(ypi*spline_o%fspl(1,j1+1,j2,j3+1)+yp*spline_o%fspl(1,j1+1,j2+1,j3+1))))
!
            somme=somme+sixth*hx2*(                                         &
     &         zpi*(                                                    &
     &           cxi*(ypi*spline_o%fspl(2,j1,j2,j3)  +yp*spline_o%fspl(2,j1,j2+1,j3))+            &
     &            cx*(ypi*spline_o%fspl(2,j1+1,j2,j3)+yp*spline_o%fspl(2,j1+1,j2+1,j3)))          &
     &        +zp*(                                                     &
     &           cxi*(ypi*spline_o%fspl(2,j1,j2,j3+1)  +yp*spline_o%fspl(2,j1,j2+1,j3+1))+        &
     &            cx*(ypi*spline_o%fspl(2,j1+1,j2,j3+1)+yp*spline_o%fspl(2,j1+1,j2+1,j3+1))))
!
            somme=somme+sixth*hy2*(                                         &
     &         zpi*(                                                    &
     &           xpi*(cyi*spline_o%fspl(3,j1,j2,j3)  +cy*spline_o%fspl(3,j1,j2+1,j3))+            &
     &            xp*(cyi*spline_o%fspl(3,j1+1,j2,j3)+cy*spline_o%fspl(3,j1+1,j2+1,j3)))          &
     &        +zp*(                                                     &
     &           xpi*(cyi*spline_o%fspl(3,j1,j2,j3+1)  +cy*spline_o%fspl(3,j1,j2+1,j3+1))+        &
     &            xp*(cyi*spline_o%fspl(3,j1+1,j2,j3+1)+cy*spline_o%fspl(3,j1+1,j2+1,j3+1))))
!
            somme=somme+sixth*hz2*(                                         &
     &         czi*(                                                    &
     &           xpi*(ypi*spline_o%fspl(4,j1,j2,j3)  +yp*spline_o%fspl(4,j1,j2+1,j3))+            &
     &            xp*(ypi*spline_o%fspl(4,j1+1,j2,j3)+yp*spline_o%fspl(4,j1+1,j2+1,j3)))          &
     &        +cz*(                                                     &
     &           xpi*(ypi*spline_o%fspl(4,j1,j2,j3+1)  +yp*spline_o%fspl(4,j1,j2+1,j3+1))+        &
     &            xp*(ypi*spline_o%fspl(4,j1+1,j2,j3+1)+yp*spline_o%fspl(4,j1+1,j2+1,j3+1))))
!
            somme=somme+z36th*hx2*hy2*(                                     &
     &         zpi*(                                                    &
     &           cxi*(cyi*spline_o%fspl(5,j1,j2,j3)  +cy*spline_o%fspl(5,j1,j2+1,j3))+            &
     &            cx*(cyi*spline_o%fspl(5,j1+1,j2,j3)+cy*spline_o%fspl(5,j1+1,j2+1,j3)))          &
     &        +zp*(                                                     &
     &           cxi*(cyi*spline_o%fspl(5,j1,j2,j3+1)  +cy*spline_o%fspl(5,j1,j2+1,j3+1))+        &
     &            cx*(cyi*spline_o%fspl(5,j1+1,j2,j3+1)+cy*spline_o%fspl(5,j1+1,j2+1,j3+1))))
!
            somme=somme+z36th*hx2*hz2*(                                     &
     &         czi*(                                                    &
     &           cxi*(ypi*spline_o%fspl(6,j1,j2,j3)  +yp*spline_o%fspl(6,j1,j2+1,j3))+            &
     &            cx*(ypi*spline_o%fspl(6,j1+1,j2,j3)+yp*spline_o%fspl(6,j1+1,j2+1,j3)))          &
     &        +cz*(                                                     &
     &           cxi*(ypi*spline_o%fspl(6,j1,j2,j3+1)  +yp*spline_o%fspl(6,j1,j2+1,j3+1))+        &
     &            cx*(ypi*spline_o%fspl(6,j1+1,j2,j3+1)+yp*spline_o%fspl(6,j1+1,j2+1,j3+1))))
!
            somme=somme+z36th*hy2*hz2*(                                     &
     &         czi*(                                                    &
     &           xpi*(cyi*spline_o%fspl(7,j1,j2,j3)  +cy*spline_o%fspl(7,j1,j2+1,j3))+            &
     &            xp*(cyi*spline_o%fspl(7,j1+1,j2,j3)+cy*spline_o%fspl(7,j1+1,j2+1,j3)))          &
     &        +cz*(                                                     &
     &           xpi*(cyi*spline_o%fspl(7,j1,j2,j3+1)  +cy*spline_o%fspl(7,j1,j2+1,j3+1))+        &
     &            xp*(cyi*spline_o%fspl(7,j1+1,j2,j3+1)+cy*spline_o%fspl(7,j1+1,j2+1,j3+1))))
!
            somme=somme+z216th*hx2*hy2*hz2*(                                &
     &         czi*(                                                    &
     &           cxi*(cyi*spline_o%fspl(8,j1,j2,j3)  +cy*spline_o%fspl(8,j1,j2+1,j3))+            &
     &            cx*(cyi*spline_o%fspl(8,j1+1,j2,j3)+cy*spline_o%fspl(8,j1+1,j2+1,j3)))          &
     &        +cz*(                                                     &
     &           cxi*(cyi*spline_o%fspl(8,j1,j2,j3+1)  +cy*spline_o%fspl(8,j1,j2+1,j3+1))+        &
     &            cx*(cyi*spline_o%fspl(8,j1+1,j2,j3+1)+cy*spline_o%fspl(8,j1+1,j2+1,j3+1))))
!
            f(v)=somme
         else
            print*,'ERROR ict(1) must be 1 (interpolation)'
         endif
!
!
      enddo                             ! vector loop
 
 
     if(ifail /= 0) ier = 97
 
  else
 
     call r8vecherm3(ict, k, p1, p2, p3, k, f, &
          & spline_o%n1, spline_o%x1pkg(1,1), &
          & spline_o%n2, spline_o%x2pkg(1,1), &
          & spline_o%n3, spline_o%x3pkg(1,1), &
          & spline_o%fspl(1,1,1,1), spline_o%n1, spline_o%n2, &
          & iwarn,ifail)
 
  endif
 
  if(ifail /= 0) ier = 97
 
end subroutine
