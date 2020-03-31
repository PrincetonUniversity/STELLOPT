
      subroutine j_inv_calc( js, num_ep_mu, NumJstar, fix_pitch,
     1     min_epmu, max_epmu, J_inv_all, epmu)
c ******************************************************************************
      use read_boozer_mod
      use B_and_J_Library
      implicit none
      external rhs
      external jac
      real(rprec), external :: sfmin
!-----------------------------------------------
!   Arguments
!-----------------------------------------------
      integer :: js, num_ep_mu, NumJstar
      logical :: fix_pitch
      real(rprec) :: min_epmu, max_epmu
      real(rprec), dimension(num_ep_mu) :: epmu
      real(rprec), dimension(NumJstar, num_ep_mu) :: J_inv_all

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(rprec), parameter :: pm4 = 1.e-4_dp, zero = 0._dp,
     >   one = 1._dp
      real(rprec), dimension(mnboz_b) :: bmn_local,
     >  abs_blocal, blocal_ordered
      integer, dimension(mnboz_b) :: m_ordered, n_ordered
      integer i, j, k, ierr, istat, nfp, iloc, numargs, iunit,
     >        iunit2, ks
      real(rprec) :: btheta, bzeta, psip, chip, an
      real(rprec) :: phi_edge, phi_center, phi_lower, phi_upper
      real(rprec), dimension(1) :: check
      real(rprec), dimension(NumJstar) :: phimin, phimax,
     >  bmin, bmax
      real(rprec), dimension(NumJstar) :: theta_J, thmn, thmx
      real(rprec), dimension(NumJstar) :: J_plus, J_minus, J_inv
      real(rprec), dimension(NumJstar) :: L_plus, L_minus, L
      integer :: neqn, itol, istate, itask, jt, istep, iopt
      integer :: liw, lrw, n_upper, jrad, istep_max
      real(rprec) :: thet, phi, bf, bf1, time_begin,
     >  time_end,time_all
      real(rprec), dimension(2) :: y, f
      real(rprec), dimension(55) :: rwork
      integer, dimension(22) :: iwork
      integer, dimension(ns_b,NumJstar) :: int_steps
      integer, dimension(NumJstar) :: int_steps_plus, int_steps_minus
      real(rprec) :: delta_phi, rtol, atol, phi_in, phi_out,
     >  maxbmin, minbmax, minbmin, maxbmax, theta_min
      real(rprec) :: theta_center, theta_lower, theta_upper,
     >  bmn_test
      real(rprec) :: J_avg, J_min, J_max, width_epmu
      integer index_dat, index_end, i_ep_mu
      integer :: imon
      logical :: lscreen, ljconf

      ep_mu = 1._dp

c
c   Get local (at flux surface js) value of quantities
c     needed for DKES run
c
      bmn_test = zero
      do i = 1,mnboz_b
         bmn_local(i) = bmnc_b(i,js)
         abs_blocal(i) = abs(bmn_local(i))
         bmn_test = bmn_test + abs_blocal(i)
      end do
      if (bmn_test .eq. zero) then
        write(*,'("   Requested surface has not been mapped")')
        write(*,'("(i.e., sum of |Bmns| at this surface is zero)")')
        stop 30
      endif
      iota_local = iota_b(js)
      length_factor = sqrt(one + iota_local*iota_local)
      btheta = buco_b(js)
      bzeta = bvco_b(js)
      phi_edge = abs(phi_b(ns_b))
      psip = phi_edge/TWOPI
      chip = psip*iota_local
      nfp = nfp_b
      an = nfp

c
c     Sort the Bmn's (along with associated m's and n's in order
c     of increasing abs(Bmn):
c
      do i=1,mnboz_b

c      check = maxval(abs_blocal)
c    Find location of this value of check in the abs_blocal array
c        do j=1,mnboz_b
c         if(abs(check - abs_blocal(j)) .lt. 1.d-6*abs_blocal(j))
c    >      iloc = j
c        end do

       check = maxloc(abs_blocal)
       iloc = check(1)
       blocal_ordered(i) = bmn_local(iloc)
       m_ordered(i) = ixm_b(iloc)
       n_ordered(i) = ixn_b(iloc)
       abs_blocal(iloc) = zero
      end do
c
c    Keep max_bmns components for use in module B_and_J_Library
c
      do i=1, max_bmns
        m_rdcd(i) = m_ordered(i)
        xm_rdcd(i) = m_ordered(i)
        n_rdcd(i) = n_ordered(i)
        xn_rdcd(i) = n_ordered(i)
        blocal_rdcd(i) = blocal_ordered(i)
c        write(*,'(i5,2x,i5,2x,e15.7)') m_rdcd(i), n_rdcd(i),
c     >     blocal_rdcd(i)
      end do
c
c    Find bmin and bmax's along adjacent field lines:
c
      do i=1,NumJstar
c
      select case (device)
c
      case("qos")
       if(J_star_opt .eq. 0) then
         theta0 = PI*(one - (iota_local/an))*
     >      (i-1)/real(NumJstar-1,rprec)
         phi_center = theta0/(an - iota_local)
c
       else if(J_star_opt .ne. 0) then
         theta0 = PI*(i-1)/real(NumJstar-1,rprec)
         phi_center = theta0/an
       endif
       if(i .eq. 1) then
        phi_lower = phi_center - (PI/(5*an))
        phi_upper = phi_center + (PI/(5*an))
       else if(i .gt. 1) then
        phi_lower = phi_center - (PI/(2*an))
        phi_upper = phi_center + (PI/(2*an))
       endif
       if(J_star_opt .eq. 0) then
        phimin(i) = sfmin(phi_lower,phi_upper,
     >      b_along_fld_line,pm4)
        bmin(i) = b_along_fld_line(phimin(i))
c
       else if(J_star_opt .ne. 0) then
        phimin(i) = sfmin(phi_lower,phi_upper,
     >      b_along_const_theta,pm4)
        bmin(i) = b_along_const_theta(phimin(i))
       endif
c
       phi_center = phi_center + (PI/an)
       phi_lower = phi_lower + (PI/an)
       phi_upper = phi_upper + (PI/an)
       if(J_star_opt .eq. 0) then
        phimax(i) = sfmin(phi_lower,phi_upper,
     >       invb_along_fld_line,pm4)
        bmax(i) = b_along_fld_line(phimax(i))
        thmn(i) = theta0 + iota_local*phimin(i)
        thmx(i) = theta0 + iota_local*phimax(i)
c
       else if(J_star_opt .ne. 0) then
        phimax(i) = sfmin(phi_lower,phi_upper,
     >      invb_along_const_theta,pm4)
        bmax(i) = b_along_const_theta(phimax(i))
        thmn(i) = theta0
        thmx(i) = theta0
       endif
c
       case("qas")
c
        phi0 = TWOPI*(i-1)/(an*(NumJstar-1))
        theta_center = zero
        theta_lower = theta_center - (PI/4)
        theta_upper = theta_center + (PI/4)

        thmn(i) = sfmin(theta_lower,theta_upper,
     >      b_along_fld_line,pm4)
        bmin(i) = b_along_fld_line(thmn(i))
        phimin(i) = thmn(i)/iota_local + phi0
c
        theta_center = PI
        theta_lower = theta_center - (PI/4)
        theta_upper = theta_center + (PI/4)
        thmx(i) = sfmin(theta_lower,theta_upper,
     >      invb_along_fld_line,pm4)
        bmax(i) = b_along_fld_line(thmx(i))
c
       case default
         write(*,'("Need to select either qo or qa device")')
         stop 23
       end select
c
      end do         ! do i=1,NumJstar
c
      maxbmin = maxval(bmin)
      minbmax = minval(bmax)
      minbmin = minval(bmin)
      maxbmax = maxval(bmax)
      width_epmu = minbmax - maxbmin

      if( .not. fix_pitch) then
         if( num_ep_mu <= 10) then
            min_epmu = maxbmin + 0.05_dp * width_epmu
            width_epmu = 0.9_dp * width_epmu
         else
            min_epmu = maxbmin + width_epmu/num_ep_mu/2
            width_epmu = width_epmu*(num_ep_mu-1)/num_ep_mu
         endif
      else
         width_epmu = max_epmu - min_epmu
      endif
c
      do i_ep_mu = 1,num_ep_mu
       if(num_ep_mu .ne. 1) then
c        ep_mu = 1.05_dp*maxbmin + ((0.95_dp*minbmax - 1.05_dp*maxbmin)
c     >    *(i_ep_mu - 1)/real(num_ep_mu - 1,rprec))
c        ep_mu = maxbmin + 0.05_dp*width_epmu + (0.9_dp*width_epmu
c    >    *(i_ep_mu - 1)/real(num_ep_mu - 1,rprec))
         ep_mu = min_epmu + width_epmu
     >    * (i_ep_mu - 1)/real(num_ep_mu - 1,rprec)

       else if(num_ep_mu .eq. 1) then
         ep_mu = min_epmu + 0.5_dp * width_epmu
       endif

c      write(*,'("ep/mu = ",f8.4)') ep_mu
c
      do i=1,NumJstar
       J_plus(i)=zero; L_plus(i)=zero; J_inv(i)=zero
       J_minus(i)=zero; L_minus(i)=zero; L(i) = zero
      end do
c
      do i=1,NumJstar
c
      if (ep_mu .ge. maxbmax) cycle     !Passing orbit
      if (ep_mu .le. minbmin) cycle     !Forbidden (i.e., v||**2 < 0) orbit
c
      select case (device)
c
      case ("qos")
       if(J_star_opt .eq. 0) then
        theta0 = PI*(one - (iota_local/an))*(i-1)
     >               /real(NumJstar-1,rprec)
        theta_min = theta0 + iota_local*phimin(i)
       else if(J_star_opt .ne. 0) then
        theta0 = PI*(i-1)/real(NumJstar-1,rprec)
       endif
        delta_phi = TWOPI/(10*an)
        istep_max = 10
      case ("qas")
       phi0 = TWOPI*(i-1)/(an*(NumJstar-1))
       delta_phi = TWOPI/(10*an)
       istep_max = int(20/iota_local)
      case default
       write(*,'("Need to select either qo or qa device")')
       stop 23
      end select
c
       y(1)=zero; y(2)=zero; istop=0; neqn=2; itol=1; iopt=0
       istate=1; itask=1; jt=10; rtol=1.d-6; atol=1.d-7
       lrw = 55; liw = 22
       rwork(5) = delta_phi
       do istep = 1,istep_max
c
      select case (device)
c
      case ("qos")
        phi_in = phimin(i) + delta_phi*(istep - 1)
        phi_out = phimin(i) + delta_phi*istep
      case ("qas")
        phi_in = phi0 + delta_phi*(istep - 1)
        phi_out = phi0 + delta_phi*istep
      case default
       write(*,'("Need to select either qo or qa device")')
       stop 23
      end select
c
        call lsode(rhs,neqn,y,phi_in,phi_out,itol,rtol,
     >  atol,itask,istate,iopt,rwork,lrw,iwork,liw,jac,jt)
        if(istop .eq. 1) exit
       end do
       int_steps_plus(i) = istep
       if(istep .lt. istep_max) J_plus(i) = abs(y(1))
       if(istep .ge. istep_max) J_plus(i) = 0
       L_plus(i) = abs(y(2))
c
      select case (device)
c
      case ("qos")
       delta_phi = TWOPI/(10*an)
       istep_max = 10
      case ("qas")
       delta_phi = TWOPI/(10*an)
       istep_max = int(20/iota_local)
      case default
       write(*,'("Need to select either qo or qa device")')
       stop 23
      end select
c
       y(1)=zero; y(2)=zero; istop=0; neqn=2; itol=1;iopt=0
       istate=1; itask=1; jt=10; rtol=1.d-6; atol=1.d-7
       lrw = 55; liw = 22
       rwork(5) = delta_phi
       do istep = 1,istep_max
c
      select case (device)
      case ("qos")
        phi_in = phimin(i) - delta_phi*(istep - 1)
        phi_out = phimin(i) - delta_phi*(istep)
      case ("qas")
        phi_in = phi0 - delta_phi*(istep - 1)
        phi_out = phi0 - delta_phi*(istep)
      case default
       write(*,'("Need to select either qo or qa device")')
       stop 23
      end select
c
        call lsode(rhs,neqn,y,phi_in,phi_out,itol,rtol,
     >  atol,itask,istate,iopt,rwork,lrw,iwork,liw,jac,jt)
        if(istop .eq. 1) exit
       end do
       int_steps_minus(i) = istep
       if(istep .lt. istep_max) J_minus(i) = abs(y(1))
       if(istep .ge. istep_max) J_minus(i) = zero
       L_minus(i) = abs(y(2))
       J_inv(i) = J_plus(i) + J_minus(i)
       L(i) = L_plus(i) + L_minus(i)
       int_steps(js,i) = int_steps_plus(i) + int_steps_minus(i)
c
      end do                                         !do i=1,NumJstar

      j_inv_all(:,i_ep_mu) = j_inv(:)
      epmu(i_ep_mu) = ep_mu
      end do                                !do i_ep_mu = 1,num_ep_mu

      end subroutine j_inv_calc
