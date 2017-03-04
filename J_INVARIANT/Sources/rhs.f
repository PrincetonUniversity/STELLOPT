
      subroutine rhs(neq, phi_loc, y, f)
      use B_and_J_Library
      integer neq
      real(rprec), parameter :: zero = 0._dp, one = 1._dp
      real(rprec), dimension(neq) :: y, f
      real(rprec) phi_loc, theta_loc, bf

      select case (device)
       case ("qos")
        if(J_star_opt .ne. 0) call b_eval(theta0,phi_loc,bf)
        if(J_star_opt .eq. 0) bf = b_along_fld_line(phi_loc)
       case ("qas")
        theta_loc = (phi_loc - phi0)*iota_local
        bf = b_along_fld_line(theta_loc)
      end select
c
      if (ep_mu .le. bf .or. istop .eq. 1) then
        f(1) = zero; f(2) = zero; istop = 1
      else if (ep_mu .gt. bf) then
        f(1) = length_factor*sqrt(one - (bf/ep_mu))
        f(2) = length_factor
        istop = 0
      end if

!     end subroutine rhs
      end
