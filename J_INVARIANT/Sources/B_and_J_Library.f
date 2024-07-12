!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
       module B_and_J_Library
       use BJdata

       contains

       subroutine b_eval(theta,zeta,bfield)
       implicit none
       real(rprec), parameter :: zero = 0._dp
       real(rprec) :: theta, zeta, bfield
       real(rprec), dimension(max_bmns) :: temp
       temp(:) = blocal_rdcd(:)*cos(theta*xm_rdcd(:)
     >   - zeta*xn_rdcd(:))
       bfield = sum(temp)

       end subroutine b_eval

       function b_along_fld_line(arg)  result(bb)
       implicit none
       real(rprec) :: zeta, theta, bb, arg
       real(rprec), dimension(max_bmns) :: temp
       select case (device)
        case ("qos")
         zeta = arg
         theta = theta0 + iota_local*zeta
        case ("qas")
         theta = arg
         zeta = theta/iota_local + phi0
        case default
         write(*,'("Need to select either qo or qa device")')
         stop 23
        end select
        temp(:) = blocal_rdcd(:)*cos(theta*xm_rdcd(:)
     >   - zeta*xn_rdcd(:))
        bb = sum(temp)
c        call b_eval(theta,zeta,bb)
       end function b_along_fld_line
c
       function invb_along_fld_line(arg)  result(bi)
       implicit none
       real(rprec) :: zeta, theta, bb, bi, arg
       real(rprec), dimension(max_bmns) :: temp
       select case (device)
        case ("qos")
         zeta = arg
         theta = theta0 + iota_local*zeta
        case ("qas")
         theta = arg
         zeta = theta/iota_local + phi0
        case default
         write(*,'("Need to select either qo or qa device")')
         stop 23
        end select
        temp(:) = blocal_rdcd(:)*cos(theta*xm_rdcd(:)
     >   - zeta*xn_rdcd(:))
        bb = sum(temp)
c        call b_eval(theta,zeta,bb)
        bi = 1./bb
       end function invb_along_fld_line
c

       function b_along_const_theta(zeta)  result(bb)
       implicit none
       real(rprec) :: zeta, theta, bb
         theta = theta0
        call b_eval(theta,zeta,bb)
       end function b_along_const_theta
c

       function invb_along_const_theta(zeta)  result(bi)
       implicit none
       real(rprec) :: zeta, theta, bb, bi
        theta = theta0
        call b_eval(theta,zeta,bb)
        bi = 1./bb
       end function invb_along_const_theta
c
       end module B_and_J_Library
