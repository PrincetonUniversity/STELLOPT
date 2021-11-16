program test_vmec_to_stella_geometry_interface

  use vmec_to_stella_geometry_interface_mod

  implicit none

  !*********************************************************************
  ! Input parameters
  !*********************************************************************

  character(len=2000) :: vmec_filename = 'equilibria/wout_w7x_standardConfig.nc'
!  character(len=2000) :: vmec_filename = 'equilibria/wout_161s1.nc'
  integer, parameter :: nalpha = 5
  integer, parameter :: nzgrid = 7
  real :: alpha0 = 0.0
  real :: zeta_center = 0.0
  real :: number_of_field_periods_to_include = 1
  real :: desired_normalized_toroidal_flux = 0.6354167d+0
  integer :: vmec_surface_option = 0
  logical :: verbose = .true.

  !*********************************************************************
  ! Output arrays
  !*********************************************************************
  
  real :: normalized_toroidal_flux_used, safety_factor_q, shat, L_reference, B_reference, nfp
  integer :: sign_toroidal_flux
  real, dimension(nalpha) :: alpha
  real, dimension(-nzgrid:nzgrid) :: zeta
  real, dimension(nalpha, -nzgrid:nzgrid) :: bmag, gradpar, gds2, gds21, gds22, gds23, gds24, gds25, gds26
  real, dimension (nalpha, -nzgrid:nzgrid) :: gbdrift, gbdrift0, cvdrift, cvdrift0
  real, dimension(nalpha, -nzgrid:nzgrid) :: theta_vmec
  real, dimension(nalpha, -nzgrid:nzgrid) :: B_sub_zeta, B_sub_theta_vmec
  ! This code uses normalizations in which kxfac is always 1, so kxfac is not presently returned.

  !*********************************************************************
  ! Variables used internally by this program
  !*********************************************************************

  integer :: j, iunit

  !*********************************************************************
  ! Beginning of executable statements
  !*********************************************************************

  call read_vmec_equilibrium (vmec_filename)

  call vmec_to_stella_geometry_interface(nalpha, alpha0, nzgrid, zeta_center, &
       number_of_field_periods_to_include, &
       desired_normalized_toroidal_flux, vmec_surface_option, verbose, &
       normalized_toroidal_flux_used, safety_factor_q, shat, L_reference, B_reference, nfp, &
       sign_toroidal_flux, &
       alpha, zeta, bmag, gradpar, gds2, gds21, gds22, gds23, gds24, gds25, gds26, &
       gbdrift, gbdrift0, cvdrift, cvdrift0, &
       theta_vmec, B_sub_zeta, B_sub_theta_vmec)

  print *,"-------------- Input parameters ------------------"
  print *,"vmec_filename: ",trim(vmec_filename)
  print *,"nalpha:",nalpha
  print *,"alpha0:",alpha0
  print *,"nzgrid:",nzgrid
  print *,"zeta_center:",zeta_center
  print *,"number_of_field_periods_to_include:",number_of_field_periods_to_include
  print *,"desired_normalized_toroidal_flux:",desired_normalized_toroidal_flux
  print *,"vmec_surface_option:",vmec_surface_option

  print *,"-------------- Output parameters -----------------"
  print *,"normalized_toroidal_flux_used:",normalized_toroidal_flux_used
  print *,"safety_factor_q:",safety_factor_q
  print *,"shat:",shat
  print *,"L_reference:",L_reference
  print *,"B_reference:",B_reference
  print *,"nfp:",nfp
  print *,"alpha:"
  print *,alpha
  print *,"zeta:"
  print *,zeta

  print *,"bmag:"
  do j=1,nalpha
     print *,bmag(j,:)
  end do

  print *,"gradpar:"
  do j=1,nalpha
     print *,gradpar(j,:)
  end do

  print *,"gds2:"
  do j=1,nalpha
     print *,gds2(j,:)
  end do

  print *,"gds21:"
  do j=1,nalpha
     print *,gds21(j,:)
  end do

  print *,"gds22:"
  do j=1,nalpha
     print *,gds22(j,:)
  end do

  print *,"gds23:"
  do j=1,nalpha
     print *,gds23(j,:)
  end do

  print *,"gds24:"
  do j=1,nalpha
     print *,gds24(j,:)
  end do

  print *,"gds25:"
  do j=1,nalpha
     print *,gds25(j,:)
  end do

  print *,"gds26:"
  do j=1,nalpha
     print *,gds26(j,:)
  end do

  print *,"gbdrift:"
  do j=1,nalpha
     print *,gbdrift(j,:)
  end do

  print *,"gbdrift0:"
  do j=1,nalpha
     print *,gbdrift0(j,:)
  end do

  print *,"cvdrfit:"
  do j=1,nalpha
     print *,cvdrift(j,:)
  end do

  print *,"cvdrift0:"
  do j=1,nalpha
     print *,cvdrift0(j,:)
  end do

  print *,"theta_vmec:"
  do j=1,nalpha
     print *,theta_vmec(j,:)
  end do

  iunit = 6
  open(file='geometry.dat',unit=iunit)
  write (iunit,*) 'nalpha nzgrid'
  write (iunit,*) nalpha, nzgrid
  write (iunit,*) 'alpha'
  write (iunit,*) alpha
  write (iunit,*) 'zeta'
  write (iunit,*) zeta

  write (iunit,*) 'bmag'
  do j=1,nalpha
     write (iunit,*) bmag(j,:)
  end do

  write (iunit,*) 'gradpar'
  do j=1,nalpha
     write (iunit,*) gradpar(j,:)
  end do

  write (iunit,*) 'gds2'
  do j=1,nalpha
     write (iunit,*) gds2(j,:)
  end do

  write (iunit,*) 'gds21'
  do j=1,nalpha
     write (iunit,*) gds21(j,:)
  end do

  write (iunit,*) 'gds22'
  do j=1,nalpha
     write (iunit,*) gds22(j,:)
  end do

  write (iunit,*) 'gds23'
  do j=1,nalpha
     write (iunit,*) gds23(j,:)
  end do

  write (iunit,*) 'gds24'
  do j=1,nalpha
     write (iunit,*) gds24(j,:)
  end do

  write (iunit,*) 'gds25'
  do j=1,nalpha
     write (iunit,*) gds25(j,:)
  end do

  write (iunit,*) 'gds26'
  do j=1,nalpha
     write (iunit,*) gds26(j,:)
  end do

  write (iunit,*) 'gbdrift'
  do j=1,nalpha
     write (iunit,*) gbdrift(j,:)
  end do

  write (iunit,*) 'gbdrift0'
  do j=1,nalpha
     write (iunit,*) gbdrift0(j,:)
  end do

  write (iunit,*) 'cvdrift'
  do j=1,nalpha
     write (iunit,*) cvdrift(j,:)
  end do

  write (iunit,*) 'cvdrift0'
  do j=1,nalpha
     write (iunit,*) cvdrift0(j,:)
  end do

  write (iunit,*) 'theta_vmec'
  do j=1,nalpha
     write (iunit,*) theta_vmec(j,:)
  end do

  close(iunit)

end program test_vmec_to_stella_geometry_interface
