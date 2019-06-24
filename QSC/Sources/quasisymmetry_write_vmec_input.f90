subroutine quasisymmetry_write_vmec_input

  use quasisymmetry_variables
  use safe_open_mod
  use vmec_input, only: vmec_nfp => nfp, lasym, ntor, raxis_cc, raxis_cs, zaxis_cc, zaxis_cs, &
       read_indata_namelist, write_indata_namelist, lfreeb, RBC, RBS, ZBC, ZBS, phiedge, AM, PMASS_TYPE, PRES_SCALE, &
       CURTOR, NCURR, PCURR_TYPE, AC, &
       nzeta,indata  !11/27,28/18.
  use vparams, only: ntord
  use quasisymmetry_Frenet_to_cylindrical_mod

  implicit none

  integer :: iunit = 40, istat
!  integer :: max_n, n (8o7)c-out
  integer :: n   !1/29/19.(8o7)
  real(dp) :: half_sum, half_difference, temp

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Read in the VMEC template namelist:
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!12/20/18.(7m12a)c-out rd of vmec_template_filename. &indata vals should come internally.
!  istat = 0
!  call safe_open(iunit,istat,trim(vmec_template_filename),'old','formatted')
!  if (istat /= 0) then
!     print "(a,a)"," Error opening vmec_template_filename: ",trim(vmec_template_filename)
!     stop
!  end if

!!  call read_indata_namelist(iunit, istat)  !11/28/18.c-out
!  read(iunit,nml=indata,iostat=istat)  !11/28/18.(7l21c) = core of rd_indata_namelist.

!  if (istat /= 0) then
!     print "(a,a)"," Error reading &indata namelist from vmec_template_filename: ",trim(vmec_template_filename)
!     stop
!  end if

!  close(iunit)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Over-write a few vmec parameters
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  vmec_nfp = nfp
  lfreeb = .false.
!4/5/19.(9q7)c-out this sec.
!  phiedge = pi * r * r * sign_psi * B0
!
!  ! Set pressure profile:
!  temp = - p2 * r * r
!  AM = 0
!  AM(0) = temp  ! Note: AM uses 0-based indicies!
!  AM(1) = -temp
!  PMASS_TYPE='power_series'
!  PRES_SCALE=1
!
!  ! Set current profile:
!  NCURR = 1
!  PCURR_TYPE = 'power_series'
!  AC = 0
!  AC(0) = 1 ! Note: AC uses 0-based indicies!
!  CURTOR = 2 * pi / mu0 * I2_over_B0 * B0 * r * r
!end (9q7)c-out.

  ! The output is not stellarator-symmetric if (1) R0s is nonzero, (2) Z0c is nonzero, or (3) sigma_initial is nonzero, or (4) B2s is nonzero:
  !lasym = (maxval(abs(R0s))>0 .or. maxval(abs(Z0c)) > 0 .or. abs(sigma_initial) > 0 .or. (order_r_squared .and. abs(B2s)>0))
  lasym = (maxval(abs(R0s))>0 .or. maxval(abs(Z0c)) > 0 .or. abs(sigma_initial) > 0 .or. ((trim(order_r_option).ne.order_r_option_r1) .and. abs(B2s)>0))

  ! We should be able to resolve (N_phi-1)/2 modes (note integer division!), but in case N_phi is very large, don't attempt more than the vmec arrays can handle.
!  ntor = min((N_phi - 1) / 2, ntord)  1/18/19.(8k11d)back out.

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Over-write vmec axis shape
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  raxis_cc = 0;  raxis_cs = 0;  zaxis_cc = 0;  zaxis_cs = 0  -out-12/21/18.ever needed?
!12/21/18.(7m12c)xfer max_n, [raxis_cc,..]=[R0c,..] initializaton to qs_rd_input.
!  max_n = min(ntord,max_axis_nmax)  !1/29/19.(8o7)c-out.

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Compute new RBC, RBS, ZBC, ZBS
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  RBC = 0
  RBS = 0
  ZBC = 0
  ZBS = 0

  select case (trim(finite_r_option))
  case (finite_r_option_nonlinear)
     call quasisymmetry_Frenet_to_cylindrical_nonlinear()

  case (finite_r_option_linear)
     if (trim(order_r_option)==order_r_option_r1) then
        mpol_nonzero = 1
     elseif (trim(order_r_option)==order_r_option_r2) then
        mpol_nonzero = 2
     else
        mpol_nonzero = 3
     end if

     ! Handle the m=0 modes of the boundary, which are the same as the axis shape through O(r), but different to O(r^2):
     RBC(0:max_n,0) = raxis_cc(0:max_n)
     RBS(0:max_n,0) = raxis_cs(0:max_n)
     ZBC(0:max_n,0) = zaxis_cc(0:max_n)
     ZBS(0:max_n,0) = zaxis_cs(0:max_n)
     
     ! Handle the n=0 m=1 modes:
     RBC(0,1) = r * sum(R1c) / N_phi
     RBS(0,1) = r * sum(R1s) / N_phi
     ZBC(0,1) = r * sum(Z1c) / N_phi
     ZBS(0,1) = r * sum(Z1s) / N_phi
     
     ! Handle the m=1 modes that have nonzero n:
     do n = 1, ntor
        ! RBC:
        half_sum        =  r * sum(R1c * cos_n_phi(:,n+1)) / N_phi
        half_difference =  r * sum(R1s * sin_n_phi(:,n+1)) / N_phi
        RBC( n,1) = half_sum + half_difference
        RBC(-n,1) = half_sum - half_difference
        
        ! ZBC:
        half_sum        =  r * sum(Z1c * cos_n_phi(:,n+1)) / N_phi
        half_difference =  r * sum(Z1s * sin_n_phi(:,n+1)) / N_phi
        ZBC( n,1) = half_sum + half_difference
        ZBC(-n,1) = half_sum - half_difference
        
        ! RBS:
        half_sum        =  r * sum(R1s * cos_n_phi(:,n+1)) / N_phi
        half_difference = -r * sum(R1c * sin_n_phi(:,n+1)) / N_phi
        RBS( n,1) = half_sum + half_difference
        RBS(-n,1) = half_sum - half_difference
        
        ! ZBS:
        half_sum        =  r * sum(Z1s * cos_n_phi(:,n+1)) / N_phi
        half_difference = -r * sum(Z1c * sin_n_phi(:,n+1)) / N_phi
        ZBS( n,1) = half_sum + half_difference
        ZBS(-n,1) = half_sum - half_difference
     end do
     
     if (trim(order_r_option) .ne. order_r_option_r1) then
        ! Add O(r^2) terms.
        ! The transformation below can be verified using
        ! m20190215_03_checkFourierTransformInQuasisymmetryCode.m

        ! Handle the n=0 m=0 modes:
        RBC(0,0) = RBC(0,0) + r * r * sum(R20) / N_phi
        ZBC(0,0) = ZBC(0,0) + r * r * sum(z20_cylindrical) / N_phi
        
        ! Handle the m=0 modes that have nonzero n:
        ! Note the sin terms have a minus sign, for the reason explained above.
        do n = 1, ntor
           RBC(n,0) = RBC(n,0) + 2 * r * r * sum(R20             * cos_n_phi(:,n+1)) / N_phi
           RBS(n,0) = RBS(n,0) - 2 * r * r * sum(R20             * sin_n_phi(:,n+1)) / N_phi
           ZBC(n,0) = ZBC(n,0) + 2 * r * r * sum(z20_cylindrical * cos_n_phi(:,n+1)) / N_phi
           ZBS(n,0) = ZBS(n,0) - 2 * r * r * sum(z20_cylindrical * sin_n_phi(:,n+1)) / N_phi
        end do

        ! Handle the n=0 m=2 modes:
        RBC(0,2) = r * r * sum(R2c) / N_phi
        RBS(0,2) = r * r * sum(R2s) / N_phi
        ZBC(0,2) = r * r * sum(z2c_cylindrical) / N_phi
        ZBS(0,2) = r * r * sum(z2s_cylindrical) / N_phi
        
        ! Handle the m=2 modes that have nonzero n:
        do n = 1, ntor
           ! RBC:
           half_sum        =  r * r * sum(R2c * cos_n_phi(:,n+1)) / N_phi
           half_difference =  r * r * sum(R2s * sin_n_phi(:,n+1)) / N_phi
           RBC( n,2) = half_sum + half_difference
           RBC(-n,2) = half_sum - half_difference
           
           ! ZBC:
           half_sum        =  r * r * sum(z2c_cylindrical * cos_n_phi(:,n+1)) / N_phi
           half_difference =  r * r * sum(z2s_cylindrical * sin_n_phi(:,n+1)) / N_phi
           ZBC( n,2) = half_sum + half_difference
           ZBC(-n,2) = half_sum - half_difference
           
           ! RBS:
           half_sum        =  r * r * sum(R2s * cos_n_phi(:,n+1)) / N_phi
           half_difference = -r * r * sum(R2c * sin_n_phi(:,n+1)) / N_phi
           RBS( n,2) = half_sum + half_difference
           RBS(-n,2) = half_sum - half_difference
           
           ! ZBS:
           half_sum        =  r * r * sum(z2s_cylindrical * cos_n_phi(:,n+1)) / N_phi
           half_difference = -r * r * sum(z2c_cylindrical * sin_n_phi(:,n+1)) / N_phi
           ZBS( n,2) = half_sum + half_difference
           ZBS(-n,2) = half_sum - half_difference
        end do
     end if

     if (trim(order_r_option).ne.order_r_option_r1 .and. trim(order_r_option).ne.order_r_option_r2) then
        ! Add O(r^3) terms.
 
        ! Handle the n=0 m=1 modes:
        RBC(0,1) = RBC(0,1) + r * r * r * sum(R3c1) / N_phi
        RBS(0,1) = RBS(0,1) + r * r * r * sum(R3s1) / N_phi
        ZBC(0,1) = ZBC(0,1) + r * r * r * sum(z3c1_cylindrical) / N_phi
        ZBS(0,1) = ZBS(0,1) + r * r * r * sum(z3s1_cylindrical) / N_phi
        
        ! Handle the m=1 modes that have nonzero n:
        do n = 1, ntor
           ! RBC:
           half_sum        =  r * r * r * sum(R3c1 * cos_n_phi(:,n+1)) / N_phi
           half_difference =  r * r * r * sum(R3s1 * sin_n_phi(:,n+1)) / N_phi
           RBC( n,1) = RBC( n,1) + half_sum + half_difference
           RBC(-n,1) = RBC(-n,1) + half_sum - half_difference
           
           ! ZBC:
           half_sum        =  r * r * r * sum(z3c1_cylindrical * cos_n_phi(:,n+1)) / N_phi
           half_difference =  r * r * r * sum(z3s1_cylindrical * sin_n_phi(:,n+1)) / N_phi
           ZBC( n,1) = ZBC( n,1) + half_sum + half_difference
           ZBC(-n,1) = ZBC(-n,1) + half_sum - half_difference
           
           ! RBS:
           half_sum        =  r * r * r * sum(R3s1 * cos_n_phi(:,n+1)) / N_phi
           half_difference = -r * r * r * sum(R3c1 * sin_n_phi(:,n+1)) / N_phi
           RBS( n,1) = RBS( n,1) + half_sum + half_difference
           RBS(-n,1) = RBS(-n,1) + half_sum - half_difference
           
           ! ZBS:
           half_sum        =  r * r * r * sum(z3s1_cylindrical * cos_n_phi(:,n+1)) / N_phi
           half_difference = -r * r * r * sum(z3c1_cylindrical * sin_n_phi(:,n+1)) / N_phi
           ZBS( n,1) = ZBS( n,1) + half_sum + half_difference
           ZBS(-n,1) = ZBS(-n,1) + half_sum - half_difference
        end do

        ! Handle the n=0 m=3 modes:
        RBC(0,3) = r * r * r * sum(R3c3) / N_phi
        RBS(0,3) = r * r * r * sum(R3s3) / N_phi
        ZBC(0,3) = r * r * r * sum(z3c3_cylindrical) / N_phi
        ZBS(0,3) = r * r * r * sum(z3s3_cylindrical) / N_phi
        
        ! Handle the m=3 modes that have nonzero n:
        do n = 1, ntor
           ! RBC:
           half_sum        =  r * r * r * sum(R3c3 * cos_n_phi(:,n+1)) / N_phi
           half_difference =  r * r * r * sum(R3s3 * sin_n_phi(:,n+1)) / N_phi
           RBC( n,3) = half_sum + half_difference
           RBC(-n,3) = half_sum - half_difference
           
           ! ZBC:
           half_sum        =  r * r * r * sum(z3c3_cylindrical * cos_n_phi(:,n+1)) / N_phi
           half_difference =  r * r * r * sum(z3s3_cylindrical * sin_n_phi(:,n+1)) / N_phi
           ZBC( n,3) = half_sum + half_difference
           ZBC(-n,3) = half_sum - half_difference
           
           ! RBS:
           half_sum        =  r * r * r * sum(R3s3 * cos_n_phi(:,n+1)) / N_phi
           half_difference = -r * r * r * sum(R3c3 * sin_n_phi(:,n+1)) / N_phi
           RBS( n,3) = half_sum + half_difference
           RBS(-n,3) = half_sum - half_difference
           
           ! ZBS:
           half_sum        =  r * r * r * sum(z3s3_cylindrical * cos_n_phi(:,n+1)) / N_phi
           half_difference = -r * r * r * sum(z3c3_cylindrical * sin_n_phi(:,n+1)) / N_phi
           ZBS( n,3) = half_sum + half_difference
           ZBS(-n,3) = half_sum - half_difference
        end do
     end if

  case default
     print *, "Invalid finite_r_option:",trim(finite_r_option)
     stop
  end select

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  new_vmec_filename = "input." // xtqsc  !11//15/18.(6l7)as in qs_rd_input().

  istat = 0
  call safe_open(iunit,istat,trim(new_vmec_filename),'unknown','formatted')
  if (istat /= 0) then
     print "(a,a)"," Error opening new_vmec_filename: ",trim(new_vmec_filename)
     stop
  end if
!  print "(a,a)"," Writing results to a vmec input file: ",trim(new_vmec_filename)!out-4/9/19(9r4)

  write(iunit,'(a)') "! This &INDATA namelist was generated by quasisymmetry.f90"
  write(iunit,'(a,a)') "! Based on template file ",trim(vmec_template_filename)
  write(iunit,'(a,es22.12)') "! r =",r

  call write_indata_namelist(iunit, istat)
  if (istat /= 0) then
     print "(a,a)"," Error writing &indata namelist to new_vmec_filename: ",trim(new_vmec_filename)
     stop
  end if

  close(iunit)

end subroutine quasisymmetry_write_vmec_input
