cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   -PENTA with Impurities-
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  To run,type: PENTA_imp [command line arguments]
c
c  where the command line arguments are:
c
c  1  unique identifier text on DKES database files
c      i.e., for (m,l,n)star_lijs_***, this would be
c      the text that replaces *** for the specific surface (ex: hsx_s10)
c  2  Er,min
c  3  Er,max  - args 2 and 3 specify the interval in Er over which
c      the search is done for the ambipolar electric field root.
c      Er,min and Er,max are dimensionless normalized electric
c      field variables, i.e.: e<a>Er,min/max/kT_e
c    - Note if the logical variable input_is_Er below is true then 
c      the search interval is set by args 2 and 3 as Er (V/cm)
c  4  flux surface number
c  5  parameter that determines whether output data should be
c      written into new files (=0) or appended to
c      existing files (=1). Typically,when a script is run for a
c      sequence of flux surfaces, this is set to 0 for the first
c      surface and to 1 for subsequent surfaces
c  6  extension of the profile_data_*** file (i.e., the *** text)
c      that comes from extracting the appropriate profiles from
c      the VMEC run using the profile_extractor.f code
c  7  unique identifier on plasma profiles file, i.e., the ***
c      text in the filename plasma_profiles***.dat. This allows
c      multiple profiles to be run without having to rename files
c  8  <B*E||> in T*V/m
c
c
c  Required input files:
c
c   [l,m,n]star_lijs_***_s# - where *** is arg1 above and # arg 4.
c       These files contain the values of efield and cmul and the
c       corresponding normalized monoenergetic coefficients.  Note
c       that L* should be L*/e**2 where e is the elementary charge.
c
c   ion_params - A file containing the namelist "ion_params" which
c       defines the variables num_ion_species, Z_ion_init and 
c       miomp_init which are the number of ion species (integer),
c       the corresponding ion charge numbers (real, array) and 
c       the ion to proton mass ratio (real, array) for each 
c       non-electron species.
c
c   profile_data_*** - where *** is arg1 above.  Contains geometry
c       variables, most likely from a VMEC run.
c
c   plasma_profilesXXX.dat - A file containing the plasma profile
c       data, where XXX is arg8 above.  This file has the format:
c       row 1: number of radial points for which data is specified.
c       all other rows: r/a, ne, Te, Ti1, ni1, Ti2, ni2 ...
c       where the densities are in units of 10**12/cc and the 
c       temperatures in eV and i1, i2... are the ion species.
c
c	Utilde2_profile - Contains the quantity <U**2>, where U is the
c		PS flow function as defined in Sugama and Nishimura.  The
c		first row is the number of points, then r/a points and the
c		corresponding <U**2> value.
c   
c  Output files:
c
c   "fluxes_vs_roa", "fluxes_vs_Er", "flows_vs_roa", "plasma_profiles_check"
c
c   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   NOTES
c    - Be sure to check parameters below
c    - based on Sugama, Nishimura PoP 9 4637 (2002) and 
c         Sugama, Nishimura PoP 15, 042502 (2008)
c    - Some code used from orignal PENTA code, written by Don Spong and 
c        modified by Jeremy Lore
c
c  7/2009 JL
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   PROGRAMMING NOTES
c   -speed optimization
c   -recursive ambipolar solver
c     -fix numargs
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    
      program penta_imp
      use penta_kind_mod
      use phys_const
      use read_input_file_mod
      use vmec_var_pass
      use coeff_var_pass
      use pprof_pass
      use Er_roots_pass
      use io_unit_spec
      use parameter_pass
      implicit none

      logical,       parameter ::  input_is_Er = .true.       !If true, Er range is (V/cm) otherwise e<a>Er/kT_e
      integer(iknd), parameter ::  num_ion_max = 20_iknd    !Maximum number of ion species to be read from namelist file
      integer(iknd), parameter ::  num_Er_test = 200_iknd    !Number of Er points in search range
      logical,    parameter ::  use_quanc8_set = .false.    !If false, rectangular approx. to convolution used
      real(rknd), parameter ::        Kmin_set = 1.e-4_rknd !Minimum K in energy convolution
      real(rknd), parameter ::        Kmax_set = 10._rknd   !Maximum K in energy convolution
      real(rknd), parameter ::      epsabs_set = 1.e-8_rknd !Absolute tolerance for quanc8  (used if use_quanc8_set=.true.)
      real(rknd), parameter ::      epsrel_set = 1.e-6_rknd !Relative tolerance for quanc8  (used if use_quanc8_set=.true.)
      integer(iknd),parameter :: numKsteps_set = 1000_iknd  !Number of K points (linear) for convolution (used if use_quanc8_set=.false.)
      logical, parameter ::     log_interp_set = .true.     !If true, logarithmic interpolation of DKES coefficients is used

      integer(iknd) :: num_species, iroot, js, i_append, 
     1  num_ion_species, inv_err, nn_inv, j, ispec, tmp_ind, ie,
     2 ispec1, ispec2, spec1_ind, spec2_ind, numargs

      real(rknd) :: vth_e, loglambda, sigma_S, Er_test, U2,
     1 lX_sum1, lX_sum2, J_E_tot, Er_min, Er_max, B_Eprl, J_E_cl, 
     2 eaEr_o_kTe

      real(rknd), dimension(:), allocatable :: Z_ion,ion_mass,vth_i,
     1  Xvec, gamma_e_ambi, q_e_ambi, B_uprl, B_qprl, u_theta, u_zeta,
     2  B_uprle, B_qprle, u_thetae, J_E_tot_parts, J_bs_ambi, L1, L2, 
     3  L3, M1, M2, M3, N1, N2, N3, charges, gamma_PS, q_PS, dens, 
     4  masses, Temps, vths

      real(rknd), dimension(:,:), allocatable :: lmat, force_mat,
     1 flux_mat, flux_mat_inv, FFmat, gamma_i, q_i, gamma_i_ambi,
     2 q_i_ambi, B_uprli, B_qprli, u_thetai

      real(rknd), dimension(num_Er_test) :: gamma_e, q_e, Er_test_vals,
     1 qg_all_i
      real(rknd), dimension(num_ion_max) :: Z_ion_init,miomp_init  

      character(1)  :: tb = char(9)
      character(10) :: pprof_char
      character(100) :: arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, 
     1  arg9, coeff_ext, run_ident, format_tmp, str_num, fstatus, fpos

      namelist / ion_params / num_ion_species,Z_ion_init,miomp_init

      !external functions
      integer(iknd), external :: iargc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   MAIN PROGRAM
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      !pass parameters to necessary functions
      use_quanc8=use_quanc8_set
      Kmin=Kmin_set
      Kmax=Kmax_set
      epsabs=epsabs_set
      epsrel=epsrel_set
      num_K_steps=numKsteps_set
      log_interp=log_interp_set

      !default values
      B_Eprl=0._rknd
     
      !Read namelist file for ion parameters
      open(iu_nl,file="ion_params",status="old",form="formatted")
      read(iu_nl,nml=ion_params)
      close(iu_nl)

      !Get command line arguments
      !numargs = iargc()
	numargs = 8
      call getarg(1, arg1)
      call getarg(2, arg2)
      call getarg(3, arg3)
      call getarg(4, arg4)
      call getarg(5, arg5)
      call getarg(6, arg6)
      call getarg(7, arg7)
      if (numargs .eq. 8) call getarg(8, arg8)

      !store command line args
      coeff_ext = trim(adjustl(arg1))
      read(arg2,*) Er_min
      read(arg3,*) Er_max
      read(arg4,*) js
      read(arg5,*) i_append
      run_ident = trim(adjustl(arg6))
      pprof_char = trim(adjustl(arg7))
      if (numargs .eq. 8) read(arg8,*) B_Eprl
  
      !allocate variables by number of species defined
      num_species=num_ion_species+1_iknd
      allocate(Z_ion(num_ion_species),ion_mass(num_ion_species))
      allocate(ni(num_ion_species),Ti(num_ion_species))
      allocate(dnidr(num_ion_species),dTidr(num_ion_species))
      allocate(vth_i(num_ion_species),J_E_tot_parts(num_species))
      allocate(charges(num_species),dens(num_species),vths(num_species))
      allocate(masses(num_species),Temps(num_species))
      allocate(q_PS(num_species),gamma_PS(num_species))
      allocate(L1(num_species),L2(num_species),L3(num_species))
      allocate(M1(num_species),M2(num_species),M3(num_species))
      allocate(N1(num_species),N2(num_species),N3(num_species))
      allocate(Xvec(num_species*2+1),lmat(num_species*2,num_species*2))
      allocate(force_mat(num_species*2,num_species*2+1))
      allocate(flux_mat(num_species*2,num_species*2))
      allocate(flux_mat_inv(num_species*2,num_species*2))
      allocate(FFmat(num_species*2,num_species*2+1))
      allocate(gamma_i(num_ion_species,num_Er_test))
      allocate(q_i(num_ion_species,num_Er_test))
      allocate(B_uprl(num_species),B_qprl(num_species))
      allocate(u_theta(num_species),u_zeta(num_species))

      !read input files
      call read_vmec_file(js,run_ident)
      call read_pprof_file(js,pprof_char,num_ion_species,roa_surf,arad)
      call read_lmn_star_files(coeff_ext)
      call read_Utilde2_file(roa_surf,U2)

      !fit lmn coefficients
      call fit_coeffs

      !change Er test range to V/cm if necessary
      if ( input_is_Er .eqv. .false.)  then
        Er_min=Er_min*Te/arad
        Er_max=Er_max*Te/arad
      endif

      !assign ion parameters
      Z_ion=Z_ion_init(1:num_ion_species)
      ion_mass=miomp_init(1:num_ion_species)*p_mass

      !Display run data to screen (if i_append==0)
      if (i_append .eq. 0) then
        write(*,*) 
        write(*,*) "Welcome to PENTA with impurities, ",
     1    "please note the following settings:"
        write(*,*)
        write(*,'(a,i3)') ' Number of ion species: ',num_ion_species
        if ( input_is_Er .eqv. .true. ) then
          write(*,*) 'Interpreting input range as Er (V/cm)'
        else
          write(*,*) 'Interpreting input range as e<a>Er/kTe'
        endif
        if ( log_interp .eqv. .true. ) then
          write(*,*) 'Performing logarithmic interpolation in Er,cmul'
        else
          write(*,*) 'Performing linear interpolation in Er,cmul'
        endif
        if ( use_quanc8 .eqv. .true. ) then
          write(*,'(a,2(a,e11.4))') 
     1     ' Using quanc8 integrator with tolerances: ',
     2      'abs: ',epsabs,' rel: ', epsrel
        else
          write(*,'(a,i6,a,a)') ' Using ',num_K_steps,' point integral',
     1       ' approximation'
        endif
        write(*,'(a,2(" ",e15.4))') ' K range on convolution integral: '
     1      ,Kmin, Kmax
        write(*,*)
        write(*,*) " <r>/<a>","   Er roots(V/cm)"
      endif

      !set specifiers for opening output files
      if(i_append .eq. 0) then
        fstatus="unknown"
        fpos="asis"
      elseif(i_append .eq. 1) then
        fstatus="old"
        fpos="append"
      endif

      !open output files
      open(unit=iu_flux_out, file="fluxes_vs_roa",
     1     position=trim(adjustl(fpos)),status=trim(adjustl(fstatus)))
      open(unit=iu_pprof_out, file="plasma_profiles_check",
     1     position=trim(adjustl(fpos)),status=trim(adjustl(fstatus)))
      open(unit=iu_fvEr_out, file="fluxes_vs_Er",
     1     position=trim(adjustl(fpos)),status=trim(adjustl(fstatus)))
      open(unit=iu_flows_out, file="flows_vs_roa",
     1     position=trim(adjustl(fpos)),status=trim(adjustl(fstatus)))

      !write legends if (i_append = 0) for most files
      if(i_append .eq. 0) then
        !fluxes vs r/a
        write(iu_flux_out,'("*",/,"r/a",a1,"Er(V/cm)",a1,"e<a>Er/kTe",  
     1    a1,"gamma_e",a1,"q_e",a1,"gamma_i",a1,"q_i",a1,"Jbs")')
     2    tb, tb, tb, tb, tb, tb, tb
        !plasma profile check
        write(iu_pprof_out,'("*",/,"r/a Te ne dnedr dTedr ",
     1      "Ti ni dnidr dTidr")')
        !flows vs r/a
        write(iu_flows_out,'("*",/,"r/a Er(V/cm) e<a>Er/kTe <B*u_||e>",
     1    " <B*u_||i> <q_||e> <q_||i> <u_pol_e> <u_pol_i> <u_tor_e> ",
     2   "<u_tor_i>")')
      endif

      !legend for fluxes vs Er is written for each surface
      write(iu_fvEr_out,'("*",/,"r/a Er(V/cm) gamma_e gamma_i",
     1    " q_e q_i")')

      !Calculate thermal velocities and loglambda
      vth_i = dsqrt(2._rknd*Ti*elem_charge/ion_mass)
      vth_e = dsqrt(2._rknd*Te*elem_charge/e_mass)
      if(Te .gt. 50._rknd) then
        loglambda = 25.3_rknd - 1.15_rknd*dlog10(ne/1.e6_rknd) + 
     1    2.3_rknd*dlog10(Te)
      else
       loglambda = 23.4_rknd - 1.15_rknd*dlog10(ne/1.e6_rknd) + 
     1    3.45_rknd*dlog10(Te)
      endif
    
      !assign arrays of parameters for all species
      charges=elem_charge*[-1._rknd,Z_ion]
      dens=[ne, ni]
      masses=[e_mass, ion_mass]
      Temps=[Te,Ti]
      vths=[vth_e,vth_i]

      !define matrix of friction coefficients (lmat)
      call define_friction_coeffs(masses,charges,vths,Temps,dens,
     1   loglambda,num_species,lmat)

      !calculate Spitzer conductivity
      sigma_S=ne**2_iknd*elem_charge**2_iknd *  
     1   lmat(2,2)/( lmat(1,1)*lmat(2,2) + (lmat(1,2))**2_iknd)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Loop over Er and evaluate fluxes(Er)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !Loop over Er to get fluxes as a function of Er
      do ie=1,num_Er_test
        !define test radial electric field in V/cm
        Er_test=(Er_min+(real(ie,rknd)-1._rknd)
     1    *(Er_max-Er_min)/(num_Er_test-1._rknd))*100._rknd
        !check for Er=0, doesn't work for log interpolation
        if (( Er_test .eq. 0._rknd ) .and. (log_interp)) then
          Er_test=(Er_min+real(ie,rknd)
     1      *(Er_max-Er_min)/(num_Er_test-1._rknd))*50._rknd  !use mean of 0 and next point
          write(*,*) 'Cannot use Er=0 with log_interp, new Er=',Er_test
        endif
        Er_test_vals(ie)=Er_test
      
        !Form thermodynamic force vector (Xvec)
        call form_Xvec(Er_test,Z_ion,B_Eprl,num_ion_species,Xvec)

        !define thermal coefficients (L1,L2,L3,M1...)
        call define_thermal_mats(Er_test,num_ion_species,vths,dens,
     1    masses,charges,loglambda,L1,L2,L3,M1,M2,M3,N1,N2,N3)

        !define flux and force matrices (flux_mat,force_mat)
        call define_flux_force_mat(num_species,lmat,
     1    L1,L2,L3,M1,M2,M3,N1,N2,N3,charges,dens,flux_mat, force_mat)   

        !invert flux mat
        nn_inv=size(flux_mat,1)
        call FINDInv(flux_mat, flux_mat_inv, nn_inv, inv_err)

        !calculate flux-force matrix
        FFmat=matmul(flux_mat_inv,force_mat)

        !calculate PS fluxes
        call calculate_PS_flux(charges,num_species,lmat,Xvec,U2,
     1    gamma_PS,q_PS)

        !calculate heat and particle fluxes
        gamma_e(ie)=sum(FFmat(1,:)*Xvec)+gamma_PS(1)            !electrons
        q_e(ie)=(sum(FFmat(2,:)*Xvec)+q_PS(1))*elem_charge*Te
        do ispec=1,num_ion_species                              !ions
          tmp_ind=3+2*(ispec-1)
          gamma_i(ispec,ie)=sum(FFmat(tmp_ind,:)*Xvec)+gamma_PS(ispec+1)
          q_i(ispec,ie)=(sum(FFmat(tmp_ind+1,:)*Xvec)+q_PS(ispec+1))*
     1        elem_charge*Ti(ispec)
        enddo


        !write fluxes vs Er
        write(str_num,*) num_ion_species*2+3
	!write(str_num,'(i)') num_ion_species*2+3
        format_tmp='(f7.4,' // trim(adjustl(str_num)) // '(" ",e15.7))'
        write(iu_fvEr_out,format_tmp) roa_surf,Er_test/100._rknd,
     1      gamma_e(ie),gamma_i(:,ie),q_e(ie),q_i(:,ie)


      enddo !efield loop
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Find roots and evaluate quantities at ambipolar roots
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !calculate classical parallel electric current
      J_E_cl=sigma_S*Xvec(num_species*2+1)

      !find the ambipolar root(s) from gamma_e = sum(Z*gamma_i)
      call find_Er_roots(gamma_e,gamma_i,Er_test_vals,Z_ion,
     1  num_Er_test,num_ion_species,qg_all_i)

      !allocate variables by number of ambipolar roots
      allocate(gamma_i_ambi(num_ion_species,num_roots))
      allocate(q_i_ambi(num_ion_species,num_roots))
      allocate(B_uprli(num_ion_species,num_roots))
      allocate(B_qprli(num_ion_species,num_roots))
      allocate(u_thetai(num_ion_species,num_roots))
      allocate(gamma_e_ambi(num_roots),q_e_ambi(num_roots))
      allocate(B_qprle(num_roots),B_uprle(num_roots))
      allocate(u_thetae(num_roots),J_bs_ambi(num_roots))

      !evaluate fluxes at ambipolar Er(s), and calculate flows
      do iroot=1,num_roots  !loop over ambipolar roots
    
        !test Er
        Er_test=Er_roots(iroot)
        eaEr_o_kTe=arad*Er_test/Te

        !Form thermodynamic force vector (Xvec)
        call form_Xvec(Er_test,Z_ion,B_Eprl,num_ion_species,Xvec)

        !define thermal coefficients (L1,L2,L3,M1...)
        call define_thermal_mats(Er_test,num_ion_species,vths,dens,
     1    masses,charges,loglambda,L1,L2,L3,M1,M2,M3,N1,N2,N3)

        !define flux and force matrices (flux_mat, force_mat)
        call define_flux_force_mat(num_species,lmat,
     1    L1,L2,L3,M1,M2,M3,N1,N2,N3,charges, dens,flux_mat, force_mat)  

        !invert flux mat
        nn_inv=size(flux_mat,1)
        call FINDInv(flux_mat, flux_mat_inv, nn_inv, inv_err)

        !calculate flux-force matrix
        FFmat=matmul(flux_mat_inv,force_mat)

        !calculate PS fluxes
        call calculate_PS_flux(charges,num_species,lmat,Xvec,U2,
     1    gamma_PS,q_PS)

        !calculate ambipolar fluxes
        gamma_e_ambi(iroot)=sum(FFmat(1,:)*Xvec)+gamma_PS(1)          !electron
        q_e_ambi(iroot)=(sum(FFmat(2,:)*Xvec)+q_PS(1))*elem_charge*Te
        do ispec=1,num_ion_species                                    !ions
          tmp_ind=3+2*(ispec-1)
          gamma_i_ambi(ispec,iroot)=sum(FFmat(tmp_ind,:)*Xvec)
     1      +gamma_PS(ispec+1)
          q_i_ambi(ispec,iroot)= (sum(FFmat(tmp_ind+1,:)*Xvec)+
     1      q_PS(ispec+1))*elem_charge*Ti(ispec)
        enddo

        !now calculate flows (<B*u_||>,<B*q_||>,<u_theta>,<u_zeta>)
        call calculate_flows(num_species,L1,L2,L3,N1,N2,N3,Xvec,
     1    gamma_e_ambi(iroot),gamma_i_ambi(:,iroot),q_e_ambi(iroot),
     2    q_i_ambi(:,iroot),charges,B_uprl,B_qprl,u_theta,u_zeta)

        !assign flows by species
        B_uprle(iroot)=B_uprl(1)
        B_qprle(iroot)=B_qprl(1)
        B_uprli(:,iroot)=B_uprl(2:num_species)
        B_qprli(:,iroot)=B_qprl(2:num_species)
        u_thetae(iroot)=u_theta(1)
        u_thetai(:,iroot)=u_theta(2:num_species)

        !calculate total parallel electric current
        J_E_tot_parts=dens*charges*B_uprl*dsqrt(Bsq)
        J_E_tot=sum(J_E_tot_parts)

        !calculate Bootstrap current
        J_bs_ambi(iroot)=J_E_tot-J_E_cl

        !write fluxes to file "flux_vs_roa"
        write(str_num,*) 2*num_ion_species+5
	 !write(str_num,'(i)') 2*num_ion_species+5
        format_tmp='(f7.4,' // trim(adjustl(str_num)) // '(" ",e15.7))'
        write(iu_flux_out,format_tmp) roa_surf,Er_test/100._rknd,
     1    eaEr_o_kTe,gamma_e_ambi(iroot),q_e_ambi(iroot),
     2    gamma_i_ambi(:,iroot),q_i_ambi(:,iroot), J_bs_ambi(iroot)

        !write flows to file "flows_vs_roa"
        write(str_num,*) 4*num_ion_species+6
	!write(str_num,'(i)') 4*num_ion_species+6
        format_tmp='(f7.4,' // trim(adjustl(str_num)) // '(" ",e15.7))'
        write(iu_flows_out,format_tmp) roa_surf,Er_test/100._rknd,
     1    eaEr_o_kTe,B_uprle(iroot),B_uprli(:,iroot),B_qprle(iroot),
     2    B_qprli(:,iroot),u_theta,u_zeta
      enddo !iroot loop


      !write output files
      !"plasma_profiles_check"
      !write(str_num,'(i)') 4*num_ion_species+5
	write(str_num,*) 4*num_ion_species+5
      format_tmp='(' // trim(adjustl(str_num)) // '(" ",e15.7))'
      write(iu_pprof_out,format_tmp) roa_surf,Te,ne,dnedr,dTedr,
     1  Ti,ni,dnidr,dTidr

      !write to screen
      !write(str_num,'(i)') num_roots
	write(str_num,*) num_roots
      format_tmp='(f7.3,' // trim(adjustl(str_num)) // '(" ",e15.4))'
      write(*,format_tmp) roa_surf,er_roots/100._rknd

      !deallocate variables
      deallocate(Z_ion,ion_mass,Xvec,lmat,gamma_i,q_i)
      deallocate(ni,Ti,dnidr,dTidr,vth_i)
      deallocate(L1,L2,L3,M1,M2,M3,N1,N2,N3)
      deallocate(flux_mat,force_mat,FFmat,charges)
      deallocate(gamma_PS,q_PS,er_roots,masses,Temps,vths)
      deallocate(gamma_i_ambi,q_i_ambi,gamma_e_ambi,q_e_ambi)
      deallocate(B_uprl,B_qprl,u_theta,u_zeta,J_bs_ambi)
      deallocate(B_uprli,B_qprli,B_uprle,B_qprle,u_thetae,u_thetai)
      deallocate(ni_prof,Ti_prof,cmul_ls,efield_ls,cmul_ms,efield_ms)
      deallocate(cmul_ns,efield_ns,coef2d_ls,coef2d_ms,coef2d_ns)

      !close output files
      close(iu_flux_out)
      close(iu_pprof_out)
      close(iu_fvEr_out)
      close(iu_flows_out)

      end program penta_imp

