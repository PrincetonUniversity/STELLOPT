cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Subroutines for reading input files 
c   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      module read_input_file_mod
      use penta_kind_mod
      use vmec_var_pass
      use pprof_pass
      use bspline
      implicit none
    
      contains
c
c----------------------------------------------------------------------------
c   Reads the file "Utilde2_profile" and returns <U**2> for the run surface
c----------------------------------------------------------------------------
c
      subroutine read_Utilde2_file(roa_test,U2)
      use penta_kind_mod
      use io_unit_spec
      implicit none
      !dummy variables
      real(rknd) :: U2, roa_test
      !local varibles
      integer(iknd) :: j, js_surf, num_pts, kord_prof
      real(rknd), dimension(:), allocatable :: roa_Utilde, U_tilde2, 
     1  U2_knot_array,spl_U2
      character(60) :: Ufile_name

      !open file and get number of defined points
      Ufile_name="Utilde2_profile"
      open(unit=iu_Ufile,file=Ufile_name,status="old")
      read(iu_Ufile,*) num_pts
      allocate(roa_Utilde(num_pts))
      allocate(U_tilde2(num_pts))
    
      !read file and assign variables
      do j = 1,num_pts
        read(iu_Ufile,*) roa_Utilde(j),U_tilde2(j)
      end do
      close(iu_Ufile)  !close file

      !spline fit profile and evaluate at r/a of current surface
      kord_prof = 3 !spline order
      allocate(U2_knot_array(num_pts+kord_prof))
      allocate(spl_U2(num_pts))
      call dbsnak(num_pts,roa_Utilde,kord_prof,U2_knot_array)
      call dbsint(num_pts,roa_Utilde,U_tilde2,kord_prof,
     1            U2_knot_array,spl_U2)
      U2 = dbsval(roa_test,kord_prof,U2_knot_array,num_pts,
     1    spl_U2)

      !deallocate variables
      deallocate(roa_Utilde,U_tilde2,U2_knot_array,spl_U2)
      end subroutine read_Utilde2_file
c
c----------------------------------------------------------------------------
c   Reads the file "profile_data_***" and assigns VMEC variables, passed
c     using module vmec_var_pass
c----------------------------------------------------------------------------
c
      subroutine read_vmec_file(js,run_ident)
      use vmec_var_pass
      use penta_kind_mod
      use io_unit_spec
      implicit none
      !dummy variables
      integer(iknd) :: js
      character(60) :: run_ident
      !local variables
      integer(iknd) :: j, js_min, js_max
	integer(iknd), dimension(:), allocatable :: js_vmec
      real(rknd), dimension(:), allocatable :: r_vmec, roa_vmec
     1 ,chip_vmec, psip_vmec, btheta_vmec, bzeta_vmec, vp_vmec, bsq_vmec
     2 ,iota_vmec
      character(60) :: vmec_fname, ch_dum, tb = char(9)

      !Read VMEC profile data file 
      vmec_fname = "profile_data_" // run_ident
      open(unit=iu_vmec,file=vmec_fname,status="old")
      read(iu_vmec,*) js_min, js_max
      allocate(js_vmec(js_max),r_vmec(js_max),roa_vmec(js_max))
      allocate(chip_vmec(js_max),psip_vmec(js_max),btheta_vmec(js_max))
      allocate(bzeta_vmec(js_max),vp_vmec(js_max),bsq_vmec(js_max))
      allocate(iota_vmec(js_max))
      read(iu_vmec,*) arad, Rmajor
      read(iu_vmec,'(a10)') ch_dum
      do j = js_min,js_max
        read(iu_vmec,'(i4,9(a1,e15.7))') js_vmec(j),tb,r_vmec(j),tb,
     1  roa_vmec(j),tb,chip_vmec(j),tb,psip_vmec(j),tb,btheta_vmec(j),
     2    tb,bzeta_vmec(j),tb,vp_vmec(j),tb,bsq_vmec(j),tb,iota_vmec(j)
      end do
      close(iu_vmec)

      ! Assign variables for the current surface
      chip = chip_vmec(js); psip = psip_vmec(js); bsq = bsq_vmec(js)
      btheta = btheta_vmec(js); bzeta = bzeta_vmec(js)
      iota = iota_vmec(js); r_surf = r_vmec(js); roa_surf = roa_vmec(js)
      vol_p = vp_vmec(js)
      !deallocate variables
      deallocate(js_vmec, r_vmec, roa_vmec, chip_vmec, psip_vmec)
      deallocate(btheta_vmec, bzeta_vmec, vp_vmec, bsq_vmec,iota_vmec)
      end subroutine read_vmec_file
c
c----------------------------------------------------------------------------
c   Reads the wout file for vmec information (SAL 05/30/13)
c----------------------------------------------------------------------------
c
      subroutine read_vmec_file_2(js,run_ident)
      use vmec_var_pass
      use penta_kind_mod
      use io_unit_spec
      use read_wout_mod, ONLY:  read_wout_file, phipf_vmec => phipf,
     1                          btheta_vmec => buco, bzeta_vmec => bvco,
     2                          iota_vmec => iotaf, vp_vmec => vp,
     3                          Rmajor_vmec => Rmajor,
     4                          Aminor_vmec => Aminor,
     5                          phi_vmec => phi, bmnc, mnmax, gmnc,
     6                          xm, xn, ns_vmec => ns
      implicit none
      !dummy variables
      integer(iknd) :: js
      character(60) :: run_ident
      !local variables
      integer(iknd) :: j, js_min, js_max, ierr, u, v, mn
	integer(iknd), dimension(:), allocatable :: js_vmec
!      real(rknd), dimension(:), allocatable :: r_vmec, roa_vmec
!     1 ,chip_vmec, psip_vmec, btheta_vmec, bzeta_vmec, vp_vmec, bsq_vmec
!     2 ,iota_vmec
      REAL(rknd) :: TWOPI, top, bottom, theta, zeta, arg, cs_arg,
     1              jacob_vmec, bbf
      character(60) :: vmec_fname, ch_dum, tb = char(9)

      TWOPI = 8*ATAN(1.0_rknd)
      !Read VMEC profile data file 
!      vmec_fname = "profile_data_" // run_ident
!      open(unit=iu_vmec,file=vmec_fname,status="old")
      CALL read_wout_file(TRIM(run_ident),ierr)
      
!      read(iu_vmec,*) js_min, js_max
!      allocate(js_vmec(js_max),r_vmec(js_max),roa_vmec(js_max))
!      allocate(chip_vmec(js_max),psip_vmec(js_max),btheta_vmec(js_max))
!      allocate(bzeta_vmec(js_max),vp_vmec(js_max),bsq_vmec(js_max))
!      allocate(iota_vmec(js_max))
!      read(iu_vmec,*) arad, Rmajor
!      read(iu_vmec,'(a10)') ch_dum
!      do j = js_min,js_max
!        read(iu_vmec,'(i4,9(a1,e15.7))') js_vmec(j),tb,r_vmec(j),tb,
!     1  roa_vmec(j),tb,chip_vmec(j),tb,psip_vmec(j),tb,btheta_vmec(j),
!     2    tb,bzeta_vmec(j),tb,vp_vmec(j),tb,bsq_vmec(j),tb,iota_vmec(j)
!      end do
!      close(iu_vmec)

      ! Assign variables for the current surface
      chip = iota_vmec(js)*phipf_vmec(js);  !note this only works for non-RFP
      psip = phipf_vmec(js)  
      ! bsq = bsq_vmec(js)
      btheta = btheta_vmec(js); bzeta = bzeta_vmec(js)
      iota = iota_vmec(js);
      roa_surf = sqrt(phi_vmec(js)/phi_vmec(ns_vmec))
      r_surf   = Aminor_vmec*roa_surf
      vol_p = vp_vmec(js)
      Rmajor = Rmajor_vmec
      arad = Aminor_vmec
      ! Now calc Bsq
      bsq = 0.0
      DO v = 1, 360
         zeta = TWOPI*REAL(v-1)/359.
         DO u = 1, 360
            theta = TWOPI*REAL(u-1)/359.
            bbf = 0.0; jacob_vmec = 0.0;
            DO mn = 1, mnmax
               arg = xm(mn)*theta - xn(mn)*zeta
               cs_arg = cos(arg)
               jacob_vmec = jacob_vmec + gmnc(mn,js)*cs_arg
               bbf = bbf + bmnc(mn,js)*cs_arg
            END DO
            top = top + jacob_vmec*bbf*bbf
            bottom = bottom + jacob_vmec
         END DO
      END DO
      bsq = top/bottom
      !deallocate variables
!      deallocate(js_vmec, r_vmec, roa_vmec, chip_vmec, psip_vmec)
!      deallocate(btheta_vmec, bzeta_vmec, vp_vmec, bsq_vmec,iota_vmec)
      end subroutine read_vmec_file_2
c
c----------------------------------------------------------------------------
c   Reads the file "plasma_profiles*.dat" and assigns plasma parameters for the
c     current surface.  Variables are passed using module pprof_pass
c----------------------------------------------------------------------------
c
      subroutine read_pprof_file(js,pprof_char,nis,roa_surf,arad)
      use penta_kind_mod
      use pprof_pass
      use bspline
      use io_unit_spec
      implicit none
      !dummy variables
      real(rknd) :: roa_surf, arad
      integer(iknd) :: js, nis
      character*10 :: pprof_char
      !local variables
      integer(iknd) :: np_prof, j, ispec, tmp_ind, kord_prof
      real(rknd), dimension(:), allocatable :: roa_prof, ne_prof, 
     1 Te_prof, ion_pprof_data, Te_knot_array, spl_Te,  ne_knot_array,
     2 spl_ne, ni_knot_array, spl_ni, Ti_knot_array, spl_Ti
      character(60) :: plasma_prof_file   

      !read r/a, ne,Te,Ti,ni from file
      plasma_prof_file="plasma_profiles"
     1  //trim(adjustl(pprof_char))//".dat"
      open(unit=iu_pprof,file=plasma_prof_file,status="old")
      read(iu_pprof,*) np_prof
      allocate(roa_prof(np_prof),ne_prof(np_prof),Te_prof(np_prof))
      allocate(Ti_prof(np_prof,nis),ni_prof(np_prof,nis))
      allocate(ion_pprof_data(nis*2_iknd))
      do j=1,np_prof
        read(iu_pprof,*) roa_prof(j),ne_prof(j),Te_prof(j),
     1             ion_pprof_data
        do ispec=1,nis
          tmp_ind=(ispec-1)*2+1
          Ti_prof(j,ispec)=ion_pprof_data(tmp_ind)
          ni_prof(j,ispec)=ion_pprof_data(tmp_ind+1)
        enddo
      enddo
      close(unit=iu_pprof)

      !spline fit ne,Te,Ti,ni profiles and evaluate at r/a of current surface
      kord_prof = 3 !spline order for profile fitting
      allocate(ne_knot_array(np_prof+kord_prof),spl_ne(np_prof))
      allocate(Te_knot_array(np_prof+kord_prof),spl_Te(np_prof))
      allocate(ni_knot_array(np_prof+kord_prof),spl_ni(np_prof))
      allocate(Ti_knot_array(np_prof+kord_prof),spl_Ti(np_prof))

      !fit electron profiles
      call dbsnak(np_prof,roa_prof,kord_prof,ne_knot_array)
      call dbsnak(np_prof,roa_prof,kord_prof,Te_knot_array)   
      call dbsint(np_prof,roa_prof,ne_prof,kord_prof,
     1            ne_knot_array,spl_ne)
      call dbsint(np_prof,roa_prof,te_prof,kord_prof,
     1            Te_knot_array,spl_Te)
    
      !loop over ion species, fit profiles and assign values
      do ispec=1,nis
        call dbsnak(np_prof,roa_prof,kord_prof,ti_knot_array)
        call dbsnak(np_prof,roa_prof,kord_prof,ni_knot_array)
        call dbsint(np_prof,roa_prof,ni_prof(:,ispec),kord_prof,
     1            ni_knot_array,spl_ni)
        call dbsint(np_prof,roa_prof,Ti_prof(:,ispec),kord_prof,
     1            Ti_knot_array,spl_ti)
        !evaluate spline fit at r/a of the test surface
        ni(ispec) = dbsval(roa_surf,kord_prof,ni_knot_array,np_prof,
     1    spl_ni)
        Ti(ispec) = dbsval(roa_surf,kord_prof,ti_knot_array,np_prof,
     1    spl_ti)
        !evaluate derivatives (d/dr) at r/a
        dnidr(ispec) = dbsder(1,roa_surf,kord_prof,ni_knot_array,
     1                   np_prof,spl_ni)/arad
        dTidr(ispec) = dbsder(1,roa_surf,kord_prof,Ti_knot_array,
     1                   np_prof,spl_Ti)/arad
      enddo
 
      !evaluate spline fit at r/a of the test surface for electrons
      ne = dbsval(roa_surf,kord_prof,ne_knot_array,np_prof,spl_ne)
      Te = dbsval(roa_surf,kord_prof,te_knot_array,np_prof,spl_te)
      !evaluate derivatives (d/dr) at r/a for electrons
      dnedr = dbsder(1,roa_surf,kord_prof,ne_knot_array,
     1                   np_prof,spl_ne)/arad
      dTedr = dbsder(1,roa_surf,kord_prof,Te_knot_array,
     1                   np_prof,spl_Te)/arad
      !convert units to mks
      ne=ne*1.e18_rknd
      ni=ni*1.e18_rknd
      dnedr=dnedr*1.e18_rknd
      dnidr=dnidr*1.e18_rknd

      !deallocate variables
      deallocate(roa_prof,ne_prof,Te_prof,ion_pprof_data,Te_knot_array)
      deallocate(spl_Te,ne_knot_array,spl_ne,ni_knot_array,spl_ni)
      deallocate(Ti_knot_array, spl_Ti)
      end subroutine read_pprof_file
c
c----------------------------------------------------------------------------
c   Reads the file *star_lijs_*** files and assigns variables for passing
c     using module coeff_var_pass
c----------------------------------------------------------------------------
c
      subroutine read_lmn_star_files(coeff_ext)
      use penta_kind_mod
      use coeff_var_pass
      use io_unit_spec
      implicit none
      !dummy variables
      character(60) :: coeff_ext
      !local variables
      real(rknd) :: junk, coeff_tmp
      integer(iknd) :: ifile, nc, ne, i, j
      real(rknd), dimension(:), allocatable :: cmul_vec, efield_vec
      real(rknd), dimension(:,:), allocatable :: coef2d
      character(60) :: fname
      character(1) :: fchar

      !loop over file type (L,M,N)
      do ifile=1,3

        !define file name
        if (ifile == 1) then      !L
          fchar = "l"
        elseif (ifile == 2) then  !M
          fchar = "m"
        elseif (ifile == 3) then  !N
          fchar = "n"
        endif
        fname=fchar//"star_lijs_"//coeff_ext

        !read file
         open(unit=iu_coeff,file=fname,status='old')
        read (iu_coeff, *) junk,junk    !wtov, dsdr no longer needed    
        read (iu_coeff, *) nc, ne       !number of cmul and efield vals 
        allocate(cmul_vec(nc),efield_vec(ne),coef2d(nc,ne))
        do i = 1, nc 
          read (iu_coeff, *) cmul_vec(i)
        enddo
        do j = 1, ne 
          read (iu_coeff, *) efield_vec(j)       
        enddo

        !Loop over efield and cmul and read coefficients into 2-D array
        do j = 1, ne
          do i = 1, nc
            read (iu_coeff, *) coeff_tmp
            if (ifile .eq. 2) then
            coef2d(i,j) = dlog(coeff_tmp)
            else
              coef2d(i,j) = coeff_tmp
            endif
          enddo
        enddo
        close(unit=iu_coeff) !close input file

        !set the output variables
        if (ifile == 1) then          !L
          allocate(cmul_ls(nc))
          allocate(efield_ls(ne))
          allocate(coef2d_ls(nc,ne))
          cmul_ls=cmul_vec
          efield_ls=efield_vec
          coef2d_ls=coef2d
        elseif (ifile == 2) then      !M
          allocate(cmul_ms(nc))
          allocate(efield_ms(ne))
          allocate(coef2d_ms(nc,ne))
          cmul_ms=cmul_vec
          efield_ms=efield_vec
          coef2d_ms=coef2d
        elseif (ifile == 3) then      !N
          allocate(cmul_ns(nc))
          allocate(efield_ns(ne))
          allocate(coef2d_ns(nc,ne))
          cmul_ns=cmul_vec
          efield_ns=efield_vec
          coef2d_ns=coef2d
        endif

        deallocate(cmul_vec,efield_vec,coef2d)
      enddo   !ifile loop


      !size of arrays
      num_cl=size(cmul_ls); num_el=size(efield_ls)
      num_cm=size(cmul_ms); num_em=size(efield_ms)
      num_cn=size(cmul_ns); num_en=size(efield_ns)
      end subroutine read_lmn_star_files

      end module read_input_file_mod