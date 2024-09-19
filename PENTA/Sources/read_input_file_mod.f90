!-----------------------------------------------------------------------------
!+ Contains subroutines for reading input files
!---------------------------------------------------------------------------
Module read_input_file_mod

!
! Description:
!   This module contains the subroutines for reading the <U**2> file,
!   the VMEC data file, the plasma profile file, and the DKES coefficient
!   files.
!
! History:
! Version   Date      Comment
! -------   ----      -------
!  1.0     01/07/2009   Original Code.  JL
!  1.1     05/25/2010   Updated for PENTA3. JL 
!  1.2     01/18/2011   Added error check in reading D* files. JL
!  1.3     09/26/2011   Added cmul,efield check in reading D* files, 
!                       <B**2> to D33_Spitzer. JL

Implicit none

Contains

  Subroutine read_beam_file(roa_surf,kord_pprof)
    Use penta_kind_mod
    Use io_unit_spec, Only: iu_beam
    Use pprof_pass, only :  &
      ! Imported scalar variables:
         beam_force
    Use bspline, only : &
      ! Imported routines
         dbsnak,           &             ! Computes spline knots
         dbsint,           &             ! Computes spline coefficients
         dbsval,           &             ! Evaluates spline
         dbsder                          ! Evaluates spline derivative
    
    Implicit None
    Real(rknd),    Intent(in) :: roa_surf
    Integer(iknd), Intent(in) :: kord_pprof       
    ! Local scalars
    Integer(iknd) :: np_prof                ! Number of radial points in file
    Integer(iknd) :: iocheck, jind
    ! Local arrays:
    Character(120) :: beam_file   
    Real(rknd), Allocatable :: roa_prof(:)        ! r/a values in file
    Real(rknd), Allocatable :: beam_prof(:)       ! beam values in file
    Real(rknd), Allocatable :: beam_knot_array(:)   ! spline knot array
    Real(rknd), Allocatable :: spl_beam(:)          ! spline coefficients
    
    beam_file="beam_force.dat"
    Open(UNIT=iu_beam,FILE=beam_file,STATUS='old',IOSTAT=iocheck)
    ! Check for success
    If ( iocheck /= 0 ) Then
      Write(*,'(a20 a40)') 'Error opening file: ', beam_file
      Stop 'Exiting: I/O Error in subroutine read_beam_file'
    Endif
    Read(iu_beam,*) np_prof
    Allocate(roa_prof(np_prof),beam_prof(np_prof))
    ! Loop over surfaces
    Do jind=1,np_prof
      ! Read each line
      Read(iu_beam,*) roa_prof(jind),beam_prof(jind)
    Enddo
    Close(unit=iu_beam)
  
    ! Spline fit profile and evaluate at r/a of surfaces

    Allocate(beam_knot_array(np_prof+kord_pprof),spl_beam(np_prof))
    
    ! Fit profile
    Call dbsnak(np_prof,roa_prof,kord_pprof,beam_knot_array)
    Call dbsint(np_prof,roa_prof,beam_prof,kord_pprof,beam_knot_array,spl_beam)
    
    ! Evaluate spline fit at r/a of the test surface
    beam_force = dbsval(roa_surf,kord_pprof,beam_knot_array,np_prof,spl_beam)

    !deallocate variables
    Deallocate(roa_prof,beam_prof,beam_knot_array,spl_beam)
    
  End Subroutine read_beam_file
  
  !---------------------------------------------------------------------------
  !+ Reads the VMEC data from the profile_data_*** file
  !---------------------------------------------------------------------------
  Subroutine read_vmec_file(js, run_ident)

  ! Description: 
  !  This routine reads the file "profile_data_***" where *** is the 
  !  string run_ident.  The variables are passed using the module
  !  vmec_var_pass for the surface defined by the integer js.
  !
  ! History:
  !  Version   Date      Comment
  !  -------   ----      -------
  !  1.0     01/7/2009   Original Code.  JL
  !  1.1     5/24/2010   Updated to add B0 input for PENTA3. JL 
  ! 
  ! Author(s): J. Lore 7/2009 - 5/24/2010 
  !
  ! Parent module: read_input_file_mod
  !
  ! Declarations:
  !
  ! Modules used:
  Use penta_kind_mod                ! Import rknd, iknd specifications
  Use io_unit_spec, Only: iu_vmec   ! VMEC file i/o unit number
  Use vmec_var_pass, Only: &        ! Used to pass the variables read
    ! Imported scalar variables
    arad, Rmajor, r_surf, roa_surf, chip,  &
    psip, btheta, bzeta, bsq, iota, vol_p, B0

  Implicit none

  ! Subroutine arguments:
  Integer(iknd), Intent(in) :: js         !see above
  Character(60), Intent(in) :: run_ident

  ! Local scalars:
  Integer(iknd) :: js_min, js_max   ! Surface range in file
  Integer(iknd) :: jind             ! Loop index
  Integer(iknd) :: iocheck                ! used to check file open
  
  ! Local arrays:
  Real(rknd),    Allocatable :: &
  r_vmec(:),                    & ! Minor radius values 
  roa_vmec(:),                  & ! Normalized minor radius values
  chip_vmec(:),                 & ! Radial deriv. of pol. flux values 
  psip_vmec(:),                 & ! Radial deriv. of tor. flux values 
  btheta_vmec(:),               & ! Poloidal B values 
  bzeta_vmec(:),                & ! Toroidal B values 
  vp_vmec(:),                   & ! Radial deriv. of Volume values 
  Bsq_vmec(:),                  & ! <B**2> values 
  iota_vmec(:),                 & ! iota/2/pi values 
  B0_vmec(:)                      ! <B> values
  Integer(iknd), Allocatable :: js_vmec(:)    ! Surface index values
  Character(120)              :: vmec_fname    ! Filename string
  Character(60)              :: ch_dum        ! Used to read legend
  Character(60)              :: tb = char(9)  ! Used to read spaces                              

  !- End of header -----------------------------------------------------------

  ! Read VMEC profile data file 
  vmec_fname = "profile_data_" // run_ident
  Open(UNIT=iu_vmec,FILE=vmec_fname,STATUS="old",IOSTAT=iocheck)
  ! Check for success
  If ( iocheck /= 0 ) then
    Write(*,'(2(a20))') 'Error opening file: ', vmec_fname
    Stop 'Exiting: I/O Error in subroutine read_vmec_file'
  Endif
  Read(iu_vmec,*) js_min, js_max
  Allocate(js_vmec(js_max),r_vmec(js_max),roa_vmec(js_max))
  Allocate(chip_vmec(js_max),psip_vmec(js_max),btheta_vmec(js_max))
  Allocate(bzeta_vmec(js_max),vp_vmec(js_max),Bsq_vmec(js_max))
  Allocate(iota_vmec(js_max),B0_vmec(js_max))
  Read(iu_vmec,*) arad, Rmajor
  Read(iu_vmec,'(a10)') ch_dum
  Do jind = js_min,js_max
    Read(iu_vmec,'(i4,10(a1,e15.7))') js_vmec(jind),tb,r_vmec(jind),tb,       &
    roa_vmec(jind),tb,chip_vmec(jind),tb,psip_vmec(jind),tb,btheta_vmec(jind),&
    tb,bzeta_vmec(jind),tb,vp_vmec(jind),tb,Bsq_vmec(jind),tb,iota_vmec(jind),&
    tb,B0_vmec(jind)
  End do
  Close(iu_vmec)

  ! Assign variables for the current surface
  r_surf = r_vmec(js)
  roa_surf = roa_vmec(js)
  chip = chip_vmec(js)
  psip = psip_vmec(js)
  btheta = btheta_vmec(js)
  bzeta = bzeta_vmec(js)
  vol_p = vp_vmec(js)
  Bsq = Bsq_vmec(js)
  iota = iota_vmec(js)
  B0 = B0_vmec(js)

  ! Deallocate variables
  Deallocate(js_vmec, r_vmec, roa_vmec, chip_vmec, psip_vmec)
  Deallocate(btheta_vmec, bzeta_vmec, vp_vmec, Bsq_vmec, iota_vmec, B0_vmec)

  EndSubroutine read_vmec_file


  !---------------------------------------------------------------------------
  !+ Reads the plasma profile data from the plasma_profiles_xxx.dat file
  !---------------------------------------------------------------------------
  Subroutine read_pprof_file(pprof_char,nis,roa_surf,arad,kord_pprof)

  ! Description: 
  !  This subroutine reads the file "plasma_profiles_xxx.dat" where xxx is the
  !  string pprof_char.  The variables are passed using the module
  !  pprof_pass for the surface defined by roa_surf.
  !
  !  Inputs: 
  !     pprof_char: (see above)
  !     nis: number of ion species
  !     roa_surf: r/a for this surface
  !     arad: minor radius (in meters)
  !     kord_pprof: Spline order for fitting profiles
  !
  ! History:
  !  Version   Date      Comment
  !  -------   ----      -------
  !  1.0     01/07/2009  Original Code.  JL
  !  1.1     05/24/2010  Updated for PENTA3. JL
  !  1.2     15/07/2010  Removed unused input js. JL
  ! 
  ! Author(s): J. Lore 7/2009 - 5/24/2010 
  !
  ! Parent module: read_input_file_mod
  !
  ! Declarations:
  !
  ! Modules used:
  Use penta_kind_mod                ! Import rknd, iknd specifications
  Use io_unit_spec, only: iu_pprof  ! plasma profile file i/o unit number
  Use pprof_pass, only :  &
  ! Imported scalar variables:
    ne, Te, dnedr, dTedr, &
  ! Imported array variables (1D): 
    ni, Ti, dnidr, dTidr, &
  ! Imported array variables (2D): 
    ni_prof, Ti_prof
  Use bspline, only : &
  ! Imported routines
    dbsnak,           &             ! Computes spline knots
    dbsint,           &             ! Computes spline coefficients
    dbsval,           &             ! Evaluates spline
    dbsder                          ! Evaluates spline derivative

  Implicit none
  
  ! Subroutine arguments:        
  Character(20), Intent(in) :: pprof_char  ! See above
  Integer(iknd), Intent(in) :: nis
  Real(rknd),    Intent(in) :: roa_surf
  Real(rknd),    Intent(in) :: arad
  Integer(iknd), Intent(in) :: kord_pprof       

  ! Local scalars:
  Integer(iknd) :: np_prof                ! Number of radial points in file
  Integer(iknd) :: jind                   ! Loop index (for VMEC surf)
  Integer(iknd) :: ispec                  ! Loop index (for ion species)
  Integer(iknd) :: tmp_ind                ! ditto
  Integer(iknd) :: iocheck                ! used to check file open

  ! Local arrays:
  Character(60) :: plasma_prof_file   
  Real(rknd), Allocatable :: roa_prof(:)        ! r/a values in file
  Real(rknd), Allocatable :: ne_prof(:)         ! n_e values in file
  Real(rknd), Allocatable :: Te_prof(:)         ! T_e values in file
  Real(rknd), Allocatable :: ion_pprof_data(:)  ! Array used to read ion data
  Real(rknd), Allocatable :: Te_knot_array(:)   ! T_e spline knot array
  Real(rknd), Allocatable :: ne_knot_array(:)   ! n_e spline knot array
  Real(rknd), Allocatable :: Ti_knot_array(:)   ! T_i spline knot array
  Real(rknd), Allocatable :: ni_knot_array(:)   ! n_i spline knot array
  Real(rknd), Allocatable :: spl_Te(:)          ! T_e spline coefficients
  Real(rknd), Allocatable :: spl_ne(:)          ! n_e spline coefficients
  Real(rknd), Allocatable :: spl_Ti(:)          ! T_i spline coefficients
  Real(rknd), Allocatable :: spl_ni(:)          ! n_i spline coefficients

  !- End of header ---------------------------------------------------------

  !
  ! Read r/a, ne, Te, Ti, ni from file
  !
  plasma_prof_file="plasma_profiles_"//Trim(Adjustl(pprof_char))//".dat"
  Open(UNIT=iu_pprof,FILE=plasma_prof_file,STATUS='old',IOSTAT=iocheck)
  ! Check for success
  If ( iocheck /= 0 ) Then
    Write(*,'(a20 a40)') 'Error opening file: ', plasma_prof_file
    Stop 'Exiting: I/O Error in subroutine read_pprof_file'
  Endif

  Read(iu_pprof,*) np_prof
  Allocate(roa_prof(np_prof),ne_prof(np_prof),Te_prof(np_prof))
  Allocate(Ti_prof(np_prof,nis),ni_prof(np_prof,nis))
  Allocate(ion_pprof_data(nis*2))

  ! Loop over surfaces
  Do jind=1,np_prof
    ! Read each line
    Read(iu_pprof,*) roa_prof(jind),ne_prof(jind),Te_prof(jind),ion_pprof_data
    ! Break apart ion data
    Do ispec=1,nis
      tmp_ind=(ispec-1)*2+1
      ni_prof(jind,ispec)=ion_pprof_data(tmp_ind)
      Ti_prof(jind,ispec)=ion_pprof_data(tmp_ind+1)
    Enddo
  Enddo
  Close(unit=iu_pprof)
  
  !
  ! Spline fit ne,Te,Ti,ni profiles and evaluate at r/a of current surface
  !

  Allocate(ne_knot_array(np_prof+kord_pprof),spl_ne(np_prof))
  Allocate(Te_knot_array(np_prof+kord_pprof),spl_Te(np_prof))
  Allocate(ni_knot_array(np_prof+kord_pprof),spl_ni(np_prof))
  Allocate(Ti_knot_array(np_prof+kord_pprof),spl_Ti(np_prof))

  ! Fit electron profiles
  Call dbsnak(np_prof,roa_prof,kord_pprof,ne_knot_array)
  Call dbsnak(np_prof,roa_prof,kord_pprof,Te_knot_array)   
  Call dbsint(np_prof,roa_prof,ne_prof,kord_pprof,ne_knot_array,spl_ne)
  Call dbsint(np_prof,roa_prof,te_prof,kord_pprof,Te_knot_array,spl_Te)
    
  ! Loop over ion species, fit profiles and assign values
  Do ispec=1,nis
    Call dbsnak(np_prof,roa_prof,kord_pprof,ti_knot_array)
    Call dbsnak(np_prof,roa_prof,kord_pprof,ni_knot_array)
    Call dbsint(np_prof,roa_prof,ni_prof(:,ispec),kord_pprof, &
                ni_knot_array,spl_ni)
    Call dbsint(np_prof,roa_prof,Ti_prof(:,ispec),kord_pprof, &
                Ti_knot_array,spl_ti)
    ! Evaluate spline fit at r/a of the test surface
    ni(ispec) = dbsval(roa_surf,kord_pprof,ni_knot_array,np_prof,spl_ni)
    Ti(ispec) = dbsval(roa_surf,kord_pprof,ti_knot_array,np_prof,spl_ti)
    ! Evaluate derivatives (d/dr) at r/a
    dnidr(ispec) = dbsder(1,roa_surf,kord_pprof,              &
                          ni_knot_array,np_prof,spl_ni)/arad
    dTidr(ispec) = dbsder(1,roa_surf,kord_pprof,              &
                          Ti_knot_array,np_prof,spl_Ti)/arad
  Enddo
 
  ! Evaluate spline fit at r/a of the test surface for electrons
  ne = dbsval(roa_surf,kord_pprof,ne_knot_array,np_prof,spl_ne)
  Te = dbsval(roa_surf,kord_pprof,te_knot_array,np_prof,spl_Te)
  ! Evaluate derivatives (d/dr) at r/a for electrons
  dnedr = dbsder(1,roa_surf,kord_pprof,ne_knot_array,np_prof,spl_ne)/arad
  dTedr = dbsder(1,roa_surf,kord_pprof,Te_knot_array,np_prof,spl_Te)/arad

  !convert units to mks
  ne=ne*1.e18_rknd
  ni=ni*1.e18_rknd
  dnedr=dnedr*1.e18_rknd
  dnidr=dnidr*1.e18_rknd

  !deallocate variables
  Deallocate(roa_prof,ne_prof,Te_prof,ion_pprof_data,Te_knot_array)
  Deallocate(Ti_prof,ni_prof)
  Deallocate(spl_Te,ne_knot_array,spl_ne,ni_knot_array,spl_ni)
  Deallocate(Ti_knot_array, spl_Ti)

  End subroutine read_pprof_file



  !---------------------------------------------------------------------------
  !+ Reads the DKES coeff. data from the D##_star_**** file
  !---------------------------------------------------------------------------
  Subroutine read_dkes_star_files(coeff_ext,Add_Spitzer_to_D33,Bsq)

  ! Description: 
  !  This subroutine reads the files "D##_star_****" where ## corresponds 
  !  to 11, 13, or 33, *** is the identifier for the configuration and the
  !  the Booz_xform (VMEC) surface index (ex: D11_star_hsx_s20).  
  !  The variables are passed using the module coeff_var_pass
  !
  !  Inputs: 
  !     coeff_ext: [Character] (see above)
  !     Add_Spitzer_to_D33: [Logical] If true the collisional or "Spitzer" 
  !                         contribution to D33* is added.  Note that
  !                         this contribution MUST be included, but this
  !                         option allows it to be included in the file ahead
  !                         of time. (Really this is to allow for the direct 
  !                         output of different versions of DKES).
  !
  ! History:
  !  Version   Date      Comment
  !  -------   ----      -------
  !  1.0     01/07/2009   Original Code.  JL
  !  1.1     05/24/2010   Updated for PENTA3. JL 
  !  1.2     01/18/2011   Check for monotonically increasing cmul, efield.JL
  !  1.3     09/26/2011   Added <B**2> term to D33 Spitzer and 
  !                       check for efield, cmul arrays. JL
  ! 
  ! Author(s): J. Lore 7/2009 - 09/26/2011
  !
  ! Parent module: read_input_file_mod
  !
  ! Declarations:
  !
  ! Modules used:
  Use penta_kind_mod                ! Import rknd, iknd specifications
  Use io_unit_spec, Only: iu_coeff  ! DKES coeff file i/o unit number
  Use coeff_var_pass, Only :                 &
    ! Imported arrays
    cmul, efield,                            & ! (1D)
    D11_mat, D13_mat, D31_mat, D33_mat,      & ! coefficient arrays (2D)
    ! Imported scalars
    num_c, num_e

  Implicit none

  ! Subroutine arguments:
  Character(60), Intent(in) :: coeff_ext     ! See above
  Logical, Intent(in) :: Add_Spitzer_to_D33  
  Real(rknd), Intent(in) :: Bsq
  ! Local scalars:
  Integer(iknd) :: ifile       ! File name index
  Integer(iknd) :: ic          ! cmul (nu/v) loop index
  Integer(iknd) :: je          ! efield (Er/v) loop index
  Integer(iknd) :: nc,ne       ! number of cmul (efield) values
  Integer(iknd) :: iocheck     ! used to check file open
  Real(rknd)    :: coeff_tmp   ! Temporary var. for reading file
  Real(rknd)    :: D33_Spitzer ! Used when adding collisional contribution
  
  ! Local arrays:
  Real(rknd), Allocatable :: cmul_vec(:)    ! cmul vals read in
  Real(rknd), Allocatable :: efield_vec(:)  ! efield vals read in
  Real(rknd), Allocatable :: coef2d(:,:)    ! temporary coefficient array
  Real(rknd), allocatable :: efield_D11(:),   & ! efield arrays for coeffs
    efield_D13(:),          &
    efield_D33(:)
  Real(rknd), allocatable :: cmul_D11(:),   & ! efield arrays for coeffs
    cmul_D13(:),          &
    cmul_D33(:)
  Character(120) :: fname      ! File name
  Character(3)  :: fchar      ! String denoting file type

  !- End of header -----------------------------------------------------------

  ! Loop over file type (11,13,33)
  Do ifile=1,3

    ! Define file name
    If (ifile == 1)     then  !D11*
      fchar = "D11"
    Elseif (ifile == 2) then  !D13*
      fchar = "D13"
    Elseif (ifile == 3) then  !D33*
      fchar = "D33"
    Endif
    fname=fchar//"_star_"//coeff_ext
    ! Read file
    Open(UNIT=iu_coeff,FILE=fname,STATUS='old',IOSTAT=iocheck)
    ! Check for success
    If ( iocheck /= 0 ) then
      Write(*,'(a20 a40)') 'Error opening file: ', fname
      Stop 'Exiting: I/O Error in subroutine read_dkes_star_files'
    Endif
    Read (iu_coeff, *) nc, ne     !number of cmul and efield vals 
    Allocate(cmul_vec(nc),efield_vec(ne),coef2d(nc,ne))

    ! Read cmul values
    Read (iu_coeff, *) cmul_vec(1)
    Do ic = 2, nc 
      Read (iu_coeff, *) cmul_vec(ic)
      If ( cmul_vec(ic) .le. cmul_vec(ic-1) ) Then
        Write(*,*) 'Cmul array must be monotonically increasing in D* files'
        Stop 'Exiting: Error in subroutine read_dkes_star_files'
      Endif
    Enddo
    
    ! Read efield values  
    Read (iu_coeff, *) efield_vec(1)     
    Do je = 2, ne 
      Read (iu_coeff, *) efield_vec(je)     
      If ( efield_vec(je) .le. efield_vec(je-1) ) Then
        Write(*,*) 'Efield array must be monotonically increasing in D* files'
        Stop 'Exiting: Error in subroutine read_dkes_star_files'
      Endif
    Enddo

    ! Loop over efield and cmul and read coefficients into 2-D array
    Do je = 1, ne
      Do ic = 1, nc
        Read (iu_coeff, *) coeff_tmp
        coef2d(ic,je) = coeff_tmp
      Enddo
    Enddo
    Close(unit=iu_coeff) !close input file

    ! Set the output variables
    If (ifile == 1) then          !D11*
      Allocate(cmul_D11(nc))
      Allocate(efield_D11(ne))
      Allocate(D11_mat(nc,ne))
      cmul_D11   = cmul_vec
      efield_D11 = efield_vec
      D11_mat = coef2d
    Elseif (ifile == 2) then      !D13*
      Allocate(cmul_D13(nc))
      Allocate(efield_D13(ne))
      Allocate(D13_mat(nc,ne))
      cmul_D13   = cmul_vec
      efield_D13 = efield_vec
      D13_mat = coef2d
    Elseif (ifile == 3) then      !D33*
      Allocate(cmul_D33(nc))
      Allocate(efield_D33(ne))
      Allocate(D33_mat(nc,ne))
      cmul_D33   = cmul_vec
      efield_D33 = efield_vec
      D33_mat = coef2d
    Endif

    Deallocate(cmul_vec,efield_vec,coef2d)
  Enddo   !ifile loop

  ! Verify the efield and cmul arrays are identical
  Do ic = 1, nc 
    If (( cmul_D11(ic) .ne. cmul_D13(ic) ) .or. ( cmul_D11(ic) .ne. cmul_D33(ic) )) Then
      Write(6,*) 'Error: cmul arrays must be equal for all D## coefficients!'
      Stop
    Endif
  Enddo
  Do je = 1, ne
    If (( efield_D11(je) .ne. efield_D13(je) ) .or. ( efield_D11(je) .ne. efield_D33(je) )) Then
      Write(6,*) 'Error: efield arrays must be equal for all D## coefficients!'
      Stop
    Endif
  Enddo

  Allocate(cmul(nc),efield(ne))
  cmul = cmul_D11
  efield = efield_D11

  Deallocate(cmul_D11,cmul_D13,cmul_D33)
  Deallocate(efield_D11,efield_D13,efield_D33)
   
  ! Get the size of the arrays
  num_c = Size(cmul)
  num_e = Size(efield)

  ! Create D31 using Onsager symmetry
  Allocate(D31_mat(num_c,num_e))
  D31_mat    = - D13_mat

  ! Optionally add collisional "Sptizer" contribution to D33*
  If ( Add_Spitzer_to_D33 ) Then
    Do je = 1, num_e
      Do ic = 1, num_c
        
        ! Calculate D33* Spitzer = (2/3)*<B**2>/cmul
        D33_Spitzer = (2._rknd/3._rknd)*Bsq/cmul(ic)

        ! Calculate D33*(Physical)=D33*(Spitzer)-D33*
        D33_mat(ic,je) = D33_Spitzer - D33_mat(ic,je)
      Enddo
    Enddo
  Endif

  End subroutine read_dkes_star_files


  !---------------------------------------------------------------------------
  !+ Reads the <U**2> data from the Utilde2_profile file
  !---------------------------------------------------------------------------
  Subroutine read_Utilde2_file(roa_test,U2,kord_pprof)

  ! Description: 
  !  This subroutine reads the file "Utilde2_profile".
  !  The value at the test surface is returned.
  !  See the file penta.f90 for information on how this file is to be 
  !  formatted.
  !
  !  Inputs: 
  !     roa_test: r/a for the current surface
  !     kord_pprof: Spline order for fitting <U**2> vs r/a
  !  Outputs:
  !     U2: Value of <U**2> at roa_test
  !
  ! History:
  !  Version   Date      Comment
  !  -------   ----      -------
  !  1.0     01/7/2009   Original Code.  JL
  !  1.1     5/25/2010   Updated for PENTA3. JL 
  ! 
  ! Author(s): J. Lore 7/2009 - 5/25/2010 
  !
  ! Parent module: read_input_file_mod
  !
  ! Declarations:
  !
  ! Modules used:
  Use penta_kind_mod                ! Import rknd, iknd specifications
  Use io_unit_spec, Only: iu_Ufile  ! file i/o unit number
  Use bspline, Only : &
  ! Imported routines
    dbsnak,           &             ! Computes spline knots
    dbsint,           &             ! Computes spline coefficients
    dbsval                          ! Evaluates spline

  Implicit none

  ! Subroutine arguments:
  Real(rknd), Intent(in)    :: roa_test   ! r/a of current surface
  Real(rknd), Intent(out)   :: U2         ! <U**2> at this r/a
  Integer(iknd), Intent(in) :: kord_pprof ! Spline order for fitting

  ! Local scalars:
  Integer(iknd) :: iocheck    ! used to check file open
  Integer(iknd) :: num_pts    ! Number of r/a and U2 points specified
  Integer(iknd) :: jind       ! loop index
  
  ! Local arrays:
  Real(rknd), Allocatable  :: roa_Utilde(:)    ! r/a points in file
  Real(rknd), Allocatable  :: U_tilde2(:)      ! <U**2> points in file
  Real(rknd), Allocatable  :: U2_knot_array(:) ! spline knot array
  Real(rknd), Allocatable  :: spl_U2(:)        ! spline coefficients
  Character(60) :: Ufile_name                  ! File name

  !- End of header -----------------------------------------------------------

  ! Open file and get number of defined points
  Ufile_name="Utilde2_profile"
  Open(UNIT=iu_Ufile,FILE=Ufile_name,STATUS='old',IOSTAT=iocheck)

  ! Check for success
  If ( iocheck /= 0 ) Then
    Write(*,'(a20 a40)') 'Error opening file: ', Ufile_name
    Stop 'Exiting: I/O Error in subroutine read_Utilde2_file'
  Endif

  Read(iu_Ufile,*) num_pts

  ! Allocate memory
  Allocate(roa_Utilde(num_pts))
  Allocate(U_tilde2(num_pts))
    
  ! Read file and assign variables
  Do jind = 1,num_pts
    Read(iu_Ufile,*) roa_Utilde(jind),U_tilde2(jind)
  Enddo
  Close(iu_Ufile)  

  ! Spline fit profile and evaluate at r/a of current surface  
  Allocate(U2_knot_array(num_pts+kord_pprof))
  Allocate(spl_U2(num_pts))
  ! Calculate spline knots
  Call dbsnak(num_pts,roa_Utilde,kord_pprof,U2_knot_array)
  ! Calculate spline coefficients
  Call dbsint(num_pts,roa_Utilde,U_tilde2,kord_pprof,U2_knot_array,spl_U2)
  ! Evaluate spline at test r/a
  U2 = dbsval(roa_test,kord_pprof,U2_knot_array,num_pts,spl_U2)

  ! Deallocate variables
  Deallocate(roa_Utilde,U_tilde2,U2_knot_array,spl_U2)
  EndSubroutine read_Utilde2_file

End module read_input_file_mod
!- End of module header --------------------------------------------------
