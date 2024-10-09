!-----------------------------------------------------------------------------
!
!   Miscellaneous subroutines 
!   
!-----------------------------------------------------------------------------
Module PENTA_subroutines
  Implicit None
Contains
  

!-----------------------------------------------------------------------------
!+ Finds the ambipolar Er roots from Sum(gamma_s * q_s) = 0
!-----------------------------------------------------------------------------
Subroutine find_Er_roots(gamma_e,gamma_i,Er_test_vals,Z_ion, &
  num_Er_test,num_ion_species,Er_roots,num_roots)
!
! Description: 
!  This subroutine QQ
!   This subroutine determines the ambipolar roots as gamma_e=sum(Z*gammai)
!     where the sum is over ion species.  The zeros are found from linear
!     interpolation and using the functions zeroin and rad_flux.  
!   The roots (er_roots(:)) and number of roots (num_roots) are passed 
!!   using the module Er_roots_pass.
!   If an even number of roots is found, the first root only is chosen.
!
! Function arguments:
!
!  QQ update this -- talk about method
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     07/01/2009  Original Code.  JL
!  1.1     09/01/2010  Adapted for PENTA3. JL
!  1.2     09/28/2011  Fixed bug that prevented choosing a single
!                      root in the case of even number of roots. JL
!
! Author(s): J. Lore 2007 - 09/28/2011
!
!
! Modules used:
Use penta_kind_mod                ! Import rknd, iknd specifications
Use penta_functions_mod, Only : &
! Imported functions
zeroin,   &     ! Finds zero of a function over a range
rad_flux        ! Returns function to zero at Er

Implicit None

! Input/output                      !See above for descriptions
Real(rknd),    Intent(in)  :: gamma_e(num_Er_test)
Real(rknd),    Intent(in)  :: gamma_i(num_Er_test,num_ion_species)
Real(rknd),    Intent(in)  :: Er_test_vals(num_Er_test)
Real(rknd),    Intent(in)  :: Z_ion(num_ion_species)
Integer(iknd), Intent(in)  :: num_Er_test
Integer(iknd), Intent(in)  :: num_ion_species
Real(rknd),    Intent(out) :: Er_roots(:) 
Integer(iknd), Intent(out) :: num_roots 

! Local scalars
Real(rknd)    :: a_test,b_test  ! Used for checking for zero crossing
Integer(iknd) :: ispec,ie,iroot ! Loop indices
Integer(iknd) :: ind            ! Current root index
Integer(iknd) :: root_diff      ! Used to check for multiple roots in fit range
! Local arrays (1D)
Real(rknd)    :: qg_all_i(num_Er_test)       ! Gamma_i * Z_i
Real(rknd)    :: I_tmp(num_Er_test)          ! Indices of sign changes

! Allocatable arrays
Real(rknd), Allocatable :: I_roots(:)   ! Root indices
Real(rknd), Allocatable :: Er_fit(:)    ! Er values used to fit
Real(rknd), Allocatable :: diff_qg(:)   ! Points of function to zero (Gamma_e -
                                        !  Zi*Gamma_i), which are fit
! Parameters for root finding
Real(rknd), Parameter    :: tol_er  = 1.e-10_rknd ! Tolerance on zeroin function
Integer(iknd), Parameter :: num_pts = 2_iknd      ! Number of points around root
                                                  !  used in fit, i.e., order of 
                                                  !  polynomial fit + 1.

!- End of header -------------------------------------------------------------

! Sum ion fluxes times charge numbers
qg_all_i = 0._rknd
Do ispec = 1_iknd, num_ion_species
  qg_all_i = qg_all_i + Z_ion(ispec)*gamma_i(:,ispec)
EndDo

! Find where the sign switches QQ root changes at end of array?
num_roots = 0_iknd
Do ie = 1_iknd, num_Er_test-1_iknd
  a_test = Dsign(1._rknd,gamma_e(ie)-qg_all_i(ie))
  b_test = Dsign(1._rknd,gamma_e(ie+1)-qg_all_i(ie+1))
  If ( Int(a_test) /= Int(b_test) ) Then
    num_roots = num_roots + 1_iknd
    I_tmp(num_roots) = ie
  EndIf
EndDo

! Check for zero or even number of roots
If ( num_roots == 0_iknd ) Then
  Write(*,*) 'No roots found in search range'
  Stop 'Error from find_Er_roots: Please modify the Er search range'
ElseIf ( Mod(num_roots,2_iknd) == 0_iknd) Then
  Write(*,*) 'Even number of roots found, choosing first root only.'
  num_roots = 1_iknd
Endif

! Allocate variables by number of roots
Allocate(I_roots(num_roots))

! Initial guess of the root indices
I_roots = I_tmp(1:num_roots)

! Find zeros using interpolation (linear if num_pts = 2)

! Allocate fit variables
Allocate(diff_qg(num_pts))
Allocate(Er_fit(num_pts))

! Loop over roots and interpolate to find fine root
Do iroot=1,num_roots
  ! This index
  ind = I_roots(iroot)

  ! Make sure there are not multiple roots in num_pts range
  If ( (num_roots >  1_iknd) .AND. (iroot /= num_roots) ) Then 
    root_diff = I_roots(iroot+1) - ind
    If (root_diff <=  num_pts-1) Then
      Write(*,*) 'Distance between roots too small for ', &
        'interpolation.  Increase num_Er_test.'
      Write(*,*) I_roots
      Stop 'Exiting from find_Er_roots'
    EndIf
  EndIf

  ! Define difference in flux and Er over fit range
  diff_qg = gamma_e(ind:ind+1) - qg_all_i(ind:ind+1)
  Er_fit  = Er_test_vals(ind:ind+1)
    
  ! Fit the values over range and use zeroin and rad_flux to find a root
  Er_roots(iroot) = zeroin(Er_fit(1),Er_fit(num_pts),rad_flux,tol_er, &
   diff_qg,Er_fit,num_pts)

EndDo

! Deallocate variables
Deallocate(I_roots,diff_qg,Er_fit)

EndSubroutine find_Er_roots

!-----------------------------------------------------------------------------
!+ Calculates the classical friction coefficients (lmat)
!-----------------------------------------------------------------------------
Subroutine define_friction_coeffs(masses,charges,v_ths,Temps,dens, &
                                   loglambda,num_species,Smax,lmat)
! Description: 
!  This subroutine calculates the classical friction coefficients using the 
!  general formulae from [Ji and Held, PoP 13, 102103, (2006)] 
!  and the PhD thesis of J. Lore.
!  Expressions for low orders can also be found in [Hirshman & Sigmar, '81]
!  [Helander and Sigmar] etc.
!
!  Inputs: 
!     masses:    [Real array] Array of masses for each plasma species starting
!                              with electrons and then each ion species (kg)
!     charges:   [Real array] Same as masses except charges (C)
!     v_ths:     [Real array] Same as above except thermal velocities (m/s)
!     Temps:     [Real array] Same as above except temperatures (eV)
!     dens:      [Real array] Same as above except densities (m^-3)
!     loglambda:       [Real] Coulomb logarithm (assumed same for all spec.)
!     num_species:  [Integer] Total number of plasma species
!     Smax:         [Integer] Order of Sonine polynomial expansion
!
!   Outputs:
!     lmat(:,:,:,:): [Real array]  Array of classical friction coefficients
!                     First two indices are the plasma species (ex: l_ei)
!                     Second two indices are the order (ex: l_ee^1,1)
!  QQ update this
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     01/7/2009   Original Code.  JL
!  1.1     5/25/2010   Updated for PENTA3. JL 
! 
! Author(s): J. Lore 7/2009 - 5/25/2010 
!
!
! Declarations:
!
! Modules used:
Use penta_kind_mod                ! Import rknd, iknd specifications
Use phys_const, Only :                 &
  ! Imported Parameters
    eps0,                              &  ! Electric constant
    pi                 
Use penta_functions_mod, Only :        &
  ! Imported functions
  calc_Mak_Ji,                         &
  calc_Nab_Ji

Implicit None

! Subroutine arguments:                   !See above for descriptions
Integer(iknd), Intent(in)   :: num_species  
Integer(iknd), Intent(in)   :: Smax  
Real(rknd),    Intent(in)   :: loglambda
Real(rknd),    Intent(in)   :: masses(num_species)
Real(rknd),    Intent(in)   :: charges(num_species)
Real(rknd),    Intent(in)   :: v_ths(num_species)
Real(rknd),    Intent(in)   :: Temps(num_species)
Real(rknd),    Intent(in)   :: dens(num_species)
Real(rknd),    Intent(out)  :: lmat((Smax+1)*num_species,(Smax+1)*num_species)

! Local scalars:
Real(rknd)    :: tau_coeff        ! Coefficient in front of collision times
Real(rknd)    ::                & ! Primary species (a) paramaters 
  ma, Ta, vth_a, charge_a, na     
Real(rknd)    ::                & ! Secondary species (b) paramaters 
  mb, Tb, vth_b, charge_b, nb     
Real(rknd)    ::                & ! Dimensionless ratios of species a on k
  chi_ab, theta_ab, mu_ab, Q_ab
Real(rknd)    :: m_ovr_tau_sum    ! Running sum of M element over tau
Real(rknd)    ::                & ! Parameters for sum over species (spec. k)
  mk, Tk, vth_k, charge_k, nk,  &
  tau_ak 
Real(rknd)    ::                & ! Dimensionless ratios of species a on k
  chi_ak, theta_ak, mu_ak, Q_ak
Real(rknd)    :: Mak              ! M matrix element of species a on k
Real(rknd)    :: Nab              ! N matrix element of species a on b
Real(rknd)    :: tau_ab           ! Collision time of a on b
Integer(iknd) :: ispec1, ispec2   ! Species loop indices
Integer(iknd) :: irow, icol       ! Submatrix loop indices
Integer(iknd) :: spec_k           ! Sum over species loop index
Integer(iknd) :: ind1, ind2       ! Indices for full lmat


! Local arrays
Real(rknd) :: lab(Smax+1,Smax+1)  ! Submatrix of coefficients

!- End of header -------------------------------------------------------------

! Coefficient for collision times
tau_coeff = 3._rknd * Sqrt(pi) * pi * eps0**2_iknd

!
! Define main friction coefficient matrix (lmat)
!

! Loop over species a (row)
Do ispec1 = 1, num_species

  ! Assign parameters for species a
  ma       = masses(ispec1)
  Ta       = Temps(ispec1)
  vth_a    = v_ths(ispec1)
  charge_a = charges(ispec1)
  na       = dens(ispec1)
  
  ! Loop over second species (columns)
  Do ispec2 = 1, num_species

    ! Loop over row and column of submatrix lab
    Do irow = 0, Smax
      Do icol = 0, Smax
        
        !
        ! Define submatrix of friction coeffs. for species Sonine order
        !

        ! Initialize current element
        lab(irow+1, icol+1) = 0._rknd

        ! Diagonal terms (Otherwise M sum is zero)
        If ( ispec1 == ispec2 ) Then

          M_ovr_tau_sum = 0._rknd

          ! Sum over all species (indicated "k")
          Do spec_k = 1, num_species
            
            ! Assign parameters for species k
            vth_k    = v_ths(spec_k)
            Tk       = Temps(spec_k)
            mk       = masses(spec_k)
            nk       = dens(spec_k)
            charge_k = charges(spec_k)

            ! Calculate the collision time for species a on k for this row, col
            tau_ak = tau_coeff * ma**2_iknd * vth_a**3_iknd   &
                     / ( nk * charge_a**2_iknd * charge_k**2_iknd * loglambda )

            ! Define dimensionless ratios
            chi_ak   = vth_k/vth_a
            theta_ak = Tk/Ta
            mu_ak    = mk/ma
            Q_ak     = 1._rknd + chi_ak**2_iknd

            ! Calculate Mak (Bragniskii matrix element) for this row, col
            Mak = calc_Mak_Ji(irow, icol, theta_ak, mu_ak, chi_ak, Q_ak)

            ! Add this Mak to the sum
            M_ovr_tau_sum = M_ovr_tau_sum + Mak / tau_ak

          Enddo  ! species k loop

          ! For diagonal elements we have a M term and a N term (added below)
          lab(irow+1,icol+1) = na * ma * M_ovr_tau_sum

        Endif  ! Diagonal element if

        ! Define parameters for species b
        vth_b     = v_ths(ispec2)
        Tb        = Temps(ispec2)
        mb        = masses(ispec2)
        charge_b  = charges(ispec2)
        nb        = dens(ispec2)

        ! Define dimensionless ratios
        chi_ab   = vth_b / vth_a 
        theta_ab = Tb / Ta
        mu_ab    = mb / ma
        Q_ab     = 1._rknd + chi_ab**2_iknd

        ! All elements of lab have the N term
        tau_ab = tau_coeff * ma ** 2_iknd * vth_a**3_iknd & 
                 /( nb * charge_a**2_iknd * charge_b**2_iknd * loglambda )
        Nab = calc_Nab_Ji(irow, icol, theta_ab, mu_ab, chi_ab, Q_ab)

        lab(irow+1,icol+1)=lab(irow+1,icol+1)+na*ma*( Nab/tau_ab )

      Enddo ! icol loop
    Enddo ! irow loop

    ! Define large friction matrix indices

    ind1 = ( ispec1 - 1 ) * ( Smax + 1 ) + 1
    ind2 = ( ispec2 - 1 ) * ( Smax + 1 ) + 1
    lmat( ind1:(ind1 + Smax) ,ind2:(ind2 + Smax) ) = lab

  Enddo !species b loop
Enddo  !species a loop
  
End Subroutine define_friction_coeffs

!-----------------------------------------------------------------------------
!+ Calculates the thermodynamic force vector (Xvec)
!-----------------------------------------------------------------------------
Subroutine form_Xvec(Er_test,Z_ion,B_Eprl,nis,Xvec)
! Description: 
!  This subroutine calculates the thermodynamic force vector (Xvec).
!  See [Sugama, Nishimura] for definitions of Xa1,Xa2,XE.  The Xvec is
!  arranged as [Xe1,Xe2,X_i1_1,X_i1_2,X_i2_1...XE]
!  Inputs: 
!     Er_test: [Real]       Er [V/m]
!     Z_ion:   [Real array] Ion charge numbers (corresponding to i1,i2...)
!     B_Eprl:  [Real]       <B dot E_||>
!     nis:      [Integer]    Number of ion species
!
!   Outputs:
!     Xvec:   [Real array]  Thermodynamic force vector (see above)
!
!  QQ update this
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     01/07/2009  Original Code.  JL
!  1.1     07/15/2010  Updated for PENTA3. JL 
! 
! Author(s): J. Lore 7/2009 - 7/15/2010 
!
!
! Declarations:
!
! Modules used:
Use penta_kind_mod                ! Import rknd, iknd specifications
Use phys_const, Only :                 &
  ! Imported Parameters
  elem_charge                     ! unit charge [C]
Use pprof_pass                    ! Import plasma profile information
Use vmec_var_pass                 ! Import geometry info (VMEC)

Implicit none

! Subroutine arguments
Integer(iknd), intent(in)  :: nis                 ! See above for definitions
Real(rknd),    intent(in)  :: Er_test
Real(rknd),    intent(in)  :: B_Eprl
Real(rknd),    intent(in)  :: Z_ion(nis)
Real(rknd),    intent(out) :: Xvec(nis*2+3)

!- End of header -------------------------------------------------------------
  
! Note: Species charge sign has been included, elem_charge has been
!  factored out of Er term because Ta is in [eV].

! Electrons
Xvec(1) = -Te*elem_charge*(dnedr/ne + dTedr/Te + Er_test/Te)
Xvec(2) = -elem_charge*dTedr

! Ions
Xvec(3:(nis-1)*2+3:2) = -elem_charge*Ti*(dnidr/ni + dTidr/Ti - Z_ion*Er_test/Ti)
Xvec(4:(nis-1)*2+4:2) = -elem_charge*dTidr 

! E_|| term
Xvec(nis*2+3) = B_Eprl/Sqrt(Bsq)

End Subroutine form_Xvec

!-----------------------------------------------------------------------------
!+  -This subroutine calculates the spline fit coefficients to the D##* data
!-----------------------------------------------------------------------------
Subroutine fit_coeffs(cmul,efield,num_c,num_e,Dstar,log_interp,kcord,keord, &
 xt_c,xt_e,Dspl,cmin,cmax,emin,emax )
! Description: 
!  This subroutine calculates the spline knots and coefficients to the 
!   "normalized" monoenergetic coefficient data D11*, D31* and D33*.
!
!  QQ update this
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     01/07/2009  Original Code.  JL
!  1.1     07/19/2010  Updated for PENTA3. JL 
!  1.2     01/18/2011  Added check for very small efield. JL
! 
! Author(s): J. Lore 7/2009 - 01/18/2011
!
!
! Declarations:
!
! Modules used:
Use penta_kind_mod                ! Import rknd, iknd specifications
Use phys_const, Only :                 &
  ! Imported Parameters
  elem_charge                     ! unit charge [C]
Use pprof_pass                    ! Import plasma profile information
Use bspline, Only :   &
    dbsnak,           &             ! Computes spline knots
    dbs2in                          ! Computes spline coefficients
Implicit none

! Local parameters
Real(rknd), Parameter :: esmall = 1.e-20_rknd ! Substituted for efield=0
                                              ! if log interp is performed.
                                              ! This should be smaller than
                                              ! min(efield).

! Subroutine arguments
Real(rknd),    Intent(in)  :: cmul(num_c)     ! See above for definitions
Real(rknd),    Intent(in)  :: efield(num_e)
Integer(iknd), Intent(in)  :: num_c
Integer(iknd), Intent(in)  :: num_e
Real(rknd),    Intent(in)  :: Dstar(num_c,num_e)
Logical,       Intent(in)  :: log_interp
Integer(iknd), Intent(in)  :: kcord
Integer(iknd), Intent(in)  :: keord
Real(rknd),    Intent(out) :: xt_c(num_c + kcord)
Real(rknd),    Intent(out) :: xt_e(num_e + keord)
Real(rknd),    Intent(out) :: Dspl(num_c,num_e)
Real(rknd),    Intent(out) :: cmin
Real(rknd),    Intent(out) :: cmax
Real(rknd),    Intent(out) :: emin
Real(rknd),    Intent(out) :: emax

! Local arrays
Real(rknd)   :: ctmp(num_c)
Real(rknd)   :: etmp(num_e)
Real(rknd)   :: enrm(num_e)

!- End of header -------------------------------------------------------------

ctmp = cmul
etmp = efield

! Take log. of x-y arrays if necessary
!  Also check for log(0)
If ( log_interp .EQV. .true. ) Then
  Where ( efield == 0._rknd ) etmp = esmall
  ! Check for very small efield
  If ( Minval(etmp(2:num_e)) .le. esmall ) Then
    Write(*,*) 'Minval(efield) is smaller than parameter esmall'
    Write(*,*) ' after efield = 0. substitution. Decrease esmall'
    Write(*,*) ' or exclude such small efield values'
    Write(*,*) 'Minval(efield), esmall =',Minval(etmp),esmall
    Stop 'Error: Exiting from subroutine fit_coeffs'
  Endif
  ctmp = Log10(ctmp)
  etmp = Log10(etmp)
EndIf

! Define limits of cmul, efield arrays
cmin = Minval(ctmp)
cmax = Maxval(ctmp)
emin = Minval(etmp)
emax = Maxval(etmp)

! Define normalized array for efield for spline fitting
enrm = ( etmp - emin ) / ( emax - emin ) 

! Calculate spline knots and coefficients
Call dbsnak(num_c,ctmp,kcord,xt_c)    !compute 'not-a-knot' sequence
Call dbsnak(num_e,enrm,keord,xt_e)    !compute 'not-a-knot' sequence
Call dbs2in(num_c,ctmp,num_e,enrm,Dstar,num_c,kcord,keord,xt_c,xt_e,Dspl)

EndSubroutine fit_coeffs


End Module PENTA_Subroutines
