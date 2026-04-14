module transport_mod
    implicit none

    ! Integer, parameter :: rknd = selected_real_kind(12,300) 
    ! Integer, parameter :: iknd = selected_int_kind(8) 
    

contains

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

Implicit None

Integer, parameter :: rknd = selected_real_kind(12,300) 
Integer, parameter :: iknd = selected_int_kind(8) 
Real(rknd),parameter :: eps0 = 8.854187817e-12_rknd        !electric const
Real(rknd),parameter :: pi = 3.14159265358979323846_rknd   !pi

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


Subroutine collision_frequency_penta(vparticles,masses,Zcharges,Temps,dens,loglambda,Nvparticles,num_species,nu)

Implicit None

Integer, parameter :: rknd = selected_real_kind(12,300) 
Integer, parameter :: iknd = selected_int_kind(8)
Real(rknd),parameter :: eps0 = 8.854187817e-12_rknd        !electric const
Real(rknd),parameter :: pi = 3.14159265358979323846_rknd   !pi
Real(rknd),parameter :: EC = 1.602176487e-19_rknd !elem. charge 

! Input/output
Integer(iknd), Intent(in) :: num_species
Integer(iknd), Intent(in) :: Nvparticles
Real(rknd),    Intent(in)   :: loglambda
Real(rknd),    Intent(in)   :: masses(num_species)
Real(rknd),    Intent(in)   :: Zcharges(num_species)
Real(rknd),    Intent(in)   :: Temps(num_species)
Real(rknd),    Intent(in)   :: dens(num_species)
Real(rknd),    Intent(in)   :: vparticles(Nvparticles)
Real(rknd),    Intent(out)  :: nu(Nvparticles)

! Local
Integer(iknd) :: j
Real(rknd) :: m1, Z1, vp
Real(rknd) :: x(num_species), aux(num_species), prefactor(num_species)

m1 = masses(1)
Z1 = Zcharges(1)

Do j=1,Nvparticles
    vp = vparticles(j)
    x = vp*vp * masses/(2*EC*Temps)
    aux = (1 - 0.5/x) * erf(sqrt(x)) + exp(-x) / sqrt(x*pi)
    prefactor = Z1**2 * Zcharges**2 * EC**4 * loglambda * dens / (m1**2 * vp**3 * 4*pi * eps0**2)
    nu(j) = SUM(prefactor * aux)
End Do

End Subroutine collision_frequency_penta


!-----------------------------------------------------------------------------
!+ Calculates the Nab Braginskii matrix element using Ji and Held formulae
!-----------------------------------------------------------------------------
Function calc_Nab_Ji(ptest, ktest, theta, mu, chi, Q) &
Result(Nab)
!
! Description: 
!  This function calculates the Mak Braginskii matrix element using the 
!  general formulae from Ji and Held.
!
! Function arguments:
!
!  QQ update this -- talk about method
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     05/28/2010  Original Code.  JL
! 
! Author(s): J. Lore 5/28/2010 
Implicit None

Integer, parameter :: rknd = selected_real_kind(12,300) 
Integer, parameter :: iknd = selected_int_kind(8)

! Input/output                      !See above for descriptions
Integer(iknd), Intent(in)  :: ptest   ! QQ check subscripts
Integer(iknd), Intent(in)  :: ktest
Real(rknd),    Intent(in)  :: theta
Real(rknd),    Intent(in)  :: mu
Real(rknd),    Intent(in)  :: chi
Real(rknd),    Intent(in)  :: Q
Real(rknd)                 :: Nab

! Local scalars
Real(rknd) :: Bpk, Gamma_q, Gamma_m, cpq, ckm
Real(rknd) :: Gamma_p, Gamma_k
Real(rknd) :: Bstar_qm, Btmp
Integer(iknd) :: qtest,mtest

!- End of header -------------------------------------------------------------


Gamma_p = Gamma_aux(ptest + 2.5_rknd)
Gamma_k = Gamma_aux(ktest + 2.5_rknd)

Bpk = 0._rknd
Do qtest = 0_iknd, ptest
  Do mtest = 0_iknd, ktest
        
    ! Calculate c coefficients
    Gamma_q = Gamma_aux(qtest + 2.5_rknd)
    Gamma_m = Gamma_aux(mtest + 2.5_rknd)

    cpq = (-1._rknd)**qtest * gamma_p / &
          ( gamma_q * ifactorial( ptest - qtest) * ifactorial(qtest))
    ckm = (-1._rknd)**mtest * gamma_k / &
          ( gamma_m * ifactorial( ktest - mtest) * ifactorial(mtest))

    ! Calculate Bstar
    Bstar_qm = calc_Bstar(qtest,mtest,theta,mu,chi,Q)
        
    ! Calculate portion of Bpk
    Btmp=cpq*ckm*Bstar_qm
        
    ! Total Bpk
    Bpk = Bpk + Btmp
  Enddo
Enddo

Nab = (2._rknd / chi) * Bpk

End Function calc_Nab_Ji

!-----------------------------------------------------------------------------
!+ Calculates the Mak Braginskii matrix element using Ji and Held formulae
!-----------------------------------------------------------------------------
Function calc_Mak_Ji(irow, icol, theta_ak, mu_ak, chi_ak, Q_ak) &
Result(Mak)
!
! Description: 
!  This function calculates the Mak Braginskii matrix element using the 
!  general formulae from Ji and Held.
!
! Function arguments:
!
!  QQ update this -- talk about method
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     05/25/2010  Original Code.  JL
! 
! Author(s): J. Lore 5/26/2010 

Implicit None

Integer, parameter :: rknd = selected_real_kind(12,300) 
Integer, parameter :: iknd = selected_int_kind(8)


! Input/output                      !See above for descriptions
Integer(iknd), Intent(in)  :: irow
Integer(iknd), Intent(in)  :: icol
Real(rknd),    Intent(in)  :: theta_ak
Real(rknd),    Intent(in)  :: mu_ak
Real(rknd),    Intent(in)  :: chi_ak
Real(rknd),    Intent(in)  :: Q_ak
Real(rknd)                 :: Mak

! Local scalars
Real(rknd)    :: Apk                 ! "A" coefficient 
Real(rknd)    :: cpq, ckm            ! "c" coefficients
Real(rknd)    ::                   & ! Gamma function evaluations 
  Gamma_p, Gamma_k, Gamma_q, Gamma_m
Integer(iknd) :: qtest, mtest        ! Loop indices


Real(rknd) :: Astar_qm, Atmp

!- End of header -------------------------------------------------------------

! Calculate the value of "A" for arbitrary p,k (irow, icol)
Apk = 0._rknd

Gamma_p = Gamma_aux(irow + 2.5_rknd)
Gamma_k = Gamma_aux(icol + 2.5_rknd)

Do qtest = 0_iknd, irow   ! (0:p)
  Do mtest = 0_iknd, icol ! (0:k)
    
    ! Calculate c coefficients
    Gamma_q = Gamma_aux(qtest + 2.5_rknd)
    Gamma_m = Gamma_aux(mtest + 2.5_rknd)

    cpq = (-1._rknd)**qtest * gamma_p / &
          ( gamma_q * ifactorial( irow - qtest) * ifactorial(qtest))
    ckm = (-1._rknd)**mtest * gamma_k / &
          ( gamma_m * ifactorial( icol - mtest) * ifactorial(mtest))

    ! Calculate A*_qm
    Astar_qm = calc_Astar(qtest,mtest,theta_ak,mu_ak,chi_ak,Q_ak)

    ! Calculate portion of  Apk
    Atmp = cpq * ckm * Astar_qm
        
    ! Total Apk
    Apk = Apk + Atmp
    

  Enddo ! m loop
Enddo ! q loop

Mak = 2._rknd * Apk

End Function calc_Mak_Ji


!-----------------------------------------------------------------------------
!+ Calculates the A*_qm element using ... QQ
!-----------------------------------------------------------------------------
Function calc_Astar(qtest,mtest,theta,mu,chi,Q) &
Result(Astar)
! Note this actually calculates Astar * (tau/3/n)

!
! Description: 
!
! Function arguments:
!
!  QQ update this -- talk about method
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     05/28/2010  Original Code.  JL
!  1.1     08/20/2010  Fixed bug due to case insensitivity JL
! 
! Author(s): J. Lore 5/26/2010-8/20/2010                           

Implicit None

Integer, parameter :: rknd = selected_real_kind(12,300) 
Integer, parameter :: iknd = selected_int_kind(8)

! Input/output                      !See above for descriptions
Integer(iknd), Intent(in)  :: qtest
Integer(iknd), Intent(in)  :: mtest
Real(rknd),    Intent(in)  :: theta
Real(rknd),    Intent(in)  :: mu
Real(rknd),    Intent(in)  :: chi   ! QQ check subscripts
Real(rknd),    Intent(in)  :: Q
Real(rknd)                 :: Astar

! Local parameters 
Integer(iknd), parameter :: iOne = 1_iknd
Integer(iknd), parameter :: iTwo = 2_iknd
Real(rknd), parameter :: One = 1._rknd
Real(rknd), parameter :: Two = 2._rknd

! Local Scalars
Real(rknd) :: a_big_E0, a_big_E1, a_big_E2
Real(rknd) :: a_lit_E0, a_lit_E1, a_lit_E2
Real(rknd) :: alpha_big_E0, alpha_big_E1, alpha_big_E2
Real(rknd) :: alpha_lit_e0, alpha_lit_e1, alpha_lit_e2
!- End of header -------------------------------------------------------------


! Calculate lower case a's
a_big_E0 = -Two * ( One - theta ) / mu
a_big_E1 = -One + ( One + Two * mtest ) * ( One - Two * theta ) / mu
a_big_E2 = Two * theta * mtest**iTwo / mu
a_lit_e0 = Two * ( One - theta ) * ( One / theta + One / mu )
a_lit_e1 = ( One + Two * mtest ) * ( One - One / mu + Two * theta / mu )
a_lit_e2 = -Two * theta * mtest**iTwo / mu

! Calc. alphas - actually alpha * (tau/3/n)

alpha_big_E0 = calc_Eba( iOne + qtest + mtest, chi, Q)
alpha_big_E1 = calc_Eba( qtest + mtest, chi, Q)

! alpha_E not defined for n=2, l=1, q=m=0 (aE2^00 = 0)
If ( qtest==0_iknd .and. mtest==0_iknd ) Then    
    alpha_big_E2 = 0._rknd
Else
    alpha_big_E2 = calc_Eba(qtest + mtest - iOne, chi, Q)
Endif

! chi**(1+2*k)*calc_eab(k,Q) gives eba

alpha_lit_e0 = chi**( iOne + iTwo * ( iTwo + qtest + mtest ) )  & 
                * calc_eab( iTwo + qtest + mtest, Q ) / chi
alpha_lit_e1 = chi**( iOne + iTwo * ( iOne + qtest + mtest ) )  &
                * calc_eab( iOne + qtest + mtest, Q ) / chi
alpha_lit_e2 = chi**( iOne + iTwo * ( qtest + mtest ) )         &
                * calc_eab( qtest + mtest , Q ) / chi

! actually Astar * (tau/3/n)
Astar = a_big_E0*alpha_big_E0 + a_big_E1*alpha_big_E1 + a_big_E2*alpha_big_E2 &
       + a_lit_e0*alpha_lit_e0 + a_lit_e1*alpha_lit_e1 + a_lit_e2*alpha_lit_e2

End Function calc_Astar


!-----------------------------------------------------------------------------
!+ Calculates the B*_qm element using ... QQ
!-----------------------------------------------------------------------------
Function calc_Bstar(qtest,mtest,theta,mu,chi,Q) &
Result(Bstar)
! Note this actually calculates Bstar * (tau/3/n)

!
! Description: 
!
! Function arguments:
!
!  QQ update this -- talk about method
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     05/28/2010  Original Code.  JL
!  1.1     08/20/2010  Chaged leading 5 in bp0 to 4 and 2 to 3 in beta_e1
! 
! Author(s): J. Lore 5/28/2010 - 8/20/2010
!
Implicit None

Integer, parameter :: rknd = selected_real_kind(12,300) 
Integer, parameter :: iknd = selected_int_kind(8)

! Input/output                      !See above for descriptions
Integer(iknd), Intent(in)  :: qtest
Integer(iknd), Intent(in)  :: mtest
Real(rknd),    Intent(in)  :: theta
Real(rknd),    Intent(in)  :: mu
Real(rknd),    Intent(in)  :: chi   ! QQ check subscripts
Real(rknd),    Intent(in)  :: Q
Real(rknd)                 :: bstar

! Local parameters
Integer(iknd), parameter :: iOne   = 1_iknd
Integer(iknd), parameter :: iTwo   = 2_iknd
Integer(iknd), parameter :: iThree = 3_iknd
Integer(iknd), parameter :: iFive  = 5_iknd
Real(rknd), parameter :: One    = 1._rknd
Real(rknd), parameter :: Two    = 2._rknd
Real(rknd), parameter :: Three  = 3._rknd
Real(rknd), parameter :: Four   = 4._rknd
Real(rknd), parameter :: Five   = 5._rknd
Real(rknd), parameter :: Eight  = 8._rknd

! Local Scalars
Integer(iknd) :: j
Real(rknd) :: b_e0, b_e1, b_p0, b_m0, b_m1
Real(rknd) :: beta_e0, beta_e1, beta_p0, beta_m0, beta_m1
Real(rknd) :: back_p0, back_m0, back_m1
!- End of header -------------------------------------------------------------


! Calculate lower case b's
b_e0 = Two * mu / theta**iTwo;
b_e1 = - Four / Five
b_p0 = Four + Eight / Five * mtest + Four/Three * mu /theta - Eight/Three/theta
b_m0 = -Eight/Three*mu/theta + Four/Three/theta
b_m1 = Eight/Five

!betas  - actually beta * (tau/3/n)
beta_e0 = chi**(iTwo * qtest + iFive) *calc_eab(iTwo+qtest+mtest,Q)
beta_e1 = chi**(iTwo * qtest + iFive) *calc_eab(iThree+qtest+mtest,Q)

back_p0=0._rknd
Do j = 0, qtest 
    back_p0 = back_p0 + ifactorial( qtest ) &
                       /ifactorial(j)*chi**(iTwo*j)*calc_eab(iTwo+mtest+j,Q)
Enddo

back_m0=0._rknd
back_m1=0._rknd
Do j=0_iknd, mtest
    back_m0 = back_m0 + One/ifactorial(j)*calc_eab(iTwo + qtest + j, Q)
    back_m1 = back_m1 + One/ifactorial(j)*calc_eab(iThree + qtest + j, Q)
Enddo

beta_p0=0.5_rknd*chi**iThree * back_p0
beta_m0=0.5_rknd*chi**(iTwo*qtest+iFive)*ifactorial(mtest)*back_m0
beta_m1=0.5_rknd*chi**(iTwo*qtest+iFive)*ifactorial(mtest)*back_m1

! actually Bstar * (tau/3/n)
Bstar=b_e0*beta_e0 + b_e1*beta_e1 + b_p0*beta_p0 + b_m0*beta_m0 + b_m1*beta_m1

End Function calc_Bstar

!-----------------------------------------------------------------------------
!+ Calculates ... QQ
!-----------------------------------------------------------------------------
Function calc_Eba(k,chi,Q) &
Result(Eba)

!
! Description: 
!
! Function arguments:
!
!  QQ update this -- talk about method
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     05/28/2010  Original Code.  JL
! 
! Author(s): J. Lore 5/28/2010 
!
!
Implicit None

Integer, parameter :: rknd = selected_real_kind(12,300) 
Integer, parameter :: iknd = selected_int_kind(8)

! Input/output                      !See above for descriptions
Integer(iknd), Intent(in)  :: k
Real(rknd),    Intent(in)  :: chi   ! QQ check subscripts
Real(rknd),    Intent(in)  :: Q
Real(rknd)                 :: Eba

! Local scalars
Real(rknd) :: coeff, Etmp
Integer(iknd) :: jtest

!- End of header -------------------------------------------------------------

coeff = ifactorial(k) / ( 2._rknd * chi )
Etmp = 0._rknd
Do jtest = 0_iknd, k 
    Etmp = Etmp + chi**( 1_iknd + 2_iknd * jtest ) &
                   * calc_eab (jtest, Q ) / ifactorial(jtest)
Enddo
Eba = coeff * Etmp


End Function calc_Eba


!-----------------------------------------------------------------------------
!+ Calculates ... QQ
!-----------------------------------------------------------------------------
Function calc_eab(k, Q) &
Result(eab)

!
! Description: 
!
! Function arguments:
!
!  QQ update this -- talk about method
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     05/28/2010  Original Code.  JL
! 
! Author(s): J. Lore 5/28/2010 
!
!

Implicit None

Integer, parameter :: rknd = selected_real_kind(12,300) 
Integer, parameter :: iknd = selected_int_kind(8)

! Input/output                      !See above for descriptions
Integer(iknd), Intent(in)  :: k ! QQ check subscripts 
Real(rknd),    Intent(in)  :: Q
Real(rknd)                 :: eab

! Local scalars
Real(rknd) ::  Xab, Gamma_k

! Local parameters
Real(rknd), parameter :: sqrtpi = 1.77245385090552_rknd

!- End of header -------------------------------------------------------------

Xab = Q**( -( k + 0.5_rknd ) )

Gamma_k = Gamma_aux(k + 0.5_rknd)

eab = Gamma_k * Xab / sqrtpi

End Function calc_eab


!-----------------------------------------------------------------------------
!+ Calculates the factorial of a positive integer
!-----------------------------------------------------------------------------
Function ifactorial(Nval) & 
Result(ifact)
!
! Description: 
!   This function simply calculates the factorial of an integer.
!
! Inputs:
!  Nval: the input [integer]
! Outputs:
!  ifact: The factorial [integer]
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     07/22/2009  Original Code.  JL
!  1.1     09/01/2010  Updated for PENTA3. JL
!
! Author(s): J. Lore 07/2009 - 9/1/2010 

Implicit None

Integer, parameter :: rknd = selected_real_kind(12,300) 
Integer, parameter :: iknd = selected_int_kind(8)

! Input/output                    ! See above for descriptions
Integer(iknd), intent(in)   :: Nval
Integer(iknd)               :: ifact

! Local scalars
Integer(iknd)             :: icount  ! loop index
!- End of header -------------------------------------------------------------

If ( Nval < 0 ) Then
  Stop 'Error: function ifactorial called for a negative number'
Elseif ( Nval > 20 ) Then
  Stop 'Error: function ifactorial called for arg > 20'
Endif

! QQ check for max value, zero

ifact = 1_iknd
Do icount = 1, Nval
  ifact = ifact*icount
Enddo

End Function ifactorial

!-----------------------------------------------------------------------------
!+ Calculates the gamma function
!-----------------------------------------------------------------------------
Function Gamma_aux(X)  &
Result(GA)
!
! Description: 
!       ==================================================
!       Purpose: Compute the gamma function �(x)
!       Input :  x  --- Argument of �(x)
!                       ( x is not equal to 0,-1,-2,��� )
!       Output:  GA --- �(x)
!       ==================================================
!
!   Source (copyrighted):
!     http://jin.ece.uiuc.edu/routines/routines.html
!
!     Updated to Fortran 90 with kind precision Jl 5/26/2010
!       Also added stop for negative integers
!

Implicit None

Integer, parameter :: rknd = selected_real_kind(12,300) 
Integer, parameter :: iknd = selected_int_kind(8)

Real(rknd), Intent(in) :: X
Real(rknd)             :: GA
Real(rknd) :: G(26),R,Z
Integer(iknd) :: K, M
Real(rknd) :: Pi, Gr
!- End of header -------------------------------------------------------------

Pi = 3.141592653589793_rknd

If ( X == Int(X,iknd) ) Then  ! If X is integer

  If ( X > 0.0_rknd ) Then  ! If X integer and positive (perform (X-1)!)
    GA = ifactorial( Int(X - 1._rknd,iknd) )
  Else                  ! If X is integer and negative
    Stop ' Subroutine Gamma called for negative integer (Inf)'
  Endif

Else  ! if X not integer

  If ( Dabs(X) > 1.0_rknd ) Then  ! If |X| > 1
    Z = Dabs(X)
    M=Int(Z,iknd)
    R=1.0_rknd
    Do K = 1, M
      R = R * (Z - K)
    Enddo
    Z = Z - M
  Else
    Z=X
  Endif

  Data  G/1.0D0              , 0.5772156649015329D0,  &
        -0.6558780715202538D0, -0.420026350340952D-1, &
        0.1665386113822915D0 ,  -.421977345555443D-1, &
        -.96219715278770D-2  , .72189432466630D-2,    &
        -.11651675918591D-2  , -.2152416741149D-3,    &
        .1280502823882D-3    , -.201348547807D-4,     &
        -.12504934821D-5     , .11330272320D-5,       &
        -.2056338417D-6      , .61160950D-8,          &
        .50020075D-8         , -.11812746D-8,         &
        .1043427D-9          , .77823D-11,            &
        -.36968D-11          , .51D-12,               &
        -.206D-13            , -.54D-14, .14D-14, .1D-15/

  GR=G(26)
  
  Do K = 25, 1, -1
    GR = GR * Z + G(K)
  Enddo

  GA = 1.0_rknd / (GR * Z)
  If ( Dabs(X) > 1.0_rknd ) Then
    GA = GA * R
    If (X < 0.0_rknd) Then
      GA = - PI / ( X*GA*Dsin(PI*X) )
    Endif
  Endif
Endif

Return
EndFunction Gamma_aux

end module transport_mod