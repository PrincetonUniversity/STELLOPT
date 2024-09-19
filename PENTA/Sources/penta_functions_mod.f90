!-----------------------------------------------------------------------------
!+ Contains various functions for PENTA
!-----------------------------------------------------------------------------
Module penta_functions_mod
! Description: 
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     05/25/2010  Original Code.  JL
!  1.1     09/23/2011  Overhaul to match new documentation.
! 
! Author(s): J. Lore 05/25/2010 - 09/23/2011
!
Contains

!-----------------------------------------------------------------------------
!+ Calculates the radial energy fluxes using direct energy convolution
!-----------------------------------------------------------------------------
Function calc_QoTs_DKES(num_species,abs_Er,Temps,dens,vths,charges,      &
     masses,loglambda,use_quanc8,Kmin,Kmax,numKsteps,log_interp,         &
     cmin,cmax,emin,emax,xt_c,xt_e,Dspl_logD11,Dspl_D31,num_c,num_e,     &
     kcord,keord,Avec,B0)                                                &
Result(QoTs)
!
! Description: 
! This function calculates the radial energy fluxes using direct energy 
!  convolution.
!
! Function arguments:
!  num_species: Total number of particle species
!  abs_Er: |Er| in V/cm
!  Temps: Array of species' temperatures in eV [Te, Ti1, Ti2...]
!  dens: Same as above but densitiy in [m**-3]
!  vths: Same as above but thermal velocities in m/s
!  charges: Same as above but charge in [C]
!  masses: Same as above but mass in kg
!  loglambda: Coulomb logarithm (assumed equal for each species)
!  use_quanc8: Logical switch for use of adapative integrator
!  Kmin, Kmax: Range of normalized energy over which to integrate
!  numKsteps:  Number of integration points for rectangular approx.
!  log_interp: Logical switch to activate 2D log. interp. of DKES coeffs.
!  cmin, cmax: Max, min value of cmul for each DKES matrix
!  emin, emax: Ditto for efield arrays
!  xt_c,Dspl_D11 etc: Spline knots and coefficients for DKES coeffs.
!  num_c,num_e: Number of cmul, efield values in array
!  kcord, keord: Spline orders for interpolation (cmul, efield)
!  Avec: Array of "A-type" thermodynamic forces
!  B0: B(m=0,n=0) in Boozer coordinates
!
! Output:
!  Total energy flux of each species divided by eT [m**-2]
!   [Q_e/eT_e, Q_i1/eT_i1, ...]
!
! See J. Lore's PhD thesis for expressions.  See also DKES references:
!  Hirshman, et al, Phys. Fluids 29, 2951 (1986)
!  van Rij and Hirshman, Phys. Fluids B 1, 563 (1989)
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     09/29/2010  Original Code.  J
!  1.1     09/23/2011  Update to match new doc. JL
!  
! Author(s): J. Lore 9/29/2010 - 09/23/2011
!
! Modules used:
Use penta_kind_mod                  ! Import rknd, iknd specifications

Implicit None

! Input/output                      !See above for descriptions
Integer(iknd), Intent(in)  :: num_species
Real(rknd),    Intent(in)  :: abs_Er
Real(rknd),    Intent(in)  :: Temps(num_species)
Real(rknd),    Intent(in)  :: dens(num_species)
Real(rknd),    Intent(in)  :: vths(num_species)
Real(rknd),    Intent(in)  :: charges(num_species)
Real(rknd),    Intent(in)  :: masses(num_species)
Real(rknd),    Intent(in)  :: loglambda
Logical,       Intent(in)  :: use_quanc8
Real(rknd),    Intent(in)  :: Kmin
Real(rknd),    Intent(in)  :: Kmax
Integer(iknd), Intent(in)  :: numKsteps
Logical,       Intent(in)  :: log_interp
Real(rknd),    Intent(in)  :: cmin
Real(rknd),    Intent(in)  :: cmax
Real(rknd),    Intent(in)  :: emin
Real(rknd),    Intent(in)  :: emax
Real(rknd),    Intent(in)  :: xt_c(num_c + kcord)
Real(rknd),    Intent(in)  :: xt_e(num_e + keord)
Real(rknd),    Intent(in)  :: Dspl_logD11(num_c,num_e)
Real(rknd),    Intent(in)  :: Dspl_D31(num_c,num_e)
Integer(iknd), Intent(in)  :: num_c
Integer(iknd), Intent(in)  :: num_e
Integer(iknd), Intent(in)  :: kcord
Integer(iknd), Intent(in)  :: keord
Real(rknd),    Intent(in)  :: Avec(num_species*3)
Real(rknd),    Intent(in)  :: B0
Real(rknd)                 :: QoTs(num_species)

! Local Scalars
Integer(iknd) :: ispec1,ind_A     ! Loop indices
Integer(iknd) ::  I               ! Index used for array constructors
Integer(iknd) ::  nu_exp          ! Exponent on collision freq. for conv.
Real(rknd)    ::                & ! Primary species (a) paramaters 
  ma,Ta,vta,qa,na
Real(rknd)    :: L21,L22,L23  ! Thermal diffusion coefficients
Real(rknd)    ::  norm_factor     ! Constants (wrt K) for convolution terms
Real(rknd)    ::  K_exp           ! Exponent on K for convolution (ignoring 
                                  !  Maxwellian contribution)

! Local arrays 
Integer(iknd) :: ikeep(num_species-1)   ! Used to index species 'b' in arrays
Real(rknd)    :: qb(num_species-1),  &  ! Paremeters for species 'b' /= 'a'
  nb(num_species-1), vtb(num_species-1) 


! Local parameters                 
Integer(iknd), Parameter ::  &
iZERO  = 0_iknd,      &
iTWO   = 2_iknd,      &
iTHREE = 3_iknd

Real(rknd), Parameter :: &
ZERO      = 0._rknd,  &
ONE       = 1._rknd,  &
THREEHALF = 1.5_rknd, &  
FIVEHALF  = 2.5_rknd, &  
SEVENHALF = 3.5_rknd, &  
TWO       = 2._rknd

!- End of header -------------------------------------------------------------
Do ispec1 = 1_iknd, num_species

  ! Assign species 'a' parameters
  ma  = masses(ispec1)
  qa  = charges(ispec1)
  Ta  = Temps(ispec1)
  na  = dens(ispec1)
  vta = vths(ispec1)

  ! Get indices of all species b =/ a
  ikeep=(/(I,I=1,ispec1-1), (I,I=ispec1+1,num_species)/)

  ! Assign species 'b' parameters
  qb  = charges(ikeep)
  nb  = dens(ikeep)
  vtb = vths(ikeep)  

  ! Integrate the radial diffusion coefficient to get the 
  !  L11, L12, etc coefficients

  norm_factor = na*ma*ma*vta*vta*vta/(TWO*qa*qa)
  nu_exp = iZERO  ! Integer

  K_exp = FIVEHALF
  L21 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps,       &
    vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
    Dspl_logD11,xt_c,xt_e,cmin,cmax,emin,emax, &
    num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false., &
    .true.)

  K_exp = SEVENHALF
  L22 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps,       &
    vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
    Dspl_logD11,xt_c,xt_e,cmin,cmax,emin,emax, &
    num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false., &
    .true.)

  ! If E_prl /= 0, then calculate L23
  ! Strictly shouldn't compare floating points, but the price is only
  !  calculating the coefficient anyway... 
  ind_A = (ispec1 - 1)*3 + 1
  If ( Avec(ind_A + 2) == ZERO ) Then 
    L23 = ZERO
  Else
    norm_factor = na*ma*vta*vta/(TWO*B0*qa)
    K_exp  = TWO    ! Real
    nu_exp = iZERO  ! Integer
    L23 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps,       &
      vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
      Dspl_D31,xt_c,xt_e,cmin,cmax,emin,emax, &
      num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false., &
      .false.)
  Endif

  ! Calculate flux due to direct energy convolution of DKES coefficients
  ind_A = (ispec1-1)*3+1
  QoTs(ispec1) = - L21*Avec(ind_A) - L22*Avec(ind_A+1) + L23*Avec(ind_A+2)

EndDo ! Species 1 loop

EndFunction calc_QoTs_DKES

!-----------------------------------------------------------------------------
!+ Calculates the radial particle fluxes from direct energy convolution.
!-----------------------------------------------------------------------------
Function calc_fluxes_DKES(num_species,abs_Er,Temps,dens,vths,charges,  &
     masses,loglambda,use_quanc8,Kmin,Kmax,numKsteps,log_interp,       &
     cmin,cmax,emin,emax,xt_c,xt_e,Dspl_logD11,Dspl_D31,num_c,num_e,   &
     kcord,keord,Avec,B0)                                              &
Result(Gammas)
!
! Description: 
! This function calculates the radial energy fluxes using direct energy 
!  convolution.
!
! Function arguments:
!  num_species: Total number of particle species
!  abs_Er: |Er| in V/cm
!  Temps: Array of species' temperatures in eV [Te, Ti1, Ti2...]
!  dens: Same as above but densitiy in [m**-3]
!  vths: Same as above but thermal velocities in m/s
!  charges: Same as above but charge in [C]
!  masses: Same as above but mass in kg
!  loglambda: Coulomb logarithm (assumed equal for each species)
!  use_quanc8: Logical switch for use of adapative integrator
!  Kmin, Kmax: Range of normalized energy over which to integrate
!  numKsteps:  Number of integration points for rectangular approx.
!  log_interp: Logical switch to activate 2D log. interp. of DKES coeffs.
!  cmin, cmax: Max, min value of cmul for each DKES matrix
!  emin, emax: Ditto for efield arrays
!  xt_c,Dspl_D11 etc: Spline knots and coefficients for DKES coeffs.
!  num_c,num_e: Number of cmul, efield values in array
!  kcord, keord: Spline orders for interpolation (cmul, efield)
!  Avec: Array of "A-type" thermodynamic forces
!  B0: B(m=0,n=0) in Boozer coordinates

! Output:
!  Radial particle flux of each species [m**-2s**-1]
!   [Gamma_e, Gamma_i1, ...]
!  
! See J. Lore's PhD thesis for expressions.  See also DKES references:
!  Hirshman, et al, Phys. Fluids 29, 2951 (1986)
!  van Rij and Hirshman, Phys. Fluids B 1, 563 (1989)
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     09/10/2010  Original Code.  JL
!  1.1     09/23/2011  Update to match new doc. JL
! 
! Author(s): J. Lore 9/10/2010 - 9/23/2011
!
!
! Modules used:
Use penta_kind_mod                  ! Import rknd, iknd specifications

Implicit None

! Input/output                      !See above for descriptions
Integer(iknd), Intent(in)  :: num_species
Real(rknd),    Intent(in)  :: abs_Er
Real(rknd),    Intent(in)  :: Temps(num_species)
Real(rknd),    Intent(in)  :: dens(num_species)
Real(rknd),    Intent(in)  :: vths(num_species)
Real(rknd),    Intent(in)  :: charges(num_species)
Real(rknd),    Intent(in)  :: masses(num_species)
Real(rknd),    Intent(in)  :: loglambda
Logical,       Intent(in)  :: use_quanc8
Real(rknd),    Intent(in)  :: Kmin
Real(rknd),    Intent(in)  :: Kmax
Integer(iknd), Intent(in)  :: numKsteps
Logical,       Intent(in)  :: log_interp
Real(rknd),    Intent(in)  :: cmin
Real(rknd),    Intent(in)  :: cmax
Real(rknd),    Intent(in)  :: emin
Real(rknd),    Intent(in)  :: emax
Real(rknd),    Intent(in)  :: xt_c(num_c + kcord)
Real(rknd),    Intent(in)  :: xt_e(num_e + keord)
Real(rknd),    Intent(in)  :: Dspl_logD11(num_c,num_e)
Real(rknd),    Intent(in)  :: Dspl_D31(num_c,num_e)
Integer(iknd), Intent(in)  :: num_c
Integer(iknd), Intent(in)  :: num_e
Integer(iknd), Intent(in)  :: kcord
Integer(iknd), Intent(in)  :: keord
Real(rknd),    Intent(in)  :: Avec(num_species*3)
Real(rknd),    Intent(in)  :: B0
Real(rknd)                 :: Gammas(num_species)

! Local Scalars
Integer(iknd) :: ispec1,ind_A     ! Loop indices
Integer(iknd) ::  I               ! Index used for array constructors
Integer(iknd) ::  nu_exp          ! Exponent on collision freq. for conv.
Real(rknd)    ::                & ! Primary species (a) paramaters 
  ma,Ta,vta,qa,na
Real(rknd)    :: L11,L12,L13      ! Thermal diffusion coefficients
Real(rknd)    ::  norm_factor     ! Constants (wrt K) for convolution terms
Real(rknd)    ::  K_exp           ! Exponent on K for convolution (ignoring 
                                  !  Maxwellian contribution)

! Local arrays 
Integer(iknd) :: ikeep(num_species-1)   ! Used to index species 'b' in arrays
Real(rknd)    :: qb(num_species-1),  &  ! Paremeters for species 'b' /= 'a'
  nb(num_species-1), vtb(num_species-1) 


! Local parameters                 
Integer(iknd), Parameter ::  &
iZERO  = 0_iknd,      &
iONE   = 1_iknd,      &
iTWO   = 2_iknd,      &
iTHREE = 3_iknd

Real(rknd), Parameter :: &
ZERO      = 0._rknd,  &
ONE       = 1._rknd,  &
THREEHALF = 1.5_rknd, &  
TWO       = 2._rknd,  &
FIVEHALF  = 2.5_rknd, &
THREE     = 3._rknd,  &  
FIVE      = 5._rknd,  &
SEVENHALF = 3.5_rknd  

!- End of header -------------------------------------------------------------
Do ispec1 = 1_iknd, num_species

  ! Assign species 'a' parameters
  ma  = masses(ispec1)
  qa  = charges(ispec1)
  Ta  = Temps(ispec1)
  na  = dens(ispec1)
  vta = vths(ispec1)

  ! Get indices of all species b =/ a
  ikeep=(/(I,I=1,ispec1-1), (I,I=ispec1+1,num_species)/)

  ! Assign species 'b' parameters
  qb  = charges(ikeep)
  nb  = dens(ikeep)
  vtb = vths(ikeep)  

  ! Integrate the radial diffusion coefficient to get the 
  !  L11, L12, etc coefficients

  norm_factor = na*ma*ma*vta*vta*vta/(TWO*qa*qa)
  K_exp = THREEHALF     ! Real
  nu_exp = iZERO  ! Integer
  L11 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps,       &
    vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
    Dspl_logD11,xt_c,xt_e,cmin,cmax,emin,emax, &
    num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.true.)

  K_exp = FIVEHALF
  L12 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps,       &
    vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
    Dspl_logD11,xt_c,xt_e,cmin,cmax,emin,emax, &
    num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.true.)


  ! If E_prl /= 0, then calculate L13
  ! Strictly shouldn't compare floating points, but the price is simply
  !  calculating the coefficient anyway... QQ
  ind_A = (ispec1 - 1)*3 + 1
  If ( Avec(ind_A + 2) == 0._rknd ) Then 
    L13 = ZERO
  Else
    norm_factor = na*vta*vta*ma/(TWO*B0*qa)  
    K_exp  = ONE    ! Real
    nu_exp = iZERO  ! Integer
    L13 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps,      &
      vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
      Dspl_D31,xt_c,xt_e,cmin,cmax,emin,emax, &
      num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.false.)
  Endif

  ! Calculate flux due to direct energy convolution of DKES coefficients
  ind_A = (ispec1-1)*3+1
  Gammas(ispec1) = - L11*Avec(ind_A) - L12*Avec(ind_A+1) + L13*Avec(ind_A+2)

EndDo ! Species 1 loop

EndFunction calc_fluxes_DKES

!------------------------------------------------------------------------------
!+ Calculates the parallel flow moments only considering primary species drives
!------------------------------------------------------------------------------
Function calc_flows_DKES(num_species,Smax,abs_Er,Temps,dens,vths,charges,  &
     masses,loglambda,B0,use_quanc8,Kmin,Kmax,numKsteps,log_interp,        &
     cmin,cmax,emin,emax,xt_c,xt_e,Dspl_D31,Dspl_logD33,num_c,num_e,       &
     kcord,keord,Avec)                                                     &
Result(Flows)
!
! Description: 
! This function calculates the parallel flow moments using only the uncorrected
!  DKES coeffs expanded to an arbitrary number of Sonine polynomials.
!
! Function arguments:
!  num_species: Total number of particle species
!  Smax: Upper index on Sonine summation
!  abs_Er: |Er| in V/cm
!  Temps: Array of species' temperatures in eV [Te, Ti1, Ti2...]
!  dens: Same as above but densitiy in [m**-3]
!  vths: Same as above but thermal velocities in m/s
!  charges: Same as above but charge in [C]
!  masses: Same as above but mass in kg
!  loglambda: Coulomb logarithm (assumed equal for each species)
!  B0: B(m=0,n=0) in Boozer coordinates
!  use_quanc8: Logical switch for use of adapative integrator
!  Kmin, Kmax: Range of normalized energy over which to integrate
!  numKsteps:  Number of integration points for rectangular approx.
!  log_interp: Logical switch to activate 2D log. interp. of DKES coeffs.
!  cmin, cmax: Max, min value of cmul for each DKES matrix
!  emin, emax: Ditto for efield arrays
!  xt_c,Dspl_D11 etc: Spline knots and coefficients for DKES coeffs.
!  num_c,num_e: Number of cmul, efield values in array
!  kcord, keord: Spline orders for interpolation (cmul, efield)
!  Avec: Array of "A-type" thermodynamic forces
!  B0: B(m=0,n=0) in Boozer coordinates
!
! Output:
!  Parallel flow moments <Bu_ks>/<B**2> where k is the Sonine index
!   and s is the species index.  [<Bu_0e> <Bu_1e> ...<Bu_0i> ...]/<B**2>
!
! See J. Lore's PhD thesis for expressions.  See also DKES references:
!  Hirshman, et al, Phys. Fluids 29, 2951 (1986)
!  van Rij and Hirshman, Phys. Fluids B 1, 563 (1989)
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     09/10/2010  Original Code.  JL
!  1.1     09/23/2011  Update to match new doc. JL
! 
! Author(s): J. Lore 9/10/2010 - 09/23/2011
!
! Modules used:
Use penta_kind_mod                  ! Import rknd, iknd specifications

Implicit None

! Input/output                      !See above for descriptions
Integer(iknd), Intent(in)  :: num_species
Integer(iknd), Intent(in)  :: Smax
Real(rknd),    Intent(in)  :: abs_Er
Real(rknd),    Intent(in)  :: Temps(num_species)
Real(rknd),    Intent(in)  :: dens(num_species)
Real(rknd),    Intent(in)  :: vths(num_species)
Real(rknd),    Intent(in)  :: charges(num_species)
Real(rknd),    Intent(in)  :: masses(num_species)
Real(rknd),    Intent(in)  :: loglambda
Real(rknd),    Intent(in)  :: B0
Logical,       Intent(in)  :: use_quanc8
Real(rknd),    Intent(in)  :: Kmin
Real(rknd),    Intent(in)  :: Kmax
Integer(iknd), Intent(in)  :: numKsteps
Logical,       Intent(in)  :: log_interp
Real(rknd),    Intent(in)  :: cmin
Real(rknd),    Intent(in)  :: cmax
Real(rknd),    Intent(in)  :: emin
Real(rknd),    Intent(in)  :: emax
Real(rknd),    Intent(in)  :: xt_c(num_c + kcord)
Real(rknd),    Intent(in)  :: xt_e(num_e + keord)
Real(rknd),    Intent(in)  :: Dspl_D31(num_c,num_e)
Real(rknd),    Intent(in)  :: Dspl_logD33(num_c,num_e)
Integer(iknd), Intent(in)  :: num_c
Integer(iknd), Intent(in)  :: num_e
Integer(iknd), Intent(in)  :: kcord
Integer(iknd), Intent(in)  :: keord
Real(rknd),    Intent(in)  :: Avec(num_species*3)
Real(rknd)                 :: Flows(num_species*(Smax+1))

! Local scalars
Integer(iknd) ::  ispec1, jval, ind_A, ind_RHS ! Loop indices
Integer(iknd) ::  I               ! Index used for array constructors
Integer(iknd) ::  nu_exp          ! Exponent on collision freq. for conv.

Real(rknd)    ::                & ! Primary species (a) paramaters 
  ma, Ta, vta, qa, na     
Real(rknd)    ::  norm_factor     ! Constants (wrt K) for convolution terms
Real(rknd)    ::  K_exp           ! Exponent on K for convolution (ignoring 
                                  !  Maxwellian contribution)
Real(rknd)    ::  RHS_1,RHS_2,  & ! RHS elements of flow equation
  RHS_3

! Local arrays 
Integer(iknd) :: ikeep(num_species-1)   ! Used to index species 'b' in arrays
Real(rknd)    :: qb(num_species-1),  &  ! Paremeters for species 'b' /= 'a'
  nb(num_species-1), vtb(num_species-1) 
Real(rknd)   :: RHS(num_species*(Smax+1)) ! 1D array for RHS of eq. sys

! Local parameters                 
Integer(iknd), Parameter ::  &
iZERO = 0_iknd,      &
iONE  = 1_iknd,      &
iTWO  = 2_iknd

Real(rknd), Parameter :: &
HALF = 0.5_rknd,    &
ZERO = 0._rknd,     &
ONE  = 1._rknd,     &
TWO  = 2._rknd 

!- End of header -------------------------------------------------------------

! Loop over species and calculate the parallel flows
Do ispec1 = 1_iknd,num_species

  ! Assign species 'a' parameters
  ma  = masses(ispec1)
  qa  = charges(ispec1)
  Ta  = Temps(ispec1)
  na  = dens(ispec1)
  vta = vths(ispec1)

  ! Get indices of all species b =/ a
  ikeep=(/(I,I=1,ispec1-1), (I,I=ispec1+1,num_species)/)

  ! Assign species 'b' parameters
  qb  = charges(ikeep)
  nb  = dens(ikeep)
  vtb = vths(ikeep)

  ! The 'j' loop defines each equation for the current species in the system
  Do jval = 0_iknd,Smax
  
    ! Define the RHS of the equation system

    ! First RHS term
    norm_factor = vta*vta*ma/(TWO*qa*B0)
    K_exp = ONE     ! Real
    nu_exp = iZERO  ! Integer
    RHS_1 = energy_conv(jval,iZERO,use_quanc8,Kmin,Kmax,numKsteps,      &
      vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
      Dspl_D31,xt_c,xt_e,cmin,cmax,emin,emax, &
      num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.false.)

    ! 2nd RHS term
    K_exp  = TWO    ! Real
    RHS_2 = energy_conv(jval,iZERO,use_quanc8,Kmin,Kmax,numKsteps,      &
      vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
      Dspl_D31,xt_c,xt_e,cmin,cmax,emin,emax, &
      num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.false.)

    ! If E_prl /= 0, then calculate 3rd RHS Term
    ! Strictly shouldn't compare floating points, but the price is simply
    !  calculating the coefficient anyway... QQ
    ind_A = (ispec1 - 1)*3 + 1
    If ( Avec(ind_A + 2) == 0._rknd ) Then 
      RHS_3 = ZERO
    Else
      norm_factor = vta/(TWO*B0*B0)
      K_exp  = HALF    ! Real
      nu_exp = iZERO  ! Integer
      RHS_3 = energy_conv(jval,iZERO,use_quanc8,Kmin,Kmax,numKsteps,      &
        vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
        Dspl_logD33,xt_c,xt_e,cmin,cmax,emin,emax, &
        num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.true.)
    Endif

    ! Sum the RHS terms multiplied by the forces 
    ind_RHS = (ispec1 - 1)*(Smax+1)+jval+1
    RHS(ind_RHS) = - RHS_1*Avec(ind_A) - RHS_2*Avec(ind_A+1) &
      - RHS_3*Avec(ind_A+2)

  EndDo ! jval loop
EndDo ! Primary species loop

Flows = RHS

EndFunction calc_flows_DKES

!-----------------------------------------------------------------------------
!+ Calculates the radial energy fluxes using the SN method
!-----------------------------------------------------------------------------
Function calc_QoTs_SN(num_species,Smax,abs_Er,Temps,dens,vths,charges,      &
     masses,loglambda,use_quanc8,Kmin,Kmax,numKsteps,log_interp,cmin,       &
     cmax,emin,emax,xt_c,xt_e,Dspl_Drat,Dspl_Drat2,Dspl_Dex,Dspl_logD11,    &
     Dspl_D31,num_c,num_e,kcord,keord,Avec,Bsq,lmat,Flows,U2,dTdrs,      &
     dndrs,flux_cap)                                                        &
Result(QoTs)
!
! Description: 
! This function calculates the radial energy fluxes using the method of 
! SN expanded to an arbitrary number of Sonine polynomials.
!
! Function arguments:
!  num_species: Total number of particle species
!  Smax: Upper index on Sonine summation
!  abs_Er: |Er| in V/cm
!  Temps: Array of species' temperatures in eV [Te, Ti1, Ti2...]
!  dens: Same as above but densitiy in [m**-3]
!  vths: Same as above but thermal velocities in m/s
!  charges: Same as above but charge in [C]
!  masses: Same as above but mass in kg
!  loglambda: Coulomb logarithm (assumed equal for each species)
!  use_quanc8: Logical switch for use of adapative integrator
!  Kmin, Kmax: Range of normalized energy over which to integrate
!  numKsteps:  Number of integration points for rectangular approx.
!  log_interp: Logical switch to activate 2D log. interp. of DKES coeffs.
!  cmin, cmax: Max, min value of cmul for each DKES matrix
!  emin, emax: Ditto for efield arrays
!  xt_c,Dspl_D11 etc: Spline knots and coefficients for DKES coeffs.
!  num_c,num_e: Number of cmul, efield values in array
!  kcord, keord: Spline orders for interpolation (cmul, efield)
!  Avec: Array of "A-type" thermodynamic forces
!  Bsq: The quantity <B**2> on this surface
!  B0: B(m=0,n=0) in Boozer coordinates
!  lmat: Matrix of classical friction coefficients
!  Flows: Array of parallel flow moments
!  U2: Related to the PS factor <U**2>
!
! Output:
!  Total energy flux of each species divided by eT [m**-2]
!   [Q_e/eT_e, Q_i1/eT_i1, ...]
!
! See J. Lore's PhD thesis for expressions.  See also:
!         Sugama, Nishimura PoP 9 4637 (2002), 
!         Sugama, Nishimura PoP 15, 042502 (2008),
!      
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     09/10/2010  Original Code.  JL
!  1.1     01/06/2011  Added convective part of PS QoT. JL
!  1.2     09/23/2011  Update to match new documentation. JL
!
! Author(s): J. Lore 09/10/2010  - 09/23/2011
!
!
! Modules used:
Use penta_kind_mod                  ! Import rknd, iknd specifications
Use phys_const, Only :  &           ! Import physical constants
! Imported parameters
elem_charge

Implicit None

! Input/output                      !See above for descriptions
Integer(iknd), Intent(in)  :: num_species
Integer(iknd), Intent(in)  :: Smax
Real(rknd),    Intent(in)  :: abs_Er
Real(rknd),    Intent(in)  :: Temps(num_species)
Real(rknd),    Intent(in)  :: dens(num_species)
Real(rknd),    Intent(in)  :: vths(num_species)
Real(rknd),    Intent(in)  :: charges(num_species)
Real(rknd),    Intent(in)  :: masses(num_species)
Real(rknd),    Intent(in)  :: loglambda
Logical,       Intent(in)  :: use_quanc8
Real(rknd),    Intent(in)  :: Kmin
Real(rknd),    Intent(in)  :: Kmax
Integer(iknd), Intent(in)  :: numKsteps
Logical,       Intent(in)  :: log_interp
Real(rknd),    Intent(in)  :: cmin
Real(rknd),    Intent(in)  :: cmax
Real(rknd),    Intent(in)  :: emin
Real(rknd),    Intent(in)  :: emax
Real(rknd),    Intent(in)  :: xt_c(num_c + kcord)
Real(rknd),    Intent(in)  :: xt_e(num_e + keord)
Real(rknd),    Intent(in)  :: Dspl_Drat(num_c,num_e)
Real(rknd),    Intent(in)  :: Dspl_Drat2(num_c,num_e)
Real(rknd),    Intent(in)  :: Dspl_Dex(num_c,num_e)
Real(rknd),    Intent(in)  :: Dspl_logD11(num_c,num_e)
Real(rknd),    Intent(in)  :: Dspl_D31(num_c,num_e)
Integer(iknd), Intent(in)  :: num_c
Integer(iknd), Intent(in)  :: num_e
Integer(iknd), Intent(in)  :: kcord
Integer(iknd), Intent(in)  :: keord
Real(rknd),    Intent(in)  :: Avec(num_species*3)
Real(rknd),    Intent(in)  :: Bsq
Real(rknd),    Intent(in)  :: lmat(num_species*(Smax+1),num_species*(Smax+1))
Real(rknd),    Intent(in)  :: Flows(num_species*(Smax+1))
Real(rknd),    Intent(in)  :: U2
Real(rknd),    Intent(in)  :: dTdrs(num_species)
Real(rknd),    Intent(in)  :: dndrs(num_species)
Logical,       Intent(in)  :: flux_cap
Real(rknd)                 :: QoTs(num_species)

! Local Scalars
Integer(iknd) :: ispec1,kval,   & ! Loop indices
  ispec2,lmat_ind1,        &
  lmat_ind2,          &
  flow_ind1,ind_A 
Integer(iknd) ::  I               ! Index used for array constructors
Integer(iknd) ::  nu_exp          ! Exponent on collision freq. for conv.
Real(rknd)    :: dn_betadr,     & ! Species beta parameters (for PS flux)
  dT_betadr,n_beta,T_beta,      &
  q_beta     
Real(rknd)    ::                & ! Primary species (a) paramaters 
  ma, Ta, vta, qa, na     
Real(rknd)    :: L21,L22          ! Thermal diffusion coefficients
Real(rknd)    :: L21_1, L21_2          ! Thermal diffusion coefficients
Real(rknd)    :: L22_1, L22_2          ! Thermal diffusion coefficients
Real(rknd)    :: norm_factor      ! Constants (wrt K) for convolution terms
Real(rknd)    :: K_exp            ! Exponent on K for convolution (ignoring 
                                  !  Maxwellian contribution)
Real(rknd)    :: lab21,lab22,   & ! Friction coefficient values for PS flux
  lab11,lab12      
Real(rknd)    :: I_1,I_2          ! Integral factors used in calculating PS flux
Real(rknd)    :: QoT_PS_std,    & ! PS fluxes
  QoT_PS_flow
Real(rknd)    :: C_PS1, C_PS2

! Local arrays 
Integer(iknd) :: ikeep(num_species-1)   ! Used to index species 'b' in arrays
Real(rknd)    :: Na_2k(Smax+1)          ! QQ
Real(rknd)    :: qb(num_species-1),  &  ! Paremeters for species 'b' /= 'a'
  nb(num_species-1), vtb(num_species-1) 
Real(rknd)    :: QoT_PS(num_species)
Real(rknd)    :: mono_QoT(num_species)
Real(rknd)    :: QoT_Ua(num_species)

! Local parameters                 
Integer(iknd), Parameter ::  &
iZERO  = 0_iknd,      &
iONE   = 1_iknd,      &
iTWO   = 2_iknd,      &
iTHREE = 3_iknd

Real(rknd), Parameter :: &
ONE       = 1._rknd,  &
THREEHALF = 1.5_rknd, &  
TWO       = 2._rknd,  &
FIVEHALF  = 2.5_rknd, &  
SEVENHALF  = 3.5_rknd, &  
THREE     = 3._rknd,  &
FIVE      = 5._rknd

!- End of header -------------------------------------------------------------
Do ispec1 = 1_iknd, num_species

  ! Assign species 'a' parameters
  ma  = masses(ispec1)
  qa  = charges(ispec1)
  Ta  = Temps(ispec1)
  na  = dens(ispec1)
  vta = vths(ispec1)

  ! Get indices of all species b =/ a
  ikeep=(/(I,I=1,ispec1-1), (I,I=ispec1+1,num_species)/)

  ! Assign species 'b' parameters
  qb  = charges(ikeep)
  nb  = dens(ikeep)
  vtb = vths(ikeep)  

  ! Integrate the radial diffusion coefficient to get the 
  !  L11, L12, etc coefficients

  norm_factor = na * vta*vta*vta*ma*ma/(TWO*qa*qa)
  nu_exp = iZERO  ! Integer

  If ( flux_cap .EQV. .true. ) Then

    K_exp = FIVEHALF
    L21 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps,       &
       vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
       Dspl_Dex,xt_c,xt_e,cmin,cmax,emin,emax, &
       num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.false.)
    
    K_exp = SEVENHALF
    L22 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps,       &
       vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
       Dspl_Dex,xt_c,xt_e,cmin,cmax,emin,emax, &
       num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.false.)

  Else

    K_exp = FIVEHALF
    L21_1 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps,       &
       vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
       Dspl_logD11,xt_c,xt_e,cmin,cmax,emin,emax, &
       num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.true.)
    K_exp = SEVENHALF
    L22_1 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps,       &
       vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
       Dspl_logD11,xt_c,xt_e,cmin,cmax,emin,emax, &
       num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.true.)

    K_exp = FIVEHALF
    L21_2 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps,       &
       vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
       Dspl_Drat2,xt_c,xt_e,cmin,cmax,emin,emax, &
       num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.false.)
    
    K_exp = SEVENHALF
    L22_2 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps,       &
       vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
       Dspl_Drat2,xt_c,xt_e,cmin,cmax,emin,emax, &
       num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.false.)

    L21 = L21_1 + L21_2
    L22 = L22_1 + L22_2

  Endif


  ! Loop over primary Sonine index (k)
  Do kval = 0_iknd,Smax

    ! Integrate QQ
    norm_factor = na * (ma*vta/qa) * (TWO*Bsq/THREE)
    nu_exp = iZERO ! Integer
    K_exp = FIVEHALF  ! Real
    Na_2k(kval+1) = energy_conv(iZERO,kval,use_quanc8,Kmin,Kmax,numKsteps, &
    vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,        &
    Dspl_Drat,xt_c,xt_e,cmin,cmax,emin,emax,   &
    num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.false.)
  EndDo

  ! Calculate PS fluxes
  QoT_PS_std = 0._rknd
  Do ispec2=1,num_species
        
    ! Specific species beta
    n_beta = dens(ispec2)
    T_beta = Temps(ispec2)
    q_beta  = charges(ispec2)
    dT_betadr = dTdrs(ispec2)
    dn_betadr = dndrs(ispec2)

    C_PS1 = (n_beta*dT_betadr+T_beta*dn_betadr)*elem_charge/(n_beta*q_beta)
    C_PS2 = dT_betadr*elem_charge/q_beta

    lmat_ind1 = ( ispec1 - 1 ) * ( Smax + 1 ) + 1
    lmat_ind2 = ( ispec2 - 1 ) * ( Smax + 1 ) + 1
    lab11 = lmat(lmat_ind1    , lmat_ind2)
    lab12 = lmat(lmat_ind1    , lmat_ind2+1)
    lab21 = lmat(lmat_ind1+1  , lmat_ind2)
    lab22 = lmat(lmat_ind1+1  , lmat_ind2+1)

    QoT_PS_std = QoT_PS_std + C_PS1*(FIVEHALF*lab11- lab21) - &
      C_PS2*(FIVEHALF*lab12 - lab22)

  EndDo


  ind_A = (ispec1-1)*3+1
  If ( flux_cap .EQV. .true. ) Then
    QoT_PS(ispec1)   = (U2/qa) * QoT_PS_std
  Else
    norm_factor = na
    nu_exp = iONE ! Integer
    K_exp = TWO  ! Real
    I_1 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps, &
      vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,        &
      Dspl_D31,xt_c,xt_e,cmin,cmax,emin,emax,      &
      num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.true.,.false.)

    K_exp = THREE  ! Real
    I_2 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps, &
      vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,        &
      Dspl_D31,xt_c,xt_e,cmin,cmax,emin,emax,      &
      num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.true.,.false.)
    QoT_PS_flow = (2._rknd/3._rknd)*(ma*Ta*elem_charge/qa)*(I_1*Avec(ind_A) + I_2*Avec(ind_A+1))
    QoT_PS(ispec1)   = (U2/qa) * (QoT_PS_std + QoT_PS_flow)
  Endif


  flow_ind1 = (ispec1-1)*(Smax+1)+1
  mono_QoT(ispec1)  = - L21*Avec(ind_A) - L22*Avec(ind_A+1)
  QoT_Ua(ispec1) = -Sum(Na_2k*Flows(flow_ind1:flow_ind1+Smax))

  ! Total QoT
  QoTs(ispec1)   = QoT_Ua(ispec1) +mono_QoT(ispec1)+QoT_PS(ispec1)

EndDo ! Species 1 loop

EndFunction calc_QoTs_SN

!-----------------------------------------------------------------------------
!+ Calculates the radial energy fluxes using the MBT method 
!-----------------------------------------------------------------------------
Function calc_QoTs_MBT(num_species,Smax,abs_Er,Temps,dens,vths,charges,      &
     masses,dTdrs,dndrs,loglambda,use_quanc8,Kmin,Kmax,numKsteps,log_interp, &
     cmin,cmax,emin,emax,xt_c,xt_e,Dspl_logD11,Dspl_D31,Dspl_Dex,     &
     num_c,num_e,kcord,keord,Avec,lmat,Flows,U2,B0,flux_cap)                 &
Result(QoTs)
!
! Description: 
! This function calculates the radial energy fluxes using the method of 
! MBT expanded to an arbitrary number of Sonine polynomials.
!
! Function arguments:
!  num_species: Total number of particle species
!  Smax: Upper index on Sonine summation
!  abs_Er: |Er| in V/cm
!  Temps: Array of species' temperatures in eV [Te, Ti1, Ti2...]
!  dens: Same as above but densitiy in [m**-3]
!  vths: Same as above but thermal velocities in m/s
!  charges: Same as above but charge in [C]
!  masses: Same as above but mass in kg
!  loglambda: Coulomb logarithm (assumed equal for each species)
!  use_quanc8: Logical switch for use of adapative integrator
!  Kmin, Kmax: Range of normalized energy over which to integrate
!  numKsteps:  Number of integration points for rectangular approx.
!  log_interp: Logical switch to activate 2D log. interp. of DKES coeffs.
!  cmin, cmax: Max, min value of cmul for each DKES matrix
!  emin, emax: Ditto for efield arrays
!  xt_c,Dspl_D11 etc: Spline knots and coefficients for DKES coeffs.
!  num_c,num_e: Number of cmul, efield values in array
!  kcord, keord: Spline orders for interpolation (cmul, efield)
!  Avec: Array of "A-type" thermodynamic forces
!  lmat: Matrix of classical friction coefficients
!  Flows: Array of parallel flow moments
!  G2: Related to the PS factor <G**2>
!  B0: B(m=0,n=0) in Boozer coordinates
!
! Output:
!  Total energy flux of each species divided by eT [m**-2]
!   [Q_e/eT_e, Q_i1/eT_i1, ...]
!
! See J. Lore's PhD thesis for expressions.  See also:
!         Maassberg, et al, PoP 15, 072504 (2009)
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     09/10/2010  Original Code.  JL
!  1.1     09/23/2011  Update to match new documentation. JL
! 
! Author(s): J. Lore 9/10/2010 - 09/23/2011
!
!
! Modules used:
Use penta_kind_mod                  ! Import rknd, iknd specifications
Use phys_const, Only :  &           ! Import physical constants
! Imported parameters
elem_charge, pi
Use penta_math_routines_mod, Only : &
! Imported functions
Gamma_aux,       &
ifactorial        

Implicit None

! Input/output                      !See above for descriptions
Integer(iknd), Intent(in)  :: num_species
Integer(iknd), Intent(in)  :: Smax
Real(rknd),    Intent(in)  :: abs_Er
Real(rknd),    Intent(in)  :: Temps(num_species)
Real(rknd),    Intent(in)  :: dens(num_species)
Real(rknd),    Intent(in)  :: vths(num_species)
Real(rknd),    Intent(in)  :: charges(num_species)
Real(rknd),    Intent(in)  :: masses(num_species)
Real(rknd),    Intent(in)  :: dTdrs(num_species)
Real(rknd),    Intent(in)  :: dndrs(num_species)
Real(rknd),    Intent(in)  :: loglambda
Logical,       Intent(in)  :: use_quanc8
Real(rknd),    Intent(in)  :: Kmin
Real(rknd),    Intent(in)  :: Kmax
Integer(iknd), Intent(in)  :: numKsteps
Logical,       Intent(in)  :: log_interp
Real(rknd),    Intent(in)  :: cmin
Real(rknd),    Intent(in)  :: cmax
Real(rknd),    Intent(in)  :: emin
Real(rknd),    Intent(in)  :: emax
Real(rknd),    Intent(in)  :: xt_c(num_c + kcord)
Real(rknd),    Intent(in)  :: xt_e(num_e + keord)
Real(rknd),    Intent(in)  :: Dspl_logD11(num_c,num_e)
Real(rknd),    Intent(in)  :: Dspl_D31(num_c,num_e)
Real(rknd),    Intent(in)  :: Dspl_Dex(num_c,num_e)
Integer(iknd), Intent(in)  :: num_c
Integer(iknd), Intent(in)  :: num_e
Integer(iknd), Intent(in)  :: kcord
Integer(iknd), Intent(in)  :: keord
Real(rknd),    Intent(in)  :: Avec(num_species*3)
Real(rknd),    Intent(in)  :: lmat(num_species*(Smax+1),num_species*(Smax+1))
Real(rknd),    Intent(in)  :: Flows(num_species*(Smax+1))
Real(rknd),    Intent(in)  :: U2
Real(rknd),    Intent(in)  :: B0
Logical,       Intent(in)  :: flux_cap
Real(rknd)                 :: QoTs(num_species)

! Local Scalars
Integer(iknd) :: ispec1,kval,   & ! Loop indices
  lval,ispec2,lmat_ind1,        &
  lmat_ind2,flow_ind2,          &
  flow_ind1,ind_A 
Integer(iknd) ::  I               ! Index used for array constructors
Integer(iknd) ::  nu_exp          ! Exponent on collision freq. for conv.
Real(rknd)    ::                & ! Primary species (a) paramaters 
  ma,Ta,vta,qa,na
Real(rknd)    :: dn_betadr,     & ! Species beta parameters (for PS flux)
  dT_betadr,n_beta,T_beta,      &
  q_beta     
Real(rknd)    :: L21,L22,L23      ! Thermal diffusion coefficients
Real(rknd)    :: norm_factor      ! Constants (wrt K) for convolution terms
Real(rknd)    :: K_exp            ! Exponent on K for convolution (ignoring 
                                  !  Maxwellian contribution)
Real(rknd)    :: lab11,lab12,   & ! Friction coefficient values for PS flux
  lab21,lab22
Real(rknd)    :: QoT_l_sum        ! Portion due to friction (sum)
Real(rknd)    :: I_1,I_2          ! Integral factors used in calculating PS flux
Real(rknd)    :: QoT_PS_std,    & ! PS fluxes
  QoT_PS_flow
Real(rknd)    :: C_PS1,C_PS2      ! PS flux coefficients

! Local arrays 
Integer(iknd) :: ikeep(num_species-1)     ! Used to index species 'b' in arrays
Real(rknd)    :: Na_2k(Smax+1)            ! First Sonine weighted coefficient
Real(rknd)    :: Ja_2(Smax+1)             ! Second Sonine weighted coefficient
Real(rknd)    :: c_l_vals(Smax+1)         ! Norm. coefficients for Sonine polys
Real(rknd)    :: friction_QoT(num_species)! Portion of Q/T due to friction
Real(rknd)    :: qb(num_species-1),  &    ! Paremeters for species 'b' /= 'a'
  nb(num_species-1), vtb(num_species-1) 
Real(rknd)    :: Ufric_QoT(num_species)   ! Portion of Q/T due to friction
Real(rknd)    :: QoT_PS(num_species)      ! PS portion of Q/T
Real(rknd)    :: mono_QoT(num_species)    ! Portion due to asymmetry
Real(rknd)    :: Uflux_QoT(num_species)   ! Portion from viscous coupling

! Local parameters                 
Integer(iknd), Parameter ::  &
iZERO  = 0_iknd,      &
iONE   = 1_iknd,      &
iTWO   = 2_iknd,      &
iTHREE = 3_iknd

Real(rknd), Parameter :: &
ZERO      = 0._rknd,  &
ONE       = 1._rknd,  &
THREEHALF = 1.5_rknd, &  
TWO       = 2._rknd,  &
FIVEHALF  = 2.5_rknd, &
THREE     = 3._rknd,  &  
FIVE      = 5._rknd,  &
SEVENHALF = 3.5_rknd  

!- End of header -------------------------------------------------------------

! Calculate terms for each species
Do ispec1 = 1_iknd, num_species

  ! Assign species 'a' parameters
  ma  = masses(ispec1)
  qa  = charges(ispec1)
  Ta  = Temps(ispec1)
  na  = dens(ispec1)
  vta = vths(ispec1)

  ! Get indices of all species b =/ a
  ikeep=(/(I,I=1,ispec1-1), (I,I=ispec1+1,num_species)/)

  ! Assign species 'b' parameters
  qb  = charges(ikeep)
  nb  = dens(ikeep)
  vtb = vths(ikeep)  

  ! Integrate the radial diffusion coefficient to get the 
  !  L11, L12, etc coefficients

  norm_factor = na*vta*vta*vta*ma*ma/(TWO*qa*qa)
  nu_exp = iZERO  ! Integer

  If ( flux_cap .EQV. .true. ) Then

    K_exp = FIVEHALF
    L21 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps,       &
       vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
       Dspl_Dex,xt_c,xt_e,cmin,cmax,emin,emax, &
       num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.false.)
    
    K_exp = SEVENHALF
    L22 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps,       &
       vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
       Dspl_Dex,xt_c,xt_e,cmin,cmax,emin,emax, &
       num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.false.)

  Else

    K_exp = FIVEHALF
    L21 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps,       &
       vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
       Dspl_logD11,xt_c,xt_e,cmin,cmax,emin,emax, &
       num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.true.)
    K_exp = SEVENHALF
    L22 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps,       &
       vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
       Dspl_logD11,xt_c,xt_e,cmin,cmax,emin,emax, &
       num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.true.)
  Endif

  ! If E_prl /= 0, then calculate L13
  ! Strictly shouldn't compare floating points, but the price is only
  !  calculating the coefficient anyway...
  ind_A = (ispec1 - 1)*3 + 1
  If ( Avec(ind_A + 2) == 0._rknd ) Then 
    L23 = ZERO
  Else
    norm_factor = na*vta*vta*ma/(TWO*B0*qa)  
    K_exp  = TWO   ! Real
    nu_exp = iZERO  ! Integer
    L23 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps,      &
      vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
      Dspl_D31,xt_c,xt_e,cmin,cmax,emin,emax, &
      num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.false.)
  Endif

  ! Calculate flux due to direct energy convolution of DKES coefficients
  mono_QoT(ispec1) = - L21*Avec(ind_A) - L22*Avec(ind_A+1) + L23*Avec(ind_A+2)

  ! Initialize sum for friction term
  friction_QoT(1:num_species) = ZERO

  ! Loop over primary Sonine index (k)
  Do kval = 0_iknd,Smax

    ! Integrate first Sonine weighted coefficient (viscous term)
    norm_factor = na * ma / qa
    nu_exp = iONE ! Integer

    K_exp = TWO ! Real
    Na_2k(kval+1) = energy_conv(iZERO,kval,use_quanc8,Kmin,Kmax,numKsteps, &
    vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,        &
    Dspl_D31,xt_c,xt_e,cmin,cmax,emin,emax,      &
    num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.false.)

    ! Loop over secondary Sonine index (l) and sum
    Do lval = 0,Smax
       
      ! Calculate the c_l coefficient
      c_l_vals(lval+1) = (0.75_rknd)*Dsqrt(pi)*  &
        ifactorial(lval)/gamma_aux(lval+2.5_rknd) ! Easy place for some optimization QQ

      ! Integrate second coefficient (viscous term)
      norm_factor = ONE / qa
      nu_exp = iZERO ! Integer        

      K_exp = TWO ! Real
      Ja_2(lval+1) = energy_conv(iZERO,lval,use_quanc8,Kmin,Kmax,numKsteps, &
        vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,      &
        Dspl_D31,xt_c,xt_e,cmin,cmax,emin,emax,    &
        num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.false.)

    EndDo ! l loop

    ! Loop over all species (including a == b)
    Do ispec2 = 1_iknd,num_species

      ! Calculate frictional contribution and sum
      lmat_ind1 = ( ispec1 - 1 ) * ( Smax + 1 ) + 1 ! j = 0:Smax is indexed below
      lmat_ind2 = ( ispec2 - 1 ) * ( Smax + 1 ) + 1 + kval

      QoT_l_sum = Sum(c_l_vals*Ja_2*lmat(lmat_ind1:lmat_ind1+Smax,lmat_ind2))

      flow_ind2 = (ispec2 - 1)*(Smax+1)+1+kval
      friction_QoT(ispec2) = friction_QoT(ispec2) + Flows(flow_ind2)*QoT_l_sum
    EndDo ! ispec2 loop
  EndDo ! k loop

  ! Calculate flux due to flows
  flow_ind1 = (ispec1 - 1)*(Smax+1)+1
  Uflux_QoT(ispec1) = -Sum(Na_2k*Flows(flow_ind1:flow_ind1+Smax))
  
  ! Calculate flux due to friction
  Ufric_QoT(ispec1) = -Sum(friction_QoT)

  ! Calculate PS fluxes
  QoT_PS_std = 0._rknd
  Do ispec2=1,num_species
        
    ! Specific species beta
    n_beta = dens(ispec2)
    T_beta = Temps(ispec2)
    q_beta  = charges(ispec2)
    dT_betadr = dTdrs(ispec2)
    dn_betadr = dndrs(ispec2)

    C_PS1 = (n_beta*dT_betadr+T_beta*dn_betadr)*elem_charge/(n_beta*q_beta)
    C_PS2 = dT_betadr*elem_charge/q_beta

    lmat_ind1 = ( ispec1 - 1 ) * ( Smax + 1 ) + 1
    lmat_ind2 = ( ispec2 - 1 ) * ( Smax + 1 ) + 1
    lab11 = lmat(lmat_ind1    , lmat_ind2)
    lab12 = lmat(lmat_ind1    , lmat_ind2+1)
    lab21 = lmat(lmat_ind1+1  , lmat_ind2)
    lab22 = lmat(lmat_ind1+1  , lmat_ind2+1)

    QoT_PS_std = QoT_PS_std + C_PS1*(FIVEHALF*lab11- lab21) - &
      C_PS2*(FIVEHALF*lab12 - lab22)

  EndDo

  ind_A = (ispec1-1)*3+1
  If ( flux_cap .EQV. .true. ) Then
    QoT_PS(ispec1)   = (U2/qa) * QoT_PS_std
  Else
    norm_factor = na
    nu_exp = iONE ! Integer
    K_exp = TWO  ! Real
    I_1 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps, &
      vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,        &
      Dspl_D31,xt_c,xt_e,cmin,cmax,emin,emax,      &
      num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.true.,.false.)

    K_exp = THREE  ! Real
    I_2 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps, &
      vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,        &
      Dspl_D31,xt_c,xt_e,cmin,cmax,emin,emax,      &
      num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.true.,.false.)
    QoT_PS_flow = (2._rknd/3._rknd)*(ma*Ta*elem_charge/qa)*(I_1*Avec(ind_A) + I_2*Avec(ind_A+1))
    QoT_PS(ispec1)   = (U2/qa) * (QoT_PS_std + QoT_PS_flow)
  Endif

  ! Total flux for this species
  QoTs(ispec1) = Uflux_QoT(ispec1)+Ufric_QoT(ispec1)+mono_QoT(ispec1)+QoT_PS(ispec1)

EndDo ! Species 1 loop

EndFunction calc_QoTs_MBT


!-----------------------------------------------------------------------------
!+ Calculates the radial particle fluxes using the MBT method 
!-----------------------------------------------------------------------------
Function calc_fluxes_MBT(num_species,Smax,abs_Er,Temps,dens,vths,charges,    &
     masses,dTdrs,dndrs,loglambda,use_quanc8,Kmin,Kmax,numKsteps,log_interp, &
     cmin,cmax,emin,emax,xt_c,xt_e,Dspl_logD11,Dspl_D31,Dspl_Dex,num_c,num_e,&
     kcord,keord,Avec,lmat,Flows,U2,B0,flux_cap) &
Result(Gammas)
!
! Description: 
! This function calculates the radial particle fluxes using the MBT method
!  for an arbitrary number of Sonine polynomial terms.
!
! Function arguments:
!  num_species: Total number of particle species
!  Smax: Upper index on Sonine summation
!  abs_Er: |Er| in V/cm
!  Temps: Array of species' temperatures in eV [Te, Ti1, Ti2...]
!  dens: Same as above but densitiy in [m**-3]
!  vths: Same as above but thermal velocities in m/s
!  charges: Same as above but charge in [C]
!  masses: Same as above but mass in kg
!  loglambda: Coulomb logarithm (assumed equal for each species)
!  use_quanc8: Logical switch for use of adapative integrator
!  Kmin, Kmax: Range of normalized energy over which to integrate
!  numKsteps:  Number of integration points for rectangular approx.
!  log_interp: Logical switch to activate 2D log. interp. of DKES coeffs.
!  cmin, cmax: Max, min value of cmul for each DKES matrix
!  emin, emax: Ditto for efield arrays
!  xt_c,Dspl_D11 etc: Spline knots and coefficients for DKES coeffs.
!  num_c,num_e: Number of cmul, efield values in array
!  kcord, keord: Spline orders for interpolation (cmul, efield)
!  Avec: Array of "A-type" thermodynamic forces
!  lmat: Matrix of classical friction coefficients
!  Flows: Array of parallel flow moments
!  G2: Related to the PS factor <G**2>
!  B0: B(m=0,n=0) in Boozer coordinates
!
! Output:
!  Radial particle flux of each species [m**-2s**-1]
!   [Gamma_e, Gamma_i1, ...]
!
! See J. Lore's PhD thesis for expressions.  See also:
!         Maassberg, et al, PoP 15, 072504 (2009)
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     09/10/2010  Original Code.  JL
!  1.1     09/23/2011  Update to match new documentation. JL
! 
! Author(s): J. Lore 9/10/2010 - 09/23/2011
!
!
! Modules used:
Use penta_kind_mod                  ! Import rknd, iknd specifications
Use phys_const, Only :  &           ! Import physical constants
! Imported parameters
elem_charge, pi
Use penta_math_routines_mod, Only : &
! Imported functions
Gamma_aux,       &
ifactorial        

Implicit None

! Input/output                      !See above for descriptions
Integer(iknd), Intent(in)  :: num_species
Integer(iknd), Intent(in)  :: Smax
Real(rknd),    Intent(in)  :: abs_Er
Real(rknd),    Intent(in)  :: Temps(num_species)
Real(rknd),    Intent(in)  :: dens(num_species)
Real(rknd),    Intent(in)  :: vths(num_species)
Real(rknd),    Intent(in)  :: charges(num_species)
Real(rknd),    Intent(in)  :: masses(num_species)
Real(rknd),    Intent(in)  :: dTdrs(num_species)
Real(rknd),    Intent(in)  :: dndrs(num_species)
Real(rknd),    Intent(in)  :: loglambda
Logical,       Intent(in)  :: use_quanc8
Real(rknd),    Intent(in)  :: Kmin
Real(rknd),    Intent(in)  :: Kmax
Integer(iknd), Intent(in)  :: numKsteps
Logical,       Intent(in)  :: log_interp
Real(rknd),    Intent(in)  :: cmin
Real(rknd),    Intent(in)  :: cmax
Real(rknd),    Intent(in)  :: emin
Real(rknd),    Intent(in)  :: emax
Real(rknd),    Intent(in)  :: xt_c(num_c + kcord)
Real(rknd),    Intent(in)  :: xt_e(num_e + keord)
Real(rknd),    Intent(in)  :: Dspl_logD11(num_c,num_e)
Real(rknd),    Intent(in)  :: Dspl_D31(num_c,num_e)
Real(rknd),    Intent(in)  :: Dspl_Dex(num_c,num_e)
Integer(iknd), Intent(in)  :: num_c
Integer(iknd), Intent(in)  :: num_e
Integer(iknd), Intent(in)  :: kcord
Integer(iknd), Intent(in)  :: keord
Real(rknd),    Intent(in)  :: Avec(num_species*3)
Real(rknd),    Intent(in)  :: lmat(num_species*(Smax+1),num_species*(Smax+1))
Real(rknd),    Intent(in)  :: Flows(num_species*(Smax+1))
Real(rknd),    Intent(in)  :: U2
Real(rknd),    Intent(in)  :: B0
logical,       Intent(in)  :: flux_cap
Real(rknd)                 :: Gammas(num_species)

! Local Scalars
Integer(iknd) :: ispec1,kval,   & ! Loop indices
  lval,ispec2,lmat_ind1,        &
  lmat_ind2,flow_ind2,          &
  flow_ind1,ind_A 
Integer(iknd) ::  I               ! Index used for array constructors
Integer(iknd) ::  nu_exp          ! Exponent on collision freq. for conv.
Real(rknd)    ::                & ! Primary species (a) paramaters 
  ma,Ta,vta,qa,na
Real(rknd)    :: dn_betadr,     & ! Species beta parameters (for PS flux)
  dT_betadr,n_beta,T_beta,      &
  q_beta     
Real(rknd)    :: L11,L12,L13      ! Thermal diffusion coefficients
Real(rknd)    ::  norm_factor     ! Constants (wrt K) for convolution terms
Real(rknd)    ::  K_exp           ! Exponent on K for convolution (ignoring 
                                  !  Maxwellian contribution)
Real(rknd)    :: lab11,lab12      ! Friction coefficient values for PS flux
Real(rknd)    :: flux_l_sum       ! Friction flux summation term
Real(rknd)    :: I_0,I_1      ! Coefficients used in calculating PS flux
Real(rknd)    :: C_PS1,C_PS2      ! PS coefficients
Real(rknd)    :: Gamma_PS_std,    & ! PS fluxes
  Gamma_PS_flow

! Local arrays 
Integer(iknd) :: ikeep(num_species-1)   ! Used to index species 'b' in arrays
Real(rknd)    :: Na_1k(Smax+1)          ! Viscous term (1st Sonine Coeff)
Real(rknd)    :: Ja_1(Smax+1)           ! Viscous term (2nd Sonine Coeff)
Real(rknd)    :: c_l_vals(Smax+1)          ! Norm. coefficients for Sonine polys
Real(rknd)    :: friction_flux(num_species)! Used in calculating friction flux
Real(rknd)    :: qb(num_species-1),  &     ! Paremeters for species 'b' /= 'a'
  nb(num_species-1), vtb(num_species-1) 
Real(rknd)    :: Ufric_Gamma(num_species)  ! Flux due to friction coupling
Real(rknd)    :: Gamma_PS(num_species)     ! PS flux
Real(rknd)    :: mono_flux(num_species)    ! Flux due to asymmetric terms
Real(rknd)    :: Uflux(num_species)        ! Flux due to flows

! Local parameters                 
Integer(iknd), Parameter ::  &
iZERO  = 0_iknd,      &
iONE   = 1_iknd,      &
iTWO   = 2_iknd,      &
iTHREE = 3_iknd

Real(rknd), Parameter :: &
ZERO      = 0._rknd,  &
ONE       = 1._rknd,  &
THREEHALF = 1.5_rknd, &  
TWO       = 2._rknd,  &
FIVEHALF  = 2.5_rknd, &
THREE     = 3._rknd,  &  
FIVE      = 5._rknd,  &
SEVENHALF = 3.5_rknd  

!- End of header -------------------------------------------------------------
Do ispec1 = 1_iknd, num_species

  ! Assign species 'a' parameters
  ma  = masses(ispec1)
  qa  = charges(ispec1)
  Ta  = Temps(ispec1)
  na  = dens(ispec1)
  vta = vths(ispec1)

  ! Get indices of all species b =/ a
  ikeep=(/(I,I=1,ispec1-1), (I,I=ispec1+1,num_species)/)

  ! Assign species 'b' parameters
  qb  = charges(ikeep)
  nb  = dens(ikeep)
  vtb = vths(ikeep)  

  norm_factor = na*vta*vta*vta*ma*ma/(TWO*qa*qa)
  nu_exp = iZERO  ! Integer

  If ( flux_cap .EQV. .true. ) Then

    K_exp = THREEHALF
    L11 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps,       &
       vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
       Dspl_Dex,xt_c,xt_e,cmin,cmax,emin,emax, &
       num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.false.)
    
    K_exp = FIVEHALF
    L12 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps,       &
       vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
       Dspl_Dex,xt_c,xt_e,cmin,cmax,emin,emax, &
       num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.false.)

  Else

    K_exp = THREEHALF
    L11 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps,       &
       vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
       Dspl_logD11,xt_c,xt_e,cmin,cmax,emin,emax, &
       num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.true.)
    K_exp = FIVEHALF
    L12 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps,       &
       vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
       Dspl_logD11,xt_c,xt_e,cmin,cmax,emin,emax, &
       num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.true.)
  Endif


  ! If E_prl /= 0, then calculate L13
  ! Strictly shouldn't compare floating points, but the price is only
  !  calculating the coefficient anyway...
  ind_A = (ispec1 - 1)*3 + 1
  If ( Avec(ind_A + 2) == 0._rknd ) Then 
    L13 = ZERO
  Else
    norm_factor = na*vta*vta*ma/(TWO*B0*qa)  
    K_exp  = ONE   ! Real
    nu_exp = iZERO  ! Integer
    L13 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps,      &
      vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
      Dspl_D31,xt_c,xt_e,cmin,cmax,emin,emax, &
      num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.false.)
  Endif

  ! Calculate flux due to direct energy convolution of DKES coefficients
  mono_flux(ispec1) = - L11*Avec(ind_A) - L12*Avec(ind_A+1) + L13*Avec(ind_A+2)

  ! Initialize sum for friction term
  friction_flux(1:num_species) = ZERO

  ! Loop over primary Sonine index (k)
  Do kval = 0_iknd,Smax

    ! Integrate first Sonine weighted coefficient (viscous term)
    norm_factor = na * ma / qa
    nu_exp = iONE ! Integer

    K_exp = ONE  ! Real
    Na_1k(kval+1) = energy_conv(iZERO,kval,use_quanc8,Kmin,Kmax,numKsteps, &
    vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,        &
    Dspl_D31,xt_c,xt_e,cmin,cmax,emin,emax,      &
    num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.false.)

    ! Loop over secondary Sonine index (l) and sum
    Do lval = 0,Smax
       
      ! Calculate the c_l coefficient
      c_l_vals(lval+1) = (0.75_rknd)*Dsqrt(pi)*  &
        ifactorial(lval)/gamma_aux(lval+2.5_rknd) ! Easy place for some optimization QQ

      ! Integrate second coefficient (viscous term)
      norm_factor = ONE / qa
      nu_exp = iZERO ! Integer        

      K_exp = ONE  ! Real
      Ja_1(lval+1) = energy_conv(iZERO,lval,use_quanc8,Kmin,Kmax,numKsteps, &
        vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,      &
        Dspl_D31,xt_c,xt_e,cmin,cmax,emin,emax,    &
        num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.false.)

    EndDo ! l loop

    ! Loop over all species (including a == b)
    Do ispec2 = 1_iknd,num_species

      lmat_ind1 = ( ispec1 - 1 ) * ( Smax + 1 ) + 1 ! j = 0:Smax is indexed below
      lmat_ind2 = ( ispec2 - 1 ) * ( Smax + 1 ) + 1 + kval

      flux_l_sum = Sum(c_l_vals*Ja_1*lmat(lmat_ind1:lmat_ind1+Smax,lmat_ind2))

      flow_ind2 = (ispec2 - 1)*(Smax+1)+1+kval
      friction_flux(ispec2) = friction_flux(ispec2) + Flows(flow_ind2)*flux_l_sum
    EndDo ! ispec2 loop
  EndDo ! k loop

  ! Calculate flux due to flows
  flow_ind1 = (ispec1 - 1)*(Smax+1)+1
  Uflux(ispec1) = -Sum(Na_1k*Flows(flow_ind1:flow_ind1+Smax))

  ! Calculate flux due to friction
  Ufric_Gamma(ispec1) = -Sum(friction_flux)

  ! Calculate PS fluxes
  Gamma_PS_std = 0._rknd
  Do ispec2=1,num_species
        
    ! Specific species beta
    n_beta = dens(ispec2)
    T_beta = Temps(ispec2)
    q_beta  = charges(ispec2)
    dT_betadr = dTdrs(ispec2)
    dn_betadr = dndrs(ispec2)

    C_PS1 = (n_beta*dT_betadr+T_beta*dn_betadr)*elem_charge/(n_beta*q_beta)
    C_PS2 = dT_betadr*elem_charge/q_beta

    lmat_ind1 = ( ispec1 - 1 ) * ( Smax + 1 ) + 1
    lmat_ind2 = ( ispec2 - 1 ) * ( Smax + 1 ) + 1
    lab11 = lmat(lmat_ind1    , lmat_ind2)
    lab12 = lmat(lmat_ind1    , lmat_ind2+1)

    Gamma_PS_std = Gamma_PS_std + C_PS1*lab11 - C_PS2*lab12

  EndDo

  ind_A = (ispec1-1)*3+1
  If ( flux_cap .EQV. .true. ) Then
    Gamma_PS(ispec1)   = (U2/qa) * Gamma_PS_std
  Else
    norm_factor = na
    nu_exp = iONE ! Integer
    K_exp = ONE  ! Real
    I_0 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps, &
      vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,        &
      Dspl_D31,xt_c,xt_e,cmin,cmax,emin,emax,      &
      num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.true.,.false.)

    K_exp = TWO  ! Real
    I_1 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps, &
      vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,        &
      Dspl_D31,xt_c,xt_e,cmin,cmax,emin,emax,      &
      num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.true.,.false.)
    Gamma_PS_flow = (2._rknd/3._rknd)*(ma*Ta*elem_charge/qa)*(I_0*Avec(ind_A) + I_1*Avec(ind_A+1))
    Gamma_PS(ispec1)   = (U2/qa) * (Gamma_PS_std + Gamma_PS_flow)
  Endif

  ! Total flux
  Gammas(ispec1) = Uflux(ispec1)+Ufric_Gamma(ispec1)+mono_flux(ispec1)+gamma_PS(ispec1)

EndDo ! Species 1 loop

EndFunction calc_fluxes_MBT


!-----------------------------------------------------------------------------
!+ Calculates the radial particle fluxes using the SN method 
!-----------------------------------------------------------------------------
Function calc_fluxes_SN(num_species,Smax,abs_Er,Temps,dens,vths,charges,masses,  &
     loglambda,use_quanc8,Kmin,Kmax,numKsteps,log_interp,cmin,cmax,emin,emax,    &
     xt_c,xt_e,Dspl_Drat,Dspl_Drat2,Dspl_Dex,Dspl_logD11,Dspl_D31,num_c,num_e,   &
     kcord,keord,Avec,Bsq,lmat,Flows,U2,dTdrs,dndrs,flux_cap)                 &
Result(Gammas)
!
! Description: 
! This function calculates the radial particle fluxes using the SN method
!  for an arbitrary number of Sonine polynomial terms.
!
! Function arguments:
!  num_species: Total number of particle species
!  Smax: Upper index on Sonine summation
!  abs_Er: |Er| in V/cm
!  Temps: Array of species' temperatures in eV [Te, Ti1, Ti2...]
!  dens: Same as above but densitiy in [m**-3]
!  vths: Same as above but thermal velocities in m/s
!  charges: Same as above but charge in [C]
!  masses: Same as above but mass in kg
!  loglambda: Coulomb logarithm (assumed equal for each species)
!  use_quanc8: Logical switch for use of adapative integrator
!  Kmin, Kmax: Range of normalized energy over which to integrate
!  numKsteps:  Number of integration points for rectangular approx.
!  log_interp: Logical switch to activate 2D log. interp. of DKES coeffs.
!  cmin, cmax: Max, min value of cmul for each DKES matrix
!  emin, emax: Ditto for efield arrays
!  xt_c,Dspl_D11 etc: Spline knots and coefficients for DKES coeffs.
!  num_c,num_e: Number of cmul, efield values in array
!  kcord, keord: Spline orders for interpolation (cmul, efield)
!  Avec: Array of "A-type" thermodynamic forces
!  Bsq: The quantity <B**2> on this surface
!  B0: B(m=0,n=0) in Boozer coordinates
!  lmat: Matrix of classical friction coefficients
!  Flows: Array of parallel flow moments
!  U2: Related to the PS factor <U**2>
!
! Output:
!  Radial particle flux of each species [m**-2s**-1]
!   [Gamma_e, Gamma_i1, ...]
!
! See J. Lore's PhD thesis for expressions.  See also:
!         Sugama, Nishimura PoP 9 4637 (2002), 
!         Sugama, Nishimura PoP 15, 042502 (2008),
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     09/10/2010  Original Code.  J
!  1.1     09/23/2011  Update to match new documentation. JL  
!
! Author(s): J. Lore 09/10/2010 - 09/23/2011
!
!
! Modules used:
Use penta_kind_mod                  ! Import rknd, iknd specifications
Use phys_const, Only :  &           ! Import physical constants
! Imported parameters
elem_charge
   
Implicit None

! Input/output                      !See above for descriptions
Integer(iknd), Intent(in)  :: num_species
Integer(iknd), Intent(in)  :: Smax
Real(rknd),    Intent(in)  :: abs_Er
Real(rknd),    Intent(in)  :: Temps(num_species)
Real(rknd),    Intent(in)  :: dens(num_species)
Real(rknd),    Intent(in)  :: vths(num_species)
Real(rknd),    Intent(in)  :: charges(num_species)
Real(rknd),    Intent(in)  :: masses(num_species)
Real(rknd),    Intent(in)  :: loglambda
Logical,       Intent(in)  :: use_quanc8
Real(rknd),    Intent(in)  :: Kmin
Real(rknd),    Intent(in)  :: Kmax
Integer(iknd), Intent(in)  :: numKsteps
Logical,       Intent(in)  :: log_interp
Real(rknd),    Intent(in)  :: cmin
Real(rknd),    Intent(in)  :: cmax
Real(rknd),    Intent(in)  :: emin
Real(rknd),    Intent(in)  :: emax
Real(rknd),    Intent(in)  :: xt_c(num_c + kcord)
Real(rknd),    Intent(in)  :: xt_e(num_e + keord)
Real(rknd),    Intent(in)  :: Dspl_Drat(num_c,num_e)
Real(rknd),    Intent(in)  :: Dspl_Drat2(num_c,num_e)
Real(rknd),    Intent(in)  :: Dspl_Dex(num_c,num_e)
Real(rknd),    Intent(in)  :: Dspl_logD11(num_c,num_e)
Real(rknd),    Intent(in)  :: Dspl_D31(num_c,num_e)
Integer(iknd), Intent(in)  :: num_c 
Integer(iknd), Intent(in)  :: num_e
Integer(iknd), Intent(in)  :: kcord
Integer(iknd), Intent(in)  :: keord
Real(rknd),    Intent(in)  :: Avec(num_species*3)
Real(rknd),    Intent(in)  :: Bsq
Real(rknd),    Intent(in)  :: lmat(num_species*(Smax+1),num_species*(Smax+1))
Real(rknd),    Intent(in)  :: Flows(num_species*(Smax+1))
Real(rknd),    Intent(in)  :: U2
Real(rknd),    Intent(in)  :: dTdrs(num_species)
Real(rknd),    Intent(in)  :: dndrs(num_species)
Logical,       Intent(in)  :: flux_cap
Real(rknd)                 :: Gammas(num_species)

! Local Scalars
Integer(iknd) :: ispec1,kval,   & ! Loop indices
  ispec2,lmat_ind1,        &
  lmat_ind2,          &
  flow_ind1,ind_A 
Integer(iknd) ::  i               ! Index used for array constructors
Integer(iknd) ::  nu_exp          ! Exponent on collision freq. for conv.
Real(rknd)    :: dn_betadr,     & ! Species beta parameters (for PS flux)
  dT_betadr,n_beta,T_beta,      &
  q_beta     
Real(rknd)    ::                & ! Primary species (a) paramaters 
  ma, Ta, vta, qa, na     
Real(rknd)    :: L11,L12  ! Thermal diffusion coefficients
Real(rknd)    :: L11_1, L11_2          ! Thermal diffusion coefficients
Real(rknd)    :: L12_1, L12_2          ! Thermal diffusion coefficients
Real(rknd)    ::  norm_factor     ! Constants (wrt K) for convolution terms
Real(rknd)    ::  K_exp           ! Exponent on K for convolution (ignoring 
                                  !  Maxwellian contribution)
Real(rknd)    :: lab11,lab12    ! Friction coefficient values for PS flux
Real(rknd)    :: I_0,I_1          ! Integral factors used in calculating PS flux
Real(rknd)    :: Gamma_PS_std,    & ! PS fluxes
  Gamma_PS_flow
Real(rknd)    :: C_PS1, C_PS2

! Local arrays 
Integer(iknd) :: ikeep(num_species-1)   ! Used to index species 'b' in arrays
Real(rknd)    :: Na_1k(Smax+1)          ! QQ
Real(rknd)    :: qb(num_species-1),  &  ! Paremeters for species 'b' /= 'a'
  nb(num_species-1), vtb(num_species-1) 
Real(rknd)    :: Gamma_PS(num_species)
Real(rknd)    :: mono_flux(num_species)
Real(rknd)    :: flux_Ua(num_species)

! Local parameters                 
Integer(iknd), Parameter ::  &
iZERO  = 0_iknd,      &
iONE   = 1_iknd,      &
iTWO   = 2_iknd,      &
iTHREE = 3_iknd

Real(rknd), Parameter :: &
THREEHALF = 1.5_rknd, &  
ONE       = 1._rknd,  &
TWO       = 2._rknd,  &
FIVEHALF  = 2.5_rknd, &  
THREE     = 3._rknd,  &
FIVE      = 5._rknd

!- End of header -------------------------------------------------------------
Do ispec1 = 1_iknd, num_species

  ! Assign species 'a' parameters
  ma  = masses(ispec1)
  qa  = charges(ispec1)
  Ta  = Temps(ispec1)
  na  = dens(ispec1)
  vta = vths(ispec1)

  ! Get indices of all species b =/ a
  ikeep=(/(I,I=1,ispec1-1), (I,I=ispec1+1,num_species)/)

  ! Assign species 'b' parameters
  qb  = charges(ikeep)
  nb  = dens(ikeep)
  vtb = vths(ikeep)  

  ! Integrate the radial diffusion coefficient to get the 
  !  L11, L12, etc coefficients

  norm_factor = na * vta*vta*vta*ma*ma/(TWO*qa*qa)
  nu_exp = iZERO  ! Integer
  If ( flux_cap .EQV. .true. ) Then

    K_exp = THREEHALF
    L11 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps,       &
       vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
       Dspl_Dex,xt_c,xt_e,cmin,cmax,emin,emax, &
       num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.false.)
    
    K_exp = FIVEHALF
    L12 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps,       &
       vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
       Dspl_Dex,xt_c,xt_e,cmin,cmax,emin,emax, &
       num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.false.)

  Else

    K_exp = THREEHALF
    L11_1 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps,       &
       vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
       Dspl_logD11,xt_c,xt_e,cmin,cmax,emin,emax, &
       num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.true.)
    K_exp = FIVEHALF
    L12_1 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps,       &
       vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
       Dspl_logD11,xt_c,xt_e,cmin,cmax,emin,emax, &
       num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.true.)

    K_exp = THREEHALF
    L11_2 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps,       &
       vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
       Dspl_Drat2,xt_c,xt_e,cmin,cmax,emin,emax, &
       num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.false.)
    
    K_exp = FIVEHALF
    L12_2 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps,       &
       vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
       Dspl_Drat2,xt_c,xt_e,cmin,cmax,emin,emax, &
       num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.false.)

    L11 = L11_1 + L11_2
    L12 = L12_1 + L12_2

  Endif

  ! Loop over primary Sonine index (k)
  Do kval = 0_iknd,Smax

    ! Integrate radial coefficient to get viscous term
    norm_factor = na * (ma*vta/qa) * (TWO*Bsq/THREE)
    nu_exp = iZERO ! Integer
    K_exp = THREEHALF  ! Real
    Na_1k(kval+1) = energy_conv(iZERO,kval,use_quanc8,Kmin,Kmax,numKsteps, &
    vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,        &
    Dspl_Drat,xt_c,xt_e,cmin,cmax,emin,emax,      &
    num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.false.)
  EndDo

  ! Calculate PS fluxes
  Gamma_PS_std = 0._rknd
  Do ispec2=1,num_species
        
    ! Specific species beta
    n_beta = dens(ispec2)
    T_beta = Temps(ispec2)
    q_beta  = charges(ispec2)
    dT_betadr = dTdrs(ispec2)
    dn_betadr = dndrs(ispec2)

    C_PS1 = (n_beta*dT_betadr+T_beta*dn_betadr)*elem_charge/(n_beta*q_beta)
    C_PS2 = dT_betadr*elem_charge/q_beta

    lmat_ind1 = ( ispec1 - 1 ) * ( Smax + 1 ) + 1
    lmat_ind2 = ( ispec2 - 1 ) * ( Smax + 1 ) + 1
    lab11 = lmat(lmat_ind1    , lmat_ind2)
    lab12 = lmat(lmat_ind1    , lmat_ind2+1)

    Gamma_PS_std = Gamma_PS_std + C_PS1*lab11 - C_PS2*lab12

  EndDo

  ind_A = (ispec1-1)*3+1
  If ( flux_cap .EQV. .true. ) Then
    Gamma_PS(ispec1)   = (U2/qa) * Gamma_PS_std
  Else
    norm_factor = na
    nu_exp = iONE ! Integer
    K_exp = ONE  ! Real
    I_0 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps, &
      vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,        &
      Dspl_D31,xt_c,xt_e,cmin,cmax,emin,emax,      &
      num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.true.,.false.)

    K_exp = TWO  ! Real
    I_1 = energy_conv(iZERO,iZERO,use_quanc8,Kmin,Kmax,numKsteps, &
      vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,        &
      Dspl_D31,xt_c,xt_e,cmin,cmax,emin,emax,      &
      num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.true.,.false.)
    Gamma_PS_flow = (2._rknd/3._rknd)*(ma*Ta*elem_charge/qa)*(I_0*Avec(ind_A) + I_1*Avec(ind_A+1))
    Gamma_PS(ispec1)   = (U2/qa) * (Gamma_PS_std + Gamma_PS_flow)
  Endif

  flow_ind1 = (ispec1-1)*(Smax+1)+1
  mono_flux(ispec1) = - L11*Avec(ind_A) - L12*Avec(ind_A+1)
  flux_Ua(ispec1) = -Sum(Na_1k*Flows(flow_ind1:flow_ind1+Smax))

  ! Total flux
  Gammas(ispec1) = flux_Ua(ispec1)+mono_flux(ispec1)+gamma_PS(ispec1)

EndDo ! Species 1 loop

EndFunction calc_fluxes_SN


!-----------------------------------------------------------------------------
!+ Calculates the parallel flow moments using the Sugama-Nishimura method 
!-----------------------------------------------------------------------------
Function calc_flows_SN(num_species,Smax,abs_Er,Temps,dens,vths,charges,  &
     masses,loglambda,B0,use_quanc8,Kmin,Kmax,numKsteps,log_interp,      &
     cmin,cmax,emin,emax,xt_c,xt_e,Dspl_Drat,Dspl_DUa,num_c,num_e,kcord, &
     keord,Avec,lmat)                                                    &
Result(Flows)
!
! Description: 
! This function calculates the parallel flow moments using only the SN
!  method expanded to an arbitrary number of Sonine polynomials.
!
! Function arguments:
!  num_species: Total number of particle species
!  Smax: Upper index on Sonine summation
!  abs_Er: |Er| in V/cm
!  Temps: Array of species' temperatures in eV [Te, Ti1, Ti2...]
!  dens: Same as above but densitiy in [m**-3]
!  vths: Same as above but thermal velocities in m/s
!  charges: Same as above but charge in [C]
!  masses: Same as above but mass in kg
!  loglambda: Coulomb logarithm (assumed equal for each species)
!  B0: B(m=0,n=0) in Boozer coordinates
!  use_quanc8: Logical switch for use of adapative integrator
!  Kmin, Kmax: Range of normalized energy over which to integrate
!  numKsteps:  Number of integration points for rectangular approx.
!  log_interp: Logical switch to activate 2D log. interp. of DKES coeffs.
!  cmin, cmax: Max, min value of cmul for each DKES matrix
!  emin, emax: Ditto for efield arrays
!  xt_c,Dspl_D11 etc: Spline knots and coefficients for DKES coeffs.
!  num_c,num_e: Number of cmul, efield values in array
!  kcord, keord: Spline orders for interpolation (cmul, efield)
!  Avec: Array of "A-type" thermodynamic forces
!  lmat: Matrix of classical friction coefficients
!
! Output:
!  Parallel flow moments <Bu_ks>/<B**2> where k is the Sonine index
!   and s is the species index.  [<Bu_0e> <Bu_1e> ...<Bu_0i> ...]/<B**2>
!
! See J. Lore's PhD thesis for expressions.  See also:
!         Sugama, Nishimura PoP 9 4637 (2002), 
!         Sugama, Nishimura PoP 15, 042502 (2008),
!                  
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     09/09/2010  Original Code.  JL
!  1.1     09/23/2011  Update to match new documentation. JL  
! 
! Author(s): J. Lore 09/09/2010 - 09/23/2011
!
!
! Modules used:
Use penta_kind_mod                  ! Import rknd, iknd specifications
Use phys_const, Only :  &           ! Import physical constants
! Imported parameters
elem_charge, pi
Use penta_math_routines_mod, Only : &
! Imported functions
Gamma_aux,       &
ifactorial,      &
idelta,          &
FINDInv,inversion_lu                             ! Subroutine to invert square matrix QQ

Implicit None

! Input/output                      !See above for descriptions
Integer(iknd), Intent(in)  :: num_species
Integer(iknd), Intent(in)  :: Smax
Real(rknd),    Intent(in)  :: abs_Er
Real(rknd),    Intent(in)  :: Temps(num_species)
Real(rknd),    Intent(in)  :: dens(num_species)
Real(rknd),    Intent(in)  :: vths(num_species)
Real(rknd),    Intent(in)  :: charges(num_species)
Real(rknd),    Intent(in)  :: masses(num_species)
Real(rknd),    Intent(in)  :: loglambda
Real(rknd),    Intent(in)  :: B0
Logical,       Intent(in)  :: use_quanc8
Real(rknd),    Intent(in)  :: Kmin
Real(rknd),    Intent(in)  :: Kmax
Integer(iknd), Intent(in)  :: numKsteps
Logical,       Intent(in)  :: log_interp
Real(rknd),    Intent(in)  :: cmin
Real(rknd),    Intent(in)  :: cmax
Real(rknd),    Intent(in)  :: emin
Real(rknd),    Intent(in)  :: emax
Real(rknd),    Intent(in)  :: xt_c(num_c + kcord)
Real(rknd),    Intent(in)  :: xt_e(num_e + keord)
Real(rknd),    Intent(in)  :: Dspl_Drat(num_c,num_e)
Real(rknd),    Intent(in)  :: Dspl_DUa(num_c,num_e)
Integer(iknd), Intent(in)  :: num_c
Integer(iknd), Intent(in)  :: num_e
Integer(iknd), Intent(in)  :: kcord
Integer(iknd), Intent(in)  :: keord
Real(rknd),    Intent(in)  :: Avec(num_species*3)
Real(rknd),    Intent(in)  :: lmat(num_species*(Smax+1),num_species*(Smax+1))
Real(rknd)                 :: Flows(num_species*(Smax+1))

! Local scalars
Integer(iknd) ::  ispec1, jval, ind_A, ind_RHS, kval, & ! Loop indices
   ind1_LHS1,   &
   ind1_LHS2, ind2_LHS2, ispec2
Integer(iknd) ::  I               ! Index used for array constructors
Integer(iknd) ::  nu_exp          ! Exponent on collision freq. for conv.
Integer(iknd) :: inv_err          ! Error flag for inversion
Real(rknd)    ::                & ! Primary species (a) paramaters 
  ma, Ta, vta, qa, na     
Real(rknd)    ::  norm_factor     ! Constants (wrt K) for convolution terms
Real(rknd)    ::  K_exp           ! Exponent on K for convolution (ignoring 
                                  !  Maxwellian contribution)
Real(rknd)    ::  RHS_1,RHS_2     ! RHS elements of flow equation

Real(rknd)    ::  LHS_tmp1        ! LHS elements of flow eq.

! Local arrays 
Integer(iknd) :: ikeep(num_species-1)   ! Used to index species 'b' in arrays
Real(rknd)    :: qb(num_species-1),  &  ! Paremeters for species 'b' /= 'a'
  nb(num_species-1), vtb(num_species-1) 
Real(rknd)    ::  &                     ! 2D arrays for the LHS of equation sys
  LHS_mat_1(num_species*(Smax+1),num_species*(Smax+1)), &
  flow_mat(num_species*(Smax+1),num_species*(Smax+1)),  & ! Flow matrix
  flow_mat_inv(num_species*(Smax+1),num_species*(Smax+1)) ! Inverted flow matrix
Real(rknd)   :: RHS(num_species*(Smax+1)) ! 1D array for RHS of eq. sys

! Local parameters                 
Integer(iknd), Parameter ::  &
iZERO = 0_iknd,      &
iONE  = 1_iknd,      &
iTWO  = 2_iknd

Real(rknd), Parameter :: &
ZERO      = 0._rknd,  &
HALF      = 0.5_rknd, &
ONE       = 1._rknd,  &
THREEHALF = 1.5_rknd, &
FIVEHALF  = 2.5_rknd, &
TWO       = 2._rknd 

!- End of header -------------------------------------------------------------

 LHS_mat_1 = 0._rknd ! Initialize LHS term for summation

! Loop over species and calculate the parallel flows
Do ispec1 = 1_iknd,num_species

  ! Assign species 'a' parameters
  ma  = masses(ispec1)
  qa  = charges(ispec1)
  Ta  = Temps(ispec1)
  na  = dens(ispec1)
  vta = vths(ispec1)

  ! Get indices of all species b =/ a
  ikeep=(/(I,I=1,ispec1-1), (I,I=ispec1+1,num_species)/)

  ! Assign species 'b' parameters
  qb  = charges(ikeep)
  nb  = dens(ikeep)
  vtb = vths(ikeep)

  ! The 'j' loop defines each equation for the current species in the system
  Do jval = 0_iknd,Smax
  
    ! Define the RHS of the equation system

    ! First RHS term
    norm_factor = na
    K_exp = THREEHALF     ! Real
    nu_exp = iZERO  ! Integer

    RHS_1 = energy_conv(jval,iZERO,use_quanc8,Kmin,Kmax,numKsteps,      &
      vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
      Dspl_Drat,xt_c,xt_e,cmin,cmax,emin,emax, &
      num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.false.)

    ! 2nd RHS term
    K_exp  = FIVEHALF    ! Real
    RHS_2 = energy_conv(jval,iZERO,use_quanc8,Kmin,Kmax,numKsteps,      &
      vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
      Dspl_Drat,xt_c,xt_e,cmin,cmax,emin,emax, &
      num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.false.)

    ! Sum the RHS terms multiplied by the forces
    ind_A = (ispec1 - 1)*3 + 1
    ind_RHS = (ispec1 - 1)*(Smax+1)+jval+1
    RHS(ind_RHS) = -RHS_1*Avec(ind_A) - RHS_2*Avec(ind_A+1) &
      + idelta(jval,0_iknd)*THREEHALF*qa*na*Avec(ind_A+2)/(ma*vta*B0) 

    ! Loop over primary Sonine index (k)
    Do kval = 0_iknd, Smax

      ! Integrate the first Ua term 
      norm_factor = na * qa/(Ta*elem_charge)
      K_exp  = THREEHALF    ! Real
      nu_exp = iZERO    ! Integer
 
      LHS_tmp1 = energy_conv(jval,kval,use_quanc8,Kmin,Kmax,numKsteps,    &
        vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
        Dspl_DUa,xt_c,xt_e,cmin,cmax,emin,emax, &
        num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.false.)
  
      ! Calculate the entire term multiplied by <U_ak>/<B**2>
      ind1_LHS1 = ( ispec1 - 1 ) * ( Smax + 1 ) + 1

      LHS_mat_1(ind1_LHS1+jval,ind1_LHS1+kval) = LHS_tmp1

      ! Loop over all species b (including b==a)
      Do ispec2 = 1,num_species

        ind1_LHS2 = ( ispec1 - 1 ) * ( Smax + 1 ) + jval +1
        ind2_LHS2 = ( ispec2 - 1 ) * ( Smax + 1 ) + kval +1 

        LHS_mat_1(ind1_LHS2,ind2_LHS2) = LHS_mat_1(ind1_LHS2,ind2_LHS2) - &
          THREEHALF*qa*lmat(ind1_LHS2,ind2_LHS2)/(ma*vta*Ta*elem_charge)
 
      EndDo ! ispec2 loop
    EndDo ! kval loop
  EndDo ! jval loop
EndDo ! Primary species loop

! Form total matrix to solve for flows
flow_mat = LHS_mat_1 

! Solve system to get flows

!Call FINDInv(flow_mat,flow_mat_inv,Size(flow_mat,1),inv_err) 
Call Inversion_lu(flow_mat,flow_mat_inv,Size(flow_mat,1),inv_err)

! Check for Inversion error
If ( inv_err /= 0 ) Then
  Write(*,*) 'Error flag returned during flow_mat inversion:',inv_err
  Stop 'Error: Exiting on inversion error in function calc_flows_SN'
EndIf

Flows = Matmul(flow_mat_inv,RHS)

EndFunction calc_flows_SN

!-----------------------------------------------------------------------------
!+ Performs the energy convolution over the local Maxwellian 
!-----------------------------------------------------------------------------
Function energy_conv(jval,kval,use_quanc8,Kmin,Kmax,numKsteps,vta,qa,qb,ma,  &
  na,nb,loglambda,abs_Er,num_species,vtb,log_interp,Dspl,xt_c,xt_e,cmin,     &
  cmax,emin,emax,nc,ne,kcord,keord,K_exp,nu_exp,norm_factor,Unity_coeff,     &
  logopt) &
Result(Coeff)
!
! Description: 
! This function convolves the specified coefficients (functions of 
! normalized energy K) over the local Maxwellian. 
!
! Some parts of this code are based on the routine DKES_coef from the original
! version of PENTA written by D. Spong (with code from van Rij?).
!
! Function arguments:
!
!  QQ update this -- talk about method
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     07/16/2010  Original Code.  JL
!  1.1     09/10/2010  Added unity_coeff option.  JL
!  1.2     09/22/2011  Added logopt option and changed to log spacing. JL
! 
! Author(s): J. Lore 7/16/2010 - 09/22/2011
!
! Modules used:
Use penta_kind_mod                    ! Import rknd, iknd specifications
Use phys_const, Only : pi             ! Import physical constants
Use penta_math_routines_mod, Only : & ! Import math routines
! Imported functions
ifactorial, &
Gamma_aux,  &
rlinspace

Implicit None

! Input/output                      ! See above for descriptions
Integer(iknd), Intent(in)  :: jval
Integer(iknd), Intent(in)  :: kval
Logical,       Intent(in)  :: use_quanc8
Real(rknd),    Intent(in)  :: Kmin
Real(rknd),    Intent(in)  :: Kmax
Integer(iknd), Intent(in)  :: numKsteps
Real(rknd),    Intent(in)  :: vta
Real(rknd),    Intent(in)  :: qa
Real(rknd),    Intent(in)  :: qb(num_species-1)
Real(rknd),    Intent(in)  :: ma
Real(rknd),    Intent(in)  :: na
Real(rknd),    Intent(in)  :: nb(num_species-1)
Real(rknd),    Intent(in)  :: loglambda
Real(rknd),    Intent(in)  :: abs_Er
Integer(iknd), Intent(in)  :: num_species
Real(rknd),    Intent(in)  :: vtb(num_species-1)
Logical,       Intent(in)  :: log_interp
Real(rknd),    Intent(in)  :: Dspl(nc,ne)
Real(rknd),    Intent(in)  :: xt_c(nc + kcord)
Real(rknd),    Intent(in)  :: xt_e(ne + keord)
Real(rknd),    Intent(in)  :: cmin
Real(rknd),    Intent(in)  :: cmax
Real(rknd),    Intent(in)  :: emin
Real(rknd),    Intent(in)  :: emax
Integer(iknd), Intent(in)  :: nc
Integer(iknd), Intent(in)  :: ne
Integer(iknd), Intent(in)  :: kcord
Integer(iknd), Intent(in)  :: keord
Real(rknd),    Intent(in)  :: K_exp
Integer(iknd), Intent(in)  :: nu_exp
Real(rknd),    Intent(in)  :: norm_factor
Logical,       Intent(in)  :: Unity_coeff
Logical,       Intent(in)  :: logopt
Real(rknd)                 :: Coeff


! Local scalars
Integer(iknd) :: mtest, iK          ! loop indices
Real(rknd) :: Gamma_k, Gamma_m, &   ! Gamma function evals
  Gamma_j
Real(rknd) :: ckm                   ! Temporary var for Sonine generation
Real(rknd) :: coeff_tmp             ! Temp. coefficient used in K loop
Real(rknd) :: K0, K1
! Local arrays
Real(rknd) :: coeffs_k(kval+1)      ! Arrays of Sonine poly coefficients
Real(rknd) :: coeffs_j(jval+1)   
Real(rknd) :: Kvals(numKsteps)      ! Energy integration values

!- End of header -------------------------------------------------------------

! Determine the coefficients of the Sonine polynomial of order kval
!  array is ordered as [K**0 K**1 ... K**kval].
Gamma_k = Gamma_aux(kval + 2.5_rknd)
Do mtest = 0,kval
  Gamma_m = Gamma_aux(mtest + 2.5_rknd)
  ckm = (-1._rknd)**mtest * Gamma_k &
    / ( Gamma_m * ifactorial( kval - mtest) * ifactorial(mtest) )
  coeffs_k(mtest+1) = ckm
EndDo

! Determine the coefficients of the Sonine polynomial of order jval
!  array is ordered as [K**0 K**1 ... K**jval].
Gamma_j = Gamma_aux(jval + 2.5_rknd)
Do mtest = 0,jval
  Gamma_m = Gamma_aux(mtest + 2.5_rknd)
  ckm = (-1._rknd)**mtest * Gamma_j &
    / ( Gamma_m * ifactorial( jval - mtest) * ifactorial(mtest) )
  coeffs_j(mtest+1) = ckm
EndDo

! Perform energy integral QQ

! Use either quanc8 or rectangular approximation
If (use_quanc8 .EQV. .true.) Then      ! quanc8
!      call quanc8(intfun, Kmin, Kmax,epsabs, epsrel,
!     1 ! !     coeff_tmp, err, nofun, int_flag)
! ! !check for integrator error
!      if(abs(int_flag) .gt. 1.e-16_rknd) then
! !   write(*,*) 'warning: quanc8 error, flag = ',int_flag
! !   write(*,*) 'May be ok, proceeding'
!      endif
  write(*,*) 'using quanc8 -- not yet implemented'
  stop
Else                      ! rectangular

  ! Log. spaced energy evaluation points
  Kvals = 10._rknd**rlinspace(log10(Kmin),log10(Kmax),numKsteps)

  ! Loop over energy values, evaluate integrand and sum 
  !  sum*dK gives integral
  ! This integrand evaluates QQ
  !  QQ replace this with vectorized integrand in future
  Coeff_tmp = 0._rknd
  Do iK = 1, numKsteps - 1
      K0 = Kvals(ik)
      K1 = Kvals(ik+1)
    Coeff_tmp =  Coeff_tmp +                                                &
      intfun(Kvals(iK),vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb, &
      log_interp,Dspl,xt_c,xt_e,cmin,cmax,emin,emax,nc,ne,kcord,keord,jval, &
      kval,K_exp,nu_exp,coeffs_j,coeffs_k,Unity_coeff,logopt)*(K1-K0)
  EndDo

EndIf ! integrator choice

Coeff = 2._rknd * Coeff_tmp * norm_factor / Dsqrt(pi)

EndFunction energy_conv

!-----------------------------------------------------------------------------
!+ Evaluates the integrand of the energy integral QQ
!-----------------------------------------------------------------------------
Function intfun(Ka,vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,       &
  log_interp,Dspl,xt_c,xt_e,cmin,cmax,emin,emax,nc,ne,kcord,keord,juse,kuse,  &
  K_exp,nu_exp,coeffs_j,coeffs_k,Unity_coeff,logopt)                                 &
Result(integrand)
!
! Description: 
! This function evaluates the integrand of the energy convolution integral
!  See QQ
!
! Talk about out of range coeffs
!
! Some parts of this code are based on the routine fun123 by D. Spong.
!
! Function arguments:
!
!  QQ update this -- talk about method
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     07/19/2010  Original Code.  JL
!  1.1     09/10/2010  Added unity_coeff option.  JL
! 
! Author(s): J. Lore 7/19/2010 - 09/10/2010
!
! Modules used:
Use penta_kind_mod                  ! Import rknd, iknd specifications
Use bspline, Only : & 
! Imported functions
dbs2vl                              ! 2D interpolation

Implicit None

! Input/output                      ! See above for descriptions
Real(rknd),    Intent(in)  :: Ka
Real(rknd),    Intent(in)  :: vta
Real(rknd),    Intent(in)  :: qa
Real(rknd),    Intent(in)  :: qb(num_species-1)
Real(rknd),    Intent(in)  :: ma
Real(rknd),    Intent(in)  :: na
Real(rknd),    Intent(in)  :: nb(num_species-1)
Real(rknd),    Intent(in)  :: loglambda
Real(rknd),    Intent(in)  :: abs_Er
Integer(iknd), Intent(in)  :: num_species
Real(rknd),    Intent(in)  :: vtb(num_species-1)
Logical,       Intent(in)  :: log_interp
Real(rknd),    Intent(in)  :: Dspl(nc,ne)
Real(rknd),    Intent(in)  :: xt_c(nc + kcord)
Real(rknd),    Intent(in)  :: xt_e(ne + keord)
Real(rknd),    Intent(in)  :: cmin
Real(rknd),    Intent(in)  :: cmax
Real(rknd),    Intent(in)  :: emin
Real(rknd),    Intent(in)  :: emax
Integer(iknd), Intent(in)  :: nc
Integer(iknd), Intent(in)  :: ne
Integer(iknd), Intent(in)  :: kcord
Integer(iknd), Intent(in)  :: keord
Integer(iknd), Intent(in)  :: juse
Integer(iknd), Intent(in)  :: kuse
Real(rknd),    Intent(in)  :: K_exp
Integer(iknd), Intent(in)  :: nu_exp
Real(rknd), Intent(in)     :: coeffs_j(juse+1)
Real(rknd), Intent(in)     :: coeffs_k(kuse+1)
Logical, Intent(in)        :: Unity_coeff
Logical,       Intent(in)  :: logopt
Real(rknd)                 :: integrand

! Local scalars
Real(rknd)    :: va              ! Test particle velocity
Real(rknd)    :: xa,xb           ! Normalized velocity (v/vth)
Real(rknd)    :: Kb              ! Normalized energy (x**2)
Real(rknd)    :: efield          ! efield parameter (|Er|/va)
Real(rknd)    :: enrm            ! normalized efield for interpolation
Real(rknd)    :: nu_perp_aa, &   ! E. dependent coll. freqs
  nu_perp_a
Real(rknd)    :: cmul_K          ! Collisionality (nu_a/va)
Real(rknd)    :: Dstar_val       ! interpolated D* value
Real(rknd)    :: kfun_tmp1,  &   ! Sonine polynomial evaluations
  kfun_tmp2, kfun, kfun2
Integer(iknd) :: ispec, mtest    ! loop indices

! Local arrays
Real(rknd)    :: nu_perp_ab(num_species - 1)      ! E. dependent coll. freqs

!- End of header -------------------------------------------------------------

! Define test species normalized velocity and real velocity
xa = Dsqrt(Ka)
va = xa * vta

! Set efield parameter (|Er|/va)
efield = abs_Er / va

! Define energy dependent collision frequency
nu_perp_aa = calc_perp_coll_freq(va,xa,qa,qa,ma,na,Ka,loglambda)

! Loop over field species and get nu_ab
Do ispec = 1, (num_species - 1)
  ! Field species norm. energy and velocity
  xb = xa*(vta/vtb(ispec))
  Kb = xb*xb
  nu_perp_ab(ispec) = calc_perp_coll_freq(va,xb,qa,qb(ispec), & 
    ma,nb(ispec),Kb,loglambda)
EndDo

! Sum over all species components
nu_perp_a = nu_perp_aa + Sum(nu_perp_ab)

! Calculate collisionality "cmul" = nu/v
cmul_K = nu_perp_a/va

If ( Unity_coeff .EQV. .true. ) Then
  Dstar_val = 1._rknd
Else
  ! Set interpolation point for logarithmic or linear interpolation 
  If ( log_interp .EQV. .true. ) Then
    efield = Dlog10(efield)
    cmul_K = Dlog10(cmul_K)
  EndIf

  ! Check for and handle out of range requests
  If ( efield < emin ) Then
    efield = emin
  ElseIf (efield > emax) Then
    efield = emax
  EndIf

  If ( cmul_K < cmin ) Then
    cmul_K = cmin
  ElseIf (cmul_K > cmax) Then
    cmul_K = cmax
  EndIf

  ! Define normalized efield for interpolation
  enrm = (efield - emin)/(emax - emin)

  ! Interpolate coefficient database
  Dstar_val = dbs2vl(cmul_K,enrm,kcord,keord,xt_c,xt_e,nc,ne,Dspl)
Endif 

! Calculate the Sonine polynomial product
kfun_tmp1 = 0._rknd
Do mtest = 0_iknd, kuse
  kfun_tmp1 = kfun_tmp1 + coeffs_k(mtest + 1)*Ka**mtest
EndDo

kfun_tmp2 = 0._rknd
Do mtest = 0_iknd, juse
  kfun_tmp2 = kfun_tmp2 + coeffs_j(mtest + 1)*Ka**mtest
EndDo

kfun = kfun_tmp1*kfun_tmp2

If (logopt .EQV. .true.) Then
  kfun2 = Dexp(Dstar_val-Ka)
Else
  kfun2 = Dexp(-Ka)*Dstar_val
Endif

integrand = kfun * kfun2 * Ka**(K_exp + 0.5_rknd) &
  * nu_perp_a**nu_exp

EndFunction intfun


!-----------------------------------------------------------------------------
!+ Returns the energy dependent collision frequency QQ
!-----------------------------------------------------------------------------
Function calc_perp_coll_freq(va,xb,qa,qb,ma,nb,Kb,loglambda)  &
Result(nu_perp_a)
!
! Description: 
!  This function returns the energy dependent collision frequency for the 
!   species a on all species b (including a == b).
!  See QQ
!
! Talk about derf? QQ
!
! INPUT QQ
!
! Function arguments:
!
!  QQ update this -- talk about method
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     07/19/2010  Original Code.  JL
! 
! Author(s): J. Lore 7/19/2010 
!
! Modules used:
Use penta_kind_mod                  ! Import rknd, iknd specifications
Use phys_const, Only :  &           ! Import physical constants
! Imported parameters
pi, eps0

Implicit None

! Input/output                      ! See above for descriptions
Real(rknd),    Intent(in)  :: va
Real(rknd),    Intent(in)  :: xb
Real(rknd),    Intent(in)  :: qa
Real(rknd),    Intent(in)  :: qb
Real(rknd),    Intent(in)  :: ma
Real(rknd),    Intent(in)  :: nb
Real(rknd),    Intent(in)  :: Kb
Real(rknd),    Intent(in)  :: loglambda
Real(rknd)                 :: nu_perp_a

! Local scalars
Real(rknd) :: va3          ! va**3
Real(rknd) :: qa2, qb2     ! qa**2, qb**2
Real(rknd) :: nu_ab        ! ref. coll. freq
Real(rknd) :: H_ab         ! related to rosenbluth potential

! Local arrays

! Local parameters
Integer(iknd), Parameter :: iTWO   = 2_iknd
Integer(iknd), Parameter :: iTHREE = 3_iknd
Real(rknd),    Parameter :: ONE    = 1._rknd
Real(rknd),    Parameter :: TWO    = 2._rknd
Real(rknd),    Parameter :: FOUR   = 4_iknd

! External Functions  (For pgf90, uncomment the following line)
!Real(rknd), External :: Derf ! QQ

!- End of header -------------------------------------------------------------

! Predefine powers for 'speed'
va3 = va*va*va
qa2 = qa*qa
qb2 = qb*qb

! Calculate reference collision freq. of a on b 
nu_ab=nb*qa2*qb2*loglambda/(ma*ma*va3*FOUR*pi*eps0*eps0)
! Calculate potential
H_ab=(one-one/(TWO*Kb))*Derf(xb) + Dexp(-Kb)/(xb*Dsqrt(pi))
! Total Collision frequency
nu_perp_a = H_ab * nu_ab

EndFunction calc_perp_coll_freq

!-----------------------------------------------------------------------------
!+ Calculates the parallel flow moments using the Taguchi method (QQ)
!-----------------------------------------------------------------------------
Function calc_flows_T(num_species,Smax,abs_Er,Temps,dens,vths,charges,  &
     masses,loglambda,B0,use_quanc8,Kmin,Kmax,numKsteps,log_interp,     &
     cmin,cmax,emin,emax,xt_c,xt_e,Dspl_D31,     &
     Dspl_logD33,num_c,num_e,kcord,   &
     keord,Avec,Bsq,lmat)                                               &
Result(Flows)
!
! Description: 
! This function calculates the parallel flow moments using the method of 
! Taguchi expanded to an arbitrary number of Sonine polynomials.
!
! Function arguments:
!
!  QQ update this -- talk about method
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     07/16/2010  Original Code.  JL
!  1.1     09/23/2011  Updated to match new documentation. JL
! 
! Author(s): J. Lore 7/16/2010 - 09/23/2011
!
!
! Modules used:
Use penta_kind_mod                  ! Import rknd, iknd specifications
Use phys_const, Only :  &           ! Import physical constants
! Imported parameters
elem_charge, pi
Use penta_math_routines_mod, Only : &
! Imported functions
Gamma_aux,       &
ifactorial,      &
idelta,          &
FINDInv,inversion_lu                             ! Subroutine to invert square matrix QQ

Implicit None

! Input/output                      !See above for descriptions
Integer(iknd), Intent(in)  :: num_species
Integer(iknd), Intent(in)  :: Smax
Real(rknd),    Intent(in)  :: abs_Er
Real(rknd),    Intent(in)  :: Temps(num_species)
Real(rknd),    Intent(in)  :: dens(num_species)
Real(rknd),    Intent(in)  :: vths(num_species)
Real(rknd),    Intent(in)  :: charges(num_species)
Real(rknd),    Intent(in)  :: masses(num_species)
Real(rknd),    Intent(in)  :: loglambda
Real(rknd),    Intent(in)  :: B0
Logical,       Intent(in)  :: use_quanc8
Real(rknd),    Intent(in)  :: Kmin
Real(rknd),    Intent(in)  :: Kmax
Integer(iknd), Intent(in)  :: numKsteps
Logical,       Intent(in)  :: log_interp
Real(rknd),    Intent(in)  :: cmin
Real(rknd),    Intent(in)  :: cmax
Real(rknd),    Intent(in)  :: emin
Real(rknd),    Intent(in)  :: emax
Real(rknd),    Intent(in)  :: xt_c(num_c + kcord)
Real(rknd),    Intent(in)  :: xt_e(num_e + keord)
Real(rknd),    Intent(in)  :: Dspl_D31(num_c,num_e)
Real(rknd),    Intent(in)  :: Dspl_logD33(num_c,num_e)
Integer(iknd), Intent(in)  :: num_c
Integer(iknd), Intent(in)  :: num_e
Integer(iknd), Intent(in)  :: kcord
Integer(iknd), Intent(in)  :: keord
Real(rknd),    Intent(in)  :: Avec(num_species*3)
Real(rknd),    Intent(in)  :: Bsq
Real(rknd),    Intent(in)  :: lmat(num_species*(Smax+1),num_species*(Smax+1))
Real(rknd)                 :: Flows(num_species*(Smax+1))

! Local scalars
Integer(iknd) ::  ispec1, jval, ind_A, ind_RHS, kval, & ! Loop indices
  mval, ind1_LHS1, lmat_ind1, lmat_ind2,   &
   ind1_LHS2, ind2_LHS2, ispec2
Integer(iknd) ::  I               ! Index used for array constructors
Integer(iknd) ::  nu_exp          ! Exponent on collision freq. for conv.
Integer(iknd) :: inv_err          ! Error flag for FINDInv
Real(rknd)    ::                & ! Primary species (a) paramaters 
  ma, Ta, vta, qa, na     
Real(rknd)    ::  norm_factor     ! Constants (wrt K) for convolution terms
Real(rknd)    ::  K_exp           ! Exponent on K for convolution (ignoring 
                                  !  Maxwellian contribution)
Real(rknd)    ::  RHS_1,RHS_2,  & ! RHS elements of flow equation
  RHS_3
Real(rknd)    ::  LHS_tmp1,     & ! LHS elements of flow eq.
  LHS2_sum, LHS_tmp2
Real(rknd)    ::  c_k, c_m        ! Norm. coefficients for Sonine polys

! Local arrays 
Integer(iknd) :: ikeep(num_species-1)   ! Used to index species 'b' in arrays
Real(rknd)    :: qb(num_species-1),  &  ! Paremeters for species 'b' /= 'a'
  nb(num_species-1), vtb(num_species-1) 
Real(rknd)    ::  &                     ! 2D arrays for the LHS of equation sys
  LHS_mat_1(num_species*(Smax+1),num_species*(Smax+1)), &
  LHS_mat_2(num_species*(Smax+1),num_species*(Smax+1)), &  
  flow_mat(num_species*(Smax+1),num_species*(Smax+1)),  & ! Flow matrix
  flow_mat_inv(num_species*(Smax+1),num_species*(Smax+1)) ! Inverted flow matrix
Real(rknd)   :: RHS(num_species*(Smax+1)) ! 1D array for RHS of eq. sys

! Local parameters                 
Integer(iknd), Parameter ::  &
iZERO = 0_iknd,      &
iONE  = 1_iknd,      &
iTWO  = 2_iknd

Real(rknd), Parameter :: &
HALF = 0.5_rknd,    &
ZERO = 0._rknd,     &
ONE  = 1._rknd,     &
TWO  = 2._rknd 

!- End of header -------------------------------------------------------------

LHS_mat_1 = 0._rknd  !initialize matrix for summation

! Loop over species and calculate the parallel flows
Do ispec1 = 1_iknd,num_species

  ! Assign species 'a' parameters
  ma  = masses(ispec1)
  qa  = charges(ispec1)
  Ta  = Temps(ispec1)
  na  = dens(ispec1)
  vta = vths(ispec1)

  ! Get indices of all species b =/ a
  ikeep=(/(I,I=1,ispec1-1), (I,I=ispec1+1,num_species)/)

  ! Assign species 'b' parameters
  qb  = charges(ikeep)
  nb  = dens(ikeep)
  vtb = vths(ikeep)

  ! The 'j' loop defines each equation for the current species in the system
  Do jval = 0_iknd,Smax
  
    ! Define the RHS of the equation system

    ! First RHS term
    norm_factor = na
    K_exp = ONE     ! Real
    nu_exp = iZERO  ! Integer
    RHS_1 = energy_conv(jval,iZERO,use_quanc8,Kmin,Kmax,numKsteps,      &
      vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
      Dspl_D31,xt_c,xt_e,cmin,cmax,emin,emax, &
      num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.false.)

    ! 2nd RHS term
    K_exp  = TWO    ! Real
    RHS_2 = energy_conv(jval,iZERO,use_quanc8,Kmin,Kmax,numKsteps,      &
      vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
      Dspl_D31,xt_c,xt_e,cmin,cmax,emin,emax, &
      num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.false.)

    ! If E_prl /= 0, then calculate 3rd RHS Term
    ! Strictly shouldn't compare floating points, but the price is simply
    !  calculating the coefficient anyway... QQ
    ind_A = (ispec1 - 1)*3 + 1
    If ( Avec(ind_A + 2) == 0._rknd ) Then 
      RHS_3 = ZERO
    Else
      norm_factor = na * qa/(B0*ma*vta)
      K_exp  = HALF    ! Real
      nu_exp = iZERO  ! Integer
      RHS_3 = energy_conv(jval,iZERO,use_quanc8,Kmin,Kmax,numKsteps,      &
        vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
        Dspl_logD33,xt_c,xt_e,cmin,cmax,emin,emax, &
        num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.true.)
    Endif

    ! Sum the RHS terms multiplied by the forces (QQ check mags)
    !  QQ -- some question over sign of 3rd term
    ind_RHS = (ispec1 - 1)*(Smax+1)+jval+1
    RHS(ind_RHS) = - RHS_1*Avec(ind_A) - RHS_2*Avec(ind_A+1) &
      + RHS_3*Avec(ind_A+2)

    ! Loop over primary Sonine index (k)
    Do kval = 0_iknd, Smax

      ! Calculate the c_k coefficient
      c_k = (0.75_rknd)*Dsqrt(pi)*ifactorial(kval)/gamma_aux(kval+2.5_rknd) ! Easy place for some optimization QQ

      ! Integrate the ||nu*sqrt(K)*L_k*D33^*||_j term
      norm_factor = na * qa/(vta*Ta*elem_charge)
      K_exp  = HALF    ! Real
      nu_exp = iONE    ! Integer
 
      LHS_tmp1 = energy_conv(jval,kval,use_quanc8,Kmin,Kmax,numKsteps,    &
        vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
        Dspl_logD33,xt_c,xt_e,cmin,cmax,emin,emax, &
        num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.true.)

      ! QQ only have to do diagonals above

      ! Calculate the entire term multiplied by <U_ak>/<B**2>
      ind1_LHS1 = ( ispec1 - 1 ) * ( Smax + 1 ) + 1

      LHS_mat_1(ind1_LHS1+jval,ind1_LHS1+kval) =   &
        idelta(jval,kval)*na*qa*Bsq/(Ta*elem_charge*c_k) - LHS_tmp1

      ! Loop over all species b (including b==a)
      Do ispec2 = 1,num_species
        ! Loop over the secondary Sonine index and sum over m
        LHS2_sum = 0._rknd
        Do mval = 0_iknd,Smax
          ! Calculate the c_m coefficient
          c_m = (0.75_rknd)*Dsqrt(pi)*ifactorial(mval)/gamma_aux(mval+2.5_rknd)
          
          ! Integrate the ||sqrt(K)*L_m*D33^*||_j Term
          norm_factor = qa/(Ta*elem_charge*ma*vta)
          K_exp  = HALF    ! Real
          nu_exp = iZERO    ! Integer
          LHS_tmp2 = energy_conv(jval,mval,use_quanc8,Kmin,Kmax,numKsteps,    &
            vta,qa,qb,ma,na,nb,loglambda,abs_Er,num_species,vtb,log_interp,   &
            Dspl_logD33,xt_c,xt_e,cmin,cmax,emin,emax, &
            num_c,num_e,kcord,keord,K_exp,nu_exp,norm_factor,.false.,.true.)

          ! Check for B0 factor here QQ and compare to matlab
          lmat_ind1 = ( ispec1 - 1 ) * ( Smax + 1 ) + 1 + mval
          lmat_ind2 = ( ispec2 - 1 ) * ( Smax + 1 ) + 1 + kval

          LHS2_sum = LHS2_sum        &
            - c_m*LHS_tmp2*lmat(lmat_ind1,lmat_ind2)
 
        EndDo ! mval loop
        
        ! Assign elements of LHS_mat_2 QQ

        ind1_LHS2 = ( ispec1 - 1 ) * ( Smax + 1 ) + jval +1
        ind2_LHS2 = ( ispec2 - 1 ) * ( Smax + 1 ) + kval +1 
 
        LHS_mat_2(ind1_LHS2,ind2_LHS2) = LHS2_sum

      EndDo ! ispec2 loop
    EndDo ! kval loop
  EndDo ! jval loop
EndDo ! Primary species loop

! Form total matrix to solve for flows
flow_mat = LHS_mat_1 + LHS_mat_2

! Solve system to get flows

!Call FINDInv(flow_mat,flow_mat_inv,Size(flow_mat,1),inv_err) !QQ
Call Inversion_lu(flow_mat,flow_mat_inv,Size(flow_mat,1),inv_err)

! Check for Inversion error
If ( inv_err /= 0 ) Then
  Write(*,*) 'Error flag returned during flow_mat inversion:',inv_err
  Stop 'Error: Exiting on inversion error in function calc_flows_T'
EndIf

Flows = Matmul(flow_mat_inv,RHS)

EndFunction calc_flows_T



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
!
!
! Modules used:
Use penta_kind_mod                ! Import rknd, iknd specifications  
Use penta_math_routines_mod, Only : &
! Imported functions
ifactorial, &
Gamma_aux

Implicit None

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
!
!
! Modules used:
Use penta_kind_mod                ! Import rknd, iknd specifications
Use penta_math_routines_mod, Only : &
! Imported functions
ifactorial, &
Gamma_aux

Implicit None

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
!
!
! Modules used:
Use penta_kind_mod                ! Import rknd, iknd specifications                            

Implicit None

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
!
! Modules used:
Use penta_kind_mod                ! Import rknd, iknd specifications                            
Use penta_math_routines_mod, Only : &
! Imported functions
  ifactorial

Implicit None

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
! Modules used:
Use penta_kind_mod                ! Import rknd, iknd specifications                            
Use penta_math_routines_mod, Only : &
! Imported functions
ifactorial
Implicit None

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
! Modules used:
Use penta_kind_mod                ! Import rknd, iknd specifications
Use penta_math_routines_mod, Only : &
! Imported functions
Gamma_aux

Implicit None

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
!+ Returns gamma_e-sum(Z*gamma_i) at a given Er using polynomial interpolation
!-----------------------------------------------------------------------------
Function rad_flux(Er_test,diff_qg,Er_fit,num_pts) &
Result(rad_flux_out)
!
! Description: 
!   This function is used in conjunction with 'zeroin' to find
!    the solution to gamma_e=sum(Z*gamma_i).  Returns difference
!    in fluxes at Er_test, using spline parameters from Er_fit_pass.
! Function arguments:
!
!  QQ update this -- talk about method
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     07/xx/2009  Original Code.  JL
!  1.1     09/01/2010  Updated for PENTA3, added arguments over module passing. JL
!
! Author(s): J. Lore 07/2009 - 9/1/2010 
!
!
! Modules used:
Use penta_kind_mod                   ! Import rknd, iknd specifications
Use penta_math_routines_mod, Only : &
! Imported routines
pol_int                              ! Polynomial interpolation routine

Implicit None

! Input/output                      !See above for descriptions
Real(rknd), Intent(in)    :: Er_test
Real(rknd), Intent(in)    :: diff_qg(num_pts)
Real(rknd), Intent(in)    :: Er_fit(num_pts)
Integer(iknd), Intent(in) :: num_pts
Real(rknd)                :: rad_flux_out

! Local scalars
Real(rknd) :: Er_tmp, err

!- End of header -------------------------------------------------------------

! Use a polynomial fit to evaluate function at Er
Call pol_int(Er_fit,diff_qg,Size(Er_fit),Er_test,Er_tmp,err)
rad_flux_out = Er_tmp

EndFunction rad_flux

!-----------------------------------------------------------------------------
!+ Finds the zero of a function over a given interval
!-----------------------------------------------------------------------------
Function zeroin(ax,bx,f,tol,diff_qg,Er_fit,num_pts) &
Result(zeroin_out)
!
! Description: 
!  a zero of the function  f(x)  is computed in the interval ax,bx
!
!  input..
!
!  ax     left endpoint of initial interval
!  bx     right endpoint of initial interval
!  f      function subprogram which evaluates f(x) for any x in
!         the interval  ax,bx
!  tol    desired length of the interval of uncertainty of the
!         final result ( .ge. 0.0)
!
!
!  output..
!
!  zeroin abcissa approximating a zero of  f  in the interval ax,bx
!
!
!      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
!  without  a  check.  zeroin  returns a zero  x  in the given interval
!  ax,bx  to within a tolerance  4*macheps*abs(x) + tol, where macheps
!  is the relative machine precision.
!      this function subprogram is a slightly  modified  translation  of
!  the algol 60 procedure  zero  given in  richard brent, algorithms for
!  minimization without derivatives, prentice - hall, inc. (1973).
!
! Function arguments:
!
!  QQ update this -- talk about method
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     xx/xx/2007  Original Code taken from Don's PENTA.
!  1.1     07/xx/2009  Modified for kind precision, iteration limit added. JL
!  1.2     09/01/2010  Updated for free form, etc for PENTA3. JL
!  1.3     09/01/2010  Added inputs, now for specific use. JL
!
! Author(s): Unknown QQ, J. Lore 5/28/2010 
!
!

! Modules used:
Use penta_kind_mod                ! Import rknd, iknd specifications

Implicit None

! Input/output                      !See above for descriptions
Real(rknd),    Intent(in)  :: ax
Real(rknd),    Intent(in)  :: bx
Real(rknd),    Intent(in)  :: diff_qg(num_pts)
Real(rknd),    Intent(in)  :: Er_fit(num_pts)
Integer(iknd), Intent(in)  :: num_pts
Real(rknd),    Intent(in)  :: tol
Real(rknd)                 :: zeroin_out

! Local scalars
Integer(iknd)        :: iter
Real(rknd) :: a, b, c, d, e, eps, fa, fb, fc, tol1,  &
  xm, p, q, r, s

! Local parameters
Integer(iknd), Parameter :: max_iter = 1000_iknd

! External functions
Real(rknd), External :: f

!- End of header -------------------------------------------------------------

! Compute eps, the relative machine precision
iter=0_iknd
eps = Epsilon(1._rknd)
tol1 = 1._rknd + eps

! QQ check need for line numbers

! Initialization
a = ax
b = bx
fa = f(a,diff_qg,Er_fit,num_pts)
fb = f(b,diff_qg,Er_fit,num_pts)
 
! Begin step
20 c = a
fc = fa
d = b - a
e = d
30 if (dabs(fc) .ge. dabs(fb)) go to 40
a = b
b = c
c = a
fa = fb
fb = fc
fc = fa

! Convergence test
40 tol1 = 2._rknd*eps*dabs(b) + 0.5_rknd*tol
xm = 0.5_rknd*(c - b)
iter=iter+1_iknd
if (iter .gt. max_iter) then
  write(*,*) 'Maximum iterations exceeded in function zeroin:', max_iter
  Stop
  Endif
if (dabs(xm) .le. tol1) go to 90
if (fb .eq. 0._rknd) go to 90
 
! Is bisection necessary
if (dabs(e) .lt. tol1) go to 70
if (dabs(fa) .le. dabs(fb)) go to 70

! is quadratic interpolation possible
if (a .ne. c) go to 50

! linear interpolation
s = fb/fa
p = 2._rknd*xm*s
q = 1._rknd - s
go to 60

! inverse quadratic interpolation
50 q = fa/fc
r = fb/fc
s = fb/fa
p = s*(2._rknd*xm*q*(q - r) - (b - a)*(r - 1._rknd))
q = (q - 1._rknd)*(r - 1._rknd)*(s - 1._rknd)
 
! adjust signs
60 if (p .gt. 0._rknd) q = -q
p = dabs(p)

! is interpolation acceptable
if ((2._rknd*p) .ge. (3._rknd*xm*q - dabs(tol1*q))) go to 70
if (p .ge. dabs(0.5_rknd*e*q)) go to 70
e = d
d = p/q
Go to 80

! bisection
70 d = xm
e = d

!complete step
80 a = b
fa = fb
if (dabs(d) .gt. tol1) b = b + d
if (dabs(d) .le. tol1) b = b + dsign(tol1, xm)
fb = f(b,diff_qg,Er_fit,num_pts)
if ((fb*(fc/dabs(fc))) .gt. 0._rknd) go to 20
go to 30

! done
90 zeroin_out = b
return
end function zeroin

End Module penta_functions_mod
!- End of module header -------------------------------------------------------
