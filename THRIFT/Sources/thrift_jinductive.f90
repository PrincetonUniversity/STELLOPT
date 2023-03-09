!-----------------------------------------------------------------------
!     Subroutine:    thrift_jinductive
!     Authors:       L. van Ham
!     Date:          01/23/2023
!     Description:   This subroutine updates the plasma response to
!                    source currents.  
!-----------------------------------------------------------------------
      SUBROUTINE thrift_jinductive
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE thrift_runtime
      USE thrift_vars
      USE thrift_funcs
      USE thrift_equil
      USE thrift_profiles_mod
      USE stel_tools
      USE EZspline_obj
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: i,j,prevtimestep, ier
      REAL(rprec) :: rho,s,drho,dt,mytime,s11,s12,etapara,pprime,&
                     temp1,temp2,source_axis,source_edge,Aminor,Rmajor,&
                     A_temp,B_temp,C_temp,D_temp,&
                     BP_temp, CP_temp, DP_temp,&
                     a1,a2,a3,a4
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: j_full, jsource_full, jplasma_full

!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
!     PRELIMINARIES
!----------------------------------------------------------------------

      ! If at zero beta, copy previous value of JPLASMA onto this timestep and skip
      IF (eq_beta == 0) THEN
         IF (mytimestep /= 1) THRIFT_JPLASMA(:,mytimestep) = THRIFT_JPLASMA(:,mytimestep-1)
         GOTO 1000 
      END IF

      ! Allocate
      ALLOCATE(j_full(nrho+2),jplasma_full(nrho+2),jsource_full(nrho+2))

      A_temp  = 0; B_temp  = 0; C_temp  = 0; D_temp  = 0
                   BP_temp = 0; CP_temp = 0; DP_temp = 0
      a1      = 0; a2      = 0; a3      = 0; a4      = 0
   
      j_full = 0; 
      jplasma_full = 0;
      jsource_full = 0; 


      ! Extrapolate source current density to magnetic axis and edge
      CALL extrapolate_arr(THRIFT_JSOURCE(:,mytimestep), jsource_full)     
      IF (lverbj) THEN
         WRITE(6,*)'==============================================================================='
         WRITE(6,*)' CALCULATING MAGNETIC VARIABLES'
         WRITE(6,*) ' RHO    DV/DPHI        <B>      <B^2>     RMAJOR     AMINOR        S11 '
      END IF

      !! ONLY NECESSARY AT CURRENT TIMESTEP !! REWRITE PLS 
      ! Store magnetic variables; both preceding and current timestep vals are necessary
      THRIFT_PHIEDGE(1,mytimestep) = eq_phiedge
      DO i = 1, nrho+2
         rho = THRIFT_RHOFULL(i)
         s = rho*rho
         ier = 0
         CALL EZspline_interp(vp_spl, rho, THRIFT_VP(i,mytimestep), ier)
         CALL get_equil_Bav(s, THRIFT_BAV(i,mytimestep), THRIFT_BSQAV(i,mytimestep), ier)
         CALL get_equil_sus(s, THRIFT_S11(i,mytimestep), temp1, temp1, temp1, ier)
         CALL get_equil_Rmajor(s, THRIFT_RMAJOR(i,mytimestep), temp1, THRIFT_AMINOR(i,mytimestep),ier)
         IF (lverbj) WRITE(6,'(F5.3,5(1X,F10.6),1X,ES10.3)') &
           rho, THRIFT_VP(i,mytimestep), THRIFT_BAV(i,mytimestep), THRIFT_BSQAV(i,mytimestep), &
           ABS(THRIFT_S11(i,mytimestep)), THRIFT_RMAJOR(i,mytimestep), THRIFT_AMINOR(i,mytimestep)
      END DO
      THRIFT_S11 = ABS(THRIFT_S11)

      ! If mytimestep = 1 ITOT=0 and continue to next iteration
      IF (mytimestep==1) THEN
         jplasma_full = -jsource_full
         THRIFT_JPLASMA(:,mytimestep) = jplasma_full(2:nrho+1)
         GOTO 1000 ! skip iteration. 
      END IF

      mytime = THRIFT_T(mytimestep) ! mytime = current sim time
      prevtimestep = mytimestep-1   ! previous time step index
      dt = THRIFT_T(mytimestep)-THRIFT_T(prevtimestep) ! dt = delta t this iter
      IF (mytime>tmax.and.THRIFT_T(prevtimestep)<=tmax.and.(nsubsteps==1)) WRITE(6,*) &
         '! THRIFT has exceeded end time of profiles file. Proceeding with profiles at t=tmax !' 
         
!----------------------------------------------------------------------
!     CALCULATING COEFFICIENTS ABCD
!----------------------------------------------------------------------
!     Calculate ABCD (everything evaluated at rho_j)
!     > A(j) = S11/(4*rho*Phi_edge)
!     > B(j) = 2*etapara*dV/dPhi*<B^2>/mu_0
!     > C(j) = 2*etapara*dV/dPhi*dp/drho
!     > D(j) = -2*etapara*dV/dPhi*<Js.B>
!----------------------------------------------------------------------
      IF (lverbj) THEN
         WRITE(6,*)'==============================================================================='
         WRITE(6,*)' CALCULATING COEFFICIENTS A,B,C,D'
         WRITE(6,*) ' RHO  ETAPARA     DV/DPHI      DP/DRHO     <J.B>      BSQAV        S11'
      END IF

      DO i = 1, nrho+2
         rho = THRIFT_RHOFULL(i)
         CALL get_prof_etapara(MIN(rho,THRIFT_RHO(nrho)),mytime,etapara)
         CALL get_prof_pprime(rho,mytime,pprime)
         temp1 = 2*etapara*THRIFT_VP(i,mytimestep) ! temp1 <- 2 eta dV/dPhi 
         IF (i > 1) &
            THRIFT_COEFF_A(i,mytimestep) = THRIFT_S11(i,mytimestep)/(4*rho*THRIFT_PHIEDGE(1,mytimestep))
         THRIFT_COEFF_B(i,mytimestep) = temp1*THRIFT_BSQAV(i,mytimestep)/mu0
         THRIFT_COEFF_C(i,mytimestep) = temp1*pprime               
         THRIFT_COEFF_D(i,mytimestep) = -temp1*jsource_full(i)*THRIFT_BAV(i,mytimestep)  
         IF (lverbj) WRITE(6,'(F5.3,6(1X,ES10.3))') &
           rho, etapara, THRIFT_VP(i,mytimestep), pprime, jsource_full(i)*THRIFT_BAV(i,mytimestep), THRIFT_BSQAV(i,mytimestep), THRIFT_S11(i,mytimestep)
      END DO
!----------------------------------------------------------------------
!     Calculate derivatives of ABCD here (see deriv1_rho_o2)
!----------------------------------------------------------------------
      
      CALL deriv1_rho_o2(THRIFT_COEFF_B(:,mytimestep), THRIFT_COEFF_BP(:,mytimestep))
      CALL deriv1_rho_o2(THRIFT_COEFF_C(:,mytimestep), THRIFT_COEFF_CP(:,mytimestep))
      CALL deriv1_rho_o2(THRIFT_COEFF_D(:,mytimestep), THRIFT_COEFF_DP(:,mytimestep))

      IF (lverbj) THEN
         WRITE(6,*)'==============================================================================='
         WRITE(6,*)' COEFFICIENTS ABCD'
         WRITE(6,*)' RHO         A         B          C          D       BDER       CDER       DDER'
         WRITE(6,*)''
         DO i = 1, nrho+2
            WRITE(6,'(F5.3, 1X, 7(ES10.2,1X))') THRIFT_RHOFULL(i), &
            THRIFT_COEFF_A(i,mytimestep), THRIFT_COEFF_B(i,mytimestep),& 
            THRIFT_COEFF_C(i,mytimestep), THRIFT_COEFF_D(i,mytimestep),&
            THRIFT_COEFF_BP(i,mytimestep), THRIFT_COEFF_CP(i,mytimestep),&
            THRIFT_COEFF_DP(i,mytimestep)
         END DO
      END IF
!----------------------------------------------------------------------
!     CALCULATING COEFFICIENTS ALPHA
!----------------------------------------------------------------------
!     Calculate alpha_1234 (everything evaluated at rho_j)
!     > a1(j) = A*dD/drho
!     > a2(j) = A*((dB/drho)/rho - B/rho^2 + dC/drho)
!     > a3(j) = A*(dB/drho + B/rho + C)
!     > a4(j) = A*B
!
!     NOTE: (Derivatives of) ABCD are required on the full [0,1] grid,
!     so their arrays are of size nrho+2. Alphas are only needed on the
!     (0,1) (THRIFT_RHO) grid, so have arrays of size nrho. Hence indices
!     for the RHS of these equations are shifted by 1 (j = i+1)
!
!         0                                1        rho
!         |  |     |     | ... |     |     |  |     ABCD:  j
!        j=1 2     3                 n    n+1 n+2
!            |     |     | ... |     |     |        a1234: i
!           i=1    2                n-1    n        => j = i + 1
!
!----------------------------------------------------------------------
      IF (lverbj) THEN
         WRITE(6,*)'==============================================================================='
         WRITE(6,*)' ALPHAS'
         WRITE(6,*)'RHO       ALPHA 1        ALPHA 2        ALPHA 3        ALPHA 4'
         WRITE(6,*)''
      END IF

      DO i = 1, nrho
         rho = THRIFT_RHO(i)
         j = i + 1
         ! Temp variables for legibility
         A_temp = THRIFT_COEFF_A(j,mytimestep) 
         B_temp = THRIFT_COEFF_B(j,mytimestep)
         C_temp = THRIFT_COEFF_C(j,mytimestep) 
         D_temp = THRIFT_COEFF_D(j,mytimestep)
         BP_temp = THRIFT_COEFF_BP(j,mytimestep)
         CP_temp = THRIFT_COEFF_CP(j,mytimestep)
         DP_temp = THRIFT_COEFF_DP(j,mytimestep)

         THRIFT_ALPHA1(i,mytimestep) = A_temp*DP_temp
         THRIFT_ALPHA2(i,mytimestep) = A_temp*(BP_temp/rho - B_temp/(rho**2) + CP_temp)  
         THRIFT_ALPHA3(i,mytimestep) = A_temp*(BP_temp + B_temp/rho + C_temp)            
         THRIFT_ALPHA4(i,mytimestep) = A_temp*B_temp   

         IF (lverbj) WRITE(6,'(F5.3, 1X, 4(ES13.5,2X))')  rho, &
            THRIFT_ALPHA1(i,mytimestep), THRIFT_ALPHA2(i,mytimestep),&
            THRIFT_ALPHA3(i,mytimestep), THRIFT_ALPHA4(i,mytimestep)
      END DO

!----------------------------------------------------------------------
!     SOLVING SYSTEM OF EQUATIONS
!----------------------------------------------------------------------
!     Populate matrix and RHS of the system of equations (N=nrho)
!
!     | B1   C1    0   0 .  0   0   0   0  | | u1  |   | D1  |  BC AXIS
!     | A2   B2   C2  X2 .  0   0   0   0  | | u2  |   | D2  |     \
!     |  0   A3   B3  C3 .  0   0   0   0  | | u3  |   | D3  |     |
!     |  .    .    .    ... .   .   .   .  | | .   | = |  .  |  evol. eq.
!     |  0    0    0   0 . AN  BN  CN   0  | | uN  |   | DN  |     |
!     |  0    0    0   0 . XN1 AN1 BN1 CN1 | | uN1 |   | DN1 |     /
!     |  0    0    0   0 .  0   0  AN2 BN2 | | UN2 |   | DN2 |  BC EDGE
!                    
!     The {uj} contain the current: u(j) = mu0*Itotal/(2*rho*phip)
!     First equation encapsulates the BC for the magnetic axis:
!        I(x=0,t) = 0 
!
!     Last equation encapsulates the BC for the plasma edge:
!        f(x=1,t) = <B^2>/mu0 d(rho*u)/drho + dp/drho*u - <Js.B> = 0
!
!     Remaining equations are the evolution equation on THRIFT_RHO grid: 
!        du/dt = a1 + a2*u + a3*du/drho + a4*d2u/drho2
!     LD = Ai   MD = Bi  UD = Ci  RHS = Di
!     NOTE: Since alphas exist on THRIFT_RHO grid but system of eqs.
!     calculates {uj} on [0,1], indices on RHS are shifted by -1.
!----------------------------------------------------------------------

      ! Magnetic axis (rho=0)
      THRIFT_MATLD(1,  mytimestep) = 0  
      THRIFT_MATMD(1,  mytimestep) = 1
      THRIFT_MATUD(1,  mytimestep) = 0
      THRIFT_MATRHS(1, mytimestep) = 0   

      drho = THRIFT_RHO(2)-THRIFT_RHO(1) 
      ! THRIFT_RHO grid (rho in (0,1))
      DO i = 2, nrho+1 
         j = i - 1 

         ! Temp variables for legibility
         a1 = THRIFT_ALPHA1(j, mytimestep) 
         a2 = THRIFT_ALPHA2(j, mytimestep) 
         a3 = THRIFT_ALPHA3(j, mytimestep) 
         a4 = THRIFT_ALPHA4(j, mytimestep)

         IF (j==1) THEN 
            THRIFT_MATLD(i, mytimestep) = -4*a3/(3*drho) + 16*a4/(5*drho**2)     ! ai, j = 1        ## i = 2
            THRIFT_MATMD(i, mytimestep) = a2 + a3/drho - 5*a4/(drho**2) - 1.0/dt ! bi, j = 1        ## i = 2
            THRIFT_MATUD(i, mytimestep) = a3/(3*drho) + 2*a4/(drho**2)           ! ci, j = 1        ## i = 2
         ELSEIF ((j>1).and.(j<nrho))
            THRIFT_MATLD(i, mytimestep) = -a3/drho + a4/(drho**2)                ! ai, 1 < j < nrho ## 2 < i < nrho+1
            THRIFT_MATMD(i, mytimestep) = a2 - 2*a4/(drho**2) - 1.0/dt           ! bi, 1 < j < nrho ## 2 < i < nrho+1
            THRIFT_MATUD(i, mytimestep) = a3/(2*drho) + a4/(drho**2)             ! ci, 1 < j < nrho ## 2 < i < nrho+1
         ELSEIF (j==nrho) THEN
            THRIFT_MATLD(i, mytimestep) = -a3/(3*drho) + 2*a4/(drho**2)          ! ai, j = nrho     ## i = nrho+1
            THRIFT_MATMD(i, mytimestep) = a2 - a3/drho - 5*a4/(drho**2) - 1.0/dt ! bi, j = nrho     ## i = nrho+1
            THRIFT_MATUD(i, mytimestep) = 4*a3/(3*drho) + 16*a4/(5*drho**2)      ! ci, j = nrho     ## i = nrho+1
         END IF
         THRIFT_MATRHS(i,mytimestep) = -THRIFT_UGRID(i, prevtimestep)/dt - a1 

      END DO 

      ! Plasma edge (rho=1)
      rho = 1
      CALL get_prof_pprime(rho, mytime, pprime)
      temp1 = THRIFT_BSQAV(nrho+2,mytimestep)/mu0
      THRIFT_MATLD(nrho+2,mytimestep)  = -1.0/drho
      THRIFT_MATMD(nrho+2,mytimestep)  = 1 + pprime/temp1 + (rho/(rho+0.5*drho))/drho
      THRIFT_MATUD(nrho+2,mytimestep)  = 0
      THRIFT_MATRHS(nrho+2,mytimestep) = jsource_full(nrho+2)*THRIFT_BAV(nrho+2,mytimestep)/temp1

      ! Nonzero element at M(2,4)
      temp1 = -THRIFT_ALPHA4(1,mytimestep)/(5*drho**2)
      temp1 = temp1/THRIFT_MATUD(3,mytimestep)
      THRIFT_MATUD( 2,mytimestep) = THRIFT_MATUD( 2,mytimestep) - temp1*THRIFT_MATMD( 3,mytimestep)
      THRIFT_MATMD( 2,mytimestep) = THRIFT_MATMD( 2,mytimestep) - temp1*THRIFT_MATLD( 3,mytimestep)
      THRIFT_MATRHS(2,mytimestep) = THRIFT_MATRHS(2,mytimestep) - temp1*THRIFT_MATRHS(3,mytimestep)

      ! Nonzero element at M(nrho+1,nrho-2)
      temp2 = -THRIFT_ALPHA4(nrho,mytimestep)/(5*drho**2)  
      temp2 = temp2/THRIFT_MATLD(nrho,mytimestep)
      THRIFT_MATLD( nrho+1,mytimestep) = THRIFT_MATLD( nrho+1,mytimestep) - temp2*THRIFT_MATMD( nrho,mytimestep)
      THRIFT_MATMD( nrho+1,mytimestep) = THRIFT_MATMD( nrho+1,mytimestep) - temp2*THRIFT_MATUD( nrho,mytimestep)
      THRIFT_MATRHS(nrho+1,mytimestep) = THRIFT_MATRHS(nrho+1,mytimestep) - temp2*THRIFT_MATRHS(nrho,mytimestep)

      ! Solve system of equations
      CALL solve_tdm( THRIFT_MATLD( :,mytimestep),&
                      THRIFT_MATMD( :,mytimestep),&
                      THRIFT_MATUD( :,mytimestep),&
                      THRIFT_MATRHS(:,mytimestep),&
                      THRIFT_UGRID( :,mytimestep))

      IF (lverbj) THEN
         WRITE(6,*) '==============================================================================='
         WRITE(6,*)'  i         LOWER           MAIN          UPPER            RHS       SOLUTION'
         WRITE(6,*)''
         DO i = 1, nrho+2
            WRITE(6,'(I4, 1X, 5(ES13.5,2X))') i, &
            THRIFT_MATLD(i,mytimestep), THRIFT_MATMD(i,mytimestep), &
            THRIFT_MATUD(i,mytimestep), THRIFT_MATRHS(i,mytimestep),THRIFT_UGRID(i,mytimestep)
         END DO
      END IF

!----------------------------------------------------------------------
!     POST SOLVING EQUATIONS
!----------------------------------------------------------------------
!     Solving system of equations yields {uj} at this timestep.
!     Get ITOTAL from u(j) = mu0*I/(2*rho*Phi_a)
!     Obtain ISOURCE from JSOURCE with curden_to_curtot subroutine.
!     Plasma current: IPLASMA = ITOTAL - ISOURCE
!     Obtain JPLASMA from IPLASMA with curtot_to_curden subroutine.
!----------------------------------------------------------------------

      THRIFT_I(:,mytimestep) = 2*THRIFT_PHIEDGE(1,mytimestep)/mu0*( THRIFT_RHOFULL*THRIFT_UGRID(:,mytimestep) )
      CALL curden_to_curtot(jsource_full,THRIFT_ISOURCE(:,mytimestep))
      THRIFT_IPLASMA(:,mytimestep) = THRIFT_I(:,mytimestep)-THRIFT_ISOURCE(:,mytimestep)
      CALL curtot_to_curden(THRIFT_IPLASMA(:,mytimestep),jplasma_full)
      THRIFT_JPLASMA(:,mytimestep) = jplasma_full(2:nrho+1)
      
      IF (lverbj) THEN
         WRITE(6,*) '==============================================================================='
         WRITE(6,*)' POST MATRIX ALGORITHM'
         WRITE(6,*)'  i        ITOTAL        ISOURCE        IPLASMA        JPLASMA        JSOURCE'
         WRITE(6,*)''
         DO i = 1, nrho+2
            WRITE(6,'(I4, 1X, 5(ES13.5,2X))') &
               i, THRIFT_I(i,mytimestep), THRIFT_ISOURCE(i,mytimestep), THRIFT_IPLASMA(i,mytimestep),&
               jplasma_full(i), jsource_full(i)
         END DO
         WRITE(6,*) '==============================================================================='
      END IF     

!----------------------------------------------------------------------
!     TRACKING TOTAL CURRENTS
!----------------------------------------------------------------------
!     The quantities THRIFT_IXXXXX are not necessary to evolve the 
!     current density profile, but are nice to have.
!----------------------------------------------------------------------
      1000  CONTINUE
      CALL extrapolate_arr(THRIFT_JBOOT(:,mytimestep),  j_full)
      CALL curden_to_curtot(j_full, THRIFT_IBOOT(:,mytimestep))
      CALL extrapolate_arr(THRIFT_JECCD(:,mytimestep),  j_full)
      CALL curden_to_curtot(j_full, THRIFT_IECCD(:,mytimestep))
      CALL extrapolate_arr(THRIFT_JNBCD(:,mytimestep),  j_full)
      CALL curden_to_curtot(j_full, THRIFT_INBCD(:,mytimestep))
      CALL extrapolate_arr(THRIFT_JOHMIC(:,mytimestep), j_full)
      CALL curden_to_curtot(j_full,THRIFT_IOHMIC(:,mytimestep))
      CALL curden_to_curtot(jplasma_full,THRIFT_IPLASMA(:,mytimestep))

      THRIFT_ISOURCE(:,mytimestep)  = THRIFT_IBOOT(:,mytimestep)&
                                    + THRIFT_IECCD(:,mytimestep)&
                                    + THRIFT_INBCD(:,mytimestep)&
                                    + THRIFT_IOHMIC(:,mytimestep)

      THRIFT_I(:,mytimestep)  = THRIFT_IPLASMA(:,mytimestep)&
                              + THRIFT_ISOURCE(:,mytimestep)

      DEALLOCATE( j_full, jplasma_full, jsource_full)

      RETURN

!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_jinductive