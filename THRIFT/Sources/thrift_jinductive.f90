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
                     a1,a2,a3,a4
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: A_temp,B_temp,C_temp,D_temp,&
                                                B_der, C_der, D_der, &
                                                alpha_1, alpha_2, alpha_3, alpha_4, &
                                                AI, BI, CI, DI, &
                                                rho_full, & 
                                                j_full, jsource_full, jplasma_full, jsourceprev_full

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
      ALLOCATE(A_temp(nrho+2), B_temp(nrho+2), C_temp(nrho+2), D_temp(nrho+2), &
               B_der(nrho+2),  C_der(nrho+2),  D_der(nrho+2), &
               alpha_1(nrho),  alpha_2(nrho),  alpha_3(nrho),  alpha_4(nrho), &
               AI(nrho+2),     BI(nrho+2),     CI(nrho+2),     DI(nrho+2), &
               rho_full(nrho+2),j_full(nrho+2),jplasma_full(nrho+2), &
               jsource_full(nrho+2),           jsourceprev_full(nrho+2))

      A_temp = 0; B_temp = 0; C_temp = 0; D_temp = 0;
      B_der  = 0; C_der  = 0; D_der  = 0; 
      alpha_1= 0; alpha_2= 0; alpha_3= 0; alpha_4= 0;
      AI     = 0; BI     = 0; CI     = 0; DI     = 0;
      
      rho_full=0;
      j_full = 0; 
      jplasma_full = 0;
      jsource_full = 0; 
      jsourceprev_full = 0;

      ! Define temp variable for the full grid
      rho_full(1) = 0.0
      rho_full(2:nrho+1) = THRIFT_RHO
      rho_full(nrho+2) = 1.0

      ! Extrapolate source current density to magnetic axis and edge
      CALL extrapolate_arr(THRIFT_JSOURCE(:,mytimestep), jsource_full)     

      !! ONLY NECESSARY AT CURRENT TIMESTEP !! REWRITE PLS 
      ! Store magnetic variables; both preceding and current timestep vals are necessary
      THRIFT_PHIEDGE(2) = eq_phiedge
      DO i = 1, nrho+2
         rho = rho_full(i)
         s = rho*rho
         ier = 0
         CALL EZspline_interp(vp_spl,rho,THRIFT_VP(i,2),ier)
         CALL get_equil_Bav(s,THRIFT_BAV(i,2), THRIFT_BSQAV(i,2), ier)
         CALL get_equil_sus(s,THRIFT_S11(i,2), temp1,temp1,temp1, ier)
         CALL get_equil_Rmajor(s,THRIFT_RMAJOR(i,2), temp1, THRIFT_AMINOR(i,2),ier)
      END DO
      THRIFT_S11 = ABS(THRIFT_S11)

      ! If mytimestep = 1 ITOT=0 and continue to next iteration
      IF (mytimestep==1) THEN
         jplasma_full = -jsource_full
         THRIFT_JPLASMA(:,mytimestep) = jplasma_full(2:nrho+1)
         GOTO 1000 ! skip iteration. 
      END IF

      mytime = THRIFT_T(mytimestep) ! mytime = current sim time
      prevtimestep = mytimestep-1 ! previous time step index
      dt = THRIFT_T(mytimestep)-THRIFT_T(prevtimestep) ! dt = delta t this iter
      IF (mytime>tmax.and.THRIFT_T(prevtimestep)<=tmax.and.(nsubsteps==1)) WRITE(6,*) &
         '! THRIFT has exceeded end time of profiles file. Proceeding with profiles at t=tmax !' 
         
!----------------------------------------------------------------------
!     CALCULATING COEFFICIENTS ABCD
!----------------------------------------------------------------------
!     Calculate ABCD (everything evaluated at rho_j)
!     > A(j) = S11/(4*rho*Phi_edge^2)
!     > B(j) = 2*etapara*dV/dPhi*<B^2>/mu_0
!     > C(j) = 2*etapara*dV/dPhi*dp/drho
!     > D(j) =-2*etapara*dV/dPhi*<Js.B>
!----------------------------------------------------------------------

      IF (lverbj) THEN
         WRITE(6,*)'==============================================================================='
         WRITE(6,*)' CALCULATING COEFFICIENTS A,B,C,D'
         WRITE(6,*) ' RHO  ETAPARA     DV/DPHI      DP/DRHO     <J.B>      BSQAV        S11'
      END IF
      DO i = 1, nrho+2
         rho = rho_full(i)
         temp2 = jsource_full(i)
         CALL get_prof_etapara(MIN(rho,THRIFT_RHO(nrho)),mytime,etapara)
         CALL get_prof_pprime(rho,mytime,pprime)
         temp1 = 2*etapara*THRIFT_VP(i,2) ! temp1 <- 2 eta dV/dPhi 
         IF (i > 1) &
            A_temp(i) = THRIFT_S11(i,2)/(4*rho*THRIFT_PHIEDGE(2)**2) ! S11/(4 rho phi_a^2)
         B_temp(i) = temp1*THRIFT_BSQAV(i,2)/mu0! 2 eta dV/dPhi <B^2>/mu_0
         C_temp(i) = temp1*pprime               ! 2 eta dV/dPhi dp/drho
         D_temp(i) = -temp1*temp2*THRIFT_BAV(i,2)   ! -2 eta dV/dPhi <J.B>
         IF (lverbj) WRITE(6,'(F5.3,6(1X,ES10.3))') &
           rho, etapara, THRIFT_VP(i,2), pprime, temp2*THRIFT_BAV(i,2), THRIFT_BSQAV(i,2), THRIFT_S11(i,2)
      END DO
      
      !! REWRITE: GET RID OF ABCD_TEMP AND USE THRIFT_COEFF INSTEAD
      THRIFT_COEFF_A(:,mytimestep) = A_temp
      THRIFT_COEFF_B(:,mytimestep) = B_temp
      THRIFT_COEFF_C(:,mytimestep) = C_temp
      THRIFT_COEFF_D(:,mytimestep) = D_temp
!----------------------------------------------------------------------
!     Calculate derivatives of ABCD here (see deriv1_rho_o2)
!     > drho = grid spacing
!----------------------------------------------------------------------
      drho = THRIFT_RHO(2)-THRIFT_RHO(1) 
      CALL deriv1_rho_o2(B_temp, drho, B_der)
      CALL deriv1_rho_o2(C_temp, drho, C_der)
      CALL deriv1_rho_o2(D_temp, drho, D_der)

      IF (lverbj) THEN
         WRITE(6,*)'==============================================================================='
         WRITE(6,*)' COEFFICIENTS ABCD'
         WRITE(6,*)' RHO         A         B          C          D       BDER       CDER       DDER'
         WRITE(6,*)''
         DO i = 1, nrho+2
            WRITE(6,'(F5.3, 1X, 7(ES10.2,1X))') rho_full(i), A_temp(i), B_temp(i), C_temp(i), D_temp(i),&
             B_der(i), C_der(i), D_der(i)
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

      DO i = 1, nrho
         rho = THRIFT_RHO(i)
         j = i + 1
         alpha_1(i) = A_temp(j)*D_der(j) 
         alpha_2(i) = A_temp(j)*(B_der(j)/rho - B_temp(j)/(rho**2) + C_der(j))  
         alpha_3(i) = A_temp(j)*(B_der(j) + B_temp(j)/rho + C_temp(j))            
         alpha_4(i) = A_temp(j)*B_temp(j)            
      END DO

      !! REWRITE: Same as coeff
      THRIFT_ALPHA1(:,mytimestep) = alpha_1; 
      THRIFT_ALPHA2(:,mytimestep) = alpha_2;
      THRIFT_ALPHA3(:,mytimestep) = alpha_3;
      THRIFT_ALPHA4(:,mytimestep) = alpha_4;

      IF (lverbj) THEN
         WRITE(6,*)'==============================================================================='
         WRITE(6,*)' ALPHAS'
         WRITE(6,*)'RHO       ALPHA 1        ALPHA 2        ALPHA 3        ALPHA 4'
         WRITE(6,*)''
         DO i = 1, nrho
            WRITE(6,'(F5.3, 1X, 4(ES13.5,2X))') THRIFT_RHO(i), alpha_1(i), alpha_2(i), alpha_3(i), alpha_4(i)
         END DO
      END IF
      
!----------------------------------------------------------------------
!     SOLVING SYSTEM OF EQUATIONS
!----------------------------------------------------------------------
!     Populate matrix and RHS of the system of equations (N=nrho)
!
!     | B1   C1    0   0 .  0   0   0   0  | | u1  |   | D1  |  BC AXIS
!     | A2   B2   C2   0 .  0   0   0   0  | | u2  |   | D2  |     \
!     |  0   A3   B3  C3 .  0   0   0   0  | | u3  |   | D3  |     |
!     |  .    .    .    ... .   .   .   .  | | .   | = |  .  |  evol. eq.
!     |  0    0    0   0 . AN  BN  CN   0  | | uN  |   | DN  |     |
!     |  0    0    0   0 .  0  AN1 BN1 CN1 | | uN1 |   | DN1 |     /
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
!
!     NOTE: Since alphas exist on THRIFT_RHO grid but system of eqs.
!     calculates {uj} on [0,1], indices on RHS are shifted by -1.
!----------------------------------------------------------------------

      ! Magnetic axis (rho=0)
      AI(1) = 0  
      BI(1) = 1
      CI(1) = 0
      DI(1) = 0   
      ! THRIFT_RHO grid (rho in (0,1))
      DO i = 2, nrho+1 
         j = i - 1 
         a1 = alpha_1(j); a2 = alpha_2(j); a3 = alpha_3(j); a4 = alpha_4(j)
         IF (i==2) THEN
            AI(i) = -a3/drho + 4*a4/(drho**2)
            BI(i) = a2 - a3/(2*drho) - 6*a4/(drho**2) - 1/dt
            CI(i) = a3/(2*drho) + 2*a4/(drho**2)
         ELSE IF (i==nrho+1) THEN
            AI(i) = -a3/(2*drho) + 2*a4/(drho**2)
            BI(i) = a2 - a3/(2*drho) - 6*a4/(drho**2) - 1/dt
            CI(i) = a3/drho + 4*a4/(drho**2)
         ELSE
            AI(i) = -a3/(2*drho) + a4/(drho**2)  
            BI(i) = a2 - 2*a4/(drho**2) - 1/dt     
            CI(i) = a3/(2*drho) + a4/(drho**2)  
         END IF
         DI(i) = - THRIFT_UGRID(j,1)/dt - a1 
      END DO
      ! Plasma edge
      rho = 1
      CALL get_prof_pprime(rho, mytime, pprime)
      temp1      = THRIFT_BSQAV(nrho+2,2)/mu0
      AI(nrho+2) = -THRIFT_RHO(nrho)/drho
      BI(nrho+2) = 1/drho + pprime/temp1
      CI(nrho+2) = 0
      DI(nrho+2) = jsource_full(nrho+2)*THRIFT_BAV(nrho+2,2)/temp1
      
      THRIFT_MATLD(:,mytimestep) = AI; 
      THRIFT_MATMD(:,mytimestep) = BI;
      THRIFT_MATUD(:,mytimestep) = CI;
      THRIFT_MATRHS(:,mytimestep)= DI;

      ! Solve system of equations
      CALL solve_tdm(AI,BI,CI,DI,THRIFT_UGRID(:,2))

      IF (lverbj) THEN
         WRITE(6,*) '==============================================================================='
         WRITE(6,*)'  i         LOWER           MAIN          UPPER            RHS       SOLUTION'
         WRITE(6,*)''
         DO i = 1, nrho+2
            WRITE(6,'(I4, 1X, 5(ES13.5,2X))') i, AI(i), BI(i), CI(i), DI(i),THRIFT_UGRID(i,2)
         END DO
      END IF

!----------------------------------------------------------------------
!     POST SOLVING EQUATIONS
!----------------------------------------------------------------------
!     Solution to previous system of equations yields {uj} at this step.
!     Get ITOTAL from u(j) = mu0*I/(2*rho*Phi_a)
!     Obtain ISOURCE from JSOURCE with curden_to_curtot subroutine.
!     Plasma current: IPLASMA = ITOTAL - ISOURCE
!     Obtain JPLASMA from IPLASMA with curtot_to_curden subroutine.
!----------------------------------------------------------------------

      THRIFT_I(:,mytimestep) = 2*THRIFT_PHIEDGE(2)/mu0*( rho_full*THRIFT_UGRID(:,2) )
      CALL curden_to_curtot(jsource_full,THRIFT_AMINOR,THRIFT_ISOURCE(:,mytimestep))
      THRIFT_IPLASMA(:,mytimestep) = THRIFT_I(:,mytimestep)-THRIFT_ISOURCE(:,mytimestep)
      CALL curtot_to_curden(THRIFT_IPLASMA(:,mytimestep),THRIFT_AMINOR,jplasma_full)
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
      CALL curden_to_curtot(j_full,THRIFT_AMINOR, THRIFT_IBOOT(:,mytimestep))
      CALL extrapolate_arr(THRIFT_JECCD(:,mytimestep),  j_full)
      CALL curden_to_curtot(j_full,THRIFT_AMINOR, THRIFT_IECCD(:,mytimestep))
      CALL extrapolate_arr(THRIFT_JNBCD(:,mytimestep),  j_full)
      CALL curden_to_curtot(j_full,THRIFT_AMINOR, THRIFT_INBCD(:,mytimestep))
      CALL extrapolate_arr(THRIFT_JOHMIC(:,mytimestep), j_full)
      CALL curden_to_curtot(j_full,THRIFT_AMINOR,THRIFT_IOHMIC(:,mytimestep))
      CALL curden_to_curtot(jplasma_full,THRIFT_AMINOR,THRIFT_IPLASMA(:,mytimestep))

      THRIFT_ISOURCE(:,mytimestep) = THRIFT_IBOOT(:,mytimestep)+THRIFT_IECCD(:,mytimestep)&
         +THRIFT_INBCD(:,mytimestep)+THRIFT_IOHMIC(:,mytimestep)
      THRIFT_I(:,mytimestep) = THRIFT_IPLASMA(:,mytimestep)+THRIFT_ISOURCE(:,mytimestep)

      DEALLOCATE( A_temp,  B_temp,  C_temp,  D_temp,  B_der, C_der, D_der, &
                  alpha_1, alpha_2, alpha_3, alpha_4, &
                  AI,      BI,      CI,      DI, &
                  rho_full, j_full, jplasma_full, jsource_full, jsourceprev_full)

      RETURN

!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_jinductive