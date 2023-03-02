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
      INTEGER :: i,itime, ier
      REAL(rprec) :: rho,s,drho,dt,t,s11,s12,etapara,pprime,&
                     temp1,temp2,source_axis,source_edge,Aminor,Rmajor
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: A_temp,B_temp,C_temp,D_temp,&
                                                B_der, C_der, D_der, &
                                                a1, a2, a3, a4, &
                                                AI, BI, CI, DI
!                                                rho_temp
!      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: datadiff
!
!
! STUPID DIFFUSION SUBROUTINE
!
! Allocations

                             
!DO i=1,nrho !init
!   rho_temp(i) = (i*1.0-1)/(nrho-1)
!   datadiff(i,1) = SIN(pi*rho_temp(i))
!END DO
!dt = THRIFT_T(5)-THRIFT_T(4)
!drho = rho_temp(2)-rho_temp(1)
!DO itime=2,20!
!
!   DO i = 1, nrho-1
!      AI(i) = 1/drho**2
!      BI(i) = -2/drho**2-1/dt
!      CI(i) = 1/drho**2
!      DI(i) = -datadiff(i,itime-1)/dt
!   END DO!
!
!   BI(1) = 1; CI(1) = 0; DI(1) = 0;
!   BI(nrho) = 1; AI(nrho-1) = 0; DI(nrho) = 0;!
!
!
 !  CALL solve_tdm(AI,BI,CI,DI,datadiff(:,itime))
!END DO
!WRITE(6,*) DATADIFF


!RETURN






!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------



      ! If at zero beta, copy previous value of JPLASMA onto this timestep
      IF (eq_beta == 0) THEN
         IF (mytimestep /= 1) THRIFT_JPLASMA(:,mytimestep) = THRIFT_JPLASMA(:,mytimestep-1)
         GOTO 1000 
      END IF

      ! Store magnetic variables; both preceding and current timestep vals are necessary
      THRIFT_PHIEDGE(2) = eq_phiedge
      DO i = 1, nrho+2
         IF (i==1) THEN
            rho = 0
         ELSE IF (i==nrho+2) THEN
            rho = 1
         ELSE
            rho = THRIFT_RHO(i-1)
         END IF
         s = rho*rho
         ier = 0
         CALL EZspline_interp(vp_spl,rho,THRIFT_VP(i,2),ier)
         CALL get_equil_Bav(s,THRIFT_BAV(i,2), THRIFT_BSQAV(i,2), ier)
         CALL get_equil_sus(s,THRIFT_S11(i,2), temp1,temp1,temp1, ier)
         CALL get_equil_Rmajor(s,THRIFT_RMAJOR(i,2), temp1, THRIFT_AMINOR(i,2),ier)
      END DO

      ! Calculate enclosed currents for source terms
      CALL curden_to_curtot(THRIFT_JBOOT,THRIFT_AMINOR,THRIFT_IBOOT,mytimestep)
      CALL curden_to_curtot(THRIFT_JECCD,THRIFT_AMINOR,THRIFT_IECCD,mytimestep)
      CALL curden_to_curtot(THRIFT_JNBCD,THRIFT_AMINOR,THRIFT_INBCD,mytimestep)
      CALL curden_to_curtot(THRIFT_JOHMIC,THRIFT_AMINOR,THRIFT_IOHMIC,mytimestep)
      THRIFT_ISOURCE(:,mytimestep) = THRIFT_IBOOT(:,mytimestep)+THRIFT_IECCD(:,mytimestep)&
         +THRIFT_INBCD(:,mytimestep)+THRIFT_IOHMIC(:,mytimestep)

      ! If mytimestep = 1 & tstart = 0, ITOT=0 and continue to next iteration
      ! If mytimestep = 1 & tstart > 0, ITOT=0 and calculate change in IPLASMA between tstart&t=0
      IF (mytimestep==1) THEN
         THRIFT_JPLASMA(:,mytimestep)=-THRIFT_JSOURCE(:,mytimestep) 
         CALL curden_to_curtot(THRIFT_JPLASMA,THRIFT_AMINOR,THRIFT_IPLASMA,mytimestep)
         IF (tstart==0) GOTO 1000
      END IF

      ! t = current simulation time
      ! dt = delta t between this and previous iteration (on first iter, dt=tstart)
      ! itime = previous timestep (on first iter, itime=mytimestep)
      t = THRIFT_T(mytimestep) ! t = current sim time
      IF (mytimestep==1) THEN
         itime = mytimestep
         dt = tstart
      ELSE
         itime = mytimestep-1
         dt = THRIFT_T(mytimestep)-THRIFT_T(itime)
      END IF
      IF (t>tmax.and.THRIFT_T(mytimestep-1)<=tmax.and.(nsubsteps==1)) WRITE(6,*) &

         '! THRIFT has exceeded end time of profiles file. Proceeding with profiles at t=tmax !' 
      ! Allocate
         ALLOCATE(A_temp(nrho+2), B_temp(nrho+2), C_temp(nrho+2), D_temp(nrho+2), &
         B_der(nrho),    C_der(nrho),    D_der(nrho), &
         a1(nrho  ), a2(nrho), a3(nrho  ), a4(nrho), &
         AI(nrho-1), BI(nrho), CI(nrho-1), DI(nrho))
         
         A_temp = 0; B_temp = 0; C_temp = 0; D_temp = 0;
         B_der  = 0; C_der  = 0; D_der  = 0;
         a1     = 0; a2     = 0; a3     = 0; a4     = 0;
         AI     = 0; BI     = 0; CI     = 0; DI     = 0;

      ! Calculate ABCD (values at **current** timestep)
      IF (lverbj) THEN
         WRITE(6,*)'==============================================================================='
         WRITE(6,*)' CALCULATING COEFFICIENTS A,B,C,D'
         WRITE(6,*) ' RHO  ETAPARA     DV/DPHI      DP/DRHO     <J.B>      BSQAV        S11'
      END IF

      ! Extrapolate source current density to magnetic axis and edge
      source_axis = (3*THRIFT_JSOURCE(1,mytimestep)-THRIFT_JSOURCE(2,mytimestep))/2
      source_edge = (-THRIFT_JSOURCE(nrho-1,mytimestep)+3*THRIFT_JSOURCE(nrho,mytimestep))/2

      DO i = 1, nrho+2
         ! calculate rho, etapara, current
         IF (i==1) THEN 
            rho = 0
            !CALL get_prof_etapara(rho,t,etapara)
            temp2 = source_axis
            !CALL get_prof_pprime(rho,t,pprime)
         ELSE IF (i==nrho+2) THEN ! etapara breaks when rho=1
            rho = 1
            !CALL get_prof_etapara(THRIFT_RHO(nrho),t,etapara)
            temp2 = source_edge
            !CALL get_prof_pprime(THRIFT_RHO(nrho),t,pprime)
         ELSE
            rho = THRIFT_RHO(i-1)
            !CALL get_prof_etapara(rho,t,etapara)
            temp2 = THRIFT_JSOURCE(i-1,mytimestep)
            !CALL get_prof_pprime(rho,t,pprime)
         END IF
         CALL get_prof_etapara(MIN(rho,THRIFT_RHO(nrho)),t,etapara)
         CALL get_prof_pprime(rho,t,pprime)
         temp1 = 2*etapara*THRIFT_VP(i,2) ! temp1 <- 2 eta dV/dPhi 
         IF (i /= 1) &
         A_temp(i) = THRIFT_S11(i,2)/(4*rho*THRIFT_PHIEDGE(2)**2) ! S11/(4 rho phi_a^2)
         B_temp(i) = temp1*THRIFT_BSQAV(i,2)/mu0! 2 eta dV/dPhi <B^2>/mu_0
         C_temp(i) = temp1*pprime               ! 2 eta dV/dPhi dp/drho
         D_temp(i) = -temp1*temp2*THRIFT_BAV(i,2)   ! -2 eta dV/dPhi <J.B>
         IF (lverbj) WRITE(6,'(F5.3,6(1X,ES10.3))') &
         rho, etapara, THRIFT_VP(i,2), pprime, temp2*THRIFT_BAV(i,2), THRIFT_BSQAV(i,2), THRIFT_S11(i,2)
      END DO
      THRIFT_COEFF_A(:,mytimestep) = A_temp
      THRIFT_COEFF_B(:,mytimestep) = B_temp
      THRIFT_COEFF_C(:,mytimestep) = C_temp
      THRIFT_COEFF_D(:,mytimestep) = D_temp
      ! drho = grid step (assuming constant spacing)
      drho = THRIFT_RHO(2)-THRIFT_RHO(1) 
      ! Calculate derivatives of ABCD
      CALL deriv1_rho_o2(B_temp, drho, B_der)
      CALL deriv1_rho_o2(C_temp, drho, C_der)
      CALL deriv1_rho_o2(D_temp, drho, D_der)

      IF (lverbj) THEN
         WRITE(6,*)'==============================================================================='
         WRITE(6,*)' COEFFICIENTS ABCD'
         WRITE(6,*)' RHO         A         B          C          D       BDER       CDER       DDER'
         WRITE(6,*)''
         WRITE(6,'(F5.3, 1X, 7(ES10.2,1X))') 0.0, A_temp(1), B_temp(1), C_temp(1), D_temp(1), &
             0.0, 0.0, 0.0
         DO i = 2, nrho+1
            WRITE(6,'(F5.3, 1X, 7(ES10.2,1X))') THRIFT_RHO(i-1), A_temp(i), B_temp(i), C_temp(i), D_temp(i),&
             B_der(i-1), C_der(i-1), D_der(i-1)
         END DO
         WRITE(6,'(F5.3, 1X, 7(ES10.2,1X))') 1.0, A_temp(nrho+2), B_temp(nrho+2), C_temp(nrho+2), D_temp(nrho+2),&
             0.0, 0.0, 0.0
      END IF

      DO i = 1, nrho
         rho = THRIFT_RHO(i)
         a1(i) = A_temp(i+1)*D_der(i) ! a1 = A dD/drho
         a2(i) = A_temp(i+1)*(B_der(i)/rho - B_temp(i+1)/(rho**2) + C_der(i))  ! a2 = A (1/rho dB/drho - B/rho^2 + dC/drho)  
         a3(i) = A_temp(i+1)*(B_der(i) + B_temp(i+1)/rho + C_temp(i+1))      ! a3 = A (dB/drho + B/rho + C)             
         a4(i) = A_temp(i+1)*B_temp(i+1)         ! a4 = A B                    
      END DO

      a2 = 0; a3 = 0;
      THRIFT_ALPHA1(:,mytimestep) = a1; 
      THRIFT_ALPHA2(:,mytimestep) = a2;
      THRIFT_ALPHA3(:,mytimestep) = a3;
      THRIFT_ALPHA4(:,mytimestep) = a4;

      IF (lverbj) THEN
         WRITE(6,*)'==============================================================================='
         WRITE(6,*)' ALPHAS'
         WRITE(6,*)'RHO       ALPHA 1        ALPHA 2        ALPHA 3        ALPHA 4'
         WRITE(6,*)''
         DO i = 1, nrho
            WRITE(6,'(F5.3, 1X, 4(ES13.5,2X))') THRIFT_RHO(i), a1(i), a2(i), a3(i), a4(i)
         END DO
      END IF
      
      ! Calculate L_ext = mu0*R_eff( ln( 8 R_eff/r_eff )-2 + F_shaping)      
      temp1 = mu0*THRIFT_RMAJOR(nrho+2,1)*&
         (log(8*THRIFT_RMAJOR(nrho+2,1)/THRIFT_AMINOR(nrho+2,1)) - 2 + 0.25) ! temp1 <- L_ext
      ! L/R time
      t = THRIFT_T(itime)
      CALL get_prof_etapara(THRIFT_RHO(nrho),t,etapara) 
      temp2 = temp1*THRIFT_AMINOR(nrho+2,1)**2/(2*etapara*THRIFT_RMAJOR(nrho+2,1)) ! temp2 <- tau_L/R
      IF (nsubsteps==1.and.(mytimestep==1.or.(mytimestep==2.and.tstart==0))) &
         WRITE(6,'(A25,F8.6)') 'Estimated tau_L/R    ',temp2
      ! Decay I_plasma at edge
      THRIFT_IPLASMA(nrho,mytimestep) = THRIFT_IPLASMA(nrho,itime)*exp(-dt/temp2)
      !WRITE(6,*) THRIFT_IPLASMA(nrho,itime)
      !WRITE(6,*) dt
      !WRITE(6,*) temp2
      !WRITE(6,*) EXP(-dt/temp2)
      !WRITE(6,*) THRIFT_IPLASMA(nrho,mytimestep)
      ! I_total at edge
      temp1 = THRIFT_IPLASMA(nrho,mytimestep)+THRIFT_ISOURCE(nrho,mytimestep)
      WRITE(6,*) 'I EDGE PREV STEP'
      WRITE(6,*) THRIFT_I(nrho,itime)
      WRITE(6,*) 'I EDGE THIS STEP (PREDICTED)'
      WRITE(6,*) temp1

       !! Calculate uedge for this timestep
      !t = THRIFT_T(itime) ! t = previous sim time (or current sim time if mytimestep=1)
      !CALL get_prof_etapara(THRIFT_RHO(INT(nrho/2)),t,etapara)
      !CALL get_prof_etapara(THRIFT_RHO(nrho),t,etapara)
      !CALL get_prof_pprime(THRIFT_RHO(nrho),t,pprime) 
      !temp2 = (-8*THRIFT_UEDGE(1)+9*THRIFT_UGRID(nrho,1)-THRIFT_UGRID(nrho-1,1))/(3*drho) ! du/drho at edge
      ! uedge [ current timestep ] =  u - dt*mu0/(2*Phi_a)*eta/Lext*dV/dphi* ...
      ! [(p' + <B^2>/mu0)*u + <B^2>/mu0*du/drho - <J.B>]  [ edge, preceding timestep ]
      !THRIFT_UEDGE(2) = THRIFT_UEDGE(1)-dt*mu0/(2*THRIFT_PHIEDGE(1))*etapara/temp1*THRIFT_VP(nrho+2,1) * &  
      !      (((pprime + THRIFT_BSQAV(nrho+2,1)/mu0)*THRIFT_UEDGE(1)) + THRIFT_BSQAV(nrho+2,1)/mu0*temp2 &               
      !      - THRIFT_JSOURCE(nrho,itime)*THRIFT_BAV(nrho+2,1))   



      ! Populate tridiagonal matrix and RHS
      ! BC1: Enclosed current at magnetic axis = 0 always
      BI(1) = 1; CI(1) = 0; DI(1) = 0   
      DO i = 2, nrho-1 ! Between boundaries
         IF (i/=nrho-1) & 
         AI(i) = -a3(i)/(2*drho)+a4(i)/(drho**2)  
         BI(i) = a2(i)-2*a4(i)/(drho**2)-1/dt     
         CI(i) = a3(i)/(2*drho) +a4(i)/(drho**2)   
         DI(i) = -THRIFT_UGRID(i,1)/dt-a1(i)  
      END DO
      ! BC2: Enclosed current at edge must equal temp1 next timestep
      BI(nrho) = 1; AI(nrho-1) = 0; DI (nrho) = mu0*temp1/(2*THRIFT_RHO(nrho)*THRIFT_PHIEDGE(2)) 

      THRIFT_MATLD(:,mytimestep) = AI; 
      THRIFT_MATMD(:,mytimestep) = BI;
      THRIFT_MATUD(:,mytimestep) = CI;
      THRIFT_MATRHS(:,mytimestep)= DI;

      ! Solve system of equations
      CALL solve_tdm(AI,BI,CI,DI,THRIFT_UGRID(:,2))
      !CALL check_sol(AI,BI,CI,DI,THRIFT_UGRID(:,2),B_der)
      !WRITE(6,*) 'RESIDUES'
      !WRITE(6,*) B_der
      !WRITE(6,*) 'END RESIDUES'

      IF (lverbj) THEN
         WRITE(6,*) '==============================================================================='
         WRITE(6,*)'  i         LOWER           MAIN          UPPER            RHS       SOLUTION'
         WRITE(6,*)''
         WRITE(6,'(I4, 1X,5(2X,ES13.5))') 1, 0.0, BI(1), CI(1), DI(1),THRIFT_UGRID(1,2)
         DO i = 2, nrho-1
            WRITE(6,'(I4, 1X, 5(ES13.5,2X))') i, AI(i-1), BI(i), CI(i), DI(i),THRIFT_UGRID(i,2)
         END DO
         WRITE(6,'(I4, 1X, 5(ES13.5,2X))') nrho, AI(nrho-1),BI(nrho), 0.0,DI(nrho), THRIFT_UGRID(nrho,2)
      END IF

      ! ITOTAL = phip*u/mu0 = 2*phi_a*rho*u/mu0
      THRIFT_I(:,mytimestep) = 2*THRIFT_PHIEDGE(2)/mu0*(THRIFT_RHO*THRIFT_UGRID(:,2))
      WRITE(6,*) 'I EDGE THIS STEP (CALCULATED)'
      WRITE(6,*) THRIFT_I(nrho,mytimestep)
       
      CALL curden_to_curtot(THRIFT_JSOURCE,THRIFT_AMINOR,THRIFT_ISOURCE,mytimestep)
      ! IPLASMA = ITOTAL - ISOURCE
      THRIFT_IPLASMA(:,mytimestep) = THRIFT_I(:,mytimestep) - THRIFT_ISOURCE(:,mytimestep)
      ! JPLASMA
      CALL curtot_to_curden(THRIFT_IPLASMA,THRIFT_AMINOR,THRIFT_JPLASMA,mytimestep)
      ! Subtract change in JSOURCE
      THRIFT_JPLASMA(:,mytimestep)=THRIFT_JPLASMA(:,mytimestep)-(THRIFT_JSOURCE(:,mytimestep)-THRIFT_JSOURCE(:,itime))

      IF (lverbj) THEN
         WRITE(6,*) '==============================================================================='
         WRITE(6,*)' POST MATRIX ALGORITHM'
         WRITE(6,*)'  i        ITOTAL        ISOURCE        IPLASMA        JPLASMA        JSOURCE'
         WRITE(6,*)''
         DO i = 1, nrho
            WRITE(6,'(I4, 1X, 5(ES13.5,2X))') &
               i, THRIFT_I(i,mytimestep), THRIFT_ISOURCE(i,mytimestep), THRIFT_IPLASMA(i,mytimestep),&
               THRIFT_JPLASMA(i,mytimestep), THRIFT_JSOURCE(i,mytimestep)
         END DO
         WRITE(6,*) '==============================================================================='
      END IF     
      DEALLOCATE(A_temp, B_temp, C_temp, D_temp, &
               B_der, C_der, D_der, &
               a1, a2, a3, a4, &
               AI, BI, CI, DI)

1000  CONTINUE
      ! Calculate enclosed plasma current
      CALL curden_to_curtot(THRIFT_JPLASMA,THRIFT_AMINOR,THRIFT_IPLASMA,mytimestep)
      THRIFT_I(:,mytimestep) = THRIFT_IPLASMA(:,mytimestep)+THRIFT_ISOURCE(:,mytimestep)
      RETURN

!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_jinductive