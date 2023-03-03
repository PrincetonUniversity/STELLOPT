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
      REAL(rprec) :: rho,s,drho,dt,t,s11,s12,etapara,pprime,&
                     temp1,temp2,source_axis,source_edge,Aminor,Rmajor
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: A_temp,B_temp,C_temp,D_temp,&
                                                B_der, C_der, D_der, &
                                                a1, a2, a3, a4, &
                                                AI, BI, CI, DI, &
                                                rho_full, & 
                                                j_full, jsource_full, jplasma_full, jsourceprev_full
!      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: datadiff
!
! STUPID DIFFUSION SUBROUTINE
!                          
!DO i=1,nrho !init
!   rho_temp(i) = (i*1.0-1)/(nrho-1)
!   datadiff(i,1) = SIN(pi*rho_temp(i))
!END DO
!dt = THRIFT_T(5)-THRIFT_T(4)
!drho = rho_temp(2)-rho_temp(1)
!DO itime=2,20!
!   DO i = 1, nrho-1
!      AI(i) = 1/drho**2
!      BI(i) = -2/drho**2-1/dt
!      CI(i) = 1/drho**2
!      DI(i) = -datadiff(i,itime-1)/dt
!   END DO!
!   BI(1) = 1; CI(1) = 0; DI(1) = 0;
!   BI(nrho) = 1; AI(nrho-1) = 0; DI(nrho) = 0;!
 !  CALL solve_tdm(AI,BI,CI,DI,datadiff(:,itime))
!END DO
!WRITE(6,*) DATADIFF
!RETURN

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

      ! Allocation (part 1)
      ALLOCATE(rho_full(nrho+2),j_full(nrho+2), jplasma_full(nrho+2), jsource_full(nrho+2), jsourceprev_full(nrho+2))

      rho_full(1) = 0.0
      rho_full(2:nrho+1) = THRIFT_RHO
      rho_full(nrho+2) = 1.0
      
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

      ! Calculate enclosed currents for source terms (tracking only)
      CALL extrapolate_arr( THRIFT_JBOOT(:,mytimestep), j_full)
      CALL curden_to_curtot(j_full,THRIFT_AMINOR, THRIFT_IBOOT(:,mytimestep))
      CALL extrapolate_arr( THRIFT_JECCD(:,mytimestep), j_full)
      CALL curden_to_curtot(j_full,THRIFT_AMINOR, THRIFT_IECCD(:,mytimestep))
      CALL extrapolate_arr( THRIFT_JNBCD(:,mytimestep), j_full)
      CALL curden_to_curtot(j_full,THRIFT_AMINOR, THRIFT_INBCD(:,mytimestep))
      CALL extrapolate_arr(THRIFT_JOHMIC(:,mytimestep), j_full)
      CALL curden_to_curtot(j_full,THRIFT_AMINOR,THRIFT_IOHMIC(:,mytimestep))
      THRIFT_ISOURCE(:,mytimestep) = THRIFT_IBOOT(:,mytimestep)+THRIFT_IECCD(:,mytimestep)&
         +THRIFT_INBCD(:,mytimestep)+THRIFT_IOHMIC(:,mytimestep)

      ! Extrapolate source current density to magnetic axis and edge
      CALL extrapolate_arr(THRIFT_JSOURCE(:,mytimestep), jsource_full)
           
      ! If mytimestep = 1 & tstart = 0, ITOT=0 and continue to next iteration
      ! If mytimestep = 1 & tstart > 0, ITOT=0 and calculate change in IPLASMA between tstart&t=0
      IF (mytimestep==1.and.tstart==0) THEN
         jplasma_full = -jsource_full
         THRIFT_JPLASMA(:,mytimestep) = jplasma_full(2:nrho+1)
         GOTO 1000 ! nothing 
      END IF

      ! t = current simulation time
      ! dt = delta t between this and previous iteration (on first iter, dt=tstart)
      ! prevtimestep = previous timestep (on first iter, prevtimestep=mytimestep)
      t = THRIFT_T(mytimestep) ! t = current sim time
      IF (mytimestep==1) THEN
         prevtimestep = mytimestep
         dt = tstart
      ELSE
         prevtimestep = mytimestep-1
         dt = THRIFT_T(mytimestep)-THRIFT_T(prevtimestep)
      END IF
      IF (t>tmax.and.THRIFT_T(mytimestep-1)<=tmax.and.(nsubsteps==1)) WRITE(6,*) &
         '! THRIFT has exceeded end time of profiles file. Proceeding with profiles at t=tmax !' 

      ! Allocate
      ALLOCATE(A_temp(nrho+2), B_temp(nrho+2), C_temp(nrho+2), D_temp(nrho+2), &
        B_der(nrho+2),    C_der(nrho+2),    D_der(nrho+2), &
        a1(nrho  ), a2(nrho), a3(nrho  ), a4(nrho), &
        AI(nrho+2), BI(nrho+2), CI(nrho+2), DI(nrho+2))
      
      A_temp = 0; B_temp = 0; C_temp = 0; D_temp = 0;
      B_der  = 0; C_der  = 0; D_der  = 0;
      a1     = 0; a2     = 0; a3     = 0; a4     = 0;
      AI     = 0; BI     = 0; CI     = 0; DI     = 0;

!----------------------------------------------------------------------
!     CALCULATING COEFFICIENTS
!----------------------------------------------------------------------

      ! Calculate ABCD (values at **current** timestep)
      IF (lverbj) THEN
         WRITE(6,*)'==============================================================================='
         WRITE(6,*)' CALCULATING COEFFICIENTS A,B,C,D'
         WRITE(6,*) ' RHO  ETAPARA     DV/DPHI      DP/DRHO     <J.B>      BSQAV        S11'
      END IF

      DO i = 1, nrho+2
         rho = rho_full(i)
         temp2 = jsource_full(i)
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
         DO i = 1, nrho+2
            WRITE(6,'(F5.3, 1X, 7(ES10.2,1X))') rho_full(i), A_temp(i), B_temp(i), C_temp(i), D_temp(i),&
             B_der(i), C_der(i), D_der(i)
         END DO
      END IF

      DO i = 1, nrho
         rho = THRIFT_RHO(i)
         j = i + 1
         a1(i) = A_temp(j)*D_der(j) ! a1 = A dD/drho
         a2(i) = A_temp(j)*(B_der(j)/rho - B_temp(j)/(rho**2) + C_der(j))  ! a2 = A (1/rho dB/drho - B/rho^2 + dC/drho)  
         a3(i) = A_temp(j)*(B_der(j) + B_temp(j)/rho + C_temp(j))      ! a3 = A (dB/drho + B/rho + C)             
         a4(i) = A_temp(j)*B_temp(j)         ! a4 = A B                    
      END DO

      !a1 = 0; a2 = 0; a3 = 0;
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
      
      !! OLD EXPONENTIAL DECAY CODE
      !! Calculate L_ext = mu0*R_eff( ln( 8 R_eff/r_eff )-2 + F_shaping)      
      !temp1 = mu0*THRIFT_RMAJOR(nrho+2,1)*&
      !   (log(8*THRIFT_RMAJOR(nrho+2,1)/THRIFT_AMINOR(nrho+2,1)) - 2 + 0.25) ! temp1 <- L_ext
      ! L/R time
      !t = THRIFT_T(itime)
      !CALL get_prof_etapara(THRIFT_RHO(nrho),t,etapara) 
      !temp2 = temp1*THRIFT_AMINOR(nrho+2,1)**2/(2*etapara*THRIFT_RMAJOR(nrho+2,1)) ! temp2 <- tau_L/R
      !IF (nsubsteps==1.and.(mytimestep==1.or.(mytimestep==2.and.tstart==0))) &
      !   WRITE(6,'(A25,F8.6)') 'Estimated tau_L/R    ',temp2
      ! Decay I_plasma at edge
      !THRIFT_IPLASMA(nrho,mytimestep) = THRIFT_IPLASMA(nrho,itime)*exp(-dt/temp2)
      ! I_total at edge
      !temp1 = THRIFT_IPLASMA(nrho,mytimestep)+THRIFT_ISOURCE(nrho,mytimestep)

      !! OLD UEDGE CODE
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

      !! OLD UEDGE CODE
      ! No gradient at edge
      !rho = 1
      !CALL get_prof_pprime(rho, t, pprime)
      !temp1 = THRIFT_BSQAV(nrho+2,2)/mu0+pprime
      !temp2 = 1+(1/(THRIFT_RHO(nrho)+drho))*(1-pprime/temp1)
      !temp1 = (source_edge*THRIFT_BAV(nrho+2,2)/temp1)/temp2

      ! Populate tridiagonal matrix and RHS
      ! BC1: Enclosed current at magnetic axis = 0 always
      AI(1) = 0; BI(1) = 1; CI(1) = 0; DI(1) = 0   
      DO i = 2, nrho+1 ! On the proper grid
         j = i - 1 
         IF (i==2) THEN
            AI(i) = -a3(j)/drho+4*a4(j)/(drho**2)
            BI(i) = a2(j)-a3(j)/(2*drho)-6*a4(1)/(drho**2)-1/dt
            CI(i) = a3(j)/(2*drho)+2*a4(j)/(drho**2)
         ELSE IF (i==nrho+1) THEN
            AI(i) = -a3(j)/(2*drho)+2*a4(j)/(drho**2)
            BI(i) = a2(j)-a3(j)/(2*drho)-6*a4(j)/(drho**2)-1/dt
            CI(i) = a3(j)/drho+4*a4(j)/(drho**2)
         ELSE
            AI(i) = -a3(j)/(2*drho)+a4(j)/(drho**2)  
            BI(i) = a2(j)-2*a4(j)/(drho**2)-1/dt     
            CI(i) = a3(j)/(2*drho) +a4(j)/(drho**2)  
         END IF
         DI(i) = -THRIFT_UGRID(j,1)/dt-a1(j)  
      END DO
      !! New BC2
      rho = 1
      CALL get_prof_pprime(rho,t, pprime)
      AI(nrho+2) = -THRIFT_RHO(nrho)/drho
      CI(nrho+2) = 0
      temp1 = THRIFT_BSQAV(nrho+2,2)/mu0
      BI(nrho+2) = 1/drho+mu0*pprime/temp1
      DI(nrho+2) = jsource_full(nrho+2)*THRIFT_BAV(nrho+2,2)/temp1
      
      !AI(nrho-1) = 0; BI(nrho) = 1; DI(nrho) = THRIFT_RHO(nrho)/(THRIFT_RHO(nrho)+drho)*temp1


      !! Old BC2: Enclosed current at edge must equal temp1 next timestep
      !BI(nrho) = 1; AI(nrho-1) = 0; DI (nrho) = mu0*temp1/(2*THRIFT_RHO(nrho)*THRIFT_PHIEDGE(2)) 
      !! Old BC2: No gradient at edge
      !AI(nrho-2)=0; BI(nrho-1)=1; CI(nrho-1)=-1; DI(nrho)=0
      !AI(nrho-1) = 0
      !BI(nrho) = a2(nrho)-2*a4(nrho)/(drho**2)-1/dt
      !DI(nrho) = -THRIFT_UGRID(nrho,1)/dt-a1(nrho)

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
         DO i = 1, nrho+2
            WRITE(6,'(I4, 1X, 5(ES13.5,2X))') i, AI(i), BI(i), CI(i), DI(i),THRIFT_UGRID(i,2)
         END DO
      END IF

      ! ITOTAL = phip*u/mu0 = 2*phi_a*rho*u/mu0
      THRIFT_I(:,mytimestep) = 2*THRIFT_PHIEDGE(2)/mu0*(rho_full*THRIFT_UGRID(:,2))
      ! JTOTAL
      CALL curtot_to_curden(THRIFT_I(:,mytimestep),THRIFT_AMINOR,j_full)
      THRIFT_J(:,mytimestep) = j_full(2:nrho+1)
      ! JPLASMA = JTOTAL - JSOURCE
      jplasma_full = j_full - jsource_full
      ! Subtract change in JSOURCE
      CALL extrapolate_arr(THRIFT_JSOURCE(:,prevtimestep), jsourceprev_full)
      !jplasma_full = jplasma_full - (jsource_full-jsourceprev_full)
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

      DEALLOCATE(A_temp, B_temp, C_temp, D_temp, &
               B_der, C_der, D_der, &
               a1, a2, a3, a4, &
               AI, BI, CI, DI)

      1000  CONTINUE
      ! Calculate enclosed plasma current
      CALL curden_to_curtot(jplasma_full,THRIFT_AMINOR,THRIFT_IPLASMA(:,mytimestep))
      DEALLOCATE(rho_full, j_full, jplasma_full, jsource_full, jsourceprev_full)
      THRIFT_I(:,mytimestep) = THRIFT_IPLASMA(:,mytimestep)+THRIFT_ISOURCE(:,mytimestep)
      RETURN

!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_jinductive