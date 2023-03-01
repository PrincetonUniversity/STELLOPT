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
      REAL(rprec) :: rho,s,h,k,t,s11,s12,etapara,pprime,&
                     temp1,temp2,source_axis,source_edge,Aminor,Rmajor
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: A_temp, B_temp, C_temp, D_temp, &
                                                B_der, C_der, D_der, &
                                                a1, a2, a3, a4, &
                                                AI, BI, CI, DI
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      ! If at zero beta, copy previous value of JPLASMA
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

      ! Calculate currents from source equations
      CALL curden_to_curtot(THRIFT_JBOOT,THRIFT_AMINOR,THRIFT_IBOOT,mytimestep)
      CALL curden_to_curtot(THRIFT_JECCD,THRIFT_AMINOR,THRIFT_IECCD,mytimestep)
      CALL curden_to_curtot(THRIFT_JNBCD,THRIFT_AMINOR,THRIFT_INBCD,mytimestep)
      CALL curden_to_curtot(THRIFT_JOHMIC,THRIFT_AMINOR,THRIFT_IOHMIC,mytimestep)
      THRIFT_ISOURCE(:,mytimestep) = THRIFT_IBOOT(:,mytimestep)+THRIFT_IECCD(:,mytimestep)&
         +THRIFT_INBCD(:,mytimestep)+THRIFT_IOHMIC(:,mytimestep)

      ! If we start at t=0, declare total current density = 0 and skip timestep
      IF (tstart==0 .and. mytimestep==1) THEN
         THRIFT_JPLASMA(:,mytimestep)=-THRIFT_JSOURCE(:,mytimestep) 
         GOTO 1000
      END IF

      ! Allocations
      ALLOCATE(A_temp(nrho+2), B_temp(nrho+2), C_temp(nrho+2), D_temp(nrho+2), &
               B_der(nrho), C_der(nrho), D_der(nrho), &
               a1(nrho), a2(nrho), a3(nrho), a4(nrho), &
               AI(nrho-1), BI(nrho), CI(nrho-1), DI(nrho))

      A_temp = 0; B_temp = 0; C_temp = 0; D_temp = 0;
      B_der = 0; C_der = 0; D_der = 0;
      a1 = 0; a2 = 0; a3 = 0; a4 = 0;
      AI = 0; BI = 0; CI = 0; DI = 0;

      ! k = Delta t
      IF (mytimestep==1) THEN
         itime = mytimestep
         k = tstart
      ELSE
         itime = mytimestep-1
         k = THRIFT_T(mytimestep)-THRIFT_T(itime)
      END IF
      t = THRIFT_T(mytimestep) ! t = current sim time
      IF (t>tmax.and.THRIFT_T(mytimestep-1)<=tmax) WRITE(6,*) &
         '! THRIFT has exceeded end time of profiles file. Proceeding with profiles at t=tmax !' 

      IF (lverbj) THEN
         WRITE(6,*)'==============================================================================='
         WRITE(6,*)' CALCULATING RESISTIVITY'
         WRITE(6,*) ' RHO         TE         NE       ZEFF     COULLN    ETAPERP    ETAPARA'
      END IF     
      DO i = 1,nrho+2
         IF (i==1) THEN
            rho = 0
         ELSE IF (i==nrho+2) THEN
            rho = 1
         ELSE
            rho = THRIFT_RHO(i-1)
         ENDIF
         CALL get_prof_te(rho, t, temp1)
         CALL get_prof_ne(rho, t, temp2)
         CALL get_prof_coulln(rho, t, s11)
         CALL get_prof_zeff(rho, t, s12)
         CALL get_prof_etapara(rho,t, etapara)
         CALL get_prof_etaperp(rho,t, pprime)
         IF (lverbj) WRITE(6,'(F5.3,6(1X,ES10.3))') rho, temp1, temp2, s12, s11, pprime, etapara
      END DO


      ! Populate A,B,C,D
      IF (lverbj) THEN
         WRITE(6,*)'==============================================================================='
         WRITE(6,*)' CALCULATING COEFFICIENTS A,B,C,D'
         WRITE(6,*) ' RHO  ETAPARA     DV/DPHI      DP/DRHO     <J.B>      BSQAV        S11'
      END IF
      
      ! A,B,C,D should be evaluated at the current timestep
      ! Extrapolate source current density to magnetic axis and edge
      source_axis = (3*THRIFT_JSOURCE(1,itime)-THRIFT_JSOURCE(2,itime))/2
      !source_edge = (-THRIFT_JSOURCE(nrho-1,itime)+3*THRIFT_JSOURCE(nrho,itime))/2
      source_edge = THRIFT_JSOURCE(nrho,itime)

      DO i = 1, nrho+2
         ! calculate rho, etapara, current
         IF (i==1) THEN 
            rho = 0
            CALL get_prof_etapara(rho,t,etapara)
            temp2 = source_axis
            CALL get_prof_pprime(rho,t,pprime)
         ELSE IF (i==nrho+2) THEN ! etapara breaks when rho=1
            rho = 1
            CALL get_prof_etapara(THRIFT_RHO(nrho),t,etapara)
            temp2 = source_edge
            CALL get_prof_pprime(THRIFT_RHO(nrho),t,pprime)
         ELSE
            rho = THRIFT_RHO(i-1)
            CALL get_prof_etapara(rho,t,etapara)
            temp2 = THRIFT_JSOURCE(i-1,itime)
            CALL get_prof_pprime(rho,t,pprime)
         END IF
         !CALL get_prof_etapara(THRIFT_RHO(INT(nrho/2)),t,etapara) ! clamp etapara
         !CALL get_prof_pprime(rho,t,pprime)
         temp1 = 2*etapara*THRIFT_VP(i,2) ! temp1 <- 2 eta dV/dPhi 
         IF (i /= 1) A_temp(i) = THRIFT_S11(i,2)/(4*rho*THRIFT_PHIEDGE(2)**2) ! S11/(4 rho phi_a^2)
         B_temp(i) = temp1*THRIFT_BSQAV(i,2)/mu0! 2 eta dV/dPhi <B^2>/mu_0
         C_temp(i) = temp1*pprime               ! 2 eta dV/dPhi dp/drho
         D_temp(i) = -temp1*temp2*THRIFT_BAV(i,2)   ! -2 eta dV/dPhi <J.B>
         IF (lverbj) WRITE(6,'(F5.3,6(1X,ES10.3))') &
         rho, etapara, THRIFT_VP(i,2), pprime, temp2*THRIFT_BAV(i,2), THRIFT_BSQAV(i,2), THRIFT_S11(i,2)
      END DO

      ! Calculate derivatives of ABCD
      h = THRIFT_RHO(2)-THRIFT_RHO(1) ! h = Delta rho
      CALL deriv1_rho_o2(B_temp, h, B_der)
      CALL deriv1_rho_o2(C_temp, h, C_der)
      CALL deriv1_rho_o2(D_temp, h, D_der)

      !DO i = 2, nrho-1
      !   B_der(i) = (B_temp(i+2)-B_temp(i))/(2*h)
      !   C_der(i) = (C_temp(i+2)-C_temp(i))/(2*h)
      !   D_der(i) = (D_temp(i+2)-D_temp(i))/(2*h)
      !END DO
      
      !B_der(1) = (B_temp(3)+3*B_temp(2)-4*B_temp(1))/(3*h)
      !C_der(1) = (C_temp(3)+3*C_temp(2)-4*C_temp(1))/(3*h)
      !D_der(1) = (D_temp(3)+3*D_temp(2)-4*D_temp(1))/(3*h)

      !B_der(nrho) = (4*B_temp(nrho+2)-3*B_temp(nrho+1)-B_temp(nrho))/(3*h)
      !C_der(nrho) = (4*C_temp(nrho+2)-3*C_temp(nrho+1)-C_temp(nrho))/(3*h)
      !D_der(nrho) = (4*D_temp(nrho+2)-3*D_temp(nrho+1)-D_temp(nrho))/(3*h)
      
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

      ! a1 = A dD/drho
      ! a2 = A (1/rho dB/drho - B/rho^2 + dC/drho)
      ! a3 = A (dB/drho + B/rho + C)
      ! a4 = A B

      DO i = 1, nrho
         rho = THRIFT_RHO(i)
         a1(i) = A_temp(i+1)*D_der(i)
         a2(i) = A_temp(i+1)*(1/rho*B_der(i) - B_temp(i+1)/(rho**2) + C_der(i))    
         a3(i) = A_temp(i+1)*(B_der(i) + B_temp(i+1)/rho + C_temp(i+1))                   
         a4(i) = A_temp(i+1)*B_temp(i+1)                             
      END DO

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
      ! First decay I_plasma
      THRIFT_IPLASMA(nrho,mytimestep) = THRIFT_IPLASMA(nrho,itime)*exp(-k/temp2)
      !WRITE(6,*) "IPLASMA EDGE PREV"
      !WRITE(6,*) THRIFT_IPLASMA(nrho,itime)
      !WRITE(6,*) "DECAY"
      !WRITE(6,*) exp(-k/temp2)
      !WRITE(6,*) "IPLASMA EDGE THIS"
      !WRITE(6,*) THRIFT_IPLASMA(nrho,mytimestep)
      !WRITE(6,*) "ISOURCE EDGE THIS"
      !WRITE(6,*) THRIFT_ISOURCE(nrho,mytimestep)
      ! I_total
      temp1 = THRIFT_IPLASMA(nrho,mytimestep)+THRIFT_ISOURCE(nrho,mytimestep)
      !WRITE(6,*) "I EDGE THIS"
      !WRITE(6,*) temp1

      !! Calculate uedge for this timestep
      !t = THRIFT_T(itime) ! t = previous sim time (or current sim time if mytimestep=1)
      !CALL get_prof_etapara(THRIFT_RHO(INT(nrho/2)),t,etapara)
      !CALL get_prof_etapara(THRIFT_RHO(nrho),t,etapara)
      !CALL get_prof_pprime(THRIFT_RHO(nrho),t,pprime) 
      !temp2 = (-8*THRIFT_UEDGE(1)+9*THRIFT_UGRID(nrho,1)-THRIFT_UGRID(nrho-1,1))/(3*h) ! du/drho at edge
      ! uedge [ current timestep ] =  u - dt*mu0/(2*Phi_a)*eta/Lext*dV/dphi* ...
      ! [(p' + <B^2>/mu0)*u + <B^2>/mu0*du/drho - <J.B>]  [ edge, preceding timestep ]
      !THRIFT_UEDGE(2) = THRIFT_UEDGE(1)-k*mu0/(2*THRIFT_PHIEDGE(1))*etapara/temp1*THRIFT_VP(nrho+2,1) * &  
      !      (((pprime + THRIFT_BSQAV(nrho+2,1)/mu0)*THRIFT_UEDGE(1)) + THRIFT_BSQAV(nrho+2,1)/mu0*temp2 &               
      !      - THRIFT_JSOURCE(nrho,itime)*THRIFT_BAV(nrho+2,1))   

      ! Populate of tridiagonal matrix and RHS; AI,CI of size nrho-1, BI,DI of size nrho
      ! (Do most of the work in one loop and fix mistakes afterwards)
      DO i = 1, nrho-1
         CI(i) = a3(i)/(2*h) +a4(i)/(h**2)   ! Upper diagonal (Correct)
         AI(i) = -a3(i)/(2*h)+a4(i)/(h**2)   ! Lower diagonal ((nrho-1) wrong)
         BI(i) = a2(i)-2*a4(i)/(h**2)-1/k    ! Middle diagonal ((1,nrho) wrong)
         DI(i) = -THRIFT_UGRID(i,1)/k-a1(i)  ! Right-hand side (B(nrho) wrong)
      END DO
      !AI(nrho-1)= -a3(nrho)/(3*h)+2*a4(nrho)/(h**2)          ! Lower diagonal fixed
        BI(1)    = a2(1)-a3(1)/(2*h)-a4(1)/(h**2)-1/k  ! Middle diagonal half fixed
      !  BI(nrho) = a2(nrho)-a3(nrho)/h+5*a4(nrho)/(h**2)-1/k! Middle diagonal fixed
      !  DI(nrho) = -THRIFT_UGRID(nrho,1)/k-a1(nrho) &
      !- THRIFT_UEDGE(2)*(4*a3(nrho)/(3*h)+16*a4(nrho)/(5*h**2)) ! Fix B(nrho)

      !s11 = BI(nrho); s12 = AI(nrho-1); temp2 = DI(nrho) ! For writing to screen
      !temp1 = -a4(nrho)/(5*h**2) ! Element at (nrho,nrho-2)
      !temp1 = temp1/AI(nrho-2) ! Row operation: [NRHO] -> [NRHO]-temp1*[NRHO-1] 
      !AI(nrho-1) = AI(nrho-1) - temp1*BI(nrho-1) ! (AI is of size nrho-1)
      !BI(nrho)   = BI(nrho)   - temp1*CI(nrho-1)
      !DI(nrho)   = DI(nrho)   - temp1*DI(nrho-1)

      ! Set u_nrho^n+1 = u_nrho+1^n+1, should be fine with good enough grid resolution
      BI(nrho) = 1; AI(nrho-1) = 0; DI(nrho) = temp1*mu0/(2*THRIFT_RHO(nrho)*THRIFT_PHIEDGE(2));

      ! Solve system of equations
      CALL solve_tdm(AI,BI,CI,DI,THRIFT_UGRID(:,2))

      IF (lverbj) THEN
         WRITE(6,*) '==============================================================================='
         WRITE(6,*) 'CALCULATING UEDGE'
         WRITE(6,'(A10,2X,ES13.5)') 'RMAJOR', THRIFT_RMAJOR(nrho+2,1)
         WRITE(6,'(A10,2X,ES13.5)') 'AMINOR',THRIFT_AMINOR(nrho+2,1)
         WRITE(6,'(A10,2X,ES13.5)') 'LEXT',temp1
         WRITE(6,'(A10,2X,ES13.5)') 'PHIEDGE',THRIFT_PHIEDGE(1)
         WRITE(6,'(A10,2X,ES13.5)') 'ETAPARA',etapara
         WRITE(6,'(A10,2X,ES13.5)') 'PPRIME',pprime
         WRITE(6,'(A10,2X,ES13.5)') 'DU/DRHO',temp2
         WRITE(6,'(A10,2X,ES13.5)') 'VP',THRIFT_VP(nrho+2,1)
         WRITE(6,'(A10,2X,ES13.5)') 'BSQAV',THRIFT_BSQAV(nrho+2,1)
         WRITE(6,'(A10,2X,ES13.5)') '<J.B>',THRIFT_JSOURCE(nrho,itime)*THRIFT_BAV(nrho+2,1)
         WRITE(6,'(A10,2X,ES13.5)') 'prev.UEDGE',THRIFT_UEDGE(1)
         WRITE(6,'(A10,2X,ES13.5)') 'this.UEDGE',THRIFT_UEDGE(2)
         WRITE(6,*) '==============================================================================='
         WRITE(6,*) 'UGRID (NEXT ITERATION)'
         WRITE(6,*) 'RHO      U'
         DO i = 1, nrho
            WRITE(6,'(F5.3, 1X, ES13.5)') THRIFT_RHO(i), THRIFT_UGRID(i,2)
         END DO
         WRITE(6, '(A5,1X,ES13.5)') 'EDGE',THRIFT_UEDGE(2)
      END IF     
      
      IF (lverbj) THEN
      !   WRITE(6,*) '==============================================================================='
      !   WRITE(6,*)' MATRIX ELEMENT AT (nrho,nrho-2)'
      !   WRITE(6,*) temp1
         WRITE(6,*) '==============================================================================='
         WRITE(6,*)'  i         LOWER           MAIN          UPPER            RHS       SOLUTION'
         WRITE(6,*)''
         WRITE(6,'(I4, 1X,5(2X,ES13.5))') 1, 0.0, BI(1), CI(1), DI(1),THRIFT_UGRID(1,2)
         DO i = 2, nrho-1
            WRITE(6,'(I4, 1X, 5(ES13.5,2X))') i, AI(i-1), BI(i), CI(i), DI(i),THRIFT_UGRID(i,2)
         END DO
         WRITE(6,'(I4, 1X, 5(ES13.5,2X))') nrho, AI(nrho-1),BI(nrho), 0.0,DI(nrho), THRIFT_UGRID(nrho,2)
         !WRITE(6,'(A5,60X,ES13.5)') 'EDGE',THRIFT_UEDGE(2)
      END IF

      ! ITOTAL: u = mu0 I / phip => I = 2*phi_a*rho*u/mu0
      THRIFT_I(:,mytimestep) = 2*eq_phiedge/mu0*(THRIFT_RHO*THRIFT_UGRID(:,2))
      CALL curden_to_curtot(THRIFT_JSOURCE,THRIFT_AMINOR,THRIFT_ISOURCE,mytimestep)

      ! ISOURCE(j) = sum_i=1->j JSOURCE(i)*Delta(i)A = ISOURCE(j-1)+JSOURCE(j)*Delta(j)A
      ! A(i) = pi*aminor(i)^2, Delta(i)A = A(i) - A(i-1) (AMINOR is in [0,1])
      !C_der = 0;
      !DO i = 1, nrho
      !   C_der(i) = THRIFT_JSOURCE(i,mytimestep)*pi*(THRIFT_AMINOR(i+1,2)**2-THRIFT_AMINOR(i,2)**2)
      !   IF (i /= 1) C_der(i) = C_der(i) + C_der(i-1)
      !END DO

      ! IPLASMA = ITOTAL - ISOURCE
      THRIFT_IPLASMA(:,mytimestep) = THRIFT_I(:,mytimestep) - THRIFT_ISOURCE(:,mytimestep)
      !WRITE(6,*) "THRIFT I THIS"
      !WRITE(6,*) THRIFT_I(:,mytimestep)
      !WRITE(6,*) "THRIFT ISOURCE THIS"
      !WRITE(6,*) THRIFT_ISOURCE(:,mytimestep)
      !WRITE(6,*) "THRIFT IPLASMA THIS"
      !WRITE(6,*) THRIFT_IPLASMA(:,mytimestep)

      ! JPLASMA(i) = (IPLASMA(i)-IPLASMA(i-1))/(pi*(aminor(i)^2-aminor(i-1)^2))
      CALL curtot_to_curden(THRIFT_IPLASMA,THRIFT_AMINOR,THRIFT_JPLASMA,mytimestep)

      ! Subtract change in JSOURCE
      THRIFT_JPLASMA(:,mytimestep) = THRIFT_JPLASMA(:,mytimestep)-(THRIFT_JSOURCE(:,mytimestep)-THRIFT_JSOURCE(:,itime))
      !WRITE(6,*) "DELTA JSOURCE THIS"
      !WRITE(6,*) THRIFT_JSOURCE(:,mytimestep)-THRIFT_JSOURCE(:,itime)
      !WRITE(6,*) "THRIFT JPLASMA THIS"
      !WRITE(6,*) THRIFT_JPLASMA(:,mytimestep)

      !temp2 = 0;
      !DO i = 1, nrho
      !   IF (i /= 1) temp2 = D_der(i-1) ! IPLASMA
      !   THRIFT_JPLASMA(i,mytimestep) = (D_der(i)-temp2)/(pi*(THRIFT_AMINOR(i+1,2)**2-THRIFT_AMINOR(i,2)**2))
      !END DO

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