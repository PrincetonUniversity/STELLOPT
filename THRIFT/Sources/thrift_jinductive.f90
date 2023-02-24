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
      LOGICAL :: lverblatin, lverbalpha, lverbmat,lverbpost
      REAL(rprec) :: rho,s,h,k,t,s11,s12,etapara,pprime,&
                     temp1,temp2,source_0,source_edge,Aminor,Rmajor
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: A_temp, B_temp, C_temp, D_temp, &
                                                B_der, C_der, D_der, &
                                                a1, a2, a3, a4, &
                                                AI, BI, CI, DI
      TYPE(EZspline1_r8) :: f_spl
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
!     The general idea here is to use equation (22) in form to update
!     THRIFT_JPLASMA.  Where THRIFT_JSOURCE is the J in <J_s.B>.

      ! Choose which items to make verbose
      lverblatin  = .false.
      lverbalpha  = .false.
      lverbmat    = .true.
      lverbpost   = .false.

      THRIFT_I        = 0; THRIFT_IBOOT    = 0; THRIFT_IPLASMA  = 0
      THRIFT_IECCD    = 0; THRIFT_INBCD    = 0; THRIFT_IOHMIC   = 0

      ! Check to make sure we're not zero beta
      IF (eq_beta == 0) THEN
         IF (mytimestep /= 1) THRIFT_JPLASMA(:,mytimestep) = THRIFT_JPLASMA(:,mytimestep-1)
         GOTO 1000
      END IF

      ! Store magnetic variables
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

      ! Check to make sure delta t != 0
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
      a1 = 0, a2 = 0, a3 = 0, a4 = 0;
      AI = 0; BI = 0; CI = 0; DI = 0;

      IF (mytimestep==1) THEN
         itime = mytimestep
         k = tstart
      ELSE
         itime = mytimestep-1
         k = THRIFT_T(mytimestep)-THRIFT_T(itime)
      END IF
      t = THRIFT_T(itime)

      ! Populate A,B,C,D
      IF (lverblatin) THEN
         WRITE(6,*)'==============================================================================='
         WRITE(6,*)'CALCULATING COEFFICIENTS A,B,C,D'
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*) 'RHO      ETAPARA     DV/DPHI     DP/DRHO     BAV      BSQAV        S11'
      END IF

      
      ! Interpolate source current density at magnetic axis and edge
      source_0 = (3*THRIFT_JSOURCE(1,itime)-THRIFT_JSOURCE(2,itime))/2
      source_edge = (-THRIFT_JSOURCE(nrho-1,itime)+3*THRIFT_JSOURCE(nrho,itime))/2

      DO i = 1, nrho+2
         ! calculate rho, etapara, current
         IF (i==1) THEN 
            rho = 0
            CALL get_prof_etapara(rho,t,etapara)
            h = source_0
         ELSE IF (i==nrho+2) THEN ! etapara breaks when rho=1
            rho = 1
            CALL get_prof_etapara(THRIFT_RHO(nrho),t,etapara)
            h = source_edge
         ELSE
            rho = THRIFT_RHO(i-1)
            CALL get_prof_etapara(rho,t,etapara)
            h = THRIFT_JSOURCE(i-1,itime)
         END IF
         CALL get_prof_pprime(rho,t,pprime)
         temp1 = 2*etapara*THRIFT_VP(i,1) ! temp1 <- 2 eta dV/dPhi 
         IF (i /= 1) A_temp(i) = THRIFT_S11(i,1)/(4*rho*THRIFT_PHIEDGE(1)**2) ! S11/(4 rho phi_a^2)
         B_temp(i) = temp1*THRIFT_BSQAV(i,1)/mu0! 2 eta dV/dPhi <B^2>/mu_0
         C_temp(i) = temp1*pprime               ! 2 eta dV/dPhi dp/drho
         D_temp(i) = -temp1*h*THRIFT_BAV(i,1)   ! -2 eta dV/dPhi <J.B>
         IF (lverblatin) WRITE(6,'(F5.3,6(1X,ES10.3))') &
         rho, etapara, THRIFT_VP(i,1), pprime, THRIFT_BAV(i,1), THRIFT_BSQAV(i,1), THRIFT_S11(i,1)
      END DO

      IF (lverblatin) THEN
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)' A'
         WRITE(6,*) '' 
         WRITE(6,*) A_temp(1:9)
         WRITE(6,*) '...'
         WRITE(6,*) A_temp(nrho-6:nrho+2)
         WRITE(6,*)''
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)' B'
         WRITE(6,*) '' 
         WRITE(6,*) B_temp(1:9)
         WRITE(6,*) '...'
         WRITE(6,*) B_temp(nrho-6:nrho+2)
         WRITE(6,*)''
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)' C'
         WRITE(6,*) '' 
         WRITE(6,*) C_temp(1:9)
         WRITE(6,*) '...'
         WRITE(6,*) C_temp(nrho-6:nrho+2)
         WRITE(6,*)''
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)' D'
         WRITE(6,*) '' 
         WRITE(6,*) D_temp(1:9)
         WRITE(6,*) '...'
         WRITE(6,*) D_temp(nrho-6:nrho+2)
         WRITE(6,*)''
      END IF

      ! Grid spacing 
      h = THRIFT_RHO(2)-THRIFT_RHO(1)  ! h <- drho

      ! Visualisation of different grids
      !j=1 2    3    4    5    6    7
      ! |  |    |    |    |    |    |     ABCD_temp on j
      !    |    |    |    |    |    |     a1,2,3,4  on i
      !   i=1   2    3    4    5    6     j = i+1 

      ! Calculate derivatives of ABCD
      ! dY/drho(i) = [Y(j+1)-Y(j-1)]/2h = [Y(i+2)-Y(i)]/2h
      DO i = 2, nrho-1
         B_der(i) = (B_temp(i+2)-B_temp(i))/(2*h)
         C_der(i) = (C_temp(i+2)-C_temp(i))/(2*h)
         D_der(i) = (D_temp(i+2)-D_temp(i))/(2*h)
      END DO
      ! Near magnetic axis
      ! dY/drho(i=1) = [Y(i=3) + 3*Y(i=2) - 4*Y(i=1)]/3h
      B_der(1) = (B_temp(3)+3*B_temp(2)-4*B_temp(1))/(3*h)
      C_der(1) = (C_temp(3)+3*C_temp(2)-4*C_temp(1))/(3*h)
      D_der(1) = (D_temp(3)+3*D_temp(2)-4*D_temp(1))/(3*h)
      ! Near plasma edge
      ! dY/drho(i=nrho) = [4*Y(i=nrho+2) - 3*Y(i=nrho+1) - Y(i=nrho)]/3h
      B_der(nrho) = (4*B_temp(nrho+2)-3*B_temp(nrho+1)-B_temp(nrho))/(3*h)
      C_der(nrho) = (4*C_temp(nrho+2)-3*C_temp(nrho+1)-C_temp(nrho))/(3*h)
      D_der(nrho) = (4*D_temp(nrho+2)-3*D_temp(nrho+1)-D_temp(nrho))/(3*h)
      
      IF (lverblatin) THEN
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)' DERIV B'
         WRITE(6,*) '' 
         WRITE(6,*) B_der(1:9)
         WRITE(6,*) '...'
         WRITE(6,*) B_der(nrho-8:nrho)
         WRITE(6,*)''
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)' DERIV C'
         WRITE(6,*) '' 
         WRITE(6,*) C_der(1:9)
         WRITE(6,*) '...'
         WRITE(6,*) C_der(nrho-8:nrho)
         WRITE(6,*)''
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)' DERIV D'
         WRITE(6,*) '' 
         WRITE(6,*) D_der(1:9)
         WRITE(6,*) '...'
         WRITE(6,*) D_der(nrho-8:nrho)
         WRITE(6,*)''
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

      IF (lverbalpha) THEN
         WRITE(6,*)'==============================================================================='
         WRITE(6,*)'COEFFICIENTS ALPHA 1,2,3,4 PREVIOUS TIMESTEP'
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)' ALPHA_1'
         WRITE(6,*) '' 
         WRITE(6,*) a1(1:9)
         WRITE(6,*) '...'
         WRITE(6,*) a1(nrho-8:nrho)
         WRITE(6,*)''
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)' ALPHA_2'
         WRITE(6,*) '' 
         WRITE(6,*) a2(1:9)
         WRITE(6,*) '...'
         WRITE(6,*) a2(nrho-8:nrho)
         WRITE(6,*)''
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)' ALPHA_3'
         WRITE(6,*) '' 
         WRITE(6,*) a3(1:9)
         WRITE(6,*) '...'
         WRITE(6,*) a3(nrho-8:nrho)
         WRITE(6,*)''
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)' ALPHA_4'
         WRITE(6,*) '' 
         WRITE(6,*) a4(1:9)
         WRITE(6,*) '...'
         WRITE(6,*) a4(nrho-8:nrho)
         WRITE(6,*)''
      END IF
      !

      ! Calculate u_edge^mytimestep from u_edge^mytimestep-1, u_grid^mytimestep-1
      ! First calculate L_ext = mu0*R_eff( ln( 8 R_eff/r_eff )-2 + F_shaping)      
      temp1 = mu0*THRIFT_RMAJOR(nrho+2,1)*&
         (log(8*THRIFT_RMAJOR(nrho+2,1)/THRIFT_AMINOR(nrho+2,1)) - 2 + 0.25) ! temp1 <- L_ext
      ! Calculate u_edge^mytimestep 
      CALL get_prof_etapara(THRIFT_RHO(nrho-1),t,etapara)  
      CALL get_prof_pprime(1.0,t,pprime) 
      temp2 = (-8*THRIFT_UEDGE(1)+9*THRIFT_UGRID(nrho,1)-THRIFT_UGRID(nrho-1,1))/(3*h) ! du/drho at edge
      ! u_edge^mytimestep =  u - dt* mu0/(2Phi_a)*eta/Lext*dV/dphi*((p' + <B^2>/mu0)*u+<B^2>/mu0 * du/drho- <J.B>)
      THRIFT_UEDGE(2) = THRIFT_UEDGE(1) - k*mu0/(2*THRIFT_PHIEDGE(1))*etapara/temp1*THRIFT_VP(nrho+2,1) * &  
            (((pprime + THRIFT_BSQAV(nrho+2,1)/mu0)*THRIFT_UEDGE(1)) + THRIFT_BSQAV(nrho+2,1)/mu0*temp2 &               
            - source_edge*THRIFT_BAV(nrho+2,1))                                   

      ! Populate diagonals and RHS; AI,CI of size nrho-1, BI, DI of size nrho
      ! (Do most of the work in one loop and fix mistakes afterwards)
      DO i = 1, nrho-1
         CI(i) = a3(i)/(2*h) +a4(i)/(h**2)   ! Upper diagonal (Correct)
         AI(i) = -a3(i)/(2*h)+a4(i)/(h**2)   ! Lower diagonal ((nrho-1) wrong)
         BI(i) = a2(i)-2*a4(i)/(h**2)-1/k    ! Middle diagonal ((1,nrho) wrong)
         DI(i) = -THRIFT_UGRID(i,1)/k-a1(i)  ! Right-hand side (B(nrho) wrong)
      END DO
      AI(nrho-1)= -a3(nrho)/(3*h)+2*a4(nrho)/(h**2)         ! Lower diagonal fixed
        BI(1)    = a2(1)-a3(1)/(2*h)-a4(1)/(h**2)-1/k  ! Middle diagonal half fixed
        BI(nrho) = a2(nrho)-a3(nrho)/h+5*a4(nrho)/(h**2)-1/k! Middle diagonal fixed
        DI(nrho) = -THRIFT_UGRID(nrho,1)/k-a1(nrho) &
      - THRIFT_UEDGE(2)*(4*a3(nrho)/(3*h)+16*a4(nrho)/(5*h**2)) ! Fix B(nrho)

      s11 = BI(nrho); s12 = AI(nrho-1); temp2 = DI(nrho)
      
      temp1 = -a4(nrho)/(5*h**2) ! Element at (nrho,nrho-2)
      temp1 = temp1/AI(nrho-2) ! Row operation: [NRHO] -> [NRHO]-temp1*[NRHO-1] 
      AI(nrho-1) = AI(nrho-1) - temp1*BI(nrho-1) ! (AI is of size nrho-1)
      BI(nrho)   = BI(nrho)   - temp1*CI(nrho-1)
      DI(nrho)   = DI(nrho)   - temp1*DI(nrho-1)

      ! Thomas algorithm: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm 
      ! Forward sweep
      C_temp = 0; D_temp = 0;
      C_temp(1) = CI(1)/BI(1) !c_1' = c_1/b_1
      D_temp(1) = DI(1)/BI(1) !d_1' = d_1/b_1
      DO i = 2, nrho ! A on wiki goes from 2,n; A here goes from 1, n-1 -> shift by -1
         IF (i /= nrho) C_temp(i) =        CI(i)/(BI(i)-AI(i-1)*C_temp(i-1))! c_i' =              c_i/(b_i - a_i*c_i-1')
         D_temp(i) = (DI(i)-AI(i-1)*D_temp(i-1))/(BI(i)-AI(i-1)*C_temp(i-1))! d_i' = (d_i-a_i*d_i-1')/(b_i - a_i*c_i-1')
      END DO
      ! Back substitution
      THRIFT_UGRID(nrho,2) = D_temp(nrho) ! x_n = d_n'
      DO i = nrho-1, 1, -1
         THRIFT_UGRID(i,2) = D_temp(i)-C_temp(i)*THRIFT_UGRID(i+1,2) ! x_i = d_i' - c_i'*x_i+1
      END DO
      
      IF (lverbmat) THEN
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)' MATRIX ELEMENT AT (nrho,nrho-2)'
         WRITE(6,*) temp1
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)'  i         LOWER           MAIN          UPPER            RHS       SOLUTION'
         WRITE(6,*)''
         WRITE(6,'(I4, 1X,A13,5(2X,ES13.5))') 1, '------------', BI(1), CI(1), DI(1), THRIFT_UGRID(1,2)
         DO i = 2, nrho-1
            WRITE(6,'(I4, 1X, 5(ES13.5,2X))') i, AI(i-1), BI(i), CI(i), DI(i), THRIFT_UGRID(i,2)
         END DO
         WRITE(6,'(I4, 1X, 2(ES13.5,2X),A13,2(2X,ES13.5))') nrho, AI(nrho-1),BI(nrho), '------------',&
                                                             DI(nrho), THRIFT_UGRID(nrho,2)
         WRITE(6,'(A5,60X,ES13.5)') 'EDGE',THRIFT_UEDGE(2)
         WRITE(6,'(A4, 1X, 2(ES13.5,2X),A13,2X,ES13.5)') 'TEMP', s12, s11,'------------', temp2
         WRITE(6,*)'-------------------------------------------------------------------------------'
      END IF

      ! Calculate Itotal: u = mu0 I / phip => I = phip*u/mu0 = 2*phi_a*rho*u/mu0
      B_temp = 0
      B_temp = 2*eq_phiedge/mu0*(THRIFT_RHO*THRIFT_UGRID(:,2))

      IF (lverbpost) THEN
         WRITE(6,*)'==============================================================================='
         WRITE(6,*)'POST SOLVING EQUATIONS'
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)' I_TOTAL AT CURRENT TIMESTEP'
         WRITE(6,*) ''
         WRITE(6,*) B_temp(1:9)
         WRITE(6,*) '...'
         WRITE(6,*) B_temp(nrho-8:nrho)
         WRITE(6,*)''
      END IF    

      ! Calculate I_source (a1)
      A_temp = 0; temp1 = 0
      DO i = 1, nrho
         rho = THRIFT_RHO(i)
         s = rho * rho
         ier = 0
         CALL get_equil_Rmajor(s,h,h,Aminor,ier)  ! aminor_i
         IF (i /= 1) THEN
            rho = THRIFT_RHO(i-1)
            s = rho*rho
            ier = 0
            CALL get_equil_Rmajor(s,h,h,temp1,ier) ! temp1 = aminor_i-1
         END IF
         A_temp(i) = THRIFT_JSOURCE(i,mytimestep)*pi*(Aminor**2-temp1**2)
         IF (i /= 1) A_temp(i) = A_temp(i) + A_temp(i-1)
      END DO

      ! I_plasma = I_total - I_source
      B_temp = B_temp - A_temp 
      IF (lverbpost) THEN
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)' I_SOURCE AT CURRENT TIMESTEP'
         WRITE(6,*)''
         WRITE(6,*) A_temp(1:9)
         WRITE(6,*) '...'
         WRITE(6,*) A_temp(nrho-8:nrho)
         WRITE(6,*)''
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)' I_PLASMA AT CURRENT TIMESTEP'
         WRITE(6,*)''
         WRITE(6,*) B_temp(1:9)
         WRITE(6,*) '...'
         WRITE(6,*) B_temp(nrho-8:nrho)
         WRITE(6,*)''
      END IF    

      ! JPLASMA(i) = (IPLASMA(i)-IPLASMA(i-1))/(pi*(a(i)**2-a(i-1)**2))
      temp1 = 0; temp2 = 0;
      DO i = 1, nrho
         rho = THRIFT_RHO(i)
         s = rho*rho
         ier = 0
         CALL get_equil_Rmajor(s,h,h,Aminor,ier)
         IF (i /= 1) THEN
            rho = THRIFT_RHO(i-1)
            s = rho*rho
            ier = 0
            CALL get_equil_Rmajor(s,h,h,temp1,ier)
         END IF
         IF (i /= 1) temp2 = B_temp(i-1) ! I_plasma(i-1)
         THRIFT_JPLASMA(i,mytimestep) = (B_temp(i)-temp2)/(pi*(Aminor**2-temp1**2))
      END DO

      IF (lverbpost) THEN
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)' J_PLASMA AT CURRENT TIMESTEP'
         WRITE(6,*)''
         WRITE(6,*) THRIFT_JPLASMA(1:9,mytimestep)
         WRITE(6,*) '...'
         WRITE(6,*) THRIFT_JPLASMA(nrho-8:nrho,mytimestep)
         WRITE(6,*)''         
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)' J_SOURCE AT CURRENT TIMESTEP'
         WRITE(6,*)''
         WRITE(6,*) THRIFT_JSOURCE(1:9,mytimestep)
         WRITE(6,*) '...'
         WRITE(6,*) THRIFT_JSOURCE(nrho-8:nrho,mytimestep)
         WRITE(6,*)'' 
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)' JINDUCTIVE DONE'
         WRITE(6,*)'==============================================================================='
      END IF    
      DEALLOCATE(A_temp, B_temp, C_temp, D_temp, &
               B_der, C_der, D_der, &
               a1, a2, a3, a4, &
               AI, BI, CI, DI)

1000  CONTINUE
      ! Calculate enclosed currents for progress
      temp1 = 0
      DO i = 1, nrho
         rho = THRIFT_RHO(i)
         s = rho * rho
         ier = 0
         CALL get_equil_Rmajor(s,h,h,Aminor,ier)  ! aminor_i
         IF (i /= 1) THEN
            rho = THRIFT_RHO(i-1)
            s = rho*rho
            CALL get_equil_Rmajor(s,h,h,temp1,ier) ! aminor_i-1
         END IF
         THRIFT_IPLASMA= THRIFT_IPLASMA + THRIFT_JPLASMA(i,mytimestep)*pi*(Aminor**2-temp1**2)
         THRIFT_IBOOT  = THRIFT_IBOOT   + THRIFT_JBOOT(i,mytimestep)  *pi*(Aminor**2-temp1**2)
         THRIFT_IECCD  = THRIFT_IECCD   + THRIFT_JECCD(i,mytimestep)  *pi*(Aminor**2-temp1**2)
         THRIFT_INBCD  = THRIFT_INBCD   + THRIFT_JNBCD(i,mytimestep)  *pi*(Aminor**2-temp1**2)
         THRIFT_IOHMIC = THRIFT_IOHMIC  + THRIFT_JOHMIC(i,mytimestep) *pi*(Aminor**2-temp1**2)
      END DO
      THRIFT_I = THRIFT_IPLASMA + THRIFT_IBOOT + THRIFT_IECCD + THRIFT_INBCD + THRIFT_IOHMIC
      RETURN


!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_jinductive