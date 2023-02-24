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
      REAL(rprec) :: rho,s,h,k,t,s11,s12,iota,Bav,Bsqav,vp,etapara,pprime,&
                     temp1,temp2,source_0,source_edge,Aminor,Rmajor
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: A_temp, B_temp, C_temp, D_temp, &
                                                B_der, C_der, D_der, &
                                                a1, a2, a3, a4,         &
                                                DU, D, DL, B
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
      ! Check to make sure delta t != 0
      IF (tstart==0 .and. mytimestep==1) THEN
         THRIFT_JPLASMA(:,mytimestep)=-THRIFT_JSOURCE(:,mytimestep)
         GOTO 1000
      END IF

      ! Allocations
      ALLOCATE(A_temp(nrho+2), B_temp(nrho+2), C_temp(nrho+2), D_temp(nrho+2), &
               B_der(nrho), C_der(nrho), D_der(nrho), &
               a1(nrho), a2(nrho), a3(nrho), a4(nrho), &
               DU(nrho-1), D(nrho), DL(nrho-1), B(nrho))

      A_temp = 0; B_temp = 0; C_temp = 0; D_temp = 0;
      B_der = 0; C_der = 0; D_der = 0;
      a1 = 0; a2 = 0; a3 = 0; a4 = 0;
      DU = 0; D = 0; DL = 0; B = 0;

      t = THRIFT_T(mytimestep)

      ! Populate A,B,C,D
      IF (lverblatin) THEN
         WRITE(6,*)'==============================================================================='
         WRITE(6,*)'CALCULATING COEFFICIENTS A,B,C,D'
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*) 'RHO      ETAPARA     DV/DPHI     DP/DRHO     BAV      BSQAV        S11'
      END IF

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

         ! Interpolate source current density at magnetic axis and edge
         source_0 = (3*THRIFT_JSOURCE(1,mytimestep)-THRIFT_JSOURCE(2,mytimestep))/2
         source_edge = (-THRIFT_JSOURCE(nrho-1,mytimestep)+3*THRIFT_JSOURCE(nrho,mytimestep))/2

         IF (i==nrho+2) THEN ! etapara breaks at rho = 1
            CALL get_prof_etapara(THRIFT_RHO(nrho), t,etapara)
            h = source_edge
         ELSE 
            CALL get_prof_etapara(rho, t,etapara)
            IF (i==1) THEN 
               h = source_0
            ELSE
               h = THRIFT_JSOURCE(i-1,mytimestep)
            END IF
         END IF
         CALL EZspline_interp(vp_spl,rho,vp,ier)
         CALL get_prof_pprime(rho,t,pprime)
         CALL get_equil_Bav(s,Bav,Bsqav,ier)
         CALL get_equil_sus(s,s11,h,h,h,ier)
         temp1 = 2*etapara*vp ! temp1 <- 2 eta dV/dPhi 
         IF (i /= 1)  A_temp(i) = s11/(4*rho*eq_phiedge) ! A not necessary at axis (no derivative required)
         B_temp(i) = temp1*Bsqav/mu0! 2 eta dV/dPhi <B^2>/mu_0
         C_temp(i) = temp1*pprime   ! 2 eta dV/dPhi dp/drho
         D_temp(i) = -temp1*h*Bav   ! -2 eta dV/dPhi <J.B>
         IF (lverblatin) WRITE(6,'(F5.3,6(1X,ES10.3))') rho, etapara, vp, pprime, bav, bsqav, s11
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

      ! Timesteps and grid spacing 
      IF (mytimestep == 1) THEN
         k = tstart
      ELSE
         k = THRIFT_T(mytimestep)-THRIFT_T(mytimestep-1)
      END IF
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
         WRITE(6,*)'CALCULATING COEFFICIENTS ALPHA 1,2,3,4'
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
      rho = 1 
      s = rho*rho
      ier = 0
      CALL get_equil_Rmajor(s,Rmajor,Bav,Aminor,ier)
      temp1 = mu0*Rmajor*(log(8*Rmajor/Aminor) - 2 + 0.25) ! temp1 <- L_ext
      ! Calculate u_edge^mytimestep 
      CALL get_prof_etapara(THRIFT_RHO(nrho-1),t,etapara)  
      CALL EZspline_interp(vp_spl,rho,vp,ier)
      CALL get_equil_Bav(s,Bav,Bsqav,ier)
      CALL get_prof_pprime(rho,t,pprime) 
      temp2 = (-8*THRIFT_UEDGE(1)+9*THRIFT_UGRID(nrho,1)-THRIFT_UGRID(nrho-1,1))/(3*h) ! du/drho at edge
      ! u_edge^mytimestep =  u - dt* mu0/(2Phi_a)*eta/Lext*dV/dphi*((p' + <B^2>/mu0)*u+<B^2>/mu0 * du/drho- <J.B>)
      ! RHS : u evaluated at previous timestep
      THRIFT_UEDGE(2) = THRIFT_UEDGE(1) - k*mu0/(2*eq_phiedge)*etapara/temp1*vp * &  
            (((pprime + Bsqav/mu0)*THRIFT_UEDGE(1)) + Bsqav/mu0*temp2 &               
            - source_edge*Bav)                                   

      ! Populate diagonals and RHS; DU and DL of size nrho-1, D of size nrho
      ! (Do most of the work in one loop and fix mistakes afterwards)
      ! If on the first time iteration, we say there's no enclosed current (i.e. B=0)
      DO i = 1, nrho-1
         DU(i) = a3(i)/(2*h) +a4(i)/(h**2)   ! Upper diagonal (Correct)
         DL(i) = -a3(i)/(2*h)+a4(i)/(h**2)   ! Lower diagonal (DL(nrho-1) wrong)
          D(i) = a2(i)-2*a4(i)/(h**2)-1/k    ! Middle diagonal (D(1,nrho) wrong)
          B(i) = -THRIFT_UGRID(i,1)/k-a1(i)  ! Right-hand side (B(nrho) wrong)
      END DO
      DL(nrho-1)= -a3(nrho)/(3*h)+2*a4(nrho)/(h**2)         ! Lower diagonal fixed
        D(1)    = a2(1)-a3(1)/(2*h)-a4(1)/(h**2)-1/k  ! Middle diagonal half fixed
        D(nrho) = a2(nrho)-a3(nrho)/h+5*a4(nrho)/(h**2)-1/k! Middle diagonal fixed
        B(nrho) = -THRIFT_UGRID(nrho,1)/k-a1(nrho) &
      - THRIFT_UEDGE(2)*(4*a3(nrho)/(3*h)+16*a4(nrho)/(5*h**2)) ! Fix B(nrho)
      s11 = D(nrho)
      s12 = DL(nrho-1)
      temp2 = B(nrho)
      
      ! Matrix is not fully tridiagonal; it has following element in (nrho,nrho-2)
      temp1 = -a4(nrho)/(5*h**2) 
      ! Eliminate that extra non-zero element to get a TDM by performing following row operations
      temp1 = temp1/DL(nrho-2) ! Row operation: [NRHO] -> [NRHO]-temp1*[NRHO-1] (DL is of size nrho-1)
      DL(nrho-1) = DL(nrho-1) - temp1* D(nrho-1) 
      D(nrho)    = D(nrho)    - temp1*DU(nrho-1)
      B(nrho)    = B(nrho)    - temp1* B(nrho-1)

      ! Thomas algorithm
      ! https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm 
      ! a_i = DL, b_i = D, c_i = DU, d_i = B
      ! note that on wiki, indices from a are in [2,nrho] whereas ours are in [1,nrho-1]
      a1 = 0; a2 = 0;
      a1(1) = DU(1)/D(1) !c_1' = c_1/b_1
      a2(1) = B(1)/D(1)  !d_1' = d_1/b_1
      DO i = 2, nrho
         IF (i /= nrho) a1(i) = DU(i)/(D(i)-DL(i-1)*a1(i-1))  ! c_i' =              c_i/(b_i - a_i*c_i-1')
         a2(i) = (B(i)-DL(i-1)*a2(i-1))/(D(i)-DL(i-1)*a1(i-1))! d_i' = (d_i-a_i*d_i-1')/(b_i - a_i*c_i-1')
      END DO
      THRIFT_UGRID(nrho,2) = a2(nrho) ! x_n = d_n'
      DO i = nrho-1, 1, -1
         THRIFT_UGRID(i,2) = a2(i)-a1(i)*THRIFT_UGRID(i+1,2) ! x_i = d_i' - c_i'*x_i+1
      END DO
      
      ! Store u in variable
      IF (lverbmat) THEN
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)' MATRIX ELEMENT AT (nrho,nrho-2)'
         WRITE(6,*) temp1
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)'  i         LOWER           MAIN          UPPER            RHS       SOLUTION'
         WRITE(6,*)''
         WRITE(6,'(I4, 1X,A13,5(2X,ES13.5))') 1, '-------------', D(1), DU(1), B(1), THRIFT_UGRID(1,2)
         DO i = 2, nrho-1
            WRITE(6,'(I4, 1X, 5(ES13.5,2X))') i, DL(i-1), D(i), DU(i), B(i), THRIFT_UGRID(i,2)
         END DO
         WRITE(6,'(I4, 1X, 2(ES13.5,2X),A13,2(2X,ES13.5))') nrho, DL(nrho-1),D(nrho), '-------------',&
                                                             B(nrho), THRIFT_UGRID(nrho,2)
         WRITE(6,'(A5,60X,ES13.5)') 'EDGE',THRIFT_UEDGE(2)
         WRITE(6,'(A5, 1X, 2(ES13.5,2X),A13,2X,ES13.5)') 'TEMP', s12, s11,'-------------', temp2
         WRITE(6,*)'-------------------------------------------------------------------------------'
      END IF
      ! B = mu0 I / phip => I = phip*B/mu0 = 2*phi_a*rho*B/mu0
      B = 2*eq_phiedge/mu0*(THRIFT_RHO*THRIFT_UGRID(:,2))

      IF (lverbpost) THEN
         WRITE(6,*)'==============================================================================='
         WRITE(6,*)'POST SOLVING EQUATIONS'
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)' I_TOTAL AT CURRENT TIMESTEP'
         WRITE(6,*) ''
         WRITE(6,*) B(1:9)
         WRITE(6,*) '...'
         WRITE(6,*) B(nrho-8:nrho)
         WRITE(6,*)''
      END IF    

      ! Calculate I_source (a1)
      a1 = 0; temp1 = 0
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
         a1(i) = THRIFT_JSOURCE(i,mytimestep)*pi*(Aminor**2-temp1**2)
         IF (i /= 1) a1(i) = a1(i) + a1(i-1)
      END DO

      ! I_plasma = I_total - I_source
      B = B - a1 
      IF (lverbpost) THEN
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)' I_SOURCE AT CURRENT TIMESTEP'
         WRITE(6,*)''
         WRITE(6,*) a1(1:9)
         WRITE(6,*) '...'
         WRITE(6,*) a1(nrho-8:nrho)
         WRITE(6,*)''
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)' I_PLASMA AT CURRENT TIMESTEP'
         WRITE(6,*)''
         WRITE(6,*) B(1:9)
         WRITE(6,*) '...'
         WRITE(6,*) B(nrho-8:nrho)
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
         IF (i /= 1) temp2 = B(i-1) ! I_plasma(i)
         THRIFT_JPLASMA(i,mytimestep) = (B(i)-temp2)/(pi*(Aminor**2-temp1**2))
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
               DU, D, DL, B)

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