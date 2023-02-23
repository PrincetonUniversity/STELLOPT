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
      LOGICAL :: lverbj
      REAL(rprec) :: rho,s,h,k,t,s11,s12,iota,Bav,Bsqav,vp,etapara,pprime,temp,Aminor,Rmajor
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

      lverbj = .true. 

      ! Check to make sure we're not zero beta
      IF (eq_beta == 0) THEN
         IF (mytimestep == 1) THEN 
            THRIFT_JPLASMA(:,mytimestep) = 0
         ELSE 
            THRIFT_JPLASMA(:,mytimestep) = THRIFT_JPLASMA(:,mytimestep-1)
         END IF
         RETURN
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
      IF (lverbj) THEN
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
         
         IF (i==nrho+2) THEN ! etapara breaks at rho = 1
            CALL get_prof_etapara(THRIFT_RHO(nrho), t,etapara)
         ELSE 
            CALL get_prof_etapara(rho, t,etapara)
         END IF
         CALL EZspline_interp(vp_spl,rho,vp,ier)
         CALL get_prof_pprime(rho,t,pprime)
         CALL get_equil_Bav(s,Bav,Bsqav,ier)
         CALL get_equil_sus(s,s11,h,h,h,ier)
         IF (i /= 1)  A_temp(i) = s11/(4*rho*eq_phiedge) ! A not necessary at axis (no derivative required)
         temp = 2*etapara*vp ! temp <- 2 eta dV/dPhi 
         B_temp(i) = temp*Bsqav/mu0 ! 2 eta dV/dPhi <B^2>/mu_0
         C_temp(i) = temp*pprime  ! 2 eta dV/dPhi dp/drho
         IF (i==1) THEN ! set jsource(1)=jsource(2), jsource(nrho+1)=jsource(nrho)
            h = THRIFT_JSOURCE(1,mytimestep)
         ELSE IF (i==nrho+2) THEN
            h = THRIFT_JSOURCE(nrho,mytimestep)
         ELSE
            h = THRIFT_JSOURCE(i-1,mytimestep)
         END IF
         D_temp(i) = -temp*h*Bav ! -2 eta dV/dPhi <J.B>
         IF (lverbj) WRITE(6,'(F5.3,6(1X,ES10.3))') rho, etapara, vp, pprime, bav, bsqav, s11
      END DO
      IF (lverbj) THEN
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)'A_TEMP'
         WRITE(6,*) A_temp
         WRITE(6,*)''
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)'B_TEMP'
         WRITE(6,*) B_temp
         WRITE(6,*)''
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)'C_TEMP'
         WRITE(6,*) C_temp
         WRITE(6,*)''
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)'D_TEMP'
         WRITE(6,*) D_temp
         WRITE(6,*)''
      END IF

      ! Timesteps and grid spacing 
      k = THRIFT_T(2)-THRIFT_T(1)      ! k <- dt
      h = THRIFT_RHO(2)-THRIFT_RHO(1)  ! h <- drho
      ! NB: If grids are changed to be non-uniform, adjust this accordingly

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
      D_der(nrho) = (4*D_temp(nrho+2)-3*D_temp(nrho+1)-D_temp(nrho))
      
      IF (lverbj) THEN
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)'DERIV B'
         WRITE(6,*) B_der
         WRITE(6,*)''
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)'DERIV C'
         WRITE(6,*) C_der
         WRITE(6,*)''
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)'DERIV D'
         WRITE(6,*) D_der
         WRITE(6,*)''
      END IF

      ! a1 = A dD/drho
      ! a2 = A (1/rho dB/drho - B/rho^2 + dC/drho)
      ! a3 = A (dB/drho + B/rho + C)
      ! a4 = A B
      DO i = 1, nrho
         rho = THRIFT_RHO(i)
         a1(i) = A_temp(i+1)*D_der(i)
         a2(i) = A_temp(i+1)*(1/rho*B_der(i) - B_temp(i+1)/rho**2 + C_der(i))    
         a3(i) = A_temp(i+1)*(B_der(i) + B_temp(i+1)/rho + C_temp(i+1))                   
         a4(i) = A_temp(i+1)*B_temp(i+1)                             
      END DO

      IF (lverbj) THEN
         WRITE(6,*)'==============================================================================='
         WRITE(6,*)'CALCULATING COEFFICIENTS ALPHA 1,2,3,4'
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)'ALPHA_1'
         WRITE(6,*) a1
         WRITE(6,*)''
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)'ALPHA_2'
         WRITE(6,*) a2
         WRITE(6,*)''
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)'ALPHA_3'
         WRITE(6,*) a3
         WRITE(6,*)''
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)'ALPHA_4'
         WRITE(6,*) a4
         WRITE(6,*)''
      END IF

      ! Populate diagonals and RHS; DU and DL of size nrho-1, D of size nrho
      ! (Do most of the work in one loop and fix mistakes afterwards)
      ! If on the first time iteration, we say there's no enclosed current (i.e. B=0)
      DO i = 1, nrho-1
         DU(i) = a3(i)/(2*h) +a4(i)/(h**2)   ! Upper diagonal (Correct)
         DL(i) = -a3(i)/(2*h)+a4(i)/(h**2)   ! Lower diagonal (DL(nrho-1) wrong)
          D(i) = a2(i)-2*a4(i)/(h**2)-1/k    ! Middle diagonal (D(1,nrho) wrong)
          IF (mytimestep /= 1) B(i) = -THRIFT_UGRID(i,mytimestep-1)/k-a1(i) ! Right-hand side (B(nrho) wrong)
      END DO
      DL(nrho-1)= -a3(nrho)/(3*h)+2*a4(nrho)/(h**2)         ! Lower diagonal fixed
        D(1)    = a2(1)-a3(1)/(2*h)-a4(1)/(h**2)-1/k  ! Middle diagonal half fixed
        D(nrho) = a2(nrho)-a3(nrho)/h+5*a4(nrho)/(h**2)-1/k! Middle diagonal fixed
      IF (mytimestep /= 1) B(nrho) = THRIFT_UGRID(nrho,mytimestep-1)/k-a1(nrho) &
      - THRIFT_UEDGE(1)*(4*a3(nrho)/(3*h)+16*a4(nrho)/(5*h**2)) ! Fix B(nrho)

      ! Calculate u = S11 iota + S12 for this timestep, this picard iteration
      ! On the grid
      IF (mytimestep /= ntimesteps) THEN
         DO i = 1, nrho
            rho = THRIFT_RHO(i) 
            s = rho*rho
            ier = 0        
            CALL get_equil_sus(s,s11,s12,h,h,ier) 
            CALL EZspline_interp(iota_spl,rho,iota,ier)
            THRIFT_UGRID(i,mytimestep) = s11*iota+s12
         END DO
      END IF
      ! On the edge
      ! First calculate L_ext = mu0*R_eff( ln( 8 R_eff/r_eff )-2 + F_shaping)      
      rho = 1 
      s = rho*rho
      ier = 0
      CALL get_equil_Rmajor(s,Rmajor,Bav,Aminor,ier)
      temp = mu0*Rmajor*(log(8*Rmajor/Aminor) - 2 + 0.25) ! temp <- L_ext
      ! Now update THRIFT_UEDGE
      CALL get_prof_etapara(THRIFT_RHO(nrho-1),t,etapara)  
      CALL EZspline_interp(vp_spl,rho,vp,ier)
      CALL get_equil_Bav(s,Bav,Bsqav,ier)
      CALL get_prof_pprime(rho,t,pprime) 
      s11 = (-8*THRIFT_UEDGE(1)+9*THRIFT_UGRID(nrho,mytimestep)-THRIFT_UGRID(nrho-1,mytimestep))/(3*h) ! du/drho at edge
      THRIFT_UEDGE(2) = THRIFT_UEDGE(1) - k*mu0/(2*eq_phiedge)*etapara/temp*vp * &  ! u - dt* mu0/(2Phi_a)*eta/Lext*dV/dphi *
            (((pprime + Bsqav/mu0)*THRIFT_UEDGE(1)) + Bsqav/mu0*s11 &               ! ((p' + <B^2>/mu0)*u+<B^2>/mu0 * du/drho
            - THRIFT_JSOURCE(nrho,mytimestep)*Bav)                                  ! - <J.B>)

      IF (lverbj) THEN
         WRITE(6,*)'==============================================================================='
         WRITE(6,*)'CALCULATING U'
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)'LEXT'
         WRITE(6,*) temp
         WRITE(6,*)''
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)'THRIFT_UGRID AT CURRENT TIMESTEP'
         WRITE(6,*) THRIFT_UGRID(:,mytimestep)
         WRITE(6,*)''
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)'THRIFT_UEDGE AT CURRENT TIMESTEP'
         WRITE(6,*) THRIFT_UEDGE
         WRITE(6,*)''
         WRITE(6,*)'==============================================================================='
         WRITE(6,*)'CALCULATING MATRIX'
         WRITE(6,*) 'D(NRHO) (PRE ROW OPERATION)'
         WRITE(6,*) D(nrho)
         WRITE(6,*) '' 
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*) 'DL(NRHO-1) (PRE ROW OPERATION)'
         WRITE(6,*) DL(nrho-1)
         WRITE(6,*) '' 
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*) 'B(NRHO) (PRE ROW OPERATION)'
         WRITE(6,*) B(nrho)
         WRITE(6,*) '' 
      END IF
      
      ! Matrix is not fully tridiagonal; it has following element in (nrho,nrho-2)
      temp = -a4(nrho)/(5*h**2) 
      ! Eliminate that extra non-zero element to get a TDM by performing following row operations
      temp = temp/DL(nrho-2) ! Row operation: [NRHO] -> [NRHO]-temp*[NRHO-1] (DL is of size nrho-1)
      DL(nrho-1) = DL(nrho-1) - temp* D(nrho-1) 
      D(nrho)    = D(nrho)    - temp*DU(nrho-1)
      B(nrho)    = B(nrho)    - temp* B(nrho-1)

      IF (lverbj) THEN
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)'MATRIX ELEMENT AT (nrho,nrho-2)'
         WRITE(6,*) temp
         WRITE(6,*)''
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)'MATRIX LOWER DIAGONAL (POST ROW OPERATION)'
         WRITE(6,*) DL
         WRITE(6,*) ''  
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)'MATRIX MAIN DIAGONAL  (POST ROW OPERATION)'
         WRITE(6,*) D
         WRITE(6,*) ''   
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)'MATRIX UPPER DIAGONAL (POST ROW OPERATION)'
         WRITE(6,*) DU
         WRITE(6,*) ''
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)'RIGHT HAND SIDE       (POST ROW OPERATION)'
         WRITE(6,*) B
         WRITE(6,*) ''
         WRITE(6,*)'-------------------------------------------------------------------------------'
      END IF     

      ! LAPACK general tridiagonal matrix solver using GE with partial pivoting
      ! See also
      !  https://netlib.org/lapack/explore-html/d4/d62/group__double_g_tsolve.html
      !CALL DGTSV(nrho,1,DL,D,DU,B,nrho)

      ! Thomas algorithm
      ! https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm 
      ! a_i = DL, b_i = D, c_i = DU, d_i = B
      a1 = 0; a2 = 0;
      a1(1) = DU(1)/D(1) !c_1' = c_1/b_1
      a2(1) = B(1)/D(1)  !d_1' = d_1/b_1
      DO i = 2, nrho-1
         a1(i) = DU(i)/(D(i)-DL(i)*a1(i-1))                 ! c_i' =              c_i/(b_i - a_i*c_i-1')
         a2(i) = (B(i)-DL(i)*a2(i-1))/(D(i)-DL(i)*a1(i-1))  ! d_i' = (d_i-a_i*d_i-1')/(b_i - a_i*c_i-1')
      END DO
      a2(nrho) = (B(nrho)-DL(nrho)*a2(nrho-1))/(D(nrho)-DL(nrho)*a1(nrho-1))
      B(nrho) = a2(nrho)
      DO i = nrho-1, 1, -1
         B(i) = a2(i)-a1(i)*B(i+1)
      END DO
      IF (lverbj) THEN
         WRITE(6,*)'==============================================================================='
         WRITE(6,*)'POST SOLVING EQUATIONS'
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)'U ARRAY AT CURRENT TIMESTEP'
         WRITE(6,*) B
         WRITE(6,*)''
      END IF
      ! B = mu0 I / phip => I = phip*B/mu0 = 2*phi_a*rho*B/mu0
      B = 2*eq_phiedge/mu0*(THRIFT_RHO*B)

      IF (lverbj) THEN
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)'I_TOTAL AT CURRENT TIMESTEP'
         WRITE(6,*) B
         WRITE(6,*)''
      END IF    

      ! Calculate I_source (a1)
      a1 = 0
      rho = THRIFT_RHO(1)
      s = rho*rho
      ier = 0
      CALL get_equil_Rmajor(s,h,h,Aminor,ier)
      a1(1) = THRIFT_JSOURCE(1,mytimestep)*pi*Aminor**2 
      DO i = 2, nrho
         rho = THRIFT_RHO(i)
         s = rho * rho
         ier = 0
         CALL get_equil_Rmajor(s,h,h,Aminor,ier)  ! aminor_i
         temp = THRIFT_RHO(i-1)
         s = temp*temp
         CALL get_equil_Rmajor(s,h,h,temp,ier) ! temp = aminor_i-1
         a1(i) = a1(i-1) + THRIFT_JSOURCE(i,mytimestep)*pi*(Aminor**2-temp**2)
      END DO

      ! I_plasma = I_total - I_source
      B = B - a1 
      IF (lverbj) THEN
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)'I_SOURCE AT CURRENT TIMESTEP'
         WRITE(6,*) a1
         WRITE(6,*)''
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)'I_PLASMA AT CURRENT TIMESTEP'
         WRITE(6,*) B
         WRITE(6,*)''
      END IF    

      ! J_plasma = dI_plasma/dA = dI/ds ds/dA =1/2rho dI/drho 2pi R/V' 
      !          = pi R/(rho v') dI/drho = pi R / (rho phi_a dV/dphi) dI_plasma/drho
      DO i = 1, nrho
         ! Calculate derivative, store in temp
         IF (i == 1) THEN ! Symmetry BC
            temp = (B(2)-B(1))/(2*h)
         ELSE IF (i == nrho) THEN ! 
            temp = B(nrho) ! Set current at nrho = nrho-1 for now
            temp = (4*temp-3*B(nrho)-B(nrho-1))/(3*h) 
         ELSE
            temp = (B(i+1)-B(i-1))/(2*h)
         END IF
         rho = THRIFT_RHO(i)
         s = rho*rho
         ier = 0
         CALL EZspline_interp(vp_spl,rho,vp,ier)
         CALL get_equil_Rmajor(s,Rmajor,s11,Aminor,ier)
         THRIFT_JPLASMA(i,mytimestep) = pi*Rmajor/(eq_phiedge*rho*vp)*temp
      END DO
      IF (lverbj) THEN
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)'J_PLASMA AT CURRENT TIMESTEP'
         WRITE(6,*) THRIFT_JPLASMA(:,mytimestep)
         WRITE(6,*)''
         WRITE(6,*)'-------------------------------------------------------------------------------'
         WRITE(6,*)'JINDUCTIVE DONE'
         WRITE(6,*)'==============================================================================='
      END IF    
      DEALLOCATE(A_temp, B_temp, C_temp, D_temp, &
               B_der, C_der, D_der, &
               a1, a2, a3, a4, &
               DU, D, DL, B)
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_jinductive

       

!---------------------------------------------
!     OLD NON-FUNCTIONING CODE
!---------------------------------------------------------
      ! Temp rho grid extended beyond boundaries (nrho+4)
      ! To get on nrho+2 grid, shift index by +1
      ! To get on nrho grid, shift index by +2
      !rho_temp(1) = -THRIFT_RHO(1)
      !rho_temp(2) = 0
      !rho_temp(3:nrho+2) = THRIFT_RHO
      !rho_temp(nrho+3) = 1
      !rho_temp(nrho+4) = 2*rho_temp(nrho+3) - rho_temp(nrho+2)

      ! MU0 I / PHIP   (nrho+4)
      !WRITE(6,*) "  rho    s11    s12     iota     mu0 I/phip"
      !WRITE(6,*) "------------------------------------------"

      !I_temp(1) = 0 ! rho = -rho(1), will be overwritten
      !I_temp(2) = 0 ! rho = 0
      !WRITE(6,'(2X,F5.3,4(3X,E9.2,A)') rho_temp(2),I_temp(2),';'
      !DO i = 1, nrho+1
      !      rho = rho_temp(i+2)
      !      s   = rho*rho
      !      ier = 0
      !      CALL get_equil_sus(s,s11,s12,Bav,Bsqav,ier) 
      !      CALL EZspline_interp(iota_spl,rho,iota,ier)   
      !      I_temp(i+2) = 2*rho*(s11*iota+s12) ! mu0 I/phip|(rho)
      !      WRITE(6,'(2X,F5.3,4(3X,E9.2),A)') rho_temp(i+2),s11,s12,iota,I_temp(i+2),';'
      !END DO
      !I_temp(nrho+4)=I_temp(nrho+3) ! no current outside rho=1
      !I_temp(1) = I_temp(3) ! rho = -rho(1)
      !WRITE(6,'(2X,F5.3,3X,E8.1,A)') rho_temp(nrho+4),I_temp(nrho+4),';'
      !WRITE(6,'(2X,F5.3,3X,E8.1,A)') rho_temp(1),I_temp(1),';'

      ! D/DS(U)   (nrho+2)
      !WRITE(6,*) ""
      !WRITE(6,*) "  rho   1/(2rho)*d/drho(mu0 I/phip)"
      !WRITE(6,*) "-----------------------------------"

      !f2(1) = 0 ! rho = 0
      !WRITE(6,'(2X,F5.3,6X,E8.1,A)') rho_temp(2),f2(1),';'
      !DO i = 2, nrho+2
      !   f2(i) = 0.5/rho_temp(i+1) * & ! 1/(2 rho)
      !            0.5*(I_temp(i+2)-I_temp(i))/(rho_temp(i+2)-rho_temp(i)) ! du/drho
      !   WRITE(6,'(2X,F5.3,6X,E8.1,A)') rho_temp(i+1),f2(i),';'
      !END DO

      ! Temporary Jsource grid with boundaries (nrho+2)
      !j_temp(1)          = THRIFT_JSOURCE(1,mytimestep)
      !j_temp(2:nrho+1)   = THRIFT_JSOURCE(:,mytimestep)
      !j_temp(nrho+2)     = 0.0

      !WRITE(6,*) ""
      !WRITE(6,*) "  rho      f"
      !WRITE(6,*) "---------------------"

      ! F in DF/DX     (nrho+2)
      !DO i = 1, nrho+2
      !   rho = rho_temp(i+1)
      !   s   = rho*rho
      !   ier = 0
      !   CALL get_prof_etapara(rho,THRIFT_T(mytimestep),etapara) !eta = eta
      !   CALL EZspline_interp(vp_spl,rho,s11,ier)           !s11 = V'
      !   CALL get_prof_pprime(rho,THRIFT_T(mytimestep),s12) !s12 = p'
      !   CALL get_equil_Bav(s,Bav,Bsqav,ier)                !Bav = <B>; Bsqav = <B^2>

         ! I on nrho+4 grid, f on nrho+2, so I_temp is shifted by +1
      !   f2(i) = etapara*s11*(s12*I_temp(i+1)  &   ! eta*V'(p'u  
      !                  + Bsqav/mu0*f2(i)       &   !     + <B^2>/mu0*du/dx 
      !                  - j_temp(i)*Bav)            !         - <J.B>)
      !   WRITE(6,'(F5.3,3X,E8.1,A)') rho, f2(i),';'

      !END DO

      ! I_temp (nrho+4), j_temp(nrho+2) freed up
      !I_temp = 0; j_temp = 0;
      
      ! Construct du/dt (nrho)
      !DO i = 1, nrho
      !   rho = THRIFT_RHO(i)
      !   s   = rho*rho
      !   ier = 0
      !   CALL get_equil_sus(s,s11,s12,Bav,Bsqav,ier) 
      !   j_temp(i) = 0.5/rho * & ! 1/(2 rho)
      !            0.5*(f2(i+2)-f2(i))/(rho_temp(i+3)-rho_temp(i+1)) ! df/drho
      !   j_temp(i) = 0.5*phip/rho * & 
      !            mu0*s11/(eq_phiedge**2)*j_temp(i) ! j_temp= dI/dt, assuming d/dt(phip)=0
      !END DO

      ! Need dj/dt though, not dI/dt
      ! dj/dt = pi R/(rho V') d/drho(dI/dt)
      ! Once again, need to evaluate d/drho using CD scheme (nrho+2)
      ! BCs: dI/dt(0) = 0, dI/dt(rho) = dI/dt(prev rho)
      !I_temp(1) = 0 
      !I_temp(2:nrho) = j_temp
      !I_temp(nrho+1) = I_temp(nrho)
      ! Note that I_temp(nrho+3)=I_temp(nrho+4)=0 are never accessed
      
      ! j_temp freed up again
      !j_temp = 0

      ! Calculate dj/dt (nrho)
      !DO i = 1, nrho
      !   rho = THRIFT_RHO(i)
      !   CALL EZspline_interp(vp_spl,rho,phip,ier)
      !   j_temp(i) = 0.5/rho * & ! 1/(2 rho)
      !               0.5*(I_temp(i+2)-I_temp(i))/(rho_temp(i+3)-rho_temp(i+1)) ! d/rho(dI/dt)
      !   j_temp(i) = j_temp(i)*pi*eq_Rmajor/(rho*phip) ! dj/dt
      !END DO     

      ! update plasma current
      ! s11 = dt
      !IF (mytimestep==1) THEN 
      !   s11 = THRIFT_T(2)-THRIFT_T(1)
      !   THRIFT_JPLASMA(:,mytimestep) = j_temp*s11
      !ELSE 
      !   s11 = THRIFT_T(mytimestep)-THRIFT_T(mytimestep-1)
      !   THRIFT_JPLASMA(:,mytimestep) = THRIFT_JPLASMA(:,mytimestep-1)+j_temp*s11
      !END IF

      ! Deallocate
