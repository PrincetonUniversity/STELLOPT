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
      REAL(rprec) :: rho,s,h,k,t,s11,s12,iota,Bav,Bsqav,vp,etapara,pprime,temp
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: ucur,         &
                                                A_temp, B_temp, C_temp, D_temp &
                                                a1, a2, a3, a4,         &
                                                DU, D, DL, B
      TYPE(EZspline1_r8) :: f_spl
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
!     The general idea here is to use equation (22) in form to update
!     THRIFT_JPLASMA.  Where THRIFT_JSOURCE is the J in <J_s.B>.

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
      ALLOCATE(ucur(nrho), &
               A_temp(nrho+2), B_temp(nrho+2), C_temp(nrho+2), D_temp(nrho+2), &
               a1(nrho), a2(nhro), a3(nrho), a4(nrho), &
               DU(nrho-1), D(nrho), DL(nrho-1), B(nrho))

      ucur = 0;
      A_temp = 0; B_temp = 0; C_temp = 0; D_temp = 0;
      a1 = 0; a2 = 0; a3 = 0; a4 = 0;
      DU = 0; D = 0; DL = 0; B = 0;

      ! Calculate u=S11 iota + S12 this iteration
      DO i = 1, nrho
         rho = THRIFT_RHO(i) 
         s = rho*rho
         ier = 0        
         CALL get_equil_sus(s,s11,s12,h,h,ier) 
         CALL EZspline_interp(iota_spl,rho,iota)
         ucur(i) = s11*iota+s12
      END DO

      t = THRIFT_T(mytimestep)

      ! Populate A,B,C,D
      DO i = 1, nrho+2
         IF (i==1) THEN 
            rho = 0
         ELSE IF (i==nrho+2) THEN
            rho = 1
         ELSE
            rho = THRIFT_RHO(i)
         END IF
         s = rho*rho
         ier = 0
         CALL get_prof_etapara(rho, t,etapara)
         CALL EZspline_interp(vp_spl,rho,vp,ier)
         CALL get_prof_pprime(rho,t,pprime)
         CALL get_equil_Bav(s,Bav,Bsqav,ier)
         CALL get_equil_sus(s,s11,h,h,h,ier)
         IF (i==1) THEN 
            A_temp(i) = 0 ! A not necessary at axis (no derivative required)
         ELSE
            A_temp(i) = s11/(4*rho*eq_phiedge)
         END IF
         temp = 2*etapara*vp ! temp <- 2 eta dV/dPhi 
         B_temp(i) = temp*Bav/mu0 ! 2 eta dV/dPhi <B^2>/mu_0
         C_temp(i) = temp*pprime  ! 2 eta dV/dPhi dp/drho
         D_temp(i) = -temp*THRIFT_JSOURCE(i,mytimestep)*Bav ! 2 eta dV/dPhi <J.B>
      END DO

      ! Timesteps and grid spacing 
      k = THRIFT_T(2)-THRIFT_T(1)      ! k <- dt
      h = THRIFT_RHO(2)-THRIFT_RHO(1)  ! h <- drho
      ! NB: If grids are changed to be non-uniform, adjust this accordingly

      ! Visualisation of different grids
      !j=1 2    3    4    5    6    7
      ! |  |    |    |    |    |    |     ABCD_temp on j
      !    |    |    |    |    |    |     a1,2,3,4  on i
      !   i=1   2    3    4    5    6     j = i+1 

      ! a1 = A dD/drho
      ! a2 = A (1/rho dB/drho - B/rho^2 + dC/drho)
      ! a3 = A (dB/drho + B/rho + C)
      ! a4 = A B

      ! Calculate a1,2,3,4 on the gridpoints in [2,nrho-1]
      ! dY/drho(i) = [Y(j+1)-Y(j-1)]/2h = [Y(i+2)-Y(i)]/2h
      DO i = 2, nrho-1
         rho = THRIFT_RHO(i)
         a1(i) = A_temp(i+1)*(      (D_temp(i+2)-D_temp(i))/(2*h))   
         a2(i) = A_temp(i+1)*(1/rho*(B_temp(i+2)-B_temp(i))/(2*h)&   
               - B_temp(i+1)/rho**2+(C_temp(i+2)-C_temp(i))/(2*h))    
         a3(i) = A_temp(i+1)*(      (B_temp(i+2)-B_temp(i))/(2*h)&   
               + B_temp(i+1)/rho  +  C_temp(i+1))                   
         a4(i) = A_temp(i+1)*B_temp(i+1)                             
      END DO

      ! Calculate a1,2,3,4 for i=1 
      ! dY/drho(i=1) = [Y(i=3) + 3*Y(i=2) - 4*Y(i=1)]/3h
      rho = THRIFT_RHO(1) 
      a1(1) = A_temp(2)*(      (D_temp(3)+3*D_temp(2)-4*D_temp(1))/(3*h))
      a2(1) = A_temp(2)*(1/rho*(B_temp(3)+3*B_temp(2)-4*B_temp(1))/(3*h)&
            - B_temp(2)/rho**2+(C_temp(3)+3*C_temp(2)-4*C_temp(1))/(3*h))
      a3(1) = A_temp(2)*(      (B_temp(3)+3*B_temp(2)-4*B_temp(1))/(3*h)&
            + B_temp(2)/rho   + C_temp(2))
      a4(1) = A_temp(2)*B_temp(2)

      ! Calculate a1,2,3,4 for i=nrho
      ! dY/drho(i=nrho) = [4*Y(i=nrho+2) - 3*Y(i=nrho+1) - Y(i=nrho)]/3h
      rho = THRIFT_RHO(nrho)
      a1(nrho) = A_temp(nrho+1)*(      (4*D_temp(nrho+2)-3*D_temp(nrho+1)-D_temp(nrho))/(3*h))
      a2(nrho) = A_temp(nrho+1)*(1/rho*(4*B_temp(nrho+2)-3*B_temp(nrho+1)-B_temp(nrho))/(3*h)&
               - B_temp(nrho+1)/rho**2+(4*C_temp(nrho+2)-3*C_temp(nrho+1)-C_temp(nrho))/(3*h))
      a3(nrho) = A_temp(nrho+1)*(      (4*B_temp(nrho+2)-3*B_temp(nrho+1)-B_temp(nrho))/(3*h)&
            + B_temp(nrho+1)/rho   + C_temp(nrho+1))
      a4(nrho) = A_temp(nrho+1)*B_temp(nrho+1)

      ! Populate diagonals; DU and DL of size nrho-1, D of size nrho
      ! Do most of the work in one loop and fix mistakes afterwards
      DO i = 1, nrho-1
         DU(i) = a3(i)/(2*h) +a4(i)/(h**2)   ! Upper diagonal (Correct)
         DL(i) = -a3(i)/(2*h)+a4(i)/(h**2)   ! Lower diagonal (DL(nrho-1) wrong)
          D(i) = a2(i)-2*a4(i)/(h**2)-1/k    ! Middle diagonal (D(1,nrho) wrong)
          B(i) = -ucur(i)/k-a1(i)            ! Right-hand side (B(nrho) wrong)
      END DO
      
      DL(nrho-1)= -a3(nrho)/(3*h)+2*a4(nrho)/(h**2)         ! Lower diagonal fixed
        D(1)    = a2(1)-a3(1)/(2*h)-a4(1)/(h**2)-1/k  ! Middle diagonal half fixed
        D(nrho) = a2(nrho)-a3(nrho)/h+5*a4(nrho)/(h**2)-1/k! Middle diagonal fixed

      ! L_ext = mu0*R_eff( ln( 8 R_eff/r_eff ) -2 + F_shaping)      
      rho = 1
      s = rho*rho
      ier = 0
      CALL get_prof_etapara(rho,t,etapara)   
      temp = mu0*eq_Rmajor*(log(8*eq_Rmajor/(rho*eq_Aminor)) - 2 + 0.25) ! temp <- L_ext
      CALL EZspline_interp(vp_spl,rho,vp,ier)
      CALL get_equil_Bav(s,Bav,Bsqav,ier)
      CALL get_prof_pprime(rho,t,pprime) 

      edge_u(2) = edge_u(1)+k*(-mu0/(2*eq_phiedge)*etapara/temp*vp*((pprime+Bsqav/mu0)*edge_u(1))&
            - THRIFT_JSOURCE(nrho,mytimestep)*Bav)

      B(nrho) = ucur(nrho)/k-a1(nrho) - edge_u(2)*(4*a3(nrho)/(3*h)+16*a4(nrho)/(5*h**2)) ! Fix B(nrho)
      temp = -a4(nrho)/(5*h**2) ! Annoying non-zero element at (nrho,nrho-2)

      ! Eliminate that extra non-zero element to get a TDM
      temp = temp/DL(nrho-2) ! Row operation: [NRHO] -> [NRHO]-s11*[NRHO-1]
      DL(nrho-1) = DL(nrho-1) - s11* D(nrho-1) 
      D(nrho)    = D(nrho)    - s11*DU(nrho-1)
      B(nrho)    = B(nrho)    - s11* B(nrho-1)
      
      ! LAPACK general tridiagonal matrix solver using GE with partial pivoting
      ! See also
      !  https://netlib.org/lapack/explore-html/d4/d62/group__double_g_tsolve.html
      CALL DGTSV(nrho,1,DL,D,DU,B,nrho)
      ! B is the solution matrix = mu0 I / phip of the next iteration
      ! phip = 2 rho phi_a
      B = 2*eq_phiedge/mu0*(THRIFT_RHO*B)
     
      ! I is total enclosed current, I_plasma = I - I_source
      DO i = 1, nrho
         s11 = V(rho_j)-V(rho_j-1)! dV
         B(i) = B(i) - 1/(pi2*eq_Rmajor)*SUM(THRIFT_JSOURCE(1:i,mytimestep))*s11
      END DO

      ! J_plasma = dI/dA=dI/ds ds/dA=1/2rho dI/drho 2pi R/V' = pi R/(rho v') dI/drho
      DO i = 1, nrho
         IF (i == 1) THEN ! Symmetry BC
            s11 = (B(2)-B(1))/(2*h)
         ELSE IF (i == nrho) THEN ! No enclosed current at edge (for now)
            s11 = (B(nrho)-B(nrho-1))/(3*h) 
         ELSE
            s11 = (B(i+1)-B(i-1))/(2*h)
         END IF
         rho = THRIFT_RHO(i)
         ier = 0
         CALL EZspline_interp(vp_spl,rho,vp,ier)
         THRIFT_JPLASMA(i,mytimestep) = pi*eq_Rmajor/(rho*vp)*s11
      END DO

      DEALLOCATE(ucur, &
               A_temp, B_temp, C_temp, D_temp &
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
