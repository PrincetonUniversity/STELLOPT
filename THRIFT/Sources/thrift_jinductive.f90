!-----------------------------------------------------------------------
!     Subroutine:    thrift_jinductive
!     Authors:       L. van Ham
!     Date:          16/03/2023
!     Description:   This subroutine calculates the inductive component
!                    of the total current. 
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
      USE EZspline
      USE EZspline_obj
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: i, j, prevtimestep, ier
      INTEGER :: bcs0(2)
      REAL(rprec) :: rho,s,drho,ds,dt,mytime,js_edge,p2,p1,temp
      TYPE(EZspline1_r8) :: splinor
      REAL(rprec), DIMENSION(:), ALLOCATABLE ::j_temp,&
                     A_temp,B_temp,C_temp,D_temp,&
                     BP_temp, CP_temp, DP_temp,temp_arr,  &
                     alpha1,alpha2,alpha3,alpha4,&
                     DIAGSUB,DIAGMID,DIAGSUP,RHS
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!======================================================================
!     CALCULATE MAGNETIC VARIABLES AND OTHER PROFILE DATA
!     > Etapara, p', Js, Phiedge, Rmajor, Aminor, S11, Bav, Bsqav, VP
!======================================================================
      ! Set up Jsource spline
      ALLOCATE(j_temp(nrho+2))
      j_temp = 0
      CALL extrapolate_arr(THRIFT_JSOURCE(:,mytimestep),j_temp)
      bcs0=(/ 0, 0/)
      CALL EZspline_init(splinor,nrho+2,bcs0,ier)
      splinor%x1        = THRIFT_RHOFULL
      splinor%isHermite = 1
      CALL EZspline_setup(splinor,j_temp,ier,EXACT_DIM=.true.)
      DEALLOCATE(j_temp) 

      ! Grab vars from profiles
      mytime = THRIFT_T(mytimestep)
      THRIFT_PHIEDGE(1,mytimestep) = eq_phiedge     
      ALLOCATE(j_temp(nsj))
      j_temp = 0
      DO i = 1, nsj
         s = THRIFT_S(i)
         rho = SQRT(s)
         ier = 0
         CALL get_equil_Rmajor(s, THRIFT_RMAJOR(i,mytimestep), temp, THRIFT_AMINOR(i,mytimestep), ier)
         CALL get_equil_sus(s, THRIFT_S11(i,mytimestep),temp,temp,temp,ier)
         CALL get_equil_Bav(s, THRIFT_BAV(i,mytimestep),THRIFT_BSQAV(i,mytimestep), ier)
         CALL EZspline_interp(vp_spl, rho, temp, ier) ! temp = dV/dPhi
         ! V' = dV/ds = dV/dPhi dPhi/ds = Phi_edge * dV/dPhi
         THRIFT_VP(i,mytimestep) = THRIFT_PHIEDGE(1,mytimestep)*temp
         ! eta breaks at rho=1(s=1) so look one gridpoint back
         CALL get_prof_etapara(MIN(rho,SQRT(THRIFT_S(nsj-1))),mytime,THRIFT_ETAPARA(i,mytimestep))
         CALL get_prof_pprime(rho, mytime, THRIFT_PPRIME(i,mytimestep))
         CALL EZspline_interp(splinor, rho, j_temp(i), ier)
      END DO
      THRIFT_S11 = ABS(THRIFT_S11)
      js_edge = j_temp(nsj)
      CALL EZspline_free(splinor,ier)  

      ! Get pprime in s-space using finite difference
      ds = THRIFT_S(2)-THRIFT_S(1)
      THRIFT_PPRIME(1,mytimestep) = 0
      DO i = 2, nsj-1
         CALL get_prof_p( SQRT(THRIFT_S(i+1)), mytime, p1)
         CALL get_prof_p( SQRT(THRIFT_S(i-1)), mytime, p2)
         THRIFT_PPRIME(i,mytimestep) = (p2-p1)/(2*ds)
      END DO
      THRIFT_PPRIME(nsj,mytimestep) = 2*THRIFT_PPRIME(nsj-1,mytimestep)-THRIFT_PPRIME(nsj-2,mytimestep)
    
      IF (lverbj) CALL print_calc_magvars()
!======================================================================
!     UPDATE TRACKER VARIABLES
!======================================================================
      CALL curden_to_curtot(THRIFT_JBOOT(:,  mytimestep),THRIFT_IBOOT(:,  mytimestep))
      CALL curden_to_curtot(THRIFT_JECCD(:,  mytimestep),THRIFT_IECCD(:,  mytimestep))
      CALL curden_to_curtot(THRIFT_JNBCD(:,  mytimestep),THRIFT_INBCD(:,  mytimestep))
      CALL curden_to_curtot(THRIFT_JOHMIC(:, mytimestep),THRIFT_IOHMIC(:, mytimestep))
      CALL curden_to_curtot(THRIFT_JPLASMA(:,mytimestep),THRIFT_IPLASMA(:,mytimestep))

      THRIFT_ISOURCE(:,mytimestep)  = THRIFT_IBOOT(:,  mytimestep) &
                                    + THRIFT_IECCD(:,  mytimestep) &
                                    + THRIFT_INBCD(:,  mytimestep) &
                                    + THRIFT_IOHMIC(:, mytimestep)
!======================================================================
!     SPECIAL CASES 
!======================================================================
      ! If mytimestep = 1 ITOT=0 and continue to next iteration
      IF (mytimestep==1) THEN
            THRIFT_JPLASMA(:,mytimestep) = -THRIFT_JSOURCE(:,mytimestep)
            THRIFT_IPLASMA(:,mytimestep) = -THRIFT_ISOURCE(:,mytimestep)
            THRIFT_I(:,mytimestep) = THRIFT_IPLASMA(:,mytimestep)+THRIFT_ISOURCE(:,mytimestep)
            RETURN 
      END IF

      ! Time variables
      prevtimestep = mytimestep-1   ! previous time step index
      dt = THRIFT_T(mytimestep)-THRIFT_T(prevtimestep) ! dt = delta t this iter
      IF (mytime>tmax.and.THRIFT_T(prevtimestep)<=tmax.and.(nsubsteps==1)) WRITE(6,*) &
         '! THRIFT has exceeded end time of profiles file. Proceeding with profiles at t=tmax !' 
 
      ! If at zero beta, copy previous value of JPLASMA and skip
      IF (eq_beta == 0) THEN
         IF (mytimestep /= 1) THEN
            THRIFT_JPLASMA(:,mytimestep) = THRIFT_JPLASMA(:,prevtimestep)
            THRIFT_IPLASMA(:,mytimestep) = THRIFT_IPLASMA(:,prevtimestep)
            THRIFT_I(:,mytimestep)       = THRIFT_I(:,prevtimestep)
         RETURN
         END IF
      END IF
!======================================================================
!     CALCULATE ABCD AND DERIVATIVES
!======================================================================
!     > A(j) = S11/Phi_a^2
!     > B(j) = etapara*V'*<B^2>/mu_0
!     > C(j) = etapara*V'*p'
!     > D(j) = -etapara*V'*<Js.B>
!======================================================================

      IF (lverbj) CALL print_calc_abcd(j_temp)
      
      ! Allocations
      ALLOCATE(A_temp(nsj),B_temp(nsj),C_temp(nsj),D_temp(nsj),&
               BP_temp(nsj),CP_temp(nsj),DP_temp(nsj), temp_arr(nsj))
      A_temp  = 0; B_temp  = 0; C_temp  = 0; D_temp = 0
      BP_temp = 0; CP_temp = 0; DP_temp = 0; temp_arr= 0

      ! Coefficients
      temp_arr= THRIFT_ETAPARA(:,mytimestep)*THRIFT_VP(:,mytimestep)
      A_temp = THRIFT_S11(:,mytimestep)/THRIFT_PHIEDGE(1,mytimestep)**2
      B_temp = temp_arr*THRIFT_BSQAV(:,mytimestep)/mu0
      C_temp = temp_arr*THRIFT_PPRIME(:,mytimestep)
      D_temp = -temp_arr*j_temp*THRIFT_BAV(:,mytimestep)

      DEALLOCATE(j_temp, temp_arr)

      ! Derivatives
      DO i = 2, nsj-1
         BP_temp(i) = (B_temp(i+1)-B_temp(i-1))/(2*ds)
         CP_temp(i) = (C_temp(i+1)-C_temp(i-1))/(2*ds)
         DP_temp(i) = (D_temp(i+1)-D_temp(i-1))/(2*ds)
      END DO
!----------------------------------------------------------------------
!     Bookkeeping
!----------------------------------------------------------------------
      THRIFT_COEFF_A(:,mytimestep) = A_temp 
      THRIFT_COEFF_B(:,mytimestep) = B_temp; THRIFT_COEFF_BP(:,mytimestep) = BP_temp
      THRIFT_COEFF_C(:,mytimestep) = C_temp; THRIFT_COEFF_CP(:,mytimestep) = CP_temp
      THRIFT_COEFF_D(:,mytimestep) = D_temp; THRIFT_COEFF_DP(:,mytimestep) = DP_temp
      IF (lverbj) CALL print_abcd()
!======================================================================
!     CALCULATE ALPHA_1,2,3,4
!======================================================================
!     > alpha1 = A*D'
!     > alpha2 = A*C'
!     > alpha3 = A*(B'+ C)
!     > alpha4 = A*B
!======================================================================
      ALLOCATE(alpha1(nsj-2),alpha2(nsj-2),alpha3(nsj-2),alpha4(nsj-2))
      alpha1 = 0; alpha2 = 0; alpha3 = 0; alpha4 = 0

      DO i = 1, nsj-2
         j = i + 1 ! ABCD in [0,1], alphas in (0,1), so index shifts
         alpha1(i) = A_temp(j)* DP_temp(j)
         alpha2(i) = A_temp(j)* CP_temp(j)
         alpha3(i) = A_temp(j)*(BP_temp(j)+C_temp(j))            
         alpha4(i) = A_temp(j)*  B_temp(j)  
      END DO
      DEALLOCATE(A_temp,B_temp,C_temp,D_temp,BP_temp,CP_temp,DP_temp)
!----------------------------------------------------------------------
!     Bookkeeping
!----------------------------------------------------------------------
      THRIFT_ALPHA1(:,mytimestep) = alpha1
      THRIFT_ALPHA2(:,mytimestep) = alpha2
      THRIFT_ALPHA3(:,mytimestep) = alpha3
      THRIFT_ALPHA4(:,mytimestep) = alpha4
      IF (lverbj) CALL print_alpha()
!======================================================================
!     SYSTEM OF EQUATIONS
!======================================================================
!     Populate matrix and RHS of the system of equations (N=nsj, NX=N-X)
!
!     | B1   C1    0   0 .  0   0   0   0  | | u1  |   | D1  |  BC AXIS
!     | A2   B2   C2   0 .  0   0   0   0  | | u2  |   | D2  |     \
!     |  0   A3   B3  C3 .  0   0   0   0  | | u3  |   | D3  |     |
!     |  .    .    .    ... .   .   .   .  | | .   | = |  .  |  evol. eq.
!     |  0    0    0   0 . AN2 BN2 CN2  0  | | uN2 |   | DN2 |     |
!     |  0    0    0   0 .  0  AN1 BN1 CN1 | | uN1 |   | DN1 |     /
!     |  0    0    0   0 .  0   0  AN  BN  | | UN  |   | DN  |  BC EDGE
!                    
!     The {u} contain the current: u(s_j) = mu0/phi_edge*I(s_j)
!     First equation encapsulates the BC for the magnetic axis:
!        I(s=0,t) = 0 
!
!     Last equation encapsulates the BC for the plasma edge:
!        f(s=1,t) = <B^2>/mu0 u' + p'u - <Js.B> = 0
!
!     Remaining equations are the evolution equation on s in (0,1)
!        du/dt = a1 + a2*u + a3*u' + a4*u"
!======================================================================
      ALLOCATE(DIAGSUB(nsj-1),DIAGMID(nsj),DIAGSUP(nsj-1),RHS(nsj))
      DIAGSUB = 0; DIAGMID = 0; DIAGSUP = 0; RHS = 0;

      ! Magnetic axis (s=0)
      DIAGMID(1) = 1
      DIAGSUP(1) = 0
      RHS(1)     = 0
      ! On the (0,1) grid
      DIAGSUB(1:nsj-2) = -alpha3/(2*ds) + alpha4/(ds**2)
      DIAGMID(2:nsj-1) =  alpha2       -2*alpha4/(ds**2) - 1.0/dt
      DIAGSUP(2:nsj-1) =  alpha3/(2*ds) + alpha4/(ds**2)
      RHS(2:nsj-1)     = -alpha1-THRIFT_UGRID(2:nsj-1,prevtimestep)/dt
      ! Plasma edge (s=1)
      temp = THRIFT_BSQAV(nsj,mytimestep)/mu0
      DIAGSUB(nsj-1) = -1.0/(2*ds)
      DIAGMID(nsj)   =  1.0/(2*ds)+THRIFT_PPRIME(nsj,mytimestep)/temp
      RHS(nsj)       = js_edge*THRIFT_BAV(nsj,mytimestep)/temp
!----------------------------------------------------------------------
!     Bookkeeping
!----------------------------------------------------------------------
      THRIFT_MATMD(:,mytimestep)  = DIAGMID
      THRIFT_MATUD(:,mytimestep)  = DIAGSUP
      THRIFT_MATLD(:,mytimestep)  = DIAGSUB
      THRIFT_MATRHS(:,mytimestep) = RHS
      if (lverbj) CALL print_syseqs()
!======================================================================
!     SOLVE SYSTEM OF EQUATIONS
!======================================================================
      CALL DGTSV(nsj, 1, DIAGSUB, DIAGMID, DIAGSUP, RHS, nsj, ier)
      ! Store solution
      THRIFT_UGRID(:,mytimestep) = RHS
      DEALLOCATE(alpha1,alpha2,alpha3,alpha4,DIAGSUB,DIAGMID,DIAGSUP,RHS)
!----------------------------------------------------------------------
!     Bookkeeping
!----------------------------------------------------------------------
      IF (lverbj) CALL check_sol(THRIFT_MATLD(:,mytimestep),THRIFT_MATMD(:,mytimestep),&
      THRIFT_MATUD(:,mytimestep),THRIFT_MATRHS(:,mytimestep),THRIFT_UGRID(:,mytimestep))
!======================================================================
!     CALCULATE PLASMA CURRENT
!======================================================================
      ALLOCATE(j_temp(nrho))
      j_temp = 0

      THRIFT_I(:,mytimestep) = (THRIFT_PHIEDGE(1,mytimestep)/mu0)*THRIFT_UGRID(:,mytimestep)
      CALL curtot_to_curden(THRIFT_I(:,mytimestep),j_temp)
      THRIFT_JPLASMA(:,mytimestep) = j_temp - THRIFT_JSOURCE(:,mytimestep)
      CALL curden_to_curtot(THRIFT_JPLASMA(:,mytimestep),THRIFT_IPLASMA(:,mytimestep))

      IF (lverbj) CALL print_postevolve(j_temp)
      DEALLOCATE(j_temp)
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_jinductive