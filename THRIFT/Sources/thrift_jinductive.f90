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
      USE EZspline
      USE EZspline_obj
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: i,j,prevtimestep,ier
      INTEGER :: bcs0(2)
      REAL(rprec) :: rho,s,drho,ds,dt,mytime,s11,jsource,temp  
      TYPE(EZspline1_r8) :: splinor
      REAL(rprec), DIMENSION(:), ALLOCATABLE ::j_temp,&
                     A_temp,B_temp,C_temp,D_temp,&
                     BP_temp, CP_temp, DP_temp,&
                     alpha1,alpha2,alpha3,alpha4,&
                     DIAGSUB,DIAGMID,DIAGSUP,RHS
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
!     CALCULATE MAGNETIC VARIABLES
!----------------------------------------------------------------------

      THRIFT_PHIEDGE(1,mytimestep) = eq_phiedge     
      DO i = 1, nssize
         s = THRIFT_S(i)
         ier = 0
         CALL get_equil_Rmajor(s, THRIFT_RMAJOR(i,mytimestep), temp, THRIFT_AMINOR(i,mytimestep), ier)
         CALL get_equil_sus(s, THRIFT_S11(i,mytimestep), temp, temp, temp, ier)
         CALL get_equil_Bav(s, THRIFT_BAV(i,mytimestep), THRIFT_BSQAV(i,mytimestep), ier)
         ! V' = dV/ds = 2*pi*R*(pi*a^2)
         THRIFT_VP(i,mytimestep) = 2*pi**2*THRIFT_RMAJOR(i,mytimestep)*THRIFT_AMINOR(i,mytimestep)**2
      END DO
      THRIFT_S11 = ABS(THRIFT_S11)
      IF (lverbj) CALL print_calc_magvars()

!----------------------------------------------------------------------
!     UPDATE TRACKER VARIABLES
!----------------------------------------------------------------------

      CALL curden_to_curtot(THRIFT_JBOOT(:,  mytimestep),THRIFT_IBOOT(:,  mytimestep))
      CALL curden_to_curtot(THRIFT_JECCD(:,  mytimestep),THRIFT_IECCD(:,  mytimestep))
      CALL curden_to_curtot(THRIFT_JNBCD(:,  mytimestep),THRIFT_INBCD(:,  mytimestep))
      CALL curden_to_curtot(THRIFT_JOHMIC(:, mytimestep),THRIFT_IOHMIC(:, mytimestep))
      CALL curden_to_curtot(THRIFT_JPLASMA(:,mytimestep),THRIFT_IPLASMA(:,mytimestep))

      THRIFT_ISOURCE(:,mytimestep)  = THRIFT_IBOOT(:,  mytimestep)&
                                    + THRIFT_IECCD(:,  mytimestep)&
                                    + THRIFT_INBCD(:,  mytimestep)&
                                    + THRIFT_IOHMIC(:, mytimestep)

!----------------------------------------------------------------------
!     TIME AND BETA=0
!----------------------------------------------------------------------

      ! SPECIAL CASE >> If mytimestep = 1 ITOT=0 and continue to next iteration
      IF (mytimestep==1) THEN
            THRIFT_JPLASMA(:,mytimestep) = -THRIFT_JSOURCE(:,mytimestep)
            THRIFT_IPLASMA(:,mytimestep) = -THRIFT_ISOURCE(:,mytimestep)
            THRIFT_I(:,mytimestep) = THRIFT_IPLASMA(:,mytimestep)+THRIFT_ISOURCE(:,mytimestep)
            RETURN 
      END IF

      ! Time variables
      mytime = THRIFT_T(mytimestep) ! mytime = current sim time
      prevtimestep = mytimestep-1   ! previous time step index
      dt = THRIFT_T(mytimestep)-THRIFT_T(prevtimestep) ! dt = delta t this iter
      IF (mytime>tmax.and.THRIFT_T(prevtimestep)<=tmax.and.(nsubsteps==1)) WRITE(6,*) &
         '! THRIFT has exceeded end time of profiles file. Proceeding with profiles at t=tmax !' 
 
      ! SPECIAL CASE >> If at zero beta, copy previous value of JPLASMA and skip
      IF (eq_beta == 0) THEN
         IF (mytimestep /= 1) THEN
            THRIFT_JPLASMA(:,mytimestep) = THRIFT_JPLASMA(:,prevtimestep)
            THRIFT_IPLASMA(:,mytimestep) = THRIFT_IPLASMA(:,prevtimestep)
            THRIFT_I(:,mytimestep)       = THRIFT_I(:,prevtimestep)
         RETURN
         END IF
      END IF
  
!----------------------------------------------------------------------
!     CALCULATING COEFFICIENTS ABCD
!----------------------------------------------------------------------
!     Calculate ABCD (everything evaluated at rho_j)
!     > A(j) = S11/Phi_a^2
!     > B(j) = etapara*V'*<B^2>/mu_0
!     > C(j) = etapara*V'*p'
!     > D(j) = -etapara*V'*<Js.B>
!----------------------------------------------------------------------

      ALLOCATE(j_temp(nrho+2),&
               A_temp(nssize),B_temp(nssize),C_temp(nssize),D_temp(nssize),&
               BP_temp(nssize),CP_temp(nssize),DP_temp(nssize))

      ! Set up J spline
      CALL extrapolate_arr(THRIFT_JSOURCE(:,mytimestep),j_temp)
      bcs0=(/ 0, 0/)
      CALL EZspline_init(splinor,nrho+2,bcs0,ier)
      splinor%x1        = THRIFT_RHOFULL
      splinor%isHermite = 1
      CALL EZspline_setup(splinor,j_temp,ier,EXACT_DIM=.true.)
      DEALLOCATE(j_temp) 
      ! J spline set up

      ! Make the calculations
      ALLOCATE(j_temp(nssize))
      DO i = 1, nssize
         s = THRIFT_S(i)
         rho = SQRT(s)
         ! A
         A_temp(i) = THRIFT_S11(i,mytimestep)/THRIFT_PHIEDGE(1,mytimestep)**2
         ! Helper variable
         CALL get_prof_etapara(MIN(rho,SQRT(THRIFT_S(nssize-1))),mytime,THRIFT_ETAPARA(i,mytimestep))
         temp = THRIFT_ETAPARA(i,mytimestep)*THRIFT_VP(i,mytimestep)
         ! B
         B_temp(i) =  temp*THRIFT_BSQAV(i,mytimestep)/mu0
         ! C
         CALL get_prof_pprime(rho, mytime, THRIFT_PPRIME(i,mytimestep))
         C_temp(i) =  temp*THRIFT_PPRIME(i,mytimestep)
         ! D
         CALL EZspline_interp(splinor, rho, j_temp(i), ier)
         jsource = j_temp(i)
         D_temp(i) = -temp*jsource*THRIFT_BAV(i,mytimestep)
      END DO

      IF (lverbj) CALL print_calc_abcd(j_temp)
      DEALLOCATE(j_temp)
      CALL EZspline_free(splinor,ier)

!     Calculate derivatives of ABCD
      ds = THRIFT_S(2)-THRIFT_S(1)
      DO i = 2, nssize-1
         BP_temp(i) = (B_temp(i+1)-B_temp(i-1))/(2*ds)
         CP_temp(i) = (C_temp(i+1)-C_temp(i-1))/(2*ds)
         DP_temp(i) = (C_temp(i+1)-D_temp(i-1))/(2*ds)
      END DO

      ! Store values
      THRIFT_COEFF_A(:,mytimestep) = A_temp 
      THRIFT_COEFF_B(:,mytimestep) = B_temp
      THRIFT_COEFF_C(:,mytimestep) = C_temp
      THRIFT_COEFF_D(:,mytimestep) = D_temp
      THRIFT_COEFF_BP(:,mytimestep) = BP_temp
      THRIFT_COEFF_CP(:,mytimestep) = CP_temp
      THRIFT_COEFF_DP(:,mytimestep) = DP_temp

      IF (lverbj) CALL print_abcd()

!----------------------------------------------------------------------
!     CALCULATING COEFFICIENTS ALPHA
!----------------------------------------------------------------------
!     Calculate alpha_1234 (everything evaluated at s_j)
!     > a1(j) = A*D'
!     > a2(j) = A*C'
!     > a3(j) = A*(B'+ C)
!     > a4(j) = A*B
!----------------------------------------------------------------------

      ALLOCATE(alpha1(nssize-2),alpha2(nssize-2),alpha3(nssize-2),alpha4(nssize-2))

      DO i = 1, nssize-2
         j = i + 1 ! ABCD in [0,1], alphas in (0,1), so index shifts
         alpha1(i) = A_temp(j)* DP_temp(j)
         alpha2(i) = A_temp(j)* CP_temp(j)
         alpha3(i) = A_temp(j)*(BP_temp(j)+C_temp(j))            
         alpha4(i) = A_temp(j)*  B_temp(j)  
      END DO

      THRIFT_ALPHA1(:,mytimestep) = alpha1
      THRIFT_ALPHA2(:,mytimestep) = alpha2
      THRIFT_ALPHA3(:,mytimestep) = alpha3
      THRIFT_ALPHA4(:,mytimestep) = alpha4

      DEALLOCATE(A_temp,B_temp,C_temp,D_temp,BP_temp,CP_temp,DP_temp)

      IF (lverbj) CALL print_alpha()

!----------------------------------------------------------------------
!     SOLVING SYSTEM OF EQUATIONS
!----------------------------------------------------------------------
!     Populate matrix and RHS of the system of equations (N=nssize)
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
!
!----------------------------------------------------------------------

      ALLOCATE(DIAGSUB(nssize-1),DIAGMID(nssize),DIAGSUP(nssize-1),RHS(nssize))

      ! Magnetic axis (s=0)
     !DIAGSUB(1) in row 2
      DIAGMID(1) = 1
      DIAGSUP(1) = 0
      RHS(1)     = 0

      ! On the (0,1) grid
      DIAGSUB(1:nssize-2) = -alpha3/(2*ds) + alpha4/(ds**2)
      DIAGMID(2:nssize-1) =  alpha2       -2*alpha4/(ds**2) - 1.0/dt
      DIAGSUP(2:nssize-1) =  alpha3/(2*ds) + alpha4/(ds**2)
      RHS(2:nssize-1)     = -alpha1-THRIFT_UGRID(2:nssize-1,prevtimestep)/dt

      ! Plasma edge (s=1)
      temp = THRIFT_BSQAV(nssize,mytimestep)/mu0
      DIAGSUB(nssize-1) = -1.0/(2*ds)
      DIAGMID(nssize)   =  1.0/(2*ds)+THRIFT_PPRIME(nssize,mytimestep)/temp
      ! jsource = jsource(nssize) still
      RHS(nssize)       = jsource*THRIFT_BAV(nssize,mytimestep)/temp

      ! Store values
      THRIFT_MATMD(:,mytimestep)  = DIAGMID
      THRIFT_MATUD(:,mytimestep)  = DIAGSUP
      THRIFT_MATLD(:,mytimestep)  = DIAGSUB
      THRIFT_MATRHS(:,mytimestep) = RHS

      if (lverbj) CALL print_syseqs()

      ! Solve system of equations
      CALL DGTSV(nssize, 1, DIAGSUB, DIAGMID, DIAGSUP, RHS, nssize, ier)
      THRIFT_UGRID(:,mytimestep) = RHS

      IF (lverbj) CALL check_sol(THRIFT_MATLD(:,mytimestep),THRIFT_MATMD(:,mytimestep),&
                                 THRIFT_MATUD(:,mytimestep),THRIFT_MATRHS(:,mytimestep),THRIFT_UGRID(:,mytimestep))

      DEALLOCATE(alpha1,  alpha2,  alpha3,  alpha4,&
                 DIAGSUB, DIAGMID, DIAGSUP, RHS)

!----------------------------------------------------------------------
!     POST SOLVING EQUATIONS
!----------------------------------------------------------------------
!     Solving system of equations yields u(s_j) at this timestep.
!     ITOTAL(s_j) = Phi_a/mu0 * u(s_j)
!     Obtain ISOURCE from JSOURCE with curden_to_curtot subroutine.
!     Plasma current: IPLASMA = ITOTAL - ISOURCE
!     Obtain JPLASMA from IPLASMA with curtot_to_curden subroutine.
!----------------------------------------------------------------------
      ALLOCATE(j_temp(nrho))

      THRIFT_I(:,mytimestep) = THRIFT_PHIEDGE(1,mytimestep)/mu0*THRIFT_UGRID(:,mytimestep)
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