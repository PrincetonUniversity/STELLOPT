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
      REAL(rprec) :: rho,s,ds,dt,mytime,temp,Lext,rmaj,amin
      REAL(rprec), DIMENSION(:), ALLOCATABLE ::j_temp,&
                     A_temp,B_temp,C_temp,D_temp,&
                     BP_temp, CP_temp, DP_temp,temp_arr,  &
                     alpha1,alpha2,alpha3,alpha4,&
                     DIAGSUB,DIAGMID,DIAGSUP,RHS
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
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
      mytime = THRIFT_T(mytimestep)
 
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
      
      ! Allocations
      ALLOCATE(A_temp(nsj),B_temp(nsj),C_temp(nsj),D_temp(nsj),&
               BP_temp(nsj),CP_temp(nsj),DP_temp(nsj), temp_arr(nsj))
      A_temp  = 0; B_temp  = 0; C_temp  = 0; D_temp = 0
      BP_temp = 0; CP_temp = 0; DP_temp = 0; temp_arr= 0

      ! Coefficients
      temp_arr= THRIFT_ETAPARA(:,mytimestep)*THRIFT_VP(:,mytimestep)
      A_temp = THRIFT_S11(:,mytimestep)/THRIFT_PHIEDGE(mytimestep)**2
      B_temp = temp_arr*THRIFT_BSQAV(:,mytimestep)/mu0
      C_temp = temp_arr*THRIFT_PPRIME(:,mytimestep)
      D_temp = -temp_arr*THRIFT_JSOURCE(:,mytimestep)*THRIFT_BAV(:,mytimestep)

      DEALLOCATE(temp_arr)

      ! Derivatives
      ds = THRIFT_S(2) - THRIFT_S(1)
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
!     Populate matrix and RHS of the system of equations (n=nsj, nX=n-X)
!
!     | b1   c1    0   0 .  0   0   0   0  | | u1  |   | y1  |  BC AXIS
!     | a2   b2   c2   0 .  0   0   0   0  | | u2  |   | y2  |     \
!     |  0   a3   b3  c3 .  0   0   0   0  | | u3  |   | y3  |     |
!     |  .    .    .    ... .   .   .   .  | | .   | = |  .  |  evol. eq.
!     |  0    0    0   0 . an2 bn2 cn2  0  | | un2 |   | yn2 |     |
!     |  0    0    0   0 .  0  an1 bn1 cn1 | | un1 |   | yn1 |     /
!     |  0    0    0   0 .  0  X2  an  cn  | | un  |   | yn  |  BC EDGE
!                    
!     The {u} contain the current: u_j = mu0/phi_edge*I_j
!     First equation encapsulates the BC for the magnetic axis:
!        I(s=0,t) = 0 
!
!     Last equation encapsulates the BC for the plasma edge:
!        E(s=1,t) = -L_ext/(2*pi*R0) * dI/dt 
!     (We have to manipulate the last row to get a TDM)
!
!     Remaining equations are the evolution equation on s in (0,1)
!        du/dt = a1 + a2*u + a3*u' + a4*u"
!======================================================================
      ALLOCATE(DIAGSUB(nsj-1),DIAGMID(nsj),DIAGSUP(nsj-1),RHS(nsj))
      DIAGSUB = 0; DIAGMID = 0; DIAGSUP = 0; RHS = 0;

      ! On the (0,1) grid
      DIAGSUB(1:nsj-2) = -alpha3/(2*ds) + alpha4/(ds**2)
      DIAGMID(2:nsj-1) =  alpha2       -2*alpha4/(ds**2) - 1.0/dt
      DIAGSUP(2:nsj-1) =  alpha3/(2*ds) + alpha4/(ds**2)
      RHS(2:nsj-1)     = -alpha1-THRIFT_UGRID(2:nsj-1,prevtimestep)/dt
      ! Magnetic axis (s=0)
      DIAGMID(1)  = 1
      DIAGSUP(1)  = 0
      RHS(1)      = 0
      ! code for dI/ds = 0
      !DIAGMID(1) = 3
      !DIAGSUP(1) = -4
      !RHS(1)     = 0

      ! Plasma edge (s=1)
      rmaj = THRIFT_RMAJOR(nsj,mytimestep); amin = THRIFT_AMINOR(nsj,mytimestep) ! R,a helpers
      Lext = mu0*rmaj*(log(8*rmaj/amin)-2) ! mu0 R (log(8R/a)-2)
      temp = 2*pi*rmaj*mu0/THRIFT_PHIEDGE(mytimestep)*THRIFT_ETAPARA(nsj,mytimestep)/Lext*THRIFT_JSOURCE(nsj,mytimestep)
      RHS(nsj) = THRIFT_UGRID(nsj,prevtimestep)/dt + temp ! u/dt + 2*pi*R*(mu0/phi_edge)*(eta/Lext)*Js
      temp = rmaj/(amin**2*ds)*THRIFT_ETAPARA(nsj,mytimestep)/Lext ! X2 = R0/(a^2 ds)*eta/Lext
      DIAGSUB(nsj-1) = -4*temp
      DIAGMID( nsj ) = 3*temp+1.0/dt

      ! Row manipulations to get TDM for DGTSV
      ! Eliminate X2
      temp = temp/DIAGSUB(nsj-2) ! X2/an1
      DIAGSUB(nsj-1)    = DIAGSUB(nsj-1)  - temp*DIAGMID(nsj-1)
      DIAGMID(nsj)      = DIAGMID(nsj)    - temp*DIAGSUP(nsj-1)
      RHS(nsj)          = RHS(nsj)        - temp*RHS(nsj-1)
      
      ! code for dI/ds = 0
      !! Eliminate X1
      !temp = 1.0/DIAGSUP(2) ! X1/c2 [X1=1]
      !DIAGMID(1)        = DIAGMID(1)      - temp*DIAGSUB(1)
      !DIAGSUP(1)        = DIAGSUP(1)      - temp*DIAGMID(2)
      !RHS(1)            = RHS(1)          - temp*RHS(2)

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
      ! We cannot use THRIFT_J directly, otherwise we cannot do proper
      ! picard iterations.

      ALLOCATE(j_temp(nsj))
      j_temp = 0 ! Will store total current density

      THRIFT_I(:,mytimestep) = (THRIFT_PHIEDGE(mytimestep)/mu0)*THRIFT_UGRID(:,mytimestep)
      CALL curtot_to_curden(THRIFT_I(:,mytimestep),j_temp)
      THRIFT_JPLASMA(:,mytimestep) = j_temp - THRIFT_JSOURCE(:,mytimestep)

      IF (lverbj) CALL print_postevolve(j_temp)
      DEALLOCATE(j_temp)
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_jinductive