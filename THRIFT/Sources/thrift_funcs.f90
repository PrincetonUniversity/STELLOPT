!-----------------------------------------------------------------------
!     Module:        thrift_funcs
!     Authors:       L. van Ham
!     Date:          01/03/2023
!     Description:   This module contains various functions used in
!                    thrift_jinductive.
!-----------------------------------------------------------------------
MODULE thrift_funcs
    !-------------------------------------------------------------------
    !     Libraries
    !-------------------------------------------------------------------
    USE stel_kinds, ONLY: rprec
    USE stel_tools
    USE EZspline
    USE EZspline_obj
    USE thrift_runtime
    USE thrift_vars
    USE thrift_equil
    USE thrift_profiles_mod

    IMPLICIT NONE
!-----------------------------------------------------------------------
!     Input Namelists
!         NONE
!-----------------------------------------------------------------------
      
!-----------------------------------------------------------------------
!     Subroutines
!-----------------------------------------------------------------------

CONTAINS

SUBROUTINE update_vars()
    IMPLICIT NONE
    REAL(rprec) :: mytime, s, rho, temp, p_p1, p_m1, ds
    INTEGER :: i, ier

    ! Grab vars from profiles
    mytime = THRIFT_T(mytimestep)
    THRIFT_PHIEDGE(mytimestep) = eq_phiedge    
    DO i = 1, nsj
        s = THRIFT_S(i)
        rho = SQRT(s)
        ier = 0
        CALL get_equil_Rmajor(s, THRIFT_RMAJOR(i,mytimestep), temp, THRIFT_AMINOR(i,mytimestep), ier)
        CALL get_equil_sus(s, THRIFT_S11(i,mytimestep),THRIFT_S12(i,mytimestep),temp,temp,ier)
        CALL get_equil_Bav(s, THRIFT_BAV(i,mytimestep),THRIFT_BSQAV(i,mytimestep), ier)
        CALL EZspline_interp(vp_spl, rho, temp, ier) ! temp = dV/dPhi
        ! V' = dV/ds = dV/dPhi dPhi/ds = Phi_edge * dV/dPhi
        THRIFT_VP(i,mytimestep) = THRIFT_PHIEDGE(mytimestep)*temp
        ! eta breaks at rho=1(s=1) so look one gridpoint back
        CALL get_prof_etapara(MIN(rho,SQRT(THRIFT_S(nsj-1))),mytime,THRIFT_ETAPARA(i,mytimestep))
        CALL get_prof_p(rho, mytime, THRIFT_P(i,mytimestep))
     END DO
     THRIFT_S11 = ABS(THRIFT_S11)

     ! Handle NaN
     IF (ANY(ISNAN(THRIFT_S11(:,mytimestep)))) THEN
        CALL handle_err(THRIFT_NAN_ERR,'THRIFT_S11',mytimestep)
        THRIFT_S11(:,mytimestep) = 0
     END IF
     IF (ANY(ISNAN(THRIFT_S12(:,mytimestep)))) THEN
        CALL handle_err(THRIFT_NAN_ERR,'THRIFT_S12',mytimestep)
        THRIFT_S12(:,mytimestep) = 0
     END IF
     IF (ANY(ISNAN(THRIFT_BAV(:,mytimestep)))) THEN
        CALL handle_err(THRIFT_NAN_ERR,'THRIFT_BAV',mytimestep)
        THRIFT_BAV(:,mytimestep) = 0
     END IF
     IF (ANY(ISNAN(THRIFT_BSQAV(:,mytimestep))))       THEN   
        CALL handle_err(THRIFT_NAN_ERR,'THRIFT_BSQAV',mytimestep)
        THRIFT_BSQAV(:,mytimestep) = 0
     END IF
     IF (ANY(ISNAN(THRIFT_VP(:,mytimestep))))        THEN     
        CALL handle_err(THRIFT_NAN_ERR,'THRIFT_VP',mytimestep)
        THRIFT_VP(:,mytimestep) = 0
     END IF
     IF (ANY(ISNAN(THRIFT_ETAPARA(:,mytimestep))))      THEN  
        IF (mytimestep/=1) CALL handle_err(THRIFT_NAN_ERR,'THRIFT_ETAPARA',mytimestep)
        THRIFT_ETAPARA(:,mytimestep) = 0
     END IF
     IF (ANY(ISNAN(THRIFT_P(:,mytimestep))))            THEN  
        IF (mytimestep/=1) CALL handle_err(THRIFT_NAN_ERR,'THRIFT_P',mytimestep)
        THRIFT_P(:,mytimestep) = 0
     END IF

    ! Get pprime in s-space using finite difference
     ds = THRIFT_S(2)-THRIFT_S(1)
     THRIFT_PPRIME(1,mytimestep) = 0
     DO i = 2, nsj-1
        p_p1 = THRIFT_P(i+1, mytimestep)
        p_m1 = THRIFT_P(i-1, mytimestep)
        THRIFT_PPRIME(i,mytimestep) = (p_p1-p_m1)/(2*ds)
     END DO
     THRIFT_PPRIME(1,mytimestep)   = 2*THRIFT_PPRIME(1,mytimestep)-THRIFT_PPRIME(2,mytimestep)
     THRIFT_PPRIME(nsj,mytimestep) = 2*THRIFT_PPRIME(nsj-1,mytimestep)-THRIFT_PPRIME(nsj-2,mytimestep)
   
     IF (lverbj) CALL print_calc_vars()
     
END SUBROUTINE update_vars

SUBROUTINE calc_iota()
    IMPLICIT NONE
    REAL(rprec) :: phiedge, curtor, s11, s12, ds
    REAL(rprec) :: ds11_0,ds11_1,ds11_2,ds12_0,ds12_1,ds12_2,dI_2,dI_1,dI_0
    INTEGER :: i

    DO i = 2, nsj
        phiedge = THRIFT_PHIEDGE(mytimestep)
        curtor = THRIFT_I(i,mytimestep)
        s11 = THRIFT_S11(i,mytimestep)
        s12 = THRIFT_S12(i,mytimestep)
        ! Calculate iota if s11 and phiedge exist 
        IF ((s11>0).and.(phiedge>0)) THRIFT_IOTA(i,mytimestep) = mu0*curtor/(s11*phiedge)-s12/s11
    END DO
    ! S11, I, S12 are 0 at the magnetic axis; use l'Hôpital's rule (extrapolate derivative to axis)
    ds = THRIFT_S(2) - THRIFT_S(1)
    ! I
    dI_2 = (THRIFT_I(4,mytimestep)-THRIFT_I(2,mytimestep))/(2*ds)
    dI_1 = (THRIFT_I(3,mytimestep)-THRIFT_I(1,mytimestep))/(2*ds)
    dI_0 = 2*dI_1 - dI_2
    ! S11
    dS11_2 = (THRIFT_S11(4,mytimestep)-THRIFT_S11(2,mytimestep))/(2*ds)
    dS11_1 = (THRIFT_S11(3,mytimestep)-THRIFT_S11(1,mytimestep))/(2*ds)
    dS11_0 = 2*dS11_1 - dS11_2    
    ! S12
    dS12_2 = (THRIFT_S12(4,mytimestep)-THRIFT_S12(2,mytimestep))/(2*ds)
    dS12_1 = (THRIFT_S12(3,mytimestep)-THRIFT_S12(1,mytimestep))/(2*ds)
    dS12_0 = 2*dS12_1 - dS12_2
    ! Iota
    THRIFT_IOTA(1,mytimestep) = (mu0*dI_0-phiedge*dS12_0)/(phiedge*dS11_0)

END SUBROUTINE calc_iota

SUBROUTINE curden_to_curtot(j_arr_in, i_arr_out)
    ! Takes a J(s) array and returns an I(s) array.
    REAL(rprec), DIMENSION(:), INTENT(in) :: j_arr_in
    REAL(rprec), DIMENSION(:), INTENT(out) :: i_arr_out
    INTEGER :: i
    REAL(rprec) :: ds

    ! Grid spacing
    ds = THRIFT_S(2)-THRIFT_S(1)
    ! No enclosed current on-axis
    i_arr_out(1) = 0
    ! Elsewhere: I(i) = I(i-1) + J(i)*dA(i) ; 
    !            dA(i) = dA/ds*ds = pi*aminor^2 * ds
    DO i = 2, nsj 
        i_arr_out(i) = i_arr_out(i-1) + j_arr_in(i)*(pi*eq_Aminor**2)*ds
    END DO
    RETURN

END SUBROUTINE curden_to_curtot

SUBROUTINE curtot_to_curden(i_arr_in, j_arr_out)
    ! Takes an I(s) array and returns a J(s) array
    REAL(rprec), DIMENSION(:), INTENT(in) :: i_arr_in
    REAL(rprec), DIMENSION(:), INTENT(out) :: j_arr_out
    INTEGER :: i
    REAL(rprec) :: ds

    ! Grid spacing
    ds = THRIFT_S(2)-THRIFT_S(1)
    ! Current density not on boundaries
    !    J = dI/dA 
    !      = dI/ds ds/dA
    !      = dI/ds / (pi*aminor^2) 
    DO i = 2, nsj-1
        j_arr_out(i) = (i_arr_in(i+1)-i_arr_in(i-1))/(2*ds)*1.0/(pi*eq_Aminor**2)
    END DO

    ! Extrapolate current density to boundaries
    j_arr_out(1)   = 2*j_arr_out(2)    -j_arr_out(3)     ! s = 0
    j_arr_out(nsj) = 2*j_arr_out(nsj-1)-j_arr_out(nsj-2) ! s = 1

    RETURN

END SUBROUTINE curtot_to_curden

!!===============================================================================
!!  PRINTER FUNCTIONS 
!!===============================================================================
SUBROUTINE print_calc_vars()
    INTEGER :: i
    WRITE(6,*)'==============================================================================='
    WRITE(6,*)' CALCULATING EVOLUTION VARIABLES'
    WRITE(6,*)'   S  VPRIME     BAV   BSQAV        S11    ETAPARA     PPRIME'
    WRITE(6,*)''
    DO i = 1, nsj
        WRITE(6,'(F5.3,3(1X,F7.3),3(1X,ES10.3))') &
            THRIFT_S(i),THRIFT_VP(i,mytimestep),THRIFT_BAV(i,mytimestep),THRIFT_BSQAV(i,mytimestep), &
            THRIFT_S11(i,mytimestep),THRIFT_ETAPARA(i,mytimestep),THRIFT_PPRIME(i,mytimestep)
    END DO
END SUBROUTINE print_calc_vars

SUBROUTINE print_abcd()
    INTEGER :: i
    WRITE(6,*)'==============================================================================='
    WRITE(6,*)' COEFFICIENTS ABCD'
    WRITE(6,*)'   S         A         B          C          D       BDER       CDER       DDER'
    WRITE(6,*)''
    DO i = 1, nsj
        WRITE(6,'(F5.3, 1X, 7(ES10.2,1X))') THRIFT_S(i),&
        THRIFT_COEFF_A(i,mytimestep),THRIFT_COEFF_B(i,mytimestep),THRIFT_COEFF_C(i,mytimestep),THRIFT_COEFF_D(i,mytimestep),&
        THRIFT_COEFF_BP(i,mytimestep),THRIFT_COEFF_CP(i,mytimestep),THRIFT_COEFF_DP(i,mytimestep)
     END DO
END SUBROUTINE print_abcd


SUBROUTINE print_alpha()
    INTEGER :: i
    WRITE(6,*)'==============================================================================='
    WRITE(6,*)' ALPHAS'
    WRITE(6,*)'  S       ALPHA 1        ALPHA 2        ALPHA 3        ALPHA 4'
    WRITE(6,*)''
    DO i = 1, nsj-2
        WRITE(6,'(F5.3, 1X, 4(ES13.5,2X))')  THRIFT_S(i+1),&
        THRIFT_ALPHA1(i,mytimestep), THRIFT_ALPHA2(i,mytimestep),&
        THRIFT_ALPHA3(i,mytimestep), THRIFT_ALPHA4(i,mytimestep)
     END DO
END SUBROUTINE print_alpha

SUBROUTINE print_syseqs()
    INTEGER :: i
    WRITE(6,*) '==============================================================================='
    WRITE(6,*)' SYSTEM OF EQUATIONS'
    WRITE(6,*)'  i         LOWER           MAIN          UPPER            RHS '
    WRITE(6,*)''
    WRITE(6,'(I4, 1X, 15X, 3(ES13.5,2X))') 1, THRIFT_MATMD(1,mytimestep), THRIFT_MATUD(1,mytimestep), THRIFT_MATRHS(1,mytimestep)
    DO i = 2, nsj-1
       WRITE(6,'(I4, 1X, 4(ES13.5,2X))') i, &
       THRIFT_MATLD(i-1,mytimestep), THRIFT_MATMD(i,mytimestep),THRIFT_MATUD(i,mytimestep), THRIFT_MATRHS(i,mytimestep)
    END DO
    i = nsj
    WRITE(6,'(I4, 1X,2(ES13.5,2X),15X,ES13.5)') i, THRIFT_MATLD(i-1,mytimestep), THRIFT_MATMD(i,mytimestep),THRIFT_MATRHS(i,mytimestep)
END SUBROUTINE print_syseqs

SUBROUTINE print_postevolve(j_arr)
    REAL(rprec), DIMENSION(:), INTENT(in) :: j_arr
    INTEGER :: i
    WRITE(6,*) '==============================================================================='
    WRITE(6,*)' POST EVOLUTION' 
    WRITE(6,*)'  i     s        ITOTAL        JTOTAL       JPLASMA        JSOURCE'
    WRITE(6,*)''
    DO i = 1, nsj
        WRITE(6,'(I4, 1X, F5.3, 3(1X,ES13.5), 1X, F5.3, 3(1X,ES13.5))') &
           i, THRIFT_S(i), THRIFT_I(i,mytimestep),j_arr(i),THRIFT_JPLASMA(i,mytimestep),THRIFT_JSOURCE(i,mytimestep)
     END DO
END SUBROUTINE print_postevolve


!!===============================================================================
!!  NO LONGER IN USE 
!!===============================================================================

SUBROUTINE solve_tdm(AI,BI,CI,DI,val) ! no longer used
    ! Thomas algorithm: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm 
    IMPLICIT NONE
    REAL(rprec), DIMENSION(:), INTENT(in) :: AI
    REAL(rprec), DIMENSION(:), INTENT(in) :: BI
    REAL(rprec), DIMENSION(:), INTENT(in) :: CI
    REAL(rprec), DIMENSION(:), INTENT(in) :: DI
    REAL(rprec), DIMENSION(:), INTENT(out) :: val
    INTEGER :: i
    REAL(rprec), DIMENSION(:), ALLOCATABLE :: c_p
    REAL(rprec), DIMENSION(:), ALLOCATABLE :: d_p
    REAL(rprec) :: denom
    ALLOCATE(c_p(nsj), d_p(nsj))
    c_p = 0; d_p = 0
    !! Forward sweep
    ! c_1' = c_1/b_1   ;  c_i' =              c_i/(b_i - a_i*c_i-1') [1<i<=n-1]
    ! d_1' = d_1/b_1   ;  d_i' = (d_i-a_i*d_i-1')/(b_i - a_i*c_i-1') [1<i<=n  ]
    c_p(1) = CI(1)/BI(1) 
    d_p(1) = DI(1)/BI(1) 
    DO i = 2, nsj
        denom = BI(i)-AI(i)*c_p(i-1)
        IF (i/=nsj) & 
          c_p(i) = CI(i)/denom
        d_p(i) = (DI(i)-AI(i)*d_p(i-1))/denom
    END DO 
    !! Back substitution
    ! x_n = d_n'
    ! x_i = d_i' - c_i'*x_i+1
    val(nsj) = d_p(nsj) 
    DO i = nsj-1, 1, -1
       val(i) = d_p(i)-c_p(i)*val(i+1) 
    END DO

    DEALLOCATE(c_p, d_p)
    RETURN

END SUBROUTINE solve_tdm

SUBROUTINE check_sol(AI,BI,CI,DI,sol) ! no longer used either
    REAL(rprec), DIMENSION(:), INTENT(in) :: AI
    REAL(rprec), DIMENSION(:), INTENT(in) :: BI
    REAL(rprec), DIMENSION(:), INTENT(in) :: CI
    REAL(rprec), DIMENSION(:), INTENT(in) :: DI
    REAL(rprec), DIMENSION(:), INTENT(in) :: sol
    REAL(rprec), DIMENSION(:), ALLOCATABLE :: residue
    INTEGER :: i

    ALLOCATE(residue(nsj))
    
    residue(1) = BI(1)*sol(1)+CI(1)*sol(2)-DI(1)
    DO i = 2, nsj-1
        residue(i) = AI(i-1)*sol(i-1)+BI(i)*sol(i)+CI(i)*sol(i+1)-DI(i)
    END DO
    residue(nsj) = AI(nsj-1)*sol(nsj)+BI(nsj)*sol(nsj)-DI(nsj)
    WRITE(6,'(A19,ES13.3)') "RESIDUE THIS ITER: ",maxval(residue)

    DEALLOCATE(residue)

    RETURN

END SUBROUTINE check_sol


SUBROUTINE extrapolate_arr(j_arr, j_extr)
    REAL(rprec), DIMENSION(:), INTENT(IN)  :: j_arr
    REAL(rprec), DIMENSION(:), INTENT(OUT) :: j_extr

    j_extr(1)        = (3*j_arr(1)-j_arr(2))/2
    j_extr(2:nrho+1) = j_arr
    j_extr(nrho+2)   = (3*j_arr(nrho)-j_arr(nrho-1))/2

    RETURN

END SUBROUTINE extrapolate_arr

SUBROUTINE Js_to_Jrho(j_s_in, j_rho_out)
    ! This subroutine takes an array of J on the THRIFT_S grid and
    ! returns an array of J on THRIFT_RHO grid using interpolation
    ! j_s_in must be of size nsj
    ! j_rho_out must be of size nrho
    REAL(rprec), DIMENSION(:), INTENT(in) :: j_s_in
    REAL(rprec), DIMENSION(:), INTENT(out) :: j_rho_out
    INTEGER :: i, ier
    INTEGER :: bcs0(2)
    TYPE(EZspline1_r8) :: j_spl
    REAL(rprec) :: rho, s

    ier = 0
    ! Setup J spline (in s space)
    bcs0=(/ 0, 0/)
    CALL EZspline_init(j_spl,nsj,bcs0,ier)
    j_spl%x1        = THRIFT_S
    j_spl%isHermite = 1
    CALL EZspline_setup(j_spl,j_s_in,ier,EXACT_DIM=.true.)  

    ! Interpolate J (in rho space)
    DO i = 1, nrho
       rho = THRIFT_RHO(i)
       s = rho*rho
       ier = 0
       CALL EZspline_interp(j_spl,s,j_rho_out(i),ier)      
    END DO
    CALL EZspline_free(j_spl,ier)

END SUBROUTINE Js_to_Jrho
END MODULE thrift_funcs


