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
    IMPLICIT NONE
!-----------------------------------------------------------------------
!     Input Namelists
!         NONE
!-----------------------------------------------------------------------
      
!-----------------------------------------------------------------------
!     Subroutines
!-----------------------------------------------------------------------

CONTAINS
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

SUBROUTINE check_sol(AI,BI,CI,DI,sol)
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
    WRITE(6,*) maxval(residue)

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


SUBROUTINE curden_to_curtot(j_arr_in, i_arr)
    ! Takes a J(rho) array and returns an I(s) array.
    REAL(rprec), DIMENSION(:), INTENT(in) :: j_arr_in
    REAL(rprec), DIMENSION(:), INTENT(out) :: i_arr
    REAL(rprec), DIMENSION(:), ALLOCATABLE :: j_full
    INTEGER :: i, ier
    INTEGER :: bcs0(2)
    REAL(rprec) :: s,rho,ds,j_interp
    TYPE(EZspline1_r8) :: j_spl

    ALLOCATE(j_full(nrho+2))

    ds = THRIFT_S(2)-THRIFT_S(1)
    CALL extrapolate_arr(j_arr_in,j_full)
    ! Create J spline (in rho space)
    bcs0=(/ 0, 0/)
    CALL EZspline_init(j_spl,nrho+2,bcs0,ier)
    j_spl%x1        = THRIFT_RHOFULL
    j_spl%isHermite = 1
    CALL EZspline_setup(j_spl,j_full,ier,EXACT_DIM=.true.)
    
    ! Calculate I (in s space)
    i_arr(1) = 0
    DO i = 2, nsj
        s = THRIFT_S(i)
        rho = SQRT(s)
        ier = 0
        CALL EZspline_interp(j_spl,rho,j_interp,ier)
        i_arr(i) = i_arr(i-1) + j_interp*(pi*THRIFT_AMINOR(nsj,mytimestep)**2)*ds
    END DO
    CALL EZspline_free(j_spl,ier)
    DEALLOCATE(j_full)

    RETURN

END SUBROUTINE curden_to_curtot



SUBROUTINE curtot_to_curden(i_arr, j_arr)
    ! Takes an I(s) array and returns a J(rho) array
    REAL(rprec), DIMENSION(:), INTENT(in) :: i_arr
    REAL(rprec), DIMENSION(:), INTENT(out) :: j_arr
    REAL(rprec), DIMENSION(:), ALLOCATABLE :: j_temp
    INTEGER :: i, ier
    INTEGER :: bcs0(2)
    TYPE(EZspline1_r8) :: j_spl
    REAL(rprec) :: ds, rho, s

    ALLOCATE(j_temp(nsj))

    ! Calculate J (in s space)
    ds = THRIFT_S(2)-THRIFT_S(1)
    DO i = 2, nsj-1
        j_temp(i) = (i_arr(i+1)-i_arr(i-1))/(2*ds)*1.0/(pi*THRIFT_AMINOR(nsj,mytimestep)**2)
    END DO

    ! Extrapolate to boundaries
    j_temp(1)      = 2*j_temp(2)       -j_temp(3)        ! s = 0
    j_temp(nsj) = 2*j_temp(nsj-1)-j_temp(nsj-2) ! s = 1

    ! Convert to J rho
    CALL Js_to_Jrho(j_temp, j_arr)
    DEALLOCATE(j_temp)
    RETURN

END SUBROUTINE curtot_to_curden

SUBROUTINE Js_to_Jrho(j_s_in, j_rho_out)
    ! This subroutine takes an array of J on THRIFT_S and
    ! returns an array of J on THRIFT_RHO by interpolating
    ! J(s) on the requested values of rho
    REAL(rprec), DIMENSION(:), INTENT(in) :: j_s_in
    REAL(rprec), DIMENSION(:), INTENT(out) :: j_rho_out
    INTEGER :: i, ier
    INTEGER :: bcs0(2)
    TYPE(EZspline1_r8) :: j_spl
    REAL(rprec) :: rho, s

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
!!===============================================================================
!!  PRINTER FUNCTIONS 
!!===============================================================================
SUBROUTINE print_calc_magvars()
    INTEGER :: i
    WRITE(6,*)'==============================================================================='
    WRITE(6,*)' CALCULATING MAGNETIC VARIABLES'
    WRITE(6,*)' S    DV/DS          <B>      <B^2>     RMAJOR     AMINOR        S11 '
    WRITE(6,*)''
    DO i = 1, nsj
        WRITE(6,'(F5.3,5(1X,F10.6),1X,ES10.3)') &
              THRIFT_S(i), THRIFT_VP(i,mytimestep), THRIFT_BAV(i,mytimestep), THRIFT_BSQAV(i,mytimestep), &
              ABS(THRIFT_S11(i,mytimestep)), THRIFT_RMAJOR(i,mytimestep), THRIFT_AMINOR(i,mytimestep)
    END DO
END SUBROUTINE print_calc_magvars

SUBROUTINE print_calc_abcd(j_arr)
    REAL(rprec), DIMENSION(:), INTENT(in) :: j_arr
    INTEGER :: i
    WRITE(6,*)'==============================================================================='
    WRITE(6,*)' CALCULATING COEFFICIENTS A,B,C,D'
    WRITE(6,*)'   S  ETAPARA       DV/DS      DP/DS       <J.B>      BSQAV        S11'
    WRITE(6,*)''
    DO i = 1, nsj
        WRITE(6,'(F5.3,6(1X,ES10.3))') &
        THRIFT_S(i), THRIFT_ETAPARA(i,mytimestep), THRIFT_VP(i,mytimestep), THRIFT_PPRIME(i,mytimestep),&
        j_arr(i)*THRIFT_BAV(i,mytimestep), THRIFT_BSQAV(i,mytimestep), THRIFT_S11(i,mytimestep)
    END DO
END SUBROUTINE print_calc_abcd

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
    WRITE(6,*)'  i  s        ITOTAL  ISOURCE        IPLASMA   rho   J     JPLASMA        JSOURCE'
    WRITE(6,*)''
    DO i = 1, nsj
        WRITE(6,'(I4, 1X, F5.3, 3(1X,ES13.5), 1X, F5.3, 3(1X,ES13.5))') &
           i, THRIFT_S(i), THRIFT_I(i,mytimestep), THRIFT_ISOURCE(i,mytimestep), THRIFT_IPLASMA(i,mytimestep),&
           THRIFT_RHO(i), j_arr(i), THRIFT_JPLASMA(i,mytimestep), THRIFT_JSOURCE(i,mytimestep)
     END DO
END SUBROUTINE print_postevolve

END MODULE thrift_funcs