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
SUBROUTINE solve_tdm(AI,BI,CI,DI,val)
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
    ALLOCATE(c_p(nrho+2), d_p(nrho+2))
    c_p = 0; d_p = 0
    !! Forward sweep
    ! c_1' = c_1/b_1   ;  c_i' =              c_i/(b_i - a_i*c_i-1') [1<i<=n-1]
    ! d_1' = d_1/b_1   ;  d_i' = (d_i-a_i*d_i-1')/(b_i - a_i*c_i-1') [1<i<=n  ]
    c_p(1) = CI(1)/BI(1) 
    d_p(1) = DI(1)/BI(1) 
    DO i = 2, nrho+2
        denom = BI(i)-AI(i)*c_p(i-1)
        IF (i/=nrho+2) & 
          c_p(i) = CI(i)/denom
        d_p(i) = (DI(i)-AI(i)*d_p(i-1))/denom
    END DO 
    !! Back substitution
    ! x_n = d_n'
    ! x_i = d_i' - c_i'*x_i+1
    val(nrho+2) = d_p(nrho+2) 
    DO i = nrho+1, 1, -1
       val(i) = d_p(i)-c_p(i)*val(i+1) 
    END DO

    DEALLOCATE(c_p, d_p)
    RETURN

END SUBROUTINE solve_tdm

SUBROUTINE check_sol(AI,BI,CI,DI,sol,residue)
    REAL(rprec), DIMENSION(:), INTENT(in) :: AI
    REAL(rprec), DIMENSION(:), INTENT(in) :: BI
    REAL(rprec), DIMENSION(:), INTENT(in) :: CI
    REAL(rprec), DIMENSION(:), INTENT(in) :: DI
    REAL(rprec), DIMENSION(:), INTENT(in) :: sol
    REAL(rprec), DIMENSION(:), INTENT(out) :: residue

    residue(1) = BI(1)*sol(1)+CI(1)*sol(2)-DI(1)
    DO i = 2, nrho-1
        residue(i) = AI(i)*sol(i-1)+BI(i)*sol(i)+CI(i)*sol(i+1)-DI(i)
    END DO
    residue(nrho) = AI(nrho)*sol(nrho-1)+BI(nrho)*sol(nrho)-DI(nrho)

    RETURN

END SUBROUTINE check_sol


SUBROUTINE extrapolate_arr(j_arr, j_extr)
    REAL(rprec), DIMENSION(:), INTENT(IN)  :: j_arr
    REAL(rprec), DIMENSION(:), INTENT(OUT) :: j_extr

    j_extr(1) = (3*j_arr(1)-j_arr(2))/2
    j_extr(2:nrho+1) = j_arr
    j_extr(nrho+2) = (-j_arr(nrho-1)+3*j_arr(nrho))/2

    RETURN

END SUBROUTINE extrapolate_arr


SUBROUTINE curden_to_curtot(j_arr, i_arr)
    ! Calculates array of enclosed current from current density array
    REAL(rprec), DIMENSION(:), INTENT(in) :: j_arr
    REAL(rprec), DIMENSION(:), INTENT(out) :: i_arr
    INTEGER :: i
    REAL(rprec) :: drhoa
    REAL(rprec), DIMENSION(:), ALLOCATABLE :: rho_full
    ALLOCATE(rho_full(nrho+2))

    rho_full(1) = 0.0
    rho_full(2:nrho+1) = THRIFT_RHO
    rho_full(nrho+2) = 1.0

    ! I(1) = I(rho=0) = 0
    ! I(i) = I(i-1)+2*pi*J(i)*(rho(i)*a(i))*d(rho*a)(i)
    i_arr(1) = 0
    DO i = 2, nrho+2
        drhoa = rho_full(i)*THRIFT_AMINOR(i,2)-rho_full(i-1)*THRIFT_AMINOR(i-1,2) !d(rho*a)
        i_arr(i) = i_arr(i-1) + 2*pi*j_arr(i)*rho_full(i)*THRIFT_AMINOR(i,2)*drhoa
    END DO
    
    DEALLOCATE(rho_full)
    RETURN

END SUBROUTINE curden_to_curtot



SUBROUTINE curtot_to_curden(i_arr, j_arr)
    ! Calculates array of density from enclosed current array
    REAL(rprec), DIMENSION(:), INTENT(in) :: i_arr
    REAL(rprec), DIMENSION(:), INTENT(out) :: j_arr
    REAL(rprec), DIMENSION(:), ALLOCATABLE :: dIdrho, dAdrho, rho_full, j_temp
    ALLOCATE(dIdrho(nrho+2), dAdrho(nrho+2), rho_full(nrho+2), j_temp(nrho_2))

    rho_full(1) = 0.0
    rho_full(2:nrho+1) = THRIFT_RHO
    rho_full(nrho+2) = 1.0

    CALL deriv1_rho_o2(i_arr, dIdrho)
    ! dA/drho = 2*pi*(rho*a**2)
    dAdrho = 2*pi*rho_full*(THRIFT_AMINOR(:,2))**2

    ! J(i) = dI/dA = dI/drho*drho/dA
    j_temp = dIdrho/dAdrho
    CALL extrapolate_arr(j_temp, j_arr)
    
    DEALLOCATE(dIdrho, dAdrho, rho_full, j_temp)
    RETURN

END SUBROUTINE curtot_to_curden

SUBROUTINE deriv1_rho_o2(arr, der_arr)
    ! Calculate the first derivative (O2 accurate)
    REAL(rprec), DIMENSION(:), INTENT(in) :: arr
    REAL(rprec), DIMENSION(:), INTENT(out) :: der_arr
    INTEGER :: i
    REAL(rprec) :: step
    step = THRIFT_RHO(2) - THRIFT_RHO(1)

    ! Derivatives at 1,nrho not necessary
    der_arr(1) = 0
    der_arr(2) = (arr(3)+3*arr(2)-4*arr(1))/(3*step) ! dY/drho(1) = [Y(3) + 3*Y(2) - 4*Y(1)]/3h
    DO i = 3, nrho-2 ! dY/drho(i) = [Y(j+1)-Y(j-1)]/2h 
        der_arr(i) = (arr(i+1)-arr(i-1))/(2*step)
    END DO
    der_arr(nrho-1) = (4*arr(nrho)-3*arr(nrho-1)-arr(nrho-2))/(3*step) ! [4*Y(n) - 3*Y(n-1) - Y(n-2)]/3h
    der_arr(nrho) = 0

    RETURN

END SUBROUTINE deriv1_rho_o2

END MODULE thrift_funcs