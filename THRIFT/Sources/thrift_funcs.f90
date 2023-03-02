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
    ! AI, CI are arrays of size n-1, BI,DI of size n, val of size n
    IMPLICIT NONE
    REAL(rprec), DIMENSION(:), INTENT(in) :: AI
    REAL(rprec), DIMENSION(:), INTENT(in) :: BI
    REAL(rprec), DIMENSION(:), INTENT(in) :: CI
    REAL(rprec), DIMENSION(:), INTENT(in) :: DI
    REAL(rprec), DIMENSION(:), INTENT(out) :: val
    INTEGER :: i,n
    REAL(rprec), DIMENSION(:), ALLOCATABLE :: c_p
    REAL(rprec), DIMENSION(:), ALLOCATABLE :: d_p
    n = SIZE(DI)
    ALLOCATE(c_p(n),d_p(n))
    c_p = 0; d_p = 0
    !! Forward sweep
    ! c_1' = c_1/b_1   ;  c_i' =              c_i/(b_i - a_i*c_i-1') [1<i<=n-1]
    ! d_1' = d_1/b_1   ;  d_i' = (d_i-a_i*d_i-1')/(b_i - a_i*c_i-1') [1<i<=n  ]
    c_p(1) = CI(1)/BI(1) 
    d_p(1) = DI(1)/BI(1) 
    DO i = 2, n-1 
       c_p(i) =          CI(i)/(BI(i)-AI(i-1)*c_p(i-1))
       d_p(i) = (DI(i)-AI(i-1)*d_p(i-1))/(BI(i)-AI(i-1)*c_p(i-1))
    END DO
    c_p(n) = (DI(n)-AI(n-1)*d_p(n-1))/(BI(n)-AI(n-1)*c_p(n-1))
    !! Back substitution
    ! x_n = d_n'
    ! x_i = d_i' - c_i'*x_i+1
    val(n) = d_p(n) 
    DO i = n-1, 1, -1
       val(i) = d_p(i)-c_p(i)*val(i+1) 
    END DO

    DEALLOCATE(c_p,d_p)
    RETURN
END SUBROUTINE solve_tdm

SUBROUTINE check_sol(AI,BI,CI,DI,sol,residue)
    REAL(rprec), DIMENSION(:), INTENT(in) :: AI
    REAL(rprec), DIMENSION(:), INTENT(in) :: BI
    REAL(rprec), DIMENSION(:), INTENT(in) :: CI
    REAL(rprec), DIMENSION(:), INTENT(in) :: DI
    REAL(rprec), DIMENSION(:), INTENT(in) :: sol
    REAL(rprec), DIMENSION(:), INTENT(out) :: residue
    INTEGER :: i,n
    n = size(BI) 
    residue(1) = BI(1)*sol(1)+CI(1)*sol(2)-DI(1)
    DO i = 2, n-1
        residue(i) = AI(i-1)*sol(i-1)+BI(i)*sol(i)+CI(i)*sol(i+1)-DI(i)
    END DO
    residue(n) = AI(n-1)*sol(n-1)+BI(n)*sol(n)-DI(n)
    RETURN
END SUBROUTINE check_sol

SUBROUTINE curden_to_curtot(j_arr, aminor_arr, i_arr,timestep)
    ! Calculates array of enclosed current from current density array
    REAL(rprec), DIMENSION(:,:), INTENT(in) :: j_arr
    REAL(rprec), DIMENSION(:,:), INTENT(in) :: aminor_arr
    REAL(rprec), DIMENSION(:,:), INTENT(out) :: i_arr
    INTEGER, INTENT(in) :: timestep
    INTEGER :: i, n
    n = SIZE(j_arr, DIM=1)
    ! I(i) = I(i-1)+J(i)*dA(i)
    i_arr(1,timestep) = j_arr(1,timestep)*pi*(aminor_arr(2,2)**2-aminor_arr(1,2)**2)
    DO i = 2, n
        i_arr(i,timestep) = i_arr(i-1,timestep) + j_arr(i,timestep)*pi*(aminor_arr(i+1,2)**2-aminor_arr(i,2)**2)
    END DO
    RETURN
END SUBROUTINE curden_to_curtot

SUBROUTINE curtot_to_curden(i_arr, aminor_arr, j_arr,timestep)
    ! Calculates array of density from enclosed current array
    REAL(rprec), DIMENSION(:,:), INTENT(in) :: i_arr
    REAL(rprec), DIMENSION(:,:), INTENT(in) :: aminor_arr
    REAL(rprec), DIMENSION(:,:), INTENT(out) :: j_arr
    INTEGER, INTENT(in) :: timestep
    INTEGER :: i, n
    n = SIZE(i_arr, DIM=1)
    ! J(i) = dI(i)/dA(i)
    j_arr(1,timestep) = i_arr(1,timestep)/(pi*(aminor_arr(2,2)**2-aminor_arr(1,2)**2))
    DO i = 2, n
        j_arr(i,timestep) = (i_arr(i,timestep)-i_arr(i-1,timestep))/(pi*(aminor_arr(i+1,2)**2-aminor_arr(i,2)**2))
    END DO
    RETURN
END SUBROUTINE curtot_to_curden

SUBROUTINE deriv1_rho_o2(arr, step, der_arr)
    ! Calculate the first derivative (O2 accurate)
    ! Visualisation of different grids
    ! j=1 2    3    4    5    6    7
    !  |  |    |    |    |    |    |  ...   arr grid
    !     |    |    |    |    |    |  ...   der grid
    !    i=1   2    3    4    5    6        => j(i) = i+1

    REAL(rprec), DIMENSION(:), INTENT(in) :: arr
    REAL(rprec), INTENT(in) :: step
    REAL(rprec), DIMENSION(:), INTENT(out) :: der_arr
    INTEGER :: i,n
    n = SIZE(arr)-2
    der_arr(1) = (arr(3)+3*arr(2)-4*arr(1))/(3*step) ! dY/drho(1) = [Y(3) + 3*Y(2) - 4*Y(1)]/3h
    DO i = 2, n-1 ! dY/drho(i) = [Y(j+1)-Y(j-1)]/2h = [Y(i+2)-Y(i)]/2h
        der_arr(i) = (arr(i+2)-arr(i))/(2*step)
    END DO
    der_arr(n) = (4*arr(n+2)-3*arr(n+1)-arr(n))/(3*step) ! [4*Y(nrho+2) - 3*Y(nrho+1) - Y(nrho)]/3h
    RETURN
END SUBROUTINE deriv1_rho_o2

END MODULE thrift_funcs