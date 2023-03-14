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
    ALLOCATE(c_p(nssize), d_p(nssize))
    c_p = 0; d_p = 0
    !! Forward sweep
    ! c_1' = c_1/b_1   ;  c_i' =              c_i/(b_i - a_i*c_i-1') [1<i<=n-1]
    ! d_1' = d_1/b_1   ;  d_i' = (d_i-a_i*d_i-1')/(b_i - a_i*c_i-1') [1<i<=n  ]
    c_p(1) = CI(1)/BI(1) 
    d_p(1) = DI(1)/BI(1) 
    DO i = 2, nssize
        denom = BI(i)-AI(i)*c_p(i-1)
        IF (i/=nssize) & 
          c_p(i) = CI(i)/denom
        d_p(i) = (DI(i)-AI(i)*d_p(i-1))/denom
    END DO 
    !! Back substitution
    ! x_n = d_n'
    ! x_i = d_i' - c_i'*x_i+1
    val(nssize) = d_p(nssize) 
    DO i = nssize-1, 1, -1
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

    ALLOCATE(residue(nssize))
    
    residue(1) = BI(1)*sol(1)+CI(1)*sol(2)-DI(1)
    DO i = 2, nssize-1
        residue(i) = AI(i)*sol(i-1)+BI(i)*sol(i)+CI(i)*sol(i+1)-DI(i)
    END DO
    residue(nssize) = AI(nssize)*sol(nssize)+BI(nssize)*sol(nssize)-DI(nssize)
    WRITE(6,*) maxval(residue)

    DEALLOCATE(residue)

    RETURN

END SUBROUTINE check_sol


SUBROUTINE extrapolate_arr(j_arr, j_extr)
    REAL(rprec), DIMENSION(:), INTENT(IN)  :: j_arr
    REAL(rprec), DIMENSION(:), INTENT(OUT) :: j_extr

    j_extr(1)        = (3*j_arr(1)-j_arr(2))/2
    j_extr(2:nrho+1) = j_arr
    j_extr(nrho+2)   = (-j_arr(nrho-1)+3*j_arr(nrho))/2

    RETURN

END SUBROUTINE extrapolate_arr


SUBROUTINE curden_to_curtot(j_arr, i_arr)
    ! Takes a J(rho) array and returns an I(s) array.
    REAL(rprec), DIMENSION(:), INTENT(in) :: j_arr
    REAL(rprec), DIMENSION(:), INTENT(out) :: i_arr
    REAL(rprec), DIMENSION(:), ALLOCATABLE :: j_temp_spl
    INTEGER :: i, ier
    INTEGER :: bcs0(2)
    REAL(rprec) :: s,rho,ds,j_temp
    TYPE(EZspline1_r8) :: j_spl

    ALLOCATE(j_temp_spl(nrho+2))

    ds = THRIFT_S(2)-THRIFT_S(1)
    CALL extrapolate_arr(j_arr,j_temp_spl)
    ! Create J spline (in rho space)

    bcs0=(/ 0, 0/)
    CALL EZspline_init(j_spl,nrho+2,bcs0,ier)
    j_spl%x1        = THRIFT_RHOFULL
    j_spl%isHermite = 1
    CALL EZspline_setup(j_spl,j_temp_spl,ier,EXACT_DIM=.true.)
    
    ! Calculate I (in s space)
    i_arr(1) = 0
    DO i = 2, nssize
        s = THRIFT_S(i)
        rho = SQRT(s)
        ier = 0
        CALL EZspline_interp(j_spl, rho, j_temp, ier)
        i_arr = i_arr(i-1) + j_temp*(pi*THRIFT_AMINOR(i,mytimestep)**2)*ds
    END DO
    CALL EZspline_free(j_spl,ier)
    DEALLOCATE(j_temp_spl)

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

    ALLOCATE(j_temp(nssize-2))

    ! Calculate J (in s space)
    DO i = 2, nssize-1
       j_temp(i) = (i_arr(i+1)-i_arr(i-1))/(2*ds)*1.0/(pi*THRIFT_AMINOR(i,mytimestep)**2)
    END DO

    ! Extrapolate to boundaries
    j_temp(1)  = 2*j_temp(2)   -j_temp(3)    ! s = 0
    j_temp(nssize) = 2*j_temp(nssize-1)-j_temp(nssize-2) ! s = 1

    ! Setup J spline (in s space)
    CALL EZspline_init(j_spl,nssize,bcs0,ier)
    j_spl%x1        = THRIFT_S
    j_spl%isHermite = 1
    CALL EZspline_setup(j_spl,j_temp,ier,EXACT_DIM=.true.)  

    ! Convert J to rho space
    DO i = 1, nrho
       rho = THRIFT_RHO(i)
       s = rho*rho
       ier = 0
       CALL EZspline_interp(j_spl, s, j_arr(i), ier)
    END DO

    CALL EZspline_free(j_spl,ier)
    DEALLOCATE(j_temp)
    RETURN

END SUBROUTINE curtot_to_curden

SUBROUTINE deriv1_rho_o2(arr, der_arr)
    ! Calculate the first derivative (O2 accurate)
    REAL(rprec), DIMENSION(:), INTENT(in) :: arr
    REAL(rprec), DIMENSION(:), INTENT(out) :: der_arr
    INTEGER :: i
    REAL(rprec) :: step
    step = THRIFT_RHO(2) - THRIFT_RHO(1)

    ! Set derivatives = 0 on boundaries
    der_arr(1) = 0
    der_arr(2) = (arr(3)+3*arr(2)-4*arr(1))/(3*step) ! dY/drho(1) = [Y(3) + 3*Y(2) - 4*Y(1)]/3h
    DO i = 3, nrho ! dY/drho(i) = [Y(j+1)-Y(j-1)]/2h 
        der_arr(i) = (arr(i+1)-arr(i-1))/(2*step)
    END DO
    der_arr(nrho+1) = (4*arr(nrho+2)-3*arr(nrho+1)-arr(nrho))/(3*step) ! [4*Y(n) - 3*Y(n-1) - Y(n-2)]/3h
    der_arr(nrho+2) = 0

    RETURN

END SUBROUTINE deriv1_rho_o2

END MODULE thrift_funcs