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
    INTEGER :: i
    
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
    REAL(rprec) :: drho

    ! I(1) = I(rho=0) = 0
    i_arr(1) = 0
    DO i = 2, nrho+2
        drho = (THRIFT_RHOFULL(i)-THRIFT_RHOFULL(i-1))
        ! I(i) = I(i-1) + J*dA/drho*drho
        i_arr(i) = i_arr(i-1) + j_arr(i)*(2*pi*THRIFT_RHOFULL(i)*THRIFT_AMINOR(i,mytimestep)**2)*drho
    END DO
    
    RETURN

END SUBROUTINE curden_to_curtot



SUBROUTINE curtot_to_curden(I_arr, j_arr)
    ! Calculates array of density from enclosed current array
    REAL(rprec), DIMENSION(:), INTENT(in) :: I_arr
    REAL(rprec), DIMENSION(:), INTENT(out) :: j_arr
    REAL(rprec), DIMENSION(:), ALLOCATABLE :: js_temp,  s_temp, jrho_temp
    INTEGER :: i, ier, ns
    INTEGER :: bcs0(2)
    TYPE(EZspline1_r8) :: I_spl
    REAL(rprec) :: ds, I_temp1, I_temp2,dIds, temp, Aminor

    ! number of gridpoints in s space
    ns = 101;

    ! Allocate s grid and temp j grids
    ALLOCATE(s_temp(ns), js_temp(ns), jrho_temp(nrho))
    s_temp = 0; js_temp = 0; jrho_temp = 0;
    
    ! Setup s grid
    ds = 1.0/(ns-1)
    DO i = 1, ns
        s_temp(i) = (i-1)*ds
    END DO

     ! Create I spline (in rho space)
    bcs0=(/ 0, 0/)
    CALL EZspline_init(I_spl,nrho+2,bcs0,ier)
    I_spl%x1        = THRIFT_RHOFULL
    I_spl%isHermite = 1
    CALL EZspline_setup(I_spl,I_arr,ier,EXACT_DIM=.true.)

    ! Calculate J in s space
    js_temp(1) = 0
    DO i = 2, ns-1
        ! Spline in rho space; rho = sqrt(s)
        ! Gives I at the required s values
        CALL EZSpline_interp(I_spl,sqrt(s_temp(i+1)),I_temp1, ier )
        CALL EZSpline_interp(I_spl,sqrt(s_temp(i-1)),I_temp2, ier )
        ! dI/ds
        dIds = (I_temp1-I_temp2)/(s_temp(i+1)-s_temp(i-1))
        ! a_eff
        CALL get_equil_Rmajor(s_temp(i),temp, temp, Aminor, ier )
        ! J = dI/ds*ds/dA = dI/ds/(pi*a_eff^2)
        js_temp(i) = dIds/(pi*Aminor**2)
    END DO
    js_temp(ns) = 0
    
    ! Free spline
    CALL EZspline_free(I_spl,ier)

    ! Setup J spline (in s space)
    bcs0=(/ 0, 0/)
    CALL EZspline_init(I_spl,ns,bcs0,ier)
    I_spl%x1        = s_temp
    I_spl%isHermite = 1
    CALL EZspline_setup(I_spl,js_temp,ier,EXACT_DIM=.true.)   

    ! Interpolate to find J in rho space
    DO i = 1, nrho
        ! Spline in s space: s = rho^2
        CALL EZSpline_interp(I_spl,THRIFT_RHO(i)**2, jrho_temp(i), ier)
    END DO

    ! Free spline
    CALL EZspline_free(I_spl,ier)
    CALL extrapolate_arr(jrho_temp, j_arr)
    
    DEALLOCATE(s_temp, js_temp, jrho_temp)
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
    der_arr(nrho+1) = (4*arr(nrho+2)-3*arr(nrho+1)-arr(nrho-2))/(3*step) ! [4*Y(n) - 3*Y(n-1) - Y(n-2)]/3h
    der_arr(nrho+2) = 0

    RETURN

END SUBROUTINE deriv1_rho_o2

END MODULE thrift_funcs