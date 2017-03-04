!-----------------------------------------------------------------------
!     Subroutine:    stelltran_bal_update
!     Authors:       J. Mittelstaedt (jmittelstaedt@uchicago.edu)
!     Date:          06/21/2016
!     Description:   This subroutine updates ne, using the particle
!                    and power balance equations. Uses implicit Euler method
!                    as described in http://w3.pppl.gov/~hammett/gyrofluid/papers/2008/Jardin-JCP-corrected.pdf
!-----------------------------------------------------------------------
      SUBROUTINE stelltran_bal_update
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE stelltran_runtime
      USE EZspline_obj
      USE EZspline
      USE stelltran_equilutils
      USE stelltran_vars
!-----------------------------------------------------------------------
!     Local Variables
!     RHS               Right hand side of the differential equation
!     init_nete         Initial value of n_e*T_e for electron power balance update
!     init_niti         Initial value of n_i*T_i for ion power balance update
!     tau_e             Collisional time for electrons for ion collisional heating
!     lambda            Coulomb logarithm
!     dr                Spatial resolution of our grid.
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), DIMENSION(prof_length) :: rho_coord, gradrho, gradrhosq, Vp_rho2, Spe_rho, ne_rho, Vps
      REAL(rprec), DIMENSION(prof_length) :: perf_rho, picol_rho, te_rho, ti_rho
      REAL(rprec), DIMENSION(prof_length-1) :: Vp_rho, De_rho, gradrhosq_rho, Xe_rho, rho_grid
      REAL(rprec) :: alpha, dr, dt ! ********** dt IS ONLY TEMPORARY
      REAL(rprec), DIMENSION(prof_length) :: A_up, B_mid, C_dn, D_src, ne_sol
      INTEGER :: j, ier
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      print*, '********* in the balance updating thing ***************************'
      dt = 0.001 ! Setting this to 1ms at first, should change when this is properly incorporated
      dr = 1./63.
      alpha = dt/dr**2.


      ! Find needed geometric factors.
      DO j = 1, prof_length
         ! Given S this function returns rho,Vp,<|grad(rho)|>,<|grad(rho)|^2>
         CALL get_equil_rho(Vp(j,1),rho_coord(j),Vps(j),gradrho(j),gradrhosq(j),ier)
      END DO

      ! Convert dV/dPhi to dV/drho
      Vp(:,2) = 2*rho_coord(:)*Vps(:)

      DO j=1,prof_length-1
            rho_grid(j) = REAL(2*j-1)/REAL(2*prof_length - 2) ! midpoints of uniform grid in rho
      END DO

      ! Put parameters onto uniform grid of midpoints of rho
      ! NOTE: If this ends up actually being used, would need to change to use the things in s_rho_vars
      !       Also, from implementation, probably more efficient to have all in their own do loops.
      DO j=1,prof_length-1
            CALL eval_prof_spline(prof_length,rho_coord,Vp(:,2),rho_grid(j),Vp_rho(j),ier) ! V'
            CALL eval_prof_spline(prof_length,rho_coord,De(:,2),rho_grid(j),De_rho(j),ier) ! D_e
            CALL eval_prof_spline(prof_length,rho_coord,Xe(:,2),rho_grid(j),Xe_rho(j),ier) ! Chi_e
            CALL eval_prof_spline(prof_length,rho_coord,Xi(:,2),rho_grid(j),Xi_rho(j),ier) ! Chi_i
            CALL eval_prof_spline(prof_length,rho_coord,gradrhosq,rho_grid(j),gradrhosq_rho(j),ier) !<|grad(rho)|^2>
      END DO

      ! Put these onto regular grid in rho
      DO j=1,prof_length
            CALL eval_prof_spline(prof_length,rho_coord,Vp(:,2),Vp(j,1),Vp_rho2(j),ier) ! V'
            CALL eval_prof_spline(prof_length,rho_coord,S_pe(:,2),S_pe(j,1),Spe_rho(j),ier) ! S_e
            CALL eval_prof_spline(prof_length,rho_coord,pe_rf(:,2),pe_rf(j,1),perf_rho(j),ier) ! pe_rf
            CALL eval_prof_spline(prof_length,rho_coord,pi_col(:,2),pi_col(j,1),picol_rho(j),ier) ! pi_col
            CALL eval_prof_spline(prof_length,rho_coord,ne(:,2),ne(j,1),ne_rho(j),ier) ! n_e
            CALL eval_prof_spline(prof_length,rho_coord,te(:,2),te(j,1),te_rho(j),ier) ! T_e
            CALL eval_prof_spline(prof_length,rho_coord,ti(:,2),ti(j,1),ti_rho(j),ier) ! T_i
      END DO

      ! Calculate matrix elements, electron particle balance (Jardin)
      C_dn(1) = 0.
      A_up(prof_length) = 0.
      DO i=1,prof_length-1
            A_up(i) = -1.*alpha*(Vp_rho(i)*gradrhosq_rho(i)*De_rho(i))/Vp_rho2(i) ! from 1 - prof_len-1
            C_dn(i+1) = -1.*alpha*(Vp_rho(i)*gradrhosq_rho(i)*De_rho(i))/Vp_rho2(i+1) ! from 2 - prof_len
      END DO
      B_mid = -1. - A(:) - C(:)
      D_src = -1.*ne_rho(:) - dt*Spe_rho(:)

      ! Solve the equation for the updated ne on a uniform grid in rho
      ! NOTE: This ~wants diagonal element in row is larger than sum of other elements in row
      solve_tridiag(C_dn,B_mid,A_up,D_src,ne_sol,prof_length)

      ! Calculate matrix elements, electron power balance (Jardin)
      ! We are assuming that we're solving for DteneDt, and internal deriv is DteneDr, *not* neDteDr
      C_dn(1) = 0.
      A_up(prof_length) = 0.
      DO i=1,prof_length-1
            A_up(i) = -1.*alpha*(Vp_rho(i)*gradrhosq_rho(i)*Xe_rho(i))/Vp_rho2(i) ! from 1 - prof_len-1
            C_dn(i+1) = -1.*alpha*(Vp_rho(i)*gradrhosq_rho(i)*Xe_rho(i))/Vp_rho2(i+1) ! from 2 - prof_len
      END DO
      B_mid = -1. - A(:) - C(:)
      D_src = -1.*ne_rho(:)*te_rho(:) - dt*(perf_rho(:)-picol_rho(:))

      ! Solve the equation for the updated ne on a uniform grid in rho
      ! NOTE: This ~wants diagonal element in row is larger than sum of other elements in row
      solve_tridiag(C_dn,B_mid,A_up,D_src,ne_sol,prof_length)

      ! Update in the magnetic s coordinates
      DO j=1,prof_length
            CALL eval_prof_spline(prof_length,ne(:,1),ne_sol,rho_coord(j),ne(j,2),ier)
      END DO

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE stelltran_bal_update

      subroutine solve_tridiag(a,b,c,d,x,n)
      implicit none
!	 a - sub-diagonal (means it is the diagonal below the main diagonal)
!	 b - the main diagonal
!	 c - sup-diagonal (means it is the diagonal above the main diagonal)
!	 d - right part
!	 x - the answer
!	 n - number of equations

        integer,intent(in) :: n
        real(rprec),dimension(n),intent(in) :: a,b,c,d
        real(rprec),dimension(n),intent(out) :: x
        real(rprec),dimension(n) :: cp,dp
        real(rprec) :: m
        integer :: i

! initialize c-prime and d-prime
        cp(1) = c(1)/b(1)
        dp(1) = d(1)/b(1)
! solve for vectors c-prime and d-prime
         do i = 2,n
           m = b(i)-cp(i-1)*a(i)
           cp(i) = c(i)/m
           dp(i) = (d(i)-dp(i-1)*a(i))/m
         enddo
! initialize x
         x(n) = dp(n)
! solve for x from the vectors c-prime and d-prime
        do i = n-1, 1, -1
          x(i) = dp(i)-cp(i)*x(i+1)
        end do

    end subroutine solve_tridiag