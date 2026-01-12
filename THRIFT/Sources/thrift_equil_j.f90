!-----------------------------------------------------------------------
!     Subroutine:    thrift_equil_j
!     Authors:       L. van Ham
!     Date:          11/XX/2022
!     Description:   This subroutine updates the equilibrium dI/ds
!-----------------------------------------------------------------------
      SUBROUTINE thrift_equil_j
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE thrift_runtime
      USE thrift_vars
      USe thrift_equil
      USE EZspline
      USE EZspline_obj
      USE stel_tools
      USE vmec_input, ONLY: ac_aux_s, ac_aux_f, pcurr_type, ncurr, &
                            curtor
      USE thrift_funcs, ONLY: curden_to_curtot
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: i, ier, itime
      INTEGER :: bcs0(2)
      REAL(rprec) :: s_val, j_val, temp
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: s_temp
      TYPE(EZspline1_r8) :: j_spl, i_spl
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      ! If first pass then just set everything to previous timestep
      IF (mytimestep.eq.1) THEN
         itime = mytimestep
      ELSE IF (nsubsteps.eq.1) THEN ! On first substep (after first timestep!), look at prev timestep value
         itime = mytimestep-1
      ELSE ! On any substep after first, look at prev substep value
         itime = mytimestep
      END IF

      ! THRIFT_J is set to 0.0 at begining of thrift_evolve
      ! so need to update it here if code is in restart mode
      IF(lrestart_from_file .AND. mytimestep.eq.1 .AND. nsubsteps .eq. 1) THEN
            ! THRIFT_J(:,itime) = J_RESTART  --> using J directly can be problematic sometimes. better to use the smoothed I
            ! THRIFT_I(:,itime) = UGRID_RESTART*eq_phiedge/mu0 --> better, but leads to small discontinuities of iota at restart times
            ! since UGRID_RESTART does not know about fpicard, while j does. The following is consistent with previous and following iters:
            ! CALL curden_to_curtot(THRIFT_J(:,itime),THRIFT_I(:,itime))  ! This may give the wrong sing?
            CALL curden_to_curtot(ABS(J_RESTART)*UGRID_RESTART/ABS(UGRID_RESTART),THRIFT_I(:,itime)) ! This for sure gives the correct current's sign

      END IF

      ! Feed VMEC with dI/ds instead of J (as was done before, see commented lines below)
      bcs0=(/ 0, 0/)
      ier = 0
      CALL EZspline_init(i_spl,nsj,bcs0,ier)
      i_spl%x1        = THRIFT_S
      i_spl%isHermite = 1
      CALL EZspline_setup(i_spl,THRIFT_I(:,itime),ier,EXACT_DIM=.true.)
      ! temporary s grid 
      ALLOCATE(s_temp(n_eq))
      s_temp = 0
      FORALL(i = 1:n_eq) s_temp(i) = DBLE(i-1)/DBLE(n_eq-1)
      IF(lvmec) THEN
            PCURR_TYPE = 'akima_spline_ip'
            NCURR = 1
            CALL EZspline_derivative1_array_r8(i_spl, 1, n_eq,s_temp,AC_AUX_F(1:n_eq),ier)
            CURTOR = THRIFT_I(nsj,itime)
            AC_AUX_S(1:n_eq) = s_temp
      END IF

!       ! Create J spline
!       bcs0=(/ 0, 0/)
!       ier = 0
!       CALL EZspline_init(j_spl,nsj,bcs0,ier)
!       j_spl%x1        = THRIFT_S
!       j_spl%isHermite = 1
!       CALL EZspline_setup(j_spl,THRIFT_J(:,itime),ier,EXACT_DIM=.true.)

!       ! temporary s grid 
!       ALLOCATE(s_temp(n_eq))
!       s_temp = 0
!       FORALL(i = 1:n_eq) s_temp(i) = DBLE(i-1)/DBLE(n_eq-1)

!       ! Update equilibrium inputs
!       IF (lvmec) THEN
!          ! VMEC requires dI/ds be specified in AC_AUX_F
!          !   A     = pi*(rho*aminor)^2 = pi*s*aminor^2
!          !   dA/ds = pi*aminor^2
!          !   j     = dI/ds * ds/dA = dI/ds / (dA/ds)
!          !   dI/ds = j * dA/ds
!          !         = j*pi*aminor^2

!          PCURR_TYPE = 'akima_spline_ip'
!          NCURR = 1
! !        Don't think this is necessary anymore
! !         ! Check to make sure we have dV/ds and Aminor
! !         IF (EZspline_allocated(vp_spl)) THEN
!             DO i = 1, n_eq
!                s_val = s_temp(i)
!                CALL EZspline_interp(j_spl,s_val,j_val,ier)
!                AC_AUX_S(i) = s_val
!                AC_AUX_F(i) = j_val*pi*eq_Aminor**2
!             END DO
!             CURTOR = SUM(AC_AUX_F(1:n_eq),DIM=1)/DBLE(n_eq-1)
! !         ELSE
! !            DO i = 1, n_eq
! !               s_val = s_temp(i)
! !               AC_AUX_S(i) = s_val
! !            END DO
! !            CURTOR = 0
! !         END IF
!       END IF

      ! Deallocate
      DEALLOCATE(s_temp)
      CALL EZspline_free(j_spl,ier)
      CALL EZspline_free(i_spl,ier)

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_equil_j

