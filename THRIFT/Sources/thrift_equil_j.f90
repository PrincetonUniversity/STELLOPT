!-----------------------------------------------------------------------
!     Subroutine:    thrift_equil_j
!     Authors:       L. van Ham
!     Date:          11/XX/2022
!     Description:   This subroutine updates the equilibrium dI/ds
!-----------------------------------------------------------------------
      SUBROUTINE thrift_equil_j(lfirst_pass)
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
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: lfirst_pass
      INTEGER :: i, ier, itime
      INTEGER :: bcs0(2)
      REAL(rprec) :: s_val, rho_val, j_val, vp, Rmajor, Aminor, temp
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: rho_temp, j_temp
      TYPE(EZspline1_r8) :: j_spl
      INTEGER, PARAMETER :: n_eq = 99
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      ! If first pass and first timestep just set everything to zero
      !IF (lfirst_pass and mytimestep==1) THEN
      !   IF (lvmec) THEN
      !      NCURR  = 1
      !      CURTOR = 0
      !   END IF
      !   RETURN
      !END IF

      ! If first pass then just set everything to previous timestep
      IF (mytimestep>1) THEN
         itime = mytimestep-1
      ELSE
         itime = mytimestep
      END IF

      ! Create J array
      ALLOCATE(rho_temp(nrho+2),j_temp(nrho+2))
      rho_temp(1)        = 0.0
      !rho_temp(1)        = -THRIFT_RHO(1)
      rho_temp(2:nrho+1) = THRIFT_RHO
      rho_temp(nrho+2)   = 1.0
      !j_temp(1)          = 0.0
      !j_temp(1)          = THRIFT_J(1,itime)
      j_temp(1)          = (3*THRIFT_J(1,itime)-THRIFT_J(2,itime))/2
      j_temp(2:nrho+1)   = THRIFT_J(:,itime)
      j_temp(nrho+2)     = (3*THRIFT_J(nrho,itime)-THRIFT_J(nrho-1,itime))/2

      ! Create Splines
      bcs0=(/ 0, 0/)
      CALL EZspline_init(j_spl,nrho+2,bcs0,ier)
      j_spl%x1        = rho_temp
      j_spl%isHermite = 1
      CALL EZspline_setup(j_spl,j_temp,ier,EXACT_DIM=.true.)

      ! Deallocate helper arrays
      DEALLOCATE(rho_temp,j_temp)

      ! Update equilibrium inputs
      IF (lvmec) THEN
         ! VMEC requires dI/ds be specified in AC_AUX_F
         !   dA/ds = dV/ds/(2*pi*R)
         !   dV/ds = dV/dPhi * dPhi/ds
         !   dPhi/ds = dPhi/drho * drho/ds = dPhi/drho / (2*rho)
         !   j     = dI/ds * ds/dA = dI/ds / (dA/ds)
         !   dI/ds = j * dA/ds
         !         = j * dV/ds / (2*pi*R)
         !         = j * dV/dPhi * dPhi/ds / (2*pi*R)
         !         = j * dV/dPhi * dPhi/drho / (4*pi*rho*R)
         !         = j * dV/dPhi * Phi_edge / (2*pi*R)
         PCURR_TYPE = 'akima_spline_ip'
         NCURR = 1
         ! Check to make sure we have dV/ds and Aminor
         IF (EZspline_allocated(vp_spl)) THEN
            DO i = 1, n_eq
               s_val = DBLE(i-1)/DBLE(n_eq-1)
               rho_val = sqrt(s_val)
               CALL EZspline_interp(j_spl,rho_val,j_val,ier)
               CALL EZspline_interp(vp_spl,rho_val,vp,ier) ! dV/dPhi
               !CALL EZspline_interp(phip_spl,rho_val,phip,ier) ! dPhi/drho
               AC_AUX_S(i) = s_val
               !AC_AUX_F(i) = j_val*vp*phip/(2*pi2*rho_val*eq_Rmajor)
               CALL get_equil_Rmajor(s_val, Rmajor, temp, Aminor, ier)
               !AC_AUX_F(i) = j_val*vp*eq_phiedge/(2*pi*Rmajor)
               AC_AUX_F(i) = j_val*pi*Aminor**2
            END DO
            !AC_AUX_F(1) = 2*AC_AUX_F(2)-AC_AUX_F(3)
            CURTOR = SUM(AC_AUX_F(1:n_eq),DIM=1)/DBLE(n_eq-1)
         ELSE
            DO i = 1, n_eq
               s_val = DBLE(i-1)/DBLE(n_eq-1)
               AC_AUX_S(i) = s_val
            END DO
            CURTOR = 0
         END IF
      END IF

      ! Deallocate Spline
      CALL EZspline_free(j_spl,ier)

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_equil_j

