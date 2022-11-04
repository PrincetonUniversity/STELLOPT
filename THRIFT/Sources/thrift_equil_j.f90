!-----------------------------------------------------------------------
!     Subroutine:    thrift_equil_j
!     Authors:       L. van Ham
!     Date:          11/XX/2022
!     Description:   This subroutine updated the equilbirium pressure
!-----------------------------------------------------------------------
      SUBROUTINE thrift_equil_j
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE thrift_runtime
      USE thrift_vars
      USE EZspline
      USE EZspline_obj
      USE vmec_input, ONLY: ac_aux_s, ac_aux_f, pcurr_type, ncurr, &
                            curtor
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: i, ier
      INTEGER :: bcs0(2)
      REAL(rprec) :: s_val, rho_val, j_val
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: rho_temp, j_temp
      TYPE(EZspline1_r8) :: j_spl
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      bcs0=(/ 0, 0/)

      ALLOCATE(rho_temp(nrho+2),j_temp(nrho+2))
      rho_temp(1)        = 0.0
      rho_temp(2:nrho+1) = THRIFT_RHO
      rho_temp(nrho+2)   = 1.0
      j_temp(1)          = 0.0
      j_temp(2:nrho+1)   = THRIFT_J(:,mytimestep)
      j_temp(nrho+2)     = 0.0

      CALL EZspline_init(j_spl,nrho+2,bcs0,ier)
      j_spl%x1        = rho_temp
      j_spl%isHermite = 1
      CALL EZspline_setup(j_spl,j_temp,ier,EXACT_DIM=.true.)

      DEALLOCATE(rho_temp,j_temp)

      IF (lvmec) THEN
         PCURR_TYPE = 'akima_spline_ip'
         NCURR = 1
         DO i = 1, 99
           s_val = DBLE(i-1)/DBLE(99-1)
           rho_val = sqrt(s_val)
           CALL EZspline_interp(j_spl,rho_val,j_val,ier)
           AC_AUX_S(i) = s_val
           AC_AUX_F(i) = j_val
         END DO
         CURTOR = SUM(AC_AUX_F(1:99),DIM=1)*(99-1)
      END IF

      CALL EZspline_free(j_spl,ier)

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_equil_j

