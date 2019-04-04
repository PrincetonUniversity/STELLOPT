!-----------------------------------------------------------------------
!     Subroutine:    stellopt_analytic
!     Authors:       J.C.Schmitt (Auburn/PPPL) jcschmitt@auburn.edu
!     Date:          2018
!     Description:   This subroutine calculates an analytic function
!                    
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_analytic(lscreen, iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_input_mod
      USE stellopt_vars
      USE equil_utils

!-----------------------------------------------------------------------
!     Subroutine Parameters
!        iflag         Error flag
!----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(inout) :: iflag
      LOGICAL, INTENT(inout)        :: lscreen

!-----------------------------------------------------------------------
!     Local Variables
!        iunit       File unit number
!----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Local Variables
!        istat         Error status
!        iunit         File unit number
      ! FOR REGCOIL
      INTEGER :: istat, iunit, ii
      real(rprec) :: mysum

!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
!      IF (iflag < 0) RETURN
      !lscreen = .true.
      IF (lscreen) then
         WRITE(6,'(a)') ' -------------  Analytic CALCULATION  ---------'
      ENDIF

      verbose = lscreen 

      DO ii = 1, INT(analytic_fcnt)
       ! print *, '<---analytic: coeff:', analytic_coeff(ii), &
       !         'x:', analytic_x,',  x_off:', analytic_x_off(ii), &
       !         'x_pow:', analytic_x_pow(ii), &
       !         'y:', analytic_y,',  y_off:', analytic_y_off(ii), &
       !         'y_pow:', analytic_y_pow(ii), &
       !         'z:', analytic_z,',  z_off:', analytic_z_off(ii), &
       !         'z_pow:', analytic_z_pow(ii)

        mysum = 0.0
        mysum = mysum + analytic_coeff(ii) * &
                (analytic_x - analytic_x_off(ii))**analytic_x_pow(ii) * &
                (analytic_y - analytic_y_off(ii))**analytic_y_pow(ii) * &
                (analytic_z - analytic_z_off(ii))**analytic_z_pow(ii) 
       ! print *, '<====sum = ', mysum
        analytic_eval(ii) = mysum
      END DO


      IF (lscreen) WRITE(6,'(a)') ' ---------------------------  Analytic CALCULATION DONE  ---------------------'

      RETURN

!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_analytic
