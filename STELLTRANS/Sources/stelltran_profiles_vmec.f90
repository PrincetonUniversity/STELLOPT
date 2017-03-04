!-----------------------------------------------------------------------
!     Subroutine:    stelltran_profiles_vmec
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          08/29/2015
!     Description:   This routine updates the VMEC profiles.
!-----------------------------------------------------------------------
      SUBROUTINE stelltran_profiles_vmec
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE stelltran_runtime
      USE stelltran_vars
      USE stelltran_data, ONLY: Ip
      USE vmec_input, ONLY: pmass_type, pcurr_type, am, ac, ncurr, &
                            ac_aux_s, ac_aux_f, am_aux_s, am_aux_f, &
                            curtor, pres_scale
      USE EZspline_obj
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Default the profile types
      ncurr = 1
      pmass_type='Akima_spline'
      pcurr_type='Akima_spline_Ip'
      
      ! Now update the am_array
      !pres_scale = 1
      am_aux_s(1:prof_length) = te(:,1)
      am_aux_f(1:prof_length) = ne(:,2)*ec*(te(:,2)+ti(:,2)/zeff(:,2))
      
      ! Now update the ac_array
      ac_aux_s(1:prof_length) = jrf(:,1)
      ac_aux_f(1:prof_length) = johm(:,2)+jboot(:,2)+jrf(:,2)+jbeam(:,2)
      
      ! Handle renormalizing if no measured value
      IF (.not. ALLOCATED(Ip)) &
      curtor = SUM(ac_aux_f(2:prof_length) * &
                   ( ac_aux_s(2:prof_length) - ac_aux_s(1:prof_length-1) ))
      
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE stelltran_profiles_vmec
