!-----------------------------------------------------------------------------
!+ Module for variables read from profile_data_***
!-----------------------------------------------------------------------------
Module vmec_var_pass

!
! Description:
!   This module contains the plasma profile data variables.
!   For definitions see subroutine read_vmec_file in read_input_file_mod.f90.
!   For directions on how to generate the profile_data_*** file see penta.f90.
!
! History:
! Version   Date      Comment
! -------   ----      -------
!  1.0     01/7/2009   Original Code.  JL
!  1.1     5/24/2010   Updated for PENTA3. JL 
! 
! Author(s): J. Lore 7/2009 - 5/24/2010 
!

  ! Modules used:
  Use penta_kind_mod                ! Import rknd, iknd specifications

  Implicit none

  Real(rknd) :: arad, Rmajor, r_surf, roa_surf, chip  &
    ,psip, btheta, bzeta, Bsq, iota, vol_p, B0
End module vmec_var_pass
!- End of header -------------------------------------------------------------
