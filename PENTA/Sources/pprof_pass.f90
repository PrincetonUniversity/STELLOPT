!-----------------------------------------------------------------------------
!+ Module for variables read from plasma_profiles*.dat
!-----------------------------------------------------------------------------
Module pprof_pass

!
! Description:
!   This module contains the plasma profile variables
!   For definitions see subroutine read_pprof_file in read_input_file_mod.f90.
!   For directions on how to generate the .dat file see penta.f90
!
! History:
! Version   Date      Comment
! -------   ----      -------
!  1.0     01/7/2009   Original Code.  JL
!  1.1     5/24/2010   Updated for PENTA3. JL 
! 
! Author(s): J. Lore 7/2009 - 5/24/2010 

  ! Modules used:
  Use penta_kind_mod         ! Import rknd, iknd specifications
  
  Implicit none
  
  Real(rknd) :: ne, Te, dnedr, dTedr, beam_force
  Real(rknd), allocatable :: ni(:), Ti(:), dnidr(:), dTidr(:)
  Real(rknd), allocatable :: ni_prof(:,:), Ti_prof(:,:)
End module pprof_pass
!- End of header -------------------------------------------------------------
