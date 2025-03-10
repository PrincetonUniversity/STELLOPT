
!-----------------------------------------------------------------------------
!+ Module for variables read from the DKES coeff. files
!-----------------------------------------------------------------------------
Module coeff_var_pass

!
! Description:
!   This module contains the DKES coeff. data variables.
!   For definitions see subroutine read_dkes_star_files 
!    in read_input_file_mod.f90.
!   For directions on how to generate the coeff. file see penta.f90.
!
! History:
! Version     Date      Comment
! -------   --------    -------
!  1.0     01/07/2009   Original Code.  JL
!  1.1     05/24/2010   Updated for PENTA3. JL 
!  1.2     09/26/2011   Consolidated cmul, efield arrays. JL
! 
! Author(s): J. Lore 7/2009 - 09/26/2011 
!     
  
  ! Modules used:
  Use penta_kind_mod         ! Import rknd, iknd specifications

  Implicit none

  ! Array variables
  Real(rknd), allocatable :: cmul(:)            ! cmul arrays for coeffs 
  Real(rknd), allocatable :: efield(:)          ! efield arrays for coeffs
  Real(rknd), allocatable :: D11_mat(:,:),    & ! 2D coefficient arrays
   D13_mat(:,:), D31_mat(:,:), D33_mat(:,:) 

  ! Scalar variables
  Integer(iknd) :: num_c, num_e                 ! Number of efield and cmul vals
End module coeff_var_pass
!- End of header -------------------------------------------------------------
