!-----------------------------------------------------------------------------
!+ Module for physical constants used in several routines
!-----------------------------------------------------------------------------
Module phys_const

!
! Description:
!   This module contains physical constant parameters.
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
  Use penta_kind_mod         ! Import rknd, iknd specifications

  Implicit none

  ! Scalar parameters
  Real(rknd),parameter :: p_mass = 1.672621637e-27_rknd      !proton mass
  Real(rknd),parameter :: e_mass = 9.10938215e-31_rknd       !electron mass
  Real(rknd),parameter :: elem_charge = 1.602176487e-19_rknd !elem. charge 
  Real(rknd),parameter :: eps0 = 8.854187817e-12_rknd        !electric const
  Real(rknd),parameter :: pi = 3.14159265358979323846_rknd   !pi
!  Real(rknd),parameter :: pi = 4._rknd*Datan(1._rknd)
  Real(rknd),parameter :: r_eps = Epsilon(1._rknd)           ! Machine prec.
End module phys_const

!- End of header -------------------------------------------------------------
