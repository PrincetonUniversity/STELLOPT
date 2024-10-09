!-----------------------------------------------------------------------------
!+ Kind specifications for PENTA
!-----------------------------------------------------------------------------
Module penta_kind_mod
!
! Description:
!   This module contains the kind specifications for PENTA and all subroutines 
!   and modules.
!
! History:
! Version   Date      Comment
! -------   ----      -------
!  1.0     01/7/2009   Original Code.  JL
!  1.1     5/24/2010   Updated for PENTA3. JL 
! 
! Author(s): J. Lore 7/2009 - 5/24/2010    
!
  Implicit none
  Integer, parameter :: rknd = selected_real_kind(12,300) 
  Integer, parameter :: iknd = selected_int_kind(8)       
End module penta_kind_mod
!- End of header -------------------------------------------------------------
