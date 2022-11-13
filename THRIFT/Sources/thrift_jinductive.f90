!-----------------------------------------------------------------------
!     Subroutine:    thrift_jinductive
!     Authors:       L. van Ham
!     Date:          11/14/2022
!     Description:   This subroutine updates the plasma response to
!                    source currents.  
!-----------------------------------------------------------------------
      SUBROUTINE thrift_jinductive
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE thrift_runtime
      USE thrift_vars
      USE stel_tools
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: i, ier
      REAL(rprec) :: rho,s,s11,s12,s21,s22,iota
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: Sup_iota, Sdn_iota
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
!     The general idea here is to use equation (22) in form to update
!     THRIFT_JPLASMA.  Where THRIFT_JSOURCE is the J in <J_s.B>.

      ! Allocations
      ALLOCATE(Sup_iota(nrho),Sdn_iota(nrho))

      ! First we calculate helpers
      DO i = 1, nrho
         rho = THRIFT_RHO(i)
         s   = rho*rho
         ier = 0
         CALL get_equil_sus(s,s11,s12,s21,s22,ier)
         CALL EZspline_interp(iota_spl,s,iota,ier)
         Sdn_iota(i) = (s21*iota+s22)
         Sup_iota(i) = (s11*iota+s12)
      END DO

      ! Deallocate
      DEALLOCATE(Sup_iota,Sdn_iota)

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_jinductive

