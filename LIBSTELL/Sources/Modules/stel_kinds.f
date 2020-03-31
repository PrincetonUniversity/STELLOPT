      MODULE stel_kinds

!----------------------------------------------------------------------
!  Kind specifications
!----------------------------------------------------------------------

      INTEGER, PARAMETER :: rprec = SELECTED_REAL_KIND(15,300)
      INTEGER, PARAMETER :: iprec = SELECTED_INT_KIND(8)
      INTEGER, PARAMETER :: cprec = KIND((1.0_rprec,1.0_rprec))
      INTEGER, PARAMETER :: dp = rprec

      END MODULE stel_kinds
