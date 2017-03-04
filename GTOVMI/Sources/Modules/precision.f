      MODULE precision
        USE stel_kinds, ONLY: rprec
!        INTEGER, PARAMETER :: rprec = SELECTED_REAL_KIND(12,100)
        INTEGER, PARAMETER :: iprec = SELECTED_INT_KIND(8)
        INTEGER, PARAMETER :: cprec = KIND((1.0_rprec,1.0_rprec))
        INTEGER, PARAMETER :: dbl = rprec
      END MODULE precision

