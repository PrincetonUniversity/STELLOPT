C  utility routines for (r8)pspltsub
 
      subroutine r8psp_tolsum(refval,tolval,sum)
 
C  for R4 sum.eq.refval expected, for R8 difference will be detectable
C  refval =~ 1.0; tolval =~1.0e-6
 
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 refval, tolval  ! input
      REAL*8 sum             ! output
 
      sum = refval + tolval*tolval
 
      return
      end
