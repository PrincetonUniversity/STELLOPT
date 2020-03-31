C  utility routines for (r8)pspltsub

      subroutine psp_tolsum(refval,tolval,sum)

C  for R4 sum.eq.refval expected, for R8 difference will be detectable
C  refval =~ 1.0; tolval =~1.0e-6

      real refval, tolval  ! input
      real sum             ! output

      sum = refval + tolval*tolval

      return
      end
