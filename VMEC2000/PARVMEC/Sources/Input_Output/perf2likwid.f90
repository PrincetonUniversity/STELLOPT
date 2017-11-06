! perflib wrapper for likwid marker-API calls
! see http://code.google.com/p/likwid/wiki/LikwidPerfCtr#NEW:_Using_the_marker_API_with_Fortran_90

subroutine perfinit
!include "likwid_f90.h"
! call likwid_markerInit()
end subroutine perfinit

subroutine perfon(region)
!#include "likwid_f90.h"
  
  character(*), intent(in) :: region
!  call likwid_markerStart(region)
! call likwid_markerStartRegion(region,len_trim(region))
end subroutine perfon

subroutine perfoff(region)
!#include "likwid_f90.h"
  character(*), intent(in) :: region
!  call likwid_markerStop(region)
! call likwid_markerStopRegion(region,len_trim(region))
end subroutine perfoff

subroutine perfout(region)
!include "likwid_f90.h"
  character(*), intent(in) :: region
! call likwid_markerClose()
end subroutine perfout


subroutine perf_context_start()

end subroutine perf_context_start

subroutine perf_context_end()

end subroutine perf_context_end
