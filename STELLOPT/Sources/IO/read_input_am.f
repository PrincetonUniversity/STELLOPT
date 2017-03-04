!-----------------------------------------------------------------------
!     FUNCTION:       read_input_am
!
!     PURPOSE:        This subroutine reads the am array from an input
!                     file.
!
!     INPUTS:         filename   Name of file to read.
!
!     OUTPUTS:        istat      Return status
!                     am_out     AM array from input file.
!
!     LIBRARIES:      lib_opt.a - kind_spec
!                               - read_namelist_mod
!                               - safe_open_mod
!
!     WRITTEN BY:     S. Lazerson (lazerson@pppl.gov)
!
!     DATE:           07/06/11
!-----------------------------------------------------------------------
      subroutine read_input_am(filename,istat,am_out)
      use stel_kinds
      use safe_open_mod
      use vmec_input, ONLY: am, read_indata_namelist
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      character*(*), intent(in) :: filename
      integer, intent(out)      :: istat
      real(rprec), intent(out)  :: am_out(0:20)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer                   :: iunit=7
      real(rprec)               :: am_old(0:20)
      
      istat=0
      am_old=am
      am_out(:)=0.0
      
      call safe_open(iunit,istat,trim(filename),'old','formatted')
      ! Note that if a min file doesn't exist we just default am back to what it was.
      if (istat .ne. 0) then
         print *,'   safe_open error: ',trim(filename)
         print *,'   istat = ',istat
         close(iunit)
         am_out=am_old
         return
      end if
      call read_indata_namelist(iunit,istat)
      if (istat .ne. 0) then
         print *,'   indata namelist read error: ',trim(filename)
         print *,'   istat = ',istat
         am_out=am_old
         return
      end if
      close(iunit)
      am_out=am
      return
 
      end subroutine read_input_am
