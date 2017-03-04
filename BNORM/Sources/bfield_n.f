!----------------------------------------------------------------------
!     DRIVER ROUTINE FOR BNORM CODE
!----------------------------------------------------------------------
      program bfield_n
      use stel_kinds
      use safe_open_mod
      use neswrite
      implicit none

!      integer, parameter :: nu=64, nv=64, mf=10, nf=10, md=20, nd=20
      integer, parameter :: ibnorm=15
      integer :: iunit = ibnorm, m, n, istat
      real(rprec), allocatable, dimension(:,:) :: bnfou, bnfou_c
!      real(rprec), dimension(0:mf,-nf:nf) :: bnfou
      character*(200) :: extension, separation
c-----------------------------------------------------------------
c     bnormal reads WOUT file (from VMEC) and delivers the
c     Fourier harmonics of B dot n :  bnfou, writing these out to the bnorm file
c-----------------------------------------------------------------
!
!     CALL BNORM AS FOLLOWS FROM COMMAND LINE:
!
!     xbnorm wout.extension coil_separation(optional)
!
!     For example:
!
!     xbnorm wout.extension 0.30
!
!     corresponds to the using the wout.ipp file for vmec input, and
!     a coil separation of 0.30 (in the correct units, of course).
!
!     This codes produces two files needed by NESCOIL (and possibly other codes):
!
!     bnorm.extension           normal b-field coefficients
!     nescin.extension          input file read for NESCOIL
!
!
!     NOTE: nescin has been modified from original format: np added to front end
!
      nu=256; nv=256; mf=24; nf=14; md=24; nd=20;
!      nu=64; nv=64; mf=12; nf=12; md=18; nd=18;

!
!     READ IN ARGUMENT FROM COMMAND LINE
!
      call bn_read_commandline(extension, separation)

      allocate(bnfou(0:mf,-nf:nf),
     1         bnfou_c(0:mf,-nf:nf),STAT=istat)
      if (istat .ne. 0) stop 'Allocation Error: bfield_n bnfou)'
      
      call bnormal(nu, nv, mf, nf, md, nd, bnfou, bnfou_c, extension)

!
!     WRITE OUT BNFOU (FOURIER COEFFICIENTS OF BNORMAL) TO BNORM.EXTENSION FILE
!
      call safe_open(iunit, istat, 'bnorm.' // extension, 'replace',
     1      'formatted')

      do m = 0, mf
         do n = -nf,nf
!            write(iunit, 1100) m, n, bnfou(m,n), bnfou_c(m,n)
            write(iunit, 1100) m, n, bnfou(m,n)
         end do
      end do

 1100 format(1x,2i5,1pe24.16,1pe24.16)

      close (iunit)


!
!     Write out input file read by NESCOIL (SPH: 9/98)
!
      call bn_write_nescoil_input(extension)
      
      deallocate(bnfou,bnfou_c)

!     print *, ' Finished writing bnorm and nescin files'

      end program bfield_n
