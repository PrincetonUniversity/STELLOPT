
!----------------------------------------------------------------------
!     NESCOIL ROUTINES
!----------------------------------------------------------------------
      subroutine bn_read_commandline(extension, separation)
      use neswrite, only: coil_separation, dp, nu, nv, mf, nf, md, nd
      implicit none
      character(len=*) :: extension, separation
      integer :: numargs, ia, ib, istat
c ----------------------------------------------------------------------
c                                                            11.09.99
c     purpose:  read command line arguments (SPH)
c
c
c ----------------------------------------------------------------------

      call getcarg(1, extension, numargs)
      if (numargs .lt.1 .or. extension .eq. '-h' .or.
     1                       extension .eq. '/h') then
         print *,
     1   ' Syntax:  xbnorm wout.extension coil_separation(optional)'
         print *, ' '
         print *,
     3   ' where the coil separation (in same units as R00) is optional'
         stop 'Invalid command line'
      else if (numargs .eq. 2) then
         call getcarg(2, separation, numargs)
         coil_separation = 0
         read (separation, *, iostat=istat) coil_separation
         if (istat .ne. 0) stop 'BNORM Error reading coil_separation'
      else if (numargs .gt. 2) then
         print *,'please enter customized dimensions: nu,nv,mf,nf,md,nd'
         READ (*,*) nu
         READ (*,*) nv
         READ (*,*) mf
         READ (*,*) nf
         READ (*,*) md
         READ (*,*) nd
         call getcarg(2, separation, numargs)
         coil_separation = 0
         read (separation, *, iostat=istat) coil_separation
         if (istat .ne. 0) stop 'BNORM Error reading coil_separation'
      endif

      ia = index(extension, 'wout_')
      ib = index(extension, 'wout.')
      if (ia .gt. 0) then
         extension = extension(ia+5:len_trim(extension))
      else if (ib .gt. 0) then
         extension = extension(ib+5:len_trim(extension))
      end if

      ia = index(extension, '.nc')
      if (ia .gt. 0) extension = extension(1:ia-1)


      end subroutine bn_read_commandline
