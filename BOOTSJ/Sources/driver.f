      program driver
      use stel_kinds
      use parambs, only: lscreen
      use safe_open_mod
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat, i, numargs, iunit_in = 10
      real(rprec) :: t1, t2
      character*120 :: extension
      character*120 :: arg1, arg2
      real(rprec) :: curtor
C-----------------------------------------------
!
!     driver: reads from command line file the wout file extension and surfaces
!     (half-radial) at which the bootstrap current is  required
!     writes the bootstrap current, JdotB to a file, jBbs.extension
!
!     call this as follows:   xbootjs input.boots (T or F)
!
!     where input.boz contains the wout file extension and the jrad values (as a
!     blank-delimited list, not necessarily on a single line):
!
!     FILE_EXTENSION
!     2  3   5   10  12
!
!     The surface numbers are relative to half-mesh vmec quantities.
!     thus, 2 is the first half mesh point.
!
!     The optional (T) or (F) argument allows (default) or suppresses output to screen.
!
      lscreen = .true.

      call getcarg(1, arg1, numargs)
      if (numargs .gt. 1) call getcarg(2, arg2, numargs)

      if (numargs .lt. 1 .or.
     1   (arg1 .eq. '-h' .or. arg1 .eq. '/h')) then
         print *,
     1   ' ENTER INPUT FILE NAME ON COMMAND LINE'
         print *,' For example: xbootsj in_bootsj.tftr'
         print *
         print *,
     1   ' Optional command line argument to suppress screen output:'
         print *,' xbootsj input_file <(T or F)>'
         print *
         print *,' where F will suppress screen output'
         stop
      else if (numargs .gt. 1) then
         if (arg2(1:1).eq.'f' .or. arg2(1:1).eq.'F') lscreen = .false.
      endif

!      INPUT FILE
!      1st line:   extension WOUT file
!      2nd line:   surfaces to be computed

      call safe_open(iunit_in, istat, trim(arg1), 'old', 'formatted')
      if (istat .ne. 0) stop 'error opening input file in bootsj'

      read (iunit_in, *,iostat=istat) extension

      call bootsj(curtor, trim(extension), iunit_in)

      end program driver
