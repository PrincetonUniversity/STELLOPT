module convert
!------------------------------------------------------------------------------
!                              AstroGK, 2009
!------------------------------------------------------------------------------
! <doc>Convert from complex variable a(d1,d2,d3, ...) to a 
! real variable ar(2,d1,d2,d3,...) and back.  
! This is necessary for saving complex variables in NetCDF format</doc>
!
!     (c) Copyright 1991 to 1998 by Michael A. Beer, William D. Dorland, 
!     P. B. Snyder, Q. P. Liu, and Gregory W. Hammett. ALL RIGHTS RESERVED.
!     
  implicit none
  private
  public :: c2r, r2c

  interface c2r
     module procedure x1c2r
     module procedure x2c2r
     module procedure x3c2r
     module procedure x4c2r
     module procedure x5c2r
     module procedure x6c2r
  end interface

  interface r2c
     module procedure x1r2c
     module procedure x2r2c
     module procedure x3r2c
     module procedure x4r2c
     module procedure x5r2c
  end interface

contains
!------------------------------------------------------------------------------
!                              AstroGK, 2009
!------------------------------------------------------------------------------
! 
  subroutine x5c2r(a, a_ri)

    complex, dimension(:,:,:,:,:), intent(in) :: a
    real, dimension(:,:,:,:,:,:), intent(out) :: a_ri

    if(size(a, 1) /= size(a_ri, 2)) call aborter (6, 'x5c2r: size(a, 1) does not match size(a_ri, 2)')
    if(size(a, 2) /= size(a_ri, 3)) call aborter (6, 'x5c2r: size(a, 2) does not match size(a_ri, 3)')
    if(size(a, 3) /= size(a_ri, 4)) call aborter (6, 'x5c2r: size(a, 3) does not match size(a_ri, 4)')
    if(size(a, 4) /= size(a_ri, 5)) call aborter (6, 'x5c2r: size(a, 4) does not match size(a_ri, 5)')
    if(size(a, 5) /= size(a_ri, 6)) call aborter (6, 'x5c2r: size(a, 5) does not match size(a_ri, 6)')
    a_ri(1,:,:,:,:,:) = real(a(:,:,:,:,:))
    a_ri(2,:,:,:,:,:) = aimag(a(:,:,:,:,:))

  end subroutine x5c2r

  subroutine x6c2r(a, a_ri)

    complex, dimension(:,:,:,:,:,:), intent(in) :: a
    real, dimension(:,:,:,:,:,:,:), intent(out) :: a_ri

    if(size(a, 1) /= size(a_ri, 2)) call aborter (6, 'x6c2r: size(a, 1) does not match size(a_ri, 2)')
    if(size(a, 2) /= size(a_ri, 3)) call aborter (6, 'x6c2r: size(a, 2) does not match size(a_ri, 3)')
    if(size(a, 3) /= size(a_ri, 4)) call aborter (6, 'x6c2r: size(a, 3) does not match size(a_ri, 4)')
    if(size(a, 4) /= size(a_ri, 5)) call aborter (6, 'x6c2r: size(a, 4) does not match size(a_ri, 5)')
    if(size(a, 5) /= size(a_ri, 6)) call aborter (6, 'x6c2r: size(a, 5) does not match size(a_ri, 6)')
    if(size(a, 6) /= size(a_ri, 7)) call aborter (6, 'x6c2r: size(a, 6) does not match size(a_ri, 7)')
    a_ri(1,:,:,:,:,:,:) = real(a(:,:,:,:,:,:))
    a_ri(2,:,:,:,:,:,:) = aimag(a(:,:,:,:,:,:))

  end subroutine x6c2r
!------------------------------------------------------------------------------
!                              AstroGK, 2009
!------------------------------------------------------------------------------
! 
  subroutine x5r2c(a, a_ri)

    real, dimension(:,:,:,:,:,:), intent(in) :: a_ri
    complex, dimension(:,:,:,:,:), intent(out) :: a

    if(size(a, 1) /= size(a_ri, 2)) call aborter (6, 'x5r2c: size(a, 1) does not match size(a_ri, 2)')
    if(size(a, 2) /= size(a_ri, 3)) call aborter (6, 'x5r2c: size(a, 2) does not match size(a_ri, 3)')
    if(size(a, 3) /= size(a_ri, 4)) call aborter (6, 'x5r2c: size(a, 3) does not match size(a_ri, 4)')
    if(size(a, 4) /= size(a_ri, 5)) call aborter (6, 'x5r2c: size(a, 4) does not match size(a_ri, 5)')
    if(size(a, 5) /= size(a_ri, 6)) call aborter (6, 'x5r2c: size(a, 5) does not match size(a_ri, 6)')
    a(:,:,:,:,:) = cmplx(a_ri(1,:,:,:,:,:), a_ri(2,:,:,:,:,:))

  end subroutine x5r2c
!------------------------------------------------------------------------------
!                              AstroGK, 2009
!------------------------------------------------------------------------------
! 
  subroutine x4c2r(a, a_ri)

    complex, dimension(:,:,:,:), intent(in) :: a
    real, dimension(:,:,:,:,:), intent(out) :: a_ri

    if(size(a, 1) /= size(a_ri, 2)) call aborter (6, 'x4c2r: size(a, 1) does not match size(a_ri, 2)')
    if(size(a, 2) /= size(a_ri, 3)) call aborter (6, 'x4c2r: size(a, 2) does not match size(a_ri, 3)')
    if(size(a, 3) /= size(a_ri, 4)) call aborter (6, 'x4c2r: size(a, 3) does not match size(a_ri, 4)')
    if(size(a, 4) /= size(a_ri, 5)) call aborter (6, 'x4c2r: size(a, 4) does not match size(a_ri, 5)')
    a_ri(1,:,:,:,:) = real(a(:,:,:,:))
    a_ri(2,:,:,:,:) = aimag(a(:,:,:,:))

  end subroutine x4c2r
!------------------------------------------------------------------------------
!                              AstroGK, 2009
!------------------------------------------------------------------------------
! 
  subroutine x4r2c(a, a_ri)

    real, dimension(:,:,:,:,:), intent(in) :: a_ri
    complex, dimension(:,:,:,:), intent(out) :: a

    if(size(a, 1) /= size(a_ri, 2)) call aborter (6, 'x4r2c: size(a, 1) does not match size(a_ri, 2)')
    if(size(a, 2) /= size(a_ri, 3)) call aborter (6, 'x4r2c: size(a, 2) does not match size(a_ri, 3)')
    if(size(a, 3) /= size(a_ri, 4)) call aborter (6, 'x4r2c: size(a, 3) does not match size(a_ri, 4)')
    if(size(a, 4) /= size(a_ri, 5)) call aborter (6, 'x4r2c: size(a, 4) does not match size(a_ri, 5)')
    a(:,:,:,:) = cmplx(a_ri(1,:,:,:,:), a_ri(2,:,:,:,:))

  end subroutine x4r2c
!------------------------------------------------------------------------------
!                              AstroGK, 2009
!------------------------------------------------------------------------------
! 
  subroutine x3c2r(a, a_ri)

    complex, dimension(:,:,:), intent(in) :: a
    real, dimension(:,:,:,:), intent(out) :: a_ri

    if(size(a, 1) /= size(a_ri, 2)) call aborter (6, 'x3c2r: size(a, 1) does not match size(a_ri, 2)')
    if(size(a, 2) /= size(a_ri, 3)) call aborter (6, 'x3c2r: size(a, 2) does not match size(a_ri, 3)')
    if(size(a, 3) /= size(a_ri, 4)) call aborter (6, 'x3c2r: size(a, 3) does not match size(a_ri, 4)')
    a_ri(1,:,:,:) = real(a(:,:,:))
    a_ri(2,:,:,:) = aimag(a(:,:,:))

  end subroutine x3c2r
!------------------------------------------------------------------------------
!                              AstroGK, 2009
!------------------------------------------------------------------------------
! 
  subroutine x3r2c(a, a_ri)

    real, dimension(:,:,:,:), intent(in) :: a_ri
    complex, dimension(:,:,:), intent(out) :: a

    if(size(a, 1) /= size(a_ri, 2)) call aborter (6, 'x3r2c: size(a, 1) does not match size(a_ri, 2)')
    if(size(a, 2) /= size(a_ri, 3)) call aborter (6, 'x3r2c: size(a, 2) does not match size(a_ri, 3)')
    if(size(a, 3) /= size(a_ri, 4)) call aborter (6, 'x3r2c: size(a, 3) does not match size(a_ri, 4)')
    a(:,:,:) = cmplx(a_ri(1,:,:,:), a_ri(2,:,:,:))

  end subroutine x3r2c
!------------------------------------------------------------------------------
!                              AstroGK, 2009
!------------------------------------------------------------------------------
! 
  subroutine x2c2r(a, a_ri)

    complex, dimension(:,:), intent(in) :: a
    real, dimension(:,:,:), intent(out) :: a_ri

    if(size(a, 1) /= size(a_ri, 2)) call aborter (6, 'x2c2r: size(a, 1) does not match size(a_ri, 2)')
    if(size(a, 2) /= size(a_ri, 3)) call aborter (6, 'x2c2r: size(a, 2) does not match size(a_ri, 3)')
    a_ri(1,:,:) = real(a(:,:))
    a_ri(2,:,:) = aimag(a(:,:))

  end subroutine x2c2r
!------------------------------------------------------------------------------
!                              AstroGK, 2009
!------------------------------------------------------------------------------
! 
  subroutine x2r2c(a, a_ri)

    real, dimension(:,:,:), intent(in) :: a_ri
    complex, dimension(:,:), intent(out) :: a

    if(size(a, 1) /= size(a_ri, 2)) call aborter (6, 'x2r2c: size(a, 1) does not match size(a_ri, 2)')
    if(size(a, 2) /= size(a_ri, 3)) call aborter (6, 'x2r2c: size(a, 2) does not match size(a_ri, 3)')
    a(:,:) = cmplx(a_ri(1,:,:), a_ri(2,:,:))

  end subroutine x2r2c
!------------------------------------------------------------------------------
!                              AstroGK, 2009
!------------------------------------------------------------------------------
! 
  subroutine x1c2r(a, a_ri)

    complex, dimension(:), intent(in) :: a
    real, dimension(:,:), intent(out) :: a_ri

    if(size(a, 1) /= size(a_ri, 2)) call aborter (6, 'x2c2r: size(a, 1) does not match size(a_ri, 2)')
    a_ri(1,:) = real(a(:))
    a_ri(2,:) = aimag(a(:))

  end subroutine x1c2r
!------------------------------------------------------------------------------
!                              AstroGK, 2009
!------------------------------------------------------------------------------
! 
  subroutine x1r2c(a, a_ri)

    real, dimension(:,:), intent(in) :: a_ri
    complex, dimension(:), intent(out) :: a

    if(size(a, 1) /= size(a_ri, 2)) call aborter (6, 'x2r2c: size(a, 1) does not match size(a_ri, 2)')
    a(:) = cmplx(a_ri(1,:), a_ri(2,:))

  end subroutine x1r2c
!------------------------------------------------------------------------------
!                              AstroGK, 2009
!------------------------------------------------------------------------------
! 
  Subroutine Aborter(iunit,ierrmsg)
!----------------------------------------------------------------------
!  ABORT A PROGRAM AFTER A FATAL ERROR CONDITION IS DETECTED.
!
!
! input: iunit  Unit Number of the file to write error messages to.
!        ierrmsg  An error message to write ilunerr
!
!       The advantage of using this subroutine is that it will
!       generate a traceback showing the chain of subroutines which
!       eventually bombed, and it forces an arithmetic error which
!       the job control system can detect as a fatal error.
!
    character ierrmsg*(*)
    real :: zz0,zz1
    integer :: iunit,ilunerr
    common /abortcmn/ zz0,zz1
!
! zz0 is in a common block to prevent an optimizing compiler
! from evaluating 1.0/zz0 during compilation rather than during
! execution.
!
    write(iunit,1001)
1001 format(//' %ABORTER:  ** FATAL ERROR.  ABORT SUBROUTINE CALLED **'//)
  
    write(iunit,1002)ierrmsg
1002 format(1x,a,//)

! on CRAY's and VAXes:
! generate a divide-by-zero error:
    zz1=1.0/zz0         !^
    ilunerr=5           !^

! on the DecStation:
!^      call exit(1)

    write(ilunerr,1010) zz0
1010 format(' ?ABORTER-- ZZ0= ',1PE11.4,' AND STILL EXECUTING...')

  end subroutine aborter
!------------------------------------------------------------------------------
end module convert
