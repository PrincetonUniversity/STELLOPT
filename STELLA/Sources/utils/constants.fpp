# include "define.inc"

module constants

!
! This module must not be compiled with a padding option
! such as -qautodbl=dbl of xlf which makes type conversion
! of variables even with explicit kind statements.
!
  implicit none

!  public :: size_of
  public :: kind_is, kind_id, kind_rs, kind_rd
  public :: zi, pi, twopi, dpi, dtwopi
# ifdef NAG_PREC
  public :: nag_kind
# endif

  private

! Symbolic names for kind type of single and double-precision reals:
! (with at least 6 and 12 digits of accuracy)

  integer, parameter :: kind_i1 = selected_int_kind (2)
  integer, parameter :: kind_ih = selected_int_kind (4)
  integer, parameter :: kind_is = selected_int_kind (8)
  integer, parameter :: kind_id = selected_int_kind (15)
  integer, parameter :: kind_rs = selected_real_kind (p=6)
  integer, parameter :: kind_rd = selected_real_kind (p=12)
  ! There is a selected_real_kind bug in xlf and the following does not work
  integer, parameter :: kind_rq = selected_real_kind (p=24)

!   integer, parameter :: sizeof_i1 = 1
!   integer, parameter :: sizeof_ih = 2
!   integer, parameter :: sizeof_is = 4
!   integer, parameter :: sizeof_id = 8
!   integer, parameter :: sizeof_rs = 4
!   integer, parameter :: sizeof_rd = 8
!   integer, parameter :: sizeof_rq = 16
!   integer, parameter :: sizeof_cs = 8
!   integer, parameter :: sizeof_cd = 16
!   integer, parameter :: sizeof_cq = 32

# if NAG_PREC == _NAGDBLE_
  integer, parameter :: nag_kind=kind_rd
# elif NAG_PREC == _NAGSNGL_
  integer, parameter :: nag_kind=kind_rs
# endif

! Symbolic names for kind type of single and double-precision complex:

!  integer, parameter :: spc = kind((1.0_sp,1.0_sp))
!  integer, parameter :: dpc = kind((1.0_dp,1.0_dp))

!  complex(dp), parameter :: ii = (0._dp, 1._dp)
!  real(dp), parameter :: pi=3.141592653589793238_dp

  complex, parameter :: zi = ( 0.0 , 1.0 )
!  real, parameter :: pi = 3.1415926535897931
!  real, parameter :: pi = 3.14159265358979323846, twopi=2.*pi
  ! this is actually quad precision
  double precision, parameter :: dpi = &
       3.14159265358979323846264338327950288419716939938, dtwopi=2.*dpi
  real, parameter :: pi = dpi, twopi= dtwopi

! Note: we will use dp="double precision" for almost everything.
!
! The fortran-90 "kind" types is kind of awkward.  But the old trick of
! using a "-r8" compiler switch to promote all real variables to 64 bits 
! does not work on some fortran 90 compilers, and so the above use of 
! the standard fortran-90 routine selected_real_kind is more portable.
!
! It may not be a good idea to mimic "-r8" by making sp to be identical
! to dp, or to write single and double-precision versions of 
! generic subroutines, since on the Cray computers both single and
! "double" precision are 64 bits, and the compiler will complain that
! it cannot distinguish the two specific subroutines.  In some cases,
! the cray compiler may be able to distinguish between two real "kinds"
! for the purposes of distinguishing overloaded procedure names,
! even though the two real kinds map to the same precision (64 bits).
!
! If this ever does become a problem, then you may be able to get around it by
! commenting out the double precision function names from the list of 
! overloaded procedures (i.e., the "module procedure" statements).
!

!   interface size_of
!      module procedure size_of_i1, size_of_ih, size_of_is, size_of_id
!      module procedure size_of_rs, size_of_rd
!      module procedure size_of_cs, size_of_cd
! !!$# ifdef QUAD
! !!$     module procedure size_of_rq, size_of_cq
! !!$# endif
!   end interface

contains

!   integer function size_of_i1 (arg)
!     integer (kind_i1) :: arg
!     size_of_i1 = sizeof_i1
!   end function size_of_i1

!   integer function size_of_ih (arg)
!     integer (kind_ih) :: arg
!     size_of_ih = sizeof_ih
!   end function size_of_ih

!   integer function size_of_is (arg)
!     integer (kind_is) :: arg
!     size_of_is = sizeof_is
!   end function size_of_is

!   integer function size_of_id (arg)
!     integer (kind_id) :: arg
!     size_of_id = sizeof_id
!   end function size_of_id

!   integer function size_of_rs (arg)
!     real (kind_rs) :: arg
!     size_of_rs = sizeof_rs
!   end function size_of_rs

!   integer function size_of_rd (arg)
!     real (kind_rd) :: arg
!     size_of_rd = sizeof_rd
!   end function size_of_rd

! !!$# ifdef QUAD
! !!$  integer function size_of_rq (arg)
! !!$    real (kind_rq) :: arg
! !!$    size_of_rq = sizeof_rq
! !!$  end function size_of_rq
! !!$# endif

!   integer function size_of_cs (arg)
!     complex (kind_rs) :: arg
!     size_of_cs = sizeof_cs
!   end function size_of_cs

!   integer function size_of_cd (arg)
!     complex (kind_rd) :: arg
!     size_of_cd = sizeof_cd
!   end function size_of_cd

!!$# ifdef QUAD
!!$  integer function size_of_cq (arg)
!!$    complex (kind_rq) :: arg
!!$    size_of_cq = sizeof_cq
!!$  end function size_of_cq
!!$# endif

end module constants
