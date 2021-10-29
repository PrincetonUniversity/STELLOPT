!!$ A C-program for MT19937: Real number version
!!$   genrand() generates one pseudorandom real number (double)
!!$ which is uniformly distributed on [0,1]-interval, for each
!!$ call. sgenrand(seed) set initial values to the working area
!!$ of 624 words. Before genrand(), sgenrand(seed) must be
!!$ called once. (seed is any 32-bit integer except for 0).
!!$ Integer generator is obtained by modifying two lines.
!!$   Coded by Takuji Nishimura, considering the suggestions by
!!$ Topher Cooper and Marc Rieffel in July-Aug. 1997.
!!$
!!$ This library is free software; you can redistribute it and/or
!!$ modify it under the terms of the GNU Library General Public
!!$ License as published by the Free Software Foundation; either
!!$ version 2 of the License, or (at your option) any later
!!$ version.
!!$ This library is distributed in the hope that it will be useful,
!!$ but WITHOUT ANY WARRANTY; without even the implied warranty of
!!$ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!!$ See the GNU Library General Public License for more details.
!!$ You should have received a copy of the GNU Library General
!!$ Public License along with this library; if not, write to the
!!$ Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
!!$ 02111-1307  USA
!!$
!!$ Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
!!$ When you use this, send an email to: matumoto@math.keio.ac.jp
!!$ with an appropriate reference to your work.
!!$
!!$***********************************************************************
!!$ Fortran translation by Hiroshi Takano.  Jan. 13, 1999.
!!$
!!$ This program uses the following non-standard intrinsics.
!!$   ishft(i,n): If n>0, shifts bits in i by n positions to left.
!!$               If n<0, shifts bits in i by n positions to right.
!!$   iand (i,j): Performs logical AND on corresponding bits of i and j.
!!$   ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
!!$   ieor (i,j): Performs exclusive OR on corresponding bits of i and j.
!!$
!!$************************************************************************
!!$ Fortran 95 translation and modularized to use in gyrokinetics project
!!$ by R. Numata, June 2, 2010.
!!$  * removed statement functions (TSHFTU,TSHFTS,TSHFTT,TSHFTL)
!!$  * uses [0,1) interval
!!$  * the above bit manipulation functions became standard.
!!$  * definition of UMASK is corrected (the original value cannot be
!!$    represented as a kind=4 integer.
!!$***********************************************************************

module mt19937
  
  implicit none

  private
  public :: sgrnd, grnd

  integer, parameter :: default_seed=4357

  ! Period parameters
  integer, parameter :: N = 624
  integer, parameter :: M = 397
  integer, parameter :: MATA = -1727483681 ! constant vector a
  integer, parameter :: LMASK = 2147483647 ! least significant r bits
  integer, parameter :: UMASK = -LMASK-1 ! most significant w-r bits

  ! Tempering parameters
  integer, parameter :: TMASKB= -1658038656
  integer, parameter :: TMASKC= -272236544

  integer, save :: mt(0:N-1) ! the array for the state vector
  integer, save :: mti=N+1 ! mti==N+1 means mt[N] is not initialized
  integer, save :: mag01(0:1)=(/ 0, MATA /) ! mag01(x) = x * MATA for x=0,1

contains

  subroutine sgrnd(seed)
    implicit none
    integer, intent(in) :: seed

!!$      setting initial seeds to mt[N] using
!!$      the generator Line 25 of Table 1 in
!!$      [KNUTH 1981, The Art of Computer Programming
!!$         Vol. 2 (2nd Ed.), pp102]

    mt(0)= iand(seed,-1)
    do mti=1,N-1
       mt(mti) = iand(69069 * mt(mti-1),-1)
    end do

    return
  end subroutine sgrnd

  function grnd ()
    implicit none
    real :: grnd
    real, parameter :: pow=4294967296.0 ! 2**32
    real, parameter :: div=1./pow       ! devided by 2**32 [0,1)-real-interval
    ! real, parameter :: div=1./(pow-1.)  ! devided by 2**32-1 [0,1]-real-interval
    integer :: y, kk

    if(mti >= N) then ! generate N words at one time
       if(mti == N+1) then ! if sgrnd() has not been called,
          call sgrnd(default_seed) ! a default initial seed is used
       endif

       do kk=0,N-M-1
          y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
          mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
       end do

       do kk=N-M,N-2
          y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
          mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
       end do

       y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
       mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
       mti = 0
    endif

    y=mt(mti)
    mti=mti+1
    y=ieor(y,ishft(y,-11))
    y=ieor(y,iand(ishft(y,7),TMASKB))
    y=ieor(y,iand(ishft(y,15),TMASKC))
    y=ieor(y,ishft(y,-18))

    grnd = real(y) 
    if(grnd < 0.) grnd=grnd + pow
    grnd = grnd * div

    return
  end function grnd

end module mt19937

!!$ this main outputs first 1000 generated numbers
!!$program main
!!$  use mt19937, only: grnd
!!$  ! use mt19937, only: sgrnd
!!$  implicit none
!!$  integer, parameter :: no=1000
!!$  real :: r(0:7)
!!$  integer :: j,k
!!$  ! call sgrnd(4357) ! any nonzero integer can be used as a seed
!!$
!!$  do j=0,no-1
!!$     r(mod(j,8))=grnd()
!!$     if(mod(j,8) == 7) then
!!$        write(*,'(8(f8.6,'' ''))') (r(k),k=0,7)
!!$     else if(j == no-1) then
!!$        write(*,'(8(f8.6,'' ''))') (r(k),k=0,mod(no-1,8))
!!$     endif
!!$  end do
!!$
!!$  stop
!!$end program main
