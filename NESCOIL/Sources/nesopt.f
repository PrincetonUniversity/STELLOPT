! ----------------------------------------------------------------------
!     From command line, type:
!     xnescoil inputfile [ output goes to nescout.ext file ]
!
      program nesopt
      USE stel_kinds

      integer :: i
c......................
c     Single call to nescoil (no loop) with 0:
      call nescoil(0)
c..................................................

c     To use nescoil in a loop (say if embedded in an optimizer)

c     First call of nescoil (with 1):
c      This call does everything except deallocations
c      1. reads nesinput file,
c      2. allocates all arrays
c      3. performs one full nescoil calculation,
c      but does not deallocate anything

c     Intermediate call to nescoil (with > 1):
c      1. no allocations or deallocations
c      2. only calculates coil surface quantities
c      3. Does full gf and solver calculations

c     Last call to nescoil (with -1):
c     This one deallocates everything

c     Sample with 3 nescoil calls total:
c     do i = 1, 2        !your loop goes here
c         call nescoil (i)
c     enddo
c     call nescoil (-1)   !last call to nescoil

c     Strategy based on iloop:
c      1. allocate   only if 0 or  1 : single run or first call
c      2. deallocate only if 0 or -1 : single run or last call
c..................................................

      end program nesopt
