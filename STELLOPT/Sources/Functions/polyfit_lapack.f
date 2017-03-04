!-----------------------------------------------------------------------
!     FUNCTION:       polyfit_lapack
!
!     PURPOSE:        This subroutine calculates a polynomial
!                     regression.  Code source was provided by the
!                     Rosetta Code project at:
!                     http://rosettacode.org/wiki/Polynomial_regression
!
!     INPUTS:         vx         X vector
!                     vy         Y vector
!                     npts       Number of points in X and Y vectors
!                     d          Order of polynomial
!
!     OUTPUTS:        a          Output polynomial array.
!
!     LIBRARIES:      lib_opt.a - kind_spec
!                     LAPACK    - DGETRF
!                               - DGETRI
!
!     WRITTEN BY:     Unknown
!                     ported by S. Lazerson
!
!     DATE:           07/06/11
!-----------------------------------------------------------------------
      subroutine polyfit_lapack(vx, vy, npts, d, a)
      use stel_kinds
      implicit none
      real(rprec), intent(in)    :: vx(npts), vy(npts)
      integer, intent(in)        :: d, npts
      real(rprec), intent(out)   :: a(d+1)
 
      real(rprec), dimension(:,:), allocatable :: X
      real(rprec), dimension(:,:), allocatable :: XT
      real(rprec), dimension(:,:), allocatable :: XTX
 
      integer :: i, j
 
      integer     :: n, lda, lwork
      integer :: info
      integer, dimension(:), allocatable :: ipiv
      real(rprec), dimension(:), allocatable :: work
 

      n = d+1
      lda = n
      lwork = n
      
      allocate(ipiv(n))
      allocate(work(lwork))
      allocate(XT(n, npts))
      allocate(X(npts, n))
      allocate(XTX(n, n))
 
      ! prepare the matrix
       do i = 0, d
         do j = 1, size(vx)
            X(j, i+1) = vx(j)**i
         end do
       end do
 
       XT  = transpose(X)
       XTX = matmul(XT, X)
 
       ! calls to LAPACK subs DGETRF and DGETRI
       call DGETRF(n, n, XTX, lda, ipiv, info)
       if ( info /= 0 ) then
          print *, "problem"
          return
       end if
       call DGETRI(n, XTX, lda, ipiv, work, lwork, info)
       if ( info /= 0 ) then
          print *, "problem"
          return
       end if
 
       a = matmul( matmul(XTX, XT), vy)
 
       deallocate(ipiv)
       deallocate(work)
       deallocate(X)
       deallocate(XT)
       deallocate(XTX)
 
       end subroutine polyfit_lapack
