      SUBROUTINE fft_local(dat)
!  fft_local
!    Subroutine to interface with 'cftfax_g' and 'cfft99', Fourier transform
!    subroutines found in LIBSTELL/FFTpack.  'cftfax_g' conditions the vectors 'trigs' and
!    'ifax', while 'cfft99' does the actual transform.  See those subroutines for more
!    information.

!    Performs a 2D (two 1D) Fourier transform on the input data.  The data
!    should have the form dat(N1, N2, 2) where N1 and N2 are the number of
!    data points in the first and second direction, respectively, and
!    the final index separates the real (1) and imaginary (2) parts of the
!    data.

      USE stel_kinds, only: rprec
      
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      
      REAL(rprec), DIMENSION(:,:,:) :: dat
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: a, trigs, work
      INTEGER :: N, i, inc, jump, lot, isign, j
      INTEGER, DIMENSION(13) :: ifax
      
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

      N=SIZE(dat,1)
      isign=-1  !Sign of exponent in transform
      lot=SIZE(dat,2)
      jump=N  !Number to skip between data vectors.  Since we only do
              !one transform at a time, it shouldn't matter, but it
              !seems to have problems for values larger than N
      inc=1   !Number of array spaces between value pairs
      
      ALLOCATE(a(2*N))
      ALLOCATE(trigs(2*N))
      ALLOCATE(work(lot*N))
      
      DO j=1, lot
         DO i=1, N
!  The values are stored in a 1-d array with alternating real and imaginary
!  parts.
            a(2*i-1)=dat(i,j,1)
            a(2*i)=dat(i,j,2)
         END DO

!  Do the transform
         CALL cftfax_g(N, ifax, trigs)
         CALL cfft99(a,work,trigs, ifax, inc, jump, N, 1, isign)

!  Put the new values back into the input/output array
         DO i=1, N
            dat(i,j,1)=a(2*i-1)
            dat(i,j,2)=a(2*i)
         END DO
      END DO

!  Reallocate to do the transform in the other index
      IF(ALLOCATED(a)) DEALLOCATE(a)
      ALLOCATE(a(2*lot))
      IF(ALLOCATED(trigs)) DEALLOCATE(trigs)
      ALLOCATE(trigs(2*lot))
      IF(ALLOCATED(work)) DEALLOCATE(work)
      ALLOCATE(work(lot*N))

      DO i=1, N
         DO j=1, lot
            a(2*j-1)=dat(i,j,1)
            a(2*j)=dat(i,j,2)
         END DO

         CALL cftfax_g(lot, ifax, trigs)
         CALL cfft99(a,work,trigs, ifax, inc, jump, lot, 1, isign)

         DO j=1, lot
            dat(i,j,1)=a(2*j-1)
            dat(i,j,2)=a(2*j)
         END DO
      END DO
      
      END SUBROUTINE fft_local