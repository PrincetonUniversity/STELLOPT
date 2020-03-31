      SUBROUTINE point_wise_to_power(npts, x, f, n, ac, i_full,
     1                               nsp_pts)
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: npts, i_full, nsp_pts, n
      REAL(rprec), DIMENSION(npts), INTENT(IN) :: x, f
      REAL(rprec), DIMENSION(0:n), INTENT(OUT) :: ac
!-----------------------------------------------
      REAL, PARAMETER :: tol_res = 0.005_dp
      INTEGER :: i, j, k , n_count, n_int
      CHARACTER(len=2) :: cdum
      CHARACTER(len=15) :: form
      REAL(rprec) :: res, ohs_sp
      REAL(rprec), DIMENSION(:), ALLOCATABLE   :: tc, norm, y,
     1                                            fy, y_sp, fy_sp
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: a_inv, b_inv,
     1                                            a, b, inorm, tleg
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      REAL(rprec) , EXTERNAL :: integral
!--------------------------------------------------------------------------------
!      MEANING OF ARGUMENTS:
!      NPTS = number of points where the point-wise function is given
!      X = independent variable values for these points in [0,1]
!      F = dependent variable values on the X points
!      N = input number, 0:n, of Legendre polynomials used in the expansion
!      N_INT = actual computed number, 0:N_INT
!      AC = final power series coefficients in x ([0,1])
!      I_FULL =0, if input data on full mesh; =1, if on half mesh
!      NSP_PTS = number of points for the spline fit of F on the half mesh.
!--------------------------------------------------------------------------------
!
!      ABSTRACT:
!      Conversion of a point-wise solution in [0,1] given on a discrete set
!      of NPTS radial points (either on the half or the full mesh, see below)
!      first to an expansion in Legendre polynomials in [-1,1] via:
!
!          F = sum_m=0^N  <F,P_m>/<P_m,P_m>   <f,g> == int_[-1,1] f(x)g(x)dx
!
!      and defining the Legendre coefficient vector TC as:
!
!          TC = {<F,P_m>/<P_m,P_m>, m=1, N}
!
!      Then calculates the equivalent power series coefficients vector in
!      [0,1], AC:
!
!                AC = TC*A_inv*B_inv
!
!      The value of Legendre Polynomials used "N" depends CRITICALLY on the
!      number of points NPTS in the discrete data set of the incoming point-wise
!      solution. Because if the integrals are evaluated using these points, the
!      wiggles in the polynomial may not be well sampled by the existing set of
!      points and may cause deviations of the ortogonality condition:
!
!                        <P_m,P_n> = 2/(2*n+1)* delta_{mn}
!
!      This makes the reconstruction deteriorate quickly for increasing "N".
!      To avoid this, the incoming data values are splined to a number of points
!      NSP_PTS on a half-like mesh using Hermite cubic splines. Then, the
!      value of "N" is chosen requiring:
!
!               | 1- .5*(2*n+1)*<P_N,P_n> | < TOL_RES
!
!      TOL_RES is prescribed to be <= .5%. Until this condition is achieved, "N" is
!      subsequently reduced. If "N" gets equal to 2 in this reduction, the
!      calculation is stopped and an error message issued.
!
!      WARNING: The point-wise solution MAY BE given on the half mesh or full mesh
!      of a equally spaced radial grid, that is:
!
!      ***I_FULL = 1 beginning with x(3/2) and finishing with x(M-1/2), where
!         x(1) = 0. and x(M) = 1.
!      ***I_FULL = 0 beginning with x(1) = 0. and finishing with x(M) = 1.
!
!--------------------------------------------------------------------------------

      ac = 0
      res = 1

!...   Map [0,1] to [-1,1]

      ALLOCATE (y(npts), fy(npts), 
     1          y_sp(nsp_pts), fy_sp(nsp_pts), stat=i)
      IF (i .ne. 0) STOP 'Allocation error in point_wise_to_power!'

      y(:npts) = 2*x(:npts) - 1
      fy(:npts) = f(:npts)
      IF (i_full .eq. 1) fy(1)= fy(2)

!...   i_full = 0 IF FY given on the full mesh, and =1 IF on the half mesh.
!      FY is splined on nsp_pts on the half_mesh

      ohs_sp = 2._dp/nsp_pts
      y_sp= (/( -1 + (i -0.5_dp)*ohs_sp, i = 1,nsp_pts)/)

      CALL spline_it(npts, y, fy, nsp_pts, y_sp, fy_sp, i_full)

!...   Determine number of Legendre Polynomials to be used

      n_count = 0
      DO WHILE(res .gt. tol_res)
         n_int = n - n_count
         IF(n_int < n) DEALLOCATE (a, b, a_inv, b_inv, tleg,
     1                             inorm, norm)
         ALLOCATE(a(0:n_int,0:n_int), b(0:n_int,0:n_int),
     1            a_inv(0:n_int,0:n_int), b_inv(0:n_int,0:n_int),
     2            tleg(0:n_int, nsp_pts), inorm(0:n_int,0:n_int),
     3            norm(0:n_int), stat=i)
         IF (i .ne. 0) STOP 'Allocation error in point_wise_to_power!'

         norm(0:n_int) = (/(2._dp/(2*i+1), i=0, n_int)/)

         CALL build_matrices_legendre (n_int, a, b, a_inv, b_inv)

!...   Build Initial Legendre polynomials on grid

         DO i = 0, n_int
            tleg(i,:nsp_pts) = 0
            DO k = n_int, 0, -1
               tleg(i,:nsp_pts) = a_inv(i,k) + y_sp(:nsp_pts)
     1                          * tleg(i,:nsp_pts)
            END DO
         END DO

!...   Compute orthogonality and norms

         DO j = 0, n_int
            DO i = 0, n_int
               inorm(i,j) = integral(nsp_pts, y_sp, tleg(j,1:nsp_pts),
     1                               tleg(i,1:nsp_pts))
            END DO
            WRITE(cdum,'(i2)')n_int+1
            cdum = ADJUSTL(cdum)
            form = '('//TRIM(cdum)//'(x,f8.3))'
            form = ADJUSTL(form)
!           WRITE(20,TRIM(form)) (inorm(i,j), i=0, n)
         END DO

!...   Compute residual. If res > tol_res, reduce N and iterate while loop

         res = ABS((norm(n_int) - inorm(n_int,n_int))/(norm(n_int)))
         n_count = n_count + 1
         IF(n_count == n-1) STOP 'Precision bad for even 2 Polynomial'

      END DO

      ALLOCATE(tc(0:n_int))
      DO j = 0, n_int
         tc(j) = integral(nsp_pts, y_sp, fy_sp, tleg(j,1:nsp_pts))
     1              /norm(j)
      END DO

!
!     FIND ac COEFFICIENTS, 0:n_int; IF n_int < n, ac(n_int:n) = 0
!
      CALL legendre_to_power(n_int, a_inv, b_inv, tc, ac(0))

      DEALLOCATE(a_inv, b_inv, a, b, tc, tleg, inorm, norm, y,
     1           fy, y_sp, fy_sp)

      END SUBROUTINE point_wise_to_power
