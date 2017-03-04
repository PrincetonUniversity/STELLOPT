!-----------------------------------------------------------------------
!     Module:        stelltran_equilutils
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          08/30/2015
!     Description:   This module contains arrays and subroutines for
!                    handling the equilibrium in a portable format.
!-----------------------------------------------------------------------
      MODULE stelltran_equilutils
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE stelltran_runtime, ONLY: pi2
      USE stel_tools
      USE EZspline
      USE EZspline_obj
!-----------------------------------------------------------------------
!     Module Variables
!         nrad        Number of radial grid points
!         nfp         Field periodicity of the equilibirum
!         lasym       Non-stellarator symmetric logical
!         R_spl       Spline of the equilibirum R values
!         Z_spl       Spline of the equilibrium Z values
!         Pow_e_spl   Spline of electron power profile
!         Pow_i_spl   Spline of ion power profile
!         Part_spl    Spline of particle influx profile
!         deriv_spl   Spline used to evaluate derivatives
!         grid_spl    Spline used to interpolate points onto various grids
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, PARAMETER :: nsys   = 16

      LOGICAL     ::  lasym
      INTEGER     ::  domain_flag
      INTEGER     ::  nrad, nfp
      INTEGER     :: npopulation, nra_ecrh, nphi_ecrh
      REAL(rprec) :: R_target, PHI_target, Z_target
      REAL(rprec) :: R0, aspect, Baxis, Aminor
      REAL(rprec), ALLOCATABLE :: rho(:), iotaf(:)
      REAL(rprec), DIMENSION(nsys) :: freq_ecrh
      REAL(rprec), DIMENSION(nsys,3)     :: antennaposition_ecrh, targetposition_ecrh,rbeam_ecrh,rfocus_ecrh
      CHARACTER, DIMENSION(nsys,5)     :: wmode_ecrh
      TYPE(EZspline1_r8) :: prof_spl, Pow_e_spl, Pow_i_spl, Part_sple, Er_spl
      
      INTEGER     ::  bcs0(2) = (/ 0, 0/)
      INTEGER     ::  bcs1(2) = (/-1,-1/)
      
!-----------------------------------------------------------------------
!     Subroutines
!         eval_prof_spline  Spline and evaluate values
!         mntouv            Transforms to real space
!         rzfunct           Returns R, Z given s,u,v (for get_equil_s)
!         get_equil_s       Returns s,u,v given R,phi,Z
!         get_equil_RZ      Returns R, Z given s,u,v
!         sort_arrays       Sorts array s and f by s in ascending order
!         get_equil_te      Returns te at s
!         get_equil_ti      Returns ti at s
!         get_equil_ne      Returns ne at s
!         get_equil_zeff    Returns zeff at s
!         *_rhs             Evaluates the right hand side of a balance equation for COLNEW
!         d*_rhs            Jacobian of right hand side of a balance equation for COLNEW
!         *_bound_cond      Boundary conditions of a balance equation for COLNEW
!         d*_bound_cond     Jacobian of boundary condition of a balance equation for COLNEW
!         fin_deriv4        Fourth order finite difference derivative
!-----------------------------------------------------------------------

      CONTAINS 
      
      SUBROUTINE eval_prof_spline(nsp,xknots,yvals,xval,fval,ier,fpval)
      IMPLICIT NONE
      INTEGER, INTENT(in)        :: nsp
      REAL(rprec), INTENT(in)    :: xknots(nsp), yvals(nsp)
      REAL(rprec), INTENT(in)    :: xval
      REAL(rprec), INTENT(out)   :: fval
      REAL(rprec), INTENT(out), OPTIONAL   :: fpval
      INTEGER, INTENT(inout)     :: ier
      INTEGER, SAVE :: nsp_old
      REAL(rprec),ALLOCATABLE, SAVE :: xknots_old(:), yvals_old(:)
      fval = 0
      IF (ier < 0) RETURN
      IF (.not. ALLOCATED(xknots_old)) THEN
         ALLOCATE(xknots_old(nsp))
         xknots_old(1:nsp) = xknots(1:nsp)
      END IF
      IF (.not. ALLOCATED(yvals_old)) THEN
         ALLOCATE(yvals_old(nsp))
         yvals_old(1:nsp) = yvals(1:nsp)
      END IF
      IF (nsp_old /= nsp .or. ANY(xknots .ne. xknots_old) .or. ANY(yvals .ne. yvals_old)) THEN
         nsp_old = nsp
         DEALLOCATE(xknots_old); ALLOCATE(xknots_old(nsp))
         DEALLOCATE(yvals_old); ALLOCATE(yvals_old(nsp))
         xknots_old(1:nsp) = xknots(1:nsp)
         yvals_old(1:nsp)  = yvals(1:nsp)
         IF (EZspline_allocated(prof_spl)) CALL EZspline_free(prof_spl,ier)
         CALL EZspline_init(prof_spl,nsp,bcs0,ier)
         IF (ier .ne. 0) RETURN
         prof_spl%x1 = xknots(1:nsp)
         prof_spl%isHermite = 1
         CALL EZspline_setup(prof_spl,yvals(1:nsp),ier)
         IF (ier .ne. 0) RETURN
         CALL EZspline_isInDomain(prof_spl,xval,ier)
         IF (ier .ne. 0) RETURN
         CALL EZspline_interp(prof_spl,xval,fval,ier)
      ELSE
         CALL EZspline_isInDomain(prof_spl,xval,ier)
         IF (ier .ne. 0) RETURN
         CALL EZspline_interp(prof_spl,xval,fval,ier)
      END IF
      IF (PRESENT(fpval)) THEN
         CALL EZspline_derivative(prof_spl,1,xval,fpval,ier)
         IF (ier .ne. 0) RETURN
      END IF
      RETURN
      END SUBROUTINE eval_prof_spline
      
      SUBROUTINE sort_array(k,sarr,farr)
      IMPLICIT NONE
      INTEGER, INTENT(in)        ::  k
      REAL(rprec), INTENT(inout)   ::  sarr(k)
      REAL(rprec), INTENT(inout)   ::  farr(k)
      INTEGER :: ik, i1
      LOGICAL :: lmask(k)
      REAL(rprec) :: tval
      lmask = .FALSE.
      DO ik = 1, k
         i1 = MINLOC(sarr,DIM=1,MASK=lmask)
         IF (i1 == 0) CYCLE
         tval = sarr(i1)
         sarr(i1) = sarr(ik)
         sarr(ik) = tval
         tval = farr(i1)
         farr(i1) = farr(ik)
         farr(ik) = tval
         lmask(ik) = .TRUE.
      END DO
      RETURN
      END SUBROUTINE sort_array
      
      function polyfit(vx, vy, d)
      !  From http://rosettacode.org/wiki/Polynomial_regression
      implicit none
      integer, intent(in)                   :: d
      integer, parameter                    :: dp = selected_real_kind(15, 307)
      real(dp), dimension(d+1)              :: polyfit
      real(dp), dimension(:), intent(in)    :: vx, vy
      
      real(dp), dimension(:,:), allocatable :: X
      real(dp), dimension(:,:), allocatable :: XT
      real(dp), dimension(:,:), allocatable :: XTX
      
      integer :: i, j
      
      integer     :: n, lda, lwork
      integer :: info
      integer, dimension(:), allocatable :: ipiv
      real(dp), dimension(:), allocatable :: work
      
      n = d+1
      lda = n
      lwork = n
      
      allocate(ipiv(n))
      allocate(work(lwork))
      allocate(XT(n, size(vx)))
      allocate(X(size(vx), n))
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
      
      polyfit = matmul( matmul(XTX, XT), vy)
      
      deallocate(ipiv)
      deallocate(work)
      deallocate(X)
      deallocate(XT)
      deallocate(XTX)
      RETURN
      
      end function polyfit

      SUBROUTINE smoothg(A,N,DLZ,B)
      ! Smoothing of A with gaussian width DLZ (H. Mynick 12/06/09)
      IMPLICIT NONE
      REAL(rprec), INTENT(IN) :: A(N)
      REAL(rprec), INTENT(IN) :: DLZ
      INTEGER, INTENT(IN) :: N
      REAL(rprec), INTENT(OUT) :: B(N)
      INTEGER :: i, j
      REAL(rprec) :: arg,gaussn, dlth, gsum
      dlth=pi2/(N-1)  !=th incrmt btwn successive meshpts.
      B=0.0_rprec 
      DO i=1,N
         gsum=0.0_rprec
         DO j=1,N
            arg=(i-j)*dlth/DLZ
            gaussn=EXP(-arg*arg)
            B(i)=B(i)+gaussn*A(j)
            gsum=gsum+gaussn
         END DO                  !j
         B(i)=B(i)/gsum
      END DO                     !i
      RETURN
      END SUBROUTINE smoothg
      
      function polyval(c,s,d)
      integer, intent(in)                   :: d
      integer, parameter                    :: dp = selected_real_kind(15, 307)
      real(dp), dimension(d), intent(in)  :: c
      real(dp), intent(in)    :: s
      real(dp)                              :: polyval
      INTEGER :: j
      polyval = 0
      DO j = d, 1, -1
         polyval = polyval*s + c(j)
      END DO
      RETURN
      end function polyval

      SUBROUTINE fin_deriv2(f, dr, qkind, df)
            ! second order finite difference derivative
            IMPLICIT NONE
            CHARACTER, INTENT(IN) :: qkind ! Left/right endpoint or midpoint
            REAL(rprec), INTENT(IN) :: f(3), dr ! input points
            REAL(rprec), INTENT(OUT) :: df ! output derivative
            INTEGER :: i
            SELECT CASE (qkind)
                  CASE('l')
                        df = (-3.*f(1)/2. + 2.*f(2) - f(3)/2.)/dr
                  CASE('m')
                        df = (f(3) - f(1))/(dr*2.)
                  CASE('r')
                        df = (3.*f(3)/2. - 2.*f(2) + f(1)/2.)/dr
            END SELECT
      END SUBROUTINE fin_deriv2

      SUBROUTINE arr_fderiv2(npts, f, dr, df)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: npts
            REAL(rprec), INTENT(IN) :: f(npts), dr
            REAL(rprec), INTENT(OUT) :: df(npts)
            INTEGER :: i
            CALL fin_deriv2(f(1:3), dr, 'l', df(1))
            DO i=2,npts-1
                  CALL fin_deriv2(f(i-1:i+1),dr,'m',df(i))
            END DO
            CALL fin_deriv2(f(npts-2:npts),dr,'r',df(npts))
      END SUBROUTINE arr_fderiv2

      SUBROUTINE fin_2deriv2(f, dr, qkind, df)
            ! second order finite difference second derivative
            IMPLICIT NONE
            CHARACTER, INTENT(IN) :: qkind ! Left/right endpoint or midpoint
            REAL(rprec), INTENT(IN) :: f(4), dr ! input points
            REAL(rprec), INTENT(OUT) :: df ! output derivative
            INTEGER :: i
            SELECT CASE (qkind)
                  CASE('l')
                        df = (2.*f(1) - 5.*f(2) + 4.*f(3) - f(4))/(dr**2.)
                  CASE('m')
                        df = (f(3) - 2.*f(2) + f(1))/(dr**2.)
                  CASE('r')
                        df = (2.*f(4) - 5.*f(3) + 4.*f(2) - f(1))/(dr**2.)
            END SELECT
      END SUBROUTINE fin_2deriv2

      SUBROUTINE arr_2fderiv2(npts, f, dr, df)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: npts
            REAL(rprec), INTENT(IN) :: f(npts), dr
            REAL(rprec), INTENT(OUT) :: df(npts)
            INTEGER :: i
            CALL fin_2deriv2(f(1:4), dr, 'l', df(1))
            DO i=2,npts-1
                  CALL fin_2deriv2(f(i-1:i+1),dr,'m',df(i))
            END DO
            CALL fin_2deriv2(f(npts-3:npts),dr,'r',df(npts))
      END SUBROUTINE arr_2fderiv2

      SUBROUTINE fin_deriv4(f, dr, qkind, df)
            ! fourth order finite difference derivative
            IMPLICIT NONE
            CHARACTER, INTENT(IN) :: qkind ! Left/right endpoint or midpoint
            REAL(rprec), INTENT(IN) :: f(5), dr ! input points
            REAL(rprec), INTENT(OUT) :: df ! output derivative
            INTEGER :: i
            SELECT CASE (qkind)
                  CASE('l')
                        df = (-25.*f(1)/12. + 4.*f(2) - 3.*f(3) + 4.*f(4)/3. - f(5)/4.)/dr
                  CASE('m')
                        df = (f(1)/12. - 2.*f(2)/3. + 2.*f(4)/3. - f(5)/12.)/dr
                  CASE('r')
                        df = (25.*f(5)/12. - 4.*f(4) + 3.*f(3) - 4.*f(2)/3. + f(1)/4.)/dr
            END SELECT
      END SUBROUTINE fin_deriv4

      SUBROUTINE arr_fderiv4(npts, f, dr, df)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: npts
            REAL(rprec), INTENT(IN) :: f(npts), dr
            REAL(rprec), INTENT(OUT) :: df(npts)
            INTEGER :: i
            CALL fin_deriv4(f(1:5), dr, 'l', df(1))
            CALL fin_deriv4(f(2:6), dr, 'l', df(2))
            DO i=3,npts-2
                  CALL fin_deriv4(f(i-2:i+2),dr,'m',df(i))
            END DO
            CALL fin_deriv4(f(npts-5:npts-1),dr,'r',df(npts-1))
            CALL fin_deriv4(f(npts-4:npts),dr,'r',df(npts))
      END SUBROUTINE arr_fderiv4

      SUBROUTINE bal_rhs(x,z,f)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: z(3)
            REAL*8, INTENT(IN) :: x
            REAL*8, INTENT(OUT) :: f(3)
            REAL*8 :: fval_powe, fval_er
            INTEGER :: ier

            ! particle blance equation part
            CALL EZspline_interp(part_sple,x,f(1),ier)

            ! Electron power balance equation part
            CALL EZspline_interp(pow_e_spl,x,fval_powe,ier)
            CALL EZspline_interp(Er_spl,x,fval_er,ier)
            f(2) = fval_powe - z(1)*fval_er

            ! Ion power balance equation part. Should eventually include z*Gi*Er part.
            CALL EZspline_interp(pow_i_spl,x,f(3),ier)
      END SUBROUTINE bal_rhs

      SUBROUTINE dbal_rhs(x,z,df)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: z(3)
            REAL*8, INTENT(IN) :: x
            REAL*8, INTENT(OUT) :: df(3,3)
            INTEGER :: ier

            df = 0.D0
            CALL EZspline_interp(Er_spl,x,df(2,1),ier) ! e pow bal depends on Ge
      END SUBROUTINE dbal_rhs

      SUBROUTINE bal_bound_cond(i,z,g)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: z(3)
            INTEGER, INTENT(IN) :: i
            REAL*8, INTENT(OUT) :: g

            g = z(i)
      END SUBROUTINE bal_bound_cond

      SUBROUTINE dbal_bound_cond(i,z,dg)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: i
            REAL*8, INTENT(IN) :: z(3)
            REAL*8, INTENT(OUT) :: dg(3)

            dg = 0.D0
            dg(i) = 1.D0
      END SUBROUTINE dbal_bound_cond

!       SUBROUTINE pow_e_rhs(x,z,f)
!             IMPLICIT NONE
!             REAL*8, INTENT(IN) :: z(1)
!             REAL(rprec), INTENT(IN) :: x
!             REAL*8, INTENT(OUT) :: f(1)
!             REAL*8 :: func_val
!             INTEGER :: ier
!             call EZspline_interp(Pow_e_spl,x,func_val,ier)
!             f(1) = func_val
!       END SUBROUTINE pow_e_rhs
! 
!       SUBROUTINE dpow_e_rhs(x,z,df)
!             IMPLICIT NONE
!             REAL*8, INTENT(IN) :: z(1)
!             REAL(rprec), INTENT(IN) :: x
!             REAL*8, INTENT(OUT) :: df(1,1) ! Jacobian of f wrt "y"
!             df(1,1) = 0.D0 ! We assume it does not depend on "y"
!       END SUBROUTINE dpow_e_rhs
! 
!       SUBROUTINE powe_bound_cond(i,z,g)
!             REAL*8, INTENT(IN) :: z(1)
!             INTEGER, INTENT(IN) :: i
!             REAL*8, INTENT(OUT) :: g
!             g = z(1) ! The boundary condition, basically what should be zero.
!       END SUBROUTINE powe_bound_cond
! 
!       SUBROUTINE dpowe_bound_cond(i,z,dg)
!             REAL*8, INTENT(IN) :: z(1)
!             INTEGER, INTENT(IN) :: i
!             REAL*8, INTENT(OUT) :: dg(1,1) ! Jacobian of g wrt "y"
!             dg = 1.D0 
!       END SUBROUTINE dpowe_bound_cond
! 
!       SUBROUTINE pow_i_rhs(x,z,f)
!             IMPLICIT NONE
!             REAL*8, INTENT(IN) :: z(1)
!             REAL(rprec), INTENT(IN) :: x
!             REAL*8, INTENT(OUT) :: f(1)
!             REAL*8 :: func_val
!             INTEGER :: ier
!             call EZspline_interp(Pow_i_spl,x,func_val,ier)
!             f(1) = func_val
!       END SUBROUTINE pow_i_rhs
! 
!       SUBROUTINE dpow_i_rhs(x,z,df)
!             IMPLICIT NONE
!             REAL*8, INTENT(IN) :: z(1)
!             REAL(rprec), INTENT(IN) :: x
!             REAL*8, INTENT(OUT) :: df(1,1) ! Jacobian of f wrt "y"
!             df(1,1) = 0.D0 ! We assume it does not depend on "y"
!       END SUBROUTINE dpow_i_rhs
! 
!       SUBROUTINE powi_bound_cond(i,z,g)
!             REAL*8, INTENT(IN) :: z(1)
!             INTEGER, INTENT(IN) :: i
!             REAL*8, INTENT(OUT) :: g
!             g = z(1) ! The boundary condition, basically what should be zero.
!       END SUBROUTINE powi_bound_cond
! 
!       SUBROUTINE dpowi_bound_cond(i,z,dg)
!             REAL*8, INTENT(IN) :: z(1)
!             INTEGER, INTENT(IN) :: i
!             REAL*8, INTENT(OUT) :: dg(1,1) ! Jacobian of g wrt "y"
!             dg = 1.D0 
!       END SUBROUTINE dpowi_bound_cond
! 
!       SUBROUTINE parte_rhs(x,z,f)
!             IMPLICIT NONE
!             REAL*8, INTENT(IN) :: z(1)
!             REAL(rprec), INTENT(IN) :: x
!             REAL*8, INTENT(OUT) :: f(1)
!             REAL*8 :: func_val
!             INTEGER :: ier
!             call EZspline_interp(Part_sple,x,func_val,ier)
!             f(1) = func_val
!       END SUBROUTINE parte_rhs
! 
!       SUBROUTINE dparte_rhs(x,z,df)
!             IMPLICIT NONE
!             REAL*8, INTENT(IN) :: z(1)
!             REAL(rprec), INTENT(IN) :: x
!             REAL*8, INTENT(OUT) :: df(1,1) ! Jacobian of f wrt "y"
!             df(1,1) = 0.D0 ! We assume it does not depend on "y"
!       END SUBROUTINE dparte_rhs
! 
!       SUBROUTINE parte_bound_cond(i,z,g) ! NEED TO THINK ABOUT THIS BOUNDARY CONDITION
!             REAL*8, INTENT(IN) :: z(1)
!             INTEGER, INTENT(IN) :: i
!             REAL*8, INTENT(OUT) :: g
!             g = z(1) ! The boundary condition, basically what should be zero.
!       END SUBROUTINE parte_bound_cond
! 
!       SUBROUTINE dparte_bound_cond(i,z,dg) ! NEED TO THINK ABOUT THIS BOUNDARY CONDITION
!             REAL*8, INTENT(IN) :: z(1)
!             INTEGER, INTENT(IN) :: i
!             REAL*8, INTENT(OUT) :: dg(1,1) ! Jacobian of g wrt "y"
!             dg = 1.D0 
!       END SUBROUTINE dparte_bound_cond
! 
!       SUBROUTINE parti_rhs(x,z,f)
!             IMPLICIT NONE
!             REAL*8, INTENT(IN) :: z(1)
!             REAL(rprec), INTENT(IN) :: x
!             REAL*8, INTENT(OUT) :: f(1)
!             REAL*8 :: func_val
!             INTEGER :: ier
!             call EZspline_interp(Part_spli,x,func_val,ier)
!             f(1) = func_val
!       END SUBROUTINE parti_rhs
! 
!       SUBROUTINE dparti_rhs(x,z,df)
!             IMPLICIT NONE
!             REAL*8, INTENT(IN) :: z(1)
!             REAL(rprec), INTENT(IN) :: x
!             REAL*8, INTENT(OUT) :: df(1,1) ! Jacobian of f wrt "y"
!             df(1,1) = 0.D0 ! We assume it does not depend on "y"
!       END SUBROUTINE dparti_rhs
! 
!       SUBROUTINE parti_bound_cond(i,z,g) ! NEED TO THINK ABOUT THIS BOUNDARY CONDITION
!             REAL*8, INTENT(IN) :: z(1)
!             INTEGER, INTENT(IN) :: i
!             REAL*8, INTENT(OUT) :: g
!             g = z(1) ! The boundary condition, basically what should be zero.
!       END SUBROUTINE parti_bound_cond
! 
!       SUBROUTINE dparti_bound_cond(i,z,dg) ! NEED TO THINK ABOUT THIS BOUNDARY CONDITION
!             REAL*8, INTENT(IN) :: z(1)
!             INTEGER, INTENT(IN) :: i
!             REAL*8, INTENT(OUT) :: dg(1,1) ! Jacobian of g wrt "y"
!             dg = 1.D0 
!       END SUBROUTINE dparti_bound_cond
      
      END MODULE stelltran_equilutils
