!-----------------------------------------------------------------------
!     Module:        fieldlines_calc_surface_fit
!     Authors:       S. Lazerson (samuel.lazeson@ipp.mpg.de)
!     Date:          11/06/2019
!     Description:   This subroutine 
!-----------------------------------------------------------------------
      SUBROUTINE fieldlines_calc_surface_fit(fit_line)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE, intrinsic :: iso_c_binding
      USE fieldlines_runtime
      USE fieldlines_lines
!-----------------------------------------------------------------------
!     Input Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(in) :: fit_line
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      LOGICAL, DIMENSION(:), ALLOCATABLE :: lmask
      INTEGER ::                  d1, d2, n, m, mfft, &
                                  nfft, mpol, ntor, mn
      INTEGER, DIMENSION(:), ALLOCATABLE :: idex
      REAL(rprec) :: x_old, y_old, norm
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: r, z, p, x, y, r0, z0, &
                                                th, dth, nth, tth, &
                                                xth, yth, dr, arr1, arr2, &
                                                raxis, zaxis
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: dx2d, dy2d, &
                                                  x_new, y_new, rmnc, zmns, &
                                                  r_temp,z_temp

      INTEGER :: i,j,n1,n2,nHull,nPoints,k
      REAL, DIMENSION(:,:), ALLOCATABLE :: Points, Hull
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: x2d, y2d

#if defined(FFTW3)
      ! FFTW3
      INCLUDE 'fftw3.f03'
      REAL(C_DOUBLE), POINTER, DIMENSION(:) :: func_1d
      COMPLEX(C_DOUBLE_COMPLEX), POINTER, DIMENSION(:) :: fmn_1d
      REAL(C_DOUBLE), POINTER, DIMENSION(:,:) :: func_2d
      COMPLEX(C_DOUBLE_COMPLEX), POINTER, DIMENSION(:,:) :: fmn_2d
      TYPE(C_PTR) :: plan
      TYPE(C_PTR) :: real_1d, complex_1d, real_2d, complex_2d
#endif

      ! Extract the surface and axis data
      i = nsteps
      ALLOCATE(r(i),z(i),p(i),r0(i),z0(i),x(i),y(i),th(i),dth(i))
      FORALL(j=0:i-1) r(j+1) = R_lines(fit_line,j)
      FORALL(j=0:i-1) z(j+1) = Z_lines(fit_line,j)
      FORALL(j=0:i-1) p(j+1) = PHI_lines(fit_line,j)
      FORALL(j=0:i-1) r0(j+1) = R_lines(1,j)
      FORALL(j=0:i-1) z0(j+1) = Z_lines(1,j)

      ! Sort into toroidal cutplanes
      n1 = nsteps/npoinc
      n2 = npoinc
      ALLOCATE(x2d(n1,n2),y2d(n1,n2))
      k=1
      DO i=1,n1
         DO j=1,n2
            x2d(i,j) = r(k)
            y2d(i,j) = z(k)
            k=k+1
         END DO
      END DO

      ! Now extract the convex hull for each
      nPoints = n1
      ALLOCATE(Points(2,0:n1-1),Hull(2,0:n1))
      DO j=1,n2
         Points(1,0:n1-1) = x2d(1:n1,j)
         Points(2,0:n1-1) = y2d(1:n1,j)
         CALL ConvHull(n1,Points,nHull,Hull)
         x2d(1:n1,j) = Hull(1,0:n1-1)
         y2d(1:n1,j) = Hull(2,0:n1-1)
      END DO
      DEALLOCATE(Points,Hull)

      DO j =1, n1
         WRITE(321,*) x2d(j,:)
      END DO

      DO j =1, n1
         WRITE(322,*) y2d(j,:)
      END DO

      DO j =1, n1
         WRITE(327,*) x2d(j,1),y2d(j,1)
      END DO

      ! Try a Manual DFT
      mfft = m/2+1
      nfft = 2*(DBLE(n)/2+1)
      ALLOCATE(r_temp(n,0:mfft-1),z_temp(n,0:mfft-1))
      ALLOCATE(rmnc(-n/2:n/2,0:mfft-1),zmns(-n/2:n/2,0:mfft-1))
      r_temp = 0; z_temp=0;
      DO i = 1, n1
         DO j=1, n2
            DO mn = 0, mfft-1
               r_temp(j,mn) = r_temp(j,mn) + x2d(i,j)*COS(mn*pi2*DBLE(i-1)/(n1))
               z_temp(j,mn) = z_temp(j,mn) + y2d(i,j)*SIN(mn*pi2*DBLE(i-1)/(n1))
            END DO
         END DO
      END DO
      DO mn = -n/2, n/2
         DO i=0, mfft-1
            DO j = 1, n
               rmnc(mn,i) = rmnc(mn,i) + r_temp(j,i)*COS(-mn*pi2*DBLE(j-1)/(n))
               zmns(mn,i) = zmns(mn,i) + r_temp(j,i)*SIN(-mn*pi2*DBLE(j-1)/(n))
            END DO
         END DO
      END DO

      DO j = 0, mfft-1
         DO i = -n/2, n/2
            IF ((j.eq.0) .and. (i<0)) CYCLE
            WRITE(555,'(2(2X,A,I4,A,I3,A,EN18.10))') 'RBC(',i,',',j,') = ', rmnc(i,j),'  ZBS(',i,',',j,') = ', zmns(i,j)
         END DO
      END DO


#if defined(FFTW3a)
      ! axis
      nfft = (n)/2+1
      PRINT *,'------------'
      PRINT *,'n = ',n
      PRINT *,'n/2 = ',n/2
      PRINT *,'nfft = ',nfft
      PRINT *,'phi = ',p(1:n)
      PRINT *,'------------'
      ALLOCATE(raxis(0:n/2),zaxis(0:n/2))
      real_1d = FFTW_ALLOC_REAL(INT(n,C_SIZE_T))
      complex_1d = FFTW_ALLOC_COMPLEX(INT(nfft,C_SIZE_T))
      CALL C_F_POINTER(real_1d, func_1d, [n])
      CALL C_F_POINTER(complex_1d, fmn_1d, [nfft])
      norm = DBLE(n)/2
      !R
      func_1d = r0(1:n)
      plan = FFTW_PLAN_DFT_R2C_1D(n, func_1d, fmn_1d, FFTW_ESTIMATE)
      CALL FFTW_EXECUTE_DFT_R2C(plan, func_1d, fmn_1d)
      CALL FFTW_DESTROY_PLAN(plan)
      WRITE(500,*) fmn_1d
      DO j = 1, n/2
         ntor = j-1
         raxis(ntor) = REAL(fmn_1d(j))/norm
      END DO
      raxis(0) = raxis(0)/2
      !Z
      func_1d = z0(1:n)
      plan = FFTW_PLAN_DFT_R2C_1D(n, func_1d, fmn_1d, FFTW_ESTIMATE)
      CALL FFTW_EXECUTE_DFT_R2C(plan, func_1d, fmn_1d)
      CALL FFTW_DESTROY_PLAN(plan)
      DO j = 1, n/2
         ntor = j-1
         zaxis(ntor) = -AIMAG(fmn_1d(j))/norm
      END DO
      DO j = 0, nfft-1
         WRITE(555,'(2(2X,A,I3,A,ES18.10))') 'RAXIS(',j,') = ', raxis(j),'  ZAXIS(',j,') = ', zaxis(j)
      END DO
      DEALLOCATE(raxis, zaxis)



      ! Surface
      mfft = m/2+1
      nfft = 2*(DBLE(n)/2+1)
      ALLOCATE(rmnc(-n/2:n/2,0:mfft-1),zmns(-n/2:n/2,0:mfft-1))
      real_2d = FFTW_ALLOC_REAL(INT(m*n,C_SIZE_T))
      complex_2d = FFTW_ALLOC_COMPLEX(INT(mfft*nfft,C_SIZE_T))
      CALL C_F_POINTER(real_2d, func_2d, [m,n])
      CALL C_F_POINTER(complex_2d, fmn_2d, [mfft,nfft])
      norm = DBLE(m*n)/2

      ! Do x
      func_2d = x_new
      DO i = 1, m
         WRITE(601,*) func_2d(i,1:n); CALL FLUSH(601)
      END DO
      plan = FFTW_PLAN_DFT_R2C_2D(n, m, func_2d, fmn_2d, FFTW_ESTIMATE)
      CALL FFTW_EXECUTE_DFT_R2C(plan, func_2d, fmn_2d)
      CALL FFTW_DESTROY_PLAN(plan)
      ! Transform to VMEC format
      DO i = 1, mfft
         DO j = 1 , n/2+1
            ntor = j-1
            mpol = i-1
            rmnc(ntor,mpol) = REAL(fmn_2d(i,j))/norm
            IF (i==1) rmnc(ntor,mpol) = rmnc(ntor,mpol) / 2
            IF (i==mfft) rmnc(ntor,mpol) = rmnc(ntor,mpol) / 2
         END DO
         DO j = n/2+2, nfft-2
            mpol = i-1
            ntor = -n/2+j-(n/2+2)+1
            rmnc(ntor,mpol) = REAL(fmn_2d(i,j))/norm
            IF (i==1) rmnc(ntor,mpol) = rmnc(ntor,mpol) / 2
            IF (i==mfft) rmnc(ntor,mpol) = rmnc(ntor,mpol) / 2
         END DO
      END DO
      IF(MOD(m,2) == 1) rmnc(:,mfft-1) = rmnc(:,mfft-1)*2

      ! Do y
      func_2d = y_new
      DO i = 1, m
         WRITE(701,*) func_2d(i,1:n); CALL FLUSH(601)
      END DO
      plan = FFTW_PLAN_DFT_R2C_2D(n, m, func_2d, fmn_2d, FFTW_ESTIMATE)
      CALL FFTW_EXECUTE_DFT_R2C(plan, func_2d, fmn_2d)
      CALL FFTW_DESTROY_PLAN(plan)
      ! Transform to VMEC format
      DO i = 1, mfft
         DO j = 1 , n/2+1
            ntor = j-1
            mpol = i-1
            zmns(ntor,mpol) =-AIMAG(fmn_2d(i,j))/norm
            IF (i==1) zmns(ntor,mpol) = zmns(ntor,mpol) / 2
            IF (i==mfft) zmns(ntor,mpol) = zmns(ntor,mpol) / 2
         END DO
         DO j = n/2+2, nfft-2
            mpol = i-1
            ntor = -n/2+j-(n/2+2)+1
            zmns(ntor,mpol) =-AIMAG(fmn_2d(i,j))/norm
            IF (i==1) zmns(ntor,mpol) = zmns(ntor,mpol) / 2
            IF (i==mfft) zmns(ntor,mpol) = zmns(ntor,mpol) / 2
         END DO
      END DO
      IF(MOD(m,2) == 1) zmns(:,mfft-1) = zmns(:,mfft-1)*2

      DO j = 0, 24
         DO i = -n/2, n/2
            IF ((j.eq.0) .and. (i<0)) CYCLE
            WRITE(555,'(2(2X,A,I4,A,I3,A,EN18.10))') 'RBC(',-i,',',j,') = ', rmnc(i,j),'  ZBS(',-i,',',j,') = ', zmns(i,j)
         END DO
      END DO
      CALL FFTW_FREE(real_2d)
      CALL FFTW_FREE(complex_2d)
#endif






      ! DEALLOCATION
      DEALLOCATE(rmnc,zmns)
      DEALLOCATE(x_new,y_new)
      DEALLOCATE(nth,arr2)
      DEALLOCATE(dx2d,dy2d)
      DEALLOCATE(x2d,y2d)
      DEALLOCATE(tth,xth,yth,dr,arr1)
      DEALLOCATE(idex)
      DEALLOCATE(lmask)
      DEALLOCATE(r,z,p,r0,z0,x,y,th,dth)      
      
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE fieldlines_calc_surface_fit



SUBROUTINE ConvHull(nPoints,Points,nHull,Hull)
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE 
!------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: nPoints
REAL,INTENT(IN)     :: Points(2,0:nPoints-1)
!------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT) :: nHull
! NOTE: allocate Hull always one point greater than Points, because we save the first value twice
REAL,INTENT(OUT)    :: Hull(2,0:nPoints)
!------------------------------------------------
! LOCAL VARIABLES
REAL                :: Lower(2,0:nPoints-1)
REAL                :: Upper(2,0:nPoints-1)
INTEGER             :: i,iLower,iUpper
!================================================
IF(nPoints.LE.1)THEN
  Hull  = Points
  nHull = nPoints
ELSE
  iLower = 0
  Lower  = -HUGE(1.)
  DO i=0,nPoints-1
    DO WHILE(iLower.GE.2.AND.Cross(Lower(:,iLower-2),Lower(:,iLower-1),Points(:,i)).LE.0.)
      Lower(:,iLower) = -HUGE(1.)
      iLower          = iLower - 1
    END DO
    Lower(:,iLower) = Points(:,i)
    iLower = iLower + 1
  END DO

  iUpper = 0
  Upper  = HUGE(1.)
  DO i=nPoints-1,0,-1
    DO WHILE(iUpper.GE.2.AND.Cross(Upper(:,iUpper-2),Upper(:,iUpper-1),Points(:,i)).LE.0.)
      Upper(:,iUpper) = HUGE(1.)
      iUpper          = iUpper - 1
    END DO
    Upper(:,iUpper) = Points(:,i)
    iUpper = iUpper + 1
  END DO

  iLower = iLower-1
  iUpper = iUpper-1
  nHull  = iLower+iUpper+1
  
  ! NOTE: Initialize Hull with zeros
  Hull   = 0.

  ! NOTE: save values in Hull
  Hull(:,0     :iLower       -1) = Lower(:,0:iLower-1)
  Hull(:,iLower:iLower+iUpper-1) = Upper(:,0:iUpper-1)

  ! NOTE: save first value twice
  Hull(:,       iLower+iUpper  ) = Hull (:,0         )
END IF


CONTAINS
FUNCTION Cross(v1,v2,v3)
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: v1(2)    !< input vector 1
REAL,INTENT(IN) :: v2(2)    !< input vector 2
REAL,INTENT(IN) :: v3(2)    !< input vector 3
!-----------------------------------------------
! OUTPUT VARIABLES
REAL            :: Cross    !< cross product
!-----------------------------------------------
! LOCAL VARIABLES
!===============================================
Cross=(v2(1)-v1(1))*(v3(2)-v1(2))-(v2(2)-v1(2))*(v3(1)-v1(1))
END FUNCTION Cross

END SUBROUTINE
