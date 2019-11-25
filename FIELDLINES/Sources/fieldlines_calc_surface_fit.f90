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
      INTEGER ::                  i, j, n1, n2, d1, d2, n, m, mfft, &
                                  nfft, mpol, ntor
      INTEGER, DIMENSION(:), ALLOCATABLE :: idex
      REAL(rprec) :: x_old, y_old, norm
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: r, z, p, x, y, r0, z0, &
                                                th, dth, nth, tth, &
                                                xth, yth, dr, arr1, arr2
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: x2d, y2d, dx2d, dy2d, &
                                                  x_new, y_new, rmnc, zmns

#if defined(FFTW3)
      ! FFTW3
      INCLUDE 'fftw3.f03'
      REAL(C_DOUBLE), POINTER, DIMENSION(:,:) :: func_2d
      COMPLEX(C_DOUBLE_COMPLEX), POINTER, DIMENSION(:,:) :: fmn_2d
      TYPE(C_PTR) :: plan
      TYPE(C_PTR) :: real_2d, complex_2d
#endif

      ! Extract the surface and axis data
      i = nsteps
      ALLOCATE(r(i),z(i),p(i),r0(i),z0(i),x(i),y(i),th(i),dth(i))
      FORALL(j=0:i-1) r(j+1) = R_lines(fit_line,j)
      FORALL(j=0:i-1) z(j+1) = Z_lines(fit_line,j)
      FORALL(j=0:i-1) p(j+1) = PHI_lines(fit_line,j)
      FORALL(j=0:i-1) r0(j+1) = R_lines(1,j)
      FORALL(j=0:i-1) z0(j+1) = Z_lines(1,j)

!      DO j =1, nsteps
!         WRITE(325,*) r(j),z(j)
!      END DO

      ! Calc helpers
      x = r-r0
      y = z-z0
      th = atan2(y,x)
      dth(1:i-1)=th(2:i)-th(1:i-1)
      th(1) = 0
      DO i = 2, nsteps
         th(i) = th(i-1)+dth(i-1)
      END DO

!      DO j =1, nsteps
!         WRITE(326,*) x(j),y(j)
!      END DO

      ! Sort the data in 2D grid
      n1 = nsteps/npoinc
      n2 = npoinc
      ALLOCATE(lmask(n1))
      ALLOCATE(idex(n1))
      ALLOCATE(tth(n1),xth(n1),yth(n1),dr(n1),arr1(n1))
      ALLOCATE(x2d(n1,n2),y2d(n1,n2))
      DO i = 1, n2
         d1 = 1
         DO j = i, nsteps-1, npoinc
            xth(d1) = x(j)
            yth(d1) = y(j)
            d1 = d1 + 1
            IF (d1 > n1) EXIT
         END DO
         tth=atan2(yth,xth)
         WHERE(tth<0) tth = tth + pi2
         CALL SORT(n1,tth,idex)
         DO j = 1, n1
            arr1(j) = xth(idex(j))
         END DO
         xth = arr1
         DO j = 1, n1
            arr1(j) = yth(idex(j))
         END DO
         yth = arr1
         x2d(1,i) = xth(1)
         y2d(1,i) = yth(1)
         x_old = xth(1)
         y_old = yth(1)
         lmask = .TRUE.
         lmask(1) = .FALSE.
         DO j = 2, n1
            dr = (xth-x_old)**2+(yth-y_old)**2
            d1 = MINLOC(dr,1,MASK=lmask)
            x2d(j,i) = xth(d1)
            y2d(j,i) = yth(d1)
            x_old = xth(d1)
            y_old = yth(d1)
            lmask(d1) = .false.
         END DO
      END DO

      DO j =1, n1
         WRITE(321,*) x2d(j,:)
      END DO

      DO j =1, n1
         WRITE(322,*) y2d(j,:)
      END DO

      DO j =1, n1
         WRITE(327,*) x2d(j,1),y2d(j,1)
      END DO

      ! Now we calculate a poloidal like angle
      ALLOCATE(dx2d(n1,n2),dy2d(n1,n2))
      dx2d=0; dy2d=0;
      DO i = 1, n1-1
         dx2d(i,:)=x2d(i+1,:)-x2d(i,:)
         dy2d(i,:)=y2d(i+1,:)-y2d(i,:)
      END DO
      dx2d = SQRT(dx2d*dx2d+dy2d*dy2d)
      dy2d = 0
      DO i = 2, n1
         dy2d(i,:) = dy2d(i-1,:)+dx2d(i-1,:)
      END DO
      DO i = 1, n2
         dy2d(:,i) = dy2d(:,i)/MAXVAL(dy2d(:,i),1)
      END DO

      DO j =1, n1
         WRITE(328,*) dx2d(j,1),dy2d(j,1)
      END DO

      ! Now calculated new 2D grid
      m = 256
      n = npoinc
      ALLOCATE(nth(m),arr2(m))
      ALLOCATE(x_new(m,n),y_new(m,n))
      FORALL(i=1:m) nth(i) = DBLE(i-1)/DBLE(256)
      DO i = 1, n-1
         xth = dy2d(:,i)
         yth = x2d(:,i)
         CALL spline_it(n1,xth,yth,m,nth,arr2,0)
         x_new(:,i) = arr2
         yth = y2d(:,i)
         CALL spline_it(n1,xth,yth,m,nth,arr2,0)
         y_new(:,i) = arr2
      END DO

      ! Make sure we start at max x
      DO i = 1, n-1
         arr2 = x_new(:,i)
         j = MAXLOC(arr2,1)
         IF (j >1) THEN
            d1 = m-j+1
            d2 = d1 + 1
            arr2(1:d1)   = x_new(j:m,i)
            arr2(d2:m) = x_new(1:j-1,i)
            x_new(:,i)   = arr2
            arr2(1:d1)   = y_new(j:m,i)
            arr2(d2:m) = y_new(1:j-1,i)
            y_new(:,i)    = arr2
         END IF
      END DO
      x_new(:,n2) = x_new(:,1)
      y_new(:,n2) = y_new(:,1)

      DO j =1, m
         WRITE(329,*) x_new(j,1),y_new(j,1)
      END DO

      DO j =1, m
         WRITE(401,*) x_new(j,:)
      END DO

      DO j =1, m
         WRITE(402,*) y_new(j,:)
      END DO

#if defined(FFTW3)
      mfft = m/2+1
      nfft = 2*(DBLE(n)/2+1)
      ALLOCATE(rmnc(-n/2:n/2,0:mfft-1),zmns(-n/2:n/2,0:mfft-1))
      real_2d = FFTW_ALLOC_REAL(INT(m*n,C_SIZE_T))
      complex_2d = FFTW_ALLOC_COMPLEX(INT(mfft*nfft,C_SIZE_T))
      CALL C_F_POINTER(real_2d, func_2d, [m,n])
      CALL C_F_POINTER(complex_2d, fmn_2d, [mfft,nfft])
      norm = m*n/2

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

      DO i = -n/2, n/2
         DO j = 0, 24
            WRITE(555,*) 'RMNC(',i,',',j,') = ', rmnc(i,j),'  ZMNS(',i,',',j,') = ', zmns(i,j)
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
