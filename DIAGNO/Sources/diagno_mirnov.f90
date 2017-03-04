!-----------------------------------------------------------------------
!     Module:        diagno_mirnov
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          04/16/2012
!     Description:   This subroutine calculates the response of the
!                    Mirnov arrays.  Here the array is assumed to be
!                    calculated via:
!                    signal[T]=\frac{\int\vec{b}\cdot d\vec{l}}{\int |d\vec{l}|}
!                    which is just the average magnetic field integrated
!                    along each array.
!-----------------------------------------------------------------------
      SUBROUTINE diagno_mirnov
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE diagno_runtime, pi2_diag => pi2
      USE biotsavart
      USE safe_open_mod
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, iunit, nseg, npts,i,j,k,ig
      REAL(rprec) :: nx, ny, nz, int_fac, xp1, yp1, zp1, xp2, yp2, zp2
      REAL(rprec) :: xvec(3),bvec(3)
      REAL(rprec), ALLOCATABLE :: flux(:)
      REAL(rprec), ALLOCATABLE :: xp(:,:), yp(:,:), rp(:,:), phip(:,:),&
                                  zp(:,:), dx(:,:), dy(:,:), dz(:,:), &
                                  xmid(:,:), ymid(:,:), zmid(:,:), &
                                  eff_area(:,:)
      
!-----------------------------------------------------------------------
!     External Functions
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      INTERFACE
         REAL FUNCTION db_bode(x0,y0,z0,x1,y1,z1,cg,ldb)
         USE stel_kinds, ONLY: rprec
         INTEGER, intent(in)     :: cg
         REAL(rprec), intent(in) :: x0,y0,z0,x1,y1,z1
         LOGICAL, intent(in), optional :: ldb
         END FUNCTION db_bode
         
         REAL FUNCTION db_midpoint(x0,y0,z0,x1,y1,z1,cg,ldb)
         USE stel_kinds, ONLY: rprec
         INTEGER, intent(in)     :: cg
         real(rprec), intent(in) :: x0,y0,z0,x1,y1,z1
         LOGICAL, intent(in), optional :: ldb
         END FUNCTION db_midpoint
      
         REAL FUNCTION db_simpson(x0,y0,z0,x1,y1,z1,cg,ldb)
         USE stel_kinds, ONLY: rprec
         INTEGER, intent(in)     :: cg
         REAL(rprec), intent(in) :: x0,y0,z0,x1,y1,z1
         LOGICAL, intent(in), optional :: ldb
         END FUNCTION db_simpson
      END INTERFACE
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      if(lverb) write(6,*)' --Calculating Mirnov Array Values'
      int_fac=1./int_step
      IF (lverb) THEN
         WRITE(6,'(5X,A,2(9X,A,9X))') 'Loop','nseg','flux'
      END IF
      
      ! Read in diagnostic file
      iunit = 36
      CALL safe_open(iunit,ier,TRIM(mirnov_file),'old','formatted')
      READ(iunit,*) nseg,npts
      ALLOCATE(xp(nseg,1:npts),yp(nseg,1:npts),zp(nseg,1:npts),rp(nseg,1:npts),phip(nseg,1:npts))
      ALLOCATE(dx(nseg,1:npts),dy(nseg,1:npts),dz(nseg,1:npts))
      ALLOCATE(xmid(nseg,0:npts),ymid(nseg,0:npts),zmid(nseg,0:npts))
      ALLOCATE(flux(nseg),eff_area(nseg,npts))
      DO i = 1, nseg
         DO j = 1, npts
            READ(iunit,*) rp(i,j), zp(i,j), phip(i,j), eff_area(i,j)
         END DO
      END DO
      CLOSE(iunit)
      phip = phip * onerad ! Input in theta
      xp = rp * cos(phip)
      yp = rp * sin(phip)
      
      ! Now create calculate helper arrays
      dx(:,2:npts-1) = xp(:,3:npts) - xp(:,1:npts-2)
      dy(:,2:npts-1) = yp(:,3:npts) - yp(:,1:npts-2)
      dz(:,2:npts-1) = zp(:,3:npts) - zp(:,1:npts-2)
      xmid(:,2:npts-1) = 0.5 * ( xp(:,3:npts) + xp(:,1:npts-2) )
      ymid(:,2:npts-1) = 0.5 * ( yp(:,3:npts) + yp(:,1:npts-2) )
      zmid(:,2:npts-1) = 0.5 * ( zp(:,3:npts) + zp(:,1:npts-2) )
      
      ! Calculate Field
      DO i = 1, nseg
         IF (lmut) THEN
            
         ELSE
            flux(i) = 0.0
            SELECT CASE (int_type)
            CASE('midpoint')
              DO j = 2, npts-1
                 DO k = 1, int_step
                    xp1 = xp(i,j) + (k-1)*int_fac*dx(i,j)
                    yp1 = yp(i,j) + (k-1)*int_fac*dy(i,j)
                    zp1 = zp(i,j) + (k-1)*int_fac*dz(i,j)
                    xp2 = xp(i,j) + k*int_fac*dx(i,j)
                    yp2 = yp(i,j) + k*int_fac*dy(i,j)
                    zp2 = zp(i,j) + k*int_fac*dz(i,j)
                    ier = 0
                    IF (.not.lvac) flux(i)=flux(i) + db_midpoint(xp1,yp1,zp1,xp2,yp2,zp2,ier,LDB = .TRUE.)
                    IF (lcoil) THEN
                       DO ig = 1, SIZE(coil_group)
                         IF (.not.luse_extcur(ig)) CYCLE
                         flux(i)=flux(i) + db_midpoint(xp1,yp1,zp1,xp2,yp2,zp2,ig,LDB = .TRUE.)
                       END DO
                    END IF
                 END DO
              END DO
            CASE('simpson')
              DO j = 2, npts-1
                 DO k = 1, int_step
                    xp1 = xp(i,j) + (k-1)*int_fac*dx(i,j)
                    yp1 = yp(i,j) + (k-1)*int_fac*dy(i,j)
                    zp1 = zp(i,j) + (k-1)*int_fac*dz(i,j)
                    xp2 = xp(i,j) + k*int_fac*dx(i,j)
                    yp2 = yp(i,j) + k*int_fac*dy(i,j)
                    zp2 = zp(i,j) + k*int_fac*dz(i,j)
                    ier = 0
                    IF (.not.lvac) flux(i)=flux(i) + db_simpson(xp1,yp1,zp1,xp2,yp2,zp2,ier,LDB = .TRUE.)
                    IF (lcoil) THEN
                       DO ig = 1, SIZE(coil_group)
                         IF (.not.luse_extcur(ig)) CYCLE
                         flux(i)=flux(i) + db_simpson(xp1,yp1,zp1,xp2,yp2,zp2,ig,LDB = .TRUE.)
                       END DO
                    END IF
                 END DO
              END DO
            CASE('bode')
              DO j = 2, npts-1
                 DO k = 1, int_step
                    xp1 = xp(i,j) + (k-1)*int_fac*dx(i,j)
                    yp1 = yp(i,j) + (k-1)*int_fac*dy(i,j)
                    zp1 = zp(i,j) + (k-1)*int_fac*dz(i,j)
                    xp2 = xp(i,j) + k*int_fac*dx(i,j)
                    yp2 = yp(i,j) + k*int_fac*dy(i,j)
                    zp2 = zp(i,j) + k*int_fac*dz(i,j)
                    ier = 0
                    IF (.not.lvac) flux(i)=flux(i) + db_bode(xp1,yp1,zp1,xp2,yp2,zp2,ier,LDB = .TRUE.)
                    IF (lcoil) THEN
                       DO ig = 1, SIZE(coil_group)
                         IF (.not.luse_extcur(ig)) CYCLE
                         flux(i)=flux(i) + db_bode(xp1,yp1,zp1,xp2,yp2,zp2,ig,LDB = .TRUE.)
                       END DO
                    END IF
                 END DO
              END DO
            CASE('adaptive')
              
            END SELECT
            IF (lverb) THEN
               WRITE(6,'(2X,i6,10X,I5,10X,E18.8)') i, npts, flux(i)
            END IF
         END IF
      END DO
      ier = 0
      
      ! Output diagnostic file
      iunit = 36
      CALL safe_open(iunit,ier,'diagno_mirnov.'//TRIM(id_string),'replace','formatted')
      WRITE(iunit,'(i6,1p,1ES22.12E3)')(i,flux(i),i=1,nseg)
      WRITE(iunit,'(A)')'   #   xp[m]      yp      zp      B.n [T]      flux [Wb]'
      CLOSE(iunit)
      
      ! Cleanup
      DEALLOCATE(xp,yp,zp,rp,phip)
      DEALLOCATE(dx,dy,dz,xmid,ymid,zmid)
      DEALLOCATE(flux,eff_area)
      
      
 
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE diagno_mirnov
