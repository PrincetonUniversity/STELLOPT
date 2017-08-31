!-----------------------------------------------------------------------
!     Subroutine:    stellopt_spline_to_coil
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/14/2016
!     Description:   This subroutine generates a coil data structure from
!                    the stellopt splines and writes it to a file.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_spline_to_coil(lscreen)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE vmec_input, ONLY: nfp, lasym, extcur
      USE stellopt_vars
      USE equil_utils, nfp_stel => nfp, lasym_stel=>lasym, extcur_stel=>extcur
      USE safe_open_mod
      IMPLICIT NONE
!-----------------------------------------------------------------------
!     Subroutine Parameters
!----------------------------------------------------------------------
      LOGICAL, INTENT(inout)        :: lscreen
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      LOGICAL ::  lmodular
      INTEGER ::  i, j, k, ier, n_coilgroups, iunit
      REAL(rprec) :: ph, rc_max, rc_min,zc_max,zc_min
      REAL(rprec), ALLOCATABLE :: x_coil(:),y_coil(:),z_coil(:),c_coil(:), &
                                  p_coil(:), r_coil(:)
      CHARACTER(LEN=8) :: coil_name

      INTEGER, PARAMETER :: nseg = 128
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (lscreen) WRITE(6,*) '---------------------------  WRITING COIL  ------------------------'
      ! Allocate coil helpers
      ALLOCATE(x_coil(nseg),y_coil(nseg),z_coil(nseg),c_coil(nseg))
      IF (.not. lasym) ALLOCATE(r_coil(nseg),p_coil(nseg))
      rc_max = 0; rc_min = 1000; zc_min = 100; zc_max = -100

      ! Count coil groups
      n_coilgroups = 0
      DO i = 1, nigroup
         IF (ANY(coil_splinesx(i,:) >-1)) n_coilgroups=n_coilgroups+1
      END DO

      ! Open files and print header
      CALL safe_open(iunit,ier,TRIM('coils.'//TRIM(proc_string)),'unknown','formatted')
      WRITE(iunit,'(A,I3)') 'periods',nfp
      WRITE(iunit,'(A)') 'begin filament'
      WRITE(iunit,'(A)') 'mirror NIL'
      CALL FLUSH(iunit)

      ! Print out some info
      IF (lscreen) THEN
          WRITE(6,'(2X,A,A)')  'COIL NAME: ',TRIM('coils.'//TRIM(proc_string))
          IF (lwindsurf) WRITE(6,'(2X,A,A)')  '      CWS: ',TRIM(windsurfname)
          WRITE(6,'(2X,A,I2)') '      NFP: ',nfp
          WRITE(6,'(2X,A,I2)') '   GROUPS: ',n_coilgroups
          CALL FLUSH(6)
      END IF

      ! Write Coils
      DO i = 1, n_coilgroups
         ! Setup the Coil name
         WRITE(coil_name,'(A,I3.3)') 'COIL_',i

         ! Reconstruct coil geometry from spline data
         CALL spline_to_coil(i, nseg, x_coil, y_coil, z_coil, lmodular)
         c_coil(1:nseg) = extcur(i)
         IF (extcur(i) == 0) c_coil = 1

         ! Write Coil
         IF (lmodular) THEN ! Coil closes in a field period
            ! Screen output
            IF (lscreen) WRITE(6,'(5X,A,2X,I3,2X,A)') coil_name,nseg,'MODULAR'

            ! Setup the coil
            x_coil(nseg) = x_coil(1); y_coil(nseg)=y_coil(1); z_coil(nseg)=z_coil(1)
            r_coil = sqrt(x_coil*x_coil + y_coil*y_coil)
            p_coil = atan2(y_coil,x_coil)

            ! Update Helpers
            rc_max = MAX(rc_max,MAXVAL(r_coil))
            rc_min = MIN(rc_min,MINVAL(r_coil))
            zc_max = MAX(zc_max,MAXVAL(z_coil))
            zc_min = MIN(zc_min,MINVAL(z_coil))

            ! Write the coil each field period
            DO k = 1, nfp
               ph = pi2*(k-1)/nfp
               DO j = 1, nseg-1
                  WRITE(iunit,'(4ES20.12)') r_coil(j)*cos(p_coil(j)+ph),&
                                            r_coil(j)*sin(p_coil(j)+ph),&
                                            z_coil(j),c_coil(j)
               END DO
               WRITE(iunit,'(4ES20.12,I5,1X,A)') r_coil(1)*cos(p_coil(1)+ph),&
                                                 r_coil(1)*sin(p_coil(1)+ph),&
                                                 z_coil(1),0.,i,coil_name
            END DO
            CALL FLUSH(6)

            ! Handle stellarator symmetry
            IF (.not. lasym) THEN
               ! Screen output
               ! Mirror Coil
               p_coil = pi2/nfp - p_coil
               z_coil = -z_coil
               DO k = 1, nfp
                  ph = pi2*(k-1)/nfp
                  DO j = nseg,2,-1
                     WRITE(iunit,'(4ES20.12)') r_coil(j)*cos(p_coil(j)+ph),&
                                               r_coil(j)*sin(p_coil(j)+ph),&
                                               z_coil(j),c_coil(j)
                  END DO
                  WRITE(iunit,'(4ES20.12,I5,1X,A)') r_coil(nseg)*cos(p_coil(nseg)+ph),&
                                                    r_coil(nseg)*sin(p_coil(nseg)+ph),&
                                                    z_coil(nseg),0.,i,coil_name
               END DO
            END IF
         ELSE ! Coil is toroidal segment
            IF (lscreen) WRITE(6,'(5X,A,2X,I3,2X,A)') coil_name,nseg,'REGULAR'
            r_coil = sqrt(x_coil*x_coil + y_coil*y_coil)
            p_coil = atan2(y_coil,x_coil)

            ! Update Helpers
            rc_max = MAX(rc_max,MAXVAL(r_coil))
            rc_min = MIN(rc_min,MINVAL(r_coil))
            zc_max = MAX(zc_max,MAXVAL(z_coil))
            zc_min = MIN(zc_min,MINVAL(z_coil))

            ! Write the coil
            DO k = 1, nfp
               ph = pi2*(k-1)/nfp
               DO j = 1, nseg-1
                  WRITE(iunit,'(4ES20.12)') r_coil(j)*cos(p_coil(j)+ph),&
                                            r_coil(j)*sin(p_coil(j)+ph),&
                                            z_coil(j),c_coil(j)
               END DO
            END DO
            WRITE(iunit,'(4ES20.12,I5,1X,A)') r_coil(1)*cos(p_coil(1)),&
                                              r_coil(1)*sin(p_coil(1)),&
                                              z_coil(1),0.,i,coil_name
         END IF
      END DO
      WRITE(iunit,'(A)') 'end'

      ! Close file
      CLOSE(iunit)

      ! General Info
      IF (lscreen) THEN
          WRITE(6,'(A,F8.5,A,F8.5,A)') '   R   = [',rc_min,',',rc_max,']'
          WRITE(6,'(A,F8.5,A,F8.5,A)') '   Z   = [',zc_min,',',zc_max,']'
          CALL FLUSH(6)
      END IF

      ! Clean up
      DEALLOCATE(x_coil,y_coil,z_coil,c_coil)
      IF (ALLOCATED(r_coil)) DEALLOCATE(r_coil)
      IF (ALLOCATED(p_coil)) DEALLOCATE(p_coil)

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_spline_to_coil

!----------------------------------------------------------------------
!     Subroutine:    spline_to_coil
!     Author:        J. Breslau (jbreslau@pppl.gov)
!     Date:          8/22/2017
!     Description:   This subroutine performs spline reconstruction of
!                    a single coil into a list of (x,y,z) vertices.
!-----------------------------------------------------------------------
      SUBROUTINE spline_to_coil(icoil, nseg, xarr, yarr, zarr, lmod)
        USE stellopt_vars
        USE windingsurface
        USE vmec_input,  ONLY : nfp
        USE equil_utils, ONLY : bcs1
        USE EZspline_obj
        USE EZspline
        IMPLICIT NONE
        INTRINSIC ABS, MODULO, EPSILON

        ! Arguments
        INTEGER, INTENT(IN)                       :: icoil, nseg
        REAL(rprec), DIMENSION(nseg), INTENT(OUT) :: xarr, yarr, zarr
        LOGICAL, INTENT(OUT)                      :: lmod  ! Does the coil close on itself within a field period?
                                                           ! Note: this is true for saddle coils as well.
        ! Local variables
        TYPE(EZspline1_r8) :: XC_spl, YC_spl, ZC_spl
        REAL(rprec)        :: s_val, u, v, xlast, ylast
        INTEGER            :: nknots, ier, j

        lmod = .FALSE.

        ! Handle spline boundary conditions
        nknots = COUNT(coil_splinesx(icoil,:) >= 0)
        IF (coil_type(icoil).eq.'M') THEN !Enforce periodicity for modular coils
           lmod = .TRUE.
           IF (lwindsurf) THEN !u0=0,uf=1,vf=v0
              coil_splinefx(icoil,     1) = 0.0d0
              coil_splinefx(icoil,nknots) = 1.0d0
              coil_splinefy(icoil,nknots) = coil_splinefy(icoil,1)
           ELSE  !z0=z(t=0), zf=z0, xf=x0,yf=y0
              coil_splinefx(icoil,nknots) = coil_splinefx(icoil,1)
              coil_splinefy(icoil,nknots) = coil_splinefy(icoil,1)
              coil_splinefz(icoil,nknots) = coil_splinefz(icoil,1)
           END IF
        END IF

        ! Set up splines
        CALL EZspline_init(XC_spl, nknots, bcs1, ier)
        XC_spl%x1 = coil_splinesx(icoil,1:nknots)
        CALL EZspline_setup(XC_spl, coil_splinefx(icoil,1:nknots), ier)
        xlast = coil_splinefx(icoil,nknots)

        CALL EZspline_init(YC_spl, nknots, bcs1, ier)
        YC_spl%x1 = coil_splinesy(icoil,1:nknots)
        CALL EZspline_setup(YC_spl, coil_splinefy(icoil,1:nknots), ier)
        ylast = coil_splinefy(icoil,nknots)

        IF (lwindsurf) THEN !Interpret coil_splinefx as u, coil_splinefy as v
           windsurf%nfp = nfp

           ! Evaluate the coil
           DO j = 1, nseg
              s_val = REAL(j-1)/REAL(nseg-1)
              CALL EZspline_interp(XC_spl, s_val, u, ier)
              CALL EZspline_interp(YC_spl, s_val, v, ier)
              CALL stellopt_uv_to_xyz(u, v, xarr(j), yarr(j), zarr(j))
           END DO

           IF ((ABS(MODULO(coil_splinefx(icoil,1),1.0d0) - MODULO(xlast,1.0D0)).le.2.0*EPSILON(u)) &
                .and. (coil_splinefy(icoil,1) == ylast) ) lmod = .TRUE.
        ELSE                !Interpret coil_splinefx as x, coil_splinefy as y, coil_splinefz as z.
           CALL EZspline_init(ZC_spl, nknots, bcs1, ier)
           ZC_spl%x1 = coil_splinesz(icoil,1:nknots)
           CALL EZspline_setup(ZC_spl, coil_splinefz(icoil,1:nknots), ier)

           ! Evaluate the coil
           DO j = 1, nseg
              s_val = REAL(j-1)/REAL(nseg-1)
              CALL EZspline_interp(XC_spl, s_val, xarr(j), ier)
              CALL EZspline_interp(YC_spl, s_val, yarr(j), ier)
              CALL EZspline_interp(ZC_spl, s_val, zarr(j), ier)
           END DO

           CALL EZspline_free(ZC_spl, ier)

           IF ((coil_splinefx(icoil,1) == xlast) &
                .and. (coil_splinefy(icoil,1) == ylast) &
                .and. (coil_splinefz(icoil,1) == coil_splinefz(icoil,nknots))) &
                lmod = .TRUE.
        END IF

        CALL EZspline_free(XC_spl, ier);  CALL EZspline_free(YC_spl, ier)
      END SUBROUTINE spline_to_coil

!-----------------------------------------------------------------------
!     Subroutine:    get_coil_length
!     Author:        J. Breslau (jbreslau@pppl.gov)
!     Date:          8/22/2017
!     Description:   Computes the length of a coil.
!     Input:         Coil index icoil
!     Output:        Length of coil
!-----------------------------------------------------------------------
      SUBROUTINE get_coil_length(icoil, length)
        USE stel_kinds, ONLY : rprec
        IMPLICIT NONE
        INTRINSIC SQRT

        ! Arguments
        INTEGER, INTENT(IN)      :: icoil
        REAL(rprec), INTENT(OUT) :: length

        ! Constants
        INTEGER, PARAMETER :: nseg_len = 360

        ! Local variables
        REAL(rprec), DIMENSION(nseg_len) :: xc, yc, zc
        INTEGER                          :: iseg
        LOGICAL                          :: ldum

        ! Get x,y,z coords along coil
        CALL spline_to_coil(icoil, nseg_len, xc, yc, zc, ldum)

        ! Compute length
        length = 0.0
        DO iseg=1,nseg_len-1  ! Assumes xc(nseg) == xc(1), etc.
           length = length + SQRT((xc(iseg+1) - xc(iseg))**2 + &
                                  (yc(iseg+1) - yc(iseg))**2 + &
                                  (zc(iseg+1) - zc(iseg))**2)
        END DO !iseg
      END SUBROUTINE get_coil_length

!-------------------------------------------------------------------------------
!     Subroutine:    get_coil_sep
!     Author:        J. Breslau (jbreslau@pppl.gov)
!     Date:          8/25/2017
!     Description:   Estimates the distance of closest approach of two coils.
!
!     Inputs:       Coil coordinate arrays and sizes
!     Output:       Shortest distance between any two coordinate pairs.
!-------------------------------------------------------------------------------
      SUBROUTINE get_coil_sep(x1, y1, z1, n1, x2, y2, z2, n2, minsep)
        USE stel_kinds, ONLY : rprec        
        IMPLICIT NONE
        INTRINSIC SQRT

        INTEGER, INTENT(IN)                    :: n1, n2
        REAL(rprec), DIMENSION(n1), INTENT(IN) :: x1, y1, z1
        REAL(rprec), DIMENSION(n2), INTENT(IN) :: x2, y2, z2
        REAL(rprec), INTENT(OUT)               :: minsep

        REAL(rprec) :: dsq
        INTEGER     :: j1, j2

        minsep = 1.0D+60

        DO j1=1,n1-1     ! Assumes 1st and last pts are identical.
           DO j2=1,n2-1
              dsq = (x1(j1) - x2(j2))**2 + (y1(j1) - y2(j2))**2 + (z1(j1) - z2(j2))**2
              IF (dsq < minsep) minsep = dsq
           END DO !j2
        END DO !j1

        minsep = SQRT(minsep)
      END SUBROUTINE get_coil_sep

!-------------------------------------------------------------------------------
!     Subroutine:    get_coil_maxcurv
!     Author:        J. Breslau (jbreslau@pppl.gov)
!     Date:          8/24/2017
!     Description:   Computes the curvature at several points along a selected
!                    coil and returns the maximum value and its location.
!     Input:         Coil index icoil
!     Outputs:       Max. curvature, parametric position of max. curvature.
!-------------------------------------------------------------------------------
      SUBROUTINE get_coil_maxcurv(icoil, maxcurv, s_max)
        USE stellopt_vars
        USE windingsurface
        USE stellopt_targets, ONLY : npts_curv
        USE vmec_input,  ONLY : nfp
        USE EZspline_obj
        USE EZspline
        IMPLICIT NONE
        INTRINSIC COUNT, SQRT

        ! Arguments
        INTEGER, INTENT(IN)      :: icoil
        REAL(rprec), INTENT(OUT) :: maxcurv, s_max

        ! Constants
        INTEGER, DIMENSION(2), PARAMETER :: ezpbc = (/-1,-1/) ! Periodic boundary condition in EZspline

        ! Local variables
        TYPE(EZspline1_r8) :: XC_spl, YC_spl, ZC_spl
        REAL(rprec)        :: xprime(3,5)
        REAL(rprec)        :: s_val, u, duds, d2uds2, v, dvds, d2vds2
        REAL(rprec)        :: dxds, dyds, dzds, d2xds2, d2yds2, d2zds2
        REAL(rprec)        :: cx, cy, cz, norm, kappa
        INTEGER            :: ier, ipt, nknots

        ! Handle spline boundary conditions
        nknots = COUNT(coil_splinesx(icoil,:) >= 0)
        IF (coil_type(icoil).eq.'M') THEN !Enforce periodicity for modular coils
           IF (lwindsurf) THEN !u0=0,uf=1,vf=v0
              coil_splinefx(icoil,     1) = 0.0d0
              coil_splinefx(icoil,nknots) = 1.0d0
              coil_splinefy(icoil,nknots) = coil_splinefy(icoil,1)
           ELSE  !z0=z(t=0), zf=z0, xf=x0,yf=y0
              coil_splinefx(icoil,nknots) = coil_splinefx(icoil,1)
              coil_splinefy(icoil,nknots) = coil_splinefy(icoil,1)
              coil_splinefz(icoil,nknots) = coil_splinefz(icoil,1)
           END IF
        END IF

        ! Set up splines
        CALL EZspline_init(XC_spl, nknots, ezpbc, ier)
        XC_spl%x1 = coil_splinesx(icoil,1:nknots)
        CALL EZspline_setup(XC_spl, coil_splinefx(icoil,1:nknots), ier)

        CALL EZspline_init(YC_spl, nknots, ezpbc, ier)
        YC_spl%x1 = coil_splinesy(icoil,1:nknots)
        CALL EZspline_setup(YC_spl, coil_splinefy(icoil,1:nknots), ier)

        IF (lwindsurf) THEN
           windsurf%nfp = nfp
        ELSE
           CALL EZspline_init(ZC_spl, nknots, ezpbc, ier)
           ZC_spl%x1 = coil_splinesz(icoil,1:nknots)
           CALL EZspline_setup(ZC_spl, coil_splinefz(icoil,1:nknots), ier)
        END IF

        ! Find maximum curvature along coil
        maxcurv = -1.0;  s_max = -1.0
        DO ipt=1,npts_curv
           s_val = REAL(ipt-1)/REAL(npts_curv)

           ! We need 1st and 2nd derivatives of (x,y,z) w/r/t s to compute curvature:
           IF (lwindsurf) THEN ! Interpret splinefx as u, splinefy as v
              ! First get 1st and 2nd derivs of u,v w/r/t s...
              CALL EZspline_interp(XC_spl, s_val, u, ier)
              CALL EZspline_derivative(XC_spl, 1, s_val, duds, ier)
              CALL EZspline_derivative(XC_spl, 2, s_val, d2uds2, ier)

              CALL EZspline_interp(YC_spl, s_val, v, ier)
              CALL EZspline_derivative(YC_spl, 1, s_val, dvds, ier)
              CALL EZspline_derivative(YC_spl, 2, s_val, d2vds2, ier)

              ! Next get 1st & 2nd derivs of (x,y,z) w/r/t (u,v)...
              CALL stellopt_uv_to_xyz_prime(u, v, xprime)

              ! Finally, apply the chain rule...
              dxds = xprime(1,1)*duds + xprime(1,2)*dvds
              dyds = xprime(2,1)*duds + xprime(2,2)*dvds
              dzds = xprime(3,1)*duds + xprime(3,2)*dvds

              d2xds2 = xprime(1,1)*d2uds2 + xprime(1,2)*d2vds2 + xprime(1,3)*(duds**2) + &
                   2.0*xprime(1,4)*duds*dvds + xprime(1,5)*(dvds**2)
              d2yds2 = xprime(2,1)*d2uds2 + xprime(2,2)*d2vds2 + xprime(2,3)*(duds**2) + &
                   2.0*xprime(2,4)*duds*dvds + xprime(2,5)*(dvds**2)
              d2zds2 = xprime(3,1)*d2uds2 + xprime(3,2)*d2vds2 + xprime(3,3)*(duds**2) + &
                   2.0*xprime(3,4)*duds*dvds + xprime(3,5)*(dvds**2)
           ELSE                ! Interpret splinefx as x, splinefy as y, splinefz as z
              CALL EZspline_derivative(XC_spl, 1, s_val,   dxds, ier)
              CALL EZspline_derivative(XC_spl, 2, s_val, d2xds2, ier)
              CALL EZspline_derivative(YC_spl, 1, s_val,   dyds, ier)
              CALL EZspline_derivative(YC_spl, 2, s_val, d2yds2, ier)
              CALL EZspline_derivative(ZC_spl, 1, s_val,   dzds, ier)
              CALL EZspline_derivative(ZC_spl, 2, s_val, d2zds2, ier)
           END IF ! winding surface?

           ! c is the cross-product (dx/ds) x (d2x/ds2)
           cx = dyds*d2zds2 - dzds*d2yds2
           cy = dzds*d2xds2 - dxds*d2zds2
           cz = dxds*d2yds2 - dyds*d2xds2
           norm = dxds**2 + dyds**2 + dzds**2
           kappa = SQRT((cx**2 + cy**2 + cz**2)/(norm**3))
           IF (kappa > maxcurv) THEN
              maxcurv = kappa;  s_max = s_val
           END IF
        END DO !ipt

        ! Clean up
        CALL EZspline_free(XC_spl, ier);  CALL EZspline_free(YC_spl, ier)
        IF (.NOT.lwindsurf) CALL EZspline_free(ZC_spl, ier)
      END SUBROUTINE get_coil_maxcurv

!-----------------------------------------------------------------------
!     Subroutine:    get_coil_self_int
!     Author:        J. Breslau (jbreslau@pppl.gov)
!     Date:          8/28/2017
!     Description:   Computes the number of self-intersections of a 2D curve.
!     Input:         Coil index icoil
!     Output:        Number of intersections
!-----------------------------------------------------------------------
      SUBROUTINE get_coil_self_int(icoil, nselfint)
        USE stellopt_vars
        USE stellopt_targets, ONLY : npts_cself
        USE EZspline_obj
        USE EZspline
        IMPLICIT NONE
        INTRINSIC COUNT
        LOGICAL, EXTERNAL :: uvintersect

        ! Arguments
        INTEGER, INTENT(IN)  :: icoil
        INTEGER, INTENT(OUT) :: nselfint

        ! Constants
        INTEGER, DIMENSION(2), PARAMETER :: ezpbc = (/-1,-1/) ! Periodic boundary condition in EZspline

        ! Local variables
        TYPE(EZspline1_r8)                     :: C_spl
        REAL(rprec), DIMENSION(:), ALLOCATABLE :: uc, vc
        REAL(rprec)                            :: s_val
        INTEGER                                :: ier, iseg, iseg2, nknots

        IF (.NOT.lwindsurf) THEN
           nselfint = -1;  RETURN
        END IF

        ! Handle spline boundary conditions
        nknots = COUNT(coil_splinesx(icoil,:) >= 0)
        IF (coil_type(icoil).eq.'M') THEN !Enforce periodicity for modular coils
           coil_splinefx(icoil,     1) = 0.0d0
           coil_splinefx(icoil,nknots) = 1.0d0
           coil_splinefy(icoil,nknots) = coil_splinefy(icoil,1)
        END IF

        ! Get u,v coords along coil
        ALLOCATE(uc(0:npts_cself-1), vc(0:npts_cself-1))

        CALL EZspline_init(C_spl, nknots, ezpbc, ier)
        C_spl%x1 = coil_splinesx(icoil,1:nknots)
        CALL EZspline_setup(C_spl, coil_splinefx(icoil,1:nknots), ier)
        DO iseg=1,npts_cself
           s_val = REAL(iseg-1)/REAL(npts_cself-1)
           CALL EZspline_interp(C_spl, s_val, uc(iseg-1), ier)
        END DO !iseg
        CALL EZspline_free(C_spl, ier)

        CALL EZspline_init(C_spl, nknots, ezpbc, ier)
        C_spl%x1 = coil_splinesy(icoil,1:nknots)
        CALL EZspline_setup(C_spl, coil_splinefy(icoil,1:nknots), ier)
        DO iseg=1,npts_cself
           s_val = REAL(iseg-1)/REAL(npts_cself-1)
           CALL EZspline_interp(C_spl, s_val, vc(iseg-1), ier)
        END DO !iseg
        CALL EZspline_free(C_spl, ier)

        ! Check for self-intersections
        nselfint = 0
        DO iseg=0,npts_cself-2
           DO iseg2=iseg+2, npts_cself-2
              IF (uvintersect(uc(iseg),  vc(iseg),  uc(iseg+1),  vc(iseg+1), &
                              uc(iseg2), vc(iseg2), uc(iseg2+1), vc(iseg2+1))) &
                nselfint = nselfint + 1
           END DO !iseg2
        END DO !iseg

        DEALLOCATE(uc, vc)
      END SUBROUTINE get_coil_self_int

!-------------------------------------------------------------------------------
! Return TRUE if the two specified line segments in a plane intersect.
! Segments intersect at (u,v) = (u11+alpha*(u12-u11), v11+alpha*(v12-v11))
! Segments intersect at (u,v) = (u21+beta*(u22-u21),  v21+beta*(v22-v21))
! Invert matrix to solve for alpha, beta.
!-------------------------------------------------------------------------------
      LOGICAL FUNCTION uvintersect(u11, v11, u12, v12, u21, v21, u22, v22)
        USE stel_kinds, ONLY : rprec
        IMPLICIT NONE

        REAL(rprec), INTENT(IN) :: u11, v11, u12, v12, u21, v21, u22, v22
        REAL(rprec), PARAMETER :: zero=0.0D0, one=1.0D0
        REAL(rprec) :: det, alpha, beta

        uvintersect = .FALSE.
        det = (u12 - u11)*(v21 - v22) - (v12 - v11)*(u21 - u22)
        IF (det == 0.0) RETURN  ! Parallel segments do not intersect.

        alpha = ((v21 - v22)*(u21 - u11) + (u22 - u21)*(v21 - v11))/det
        IF ((alpha < zero) .OR. (alpha > one)) RETURN

        beta  = ((v11 - v12)*(u21 - u11) + (u12 - u11)*(v21 - v11))/det
        IF ((beta.ge.zero).AND.(beta.le.one)) uvintersect = .TRUE.
      END FUNCTION uvintersect
