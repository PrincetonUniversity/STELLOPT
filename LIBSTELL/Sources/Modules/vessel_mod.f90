!-----------------------------------------------------------------------
!     Module:        vessel_mod
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/03/2012
!     Description:   This module contains routines for working with
!                    vessel files.  It is designed to load a vessel
!                    into memory, calculate the minimum distance between
!                    a point in space and the vessel, and free the
!                    memory.  It also determines if the point is inside
!                    or outside the vessel.  The sequence of calls is:
!                    CALL vessel_load_txt(file,ier)
!                    CALL vessel_dist_cyl(r,phi,z,dist,outside)
!                              - or -
!                    CALL vessel_dist_cart(x,y,z,dist,outside)
!                    CALL vessel_free
!                    The dist variable is the minimum distance between
!                    the point and the vessel.  The outside variable is
!                    an integer which will be positive if the point lies
!                    on the outside of the vessel.
!-----------------------------------------------------------------------
      MODULE vessel_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE safe_open_mod
      USE EZspline_obj
      USE EZspline
      
!-----------------------------------------------------------------------
!     Module Variables
!         
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER   :: nvertex_ves, nu_ves, nv_ves
      REAL ::  phives
      REAL, ALLOCATABLE :: R_ves(:), PHI_ves(:), Z_ves(:), retmap(:), &
                           X_ves(:), Y_ves(:), SNX_ves(:), SNY_ves(:), &
                           SNZ_ves(:)
      CHARACTER(LEN=256) :: machine_string
      CHARACTER(LEN=256) :: date
      TYPE(EZspline2_r8) :: R_ves_spl, Z_ves_spl
      
!-----------------------------------------------------------------------
!     Subroutines
!         vessel_load_txt:   Reads a vessel file and loads the vessel.
!         vessel_dist_cart:  Calculates the distance to the vessel.
!         vessel_dist_cyl:   Calculates the distance to the vessel.
!         vessel_free:   Free allocated variables.
!-----------------------------------------------------------------------
      INTERFACE vessel_dist_cyl
         MODULE PROCEDURE vessel_dist_cyl_flt, vessel_dist_cyl_dbl
      END INTERFACE
      CONTAINS
      
      SUBROUTINE vessel_load_txt(filename,istat)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER, INTENT(inout)       :: istat
      INTEGER :: iunit, u, v, uv, ntheta, nphi, ier
      INTEGER :: bcs0(2) = (/-1,-1/)
      REAL    :: temp1, temp2, temp3, pi2, r1, z1, r2, z2, theta
      DOUBLE PRECISION, ALLOCATABLE :: r_temp(:,:),z_temp(:,:)
      CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
      READ(iunit,'(A)') machine_string
      READ(iunit,'(A)') date
      nvertex_ves = 0
      istat   = 0
      pi2 = 8.0 * ATAN(1.0)
      DO WHILE (istat == 0)
         nvertex_ves = nvertex_ves + 1
         READ(iunit,*,IOSTAT=istat) temp1,temp2,temp3
      END DO
      nvertex_ves = nvertex_ves - 1
      CLOSE(iunit)
      CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
      ALLOCATE(R_ves(nvertex_ves),PHI_ves(nvertex_ves),Z_ves(nvertex_ves), &
               X_ves(nvertex_ves),Y_ves(nvertex_ves),retmap(nvertex_ves), &
               SNX_ves(nvertex_ves), SNY_ves(nvertex_ves), SNZ_ves(nvertex_ves))
      READ(iunit,'(A)') machine_string
      READ(iunit,'(A)') date
      istat   = 0
      DO u = 1, nvertex_ves
         READ(iunit,*,IOSTAT=istat) R_ves(u), Z_ves(u), PHI_ves(u)
      END DO
      CLOSE(iunit)
      X_ves = R_ves * cos (PHI_ves)
      Y_ves = R_ves * sin (PHI_ves)
      phives = MAXVAL(PHI_ves)
      nu_ves = MINLOC(PHI_ves, DIM=1, MASK = PHI_ves > 0.0 ) - 1
      nv_ves = nvertex_ves / nu_ves
      uv = 1
      DO v = 1, nv_ves - 1
         DO u = 1, nu_ves - 1
            r1 = R_ves(uv+1) - R_ves(uv)
            z1 = Z_ves(uv+1) - Z_ves(uv)
            theta = ATAN2(z1,r1)
            IF (theta < 0) theta = theta + pi2
            theta = theta - 0.25*pi2
            r2 = SQRT(r1*r1+z1*z1) * cos(theta)
            z2 = SQRT(r1*r1+z1*z1) * sin(theta)
            SNX_ves(uv) = r2 * cos(PHI_ves(uv))
            SNY_ves(uv) = r2 * sin(PHI_ves(uv))
            SNZ_ves(uv) = z2
            uv = uv + 1
         END DO
         ! Handles last poloidal point (mirrored point)
         SNX_ves(uv) = SNX_ves(uv-nu_ves)
         SNY_ves(uv) = SNY_ves(uv-nu_ves)
         SNZ_ves(uv) = SNZ_ves(uv-nu_ves)
         uv = uv + 1       
      END DO
      ! Handle the last toroidal point
      DO u = 1, nu_ves - 1
         r1 = R_ves(uv+1) - R_ves(uv)
         z1 = Z_ves(uv+1) - Z_ves(uv)
         theta = ATAN2(z1,r1)
         IF (theta < 0) theta = theta + pi2
         theta = theta - 0.25*pi2
         r2 = SQRT(r1*r1+z1*z1) * cos(theta)
         z2 = SQRT(r1*r1+z1*z1) * sin(theta)
         SNX_ves(uv) = r2 * cos(PHI_ves(uv))
         SNY_ves(uv) = r2 * sin(PHI_ves(uv))
         SNZ_ves(uv) = z2
         uv = uv + 1
      END DO
      ! Handle the very last point
      SNX_ves(uv) = SNX_ves(uv-nu_ves)
      SNY_ves(uv) = SNY_ves(uv-nu_ves)
      SNZ_ves(uv) = SNZ_ves(uv-nu_ves)
      ! Construct Splines
      !ntheta=COUNT(phi == 0.0)
      !nphi=nvertex_ves/ntheta
      ! Construct Spline represenation
      ALLOCATE(r_temp(nu_ves,nv_ves),z_temp(nu_ves,nv_ves))
      uv = 1
      DO v = 1, nv_ves
         DO u = 1, nu_ves
            r_temp(u,v) = R_ves(uv)
            z_temp(u,v) = Z_ves(uv)
            uv = uv + 1
         END DO
      END DO
      ier = 0
      IF (EZspline_allocated(R_ves_spl)) CALL EZspline_free(R_ves_spl,ier)
      IF (EZspline_allocated(Z_ves_spl)) CALL EZspline_free(Z_ves_spl,ier)
      CALL EZspline_init(R_ves_spl,nu_ves,nv_ves,bcs0,bcs0,ier)
      CALL EZspline_init(Z_ves_spl,nu_ves,nv_ves,bcs0,bcs0,ier)
      R_ves_spl%isHermite = 1
      Z_ves_spl%isHermite = 1
      CALL EZspline_setup(R_ves_spl,r_temp,ier)
      CALL EZspline_setup(Z_ves_spl,z_temp,ier)
      DEALLOCATE(r_temp,z_temp)
      
      RETURN
      END SUBROUTINE vessel_load_txt
      
      SUBROUTINE vessel_dist_cyl_flt(rp,phip,zp,dist,out)
      IMPLICIT NONE
      REAL, INTENT(in) :: rp, phip, zp
      REAL, INTENT(out) :: dist
      INTEGER, INTENT(out) :: out
      REAL :: xp, yp, phi_temp
      phi_temp = MOD(phip,phives)
      xp = rp * cos(phi_temp)
      yp = rp * sin(phi_temp)
      CALL vessel_dist_cart(xp,yp,zp,dist,out)
      RETURN
      END SUBROUTINE vessel_dist_cyl_flt
      
      SUBROUTINE vessel_dist_cyl_dbl(rp,phip,zp,dist,out)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: rp, phip, zp
      DOUBLE PRECISION, INTENT(out) :: dist
      INTEGER, INTENT(out) :: out
      REAL :: xp_flt, yp_flt, zp_flt, dist_flt, phi_temp
      phi_temp = MOD(REAL(phip),phives)
      xp_flt = rp * cos(phi_temp)
      yp_flt = rp * sin(phi_temp)
      zp_flt = zp
      CALL vessel_dist_cart(xp_flt,yp_flt,zp_flt,dist_flt,out)
      dist   = dist_flt
      RETURN
      END SUBROUTINE vessel_dist_cyl_dbl
      
      SUBROUTINE vessel_dist_cart(xp,yp,zp,dist,out)
      IMPLICIT NONE
      REAL, INTENT(in) :: xp, yp, zp
      REAL, INTENT(out) :: dist
      INTEGER, INTENT(out) :: out
      INTEGER :: loc
      REAL :: dx, dy, dz
      out = 0
      retmap(:) = SQRT( (xp-X_ves(:))**2 + (yp-Y_ves(:))**2 + (zp-Z_ves(:))**2 )
      loc = MINLOC(retmap,DIM=1)
      dist = retmap(loc)
      dx   = xp - X_ves(loc)
      dy   = yp - Y_ves(loc)
      dz   = zp - Z_ves(loc)
      out = SIGN(1.0,(dx*SNX_ves(loc))+(dy*SNY_ves(loc))+(dz*SNZ_ves(loc)))
      RETURN
      END SUBROUTINE vessel_dist_cart
      
      SUBROUTINE vessel_info(iunit)
      IMPLICIT NONE
      INTEGER, INTENT(in)    :: iunit
      WRITE(iunit,'(A)')         ' -----  Vessel Information  -----'
      WRITE(iunit,'(3X,A,A)')    'Vessel Name: ',TRIM(machine_string(10:))
      WRITE(iunit,'(3X,A,A)')    'Date       : ',TRIM(date(6:))
      WRITE(iunit,'(3X,A,I7)')   'Elements   : ',nvertex_ves
      WRITE(iunit,'(3X,A,F7.3)') 'Phi_ves    : ',phives
      WRITE(iunit,'(3X,A,I7)')   'NU         : ',nu_ves
      WRITE(iunit,'(3X,A,I7)')   'NV         : ',nv_ves
      RETURN
      END SUBROUTINE vessel_info
      
      SUBROUTINE vessel_free(istat)
      IMPLICIT NONE
      INTEGER, INTENT(inout) :: istat
      INTEGER :: ier
      DEALLOCATE(R_ves, PHI_ves, Z_ves, X_ves, Y_ves, retmap, STAT=istat)
      IF (EZspline_allocated(R_ves_spl)) CALL EZspline_free(R_ves_spl,ier)
      IF (EZspline_allocated(Z_ves_spl)) CALL EZspline_free(Z_ves_spl,ier)
      machine_string=''
      date=''
      nvertex_ves = -1
      RETURN
      END SUBROUTINE vessel_free

!-----------------------------------------------------------------------
!     End Module
!-----------------------------------------------------------------------
      END MODULE vessel_mod
