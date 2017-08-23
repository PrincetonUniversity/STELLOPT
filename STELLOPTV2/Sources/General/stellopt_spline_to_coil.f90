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

      ! Make sure winding surface structure is set up correctly
      IF (lwindsurf) THEN
         IF ((windsurf%mmax.lt.0).or.(windsurf%nmax.lt.0).or.&
              (.not.allocated(windsurf%rctab)).or.(.not.allocated(windsurf%zstab))) &
              CALL handle_err(BAD_CWS_ERR, 'stellopt_spline_to_coil', 1)
         windsurf%nfp = nfp
      END IF

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

        ! Local variables
        TYPE(EZspline1_r8) :: XC_spl, YC_spl, ZC_spl
        REAL(rprec)        :: s_val, u, v, xlast, ylast
        INTEGER            :: nknots, ier, j

        lmod = .FALSE.

        ! Set up splines
        nknots = COUNT(coil_splinesx(icoil,:) >= 0)
        CALL EZspline_init(XC_spl, nknots, bcs1, ier)
        XC_spl%x1 = coil_splinesx(icoil,1:nknots)
        XC_spl%isHermite = 1
        CALL EZspline_setup(XC_spl, coil_splinefx(icoil,1:nknots), ier)
        xlast = coil_splinefx(icoil,nknots)

        nknots = COUNT(coil_splinesy(icoil,:) >= 0)
        CALL EZspline_init(YC_spl, nknots, bcs1, ier)
        YC_spl%x1 = coil_splinesy(icoil,1:nknots)
        YC_spl%isHermite = 1
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
           nknots = COUNT(coil_splinesz(icoil,:) >= 0)
           CALL EZspline_init(ZC_spl, nknots, bcs1, ier)
           ZC_spl%x1 = coil_splinesz(icoil,1:nknots)
           ZC_spl%isHermite = 1
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
!     Subroutine:    stellopt_uv_to_xyz
!     Author:        J. Breslau (jbreslau@pppl.gov)
!     Date:          8/21/2017
!     Description:   This subroutine converts a u,v coordinate pair
!                    on a winding surface to Cartesian coordinates.
!     Inputs:
!                    Poloidal coordinate 0 < u < 1.
!                    Toroidal coordinate 0 < v < nfp.
!     Outputs:
!                    Cartesian values (x, y, z).
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_uv_to_xyz(u, v, x, y, z)
        USE stel_kinds, ONLY    : rprec
        USE stellopt_vars, ONLY : windsurf
        IMPLICIT NONE
        INTRINSIC COS, SIN

        ! Arguments
        REAL(rprec), INTENT(IN)  :: u, v
        REAL(rprec), INTENT(OUT) :: x, y, z

        ! Constants
        REAL(rprec), PARAMETER :: twopi = 6.283185307179586476925286766559D0

        ! Local variables
        REAL(rprec) :: theta, Nphi         ! Winding surface angular coordinates, in radians
        REAL(rprec) :: cmt, smt            ! cos(m theta), sin(m theta)
        REAL(rprec) :: cap, cam, sap, sam  ! cos(m theta +/- n N phi), sin(m theta +/- n N phi)
        REAL(rprec) :: R                   ! Cylindrical R coordinate
        REAL(rprec) :: Cu, Su, Cv, Sv, tmp
        INTEGER     :: m, n                ! Poloidal, toroidal mode numbers

        ! Initialize N phi trig. functions
        Nphi = twopi*v
        Cv = COS(Nphi);  Sv = SIN(Nphi)

        ! Initialize m theta trig. fcns, recurrence relation
        theta = twopi*u
        Cu = COS(theta);  Su = SIN(theta)
        cmt = 1.0;  smt = 0.0

        ! Outer loop over poloidal mode numbers
        R = 0.0;  z = 0.0
        DO m=0,windsurf%mmax
           ! Initialize +/- N phi recurrence relations
           cam = cmt;  cap = cmt;  sam = smt;  sap = smt

           ! n=0 mode contributions
           R = R + windsurf%rctab(0,m)*cap
           z = z + windsurf%zstab(0,m)*sap

           ! Inner loop over toroidal mode numbers
           DO n=1,windsurf%nmax
              ! Increment/decrement sin(arg), cos(arg)
              tmp = cap*Cv - sap*Sv;  sap = sap*Cv + cap*Sv;  cap = tmp
              tmp = cam*Cv + sam*Sv;  sam = sam*Cv - cam*Sv;  cam = tmp

              ! Add terms to series
              R = R + windsurf%rctab(n,m)*cap + windsurf%rctab(-n,m)*cam
              z = z + windsurf%zstab(n,m)*sap + windsurf%zstab(-n,m)*sam
           END DO !n

           ! Increment sin(m theta), cos(m theta)
           tmp = cmt*Cu - smt*Su;  smt = smt*Cu + cmt*Su;  cmt = tmp
        END DO !m

        ! Convert from cylindrical to Cartesian coordinates
        x = R*COS(Nphi/windsurf%nfp);  y = R*SIN(Nphi/windsurf%nfp)
      END SUBROUTINE stellopt_uv_to_xyz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE get_coil_length(icoil, length)
        USE stel_kinds, ONLY    : rprec !TMP
        IMPLICIT NONE
        INTRINSIC SQRT

        ! Arguments
        INTEGER, INTENT(IN)      :: icoil
        REAL(rprec), INTENT(OUT) :: length

        ! Constants
        INTEGER, PARAMETER :: nseg_len = 128

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
