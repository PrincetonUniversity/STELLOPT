!-----------------------------------------------------------------------
!     Subroutine:    stellopt_spline_to_coil
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/14/2016
!     Description:   This subroutine generates a coil data structure from
!                    the stellopt splines and writes it to a file.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_spline_to_coil(nseg, fixedcoils, lscreen)
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
      INTEGER, INTENT(in)           :: nseg
      CHARACTER(*), INTENT(in)      :: fixedcoils
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
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (lscreen) WRITE(6,*) '---------------------------  WRITING COIL  ------------------------'
      IF (nseg.lt.3) CALL handle_err(0, 'stellopt_spline_to_coil: nseg', nseg)
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

      ! Append fixed coil data, if present
      IF (LEN_TRIM(fixedcoils).GT.0) &
         CALL append_fixed_coils(iunit, fixedcoils, n_coilgroups, lscreen)

      ! Close file
      WRITE(iunit,'(A)') 'end'
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
        USE stellopt_cbspline
        IMPLICIT NONE

        ! Arguments
        INTEGER, INTENT(IN)                       :: icoil, nseg
        REAL(rprec), DIMENSION(nseg), INTENT(OUT) :: xarr, yarr, zarr
        LOGICAL, INTENT(OUT)                      :: lmod  ! Does the coil close on itself within a field period?
                                                           ! Note: this is true for saddle coils as well.
        ! Local variables
        TYPE(cbspline)     :: XC_spl, YC_spl, ZC_spl
        REAL(rprec)        :: s_val, u, v, v0
        INTEGER            :: nknots, ncoefs, ier, j, js1, js2

        ! Handle spline boundary conditions
        CALL enforce_spline_bcs(icoil, lmod)
        ncoefs = coil_nctrl(icoil)
        nknots = ncoefs + 4

        ! Set up splines
        CALL cbspline_init(XC_spl, ncoefs, ier)
        CALL cbspline_setup(XC_spl, ncoefs, &
             coil_splinesx(icoil,1:nknots), coil_splinefx(icoil,1:ncoefs), ier)

        CALL cbspline_init(YC_spl, ncoefs, ier)
        CALL cbspline_setup(YC_spl, ncoefs, &
             coil_splinesy(icoil,1:nknots), coil_splinefy(icoil,1:ncoefs), ier)

        IF (lwindsurf) THEN !Interpret coil_splinefx as u, coil_splinefy as v
           windsurf%nfp = nfp

           ! Evaluate the coil
           v0 = coil_splinefy(icoil,1)
           js1 =   FLOOR(1.0 + coil_splinesx(icoil,1)*(nseg - 1))
           js2 = CEILING(1.0 + coil_splinesx(icoil,nknots)*(nseg - 1))
           DO j = 1, js1  ! Straight section above midplane
              u = REAL(j-1)/REAL(nseg-1)
              CALL stellopt_uv_to_xyz(u, v0, xarr(j), yarr(j), zarr(j))
           END DO
           DO j = js1+1, js2-1
              s_val = REAL(j-1)/REAL(nseg-1)
              CALL cbspline_eval(XC_spl, s_val, u, ier)
              CALL cbspline_eval(YC_spl, s_val, v, ier)
              CALL stellopt_uv_to_xyz(u, v, xarr(j), yarr(j), zarr(j))
           END DO
           DO j = js2, nseg  ! Straight section below midplane
              u = REAL(j-1)/REAL(nseg-1)
              CALL stellopt_uv_to_xyz(u, v0, xarr(j), yarr(j), zarr(j))
           END DO
        ELSE  !Interpret coil_splinefx as x, coil_splinefy as y, coil_splinefz as z.
           CALL cbspline_init(ZC_spl, ncoefs, ier)
           CALL cbspline_setup(ZC_spl, ncoefs, &
                coil_splinesz(icoil,1:nknots), coil_splinefz(icoil,1:ncoefs), ier)

           ! Evaluate the coil
           DO j = 1, nseg
              s_val = REAL(j-1)/REAL(nseg-1)
              CALL cbspline_eval(XC_spl, s_val, xarr(j), ier)
              CALL cbspline_eval(YC_spl, s_val, yarr(j), ier)
              CALL cbspline_eval(ZC_spl, s_val, zarr(j), ier)
           END DO

           CALL cbspline_delete(ZC_spl)
        END IF

        CALL cbspline_delete(XC_spl);  CALL cbspline_delete(YC_spl)
      END SUBROUTINE spline_to_coil

!-----------------------------------------------------------------------
      SUBROUTINE enforce_spline_bcs(icoil, isper)
        USE stel_kinds, ONLY: rprec
        USE stellopt_vars, ONLY : coil_splinefx, coil_splinefy, coil_splinefz, &
             coil_splinesx, coil_nctrl, coil_type, lwindsurf
        USE vmec_input,  ONLY : nfp
        IMPLICIT NONE
        INTRINSIC ABS, MODULO

        INTEGER, INTENT(IN)  :: icoil
        LOGICAL, INTENT(OUT) :: isper  !Does the coil repeat each field period?

        REAL(rprec), PARAMETER :: zero=0.0d0, one=1.0d0, abstol=2.0d-15
        REAL(rprec) krat
        INTEGER ncoefs
        ncoefs = coil_nctrl(icoil)

        SELECT CASE (coil_type(icoil))
        CASE ('M') ! Modular coil
           isper = .TRUE. ! Modular coils repeat each field period.
           coil_splinefy(icoil,ncoefs) = coil_splinefy(icoil,1)
           IF (lwindsurf) THEN !u0=k0,uf=kf,vf=v0
              coil_splinefx(icoil,1)      = coil_splinesx(icoil,1)
              coil_splinefx(icoil,ncoefs) = coil_splinesx(icoil,ncoefs+4)
           ELSE  !z0=z(t=0), zf=z0,xf=x0,yf=y0
              coil_splinefx(icoil,ncoefs) = coil_splinefx(icoil,1)
              coil_splinefz(icoil,ncoefs) = coil_splinefz(icoil,1)
           END IF
        CASE ('A') ! Modular coil with straight section represented by spline
           isper = .TRUE. ! Modular coils repeat each field period.
           IF (lwindsurf) THEN
              coil_splinefy(icoil,1) = coil_splinefy(icoil,4)
              coil_splinefy(icoil,2) = coil_splinefy(icoil,4)
              coil_splinefy(icoil,3) = coil_splinefy(icoil,4)
              coil_splinefy(icoil,ncoefs) = coil_splinefy(icoil,4)
              coil_splinefy(icoil,ncoefs-1) = coil_splinefy(icoil,4)
              coil_splinefy(icoil,ncoefs-2) = coil_splinefy(icoil,4)
              coil_splinefy(icoil,ncoefs-3) = coil_splinefy(icoil,4)
              RETURN
           ELSE  !Just treat this the same as type 'M' for now...
              coil_splinefx(icoil,ncoefs) = coil_splinefx(icoil,1)
              coil_splinefy(icoil,ncoefs) = coil_splinefy(icoil,1)
              coil_splinefz(icoil,ncoefs) = coil_splinefz(icoil,1)
          END IF
        CASE ('S') ! Saddle coil
           isper = .TRUE. ! Saddle coils repeat each field period like modular coils.
           coil_splinefx(icoil,ncoefs) = coil_splinefx(icoil,1)
           coil_splinefy(icoil,ncoefs) = coil_splinefy(icoil,1)
           IF (.NOT.lwindsurf) coil_splinefz(icoil,ncoefs) = coil_splinefz(icoil,1)
        CASE ('P') ! Wavy PF coil
           isper = .FALSE. ! PF coils wrap around, do not repeat.
           coil_splinefx(icoil,ncoefs) = coil_splinefx(icoil,1)
           IF (lwindsurf) THEN !v0=0,vf=N,uf=u0
              coil_splinefy(icoil,1) = zero;  coil_splinefy(icoil,ncoefs) = nfp
           ELSE
              coil_splinefy(icoil,ncoefs) = coil_splinefy(icoil,1)
              coil_splinefz(icoil,ncoefs) = coil_splinefz(icoil,1)
           END IF
        CASE DEFAULT
           IF (lwindsurf) THEN
              isper = ((ABS(MODULO(coil_splinefx(icoil,1),one) - &
                   MODULO(coil_splinefx(icoil,ncoefs),one)).LE.abstol) &
                    .AND.(coil_splinefy(icoil,1) == coil_splinefy(icoil,ncoefs)) )
           ELSE
              isper = ((coil_splinefx(icoil,1) == coil_splinefx(icoil,ncoefs)) &
                  .AND.(coil_splinefy(icoil,1) == coil_splinefy(icoil,ncoefs)) &
                  .AND.(coil_splinefz(icoil,1) == coil_splinefz(icoil,ncoefs)))
           END IF
        END SELECT

        !Enforce continuity of 1st derivatives at s=0
        krat = (coil_splinesx(icoil,ncoefs+1) - coil_splinesx(icoil,ncoefs))/&
             (coil_splinesx(icoil,5) - coil_splinesx(icoil,4))
        coil_splinefx(icoil,ncoefs-1) = coil_splinefx(icoil,ncoefs) - &
             krat*(coil_splinefx(icoil,2) - coil_splinefx(icoil,1))
        coil_splinefy(icoil,ncoefs-1) = coil_splinefy(icoil,ncoefs) - &
             krat*(coil_splinefy(icoil,2) - coil_splinefy(icoil,1))
        IF (.NOT.lwindsurf) &
           coil_splinefz(icoil,ncoefs-1) = coil_splinefz(icoil,ncoefs) - &
                krat*(coil_splinefz(icoil,2) - coil_splinefz(icoil,1))
      END SUBROUTINE enforce_spline_bcs

!-----------------------------------------------------------------------
!     Subroutine:    get_coil_length
!     Author:        J. Breslau (jbreslau@pppl.gov)
!     Date:          8/22/2017
!     Description:   Computes the length of a coil.
!                    Also computes normalized standard dev. of seg length.
!     Input:         Coil index icoil
!     Output:        Length of coil, (segment length s.d.)/(mean length)
!-----------------------------------------------------------------------
      SUBROUTINE get_coil_length(icoil, length, relvar)
        USE stel_kinds, ONLY : rprec
        USE stellopt_targets, ONLY : npts_clen
        IMPLICIT NONE
        INTRINSIC SQRT

        ! Arguments
        INTEGER, INTENT(IN)      :: icoil
        REAL(rprec), INTENT(OUT) :: length, relvar

        ! Constants
        REAL(rprec), PARAMETER :: one=1.0d0

        ! Local variables
        REAL(rprec), DIMENSION(:), ALLOCATABLE :: xc, yc, zc
        REAL(rprec)                            :: seglen
        INTEGER                                :: iseg
        LOGICAL                                :: ldum

        ! Allocate x,y,z storage
        ALLOCATE(xc(npts_clen), yc(npts_clen), zc(npts_clen))

        ! Get x,y,z coords along coil
        CALL spline_to_coil(icoil, npts_clen, xc, yc, zc, ldum)

        ! Compute length
        length = 0.0;  relvar = 0.0
        DO iseg=1,npts_clen-1  ! Assumes xc(nseg) == xc(1), etc.
           seglen = SQRT((xc(iseg+1) - xc(iseg))**2 + &
                         (yc(iseg+1) - yc(iseg))**2 + &
                         (zc(iseg+1) - zc(iseg))**2)
           length = length + seglen
           relvar = relvar + seglen**2
        END DO !iseg
        relvar = SQRT((npts_clen-1)*relvar/(length**2) - one)

        DEALLOCATE(xc, yc, zc)
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
      SUBROUTINE get_coil_sep(x1, y1, z1, n1, x2, y2, z2, n2, minsep, s1, s2)
        USE stel_kinds, ONLY : rprec        
        IMPLICIT NONE
        INTRINSIC SQRT

        INTEGER, INTENT(IN)                    :: n1, n2
        REAL(rprec), DIMENSION(n1), INTENT(IN) :: x1, y1, z1
        REAL(rprec), DIMENSION(n2), INTENT(IN) :: x2, y2, z2
        REAL(rprec), INTENT(OUT)               :: minsep, s1, s2

        REAL(rprec) :: dsq
        INTEGER     :: j1, j2, sl1, sl2

        minsep = 1.0D+60
        sl1 = 1;  sl2 = 1

        DO j1=1,n1-1     ! Assumes 1st and last pts are identical.
           DO j2=1,n2-1
              dsq = (x1(j1) - x2(j2))**2 + (y1(j1) - y2(j2))**2 + (z1(j1) - z2(j2))**2
              IF (dsq < minsep) THEN
                 minsep = dsq
                 sl1 = j1;  sl2 = j2
              END IF
           END DO !j2
        END DO !j1

        minsep = SQRT(minsep)
        s1 = REAL(sl1 - 1)/REAL(n1 - 1)
        s2 = REAL(sl2 - 1)/REAL(n2 - 1)
      END SUBROUTINE get_coil_sep

!-------------------------------------------------------------------------------
!     Subroutine:    get_coil_maxcurv
!     Author:        J. Breslau (jbreslau@pppl.gov)
!     Date:          8/24/2017
!     Description:   Computes the curvature at several points along a selected
!                    coil and returns the maximum value and its location.
!     Input:         Coil index icoil
!     Outputs:       Max. curvature, parametric position & u,v of max. curvature.
!-------------------------------------------------------------------------------
      SUBROUTINE get_coil_maxcurv(icoil, maxcurv, s_max, u_max, v_max, uout)
        USE stellopt_vars
        USE windingsurface
        USE stellopt_targets, ONLY : npts_curv
        USE vmec_input,  ONLY : nfp
        USE stellopt_cbspline
        IMPLICIT NONE
        INTRINSIC SQRT

        ! Arguments
        INTEGER, INTENT(IN)           :: icoil
        REAL(rprec), INTENT(OUT)      :: maxcurv, s_max, u_max, v_max
        INTEGER, INTENT(IN), OPTIONAL :: uout

        ! Local variables
        TYPE(cbspline)     :: XC_spl, YC_spl, ZC_spl
        REAL(rprec)        :: xprime(3,5)
        REAL(rprec)        :: s_val, u, duds, d2uds2, v, dvds, d2vds2
        REAL(rprec)        :: dxds, dyds, dzds, d2xds2, d2yds2, d2zds2
        REAL(rprec)        :: cx, cy, cz, norm, kappa
        INTEGER            :: ier, ipt, nknots, ncoefs
        LOGICAL            :: ldum

        ! Handle spline boundary conditions
        CALL enforce_spline_bcs(icoil, ldum)
        ncoefs = coil_nctrl(icoil)
        nknots = ncoefs + 4

        ! Set up splines -- cannot use spline_to_coil b/c we also need derivatives
        CALL cbspline_init(XC_spl, ncoefs, ier)
        CALL cbspline_setup(XC_spl, ncoefs, &
             coil_splinesx(icoil,1:nknots), coil_splinefx(icoil,1:ncoefs), ier)

        CALL cbspline_init(YC_spl, ncoefs, ier)
        CALL cbspline_setup(YC_spl, ncoefs, &
             coil_splinesy(icoil,1:nknots), coil_splinefy(icoil,1:ncoefs), ier)

        IF (lwindsurf) THEN
           windsurf%nfp = nfp
        ELSE
           CALL cbspline_init(ZC_spl, ncoefs, ier)
           CALL cbspline_setup(ZC_spl, ncoefs, &
                coil_splinesz(icoil,1:nknots), coil_splinefz(icoil,1:ncoefs), ier)
        END IF

        ! Write header to optional output file
        IF (PRESENT(uout)) WRITE(uout,'(A)') 's  K  u  v  x  y  z'

        ! Find maximum curvature along coil
        maxcurv = -1.0;  s_max = -1.0;  u_max = -1.0;  v_max = -1.0
        DO ipt=1,npts_curv
           s_val = REAL(ipt-1)/REAL(npts_curv)

           ! We need 1st and 2nd derivatives of (x,y,z) w/r/t s to compute curvature:
           IF (lwindsurf) THEN ! Interpret splinefx as u, splinefy as v
              ! First get 1st and 2nd derivs of u,v w/r/t s...
              IF ((s_val.LT.coil_splinesx(icoil,1)) .OR. &
                  (s_val.GT.coil_splinesx(icoil,nknots))) THEN  ! Straight section
                 u = s_val;  duds = 1.0;  d2uds2 = 0.0
                 v = coil_splinefy(icoil,1);  dvds = 0.0;  d2vds2 = 0.0
              ELSE
                 CALL cbspline_derivs(XC_spl, s_val, u, duds, d2uds2, ier)
                 CALL cbspline_derivs(YC_spl, s_val, v, dvds, d2vds2, ier)
              END IF

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
              CALL cbspline_derivs(XC_spl, s_val, u, dxds, d2xds2, ier)
              CALL cbspline_derivs(YC_spl, s_val, u, dyds, d2yds2, ier)
              CALL cbspline_derivs(ZC_spl, s_val, u, dzds, d2zds2, ier)
           END IF ! winding surface?

           ! c is the cross-product (dx/ds) x (d2x/ds2)
           cx = dyds*d2zds2 - dzds*d2yds2
           cy = dzds*d2xds2 - dxds*d2zds2
           cz = dxds*d2yds2 - dyds*d2xds2
           norm = dxds**2 + dyds**2 + dzds**2
           kappa = SQRT((cx**2 + cy**2 + cz**2)/(norm**3))
           IF (PRESENT(uout)) THEN  ! Write curvature to optional output file
              CALL stellopt_uv_to_xyz(u, v, xprime(1,1), xprime(2,1), xprime(3,1))
              WRITE(uout,'(7ES17.8E2)') s_val,kappa,u,v,xprime(1,1),xprime(2,1),xprime(3,1)
           END IF
           IF (kappa > maxcurv) THEN
              maxcurv = kappa;  s_max = s_val
              u_max = u;  v_max = v
           END IF
        END DO !ipt

        ! Clean up
        CALL cbspline_delete(XC_spl);  CALL cbspline_delete(YC_spl)
        IF (.NOT.lwindsurf) CALL cbspline_delete(ZC_spl)
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
        USE stellopt_cbspline
        IMPLICIT NONE
        LOGICAL, EXTERNAL :: uvintersect

        ! Arguments
        INTEGER, INTENT(IN)  :: icoil
        INTEGER, INTENT(OUT) :: nselfint

        ! Local variables
        TYPE(cbspline)                         :: C_spl
        REAL(rprec), DIMENSION(:), ALLOCATABLE :: uc, vc
        REAL(rprec)                            :: s_val
        INTEGER                                :: ier, iseg, iseg2, nknots, ncoefs, js1, js2
        LOGICAL                                :: ldum

        IF (.NOT.lwindsurf) THEN
           nselfint = -1;  RETURN
        END IF

        ! Handle spline boundary conditions
        CALL enforce_spline_bcs(icoil, ldum)
        ncoefs = coil_nctrl(icoil)
        nknots = ncoefs + 4

        ! Get u,v coords along coil -- don't use spline_to_coil b/c we don't need x,y,z coords
        ALLOCATE(uc(npts_cself), vc(npts_cself))

        js1 =   FLOOR(1.0 + coil_splinesx(icoil,1)*(npts_cself - 1))
        js2 = CEILING(1.0 + coil_splinesx(icoil,nknots)*(npts_cself - 1))

        CALL cbspline_init(C_spl, ncoefs, ier)
        CALL cbspline_setup(C_spl, ncoefs, &
             coil_splinesx(icoil,1:nknots), coil_splinefx(icoil,1:ncoefs), ier)
        DO iseg = 1, js1  ! Straight section above midplane
           uc(iseg) = REAL(iseg-1)/REAL(npts_cself-1)
        END DO
        DO iseg = js1+1, js2-1
           s_val = REAL(iseg-1)/REAL(npts_cself-1)
           CALL cbspline_eval(C_spl, s_val, uc(iseg), ier)
        END DO !iseg
        DO iseg = js2, npts_cself  ! Straight section below midplane
           uc(iseg) = REAL(iseg-1)/REAL(npts_cself-1)
        END DO

        CALL cbspline_setup(C_spl, ncoefs, &
             coil_splinesy(icoil,1:nknots), coil_splinefy(icoil,1:ncoefs), ier)
        DO iseg = js1+1, js2-1
           s_val = REAL(iseg-1)/REAL(npts_cself-1)
           CALL cbspline_eval(C_spl, s_val, vc(iseg), ier)
        END DO !iseg
        vc(1:js1) = coil_splinefy(icoil,1);  vc(js2:npts_cself) = coil_splinefy(icoil,1)
        CALL cbspline_delete(C_spl)

        ! Check for self-intersections
        nselfint = 0
        DO iseg=1,npts_cself-1
           DO iseg2=iseg+2, npts_cself-1
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

!-----------------------------------------------------------------------
!     Subroutine:    get_coil_tor_excur
!     Author:        J. Breslau (jbreslau@pppl.gov)
!     Date:          1/11/2019
!     Description:   Computes the weighted RMS excursion of a coil from a
!                    constant-phi plane.
!     Input:         Coil index icoil, theta weighting coefficient
!     Output:        RMS excursion in radians
!-----------------------------------------------------------------------
      SUBROUTINE get_coil_tor_excur(icoil, nu, excur)
        USE stel_kinds, ONLY : rprec
        USE stellopt_targets, ONLY : npts_torx
        USE stellopt_vars
        USE windingsurface
        USE vmec_input,  ONLY : nfp
        USE stellopt_cbspline
        IMPLICIT NONE
        INTRINSIC ATAN2, COS, SQRT

        ! Arguments
        INTEGER, INTENT(IN)      :: icoil  ! Coil index
        REAL(rprec), INTENT(IN)  :: nu     ! Theta non-uniformity of penalty weight
        REAL(rprec), INTENT(OUT) :: excur  ! RMS excursion

        ! Constants
        REAL(rprec), PARAMETER :: zero=0.0D0, pi2r = 1.5707963267948966192313216916398_rprec

        ! Local variables
        TYPE(cbspline)                         :: XC_spl, YC_spl, ZC_spl
        REAL(rprec), DIMENSION(:), ALLOCATABLE :: dl, wt, va
        REAL(rprec)                            :: s_val, u00, u0, u1, v00, v0, v1
        REAL(rprec)                            :: x00, y00, z00, x0, y0, z0, x1, y1, z1
        REAL(rprec)                            :: hsgn, len, vbar, vvar
        INTEGER                                :: nknots, ncoefs, ier, j
        LOGICAL                                :: lmod

        IF ((coil_type(icoil).NE.'M').AND.(coil_type(icoil).NE.'A')) THEN
           excur = zero
           RETURN
        END IF

        ! Allocate storage
        ALLOCATE(dl(npts_torx), wt(npts_torx), va(npts_torx))

        ! Handle spline boundary conditions
        CALL enforce_spline_bcs(icoil, lmod)
        ncoefs = coil_nctrl(icoil)
        nknots = ncoefs + 4

        ! Set up splines
        CALL cbspline_init(XC_spl, ncoefs, ier)
        CALL cbspline_setup(XC_spl, ncoefs, &
             coil_splinesx(icoil,1:nknots), coil_splinefx(icoil,1:ncoefs), ier)

        CALL cbspline_init(YC_spl, ncoefs, ier)
        CALL cbspline_setup(YC_spl, ncoefs, &
             coil_splinesy(icoil,1:nknots), coil_splinefy(icoil,1:ncoefs), ier)

        len = zero;  vbar = zero
        IF (lwindsurf) THEN !Interpret coil_splinefx as u, coil_splinefy as v
           windsurf%nfp = nfp

           u0 = zero;  v0 = coil_splinefy(icoil,1)
           CALL stellopt_uv_to_xyz(u0, v0, x0, y0, z0)
           u00 = u0; v00 = v0;  x00 = x0; y00 = y0; z00 = z0
           DO j = 2, npts_torx
              s_val = REAL(j-1)/REAL(npts_torx-1)
              IF ((s_val.GT.coil_splinesx(icoil,1)) .AND. (s_val.LT.coil_splinesx(icoil,nknots))) THEN
                 CALL cbspline_eval(XC_spl, s_val, u1, ier)
                 CALL cbspline_eval(YC_spl, s_val, v1, ier)
              ELSE
                 u1 = s_val
                 v1 = coil_splinefy(icoil,1)
              END IF
              CALL stellopt_uv_to_xyz(u1, v1, x1, y1, z1)
              dl(j) = SQRT((x1 - x0)**2 + (y1 - y0)**2 + (z1 - z0)**2)
              x0 = x1;  y0 = y1;  z0 = z1
              len = len + dl(j)
              va(j) = 0.5d0*(v0 + v1)
              vbar = vbar + dl(j)*va(j)
              wt(j) = (1.0d0 - nu) + nu*COS(pi2r*(u0 + u1))**4
              u0 = u1;  v0 = v1
           END DO !j

           va(1) = 0.5d0*(v00 + v0)
           wt(1) = (1.0d0 - nu) + nu*COS(pi2r*(u00 + u0))**4
        ELSE  !Interpret coil_splinefx as x, coil_splinefy as y, coil_splinefz as z.
           CALL cbspline_init(ZC_spl, ncoefs, ier)
           CALL cbspline_setup(ZC_spl, ncoefs, &
                coil_splinesz(icoil,1:nknots), coil_splinefz(icoil,1:ncoefs), ier)

           CALL cbspline_eval(XC_spl, zero, x0, ier);  x00 = x0
           IF (x0 < zero) THEN
              hsgn = -0.5d0
           ELSE
              hsgn = 0.5d0
           END IF
           CALL cbspline_eval(YC_spl, zero, y0, ier);  y00 = y0
           CALL cbspline_eval(ZC_spl, zero, z0, ier);  z00 = z0
           DO j = 2, npts_torx
              s_val = REAL(j-1)/REAL(npts_torx-1)
              CALL cbspline_eval(XC_spl, s_val, x1, ier)
              CALL cbspline_eval(YC_spl, s_val, y1, ier)
              CALL cbspline_eval(ZC_spl, s_val, z1, ier)
              dl(j) = SQRT((x1 - x0)**2 + (y1 - y0)**2 + (z1 - z0)**2)
              va(j) = ATAN2(0.5*(y0 + y1), hsgn*(x0 + x1))
              x0 = x1;  y0 = y1;  z0 = z1
              len = len + dl(j)
              vbar = vbar + dl(j)*va(j)
           END DO !j

           CALL cbspline_delete(ZC_spl)

           va(1) = ATAN2(0.5*(y00 + y0), hsgn*(x00 + x0))
           wt = 1.0d0  ! Theta is not well defined for general 3D coil curves
        END IF !lwindsurf

        !Last segment
        dl(1) = SQRT((x00 - x0)**2 + (y00 - y0)**2 + (z00 - z0)**2)
        len = len + dl(1)
        vbar = (vbar + dl(1)*va(1))/len
        
        vvar = 0.0
        DO j=1, npts_torx
           vvar = vvar + dl(j)*wt(j)*(va(j) - vbar)**2
        END DO !j
        excur = SQRT(vvar/len)

        CALL cbspline_delete(XC_spl);  CALL cbspline_delete(YC_spl)
        DEALLOCATE(dl, wt, va)
      END SUBROUTINE get_coil_tor_excur

!-----------------------------------------------------------------------
!     Subroutine:    get_coil_vbounds
!     Author:        J. Breslau (jbreslau@pppl.gov)
!     Date:          6/4/2019
!     Description:   Computes the min & max v limits of a ws coil within
!                     specified du bounds.
!     Input:         Coil index icoil, u bounds
!     Output:        v_min, v_max
!-----------------------------------------------------------------------
SUBROUTINE get_coil_vbounds(icoil, duu, dul, vmin, vmax)
  USE stel_kinds, ONLY : rprec
  USE stellopt_vars
  USE stellopt_targets, ONLY : npts_crect
  USE stellopt_cbspline
  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(IN)      :: icoil
  REAL(rprec), INTENT(IN)  :: duu, dul
  REAL(rprec), INTENT(OUT) :: vmin, vmax

  ! Local variables
  TYPE(cbspline) :: XC_spl, YC_spl
  REAL(rprec)    :: s_val, u_min, u_val, v_val
  INTEGER        :: nknots, ncoefs, ier, j
  LOGICAL        :: lmod

  vmin = 1.0;  vmax = 0.0

  IF (.NOT.lwindsurf) RETURN
  IF ((coil_type(icoil).NE.'M').AND.(coil_type(icoil).NE.'A')) RETURN

  ! Handle spline boundary conditions
  CALL enforce_spline_bcs(icoil, lmod)
  ncoefs = coil_nctrl(icoil)
  nknots = ncoefs + 4

  ! Set up splines
  CALL cbspline_init(XC_spl, ncoefs, ier)
  CALL cbspline_setup(XC_spl, ncoefs, &
       coil_splinesx(icoil,1:nknots), coil_splinefx(icoil,1:ncoefs), ier)

  CALL cbspline_init(YC_spl, ncoefs, ier)
  CALL cbspline_setup(YC_spl, ncoefs, &
       coil_splinesy(icoil,1:nknots), coil_splinefy(icoil,1:ncoefs), ier)

  ! Upper segment
  DO j=1, npts_crect
     s_val = REAL(j-1,rprec)/REAL(npts_crect-1,rprec)
     CALL cbspline_eval(XC_spl, s_val, u_val, ier)
     IF (u_val > duu) exit
     CALL cbspline_eval(YC_spl, s_val, v_val, ier)
     IF (v_val < vmin) vmin = v_val
     IF (v_val > vmax) vmax = v_val
  END DO !j

  u_min = MAX(1.0d0 - dul, duu)

  ! Lower segment
  DO j=npts_crect, 1, -1
     s_val = REAL(j-1,rprec)/REAL(npts_crect-1,rprec)
     CALL cbspline_eval(XC_spl, s_val, u_val, ier)
     IF (u_val < u_min) exit
     CALL cbspline_eval(YC_spl, s_val, v_val, ier)
     IF (v_val < vmin) vmin = v_val
     IF (v_val > vmax) vmax = v_val
  END DO !j

  ! Clean up
  CALL cbspline_delete(XC_spl);  CALL cbspline_delete(YC_spl)
END SUBROUTINE get_coil_vbounds

!-----------------------------------------------------------------------
SUBROUTINE append_fixed_coils(ounit, inname, maxdex, lscreen)
  USE stel_kinds, ONLY : rprec
  USE safe_open_mod
  IMPLICIT NONE

! Arguments
  INTEGER, INTENT(IN)      :: ounit, maxdex
  CHARACTER(*), INTENT(IN) :: inname
  LOGICAL, INTENT(IN)      :: lscreen

! Local variables
  LOGICAL        :: firstgroup
  INTEGER        :: istat, iunit
  INTEGER        :: nlines, igroup, mingroup, offset
  REAL(rprec)    :: xw, yw, zw, currin
  CHARACTER(256) :: line, group_id

! Open fixed coils file
  iunit = ounit + 1
  CALL safe_open(iunit, istat, TRIM(inname), 'old', 'formatted')
  IF (istat.ne.0) THEN
     PRINT *,'Error opening fixed coil file '//TRIM(inname)
     RETURN
  ENDIF

! Determine the range of coil group indices.
! For now, assume there is no start string.
  firstgroup = .TRUE.
  nlines = 0
  DO
     READ(iunit, '(a)', IOSTAT=istat) line
     IF (istat.NE.0) EXIT !EOF condition
     nlines = nlines + 1
     READ(line,*,iostat=istat) xw, yw, zw, currin, igroup, group_id
     IF (istat.eq.0) THEN ! Group info present
        IF (firstgroup) THEN
           mingroup = igroup
           firstgroup = .FALSE.
        ELSE
           mingroup = MIN(mingroup, igroup)
        ENDIF !!firstgroup
     ENDIF ! Group info present
  ENDDO

  IF (lscreen) PRINT *,'Appending ',nlines,' lines from '//TRIM(inname)//' to coils.'
  offset = maxdex + 1 - mingroup
  REWIND(iunit)
  DO
     READ(iunit, '(a)', IOSTAT=istat) line
     IF (istat.NE.0) EXIT !EOF condition
     READ(line,*,iostat=istat) xw, yw, zw, currin, igroup, group_id
     IF (istat.ne.0) THEN ! Group info not present
        IF (LEN(TRIM(line)).GT.0) WRITE(ounit,'(A)') TRIM(line)
     ELSE !Group info present
        igroup = igroup + offset
        WRITE(ounit,'(4e25.17,i5,A)') xw, yw, zw, currin, igroup, ' '//TRIM(group_id)
     ENDIF
  ENDDO

  CALL safe_close(iunit) !Close the file when done
END SUBROUTINE append_fixed_coils
