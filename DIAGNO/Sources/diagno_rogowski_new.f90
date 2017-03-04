!-----------------------------------------------------------------------
!     Module:        diagno_rogowski_new
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/02/2012
!     Description:   This subroutine calculates the response of a
!                    simulated Rogowski loop.
!                    The file defining the Rogowski coil must have the
!                    following format
!  NLOOPS
!  NL     IFLFLG   IDIAFLG   TITLE1
!  X1     Y1     Z1    EFF_AREA1   
!  X2     Y2     Z2    EFF_AREA2           
!   .      .      .         .
!   .      .      .         .
!   .      .      .         .
!  XNL-1  YNL-1  ZNL-1 EFF_AREA NL-1
!  XNL    YNL    ZNL   EFF_AREA NL
!  NL     IFLFLG   IDIAFLG   TITLE2
!  X1     Y1     Z1    EFF_AREA1   
!  X2     Y2     Z2    EFF_AREA2           
!   .      .      .         .
!   .      .      .         .
!   .      .      .         .
!  XNL-1  YNL-1  ZNL-1 EFF_AREA NL-1
!  XNL    YNL    ZNL   EFF_AREA NL
!
!  WHERE
!  NLOOPS  : Total Number of Rogowski segements in file
!  NL      : Number of elements in a segment
!  IFLFLG  : Flux Loop Flag (set to 1 if a segment should be repeated)
!  IDIAFLG : Diamagnetic Flag (set to 1 if PHIEDGE should be subtracted)
!  TITLE   : Name of Flux Loop
!  X/Y/Z   : Coordinates of points defining loop (should not close)
!  EFF_AREA: Coefficient to multiply B.dl by (set to zero to start next
!                   segment)
!
!-----------------------------------------------------------------------
      SUBROUTINE diagno_rogowski_new
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE diagno_runtime
      USE biotsavart
      USE safe_open_mod
      USE EZspline_obj
      USE EZspline
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, i, j, k, iunit, nfl, nfl2, nsegmx, nseg, ig,&
                 i1, i2, nfl_mut
      INTEGER, ALLOCATABLE :: nseg_ar(:), iflflg(:), idia(:)
      REAL(rprec) :: br, bphi, bz, temp, xp, yp, zp, xp1, yp1, zp1,int_fac
      REAL(rprec), ALLOCATABLE :: flux(:)
      REAL(rprec), ALLOCATABLE :: xfl(:,:), yfl(:,:), zfl(:,:), &
                                  dx(:,:), dy(:,:), dz(:,:),&
                                  xmid(:,:), ymid(:,:), zmid(:,:),&
                                  flux_ind(:,:), eff_area(:,:)
      CHARACTER(256) :: format_flux_headers
      CHARACTER(48), DIMENSION(:), ALLOCATABLE :: title
      TYPE(EZspline1_r8) :: x_spl,y_spl,z_spl
      INTEGER     ::  bcs0(2) = (/ 0, 0/)
      
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
      if(lverb) write(6,*)' --Calculating Rogowski Coil Values'
      int_fac=1./int_step
      !Define format control variables
      format_flux_headers='(5X,A14,2X,I3,4X,I3,6X,I1,9X,I1,6X,F5.1)'
      IF (lverb) THEN
         WRITE(6,'(5X,A,2(9X,A,9X))') 'Loop','nseg','Signal'
      END IF
      
      ! Read the diagnostics
      iunit = 36
      CALL safe_open(iunit,ier,TRIM(seg_rog_file),'old','formatted')
      READ(iunit,'(1I6)') nfl
      ALLOCATE(nseg_ar(nfl), title(nfl), flux(nfl),&
               iflflg(nfl), idia(nfl), STAT=ier)

      nseg_ar = 0
      DO i = 1, nfl
         READ(iunit,'(3I6,A48)') nseg_ar(i), iflflg(i), idia(i), title(i)
         DO j = 1, nseg_ar(i)
            READ(iunit,*) xp,yp,zp,xp1
         END DO
      END DO
      nsegmx=MAXVAL(nseg_ar)
      ALLOCATE(xfl(nfl,1:nsegmx),yfl(nfl,1:nsegmx),zfl(nfl,1:nsegmx),&
               eff_area(nfl,1:nsegmx),STAT=ier)
      REWIND(UNIT=iunit)
      READ(iunit,'(1I6)') nfl
      DO i = 1, nfl
         READ(iunit,'(3I6,A48)') nseg_ar(i), iflflg(i), idia(i), title(i)
         DO j = 1, nseg_ar(i)
            READ(iunit,*) xfl(i,j), yfl(i,j), zfl(i,j), eff_area(i,j)
         END DO
      END DO
      CLOSE(iunit)
      
      ! Calculate the endpoints
      ALLOCATE(dx(nfl,1:nsegmx),dy(nfl,1:nsegmx),dz(nfl,1:nsegmx),&
               xmid(nfl,1:nsegmx),ymid(nfl,1:nsegmx),zmid(nfl,1:nsegmx),STAT=ier)
      DO i = 1, nfl
         nseg = nseg_ar(i)
         dx(i,1:nseg-1) = xfl(i,2:nseg)-xfl(i,1:nseg)
         dy(i,1:nseg-1) = yfl(i,2:nseg)-yfl(i,1:nseg)
         dz(i,1:nseg-1) = zfl(i,2:nseg)-zfl(i,1:nseg)
         xmid(i,1:nseg-1) = 0.5 * ( xfl(i,2:nseg) + xfl(i,2:nseg-1) )
         ymid(i,1:nseg-1) = 0.5 * ( yfl(i,2:nseg) + yfl(i,2:nseg-1) )
         zmid(i,1:nseg-1) = 0.5 * ( zfl(i,2:nseg) + zfl(i,2:nseg-1) )
      END DO
        
      ! Calculate the fields
      DO i = 1, nfl
         flux(i) = 0.0
         IF (lskip_rogo(i)) CYCLE
         nseg = nseg_ar(i)
         IF (lmut) THEN
         	 
         ELSE
            SELECT CASE (int_type)
               CASE('midpoint')
                 DO j = 1, nseg-1
                    DO k = 1, int_step
                       xp = xfl(i,j) + (k-1)*int_fac*dx(i,j)
                       yp = yfl(i,j) + (k-1)*int_fac*dy(i,j)
                       zp = zfl(i,j) + (k-1)*int_fac*dz(i,j)
                       xp1 = xfl(i,j) + k*int_fac*dx(i,j)
                       yp1 = yfl(i,j) + k*int_fac*dy(i,j)
                       zp1 = zfl(i,j) + k*int_fac*dz(i,j)
                       ier = 0
                       IF (.not.lvac) flux(i)=flux(i) + db_midpoint(xp,yp,zp,xp1,yp1,zp1,ier)*eff_area(i,j)
                       IF (lcoil) THEN
                          DO ig = 1, SIZE(coil_group)
                            IF (.not.luse_extcur(ig)) CYCLE
                            flux(i)=flux(i) + db_midpoint(xp,yp,zp,xp1,yp1,zp1,ig)*eff_area(i,j)
                          END DO
                       END IF
                    END DO
                 END DO
               CASE('simpson')
                 DO j = 1, nseg-1
                    DO k = 1, int_step
                       xp = xfl(i,j) + (k-1)*int_fac*dx(i,j)
                       yp = yfl(i,j) + (k-1)*int_fac*dy(i,j)
                       zp = zfl(i,j) + (k-1)*int_fac*dz(i,j)
                       xp1 = xfl(i,j) + k*int_fac*dx(i,j)
                       yp1 = yfl(i,j) + k*int_fac*dy(i,j)
                       zp1 = zfl(i,j) + k*int_fac*dz(i,j)
                       ier = 0
                       IF (.not.lvac) flux(i)=flux(i) + db_simpson(xp,yp,zp,xp1,yp1,zp1,ier)*eff_area(i,j)
                       IF (lcoil) THEN
                          DO ig = 1, SIZE(coil_group)
                            IF (.not.luse_extcur(ig)) CYCLE
                            flux(i)=flux(i) + db_simpson(xp,yp,zp,xp1,yp1,zp1,ig)*eff_area(i,j)
                          END DO
                       END IF
                    END DO
                 END DO
               CASE('bode')
                 DO j = 1, nseg-1
                    DO k = 1, int_step
                       xp = xfl(i,j) + (k-1)*int_fac*dx(i,j)
                       yp = yfl(i,j) + (k-1)*int_fac*dy(i,j)
                       zp = zfl(i,j) + (k-1)*int_fac*dz(i,j)
                       xp1 = xfl(i,j) + k*int_fac*dx(i,j)
                       yp1 = yfl(i,j) + k*int_fac*dy(i,j)
                       zp1 = zfl(i,j) + k*int_fac*dz(i,j)
                       ier = 0
                       IF (.not.lvac) flux(i)=flux(i) + db_bode(xp,yp,zp,xp1,yp1,zp1,ier)*eff_area(i,j)
                       IF (lcoil) THEN
                          DO ig = 1, SIZE(coil_group)
                            IF (.not.luse_extcur(ig)) CYCLE
                            flux(i)=flux(i) + db_bode(xp,yp,zp,xp1,yp1,zp1,ig)*eff_area(i,j)
                          END DO
                       END IF
                    END DO
                 END DO
               CASE('simpson_spline')
                  !IF (EZspline_allocated(x_spl)) CALL EZspline_free(x_spl,ier)
                  !IF (EZspline_allocated(y_spl)) CALL EZspline_free(y_spl,ier)
                  !IF (EZspline_allocated(z_spl)) CALL EZspline_free(z_spl,ier)
                  !CALL EZspline_init(x_spl,nseg+2,bcs0,ier)
                  !CALL EZspline_init(y_spl,nseg+2,bcs0,ier)
                  !CALL EZspline_init(z_spl,nseg+2,bcs0,ier)
                  !x_spl%isHermite = 1
                  !y_spl%isHermite = 1
                  !z_spl%isHermite = 1
                  !CALL EZspline_setup(x_spl,xfl(0:nseg+1),ier)
                  !CALL EZspline_setup(y_spl,yp(0:nseg+1),ier)
                  !CALL EZspline_setup(z_spl,zp(0:nseg+1),ier)
            END SELECT
            IF (lverb) WRITE(6,'(2X,i6,10X,I5,10X,E18.8)') i, nseg, flux(i)
         END IF
      END DO
      
      ! Output diagnostic file
      iunit = 36
      CALL safe_open(iunit,ier,'diagno_seg.'//TRIM(id_string),'replace','formatted')
      WRITE(iunit,'(1X,I4.4)') nfl
      DO i = 1, nfl
         WRITE(iunit,'(1X,1ES22.12E3)') flux(i)
      END DO
      DO i = 1, nfl
         WRITE(iunit,'(A)') title(i)
      END DO
      CALL FLUSH(iunit)
      CLOSE(iunit)
      
      ! Clean up
      IF (ALLOCATED(flux_ind)) DEALLOCATE(flux_ind)
      DEALLOCATE(nseg_ar, title, flux, iflflg, idia)
      DEALLOCATE(xfl,yfl,zfl,dx,dy,dz,xmid,ymid,zmid,eff_area)
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE diagno_rogowski_new
