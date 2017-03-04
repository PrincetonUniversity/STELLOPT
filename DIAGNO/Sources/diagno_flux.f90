!-----------------------------------------------------------------------
!     Module:        diagno_flux
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/02/2012
!     Description:   This subroutine calculates the flux through a
!                    simulated flux loop.
!                    The file defining the flux loop must have the
!                    following format
!  NLOOPS
!  NL     IFLFLG   IDIAFLG   TITLE1
!  X1     Y1     Z1   
!  X2     Y2     Z2           
!   .      .      .         
!   .      .      .         
!   .      .      .
!  XNL-1  YNL-1  ZNL-1
!  XNL    YNL    ZNL
!  NL     IFLFLG   IDIAFLG   TITLE2
!  X1     Y1     Z1   
!  X2     Y2     Z2           
!   .      .      .         
!   .      .      .         
!   .      .      .
!  XNL-1  YNL-1  ZNL-1
!  XNL    YNL    ZNL
!
!  WHERE
!  NLOOPS  : Total Number of flux loops in file
!  NL      : Number of elements in a loop
!  IFLFLG  : Flux Loop Flag (set to 1 if a segment should be repeated)
!  IDIAFLG : Diamagnetic Flag (set to 1 if PHIEDGE should be subtracted)
!  TITLE   : Name of Flux Loop
!  X/Y/Z   : Coordinates of points defining loop (should not close)
!
!-----------------------------------------------------------------------
      SUBROUTINE diagno_flux
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE read_wout_mod, ONLY: nfp
      USE diagno_runtime
      USE safe_open_mod
      USE biotsavart
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
                                  flux_ind(:,:)
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
         REAL FUNCTION dflux_bode(x0,y0,z0,x1,y1,z1,cg)
         USE stel_kinds, ONLY: rprec
         INTEGER, intent(in)     :: cg
         REAL(rprec), intent(in) :: x0,y0,z0,x1,y1,z1
         END FUNCTION dflux_bode
         
         REAL FUNCTION dflux_midpoint(x0,y0,z0,x1,y1,z1,cg)
         USE stel_kinds, ONLY: rprec
         INTEGER, intent(in)     :: cg
         real(rprec), intent(in) :: x0,y0,z0,x1,y1,z1
         END FUNCTION dflux_midpoint
      
         REAL FUNCTION dflux_simpson(x0,y0,z0,x1,y1,z1,cg)
         USE stel_kinds, ONLY: rprec
         INTEGER, intent(in)     :: cg
         REAL(rprec), intent(in) :: x0,y0,z0,x1,y1,z1
         END FUNCTION dflux_simpson
      END INTERFACE
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      if(lverb) write(6,*)' --Calculating Fluxloop Values'
      int_fac=1./int_step
      !Define format control variables
      format_flux_headers='(5X,A14,2X,I3,4X,I3,6X,I1,9X,I1,6X,F5.1)'
      IF (lverb) THEN
         WRITE(6,'(5X,A,2(9X,A,9X))') 'Loop','nseg','flux'
      END IF
      
      ! Read the diagnostics
      iunit = 36
      CALL safe_open(iunit,ier,TRIM(flux_diag_file),'old','formatted')
      READ(iunit,'(1I6)') nfl
      nfl2 = nfl+2
      ALLOCATE(nseg_ar(nfl), title(nfl2), flux(nfl2),&
               iflflg(nfl), idia(nfl), STAT=ier)
      IF (lmut) ALLOCATE(flux_ind(nfl,SIZE(coil_group)), STAT=ier)

      nseg_ar = 0
      DO i = 1, nfl
         READ(iunit,'(3I6,A48)') nseg_ar(i), iflflg(i), idia(i), title(i)
         DO j = 1, nseg_ar(i)
            READ(iunit,*) xp,yp,zp
         END DO
      END DO
      nsegmx=MAXVAL(nseg_ar)+1
      ALLOCATE(xfl(nfl,0:nsegmx),yfl(nfl,0:nsegmx),zfl(nfl,0:nsegmx), STAT=ier)
      REWIND(UNIT=iunit)
      READ(iunit,'(1I6)') nfl
      DO i = 1, nfl
         READ(iunit,'(3I6,A48)') nseg_ar(i), iflflg(i), idia(i), title(i)
         DO j = 1, nseg_ar(i)
            READ(iunit,*) xfl(i,j), yfl(i,j), zfl(i,j)
         END DO
      END DO
      CLOSE(iunit)
      
      ! Calculate the endpoints
      ALLOCATE(dx(nfl,0:nsegmx),dy(nfl,0:nsegmx),dz(nfl,0:nsegmx),&
               xmid(nfl,0:nsegmx),ymid(nfl,0:nsegmx),zmid(nfl,0:nsegmx),STAT=ier)
      DO i = 1, nfl
         nseg = nseg_ar(i)
         IF (iflflg(i) <= 0) THEN
            xfl(i,0) = xfl(i,nseg)
            yfl(i,0) = yfl(i,nseg)
            xfl(i,nseg+1) = xfl(i,1)
            yfl(i,nseg+1) = yfl(i,1)
         ELSE
            xfl(i,0) = xfl(i,nseg)*cos(pi2/nfp)+yfl(i,nseg)*sin(pi2/nfp)
            yfl(i,0) = yfl(i,nseg)*cos(pi2/nfp)-xfl(i,nseg)*sin(pi2/nfp)
            xfl(i,nseg+1) = xfl(i,1)*cos(pi2/nfp) - yfl(i,1)*sin(pi2/nfp)
            yfl(i,nseg+1) = yfl(i,1)*cos(pi2/nfp) + xfl(i,1)*sin(pi2/nfp)
         END IF
         zfl(i,0) = zfl(i,nseg)
         zfl(i,nseg+1) = zfl(i,1)
         dx(i,1:nseg) = xfl(i,2:nseg+1) - xfl(i,1:nseg)
         dy(i,1:nseg) = yfl(i,2:nseg+1) - yfl(i,1:nseg)
         dz(i,1:nseg) = zfl(i,2:nseg+1) - zfl(i,1:nseg)
         IF (iflflg(i) == 0) THEN
            dx(i,0) = dx(i,nseg)
            dy(i,0) = dy(i,nseg)
            dx(i,nseg+1) = dx(i,1)
            dy(i,nseg+1) = dy(i,1)
         ELSE
            dx(i,0) = dx(i,nseg)*cos(pi2/nfp)+dy(i,nseg)*sin(pi2/nfp)
            dy(i,0) = dy(i,nseg)*cos(pi2/nfp)-dx(i,nseg)*sin(pi2/nfp)
            dx(i,nseg+1) = dx(i,1)*cos(pi2/nfp) - dy(i,1)*sin(pi2/nfp)
            dy(i,nseg+1) = dy(i,1)*cos(pi2/nfp) + dx(i,1)*sin(pi2/nfp)
         END IF
         dz(i,0) = dz(i,nseg)
         dz(i,nseg+1) = dz(i,1)
         xmid(i,1:nseg) = 0.5*(xfl(i,1:nseg) + xfl(i,1:nseg))
         ymid(i,1:nseg) = 0.5*(yfl(i,1:nseg) + yfl(i,1:nseg))
         zmid(i,1:nseg) = 0.5*(zfl(i,1:nseg) + zfl(i,1:nseg))
      END DO
      
      ! Read the mutual induction file if present
      iunit = 36
      IF (luse_mut) THEN
         CALL safe_open(iunit,ier,TRIM(flux_mut_file),'old','formatted')
         READ(iunit,'(2I4)') nfl_mut, nextcur
         ALLOCATE(flux_ind(nfl,nextcur), STAT=ier)
         DO i = 1, nfl
            DO ig = 1, nextcur
               READ(iunit,'(1X,I4,1X,I4,E22.14)') i1,i2,flux_ind(i,ig)
            END DO
         END DO
         CLOSE(iunit)
      END IF
      
        
      ! Calculate the fields
      DO i = 1, nfl
         flux(i) = 0.0
         IF (lskip_flux(i)) CYCLE
         nseg = nseg_ar(i)
         IF (lmut) THEN
            SELECT CASE (int_type)
               CASE('midpoint')
                 DO j = 1, nseg
                    DO k = 1, int_step
                       xp = xfl(i,j) + (k-1)*int_fac*dx(i,j)
                       yp = yfl(i,j) + (k-1)*int_fac*dy(i,j)
                       zp = zfl(i,j) + (k-1)*int_fac*dz(i,j)
                       xp1 = xfl(i,j) + k*int_fac*dx(i,j)
                       yp1 = yfl(i,j) + k*int_fac*dy(i,j)
                       zp1 = zfl(i,j) + k*int_fac*dz(i,j)
                       ier = 0
                       IF (lcoil) THEN
                          DO ig = 1, SIZE(coil_group)
                            flux_ind(i,ig)=flux_ind(i,ig) + dflux_midpoint(xp,yp,zp,xp1,yp1,zp1,ig)
                          END DO
                       END IF
                    END DO
                 END DO
               CASE('simpson')
                 DO j = 1, nseg
                    DO k = 1, int_step
                       xp = xfl(i,j) + (k-1)*int_fac*dx(i,j)
                       yp = yfl(i,j) + (k-1)*int_fac*dy(i,j)
                       zp = zfl(i,j) + (k-1)*int_fac*dz(i,j)
                       xp1 = xfl(i,j) + k*int_fac*dx(i,j)
                       yp1 = yfl(i,j) + k*int_fac*dy(i,j)
                       zp1 = zfl(i,j) + k*int_fac*dz(i,j)
                       ier = 0
                       IF (lcoil) THEN
                          DO ig = 1, SIZE(coil_group)
                            flux_ind(i,ig)=flux_ind(i,ig) + dflux_simpson(xp,yp,zp,xp1,yp1,zp1,ig)
                          END DO
                       END IF
                    END DO
                 END DO
               CASE('bode')
                 DO j = 1, nseg
                    DO k = 1, int_step
                       xp = xfl(i,j) + (k-1)*int_fac*dx(i,j)
                       yp = yfl(i,j) + (k-1)*int_fac*dy(i,j)
                       zp = zfl(i,j) + (k-1)*int_fac*dz(i,j)
                       xp1 = xfl(i,j) + k*int_fac*dx(i,j)
                       yp1 = yfl(i,j) + k*int_fac*dy(i,j)
                       zp1 = zfl(i,j) + k*int_fac*dz(i,j)
                       ier = 0
                       IF (lcoil) THEN
                          DO ig = 1, SIZE(coil_group)
                            flux_ind(i,ig)=flux_ind(i,ig) + dflux_bode(xp,yp,zp,xp1,yp1,zp1,ig)
                          END DO
                       END IF
                    END DO
                 END DO
               CASE('adaptive')
              
            END SELECT
            IF (lverb) WRITE(6,'(2X,i6,10X,I5)') i, nseg
         ELSE
            SELECT CASE (int_type)
               CASE('midpoint')
                 DO j = 1, nseg
                    DO k = 1, int_step
                       xp = xfl(i,j) + (k-1)*int_fac*dx(i,j)
                       yp = yfl(i,j) + (k-1)*int_fac*dy(i,j)
                       zp = zfl(i,j) + (k-1)*int_fac*dz(i,j)
                       xp1 = xfl(i,j) + k*int_fac*dx(i,j)
                       yp1 = yfl(i,j) + k*int_fac*dy(i,j)
                       zp1 = zfl(i,j) + k*int_fac*dz(i,j)
                       ier = 0
                       IF (.not.lvac) flux(i)=flux(i) + dflux_midpoint(xp,yp,zp,xp1,yp1,zp1,ier)
                       IF (lcoil) THEN
                          DO ig = 1, SIZE(coil_group)
                            IF (.not.luse_extcur(ig)) CYCLE
                            flux(i)=flux(i) + dflux_midpoint(xp,yp,zp,xp1,yp1,zp1,ig)
                          END DO
                       END IF
                    END DO
                 END DO
               CASE('simpson')
                 DO j = 1, nseg
                    DO k = 1, int_step
                       xp = xfl(i,j) + (k-1)*int_fac*dx(i,j)
                       yp = yfl(i,j) + (k-1)*int_fac*dy(i,j)
                       zp = zfl(i,j) + (k-1)*int_fac*dz(i,j)
                       xp1 = xfl(i,j) + k*int_fac*dx(i,j)
                       yp1 = yfl(i,j) + k*int_fac*dy(i,j)
                       zp1 = zfl(i,j) + k*int_fac*dz(i,j)
                       ier = 0
                       IF (.not.lvac) flux(i)=flux(i) + dflux_simpson(xp,yp,zp,xp1,yp1,zp1,ier)
                       IF (lcoil) THEN
                          DO ig = 1, SIZE(coil_group)
                            IF (.not.luse_extcur(ig)) CYCLE
                            flux(i)=flux(i) + dflux_simpson(xp,yp,zp,xp1,yp1,zp1,ig)
                          END DO
                       END IF
                    END DO
                 END DO
               CASE('bode')
                 DO j = 1, nseg
                    DO k = 1, int_step
                       xp = xfl(i,j) + (k-1)*int_fac*dx(i,j)
                       yp = yfl(i,j) + (k-1)*int_fac*dy(i,j)
                       zp = zfl(i,j) + (k-1)*int_fac*dz(i,j)
                       xp1 = xfl(i,j) + k*int_fac*dx(i,j)
                       yp1 = yfl(i,j) + k*int_fac*dy(i,j)
                       zp1 = zfl(i,j) + k*int_fac*dz(i,j)
                       ier = 0
                       IF (.not.lvac) flux(i)=flux(i) + dflux_bode(xp,yp,zp,xp1,yp1,zp1,ier)
                       IF (lcoil) THEN
                          DO ig = 1, SIZE(coil_group)
                            IF (.not.luse_extcur(ig)) CYCLE
                            flux(i)=flux(i) + dflux_bode(xp,yp,zp,xp1,yp1,zp1,ig)
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
            IF (luse_mut) THEN
               DO ig = 1, nextcur
                  flux(i)=flux(i) + flux_ind(i,ig)*extcur(ig)
               END DO
            END IF
            IF (lverb) WRITE(6,'(2X,i6,10X,I5,10X,E18.8)') i, nseg, flux(i)
         END IF
      END DO
      
      iunit = 36
      IF (lmut) THEN
         CALL safe_open(iunit,ier,TRIM(flux_mut_file),'replace','formatted')
         WRITE(iunit,'(1X,I4,1X,I4)') i, SIZE(coil_group)
         DO i = 1, nfl
            DO ig = 1, SIZE(coil_group)
               WRITE(iunit,'(1X,I4,1X,I4,ES22.12E3)') i,ig,flux_ind(i,ig)
            END DO
         END DO
         DO i = 1, nfl
            WRITE(iunit,'(I4,A)') i, title(i)
         END DO
         CLOSE(iunit)
      ELSE
         DO i = 1, nfl
            IF (iflflg(i) == 1) flux(i) = flux(i) * nfp_diagno
            IF (idia(i) == 1) flux(i) = flux(i) + phiedge * eq_sgns
            IF (idia(i) < 0) flux(i) = flux(i) - flux(ABS(idia(i))) ! Do this so we can do flux subtraction
            flux(i) = flux(i) * flux_turns(i)
            IF (lskip_flux(i)) flux(i) = 0
         END DO
         !IF (lcomp_dia) THEN
         !   flux(nfl+1) = SUM(flux(1:nfl)*dia_comp_coefs(1:nfl))
         !   title(nfl+1) = 'Comp. Diamag. Loop.'
         !   flux_comb_coefs(nfl+1) = 0
         !END IF
         !IF (lflux_comb) THEN
         !   flux(nfl+2) = SUM(flux(1:nfl+1)*flux_comb_coefs(1:nfl+1))
         !   title(nfl+2) = 'Combined Flux Loops'
         !END IF
         
         ! Output the response
         CALL safe_open(iunit,ier,'diagno_flux.'//TRIM(id_string),'replace','formatted')
         WRITE(iunit,'(1X,I4.4)') nfl
         DO i = 1, nfl
            WRITE(iunit,'(1X,1ES22.12E3)') flux(i)
         END DO
         DO i = 1, nfl
            WRITE(iunit,'(A)') title(i)
         END DO
         CALL FLUSH(iunit)
         CLOSE(iunit)
      END IF
      
      ! Clean up
      IF (ALLOCATED(flux_ind)) DEALLOCATE(flux_ind)
      DEALLOCATE(nseg_ar, title, flux, iflflg, idia)
      DEALLOCATE(xfl,yfl,zfl,dx,dy,dz,xmid,ymid,zmid)
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE diagno_flux
