!-----------------------------------------------------------------------
!     Subroutine:    spec_init_wout
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          04/9/2013
!     Description:   This subroutine initizlizes the SPEC interfaces
!                    from the VMEC interfaces by iterpolating over the
!                    Fourier modes.
!-----------------------------------------------------------------------
      SUBROUTINE spec_init_wout
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE read_wout_mod, nfp_in => nfp, lasym_in => lasym, &
                      lfreeb_in => lfreeb, mnmax_in => mnmax, &
                      rmnc_in => rmnc, zmns_in => zmns, &
                      rmns_in => rmns, zmnc_in => zmnc, &
                      bsupumnc_in => bsupumnc, bsupumns_in => bsupumns,&
                      bsupvmnc_in => bsupvmnc, bsupvmns_in => bsupvmns,&
                      xn_in => xn, xm_in => xm, mgrid_file_in => mgrid_file,&
                      nextcur_in => nextcur, extcur_in => extcur, &
                      iotaf_in => iotaf
      USE spec_background
      USE spec_runtime
      USE spec_profile, ONLY: press, iprime, p_spl, ip_spl,torflux_edge,&
                              curtor, p_cubspl, ip_cubspl, iotaf, iota_spl,&
                              adiab, pscale
      USE EZspline_obj
      USE EZspline
!-----------------------------------------------------------------------
!     Local Variables
!          iunit       File Unit Number
!          ierr        Error flag
!          mn          Dummy index over modes
!          im          Dummy index over modes
!          in          Dummy index over modes
!          uv          Dummy index over real-space
!          ik          Dummy index over radial surfaces
!          i           Dummy index
!          j           Dummy index
!          bcs1        Boundary Condition Array for EZspline
!          val1        Dummy value used in boundary extension
!          val2        Dummy value used in boundary extension
!          dval        Dummy value used in boundary extension
!          rhomax      Twice the maximum VMEC rho in each v section
!          rreal_in    VMEC R in (s,u,v) space
!          zreal_in    VMEC Z in (s,u,v) space
!          bureal_in   VMEC Bsupu in (s,u,v) space
!          bvreal_in   VMEC Bsupv in (s,u,v) space
!          rho_vmec    VMEC radial grid
!          f_spl       Dummy spline object for spline to PIES grid
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: iunit, ier, mn, im, in, ik , i, j, n_int
      INTEGER :: bcs1(2)
      REAL(rprec) :: val1, val2, dval, rhomax,pi2, scale, theta,r1, z1,&
                     r2, z2, theta2, dtheta, iota_min, iota_max
      REAL(rprec) :: p1,q1,p2,q2, flux1,flux2,dflux
      REAL(rprec), ALLOCATABLE :: rho_vmec(:), ftemp(:), ptemp(:),&
                                  jtemp(:), btemp(:), itemp(:)
      REAL(rprec), ALLOCATABLE :: rmnc_temp(:,:), zmns_temp(:,:)
      REAL(rprec), ALLOCATABLE :: rmns_temp(:,:), zmnc_temp(:,:)
      TYPE(EZspline1_r8) :: f_spl, jdotb_spl, b_spl
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      bcs1=(/0,0/)
      pi2 = 8 * ATAN(1._rprec)
      CALL read_wout_file(TRIM(id_string(7:)),ierr=ier)
      IF (ier /= 0) CALL handle_err(VMEC_WOUT_ERR,id_string,ier)
      ! Now set some background values
      lasym = lasym_in
      lfreeb = lfreeb_in
      ! Initialize flux
      IF (ALL(TFLUX == 0)) CALL spec_rat_prof
      ! Do nothing for now
      IF (lfreeb) THEN
      !   mgrid_file = TRIM(mgrid_file_in)
      !   nv_vmec = 32 ! Will need to read the INDATA namelist to get this
      !   nextcur = nextcur_in
      !   ALLOCATE(extcur(1:nextcur),STAT=ier)
      !   IF (ier /= 0) CALL handle_err(ALLOC_ERR,'EXTCUR',ier)
      !   extcur(1:nextcur) = extcur_in(1:nextcur_in)
      ENDIF
      IF (isigng >0) THEN
         signgs = -1
      ELSE
         signgs = 1
      END IF
      nfp = nfp_in
      m = mpol-1
      n = ntor
      mnmax = (m+1)*(2*n+1)
      iotamn = MINVAL(iotaf_in(:))
      iotamx = MAXVAL(iotaf_in(:))
      torflux_edge = phi(ns)
      curtor = Itor
      id_string = TRIM(input_extension)
      if (lverb) THEN
         write(6,*) '-----VMEC File Parameters-----'
         write(6,'(A,A)') '    file: ',TRIM(id_string)
         write(6,'(A,I3,A,I3)') '       m: ',m
         write(6,'(A,I3,A,I3)') '       n: ',n
         write(6,'(A,I4,A,I4,A,I5)') '   mnmax: ',mnmax
         write(6,'(A,I3)') '     nfp: ',nfp
         write(6,'(A,I3)') '      ns: ',ns
         write(6,'(A,L3)') '  lfreeb: ',lfreeb
         write(6,'(A,F6.3,A,F6.3,A)')    '    iota: [',iotamn,',',iotamx,']'  
         write(6,'(A,F7.3)')             'torflux_edge: ',torflux_edge 
         IF (ABS(curtor) .ge. 1.0E6) THEN
            write(6,'(A,F8.3,A)')             'Total Current: ',curtor/1.0E6,' [MA]'
         ELSE IF (ABS(curtor) .ge. 1.0E3) THEN
            write(6,'(A,F8.3,A)')             'Total Current: ',curtor/1.0E3,' [kA]'
         ELSE
            write(6,'(A,F8.3,A)')             'Total Current: ',curtor,' [A]'
         END IF
      END IF
      ! Now allocate 1D variables
      ALLOCATE(xm(1:mnmax),xn(1:mnmax),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'XM XN',ier)
      ALLOCATE(rho(0:nvol),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'RHO',ier)
      ALLOCATE(rho_vmec(1:ns),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'RHO VMEC',ier)
      ! Now allocate 2D variables
      ALLOCATE(rmnc(1:mnmax,0:nvol),zmns(1:mnmax,0:nvol),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'R Z',ier)
      IF (lasym) THEN
         ALLOCATE(rmns(1:mnmax,0:nvol),zmnc(1:mnmax,0:nvol),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'Rs Zc',ier)
      END IF
      ! Now initizlie modes
      mn = 1
      DO i=0,m
         DO j=-n,n
           xm(mn) = i
           xn(mn) = j
           mn = mn + 1
         END DO
      END DO
      ! Initialize rho
      rho(0) = 0.0
      ! Now set the interpolation axis to normalized toroidal flux
      rho_vmec = phi / phi(ns)
      rho(0) = 0.0
      rho(1:nvol) = tflux(1:nvol)
      ! Now interpolate to SPEC Surfaces
      ALLOCATE(ftemp(1:ns),STAT=ier)
      ALLOCATE(rmnc_temp(1:mnmax_in,0:nvol),zmns_temp(1:mnmax_in,0:nvol),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'R Z',ier)
      IF (lasym) THEN
         ALLOCATE(rmns_temp(1:mnmax_in,0:nvol),zmnc_temp(1:mnmax_in,0:nvol),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'R Z',ier)
      END IF
      DO mn = 1, mnmax_in
         !RMNC
         ! We want the fourier coefficients to go to zero linearly in toroidal flux
         ! so we multiply by flux**(-m)
         ! we then interpolate in flux space
         CALL EZspline_init(f_spl,ns,bcs1,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/pies_init_wout',ier)
         f_spl%isHermite = 1
         f_spl%x1 = rho_vmec
         ftemp(:) = rmnc_in(mn,:)*rho_vmec(:)**(-xm(mn))
         IF (xm(mn) /= 0.0) ftemp(1) = 0.0
         CALL EZspline_setup(f_spl,ftemp,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/pies_init_wout',ier)
         DO ik = 0, nvol
            CALL EZspline_interp(f_spl,rho(ik),rmnc_temp(mn,ik),ier)
         END DO
         rmnc_temp(mn,:) = rmnc_temp(mn,:)*rho(:)**(xm(mn))
         CALL EZspline_free(f_spl,ier)
         ! ZMNS
         CALL EZspline_init(f_spl,ns,bcs1,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/pies_init_wout',ier)
         f_spl%isHermite = 1
         f_spl%x1 = rho_vmec
         ftemp(:) = zmns_in(mn,:)*rho_vmec(:)**(-xm(mn))
         IF (xm(mn) /= 0.0) ftemp(1) = 0.0
         CALL EZspline_setup(f_spl,ftemp,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/pies_init_wout',ier)
         DO ik = 0, nvol
            CALL EZspline_interp(f_spl,rho(ik),zmns_temp(mn,ik),ier)
         END DO
         zmns_temp(mn,:) = zmns_temp(mn,:)*rho(:)**(xm(mn))
         CALL EZspline_free(f_spl,ier)
         IF (lasym) THEN
            !RMNS
            CALL EZspline_init(f_spl,ns,bcs1,ier)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/pies_init_wout',ier)
            f_spl%isHermite = 1
            f_spl%x1 = rho_vmec
            ftemp(:) = rmns_in(mn,:)*rho_vmec(:)**(-xm(mn))
            IF (xm(mn) /= 0.0) ftemp(1) = 0.0
            CALL EZspline_setup(f_spl,ftemp,ier)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/pies_init_wout',ier)
            DO ik = 0, nvol
               CALL EZspline_interp(f_spl,rho(ik),rmns_temp(mn,ik),ier)
            END DO
            rmns_temp(mn,:) = rmns_temp(mn,:)*rho(:)**(xm(mn))
            CALL EZspline_free(f_spl,ier)
            ! ZMNC
            CALL EZspline_init(f_spl,ns,bcs1,ier)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/pies_init_wout',ier)
            f_spl%isHermite = 1
            f_spl%x1 = rho_vmec
            ftemp(:) = zmnc_in(mn,:)*rho_vmec(:)**(-xm(mn))
            IF (xm(mn) /= 0.0) ftemp(1) = 0.0
            CALL EZspline_setup(f_spl,ftemp,ier)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/pies_init_wout',ier)
            DO ik = 0, nvol
               CALL EZspline_interp(f_spl,rho(ik),zmnc_temp(mn,ik),ier)
            END DO
            zmnc_temp(mn,:) = zmnc_temp(mn,:)*rho(:)**(xm(mn))
            CALL EZspline_free(f_spl,ier)
         END IF
      END DO
      DEALLOCATE(ftemp)
      ! Now copy from rmnc_temp to rmnc  (mu-nv) -> (mn-nv)
      rmnc = 0.0
      zmns = 0.0
      xn_in = xn_in/nfp
      DO mn = 1, mnmax
         DO j = 1, mnmax_in
            IF ((xm(mn) .eq. xm_in(j)) .and. (xn(mn) .eq. xn_in(j))) THEN 
               rmnc(mn,0:nvol) =  rmnc_temp(j,0:nvol)
               zmns(mn,0:nvol) =  zmns_temp(j,0:nvol)
               IF (lasym) THEN
                  rmns(mn,0:nvol) =  rmns_temp(j,0:nvol)
                  zmnc(mn,0:nvol) =  zmnc_temp(j,0:nvol)
               END IF
            END IF
         END DO
      END DO
      DEALLOCATE(rmnc_temp,zmns_temp)
      IF (lasym) DEALLOCATE(rmns_temp,zmnc_temp)
      
      ! Adjust modes so clockwise theta coordinate
      DO mn = 1, mnmax
         IF ((xm(mn) == 1) .and. (xn(mn) == 0) .and. (zmns(mn,nvol) > 0)) THEN
            xn = -xn
            zmns = -zmns
            IF (lasym) rmns = -rmns
            lflipped = .true.
            EXIT
         END IF
      END DO
      
      ! Now initialize pressure profiles
      n_int = 100
      pflux = 0.0
      iota(0) = iotaf_in(1)
      ALLOCATE(press(1:nvol),adiab(1:nvol),STAT=ier)
      ALLOCATE(ftemp(1:n_int),ptemp(1:n_int),jtemp(1:n_int),btemp(1:n_int),itemp(1:n_int),STAT=ier)
      CALL EZspline_init(p_spl,ns,bcs1,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/spec_init_wout(p)',ier)
      CALL EZspline_init(jdotb_spl,ns,bcs1,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/spec_init_wout(jdotb)',ier)
      CALL EZspline_init(b_spl,ns,bcs1,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/spec_init_wout(b)',ier)
      CALL EZspline_init(iota_spl,ns,bcs1,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/spec_init_wout(iota)',ier)
      p_spl%isHermite = 1      
      p_spl%x1 = phi/phi(ns)       !Spline over normalized toroidal flux
      jdotb_spl%isHermite = 1      
      jdotb_spl%x1 = phi/phi(ns)       !Spline over normalized toroidal flux
      b_spl%isHermite = 1      
      b_spl%x1 = phi/phi(ns)       !Spline over normalized toroidal flux
      iota_spl%isHermite = 1      
      iota_spl%x1 = phi/phi(ns)       !Spline over normalized toroidal flux
      CALL EZspline_setup(p_spl,presf,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/spec_init_wout',ier)
      CALL EZspline_setup(jdotb_spl,jdotb,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/spec_init_wout',ier)
      CALL EZspline_setup(b_spl,bdotgradv,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/spec_init_wout',ier)
      CALL EZspline_setup(iota_spl,iotaf_in,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/spec_init_wout',ier)
      !PRINT *,'tflux ',tflux
      flux1 = 0.0
      pflux = 0.0
      DO ik = 1, nvol
         flux2 = tflux(ik)
         dflux = flux2-flux1
         DO i = 1, n_int
            ftemp(i) = dflux*(i-1)/(n_int-1)+flux1
            !ftemp(i) = tflux(ik)*i/n_int+tflux(ik-1)
         END DO
         CALL EZspline_interp(p_spl,n_int,ftemp,ptemp,ier)
         CALL EZspline_interp(jdotb_spl,n_int,ftemp,jtemp,ier)
         CALL EZspline_interp(b_spl,n_int,ftemp,btemp,ier)
         CALL EZspline_interp(iota_spl,n_int,ftemp,itemp,ier)
         press(ik) = mu0*SUM(ptemp*dflux/n_int)/dflux
         mu(ik) = mu0*SUM(jtemp/btemp/btemp)/n_int
         adiab(ik) = 5./3.
         ! Chi' = iota*Phi'
         pflux(ik) = SUM(itemp*dflux)/torflux_edge
         flux1 = flux2
      END DO
      DO ik = 2, nvol
         pflux(ik) = pflux(ik) + pflux(ik-1)
         PRINT *,tflux(ik),press(ik),mu(ik),pflux(ik)
      END DO
      CALL EZspline_free(p_spl,ier)
      ! Get Pressure and mu in magnetic units
      !press = press*mu0
      !mu    = abs(mu0*mu/phi(ns))
      !mu    = mu*mu0
      IF (ANY(press > 0.0)) pscale = 1.0
      
      ! Now calculate Current/iota Profile
      !IF (ncurr == 0) THEN
         !pflux(0) = 0.0
         !iota(0) = iotaf_in(1)
         !CALL EZspline_init(iota_spl,ns,bcs1,ier)
         !IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/spec_init_wout',ier)
         !iota_spl%isHermite = 1      
         !iota_spl%x1 = phi/phi(ns)       !Spline over normalized toroidal flux
         !CALL EZspline_setup(iota_spl,iotaf_in,ier)
         !IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/spec_init_wout',ier)
         !DO ik = 1, nvol
         !   CALL EZspline_interp(iota_spl,tflux(ik),iota(ik),ier)
         !   pflux(ik) = tflux(ik) * iota(ik)
         !END DO
         !pr(0:nvol) = 0.0
         !qr(0:nvol) = 1.0
         !pl(0:nvol) = 0.0
         !ql(0:nvol) = 1.0
      !ELSE IF (ncurr == 1) THEN
      !   curtor=Itor
      !   stop 'ncurr==1 not supported yet'
      !ELSE
      !   ! Some ncurr .ne. 1 or 0 Error
      !   stop 'ncurr .ne 1 or 0'
      !END IF
      
         
      
      
      ! Free VMEC wout memory
      CALL read_wout_deallocate
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE spec_init_wout
