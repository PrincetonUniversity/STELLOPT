!-----------------------------------------------------------------------
!     Subroutine:    thrift_load_vmec
!     Authors:       S. Lazerson
!     Date:          11/XX/2022
!     Description:   This subroutine handles the VMEC output for
!                    running stel_tools
!                    A note on VMEC definitions.  VMEC outputs a 1D
!                    array called jcurv.  This should be <j^phi> but
!                    is actually jcurv = (dI/ds)/(2*pi).  
!-----------------------------------------------------------------------
      SUBROUTINE thrift_load_vmec
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE thrift_runtime
      USE thrift_vars
      USE thrift_equil
      USE read_wout_mod
      USE stel_tools
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: nu, nv, mnmax_temp, k, iflag, u, v
      INTEGER, ALLOCATABLE :: xm_temp(:), xn_temp(:)
      REAL(rprec) :: temp, s_temp, u_temp, phi_temp
      REAL(rprec), ALLOCATABLE :: xu(:), xv(:)
      REAL(rprec), ALLOCATABLE :: rmnc_temp(:,:), zmns_temp(:,:), lmns_temp(:,:)
      REAL(rprec), ALLOCATABLE :: rmns_temp(:,:), zmnc_temp(:,:), lmnc_temp(:,:)
      REAL(rprec), ALLOCATABLE :: gmnc_temp(:,:), bmnc_temp(:,:)
      REAL(rprec), ALLOCATABLE :: bsupumnc_temp(:,:), bsupvmnc_temp(:,:)
      REAL(rprec), ALLOCATABLE :: gmns_temp(:,:), bmns_temp(:,:)
      REAL(rprec), ALLOCATABLE :: bsupumns_temp(:,:), bsupvmns_temp(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: mfact(:,:)
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      ! Get the realspace R and Z and metric elements
      nu = 8 * mpol + 1
      nu = 2 ** CEILING(log(REAL(nu))/log(2.0_rprec))
      nv = 4 * ntor + 5                                      ! Use at least 5 toroidal points
      nv = 2 ** CEILING(log(REAL(nv))/log(2.0_rprec)) + 1  ! Odd so we get nfp/2 plane

      ! Place the non-nyquist sized arrays on nyquist sized arrays
      ! Allocate helper arrays for 
      mnmax_temp = mnmax_nyq
      ALLOCATE(xm_temp(mnmax_temp),xn_temp(mnmax_temp))
      ALLOCATE(rmnc_temp(mnmax_temp,ns), zmns_temp(mnmax_temp,ns), lmns_temp(mnmax_temp,ns))
      rmnc_temp=0; zmns_temp=0; lmns_temp = 0
      IF (lasym) THEN
         ALLOCATE(rmns_temp(mnmax_temp,ns), zmnc_temp(mnmax_temp,ns), lmnc_temp(mnmax_temp,ns))
         rmns_temp=0; zmnc_temp=0; lmnc_temp = 0
      END IF
      xm_temp(1:mnmax_temp) = xm_nyq(1:mnmax_temp)
      xn_temp(1:mnmax_temp) = xn_nyq(1:mnmax_temp)
      DO u = 1,mnmax_temp
         DO v = 1, mnmax
            IF ((xm(v) .eq. xm_nyq(u)) .and. (xn(v) .eq. xn_nyq(u))) THEN
               rmnc_temp(u,:) = rmnc(v,:)
               zmns_temp(u,:) = zmns(v,:)
               lmns_temp(u,:) = lmns(v,:)
               IF (lasym) THEN
                  rmns_temp(u,:) = rmns(v,:)
                  zmnc_temp(u,:) = zmnc(v,:)
                  lmnc_temp(u,:) = lmnc(v,:)
               END IF
            END IF
         END DO
      END DO

      ! Allocate Half grid helper quantities
      ALLOCATE(gmnc_temp(mnmax_temp,ns), bmnc_temp(mnmax_temp,ns))
      ALLOCATE(bsupumnc_temp(mnmax_temp,ns), bsupvmnc_temp(mnmax_temp,ns))
      gmnc_temp = 0; bmnc_temp = 0; bsupumnc_temp = 0; bsupvmnc_temp = 0
      IF (lasym) THEN
         ALLOCATE(gmns_temp(mnmax_temp,ns), bmns_temp(mnmax_temp,ns))
         ALLOCATE(bsupumns_temp(mnmax_temp,ns), bsupvmns_temp(mnmax_temp,ns))
         gmns_temp = 0; bmns_temp = 0; bsupumns_temp = 0; bsupvmns_temp = 0
      END IF
      ALLOCATE(mfact(mnmax_temp,2))

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Extrapolate Half grid quantities 
      !    x = [0,1] between bounding gridpoints
      !    whi_even = x
      !    wlo_even = 1-x
      !    whi_odd  = whi_even*sqrt(s/s_hi)
      !    wlo_odd  = wlo_even*sqrt(s/s_lo)
      !    f = f(jlo)*wlo + f(jhi)*whi  Interpolant
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !   Second extrapolate to axis (not reference below)
      k = 1
      WHERE (MOD(NINT(REAL(xm_temp(:))),2) .eq. 0)
         mfact(:,1)= 1.0/2.0
         mfact(:,2)= 3.0/2.0
      ELSEWHERE
         mfact(:,1)= 0
         mfact(:,2)= 0
      ENDWHERE
      lmns_temp(:,k) = mfact(:,1)*lmns_temp(:,k+1)+mfact(:,2)*lmns_temp(:,k+2)
      gmnc_temp(:,k) = mfact(:,1)*gmnc(:,k+1)+mfact(:,2)*gmnc(:,k+2)
      bmnc_temp(:,k) = mfact(:,1)*bmnc(:,k+1)+mfact(:,2)*bmnc(:,k+2)
      bsupumnc_temp(:,k) = mfact(:,1)*bsupumnc(:,k+1)+mfact(:,2)*bsupumnc(:,k+2)
      bsupvmnc_temp(:,k) = mfact(:,1)*bsupvmnc(:,k+1)+mfact(:,2)*bsupvmnc(:,k+2)
      IF (lasym) THEN
         lmnc_temp(:,k) = mfact(:,1)*lmnc_temp(:,k+1)+mfact(:,2)*lmnc_temp(:,k+2)
         gmns_temp(:,k) = mfact(:,1)*gmns(:,k+1)+mfact(:,2)*gmns(:,k+2)
         bmns_temp(:,k) = mfact(:,1)*bmns(:,k+1)+mfact(:,2)*bmns(:,k+2)
         bsupumns_temp(:,k) = mfact(:,1)*bsupumns(:,k+1)+mfact(:,2)*bsupumns(:,k+2)
         bsupvmns_temp(:,k) = mfact(:,1)*bsupvmns(:,k+1)+mfact(:,2)*bsupvmns(:,k+2)
      END IF

      !   Third interpolate from half grid to full (respect overwrite indexing)
      DO k = 2, ns-1
         WHERE (MOD(NINT(REAL(xm_temp(:))),2) .eq. 0)
            mfact(:,1)= 0.5
            mfact(:,2)= 0.5
         ELSEWHERE
            mfact(:,1)= 0.5*SQRT((k-1.0)/(k-1.5)) !rho/rholo
            mfact(:,2)= 0.5*SQRT((k-1.0)/(k-0.5)) !rho/rhohi
         ENDWHERE
         lmns_temp(:,k) = mfact(:,1)*lmns_temp(:,k)+mfact(:,2)*lmns_temp(:,k+1)
         gmnc_temp(:,k) = mfact(:,1)*gmnc(:,k)+mfact(:,2)*gmnc(:,k+1)
         bmnc_temp(:,k) = mfact(:,1)*bmnc(:,k)+mfact(:,2)*bmnc(:,k+1)
         bsupumnc_temp(:,k) = mfact(:,1)*bsupumnc(:,k)+mfact(:,2)*bsupumnc(:,k+1)
         bsupvmnc_temp(:,k) = mfact(:,1)*bsupvmnc(:,k)+mfact(:,2)*bsupvmnc(:,k+1)
         IF (lasym) THEN
            lmnc_temp(:,k) = mfact(:,1)*lmnc_temp(:,k)+mfact(:,2)*lmnc_temp(:,k+1)
            gmns_temp(:,k) = mfact(:,1)*gmns(:,k)+mfact(:,2)*gmns(:,k+1)
            bmns_temp(:,k) = mfact(:,1)*bmns(:,k)+mfact(:,2)*bmns(:,k+1)
            bsupumns_temp(:,k) = mfact(:,1)*bsupumns(:,k)+mfact(:,2)*bsupumns(:,k+1)
            bsupvmns_temp(:,k) = mfact(:,1)*bsupvmns(:,k)+mfact(:,2)*bsupvmns(:,k+1)
         END IF
      END DO

      !   Fourth, extrapolate to ns
      !       note that ns-1 is full grid but ns is on half grid
      k = ns
      WHERE (MOD(NINT(REAL(xm_temp(:))),2) .eq. 0)
         mfact(:,1)= 2.0 ! ns (half grid point)
         mfact(:,2)=-1.0 ! ns-1 (full grid point)
      ELSEWHERE
         mfact(:,1)= 2.0*SQRT((ns-1)/(k-1.5))
         mfact(:,2)=-1.0*SQRT((ns-1)/(k-2.0))
      ENDWHERE
      lmns_temp(:,k) = mfact(:,1)*lmns_temp(:,k)+mfact(:,2)*lmns_temp(:,k-1)
      gmnc_temp(:,k) = mfact(:,1)*gmnc(:,k)+mfact(:,2)*gmnc(:,k-1)
      bmnc_temp(:,k) = mfact(:,1)*bmnc(:,k)+mfact(:,2)*bmnc(:,k-1)
      bsupumnc_temp(:,k) = mfact(:,1)*bsupumnc(:,k)+mfact(:,2)*bsupumnc(:,k-1)
      bsupvmnc_temp(:,k) = mfact(:,1)*bsupvmnc(:,k)+mfact(:,2)*bsupvmnc(:,k-1)
      IF (lasym) THEN
         lmnc_temp(:,k) = mfact(:,1)*lmnc_temp(:,k)+mfact(:,2)*lmnc_temp(:,k-1)
         gmns_temp(:,k) = mfact(:,1)*gmns(:,k)+mfact(:,2)*gmns(:,k-1)
         bmns_temp(:,k) = mfact(:,1)*bmns(:,k)+mfact(:,2)*bmns(:,k-1)
         bsupumns_temp(:,k) = mfact(:,1)*bsupumns(:,k)+mfact(:,2)*bsupumns(:,k-1)
         bsupvmns_temp(:,k) = mfact(:,1)*bsupvmns(:,k)+mfact(:,2)*bsupvmns(:,k-1)
      END IF
      DEALLOCATE(mfact)

      ! Load STEL_TOOLS
      IF (lasym) THEN
         CALL load_fourier_geom(1,ns,mnmax_temp,nu,nv,INT(xm_temp),INT(-xn_temp),iflag,&
                     DBLE(rmnc_temp),DBLE(zmns_temp),RMNS=DBLE(rmns_temp),ZMNC=DBLE(zmnc_temp),&
                     BUMNC=DBLE(bsupumnc_temp),BVMNC=DBLE(bsupvmnc_temp),&
                     BUMNS=DBLE(bsupumns_temp),BVMNS=DBLE(bsupvmns_temp),&
                     LMNS=DBLE(lmns_temp),LMNC=DBLE(lmnc_temp),&
                     BMNC=DBLE(bmnc_temp),BMNS=DBLE(bmns_temp),&
                     GMNC=DBLE(gmnc_temp),GMNS=DBLE(gmns_temp))
         DEALLOCATE(rmns_temp,zmnc_temp,lmnc_temp,&
            bmns_temp,gmns_temp,bsupumns_temp,bsupvmns_temp)
      ELSE
         CALL load_fourier_geom(1,ns,mnmax_temp,nu,nv,INT(xm_temp),INT(-xn_temp),iflag,&
                     DBLE(rmnc_temp),DBLE(zmns_temp),&
                     BUMNC=DBLE(bsupumnc_temp),BVMNC=DBLE(bsupvmnc_temp),&
                     LMNS=DBLE(lmns_temp),&
                     BMNC=DBLE(bmnc_temp),&
                     GMNC=DBLE(gmnc_temp))
      END IF
      DEALLOCATE(xm_temp,xn_temp,rmnc_temp,zmns_temp,lmns_temp,&
         bmnc_temp,gmnc_temp,bsupumnc_temp,bsupvmnc_temp)

      ! Helpers
      b0 = RBtor0/Rmajor
      iota0 = iotaf(1)

      ! We want dV/dPhi = dV/ds*ds/dPhi = dV/ds / (dPhi/ds)
      !     Note that Vp is missing a factor of 4*pi*pi
      vp = pi2*pi2*vp
      vp = vp/phipf

      ! Iota Spline
      bcs1=(/ 0, 0/)
      IF (EZspline_allocated(iota_spl)) CALL EZspline_free(iota_spl,iflag)
      CALL EZspline_init(iota_spl,ns,bcs1,iflag)
      IF (iflag /=0) CALL handle_err(EZSPLINE_ERR,'thrift_load_vmec: iota',iflag)
      iota_spl%isHermite   = 0
      FORALL (k=1:ns) iota_spl%x1(k) = sqrt(DBLE(k-1)/DBLE(ns-1))
      CALL EZspline_setup(iota_spl,iotaf,iflag,EXACT_DIM=.true.)
      IF (iflag /=0) CALL handle_err(EZSPLINE_ERR,'thrift_load_vmec: iota',iflag)

      ! PHI' Spline (toroidal flux derivative)
      bcs1=(/ 0, 0/)
      IF (EZspline_allocated(phip_spl)) CALL EZspline_free(phip_spl,iflag)
      CALL EZspline_init(phip_spl,ns,bcs1,iflag)
      IF (iflag /=0) CALL handle_err(EZSPLINE_ERR,'thrift_load_vmec: phip',iflag)
      phip_spl%isHermite   = 0
      FORALL (k=1:ns) phip_spl%x1(k) = sqrt(DBLE(k-1)/DBLE(ns-1))
      phipf = 2*phipf*phip_spl%x1
      CALL EZspline_setup(phip_spl,phipf,iflag,EXACT_DIM=.true.)
      IF (iflag /=0) CALL handle_err(EZSPLINE_ERR,'thrift_load_vmec: phip',iflag)

      ! dV/dPhi Spline (Volume derivative)
      bcs1=(/ 0, 0/)
      IF (EZspline_allocated(vp_spl)) CALL EZspline_free(vp_spl,iflag)
      CALL EZspline_init(vp_spl,ns,bcs1,iflag)
      IF (iflag /=0) CALL handle_err(EZSPLINE_ERR,'thrift_load_vmec: vp',iflag)
      vp_spl%isHermite   = 0
      FORALL (k=1:ns) vp_spl%x1(k) = sqrt(DBLE(k-1)/DBLE(ns-1))
      CALL EZspline_setup(vp_spl,vp,iflag,EXACT_DIM=.true.)
      IF (iflag /=0) CALL handle_err(EZSPLINE_ERR,'thrift_load_vmec: vp',iflag)


      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_load_vmec

