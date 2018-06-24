!-----------------------------------------------------------------------
!     Module:        virtual_casing_mod
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/27/2012
!     Description:   This subroutine calculates the plasma response in
!                    real space coordinates through a virtual casing
!                    principle.  The code utilizes the ezspline package
!                    which is part of the NTCC.  Note if adapt_tol is
!                    set less than zero the adaptive call will call the
!                    non-adaptive routines.
!     Usage:
!                    To load an equilibrium using virtual casing from VMEC
!
!                    ALLOCATE(xm_temp(mnmax),xn_temp(mnmax)) ! Integer Arrays
!                    ALLOCATE(rmnc_temp(mnmax,2),zmns_temp(mnmax,2)) ! DOUBLE PRECISION
!                    ALLOCATE(bumnc_temp(mnmax,1),bvmnc_temp(mnmax,1)) ! DOUBLE PRECISION
!                    xm_temp=xm
!                    xn_temp=-xn/nfp
!                    rmnc_temp(:,1)=rmnc(:,ns-1)
!                    rmnc_temp(:,2)=rmnc(:,ns)
!                    zmns_temp(:,1)=zmns(:,ns-1)
!                    zmns_temp(:,2)=zmns(:,ns)
!                    bumnc_temp(:,1) = (1.5*bsupumnc(:,ns) - 0.5*bsupumnc(:,ns-1))
!                    bvmnc_temp(:,1) = (1.5*bsupvmnc(:,ns) - 0.5*bsupvmnc(:,ns-1))                    
!                    CALL init_virtual_casing(mnmax,nu2,nv2,xm_temp,xn_temp,&
!                                         rmnc_temp,zmns_temp,nfp,&
!                                         BUMNC=bumnc_temp,BVMNC=bvmnc_temp)
!                               -or-                    
!                    CALL init_virtual_casing(mnmax,nu2,nv2,xm_temp,xn_temp,&
!                                         rmnc_temp,zmns_temp,nfp,&
!                                         RMNS=rmns_temp, ZMNC=zmnc_temp,&
!                                         BUMNC=bumnc_temp,BVMNC=bvmnc_temp,&
!                                         BUMNS=bumns_temp,BVMNS=bvmns_temp)
!
!                    To load an equilibrium using current density from VMEC
!
!                    ALLOCATE(xm_temp(mnmax),xn_temp(mnmax)) ! Integer Arrays
!                    ALLOCATE(rmnc_temp(mnmax,ns),zmns_temp(mnmax,ns)) ! DOUBLE PRECISION
!                    ALLOCATE(jumnc_temp(mnmax,ns),jvmnc_temp(mnmax,ns)) ! DOUBLE PRECISION
!                    xm_temp = xm
!                    xn_temp = -xn
!                    rmnc_temp = rmnc
!                    zmns_temp = zmns
!                    jumnc_temp = isigng*currumnc
!                    jvmnc_temp = isigng*currvmnc
!                    CALL init_volint(mnmax,nu2,nv2,ns,xm_temp,xn_temp, &
!                                     rmnc_temp,zmns_temp,nfp,&
!                                     JUMNC=jumnc_temp, JVMNC=jvmnc_temp)
!
!                    To evaluate the field at a point in space call
!                    CALL bfield_vc(xp,yp,zp,bxp,byp,bzp,ier)
!                    -or-
!                    CALL vecpot_vc(xp,yp,zp,ax,ay,az,ier)
!
!                    Note that setting adapt_tol < 0 will cause the code
!                    to use non-adaptive integration routines.  Also note
!                    that if the code is compiled without the NTCC PSLPLINE
!                    routines the calls will default to non-adaptive integration.
!                 
!
!
!     Notes:         6/15/12 - SAL
!                    Reorganized the loops for left-hand memory
!                    efficiency.  Realspace matricies are now stored in
!                    [nu,nv,k] order instead of [k,nu,nv]
!                    4/29/13 - SAL
!                    Set minimum number of calls to equal to 0.
!                    7/29/13 - SAL
!                    Minimum number of calls may be set by calling codes
!                    7/30/13 - SAL
!                    Adjusted screen output to reflect non-adaptive
!                    routine selection.
!                    Caught bug in init_virtual_casing_dbl where part of
!                    the nx, ny, nz arrays weren't being calculated
!                    properly.
!                    Added IWRK Parameter to control NAG helper array
!                    size.
!                    1/9/14 - SAL
!                    Benchmarking of code for various configurations
!                    resulted in improved code character.  Surface
!                    normals benchmarked against VMEC ouput (surface
!                    area).  B-Field benchmarked against VMEC.  These
!                    values now appear to be correct.  Although small
!                    errors in B for stellarators and big errors in A
!                    persist.  Benchmarking on flux loops indicates
!                    flux correctly calculated despite A discrepancy.
!                    9/10/14 - SAL
!                    Fixed local deallocations.
!                    12/17/15 - SAL
!                    Changed INTEGER to INTEGER(KIND=8) for gfortran and nag
!                    Change deffinitino of norm_nag to remove factor of nfp
!-----------------------------------------------------------------------
      MODULE virtual_casing_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE EZspline_obj
      USE EZspline
!-----------------------------------------------------------------------
!     Module Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER                          :: nu_vc,nv_vc,nvp,nuvp,nr_vc,nextcur_vc,&
                                          nlastcall
      INTEGER                          :: MIN_CLS = 0
      DOUBLE PRECISION                 :: norm, min_delta_x, x_nag, y_nag, z_nag, &
                                            adapt_tol, norm_nag, norm_3d, adapt_rel, &
                                          surf_area
      DOUBLE PRECISION, ALLOCATABLE    :: xsurf(:),ysurf(:),zsurf(:)
      DOUBLE PRECISION, ALLOCATABLE    :: nxsurf(:),nysurf(:),nzsurf(:)
      DOUBLE PRECISION, ALLOCATABLE    :: btopx(:),btopy(:),btopz(:),btops(:)
      DOUBLE PRECISION, ALLOCATABLE    :: jx_3d(:), jy_3d(:),jz_3d(:)
      DOUBLE PRECISION, PARAMETER      :: pi2 = 6.283185482025146D+00
!      INTEGER, PARAMETER               :: IWRK = 67108864 ! 2^26
      INTEGER, PARAMETER               :: IWRK = 16777216 ! 2^24
      CHARACTER(LEN=256)               :: vc_type_str=''
      DOUBLE PRECISION, PRIVATE, PARAMETER      :: zero = 0.0D+0
      DOUBLE PRECISION, PRIVATE, PARAMETER      :: one  = 1.0D+0

      TYPE(EZspline2_r8)   :: x_spl, y_spl, z_spl
      TYPE(EZspline2_r8)   :: nx_spl, ny_spl, nz_spl
      TYPE(EZspline2_r8)   :: kx_spl, ky_spl, kz_spl, bn_spl
      TYPE(EZspline3_r8)   :: x3d_spl, y3d_spl, z3d_spl, jx3d_spl, jy3d_spl, jz3d_spl

!-----------------------------------------------------------------------
!     Module Subroutines
!          init_virtual_casing    Initializes the atop and surf arrays
!          init_virtual_casing_realspace  Initizliaes the VC algorithm with realspace values
!          virtual_casing_info    Outputs Information to the Screen
!          virtual_casing_dist    Returns distance to surface
!          free_virtual_casing    Frees Allocated Arrays
!          bfield_virtual_casing  Calculates the field
!          vecpot_virtual_casing  Calculates the Vector Potential (A)
!          bfield_virtual_casing_adapt  Calculates the field using NAG adaptive integration
!          vecpot_virtual_casing_adapt  Calculates the Vector Potential using NAG adaptive integration (A)
!-----------------------------------------------------------------------
      INTERFACE init_virtual_casing
         MODULE PROCEDURE init_virtual_casing_flt, init_virtual_casing_dbl
      END INTERFACE
      INTERFACE init_virtual_casing_realspace
         MODULE PROCEDURE init_virtual_casing_realspace_flt, init_virtual_casing_realspace_dbl
      END INTERFACE
      INTERFACE init_volint
         MODULE PROCEDURE init_volint_flt, init_volint_dbl
      END INTERFACE
      INTERFACE bfield_virtual_casing
         MODULE PROCEDURE bfield_virtual_casing_flt, bfield_virtual_casing_dbl
      END INTERFACE
      INTERFACE bfield_virtual_casing_adapt
         MODULE PROCEDURE bfield_virtual_casing_adapt_flt, bfield_virtual_casing_adapt_dbl
      END INTERFACE
      INTERFACE vecpot_virtual_casing
         MODULE PROCEDURE vecpot_virtual_casing_flt, vecpot_virtual_casing_dbl
      END INTERFACE
      INTERFACE vecpot_virtual_casing_adapt
         MODULE PROCEDURE vecpot_virtual_casing_adapt_flt, vecpot_virtual_casing_adapt_dbl
      END INTERFACE
      INTERFACE virtual_casing_dist
         MODULE PROCEDURE virtual_casing_dist_flt, virtual_casing_dist_dbl
      END INTERFACE
      INTERFACE bfield_volint_adapt
         MODULE PROCEDURE bfield_volint_adapt_flt, bfield_volint_adapt_dbl
      END INTERFACE
      INTERFACE vecpot_volint_adapt
         MODULE PROCEDURE vecpot_volint_adapt_flt, vecpot_volint_adapt_dbl
      END INTERFACE
      INTERFACE bfield_vc
         MODULE PROCEDURE bfield_vc_flt, bfield_vc_dbl
      END INTERFACE
      INTERFACE vecpot_vc
         MODULE PROCEDURE vecpot_vc_flt, vecpot_vc_dbl
      END INTERFACE
      CONTAINS
      
      !-----------------------------------------------------------------
      SUBROUTINE init_virtual_casing_dbl(mnmax,nu,nv,xm,xn,rmnc,zmns,nfp,bumnc,bvmnc,&
                                     rmns,zmnc,bsmns,bsmnc,bumns,bvmns,dr)
      IMPLICIT NONE
      ! INPUT VARIABLES
      INTEGER, INTENT(in) :: mnmax, nu, nv, nfp
      INTEGER, INTENT(in) :: xm(1:mnmax),xn(1:mnmax)
      DOUBLE PRECISION, INTENT(in) :: rmnc(1:mnmax,2),zmns(1:mnmax,2)
      DOUBLE PRECISION, INTENT(in), OPTIONAL :: bumnc(1:mnmax,1),bvmnc(1:mnmax,1)
      DOUBLE PRECISION, INTENT(in), OPTIONAL :: bsmns(1:mnmax,1),bsmnc(1:mnmax,1)
      DOUBLE PRECISION, INTENT(in), OPTIONAL :: bumns(1:mnmax,1),bvmns(1:mnmax,1)
      DOUBLE PRECISION, INTENT(in), OPTIONAL :: rmns(1:mnmax,2),zmnc(1:mnmax,2)
      DOUBLE PRECISION, INTENT(in), OPTIONAL :: dr
      
      ! LOCAL VARIABLES
      INTEGER :: nuv, mn, uv, i, u, v, dex1, dex2, ier, nuvm
      INTEGER :: bcs1(2), bcs2(2)
      DOUBLE PRECISION :: snr, snphi,snz,snx,sny,brs,bphis,bzs,bxs,bys,&
                     factor, cop, sip, dx, dy, dz, sn, signs,xt,yt,zt, dr_temp,&
                     br_vac,bphi_vac,bz_vac, xu_temp, xv_temp, yu, yv, stotal
      DOUBLE PRECISION, ALLOCATABLE :: xu(:), xv(:)
      DOUBLE PRECISION, ALLOCATABLE :: fmn_temp(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: xreal(:,:),yreal(:,:),zreal(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: nxreal(:,:),nyreal(:,:),nzreal(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: kxreal(:,:),kyreal(:,:),kzreal(:,:), btopreal(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: r_temp(:,:,:), z_temp(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: bs_temp(:,:,:), bu_temp(:,:,:), bv_temp(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: rs(:,:,:), ru(:,:,:), rv(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: zs(:,:,:), zu(:,:,:), zv(:,:,:)
      ! BEGIN SUBROUTINE
      ! Deallocate anything globally allocated
      IF (ALLOCATED(xsurf)) DEALLOCATE(xsurf)
      IF (ALLOCATED(ysurf)) DEALLOCATE(ysurf)
      IF (ALLOCATED(zsurf)) DEALLOCATE(zsurf)
      IF (ALLOCATED(nxsurf)) DEALLOCATE(nxsurf)
      IF (ALLOCATED(nysurf)) DEALLOCATE(nysurf)
      IF (ALLOCATED(nzsurf)) DEALLOCATE(nzsurf) 
      IF (ALLOCATED(btopx)) DEALLOCATE(btopx)
      IF (ALLOCATED(btopy)) DEALLOCATE(btopy)
      IF (ALLOCATED(btopz)) DEALLOCATE(btopz)
      IF (ALLOCATED(btops)) DEALLOCATE(btops)
      IF (EZspline_allocated(x_spl)) CALL EZspline_free(x_spl,ier)
      IF (EZspline_allocated(y_spl)) CALL EZspline_free(y_spl,ier)
      IF (EZspline_allocated(z_spl)) CALL EZspline_free(z_spl,ier)
      IF (EZspline_allocated(nx_spl)) CALL EZspline_free(nx_spl,ier)
      IF (EZspline_allocated(ny_spl)) CALL EZspline_free(ny_spl,ier)
      IF (EZspline_allocated(nz_spl)) CALL EZspline_free(nz_spl,ier)
      IF (EZspline_allocated(kx_spl)) CALL EZspline_free(kx_spl,ier)
      IF (EZspline_allocated(ky_spl)) CALL EZspline_free(ky_spl,ier)
      IF (EZspline_allocated(kz_spl)) CALL EZspline_free(kz_spl,ier)
      IF (EZspline_allocated(bn_spl)) CALL EZspline_free(bn_spl,ier)

      ! Initialize varaibles
      bcs1=(/-1,-1/)
      bcs2=(/-1,-1/)
      adapt_tol = 1.0D-5
      adapt_rel = 1.0D-3
      factor = pi2 / DBLE(nfp)
      nr_vc = 1
      nu_vc = nu
      nv_vc = nv
      nuv = nu * nv
      nvp = nv * nfp
      nuvp = nu * nv * nfp
      norm = nfp/(2*pi2*nuvp)  !Checked against DIAGMAGNETIC FLUX for 5 FP axisymmetric config
      norm_nag = nfp/(pi2*2)   !CHECKED against DIAGMAGNETIC FLUX for 5 FP axisymmetric config
      vc_type_str='Surface Current'
      ! Alloocations
      ALLOCATE(xu(1:nu),xv(1:nv))
      ALLOCATE(xsurf(1:nuvp),ysurf(1:nuvp),zsurf(1:nuvp))
      ALLOCATE(nxsurf(1:nuvp),nysurf(1:nuvp),nzsurf(1:nuvp))
      ALLOCATE(r_temp(1:nu,1:nv,2),z_temp(1:nu,1:nv,2))
      ALLOCATE(bs_temp(1:nu,1:nv,1),bu_temp(1:nu,1:nv,1),bv_temp(1:nu,1:nv,1))
      ALLOCATE(rs(1:nu,1:nv,1),ru(1:nu,1:nv,1),rv(1:nu,1:nv,1))
      ALLOCATE(zs(1:nu,1:nv,1),zu(1:nu,1:nv,1),zv(1:nu,1:nv,1))
      ALLOCATE(xreal(1:nu+1,1:nvp+1),yreal(1:nu+1,1:nvp+1),zreal(1:nu+1,1:nvp+1))
      ALLOCATE(nxreal(1:nu+1,1:nvp+1),nyreal(1:nu+1,1:nvp+1),nzreal(1:nu+1,1:nvp+1))
      ALLOCATE(kxreal(1:nu+1,1:nvp+1),kyreal(1:nu+1,1:nvp+1),kzreal(1:nu+1,1:nvp+1))
      ALLOCATE(btopreal(1:nu+1,1:nvp+1))
      r_temp=zero; z_temp=zero; bs_temp=zero; bu_temp=zero; bv_temp=zero
      dr_temp=one;
      FORALL(u=1:nu) xu(u) = DBLE(u-1)/DBLE(nu)
      FORALL(v=1:nv) xv(v) = DBLE(v-1)/DBLE(nv)
      CALL mntouv_local(1,2,mnmax,nu,nv,xu,xv,rmnc,xm,xn,r_temp,0,1)
      CALL mntouv_local(1,2,mnmax,nu,nv,xu,xv,zmns,xm,xn,z_temp,1,0)
      IF (PRESENT(bsmns)) CALL mntouv_local(1,1,mnmax,nu,nv,xu,xv,bsmns,xm,xn,bs_temp,1,0)
      IF (PRESENT(bumnc)) CALL mntouv_local(1,1,mnmax,nu,nv,xu,xv,bumnc,xm,xn,bu_temp,0,0)
      IF (PRESENT(bvmnc)) CALL mntouv_local(1,1,mnmax,nu,nv,xu,xv,bvmnc,xm,xn,bv_temp,0,0)
      IF (PRESENT(rmns)) CALL mntouv_local(1,2,mnmax,nu,nv,xu,xv,rmns,xm,xn,r_temp,1,0)
      IF (PRESENT(zmnc)) CALL mntouv_local(1,2,mnmax,nu,nv,xu,xv,zmnc,xm,xn,z_temp,0,0)
      IF (PRESENT(bsmnc)) CALL mntouv_local(1,1,mnmax,nu,nv,xu,xv,bsmnc,xm,xn,bs_temp,0,0)
      IF (PRESENT(bumns)) CALL mntouv_local(1,1,mnmax,nu,nv,xu,xv,bumns,xm,xn,bu_temp,1,0)
      IF (PRESENT(bvmns)) CALL mntouv_local(1,1,mnmax,nu,nv,xu,xv,bvmns,xm,xn,bv_temp,1,0)
      IF (PRESENT(dr)) dr_temp = dr
      
      xreal = zero; yreal = zero; zreal=zero
      uv = 1
      DO v = 1, nv
         DO u = 1, nu
            xsurf(uv) = r_temp(u,v,2)*DCOS(factor*xv(v))
            ysurf(uv) = r_temp(u,v,2)*DSIN(factor*xv(v))
            zsurf(uv) = z_temp(u,v,2)
            xreal(u,v) = xsurf(uv)
            yreal(u,v) = ysurf(uv)
            zreal(u,v) = zsurf(uv)
            uv = uv + 1
         END DO
      END DO
      ! Note we need to extend to nu+1
      xreal(nu+1,:) = xreal(1,:)
      yreal(nu+1,:) = yreal(1,:)
      zreal(nu+1,:) = zreal(1,:)
      
      ! Now we calculate the edge metric elements
      ALLOCATE(fmn_temp(1:mnmax,1))
      rs = zero; zs = zero; ru = zero; zu = zero; rv = zero; zv = zero
      FORALL(mn = 1:mnmax) fmn_temp(mn,1) = -rmnc(mn,2)*xm(mn)
      CALL mntouv_local(1,1,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,ru,1,0)
      FORALL(mn = 1:mnmax) fmn_temp(mn,1) = -rmnc(mn,2)*xn(mn)
      CALL mntouv_local(1,1,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,rv,1,0)  
      FORALL(mn = 1:mnmax) fmn_temp(mn,1) = zmns(mn,2)*xm(mn)
      CALL mntouv_local(1,1,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,zu,0,0) 
      FORALL(mn = 1:mnmax) fmn_temp(mn,1) = zmns(mn,2)*xn(mn)
      CALL mntouv_local(1,1,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,zv,0,0)  
      IF (PRESENT(rmns)) THEN
         FORALL(mn = 1:mnmax) fmn_temp(mn,1) = rmns(mn,2)*xm(mn)
         CALL mntouv_local(1,1,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,ru,0,0)
         FORALL(mn = 1:mnmax) fmn_temp(mn,1) = rmns(mn,2)*xn(mn)
         CALL mntouv_local(1,1,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,rv,0,0)
      END IF
      IF (PRESENT(zmnc)) THEN
         FORALL(mn = 1:mnmax) fmn_temp(mn,1) = -zmnc(mn,2)*xm(mn)
         CALL mntouv_local(1,1,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,zu,1,0)  
         FORALL(mn = 1:mnmax) fmn_temp(mn,1) = -zmnc(mn,2)*xn(mn)
         CALL mntouv_local(1,1,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,zv,1,0)
      END IF
      DO u=1,nu
         DO v=1,nv
            rs(u,v,1)     = (r_temp(u,v,2) - r_temp(u,v,1))
            zs(u,v,1)     = (z_temp(u,v,2) - z_temp(u,v,1))
         END DO
      END DO
      
      ALLOCATE(btopx(1:nuvp), btopy(1:nuvp), btopz(1:nuvp),btops(1:nuvp))
      btopx = zero; btopy = zero; btopz = zero; btops = zero; btopreal = zero
      kxreal = zero; kyreal = zero; kzreal = zero;
      nxreal = zero; nyreal = zero; nzreal = zero;
      uv = 1
      ! Check sign of jacobian
      signs = one
      stotal = zero
      IF (zreal(2,1)-zreal(1,1) < 0) signs = -one
      ! u then v if not outputting
      DO v = 1, nv
         DO u = 1, nu
            cop   = DCOS(factor*xv(v))
            sip   = DSIN(factor*xv(v))
            xu_temp = pi2*ru(u,v,1)*cop
            yu      = pi2*ru(u,v,1)*sip
            xv_temp = pi2*rv(u,v,1)*cop - factor*yreal(u,v)
            yv      = pi2*rv(u,v,1)*sip + factor*xreal(u,v)
            snx     = -yu*zv(u,v,1)*pi2      + zu(u,v,1)*yv*pi2
            sny     = -xv_temp*zu(u,v,1)*pi2 + zv(u,v,1)*xu_temp*pi2
            snz     = -xu_temp*yv        + yu*xv_temp
            sn      = DSQRT(snx*snx+sny*sny+snz*snz)
            stotal = stotal + sn
            brs   = bs_temp(u,v,1)*rs(u,v,1)+bu_temp(u,v,1)*ru(u,v,1)+bv_temp(u,v,1)*rv(u,v,1)*nfp
            bphis = bv_temp(u,v,1)*r_temp(u,v,2)
            bzs   = bs_temp(u,v,1)*zs(u,v,1)+bu_temp(u,v,1)*zu(u,v,1)+bv_temp(u,v,1)*zv(u,v,1)*nfp
            bxs   = brs*cop - bphis*sip
            bys   = brs*sip + bphis*cop
            btopx(uv) = -( bys * snz - bzs * sny )
            btopy(uv) = -( bzs * snx - bxs * snz )
            btopz(uv) = -( bxs * sny - bys * snx )
            kxreal(u,v) = btopx(uv)
            kyreal(u,v) = btopy(uv)
            kzreal(u,v) = btopz(uv)
            nxsurf(uv)  = snx/sn
            nysurf(uv)  = sny/sn
            nzsurf(uv)  = snz/sn
            nxreal(u,v) = snx/sn
            nyreal(u,v) = sny/sn
            nzreal(u,v) = snz/sn
            uv = uv + 1
         END DO
      END DO
      surf_area = nfp*stotal/nuv
      CALL FLUSH(6)
      ! NUMAX
      kxreal(nu+1,:) = kxreal(1,:)
      kyreal(nu+1,:) = kyreal(1,:)
      kzreal(nu+1,:) = kzreal(1,:)
      nxreal(nu+1,:) = nxreal(1,:)
      nyreal(nu+1,:) = nyreal(1,:)
      nzreal(nu+1,:) = nzreal(1,:)
      btopreal(nu+1,:) = btopreal(1,:)
      ! NVMAX
      kxreal(:,nv+1) = kxreal(:,1)
      kyreal(:,nv+1) = kyreal(:,1)
      kzreal(:,nv+1) = kzreal(:,1)
      nxreal(:,nv+1) = nxreal(:,1)
      nyreal(:,nv+1) = nyreal(:,1)
      nzreal(:,nv+1) = nzreal(:,1)
      btopreal(:,nv+1) = btopreal(:,1)
      ! Last point
      kxreal(nu+1,nv+1) = kxreal(1,1)
      kyreal(nu+1,nv+1) = kyreal(1,1)
      kzreal(nu+1,nv+1) = kzreal(1,1)
      nxreal(nu+1,nv+1) = nxreal(1,1)
      nyreal(nu+1,nv+1) = nyreal(1,1)
      nzreal(nu+1,nv+1) = nzreal(1,1)
      btopreal(nu+1,nv+1) = btopreal(nu+1,nv+1)
      ! Finish off
      DO v = 2, nfp
         cop  = DCOS((v-1)*factor)
         sip  = DSIN((v-1)*factor)
         dex1 = (v-1)*nuv+1
         dex2 = v*nuv
         btopx(dex1:dex2) = btopx(1:nuv)*cop - btopy(1:nuv)*sip
         btopy(dex1:dex2) = btopy(1:nuv)*cop + btopx(1:nuv)*sip
         btopz(dex1:dex2) = btopz(1:nuv)
         btops(dex1:dex2) = btops(1:nuv)
         xsurf(dex1:dex2) = xsurf(1:nuv)*cop - ysurf(1:nuv)*sip
         ysurf(dex1:dex2) = ysurf(1:nuv)*cop + xsurf(1:nuv)*sip
         zsurf(dex1:dex2) = zsurf(1:nuv)
         nxsurf(dex1:dex2) = nxsurf(1:nuv)*cop - nysurf(1:nuv)*sip
         nysurf(dex1:dex2) = nysurf(1:nuv)*cop + nxsurf(1:nuv)*sip
         nzsurf(dex1:dex2) = nzsurf(1:nuv)
         dex1 = (v-1)*nv+1
         dex2 = v*nv
         xreal(1:nu+1,dex1:dex2) = xreal(1:nu+1,1:nv)*cop - yreal(1:nu+1,1:nv)*sip
         yreal(1:nu+1,dex1:dex2) = yreal(1:nu+1,1:nv)*cop + xreal(1:nu+1,1:nv)*sip
         zreal(1:nu+1,dex1:dex2) = zreal(1:nu+1,1:nv)
         kxreal(1:nu+1,dex1:dex2) = kxreal(1:nu+1,1:nv)*cop - kyreal(1:nu+1,1:nv)*sip
         kyreal(1:nu+1,dex1:dex2) = kyreal(1:nu+1,1:nv)*cop + kxreal(1:nu+1,1:nv)*sip
         kzreal(1:nu+1,dex1:dex2) = kzreal(1:nu+1,1:nv)
         nxreal(1:nu+1,dex1:dex2) = nxreal(1:nu+1,1:nv)*cop - nyreal(1:nu+1,1:nv)*sip
         nyreal(1:nu+1,dex1:dex2) = nyreal(1:nu+1,1:nv)*cop + nxreal(1:nu+1,1:nv)*sip
         nzreal(1:nu+1,dex1:dex2) = nzreal(1:nu+1,1:nv)
         btopreal(1:nu+1,dex1:dex2) = btopreal(1:nu+1,1:nv)
      END DO
      xreal(:,nvp+1)  = xreal(:,1)
      yreal(:,nvp+1)  = yreal(:,1)
      zreal(:,nvp+1)  = zreal(:,1)
      kxreal(:,nvp+1) = kxreal(:,1)
      kyreal(:,nvp+1) = kyreal(:,1)
      kzreal(:,nvp+1) = kzreal(:,1)
      nxreal(:,nvp+1) = nxreal(:,1)
      nyreal(:,nvp+1) = nyreal(:,1)
      nzreal(:,nvp+1) = nzreal(:,1)
      btopreal(:,nvp+1) = btopreal(:,1)
      ! Calculate Return map
      min_delta_x = 1.0D+10
      DO v = 2, nv-1
         DO u = 2, nu-1
             dx = xreal(nu+1,nv) - xreal(nu-1,nv)
             dy = yreal(nu+1,nv) - yreal(nu-1,nv)
             dz = zreal(nu+1,nv) - zreal(nu-1,nv)
             min_delta_x=MIN(min_delta_x,SQRT(dx*dx+dy*dy+dz*dz))
             dx = xreal(nu,nv+1) - xreal(nu,nv-1)
             dy = yreal(nu,nv+1) - yreal(nu,nv-1)
             dz = zreal(nu,nv+1) - zreal(nu,nv-1)
             min_delta_x=MIN(min_delta_x,SQRT(dx*dx+dy*dy+dz*dz))
         END DO
      END DO
      ! Construct Splines
      CALL EZspline_init(x_spl,nu+1,nvp+1,bcs1,bcs2,ier)
      CALL EZspline_init(y_spl,nu+1,nvp+1,bcs1,bcs2,ier)
      CALL EZspline_init(z_spl,nu+1,nvp+1,bcs1,bcs2,ier)
      CALL EZspline_init(nx_spl,nu+1,nvp+1,bcs1,bcs2,ier)
      CALL EZspline_init(ny_spl,nu+1,nvp+1,bcs1,bcs2,ier)
      CALL EZspline_init(nz_spl,nu+1,nvp+1,bcs1,bcs2,ier)
      CALL EZspline_init(kx_spl,nu+1,nvp+1,bcs1,bcs2,ier)
      CALL EZspline_init(ky_spl,nu+1,nvp+1,bcs1,bcs2,ier)
      CALL EZspline_init(kz_spl,nu+1,nvp+1,bcs1,bcs2,ier)
      CALL EZspline_init(bn_spl,nu+1,nvp+1,bcs1,bcs2,ier)
      ! Define arrays from 0 to 1 because integrand already contains
      ! dA information.
      DO  u = 1, nu+1
         x_spl%x1(u) = DBLE(u-1)/DBLE(nu)
         y_spl%x1(u) = DBLE(u-1)/DBLE(nu)
         z_spl%x1(u) = DBLE(u-1)/DBLE(nu)
         nx_spl%x1(u) = DBLE(u-1)/DBLE(nu)
         ny_spl%x1(u) = DBLE(u-1)/DBLE(nu)
         nz_spl%x1(u) = DBLE(u-1)/DBLE(nu)
         kx_spl%x1(u) = DBLE(u-1)/DBLE(nu)
         ky_spl%x1(u) = DBLE(u-1)/DBLE(nu)
         kz_spl%x1(u) = DBLE(u-1)/DBLE(nu)
         bn_spl%x1(u) = DBLE(u-1)/DBLE(nu)
      END DO
      DO  v = 1, nvp+1
         x_spl%x2(v) = DBLE(v-1)/DBLE(nvp)
         y_spl%x2(v) = DBLE(v-1)/DBLE(nvp)
         z_spl%x2(v) = DBLE(v-1)/DBLE(nvp)
         nx_spl%x2(v) = DBLE(v-1)/DBLE(nvp)
         ny_spl%x2(v) = DBLE(v-1)/DBLE(nvp)
         nz_spl%x2(v) = DBLE(v-1)/DBLE(nvp)
         kx_spl%x2(v) = DBLE(v-1)/DBLE(nvp)
         ky_spl%x2(v) = DBLE(v-1)/DBLE(nvp)
         kz_spl%x2(v) = DBLE(v-1)/DBLE(nvp)
         bn_spl%x2(v) = DBLE(v-1)/DBLE(nvp)
      END DO
      x_spl%isHermite  = 1
      y_spl%isHermite  = 1
      z_spl%isHermite  = 1
      nx_spl%isHermite = 1
      ny_spl%isHermite = 1
      nz_spl%isHermite = 1
      kx_spl%isHermite = 1
      ky_spl%isHermite = 1
      kz_spl%isHermite = 1
      bn_spl%isHermite = 1
      CALL EZspline_setup(x_spl,xreal,ier)
      CALL EZspline_setup(y_spl,yreal,ier)
      CALL EZspline_setup(z_spl,zreal,ier)
      CALL EZspline_setup(nx_spl,nxreal,ier)
      CALL EZspline_setup(ny_spl,nyreal,ier)
      CALL EZspline_setup(nz_spl,nzreal,ier)
      CALL EZspline_setup(kx_spl,kxreal,ier)
      CALL EZspline_setup(ky_spl,kyreal,ier)
      CALL EZspline_setup(kz_spl,kzreal,ier)
      CALL EZspline_setup(bn_spl,btopreal,ier)
      DEALLOCATE(xu,xv)
      DEALLOCATE(r_temp,z_temp)
      DEALLOCATE(bs_temp,bu_temp,bv_temp)
      DEALLOCATE(rs,ru,rv)
      DEALLOCATE(zs,zu,zv)
      DEALLOCATE(xreal,yreal,zreal)
      DEALLOCATE(nxreal,nyreal,nzreal)
      DEALLOCATE(kxreal,kyreal,kzreal,btopreal)
      DEALLOCATE(fmn_temp)
      ! END SUBROUTINE
      END SUBROUTINE init_virtual_casing_dbl
      !-----------------------------------------------------------------
         
      !-----------------------------------------------------------------
      SUBROUTINE init_virtual_casing_realspace_dbl(nu,nv,nfp,phi,&
                                                   rreal,zreal,&
                                                   snr,snphi,snz,&
                                                   brreal,bphireal,bzreal)
      IMPLICIT NONE
      ! INPUT VARIABLES
      INTEGER, INTENT(in) :: nu, nv, nfp
      DOUBLE PRECISION, INTENT(in) :: phi(nv)
      DOUBLE PRECISION, INTENT(in) :: rreal(nu,nv), zreal(nu,nv)
      DOUBLE PRECISION, INTENT(in) :: snr(nu,nv), snphi(nu,nv), snz(nu,nv)
      DOUBLE PRECISION, INTENT(in) :: brreal(nu,nv), bphireal(nu,nv), bzreal(nu,nv)
      ! LOCAL VARIABLES
      INTEGER :: nuv, mn, uv, i, u, v, dex1, dex2, ier, nuvm
      INTEGER :: bcs1(2), bcs2(2)
      DOUBLE PRECISION :: snx,sny,bzs,bxs,bys,&
                     factor, cop, sip, dx, dy, dz, sn, signs,xt,yt,zt, dr_temp
      DOUBLE PRECISION, ALLOCATABLE :: fmn_temp(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: xreal(:,:),yreal(:,:),zreal2(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: nxreal(:,:),nyreal(:,:),nzreal(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: kxreal(:,:),kyreal(:,:),kzreal(:,:), btopreal(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: r_temp(:,:,:), z_temp(:,:,:)
      ! BEGIN SUBROUTINE
      ! Deallocate anything globally allocated
      IF (ALLOCATED(xsurf)) DEALLOCATE(xsurf)
      IF (ALLOCATED(ysurf)) DEALLOCATE(ysurf)
      IF (ALLOCATED(zsurf)) DEALLOCATE(zsurf)
      IF (ALLOCATED(nxsurf)) DEALLOCATE(nxsurf)
      IF (ALLOCATED(nysurf)) DEALLOCATE(nysurf)
      IF (ALLOCATED(nzsurf)) DEALLOCATE(nzsurf) 
      IF (ALLOCATED(btopx)) DEALLOCATE(btopx)
      IF (ALLOCATED(btopy)) DEALLOCATE(btopy)
      IF (ALLOCATED(btopz)) DEALLOCATE(btopz)
      IF (ALLOCATED(btops)) DEALLOCATE(btops)
      IF (EZspline_allocated(x_spl)) CALL EZspline_free(x_spl,ier)
      IF (EZspline_allocated(y_spl)) CALL EZspline_free(y_spl,ier)
      IF (EZspline_allocated(z_spl)) CALL EZspline_free(z_spl,ier)
      IF (EZspline_allocated(nx_spl)) CALL EZspline_free(nx_spl,ier)
      IF (EZspline_allocated(ny_spl)) CALL EZspline_free(ny_spl,ier)
      IF (EZspline_allocated(nz_spl)) CALL EZspline_free(nz_spl,ier)
      IF (EZspline_allocated(kx_spl)) CALL EZspline_free(kx_spl,ier)
      IF (EZspline_allocated(ky_spl)) CALL EZspline_free(ky_spl,ier)
      IF (EZspline_allocated(kz_spl)) CALL EZspline_free(kz_spl,ier)
      IF (EZspline_allocated(bn_spl)) CALL EZspline_free(bn_spl,ier)
      ! Initialize varaibles
      bcs1=(/-1,-1/)
      bcs2=(/-1,-1/)
      adapt_tol = 1.0D-5
      adapt_rel = 1.0D-3
      factor = pi2 / DBLE(nfp)
      nu_vc = nu
      nv_vc = nv
      nuv = nu * nv
      nvp = nv * nfp
      nuvp = nu * nv * nfp
      norm = one/(2*pi2*nu*nv)
      norm_nag = one/(2*pi2)
      vc_type_str='Surface Current (REALSPACE)'
      ! Alloocations
      ALLOCATE(xsurf(1:nuvp),ysurf(1:nuvp),zsurf(1:nuvp))
      ALLOCATE(nxsurf(1:nuvp),nysurf(1:nuvp),nzsurf(1:nuvp))
      ALLOCATE(r_temp(2,1:nu,1:nv),z_temp(2,1:nu,1:nv))
      ALLOCATE(xreal(1:nu+1,1:nvp+1),yreal(1:nu+1,1:nvp+1),zreal2(1:nu+1,1:nvp+1))
      ALLOCATE(nxreal(1:nu+1,1:nvp+1),nyreal(1:nu+1,1:nvp+1),nzreal(1:nu+1,1:nvp+1))
      ALLOCATE(kxreal(1:nu+1,1:nvp+1),kyreal(1:nu+1,1:nvp+1),kzreal(1:nu+1,1:nvp+1))
      ALLOCATE(btopreal(1:nu+1,1:nvp+1))
      dr_temp=one;
      uv = 1
      DO v = 1, nv
         DO u = 1, nu
            xsurf(uv)   = rreal(u,v)*dcos(phi(nv))
            ysurf(uv)   = rreal(u,v)*dsin(phi(nv))
            zsurf(uv)   = zreal(u,v)
            xreal(u,v)  = xsurf(uv)
            yreal(u,v)  = ysurf(uv)
            zreal2(u,v) = zsurf(uv)
            uv = uv + 1
         END DO
      END DO
      ! Note we need to extend to nu+1
      xreal(nu+1,:) = xreal(1,:)
      yreal(nu+1,:) = yreal(1,:)
      zreal2(nu+1,:) = zreal2(1,:)
      
      ALLOCATE(btopx(1:nuvp), btopy(1:nuvp), btopz(1:nuvp),btops(1:nuvp))
      btopx = zero; btopy = zero; btopz = zero; btops = zero; btopreal = zero
      uv = 1
      ! Check sign of jacobian
      signs = -one
      IF (zreal2(2,1) < 0) signs = 1.0
      ! u then v if not outputting
      DO v = 1, nv
         DO u = 1, nu
            cop   = dcos(phi(v))
            sip   = dsin(phi(v))
            snx   = snr(u,v)*cop - snphi(u,v)*sip
            sny   = snr(u,v)*sip + snphi(u,v)*cop
            sn    = DSQRT(snx*snx+sny*sny+snz(u,v)*snz(u,v))
            bxs   = brreal(u,v)*cop - bphireal(u,v)*sip
            bys   = brreal(u,v)*sip + bphireal(u,v)*cop
            bzs   = bzreal(u,v)
            btopx(uv) = signs * ( bys * snz(u,v) - bzs * sny )
            btopy(uv) = signs * ( bzs * snx - bxs * snz(u,v) )
            btopz(uv) = signs * ( bxs * sny - bys * snx )
            kxreal(u,v) = btopx(uv)
            kyreal(u,v) = btopy(uv)
            kzreal(u,v) = btopz(uv)
            nxsurf(uv)  = snx/sn
            nysurf(uv)  = sny/sn
            nzsurf(uv)  = snz(u,v)/sn
            nxreal(u,v) = snx/sn
            nyreal(u,v) = sny/sn
            nzreal(u,v) = snz(u,v)/sn
            btops(uv) = bxs*snx+bys*sny+bzs*snz(u,v)
            btopreal(u,v) = btops(uv)
            uv = uv + 1
         END DO
      END DO
      kxreal(nu+1,:) = kxreal(1,:)
      kyreal(nu+1,:) = kyreal(1,:)
      kzreal(nu+1,:) = kzreal(1,:)
      nxreal(nu+1,:) = nxreal(1,:)
      nyreal(nu+1,:) = nyreal(1,:)
      nzreal(nu+1,:) = nzreal(1,:)
      DO v = 2, nfp
         cop  = dcos((v-1)*factor)
         sip  = dsin((v-1)*factor)
         dex1 = (v-1)*nuv+1
         dex2 = v*nuv
         btopx(dex1:dex2) = btopx(1:nuv)*cop - btopy(1:nuv)*sip
         btopy(dex1:dex2) = btopy(1:nuv)*cop + btopx(1:nuv)*sip
         btopz(dex1:dex2) = btopz(1:nuv)
         btops(dex1:dex2) = btops(1:nuv)
         xsurf(dex1:dex2) = xsurf(1:nuv)*cop - ysurf(1:nuv)*sip
         ysurf(dex1:dex2) = ysurf(1:nuv)*cop + xsurf(1:nuv)*sip
         zsurf(dex1:dex2) = zsurf(1:nuv)
         nxsurf(dex1:dex2) = nxsurf(1:nuv)*cop - nysurf(1:nuv)*sip
         nysurf(dex1:dex2) = nysurf(1:nuv)*cop + nxsurf(1:nuv)*sip
         nzsurf(dex1:dex2) = nzsurf(1:nuv)
         dex1 = (v-1)*nv+1
         dex2 = v*nv
         xreal(1:nu+1,dex1:dex2) = xreal(1:nu+1,1:nv)*cop - yreal(1:nu+1,1:nv)*sip
         yreal(1:nu+1,dex1:dex2) = yreal(1:nu+1,1:nv)*cop + xreal(1:nu+1,1:nv)*sip
         zreal2(1:nu+1,dex1:dex2) = zreal2(1:nu+1,1:nv)
         kxreal(1:nu+1,dex1:dex2) = kxreal(1:nu+1,1:nv)*cop - kyreal(1:nu+1,1:nv)*sip
         kyreal(1:nu+1,dex1:dex2) = kyreal(1:nu+1,1:nv)*cop + kxreal(1:nu+1,1:nv)*sip
         kzreal(1:nu+1,dex1:dex2) = kzreal(1:nu+1,1:nv)
         nxreal(1:nu+1,dex1:dex2) = nxreal(1:nu+1,1:nv)*cop - nyreal(1:nu+1,1:nv)*sip
         nyreal(1:nu+1,dex1:dex2) = nyreal(1:nu+1,1:nv)*cop + nxreal(1:nu+1,1:nv)*sip
         nzreal(1:nu+1,dex1:dex2) = nzreal(1:nu+1,1:nv)
         btopreal(1:nu+1,dex1:dex2) = btopreal(1:nu+1,1:nv)
      END DO
      xreal(:,nvp+1)  = xreal(:,1)
      yreal(:,nvp+1)  = yreal(:,1)
      zreal2(:,nvp+1)  = zreal2(:,1)
      kxreal(:,nvp+1) = kxreal(:,1)
      kyreal(:,nvp+1) = kyreal(:,1)
      kzreal(:,nvp+1) = kzreal(:,1)
      nxreal(:,nvp+1) = nxreal(:,1)
      nyreal(:,nvp+1) = nyreal(:,1)
      nzreal(:,nvp+1) = nzreal(:,1)
      btopreal(:,nvp+1) = btopreal(:,1)
      ! Calculate Return map
      min_delta_x = 1.0D+10
      DO v = 2, nv-1
         DO u = 2, nu-1
             dx = xreal(nu+1,nv) - xreal(nu-1,nv)
             dy = yreal(nu+1,nv) - yreal(nu-1,nv)
             dz = zreal2(nu+1,nv) - zreal2(nu-1,nv)
             min_delta_x=MIN(min_delta_x,SQRT(dx*dx+dy*dy+dz*dz))
             dx = xreal(nu,nv+1) - xreal(nu,nv-1)
             dy = yreal(nu,nv+1) - yreal(nu,nv-1)
             dz = zreal2(nu,nv+1) - zreal2(nu,nv-1)
             min_delta_x=MIN(min_delta_x,SQRT(dx*dx+dy*dy+dz*dz))
         END DO
      END DO
      ! Construct Splines
      CALL EZspline_init(x_spl,nu+1,nvp+1,bcs1,bcs2,ier)
      CALL EZspline_init(y_spl,nu+1,nvp+1,bcs1,bcs2,ier)
      CALL EZspline_init(z_spl,nu+1,nvp+1,bcs1,bcs2,ier)
      CALL EZspline_init(nx_spl,nu+1,nvp+1,bcs1,bcs2,ier)
      CALL EZspline_init(ny_spl,nu+1,nvp+1,bcs1,bcs2,ier)
      CALL EZspline_init(nz_spl,nu+1,nvp+1,bcs1,bcs2,ier)
      CALL EZspline_init(kx_spl,nu+1,nvp+1,bcs1,bcs2,ier)
      CALL EZspline_init(ky_spl,nu+1,nvp+1,bcs1,bcs2,ier)
      CALL EZspline_init(kz_spl,nu+1,nvp+1,bcs1,bcs2,ier)
      CALL EZspline_init(bn_spl,nu+1,nvp+1,bcs1,bcs2,ier)
      ! Define arrays from 0 to 1 because integrand already contains
      ! dA information.
      DO  u = 1, nu+1
         x_spl%x1(u) = DBLE(u-1)/DBLE(nu)
         y_spl%x1(u) = DBLE(u-1)/DBLE(nu)
         z_spl%x1(u) = DBLE(u-1)/DBLE(nu)
         nx_spl%x1(u) = DBLE(u-1)/DBLE(nu)
         ny_spl%x1(u) = DBLE(u-1)/DBLE(nu)
         nz_spl%x1(u) = DBLE(u-1)/DBLE(nu)
         kx_spl%x1(u) = DBLE(u-1)/DBLE(nu)
         ky_spl%x1(u) = DBLE(u-1)/DBLE(nu)
         kz_spl%x1(u) = DBLE(u-1)/DBLE(nu)
         bn_spl%x1(u) = DBLE(u-1)/DBLE(nu)
      END DO
      DO  v = 1, nvp+1
         x_spl%x2(v) = DBLE(v-1)/DBLE(nvp)
         y_spl%x2(v) = DBLE(v-1)/DBLE(nvp)
         z_spl%x2(v) = DBLE(v-1)/DBLE(nvp)
         nx_spl%x2(v) = DBLE(v-1)/DBLE(nvp)
         ny_spl%x2(v) = DBLE(v-1)/DBLE(nvp)
         nz_spl%x2(v) = DBLE(v-1)/DBLE(nvp)
         kx_spl%x2(v) = DBLE(v-1)/DBLE(nvp)
         ky_spl%x2(v) = DBLE(v-1)/DBLE(nvp)
         kz_spl%x2(v) = DBLE(v-1)/DBLE(nvp)
         bn_spl%x2(v) = DBLE(v-1)/DBLE(nvp)
      END DO
      x_spl%isHermite  = 0
      y_spl%isHermite  = 0
      z_spl%isHermite  = 0
      nx_spl%isHermite = 0
      ny_spl%isHermite = 0
      nz_spl%isHermite = 0
      kx_spl%isHermite = 0
      ky_spl%isHermite = 0
      kz_spl%isHermite = 0
      bn_spl%isHermite = 0
      CALL EZspline_setup(x_spl,xreal,ier)
      CALL EZspline_setup(y_spl,yreal,ier)
      CALL EZspline_setup(z_spl,zreal2,ier)
      CALL EZspline_setup(kx_spl,kxreal,ier)
      CALL EZspline_setup(ky_spl,kyreal,ier)
      CALL EZspline_setup(kz_spl,kzreal,ier)
      CALL EZspline_setup(nx_spl,nxreal,ier)
      CALL EZspline_setup(ny_spl,nyreal,ier)
      CALL EZspline_setup(nz_spl,nzreal,ier)
      CALL EZspline_setup(bn_spl,btopreal,ier)

      DEALLOCATE(r_temp,z_temp)
      DEALLOCATE(xreal,yreal,zreal2)
      DEALLOCATE(kxreal,kyreal,kzreal,btopreal)
      DEALLOCATE(nxreal,nyreal,nzreal)
      ! END SUBROUTINE
      END SUBROUTINE init_virtual_casing_realspace_dbl
      !-----------------------------------------------------------------        
      
      !-----------------------------------------------------------------
      !  Note optional arguments must have a different name so the
      !  Interface will pick the proper variables.  Here:  _FLT
      SUBROUTINE init_virtual_casing_flt(mnmax_flt,nu_flt,nv_flt,xm_flt,&
                                         xn_flt,rmnc_flt,zmns_flt,nfp_flt,&
                                         bumnc_flt,bvmnc_flt,&
                                         rmns_flt,zmnc_flt,bsmns_flt,&
                                         bsmnc_flt,bumns_flt,bvmns_flt,&
                                         dr_flt)
      IMPLICIT NONE
      ! INPUT VARIABLES
      INTEGER, INTENT(in) :: mnmax_flt, nu_flt, nv_flt, nfp_flt
      INTEGER, INTENT(in) :: xm_flt(1:mnmax_flt),xn_flt(1:mnmax_flt)
      REAL, INTENT(in) :: rmnc_flt(1:mnmax_flt,2),zmns_flt(1:mnmax_flt,2)
      REAL, INTENT(in), OPTIONAL :: bsmns_flt(1:mnmax_flt,1),bsmnc_flt(1:mnmax_flt,1)
      REAL, INTENT(in), OPTIONAL :: bumnc_flt(1:mnmax_flt,1),bvmnc_flt(1:mnmax_flt,1)
      REAL, INTENT(in), OPTIONAL :: bumns_flt(1:mnmax_flt,1),bvmns_flt(1:mnmax_flt,1)
      REAL, INTENT(in), OPTIONAL :: rmns_flt(1:mnmax_flt,2),zmnc_flt(1:mnmax_flt,2)
      REAL, INTENT(in), OPTIONAL :: dr_flt
      ! LOCAL VARIABLES
      DOUBLE PRECISION :: rmnct(1:mnmax_flt,2),zmnst(1:mnmax_flt,2)
      DOUBLE PRECISION :: rmnst(1:mnmax_flt,2),zmnct(1:mnmax_flt,2)
      DOUBLE PRECISION :: bumnct(1:mnmax_flt,1),bvmnct(1:mnmax_flt,1)
      DOUBLE PRECISION :: bsmnst(1:mnmax_flt,1),bsmnct(1:mnmax_flt,1)
      DOUBLE PRECISION :: bumnst(1:mnmax_flt,1),bvmnst(1:mnmax_flt,1)
      DOUBLE PRECISION :: drt
      ! BEGIN SUBROUTINE
      rmnct = zero; rmnst = zero; zmnct = zero; zmnst = zero
      bsmnct = zero; bsmnst = zero
      bumnct = zero; bumnst = zero
      bvmnct = zero; bvmnst = zero
      drt    = one
      rmnct  = rmnc_flt
      zmnst  = zmns_flt
      IF (PRESENT(rmns_flt))  rmnst = rmns_flt
      IF (PRESENT(zmnc_flt))  zmnct = zmnc_flt
      IF (PRESENT(bsmnc_flt)) bsmnct = bsmnc_flt
      IF (PRESENT(bsmns_flt)) bsmnst = bsmns_flt
      IF (PRESENT(bumnc_flt)) bumnct = bumnc_flt
      IF (PRESENT(bumns_flt)) bumnst = bumns_flt
      IF (PRESENT(bvmnc_flt)) bvmnct = bvmnc_flt
      IF (PRESENT(bvmns_flt)) bvmnst = bvmns_flt
      IF (PRESENT(dr_flt))    drt    = dr_flt
      CALL init_virtual_casing_dbl(mnmax_flt,nu_flt,nv_flt,xm_flt,xn_flt, &
                                   rmnct,zmnst,nfp_flt, &
                                   RMNS=rmnst,ZMNC=zmnct, &
                                   BSMNS=bsmnst, BSMNC=bsmnct, &
                                   BUMNC=bumnct, BUMNS=bumnst, &
                                   BVMNC=bvmnct, BVMNS=bvmnst, &
                                   DR=drt)
      ! END SUBROUTINE
      END SUBROUTINE init_virtual_casing_flt
      !-----------------------------------------------------------------        
      
      !-----------------------------------------------------------------
      !  Note optional arguments must have a different name so the
      !  Interface will pick the proper variables.  Here:  _FLT
      SUBROUTINE init_virtual_casing_realspace_flt(nu_flt,nv_flt,nfp_flt,phi_flt,&
                                                   rreal_flt,zreal_flt,&
                                                   snr_flt,snphi_flt,snz_flt,&
                                                   brreal_flt,bphireal_flt,bzreal_flt)
      IMPLICIT NONE
      ! INPUT VARIABLES
      INTEGER, INTENT(in) :: nu_flt, nv_flt, nfp_flt
      REAL, INTENT(in) :: phi_flt(nv_flt)
      REAL, INTENT(in) :: rreal_flt(nu_flt,nv_flt), zreal_flt(nu_flt,nv_flt)
      REAL, INTENT(in) :: snr_flt(nu_flt,nv_flt), snphi_flt(nu_flt,nv_flt), snz_flt(nu_flt,nv_flt)
      REAL, INTENT(in) :: brreal_flt(nu_flt,nv_flt), bphireal_flt(nu_flt,nv_flt), bzreal_flt(nu_flt,nv_flt)
      ! LOCAL VARIABLES
      DOUBLE PRECISION :: phi_dbl(nv_flt)
      DOUBLE PRECISION :: rreal_dbl(nu_flt,nv_flt), zreal_dbl(nu_flt,nv_flt)
      DOUBLE PRECISION :: snr_dbl(nu_flt,nv_flt), snphi_dbl(nu_flt,nv_flt), snz_dbl(nu_flt,nv_flt)
      DOUBLE PRECISION :: brreal_dbl(nu_flt,nv_flt), bphireal_dbl(nu_flt,nv_flt), bzreal_dbl(nu_flt,nv_flt)
      ! BEGIN SUBROUTINE
      phi_dbl      = phi_flt
      rreal_dbl    = rreal_flt
      zreal_dbl    = zreal_flt
      snr_dbl      = snr_flt
      snphi_dbl    = snphi_flt
      snz_dbl      = snz_flt
      brreal_dbl   = brreal_flt
      bphireal_dbl = bphireal_flt
      bzreal_dbl   = bzreal_flt
      CALL init_virtual_casing_realspace_dbl(nu_flt,nv_flt,nfp_flt,phi_dbl,&
                                             rreal_dbl,zreal_dbl,&
                                             snr_dbl,snphi_dbl,snz_dbl,&
                                             brreal_dbl,bphireal_dbl,bzreal_dbl)
      ! END SUBROUTINE
      END SUBROUTINE init_virtual_casing_realspace_flt
      !-----------------------------------------------------------------
      
         
      !-----------------------------------------------------------------
      SUBROUTINE free_virtual_casing
      IMPLICIT NONE
      INTEGER :: ier
      ! BEGIN SUBROUTINE
      nuvp = 0
      norm = -one
      IF (ALLOCATED(xsurf))  DEALLOCATE(xsurf, ysurf, zsurf)
      IF (ALLOCATED(nxsurf)) DEALLOCATE(nxsurf, nysurf, nzsurf)
      IF (ALLOCATED(btopx))  DEALLOCATE(btopx, btopy, btopz, btops)
      IF (ALLOCATED(jx_3d))  DEALLOCATE(jx_3d, jy_3d, jz_3d)

      IF (EZspline_allocated(x_spl)) CALL EZspline_free(x_spl,ier)
      IF (EZspline_allocated(y_spl)) CALL EZspline_free(y_spl,ier)
      IF (EZspline_allocated(z_spl)) CALL EZspline_free(z_spl,ier)
      IF (EZspline_allocated(nx_spl)) CALL EZspline_free(nx_spl,ier)
      IF (EZspline_allocated(ny_spl)) CALL EZspline_free(ny_spl,ier)
      IF (EZspline_allocated(nz_spl)) CALL EZspline_free(nz_spl,ier)
      IF (EZspline_allocated(kx_spl)) CALL EZspline_free(kx_spl,ier)
      IF (EZspline_allocated(ky_spl)) CALL EZspline_free(ky_spl,ier)
      IF (EZspline_allocated(kz_spl)) CALL EZspline_free(kz_spl,ier)
      IF (EZspline_allocated(bn_spl)) CALL EZspline_free(bn_spl,ier)
      IF (EZspline_allocated(x3d_spl)) CALL EZspline_free(x3d_spl,ier)
      IF (EZspline_allocated(y3d_spl)) CALL EZspline_free(y3d_spl,ier)
      IF (EZspline_allocated(z3d_spl)) CALL EZspline_free(z3d_spl,ier)
      IF (EZspline_allocated(jx3d_spl)) CALL EZspline_free(jx3d_spl,ier)
      IF (EZspline_allocated(jy3d_spl)) CALL EZspline_free(jy3d_spl,ier)
      IF (EZspline_allocated(jz3d_spl)) CALL EZspline_free(jz3d_spl,ier)

      ! END SUBROUTINE
      END SUBROUTINE free_virtual_casing                        
      !-----------------------------------------------------------------
         
         
      !-----------------------------------------------------------------
      SUBROUTINE bfield_virtual_casing_dbl(x,y,z,bx,by,bz)
      IMPLICIT NONE
      ! INPUT VARIABLES
      DOUBLE PRECISION, INTENT(in)  :: x, y, z
      DOUBLE PRECISION, INTENT(out) :: bx, by, bz
      ! LOCAL VARIABLES
      DOUBLE PRECISION ::  gf(1:nuvp), gf3(1:nuvp)
      ! BEGIN SUBROUTINE
      gf(1:nuvp) = one/dsqrt(  (x-xsurf(1:nuvp))**2 &
                                 + (y-ysurf(1:nuvp))**2 &
                                 + (z-zsurf(1:nuvp))**2)
      gf3 = gf*gf*gf
      bx  = norm*sum( (  btopy(1:nuvp)*(z-zsurf(1:nuvp)) &
                       - btopz(1:nuvp)*(y-ysurf(1:nuvp)) ) * gf3(1:nuvp))
      by  = norm*sum( (  btopz(1:nuvp)*(x-xsurf(1:nuvp)) &
                       - btopx(1:nuvp)*(z-zsurf(1:nuvp)) ) * gf3(1:nuvp))
      bz  = norm*sum( (  btopx(1:nuvp)*(y-ysurf(1:nuvp)) &
                       - btopy(1:nuvp)*(x-xsurf(1:nuvp)) ) * gf3(1:nuvp))
      bx  = bx + norm*sum((btops(1:nuvp)*(x-xsurf(1:nuvp)))*gf3(1:nuvp))
      by  = by + norm*sum((btops(1:nuvp)*(y-ysurf(1:nuvp)))*gf3(1:nuvp))
      bz  = bz + norm*sum((btops(1:nuvp)*(z-zsurf(1:nuvp)))*gf3(1:nuvp))
      RETURN
      ! END SUBROUTINE
      END SUBROUTINE bfield_virtual_casing_dbl
      !-----------------------------------------------------------------
         
         
      !-----------------------------------------------------------------
      SUBROUTINE bfield_virtual_casing_flt(x_flt,y_flt,z_flt,bx_flt,by_flt,bz_flt)
      IMPLICIT NONE
      ! INPUT VARIABLES
      REAL, INTENT(in)  :: x_flt, y_flt, z_flt
      REAL, INTENT(out) :: bx_flt, by_flt, bz_flt
      ! LOCAL VARIABLES
      DOUBLE PRECISION :: xt,yt,zt,bxt,byt,bzt
      ! BEGIN SUBROUTINE
      xt  = x_flt
      yt  = y_flt
      zt  = z_flt
      bxt = zero
      byt = zero
      bzt = zero
      CALL bfield_virtual_casing_dbl(xt,yt,zt,bxt,byt,bzt)
      bx_flt = bxt
      by_flt = byt
      bz_flt = bzt
      RETURN
      ! END SUBROUTINE
      END SUBROUTINE bfield_virtual_casing_flt
      !-----------------------------------------------------------------
         
         
      !-----------------------------------------------------------------
      SUBROUTINE bfield_virtual_casing_adapt_dbl(x,y,z,bx,by,bz,istat)
      IMPLICIT NONE
      ! INPUT VARIABLES
      DOUBLE PRECISION, INTENT(in)  :: x, y, z
      DOUBLE PRECISION, INTENT(out) :: bx, by, bz
      INTEGER, INTENT(inout) :: istat
      ! LOCAL VARIABLES
      LOGICAL            :: adapt_rerun
      INTEGER(KIND=8), PARAMETER :: ndim_nag = 2 ! theta,zeta
      INTEGER(KIND=8), PARAMETER :: nfun_nag = 3 ! Bx, By, Bz
      !INTEGER, PARAMETER :: lenwrk_nag = 6*ndim_nag+9*nfun_nag+(ndim_nag+nfun_nag+2)*(1+1)
      INTEGER(KIND=8), PARAMETER :: lenwrk_nag = IWRK
      INTEGER(KIND=8) :: maxcls_nag,mincls_nag, subs, restar, wrklen, rulcls, wrksbs, n, m, funcls
      !DOUBLE PRECISION :: gf(1:nuvp), gf3(1:nuvp)
      DOUBLE PRECISION :: absreq_nag, relreq_nag
      DOUBLE PRECISION :: wrkstr_nag(lenwrk_nag)
      DOUBLE PRECISION :: a_nag(ndim_nag), b_nag(ndim_nag), &
                          finest_nag(nfun_nag), absest_nag(nfun_nag)
      DOUBLE PRECISION, ALLOCATABLE :: vrtwrk(:)

#ifdef NAG
      EXTERNAL :: D01EAF

#else
      EXTERNAL :: dcuhre

#endif

      ! BEGIN SUBROUTINE
      IF (adapt_tol < 0) THEN
         CALL bfield_virtual_casing_dbl(x,y,z,bx,by,bz)
         RETURN
      END IF

      a_nag(1:2) = zero
      b_nag(1:2) = one
      mincls_nag = MIN_CLS
      maxcls_nag = IWRK
      !absreq_nag = zero
      absreq_nag = adapt_tol       ! Talk to Stuart about proper values
      relreq_nag = adapt_rel ! Talk to Stuart about proper values
      finest_nag = zero
      absest_nag = zero
      x_nag      = x
      y_nag      = y
      z_nag      = z
      adapt_rerun = .true.
      subs = 1
      restar = 0
      DO WHILE (adapt_rerun) 

#ifdef NAG
         CALL D01EAF(ndim_nag,a_nag,b_nag,mincls_nag,maxcls_nag,nfun_nag,funsub_nag_b,absreq_nag,&
                   relreq_nag,lenwrk_nag,wrkstr_nag,finest_nag,absest_nag,istat)
         IF (istat == 1 .or. istat == 3) THEN
            maxcls_nag = maxcls_nag*10
            mincls_nag = -1
            restar = 1
            WRITE(6,*) '!!!!!  WARNING Could not reach desired tollerance  !!!!!'
            WRITE(6,*) '  BX = ',finest_nag(1),' +/- ',absest_nag(1)
            WRITE(6,*) '  BY = ',finest_nag(2),' +/- ',absest_nag(2)
            WRITE(6,*) '  BZ = ',finest_nag(3),' +/- ',absest_nag(3)
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         ELSE IF (mincls_nag <= 32) THEN
            mincls_nag = 64
            restar = 1
         ELSE IF (istat < 0) THEN
            bx = zero
            by = zero
            bz = zero
            adapt_rerun=.false.
         ELSE
            bx = finest_nag(1)
            by = finest_nag(2)
            bz = finest_nag(3)
            adapt_rerun=.false.
         END IF

#else
         IF (.not.ALLOCATED(vrtwrk)) THEN
            wrklen = ((maxcls_nag-ndim_nag)/(2*ndim_nag) + 1)*(2*ndim_nag+2*nfun_nag+2) + 17*nfun_nag + 256
            ALLOCATE(vrtwrk(wrklen),STAT=istat)
            IF (istat .ne. 0) THEN
               WRITE(6,*) ' ALLOCATION ERROR IN: bfield_virtual_casing_adapt_dbl'
               WRITE(6,*) '   VARIABLE: VRTWRK, SIZE: ',wrklen
               WRITE(6,*) '   ISTAT: ',istat
               RETURN
            END IF
         END IF
         CALL dcuhre(ndim_nag,nfun_nag,a_nag,b_nag,mincls_nag,maxcls_nag,funsub_nag_b,absreq_nag,&
                     relreq_nag,0,wrklen,restar,finest_nag,absest_nag,funcls,istat,vrtwrk)
         !PRINT *,istat,mincls_nag,funcls,maxcls_nag,finest_nag(1),finest_nag(2),finest_nag(3)
         IF (istat == 1) THEN
            maxcls_nag = maxcls_nag*10
            mincls_nag = funcls
            restar = 1
            istat = 0
         ELSE IF (istat > 1) THEN
            bx = zero
            by = zero
            bz = zero
            adapt_rerun=.false.
            DEALLOCATE(vrtwrk)
         ELSE
            bx = finest_nag(1)
            by = finest_nag(2)
            bz = finest_nag(3)
            adapt_rerun=.false.
            DEALLOCATE(vrtwrk)
         END IF

#endif
      END DO
      nlastcall=mincls_nag
      RETURN
      ! END SUBROUTINE
      END SUBROUTINE bfield_virtual_casing_adapt_dbl
      !-----------------------------------------------------------------
         
         
      !-----------------------------------------------------------------
      SUBROUTINE bfield_virtual_casing_adapt_flt(x_flt,y_flt,z_flt,bx_flt,by_flt,bz_flt,istat)
      IMPLICIT NONE
      ! INPUT VARIABLES
      REAL, INTENT(in)  :: x_flt, y_flt, z_flt
      REAL, INTENT(out) :: bx_flt, by_flt, bz_flt
      INTEGER, INTENT(inout) :: istat
      ! LOCAL VARIABLES
      DOUBLE PRECISION :: xt,yt,zt,bxt,byt,bzt
      ! BEGIN SUBROUTINE
      xt = x_flt
      yt = y_flt
      zt = z_flt
      bxt = zero
      byt = zero
      bzt = zero
      CALL bfield_virtual_casing_adapt_dbl(xt,yt,zt,bxt,byt,bzt,istat)
      bx_flt = bxt
      by_flt = byt
      bz_flt = bzt
      RETURN
      ! END SUBROUTINE
      END SUBROUTINE bfield_virtual_casing_adapt_flt
      !-----------------------------------------------------------------
         
         
      !-----------------------------------------------------------------
      SUBROUTINE funsub_nag_b(ndim, vec, nfun, f)
      IMPLICIT NONE
      ! INPUT VARIABLES
      INTEGER, INTENT(in) :: ndim, nfun
      DOUBLE PRECISION, INTENT(in) :: vec(ndim)
      DOUBLE PRECISION, INTENT(out) :: f(nfun)
      ! LOCAL VARIABLES
      INTEGER :: ier
      DOUBLE PRECISION :: bn, bx, by, bz , xs, ys, zs, gf, gf3, nx, ny, nz, ax, ay, az
      ! BEGIN SUBROUTINE

      xs = zero; ys = zero; zs = zero
      CALL EZspline_interp(x_spl,vec(1),vec(2),xs,ier)
      CALL EZspline_interp(y_spl,vec(1),vec(2),ys,ier)
      CALL EZspline_interp(z_spl,vec(1),vec(2),zs,ier)
      CALL EZspline_interp(kx_spl,vec(1),vec(2),ax,ier)
      CALL EZspline_interp(ky_spl,vec(1),vec(2),ay,ier)
      CALL EZspline_interp(kz_spl,vec(1),vec(2),az,ier)

      !CALL EZspline_interp(bn_spl,vec(1),vec(2),bn,ier)
      bn = zero
      gf   = one/DSQRT((x_nag-xs)*(x_nag-xs)+(y_nag-ys)*(y_nag-ys)+(z_nag-zs)*(z_nag-zs))
      gf3  = gf*gf*gf
      f(1) = norm_nag*(ay*(z_nag-zs)-az*(y_nag-ys)+bn*(x_nag-xs))*gf3
      f(2) = norm_nag*(az*(x_nag-xs)-ax*(z_nag-zs)+bn*(y_nag-ys))*gf3
      f(3) = norm_nag*(ax*(y_nag-ys)-ay*(x_nag-xs)+bn*(z_nag-zs))*gf3
      !WRITE(427,*) vec(1),vec(2),xs,ys,zs
      RETURN
      ! END SUBROUTINE
      END SUBROUTINE funsub_nag_b
      !-----------------------------------------------------------------
         
         
      !-----------------------------------------------------------------
      SUBROUTINE vecpot_virtual_casing_dbl(x,y,z,ax,ay,az)
      IMPLICIT NONE
      ! INPUT VARIABLES
      DOUBLE PRECISION, INTENT(in)  :: x, y, z
      DOUBLE PRECISION, INTENT(out) :: ax, ay,az
      ! LOCAL VARIABLES
      DOUBLE PRECISION ::  gf(1:nuvp)
      ! BEGIN SUBROUTINE
      gf(1:nuvp) = one/dsqrt(  (x-xsurf(1:nuvp))**2 &
                                 + (y-ysurf(1:nuvp))**2 &
                                 + (z-zsurf(1:nuvp))**2)
      ax = norm*sum((btopx(1:nuvp)+btops(1:nuvp)*nxsurf(1:nuvp))*gf(1:nuvp))
      ay = norm*sum((btopy(1:nuvp)+btops(1:nuvp)*nysurf(1:nuvp))*gf(1:nuvp))
      az = norm*sum((btopz(1:nuvp)+btops(1:nuvp)*nzsurf(1:nuvp))*gf(1:nuvp))
      ! END SUBROUTINE
      END SUBROUTINE vecpot_virtual_casing_dbl
      !-----------------------------------------------------------------
         
         
      !-----------------------------------------------------------------
      SUBROUTINE vecpot_virtual_casing_flt(x_flt,y_flt,z_flt,ax_flt,ay_flt,az_flt)
      IMPLICIT NONE
      ! INPUT VARIABLES
      REAL, INTENT(in)  :: x_flt, y_flt, z_flt
      REAL, INTENT(out) :: ax_flt, ay_flt,az_flt
      ! LOCAL VARIABLES
      DOUBLE PRECISION  :: xt, yt, zt
      DOUBLE PRECISION  :: axt, ayt,azt
      ! BEGIN SUBROUTINE
      xt   = x_flt
      yt   = y_flt
      zt   = z_flt
      axt = zero
      ayt = zero
      azt = zero
      CALL vecpot_virtual_casing_dbl(xt,yt,zt,axt,ayt,azt)
      ax_flt = axt
      ay_flt = ayt
      az_flt = azt
      RETURN
      ! END SUBROUTINE
      END SUBROUTINE vecpot_virtual_casing_flt
      !-----------------------------------------------------------------
         
         
      !-----------------------------------------------------------------
      SUBROUTINE vecpot_virtual_casing_adapt_dbl(x,y,z,ax,ay,az,istat)
      IMPLICIT NONE
      ! INPUT VARIABLES
      DOUBLE PRECISION, INTENT(in)  :: x, y, z
      DOUBLE PRECISION, INTENT(out) :: ax, ay, az
      INTEGER, INTENT(inout) :: istat
      ! LOCAL VARIABLES
      LOGICAL            :: adapt_rerun
      INTEGER(KIND=8), PARAMETER :: ndim_nag = 2 ! theta,zeta
      INTEGER(KIND=8), PARAMETER :: nfun_nag = 3 ! Bx, By, Bz
      !INTEGER, PARAMETER :: lenwrk_nag = 6*ndim_nag+9*nfun_nag+(ndim_nag+nfun_nag+2)*(1+1)
      INTEGER(KIND=8), PARAMETER :: lenwrk_nag = IWRK
      INTEGER(KIND=8) :: maxcls_nag,mincls_nag, subs, restar, wrklen, rulcls, wrksbs, n, m, funcls
      !DOUBLE PRECISION ::  gf(1:nuvp), gf3(1:nuvp)
      DOUBLE PRECISION :: absreq_nag, relreq_nag
      DOUBLE PRECISION :: wrkstr_nag(lenwrk_nag)
      DOUBLE PRECISION :: a_nag(ndim_nag), b_nag(ndim_nag), &
                          finest_nag(nfun_nag), absest_nag(nfun_nag)
      DOUBLE PRECISION, ALLOCATABLE :: vrtwrk(:)

#ifdef NAG
      EXTERNAL :: D01EAF

#else
      EXTERNAL :: dcuhre

#endif
      ! BEGIN SUBROUTINE 
      IF (adapt_tol < 0) THEN
         CALL vecpot_virtual_casing_dbl(x,y,z,ax,ay,az)
         RETURN
      END IF

      a_nag(1:2) = zero
      b_nag(1:2) = one
      mincls_nag = MIN_CLS
      maxcls_nag = IWRK
      absreq_nag = adapt_tol       ! Talk to Stuart about proper values
      relreq_nag = adapt_rel ! Talk to Stuart about proper values
      finest_nag = zero
      absest_nag = zero
      x_nag      = x
      y_nag      = y
      z_nag      = z
      adapt_rerun = .true.
      subs = 1
      restar = 0
      DO WHILE (adapt_rerun)

#ifdef NAG
         CALL D01EAF(ndim_nag,a_nag,b_nag,mincls_nag,maxcls_nag,nfun_nag,funsub_nag_a,absreq_nag,&
                   relreq_nag,lenwrk_nag,wrkstr_nag,finest_nag,absest_nag,istat)
         IF (istat == 1 .or. istat == 3) THEN
            maxcls_nag = maxcls_nag*10
            mincls_nag = -1
            restar = 1
            WRITE(6,*) '!!!!!  WARNING Could not reach desired tollerance  !!!!!'
            WRITE(6,*) '  AX = ',finest_nag(1),' +/- ',absest_nag(1)
            WRITE(6,*) '  AY = ',finest_nag(2),' +/- ',absest_nag(2)
            WRITE(6,*) '  AZ = ',finest_nag(3),' +/- ',absest_nag(3)
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         ELSE IF (mincls_nag <= 32) THEN
            mincls_nag = 64
            restar = 1
         ELSE IF (istat < 0) THEN
            ax = zero
            ay = zero
            az = zero
            adapt_rerun=.false.
            print *,'error',istat
         ELSE
            ax = finest_nag(1)
            ay = finest_nag(2)
            az = finest_nag(3)
            adapt_rerun=.false.
         END IF

#else
         IF (.not.ALLOCATED(vrtwrk)) THEN
            wrklen = ((IWRK-ndim_nag)/(2*ndim_nag) + 1)*(2*ndim_nag+2*nfun_nag+2) + 17*nfun_nag + 256
            ALLOCATE(vrtwrk(wrklen),STAT=istat)
            IF (istat .ne. 0) THEN
               WRITE(6,*) ' ALLOCATION ERROR IN: vecpot_virtual_casing_adapt_dbl'
               WRITE(6,*) '   VARIABLE: VRTWRK, SIZE: ',wrklen
               WRITE(6,*) '   ISTAT: ',istat
               RETURN
            END IF
         END IF
         CALL dcuhre(ndim_nag,nfun_nag,a_nag,b_nag,mincls_nag,maxcls_nag,funsub_nag_a,absreq_nag,&
                     relreq_nag,0,wrklen,restar,finest_nag,absest_nag,funcls,istat,vrtwrk)
         !DEALLOCATE(vrtwrk)
         IF (istat == 1) THEN
            maxcls_nag = maxcls_nag*10
            mincls_nag = -1
            restar = 1
         ELSE IF (istat > 1) THEN
            ax = zero
            ay = zero
            az = zero
            adapt_rerun=.false.
            DEALLOCATE(vrtwrk)
         ELSE
            ax = finest_nag(1)
            ay = finest_nag(2)
            az = finest_nag(3)
            adapt_rerun=.false.
            DEALLOCATE(vrtwrk)
         END IF

#endif
      END DO
      nlastcall=mincls_nag
      RETURN
      ! END SUBROUTINE
      END SUBROUTINE vecpot_virtual_casing_adapt_dbl
      !-----------------------------------------------------------------
         
         
      !-----------------------------------------------------------------
      SUBROUTINE vecpot_virtual_casing_adapt_flt(x_flt,y_flt,z_flt,ax_flt,ay_flt,az_flt,istat)
      IMPLICIT NONE
      ! INPUT VARIABLES
      REAL, INTENT(in)  :: x_flt, y_flt, z_flt
      REAL, INTENT(out) :: ax_flt, ay_flt,az_flt
      INTEGER, INTENT(inout) :: istat
      ! LOCAL VARIABLES
      DOUBLE PRECISION  :: xt, yt, zt
      DOUBLE PRECISION  :: axt, ayt,azt
      ! BEGIN SUBROUTINE 
      xt   = x_flt
      yt   = y_flt
      zt   = z_flt
      axt = zero
      ayt = zero
      azt = zero
      CALL vecpot_virtual_casing_adapt_dbl(xt,yt,zt,axt,ayt,azt,istat)
      ax_flt = axt
      ay_flt = ayt
      az_flt = azt
      RETURN
      ! END SUBROUTINE
      END SUBROUTINE vecpot_virtual_casing_adapt_flt
      !-----------------------------------------------------------------
         
         
      !-----------------------------------------------------------------
      SUBROUTINE funsub_nag_a(ndim, vec, nfun, f)
      IMPLICIT NONE
      ! INPUT VARIABLES
      INTEGER, INTENT(in) :: ndim, nfun
      DOUBLE PRECISION, INTENT(in) :: vec(ndim)
      DOUBLE PRECISION, INTENT(out) :: f(nfun)
      ! LOCAL VARIABLES
      INTEGER :: ier
      DOUBLE PRECISION :: bn, bx, by, bz , xs, ys, zs, gf, nx, ny, nz, ax, ay, az
      ! BEGIN SUBROUTINE

      xs = zero; ys = zero; zs = zero; ax = zero; ay = zero; az = zero
      nx = zero; ny = zero; nz = zero; bn = zero
      CALL EZspline_interp(x_spl,vec(1),vec(2),xs,ier)
      CALL EZspline_interp(y_spl,vec(1),vec(2),ys,ier)
      CALL EZspline_interp(z_spl,vec(1),vec(2),zs,ier)
      CALL EZspline_interp(kx_spl,vec(1),vec(2),ax,ier)
      CALL EZspline_interp(ky_spl,vec(1),vec(2),ay,ier)
      CALL EZspline_interp(kz_spl,vec(1),vec(2),az,ier)
      CALL EZspline_interp(nx_spl,vec(1),vec(2),nx,ier)
      CALL EZspline_interp(ny_spl,vec(1),vec(2),ny,ier)
      CALL EZspline_interp(nz_spl,vec(1),vec(2),nz,ier)
      CALL EZspline_interp(bn_spl,vec(1),vec(2),bn,ier)

      gf   = one/DSQRT((x_nag-xs)*(x_nag-xs)+(y_nag-ys)*(y_nag-ys)+(z_nag-zs)*(z_nag-zs))
      f(1) = norm_nag*(ax+bn*nx)*gf
      f(2) = norm_nag*(ay+bn*ny)*gf
      f(3) = norm_nag*(az+bn*nz)*gf    
      RETURN
      ! END SUBROUTINE
      END SUBROUTINE funsub_nag_a
      !-----------------------------------------------------------------
         
         
      !-----------------------------------------------------------------
      SUBROUTINE virtual_casing_info(iunit)
      IMPLICIT NONE
      ! INPUT VARIABLES
      INTEGER, INTENT(in)  :: iunit
      ! LOCAL VARIABLES
      ! BEGIN SUBROUTINE
      WRITE(iunit,'(A)')                     '----- Virtual Casing Information -----'
#ifdef NAG 
      WRITE(iunit,'(A,A,A)')                   '   INTEGRAL TYPE: ',TRIM(vc_type_str),' (NAG) '
#else
      WRITE(iunit,'(A,A,A)')                   '   INTEGRAL TYPE: ',TRIM(vc_type_str),' (DCUHRE) '
#endif
      WRITE(iunit,'(A,ES11.4)')              '   MIN_GRID_DISTANCE = ',min_delta_x
      WRITE(iunit,'(A,ES11.4)')              '   NORMAL_AREA = ',surf_area
      WRITE(iunit,'(A,I4,A,I4,A,I4,A,I3)')   '   NR = ',nr_vc,';   NU = ',nu_vc,';  NV = ',nv_vc,';  NFP = ',nvp/nv_vc
      WRITE(iunit,'(A,I6)')                  '   NUVP = ',nuvp
      IF (adapt_tol > 0 .or. adapt_rel >0) THEN
         WRITE(iunit,'(A,ES11.4,A,ES11.4)')  '   ABS_TOL = ',adapt_tol,';   REL_TOL = ',adapt_rel
      ELSE
         WRITE(iunit,'(A)')                  '   !!!!! Using discrete integration !!!!!'
      END IF
      WRITE(iunit,'(A,I6,A,I8,A)')                  '   MIN_CLS = ',MIN_CLS,'   (',IWRK,')'
      CALL FLUSH(iunit)
      ! END SUBROUTINE
      END SUBROUTINE virtual_casing_info
      !-----------------------------------------------------------------
         
         
      !-----------------------------------------------------------------
      SUBROUTINE virtual_casing_surf_dump(iunit)
      IMPLICIT NONE
      ! INPUT VARIABLES
      INTEGER, INTENT(in)  :: iunit
      ! LOCAL VARIABLES
      INTEGER :: i, j, k, ik
      ! BEGIN SUBROUTINE
      ik = 1
      WRITE(iunit,'(A)') ' n  nfp  nu  nv  x  y  z  nx  ny  nz  btopx  btopy  btopz  btops'
      DO i = 1, nuvp/(nu_vc*nv_vc)
         DO k = 1, nu_vc
            DO j = 1, nv_vc
               WRITE(iunit,'(4I6,10ES23.16)') ik, i, k, j, &
                                 xsurf(ik), ysurf(ik), zsurf(ik), &
                                 nxsurf(ik), nysurf(ik), nzsurf(ik), &
                                 btopx(ik), btopy(ik), btopz(ik), &
                                 btops(ik)
               ik = ik + 1
            END DO
         END DO
      END DO
      CALL FLUSH(iunit)
      ! END SUBROUTINE
      END SUBROUTINE virtual_casing_surf_dump
      !-----------------------------------------------------------------
         
         
      !-----------------------------------------------------------------
      DOUBLE PRECISION FUNCTION virtual_casing_dist_dbl(xp,yp,zp)
      IMPLICIT NONE
      ! INPUT VARIABLES
      DOUBLE PRECISION, INTENT(in)  :: xp,yp,zp
      ! LOCAL VARIABLES
      ! BEGIN FUNCTION
      virtual_casing_dist_dbl = MINVAL(DSQRT((xp-xsurf)**2+(yp-ysurf)**2+(zp-zsurf)**2))
      ! END FUNCTION
      END FUNCTION virtual_casing_dist_dbl
      !-----------------------------------------------------------------
         
         
      !-----------------------------------------------------------------
      REAL FUNCTION virtual_casing_dist_flt(xp_flt,yp_flt,zp_flt)
      IMPLICIT NONE
      ! INPUT VARIABLES
      REAL, INTENT(in)  :: xp_flt,yp_flt,zp_flt
      ! LOCAL VARIABLES
      ! BEGIN FUNCTION
      virtual_casing_dist_flt = MINVAL(DSQRT((xp_flt-xsurf)**2+(yp_flt-ysurf)**2+(zp_flt-zsurf)**2))
      ! END FUNCTION
      END FUNCTION virtual_casing_dist_flt
      !-----------------------------------------------------------------
         
         
      !-----------------------------------------------------------------
      SUBROUTINE init_volint_dbl(mnmax,nu,nv,ns,xm,xn,rmnc,zmns,nfp,jumnc,jvmnc,&
                                     rmns,zmnc,jumns,jvmns)
      IMPLICIT NONE
      ! INPUT VARIABLES
      INTEGER, INTENT(in) :: mnmax, nu, nv, nfp,ns
      INTEGER, INTENT(in) :: xm(1:mnmax),xn(1:mnmax)
      DOUBLE PRECISION, INTENT(in) :: rmnc(1:mnmax,ns),zmns(1:mnmax,ns)
      DOUBLE PRECISION, INTENT(in) :: jumnc(1:mnmax,ns),jvmnc(1:mnmax,ns)
      DOUBLE PRECISION, INTENT(in), OPTIONAL :: jumns(1:mnmax,ns),jvmns(1:mnmax,ns)
      DOUBLE PRECISION, INTENT(in), OPTIONAL :: rmns(1:mnmax,ns),zmnc(1:mnmax,ns)
      ! LOCAL VARIABLES
      INTEGER :: mn, uv, i, u, v, dex1, dex2, ier, nuv, k
      INTEGER :: bcs1(2), bcs2(2), bcs3(2)
      DOUBLE PRECISION :: cop, sip
      DOUBLE PRECISION, ALLOCATABLE :: xu(:), xv(:)
      DOUBLE PRECISION, ALLOCATABLE :: fmn_temp(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: r_temp(:,:,:), z_temp(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: x_temp(:,:,:), y_temp(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: jv_temp(:,:,:), ju_temp(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: ru(:,:,:), rv(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: zu(:,:,:), zv(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: jr(:,:,:), jphi(:,:,:), jz(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: jx(:,:,:), jy(:,:,:)
      ! BEGIN SUBROUTINE
      ! Deallocate anything globally allocated

      IF (EZspline_allocated(x3d_spl)) CALL EZspline_free(x3d_spl,ier)
      IF (EZspline_allocated(y3d_spl)) CALL EZspline_free(y3d_spl,ier)
      IF (EZspline_allocated(z3d_spl)) CALL EZspline_free(z3d_spl,ier)
      IF (EZspline_allocated(jx3d_spl)) CALL EZspline_free(jx3d_spl,ier)
      IF (EZspline_allocated(jy3d_spl)) CALL EZspline_free(jy3d_spl,ier)
      IF (EZspline_allocated(jz3d_spl)) CALL EZspline_free(jz3d_spl,ier)

      ! Initialize varaibles
      bcs1=(/-1,-1/)
      bcs2=(/-1,-1/)
      bcs3=(/0,0/)
      adapt_tol = 1.0D-5
      adapt_rel = 1.0D-3
      nr_vc = ns
      nu_vc = nu
      nv_vc = nv
      nuv = nu * nv
      nvp = nv * nfp
      nuvp = nu * nv * nfp
      norm_3d = 1.0D-7 ! Worked for integration 0 to 2pi
      norm_3d = pi2*pi2*1.0D-7
      vc_type_str='Volume Integral'
      ! Alloocations
      ALLOCATE(xu(nu),xv(nvp))
      ALLOCATE(fmn_temp(mnmax,ns))
      ALLOCATE(r_temp(nu,nvp,ns),z_temp(nu,nvp,ns))
      ALLOCATE(x_temp(nu,nvp,ns),y_temp(nu,nvp,ns))
      ALLOCATE(ju_temp(nu,nvp,ns),jv_temp(nu,nvp,ns))
      ALLOCATE(ru(nu,nvp,ns),rv(nu,nvp,ns))
      ALLOCATE(zu(nu,nvp,ns),zv(nu,nvp,ns))
      ALLOCATE(jx(nu,nvp,ns),jy(nu,nvp,ns))
      ALLOCATE(jr(nu,nvp,ns),jphi(nu,nvp,ns),jz(nu,nvp,ns))
      ALLOCATE(jx_3d(nu*nvp*ns),jy_3d(nu*nvp*ns),jz_3d(nu*nvp*ns))
      ALLOCATE(xsurf(nu*nvp*ns),ysurf(nu*nvp*ns),zsurf(nu*nvp*ns))
      
      ! Get the major quantities (not require use the VMEC convenction of j of j=j*jac
      r_temp=zero; z_temp=zero; ju_temp=zero; jv_temp=zero
      FORALL(u=1:nu) xu(u) = DBLE(u-1)/DBLE(nu-1)
      FORALL(v=1:nvp) xv(v) = DBLE(v-1)/DBLE(nvp-1)
      CALL mntouv_local(1,ns,mnmax,nu,nvp,xu,xv,rmnc,xm,xn,r_temp,0,1)
      CALL mntouv_local(1,ns,mnmax,nu,nvp,xu,xv,zmns,xm,xn,z_temp,1,0)
      CALL mntouv_local(1,ns,mnmax,nu,nvp,xu,xv,jumnc,xm,xn,ju_temp,0,0)
      CALL mntouv_local(1,ns,mnmax,nu,nvp,xu,xv,jvmnc,xm,xn,jv_temp,0,0)
      IF (PRESENT(rmns)) CALL mntouv_local(1,ns,mnmax,nu,nvp,xu,xv,rmns,xm,xn,r_temp,1,0)
      IF (PRESENT(zmnc)) CALL mntouv_local(1,ns,mnmax,nu,nvp,xu,xv,zmnc,xm,xn,z_temp,0,0)
      IF (PRESENT(jumns)) CALL mntouv_local(1,ns,mnmax,nu,nvp,xu,xv,jumns,xm,xn,ju_temp,1,0)
      IF (PRESENT(jvmns)) CALL mntouv_local(1,ns,mnmax,nu,nvp,xu,xv,jvmns,xm,xn,jv_temp,1,0)
      ! Now we calculate the edge metric elements
      ru = zero; zu = zero; rv = zero; zv = zero
      FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -rmnc(mn,:)*xm(mn)
      CALL mntouv_local(1,ns,mnmax,nu,nvp,xu,xv,fmn_temp,xm,xn,ru,1,0)
      FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -rmnc(mn,:)*xn(mn)
      CALL mntouv_local(1,ns,mnmax,nu,nvp,xu,xv,fmn_temp,xm,xn,rv,1,0)  
      FORALL(mn = 1:mnmax) fmn_temp(mn,:) = zmns(mn,:)*xm(mn)
      CALL mntouv_local(1,ns,mnmax,nu,nvp,xu,xv,fmn_temp,xm,xn,zu,0,0) 
      FORALL(mn = 1:mnmax) fmn_temp(mn,:) = zmns(mn,:)*xn(mn)
      CALL mntouv_local(1,ns,mnmax,nu,nvp,xu,xv,fmn_temp,xm,xn,zv,0,0)  
      IF (PRESENT(rmns)) THEN
         FORALL(mn = 1:mnmax) fmn_temp(mn,:) = rmns(mn,:)*xm(mn)
         CALL mntouv_local(1,ns,mnmax,nu,nvp,xu,xv,fmn_temp,xm,xn,ru,0,0)
         FORALL(mn = 1:mnmax) fmn_temp(mn,:) = rmns(mn,:)*xn(mn)
         CALL mntouv_local(1,ns,mnmax,nu,nvp,xu,xv,fmn_temp,xm,xn,rv,0,0)  
      END IF 
      IF (PRESENT(zmnc)) THEN
         FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -zmnc(mn,:)*xm(mn)
         CALL mntouv_local(1,ns,mnmax,nu,nvp,xu,xv,fmn_temp,xm,xn,zu,1,0)  
         FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -zmnc(mn,:)*xn(mn)
         CALL mntouv_local(1,ns,mnmax,nu,nvp,xu,xv,fmn_temp,xm,xn,zv,1,0)
      END IF
      ! Calculate jr, jphi, jz, jx, jy
      jr   = ju_temp*ru+jv_temp*rv*nfp
      jphi = r_temp * jv_temp
      jz   = ju_temp*zu+jv_temp*zv*nfp
      DO u = 1, nu
         DO v = 1, nvp
            cop = DCOS(pi2*xv(v))
            sip = DSIN(pi2*xv(v))
            x_temp(u,v,:) = r_temp(u,v,:) * cop
            y_temp(u,v,:) = r_temp(u,v,:) * sip
            jx(u,v,:) = jr(u,v,:) * cop - jphi(u,v,:) * sip
            jy(u,v,:) = jr(u,v,:) * sip + jphi(u,v,:) * cop
         END DO
      END DO
      i=1
      DO u = 1, nu
         DO v = 1, nvp
            DO k = 1, ns
               jx_3d(i) = jx(u,v,k)
               jy_3d(i) = jy(u,v,k)
               jz_3d(i) = jz(u,v,k)
               xsurf(i) = x_temp(u,v,k)
               ysurf(i) = y_temp(u,v,k)
               zsurf(i) = z_temp(u,v,k)
               i = i + 1
            END DO
         END DO
      END DO

      ! Construct Splines
      CALL EZspline_init(x3d_spl,nu,nvp,ns,bcs1,bcs2,bcs3,ier)
      CALL EZspline_init(y3d_spl,nu,nvp,ns,bcs1,bcs2,bcs3,ier)
      CALL EZspline_init(z3d_spl,nu,nvp,ns,bcs1,bcs2,bcs3,ier)
      CALL EZspline_init(jx3d_spl,nu,nvp,ns,bcs1,bcs2,bcs3,ier)
      CALL EZspline_init(jy3d_spl,nu,nvp,ns,bcs1,bcs2,bcs3,ier)
      CALL EZspline_init(jz3d_spl,nu,nvp,ns,bcs1,bcs2,bcs3,ier)
      x3d_spl%x1=xu
      y3d_spl%x1=xu
      z3d_spl%x1=xu
      jx3d_spl%x1=xu
      jy3d_spl%x1=xu
      jz3d_spl%x1=xu
      x3d_spl%x2=xv
      y3d_spl%x2=xv
      z3d_spl%x2=xv
      jx3d_spl%x2=xv
      jy3d_spl%x2=xv
      jz3d_spl%x2=xv
      ! default [0 1] x3
      x3d_spl%isHermite = 1
      y3d_spl%isHermite = 1
      z3d_spl%isHermite = 1
      jx3d_spl%isHermite = 1
      jy3d_spl%isHermite = 1
      jz3d_spl%isHermite = 1
      CALL EZspline_setup(x3d_spl,x_temp,ier)
      CALL EZspline_setup(y3d_spl,y_temp,ier)
      CALL EZspline_setup(z3d_spl,z_temp,ier)
      CALL EZspline_setup(jx3d_spl,jx,ier)
      CALL EZspline_setup(jy3d_spl,jy,ier)
      CALL EZspline_setup(jz3d_spl,jz,ier)

      ! DEALLOCATE
      DEALLOCATE(xu,xv)
      DEALLOCATE(fmn_temp)
      DEALLOCATE(r_temp,z_temp)
      DEALLOCATE(x_temp,y_temp)
      DEALLOCATE(ju_temp,jv_temp)
      DEALLOCATE(ru,rv)
      DEALLOCATE(zu,zv)
      DEALLOCATE(jx,jy,jr,jphi,jz)
      ! END SUBROUTINE
      END SUBROUTINE init_volint_dbl
                              
      !-----------------------------------------------------------------      
      
      !-----------------------------------------------------------------
      !  Note optional arguments must have a different name so the
      !  Interface will pick the proper variables.  Here:  _FLT
      SUBROUTINE init_volint_flt(mnmax_flt,nu_flt,nv_flt,ns_flt,xm_flt,&
                                         xn_flt,rmnc_flt,zmns_flt,nfp_flt,&
                                         jumnc_flt,jvmnc_flt,&
                                         rmns_flt,zmnc_flt,&
                                         jumns_flt,jvmns_flt)
      IMPLICIT NONE
      ! INPUT VARIABLES
      INTEGER, INTENT(in) :: mnmax_flt, nu_flt, nv_flt, nfp_flt, ns_flt
      INTEGER, INTENT(in) :: xm_flt(1:mnmax_flt),xn_flt(1:mnmax_flt)
      REAL, INTENT(in) :: rmnc_flt(1:mnmax_flt,ns_flt),zmns_flt(1:mnmax_flt,ns_flt)
      REAL, INTENT(in) :: jumnc_flt(1:mnmax_flt,ns_flt),jvmnc_flt(1:mnmax_flt,ns_flt)
      REAL, INTENT(in), OPTIONAL :: jumns_flt(1:mnmax_flt,ns_flt),jvmns_flt(1:mnmax_flt,ns_flt)
      REAL, INTENT(in), OPTIONAL :: rmns_flt(1:mnmax_flt,ns_flt),zmnc_flt(1:mnmax_flt,ns_flt)
      ! LOCAL VARIABLES
      DOUBLE PRECISION :: rmnct(1:mnmax_flt,ns_flt),zmnst(1:mnmax_flt,ns_flt)
      DOUBLE PRECISION :: rmnst(1:mnmax_flt,ns_flt),zmnct(1:mnmax_flt,ns_flt)
      DOUBLE PRECISION :: jumnct(1:mnmax_flt,ns_flt),jvmnct(1:mnmax_flt,ns_flt)
      DOUBLE PRECISION :: jumnst(1:mnmax_flt,ns_flt),jvmnst(1:mnmax_flt,ns_flt)
      ! BEGIN SUBROUTINE
      rmnct = zero; rmnst = zero; zmnct = zero; zmnst = zero
      jumnct = zero; jumnst = zero
      jvmnct = zero; jvmnst = zero
      rmnct  = rmnc_flt
      zmnst  = zmns_flt
      jumnct = jumnc_flt
      jvmnct = jvmnc_flt
      IF (PRESENT(rmns_flt))  rmnst = rmns_flt
      IF (PRESENT(zmnc_flt))  zmnct = zmnc_flt
      IF (PRESENT(jumns_flt)) jumnst = jumns_flt
      IF (PRESENT(jvmns_flt)) jvmnst = jvmns_flt
      CALL init_volint_dbl(mnmax_flt,nu_flt,nv_flt,ns_flt,xm_flt,xn_flt, &
                                   rmnct,zmnst,nfp_flt, &
                                   RMNS=rmnst,ZMNC=zmnct, &
                                   JUMNC=jumnct, JUMNS=jumnst, &
                                   JVMNC=jvmnct, JVMNS=jvmnst)
      !CALL init_volint_dbl(mnmax_flt,nu_flt,nv_flt,ns_flt,xm_flt,xn_flt, &
      !                             rmnct,zmnst,nfp_flt, &
      !                             JUMNC=jumnct, JVMNC=jvmnct)
      ! END SUBROUTINE
      END SUBROUTINE init_volint_flt
      !----------------------------------------------------------------- 
         
         
      !-----------------------------------------------------------------
      SUBROUTINE funsub_nag_a3d(ndim, vec, nfun, f)
      IMPLICIT NONE
      ! INPUT VARIABLES
      INTEGER, INTENT(in) :: ndim, nfun
      DOUBLE PRECISION, INTENT(in) :: vec(ndim)
      DOUBLE PRECISION, INTENT(out) :: f(nfun)
      ! LOCAL VARIABLES
      INTEGER :: ier
      DOUBLE PRECISION :: ax, ay, az , xs, ys, zs, gf, nx, ny, nz
      ! BEGIN SUBROUTINE

      xs = zero; ys = zero; zs = zero
      ax = zero; ay = zero; az = zero
      CALL EZspline_interp(x3d_spl,vec(1),vec(2),vec(3),xs,ier)
      CALL EZspline_interp(y3d_spl,vec(1),vec(2),vec(3),ys,ier)
      CALL EZspline_interp(z3d_spl,vec(1),vec(2),vec(3),zs,ier)
      CALL EZspline_interp(jx3d_spl,vec(1),vec(2),vec(3),ax,ier)
      CALL EZspline_interp(jy3d_spl,vec(1),vec(2),vec(3),ay,ier)
      CALL EZspline_interp(jz3d_spl,vec(1),vec(2),vec(3),az,ier)

      gf   = one/DSQRT((x_nag-xs)*(x_nag-xs)+(y_nag-ys)*(y_nag-ys)+(z_nag-zs)*(z_nag-zs))
      !PRINT *,xs,ys,zs,ax,ay,ax,gf
      f(1) = norm_3d*ax*gf
      f(2) = norm_3d*ay*gf
      f(3) = norm_3d*az*gf    
      RETURN
      ! END SUBROUTINE
      END SUBROUTINE funsub_nag_a3d
      !-----------------------------------------------------------------
         
         
      !-----------------------------------------------------------------
      SUBROUTINE funsub_nag_b3d(ndim, vec, nfun, f)
      IMPLICIT NONE
      ! INPUT VARIABLES
      INTEGER, INTENT(in) :: ndim, nfun
      DOUBLE PRECISION, INTENT(in) :: vec(ndim)
      DOUBLE PRECISION, INTENT(out) :: f(nfun)
      ! LOCAL VARIABLES
      INTEGER :: ier
      DOUBLE PRECISION :: bn, ax, ay, az , xs, ys, zs, gf, gf3
      ! BEGIN SUBROUTINE

      xs = zero; ys = zero; zs = zero
      ax = zero; ay = zero; az = zero
      CALL EZspline_interp(x3d_spl,vec(1),vec(2),vec(3),xs,ier)
      CALL EZspline_interp(y3d_spl,vec(1),vec(2),vec(3),ys,ier)
      CALL EZspline_interp(z3d_spl,vec(1),vec(2),vec(3),zs,ier)
      CALL EZspline_interp(jx3d_spl,vec(1),vec(2),vec(3),ax,ier)
      CALL EZspline_interp(jy3d_spl,vec(1),vec(2),vec(3),ay,ier)
      CALL EZspline_interp(jz3d_spl,vec(1),vec(2),vec(3),az,ier)

      gf   = one/DSQRT((x_nag-xs)*(x_nag-xs)+(y_nag-ys)*(y_nag-ys)+(z_nag-zs)*(z_nag-zs))
      gf3  = gf*gf*gf
      f(1) = norm_3d*(ay*(z_nag-zs)-az*(y_nag-ys))*gf3
      f(2) = norm_3d*(az*(x_nag-xs)-ax*(z_nag-zs))*gf3
      f(3) = norm_3d*(ax*(y_nag-ys)-ay*(x_nag-xs))*gf3
      !WRITE(427,*) xs,ys,zs
      RETURN
      ! END SUBROUTINE
      END SUBROUTINE funsub_nag_b3d
      !-----------------------------------------------------------------
         
         
      !-----------------------------------------------------------------
      SUBROUTINE vecpot_volint_adapt_flt(x_flt,y_flt,z_flt,ax_flt,ay_flt,az_flt,istat)
      IMPLICIT NONE
      ! INPUT VARIABLES
      REAL, INTENT(in)  :: x_flt, y_flt, z_flt
      REAL, INTENT(out) :: ax_flt, ay_flt,az_flt
      INTEGER, INTENT(inout) :: istat
      ! LOCAL VARIABLES
      DOUBLE PRECISION  :: xt, yt, zt
      DOUBLE PRECISION  :: axt, ayt,azt
      ! BEGIN SUBROUTINE 
      xt   = x_flt
      yt   = y_flt
      zt   = z_flt
      axt = zero
      ayt = zero
      azt = zero
      CALL vecpot_volint_adapt_dbl(xt,yt,zt,axt,ayt,azt,istat)
      ax_flt = axt
      ay_flt = ayt
      az_flt = azt
      RETURN
      ! END SUBROUTINE
      END SUBROUTINE vecpot_volint_adapt_flt
      !-----------------------------------------------------------------
         
         
      !-----------------------------------------------------------------
      SUBROUTINE vecpot_volint_dbl(x,y,z,ax,ay,az)
      IMPLICIT NONE
      ! INPUT VARIABLES
      DOUBLE PRECISION, INTENT(in)  :: x, y, z
      DOUBLE PRECISION, INTENT(out) :: ax, ay, az
      ! LOCAL VARIABLES
      DOUBLE PRECISION ::  gf(nuvp*nr_vc)
      ! BEGIN SUBROUTINE
      gf = one/dsqrt((x-xsurf)**2 + (y-ysurf)**2 + (z-zsurf)**2)
      ax  = norm_3d*SUM(jx_3d*gf)/nuvp
      ay  = norm_3d*SUM(jy_3d*gf)/nuvp
      az  = norm_3d*SUM(jz_3d*gf)/nuvp
      RETURN
      ! END SUBROUTINE
      END SUBROUTINE vecpot_volint_dbl
      !-----------------------------------------------------------------
         
         
      !-----------------------------------------------------------------
      SUBROUTINE bfield_volint_dbl(x,y,z,bx,by,bz)
      IMPLICIT NONE
      ! INPUT VARIABLES
      DOUBLE PRECISION, INTENT(in)  :: x, y, z
      DOUBLE PRECISION, INTENT(out) :: bx, by, bz
      ! LOCAL VARIABLES
      DOUBLE PRECISION ::  gf(nuvp*nr_vc), gf3(nuvp*nr_vc)
      ! BEGIN SUBROUTINE
      gf = one/dsqrt((x-xsurf)**2 + (y-ysurf)**2 + (z-zsurf)**2)
      gf3  = gf*gf*gf
      bx  = norm_3d*SUM((jy_3d*(z-zsurf)-jz_3d*(y-ysurf))*gf3)/nuvp
      by  = norm_3d*SUM((jz_3d*(x-xsurf)-jx_3d*(z-zsurf))*gf3)/nuvp
      bz  = norm_3d*SUM((jx_3d*(y-ysurf)-jy_3d*(x-xsurf))*gf3)/nuvp
      RETURN
      ! END SUBROUTINE
      END SUBROUTINE bfield_volint_dbl
      !-----------------------------------------------------------------
         
         
      !-----------------------------------------------------------------
      SUBROUTINE vecpot_volint_adapt_dbl(x,y,z,ax,ay,az,istat)
      IMPLICIT NONE
      ! INPUT VARIABLES
      DOUBLE PRECISION, INTENT(in)  :: x, y, z
      DOUBLE PRECISION, INTENT(out) :: ax, ay, az
      INTEGER, INTENT(inout) :: istat
      ! LOCAL VARIABLES
      LOGICAL            :: adapt_rerun
      INTEGER(KIND=8), PARAMETER :: ndim_nag = 3 ! theta,zeta,ns
      INTEGER(KIND=8), PARAMETER :: nfun_nag = 3 ! Ax, Ay, Az
      INTEGER(KIND=8), PARAMETER :: lenwrk_nag = IWRK
      INTEGER(KIND=8) :: maxcls_nag,mincls_nag, subs, restar, wrklen, rulcls, wrksbs, n, m, funcls
      DOUBLE PRECISION :: absreq_nag, relreq_nag
      DOUBLE PRECISION :: wrkstr_nag(lenwrk_nag)
      DOUBLE PRECISION :: a_nag(ndim_nag), b_nag(ndim_nag), &
                          finest_nag(nfun_nag), absest_nag(nfun_nag)
      DOUBLE PRECISION, ALLOCATABLE :: vrtwrk(:)

#ifdef NAG
      EXTERNAL :: D01EAF

#else
      EXTERNAL :: dcuhre

#endif
      ! BEGIN SUBROUTINE 
      IF (adapt_tol < 0) THEN
         ! Not implmented
         CALL vecpot_volint_dbl(x,y,z,ax,ay,az)
         RETURN
      END IF
      a_nag(1:3) = zero
      b_nag(1:3) = one
      !b_nag(3)   = one
      mincls_nag = MIN_CLS
      maxcls_nag = IWRK
      absreq_nag = adapt_tol       ! Talk to Stuart about proper values
      relreq_nag = adapt_rel ! Talk to Stuart about proper values
      finest_nag = zero
      absest_nag = zero
      x_nag      = x
      y_nag      = y
      z_nag      = z
      adapt_rerun = .true.
      subs = 1
      restar = 0
      DO WHILE (adapt_rerun)

#ifdef NAG
         CALL D01EAF(ndim_nag,a_nag,b_nag,mincls_nag,maxcls_nag,nfun_nag,funsub_nag_a3d,absreq_nag,&
                   relreq_nag,lenwrk_nag,wrkstr_nag,finest_nag,absest_nag,istat)
         IF (istat == 1 .or. istat == 3) THEN
            maxcls_nag = maxcls_nag*10
            mincls_nag = -1
            restar = 1
            WRITE(6,*) '!!!!!  WARNING Could not reach desired tollerance  !!!!!'
            WRITE(6,*) '  AX = ',finest_nag(1),' +/- ',absest_nag(1)
            WRITE(6,*) '  AY = ',finest_nag(2),' +/- ',absest_nag(2)
            WRITE(6,*) '  AZ = ',finest_nag(3),' +/- ',absest_nag(3)
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         ELSE IF (istat < 0) THEN
            ax = zero
            ay = zero
            az = zero
            adapt_rerun=.false.
         ELSE
            ax = finest_nag(1)
            ay = finest_nag(2)
            az = finest_nag(3)
            adapt_rerun=.false.
         END IF

#else
         IF (.not.ALLOCATED(vrtwrk)) THEN
            wrklen = ((IWRK-ndim_nag)/(2*ndim_nag) + 1)*(2*ndim_nag+2*nfun_nag+2) + 17*nfun_nag + 256
            ALLOCATE(vrtwrk(wrklen),STAT=istat)
            IF (istat .ne. 0) THEN
               WRITE(6,*) ' ALLOCATION ERROR IN: vecpot_volint_adapt_dbl'
               WRITE(6,*) '   VARIABLE: VRTWRK, SIZE: ',wrklen
               WRITE(6,*) '   ISTAT: ',istat
               RETURN
            END IF
         END IF
         CALL dcuhre(ndim_nag,nfun_nag,a_nag,b_nag,mincls_nag,maxcls_nag,funsub_nag_a3d,absreq_nag,&
                     relreq_nag,0,wrklen,restar,finest_nag,absest_nag,funcls,istat,vrtwrk)
         !DEALLOCATE(vrtwrk)
         IF (istat == 1) THEN
            maxcls_nag = maxcls_nag*10
            mincls_nag = -1
            restar = 1
         ELSE IF (istat > 1) THEN
            ax = zero
            ay = zero
            az = zero
            adapt_rerun=.false.
            DEALLOCATE(vrtwrk)
         ELSE
            ax = finest_nag(1)
            ay = finest_nag(2)
            az = finest_nag(3)
            adapt_rerun=.false.
            DEALLOCATE(vrtwrk)
         END IF

#endif
      END DO
      nlastcall=mincls_nag
      RETURN
      ! END SUBROUTINE
      END SUBROUTINE vecpot_volint_adapt_dbl
      !-----------------------------------------------------------------
         
         
      !-----------------------------------------------------------------
      SUBROUTINE bfield_volint_adapt_flt(x_flt,y_flt,z_flt,bx_flt,by_flt,bz_flt,istat)
      IMPLICIT NONE
      ! INPUT VARIABLES
      REAL, INTENT(in)  :: x_flt, y_flt, z_flt
      REAL, INTENT(out) :: bx_flt, by_flt, bz_flt
      INTEGER, INTENT(inout) :: istat
      ! LOCAL VARIABLES
      DOUBLE PRECISION :: xt,yt,zt,bxt,byt,bzt
      ! BEGIN SUBROUTINE
      xt = x_flt
      yt = y_flt
      zt = z_flt
      bxt = zero
      byt = zero
      bzt = zero
      CALL bfield_volint_adapt_dbl(xt,yt,zt,bxt,byt,bzt,istat)
      bx_flt = bxt
      by_flt = byt
      bz_flt = bzt
      RETURN
      ! END SUBROUTINE
      END SUBROUTINE bfield_volint_adapt_flt
      !-----------------------------------------------------------------
         
         
      !-----------------------------------------------------------------
      SUBROUTINE bfield_vc_flt(x_flt,y_flt,z_flt,bx_flt,by_flt,bz_flt,istat)
      IMPLICIT NONE
      ! INPUT VARIABLES
      REAL, INTENT(in)  :: x_flt, y_flt, z_flt
      REAL, INTENT(out) :: bx_flt, by_flt, bz_flt
      INTEGER, INTENT(inout) :: istat
      ! LOCAL VARIABLES
      DOUBLE PRECISION :: xt,yt,zt,bxt,byt,bzt
      ! BEGIN SUBROUTINE
      xt = x_flt
      yt = y_flt
      zt = z_flt
      bxt = zero
      byt = zero
      bzt = zero
      SELECT CASE(TRIM(vc_type_str))
         CASE('Surface Current')
            CALL bfield_virtual_casing_adapt_dbl(xt,yt,zt,bxt,byt,bzt,istat)
         CASE('Volume Integral')
            CALL bfield_volint_adapt_dbl(xt,yt,zt,bxt,byt,bzt,istat)
      END SELECT
      bx_flt = bxt
      by_flt = byt
      bz_flt = bzt
      RETURN
      ! END SUBROUTINE
      END SUBROUTINE bfield_vc_flt
      !-----------------------------------------------------------------
         
         
      !-----------------------------------------------------------------
      SUBROUTINE vecpot_vc_flt(x_flt,y_flt,z_flt,bx_flt,by_flt,bz_flt,istat)
      IMPLICIT NONE
      ! INPUT VARIABLES
      REAL, INTENT(in)  :: x_flt, y_flt, z_flt
      REAL, INTENT(out) :: bx_flt, by_flt, bz_flt
      INTEGER, INTENT(inout) :: istat
      ! LOCAL VARIABLES
      DOUBLE PRECISION :: xt,yt,zt,bxt,byt,bzt
      ! BEGIN SUBROUTINE
      xt = x_flt
      yt = y_flt
      zt = z_flt
      bxt = zero
      byt = zero
      bzt = zero
      SELECT CASE(TRIM(vc_type_str))
         CASE('Surface Current')
            CALL vecpot_virtual_casing_adapt_dbl(xt,yt,zt,bxt,byt,bzt,istat)
         CASE('Volume Integral')
            CALL vecpot_volint_adapt_dbl(xt,yt,zt,bxt,byt,bzt,istat)
      END SELECT
      bx_flt = bxt
      by_flt = byt
      bz_flt = bzt
      RETURN
      ! END SUBROUTINE
      END SUBROUTINE vecpot_vc_flt
      !-----------------------------------------------------------------
         
         
      !-----------------------------------------------------------------
      SUBROUTINE bfield_vc_dbl(x_dbl,y_dbl,z_dbl,bx_dbl,by_dbl,bz_dbl,istat)
      IMPLICIT NONE
      ! INPUT VARIABLES
      DOUBLE PRECISION, INTENT(in)  :: x_dbl, y_dbl, z_dbl
      DOUBLE PRECISION, INTENT(out) :: bx_dbl, by_dbl, bz_dbl
      INTEGER, INTENT(inout) :: istat
      ! LOCAL VARIABLES
      ! BEGIN SUBROUTINE
      SELECT CASE(TRIM(vc_type_str))
         CASE('Surface Current')
            CALL bfield_virtual_casing_adapt_dbl(x_dbl,y_dbl,z_dbl,bx_dbl,by_dbl,bz_dbl,istat)
         CASE('Volume Integral')
            CALL bfield_volint_adapt_dbl(x_dbl,y_dbl,z_dbl,bx_dbl,by_dbl,bz_dbl,istat)
      END SELECT
      RETURN
      ! END SUBROUTINE
      END SUBROUTINE bfield_vc_dbl
      !-----------------------------------------------------------------
         
         
      !-----------------------------------------------------------------
      SUBROUTINE vecpot_vc_dbl(x_dbl,y_dbl,z_dbl,bx_dbl,by_dbl,bz_dbl,istat)
      IMPLICIT NONE
      ! INPUT VARIABLES
      DOUBLE PRECISION, INTENT(in)  :: x_dbl, y_dbl, z_dbl
      DOUBLE PRECISION, INTENT(out) :: bx_dbl, by_dbl, bz_dbl
      INTEGER, INTENT(inout) :: istat
      ! LOCAL VARIABLES
      ! BEGIN SUBROUTINE
      SELECT CASE(TRIM(vc_type_str))
         CASE('Surface Current')
            CALL vecpot_virtual_casing_adapt_dbl(x_dbl,y_dbl,z_dbl,bx_dbl,by_dbl,bz_dbl,istat)
         CASE('Volume Integral')
            CALL vecpot_volint_adapt_dbl(x_dbl,y_dbl,z_dbl,bx_dbl,by_dbl,bz_dbl,istat)
      END SELECT
      RETURN
      ! END SUBROUTINE
      END SUBROUTINE vecpot_vc_dbl
      !-----------------------------------------------------------------
         
         
      !-----------------------------------------------------------------
      SUBROUTINE bfield_volint_adapt_dbl(x,y,z,bx,by,bz,istat)
      IMPLICIT NONE
      ! INPUT VARIABLES
      DOUBLE PRECISION, INTENT(in)  :: x, y, z
      DOUBLE PRECISION, INTENT(out) :: bx, by, bz
      INTEGER, INTENT(inout) :: istat
      ! LOCAL VARIABLES
      LOGICAL            :: adapt_rerun
      INTEGER(KIND=8), PARAMETER :: ndim_nag = 3 ! theta,zeta,s
      INTEGER(KIND=8), PARAMETER :: nfun_nag = 3 ! Bx, By, Bz
      INTEGER(KIND=8), PARAMETER :: lenwrk_nag = IWRK
      INTEGER(KIND=8) :: maxcls_nag,mincls_nag, subs, restar, wrklen, rulcls, wrksbs, n, m, funcls
      DOUBLE PRECISION :: absreq_nag, relreq_nag
      DOUBLE PRECISION :: wrkstr_nag(lenwrk_nag)
      DOUBLE PRECISION :: a_nag(ndim_nag), b_nag(ndim_nag), &
                          finest_nag(nfun_nag), absest_nag(nfun_nag)
      DOUBLE PRECISION, ALLOCATABLE :: vrtwrk(:)

#ifdef NAG
      EXTERNAL :: D01EAF

#else
      EXTERNAL :: dcuhre

#endif
      ! BEGIN SUBROUTINE
      IF (adapt_tol < 0) THEN
         ! Not implmented
         CALL bfield_volint_dbl(x,y,z,bx,by,bz)
         RETURN
      END IF
      a_nag(1:3) = zero
      b_nag(1:3) = one
      !b_nag(3)   = one
      mincls_nag = MIN_CLS
      maxcls_nag = IWRK
      absreq_nag = adapt_tol       ! Talk to Stuart about proper values
      relreq_nag = adapt_rel ! Talk to Stuart about proper values
      finest_nag = zero
      absest_nag = zero
      x_nag      = x
      y_nag      = y
      z_nag      = z
      adapt_rerun = .true.
      subs = 1
      restar = 0
      DO WHILE (adapt_rerun)

#ifdef NAG
         CALL D01EAF(ndim_nag,a_nag,b_nag,mincls_nag,maxcls_nag,nfun_nag,funsub_nag_b3d,absreq_nag,&
                   relreq_nag,lenwrk_nag,wrkstr_nag,finest_nag,absest_nag,istat)
         IF (istat == 1 .or. istat == 3) THEN
            maxcls_nag = maxcls_nag*10
            mincls_nag = -1
            restar = 1
            WRITE(6,*) '!!!!!  WARNING Could not reach desired tollerance  !!!!!'
            WRITE(6,*) '  BX = ',finest_nag(1),' +/- ',absest_nag(1)
            WRITE(6,*) '  BY = ',finest_nag(2),' +/- ',absest_nag(2)
            WRITE(6,*) '  BZ = ',finest_nag(3),' +/- ',absest_nag(3)
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         ELSE IF (istat < 0) THEN
            bx = zero
            by = zero
            bz = zero
            adapt_rerun=.false.
         ELSE
            bx = finest_nag(1)
            by = finest_nag(2)
            bz = finest_nag(3)
            adapt_rerun=.false.
         END IF

#else
         IF (.not.ALLOCATED(vrtwrk)) THEN
            wrklen = ((IWRK-ndim_nag)/(2*ndim_nag) + 1)*(2*ndim_nag+2*nfun_nag+2) + 17*nfun_nag + 256
            ALLOCATE(vrtwrk(wrklen),STAT=istat)
            IF (istat .ne. 0) THEN
               WRITE(6,*) ' ALLOCATION ERROR IN: bfield_volint_adapt_dbl'
               WRITE(6,*) '   VARIABLE: VRTWRK, SIZE: ',wrklen
               WRITE(6,*) '   ISTAT: ',istat
               RETURN
            END IF
         END IF
         CALL dcuhre(ndim_nag,nfun_nag,a_nag,b_nag,mincls_nag,maxcls_nag,funsub_nag_b3d,absreq_nag,&
                     relreq_nag,0,wrklen,restar,finest_nag,absest_nag,funcls,istat,vrtwrk)
         !DEALLOCATE(vrtwrk)
         IF (istat == 1) THEN
            maxcls_nag = maxcls_nag*10
            mincls_nag = funcls
            restar = 1
         ELSE IF (istat > 1) THEN
            bx = zero
            by = zero
            bz = zero
            adapt_rerun=.false.
            DEALLOCATE(vrtwrk)
         ELSE
            bx = finest_nag(1)
            by = finest_nag(2)
            bz = finest_nag(3)
            adapt_rerun=.false.
            DEALLOCATE(vrtwrk)
         END IF

#endif
      END DO
      nlastcall=mincls_nag
      RETURN
      ! END SUBROUTINE
      END SUBROUTINE bfield_volint_adapt_dbl
      !-----------------------------------------------------------------
         
         
      !-----------------------------------------------------------------
      SUBROUTINE mntouv_local(k1,k,mnmax,nu,nv,xu,xv,fmn,xm,xn,f,signs,calc_trig)
      IMPLICIT NONE
      ! INPUT VARIABLES
      INTEGER, INTENT(in) :: k1
      INTEGER, INTENT(in) :: k
      INTEGER, INTENT(in) :: mnmax
      INTEGER, INTENT(in) :: nu
      INTEGER, INTENT(in) :: nv
      DOUBLE PRECISION, INTENT(in) :: xu(1:nu)
      DOUBLE PRECISION, INTENT(in) :: xv(1:nv)           
      DOUBLE PRECISION, INTENT(in) :: fmn(1:mnmax,k1:k)
      INTEGER, INTENT(in) :: xm(1:mnmax)
      INTEGER, INTENT(in) :: xn(1:mnmax)
      DOUBLE PRECISION, INTENT(inout) :: f(1:nu,1:nv,k1:k)
      INTEGER, INTENT(in) :: signs
      INTEGER, INTENT(in) :: calc_trig
      ! LOCAL VARIABLES
      INTEGER     :: mn, i, ier, ik
      DOUBLE PRECISION :: xm_temp(1:mnmax,1)
      DOUBLE PRECISION :: xn_temp(1:mnmax,1)
      DOUBLE PRECISION :: pi2_l
      DOUBLE PRECISION :: mt(1:mnmax,1:nu)
      DOUBLE PRECISION :: nz(1:mnmax,1:nv)
      DOUBLE PRECISION :: fmn_temp(1:mnmax,1:nu)
      DOUBLE PRECISION :: xu_temp(1,1:nu)
      DOUBLE PRECISION :: xv_temp(1,1:nv)
      DOUBLE PRECISION :: fmn_help(1:mnmax)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: cosmt(:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: sinmt(:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: cosnz(:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: sinnz(:,:)
      ! BEGIN SUBROUTINE
      pi2_l = 8.0D+0 * ATAN(1.)
      IF (calc_trig == 1) THEN
         IF (ALLOCATED(cosmt)) DEALLOCATE(cosmt)
         IF (ALLOCATED(sinmt)) DEALLOCATE(sinmt)
         IF (ALLOCATED(cosnz)) DEALLOCATE(cosnz)
         IF (ALLOCATED(sinnz)) DEALLOCATE(sinnz)
         ALLOCATE(cosmt(1:mnmax,1:nu),sinmt(1:mnmax,1:nu),&
                  cosnz(1:mnmax,1:nv),sinnz(1:mnmax,1:nv),STAT=ier)
         FORALL(i=1:mnmax) xm_temp(i,1)=DBLE(xm(i))
         FORALL(i=1:mnmax) xn_temp(i,1)=DBLE(xn(i))
         FORALL(i=1:nu) xu_temp(1,i)=xu(i)
         FORALL(i=1:nv) xv_temp(1,i)=xv(i)
         mt = MATMUL(xm_temp,xu_temp)
         nz = MATMUL(xn_temp,xv_temp)
         FORALL(mn=1:mnmax,i=1:nu) cosmt(mn,i) = dcos(pi2_l*mt(mn,i))
         FORALL(mn=1:mnmax,i=1:nu) sinmt(mn,i) = dsin(pi2_l*mt(mn,i))
         FORALL(mn=1:mnmax,i=1:nv) cosnz(mn,i) = dcos(pi2_l*nz(mn,i))
         FORALL(mn=1:mnmax,i=1:nv) sinnz(mn,i) = dsin(pi2_l*nz(mn,i))
      END IF
      IF (SIGNS == 0) THEN
         DO ik = k1,k
            FORALL(mn=1:mnmax) fmn_help(mn)=fmn(mn,ik)
            fmn_temp=SPREAD(fmn_help,2,nu)
            f(1:nu,1:nv,ik) = f(1:nu,1:nv,ik)  + MATMUL(TRANSPOSE((fmn_temp*cosmt)),cosnz) &
                                   - MATMUL(TRANSPOSE((fmn_temp*sinmt)),sinnz)
         END DO
      ELSE IF (SIGNS == 1) THEN
         DO ik = k1,k
            FORALL(mn=1:mnmax) fmn_help(mn)=fmn(mn,ik)
            fmn_temp=SPREAD(fmn_help,2,nu)
            f(1:nu,1:nv,ik) = f(1:nu,1:nv,ik) + MATMUL(TRANSPOSE((fmn_temp*sinmt)),cosnz) &
                                  + MATMUL(TRANSPOSE((fmn_temp*cosmt)),sinnz)
         END DO
      END IF
      ! END SUBROUTINE
      END SUBROUTINE mntouv_local
      !-----------------------------------------------------------------
         
      !-----------------------------------------------------------------
      SUBROUTINE uvtomn_local(k1,k,mnmax,nu,nv,xu,xv,fmn,xm,xn,f,signs,calc_trig)
      IMPLICIT NONE
      ! INPUT VARIABLES
      INTEGER, INTENT(in) :: k1
      INTEGER, INTENT(in) :: k
      INTEGER, INTENT(in) :: mnmax
      INTEGER, INTENT(in) :: nu
      INTEGER, INTENT(in) :: nv
      DOUBLE PRECISION, INTENT(in) :: xu(1:nu)
      DOUBLE PRECISION, INTENT(in) :: xv(1:nv)
      DOUBLE PRECISION, INTENT(inout) :: fmn(1:mnmax,k1:k)
      INTEGER, INTENT(in) :: xm(1:mnmax)
      INTEGER, INTENT(in) :: xn(1:mnmax)
      DOUBLE PRECISION, INTENT(in) :: f(1:nu,1:nv,k1:k)
      INTEGER, INTENT(in) :: signs
      INTEGER, INTENT(in) :: calc_trig
      ! LOCAL VARIABLES
      INTEGER     :: mn, i, j, ier, ik
      DOUBLE PRECISION :: xn_temp(1:mnmax,1)
      DOUBLE PRECISION :: xm_temp(1:mnmax,1)
      DOUBLE PRECISION :: pi2_l
      DOUBLE PRECISION :: mt(1:mnmax,1:nu)
      DOUBLE PRECISION :: nz(1:mnmax,1:nv)
      DOUBLE PRECISION :: xu_temp(1,1:nu)
      DOUBLE PRECISION :: xv_temp(1,1:nv)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: fnuv(:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: cosmt(:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: sinmt(:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: cosnz(:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: sinnz(:,:)
      ! BEGIN SUBROUTINE
      pi2_l = 8.0D+0 * ATAN(1.)
      IF (calc_trig == 1) THEN
         IF (ALLOCATED(cosmt)) DEALLOCATE(cosmt)
         IF (ALLOCATED(sinmt)) DEALLOCATE(sinmt)
         IF (ALLOCATED(cosnz)) DEALLOCATE(cosnz)
         IF (ALLOCATED(sinnz)) DEALLOCATE(sinnz)
         IF (ALLOCATED(fnuv)) DEALLOCATE(fnuv)
         ALLOCATE(fnuv(1:mnmax),STAT=ier)
         ALLOCATE(cosmt(1:mnmax,1:nu),sinmt(1:mnmax,1:nu),&
                  cosnz(1:mnmax,1:nv),sinnz(1:mnmax,1:nv),STAT=ier)
         FORALL(i=1:mnmax) xm_temp(i,1)=DBLE(xm(i))
         FORALL(i=1:mnmax) xn_temp(i,1)=DBLE(xn(i))
         FORALL(i=1:mnmax) fnuv(i) = 2.0D+0/DBLE(nu*nv)
         WHERE(xm == 0) fnuv = 5.0D-1*fnuv
         FORALL(i=1:nu) xu_temp(1,i)=xu(i)
         FORALL(i=1:nv) xv_temp(1,i)=xv(i)
         mt = MATMUL(xm_temp,xu_temp)
         nz = MATMUL(xn_temp,xv_temp)
         FORALL(mn=1:mnmax,i=1:nu) cosmt(mn,i) = dcos(pi2_l*mt(mn,i))
         FORALL(mn=1:mnmax,i=1:nu) sinmt(mn,i) = dsin(pi2_l*mt(mn,i))
         FORALL(mn=1:mnmax,i=1:nv) cosnz(mn,i) = dcos(pi2_l*nz(mn,i))
         FORALL(mn=1:mnmax,i=1:nv) sinnz(mn,i) = dsin(pi2_l*nz(mn,i))
      END IF
      fmn = zero
      IF (signs == 0) THEN
         DO mn = 1, mnmax
            DO i = 1, nu
               DO j = 1, nv
                  !PRINT *,cosmt(mn,i), cosnz(mn,j), sinmt(mn,i), sinnz(mn,j),fnuv(mn)
                  fmn(mn,k1:k) = fmn(mn,k1:k) + f(i,j,k1:k)*(cosmt(mn,i)*cosnz(mn,j) &
                  - sinmt(mn,i)*sinnz(mn,j))*fnuv(mn)
               END DO
            END DO
         END DO
      ELSE IF (signs == 1) THEN
         DO mn = 1, mnmax
            DO i = 1, nu
               DO j = 1, nv
                  fmn(mn,k1:k) = fmn(mn,k1:k) + f(i,j,k1:k)*(sinmt(mn,i)*cosnz(mn,j) &
                  + cosmt(mn,i)*sinnz(mn,j))*fnuv(mn)
               END DO
            END DO
         END DO
      END IF
      ! END SUBROUTINE
      END SUBROUTINE uvtomn_local
      !-----------------------------------------------------------------
      
!-----------------------------------------------------------------------
!     End Module
!-----------------------------------------------------------------------
      END MODULE virtual_casing_mod
      
      
