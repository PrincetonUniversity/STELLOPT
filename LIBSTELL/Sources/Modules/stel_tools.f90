!-----------------------------------------------------------------------
!     Module:        stel_tools
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          01/20/2016
!     Description:   This module is designed as a replacement to the 
!                    AJAX stellarator utilities.  It provides
!                    routines for the inversion of stellerator
!                    equilibira.
!
!     LOAD_FOURIER_GEOM(k1,k2,mnmax,nu,nv,xm,xn,iflag,rmnc,zmns)
!        Loads a given set of Fourier variables into memory.  The
!        B-Field Fourier Harmonics must be supplied if the user wishes
!        to evaluate the field.  The Lambda variables are only needed
!        if calling the vmec2pest routine.  Assmues Kernel of the form
!        (mu+nv).
!        K1     First radial index (usually 1)
!        K2     Last radial index
!        MNMAX  Number of Fourier modes
!        NU     Number of real space poloidal gridpoints
!        NV     Number of real space toroidal gridpoints
!        XM     Array of poloidal modes (1:MNMAX)
!        XN     Array of toroidal modes (1:MNMAX)
!        IFLAG  Error flag
!        RMNC   Array of R modes (1:MNMAX,K1:K2) (cos)
!        ZMNS   Array of Z modes (1:MNMAX,K1:K2) (sin)
!        (OPTIONAL VARS)
!        RMNS   Array of R modes (1:MNMAX,K1:K2) (sin)
!        ZMNC   Array of Z modes (1:MNMAX,K1:K2) (cos)
!        BSMNS  Array of B^s modes (1:MNMAX,K1:K2) (sin)
!        BUMNC  Array of B^u modes (1:MNMAX,K1:K2) (cos)
!        BVMNC  Array of B^v modes (1:MNMAX,K1:K2) (cos)
!        BSMNC  Array of B^s modes (1:MNMAX,K1:K2) (cos)
!        BUMNS  Array of B^u modes (1:MNMAX,K1:K2) (sin)
!        BVMNS  Array of B^v modes (1:MNMAX,K1:K2) (sin)
!        LMNS   Array of Lambda modes (1:MNMAX,K1:K2) (sin)
!        LMNC   Array of Lambda modes (1:MNMAX,K1:K2) (cos)
!
!     GET_EQUIL_S(r_val,phi_val,z_val,s_val,iflag)
!        Returns the flux coordinates given a real space cylindrical
!        coordinate.
!        R_VAL    R Cylindrical Coordinate
!        PHI_VAL  PHI Cylindrical Coordiante
!        Z_VAL    Z Clyindrical Coorindate
!        S_VAL    Radial coordiante (output)
!        IFLAG    Error Flag
!        (OPTIONAL VARS)
!        U_VAL    Poloidal coordinate (output)
!!!!!!WARNING
!        Where the code takes V as an input it is asking for a value
!        running from 0 to 2*pi over a field period.  However,
!        where the code takes PHI as an input it is asking for
!        the real toroidal angle, 0 to 2*pi over the device.  So
!        PHI = V/NFP.
!-----------------------------------------------------------------------
      MODULE stel_tools
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE EZspline_obj
!-----------------------------------------------------------------------
!     Module Variables
!     PUBLIC
!        LINTSTEPS     Number of discrete steps for line integrals (256)
!                      
!     PRIVATE
!        DOMAIN_FLAG     Used to determin if in domain durring search
!        NFP             Field periodicity determined from XN array
!        R(PHI/Z)_target Used durring search
!        BCS(0/1)        Spline Boundary Conditions
!        PI2             2*pi
!        X_SPL           Spline Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, PARAMETER ::  lintsteps=512
      INTEGER, PARAMETER, PRIVATE ::  neval_max=10000
      INTEGER, PRIVATE      ::  domain_flag, nfp
      DOUBLE PRECISION, PRIVATE ::  R_target, PHI_target, Z_target
      INTEGER, PRIVATE     ::  bcs0(2) = (/ 0, 0/)
      INTEGER, PRIVATE     ::  bcs1(2) = (/-1,-1/)
      DOUBLE PRECISION, PARAMETER,PRIVATE      :: pi2 = 6.283185482025146D+00
      DOUBLE PRECISION, PARAMETER,PRIVATE      :: one = 1.000000000000000D+00
      DOUBLE PRECISION, PARAMETER,PRIVATE      :: search_tol = 1.000000000000000D-12
      DOUBLE PRECISION, DIMENSION(:), POINTER, PRIVATE :: x1, x2, x3
      DOUBLE PRECISION, DIMENSION(:,:), POINTER, PRIVATE ::&
                                    VP2D, GRHO2D, GRHO22D, BAV2D, BSQ2D, &
                                    S112D, S122D, S212D, S222D
      DOUBLE PRECISION, DIMENSION(:,:,:,:), POINTER, PRIVATE ::  R4D, Z4D, &
                                    G4D, RU4D, ZU4D, RV4D, ZV4D, BS4D, BU4D, &
                                    BV4D, B4D, L4D, LU4D, LV4D
      INTEGER, PRIVATE :: win_R4D, win_Z4D, win_G4D, win_RU4D, win_ZU4D, win_RV4D, &
                          win_ZV4D, win_BS4D, win_BU4D, win_BV4D, win_B4D, &   
                          win_L4D, win_LU4D, win_LV4D, win_VP2D, win_GRHO2D, win_GRHO22D, &
                          win_BAV2D, win_BSQ2D, win_S112D, win_S122D, win_S212D, win_S222D, &
                          win_x1, win_x2, win_x3, nx1, nx2, nx3
      REAL*8, PRIVATE :: eps1, eps2, eps3, x1_min, x1_max, x2_min, x2_max, x3_min, x3_max
      REAL*8, parameter, PRIVATE :: small = 1.e-10_ezspline_r8
!-----------------------------------------------------------------------
!     Private Subroutines
!-----------------------------------------------------------------------
      PRIVATE :: mntouv, isingrid, lookupgrid1d, lookupgrid3d
!-----------------------------------------------------------------------
!     INTERFACE Modules
!-----------------------------------------------------------------------
      INTERFACE load_fourier_geom
         MODULE PROCEDURE load_fourier_geom_dbl, load_fourier_geom_sgl
      END INTERFACE
      INTERFACE get_equil_s
         MODULE PROCEDURE get_equil_s_dbl, get_equil_s_sgl
      END INTERFACE get_equil_s
      INTERFACE get_equil_guess_flx
         MODULE PROCEDURE get_equil_guess_flx_dbl, get_equil_guess_flx_sgl
      END INTERFACE get_equil_guess_flx
      INTERFACE get_equil_RZ
         MODULE PROCEDURE get_equil_RZ_dbl, get_equil_RZ_sgl
      END INTERFACE
      INTERFACE get_equil_L
         MODULE PROCEDURE get_equil_L_dbl, get_equil_L_sgl
      END INTERFACE
      INTERFACE get_equil_nhat
         MODULE PROCEDURE get_equil_nhat_dbl, get_equil_nhat_sgl
      END INTERFACE
      INTERFACE get_equil_B
         MODULE PROCEDURE get_equil_B_dbl, get_equil_B_sgl
      END INTERFACE
      INTERFACE get_equil_Bcylsuv
         MODULE PROCEDURE get_equil_Bcylsuv_dbl, get_equil_Bcylsuv_sgl
      END INTERFACE
      INTERFACE get_equil_Bcyl
         MODULE PROCEDURE get_equil_Bcyl_dbl, get_equil_Bcyl_sgl
      END INTERFACE
      INTERFACE get_equil_Bflx
         MODULE PROCEDURE get_equil_Bflx_dbl, get_equil_Bflx_sgl
      END INTERFACE
      INTERFACE get_equil_Bav
         MODULE PROCEDURE get_equil_Bav_dbl, get_equil_Bav_sgl
      END INTERFACE
      INTERFACE pest2vmec
         MODULE PROCEDURE pest2vmec_dbl, pest2vmec_sgl
      END INTERFACE
      INTERFACE get_equil_kappa
         MODULE PROCEDURE get_equil_kappa_dbl, get_equil_kappa_sgl
      END INTERFACE
      INTERFACE get_equil_rho
         MODULE PROCEDURE get_equil_rho_dbl, get_equil_rho_sgl
      END INTERFACE
      INTERFACE get_equil_sus
         MODULE PROCEDURE get_equil_sus_dbl, get_equil_sus_sgl
      END INTERFACE
      INTERFACE line_int
         MODULE PROCEDURE line_int_dbl, line_int_sgl
      END INTERFACE
      INTERFACE line_modB
         MODULE PROCEDURE line_modB_dbl, line_modB_sgl
      END INTERFACE
      CONTAINS
      
      SUBROUTINE load_fourier_geom_dbl(k1,k2,mnmax,nu,nv,xm,xn_in,iflag,rmnc,zmns,&
                                   rmns,zmnc,bsmns,bumnc,bvmnc,bsmnc,bumns,bvmns,&
                                   lmns,lmnc,bmnc,bmns,gmnc,gmns,comm)
      USE EZspline
      USE mpi_sharmem
#if defined(MPI_OPT)
      USE mpi
#endif
      IMPLICIT NONE
      INTEGER, INTENT(in)        :: k1
      INTEGER, INTENT(in)        :: k2
      INTEGER, INTENT(in)        :: mnmax
      INTEGER, INTENT(in)        :: nu
      INTEGER, INTENT(in)        :: nv
      INTEGER, INTENT(in) :: xm(1:mnmax)
      INTEGER, INTENT(in) :: xn_in(1:mnmax)
      INTEGER, INTENT(inout) :: iflag
      DOUBLE PRECISION, INTENT(in) :: rmnc(1:mnmax,k1:k2), zmns(1:mnmax,k1:k2)
      DOUBLE PRECISION, INTENT(in),OPTIONAL :: rmns(1:mnmax,k1:k2), zmnc(1:mnmax,k1:k2)
      DOUBLE PRECISION, INTENT(in),OPTIONAL :: bsmns(1:mnmax,k1:k2),bumnc(1:mnmax,k1:k2),&
                                               bvmnc(1:mnmax,k1:k2)
      DOUBLE PRECISION, INTENT(in),OPTIONAL :: bsmnc(1:mnmax,k1:k2),bumns(1:mnmax,k1:k2),&
                                               bvmns(1:mnmax,k1:k2)
      DOUBLE PRECISION, INTENT(in),OPTIONAL :: lmns(1:mnmax,k1:k2),lmnc(1:mnmax,k1:k2)
      DOUBLE PRECISION, INTENT(in),OPTIONAL :: bmnc(1:mnmax,k1:k2),bmns(1:mnmax,k1:k2)
      DOUBLE PRECISION, INTENT(in),OPTIONAL :: gmnc(1:mnmax,k1:k2),gmns(1:mnmax,k1:k2)
      INTEGER, INTENT(in), OPTIONAL :: comm
      INTEGER ::  ns_t, u, mn, isherm
      INTEGER ::  shar_comm, shar_rank, shar_size
      INTEGER ::  xn(1:mnmax)
      DOUBLE PRECISION, ALLOCATABLE :: xu(:),xv(:),rho(:),vp(:),grho(:),grho2(:)
      DOUBLE PRECISION, ALLOCATABLE :: fmn_temp(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: f_temp(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: gsr(:,:,:),gsp(:,:,:),gsz(:,:,:),gs(:,:,:)
      TYPE(EZspline1_r8) :: Vp_spl, grho_spl, grho2_spl, Bav_spl, Bsq_spl
      TYPE(EZspline1_r8) :: S11_spl, S12_spl, S21_spl, S22_spl
      TYPE(EZspline3_r8) :: R_spl, Z_spl, G_spl
      TYPE(EZspline3_r8) :: Ru_spl, Zu_spl
      TYPE(EZspline3_r8) :: Rv_spl, Zv_spl
      TYPE(EZspline3_r8) :: Bs_spl, Bu_spl, Bv_spl, B_spl
      TYPE(EZspline3_r8) :: L_spl, Lu_spl, Lv_spl

      !Helper vars
      iflag = 0
      ns_t=k2-k1+1
      isherm = 0  ! Cannot change now
      ! Preform checks
      IF (ns_t < 1) iflag = -2
      IF (mnmax< 1) iflag = -3
      IF (nu < 1 .or. nv < 1) iflag = -4
      IF (PRESENT(rmns).NEQV.PRESENT(zmnc)) iflag = -5
      IF (PRESENT(bumnc).NEQV.PRESENT(bvmnc)) iflag = -6
      IF (iflag <0) RETURN
      ! Find NFP
      nfp = 1
      nfp = MINVAL(ABS(xn_in),MASK=(xn_in>0))
      IF (nfp == 0) nfp = 1
      xn = xn_in / nfp
      ! These must be consistent with splines below
      nx1    = nu;  nx2    = nv;   nx3    = ns_t
      x1_min = 0;   x2_min = 0;    x3_min = 0
      x1_max = pi2; x2_max = pi2;  x3_max = 1
      eps1 = (x1_max-x1_min)*small
      eps2 = (x2_max-x2_min)*small
      eps3 = (x3_max-x3_min)*small
      ! Handle Allocating the 4D arrays
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         ! Get rank
         CALL MPI_COMM_SPLIT_TYPE(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shar_comm, iflag)
         CALL MPI_COMM_RANK(shar_comm, shar_rank, iflag)
         CALL MPI_COMM_SIZE(shar_comm, shar_size, iflag)
         ! Free if allocated
         IF (ASSOCIATED(R4D))  CALL mpidealloc(R4D,  win_R4D)
         IF (ASSOCIATED(Z4D))  CALL mpidealloc(Z4D,  win_Z4D)
         IF (ASSOCIATED(G4D))  CALL mpidealloc(G4D,  win_G4D)
         IF (ASSOCIATED(RU4D)) CALL mpidealloc(RU4D, win_RU4D)
         IF (ASSOCIATED(RV4D)) CALL mpidealloc(RV4D, win_RV4D)
         IF (ASSOCIATED(ZU4D)) CALL mpidealloc(ZU4D, win_ZU4D)
         IF (ASSOCIATED(ZV4D)) CALL mpidealloc(ZV4D, win_ZV4D)
         IF (ASSOCIATED(BS4D)) CALL mpidealloc(BS4D, win_BS4D)
         IF (ASSOCIATED(BU4D)) CALL mpidealloc(BU4D, win_BU4D)
         IF (ASSOCIATED(BV4D)) CALL mpidealloc(BV4D, win_BV4D)
         IF (ASSOCIATED(B4D))  CALL mpidealloc(B4D,  win_B4D)
         IF (ASSOCIATED(L4D))  CALL mpidealloc(L4D,  win_L4D)
         IF (ASSOCIATED(LU4D)) CALL mpidealloc(LU4D, win_LU4D)
         IF (ASSOCIATED(LV4D)) CALL mpidealloc(LV4D, win_LV4D)
         IF (ASSOCIATED(x1)) CALL mpidealloc(x1, win_x1)
         IF (ASSOCIATED(x2)) CALL mpidealloc(x2, win_x2)
         IF (ASSOCIATED(x3)) CALL mpidealloc(x3, win_x3)
         ! ALLOCATE
         CALL mpialloc(R4D, 8, nu, nv, ns_t, shar_rank, 0, shar_comm, win_R4D)
         CALL mpialloc(Z4D, 8, nu, nv, ns_t, shar_rank, 0, shar_comm, win_Z4D)
         CALL mpialloc(G4D, 8, nu, nv, ns_t, shar_rank, 0, shar_comm, win_G4D)
         CALL mpialloc(RU4D, 8, nu, nv, ns_t, shar_rank, 0, shar_comm, win_RU4D)
         CALL mpialloc(ZU4D, 8, nu, nv, ns_t, shar_rank, 0, shar_comm, win_ZU4D)
         CALL mpialloc(RV4D, 8, nu, nv, ns_t, shar_rank, 0, shar_comm, win_RV4D)
         CALL mpialloc(ZV4D, 8, nu, nv, ns_t, shar_rank, 0, shar_comm, win_ZV4D)
         CALL mpialloc(BS4D, 8, nu, nv, ns_t, shar_rank, 0, shar_comm, win_BS4D)
         CALL mpialloc(BU4D, 8, nu, nv, ns_t, shar_rank, 0, shar_comm, win_BU4D)
         CALL mpialloc(BV4D, 8, nu, nv, ns_t, shar_rank, 0, shar_comm, win_BV4D)
         CALL mpialloc(B4D, 8, nu, nv, ns_t, shar_rank, 0, shar_comm, win_B4D)
         CALL mpialloc(L4D, 8, nu, nv, ns_t, shar_rank, 0, shar_comm, win_L4D)
         CALL mpialloc(LU4D, 8, nu, nv, ns_t, shar_rank, 0, shar_comm, win_LU4D)
         CALL mpialloc(LV4D, 8, nu, nv, ns_t, shar_rank, 0, shar_comm, win_LV4D)
         CALL mpialloc(x1, nu, shar_rank, 0, shar_comm, win_x1)
         CALL mpialloc(x2, nv, shar_rank, 0, shar_comm, win_x2)
         CALL mpialloc(x3, ns_t, shar_rank, 0, shar_comm, win_x3)
         ! Handle the 1D arrays
         IF (PRESENT(gmnc) .or. PRESENT(gmns)) THEN
            IF (ASSOCIATED(VP2D))    CALL mpidealloc(VP2D,    win_VP2D)
            IF (ASSOCIATED(GRHO2D))  CALL mpidealloc(GRHO2D,  win_GRHO2D)
            IF (ASSOCIATED(GRHO22D)) CALL mpidealloc(GRHO22D, win_GRHO22D)
            IF (ASSOCIATED(BAV2D))   CALL mpidealloc(BAV2D,   win_BAV2D)
            IF (ASSOCIATED(BSQ2D))   CALL mpidealloc(BSQ2D,   win_BSQ2D)
            IF (ASSOCIATED(S112D))   CALL mpidealloc(S112D,   win_S112D)
            IF (ASSOCIATED(S122D))   CALL mpidealloc(S122D,   win_S122D)
            IF (ASSOCIATED(S212D))   CALL mpidealloc(S212D,   win_S212D)
            IF (ASSOCIATED(S222D))   CALL mpidealloc(S222D,   win_S222D)
            CALL mpialloc(VP2D,    2, ns_t, shar_rank, 0, shar_comm, win_VP2D)
            CALL mpialloc(GRHO2D,  2, ns_t, shar_rank, 0, shar_comm, win_GRHO2D)
            CALL mpialloc(GRHO22D, 2, ns_t, shar_rank, 0, shar_comm, win_GRHO22D)
            CALL mpialloc(BAV2D,   2, ns_t, shar_rank, 0, shar_comm, win_BAV2D)
            CALL mpialloc(BSQ2D,   2, ns_t, shar_rank, 0, shar_comm, win_BSQ2D)
            CALL mpialloc(S112D,   2, ns_t, shar_rank, 0, shar_comm, win_S112D)
            CALL mpialloc(S122D,   2, ns_t, shar_rank, 0, shar_comm, win_S122D)
            CALL mpialloc(S212D,   2, ns_t, shar_rank, 0, shar_comm, win_S212D)
            CALL mpialloc(S222D,   2, ns_t, shar_rank, 0, shar_comm, win_S222D)
         END IF
      ELSE
#endif
         shar_rank = 0; shar_size = 1
         IF (ASSOCIATED(R4D))   DEALLOCATE(R4D)
         IF (ASSOCIATED(Z4D))   DEALLOCATE(Z4D)
         IF (ASSOCIATED(G4D))   DEALLOCATE(G4D)
         IF (ASSOCIATED(RU4D))  DEALLOCATE(RU4D)
         IF (ASSOCIATED(RV4D))  DEALLOCATE(RV4D)
         IF (ASSOCIATED(ZU4D))  DEALLOCATE(ZU4D)
         IF (ASSOCIATED(ZV4D))  DEALLOCATE(ZV4D)
         IF (ASSOCIATED(BS4D))  DEALLOCATE(BS4D)
         IF (ASSOCIATED(BU4D))  DEALLOCATE(BU4D)
         IF (ASSOCIATED(BV4D))  DEALLOCATE(BV4D)
         IF (ASSOCIATED(B4D))   DEALLOCATE(B4D)
         IF (ASSOCIATED(L4D))   DEALLOCATE(L4D)
         IF (ASSOCIATED(LU4D))  DEALLOCATE(LU4D)
         IF (ASSOCIATED(LV4D))  DEALLOCATE(LV4D)
         IF (ASSOCIATED(x1))  DEALLOCATE(x1)
         IF (ASSOCIATED(x2))  DEALLOCATE(x2)
         IF (ASSOCIATED(x3))  DEALLOCATE(x3)
         ALLOCATE(R4D(8,nu,nv,ns_t))
         ALLOCATE(Z4D(8,nu,nv,ns_t))
         ALLOCATE(G4D(8,nu,nv,ns_t))
         ALLOCATE(RU4D(8,nu,nv,ns_t))
         ALLOCATE(RV4D(8,nu,nv,ns_t))
         ALLOCATE(ZU4D(8,nu,nv,ns_t))
         ALLOCATE(ZV4D(8,nu,nv,ns_t))
         ALLOCATE(BS4D(8,nu,nv,ns_t))
         ALLOCATE(BU4D(8,nu,nv,ns_t))
         ALLOCATE(BV4D(8,nu,nv,ns_t))
         ALLOCATE(B4D(8,nu,nv,ns_t))
         ALLOCATE(L4D(8,nu,nv,ns_t))
         ALLOCATE(LU4D(8,nu,nv,ns_t))
         ALLOCATE(LV4D(8,nu,nv,ns_t))
         ALLOCATE(x1(nu))
         ALLOCATE(x2(nv))
         ALLOCATE(x3(ns_t))
         IF (PRESENT(gmnc) .or. PRESENT(gmns)) THEN
            IF (ASSOCIATED(VP2D))    DEALLOCATE(VP2D)
            IF (ASSOCIATED(GRHO2D))  DEALLOCATE(GRHO2D)
            IF (ASSOCIATED(GRHO22D)) DEALLOCATE(GRHO22D)
            IF (ASSOCIATED(BAV2D))   DEALLOCATE(BAV2D)
            IF (ASSOCIATED(BSQ2D))   DEALLOCATE(BSQ2D)
            IF (ASSOCIATED(S112D))   DEALLOCATE(S112D)
            IF (ASSOCIATED(S122D))   DEALLOCATE(S122D)
            IF (ASSOCIATED(S212D))   DEALLOCATE(S212D)
            IF (ASSOCIATED(S222D))   DEALLOCATE(S222D)
            ALLOCATE(VP2D(2,ns_t))
            ALLOCATE(GRHO2D(2,ns_t))
            ALLOCATE(GRHO22D(2,ns_t))
            ALLOCATE(BAV2D(2,ns_t))
            ALLOCATE(BSQ2D(2,ns_t))
            ALLOCATE(S112D(2,ns_t))
            ALLOCATE(S122D(2,ns_t))
            ALLOCATE(S212D(2,ns_t))
            ALLOCATE(S222D(2,ns_t))
         END IF
#if defined(MPI_OPT)
      END IF
#endif
      IF (shar_rank == 0) THEN
         !Allocations
         ALLOCATE(xu(nu),xv(nv),rho(k1:k2))
         ALLOCATE(fmn_temp(1:mnmax,k1:k2))
         ALLOCATE(f_temp(nu,nv,k1:k2))
         FORALL(u=k1:k2) rho(u) = REAL(u-1)/REAL(ns_t-1)
         rho = SQRT(rho) ! Improves lookup near axis
         FORALL(u=1:nu) xu(u) = REAL(u-1)/REAL(nu-1)
         FORALL(u=1:nv) xv(u) = REAL(u-1)/REAL(nv-1)
         ! Preform Init
         CALL EZspline_init(R_spl,nu,nv,ns_t,bcs1,bcs1,bcs0,iflag)
         CALL EZspline_init(Z_spl,nu,nv,ns_t,bcs1,bcs1,bcs0,iflag)
         CALL EZspline_init(G_spl,nu,nv,ns_t,bcs1,bcs1,bcs0,iflag)
         CALL EZspline_init(Ru_spl,nu,nv,ns_t,bcs1,bcs1,bcs0,iflag)
         CALL EZspline_init(Rv_spl,nu,nv,ns_t,bcs1,bcs1,bcs0,iflag)
         CALL EZspline_init(Zu_spl,nu,nv,ns_t,bcs1,bcs1,bcs0,iflag)
         CALL EZspline_init(Zv_spl,nu,nv,ns_t,bcs1,bcs1,bcs0,iflag)
         CALL EZspline_init(Bs_spl,nu,nv,ns_t,bcs1,bcs1,bcs0,iflag)
         CALL EZspline_init(Bu_spl,nu,nv,ns_t,bcs1,bcs1,bcs0,iflag)
         CALL EZspline_init(Bv_spl,nu,nv,ns_t,bcs1,bcs1,bcs0,iflag)
         CALL EZspline_init(B_spl,nu,nv,ns_t,bcs1,bcs1,bcs0,iflag)
         CALL EZspline_init(L_spl,nu,nv,ns_t,bcs1,bcs1,bcs0,iflag)
         CALL EZspline_init(Lu_spl,nu,nv,ns_t,bcs1,bcs1,bcs0,iflag)
         CALL EZspline_init(Lv_spl,nu,nv,ns_t,bcs1,bcs1,bcs0,iflag)
         ! R/Z
         f_temp =0
         R_spl%x1 = xu*pi2; R_spl%x2 = xv*pi2; R_spl%x3 = rho; R_spl%isHermite = isherm
         CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,rmnc,xm,xn,f_temp,0,1)
         IF (PRESENT(rmns)) CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,rmns,xm,xn,f_temp,1,0)
         CALL EZspline_setup(R_spl,f_temp,iflag); f_temp=0
         Z_spl%x1 = xu*pi2; Z_spl%x2 = xv*pi2; Z_spl%x3 = rho; Z_spl%isHermite = isherm
         CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,zmns,xm,xn,f_temp,1,0)
         IF (PRESENT(zmnc)) CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,zmnc,xm,xn,f_temp,0,0)
         CALL EZspline_setup(Z_spl,f_temp,iflag); f_temp = 0
         ! dR/Du Derivatives
         Ru_spl%x1 = xu*pi2; Ru_spl%x2 = xv*pi2; Ru_spl%x3 = rho; Ru_spl%isHermite = isherm
         FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -rmnc(mn,:)*xm(mn)
         CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,f_temp,1,0)
         IF (PRESENT(rmns)) THEN
            FORALL(mn = 1:mnmax) fmn_temp(mn,:) = rmns(mn,:)*xm(mn)
            CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,f_temp,0,0)
         END IF
         CALL EZspline_setup(Ru_spl,f_temp,iflag); f_temp = 0
         ! dR/Dv Derivatives
         Rv_spl%x1 = xu*pi2; Rv_spl%x2 = xv*pi2; Rv_spl%x3 = rho; Rv_spl%isHermite = isherm
         FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -rmnc(mn,:)*xn(mn)
         CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,f_temp,1,0)
         IF (PRESENT(rmns)) THEN
            FORALL(mn = 1:mnmax) fmn_temp(mn,:) = rmns(mn,:)*xn(mn)
            CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,f_temp,0,0)
         END IF
         CALL EZspline_setup(Rv_spl,f_temp,iflag); f_temp = 0
         ! dZ/Du Derivatives
         Zu_spl%x1 = xu*pi2; Zu_spl%x2 = xv*pi2; Zu_spl%x3 = rho; Zu_spl%isHermite = isherm
         FORALL(mn = 1:mnmax) fmn_temp(mn,:) = zmns(mn,:)*xm(mn)
         CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,f_temp,0,0)
         IF (PRESENT(zmnc)) THEN
            FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -zmnc(mn,:)*xm(mn)
            CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,f_temp,1,0)
         END IF
         CALL EZspline_setup(Zu_spl,f_temp,iflag); f_temp = 0
         ! dZ/Dv Derivatives
         Zv_spl%x1 = xu*pi2; Zv_spl%x2 = xv*pi2; Zv_spl%x3 = rho; Zv_spl%isHermite = isherm
         FORALL(mn = 1:mnmax) fmn_temp(mn,:) = zmns(mn,:)*xn(mn)
         CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,f_temp,0,0)
         IF (PRESENT(zmnc)) THEN
            FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -zmnc(mn,:)*xn(mn)
            CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,f_temp,1,0)
         END IF
         CALL EZspline_setup(Zv_spl,f_temp,iflag); f_temp = 0
         ! B^s
         Bs_spl%x1 = xu*pi2; Bs_spl%x2 = xv*pi2; Bs_spl%x3 = rho; Bs_spl%isHermite = isherm
         IF (PRESENT(bsmns)) CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,bsmns,xm,xn,f_temp,1,0)
         IF (PRESENT(bsmnc)) CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,bsmnc,xm,xn,f_temp,0,0)
         CALL EZspline_setup(Bs_spl,f_temp,iflag); f_temp = 0
         ! B^u
         Bu_spl%x1 = xu*pi2; Bu_spl%x2 = xv*pi2; Bu_spl%x3 = rho; Bu_spl%isHermite = isherm
         IF (PRESENT(bumnc)) CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,bumnc,xm,xn,f_temp,0,0)
         IF (PRESENT(bumns)) CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,bumns,xm,xn,f_temp,1,0)
         CALL EZspline_setup(Bu_spl,f_temp,iflag); f_temp = 0
         ! B^v
         Bv_spl%x1 = xu*pi2; Bv_spl%x2 = xv*pi2; Bv_spl%x3 = rho; Bv_spl%isHermite = isherm
         IF (PRESENT(bvmnc)) CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,bvmnc,xm,xn,f_temp,0,0)
         IF (PRESENT(bvmns)) CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,bvmns,xm,xn,f_temp,1,0)
         CALL EZspline_setup(Bv_spl,f_temp,iflag); f_temp = 0
         ! ModB
         B_spl%x1 = xu*pi2; B_spl%x2 = xv*pi2; B_spl%x3 = rho; B_spl%isHermite = isherm
         IF (PRESENT(bmnc)) CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,bmnc,xm,xn,f_temp,0,0)
         IF (PRESENT(bmns)) CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,bmns,xm,xn,f_temp,1,0)
         CALL EZspline_setup(B_spl,f_temp,iflag); f_temp = 0
         ! Jacobian sqrt(g)
         G_spl%x1 = xu*pi2; G_spl%x2 = xv*pi2; G_spl%x3 = rho; G_spl%isHermite = isherm
         IF (PRESENT(gmnc)) CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,gmnc,xm,xn,f_temp,0,0)
         IF (PRESENT(gmns)) CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,gmns,xm,xn,f_temp,1,0)
         CALL EZspline_setup(G_spl,f_temp,iflag); f_temp = 0
         ! Lambda
         L_spl%x1 = xu*pi2; L_spl%x2 = xv*pi2; L_spl%x3 = rho; L_spl%isHermite = isherm
         IF (PRESENT(lmns)) CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,lmns,xm,xn,f_temp,1,0)
         IF (PRESENT(lmnc)) CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,lmnc,xm,xn,f_temp,0,0)
         CALL EZspline_setup(L_spl,f_temp,iflag); f_temp = 0
         ! Lambda/u
         Lu_spl%x1 = xu*pi2; Lu_spl%x2 = xv*pi2; Lu_spl%x3 = rho; Lu_spl%isHermite = isherm
         IF (PRESENT(lmns)) THEN
            FORALL(mn = 1:mnmax) fmn_temp(mn,:) = lmns(mn,:)*xm(mn)
            CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,f_temp,0,0)
         END IF
         IF (PRESENT(lmnc)) THEN
            FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -lmnc(mn,:)*xm(mn)
            CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,f_temp,1,0)
         END IF
         CALL EZspline_setup(Lu_spl,f_temp,iflag); f_temp = 0
         ! Lambda/v
         Lv_spl%x1 = xu*pi2; Lv_spl%x2 = xv*pi2; Lv_spl%x3 = rho; Lv_spl%isHermite = isherm
         IF (PRESENT(lmns)) THEN
            FORALL(mn = 1:mnmax) fmn_temp(mn,:) = lmns(mn,:)*xn(mn)
            CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,f_temp,0,0)
         END IF
         IF (PRESENT(lmnc)) THEN
            FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -lmnc(mn,:)*xn(mn)
            CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,f_temp,1,0)
         END IF
         CALL EZspline_setup(Lv_spl,f_temp,iflag); f_temp = 0

         ! Here the part where we copy and delete everything
         x1   = R_SPL%x1
         x2   = R_SPL%x2
         x3   = R_SPL%x3
         R4D  = R_SPL%fspl
         Z4D  = Z_SPL%fspl
         G4D  = G_SPL%fspl
         RU4D = RU_SPL%fspl
         ZU4D = ZU_SPL%fspl
         RV4D = RV_SPL%fspl
         ZV4D = ZV_SPL%fspl
         BS4D = BS_SPL%fspl
         BU4D = BU_SPL%fspl
         BV4D = BV_SPL%fspl
         B4D  = B_SPL%fspl
         L4D  = L_SPL%fspl
         LU4D  = LU_SPL%fspl
         LV4D  = LV_SPL%fspl
         CALL EZspline_free(R_spl,iflag)
         CALL EZspline_free(Z_spl,iflag)
         CALL EZspline_free(G_spl,iflag)
         CALL EZspline_free(RU_spl,iflag)
         CALL EZspline_free(RV_spl,iflag)
         CALL EZspline_free(ZU_spl,iflag)
         CALL EZspline_free(ZV_spl,iflag)
         CALL EZspline_free(BS_spl,iflag)
         CALL EZspline_free(BU_spl,iflag)
         CALL EZspline_free(BV_spl,iflag)
         CALL EZspline_free(B_spl,iflag)
         CALL EZspline_free(L_spl,iflag)
         CALL EZspline_free(LU_spl,iflag)
         CALL EZspline_free(LV_spl,iflag)
         
         ! Calculate rho_s
         IF (PRESENT(gmnc) .or. PRESENT(gmns)) THEN
            ! Allocations
            IF (EZspline_allocated(Vp_spl)) CALL EZspline_free(Vp_spl,iflag)
            IF (EZspline_allocated(grho_spl)) CALL EZspline_free(grho_spl,iflag)
            IF (EZspline_allocated(grho2_spl)) CALL EZspline_free(grho2_spl,iflag)
            IF (EZspline_allocated(S11_spl)) CALL EZspline_free(S11_spl,iflag)
            IF (EZspline_allocated(S12_spl)) CALL EZspline_free(S12_spl,iflag)
            IF (EZspline_allocated(S21_spl)) CALL EZspline_free(S21_spl,iflag)
            IF (EZspline_allocated(S22_spl)) CALL EZspline_free(S22_spl,iflag)
            IF (EZspline_allocated(Bav_spl)) CALL EZspline_free(Bav_spl,iflag)
            IF (EZspline_allocated(Bsq_spl)) CALL EZspline_free(Bsq_spl,iflag)
            CALL EZspline_init(Vp_spl,ns_t,bcs0,iflag)
            CALL EZspline_init(grho_spl,ns_t,bcs0,iflag)
            CALL EZspline_init(grho2_spl,ns_t,bcs0,iflag)
            CALL EZspline_init(S11_spl,ns_t,bcs0,iflag)
            CALL EZspline_init(S12_spl,ns_t,bcs0,iflag)
            CALL EZspline_init(S21_spl,ns_t,bcs0,iflag)
            CALL EZspline_init(S22_spl,ns_t,bcs0,iflag)
            CALL EZspline_init(Bav_spl,ns_t,bcs0,iflag)
            CALL EZspline_init(Bsq_spl,ns_t,bcs0,iflag)
            Vp_spl%x1 = rho; grho_spl%x1 = rho; grho2_spl%x1 = rho
            S11_spl%x1 = rho; S12_spl%x1 = rho; S21_spl%x1 = rho; S22_spl%x1=rho
            Bav_spl%x1 = rho; Bsq_spl%x1 = rho
            ALLOCATE(Vp(k1:k2),grho(k1:k2),grho2(k1:k2))
            ALLOCATE(gsr(nu,nv,k1:k2),gsp(nu,nv,k1:k2),gsz(nu,nv,k1:k2),&
                     gs(nu,nv,k1:k2))
            ! Calc grad(s) components dR/du X dR/dv / sqrt(g)
            !    Note component of R_spl comes from dphi/dphi and cyl coordinates
            gsr = - ZU4D(1,:,:,:)*R4D(1,:,:,:)
            gsp = (ZU4D(1,:,:,:)*RV4D(1,:,:,:) - RU4D(1,:,:,:)*ZV4D(1,:,:,:))*nfp
            gsz =   RU4D(1,:,:,:)*R4D(1,:,:,:)
            f_temp   = G4D(1,:,:,:)
            gs  = (gsr*gsr+gsp*gsp+gsz*gsz)/(f_temp*f_temp)  !|grad(s)|^2
            FORALL(u=k1:k2) gs(:,:,u) = gs(:,:,u)/(4*rho(u)*rho(u)) !|grad(rho)|^2
            ! dV/ds
            Vp = SUM(SUM(f_temp,DIM=1),DIM=1)
            !Vp(1) = 2*Vp(2) - Vp(3)
            ! <|grad(rho)|^2>
            grho2 = SUM(SUM(gs*f_temp,DIM=1),DIM=1)
            grho2 = grho2 / Vp
            grho2(1) = 2*grho2(2) - grho2(3)
            ! <|grad(rho|>|
            gs = sqrt(gs) !|grad(rho)|
            grho = SUM(SUM(gs*f_temp,DIM=1),DIM=1)
            grho = grho / Vp
            grho(1) = 2*grho(2) - grho(3)
            ! Construct splines
            CALL EZspline_setup(Vp_spl,ABS(Vp*pi2*pi2/(nu*nv)),iflag) ! ABS because of negative Jacobian
            CALL EZspline_setup(grho_spl,grho,iflag)
            CALL EZspline_setup(grho2_spl,grho2,iflag)
            f_temp = 0; grho = 0
            ! Calc S11
            f_temp = (RU4D(1,:,:,:)*RU4D(1,:,:,:)+ &
                      ZU4D(1,:,:,:)*ZU4D(1,:,:,:))
            f_temp = f_temp / G4D(1,:,:,:)
            grho   = SUM(SUM(f_temp,DIM=1),DIM=1)/(nu*nv)
            grho(1) = 2*grho(2) - grho(3)
            CALL EZspline_setup(S11_spl,grho,iflag); f_temp = 0; grho = 0
            ! Calc S21
            f_temp = (RU4D(1,:,:,:)*RV4D(1,:,:,:)+ &
                      ZU4D(1,:,:,:)*ZV4D(1,:,:,:))*nfp
            f_temp = f_temp / G4D(1,:,:,:)
            grho   = SUM(SUM(f_temp,DIM=1),DIM=1)/(nu*nv)
            grho(1) = 2*grho(2) - grho(3)
            CALL EZspline_setup(S21_spl,grho,iflag); f_temp = 0; grho = 0
            ! Calc S12
            f_temp = (RU4D(1,:,:,:)*RV4D(1,:,:,:)+ &
                      ZU4D(1,:,:,:)*ZV4D(1,:,:,:))* &
                     (one+LU4D(1,:,:,:))*nfp
            f_temp = f_temp - (RU4D(1,:,:,:)*RU4D(1,:,:,:)+ &
                               ZU4D(1,:,:,:)*ZU4D(1,:,:,:))*&
                              LV4D(1,:,:,:)
            f_temp = f_temp / G4D(1,:,:,:)
            grho   = SUM(SUM(f_temp,DIM=1),DIM=1)/(nu*nv)
            grho(1) = 2*grho(2) - grho(3)
            CALL EZspline_setup(S12_spl,grho,iflag); f_temp = 0; grho = 0
            ! Calc S22
            f_temp = (RV4D(1,:,:,:)*RV4D(1,:,:,:)*nfp*nfp+ &
                      ZV4D(1,:,:,:)*ZV4D(1,:,:,:)*nfp*nfp+ &
                      R4D(1,:,:,:)* R4D(1,:,:,:))* &
                     (one+LU4D(1,:,:,:))
            f_temp = f_temp - (RU4D(1,:,:,:)*RV4D(1,:,:,:)+ &
                               ZU4D(1,:,:,:)*ZV4D(1,:,:,:))*&
                              LV4D(1,:,:,:)*nfp
            f_temp = f_temp / G4D(1,:,:,:)
            grho   = SUM(SUM(f_temp,DIM=1),DIM=1)/(nu*nv)
            !grho(1) = 2*grho(2) - grho(3)
            CALL EZspline_setup(S22_spl,grho,iflag); f_temp = 0; grho = 0
            ! Bav
            f_temp = B4D(1,:,:,:)*G4D(1,:,:,:)
            grho   = SUM(SUM(f_temp,DIM=1),DIM=1)/(nu*nv)
            grho2 = grho2 / Vp
            grho2(1) = 2*grho2(2) - grho2(3)
            CALL EZspline_setup(Bav_spl,grho,iflag); f_temp = 0; grho = 0
            ! Bsq
            f_temp = B4D(1,:,:,:)*B4D(1,:,:,:)*G4D(1,:,:,:)
            grho   = SUM(SUM(f_temp,DIM=1),DIM=1)/(nu*nv)
            grho2 = grho2 / Vp
            grho2(1) = 2*grho2(2) - grho2(3)
            CALL EZspline_setup(Bsq_spl,grho,iflag); f_temp = 0; grho = 0
            ! Deallocate arrays
            DEALLOCATE(gsr,gsp,gsz,gs,Vp,grho,grho2)
            f_temp = 0
            VP2D    = VP_SPL%fspl
            GRHO2D  = GRHO_SPL%fspl
            GRHO22D = GRHO2_SPL%fspl
            BAV2D   = BAV_SPL%fspl
            BSQ2D   = BSQ_SPL%fspl
            S112D   = S11_SPL%fspl
            S122D   = S12_SPL%fspl
            S212D   = S21_SPL%fspl
            S222D   = S22_SPL%fspl
            CALL EZspline_free(VP_spl,iflag)
            CALL EZspline_free(GRHO_spl,iflag)
            CALL EZspline_free(GRHO2_spl,iflag)
            CALL EZspline_free(BAV_spl,iflag)
            CALL EZspline_free(BSQ_spl,iflag)
            CALL EZspline_free(S11_spl,iflag)
            CALL EZspline_free(S12_spl,iflag)
            CALL EZspline_free(S21_spl,iflag)
            CALL EZspline_free(S22_spl,iflag)
         END IF
         ! DEALLOCATIONS
         DEALLOCATE(xu,xv,rho)
         DEALLOCATE(f_temp)
      END IF !So shared memory doesn't do work
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL MPI_BARRIER(shar_comm,iflag)
         CALL MPI_COMM_FREE(shar_comm,iflag)
      END IF
#endif
      RETURN
      END SUBROUTINE load_fourier_geom_dbl

      SUBROUTINE load_fourier_geom_sgl(k1,k2,mnmax,nu,nv,xm,xn_in,iflag,rmnc,zmns,&
                                   rmns,zmnc,bsmns,bumnc,bvmnc,&
                                   bsmnc,bumns,bvmns,&
                                   lmns,lmnc,bmnc,bmns,gmnc,gmns,comm)
      USE EZspline
      IMPLICIT NONE
      INTEGER, INTENT(in)        :: k1
      INTEGER, INTENT(in)        :: k2
      INTEGER, INTENT(in)        :: mnmax
      INTEGER, INTENT(in)        :: nu
      INTEGER, INTENT(in)        :: nv
      INTEGER, INTENT(in) :: xm(1:mnmax)
      INTEGER, INTENT(in) :: xn_in(1:mnmax)
      INTEGER, INTENT(inout) :: iflag
      REAL, INTENT(in) :: rmnc(1:mnmax,k1:k2), zmns(1:mnmax,k1:k2)
      REAL, INTENT(in),OPTIONAL :: rmns(1:mnmax,k1:k2), zmnc(1:mnmax,k1:k2)
      REAL, INTENT(in),OPTIONAL :: bsmns(1:mnmax,k1:k2),bumnc(1:mnmax,k1:k2),&
                                               bvmnc(1:mnmax,k1:k2)
      REAL, INTENT(in),OPTIONAL :: bsmnc(1:mnmax,k1:k2),bumns(1:mnmax,k1:k2),&
                                               bvmns(1:mnmax,k1:k2)
      REAL, INTENT(in),OPTIONAL :: lmns(1:mnmax,k1:k2), lmnc(1:mnmax,k1:k2)
      REAL, INTENT(in),OPTIONAL :: bmnc(1:mnmax,k1:k2), bmns(1:mnmax,k1:k2)
      REAL, INTENT(in),OPTIONAL :: gmnc(1:mnmax,k1:k2), gmns(1:mnmax,k1:k2)
      INTEGER, INTENT(in), OPTIONAL :: comm
      DOUBLE PRECISION :: rmnc_dbl(1:mnmax,k1:k2), zmns_dbl(1:mnmax,k1:k2)
      DOUBLE PRECISION :: rmns_dbl(1:mnmax,k1:k2), zmnc_dbl(1:mnmax,k1:k2)
      DOUBLE PRECISION :: bsmns_dbl(1:mnmax,k1:k2),bumnc_dbl(1:mnmax,k1:k2),&
                                               bvmnc_dbl(1:mnmax,k1:k2)
      DOUBLE PRECISION :: bsmnc_dbl(1:mnmax,k1:k2),bumns_dbl(1:mnmax,k1:k2),&
                                               bvmns_dbl(1:mnmax,k1:k2)
      DOUBLE PRECISION :: lmns_dbl(1:mnmax,k1:k2), lmnc_dbl(1:mnmax,k1:k2)
      DOUBLE PRECISION :: bmnc_dbl(1:mnmax,k1:k2), bmns_dbl(1:mnmax,k1:k2)
      DOUBLE PRECISION :: gmnc_dbl(1:mnmax,k1:k2), gmns_dbl(1:mnmax,k1:k2)
      rmnc_dbl = rmnc
      zmns_dbl = zmns
      rmns_dbl = 0; zmnc_dbl = 0
      bsmns_dbl = 0; bumnc_dbl = 0; bvmnc_dbl = 0
      bsmnc_dbl = 0; bumns_dbl = 0; bvmns_dbl = 0
      lmns_dbl = 0; lmnc_dbl = 0
      bmnc_dbl = 0; bmns_dbl = 0
      gmnc_dbl = 0; gmns_dbl = 0
      IF (PRESENT(rmns)) rmns_dbl = rmns
      IF (PRESENT(zmnc)) zmnc_dbl = zmnc
      IF (PRESENT(bsmns)) bsmns_dbl = bsmns
      IF (PRESENT(bsmnc)) bsmnc_dbl = bsmnc
      IF (PRESENT(bumns)) bumns_dbl = bumns
      IF (PRESENT(bumnc)) bumnc_dbl = bumnc
      IF (PRESENT(bvmns)) bvmns_dbl = bvmns
      IF (PRESENT(bvmnc)) bvmnc_dbl = bvmnc
      IF (PRESENT(lmns)) lmns_dbl = lmns
      IF (PRESENT(lmnc)) lmnc_dbl = lmnc
      IF (PRESENT(bmns)) bmns_dbl = bmns
      IF (PRESENT(bmnc)) bmnc_dbl = bmnc
      IF (PRESENT(comm)) THEN
         CALL load_fourier_geom_dbl(k1,k2,mnmax,nu,nv,xm,xn_in,iflag,rmnc_dbl,zmns_dbl,&
           RMNS=rmns_dbl,ZMNC=zmnc_dbl,&
           BSMNS=bsmns_dbl,BUMNC=bumnc_dbl,BVMNC=bvmnc_dbl,&
           BSMNC=bsmnc_dbl,BUMNS=bumns_dbl,BVMNS=bvmns_dbl,&
           LMNS=lmns_dbl,LMNC=lmnc_dbl,BMNC=bmnc_dbl,BMNS=bmns_dbl,&
           GMNC=gmnc_dbl,GMNS=gmns_dbl,COMM=comm)
      ELSE
         CALL load_fourier_geom_dbl(k1,k2,mnmax,nu,nv,xm,xn_in,iflag,rmnc_dbl,zmns_dbl,&
           RMNS=rmns_dbl,ZMNC=zmnc_dbl,&
           BSMNS=bsmns_dbl,BUMNC=bumnc_dbl,BVMNC=bvmnc_dbl,&
           BSMNC=bsmnc_dbl,BUMNS=bumns_dbl,BVMNS=bvmns_dbl,&
           LMNS=lmns_dbl,LMNC=lmnc_dbl,BMNC=bmnc_dbl,BMNS=bmns_dbl,&
           GMNC=gmnc_dbl,GMNS=gmns_dbl)
      END IF
      RETURN
      END SUBROUTINE load_fourier_geom_sgl
      
      SUBROUTINE rzfunct_stel_tool(m,n,x,fvec,fjac,ldfjac,iflag)
      USE EZspline
      IMPLICIT NONE
      INTEGER :: m,n,ldfjac,iflag, ier
      DOUBLE PRECISION :: x(n),fvec(m),fjac(ldfjac,n)
      DOUBLE PRECISION :: R_temp, Z_temp
      DOUBLE PRECISION :: R_grad(3), Z_grad(3)
      INTEGER :: i,j,k
      REAL*8 :: xparam, yparam, zparam, hx, hy, hz, hxi, hyi, hzi
      REAL*8 :: fval(4)
      INTEGER, parameter :: ict(10)=(/1,1,1,1,0,0,0,0,0,0/)
      IF (x(2) < 0.0) x(2) = x(2) + pi2
      x(2) = MOD(x(2),pi2)
      IF (x(1) < 0) THEN
         x(1) = ABS(x(1))
         x(2) = x(2)+pi2*0.5
         x(2) = MOD(x(2),pi2)
      END IF
      ier = 0; domain_flag = 0
      !CALL EZspline_isInDomain(R_spl,x(2),PHI_Target,x(1),ier)
      !IF (ier .ne. 0) THEN ! Outside domain, extrapolate
      IF (.not. isingrid(x(2),PHI_Target,x(1))) THEN
         !CALL EZspline_interp(R_spl,x(2),PHI_Target,one,R_temp,iflag)
         !CALL EZspline_interp(Z_spl,x(2),PHI_Target,one,Z_temp,iflag)
         !CALL EZspline_gradient(R_spl,x(2),PHI_Target,one,R_grad,iflag)
         !CALL EZspline_gradient(Z_spl,x(2),PHI_Target,one,Z_grad,iflag)
         CALL lookupgrid3d(x(2),PHI_Target,one,i,j,k,hx,hy,hz,hxi,hyi,hzi,xparam,yparam,zparam)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         R4D(1,1,1,1), nx1, nx2, nx3)
         R_temp = fval(1); R_grad(1:3) = fval(2:4)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         Z4D(1,1,1,1), nx1, nx2, nx3)
         Z_temp = fval(1); Z_grad(1:3) = fval(2:4)
         IF (iflag == 1) THEN
            R_temp = R_temp + R_grad(3)*(x(1)-one)
            Z_temp = Z_temp + Z_grad(3)*(x(1)-one)
            fvec(1) = (R_temp - R_target)
            fvec(2) = (Z_temp - Z_target)
         ELSE IF (iflag == 2) THEN
            fjac(1,1) = R_grad(3) !dR/ds
            fjac(1,2) = R_grad(1) !dR/du
            fjac(2,1) = Z_grad(3) !dZ/ds
            fjac(2,2) = Z_grad(1) !dZ/du
         END IF
         domain_flag = -1
         RETURN
      END IF
      ! Inside domain
      IF (iflag == 1) THEN
         !CALL EZspline_interp(R_spl,x(2),PHI_Target,x(1),R_temp,iflag)
         !CALL EZspline_interp(Z_spl,x(2),PHI_Target,x(1),Z_temp,iflag)
         CALL lookupgrid3d(x(2),PHI_Target,x(1),i,j,k,hx,hy,hz,hxi,hyi,hzi,xparam,yparam,zparam)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         R4D(1,1,1,1), nx1, nx2, nx3)
         R_temp = fval(1)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         Z4D(1,1,1,1), nx1, nx2, nx3)
         Z_temp = fval(1)
         fvec(1) = (R_temp - R_target)
         fvec(2) = (Z_temp - Z_target)
      ELSE IF (iflag == 2) THEN
         !CALL EZspline_gradient(R_spl,x(2),PHI_Target,x(1),R_grad,iflag)
         !CALL EZspline_gradient(Z_spl,x(2),PHI_Target,x(1),Z_grad,iflag)
         CALL lookupgrid3d(x(2),PHI_Target,x(1),i,j,k,hx,hy,hz,hxi,hyi,hzi,xparam,yparam,zparam)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         R4D(1,1,1,1), nx1, nx2, nx3)
         R_temp = fval(1); R_grad(1:3) = fval(2:4)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         Z4D(1,1,1,1), nx1, nx2, nx3)
         Z_temp = fval(1); Z_grad(1:3) = fval(2:4)
         fvec(1) = (R_temp - R_target)
         fvec(2) = (Z_temp - Z_target)
         fjac(1,1) = R_grad(3) !dR/ds
         fjac(1,2) = R_grad(1) !dR/du
         fjac(2,1) = Z_grad(3) !dZ/ds
         fjac(2,2) = Z_grad(1) !dZ/du
      END IF
      RETURN
      END SUBROUTINE rzfunct_stel_tool
      
      SUBROUTINE get_equil_s_dbl(r_val,phi_val,z_val,s_val,ier,u_val)
      ! Purpose:
      !   Given a set of real space (cylindrical) coordinates return
      !   a set of flux coordinates.
      !
      ! Error codes:
      !   ier = 9  Value is outside of flux space domain.
      !  
      USE EZspline
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in)    ::  r_val
      DOUBLE PRECISION, INTENT(in)    ::  phi_val
      DOUBLE PRECISION, INTENT(in)    ::  z_val
      DOUBLE PRECISION, INTENT(out)   ::  s_val
      DOUBLE PRECISION, INTENT(out), OPTIONAL   ::  u_val
      INTEGER, INTENT(inout)     ::  ier
      INTEGER, PARAMETER :: mfunct=2
      INTEGER, PARAMETER :: nvars=2
      INTEGER, PARAMETER :: ldfjac=2
      INTEGER :: ik, maxfev_local, nfev, info, njev, maxfev, nprint, mode
      INTEGER, DIMENSION(nvars) :: ipvt
      DOUBLE PRECISION :: ftol,xtol,gtol,factor
      DOUBLE PRECISION, DIMENSION(nvars)  :: xc_opt, diag, qtf, wa1, wa2, wa3
      DOUBLE PRECISION, DIMENSION(mfunct) :: fval,wa4
      DOUBLE PRECISION, DIMENSION(ldfjac,nvars) :: fjac
      !DOUBLE PRECISION :: guess_flx(3)
      !LOGICAL :: l_make_guess

      DOUBLE PRECISION :: enorm
      
      !l_make_guess = .true.
      !
      !IF (l_make_guess) THEN
      !   CALL get_equil_guess_flx(r_val,phi_val,z_val,guess_flx(1),guess_flx(2),guess_flx(3),ier)
      !ELSE
      !   guess_flx(1) = 0.5
      !   guess_flx(2) = pi2*0.6
      !   guess_flx(3) = MOD(phi_val,pi2/nfp)*nfp
      !END IF
      
         
      IF (s_val > 1 .or. s_val < 0) s_val = 0
      IF (ier < 0) RETURN
      ier = 0
      IF (ASSOCIATED(R4D) .and. ASSOCIATED(Z4D)) THEN
         ! These target variables are modules level variables that are used
         ! in rzfunct_stel_tool as part of the optimization loop.
         R_target = r_val
         Z_target = z_val
         PHI_target=phi_val
         IF (PHI_target < 0) PHI_target = PHI_target + pi2
         PHI_target = MOD(PHI_target,pi2/nfp)*nfp
         xc_opt(1) = SQRT(s_val)
         xc_opt(2) = pi2*0.5
         IF (PRESENT(u_val)) xc_opt(2) = u_val
         fval = 1.0E-30
         fjac = 0
         ftol = search_tol
         xtol = search_tol
         gtol = 1.0E-30
         maxfev_local = neval_max
         diag = (/1.0,4.0/)
         mode = 2
         factor = 0.1
         nprint = 0
         info   = 9
         nfev   = 0
         njev   = 0
         ipvt   = 0
         qtf    = 0.0
         wa1 = 0; wa2 = 0; wa3 = 0; wa4 = 0
         domain_flag = 0
         DO ik = 1, 4
            info = 0
            nfev = 0
            njev = 0
            CALL lmder_serial(rzfunct_stel_tool,mfunct,nvars,xc_opt,fval,fjac,ldfjac,ftol,xtol,gtol,&
                    maxfev_local,diag,mode,factor,nprint,info,nfev,njev,ipvt,qtf,&
                    wa1,wa2,wa3,wa4)
            ! Check for absolute chi-sq convergence. This test is not done within the L-M routine.
            IF (info > 0) THEN
               IF (enorm(mfunct,fval)**2 < ftol) THEN
                  info = 1
               ENDIF
            ENDIF
            IF (info < 4) EXIT
            ! If we have not achived success, try to bump the fitter into a better space.
            IF (ik < 4) THEN 
               xc_opt(1) = xc_opt(1) + 1e-4 * 10**(ik-1)
               xc_opt(2) = xc_opt(2) + 1e-4*pi2 * 10**(ik-1)
            END IF
         END DO
         ier = info
         ! Any values of info greater than 0 may indicate success, however
         ! values greater than 3 may indicate a problem in convergence.
         IF ((info > 0) .AND. (info < 4)) ier = 0
         IF (domain_flag .ne. 0) THEN
            ier = 9
            s_val = 1.5
            IF (PRESENT(u_val)) u_val = 2*pi2
            RETURN
         END IF
         s_val = xc_opt(1)*xc_opt(1)
         IF (PRESENT(u_val)) THEN
            u_val = xc_opt(2)
            u_val = MOD(u_val,pi2)
            DO WHILE (u_val < 0)
               u_val = u_val + pi2
            END DO
         END IF
      ELSE
         ier = -1
      END IF
      RETURN
      END SUBROUTINE get_equil_s_dbl
      
      SUBROUTINE get_equil_s_sgl(r_val,phi_val,z_val,s_val,ier,u_val)
      IMPLICIT NONE
      REAL, INTENT(in)    ::  r_val
      REAL, INTENT(in)    ::  phi_val
      REAL, INTENT(in)    ::  z_val
      REAL, INTENT(out)   ::  s_val
      REAL, INTENT(out), OPTIONAL   ::  u_val
      INTEGER, INTENT(inout)     ::  ier
      DOUBLE PRECISION    ::  r_dbl
      DOUBLE PRECISION    ::  phi_dbl
      DOUBLE PRECISION    ::  z_dbl
      DOUBLE PRECISION    ::  s_dbl
      DOUBLE PRECISION    ::  u_dbl
      r_dbl = r_val; phi_dbl = phi_val; z_dbl = z_val
      CALL get_equil_s_dbl(r_dbl,phi_dbl,z_dbl,s_dbl,ier,u_dbl)
      IF (PRESENT(u_val)) u_val = u_val
      s_val = s_dbl
      RETURN
      END SUBROUTINE get_equil_s_sgl

      SUBROUTINE get_equil_guess_flx_dbl(r_val,phi_val,z_val,s_val,u_val,v_val,ier)
      ! Purpose:
      !   Make a guess of flux space coordinates given a set of real space
      !   (cylindrical) coordinates.  This is only valid for toroidal geometries.
      !
      !   For certain geometrys this guess will be reasonably good, however for
      !   many stellarator geometrys this will only be a crude guess.
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in)  ::  r_val
      DOUBLE PRECISION, INTENT(in)  ::  phi_val
      DOUBLE PRECISION, INTENT(in)  ::  z_val
      DOUBLE PRECISION, INTENT(out) ::  s_val
      DOUBLE PRECISION, INTENT(out) ::  u_val
      DOUBLE PRECISION, INTENT(out) ::  v_val
      INTEGER, INTENT(inout) ::  ier

      DOUBLE PRECISION :: axis(3)
      DOUBLE PRECISION :: edge(3)
      DOUBLE PRECISION :: s_temp, u_temp
      IF (ier < 0) RETURN
      
      v_val = MOD(phi_val, pi2/nfp)*nfp
      
      ! Find the axis at the given phi_value.
      s_temp = 0.0
      u_temp = 0.0
      CALL get_equil_rz(s_temp,u_temp,v_val,axis(1),axis(3),ier)
      axis(2) = phi_val
      IF (ier < 0) RETURN
      
      ! Calculate the u_val guess.
      u_val = atan2(z_val-axis(3), r_val-axis(1))
      IF (u_val < 0D0) THEN
         u_val = u_val+pi2
      END IF

      ! Find the edge at the given u and v values.
      s_temp = 1.0
      CALL get_equil_rz(s_temp,u_val,v_val,edge(1),edge(3),ier)
      edge(2) = phi_val
      IF (ier < 0) RETURN
      
      ! Make a guess for s. This assumes that sqrt(s) is approx proportional
      ! to real space distance. 
      s_val = ((r_val-axis(1))**2 + (z_val-axis(3))**2)/((edge(1)-axis(1))**2+(edge(3)-axis(3))**2)

      IF (s_val > 1) THEN
         s_val = 1.0
      ENDIF
      RETURN
      END SUBROUTINE get_equil_guess_flx_dbl
    
      SUBROUTINE get_equil_guess_flx_sgl(r_val,phi_val,z_val,s_val,u_val,v_val,ier)
      IMPLICIT NONE
      REAL, INTENT(in)    ::  r_val
      REAL, INTENT(in)    ::  phi_val
      REAL, INTENT(in)    ::  z_val
      REAL, INTENT(out)   ::  s_val
      REAL, INTENT(out)   ::  u_val
      REAL, INTENT(out)   ::  v_val
      INTEGER, INTENT(inout)     ::  ier
      DOUBLE PRECISION    ::  r_dbl
      DOUBLE PRECISION    ::  phi_dbl
      DOUBLE PRECISION    ::  z_dbl
      DOUBLE PRECISION    ::  s_dbl
      DOUBLE PRECISION    ::  u_dbl
      DOUBLE PRECISION    ::  v_dbl
      r_dbl = r_val; phi_dbl = phi_val; z_dbl = z_val
      CALL get_equil_guess_flx_dbl(r_dbl,phi_dbl,z_dbl,s_dbl,u_dbl,v_dbl,ier)
      s_val = s_dbl; u_val = u_dbl; v_val = v_dbl
      RETURN
      END SUBROUTINE get_equil_guess_flx_sgl
      
      SUBROUTINE get_equil_RZ_dbl(s_val,u_val,v_val,R_val,Z_val,ier,R_grad,Z_grad)
      USE EZspline
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in)    ::  s_val
      DOUBLE PRECISION, INTENT(in)    ::  u_val
      DOUBLE PRECISION, INTENT(in)    ::  v_val
      DOUBLE PRECISION, INTENT(out)   ::  R_val
      DOUBLE PRECISION, INTENT(out)   ::  Z_val
      DOUBLE PRECISION, INTENT(out),OPTIONAL   ::  R_grad(3)
      DOUBLE PRECISION, INTENT(out),OPTIONAL   ::  Z_grad(3)
      DOUBLE PRECISION :: rho_val
      INTEGER, INTENT(inout)     ::  ier
      INTEGER :: i,j,k
      REAL*8 :: xparam, yparam, zparam, hx, hy, hz, hxi, hyi, hzi
      REAL*8 :: fval(4)
      INTEGER, parameter :: ict(10)=(/1,1,1,1,0,0,0,0,0,0/)
      R_val = 0; Z_val = 0
      IF (ier < 0) RETURN
      rho_val = SQRT(s_val)
      !CALL EZspline_isInDomain(R_spl,u_val,v_val,rho_val,ier)
      !IF (ier == 0) THEN
      IF (isingrid(u_val,v_val,rho_val)) THEN
         !CALL EZspline_interp(R_spl,u_val,v_val,rho_val,R_val,ier)
         !CALL EZspline_interp(Z_spl,u_val,v_val,rho_val,Z_val,ier)
         !IF (PRESENT(R_grad)) CALL EZspline_gradient(R_spl,u_val,v_val,rho_val,R_grad,ier)
         !IF (PRESENT(Z_grad)) CALL EZspline_gradient(Z_spl,u_val,v_val,rho_val,Z_grad,ier)
         CALL lookupgrid3d(u_val,v_val,rho_val,i,j,k,hx,hy,hz,hxi,hyi,hzi,xparam,yparam,zparam)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         R4D(1,1,1,1), nx1, nx2, nx3)
         R_val = fval(1)
         IF (PRESENT(R_grad)) R_grad = fval(2:4)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         Z4D(1,1,1,1), nx1, nx2, nx3)
         Z_val = fval(1)
         IF (PRESENT(Z_grad)) Z_grad = fval(2:4)
      ELSE
         ier=9
      END IF
      RETURN
      END SUBROUTINE get_equil_RZ_dbl
      
      SUBROUTINE get_equil_RZ_sgl(s_val,u_val,v_val,R_val,Z_val,ier,R_grad,Z_grad)
      IMPLICIT NONE
      REAL, INTENT(in)    ::  s_val
      REAL, INTENT(in)    ::  u_val
      REAL, INTENT(in)    ::  v_val
      REAL, INTENT(out)   ::  R_val
      REAL, INTENT(out)   ::  Z_val
      INTEGER, INTENT(inout)     ::  ier
      REAL, INTENT(out),OPTIONAL   ::  R_grad(3)
      REAL, INTENT(out),OPTIONAL   ::  Z_grad(3)
      DOUBLE PRECISION    ::  s_dbl
      DOUBLE PRECISION    ::  u_dbl
      DOUBLE PRECISION    ::  v_dbl
      DOUBLE PRECISION    ::  R_dbl
      DOUBLE PRECISION    ::  Z_dbl
      DOUBLE PRECISION    ::  R_grad_dbl(3)
      DOUBLE PRECISION    ::  Z_grad_dbl(3)
      s_dbl = s_val; u_dbl = u_val; v_dbl = v_val
      R_grad_dbl = 0; Z_grad_dbl = 0
      CALL get_equil_RZ_dbl(s_dbl,u_dbl,v_dbl,R_dbl,Z_dbl,ier,&
            R_GRAD=R_grad_dbl,Z_GRAD=Z_grad_dbl)
      R_val = R_dbl; Z_val = Z_dbl
      IF(PRESENT(R_grad)) R_grad = R_grad_dbl
      IF(PRESENT(Z_grad)) Z_grad = Z_grad_dbl
      RETURN
      END SUBROUTINE get_equil_RZ_sgl
      
      SUBROUTINE get_equil_L_dbl(s_val,u_val,v_val,L_val,ier,L_grad)
      USE EZspline
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in)    ::  s_val
      DOUBLE PRECISION, INTENT(in)    ::  u_val
      DOUBLE PRECISION, INTENT(in)    ::  v_val
      DOUBLE PRECISION, INTENT(out)   ::  L_val
      DOUBLE PRECISION, INTENT(out),OPTIONAL   ::  L_grad(3)
      DOUBLE PRECISION :: rho_val
      INTEGER, INTENT(inout)     ::  ier
      INTEGER :: i,j,k
      REAL*8 :: xparam, yparam, zparam, hx, hy, hz, hxi, hyi, hzi
      REAL*8 :: fval(4)
      INTEGER, parameter :: ict(10)=(/1,1,1,1,0,0,0,0,0,0/)
      L_val = 0
      IF (ier < 0) RETURN
      rho_val = SQRT(s_val)
      !CALL EZspline_isInDomain(L_spl,u_val,v_val,rho_val,ier)
      !IF (ier == 0) THEN
      !   CALL EZspline_interp(L_spl,u_val,v_val,rho_val,L_val,ier)
      !   IF (PRESENT(L_grad)) CALL EZspline_gradient(L_spl,u_val,v_val,rho_val,L_grad,ier)
      IF (isingrid(u_val,v_val,rho_val)) THEN
         CALL lookupgrid3d(u_val,v_val,rho_val,i,j,k,hx,hy,hz,hxi,hyi,hzi,xparam,yparam,zparam)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         L4D(1,1,1,1), nx1, nx2, nx3)
         L_val = fval(1)
         IF (PRESENT(L_grad)) L_grad = fval(2:4)
      ELSE
         ier=9
      END IF
      RETURN
      END SUBROUTINE get_equil_L_dbl
      
      SUBROUTINE get_equil_L_sgl(s_val,u_val,v_val,L_val,ier,L_grad)
      IMPLICIT NONE
      REAL, INTENT(in)    ::  s_val
      REAL, INTENT(in)    ::  u_val
      REAL, INTENT(in)    ::  v_val
      REAL, INTENT(out)   ::  L_val
      INTEGER, INTENT(inout)     ::  ier
      REAL, INTENT(out),OPTIONAL   ::  L_grad(3)
      DOUBLE PRECISION    ::  s_dbl
      DOUBLE PRECISION    ::  u_dbl
      DOUBLE PRECISION    ::  v_dbl
      DOUBLE PRECISION    ::  L_dbl
      DOUBLE PRECISION    ::  L_grad_dbl(3)
      s_dbl = s_val; u_dbl = u_val; v_dbl = v_val
      L_grad_dbl = 0
      CALL get_equil_L_dbl(s_dbl,u_dbl,v_dbl,L_dbl,ier,L_GRAD=L_grad_dbl)
      L_val = L_dbl
      IF(PRESENT(L_grad)) L_grad = L_grad_dbl
      RETURN
      END SUBROUTINE get_equil_L_sgl
      
      SUBROUTINE get_equil_rho_dbl(s_val,rho,vp,gradrho,gradrho2,ier)
      USE EZspline
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in)    ::  s_val
      DOUBLE PRECISION, INTENT(out)   ::  rho
      DOUBLE PRECISION, INTENT(out)   ::  vp
      DOUBLE PRECISION, INTENT(out)   ::  gradrho
      DOUBLE PRECISION, INTENT(out)   ::  gradrho2
      INTEGER, INTENT(inout)     ::  ier
      DOUBLE PRECISION :: rho_val
      INTEGER :: k
      REAL*8 :: zparam, hz, hzi
      REAL*8 :: fval(1)
      INTEGER, parameter :: ict(3)=(/1,0,0/)
      rho = -1; vp = 0; gradrho=0; gradrho2=0
      IF (ier < 0) RETURN
      rho_val = SQRT(s_val)
      IF (s_val >= 0 .and. s_val <= 1) THEN
         rho = sqrt(s_val)
         !CALL EZspline_interp(Vp_spl,rho_val,vp,ier)
         !CALL EZspline_interp(grho_spl,rho_val,gradrho,ier)
         !CALL EZspline_interp(grho2_spl,rho_val,gradrho2,ier)
         CALL lookupgrid1d(rho_val,k,hz,hzi,zparam)
         CALL r8fvspline(ict,1,1,fval,k,zparam,hz,hzi,VP2D(1,1),nx3)
         vp = fval(1);
         CALL r8fvspline(ict,1,1,fval,k,zparam,hz,hzi,GRHO2D(1,1),nx3)
         gradrho = fval(1);
         CALL r8fvspline(ict,1,1,fval,k,zparam,hz,hzi,GRHO22D(1,1),nx3)
         gradrho2 = fval(1);
      ELSE
         ier=-1
      END IF
      RETURN
      END SUBROUTINE get_equil_rho_dbl
      
      SUBROUTINE get_equil_rho_sgl(s_val,rho,vp,gradrho,gradrho2,ier)
      USE EZspline
      IMPLICIT NONE
      REAL, INTENT(in)    ::  s_val
      REAL, INTENT(out)   ::  rho
      REAL, INTENT(out)   ::  vp
      REAL, INTENT(out)   ::  gradrho
      REAL, INTENT(out)   ::  gradrho2
      INTEGER, INTENT(inout)     ::  ier
      DOUBLE PRECISION    ::  s_dbl
      DOUBLE PRECISION    ::  rho_dbl
      DOUBLE PRECISION    ::  vp_dbl
      DOUBLE PRECISION    ::  gradrho_dbl
      DOUBLE PRECISION    ::  gradrho2_dbl
      s_dbl = s_val
      CALL get_equil_rho_dbl(s_dbl,rho_dbl,vp_dbl,gradrho_dbl,gradrho2_dbl,ier)
      rho = rho_dbl; vp = vp_dbl; gradrho = gradrho_dbl; gradrho2 = gradrho2_dbl
      RETURN
      END SUBROUTINE get_equil_rho_sgl
      
      SUBROUTINE get_equil_nhat_dbl(s_val,u_val,v_val,nhat,ier)
      USE EZspline
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in)    ::  s_val
      DOUBLE PRECISION, INTENT(in)    ::  u_val
      DOUBLE PRECISION, INTENT(in)    ::  v_val
      DOUBLE PRECISION, INTENT(out)   ::  nhat(3)
      INTEGER, INTENT(inout)     ::  ier
      DOUBLE PRECISION :: rho_val, x_val, y_val, xu,yu, xv, yv, R_val
      DOUBLE PRECISION ::  R_grad(3),Z_grad(3)
      INTEGER :: i,j,k
      REAL*8 :: xparam, yparam, zparam, hx, hy, hz, hxi, hyi, hzi
      REAL*8 :: fval(4)
      INTEGER, parameter :: ict(10)=(/1,1,1,1,0,0,0,0,0,0/)
      R_val = 0; nhat = 0
      IF (ier < 0) RETURN
      rho_val = SQRT(s_val)
      !CALL EZspline_isInDomain(R_spl,u_val,v_val,rho_val,ier)
      !IF (ier == 0) THEN
      IF (isingrid(u_val,v_val,rho_val)) THEN
         CALL lookupgrid3d(u_val,v_val,rho_val,i,j,k,hx,hy,hz,hxi,hyi,hzi,xparam,yparam,zparam)
         !CALL EZspline_interp(R_spl,u_val,v_val,rho_val,R_val,ier)
         !CALL EZspline_gradient(R_spl,u_val,v_val,rho_val,R_grad,ier)
         !CALL EZspline_gradient(Z_spl,u_val,v_val,rho_val,Z_grad,ier)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         R4D(1,1,1,1), nx1, nx2, nx3)
         R_val = fval(1); R_grad = fval(2:4)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         Z4D(1,1,1,1), nx1, nx2, nx3)
         Z_grad = fval(2:4)
         x_val = R_val * DCOS(v_val)
         y_val = R_val * DSIN(v_val)
         xu    = R_grad(1)*DCOS(v_val)
         yu    = R_grad(1)*DSIN(v_val)
         xv    = R_grad(2)*DCOS(v_val) - y_val*pi2/nfp
         yv    = R_grad(2)*DSIN(v_val) - x_val*pi2/nfp
         ! returns snx,sny,snz
         nhat(1) = -yu*Z_grad(2) + yv*Z_grad(1)
         nhat(2) = -xv*Z_grad(1) + xu*Z_grad(2)
         nhat(3) = -xu*yv        + yu*xv
      ELSE
         ier=9
      END IF
      RETURN
      END SUBROUTINE get_equil_nhat_dbl
      
      SUBROUTINE get_equil_nhat_sgl(s_val,u_val,v_val,nhat,ier)
      IMPLICIT NONE
      REAL,    INTENT(in)    ::  s_val
      REAL,    INTENT(in)    ::  u_val
      REAL,    INTENT(in)    ::  v_val
      REAL,    INTENT(out)   ::  nhat(3)
      INTEGER, INTENT(inout) ::  ier
      DOUBLE PRECISION    ::  s_dbl
      DOUBLE PRECISION    ::  u_dbl
      DOUBLE PRECISION    ::  v_dbl
      DOUBLE PRECISION    ::  nhat_dbl(3)
      s_dbl = s_val; u_dbl = u_val; v_dbl = v_val
      nhat_dbl = 0; nhat = 0
      CALL get_equil_nhat_dbl(s_dbl,u_dbl,v_dbl,nhat_dbl,ier)
      nhat = nhat_dbl
      RETURN
      END SUBROUTINE get_equil_nhat_sgl
      
      SUBROUTINE get_equil_sus_dbl(s_val,s11,s12,s21,s22,ier)
      USE EZspline
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in)    ::  s_val
      DOUBLE PRECISION, INTENT(out)   ::  s11
      DOUBLE PRECISION, INTENT(out)   ::  s12
      DOUBLE PRECISION, INTENT(out)   ::  s21
      DOUBLE PRECISION, INTENT(out)   ::  s22
      INTEGER, INTENT(inout)     ::  ier
      DOUBLE PRECISION :: rho_val
      INTEGER :: k
      REAL*8 :: zparam, hz, hzi
      REAL*8 :: fval(1)
      INTEGER, parameter :: ict(3)=(/1,0,0/)
      s11 = 0; s12 = 0; s21=0; s22=0
      IF (ier < 0) RETURN
      rho_val = SQRT(s_val)
      IF (s_val >= 0 .and. s_val <= 1) THEN
         CALL lookupgrid1d(rho_val,k,hz,hzi,zparam)
         CALL r8fvspline(ict,1,1,fval,k,zparam,hz,hzi,S112D(1,1),nx3)
         s11 = fval(1);
         CALL r8fvspline(ict,1,1,fval,k,zparam,hz,hzi,S122D(1,1),nx3)
         s12 = fval(1);
         CALL r8fvspline(ict,1,1,fval,k,zparam,hz,hzi,S212D(1,1),nx3)
         s21 = fval(1);
         CALL r8fvspline(ict,1,1,fval,k,zparam,hz,hzi,S222D(1,1),nx3)
         s22 = fval(1);
         !CALL EZspline_interp(S11_spl,rho_val,s11,ier)
         !CALL EZspline_interp(S12_spl,rho_val,s12,ier)
         !CALL EZspline_interp(S21_spl,rho_val,s21,ier)
         !CALL EZspline_interp(S22_spl,rho_val,s22,ier)
      ELSE
         ier=-1
      END IF
      RETURN
      END SUBROUTINE get_equil_sus_dbl
      
      SUBROUTINE get_equil_sus_sgl(s_val,s11,s12,s21,s22,ier)
      USE EZspline
      IMPLICIT NONE
      REAL, INTENT(in)    ::  s_val
      REAL, INTENT(out)   ::  s11
      REAL, INTENT(out)   ::  s12
      REAL, INTENT(out)   ::  s21
      REAL, INTENT(out)   ::  s22
      INTEGER, INTENT(inout)     ::  ier
      DOUBLE PRECISION    ::  s_dbl
      DOUBLE PRECISION    ::  s11_dbl
      DOUBLE PRECISION    ::  s12_dbl
      DOUBLE PRECISION    ::  s21_dbl
      DOUBLE PRECISION    ::  s22_dbl
      s_dbl = s_val
      CALL get_equil_sus_dbl(s_dbl,s11_dbl,s12_dbl,s21_dbl,s22_dbl,ier)
      s11 = s11_dbl; s12 = s12_dbl; s21 = s21_dbl; s22 = s22_dbl
      RETURN
      END SUBROUTINE get_equil_sus_sgl
      
      SUBROUTINE get_equil_kappa_dbl(s_val,u_val,v_val,kappa,ier)
      USE EZspline
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in)    ::  s_val
      DOUBLE PRECISION, INTENT(in)    ::  u_val
      DOUBLE PRECISION, INTENT(in)    ::  v_val
      DOUBLE PRECISION, INTENT(out)   ::  kappa
      INTEGER, INTENT(inout)     ::  ier
      DOUBLE PRECISION :: rho_val
      DOUBLE PRECISION :: xp, xpp, zp, zpp, denom
      INTEGER :: i,j,k
      REAL*8 :: xparam, yparam, zparam, hx, hy, hz, hxi, hyi, hzi
      REAL*8 :: fval(2)
      INTEGER, parameter :: ict(10)=(/1,1,0,0,0,0,0,0,0,0/)
      kappa = 0
      IF (ier < 0) RETURN
      rho_val = SQRT(s_val)
!      CALL EZspline_isInDomain(R_spl,u_val,v_val,rho_val,ier)
!      IF (ier == 0) THEN
!         !CALL EZspline_derivative(R_spl,1,0,0,u_val,v_val,rho_val,xp,ier)
!         CALL EZspline_interp(Ru_spl,u_val,v_val,rho_val,xp,ier)
!         CALL EZspline_derivative(Ru_spl,1,0,0,u_val,v_val,rho_val,xpp,ier)
!         !CALL EZspline_derivative(Z_spl,1,0,0,u_val,v_val,rho_val,zp,ier)
!         CALL EZspline_interp(Zu_spl,u_val,v_val,rho_val,zp,ier)
!         CALL EZspline_derivative(Zu_spl,1,0,0,u_val,v_val,rho_val,zpp,ier)
      IF (isingrid(u_val,v_val,rho_val)) THEN
         CALL lookupgrid3d(u_val,v_val,rho_val,i,j,k,hx,hy,hz,hxi,hyi,hzi,xparam,yparam,zparam)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         RU4D(1,1,1,1), nx1, nx2, nx3)
         xp = fval(1); xpp = fval(2)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         ZU4D(1,1,1,1), nx1, nx2, nx3)
         zp = fval(1); zpp = fval(2)
         denom = (xp*xp+zp*zp)**1.5
         IF (ABS(denom) > 0) THEN
            kappa = ABS(xp*zpp-zp*xpp)/denom
         END IF
      ELSE
         ier=9
      END IF
      RETURN
      END SUBROUTINE get_equil_kappa_dbl
      
      SUBROUTINE get_equil_kappa_sgl(s_val,u_val,v_val,kappa,ier)
      USE EZspline
      IMPLICIT NONE
      REAL, INTENT(in)    ::  s_val
      REAL, INTENT(in)    ::  u_val
      REAL, INTENT(in)    ::  v_val
      REAL, INTENT(out)   ::  kappa
      INTEGER, INTENT(inout)     ::  ier
      DOUBLE PRECISION    ::  s_dbl
      DOUBLE PRECISION    ::  u_dbl
      DOUBLE PRECISION    ::  v_dbl
      DOUBLE PRECISION    ::  kappa_dbl
      s_dbl = s_val; u_dbl = u_val; v_dbl = v_val
      CALL get_equil_kappa_dbl(s_dbl,u_dbl,v_dbl,kappa_dbl,ier)
      kappa = kappa_dbl
      RETURN
      END SUBROUTINE get_equil_kappa_sgl
      
      SUBROUTINE get_equil_Bcylsuv_dbl(s_val,u_val,v_val,br,bphi,bz,ier,modb_val,B_grad)
      USE EZspline
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in)    ::  s_val
      DOUBLE PRECISION, INTENT(in)    ::  u_val
      DOUBLE PRECISION, INTENT(in)    ::  v_val
      DOUBLE PRECISION, INTENT(out)   ::  br
      DOUBLE PRECISION, INTENT(out)   ::  bphi
      DOUBLE PRECISION, INTENT(out)   ::  bz
      DOUBLE PRECISION, INTENT(out), OPTIONAL   ::  modb_val
      DOUBLE PRECISION, INTENT(out), OPTIONAL   ::  B_grad(3)
      INTEGER, INTENT(inout)     ::  ier
      DOUBLE PRECISION :: rho_val
      DOUBLE PRECISION :: Bs, Bu, Bv, r_val
      DOUBLE PRECISION :: R_grad(3), Z_grad(3)
      INTEGER :: i,j,k
      REAL*8 :: xparam, yparam, zparam, hx, hy, hz, hxi, hyi, hzi
      REAL*8 :: fval(4)
      INTEGER, parameter :: ict(10)=(/1,1,1,1,0,0,0,0,0,0/)
      IF (ier < 0) RETURN
      rho_val = SQRT(s_val)
!      CALL EZSPLINE_isInDomain(R_spl,u_val,v_val,rho_val,ier)
!      IF (ier == 0) THEN
!         R_grad = 0; Z_grad = 0
!         CALL EZspline_interp(R_spl,u_val,v_val,rho_val,r_val,ier)
!         CALL EZspline_interp(Bs_spl,u_val,v_val,rho_val,Bs,ier)
!         CALL EZspline_interp(Bu_spl,u_val,v_val,rho_val,Bu,ier)
!         CALL EZspline_interp(Bv_spl,u_val,v_val,rho_val,Bv,ier)
!         CALL EZspline_interp(Ru_spl,u_val,v_val,rho_val,R_grad(1),ier)
!         CALL EZspline_interp(Rv_spl,u_val,v_val,rho_val,R_grad(2),ier)
!         CALL EZspline_interp(Zu_spl,u_val,v_val,rho_val,Z_grad(1),ier)
!         CALL EZspline_interp(Zv_spl,u_val,v_val,rho_val,Z_grad(2),ier)
      IF (isingrid(u_val,v_val,rho_val)) THEN
         CALL lookupgrid3d(u_val,v_val,rho_val,i,j,k,hx,hy,hz,hxi,hyi,hzi,xparam,yparam,zparam)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         R4D(1,1,1,1), nx1, nx2, nx3)
         r_val = fval(1)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         BS4D(1,1,1,1), nx1, nx2, nx3)
         Bs = fval(1)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         BU4D(1,1,1,1), nx1, nx2, nx3)
         Bu = fval(1)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         BV4D(1,1,1,1), nx1, nx2, nx3)
         Bv = fval(1)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         RU4D(1,1,1,1), nx1, nx2, nx3)
         R_grad(1) = fval(1)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         RV4D(1,1,1,1), nx1, nx2, nx3)
         R_grad(2) = fval(1)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         ZU4D(1,1,1,1), nx1, nx2, nx3)
         Z_grad(1) = fval(1)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         ZV4D(1,1,1,1), nx1, nx2, nx3)
         Z_grad(2) = fval(1)
         br = R_grad(3)*Bs + R_grad(1)*Bu + R_grad(2)*Bv*nfp
         bphi = r_val * Bv
         bz = Z_grad(3)*Bs + Z_grad(1)*Bu + Z_grad(2)*Bv*nfp
!         IF (PRESENT(B_grad)) CALL EZspline_gradient(B_spl,u_val,v_val,s_val,B_grad,ier)
         IF (PRESENT(B_grad))  THEN
            CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                           hx, hxi, hy, hyi, hz, hzi, &
                           B4D(1,1,1,1), nx1, nx2, nx3)
            B_grad = fval(2:4)
         END IF
      ELSE
         ier   = 9
         br    = 0
         bphi  = 0
         bz    = 0
         IF (PRESENT(B_grad)) B_grad = 0
      END IF
      IF (PRESENT(modb_val)) modb_val = sqrt(br*br+bphi*bphi+bz*bz)
      RETURN
      END SUBROUTINE get_equil_Bcylsuv_dbl
      
      SUBROUTINE get_equil_Bcylsuv_sgl(s_val,u_val,v_val,br,bphi,bz,ier,modb_val,B_grad)
      IMPLICIT NONE
      REAL, INTENT(in)    ::  s_val
      REAL, INTENT(in)    ::  u_val
      REAL, INTENT(in)    ::  v_val
      REAL, INTENT(out)   ::  br
      REAL, INTENT(out)   ::  bphi
      REAL, INTENT(out)   ::  bz
      REAL, INTENT(out), OPTIONAL   ::  modb_val
      REAL, INTENT(out), OPTIONAL   ::  B_grad(3)
      INTEGER, INTENT(inout)     ::  ier
      DOUBLE PRECISION    ::  s_dbl
      DOUBLE PRECISION    ::  u_dbl
      DOUBLE PRECISION    ::  v_dbl
      DOUBLE PRECISION    ::  br_dbl
      DOUBLE PRECISION    ::  bphi_dbl
      DOUBLE PRECISION    ::  bz_dbl
      DOUBLE PRECISION    ::  modb_dbl
      DOUBLE PRECISION    ::  B_grad_dbl(3)
      s_dbl = s_val; u_dbl = u_val; v_dbl = v_val
      CALL get_equil_Bcylsuv_dbl(s_dbl,u_dbl,v_dbl,br_dbl,bphi_dbl,bz_dbl,ier,modb_dbl,B_grad_dbl)
      br = br_dbl; bphi = bphi_dbl; bz = bz_dbl
      IF(PRESENT(modb_val)) modb_val = modb_dbl
      IF(PRESENT(B_grad)) B_grad = B_grad_dbl
      RETURN
      END SUBROUTINE get_equil_Bcylsuv_sgl
      
      SUBROUTINE get_equil_B_dbl(r_val,phi_val,z_val,bx,by,bz,ier,modb_val,B_grad)
      USE EZspline
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in)    ::  r_val
      DOUBLE PRECISION, INTENT(in)    ::  phi_val
      DOUBLE PRECISION, INTENT(in)    ::  z_val
      DOUBLE PRECISION, INTENT(out)   ::  bx
      DOUBLE PRECISION, INTENT(out)   ::  by
      DOUBLE PRECISION, INTENT(out)   ::  bz
      DOUBLE PRECISION, INTENT(out), OPTIONAL   ::  modb_val
      DOUBLE PRECISION, INTENT(out), OPTIONAL   ::  B_grad(3)
      INTEGER, INTENT(inout)     ::  ier
      DOUBLE PRECISION :: rho_val
      DOUBLE PRECISION :: s_val, u_val, v_val, br, bphi
      DOUBLE PRECISION :: R_grad(3), Z_grad(3)
      IF (ier < 0) RETURN
      s_val = 0.5
      CALL get_equil_s(r_val,phi_val,z_val,s_val,ier,u_val)
      IF (ier < 0) RETURN
      v_val = PHI_target
      CALL get_equil_Bcylsuv_dbl(s_val,u_val,v_val,br,bphi,bz,ier,modb_val,B_grad)
      bx = br * cos(v_val) - bphi * sin(phi_val)
      by = br * sin(v_val) + bphi * cos(phi_val)
      RETURN
      END SUBROUTINE get_equil_B_dbl
      
      SUBROUTINE get_equil_B_sgl(r_val,phi_val,z_val,bx,by,bz,ier,modb_val,B_grad)
      IMPLICIT NONE
      REAL, INTENT(in)    ::  r_val
      REAL, INTENT(in)    ::  phi_val
      REAL, INTENT(in)    ::  z_val
      REAL, INTENT(out)   ::  bx
      REAL, INTENT(out)   ::  by
      REAL, INTENT(out)   ::  bz
      REAL, INTENT(out), OPTIONAL   ::  modb_val
      REAL, INTENT(out), OPTIONAL   ::  B_grad(3)
      INTEGER, INTENT(inout)     ::  ier
      DOUBLE PRECISION    ::  r_dbl
      DOUBLE PRECISION    ::  phi_dbl
      DOUBLE PRECISION    ::  z_dbl
      DOUBLE PRECISION    ::  bx_dbl
      DOUBLE PRECISION    ::  by_dbl
      DOUBLE PRECISION    ::  bz_dbl
      DOUBLE PRECISION    ::  modb_dbl
      DOUBLE PRECISION    ::  B_grad_dbl(3)
      r_dbl = r_val; phi_dbl = phi_val; z_dbl = z_val
      CALL get_equil_B_dbl(r_dbl,phi_dbl,z_dbl,bx_dbl,by_dbl,bz_dbl,ier,modb_dbl,B_grad_dbl)
      bx = bx_dbl; by = by_dbl; bz = bz_dbl
      IF(PRESENT(modb_val)) modb_val = modb_dbl
      IF(PRESENT(B_grad)) B_grad = B_grad_dbl
      RETURN
      END SUBROUTINE get_equil_B_sgl
      
      SUBROUTINE get_equil_Bcyl_dbl(r_val,phi_val,z_val,Br,Bphi,Bz,ier)
      USE EZspline
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in)    ::  r_val
      DOUBLE PRECISION, INTENT(in)    ::  phi_val
      DOUBLE PRECISION, INTENT(in)    ::  z_val
      DOUBLE PRECISION, INTENT(out)   ::  Br
      DOUBLE PRECISION, INTENT(out)   ::  Bphi
      DOUBLE PRECISION, INTENT(out)   ::  Bz
      INTEGER, INTENT(inout)     ::  ier
      DOUBLE PRECISION :: rho_val
      DOUBLE PRECISION :: s_val, u_val, v_val, Bs, Bu, Bv
      DOUBLE PRECISION :: R_grad(3), Z_grad(3)
      INTEGER :: i,j,k
      REAL*8 :: xparam, yparam, zparam, hx, hy, hz, hxi, hyi, hzi
      REAL*8 :: fval(1)
      INTEGER, parameter :: ict(10)=(/1,0,0,0,0,0,0,0,0,0/)
      IF (ier < 0) RETURN
      CALL get_equil_s(r_val,phi_val,z_val,s_val,ier,u_val)
      IF (ier < 0) RETURN
      v_val = PHI_target
      rho_val = SQRT(s_val)
!      CALL EZSPLINE_isInDomain(R_spl,u_val,v_val,rho_val,ier)
!      IF (ier == 0) THEN
!         CALL EZspline_interp(Bs_spl,u_val,v_val,rho_val,Bs,ier)
!         CALL EZspline_interp(Bu_spl,u_val,v_val,rho_val,Bu,ier)
!         CALL EZspline_interp(Bv_spl,u_val,v_val,rho_val,Bv,ier)
!         CALL EZspline_interp(Ru_spl,u_val,v_val,rho_val,R_grad(1),ier)
!         CALL EZspline_interp(Rv_spl,u_val,v_val,rho_val,R_grad(2),ier)
!         CALL EZspline_interp(Zu_spl,u_val,v_val,rho_val,Z_grad(1),ier)
!         CALL EZspline_interp(Zv_spl,u_val,v_val,rho_val,Z_grad(2),ier)
      IF (isingrid(u_val,v_val,rho_val)) THEN
         R_grad = 0; Z_grad = 0
         CALL lookupgrid3d(u_val,v_val,rho_val,i,j,k,hx,hy,hz,hxi,hyi,hzi,xparam,yparam,zparam)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         BS4D(1,1,1,1), nx1, nx2, nx3)
         Bs = fval(1)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         BU4D(1,1,1,1), nx1, nx2, nx3)
         Bu = fval(1)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         BV4D(1,1,1,1), nx1, nx2, nx3)
         Bv = fval(1)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         RU4D(1,1,1,1), nx1, nx2, nx3)
         R_grad(1) = fval(1)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         RV4D(1,1,1,1), nx1, nx2, nx3)
         R_grad(2) = fval(1)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         ZU4D(1,1,1,1), nx1, nx2, nx3)
         Z_grad(1) = fval(1)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         ZV4D(1,1,1,1), nx1, nx2, nx3)
         Z_grad(2) = fval(1)
         Br = R_grad(3)*Bs + R_grad(1)*Bu + R_grad(2)*Bv*nfp
         Bphi = r_val * Bv
         bz = Z_grad(3)*Bs + Z_grad(1)*Bu + Z_grad(2)*Bv*nfp
      ELSE
         ier   = 9
         Br    = 0
         Bphi  = 0
         Bz    = 0
      END IF
      RETURN
      END SUBROUTINE get_equil_Bcyl_dbl
      
      SUBROUTINE get_equil_Bcyl_sgl(r_val,phi_val,z_val,Br,Bphi,Bz,ier)
      IMPLICIT NONE
      REAL, INTENT(in)    ::  r_val
      REAL, INTENT(in)    ::  phi_val
      REAL, INTENT(in)    ::  z_val
      REAL, INTENT(out)   ::  Br
      REAL, INTENT(out)   ::  Bphi
      REAL, INTENT(out)   ::  Bz
      INTEGER, INTENT(inout)     ::  ier
      DOUBLE PRECISION    ::  r_dbl
      DOUBLE PRECISION    ::  phi_dbl
      DOUBLE PRECISION    ::  z_dbl
      DOUBLE PRECISION    ::  Br_dbl
      DOUBLE PRECISION    ::  Bphi_dbl
      DOUBLE PRECISION    ::  Bz_dbl
      r_dbl = r_val; phi_dbl = phi_val; z_dbl = z_val
      CALL get_equil_Bcyl_dbl(r_dbl,phi_dbl,z_dbl,Br_dbl,Bphi_dbl,Bz_dbl,ier)
      Br = Br_dbl; Bphi = Bphi_dbl; Bz = Bz_dbl
      RETURN
      END SUBROUTINE get_equil_Bcyl_sgl
      
      SUBROUTINE get_equil_Bflx_dbl(s_val,u_val,v_val,bs,bu,bv,ier,modb_val,B_grad)
      USE EZspline
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in)    ::  s_val
      DOUBLE PRECISION, INTENT(in)    ::  u_val
      DOUBLE PRECISION, INTENT(in)    ::  v_val
      DOUBLE PRECISION, INTENT(out)   ::  bs
      DOUBLE PRECISION, INTENT(out)   ::  bu
      DOUBLE PRECISION, INTENT(out)   ::  bv
      DOUBLE PRECISION, INTENT(out), OPTIONAL   ::  modb_val
      DOUBLE PRECISION, INTENT(out), OPTIONAL   ::  B_grad(3)
      INTEGER, INTENT(inout)     ::  ier
      DOUBLE PRECISION :: rho_val
      INTEGER :: i,j,k
      REAL*8 :: xparam, yparam, zparam, hx, hy, hz, hxi, hyi, hzi
      REAL*8 :: fval(4)
      INTEGER, parameter :: ict(10)=(/1,1,1,1,0,0,0,0,0,0/)
      IF (ier < 0) RETURN
      rho_val = SQRT(s_val)
!      CALL EZspline_interp(Bs_spl,u_val,v_val,rho_val,bs,ier)
!      CALL EZspline_interp(Bu_spl,u_val,v_val,rho_val,bu,ier)
!      CALL EZspline_interp(Bv_spl,u_val,v_val,rho_val,bv,ier)
!      IF (PRESENT(modb_val)) CALL EZspline_interp(B_spl,u_val,v_val,rho_val,modb_val,ier)
!      IF (PRESENT(B_grad)) CALL EZspline_gradient(B_spl,u_val,v_val,rho_val,B_grad,ier)
      CALL lookupgrid3d(u_val,v_val,rho_val,i,j,k,hx,hy,hz,hxi,hyi,hzi,xparam,yparam,zparam)
      CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                      hx, hxi, hy, hyi, hz, hzi, &
                      BS4D(1,1,1,1), nx1, nx2, nx3)
      bs = fval(1)
      CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                      hx, hxi, hy, hyi, hz, hzi, &
                      BU4D(1,1,1,1), nx1, nx2, nx3)
      bu = fval(1)
      CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                      hx, hxi, hy, hyi, hz, hzi, &
                      BV4D(1,1,1,1), nx1, nx2, nx3)
      bv = fval(1)
      IF (PRESENT(modb_val) .or. PRESENT(B_GRAD)) &
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         B4D(1,1,1,1), nx1, nx2, nx3)
      IF (PRESENT(modb_val)) modb_val = fval(1)
      IF (PRESENT(B_grad)) B_grad = fval(2:4)
      RETURN
      END SUBROUTINE get_equil_Bflx_dbl
      
      SUBROUTINE get_equil_Bflx_sgl(s_val,u_val,v_val,bs,bu,bv,ier,modb_val,B_grad)
      USE EZspline
      IMPLICIT NONE
      REAL, INTENT(in)    ::  s_val
      REAL, INTENT(in)    ::  u_val
      REAL, INTENT(in)    ::  v_val
      REAL, INTENT(out)   ::  bs
      REAL, INTENT(out)   ::  bu
      REAL, INTENT(out)   ::  bv
      REAL, INTENT(out), OPTIONAL   ::  modb_val
      REAL, INTENT(out), OPTIONAL   ::  B_grad(3)
      INTEGER, INTENT(inout)     ::  ier
      DOUBLE PRECISION    ::  s_dbl
      DOUBLE PRECISION    ::  u_dbl
      DOUBLE PRECISION    ::  v_dbl
      DOUBLE PRECISION   ::  bs_dbl
      DOUBLE PRECISION   ::  bu_dbl
      DOUBLE PRECISION   ::  bv_dbl
      DOUBLE PRECISION   ::  modb_dbl
      DOUBLE PRECISION   ::  B_grad_dbl(3)
      s_dbl = s_val; u_dbl = u_val; v_dbl = v_val;
      CALL get_equil_Bflx_dbl(s_dbl,u_dbl,v_dbl,bs_dbl,bu_dbl,bv_dbl,ier,&
            MODB_VAL=modb_dbl,B_GRAD=B_grad_dbl)
      bs = bs_dbl; bu = bu_dbl; bv = bv_dbl;
      IF (PRESENT(modb_val)) modb_val = modb_dbl
      IF (PRESENT(B_grad)) B_grad = B_grad_dbl
      RETURN
      END SUBROUTINE get_equil_Bflx_sgl

      SUBROUTINE get_equil_Bav_dbl(s_val,Bav,Bsqav,ier,Bsqavp_val)
      USE EZspline
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in)    ::  s_val
      DOUBLE PRECISION, INTENT(out)   ::  Bav
      DOUBLE PRECISION, INTENT(out)   ::  Bsqav
      DOUBLE PRECISION, INTENT(out), OPTIONAL   ::  Bsqavp_val
      INTEGER, INTENT(inout)     ::  ier
      DOUBLE PRECISION :: rho_val, vp_val
      INTEGER :: k
      REAL*8 :: zparam, hz, hzi
      REAL*8 :: fval(2)
      INTEGER, parameter :: ict(3)=(/1,1,0/)
      IF (ier < 0) RETURN
      rho_val = SQRT(s_val)
      !CALL EZspline_interp(Bav_spl,rho_val,Bav,ier)
      !CALL EZspline_interp(Bsq_spl,rho_val,Bsqav,ier)
      CALL lookupgrid1d(rho_val,k,hz,hzi,zparam)
      CALL r8fvspline(ict,1,1,fval,k,zparam,hz,hzi,BAV2D(1,1),nx3)
      Bav = fval(1)
      CALL r8fvspline(ict,1,1,fval,k,zparam,hz,hzi,BSQ2D(1,1),nx3)
      Bsqav = fval(1)
      IF (PRESENT(Bsqavp_val))  THEN
         !CALL EZspline_derivative(Bsq_spl,1,rho_val,Bsqavp_val,ier)
         !CALL EZspline_interp(Vp_spl,rho_val,vp_val,ier)
         Bsqavp_val = fval(2)
         CALL r8fvspline(ict,1,1,fval,k,zparam,hz,hzi,VP2D(1,1),nx3)
         vp_val = fval(1)
         Bsqavp_val = 2*rho_val*Bsqavp_val/vp_val  ! d/dV = (dPhi/drho)*(dV/dPhi)^-1 * d/drho
      END IF
      RETURN
      END SUBROUTINE get_equil_Bav_dbl

      SUBROUTINE get_equil_Bav_sgl(s_val,Bav,Bsqav,ier,Bsqavp_val)
      USE EZspline
      IMPLICIT NONE
      REAL, INTENT(in)    ::  s_val
      REAL, INTENT(out)   ::  Bav
      REAL, INTENT(out)   ::  Bsqav
      REAL, INTENT(out), OPTIONAL   ::  Bsqavp_val
      INTEGER, INTENT(inout)     ::  ier
      DOUBLE PRECISION    ::  s_dbl
      DOUBLE PRECISION   ::  Bav_dbl
      DOUBLE PRECISION   ::  Bsqav_dbl
      DOUBLE PRECISION   ::  Bsqavp_dbl
      s_dbl = s_val
      CALL get_equil_Bav_dbl(s_dbl,Bav_dbl,Bsqav_dbl,ier,BSQAVP_VAL=Bsqavp_dbl)
      Bav = Bav_dbl
      Bsqav = Bsqav_dbl
      IF (PRESENT(Bsqavp_val)) Bsqavp_val = Bsqavp_dbl
      RETURN
      END SUBROUTINE get_equil_Bav_sgl

      SUBROUTINE pest2vmec_dbl(coord)
      USE EZspline
      IMPLICIT none 
      DOUBLE PRECISION,INTENT(inout) :: coord(3)
      INTEGER :: ier, n1, n2
      DOUBLE PRECISION :: rho_val
      DOUBLE PRECISION :: th, th1, dth, phi, s, lam, dlam
      DOUBLE PRECISION,PARAMETER :: eps_newt = 1.0D-12
      INTEGER :: i,j,k
      REAL*8 :: xparam, yparam, zparam, hx, hy, hz, hxi, hyi, hzi
      REAL*8 :: fval(1)
      INTEGER, parameter :: ict(10)=(/1,0,0,0,0,0,0,0,0,0/)
      s = coord(1)
      th = coord(2)
      phi = coord(3)
      phi = MOD(phi,pi2/nfp)*nfp
      IF (phi < 0) THEN
         phi = -MOD(ABS(phi),pi2)
         phi = phi + pi2
      END IF
      IF (th < 0) THEN
         th = -MOD(ABS(th),pi2)
         th = th + pi2
      END IF
      IF (phi > pi2) phi = DMOD(phi,pi2)
      IF (th > pi2) th = DMOD(th,pi2)
      th1 = th
      dth = one
      n1 = 0
      rho_val = SQRT(s)
      DO WHILE(ABS(dth) >= search_tol .and. n1 < 500)
         IF (th < 0) th = th + pi2
         IF (th > pi2) th = MOD(th,pi2)
         !CALL EZSpline_interp(L_spl,th,phi,rho_val,lam,ier)
         !CALL EZSpline_interp(Lu_spl,th,phi,rho_val,dlam,ier)
         CALL lookupgrid3d(th,phi,rho_val,i,j,k,hx,hy,hz,hxi,hyi,hzi,xparam,yparam,zparam)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         L4D(1,1,1,1), nx1, nx2, nx3)
         lam = fval(1)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         LU4D(1,1,1,1), nx1, nx2, nx3)
         dlam = fval(1)
         dth = -(th + lam - th1)/(one+dlam)
         n1 = n1 + 1
         th = th + 0.5*dth
      END DO
      coord(1) = s
      coord(2) = th
      coord(3) = phi
      END SUBROUTINE pest2vmec_dbl

      SUBROUTINE pest2vmec_sgl(coord)
      IMPLICIT none 
      REAL,INTENT(inout) :: coord(3)
      DOUBLE PRECISION :: coord_dbl(3)
      coord_dbl = coord
      CALL pest2vmec_dbl(coord_dbl)
      coord = coord_dbl
      RETURN
      END SUBROUTINE pest2vmec_sgl
      
      SUBROUTINE line_modb_dbl(r1,r2,target_B,s,length)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in)   :: r1(3)
      DOUBLE PRECISION, INTENT(inout)   :: r2(3)
      DOUBLE PRECISION, INTENT(in)  :: target_B
      DOUBLE PRECISION, INTENT(out) :: s
      DOUBLE PRECISION, INTENT(out) :: length
      INTEGER     :: i, ik, ier
      DOUBLE PRECISION :: x1, y1, z1, x2, y2, z2, dx, dy, dz, xp, yp, zp,s1
      DOUBLE PRECISION :: Br, Bphi, Bz, modB1, modB2, deltaB, dl, dldB
      DOUBLE PRECISION :: phip, rp
      DOUBLE PRECISION, PARAMETER :: eps = 1.0E-3
      DOUBLE PRECISION, PARAMETER :: alpha = 0.7
      INTEGER, PARAMETER :: max_iter = 100
      phip = r1(2)
      s1    = r2(2)
      ier = 0
      x1 = r1(1)*cos(phip); x2 = r2(1)*cos(s1)
      y1 = r1(1)*sin(phip); y2 = r2(1)*sin(s1)
      z1 = r1(3);           z2 = r2(3);
      dx = (x2-x1)/lintsteps
      dy = (y2-y1)/lintsteps
      dz = (z2-z1)/lintsteps
      DO i = 1, lintsteps ! Search to edge
         xp = x1+dx*(i-1)
         yp = y1+dy*(i-1)
         zp = z1+dz*(i-1)
         rp = sqrt(xp*xp+yp*yp)
         phip = ATAN2(yp,xp)
         IF (phip < 0) phip = phip+pi2
         CALL get_equil_s(rp,phip,zp,s1,ier)
         IF (s1 <= 1) THEN
            x2 = xp+dx; x1 = xp
            y2 = yp+dy; y1 = yp
            z2 = zp+dz; z1 = zp
            EXIT
         END IF
      END DO
      IF (i== lintsteps) THEN ! In this case we miss the plasma
         r2(1) = -one; r2(2) = 0.0; r2(3) = 0.0
         s = 1.5
         length = 0
         RETURN
      END IF
      ik = 0
      DO WHILE (ik < max_iter)
         !WRITE(327,*) x1,y1,z1
         rp = sqrt(x1*x1+y1*y1)
         phip = ATAN2(y1,x1)
         IF (phip < 0) phip = phip+pi2
         CALL get_equil_Bcyl(rp,phip,z1,Br,Bphi,Bz,ier)
         modB1 = SQRT(Br*Br+Bphi*Bphi+Bz*Bz)
         !WRITE(329,*) Br,Bphi,Bz
         deltaB = target_B-modb1
      	 IF (ABS(deltaB) < eps) THEN
      	    CALL get_equil_s(rp,phip,zp,s1,ier)
      	    r2(1) = rp; r2(2) = phip; r2(3) = zp; s = s1
      	    x2 = r1(1)*cos(r1(2));
      	    y2 = r1(1)*sin(r1(2));
      	    z2 = r1(3);
      	    dx = x1-x2; dy = y1-y2; dz = z1-z2
      	    length = SQRT(dx*dx+dy*dy+dz*dz)
      	    EXIT
      	 END IF
         rp = sqrt(x2*x2+y2*y2)
         phip = ATAN2(y2,x2)
         IF (phip < 0) phip = phip+pi2
         CALL get_equil_Bcyl(rp,phip,z2,Br,Bphi,Bz,ier)
         modB2 = SQRT(Br*Br+Bphi*Bphi+Bz*Bz)
         !WRITE(328,*) modB1,modB2
         IF (modB1==0 .or. modB2==0) THEN
            r2(1) = -one; r2(2) = 0.0; r2(3) = 0.0
            s = 1.5
            length = 0
            RETURN
         END IF
         IF (modB1 < target_B .and. target_B < modB2) THEN
      	    dl = SQRT(dx*dx+dy*dy+dz*dz)
      	    dldB = dl/(modB2-modB1)
      	    dl = dldB*deltaB
            x1 = x1 + alpha*dl*dx; x2 = x2 - alpha*dl*dx
            y1 = y1 + alpha*dl*dy; y2 = y2 - alpha*dl*dy
            z1 = z1 + alpha*dl*dz; z2 = z2 - alpha*dl*dz
         ELSE
            x1 = x2; y1 = y2; z1 = z2;
            x2 = x2 + dx
            y2 = y2 + dy
            z2 = z2 + dz
         END IF
      	 ik = ik + 1
      END DO
      RETURN
      END SUBROUTINE line_modb_dbl
      
      SUBROUTINE line_modb_sgl(r1,r2,target_B,s,length)
      IMPLICIT NONE
      REAL, INTENT(in)   :: r1(3)
      REAL, INTENT(inout)   :: r2(3)
      REAL, INTENT(in)  :: target_B
      REAL, INTENT(out) :: s
      REAL, INTENT(out) :: length
      DOUBLE PRECISION :: r1_dbl(3)
      DOUBLE PRECISION :: r2_dbl(3)
      DOUBLE PRECISION :: target_B_dbl
      DOUBLE PRECISION :: s_dbl
      DOUBLE PRECISION :: length_dbl
      r1_dbl = r1; r2_dbl = r2; target_B_dbl=target_B
      CALL line_modb_dbl(r1_dbl,r2_dbl,target_B_dbl,s_dbl,length_dbl)
      r2 = r2_dbl
      s  = s_dbl
      length = length_dbl
      RETURN
      END SUBROUTINE line_modb_sgl
      
      SUBROUTINE line_int_dbl(fcn,r1,r2,val,length)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in)   :: r1(3), r2(3)
      DOUBLE PRECISION, INTENT(out)  :: val
      DOUBLE PRECISION, INTENT(out), OPTIONAL :: length
      INTEGER     :: i, j, k, ier
      DOUBLE PRECISION :: x1, y1, z1, x2, y2, z2, dx, dy, dz, xp, yp, zp, int_fac
      DOUBLE PRECISION :: xp1, yp1, zp1, delf, delt, rp, phip,s, dx2, dy2, dz2,uval
      DOUBLE PRECISION :: f_val
      INTEGER, PARAMETER :: nop=3
      INTEGER, PARAMETER :: int_step=2
      DOUBLE PRECISION, dimension(nop), parameter :: ci=(/1./6.,2./3.,1./6./)
      EXTERNAL fcn
      phip = r1(2)
      s    = r2(2)
      ier = 0; f_val = 0
      x1 = r1(1)*cos(phip); x2 = r2(1)*cos(s)
      y1 = r1(1)*sin(phip); y2 = r2(1)*sin(s)
      z1 = r1(3);           z2 = r2(3);
      dx = (x2-x1)/lintsteps
      dy = (y2-y1)/lintsteps
      dz = (z2-z1)/lintsteps
      s=0; phip=0
      DO i = 1, lintsteps ! Get first boundary point
         xp = x1+dx*(i-1)
         yp = y1+dy*(i-1)
         zp = z1+dz*(i-1)
         rp = sqrt(xp*xp+yp*yp)
         phip = ATAN2(yp,xp)
         IF (phip < 0) phip = phip+pi2
         CALL get_equil_s(rp,phip,zp,s,ier)
         IF (s <= 1) THEN
            x1 = xp-dx
            y1 = yp-dy
            z1 = zp-dz
            EXIT
         END IF
      END DO
      IF (i== lintsteps) THEN
         val = 0
         IF (PRESENT(length)) length = 0
         RETURN
      END IF
      DO i = 1, lintsteps ! Get second boundary point
         xp = x2-dx*(i-1)
         yp = y2-dy*(i-1)
         zp = z2-dz*(i-1)
         IF ((xp == x1) .and. (yp == y1) .and. (zp == z1)) THEN
            val = 0
            RETURN
         END IF
         rp = sqrt(xp*xp+yp*yp)
         phip = ATAN2(yp,xp)
         IF (phip < 0) phip = phip+pi2
         CALL get_equil_s(rp,phip,zp,s,ier)
         IF (s <= 1) THEN
            x2 = xp+dx
            y2 = yp+dy
            z2 = zp+dz
            dx = (x2-x1)/lintsteps
            dy = (y2-y1)/lintsteps
            dz = (z2-z1)/lintsteps
            EXIT
         END IF
      END DO
      IF (i== lintsteps) THEN
         val = 0
         IF (PRESENT(length)) length = 0
         RETURN
      END IF
      val = 0
      delt = one/DBLE(nop-1)
      int_fac = one/DBLE(int_step)
      IF (PRESENT(length)) length = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1)) 
      DO i = 1, lintsteps
         DO j = 1, int_step
            xp = x1+dx*(i-1)+(j-1)*int_fac*dx
            yp = y1+dy*(i-1)+(j-1)*int_fac*dy
            zp = z1+dz*(i-1)+(j-1)*int_fac*dz
            delf = 0
            dx2 = x1+dx*(i-1)+(j)*int_fac*dx - xp
            dy2 = y1+dy*(i-1)+(j)*int_fac*dy - yp
            dz2 = z1+dz*(i-1)+(j)*int_fac*dz - zp
            DO k = 1, nop
               rp = sqrt(xp*xp+yp*yp)
               phip = ATAN2(yp,xp)
               IF (phip < 0) phip = phip+pi2
               CALL get_equil_s(rp,phip,zp,s,ier,uval)
               IF (s <= 1 .and. s >= 0) THEN
                  ier = 0
                  CALL fcn(s,uval,phip,dx2,dy2,dz2,f_val,ier)
                  IF (ier /= 0) f_val = 0
               ELSE
                  f_val = 0
               END IF
               !WRITE(327,*) xp,yp,zp,s,f_val,ier
               delf = delf + ci(k)*f_val
               xp   = xp + k*dx2*delt
               yp   = yp + k*dy2*delt
               zp   = zp + k*dz2*delt
            END DO
            val = val + delf
         END DO
      END DO
      RETURN
      END SUBROUTINE line_int_dbl
      
      SUBROUTINE line_int_sgl(fcn,r1,r2,val,length)
      IMPLICIT NONE
      REAL, INTENT(in)   :: r1(3), r2(3)
      REAL, INTENT(out)  :: val
      REAL, INTENT(out), OPTIONAL :: length
      DOUBLE PRECISION   :: r1_dbl(3), r2_dbl(3)
      DOUBLE PRECISION   :: val_dbl
      DOUBLE PRECISION   :: length_dbl
      EXTERNAL fcn
      r1_dbl = r1; r2_dbl = r2
      CALL line_int_dbl(fcn,r1_dbl,r2_dbl,val_dbl,length_dbl)
      val = val_dbl
      IF (PRESENT(length)) length = length_dbl
      RETURN
      END SUBROUTINE line_int_sgl

      !-----------------------------------------------------------------
      SUBROUTINE mntouv(k1,k,mnmax,nu,nv,xu,xv,fmn,xm,xn,f,signs,calc_trig)
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
         FORALL(mn=1:mnmax,i=1:nu) cosmt(mn,i) = dcos(pi2*mt(mn,i))
         FORALL(mn=1:mnmax,i=1:nu) sinmt(mn,i) = dsin(pi2*mt(mn,i))
         FORALL(mn=1:mnmax,i=1:nv) cosnz(mn,i) = dcos(pi2*nz(mn,i))
         FORALL(mn=1:mnmax,i=1:nv) sinnz(mn,i) = dsin(pi2*nz(mn,i))
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
      END SUBROUTINE mntouv

      LOGICAL*2 FUNCTION isingrid(x1_in,x2_in,x3_in)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: x1_in
      DOUBLE PRECISION, INTENT(in), OPTIONAL :: x2_in, x3_in
      isingrid = .FALSE.
      isingrid = ((x1_in >= x1_min-eps1) .and. (x1_in <= x1_max+eps1))
      IF (PRESENT(x2_in)) isingrid = (isingrid .and. (x2_in >= x2_min-eps2) .and. (x2_in <= x2_max+eps2))
      IF (PRESENT(x3_in)) isingrid = (isingrid .and. (x3_in >= x3_min-eps3) .and. (x3_in <= x3_max+eps3))
      RETURN
      END FUNCTION isingrid

      SUBROUTINE lookupgrid1d(x3_in,k,hz,hzi,zparam)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: x3_in
      INTEGER, INTENT(out) :: k
      DOUBLE PRECISION, INTENT(out) :: hz, hzi, zparam
      k = MIN(MAX(COUNT(x3 < x3_in),1),nx3-1)
      hz     = x3(k+1) - x3(k)
      hzi    = 1 / hz
      zparam = (x3_in - x3(k)) * hzi
      RETURN
      END SUBROUTINE lookupgrid1d

      SUBROUTINE lookupgrid3d(x1_in,x2_in,x3_in,i,j,k,hx,hy,hz,hxi,hyi,hzi,xparam,yparam,zparam)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: x1_in, x2_in, x3_in
      INTEGER, INTENT(out) :: i,j,k
      DOUBLE PRECISION, INTENT(out) :: hx, hy, hz, hxi, hyi, hzi, xparam, yparam, zparam
      i = MIN(MAX(COUNT(x1 < x1_in),1),nx1-1)
      j = MIN(MAX(COUNT(x2 < x2_in),1),nx2-1)
      k = MIN(MAX(COUNT(x3 < x3_in),1),nx3-1)
      hx     = x1(i+1) - x1(i)
      hy     = x2(j+1) - x2(j)
      hz     = x3(k+1) - x3(k)
      hxi    = 1 / hx
      hyi    = 1 / hy
      hzi    = 1 / hz
      xparam = (x1_in - x1(i)) * hxi
      yparam = (x2_in - x2(j)) * hyi
      zparam = (x3_in - x3(k)) * hzi
      RETURN
      END SUBROUTINE lookupgrid3d
      !-----------------------------------------------------------------

!-----------------------------------------------------------------------
!     END MODULE
!-----------------------------------------------------------------------
      END MODULE stel_tools
