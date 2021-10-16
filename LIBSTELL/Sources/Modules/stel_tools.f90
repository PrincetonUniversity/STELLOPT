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
                                    BV4D, B4D, L4D, LU4D, LV4D, &
                                    B_R4D, B_Z4D, B_Phi4D
      INTEGER, PRIVATE :: win_B_R4D, win_B_Z4D, win_B_Phi4D
      INTEGER, PRIVATE :: win_R4D, win_Z4D, win_G4D, win_RU4D, win_ZU4D, win_RV4D, &
                          win_ZV4D, win_BS4D, win_BU4D, win_BV4D, win_B4D, &   
                          win_L4D, win_LU4D, win_LV4D, win_VP2D, win_GRHO2D, win_GRHO22D, &
                          win_BAV2D, win_BSQ2D, win_S112D, win_S122D, win_S212D, win_S222D, &
                          win_x1, win_x2, win_x3, nx1, nx2, nx3
      REAL*8, PRIVATE :: eps1, eps2, eps3, x1_min, x1_max, x2_min, x2_max, x3_min, x3_max
      REAL*8, parameter, PRIVATE :: small = 1.e-10_ezspline_r8
      LOGICAL, PARAMETER, PRIVATE :: lcalc_deriv = .false.
!-----------------------------------------------------------------------
!     Private Subroutines
!-----------------------------------------------------------------------
      PRIVATE :: mntouv, isingrid, lookupgrid1d, lookupgrid3d, mycross
!-----------------------------------------------------------------------
!     INTERFACE Modules
!-----------------------------------------------------------------------
      INTERFACE load_fourier_geom
         MODULE PROCEDURE load_fourier_geom_dbl, load_fourier_geom_sgl
      END INTERFACE
      INTERFACE load_vmec_geom
         MODULE PROCEDURE load_vmec_geom_dbl, load_vmec_geom_sgl
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
      INTERFACE get_equil_Bcylsuv2
         MODULE PROCEDURE get_equil_Bcylsuv2_dbl, get_equil_Bcylsuv2_sgl
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
      INTERFACE get_equil_kappa2
         MODULE PROCEDURE get_equil_kappa2_dbl, get_equil_kappa2_sgl
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
      INTEGER ::  ns_t, u, mn, isherm, nu1, nv1
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
      TYPE(EZspline3_r8) :: B_R_spl, B_Z_spl, B_Phi_spl
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
      nu1 = nu-1
      nv1 = nv-1
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
         IF (ASSOCIATED(B_R4D)) CALL mpidealloc(B_R4D, win_B_R4D)
         IF (ASSOCIATED(B_Z4D)) CALL mpidealloc(B_Z4D, win_B_Z4D)
         IF (ASSOCIATED(B_Phi4D)) CALL mpidealloc(B_Phi4D, win_B_Phi4D)
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
         CALL mpialloc(B_R4D, 8, nu, nv, ns_t, shar_rank, 0, shar_comm, win_B_R4D)
         CALL mpialloc(B_Z4D, 8, nu, nv, ns_t, shar_rank, 0, shar_comm, win_B_Z4D)
         CALL mpialloc(B_Phi4D, 8, nu, nv, ns_t, shar_rank, 0, shar_comm, win_B_Phi4D)
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
         IF (ASSOCIATED(B_R4D))  DEALLOCATE(B_R4D)
         IF (ASSOCIATED(B_Z4D))  DEALLOCATE(B_Z4D)
         IF (ASSOCIATED(B_Phi4D))  DEALLOCATE(B_Phi4D)
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
         ALLOCATE(B_R4D(8,nu,nv,ns_t))
         ALLOCATE(B_Z4D(8,nu,nv,ns_t))
         ALLOCATE(B_Phi4D(8,nu,nv,ns_t))
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
         CALL EZspline_init(B_R_spl,nu,nv,ns_t,bcs1,bcs1,bcs0,iflag)
         CALL EZspline_init(B_Z_spl,nu,nv,ns_t,bcs1,bcs1,bcs0,iflag)
         CALL EZspline_init(B_Phi_spl,nu,nv,ns_t,bcs1,bcs1,bcs0,iflag)
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

         ! B_R, B_Z, B_Phi
         B_R_spl%x1 = xu*pi2; B_R_spl%x2 = xv*pi2; B_R_spl%x3 = rho; B_R_spl%isHermite = isherm
         !br = RU_SPL%fspl * BU_SPL%fspl + RV_SPL%fspl * BV_SPL%fspl *nfp
         ! this part is 0: BS_SPL%fspl * RS_SPL%fspl
         f_temp = RU_SPL%fspl(1,:,:,:) * BU_SPL%fspl(1,:,:,:) + RV_SPL%fspl(1,:,:,:) * BV_SPL%fspl(1,:,:,:) *nfp 
         CALL EZspline_setup(B_R_spl,f_temp,iflag); f_temp = 0

         B_Z_spl%x1 = xu*pi2; B_Z_spl%x2 = xv*pi2; B_Z_spl%x3 = rho; B_Z_spl%isHermite = isherm
         !bz = Z_grad(3)*Bs + Z_grad(1)*Bu + Z_grad(2)*Bv*nfp
         f_temp = ZU_SPL%fspl(1,:,:,:) * BU_SPL%fspl(1,:,:,:) + ZV_SPL%fspl(1,:,:,:) * BV_SPL%fspl(1,:,:,:)*nfp
         CALL EZspline_setup(B_Z_spl,f_temp,iflag); f_temp = 0

         B_Phi_spl%x1 = xu*pi2; B_Phi_spl%x2 = xv*pi2; B_Phi_spl%x3 = rho; B_Phi_spl%isHermite = isherm
         !bphi = r_val * Bv
         f_temp = R_SPL%fspl(1,:,:,:) * BV_SPL%fspl(1,:,:,:)
         CALL EZspline_setup(B_Phi_spl,f_temp,iflag); f_temp = 0


         ! Now we can get rid of some stuff
         x1   = R_SPL%x1
         x2   = R_SPL%x2
         x3   = R_SPL%x3
         R4D  = R_SPL%fspl
         Z4D  = Z_SPL%fspl
         RU4D = RU_SPL%fspl
         ZU4D = ZU_SPL%fspl
         RV4D = RV_SPL%fspl
         ZV4D = ZV_SPL%fspl
         L4D  = L_SPL%fspl
         LU4D  = LU_SPL%fspl
         LV4D  = LV_SPL%fspl
         G4D  = G_SPL%fspl
         BS4D = BS_SPL%fspl
         BU4D = BU_SPL%fspl
         BV4D = BV_SPL%fspl
         B4D  = B_SPL%fspl
         B_R4D  = B_R_SPL%fspl
         B_Z4D  = B_Z_SPL%fspl
         B_Phi4D  = B_Phi_SPL%fspl
         CALL EZspline_free(R_spl,iflag)
         CALL EZspline_free(Z_spl,iflag)
         CALL EZspline_free(RU_spl,iflag)
         CALL EZspline_free(RV_spl,iflag)
         CALL EZspline_free(ZU_spl,iflag)
         CALL EZspline_free(ZV_spl,iflag)
         CALL EZspline_free(L_spl,iflag)
         CALL EZspline_free(LU_spl,iflag)
         CALL EZspline_free(LV_spl,iflag)
         CALL EZspline_free(G_spl,iflag)
         CALL EZspline_free(BS_spl,iflag)
         CALL EZspline_free(BU_spl,iflag)
         CALL EZspline_free(BV_spl,iflag)
         CALL EZspline_free(B_spl,iflag)
         CALL EZspline_free(B_R_spl,iflag)
         CALL EZspline_free(B_Z_spl,iflag)
         CALL EZspline_free(B_Phi_spl,iflag)
         
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
            grho   = SUM(SUM(f_temp(1:nu1,1:nv1,:),DIM=1),DIM=1)/(nu1*nv1)
            CALL EZspline_setup(S11_spl,grho,iflag); f_temp = 0; grho = 0
            ! Calc S21
            f_temp = (RU4D(1,:,:,:)*RV4D(1,:,:,:)+ &
                      ZU4D(1,:,:,:)*ZV4D(1,:,:,:))*nfp
            f_temp = f_temp / G4D(1,:,:,:)
            grho   = SUM(SUM(f_temp(1:nu1,1:nv1,:),DIM=1),DIM=1)/(nu1*nv1)
            CALL EZspline_setup(S21_spl,grho,iflag); f_temp = 0; grho = 0
            ! Calc S12
            f_temp = (RU4D(1,:,:,:)*RV4D(1,:,:,:)+ &
                      ZU4D(1,:,:,:)*ZV4D(1,:,:,:))* &
                     (one+LU4D(1,:,:,:))*nfp
            f_temp = f_temp - (RU4D(1,:,:,:)*RU4D(1,:,:,:)+ &
                               ZU4D(1,:,:,:)*ZU4D(1,:,:,:))*&
                              LV4D(1,:,:,:)*nfp
            f_temp = f_temp / G4D(1,:,:,:)
            grho   = SUM(SUM(f_temp(1:nu1,1:nv1,:),DIM=1),DIM=1)/(nu1*nv1)
            CALL EZspline_setup(S12_spl,grho,iflag); f_temp = 0; grho = 0
            ! Calc S22
            f_temp = (RV4D(1,:,:,:)*RV4D(1,:,:,:)*nfp*nfp+ &
                      ZV4D(1,:,:,:)*ZV4D(1,:,:,:)*nfp*nfp+ &
                      R4D(1,:,:,:)* R4D(1,:,:,:))* &
                     (one+LU4D(1,:,:,:))
            f_temp = f_temp - (RU4D(1,:,:,:)*RV4D(1,:,:,:)+ &
                               ZU4D(1,:,:,:)*ZV4D(1,:,:,:))*&
                              LV4D(1,:,:,:)*nfp*nfp
            f_temp = f_temp / G4D(1,:,:,:)
            grho   = SUM(SUM(f_temp(1:nu1,1:nv1,:),DIM=1),DIM=1)/(nu1*nv1)
            CALL EZspline_setup(S22_spl,grho,iflag); f_temp = 0; grho = 0
            ! Bav
            f_temp = B4D(1,:,:,:)*G4D(1,:,:,:)
            grho   = SUM(SUM(f_temp,DIM=1),DIM=1)
            grho2 = grho / Vp
            grho2(1) = 2*grho2(2) - grho2(3)
            CALL EZspline_setup(Bav_spl,grho2,iflag); f_temp = 0; grho = 0
            ! Bsq
            f_temp = B4D(1,:,:,:)*B4D(1,:,:,:)*G4D(1,:,:,:)
            grho   = SUM(SUM(f_temp,DIM=1),DIM=1)
            grho2 = grho / Vp
            grho2(1) = 2*grho2(2) - grho2(3)
            CALL EZspline_setup(Bsq_spl,grho2,iflag); f_temp = 0; grho = 0
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
      END IF !So shared memory doesnt do work
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

      SUBROUTINE load_vmec_geom_dbl(ns,mnmax,nu,nv,xm,xn_in,iflag,rmnc,zmns,lmns,&
                                    phiprime,iota,rmns,zmnc,lmnc,comm)
      ! Couple of notes here
      ! Lambda, phiprime, and iota are on the half mesh coming into this
      USE EZspline
      USE mpi_sharmem
#if defined(MPI_OPT)
      USE mpi
#endif
      IMPLICIT NONE
      INTEGER, INTENT(in)        :: ns
      INTEGER, INTENT(in)        :: mnmax
      INTEGER, INTENT(in)        :: nu
      INTEGER, INTENT(in)        :: nv
      INTEGER, INTENT(in) :: xm(1:mnmax)
      INTEGER, INTENT(in) :: xn_in(1:mnmax)
      INTEGER, INTENT(inout) :: iflag
      DOUBLE PRECISION, INTENT(in) :: rmnc(1:mnmax,1:ns), zmns(1:mnmax,1:ns), lmns(1:mnmax,1:ns)
      DOUBLE PRECISION, INTENT(in) :: iota(1:ns),phiprime(1:ns)
      DOUBLE PRECISION, INTENT(in),OPTIONAL :: rmns(1:mnmax,1:ns), zmnc(1:mnmax,1:ns), lmnc(1:mnmax,1:ns)
      INTEGER, INTENT(in), OPTIONAL :: comm
      INTEGER ::  ns_t, u, mn, isherm, nu1, nv1, k1p, k1, k2
      INTEGER ::  shar_comm, shar_rank, shar_size
      INTEGER ::  xn(1:mnmax)
      DOUBLE PRECISION :: ohs
      DOUBLE PRECISION, ALLOCATABLE :: xu(:),xv(:),rho(:),vp(:),grho(:),grho2(:),rhoinv(:)
      DOUBLE PRECISION, ALLOCATABLE :: fmn_temp(:,:), fmn_o(:,:),fmn_e(:,:), fumn_o(:,:), fumn_e(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: f_temp(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: r_e(:,:,:), r_o(:,:,:), z_e(:,:,:), z_o(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: ru_e(:,:,:), ru_o(:,:,:), zu_e(:,:,:), zu_o(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: rs(:,:,:), zs(:,:,:), ru12(:,:,:), zu12(:,:,:)
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
      k1 = 1
      k2 = ns
      k1p= 2
      ns_t = ns
      isherm = 0  ! Cannot change now
      ! Preform checks
      IF (mnmax< 1) iflag = -3
      IF (nu < 1 .or. nv < 1) iflag = -4
      IF (PRESENT(rmns).NEQV.PRESENT(zmnc)) iflag = -5
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
      nu1 = nu-1
      nv1 = nv-1
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
#if defined(MPI_OPT)
      END IF
#endif
      IF (shar_rank == 0) THEN
         !Allocations
         ALLOCATE(xu(nu),xv(nv),rho(k1:k2),rhoinv(k1:k2))
         ALLOCATE(fmn_temp(1:mnmax,k1:k2))
         ALLOCATE(f_temp(nu,nv,k1:k2))
         FORALL(u=k1:k2) rho(u) = REAL(u-1)/REAL(ns_t-1)
         rho = SQRT(rho) ! Improves lookup near axis
         ohs = k2-k1
         rhoinv = 1.0/rho
         rhoinv(1) = 1.0
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

         ! Define Even Odd quantities
         ALLOCATE(fmn_e(1:mnmax,k1:k2),fmn_o(1:mnmax,k1:k2))
         ALLOCATE(fumn_e(1:mnmax,k1:k2),fumn_o(1:mnmax,k1:k2))
         ALLOCATE(r_e(nu,nv,k1:k2),r_o(nu,nv,k1:k2),z_e(nu,nv,k1:k2),z_o(nu,nv,k1:k2))
         ALLOCATE(ru_e(nu,nv,k1:k2),ru_o(nu,nv,k1:k2),zu_e(nu,nv,k1:k2),zu_o(nu,nv,k1:k2))
         ALLOCATE(rs(nu,nv,k1:k2),zs(nu,nv,k1:k2),ru12(nu,nv,k1:k2),zu12(nu,nv,k1:k2))
         DO mn = 1, mnmax
            fmn_e(mn,:) = 0; fmn_o(mn,:) = 0
            fumn_e(mn,:) = 0; fumn_o(mn,:) = 0
            IF (MOD(xm(mn),2)==1) THEN
               fmn_o(mn,:)  =         rmnc(mn,:)*rhoinv
               fumn_o(mn,:) = -xm(mn)*rmnc(mn,:)*rhoinv
            ELSE
               fmn_e(mn,:)  =         rmnc(mn,:)
               fumn_e(mn,:) = -xm(mn)*rmnc(mn,:)
            END IF
         END DO
         CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fmn_e,xm,xn,r_e,0,1)
         CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fmn_o,xm,xn,r_o,0,0)
         CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fumn_e,xm,xn,ru_e,1,0)
         CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fumn_o,xm,xn,ru_o,1,0)
         IF (PRESENT(rmns)) THEN
            DO mn = 1, mnmax
               fmn_e(mn,:) = 0; fmn_o(mn,:) = 0
               fumn_e(mn,:) = 0; fumn_o(mn,:) = 0
               IF (MOD(xm(mn),2)==1) THEN
                  fmn_o(mn,:)  =         rmns(mn,:)*rhoinv
                  fumn_o(mn,:) =  xm(mn)*rmns(mn,:)*rhoinv
               ELSE
                  fmn_e(mn,:)  =         rmns(mn,:)
                  fumn_e(mn,:) =  xm(mn)*rmns(mn,:)
               END IF
            END DO
            CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fmn_e,xm,xn,r_e,1,0)
            CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fmn_o,xm,xn,r_o,1,0)
            CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fumn_e,xm,xn,ru_e,0,0)
            CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fumn_o,xm,xn,ru_o,0,0)
         END IF
         DO mn = 1, mnmax
            fmn_e(mn,:) = 0; fmn_o(mn,:) = 0
            fumn_e(mn,:) = 0; fumn_o(mn,:) = 0
            IF (MOD(xm(mn),2)==1) THEN
               fmn_o(mn,:)  =         zmns(mn,:)*rhoinv
               fumn_o(mn,:) =  xm(mn)*zmns(mn,:)*rhoinv
            ELSE
               fmn_e(mn,:)  =         zmns(mn,:)
               fumn_e(mn,:) =  xm(mn)*zmns(mn,:)
            END IF
         END DO
         CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fmn_e,xm,xn,z_e,1,0)
         CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fmn_o,xm,xn,z_o,1,0)
         CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fumn_e,xm,xn,zu_e,0,0)
         CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fumn_o,xm,xn,zu_o,0,0)
         IF (PRESENT(zmnc)) THEN
            DO mn = 1, mnmax
               fmn_e(mn,:) = 0; fmn_o(mn,:) = 0
               fumn_e(mn,:) = 0; fumn_o(mn,:) = 0
               IF (MOD(xm(mn),2)==1) THEN
                  fmn_o(mn,:)  =         zmnc(mn,:)*rhoinv
                  fumn_o(mn,:) = -xm(mn)*zmnc(mn,:)*rhoinv
               ELSE
                  fmn_e(mn,:)  =         zmnc(mn,:)
                  fumn_e(mn,:) = -xm(mn)*zmnc(mn,:)
               END IF
            END DO
            CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fmn_e,xm,xn,z_e,0,0)
            CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fmn_o,xm,xn,z_o,0,0)
            CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fumn_e,xm,xn,zu_e,1,0)
            CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fumn_o,xm,xn,zu_o,1,0)
         END IF
         rs = 0; zs = 0;
         DO mn = k1p,k2
            rs(:,:,mn) = ohs*(r_e(:,:,mn)-r_e(:,:,mn-1) &
                         + rho(mn)*(r_o(:,:,mn)-r_o(:,:,mn-1)))
            zs(:,:,mn) = ohs*(z_e(:,:,mn)-z_e(:,:,mn-1) &
                         + rho(mn)*(z_o(:,:,mn)-z_o(:,:,mn-1)))
         END DO
         rs(:,:,k1) = rs(:,:,k1p)
         zs(:,:,k1) = zs(:,:,k1p)
         DEALLOCATE(fmn_o,fmn_e)
         DEALLOCATE(fumn_o,fumn_e)

         ! R
         f_temp = 0;
         R_spl%x1 = xu*pi2; R_spl%x2 = xv*pi2; R_spl%x3 = rho; R_spl%isHermite = isherm
         FORALL(mn = k1:k2) f_temp(:,:,mn) = r_e(:,:,mn) + rho(mn)*r_o(:,:,mn) 
         CALL EZspline_setup(R_spl,f_temp,iflag); f_temp=0
         ! Z
         Z_spl%x1 = xu*pi2; Z_spl%x2 = xv*pi2; Z_spl%x3 = rho; Z_spl%isHermite = isherm
         FORALL(mn = k1:k2) f_temp(:,:,mn) = z_e(:,:,mn) + rho(mn)*z_o(:,:,mn) 
         CALL EZspline_setup(Z_spl,f_temp,iflag); f_temp=0
         ! dR/du
         Ru_spl%x1 = xu*pi2; Ru_spl%x2 = xv*pi2; Ru_spl%x3 = rho; Ru_spl%isHermite = isherm
         FORALL(mn = k1:k2) f_temp(:,:,mn) = ru_e(:,:,mn) + rho(mn)*ru_o(:,:,mn) 
         ru12 = 0;
         DO mn = k1p,k2
            ru12(:,:,mn) = (f_temp(:,:,mn)+f_temp(:,:,mn-1))*0.5
         END DO
         ru12(:,:,k1)=ru12(:,:,k1p)
         CALL EZspline_setup(Ru_spl,f_temp,iflag); f_temp = 0
         ! dZ/du
         Zu_spl%x1 = xu*pi2; Zu_spl%x2 = xv*pi2; Zu_spl%x3 = rho; Zu_spl%isHermite = isherm
         FORALL(mn = k1:k2) f_temp(:,:,mn) = zu_e(:,:,mn) + rho(mn)*zu_o(:,:,mn) 
         zu12 = 0;
         DO mn = k1p,k2
            zu12(:,:,mn) = (f_temp(:,:,mn)+f_temp(:,:,mn-1))*0.5
         END DO
         zu12(:,:,k1)=zu12(:,:,k1p)
         CALL EZspline_setup(Zu_spl,f_temp,iflag); f_temp = 0
         ! dR/Dv Derivatives
         Rv_spl%x1 = xu*pi2; Rv_spl%x2 = xv*pi2; Rv_spl%x3 = rho; Rv_spl%isHermite = isherm
         FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -rmnc(mn,:)*xn(mn)
         CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,f_temp,1,0)
         IF (PRESENT(rmns)) THEN
            FORALL(mn = 1:mnmax) fmn_temp(mn,:) = rmns(mn,:)*xn(mn)
            CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,f_temp,0,0)
         END IF
         CALL EZspline_setup(Rv_spl,f_temp,iflag); f_temp = 0
         ! dZ/Dv Derivatives
         Zv_spl%x1 = xu*pi2; Zv_spl%x2 = xv*pi2; Zv_spl%x3 = rho; Zv_spl%isHermite = isherm
         FORALL(mn = 1:mnmax) fmn_temp(mn,:) = zmns(mn,:)*xn(mn)
         CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,f_temp,0,0)
         IF (PRESENT(zmnc)) THEN
            FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -zmnc(mn,:)*xn(mn)
            CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,f_temp,1,0)
         END IF
         CALL EZspline_setup(Zv_spl,f_temp,iflag); f_temp = 0
         ! Lambda (on half grid)
         L_spl%x1 = xu*pi2; L_spl%x2 = xv*pi2; L_spl%x3 = rho; L_spl%isHermite = isherm
         CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,lmns,xm,xn,f_temp,1,0)
         IF (PRESENT(lmnc)) CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,lmnc,xm,xn,f_temp,0,0)
         CALL EZspline_setup(L_spl,f_temp,iflag); f_temp = 0
         ! Lambda/u
         Lu_spl%x1 = xu*pi2; Lu_spl%x2 = xv*pi2; Lu_spl%x3 = rho; Lu_spl%isHermite = isherm
         FORALL(mn = 1:mnmax) fmn_temp(mn,:) = lmns(mn,:)*xm(mn)
         CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,f_temp,0,0)
         IF (PRESENT(lmnc)) THEN
            FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -lmnc(mn,:)*xm(mn)
            CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,f_temp,1,0)
         END IF
         CALL EZspline_setup(Lu_spl,f_temp,iflag); f_temp = 0
         ! Lambda/v
         Lv_spl%x1 = xu*pi2; Lv_spl%x2 = xv*pi2; Lv_spl%x3 = rho; Lv_spl%isHermite = isherm
         FORALL(mn = 1:mnmax) fmn_temp(mn,:) = lmns(mn,:)*xn(mn)
         CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,f_temp,0,0)
         IF (PRESENT(lmnc)) THEN
            FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -lmnc(mn,:)*xn(mn)
            CALL mntouv(k1,k2,mnmax,nu,nv,xu,xv,fmn_temp,xm,xn,f_temp,1,0)
         END IF
         CALL EZspline_setup(Lv_spl,f_temp,iflag); f_temp = 0

         DEALLOCATE(fmn_temp)

         ! Do this here for easy access 
         x1   = R_SPL%x1
         x2   = R_SPL%x2
         x3   = R_SPL%x3
         R4D  = R_SPL%fspl
         Z4D  = Z_SPL%fspl
         RU4D = RU_SPL%fspl
         ZU4D = ZU_SPL%fspl
         RV4D = RV_SPL%fspl
         ZV4D = ZV_SPL%fspl
         L4D  = L_SPL%fspl
         LU4D  = LU_SPL%fspl
         LV4D  = LV_SPL%fspl

         ! Calc Gsqrt
         ! SQRT(G) = R(RuZs-RsZu) Eq17 Hirshman 83
         G_spl%x1 = xu*pi2; G_spl%x2 = xv*pi2; G_spl%x3 = rho; G_spl%isHermite = isherm; f_temp=0
         rhoinv(k1p:k2) = 2/(rho(1:k2-1)+rho(k1p:k2))
         DO mn = k1p, k2
            f_temp(:,:,mn) = ru12(:,:,mn)*zs(:,:,mn) &
                           + 0.25*(   ru_o(:,:,mn)*z_o(:,:,mn) + ru_o(:,:,mn-1)*z_o(:,:,mn-1) &
                                   + (ru_e(:,:,mn)*z_o(:,:,mn) + ru_e(:,:,mn-1)*z_o(:,:,mn-1))*rhoinv(mn))
            f_temp(:,:,mn) = f_temp(:,:,mn) - zu12(:,:,mn)*rs(:,:,mn)&
                           - 0.25*(   zu_o(:,:,mn)*r_o(:,:,mn) + zu_o(:,:,mn-1)*r_o(:,:,mn-1) &
                                   + (zu_e(:,:,mn)*r_o(:,:,mn) + zu_e(:,:,mn-1)*r_o(:,:,mn-1))*rhoinv(mn))
            f_temp(:,:,mn) = 0.5*(R4D(1,:,:,mn) + R4D(1,:,:,mn-1))*f_temp(:,:,mn)
         END DO
         f_temp(:,:,k1) = f_temp(:,:,k1p)

         ! Put on full grid
         ! From the half grid (2:ns) just average to (2:ns-1)
         f_temp(:,:,k1p:k2-1) = 0.5*(f_temp(:,:,k1p:k2-1)+f_temp(:,:,k1p+1:k2))
         f_temp(:,:,k2) = 2.0*f_temp(:,:,k2) - f_temp(:,:,k2-1) ! Extrapolate to edge
         !This works but is probably not the correct way to extrapolate near the axis
         f_temp(:,:,k1p) = 2*f_temp(:,:,k1p+1) - f_temp(:,:,k1p+2)
         f_temp(:,:,k1) = 2*f_temp(:,:,k1p) - f_temp(:,:,k1p+1)

         ! Create Spline
         CALL EZspline_setup(G_spl,f_temp,iflag); f_temp = 0

         ! B^s
         Bs_spl%x1 = xu*pi2; Bs_spl%x2 = xv*pi2; Bs_spl%x3 = rho; Bs_spl%isHermite = isherm
         CALL EZspline_setup(Bs_spl,f_temp,iflag); f_temp = 0

         ! B^u = phip*(iota-Lv)/sqrt(g)
         Bu_spl%x1 = xu*pi2; Bu_spl%x2 = xv*pi2; Bu_spl%x3 = rho; Bu_spl%isHermite = isherm
         f_temp = -LV4D(1,:,:,:)*nfp
         FORALL(u=k1:k2) f_temp(:,:,u) = -(f_temp(:,:,u)+iota(u))*phiprime(u)/pi2
         f_temp = f_temp / G_SPL%fspl(1,:,:,:)
         CALL EZspline_setup(Bu_spl,f_temp,iflag); f_temp = 0

         ! B^v = phip*(1+Lu)/sqrt(g)
         Bv_spl%x1 = xu*pi2; Bv_spl%x2 = xv*pi2; Bv_spl%x3 = rho; Bv_spl%isHermite = isherm
         f_temp =  LU4D(1,:,:,:)+1
         FORALL(u=k1:k2) f_temp(:,:,u) = -f_temp(:,:,u)*phiprime(u)/pi2
         f_temp = f_temp / G_SPL%fspl(1,:,:,:)
         CALL EZspline_setup(Bv_spl,f_temp,iflag); f_temp = 0

         ! |B|^2 = Bu**2*guu+2*Bu*Bv*guv+Bv**2*gvv Eq8b Hirshman 83 (Bk=B^k)
         !  guu = Ru*Ru+Zv*Zv (Ru = dR/du)
         !  guv = Ru*Rv+Zu*Zv
         !  gvv = Rv*Rv+R**2+Zv*Zv
         B_spl%x1 = xu*pi2; B_spl%x2 = xv*pi2; B_spl%x3 = rho; B_spl%isHermite = isherm
         f_temp = (RU4D(1,:,:,:)*RU4D(1,:,:,:)+ZU4D(1,:,:,:)*ZU4D(1,:,:,:))*BU_SPL%fspl(1,:,:,:)*BU_SPL%fspl(1,:,:,:)
         f_temp = f_temp + (RV4D(1,:,:,:)*RV4D(1,:,:,:)+ZV4D(1,:,:,:)*ZV4D(1,:,:,:))*BV_SPL%fspl(1,:,:,:)*BV_SPL%fspl(1,:,:,:)*nfp*nfp
         f_temp = f_temp + (RU4D(1,:,:,:)*RV4D(1,:,:,:)*nfp+ZU4D(1,:,:,:)*ZV4D(1,:,:,:)*nfp+R4D(1,:,:,:)*R4D(1,:,:,:))*BU_SPL%fspl(1,:,:,:)*BV_SPL%fspl(1,:,:,:)
         f_temp = SQRT(f_temp)
         CALL EZspline_setup(B_spl,f_temp,iflag); f_temp = 0
         DEALLOCATE(Rs,Zs,ru12,zu12)
         DEALLOCATE(ru_e,ru_o,zu_e,zu_o)
         DEALLOCATE(r_e,r_o,z_e,z_o)

         ! Now we can get rid of some stuff

         ! Here the part where we copy and delete everything
         G4D  = G_SPL%fspl
         BS4D = BS_SPL%fspl
         BU4D = BU_SPL%fspl
         BV4D = BV_SPL%fspl
         B4D  = B_SPL%fspl
         CALL EZspline_free(R_spl,iflag)
         CALL EZspline_free(Z_spl,iflag)
         CALL EZspline_free(RU_spl,iflag)
         CALL EZspline_free(RV_spl,iflag)
         CALL EZspline_free(ZU_spl,iflag)
         CALL EZspline_free(ZV_spl,iflag)
         CALL EZspline_free(L_spl,iflag)
         CALL EZspline_free(LU_spl,iflag)
         CALL EZspline_free(LV_spl,iflag)
         CALL EZspline_free(G_spl,iflag)
         CALL EZspline_free(BS_spl,iflag)
         CALL EZspline_free(BU_spl,iflag)
         CALL EZspline_free(BV_spl,iflag)
         CALL EZspline_free(B_spl,iflag)
         
         !
         ! Now we just go ahead and calculate the surface averaged quantities
         !

         ! Deallocate if allocated
         IF (EZspline_allocated(Vp_spl)) CALL EZspline_free(Vp_spl,iflag)
         IF (EZspline_allocated(grho_spl)) CALL EZspline_free(grho_spl,iflag)
         IF (EZspline_allocated(grho2_spl)) CALL EZspline_free(grho2_spl,iflag)
         IF (EZspline_allocated(S11_spl)) CALL EZspline_free(S11_spl,iflag)
         IF (EZspline_allocated(S12_spl)) CALL EZspline_free(S12_spl,iflag)
         IF (EZspline_allocated(S21_spl)) CALL EZspline_free(S21_spl,iflag)
         IF (EZspline_allocated(S22_spl)) CALL EZspline_free(S22_spl,iflag)
         IF (EZspline_allocated(Bav_spl)) CALL EZspline_free(Bav_spl,iflag)
         IF (EZspline_allocated(Bsq_spl)) CALL EZspline_free(Bsq_spl,iflag)

         ! Initialize Splines
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
         grho   = SUM(SUM(f_temp(1:nu1,1:nv1,:),DIM=1),DIM=1)/(nu1*nv1)
         CALL EZspline_setup(S11_spl,grho,iflag); f_temp = 0; grho = 0
         ! Calc S21
         f_temp = (RU4D(1,:,:,:)*RV4D(1,:,:,:)+ &
                   ZU4D(1,:,:,:)*ZV4D(1,:,:,:))*nfp
         f_temp = f_temp / G4D(1,:,:,:)
         grho   = SUM(SUM(f_temp(1:nu1,1:nv1,:),DIM=1),DIM=1)/(nu1*nv1)
         CALL EZspline_setup(S21_spl,grho,iflag); f_temp = 0; grho = 0
         ! Calc S12
         f_temp = (RU4D(1,:,:,:)*RV4D(1,:,:,:)+ &
                   ZU4D(1,:,:,:)*ZV4D(1,:,:,:))* &
                  (one+LU4D(1,:,:,:))*nfp
         f_temp = f_temp - (RU4D(1,:,:,:)*RU4D(1,:,:,:)+ &
                            ZU4D(1,:,:,:)*ZU4D(1,:,:,:))*&
                            LV4D(1,:,:,:)*nfp
         f_temp = f_temp / G4D(1,:,:,:)
         grho   = SUM(SUM(f_temp(1:nu1,1:nv1,:),DIM=1),DIM=1)/(nu1*nv1)
         CALL EZspline_setup(S12_spl,grho,iflag); f_temp = 0; grho = 0
         ! Calc S22
         f_temp = (RV4D(1,:,:,:)*RV4D(1,:,:,:)*nfp*nfp+ &
                   ZV4D(1,:,:,:)*ZV4D(1,:,:,:)*nfp*nfp+ &
                   R4D(1,:,:,:)* R4D(1,:,:,:))* &
                  (one+LU4D(1,:,:,:))
         f_temp = f_temp - (RU4D(1,:,:,:)*RV4D(1,:,:,:)+ &
                            ZU4D(1,:,:,:)*ZV4D(1,:,:,:))*&
                            LV4D(1,:,:,:)*nfp*nfp
         f_temp = f_temp / G4D(1,:,:,:)
         grho   = SUM(SUM(f_temp(1:nu1,1:nv1,:),DIM=1),DIM=1)/(nu1*nv1)
         CALL EZspline_setup(S22_spl,grho,iflag); f_temp = 0; grho = 0
         ! Bav
         f_temp = B4D(1,:,:,:)*G4D(1,:,:,:)
         grho   = SUM(SUM(f_temp,DIM=1),DIM=1)
         grho2 = grho / Vp
         grho2(1) = 2*grho2(2) - grho2(3)
         CALL EZspline_setup(Bav_spl,grho2,iflag); f_temp = 0; grho = 0
         ! Bsq
         f_temp = B4D(1,:,:,:)*B4D(1,:,:,:)*G4D(1,:,:,:)
         grho   = SUM(SUM(f_temp,DIM=1),DIM=1)
         grho2 = grho / Vp
         grho2(1) = 2*grho2(2) - grho2(3)
         CALL EZspline_setup(Bsq_spl,grho2,iflag); f_temp = 0; grho = 0
         
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
         
         ! DEALLOCATIONS
         DEALLOCATE(xu,xv,rho,rhoinv)
         DEALLOCATE(f_temp)
      END IF !So shared memory doesnt do work
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL MPI_BARRIER(shar_comm,iflag)
         CALL MPI_COMM_FREE(shar_comm,iflag)
      END IF
#endif
      RETURN
      END SUBROUTINE load_vmec_geom_dbl

      SUBROUTINE load_vmec_geom_sgl(ns,mnmax,nu,nv,xm,xn_in,iflag,rmnc,zmns,lmns,&
                                    phiprime,iota,rmns,zmnc,lmnc,comm)
      USE EZspline
      IMPLICIT NONE
      INTEGER, INTENT(in)        :: ns
      INTEGER, INTENT(in)        :: mnmax
      INTEGER, INTENT(in)        :: nu
      INTEGER, INTENT(in)        :: nv
      INTEGER, INTENT(in) :: xm(1:mnmax)
      INTEGER, INTENT(in) :: xn_in(1:mnmax)
      INTEGER, INTENT(inout) :: iflag
      REAL, INTENT(in) :: rmnc(1:mnmax,1:ns), zmns(1:mnmax,1:ns), lmns(1:mnmax,1:ns)
      REAL, INTENT(in) :: iota(1:ns),phiprime(1:ns)
      REAL, INTENT(in),OPTIONAL :: rmns(1:mnmax,1:ns), zmnc(1:mnmax,1:ns), lmnc(1:mnmax,1:ns)
      INTEGER, INTENT(in), OPTIONAL :: comm
      DOUBLE PRECISION :: rmnc_dbl(1:mnmax,1:ns), zmns_dbl(1:mnmax,1:ns), lmns_dbl(1:mnmax,1:ns)
      DOUBLE PRECISION :: rmns_dbl(1:mnmax,1:ns), zmnc_dbl(1:mnmax,1:ns), lmnc_dbl(1:mnmax,1:ns)
      DOUBLE PRECISION :: iota_dbl(1:ns),phiprime_dbl(1:ns)
      rmnc_dbl = rmnc
      zmns_dbl = zmns
      lmns_dbl = lmns
      phiprime_dbl = phiprime
      iota_dbl = iota
      rmns_dbl = 0; zmnc_dbl = 0; lmnc_dbl = 0
      IF (PRESENT(rmns)) rmns_dbl = rmns
      IF (PRESENT(zmnc)) zmnc_dbl = zmnc
      IF (PRESENT(lmnc)) lmnc_dbl = lmnc
      IF (PRESENT(comm)) THEN
         CALL load_vmec_geom_dbl(ns,mnmax,nu,nv,xm,xn_in,iflag,rmnc_dbl,zmns_dbl,lmns_dbl,&
            phiprime_dbl, iota_dbl, RMNS=rmns_dbl,ZMNC=zmnc_dbl,LMNC=lmnc_dbl, COMM=comm)
      ELSE
         CALL load_vmec_geom_dbl(ns,mnmax,nu,nv,xm,xn_in,iflag,rmnc_dbl,zmns_dbl,lmns_dbl,&
            phiprime_dbl, iota_dbl, RMNS=rmns_dbl,ZMNC=zmnc_dbl,LMNC=lmnc_dbl)
      END IF
      RETURN
      END SUBROUTINE load_vmec_geom_sgl
      
      SUBROUTINE rzfunct_stel_tool(m,n,x,fvec,fjac,ldfjac,iflag)
      USE EZspline
      IMPLICIT NONE
      INTEGER :: m,n,ldfjac,iflag, ier
      DOUBLE PRECISION :: x(n),fvec(m),fjac(ldfjac,n)
      DOUBLE PRECISION :: R_temp, Z_temp
      DOUBLE PRECISION :: R_grad(2), Z_grad(2)
      INTEGER :: i,j,k
      REAL*8 :: xparam, yparam, zparam, hx, hy, hz, hxi, hyi, hzi
      REAL*8 :: fval1(1,1), fval3(1,3), fval4(1,4)
      INTEGER, parameter :: ict1(10)=(/1,0,0,0,0,0,0,0,0,0/)
      INTEGER, parameter :: ict3(10)=(/1,1,0,1,0,0,0,0,0,0/)
      INTEGER, parameter :: ict4(10)=(/1,1,1,1,0,0,0,0,0,0/)
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
         CALL r8fvtricub(ict3, 1, 1, fval3, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         R4D(1,1,1,1), nx1, nx2, nx3)
         R_temp = fval3(1,1)
         R_grad(1) = fval3(1,2); R_grad(2) = fval3(1,3)
         CALL r8fvtricub(ict3, 1, 1, fval3, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         Z4D(1,1,1,1), nx1, nx2, nx3)
         Z_temp = fval3(1,1)
         Z_grad(1) = fval3(1,2); Z_grad(2) = fval3(1,3)
         IF (iflag == 1) THEN
            R_temp = R_temp + R_grad(2)*(x(1)-one)
            Z_temp = Z_temp + Z_grad(2)*(x(1)-one)
            fvec(1) = (R_temp - R_target)
            fvec(2) = (Z_temp - Z_target)
         ELSE IF (iflag == 2) THEN
            fjac(1,1) = R_grad(2) !dR/ds
            fjac(1,2) = R_grad(1) !dR/du
            fjac(2,1) = Z_grad(2) !dZ/ds
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
         CALL r8fvtricub(ict1, 1, 1, fval1(1,1), i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         R4D(1,1,1,1), nx1, nx2, nx3)
         R_temp = fval1(1,1)
         CALL r8fvtricub(ict1, 1, 1, fval1(1,1), i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         Z4D(1,1,1,1), nx1, nx2, nx3)
         Z_temp = fval1(1,1)
         fvec(1) = (R_temp - R_target)
         fvec(2) = (Z_temp - Z_target)
      ELSE IF (iflag == 2) THEN
         !CALL EZspline_gradient(R_spl,x(2),PHI_Target,x(1),R_grad,iflag)
         !CALL EZspline_gradient(Z_spl,x(2),PHI_Target,x(1),Z_grad,iflag)
         CALL lookupgrid3d(x(2),PHI_Target,x(1),i,j,k,hx,hy,hz,hxi,hyi,hzi,xparam,yparam,zparam)
         CALL r8fvtricub(ict3, 1, 1, fval3, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         R4D(1,1,1,1), nx1, nx2, nx3)
         R_temp = fval3(1,1)
         R_grad(1) = fval3(1,2); R_grad(2) = fval3(1,3)
         CALL r8fvtricub(ict3, 1, 1, fval3, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         Z4D(1,1,1,1), nx1, nx2, nx3)
         Z_temp = fval3(1,1)
         Z_grad(1) = fval3(1,2); Z_grad(2) = fval3(1,3)
         fvec(1) = (R_temp - R_target)
         fvec(2) = (Z_temp - Z_target)
         fjac(1,1) = R_grad(2) !dR/ds
         fjac(1,2) = R_grad(1) !dR/du
         fjac(2,1) = Z_grad(2) !dZ/ds
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
      REAL*8 :: fval(1,4)
      INTEGER, parameter :: ict1(10)=(/1,0,0,0,0,0,0,0,0,0/)
      INTEGER, parameter :: ict4(10)=(/1,1,1,1,0,0,0,0,0,0/)
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
         IF (PRESENT(R_grad)) THEN
            CALL r8fvtricub(ict4, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                            hx, hxi, hy, hyi, hz, hzi, &
                            R4D(1,1,1,1), nx1, nx2, nx3)
            R_grad(1) = fval(1,2); R_grad(2) = fval(1,3); R_grad(3) = fval(1,4)
         ELSE
            CALL r8fvtricub(ict1, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                            hx, hxi, hy, hyi, hz, hzi, &
                            R4D(1,1,1,1), nx1, nx2, nx3)
         END IF
         R_val = fval(1,1)
         IF (PRESENT(Z_grad)) THEN
            CALL r8fvtricub(ict4, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                            hx, hxi, hy, hyi, hz, hzi, &
                            Z4D(1,1,1,1), nx1, nx2, nx3)
            Z_grad(1) = fval(1,2); Z_grad(2) = fval(1,3); Z_grad(3) = fval(1,4)
         ELSE
            CALL r8fvtricub(ict1, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                            hx, hxi, hy, hyi, hz, hzi, &
                            Z4D(1,1,1,1), nx1, nx2, nx3)
         END IF
         Z_val = fval(1,1)
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
      REAL*8 :: fval(1,4)
      INTEGER, parameter :: ict1(10)=(/1,0,0,0,0,0,0,0,0,0/)
      INTEGER, parameter :: ict4(10)=(/1,1,1,1,0,0,0,0,0,0/)
      L_val = 0
      IF (ier < 0) RETURN
      rho_val = SQRT(s_val)
      !CALL EZspline_isInDomain(L_spl,u_val,v_val,rho_val,ier)
      !IF (ier == 0) THEN
      !   CALL EZspline_interp(L_spl,u_val,v_val,rho_val,L_val,ier)
      !   IF (PRESENT(L_grad)) CALL EZspline_gradient(L_spl,u_val,v_val,rho_val,L_grad,ier)
      IF (isingrid(u_val,v_val,rho_val)) THEN
         CALL lookupgrid3d(u_val,v_val,rho_val,i,j,k,hx,hy,hz,hxi,hyi,hzi,xparam,yparam,zparam)
         IF (PRESENT(L_grad)) THEN
            CALL r8fvtricub(ict4, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                            hx, hxi, hy, hyi, hz, hzi, &
                            L4D(1,1,1,1), nx1, nx2, nx3)
            L_grad(1) = fval(1,2); L_grad(2) = fval(1,3); L_grad(3) = fval(1,4)
         ELSE
            CALL r8fvtricub(ict1, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                            hx, hxi, hy, hyi, hz, hzi, &
                            L4D(1,1,1,1), nx1, nx2, nx3)
         END IF
         L_val = fval(1,1)
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
      DOUBLE PRECISION ::  R_grad(2),Z_grad(2)
      INTEGER :: i,j,k
      REAL*8 :: xparam, yparam, zparam, hx, hy, hz, hxi, hyi, hzi
      REAL*8 :: fval3(1,3), fval2(1,2)
      INTEGER, parameter :: ict3(10)=(/0,1,1,0,0,0,0,0,0,0/)
      INTEGER, parameter :: ict4(10)=(/1,1,1,0,0,0,0,0,0,0/)
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
         CALL r8fvtricub(ict4, 1, 1, fval3, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         R4D(1,1,1,1), nx1, nx2, nx3)
         R_val = fval3(1,1);
         R_grad(1) = fval3(1,2); R_grad(2) = fval3(1,3)
         CALL r8fvtricub(ict3, 1, 1, fval2, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         Z4D(1,1,1,1), nx1, nx2, nx3)
         Z_grad(1) = fval2(1,1); Z_grad(2) = fval2(1,2)
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
      REAL*8 :: fval(1,2)
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
         xp = fval(1,1); xpp = fval(1,2)
         CALL r8fvtricub(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         ZU4D(1,1,1,1), nx1, nx2, nx3)
         zp = fval(1,1); zpp = fval(1,2)
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



      SUBROUTINE get_equil_kappa2_dbl(s_val,u_val,v_val,phiedge,zeta_p,kappa2,kappa2v2,ier,diagnostic)
      USE EZspline
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in)    ::  s_val
      DOUBLE PRECISION, INTENT(in)    ::  u_val
      DOUBLE PRECISION, INTENT(in)    ::  v_val
      DOUBLE PRECISION, INTENT(in)    ::  phiedge
      DOUBLE PRECISION, INTENT(in)    ::  zeta_p
      DOUBLE PRECISION, INTENT(out)   ::  kappa2, kappa2v2
      INTEGER, INTENT(inout)     ::  ier
      DOUBLE PRECISION, INTENT(out),OPTIONAL   :: diagnostic(1,71) 
      real*8 :: rho_val
      real*8 :: R, d2R_dudv, d2R_duds, d2R_dvds
      real*8 :: Z, d2Z_dudv, d2Z_duds, d2Z_dvds

      real*8 :: Bsupu, Bsupv, Bsups
      real*8 :: B_R, dB_R_du, dB_R_dv, B_Z, dB_Z_du, dB_Z_dv, dB_Phi_du
      real*8 :: dB_R_ds, dB_Z_ds, dB_Phi_ds
      real*8 :: B_Phi, modB, B_X_v2, B_Y_v2, B_Z_v2
      real*8 :: sqrtg, norm_binormal
      real*8 :: subpart2, subpart3, subpart4
      real*8 :: R_grad(1,3), Z_grad(1,3), B_cyl(1,3), B_xyz(1,3)
      real*8 :: B_cyl2(1,3)
      real*8 :: esubs(1,3), esubu(1,3), esubv(1,3)
      real*8 :: es(1,3), eu(1,3), ev(1,3), grad_psi(1,3), binormal(1,3)
      real*8 :: binormal_xyz(1,3), grad_psi_xyz(1,3)
      real*8 :: subpart1(1,3), subpart5(1,3), subpart6(1,3), subpart7(1,3)
      real*8 :: subpart6v2(1,3)
      INTEGER :: i,j,k
      REAL*8 :: xparam, yparam, zparam, hx, hy, hz, hxi, hyi, hzi
      REAL*8 :: fval1(1,1)
      REAL*8 :: fval2(1,4)
      REAL*8 :: fval3(1,3)
      REAL*8 :: fval4(1,7)
      INTEGER, parameter :: ict1(10)=(/1,0,0,0,0,0,0,0,0,0/)
      INTEGER, parameter :: ict2(10)=(/1,1,1,1,0,0,0,0,0,0/)
      INTEGER, parameter :: ict3(10)=(/0,1,1,1,0,0,0,0,0,0/)
      INTEGER, parameter :: ict4(10)=(/1,1,1,1,0,0,0,1,1,1/)







      kappa2 = 0
      IF (ier < 0) RETURN
      rho_val = SQRT(s_val)
      IF (isingrid(u_val,v_val,rho_val)) THEN
         ! This will return the geodesic curvature of the surface, or
         ! (B cross grad(Psi)/|B cross grad(Psi)| ) dot
         !          [ (B^u grad(v) + B^v grad(u)) (dB_v/du - dB_u/dv)
         !             + grad(s) [ B^v (dB_s/dv-dB_v/ds) + B^u (dB_s/du -dB_u/ds) ] ]
         ! Three of the terms can be expanded in vmec coordinates:
         ! dB_v/du - dB_u/dv = dB_R/du dR/dv - dB_R/dv dR/du +
         !                     dB_Z/du dZ/dv - dB_Z/dv dZ/du +
         !                     (1/NFP) dB_Phi/du
         ! dB_s/dv-dB_v/ds  =   dB_R/dv dR/ds - dB_R/ds dR/dv + 
         !                  dB_Z/dv dZ/ds - dB_Z/ds dZ/dv - (1/NFP) dB_phi/ds
         ! dB_s/du-dB_u/ds =  dB_R/du dR/ds - dB_R/ds R/du +
         !                  dB_Z/du dZ/ds - dB_Z/ds dZ/du
         ! 
         ! Variables needed:
         ! Scalars
         ! B^u, B^v : get_equil_Bflx
         ! NFP: set elsewhere
         ! dR/dv, dR/du, dR/ds :gradR
         ! dZ/dv, dZ/du, dZ/ds : gradZ
         ! use new (not get_equil_Bcylsuv ?)
         ! dB_Phi/du, dB_Phi/ds
         ! dB_R/ds, dB_R/du, dB_R/dv, 
         ! dB_Z/ds, dB_Z/du, dB_Z/dv
         
         ! Vectors
         ! grad(Psi) , grad(s), grad(v), grad(u): see lines 346-399 of
         ! chisq_gamma_c_v2_mod.f90

         ! Gather the scalars that we need

         ! B^u, B^v : get_equil_Bflx
         !CALL get_equil_Bflx(phi_N, u, v, Bsups(j), Bsupu(j), Bsupv(j), ier, B_GRAD = gradB_init)
         CALL get_equil_Bflx(s_val, u_val, v_val, Bsups, Bsupu, Bsupv, ier)
         ! NFP: set elsewhere
         !nfp


         ! use new
         ! dR/dv, dR/du, dR/ds :gradR
         ! dZ/dv, dZ/du, dZ/ds : gradZ
         ! also see: get_equil_RZ(phi_N, u, v, R, Z, ier, gradR, gradZ)
         ! d^2R/(duds), d^2R/(dvds), d^2Z/(duds), d^2Z/(dvds)
         CALL lookupgrid3d(u_val,v_val,rho_val,i,j,k,hx,hy,hz,hxi,hyi,hzi,xparam,yparam,zparam)
         CALL r8fvtricub(ict2, 1, 1, fval4, i, j, k, xparam, yparam, zparam, &
                            hx, hxi, hy, hyi, hz, hzi, &
                            R4D(1,1,1,1), nx1, nx2, nx3)
         R = fval4(1,1)
         !R = fval4(1,0)
         R_grad(1,1) = fval4(1,2); R_grad(1,2) = fval4(1,3); R_grad(1,3) = fval4(1,4)
         !d2R_dudv = fval4(1,5); d2R_duds = fval4(1,6); d2R_dvds = fval4(1,7)

         CALL r8fvtricub(ict2, 1, 1, fval4, i, j, k, xparam, yparam, zparam, &
                            hx, hxi, hy, hyi, hz, hzi, &
                            Z4D(1,1,1,1), nx1, nx2, nx3)
         Z = fval4(1,1)
         Z_grad(1,1) = fval4(1,2); Z_grad(1,2) = fval4(1,3); Z_grad(1,3) = fval4(1,4)
         !d2Z_dudv = fval4(1,5); d2Z_duds = fval4(1,6); d2Z_dvds = fval4(1,7)

         !   Convert radial gradients from d/drho to d/ds
         !   phi_Ns = rho^2; ds/drho = 2*rho; drho/ds = 1/(2*rho) = 0.5/rovera
         !   - use rovera from above
         !   gradX(1) = dX/du, gradX(2) = dX/dv, gradX(3) = dX/dsqrt(s)
         R_grad(1,3) = R_grad(1,3) / (2.0 * rho_val)
         !d2R_duds = d2R_duds / (2.0 * rho_val)
         !d2R_dvds = d2R_dvds  / (2.0 * rho_val)
         Z_grad(1,3) = Z_grad(1,3) / (2.0 * rho_val)
         !d2Z_duds = d2Z_duds / (2.0 * rho_val)
         !d2Z_dvds = d2Z_dvds  / (2.0 * rho_val)




         ! dB_Phi/du, dB_Phi/ds
         CALL r8fvtricub(ict2, 1, 1, fval2, i, j, k, xparam, yparam, zparam, &
                            hx, hxi, hy, hyi, hz, hzi, &
                            B_Phi4D(1,1,1,1), nx1, nx2, nx3)
         B_Phi = fval2(1,1)
         dB_Phi_du = fval2(1,2)
         dB_Phi_ds = fval2(1,4)
  
       ! dB_R/du, dB_R/dv, dB_R/ds
         CALL r8fvtricub(ict2, 1, 1, fval2, i, j, k, xparam, yparam, zparam, &
                            hx, hxi, hy, hyi, hz, hzi, &
                            B_R4D(1,1,1,1), nx1, nx2, nx3)
         B_R = fval2(1,1)
         dB_R_du = fval2(1,2)
         dB_R_dv = fval2(1,3)*nfp
         dB_R_ds = fval2(1,4)

         ! dB_z/du, dB_Z/dv, dB_Z/ds
         CALL r8fvtricub(ict2, 1, 1, fval2, i, j, k, xparam, yparam, zparam, &
                            hx, hxi, hy, hyi, hz, hzi, &
                            B_Z4D(1,1,1,1), nx1, nx2, nx3)
         B_Z = fval2(1,1)
         dB_Z_du = fval2(1,2)
         dB_Z_dv = fval2(1,3)*nfp
         dB_Z_ds = fval2(1,4)


      !         !  
         ! Vectors
         ! grad(Psi) , grad(s), grad(v), grad(u): see lines 346-399 of
         ! chisq_gamma_c_v2_mod.f90
            esubs(1,1) = R_grad(1,3)    ! dR/ds
            esubs(1,2) = 0.0        ! dPhi/ds
            esubs(1,3) = Z_grad(1,3)    ! dZ/ds

            esubu(1,1) = R_grad(1,1)    ! dR/du
            esubu(1,2) = 0.0        ! dPhi/du
            esubu(1,3) = Z_grad(1,1)    ! dZ/du

            esubv(1,1) = R_grad(1,2)    ! dR/dv
            !esubv(2) = one
            esubv(1,2) = one/nfp         ! dPhi/dv:  v= nfp *mod(phi,2*pi/nfp) -> dv ~ nfp * dphi
            esubv(1,3) = Z_grad(1,2)    ! dZ/dv

            !  esubv x esubs = (R) esubv(2) * esubs(3) - esubv(3) * esubs(2) +
            !                 (Phi) esubv(3) * esubs(1) - esubv(1) * esubs(3) +
            !                  (Z) esubv(1) * esubs(2) - esubv(2) * esubs(1)
            sqrtg = esubu(1,1) * (esubv(1,2) * esubs(1,3)) +  &
                    esubu(1,2) * (esubv(1,3) * esubs(1,1)) +  &
                    esubu(1,3) * (-esubv(1,2) * esubs(1,1))


            CALL mycross(esubu,esubv,es)
            CALL mycross(esubv,esubs,eu)
            CALL mycross(esubs,esubu,ev)

!           es, eu, ev = e^s,e^u, e^v =  grad(s), grad(u), and grad(v) in cylindrical coordinates
!           (R, phi, Z) for the points on the single field period.
            es = es/sqrtg
            eu = eu/sqrtg
            ev = ev/sqrtg
            ! gradS is grad psi
            grad_psi = es * phiedge
            ! grad_psi is for the points on a single field period (not full torus)

         ! Now, build combine the terms
         ! (B cross grad(Psi)/|B cross grad(Psi)| ) dot
         !          [ (B^u grad(v) + B^v grad(u)) (dB_v/du - dB_u/dv)
         !old       + grad(s) ( B^v dB_s/dv - B^u dB_s/du ) ]
         !             + grad(s) [ B^v (dB_s/dv-dB_v/ds) + B^u (dB_s/du -dB_u/ds) ] ]
         ! Three of the terms can be expanded in vmec coordinates:
         ! dB_v/du - dB_u/dv = dB_R/du dR/dv - dB_R/dv dR/du +
         !                     dB_Z/du dZ/dv - dB_Z/dv dZ/du +
         !                     (1/NFP) dB_Phi/du
         ! dB_s/dv-dB_v/ds  =   dB_R/dv dR/ds - dB_R/ds dR/dv + 
         !         dB_Z/dv dZ/ds - dB_Z/ds dZ/dv - (1/NFP) dB_phi/ds
         ! dB_s/du-dB_u/ds =  dB_R/du dR/ds - dB_R/ds R/du +
         !         dB_Z/du dZ/ds - dB_Z/ds dZ/du


         !ok subpart1 = B^u grad(v) + B^v grad(u)
         subpart1 = Bsupu * ev + Bsupv * eu
         !subpart2 = dB_v/du - dB_u/dv
         !ok dB_v/du - dB_u/dv = dB_R/du dR/dv - dB_R/dv dR/du +
         !                     dB_Z/du dZ/dv - dB_Z/dv dZ/du +
         !                     (1/NFP) dB_Phi/du
         subpart2 = dB_R_du * R_grad(1,2) - dB_R_dv * R_grad(1,1) + &
                    dB_Z_du * Z_grad(1,2) - dB_Z_dv * Z_grad(1,1) + &
                    dB_Phi_du / nfp
         !old: subpart3 = B^v dB_s/dv
         !old:  B^v dB_s/dv  = B^v * [ B_R d^2R/(dvds) + 
         !old:                   dB_R/dv dR/ds + B_z d^2Z/(dvds) +
         !old:                   dB_Z/dv dZ/ds ]
         !old:  * phiedge, to normalized correctly
      !old subpart3 = Bsupv * (B_R * d2R_dvds + dB_R_dv * R_grad(1,3) + &
      !old              B_Z * d2Z_dvds + dB_Z_dv * Z_grad(1,3) ) 

     ! new subpart3 = B^v * ( dB_s/dv-dB_v/ds)
     !             = B^v*[  dB_R/dv dR/ds - dB_R/ds dR/dv + 
     !    dB_Z/dv dZ/ds - dB_Z/ds dZ/dv - (1/NFP) dB_phi/ds]
     ! 
         subpart3 = Bsupv * (dB_R_dv * R_grad(1,3) - &
                             dB_R_ds * R_grad(1,2) + &
                             dB_Z_dv * Z_grad(1,3)  - &
                             dB_Z_ds * Z_grad(1,2) - &
                             (one/NFP) * dB_Phi_ds ) 
         !old: subpart4 = B^u dB_s/du
         !old:  B^u dB_s/du = B^u * [ B_R d^2R/(duds) + 
         !old:                   dB_R/du dR/ds + B_Z d^2Z/(duds) +
         !old:                   dB_Z/du dZ/ds ]
         !old:  * phiedge, to normalized correctly
         !oldsubpart4 = Bsupu * (B_R * d2R_duds + dB_R_du * d2Z_duds + &
         !old              B_Z * d2Z_duds + dB_Z_du * Z_grad(1,3) )

         !new: subpart4 = B^u *(dB_s/du - dB_u/ds)

         ! Bsupu * ( dB_R/du dR/ds - dB_R/ds R/du +
         !         dB_Z/du dZ/ds - dB_Z/ds dZ/du )

         subpart4 = Bsupu * (dB_R_du * R_grad(1,3) - &
                             dB_R_ds * R_grad(1,1) + &
                             dB_Z_du * Z_grad(1,3) -  &
                             dB_Z_ds * Z_grad(1,1)   )
         !old: subpart5 = grad(s)*(subpart3 - subpart4)
         !old: subpart5 = gradS*(subpart3 - subpart4)
         !olssubpart5 = es*(subpart3 - subpart4)
         !new: subpart5 = gradS*(subpart3 + subpart4)
         subpart5 = es*(subpart3 + subpart4)
         ! ok - thi is the total curvature, kappa
         ! new: I think I only need the 'subpart1 * subpart2' portion
         !  the subpart 5 portion is normal to the surface and should not contribute
         ! to the projection into the binormal direction.
         !  (the binormal and normal direction are perpendicular by construction)
         subpart6 = subpart1 * subpart2 + subpart5
         subpart6v2 = subpart1 * subpart2 

         !kappa2 = dot_product((B cross grad(Psi)/|B cross grad(Psi)| ) , subpart6)
         ! subpart7 = (B cross grad(Psi)/|B cross grad(Psi)| ) 
         !CALL cross_product(grad_psi_xyz/norm_grad_psi_xyz, bnxyz, binormal(j,:))
         B_cyl(1,1) = B_R*sign(one,sqrtg)
         B_cyl(1,2) = B_Phi*sign(one,sqrtg)
         B_cyl(1,3) = B_Z*sign(one,sqrtg) ! why is this -1??? see gc12
         modB = sqrt(B_R**2 + B_Phi**2 + B_Z**2)
         B_X_v2 = sign(one,sqrtg) * ( B_R*cos(-zeta_p) - B_Phi*sin(-zeta_p) )
         B_Y_v2 = sign(one,sqrtg) * ( B_R*sin(-zeta_p) + B_Phi*cos(-zeta_p) )
         B_Z_v2 = B_Z *sign(one,sqrtg)
         B_cyl2(1,1) =  B_X_v2 * cos(-zeta_p) + B_Y_v2 * sin(-zeta_p) ! BR
         B_cyl2(1,2) = -B_X_v2 * sin(-zeta_p) + B_Y_v2  * cos(-zeta_p) !  BPhi
         !B_X_v2 = sign(one,sqrtg) * ( B_R*cos(v_val) - B_Phi*sin(v_val) )
         !B_Y_v2 = sign(one,sqrtg) * ( B_R*sin(v_val) + B_Phi*cos(v_val) )
         !B_Z_v2 = B_Z *sign(one,sqrtg)
         !B_cyl2(1,1) =  B_X_v2 * cos(v_val) + B_Y_v2 * sin(v_val) ! BR
         !B_cyl2(1,2) = -B_X_v2 * sin(v_val) + B_Y_v2  * cos(v_val) !  BPhi
         B_cyl2(1,3) =  B_Z_v2 ! Bz
         CALL mycross(es, B_cyl, binormal)
         !CALL mycross(es, B_cyl2, binormal)
         norm_binormal = sqrt(binormal(1,1)**2 + binormal(1,2)**2 + &
                              binormal(1,3)**2)
         subpart7 = binormal /( norm_binormal )
         ! ok - geodesic curvature
         kappa2 = subpart7(1,1) * subpart6(1,1) + subpart7(1,2) * subpart6(1,2) + &
                  subpart7(1,3) * subpart6(1,3)
         kappa2v2 = subpart7(1,1) * subpart6v2(1,1) + subpart7(1,2) * subpart6v2(1,2) + &
                  subpart7(1,3) * subpart6v2(1,3)
!
         IF (PRESENT(diagnostic))  THEN
            diagnostic(1,1) = s_val
            diagnostic(1,2) = rho_val
            diagnostic(1,3) = u_val
            diagnostic(1,4) = v_val
            diagnostic(1,5) = B_R
            diagnostic(1,6) = B_Phi
            diagnostic(1,7) = B_Z
            diagnostic(1,8) = Bsups
            diagnostic(1,9) = Bsupu
            diagnostic(1,10) = Bsupv
            diagnostic(1,11) = dB_Phi_du
            diagnostic(1,12) = dB_R_du
            diagnostic(1,13) = dB_R_dv
            diagnostic(1,14) = dB_Z_du
            diagnostic(1,15) = dB_Z_dv
            diagnostic(1,16) = R
            diagnostic(1,17) = R_grad(1,1)
            diagnostic(1,18) = R_grad(1,2)
            diagnostic(1,19) = R_grad(1,3)
            diagnostic(1,20) = B_cyl2(1,1)
            diagnostic(1,21) = B_cyl2(1,2)
            diagnostic(1,22) = B_cyl2(1,3)
            diagnostic(1,23) = Z
            diagnostic(1,24) = Z_grad(1,1)
            diagnostic(1,25) = Z_grad(1,2)
            diagnostic(1,26) = Z_grad(1,3)
            diagnostic(1,27) = d2Z_dudv
            diagnostic(1,28) = d2Z_duds
            diagnostic(1,29) = d2Z_dvds
            diagnostic(1,30) = esubs(1,1)
            diagnostic(1,31) = esubs(1,2)
            diagnostic(1,32) = esubs(1,3)
            diagnostic(1,33) = esubu(1,1)
            diagnostic(1,34) = esubu(1,2)
            diagnostic(1,35) = esubu(1,3)
            diagnostic(1,36) = esubv(1,1)
            diagnostic(1,37) = esubv(1,2)
            diagnostic(1,38) = esubv(1,3)
            diagnostic(1,39) = es(1,1)
            diagnostic(1,40) = es(1,2)
            diagnostic(1,41) = es(1,3)
            diagnostic(1,42) = eu(1,1)
            diagnostic(1,43) = eu(1,2)
            diagnostic(1,44) = eu(1,3)
            diagnostic(1,45) = ev(1,1)
            diagnostic(1,46) = ev(1,2)
            diagnostic(1,47) = ev(1,3)
            diagnostic(1,48) = sqrtg
            diagnostic(1,49) = phiedge
            diagnostic(1,50) = subpart1(1,1)
            diagnostic(1,51) = subpart1(1,2)
            diagnostic(1,52) = subpart1(1,3)
            diagnostic(1,53) = subpart2
            diagnostic(1,54) = subpart3
            diagnostic(1,55) = subpart4
            diagnostic(1,56) = subpart5(1,1)
            diagnostic(1,57) = subpart5(1,2)
            diagnostic(1,58) = subpart5(1,3)
            diagnostic(1,59) = subpart6(1,1)
            diagnostic(1,60) = subpart6(1,2)
            diagnostic(1,61) = subpart6(1,3)
            diagnostic(1,62) = binormal(1,1)
            diagnostic(1,63) = binormal(1,2)
            diagnostic(1,64) = binormal(1,3)
            diagnostic(1,65) = subpart7(1,1)
            diagnostic(1,66) = subpart7(1,2)
            diagnostic(1,67) = subpart7(1,3)
            diagnostic(1,68) = subpart6v2(1,1)
            diagnostic(1,69) = subpart6v2(1,2)
            diagnostic(1,70) = subpart6v2(1,3)
            diagnostic(1,71) = zeta_p

         END IF

      ELSE
         ier=9
      END IF
      RETURN
      END SUBROUTINE get_equil_kappa2_dbl
      
      SUBROUTINE get_equil_kappa2_sgl(s_val,u_val,v_val,phiedge,zeta_p,kappa2,kappa2v2,ier)
      USE EZspline
      IMPLICIT NONE
      REAL, INTENT(in)    ::  s_val
      REAL, INTENT(in)    ::  u_val
      REAL, INTENT(in)    ::  v_val
      REAL, INTENT(in)    ::  phiedge
      REAL, INTENT(in)    ::  zeta_p
      REAL, INTENT(out)   ::  kappa2, kappa2v2 
      INTEGER, INTENT(inout)     ::  ier
      DOUBLE PRECISION    ::  s_dbl
      DOUBLE PRECISION    ::  u_dbl
      DOUBLE PRECISION    ::  v_dbl
      DOUBLE PRECISION    ::  phiedge_dbl
      DOUBLE PRECISION    ::  zeta_p_dbl
      DOUBLE PRECISION    ::  kappa2_dbl, kappa2v2_dbl
      s_dbl = s_val; u_dbl = u_val; v_dbl = v_val; phiedge_dbl = phiedge
      zeta_p_dbl = zeta_p
      CALL get_equil_kappa2_dbl(s_dbl,u_dbl,v_dbl,phiedge_dbl,zeta_p_dbl,kappa2_dbl,kappa2v2_dbl,ier)
      kappa2 = kappa2_dbl
      kappa2v2 = kappa2v2_dbl  
      RETURN
      END SUBROUTINE get_equil_kappa2_sgl


      SUBROUTINE mycross(a, b, c)
        real*8 , INTENT(IN)  :: a(1,3), b(1,3)
        real*8 , INTENT(OUT) :: c(1,3)
      
        c(1,1) = a(1,2) * b(1,3) - a(1,3) * b(1,2)
        c(1,2) = a(1,3) * b(1,1) - a(1,1) * b(1,3)
        c(1,3) = a(1,1) * b(1,2) - a(1,2) * b(1,1)
      END SUBROUTINE mycross
 

      
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
      REAL*8 :: fval1(1)
      REAL*8 :: fval2(1,3)
      INTEGER, parameter :: ict1(10)=(/1,0,0,0,0,0,0,0,0,0/)
      INTEGER, parameter :: ict3(10)=(/0,1,1,1,0,0,0,0,0,0/)
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
         CALL r8fvtricub(ict1, 1, 1, fval1, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         R4D(1,1,1,1), nx1, nx2, nx3)
         r_val = fval1(1)
         CALL r8fvtricub(ict1, 1, 1, fval1, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         BS4D(1,1,1,1), nx1, nx2, nx3)
         Bs = fval1(1)
         CALL r8fvtricub(ict1, 1, 1, fval1, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         BU4D(1,1,1,1), nx1, nx2, nx3)
         Bu = fval1(1)
         CALL r8fvtricub(ict1, 1, 1, fval1, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         BV4D(1,1,1,1), nx1, nx2, nx3)
         Bv = fval1(1)
         CALL r8fvtricub(ict1, 1, 1, fval1, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         RU4D(1,1,1,1), nx1, nx2, nx3)
         R_grad(1) = fval1(1)
         CALL r8fvtricub(ict1, 1, 1, fval1, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         RV4D(1,1,1,1), nx1, nx2, nx3)
         R_grad(2) = fval1(1)
         CALL r8fvtricub(ict1, 1, 1, fval1, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         ZU4D(1,1,1,1), nx1, nx2, nx3)
         Z_grad(1) = fval1(1)
         CALL r8fvtricub(ict1, 1, 1, fval1, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         ZV4D(1,1,1,1), nx1, nx2, nx3)
         Z_grad(2) = fval1(1)
         br = R_grad(3)*Bs + R_grad(1)*Bu + R_grad(2)*Bv*nfp
         bphi = r_val * Bv
         bz = Z_grad(3)*Bs + Z_grad(1)*Bu + Z_grad(2)*Bv*nfp
!         IF (PRESENT(B_grad)) CALL EZspline_gradient(B_spl,u_val,v_val,s_val,B_grad,ier)
         IF (PRESENT(B_grad))  THEN
            CALL r8fvtricub(ict3, 1, 1, fval2, i, j, k, xparam, yparam, zparam, &
                           hx, hxi, hy, hyi, hz, hzi, &
                           B4D(1,1,1,1), nx1, nx2, nx3)
            B_grad(1) = fval2(1,1); B_grad(2) = fval2(1,2); B_grad(3) = fval2(1,3)
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
      
      SUBROUTINE get_equil_Bcylsuv2_dbl(s_val,u_val,v_val,ier,dBphidpsi)
      USE EZspline
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in)    ::  s_val
      DOUBLE PRECISION, INTENT(in)    ::  u_val
      DOUBLE PRECISION, INTENT(in)    ::  v_val
      DOUBLE PRECISION, INTENT(out)   ::  dBphidpsi(3)
      INTEGER, INTENT(inout)     ::  ier
      DOUBLE PRECISION :: rho_val
      INTEGER :: i,j,k
      REAL*8 :: xparam, yparam, zparam, hx, hy, hz, hxi, hyi, hzi
      REAL*8 :: fval1(1)
      REAL*8 :: fval2(1,3)
      INTEGER, parameter :: ict1(10)=(/1,0,0,0,0,0,0,0,0,0/)
      INTEGER, parameter :: ict3(10)=(/0,1,1,1,0,0,0,0,0,0/)
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
         !CALL r8fvtricub(ict1, 1, 1, fval1, i, j, k, xparam, yparam, zparam, &
         !                hx, hxi, hy, hyi, hz, hzi, &
         !                BV4D(1,1,1,1), nx1, nx2, nx3)
         !Bv = fval1(1)
            CALL r8fvtricub(ict3, 1, 1, fval2, i, j, k, xparam, yparam, zparam, &
                           hx, hxi, hy, hyi, hz, hzi, &
                           BV4D(1,1,1,1), nx1, nx2, nx3)
            dBphidpsi(1) = fval2(1,1)
            dBphidpsi(2) = fval2(1,2)
            dBphidpsi(3) = fval2(1,3)
      ELSE
         ier   = 9
            dBphidpsi(1) = 0
            dBphidpsi(2) = 0
            dBphidpsi(3) = 0
      END IF
      RETURN
      END SUBROUTINE get_equil_Bcylsuv2_dbl
      
      SUBROUTINE get_equil_Bcylsuv2_sgl(s_val,u_val,v_val,ier,dBphidpsi)
      IMPLICIT NONE
      REAL, INTENT(in)    ::  s_val
      REAL, INTENT(in)    ::  u_val
      REAL, INTENT(in)    ::  v_val
      REAL, INTENT(out)   ::  dBphidpsi(3)
      INTEGER, INTENT(inout)     ::  ier
      DOUBLE PRECISION    ::  s_dbl
      DOUBLE PRECISION    ::  u_dbl
      DOUBLE PRECISION    ::  v_dbl
      DOUBLE PRECISION    ::  dBphidpsi_dbl(3)
      s_dbl = s_val; u_dbl = u_val; v_dbl = v_val
      CALL get_equil_Bcylsuv2_dbl(s_dbl,u_dbl,v_dbl,ier,dBphidpsi_dbl)
      dBphidpsi = dBphidpsi_dbl
      RETURN
      END SUBROUTINE get_equil_Bcylsuv2_sgl
      
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
      REAL*8 :: fval1(1)
      REAL*8 :: fval2(1,3)
      INTEGER, parameter :: ict1(10)=(/1,0,0,0,0,0,0,0,0,0/)
      INTEGER, parameter :: ict2(10)=(/0,1,1,1,0,0,0,0,0,0/)
      IF (ier < 0) RETURN
      rho_val = SQRT(s_val)
!      CALL EZspline_interp(Bs_spl,u_val,v_val,rho_val,bs,ier)
!      CALL EZspline_interp(Bu_spl,u_val,v_val,rho_val,bu,ier)
!      CALL EZspline_interp(Bv_spl,u_val,v_val,rho_val,bv,ier)
!      IF (PRESENT(modb_val)) CALL EZspline_interp(B_spl,u_val,v_val,rho_val,modb_val,ier)
!      IF (PRESENT(B_grad)) CALL EZspline_gradient(B_spl,u_val,v_val,rho_val,B_grad,ier)
      CALL lookupgrid3d(u_val,v_val,rho_val,i,j,k,hx,hy,hz,hxi,hyi,hzi,xparam,yparam,zparam)
      CALL r8fvtricub(ict1, 1, 1, fval1, i, j, k, xparam, yparam, zparam, &
                      hx, hxi, hy, hyi, hz, hzi, &
                      BS4D(1,1,1,1), nx1, nx2, nx3)
      bs = fval1(1)
      CALL r8fvtricub(ict1, 1, 1, fval1, i, j, k, xparam, yparam, zparam, &
                      hx, hxi, hy, hyi, hz, hzi, &
                      BU4D(1,1,1,1), nx1, nx2, nx3)
      bu = fval1(1)
      CALL r8fvtricub(ict1, 1, 1, fval1, i, j, k, xparam, yparam, zparam, &
                      hx, hxi, hy, hyi, hz, hzi, &
                      BV4D(1,1,1,1), nx1, nx2, nx3)
      bv = fval1(1)
      IF (PRESENT(modb_val)) THEN
         CALL r8fvtricub(ict1, 1, 1, fval1, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         B4D(1,1,1,1), nx1, nx2, nx3)
         modb_val = fval1(1)
      END IF
      IF (PRESENT(B_GRAD)) THEN
         CALL r8fvtricub(ict2, 1, 1, fval2, i, j, k, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, &
                         B4D(1,1,1,1), nx1, nx2, nx3)
         B_grad(1) = fval2(1,1); B_grad(2) = fval2(1,2); B_grad(3) = fval2(1,3)
      END IF
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
      REAL*8 :: fval(1,2)
      INTEGER, parameter :: ict1(3)=(/1,0,0/)
      INTEGER, parameter :: ict2(3)=(/1,1,0/)
      IF (ier < 0) RETURN
      rho_val = SQRT(s_val)
      !CALL EZspline_interp(Bav_spl,rho_val,Bav,ier)
      !CALL EZspline_interp(Bsq_spl,rho_val,Bsqav,ier)
      CALL lookupgrid1d(rho_val,k,hz,hzi,zparam)
      CALL r8fvspline(ict1,1,1,fval,k,zparam,hz,hzi,BAV2D(1,1),nx3)
      Bav = fval(1,1)
      IF (PRESENT(Bsqavp_val))  THEN
         CALL r8fvspline(ict2,1,1,fval,k,zparam,hz,hzi,BSQ2D(1,1),nx3)
         Bsqavp_val = fval(1,2)
      ELSE
         CALL r8fvspline(ict1,1,1,fval,k,zparam,hz,hzi,BSQ2D(1,1),nx3)
      END IF
      Bsqav = fval(1,1)
      IF (PRESENT(Bsqavp_val))  THEN
         !CALL EZspline_derivative(Bsq_spl,1,rho_val,Bsqavp_val,ier)
         !CALL EZspline_interp(Vp_spl,rho_val,vp_val,ier)
         CALL r8fvspline(ict1,1,1,fval,k,zparam,hz,hzi,VP2D(1,1),nx3)
         vp_val = fval(1,1)
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
      !th = coord(2)
      !phi = coord(3)
      th = coord(2)
      phi = coord(3)

      ! suppose nfp = 4, then:  \phi \in -pi/2 to pi/2: phi goes to phi*nfp
      !                         \phi \in pi/2 to pi: phi goes to (phi - pi/2)* nfp
      !                         \phi \in pi to 3*pi/2: phi goes to (phi - pi)* nfp
      !                         \phi \in 3*pi/2 to 2*pi: phi goes to (phi - 3*pi/2)* nfp
      ! at 2pi/nfp = pi/2   (-)  : phi*nfp     at pi/2   (+): 0*nfp
      !  at pi     (-)  : phi*nfp/2   at pi     (+): 0*nfp
      !  at 3*pi/2 (-)  : phi*nfp/2   at 3*pi/2 (+): 0*nfp
      !  at 2*pi   (-)  : phi*nfp/2   at 2*pi   (+): 0*nfp

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
         IF (th < 0) then
           th = th + pi2
           th1 = th1 + pi2
         end if
         IF (th > pi2) then
           th = MOD(th,pi2)
           th1 = th1 - pi2
         end if
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
         ! JCS - how does this behave near the th~th1~2pi boundary?
         dth = -(th + lam - th1)/(one+dlam)
         n1 = n1 + 1
!         write (*,'(A,I4,(2es22.12))') "<--n1,th,dth= ",n1,th,dth
         th = th + 0.5*dth
!         if (n1 .eq. 500) write (*,'(A,(3es22.12))') "<--1 = 500,s,th,phi= ",s,th,phi
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
      
      SUBROUTINE line_int_dbl(fcn,r1,r2,val,length)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in)   :: r1(3), r2(3)
      DOUBLE PRECISION, INTENT(out)  :: val
      DOUBLE PRECISION, INTENT(out), OPTIONAL :: length
      INTEGER     :: i, j, k, ier
      DOUBLE PRECISION :: x1, y1, z1, x2, y2, z2, dx, dy, dz, xp, yp, zp, int_fac
      DOUBLE PRECISION :: xp1, yp1, zp1, delf, delt, rp, phip, s, dx2, dy2, dz2, uval
      DOUBLE PRECISION :: f_val
      INTEGER, PARAMETER :: nop=3
      INTEGER, PARAMETER :: int_step=2
      DOUBLE PRECISION, dimension(nop), parameter :: ci=(/1./6.,2./3.,1./6./)
      EXTERNAL fcn
      phip = r1(2)
      s    = r2(2)
      ier = 0; f_val = 0; val = 0
      IF (PRESENT(length)) length = 0
      x1 = r1(1)*cos(phip); x2 = r2(1)*cos(s)
      y1 = r1(1)*sin(phip); y2 = r2(1)*sin(s)
      z1 = r1(3);           z2 = r2(3);
      dx = (x2-x1)/lintsteps
      dy = (y2-y1)/lintsteps
      dz = (z2-z1)/lintsteps
      s=2; phip=0; i = 1
      DO WHILE ((i <= lintsteps) .and. (s>1))
         xp = x1+dx*(i-1)
         yp = y1+dy*(i-1)
         zp = z1+dz*(i-1)
         rp = sqrt(xp*xp+yp*yp)
         phip = ATAN2(yp,xp)
         IF (phip < 0) phip = phip+pi2
         CALL get_equil_s(rp,phip,zp,s,ier)
         i = i + 1
      END DO
      IF (i == lintsteps) RETURN
      x1 = xp-dx
      y1 = yp-dy
      z1 = zp-dz
      i=1; s=2
      DO WHILE ((i <= lintsteps) .and. (s>1))
         xp = x2-dx*(i-1)
         yp = y2-dy*(i-1)
         zp = z2-dz*(i-1)
         rp = sqrt(xp*xp+yp*yp)
         phip = ATAN2(yp,xp)
         IF (phip < 0) phip = phip+pi2
         CALL get_equil_s(rp,phip,zp,s,ier)
         i = i + 1
      END DO
      IF (i== lintsteps) RETURN
      x2 = xp+dx
      y2 = yp+dy
      z2 = zp+dz
      dx = (x2-x1)/lintsteps
      dy = (y2-y1)/lintsteps
      dz = (z2-z1)/lintsteps
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
