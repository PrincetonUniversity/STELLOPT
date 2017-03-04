
      SUBROUTINE init_ajax
c-----------------------------------------------
c   M o d u l e s
c-----------------------------------------------
      USE stel_kinds
      USE AJAX_MOD
      USE vmec_input, only: phiedge, lasym
      USE optim_params, only: rgrid_min, rgrid_max, zgrid_min, zgrid_max
      USE optim
      USE boozer_params
      USE safe_open_mod
      USE system_mod
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------

c-----------------------------------------------
c   L o c a l   V a r i a b l e s
c-----------------------------------------------
      INTEGER :: i, iflag
      INTEGER, DIMENSION(mnmax_opt) :: m, n
      REAL(rprec), DIMENSION(nrad, mnmax_opt) :: rc_aj, zs_aj,
     1                                           rs_aj, zc_aj
      REAL(rprec), DIMENSION(nrad-1, mnmax_opt) :: lams_aj, lamc_aj
      REAL(rprec) :: hs
      REAL(rprec), DIMENSION(nrad) :: s_aj, sl_aj
      CHARACTER*120 :: message

c-----------------------------------------------
      rc_aj   = transpose(rmnc_opt)
      zs_aj   = transpose(zmns_opt)
      lams_aj = transpose(lmns_opt(:,2:nrad))

      m = NINT(xm_bdy)
      n = NINT(xn_bdy)  ! SAM 1/11/11

      hs = 1/REAL((nrad - 1),rprec)
      s_aj(1) = 0
      DO i=2, nrad
         s_aj(i) = (i-1)*hs
         sl_aj(i-1) = (i-1.5_rprec)*hs
      END DO

      IF (.not. lasym) then
!     2         NTHETA_AJAX=4*mpol_opt+1,
!     2         NZETA_AJAX=(4*ntor_opt)*nfp_opt+1,
         CALL ajax_load_rzlam(nrad, mnmax_opt, s_aj, m, n, rc_aj, zs_aj,
     1         iflag, message,
     1         K_GRID=1, L_MFILTER_AJAX=.false., NRHO_AJAX=nrad,
     2         NTHETA_AJAX=721,
     2         NZETA_AJAX=721, RHOMAX_AJAX=1._rprec,
     2         NR_LAM=nrad-1, NK_LAM=mnmax_opt,
     3         RHO_LAM=sl_aj, LAM=lams_aj)
      ELSE
         rs_aj   = transpose(rmns_opt)
         zc_aj   = transpose(zmnc_opt)
         lamc_aj = transpose(lmnc_opt(:,2:nrad))
         CALL ajax_load_rzlam(nrad, mnmax_opt, s_aj, m, n, rc_aj, zs_aj,
     1         iflag, message,
     2         K_GRID=1, L_MFILTER_AJAX=.false., NRHO_AJAX=nrad,
     3         NTHETA_AJAX=4*mpol_opt+1,
     4         NZETA_AJAX=(4*ntor_opt)*nfp_opt+1, RHOMAX_AJAX=1._rprec,
     5         NR_LAM=nrad-1, NK_LAM=mnmax_opt,
     6         RHO_LAM=sl_aj, LAM=lams_aj,
     7         R_S=rs_aj, Z_C=zc_aj, LAM_C=lamc_aj)
      ENDIF
         

      IF (iflag < 0) then
         PRINT *,' Init_Ajax Load_RZLam warning:'
         PRINT *,message

      ELSE IF (iflag > 0) THEN
         PRINT *,' Init_Ajax Load_RZLam Error:'
         PRINT *,message
         STOP
      ENDIF

!      CALL ajax_load_magflux(phiedge, 1, nrad, s_aj, iota_opt, iflag,
!     1                       message)
      CALL ajax_load_magflux(phiedge, 1, nrad-1, sl_aj,
     1                       iota_opt(2:nrad), iflag,
     1                       message)

      IF (iflag < 0) THEN
         print *,' Init_Ajax Load_MagFlux warning:'
         print *,message

      ELSE IF (iflag > 0) THEN
         PRINT *,' Init_Ajax Load_MagFlux Error:'
         PRINT *,message
         STOP
      END IF
      END SUBROUTINE init_ajax
