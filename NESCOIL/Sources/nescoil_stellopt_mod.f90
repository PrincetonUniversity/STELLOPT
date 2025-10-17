!-----------------------------------------------------------------------
!     Module:        nescoil_stellopt_mod
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          10/17/2025
!     Description:   This module contains routines for interfacing
!                    NESCOIL and STELLOPT
!-----------------------------------------------------------------------
      MODULE nescoil_stellopt_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE Vmeshes
      USE Vprecal1
      USE Vvacuum2, ONLY: iota_edge, phip_edge, curpol, ms1, ns1, cr1, &
                          cz1, cl1
      USE Vvacuum3, ONLY: ibex, cup, cut
      USE SvdCtrl,  ONLY: mstrt, mstep, mkeep, mdspw, curwt, trgwt, noaccuracy
      USE OutCtrl,  ONLY: w_psurf, w_csurf, w_bnuv, w_jsurf, w_xerr,   &
                          w_svd
      USE Vvacuum1, ONLY: cr, cz, ms, ns
      USE Vbiopre,  ONLY: vx, vy, vz, dx, dy, dz
      USE Vculine1, ONLY: xw, yw, zw, curre, nw
      USE Vprecal4, ONLY: cok, sik
      USE NumParams, ONLY : inesc, zero
!-----------------------------------------------------------------------
!     Module Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      
!-----------------------------------------------------------------------
!     Module SUBROUTINES/FUNCTIONS
!-----------------------------------------------------------------------
      CONTAINS

      SUBROUTINE set_nescoil_grid(nu_in, nv_in, nu1_in, nv1_in, npol_in, ntor_in)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: nu_in, nv_in, nu1_in, nv1_in, npol_in, ntor_in
      nu = nu_in
      nv = nv_in
      nu1 = nu1_in
      nv1 = nv1_in
      npol = npol_in
      ntor = ntor_in
      nmax  = 3*nv*(nu+2)+1
      nuv = nu*nv
      nuv1 = nu1*nv1
      nuvh = nuv/2 + nu
      nuvh1 = nuv1/2 + nu1
      RETURN
      END SUBROUTINE set_nescoil_grid

      SUBROUTINE set_nescoil_fourier(mf_in, nf_in, md_in, nd_in)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: mf_in, nf_in, md_in, nd_in
      mf = mf_in
      nf = nf_in
      md = md_in
      nd = nd_in
      mnd   = (md + 1)*(2*nd + 1)
      RETURN
      END SUBROUTINE set_nescoil_fourier

      SUBROUTINE set_nescoil_plasma(np_in, iota_edge_in, phip_edge_in, curpol_in)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: np_in
      REAL(rprec), INTENT(in) :: iota_edge_in, phip_edge_in, curpol_in
      np = np_in
      iota_edge = iota_edge_in
      phip_edge = phip_edge_in
      curpol = curpol_in
      RETURN
      END SUBROUTINE set_nescoil_plasma

      SUBROUTINE set_nescoil_current(cut_in, cup_in, ibex_in)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: ibex_in
      REAL(rprec), INTENT(in) :: cut_in, cup_in
      cut = cut_in
      cup = cup_in
      ibex = ibex_in
      RETURN
      END SUBROUTINE set_nescoil_current

      SUBROUTINE set_nescoil_svd(mstrt_in, mstep_in, mkeep_in, mdspw_in, curwt_in, trgwt_in, noaccuracy_in)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: mstrt_in, mstep_in, mkeep_in, mdspw_in, noaccuracy_in
      REAL(rprec), INTENT(in) :: curwt_in, trgwt_in
      mstrt = mstrt_in
      mstep = mstep_in
      mkeep = mkeep_in
      mdspw = mdspw_in
      curwt = curwt_in
      trgwt = trgwt_in
      noaccuracy = noaccuracy_in
      RETURN
      END SUBROUTINE set_nescoil_svd

      SUBROUTINE set_nescoil_output(w_psurf_in, w_csurf_in, w_bnuv_in, w_jsurf_in, w_xerr_in, w_svd_in)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: w_psurf_in, w_csurf_in, w_bnuv_in, w_jsurf_in, w_xerr_in, w_svd_in
      w_psurf = w_psurf_in
      w_csurf = w_csurf_in
      w_bnuv  = w_bnuv_in
      w_jsurf = w_jsurf_in
      w_xerr  = w_xerr_in
      w_svd   = w_svd_in
      RETURN
      END SUBROUTINE set_nescoil_output

      SUBROUTINE set_nescoil_plasma_boundary(mnmax_in,xm_in,xn_in,rbc_in,zbs_in,lbs_in)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: mnmax_in
      REAL(rprec), DIMENSION(mnmax_in), INTENT(in) :: xm_in, xn_in, rbc_in, zbs_in, lbs_in
      INTEGER :: mn, m, n
      ms1 = INT(MAXVAL(ABS(xm_in)))
      ns1 = INT(MAXVAL(ABS(xn_in)))
      IF (ALLOCATED(cr1)) DEALLOCATE(cr1)
      IF (ALLOCATED(cz1)) DEALLOCATE(cz1)
      IF (ALLOCATED(cl1)) DEALLOCATE(cl1)
      ALLOCATE(cr1(0:md,-nd:nd), cz1(0:md,-nd:nd), cl1(0:md,-nd:nd))
      cr1 = zero; cz1 = zero; cl1 = zero
      DO mn = 1, mnmax_in
         m = INT(xm_in(mn))
         n = INT(xn_in(mn))
         cr1(m,n) = rbc_in(mn)
         cz1(m,n) = zbs_in(mn)
         cl1(m,n) = lbs_in(mn)
      END DO
      RETURN
      END SUBROUTINE set_nescoil_plasma_boundary

      SUBROUTINE set_nescoil_current_surface(mnmax_in,xm_in,xn_in,rbc_in,zbs_in,lbs_in)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: mnmax_in
      REAL(rprec), DIMENSION(mnmax_in), INTENT(in) :: xm_in, xn_in, rbc_in, zbs_in, lbs_in
      INTEGER :: mn, m, n
      ms = INT(MAXVAL(ABS(xm_in)))
      ns = INT(MAXVAL(ABS(xn_in)))
      IF (ALLOCATED(cr)) DEALLOCATE(cr)
      IF (ALLOCATED(cz)) DEALLOCATE(cz)
      ALLOCATE(cr(0:md,-nd:nd), cz(0:md,-nd:nd))
      cr = zero; cz = zero
      DO mn = 1, mnmax_in
         m = INT(xm_in(mn))
         n = INT(xn_in(mn))
         cr(m,n) = rbc_in(mn)
         cz(m,n) = zbs_in(mn)
      END DO
      RETURN
      END SUBROUTINE set_nescoil_current_surface

      SUBROUTINE run_nescoil(lscreen)
      IMPLICIT NONE
      LOGICAL, INTENT(in) :: lscreen
      INTEGER :: istat
      real(rprec) :: t1, t2, tbeg, tend
      character*10 :: date, time
      IF (lscreen) THEN
         write (inesc, '(A)') '-----  Begin Nescoil run  -----'
         call date_and_time(date,time)
         write (inesc, 100) date(5:6),date(7:8),date(1:4),time(1:2),time(3:4),time(5:6)
      END IF
      call second0(tbeg)
      nw = ntor*np*(npol+1)
      istat = 0
      IF (ALLOCATED(vx)) DEALLOCATE(vx)
      IF (ALLOCATED(vy)) DEALLOCATE(vy)
      IF (ALLOCATED(vz)) DEALLOCATE(vz)
      IF (ALLOCATED(dx)) DEALLOCATE(dx)
      IF (ALLOCATED(dy)) DEALLOCATE(dy)
      IF (ALLOCATED(dz)) DEALLOCATE(dz)
      IF (ALLOCATED(xw)) DEALLOCATE(xw)
      IF (ALLOCATED(yw)) DEALLOCATE(yw)
      IF (ALLOCATED(zw)) DEALLOCATE(zw)
      IF (ALLOCATED(curre)) DEALLOCATE(curre)
      IF (ALLOCATED(cok)) DEALLOCATE(cok)
      IF (ALLOCATED(sik)) DEALLOCATE(sik)
      ALLOCATE(vx(nw), vy(nw), vz(nw), dx(nw), dy(nw), dz(nw), &
               xw(nw), yw(nw), zw(nw), curre(nw), &
               cok(0:np - 1), sik(0:np - 1), stat=istat)
      vx = zero; vy = zero; vz = zero
      dx = zero; dy = zero; dz = zero
      xw = zero; yw = zero; zw = zero
      curre = zero; cok = zero; sik = zero
      ! Do all precalculations
      CALL precal
      CALL second0(t2)
      IF(lscreen) write (inesc,"('PRECAL took ',g12.3,' sec')") t2-tbeg
      ! Compute quantities on plasma surface
      IF(lscreen) write (inesc, '(A)') '----- Calling Surface_Plasma -----'
      CALL second0(t1)
      CALL surface_plas
      CALL second0(t2)
      IF(lscreen) write (inesc,"('Time in Surface_Plasma: ',g12.3,' sec')") t2-t1
      ! Compute quantities on coil surface
      IF(lscreen) write (inesc, '(A)') '----- Calling Surface_Coil -----'
      CALL second0(t1)
      CALL surface_coil
      CALL second0(t2)
      IF(lscreen) write (inesc,"('Time in Surface_Coil: ',g12.3,' sec')") t2-t1
      ! Solve boundary value problem
      IF(lscreen) write (inesc, '(A)') '----- Calling Solver -----'
      CALL second0(t1)
      CALL solver_nescoil
      CALL second0(t2)
      IF(lscreen) write (inesc,"('Time in Solver: ',g12.3,' sec')") t2-t1
      ! Post-process solution for various answers
      IF(lscreen) write (inesc, '(A)') '----- Calling Surfcur_Diag -----'
      CALL second0(t1)
      CALL surfcur_diag
      CALL second0(t2)
      IF(lscreen) write (inesc,"('Time in Surfcur_Diag: ',g12.3,' sec')") t2-t1
      IF(noaccuracy .eq. 0) then
         IF(lscreen)write (inesc, '(A)') '----- Calling Accuracy -----'
         CALL second0(t1)
         CALL accuracy
         CALL second0(t2)
         IF(lscreen)write (inesc,"('ACCURACY took ',g12.3,' sec')") t2-t1
      ENDIF
      CALL second0(tend)
      IF(lscreen) write (inesc,"('ONE NESCOIL RUN took ',g12.3,' sec')") tend-tbeg
      RETURN
 100  format('DATE = ',a2,'-',a2,'-',a4,' ',' TIME = ',2(a2,':'),a2)
      END SUBROUTINE run_nescoil

!-----------------------------------------------------------------------
!     End Module
!-----------------------------------------------------------------------
      END MODULE nescoil_stellopt_mod