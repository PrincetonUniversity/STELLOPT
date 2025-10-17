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
      USE NumParams, ONLY : inesc, zero, one, pi2
      USE Vdiagno2, ONLY: pot
      USE bsc_T
      USE biotsavart, ONLY: coil_group, nfp => nfp_bs
      USE stel_kinds, ONLY: rprec
!-----------------------------------------------------------------------
!     Module Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, PRIVATE :: factor
      
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
      factor = pi2/np
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

      SUBROUTINE set_nescoil_current_surface(mnmax_in,xm_in,xn_in,rbc_in,zbs_in)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: mnmax_in
      REAL(rprec), DIMENSION(mnmax_in), INTENT(in) :: xm_in, xn_in, rbc_in, zbs_in
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

      SUBROUTINE cut_nescoil(ncoils)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: ncoils
      INTEGER :: m, n, i, j, k, mn
      REAL(rprec) :: u, v, arg, potmax, potmin, level, cop, sip
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: xu,xv, theta_coil, zeta_coil, &
                                                r_coil, z_coil
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: potuv, xnod_in, xnod_ss, xnod_fp
      CHARACTER(len=100) :: s_name
      CHARACTER(len=100) :: l_name
      CHARACTER(len=100) :: c_name
      TYPE(bsc_coil)     :: coil_temp
      TYPE(bsc_rs)       :: rot_mat
      ! Fourier transform pot
      ALLOCATE(xu(nu),xv(nv),potuv(nu,nv))
      potuv = 0.0
      FORALL(i=1:nu) xu(i) = DBLE(i-1)/DBLE(nu-1)
      FORALL(j=1:nv) xv(j) = 0.5*DBLE(j-1)/DBLE(nv-1)
      DO i = 1, nu
         DO j = 1, nv
            DO m = 0, mf
               DO n = -nf, nf
                  arg = pi2*(m*xu(i)+n*xv(j))
                  potuv(i,j) = potuv(i,j) + pot(m,n)*SIN(arg)
               END DO
            END DO
         END DO
      END DO
      ! Do secular parts
      FORALL(i=1:nu) potuv(i,:) = potuv(i,:) - cut*0.5*xu(i)
      FORALL(j=1:nv) potuv(:,j) = potuv(:,j) - cup*0.5*xv(i)
      ! Compute min/max
      potmax = MAXVAL(potuv)
      potmin = MINVAL(potuv)
      IF (ALLOCATED(theta_coil)) DEALLOCATE(theta_coil)
      IF (ALLOCATED(zeta_coil)) DEALLOCATE(zeta_coil)
      ALLOCATE(theta_coil(1),zeta_coil(1))
      m = 0
      ! Find longest paths
      DO i = 1, ncoils
         level = potmin + (potmax-potmin)*(i-0.5)/ncoils
         CALL conrec(potuv,1,nu,1,nv,xu,xv,level,0,n,theta_coil,zeta_coil)
         m = MAX(n,m)
      END DO
      DEALLOCATE(theta_coil,zeta_coil)
      ! Deallocated the coils if allocated
      IF (ALLOCATED(coil_group)) THEN
         DO i = 1, SIZE(coil_group)
            CALL bsc_destroy(coil_group(i))
         END DO
         DEALLOCATE(coil_group)
      END IF
      ! Allocate the coilgroups and helpers
      ALLOCATE(coil_group(ncoils))
      ALLOCATE(theta_coil(m),zeta_coil(m))
      theta_coil = 0; zeta_coil = -3.14
      ! Compute the coils
      DO i = 1, ncoils
         WRITE(l_name,*) 'i = ',i
         WRITE(c_name,'(A,I2.2)') 'MODULAR_COIL_',i
         CALL bsc_construct_coilcoll(coil_group(i),TRIM(c_name),l_name)
         level = potmin + (potmax-potmin)*(i-0.5)/ncoils
         CALL conrec(potuv,1,nu,1,nv,xu,xv,level,1,k,theta_coil,zeta_coil)
         ALLOCATE(r_coil(n),z_coil(n))
         r_coil = 0.0
         z_coil = 0.0
         ! Fourier Transform
         DO j = 1, k
            DO m = 0, mf
               DO n = -nf, nf
                  arg = pi2*(m*theta_coil(j)+n*zeta_coil(j))
                  r_coil(j) = r_coil(j) + cr(m,n)*COS(arg)
                  z_coil(j) = z_coil(j) + cz(m,n)*SIN(arg)
               END DO
            END DO
         END DO
         ALLOCATE(xnod_in(3,k),xnod_ss(3,k),xnod_fp(3,k))
         xnod_in(1,:) = r_coil * COS(zeta_coil*factor)
         xnod_in(2,:) = r_coil * SIN(zeta_coil*factor)
         xnod_in(3,:) = z_coil
         xnod_in(:,k) = xnod_in(:,1)
         ! Now create the first coil
         WRITE(s_name, '(a4,i5.5)') 'ID #', 1
         CALL bsc_construct_coil(coil_temp,'fil_loop',s_name,'',one,xnod_in(1:3,1:k))
         CALL bsc_append(coil_group(i),coil_temp)
         ! Now create the stellarator symmetric coil
         xnod_ss(1,:) = r_coil * COS((1.0-zeta_coil)*factor)
         xnod_ss(2,:) = r_coil * SIN((1.0-zeta_coil)*factor)
         xnod_ss(3,:) =-z_coil
         xnod_ss(:,2:k) = xnod_ss(:,k-1:1:-1)
         xnod_ss(:,1) = xnod_ss(:,k)
         CALL bsc_construct_coil(coil_temp,'fil_loop',s_name,'',one,xnod_ss(1:3,1:ns))
         CALL bsc_append(coil_group(i),coil_temp)
         DO j = 2, np
            cop  = cos((j-1)*factor)
            sip  = sin((j-1)*factor)
            xnod_fp(1,:) = xnod_in(1,:)*cop - xnod_in(2,:)*sip
            xnod_fp(2,:) = xnod_in(2,:)*cop + xnod_in(1,:)*sip
            xnod_fp(3,:) = xnod_in(3,:)
            CALL bsc_construct_coil(coil_temp,'fil_loop',s_name,'',one,xnod_fp(1:3,1:ns))
            CALL bsc_append(coil_group(i),coil_temp)
            xnod_fp(1,:) = xnod_ss(1,:)*cop - xnod_ss(2,:)*sip
            xnod_fp(2,:) = xnod_ss(2,:)*cop + xnod_ss(1,:)*sip
            xnod_fp(3,:) = xnod_ss(3,:)
            CALL bsc_construct_coil(coil_temp,'fil_loop',s_name,'',one,xnod_fp(1:3,1:ns))
            CALL bsc_append(coil_group(i),coil_temp)
         END DO
         DEALLOCATE(xnod_in,xnod_ss,xnod_fp)
         DEALLOCATE(r_coil,z_coil)
      END DO
      DEALLOCATE(xu,xv,potuv)
      RETURN
      END SUBROUTINE cut_nescoil

!-----------------------------------------------------------------------
!     End Module
!-----------------------------------------------------------------------
      END MODULE nescoil_stellopt_mod