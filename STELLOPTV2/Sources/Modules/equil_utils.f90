!-----------------------------------------------------------------------
!     Module:        equil_utils
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          06/08/2012
!     Description:   This module contains the various utilies STELLOPT
!                    utilizes to calculate equilibrium values.
!-----------------------------------------------------------------------
      MODULE equil_utils
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE stellopt_runtime
      USE stellopt_vars
      USE stellopt_targets
      USE equil_vals
      USE read_boozer_mod, ONLY: read_boozer_file, read_boozer_deallocate,&
                                 write_boozer_nc, write_boozer_bin
      USE EZspline_obj
      USE EZspline
      USE stel_tools
      
!-----------------------------------------------------------------------
!     Module Variables
!         
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER     :: domain_flag, nfit_targs, nfit_coefs
      REAL(rprec) :: R_target, PHI_target, Z_target, SUM_target
      REAL(rprec), ALLOCATABLE :: fit_targs(:,:), fit_coefs(:)
      CHARACTER(LEN=256) :: fit_type
      TYPE(EZspline1_r8) :: prof_spl,Bhat_spl,L2_spl
      DOUBLE PRECISION, PARAMETER :: delta_TEM = 0.9999
      DOUBLE PRECISION :: lam_TEM
!-----------------------------------------------------------------------
!     Subroutines
!         get_equil_s:     Returns s given R, PHI(degrees), Z
!         profile_norm:    Returns norm of a profile function
!         get_equil_iota:  Returns iota (and u) given R, PHI(degrees), Z
!         get_equil_jdotb: Returns <j*B> given radial coordiante
!         get_equil_p:     Returns p (and u) given R, PHI(degrees), Z
!         get_equil_phi:   Returns phi (and u) given R, PHI(degrees), Z (electrostatic potential)
!         get_equil_ne:    Returns ne (and u) given R, PHI(degrees), Z
!         get_equil_te:    Returns te (and u) given R, PHI(degrees), Z
!         get_equil_ti:    Returns ti (and u) given R, PHI(degrees), Z
!         mntouv:          Transforms to real space
!         j_star:          Computes trapped branch of jstar
!         
!-----------------------------------------------------------------------

      
      CONTAINS
      
      SUBROUTINE eval_prof_spline(nsp,xknots,yvals,xval,fval,ier,fpval)
      IMPLICIT NONE
      INTEGER, INTENT(in)        :: nsp
      REAL(rprec), INTENT(in)    :: xknots(nsp), yvals(nsp)
      REAL(rprec), INTENT(in)    :: xval
      REAL(rprec), INTENT(out)   :: fval
      REAL(rprec), INTENT(out), OPTIONAL   :: fpval
      INTEGER, INTENT(inout)     :: ier
      INTEGER, SAVE :: nsp_old
      REAL(rprec),ALLOCATABLE, SAVE :: xknots_old(:), yvals_old(:)
      fval = 0
      IF (ier < 0) RETURN
      IF (.not. ALLOCATED(xknots_old)) THEN
         ALLOCATE(xknots_old(nsp))
         xknots_old(1:nsp) = xknots(1:nsp)
      END IF
      IF (.not. ALLOCATED(yvals_old)) THEN
         ALLOCATE(yvals_old(nsp))
         yvals_old(1:nsp) = yvals(1:nsp)
      END IF
      IF (nsp_old /= nsp .or. ANY(xknots .ne. xknots_old) .or. ANY(yvals .ne. yvals_old)) THEN
         nsp_old = nsp
         DEALLOCATE(xknots_old); ALLOCATE(xknots_old(nsp))
         DEALLOCATE(yvals_old); ALLOCATE(yvals_old(nsp))
         xknots_old(1:nsp) = xknots(1:nsp)
         yvals_old(1:nsp)  = yvals(1:nsp)
         IF (EZspline_allocated(prof_spl)) CALL EZspline_free(prof_spl,ier)
         CALL EZspline_init(prof_spl,nsp,bcs0,ier)
         IF (ier .ne. 0) RETURN
         prof_spl%x1 = xknots(1:nsp)
         prof_spl%isHermite = 1
         CALL EZspline_setup(prof_spl,yvals(1:nsp),ier)
         IF (ier .ne. 0) RETURN
         CALL EZspline_isInDomain(prof_spl,xval,ier)
         IF (ier .ne. 0) RETURN
         CALL EZspline_interp(prof_spl,xval,fval,ier)
      ELSE
         CALL EZspline_isInDomain(prof_spl,xval,ier)
         IF (ier .ne. 0) RETURN
         CALL EZspline_interp(prof_spl,xval,fval,ier)
      END IF
      IF (PRESENT(fpval)) THEN
         CALL EZspline_derivative(prof_spl,1,xval,fpval,ier)
         IF (ier .ne. 0) RETURN
      END IF
      RETURN
      END SUBROUTINE eval_prof_spline
      
      SUBROUTINE get_equil_iota(s_val,val,ier)
      IMPLICIT NONE
      REAL(rprec), INTENT(inout) ::  s_val
      REAL(rprec), INTENT(out)   ::  val
      INTEGER, INTENT(inout)     ::  ier
      IF (ier < 0) RETURN
      IF (EZspline_allocated(iota_spl)) THEN
         CALL EZspline_isInDomain(iota_spl,s_val,ier)
         IF (ier .ne. 0) RETURN
         CALL EZspline_interp(iota_spl,s_val,val,ier)
      ELSE
         ier = -1
      END IF
      RETURN
      END SUBROUTINE get_equil_iota
      
      SUBROUTINE get_equil_p(s_val,val,ier,pval)
      IMPLICIT NONE
      REAL(rprec), INTENT(inout) ::  s_val
      REAL(rprec), INTENT(out)   ::  val
      REAL(rprec), INTENT(out), OPTIONAL :: pval ! P'=dp/ds
      INTEGER, INTENT(inout)     ::  ier
      IF (ier < 0) RETURN
      IF (EZspline_allocated(pres_spl)) THEN
         CALL EZspline_isInDomain(pres_spl,s_val,ier)
         IF (ier .ne. 0) RETURN
         CALL EZspline_interp(pres_spl,s_val,val,ier)
         IF (PRESENT(pval)) CALL EZspline_derivative(pres_spl,1,s_val,pval,ier)
      ELSE
         ier = -1
      END IF
      RETURN
      END SUBROUTINE get_equil_p
      
      SUBROUTINE get_equil_volume(s_val,val,ier,pval)
      IMPLICIT NONE
      REAL(rprec), INTENT(inout) ::  s_val
      REAL(rprec), INTENT(out)   ::  val
      REAL(rprec), INTENT(out), OPTIONAL :: pval ! V'=dV/ds
      INTEGER, INTENT(inout)     ::  ier
      IF (ier < 0) RETURN
      IF (EZspline_allocated(V_spl)) THEN
         CALL EZspline_isInDomain(V_spl,s_val,ier)
         IF (ier .ne. 0) RETURN
         CALL EZspline_interp(V_spl,s_val,val,ier)
         IF (PRESENT(pval)) CALL EZspline_derivative(V_spl,1,s_val,pval,ier)
      ELSE
         ier = -1
      END IF
      RETURN
      END SUBROUTINE get_equil_volume
      
      SUBROUTINE get_equil_jdotb(s_val,val,ier)
      IMPLICIT NONE
      REAL(rprec), INTENT(inout) ::  s_val
      REAL(rprec), INTENT(out)   ::  val
      INTEGER, INTENT(inout)     ::  ier
      IF (ier < 0) RETURN
      IF (EZspline_allocated(jdotb_spl)) THEN
         CALL EZspline_isInDomain(jdotb_spl,s_val,ier)
         IF (ier .ne. 0) RETURN
         CALL EZspline_interp(jdotb_spl,s_val,val,ier)
      ELSE
         ier = -1
      END IF
      RETURN
      END SUBROUTINE get_equil_jdotb
      
      SUBROUTINE get_equil_jcurv(s_val,val,ier)
      IMPLICIT NONE
      REAL(rprec), INTENT(inout) ::  s_val
      REAL(rprec), INTENT(out)   ::  val
      INTEGER, INTENT(inout)     ::  ier
      IF (ier < 0) RETURN
      IF (EZspline_allocated(jcurv_spl)) THEN
         CALL EZspline_isInDomain(jcurv_spl,s_val,ier)
         IF (ier .ne. 0) RETURN
         CALL EZspline_interp(jcurv_spl,s_val,val,ier)
      ELSE
         ier = -1
      END IF
      RETURN
      END SUBROUTINE get_equil_jcurv
      
      SUBROUTINE get_equil_omega(s_val,val,ier)
      IMPLICIT NONE
      REAL(rprec), INTENT(inout) ::  s_val
      REAL(rprec), INTENT(out)   ::  val
      INTEGER, INTENT(inout)     ::  ier
      IF (ier < 0) RETURN
      IF (EZspline_allocated(omega_spl)) THEN
         CALL EZspline_isInDomain(omega_spl,s_val,ier)
         IF (ier .ne. 0) RETURN
         CALL EZspline_interp(omega_spl,s_val,val,ier)
      ELSE
         ier = -1
      END IF
      RETURN
      END SUBROUTINE get_equil_omega
      
      SUBROUTINE get_equil_nustar(s_val,val,ier)
      IMPLICIT NONE
      REAL(rprec), INTENT(inout) ::  s_val
      REAL(rprec), INTENT(out)   ::  val
      INTEGER, INTENT(inout)     ::  ier
      IF (ier < 0) RETURN
      IF (EZspline_allocated(nustar_spl)) THEN
         CALL EZspline_isInDomain(nustar_spl,s_val,ier)
         IF (ier .ne. 0) RETURN
         CALL EZspline_interp(nustar_spl,s_val,val,ier)
      ELSE
         ier = -1
      END IF
      RETURN
      END SUBROUTINE get_equil_nustar
      
      SUBROUTINE eval_prof_stel(s_val,type,val,ncoefs,coefs,ier,spl_obj)
      IMPLICIT NONE
      REAL(rprec),      INTENT(in)             :: s_val
      CHARACTER(LEN=*), INTENT(in)             :: type
      REAL(rprec),      INTENT(inout)          :: val
      INTEGER,          INTENT(in)             :: ncoefs
      REAL(rprec),      INTENT(inout)          :: coefs(ncoefs)
      INTEGER,          INTENT(inout)          ::  ier
      TYPE(EZspline1_r8),OPTIONAL              :: spl_obj
      INTEGER :: i
      REAL(rprec) :: x0,x1, x2, h, x3, xp
      REAL(rprec), PARAMETER :: one = 1.0_rprec
      REAL(rprec), DIMENSION(10), PARAMETER :: glx = (/                       &
     &   0.01304673574141414, 0.06746831665550774, 0.1602952158504878,         &
     &   0.2833023029353764, 0.4255628305091844, 0.5744371694908156,           &
     &   0.7166976970646236, 0.8397047841495122, 0.9325316833444923,           &
     &   0.9869532642585859 /)
      REAL(rprec), DIMENSION(10), PARAMETER :: glw = (/                       &
     &   0.03333567215434407, 0.0747256745752903, 0.1095431812579910,          &
     &   0.1346333596549982, 0.1477621123573764, 0.1477621123573764,           &
     &   0.1346333596549982, 0.1095431812579910, 0.0747256745752903,           &
     &   0.03333567215434407 /)
      val = 0
      IF (ier < 0) RETURN
      CALL tolower(type)
      SELECT CASE (type)
         CASE ('gauss_trunc')
            val = (coefs(1)/(one - EXP(-(one/coefs(2))**2)))*&
                  (EXP(-(s_val/coefs(2))**2)-EXP(-(one/coefs(2))**2))
         CASE ('gauss_trunc_offset')
            val = coefs(3)+(coefs(1)/(one - EXP(-(one/coefs(2))**2)))*&
                  (EXP(-(s_val/coefs(2))**2)-EXP(-(one/coefs(2))**2))
         CASE ('two_lorentz')
            val = coefs(1)*(coefs(2)*(one/(one+(  s_val/coefs(3)**2)**coefs(4))**coefs(5) &
                                       -one/(one+(one/coefs(3)**2)**coefs(4))**coefs(5))/   &
                                   (one-one/(one+(one/coefs(3)**2)**coefs(4))**coefs(5))+   &
                   (one-coefs(2))*(one/(one+(  s_val/coefs(6)**2)**coefs(7))**coefs(8)     &
                                   -one/(one+(one/coefs(6)**2)**coefs(7))**coefs(8))/       &
                               (one-one/(one+(one/coefs(6)**2)**coefs(7))**coefs(8)))
         CASE ('two_power')
            val = coefs(1) * (one - s_val**coefs(2))**coefs(3)
         CASE ('two_power_hollow')
            val = s_val * coefs(1) * (one - s_val**coefs(2))**coefs(3)
         CASE ('two_power_offset')
            val = coefs(4) + coefs(1) * (one - s_val**coefs(2))**coefs(3)
         CASE ('power_series')
            DO i = UBOUND(coefs,DIM=1), LBOUND(coefs,DIM=1), -1
               val = s_val*val + coefs(i)
            END DO
         CASE ('power_series_0_boundaries')
            DO i = UBOUND(coefs,DIM=1), LBOUND(coefs,DIM=1), -1
               val = s_val*val + coefs(i)
            END DO
            val = val * s_val * (1 - s_val)
         CASE ('power_series_0i0')
            DO i = UBOUND(coefs,DIM=1), LBOUND(coefs,DIM=1), -1
               val = val*s_val**0.25 + coefs(i)
            END DO
            val = 4*val*s_val*(1-s_val)
         CASE ('power_series_edge0')
            i=UBOUND(coefs,DIM=1)
            !coefs(i) = -SUM(coefs(LBOUND(coefs,DIM=1):i-1))
            val = -SUM(coefs(LBOUND(coefs,DIM=1):i))
            DO i = UBOUND(coefs,DIM=1), LBOUND(coefs,DIM=1), -1
               val = s_val*val + coefs(i)
            END DO
         CASE ('power_series_i') ! dI/ds = a1+2*a2*x+a
            DO i = UBOUND(coefs,DIM=1), LBOUND(coefs,DIM=1), -1
               val = s_val*val + coefs(i)*(i-1) ! OK (1-1)=0
            END DO
            IF (s_val >0) val = val / s_val
         CASE ('power_series_i_edge0') ! dI/ds = a1+2*a2*x+a
            DO i = UBOUND(coefs,DIM=1), LBOUND(coefs,DIM=1), -1
               val = val - coefs(i)*(i-1)
            END DO
            DO i = UBOUND(coefs,DIM=1), LBOUND(coefs,DIM=1), -1
               val = s_val*val + coefs(i)*(i-1) ! OK (1-1)=0
            END DO
            IF (s_val >0) val = val / s_val
         ! For now don't handle spline
         CASE ('spline','akima_spline','akima_spline_ip')
            IF (EZspline_allocated(spl_obj)) THEN
               CALL EZspline_isInDomain(spl_obj,s_val,ier)
               IF (ier .ne. 0) RETURN
               CALL EZspline_interp(spl_obj,s_val,val,ier)
            ELSE
               ier = -1
            END IF
         CASE ('pedestal')
            DO i=15, LBOUND(coefs,1),-1
               val = s_val*val+coefs(i)
            END DO
            i=17
            IF(coefs(i+3) .le. 0.0_rprec) THEN
               coefs(i:i+4) = 0
               coefs(i+3)   = 1.0D+30
            ELSE
               coefs(i+4) = one / (TANH(2*coefs(i+2)/coefs(i+3))-TANH(2*(coefs(i+2)-1)/coefs(i+3)))
            END IF
            val = val + coefs(i+4) * coefs(i+1) * ( TANH(2*(coefs(i+2)-SQRT(s_val))/coefs(i+3) )    &
                                                     -TANH(2*(coefs(i+2)-one)/coefs(i+3)      ) )
         CASE('sum_atan')
            IF (s_val >= 1) THEN
               val = coefs(1)+coefs(2)+coefs(6)+coefs(10)+coefs(14)+coefs(18)
            ELSE
               val = coefs(1) + &
                     (4/pi2) * ( coefs(2) * ATAN(coefs(3)*s_val**coefs(4)/(1-s_val)**coefs(5)) &
                                +coefs(6) * ATAN(coefs(7)*s_val**coefs(8)/(1-s_val)**coefs(9)) &
                                +coefs(10) * ATAN(coefs(11)*s_val**coefs(12)/(1-s_val)**coefs(13)) &
                                +coefs(14) * ATAN(coefs(15)*s_val**coefs(16)/(1-s_val)**coefs(17)) &
                                +coefs(18) * ATAN(coefs(19)*s_val**coefs(20)/(1-s_val)**coefs(21)))
            END IF
         CASE ('bump')
            ! bootj_aux_f(1) : x0 center of bump
            ! bootj_aux_f(2) : height of bump
            ! bootj_aux_f(3) : height of bulk
            x0  = bootj_aux_f(1)
            x1  = 1.0_rprec-2*(1.0_rprec-x0)
            x2  = 1.0_rprec
            x3  = bootj_aux_f(3)
            h   = bootj_aux_f(2)/((x0-x1)*(x0-x2))
            if ((s_val > x1) .and. (s_val < 1.0_rprec)) val = h*(s_val-x1)*(s_val-x2)
            val = val + x3*s_val*(s_val-1)/(-0.25_rprec)
         CASE ('hollow')
            ! coefs(1) : Value at core
            ! coefs(2) : Value at peak
            ! coefs(3) : grid scaling
            xp  = s_val**coefs(3)
            x1  = 4*coefs(2)-3*coefs(1)
            x2  = 2*coefs(1)-4*coefs(2)
            val = coefs(1) + x1*xp + x2*xp*xp
         CASE ('hollow2')
            ! coefs(1) : Value at core
            ! coefs(2) : Value at peak
            ! coefs(3) : grid scaling
            xp  = s_val**coefs(3)
            x1  = 8*coefs(2)-2*coefs(1)
            val = (coefs(1)+x1*xp)*(xp-1)**2
         CASE DEFAULT
            PRINT *,"Error! Unknown profile type in subroutine eval_prof_stel:",type
            STOP
      END SELECT
      RETURN
      END SUBROUTINE eval_prof_stel

      SUBROUTINE setup_prof_spline(spl_obj,l_aux,s_aux,f_aux,ier)
      IMPLICIT NONE
      TYPE(EZspline1_r8), INTENT(inout) :: spl_obj
      INTEGER,            INTENT(inout) :: l_aux
      REAL(rprec),        INTENT(inout) :: s_aux(l_aux)
      REAL(rprec),        INTENT(inout) :: f_aux(l_aux)
      INTEGER,            INTENT(inout) ::  ier
      ier = 0
      IF (EZspline_allocated(spl_obj)) CALL EZspline_free(spl_obj,ier)
      IF (ier /= 0) RETURN
      CALL EZspline_init(spl_obj,l_aux,bcs0,ier)
      IF (ier /= 0) RETURN
      spl_obj%x1        = s_aux(1:l_aux)
      spl_obj%isHermite = 1
      CALL EZspline_setup(spl_obj,f_aux,ier)
      RETURN
      END SUBROUTINE setup_prof_spline
      
      SUBROUTINE get_equil_phi(s_val,type,val,ier,pval)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) ::  s_val
      CHARACTER(LEN=*), INTENT(in)   :: type
      REAL(rprec), INTENT(inout)   ::  val
      REAL(rprec), INTENT(inout),OPTIONAL   ::  pval
      INTEGER, INTENT(inout)     ::  ier
      INTEGER :: i
      REAL(rprec), PARAMETER :: one = 1.0_rprec
      IF (ier < 0) RETURN
      CALL tolower(type)
      SELECT CASE (type)
         CASE ('spline','akima_spline','akima_spline_ip')
            CALL eval_prof_stel(s_val,type,val,20,phi_aux_f(1:20),ier,phi_spl)
            IF (PRESENT(pval)) THEN
               CALL EZspline_derivative(phi_spl,1,s_val,pval,ier)
               IF (s_val>one) pval=0
            END IF
         CASE DEFAULT
            CALL eval_prof_stel(s_val,type,val,20,phi_aux_f(1:20),ier)
            IF (PRESENT(pval)) THEN
               CALL eval_prof_stel(s_val+1E-6,type,pval,20,phi_aux_f(1:20),ier)
               pval = (pval-val)*1E6
               IF (s_val>=one) pval=0
            END IF
      END SELECT
      RETURN
      END SUBROUTINE get_equil_phi
      
      SUBROUTINE get_equil_E(r_val,phi_val,z_val,er,ez,ier)
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  r_val
      REAL(rprec), INTENT(in)    ::  phi_val
      REAL(rprec), INTENT(in)    ::  z_val
      REAL(rprec), INTENT(out)   ::  er
      REAL(rprec), INTENT(out)   ::  ez
      INTEGER, INTENT(inout)     ::  ier
      REAL(rprec) :: s_val, u_val, v_val, phi_prime, phi2,R1,Z1
      REAL(rprec) :: R_grad(3), Z_grad(3)
      IF (ier < 0) RETURN
      CALL get_equil_s(r_val,phi_val,z_val,s_val,ier,u_val)
      IF (ier < 0) RETURN
      v_val = MOD(phi_val,pi2/nfp)*nfp
      CALL get_equil_RZ(s_val,u_val,v_val,R1,Z1,ier,&
            R_GRAD=R_grad,Z_GRAD=Z_grad)
      IF (ier == 0) THEN
         CALL get_equil_phi(s_val,phi_type,phi2,ier,phi_prime)
         er = R_grad(3)*phi_prime
         ez = Z_grad(3)*phi_prime
      ELSE
        er  = 0
        ez  = 0
      END IF
      RETURN
      END SUBROUTINE get_equil_E
      
      SUBROUTINE get_equil_ne(s_val,type,val,ier)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) ::  s_val
      CHARACTER(LEN=*), INTENT(in)   :: type
      REAL(rprec), INTENT(inout)   ::  val
      INTEGER, INTENT(inout)     ::  ier
      INTEGER :: i
      REAL(rprec), PARAMETER :: one = 1.0_rprec
      IF (ier < 0) RETURN
      CALL tolower(type)
      SELECT CASE (type)
         CASE ('spline','akima_spline','akima_spline_ip')
            CALL eval_prof_stel(s_val,type,val,21,ne_opt(0:20),ier,ne_spl)
         CASE DEFAULT
            CALL eval_prof_stel(s_val,type,val,21,ne_opt(0:20),ier)
      END SELECT
      val = val * ne_norm
      RETURN
      END SUBROUTINE get_equil_ne
      
      SUBROUTINE get_equil_te(s_val,type,val,ier)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) ::  s_val
      CHARACTER(LEN=*), INTENT(in)   :: type
      REAL(rprec), INTENT(inout)   ::  val
      INTEGER, INTENT(inout)     ::  ier
      INTEGER :: i
      REAL(rprec), PARAMETER :: one = 1.0_rprec
      IF (ier < 0) RETURN
      CALL tolower(type)
      SELECT CASE (type)
         CASE ('spline','akima_spline','akima_spline_ip')
            CALL eval_prof_stel(s_val,type,val,21,te_opt(0:20),ier,te_spl)
            !IF (EZspline_allocated(te_spl)) THEN
            !   CALL EZspline_isInDomain(te_spl,s_val,ier)
            !   IF (ier .ne. 0) RETURN
            !   CALL EZspline_interp(te_spl,s_val,val,ier)
            !ELSE
            !   ier = -1
            !END IF
         CASE DEFAULT
            CALL eval_prof_stel(s_val,type,val,21,te_opt(0:20),ier)
      END SELECT
      RETURN
      END SUBROUTINE get_equil_te
      
      SUBROUTINE get_equil_ti(s_val,type,val,ier)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) ::  s_val
      CHARACTER(LEN=*), INTENT(in)   :: type
      REAL(rprec), INTENT(inout)   ::  val
      INTEGER, INTENT(inout)     ::  ier
      INTEGER :: i
      REAL(rprec), PARAMETER :: one = 1.0_rprec
      IF (ier < 0) RETURN
      CALL tolower(type)
      SELECT CASE (type)
         CASE ('spline','akima_spline','akima_spline_ip')
            CALL eval_prof_stel(s_val,type,val,21,ti_opt(0:20),ier,ti_spl)
            !IF (EZspline_allocated(ti_spl)) THEN
            !   CALL EZspline_isInDomain(ti_spl,s_val,ier)
            !   IF (ier .ne. 0) RETURN
            !   CALL EZspline_interp(ti_spl,s_val,val,ier)
            !ELSE
            !   ier = -1
            !END IF
         CASE DEFAULT
            CALL eval_prof_stel(s_val,type,val,21,ti_opt(0:20),ier)
      END SELECT
      RETURN
      END SUBROUTINE get_equil_ti
      
      SUBROUTINE get_equil_emis_xics(s_val,type,val,ier)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) ::  s_val
      CHARACTER(LEN=*), INTENT(in)   :: type
      REAL(rprec), INTENT(inout)   ::  val
      INTEGER, INTENT(inout)     ::  ier
      INTEGER :: i
      REAL(rprec), PARAMETER :: one = 1.0_rprec
      IF (ier < 0) RETURN
      CALL tolower(type)
      SELECT CASE (type)
         CASE ('spline','akima_spline','akima_spline_ip')
            CALL eval_prof_stel(s_val,type,val,20,emis_xics_f(1:20),ier,emis_xics_spl)
         CASE DEFAULT
            CALL eval_prof_stel(s_val,type,val,20,emis_xics_f(1:20),ier)
      END SELECT
      RETURN
      END SUBROUTINE get_equil_emis_xics
      
      SUBROUTINE get_equil_zeff(s_val,type,val,ier)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) ::  s_val
      CHARACTER(LEN=*), INTENT(in)   :: type
      REAL(rprec), INTENT(inout)   ::  val
      INTEGER, INTENT(inout)     ::  ier
      INTEGER :: i
      REAL(rprec), PARAMETER :: one = 1.0_rprec
      IF (ier < 0) RETURN
      CALL tolower(type)
      SELECT CASE (type)
         CASE ('spline','akima_spline','akima_spline_ip')
            CALL eval_prof_stel(s_val,type,val,21,zeff_opt(0:20),ier,zeff_spl)
            !IF (EZspline_allocated(zeff_spl)) THEN
            !   CALL EZspline_isInDomain(zeff_spl,s_val,ier)
            !   IF (ier .ne. 0) RETURN
            !   CALL EZspline_interp(zeff_spl,s_val,val,ier)
            !ELSE
            !   ier = -1
            !END IF
         CASE DEFAULT
            CALL eval_prof_stel(s_val,type,val,21,zeff_opt(0:20),ier)
      END SELECT
      RETURN
      END SUBROUTINE get_equil_zeff

      SUBROUTINE fcn_linene(s,u,v,dx,dy,dz,fval,ier)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: s,u,v,dx,dy,dz
      REAL(rprec), INTENT(out) :: fval
      INTEGER, INTENT(inout) :: ier
      fval = 0
      IF (s>1) RETURN
      CALL get_equil_ne(s,TRIM(ne_type),fval,ier)
      fval = MAX(fval,0.0)*sqrt(dx*dx+dy*dy+dz*dz)
      RETURN
      END SUBROUTINE fcn_linene

      SUBROUTINE fcn_linete(s,u,v,dx,dy,dz,fval,ier)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: s,u,v,dx,dy,dz
      REAL(rprec), INTENT(out) :: fval
      INTEGER, INTENT(inout) :: ier
      fval = 0
      IF (s>1) RETURN
      CALL get_equil_te(s,TRIM(te_type),fval,ier)
      fval = MAX(fval,0.0)*sqrt(dx*dx+dy*dy+dz*dz)
      RETURN
      END SUBROUTINE fcn_linete

      SUBROUTINE fcn_lineti(s,u,v,dx,dy,dz,fval,ier)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: s,u,v,dx,dy,dz
      REAL(rprec), INTENT(out) :: fval
      INTEGER, INTENT(inout) :: ier
      fval = 0
      IF (s>1) RETURN
      CALL get_equil_ti(s,TRIM(ti_type),fval,ier)
      fval = MAX(fval,0.0)*sqrt(dx*dx+dy*dy+dz*dz)
      RETURN
      END SUBROUTINE fcn_lineti

      SUBROUTINE fcn_xics_bright(s,u,v,dx,dy,dz,fval,ier)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: s,u,v,dx,dy,dz
      REAL(rprec), INTENT(out) :: fval
      INTEGER, INTENT(inout) :: ier
      fval = 0
      IF (s>1) RETURN
      CALL get_equil_emis_xics(s,TRIM(emis_xics_type),fval,ier)
      fval = MAX(fval,0.0)*sqrt(dx*dx+dy*dy+dz*dz)
      RETURN
      END SUBROUTINE fcn_xics_bright

      SUBROUTINE fcn_xics(s,u,v,dx,dy,dz,fval,ier)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: s,u,v,dx,dy,dz
      REAL(rprec), INTENT(out) :: fval
      REAL(rprec) :: f1, f2
      INTEGER, INTENT(inout) :: ier
      fval = 0
      IF (s>1) RETURN
      CALL get_equil_ti(s,TRIM(ti_type),f1,ier)
      CALL get_equil_emis_xics(s,TRIM(emis_xics_type),f2,ier)
      fval = MAX(f1*f2,0.0)*sqrt(dx*dx+dy*dy+dz*dz)
      RETURN
      END SUBROUTINE fcn_xics

      SUBROUTINE fcn_xics_w3(s,u,v,dx,dy,dz,fval,ier)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: s,u,v,dx,dy,dz
      REAL(rprec), INTENT(out) :: fval
      INTEGER, INTENT(inout) :: ier
      REAL(rprec) :: w, te
      ! natural log fit from VAINSHTEIN_2011 and MARCHUCK_2004
      REAL(rprec), PARAMETER :: c3 = -7.5306337062E-02
      REAL(rprec), PARAMETER :: c2 =  1.7282006139E+00
      REAL(rprec), PARAMETER :: c1 = -1.4811899942E+01
      !REAL(rprec), PARAMETER :: c0 =  3.9207616355E+01  ! From equations
      REAL(rprec), PARAMETER :: c0 =  3.8707616355E+01  ! From Fitting
      REAL(rprec), PARAMETER :: temin =  1000.0  ! Comes from table (sort of)
      fval = 0
      IF (s>1) RETURN
      CALL get_equil_emis_xics(s,TRIM(emis_xics_type),w,ier)
      w = MAX(w,0.0)
      CALL get_equil_te(s,TRIM(te_type),te,ier)
      te = MAX(te,temin)
      te = LOG(te) ! interpolation in log of te
      fval = w*EXP((((c3*te)+c2)*te+c1)*te+c0)
      fval = MAX(fval,0.0)*SQRT(dx*dx+dy*dy+dz*dz)
      RETURN
      END SUBROUTINE fcn_xics_w3

      SUBROUTINE fcn_xics_v(s,u,v,dx,dy,dz,fval,ier)
      IMPLICIT NONE
      REAL(rprec), INTENT(inout) :: s,u,v,dx,dy,dz
      REAL(rprec), INTENT(out) :: fval
      REAL(rprec) :: phi_val, phi_prime,Bav,Bavsq, &
                     rho,vp,grho,grho2,&
                     br,bp,bx,by,bz,modb
      REAL(rprec) :: nhat(3), uperp(3)
      INTEGER, INTENT(inout) :: ier
      fval = 0
      IF (s>1) RETURN
      CALL get_equil_phi(s,phi_type,phi_val,ier,phi_prime)
      !CALL get_equil_Bav(s,Bav,Bavsq,ier) !<B>,<B^2>
      CALL get_equil_rho(s,rho,vp,grho,grho2,ier)
      ! To make these lines work we need to implement
      ! passing u and v
      CALL get_equil_nhat(s,u,v,nhat,ier)
      CALL get_equil_Bcylsuv(s,u,v,br,bp,bz,ier,modb)
      phi_prime = phi_prime * 2 * rho * grho !dphi/drho=dphi/ds*ds/drho
      nhat = nhat / SQRT(SUM(nhat*nhat))
      bx = br * cos(v) - bp * sin(v)
      by = br * sin(v) + bp * cos(v)
      uperp(1) = by*nhat(3)-bz*nhat(2)
      uperp(2) = bz*nhat(1)-bx*nhat(1)
      uperp(3) = bx*nhat(2)-by*nhat(1)
      uperp    = uperp*phi_prime/(modb*modb)
      fval = uperp(1)*dx+uperp(2)*dy+uperp(3)*dz
      RETURN
      END SUBROUTINE fcn_xics_v

      SUBROUTINE fcn_sxr(s,u,v,dx,dy,dz,fval,ier)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: s,u,v,dx,dy,dz
      REAL(rprec), INTENT(out) :: fval
      INTEGER, INTENT(inout) :: ier
      REAL(rprec) :: ne_val,te_val,zeff_val
      fval = 0
      IF (s>1) RETURN
      CALL get_equil_ne(s,TRIM(ne_type),ne_val,ier)
      CALL get_equil_te(s,TRIM(ne_type),te_val,ier)
      CALL get_equil_zeff(s,TRIM(ne_type),zeff_val,ier)
      IF (abs(te_val) > 0) THEN
         fval = zeff_val*ne_val*ne_val*sqrt(te_val)
      ELSE
         fval = 0
      END IF
      fval = fval*sqrt(dx*dx+dy*dy+dz*dz)
      RETURN
      END SUBROUTINE fcn_sxr
      
      SUBROUTINE line_int_faraday(r1,r2,val)
      IMPLICIT NONE
      REAL(rprec), INTENT(in)   :: r1(3), r2(3)
      REAL(rprec), INTENT(out)  :: val
      INTEGER     :: i, j, k, ier
      REAL(rprec) :: x1, y1, z1, x2, y2, z2, dx, dy, dz, xp, yp, zp, int_fac
      REAL(rprec) :: xp1, yp1, zp1, delf, delt, rp, phip,s, dx2, dy2, dz2
      REAL(rprec) :: ne_val, bx, by, bz
      INTEGER, PARAMETER ::  nsteps=50
      INTEGER, PARAMETER :: nop=3
      INTEGER, PARAMETER :: int_step=2
      real(rprec), dimension(nop), parameter :: ci=(/1._rprec/6._rprec,2._rprec/3._rprec,1._rprec/6._rprec/)
      phip = r1(2);
      s    = r2(2);
      !IF (phip < 0) phip = phip+2*pi
      !IF (s < 0) s = s+2*pi
      !phip = MOD(phip,pi2/nfp)*nfp
      !s    = MOD(s,pi2/nfp)*nfp
      ier = 0; ne_val = 0.0
      x1 = r1(1)*cos(phip); x2 = r2(1)*cos(s)
      y1 = r1(1)*sin(phip); y2 = r2(1)*sin(s)
      z1 = r1(3);            z2 = r2(3);
      dx = (x2-x1)/nsteps
      dy = (y2-y1)/nsteps
      dz = (z2-z1)/nsteps
      s=0; phip=0
      DO i = 1, nsteps ! Get first boundary point
         xp = x1+dx*(i-1)
         yp = y1+dy*(i-1)
         zp = z1+dz*(i-1)
         rp = sqrt(xp*xp+yp*yp)
         phip = ATAN2(yp,xp)
         IF (phip < 0) phip = phip+2*pi
         CALL get_equil_s(rp,phip,zp,s,ier)
         IF (s <= 1.0_rprec) THEN
            x1 = xp-dx
            y1 = yp-dy
            z1 = zp-dz
            EXIT
         END IF
      END DO
      DO i = 1, nsteps ! Get second boundary point
         xp = x2-dx*(i-1)
         yp = y2-dy*(i-1)
         zp = z2-dz*(i-1)
         IF ((xp == x1) .and. (yp == y1) .and. (zp == z1)) THEN
            val = 0.0_rprec
            RETURN
         END IF
         rp = sqrt(xp*xp+yp*yp)
         phip = ATAN2(yp,xp)
         IF (phip < 0) phip = phip+2*pi
         CALL get_equil_s(rp,phip,zp,s,ier)
         IF (s <= 1.0_rprec) THEN
            x2 = xp+dx
            y2 = yp+dy
            z2 = zp+dz
            dx = (x2-x1)/nsteps
            dy = (y2-y1)/nsteps
            dz = (z2-z1)/nsteps
            EXIT
         END IF
      END DO
      val =0
      delt = 1.0_rprec/REAL(nop-1)
      int_fac = 1.0_rprec/REAL(int_step)
      DO i = 1, nsteps
         DO j = 1, int_step
            xp = x1+dx*(i-1)+(j-1)*int_fac*dx
            yp = y1+dy*(i-1)+(j-1)*int_fac*dy
            zp = z1+dz*(i-1)+(j-1)*int_fac*dz
            delf = 0.0_rprec
            dx2 = x1+dx*(i-1)+(j)*int_fac*dx - xp
            dy2 = y1+dy*(i-1)+(j)*int_fac*dy - yp
            dz2 = z1+dz*(i-1)+(j)*int_fac*dz - zp
            DO k = 1, nop
               rp = sqrt(xp*xp+yp*yp)
               phip = ATAN2(yp,xp)
               IF (phip < 0) phip = phip+2*pi
               CALL get_equil_s(rp,phip,zp,s,ier)
               IF (s <= 1.0_rprec .and. s >= 0.0_rprec) THEN
                  CALL get_equil_ne(s,TRIM(ne_type),ne_val,ier)
                  IF (ier /= 0) ne_val = 0.0_rprec
                  CALL get_equil_B(rp,phip,zp,bx,by,bz,ier)
                  IF (ier /= 0) THEN
                     bx=0.0_rprec
                     by=0.0_rprec
                     bz=0.0_rprec
                  END IF
               ELSE
                  ne_val = 0.0_rprec
                  bx = 0.0_rprec
                  by = 0.0_rprec
                  bz = 0.0_rprec
               END IF
               !WRITE(27,*) rp, phip,zp,ne_val,s,delf
               delf = delf + ci(k)*ne_val*ABS(bx*dx2+by*dy2+bz*dz2)
               xp   = xp + k*dx2*delt
               yp   = yp + k*dy2*delt
               zp   = zp + k*dz2*delt
            END DO
            val = val + delf
         END DO
      END DO
      RETURN
      END SUBROUTINE line_int_faraday
         
         
      SUBROUTINE mntouv(k1,k,mnmax,nu,nv,xu,xv,fmn,xm,xn,f,signs,calc_trig)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: k1
      INTEGER, INTENT(in) :: k
      INTEGER, INTENT(in) :: mnmax
      INTEGER, INTENT(in) :: nu
      INTEGER, INTENT(in) :: nv
      REAL(rprec), INTENT(in) :: xu(1:nu)
      REAL(rprec), INTENT(in) :: xv(1:nv)           
      REAL(rprec), INTENT(in) :: fmn(1:mnmax,k1:k)
      INTEGER, INTENT(in) :: xm(1:mnmax)
      INTEGER, INTENT(in) :: xn(1:mnmax)
      REAL(rprec), INTENT(inout) :: f(1:nu,1:nv,k1:k)
      INTEGER, INTENT(in) :: signs
      INTEGER, INTENT(in) :: calc_trig
      INTEGER     :: mn, i, ier, ik
      REAL(rprec) :: xm_temp(1:mnmax,1)
      REAL(rprec) :: xn_temp(1:mnmax,1)
      REAL(rprec) :: pi2_l
      REAL(rprec) :: mt(1:mnmax,1:nu)
      REAL(rprec) :: nz(1:mnmax,1:nv)
      REAL(rprec) :: fmn_temp(1:mnmax,1:nu)
      REAL(rprec) :: xu_temp(1,1:nu)
      REAL(rprec) :: xv_temp(1,1:nv)
      REAL(rprec) :: fmn_help(1:mnmax)
      REAL(rprec), ALLOCATABLE, SAVE :: cosmt(:,:)
      REAL(rprec), ALLOCATABLE, SAVE :: sinmt(:,:)
      REAL(rprec), ALLOCATABLE, SAVE :: cosnz(:,:)
      REAL(rprec), ALLOCATABLE, SAVE :: sinnz(:,:)
      pi2_l = 8 * ATAN(1.)
      IF (calc_trig == 1) THEN
         IF (ALLOCATED(cosmt)) DEALLOCATE(cosmt)
         IF (ALLOCATED(sinmt)) DEALLOCATE(sinmt)
         IF (ALLOCATED(cosnz)) DEALLOCATE(cosnz)
         IF (ALLOCATED(sinnz)) DEALLOCATE(sinnz)
         ALLOCATE(cosmt(1:mnmax,1:nu),sinmt(1:mnmax,1:nu),&
                  cosnz(1:mnmax,1:nv),sinnz(1:mnmax,1:nv),STAT=ier)
         FORALL(i=1:mnmax) xm_temp(i,1)=REAL(xm(i))
         FORALL(i=1:mnmax) xn_temp(i,1)=REAL(xn(i))
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
            IF (.not.ANY(fmn(:,ik) /= 0)) CYCLE
            FORALL(mn=1:mnmax) fmn_help(mn)=fmn(mn,ik)
            fmn_temp=SPREAD(fmn_help,2,nu)
            f(1:nu,1:nv,ik) = f(1:nu,1:nv,ik)  + MATMUL(TRANSPOSE((fmn_temp*cosmt)),cosnz) &
                                   - MATMUL(TRANSPOSE((fmn_temp*sinmt)),sinnz)
         END DO
      ELSE IF (SIGNS == 1) THEN
         DO ik = k1,k
            IF (.not.ANY(fmn(:,ik) /= 0)) CYCLE
            FORALL(mn=1:mnmax) fmn_help(mn)=fmn(mn,ik)
            fmn_temp=SPREAD(fmn_help,2,nu)
            f(1:nu,1:nv,ik) = f(1:nu,1:nv,ik) + MATMUL(TRANSPOSE((fmn_temp*sinmt)),cosnz) &
                                  + MATMUL(TRANSPOSE((fmn_temp*cosmt)),sinnz)
         END DO
      END IF
      END SUBROUTINE mntouv

      
      REAL(rprec) FUNCTION profile_norm(x,prof_type)
      REAL(rprec), INTENT(in) :: x(:)
      CHARACTER(LEN=*), INTENT(in)   :: prof_type
      INTEGER :: ik
      CALL tolower(prof_type)
      profile_norm = 0.0_rprec
      SELECT CASE (prof_type)
         CASE ('two_power','two_power_hollow','two_power_offset','two_lorentz','gauss_trunc', &
               'gauss_trunc_offset','sum_atan','pedestal','bump','hollow','hollow2')
            profile_norm = 0.0_rprec  ! Don't normalize as we don't want to screw up our coefficients
         CASE ('power_series','power_series_edge0','power_series_0_boundaries')
            DO ik = LBOUND(x,DIM=1), UBOUND(x,DIM=1)
               profile_norm = profile_norm + x(ik)/(ik+1)
            END DO
         CASE ('spline','akima_spline','akima_spline_ip')
            DO ik = LBOUND(x,DIM=1), UBOUND(x,DIM=1)
               profile_norm = profile_norm + x(ik)
            END DO
         CASE DEFAULT
            PRINT *,"Error! Unknown profile type in subroutine profile_norm:",prof_type
            PRINT *,x(1:10)
            STOP
      END SELECT
      IF (ltriangulate) profile_norm = 0.0_rprec ! Don't use normalization in triangulation mode
      RETURN
      END FUNCTION profile_norm

      INTEGER FUNCTION count_vars(lvar,xvar,profile_type)
      LOGICAL, INTENT(inout) :: lvar(:)
      REAL(rprec), INTENT(inout) :: xvar(:)
      CHARACTER(LEN=*), INTENT(inout) :: profile_type
      REAL(rprec) :: norm
      count_vars = 0; norm = 0
      count_vars = COUNT(lvar)
      norm = profile_norm(xvar,profile_type)
      IF (norm /= 0) count_vars = count_vars + 1
      RETURN
      END FUNCTION count_vars
      

      SUBROUTINE scale_profile(profile_type, profile_coefficients, factor)
        ! This subroutine takes a profile and scales it by 'factor'.
        ! For some profiles types like 'power_series', this involves multiplying all the profile coefficients
        ! by factor. For other profile types like 'two_power', only certain profile coefficients get scaled.
        IMPLICIT NONE
        CHARACTER(LEN=*) :: profile_type
        REAL(rprec), DIMENSION(:) :: profile_coefficients
        REAL(rprec), INTENT(IN) :: factor

        CALL tolower(profile_type)
        SELECT CASE (profile_type)
        CASE ('power_series','spline','akima_spline','akima_spline_ip','power_series_0_boundaries')
           profile_coefficients = profile_coefficients * factor
        CASE ('two_power','two_power_hollow')
           profile_coefficients(1) = profile_coefficients(1) * factor
        CASE DEFAULT
           PRINT *,"Error! Unknown profile type in subroutine scale_profile:",profile_type
           STOP
        END SELECT
        
      END SUBROUTINE scale_profile



      SUBROUTINE get_equil_beamj(s_val,val,ier)
      IMPLICIT NONE
      REAL(rprec), INTENT(inout) ::  s_val
      REAL(rprec), INTENT(inout)   ::  val
      INTEGER, INTENT(inout)     ::  ier
      INTEGER :: dex, i
      REAL(rprec) :: xp
      REAL(rprec), DIMENSION(10), PARAMETER :: glx = (/                       &
     &   0.01304673574141414, 0.06746831665550774, 0.1602952158504878,         &
     &   0.2833023029353764, 0.4255628305091844, 0.5744371694908156,           &
     &   0.7166976970646236, 0.8397047841495122, 0.9325316833444923,           &
     &   0.9869532642585859 /)
      REAL(rprec), DIMENSION(10), PARAMETER :: glw = (/                       &
     &   0.03333567215434407, 0.0747256745752903, 0.1095431812579910,          &
     &   0.1346333596549982, 0.1477621123573764, 0.1477621123573764,           &
     &   0.1346333596549982, 0.1095431812579910, 0.0747256745752903,           &
     &   0.03333567215434407 /)
      IF (ier < 0) RETURN
      SELECT CASE(beamj_type)
         CASE('spline','akima_spline')
            dex = MINLOC(beamj_aux_s(2:),DIM=1)
            CALL eval_prof_spline(dex,beamj_aux_s(1:dex),beamj_aux_f(1:dex),s_val,val,ier)
         CASE DEFAULT
            CALL eval_prof_stel(s_val,beamj_type,val,21,beamj_aux_f(1:21),ier)
!         CASE ('power_series')
!            val = 0
!            DO i = UBOUND(beamj_aux_f,DIM=1), LBOUND(beamj_aux_f,DIM=1), -1
!               val = s_val*val + beamj_aux_f(i)
!            END DO
!         CASE ('two_power')
!            val = 0
!            val = beamj_aux_f(1) * (1.0 - s_val**beamj_aux_f(2))**beamj_aux_f(3)
!         CASE ('gauss_trunc')
!            val = 0
!            DO i = 1, 10
!               xp = s_val*glx(i)
!               val = val + glw(i) * beamj_aux_f(1) * ( EXP(-(xp/beamj_aux_f(2))**2) &
!                                                     -EXP(-(1/beamj_aux_f(2))**2))
!            END DO
!            val = val * s_val
!         CASE('sum_atan')
!            val = 0
!            IF (s_val >= 1) THEN
!               val = beamj_aux_f(1)+beamj_aux_f(2)+beamj_aux_f(6)+beamj_aux_f(10)+beamj_aux_f(14)+beamj_aux_f(18)
!            ELSE
!               val = beamj_aux_f(1) + &
!                     (4/pi2) * ( beamj_aux_f(2) * ATAN(beamj_aux_f(3)*s_val**beamj_aux_f(4)/(1-s_val)**beamj_aux_f(5)) &
!                                +beamj_aux_f(6) * ATAN(beamj_aux_f(7)*s_val**beamj_aux_f(8)/(1-s_val)**beamj_aux_f(9)) &
!                                +beamj_aux_f(10) * ATAN(beamj_aux_f(11)*s_val**beamj_aux_f(12)/(1-s_val)**beamj_aux_f(13)) &
!                                +beamj_aux_f(14) * ATAN(beamj_aux_f(15)*s_val**beamj_aux_f(16)/(1-s_val)**beamj_aux_f(17)) &
!                                +beamj_aux_f(18) * ATAN(beamj_aux_f(19)*s_val**beamj_aux_f(20)/(1-s_val)**beamj_aux_f(21)))
!            END IF
      END SELECT
      RETURN
      END SUBROUTINE get_equil_beamj
      
      SUBROUTINE get_equil_bootj(s_val,val,ier)
      IMPLICIT NONE
      REAL(rprec), INTENT(inout) ::  s_val
      REAL(rprec), INTENT(inout)   ::  val
      INTEGER, INTENT(inout)     ::  ier
      INTEGER :: dex, i
      REAL(rprec) :: x0,x1, x2, h, x3, xp
      REAL(rprec), DIMENSION(10), PARAMETER :: glx = (/                       &
     &   0.01304673574141414, 0.06746831665550774, 0.1602952158504878,         &
     &   0.2833023029353764, 0.4255628305091844, 0.5744371694908156,           &
     &   0.7166976970646236, 0.8397047841495122, 0.9325316833444923,           &
     &   0.9869532642585859 /)
      REAL(rprec), DIMENSION(10), PARAMETER :: glw = (/                       &
     &   0.03333567215434407, 0.0747256745752903, 0.1095431812579910,          &
     &   0.1346333596549982, 0.1477621123573764, 0.1477621123573764,           &
     &   0.1346333596549982, 0.1095431812579910, 0.0747256745752903,           &
     &   0.03333567215434407 /)
      IF (ier < 0) RETURN
      SELECT CASE(bootj_type)
         CASE('spline','akima_spline')
            dex = MINLOC(bootj_aux_s(2:),DIM=1)
            CALL eval_prof_spline(dex,bootj_aux_s(1:dex),bootj_aux_f(1:dex),s_val,val,ier)
         CASE DEFAULT
            CALL eval_prof_stel(s_val,bootj_type,val,21,bootj_aux_f(1:21),ier)
!         CASE ('power_series')
!            val = 0
!            DO i = UBOUND(bootj_aux_f,DIM=1), LBOUND(bootj_aux_f,DIM=1), -1
!               val = s_val*val + bootj_aux_f(i)
!            END DO
!         CASE ('bump')
!            ! bootj_aux_f(1) : x0 center of bump
!            ! bootj_aux_f(2) : height of bump
!            ! bootj_aux_f(3) : height of bulk
!            x0  = bootj_aux_f(1)
!            x1  = 1.0_rprec-2*(1.0_rprec-x0)
!            x2  = 1.0_rprec
!            x3  = bootj_aux_f(3)
!            h   = bootj_aux_f(2)/((x0-x1)*(x0-x2))
!            val = 0
!            if ((s_val > x1) .and. (s_val < 1.0_rprec)) val = h*(s_val-x1)*(s_val-x2)
!            val = val + x3*s_val*(s_val-1)/(-0.25_rprec)
!         CASE ('two_power')
!            val = 0
!            val = bootj_aux_f(1) * (1.0_rprec - s_val**bootj_aux_f(2))**bootj_aux_f(3)
!         CASE ('gauss_trunc')
!            val = 0
!            DO i = 1, 10
!               xp = s_val*glx(i)
!               val = val + glw(i) * bootj_aux_f(1) * ( EXP(-(xp/bootj_aux_f(2))**2) &
!                                                     -EXP(-(1/bootj_aux_f(2))**2))
!            END DO
!            val = val * s_val
!         CASE('sum_atan')
!            val = 0
!            IF (s_val >= 1) THEN
!               val = bootj_aux_f(1)+bootj_aux_f(2)+bootj_aux_f(6)+bootj_aux_f(10)+bootj_aux_f(14)+bootj_aux_f(18)
!            ELSE
!               val = bootj_aux_f(1) + &
!                     (4/pi2) * ( bootj_aux_f(2) * ATAN(bootj_aux_f(3)*s_val**bootj_aux_f(4)/(1-s_val)**bootj_aux_f(5)) &
!                                +bootj_aux_f(6) * ATAN(bootj_aux_f(7)*s_val**bootj_aux_f(8)/(1-s_val)**bootj_aux_f(9)) &
!                                +bootj_aux_f(10) * ATAN(bootj_aux_f(11)*s_val**bootj_aux_f(12)/(1-s_val)**bootj_aux_f(13)) &
!                                +bootj_aux_f(14) * ATAN(bootj_aux_f(15)*s_val**bootj_aux_f(16)/(1-s_val)**bootj_aux_f(17)) &
!                                +bootj_aux_f(18) * ATAN(bootj_aux_f(19)*s_val**bootj_aux_f(20)/(1-s_val)**bootj_aux_f(21)))
!            END IF
      END SELECT
      RETURN
      END SUBROUTINE get_equil_bootj

      SUBROUTINE eval_profile(s_val, profile_type, val, profile_aux_s, profile_aux_f, ier)
        ! s_val = Value of normalized flux s at which you wish to evaluate the profile.
        ! On return, 'val' stores the value of the profile function at the selected s.
        IMPLICIT NONE
        REAL(rprec) ::  s_val, val
        INTEGER ::  ier, dex
        CHARACTER(LEN=*) :: profile_type
        REAL(rprec), DIMENSION(21) :: profile_aux_s, profile_aux_f
        IF (ier < 0) RETURN
        CALL tolower(profile_type)
        SELECT CASE(TRIM(profile_type))
        CASE('spline','akima_spline')
           dex = MINLOC(profile_aux_s(2:),DIM=1)
           CALL eval_prof_spline(dex, profile_aux_s(1:dex), profile_aux_f(1:dex), s_val, val, ier)
        CASE DEFAULT
           CALL eval_prof_stel(s_val, profile_type, val, 21, profile_aux_f(1:21), ier)
        END SELECT
        RETURN
      END SUBROUTINE eval_profile

!     J_STAR by DA Spong
!     computes the trapped branch of jstar on a single flux surface
!     and for a single value of ep/mu.  jstar is only non-zero for
!     values of theta where trapped particle orbits can exist.  at
!     theta values which are in the passing particle regime (ep/mu > bmax)
!     or the forbidden regime (ep/mu < bmin) jstar is set to 0.  a jstar
!     topology is assumed here such that as theta runs from 0 to pi,
!     one first encounters the ep/mu = bmax point and then
!     the ep/mu = bmin point
      SUBROUTINE j_star(modb, bmin, bmax, ep_mu, jstar, nzeta, ntheta)
      IMPLICIT NONE
      REAL(rprec), DIMENSION(nzeta,ntheta), INTENT(in) :: modb
      REAL(rprec), DIMENSION(ntheta), INTENT(in) :: bmin, bmax
      REAL(rprec), INTENT(in) :: ep_mu
      REAL(rprec), DIMENSION(ntheta), INTENT(out) :: jstar
      INTEGER nzeta, ntheta
      REAL(rprec), PARAMETER :: zero = 0, one = 1
      INTEGER :: ku_lw, ku_up, ku, istat
      REAL(rprec) :: dzeta
      REAL(rprec), DIMENSION(ntheta) :: test_bmin, test_bmax
      REAL(rprec) , DIMENSION(ntheta) :: test1, test2
      REAL(rprec) , ALLOCATABLE,  DIMENSION(:,:) :: vpar1
      LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: l3v
      ku_lw = 1
      ku_up = ntheta
      jstar = zero
      dzeta = one
      IF (nzeta .gt. 1) dzeta = one/REAL((nzeta-1),rprec)
      test_bmax = ep_mu - bmax(:ntheta)
      test_bmin = ep_mu - bmin(:ntheta)
      test1(2:ntheta) = test_bmax(2:ntheta)*test_bmax(:ntheta-1)
      test2(2:ntheta) = test_bmin(2:ntheta)*test_bmin(:ntheta-1)
      DO ku = 2,ntheta
          IF(test1(ku) .le. zero) ku_lw = ku
          IF(test2(ku) .le. zero) ku_up = ku
      END DO
      IF (ku_lw .ge. ku_up) RETURN
      ALLOCATE (vpar1(nzeta,ku_up-ku_lw+1), l3v(nzeta,ku_up-ku_lw+1))
      vpar1 = one - modb(:,ku_lw:ku_up)/ep_mu
      l3v = vpar1 .gt. zero
      WHERE (l3v) vpar1 = SQRT(vpar1)/modb(:,ku_lw:ku_up)
      jstar(ku_lw:ku_up) = dzeta * SUM(vpar1, mask=l3v, dim=1)
      DEALLOCATE (vpar1, l3v)
      END SUBROUTINE j_star
      
      
      SUBROUTINE cross_product(A,B,AxB)
      IMPLICIT none 
      REAL(rprec),INTENT(in) :: A(3),B(3)
      REAL(rprec),INTENT(out) :: AxB(3)

      AxB = (/ (A(2)*B(3)-A(3)*B(2)), & 
           (A(3)*B(1)-A(1)*B(3)), & 
          (A(1)*B(2)-A(2)*B(1)) /)
 
      END SUBROUTINE cross_product

      SUBROUTINE smoothg(A,N,DLZ,B)
      ! Smoothing of A with gaussian width DLZ (H. Mynick 12/06/09)
      IMPLICIT NONE
      REAL(rprec), INTENT(IN) :: A(N)
      REAL(rprec), INTENT(IN) :: DLZ
      INTEGER, INTENT(IN) :: N
      REAL(rprec), INTENT(OUT) :: B(N)
      INTEGER :: i, j
      REAL(rprec) :: arg,gaussn, dlth, gsum
      dlth=pi2/(N-1)  !=th incrmt btwn successive meshpts.
      B=0.0_rprec 
      DO i=1,N
         gsum=0.0_rprec
         DO j=1,N
            arg=(i-j)*dlth/DLZ
            gaussn=EXP(-arg*arg)
            B(i)=B(i)+gaussn*A(j)
            gsum=gsum+gaussn
         END DO                  !j
         B(i)=B(i)/gsum
      END DO                     !i
      RETURN
      END SUBROUTINE smoothg
      
      FUNCTION polyfit(vx, vy, d)
      ! Fit a polynomial of degree d to data (vx,vy) using one of the LAPACK
      ! subroutines for least-squares fitting, DGELS.
      implicit none
      integer, intent(in)                   :: d
      integer, parameter                    :: dp = selected_real_kind(15, 307)
      real(dp), dimension(d+1)              :: polyfit
      real(dp), dimension(:), intent(in)    :: vx, vy
      
      real(dp), dimension(:,:), allocatable :: matrix
      integer :: n, m, j, lwork, info, ldb
      real(dp), dimension(:), allocatable :: work, b

      n = d+1
      m = SIZE(vx)
      IF (SIZE(vy) .ne. m) STOP "Error in polyfit: sizes of vx and vy do not match"

      ! Form the Vandermonde matrix:
      ALLOCATE(matrix(m,n))
      matrix(:,1) = 1
      DO j = 1,d
         matrix(:,j+1) = matrix(:,j) * vx
      END DO

      ldb = MAX(m,n)
      ALLOCATE(b(ldb))
      b(1:m) = vy
      ! The DGELS argument b stores the data vector b on input, and it gets over-written with the solution vector.
      ! Note that the sizes of the input and solution vector can be different.

      ! Get size of work array:
      lwork = -1
      ALLOCATE(work(1))
      CALL DGELS('N',m,n,1,matrix,m,b,ldb,work,lwork,info)
      IF (info .ne. 0) THEN
         PRINT *,"Error 1 in polyfit LAPACK DGELS: info=",info
         STOP
      END IF
      lwork = work(1)
      DEALLOCATE(work)
      ALLOCATE(work(lwork))

      ! Main LAPACK call:
      CALL DGELS('N',m,n,1,matrix,m,b,ldb,work,lwork,info)
      IF (info .ne. 0) THEN
         PRINT *,"Error 2 in polyfit LAPACK DGELS: info=",info
         STOP
      END IF

      polyfit = b(1:n)
      DEALLOCATE(matrix,work,b)

      END FUNCTION polyfit


      FUNCTION polyfit_0_boundaries(vx, vy, d)
      ! Fit a polynomial of the form p(x) * x * (1-x) to data (vx,vy), where p(s) is a polynomial
      ! of degree d, using one of the LAPACK subroutines for least-squares fitting, DGELS.
      implicit none
      integer, intent(in)                   :: d
      integer, parameter                    :: dp = selected_real_kind(15, 307)
      real(dp), dimension(d+1)              :: polyfit_0_boundaries
      real(dp), dimension(:), intent(in)    :: vx, vy
      
      real(dp), dimension(:,:), allocatable :: matrix
      integer :: n, m, j, lwork, info, ldb
      real(dp), dimension(:), allocatable :: work, b

      n = d + 1
      m = SIZE(vx)
      IF (SIZE(vy) .ne. m) STOP "Error in polyfit: sizes of vx and vy do not match"

      ! Form the Vandermonde matrix * vx * (1-vx):
      ALLOCATE(matrix(m,n))
      matrix(:,1) = vx * (1 - vx)
      DO j = 1,d
         matrix(:,j+1) = matrix(:,j) * vx
      END DO

      ldb = MAX(m,n)
      ALLOCATE(b(ldb))
      b(1:m) = vy
      ! The DGELS argument b stores the data vector b on input, and it gets over-written with the solution vector.
      ! Note that the sizes of the input and solution vector can be different.

      ! Get size of work array:
      lwork = -1
      ALLOCATE(work(1))
      CALL DGELS('N',m,n,1,matrix,m,b,ldb,work,lwork,info)
      IF (info .ne. 0) THEN
         PRINT *,"Error 1 in polyfit LAPACK DGELS: info=",info
         STOP
      END IF
      lwork = work(1)
      DEALLOCATE(work)
      ALLOCATE(work(lwork))

      ! Main LAPACK call:
      CALL DGELS('N',m,n,1,matrix,m,b,ldb,work,lwork,info)
      IF (info .ne. 0) THEN
         PRINT *,"Error 2 in polyfit LAPACK DGELS: info=",info
         STOP
      END IF

      polyfit_0_boundaries = b(1:n)
      DEALLOCATE(matrix,work,b)

      END FUNCTION polyfit_0_boundaries


      function polyval(c,s,d)
      integer, intent(in)                   :: d
      integer, parameter                    :: dp = selected_real_kind(15, 307)
      real(dp), dimension(d), intent(in)  :: c
      real(dp), intent(in)    :: s
      real(dp)                              :: polyval
      INTEGER :: j
      polyval = 0
      DO j = d, 1, -1
         polyval = polyval*s + c(j)
      END DO
      RETURN
      end function polyval

      SUBROUTINE move_txtfile(file_in, file_out)
      USE safe_open_mod
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(in) :: file_in
      CHARACTER(LEN=*), INTENT(in) :: file_out
      LOGICAL ::  lfile_check
      INTEGER ::  iunit_in, iunit_out, iflag_in, iflag_out
      CHARACTER(LEN=4096) :: temp_str
      iflag_in  = 0
      iflag_out = 0
      iunit_in = 13
      iunit_out = 14
      INQUIRE(FILE=TRIM(file_in),EXIST=lfile_check)
      IF (.not.lfile_check) RETURN
      CALL safe_open(iunit_in,iflag_in,TRIM(file_in),'OLD','formatted') 
      CALL safe_open(iunit_out,iflag_out,TRIM(file_out),'REPLACE','formatted')
      DO
         temp_str = ''
         READ(UNIT=iunit_in,FMT='(A)',IOSTAT=iflag_in) temp_str
         IF (iflag_in .lt. 0) EXIT
         WRITE(UNIT=iunit_out,FMT='(A)',IOSTAT=iflag_out) TRIM(temp_str)
         IF (iflag_in .gt. 0) PRINT *,'ipecopt_move_txtfile error IOSTAT=',iflag_in
      END DO
      CLOSE(iunit_in,STATUS='delete')
      CLOSE(iunit_out)
      RETURN
      END SUBROUTINE move_txtfile

      SUBROUTINE copy_txtfile(file_in, file_out)
      USE safe_open_mod
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(in) :: file_in
      CHARACTER(LEN=*), INTENT(in) :: file_out
      INTEGER ::  iunit_in, iunit_out, iflag_in, iflag_out
      LOGICAL ::  lfile_check
      CHARACTER(LEN=4096) :: temp_str
      iflag_in  = 0
      iflag_out = 0
      iunit_in = 13
      iunit_out = 14
      INQUIRE(FILE=TRIM(file_in),EXIST=lfile_check)
      IF (.not.lfile_check) RETURN
      CALL safe_open(iunit_in,iflag_in,TRIM(file_in),'OLD','formatted') 
      CALL safe_open(iunit_out,iflag_out,TRIM(file_out),'REPLACE','formatted')
      DO
         temp_str = ''
         READ(UNIT=iunit_in,FMT='(A)',IOSTAT=iflag_in) temp_str
         IF (iflag_in .lt. 0) EXIT
         WRITE(UNIT=iunit_out,FMT='(A)',IOSTAT=iflag_out) TRIM(temp_str)
         IF (iflag_in .gt. 0) PRINT *,'ipecopt_copy_txtfile error IOSTAT=',iflag_in
      END DO
      CLOSE(iunit_in)
      CLOSE(iunit_out)
      RETURN
      END SUBROUTINE copy_txtfile

      FUNCTION TEMfunc_proll(x) ! TEM function for DCUHRE
      IMPLICIT NONE
      DOUBLE PRECISION :: x !lambda
      DOUBLE PRECISION :: TEMfunc_proll
      DOUBLE PRECISION :: result,abserr,A,B
      INTEGER :: neval,ier, last
      INTEGER, PARAMETER :: KEY = 3, NW=1024
      DOUBLE PRECISION, PARAMETER :: EPSABS=1.0D-9, EPSREL=1.0D-3
      INTEGER, DIMENSION(nw) :: work
      DOUBLE PRECISION, DIMENSION(nw) :: iwork

      lam_TEM = x
      A = 0; B= pi2
      CALL DQAG(TEMsubfunc_proll,A,B,EPSABS,EPSREL,KEY,result,abserr,neval,ier,&
                 128,NW,last,iwork,work)
      TEMfunc_proll=result
      RETURN

      END FUNCTION TEMfunc_proll

      FUNCTION  TEMsubfunc_proll(x) ! TEM function for DCUHRE
      IMPLICIT NONE
      DOUBLE PRECISION :: x !z
      DOUBLE PRECISION :: TEMsubfunc_proll
      INTEGER :: ier
      DOUBLE PRECISION :: B0, L20
      DOUBLE PRECISION, PARAMETER ::  zero = 0.0E+00
      DOUBLE PRECISION, PARAMETER ::  one = 1.0E+00
      DOUBLE PRECISION, PARAMETER ::  half = 0.5D+00

      CALL EZspline_interp(Bhat_spl,x,B0,ier)
      CALL EZspline_interp(L2_spl,x,L20,ier)
      !PRINT *,'B0,L2',B0,L20,x
      IF ((delta_TEM/lam_TEM - B0) < 0 ) THEN
         TEMsubfunc_proll = zero
      ELSE
         TEMsubfunc_proll = L20*(one-half*lam_TEM*B0)/DSQRT(one-lam_TEM*B0)
      END IF
      !WRITE(327,*) x,B0,L20,lam_TEM,TEMsubfunc_proll
      RETURN

      END FUNCTION TEMsubfunc_proll



      FUNCTION TEMfunc_proll_omegatau(x) ! TEM function for DCUHRE
      IMPLICIT NONE
      DOUBLE PRECISION :: x !lambda
      DOUBLE PRECISION :: TEMfunc_proll_omegatau
      DOUBLE PRECISION :: result_w,result_t,abserr,A,B
      INTEGER :: neval,ier, last
      INTEGER, PARAMETER :: KEY = 3, NW=1024
      DOUBLE PRECISION, PARAMETER :: EPSABS=1.0D-9, EPSREL=1.0D-3
      INTEGER, DIMENSION(nw) :: work
      DOUBLE PRECISION, DIMENSION(nw) :: iwork

      lam_TEM = x
      A = 0; B= pi2
      CALL DQAG(TEMsubfunc_proll,A,B,EPSABS,EPSREL,KEY,result_w,abserr,neval,ier,&
                 128,NW,last,iwork,work)
      IF (result_w == 0) THEN
         TEMfunc_proll_omegatau = 0
         RETURN
      END IF
      lam_TEM = x
      A = 0; B= pi2
      CALL DQAG(TEMsubfunc_proll_tau,A,B,EPSABS,EPSREL,KEY,result_t,abserr,neval,ier,&
                 128,NW,last,iwork,work)
      !WRITE(327,*) x,result_w,result_t
      TEMfunc_proll_omegatau=result_w/result_t
      RETURN

      END FUNCTION TEMfunc_proll_omegatau

      FUNCTION  TEMsubfunc_proll_tau(x) ! TEM function for DCUHRE
      IMPLICIT NONE
      DOUBLE PRECISION :: x !z
      DOUBLE PRECISION :: TEMsubfunc_proll_tau
      INTEGER :: ier
      DOUBLE PRECISION :: B0, L20
      DOUBLE PRECISION, PARAMETER ::  zero = 0.0E+00
      DOUBLE PRECISION, PARAMETER ::  one = 1.0E+00
      DOUBLE PRECISION, PARAMETER ::  half = 0.5D+00

      CALL EZspline_interp(Bhat_spl,x,B0,ier)
      CALL EZspline_interp(L2_spl,x,L20,ier)
      !PRINT *,'B0,L2',B0,L20,x
      IF ((delta_TEM/lam_TEM - B0) < 0 ) THEN
         TEMsubfunc_proll_tau = zero
      ELSE
         TEMsubfunc_proll_tau = REAL(1)/DSQRT(one-lam_TEM*B0)
      END IF
      !WRITE(327,*) x,B0,L20,lam_TEM,TEMsubfunc_proll
      RETURN

      END FUNCTION TEMsubfunc_proll_tau

      SUBROUTINE copy_boozer_file(file_in,file_out)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(in) :: file_in
      CHARACTER(LEN=*), INTENT(in) :: file_out
      INTEGER :: ier
      LOGICAL :: lfile_check

      ier = 0
!DEC$ IF DEFINED (NETCDF)
      INQUIRE(FILE='boozmn_'//TRIM(file_in)//'.nc',EXIST=lfile_check)
!DEC$ ELSE
      INQUIRE(FILE='boozmn.'//TRIM(file_in),EXIST=lfile_check)
!DEC$ ENDIF
      IF (.not.lfile_check) RETURN
      CALL read_boozer_file(file_in,ier)
      IF (ier ==0)  THEN
!DEC$ IF DEFINED (NETCDF)
         CALL write_boozer_nc(file_out,ier)
!DEC$ ELSE
         CALL write_boozer_bin(file_out,ier)
!DEC$ ENDIF
         !CALL read_boozer_deallocate
      END IF

      RETURN
      END SUBROUTINE copy_boozer_file

      SUBROUTINE fit_profile(ptype,ntarg,sarr,farr,ncoefs,coefs)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(in)   :: ptype
      INTEGER, INTENT(in) :: ntarg,ncoefs
      REAL(rprec), INTENT(in) :: sarr(ntarg), farr(ntarg)
      REAL(rprec), INTENT(inout) :: coefs(ncoefs)
      INTEGER :: nc, ik, maxfev_local, nfev, info, njev, maxfev, nprint, mode
      INTEGER, DIMENSION(ncoefs) :: ipvt
      DOUBLE PRECISION :: ftol,xtol,gtol,factor
      DOUBLE PRECISION, ALLOCATABLE :: xc_opt(:), diag(:), qtf(:), wa1(:), wa2(:), wa3(:)
      DOUBLE PRECISION, ALLOCATABLE :: fval(:),wa4(:)
      DOUBLE PRECISION, ALLOCATABLE :: fjac(:,:)

      ! Adjust number of coefficients per fit type
      nc = ncoefs
      CALL tolower(ptype)
      SELECT CASE (ptype)
         CASE('two_power')
            nc = 3
         CASE('two_power_hollow')
            nc = 3
         CASE('two_lorentz')
            nc = 8
         CASE('gauss_trunc')
            nc = 2
         CASE('power_series','power_series_0i0','power_series_edge0','power_series_i','power_series_i_edge0','power_series_0_boundaries')
            DO ik = 1, ncoefs
               IF (coefs(ik) /=0 ) nc = ik
            END DO
         CASE('pedestal','sum_atan')
            nc = 21
         CASE('bump')
            nc = 3
         CASE DEFAULT
            PRINT *,"Error! Unknown profile type in subroutine fit_profile:",ptype
            STOP
      END SELECT
      IF (ptype == 'power_series') THEN
         ! For fitting polynomials, 'polyfit' is more robust than LMDER, so use polyfit.
         coefs(1:nc) = polyfit(sarr,farr,nc-1)
         RETURN
      ELSEIF (ptype == 'power_series_0_boundaries') THEN
         coefs(1:nc) = polyfit_0_boundaries(sarr,farr,nc-1)
         RETURN
      END IF

      nfit_targs = ntarg+1
      ALLOCATE(fit_targs(nfit_targs,2))
      fit_targs(:,1) = sarr
      fit_targs(:,2) = farr
      SUM_target = SUM(farr)
      fit_type = ptype
      ! ALLOCATE Vars
      ALLOCATE(xc_opt(nc),diag(nc),qtf(nc),wa1(nc),wa2(nc),wa3(nc),fjac(nfit_targs,nc),&
               fval(nfit_targs),wa4(nfit_targs))
      ! Setup LMDER
      xc_opt(1:nc) = coefs(1:nc)
      fval = 1.0E-30
      fjac = 0
      ftol = 1.0E-6
      xtol = 1.0E-6
      gtol = 1.0E-30
      maxfev_local = 2000
      diag(:) = 1
      mode = 1
      factor = 0.1
      nprint = 0
      info   = 9
      nfev   = 0
      njev   = 0
      ipvt   = 0
      qtf    = 0.0
      wa1 = 0; wa2 = 0; wa3 = 0; wa4 = 0
!!$      print *,"About to call lmder_serial."
!!$      print *,"  nfit_targs:",nfit_targs
!!$      print *,"  nc:",nc
!!$      print *,"  xc_opt:",xc_opt
!!$      print *,"  fval:",fval
!!$      print *,"  fjac:",fjac
!!$      print *,"  nfit_targs:",nfit_targs
!!$      print *,"  ftol:",ftol
!!$      print *,"  xtol:",xtol
!!$      print *,"  gtol:",gtol
!!$      print *,"  maxfev_local:",maxfev_local
!!$      print *,"  diag:",diag
!!$      print *,"  mode:",mode
!!$      print *,"  factor:",factor
!!$      print *,"  nprint:",nprint
!!$      print *,"  info:",info
!!$      print *,"  nfev:",nfev
!!$      print *,"  njev:",njev
!!$      print *,"  ipvt:",ipvt
!!$      print *,"  qtf:",qtf
!!$      print *,"  wa1:",wa1
!!$      print *,"  wa2:",wa2
!!$      print *,"  wa3:",wa3
!!$      print *,"  wa4:",wa4
      CALL lmder_serial(fit_prof_fcn,nfit_targs,nc,xc_opt,fval,fjac,nfit_targs,ftol,xtol,gtol,&
                    maxfev_local,diag,mode,factor,nprint,info,nfev,njev,ipvt,qtf,&
                    wa1,wa2,wa3,wa4)
!!$      print *,"After call to lmder_serial:"
!!$      print *,"  nfit_targs:",nfit_targs
!!$      print *,"  nc:",nc
!!$      print *,"  xc_opt:",xc_opt
!!$      print *,"  fval:",fval
!!$      print *,"  fjac:",fjac
!!$      print *,"  nfit_targs:",nfit_targs
!!$      print *,"  ftol:",ftol
!!$      print *,"  xtol:",xtol
!!$      print *,"  gtol:",gtol
!!$      print *,"  maxfev_local:",maxfev_local
!!$      print *,"  diag:",diag
!!$      print *,"  mode:",mode
!!$      print *,"  factor:",factor
!!$      print *,"  nprint:",nprint
!!$      print *,"  info:",info
!!$      print *,"  nfev:",nfev
!!$      print *,"  njev:",njev
!!$      print *,"  ipvt:",ipvt
!!$      print *,"  qtf:",qtf
!!$      print *,"  wa1:",wa1
!!$      print *,"  wa2:",wa2
!!$      print *,"  wa3:",wa3
!!$      print *,"  wa4:",wa4
      coefs(1:nc) = xc_opt(1:nc)
      DEALLOCATE(fit_targs,xc_opt,diag,qtf,wa1,wa2,wa3,fjac,fval,wa4)
      RETURN
      END SUBROUTINE fit_profile
      
      SUBROUTINE  fit_prof_fcn(m,n,x,fvec,fjac,ldfjac,iflag)
      IMPLICIT NONE
      INTEGER m,n,ldfjac,iflag, ier
      DOUBLE PRECISION x(n),fvec(m),fjac(ldfjac,n), x_temp(n)
      INTEGER :: j,k
      DOUBLE PRECISION :: val, val2, dx, vals, vals2
      DOUBLE PRECISION :: x_temp2(n)
      DOUBLE PRECISION, PARAMETER :: delta = 1.000001
      IF (iflag == 1) THEN
         vals =0
         DO k = 1, m-1
            CALL eval_prof_stel(fit_targs(k,1),fit_type,val,n,x,ier)
            vals = vals + val
            fvec(k) = val - fit_targs(k,2)
         END DO
         fvec(m) = vals - SUM_target
      ELSE IF (iflag == 2) THEN
         x_temp  = x
         x_temp2 = x
         DO j = 1, n
            vals = 0
            vals2 = 0
            x_temp2(j) = x_temp(j)*delta
            dx = x_temp2(j)-x_temp(j)
            IF (dx == 0) THEN
               dx = delta-1
               x_temp2(j) = dx
            END IF
            DO k = 1, m-1
               CALL eval_prof_stel(fit_targs(k,1),fit_type,val,n,x_temp,ier)
               CALL eval_prof_stel(fit_targs(k,1),fit_type,val2,n,x_temp2,ier)
               fjac(k,j) = (val2-val)/dx
               vals = vals + val
               vals2 = vals2 + val2
            END DO
            fjac(m,j) = (vals2-vals)/dx
            x_temp2(j) = x_temp(j)
         END DO
      END IF
      RETURN
      END SUBROUTINE fit_prof_fcn
      
      END MODULE equil_utils
