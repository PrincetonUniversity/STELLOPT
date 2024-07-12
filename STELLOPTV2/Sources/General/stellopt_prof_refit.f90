!-----------------------------------------------------------------------
!     Subroutine:    stellopt_prof_refit
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/01/2012
!     Description:   This subroutine refits the existing profiles to
!                    the data in order to speed reconstruction.
!                    Note we use polyfit from LAPACK.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_prof_refit(n,x,lscreen)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE stellopt_vars
      USE equil_utils
      USE parambs, ONLY: ajBbs, bsnorm, rhoar, l_boot_all, aibs
      USE EZspline
      
      IMPLICIT NONE
      
!-----------------------------------------------------------------------
!     Input Variables
!        n       Number of function variables
!        x       Vector of function variables
!     lscreen    Controls printing to the screen
!----------------------------------------------------------------------
      LOGICAL, INTENT(in)      :: lscreen
      INTEGER, INTENT(in)      ::  n
      REAL(rprec), INTENT(inout)  ::  x(n)
      
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      INTEGER :: n_ne, n_te, n_ti, n_ne_line, n_boot, i, j, dex, ier
      INTEGER, ALLOCATABLE :: idx(:)
      REAL(rprec) :: temp, norm_ne, norm_te, norm_ti, avg_jdotb, alpha,&
                     j_curv, s_temp, ne_err, new_area
      REAL(rprec) :: x0(3), x1(3)
      REAL(rprec), ALLOCATABLE :: s_val(:), f_val(:)
      LOGICAL, PARAMETER :: LBOOT_REFIT =.false.
      INTEGER, PARAMETER :: NPRESS = 512
      INTEGER, PARAMETER :: A_COEFF = 5
      REAL(rprec) :: polynom(A_COEFF + 1)
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      n_ne = COUNT(sigma_ne < bigno_ne)
      n_te = COUNT(sigma_te < bigno)
      n_ti = COUNT(sigma_ti < bigno)
      n_ne_line = COUNT(sigma_ne_line < bigno_ne)
      n_boot = COUNT(sigma_bootstrap < bigno)
      !----------------------------------------------------------------
      !     BOOTSTRAP CURRENT
      !----------------------------------------------------------------
      IF (n_boot > 0 .and. LBOOT_REFIT) THEN
         dex = n_boot+1
         ier=0
         CALL stellopt_toboozer(lscreen,ier)
         IF (ier == 0) THEN
            CALL stellopt_bootsj(lscreen,ier)
            IF (ier == 0) THEN
               ALLOCATE(s_val(dex),f_val(dex))
               s_val(1) = 0.0; f_val(1) = 0.0
               j=2
               DO i = 1, nprof
                  IF (sigma_bootstrap(i) >= bigno) CYCLE
                  s_val(j) = rhoar(i)
                  f_val(j) = ajBbs(i)
                  j=j+1
               END DO
               f_val(1) = f_val(2)
               s_val(dex) = 1.0; f_val(dex) = 0.0
               DO i = 1, ndatafmax
                  IF ((bootj_aux_s(i) > 1.0) .or. (bootj_aux_s(i) < 0.0)) CYCLE
                  ier = 0
                  CALL get_equil_jdotb(bootj_aux_s(i),avg_jdotb,ier)
                  CALL get_equil_jcurv(bootj_aux_s(i),j_curv,ier)
                  CALL eval_prof_spline(dex,s_val(1:dex),f_val(1:dex),bootj_aux_s(i),temp,ier)
                  alpha  = j_curv*temp/avg_jdotb
                  !WRITE(6,'(6F12.5)'),i,bootj_aux_s(i),avg_jdotb,j_curv,temp,alpha
                  bootj_aux_f(i) = bootj_aux_f(i)-refit_param*(bootj_aux_f(i)-alpha)
               END DO
            ELSE
               IF (lscreen) WRITE(6,*) ' BOOTSTRAP FITTING FAILED: BOOTSJ'
            END IF
         ELSE
            IF (lscreen) WRITE(6,*) ' BOOTSTRAP FITTING FAILED: BOOZ_XFORM'
         END IF
      END IF
      IF (ALLOCATED(s_val)) DEALLOCATE(s_val)
      IF (ALLOCATED(f_val)) DEALLOCATE(f_val)
      !----------------------------------------------------------------
      !     ELECTRON DENSITY (NE)
      !----------------------------------------------------------------
      IF (n_ne > 3)  THEN
         !PRINT *,'---------S - Calc--------'
         ALLOCATE(s_val(n_ne),idx(n_ne),f_val(n_ne))
         j = 1
         s_val = 0.0; f_val = 0.0
         DO i = 1, nprof
            IF (sigma_ne(i) >= bigno_ne) CYCLE
            ier = 0
            f_val(j) = target_ne(i)/ne_norm    
            CALL get_equil_s(r_ne(i),phi_ne(i),z_ne(i),s_val(j),ier)
            IF (ier < 0 .or. s_val(j) > 1.0 .or. s_val(j) < 0.0) THEN 
               s_val(j) = 1.5+0.1*j
               f_val(j) = 0.0
            END IF
            j = j + 1
         END DO
         DO i = 1, n_ne-1  ! Replace duplicates
            DO j = i+1, n_ne
               IF (s_val(i) == s_val(j)) THEN
                  s_val(j) = s_val(j) + 1.5+0.1*j
                  f_val(j) = 0.0
               END IF
            END DO
         END DO
         idx = 1
         DO i = 1, n_ne ! Sort
            dex = MINLOC(s_val,DIM=1,MASK=idx > 0)
            temp = s_val(i)
            s_val(i) = s_val(dex)
            s_val(dex) = temp
            temp = f_val(i)
            f_val(i) = f_val(dex)
            f_val(dex) = temp
            idx(i) = 0
         END DO
         dex = MINLOC(s_val, DIM=1, MASK = s_val >= 1)
         IF (dex > 2) THEN
            DO WHILE(f_val(dex-1) == 0.0)
              dex = dex - 1
            END DO
            s_val(dex) = 0.5*(1.0-s_val(dex-1))+s_val(dex-1)
            f_val(dex) = 0.5*(f_val(dex-1))
            s_val(dex) = 1.0
            f_val(dex) = 0.0
            polynom = polyfit(s_val,f_val,A_COEFF)
            polynom(A_COEFF+1) = -SUM(polynom(1:A_COEFF))
            DO i = 2, ndatafmax  ! Adjust all but central spline knot
               IF (ne_aux_s(i) > 1.0 .or. ne_aux_s(i) < 0.0) CYCLE
               ier = 0
               s_temp = ne_aux_s(i)
               !CALL eval_prof_spline(dex,s_val(1:dex),f_val(1:dex),s_temp,temp,ier)
               temp = polyval(polynom,s_temp,A_COEFF+1)
               ne_aux_f(i) = ne_aux_f(i) - refit_param*(ne_aux_f(i) - temp)
               IF (ne_aux_f(i) < 0.0) ne_aux_f(i) = 0.0
            END DO
            ne_aux_f(1) = ne_aux_f(1) - refit_param*(ne_aux_f(1) - f_val(1))
            IF (ne_aux_s(1) < 0.0) ne_aux_f(1) = ne_aux_f(2)
            DO i = 2, ndatafmax
              IF (ne_aux_s(i) >= 1.0 .or. ne_aux_s(i) < 0.0) CYCLE
              IF (ne_aux_f(i) == 0.0) ne_aux_f(i) = ne_aux_f(i-1)
            END DO
         END IF
      END IF
      IF (ALLOCATED(s_val)) DEALLOCATE(s_val)
      IF (ALLOCATED(idx)) DEALLOCATE(idx)
      IF (ALLOCATED(f_val)) DEALLOCATE(f_val)
      !----------------------------------------------------------------
      !     ELECTRON TEMPERATURE (TE)
      !----------------------------------------------------------------
      IF (n_te > 3)  THEN
         ALLOCATE(s_val(n_te),idx(n_te),f_val(n_te))
         j = 1
         DO i = 1, nprof
            IF (sigma_te(i) >= bigno) CYCLE
            ier = 0
            f_val(j) = target_te(i)
            CALL get_equil_s(r_te(i),phi_te(i),z_te(i),s_val(j),ier)
            IF (ier < 0 .or. s_val(j) > 1.0 .or. s_val(j) < 0.0) THEN 
               s_val(j) = 1.5+0.1*j
               f_val(j) = 0.0
            END IF
            j = j + 1
         END DO
         DO i = 1, n_te-1  ! Replace duplicates
            DO j = i+1, n_te
               IF (s_val(i) == s_val(j)) THEN
                  s_val(j) = s_val(j) + 1.5+0.1*j
                  f_val(j) = 0.0
               END IF
            END DO
         END DO
         idx = 1
         DO i = 1, n_te ! Sort
            dex = MINLOC(s_val,DIM=1,MASK=idx > 0)
            temp = s_val(i)
            s_val(i) = s_val(dex)
            s_val(dex) = temp
            temp = f_val(i)
            f_val(i) = f_val(dex)
            f_val(dex) = temp
            idx(i) = 0
         END DO
         dex = MINLOC(s_val, DIM=1, MASK = s_val >= 1)
         IF (dex > 2) THEN
            DO WHILE(f_val(dex-1) == 0.0)
              dex = dex - 1
            END DO
            s_val(dex) = 0.5*(1.0-s_val(dex-1))+s_val(dex-1)
            f_val(dex) = 0.5*(f_val(dex-1))
            s_val(dex+1) = 1.0
            f_val(dex+1) = 0.0
            polynom = polyfit(s_val,f_val,A_COEFF)
            polynom(A_COEFF+1) = -SUM(polynom(1:A_COEFF))
            DO i = 2, ndatafmax  ! Adjust spline knots
               IF (te_aux_s(i) > 1.0 .or. te_aux_s(i) < 0.0) CYCLE
               ier = 0
               s_temp = te_aux_s(i)
               !CALL eval_prof_spline(dex,s_val(1:dex),f_val(1:dex),s_temp,temp,ier)
               temp = polyval(polynom,s_temp,A_COEFF+1)
               te_aux_f(i) = te_aux_f(i) - refit_param*(te_aux_f(i) - temp)
               IF (te_aux_s(i) == 1.0) te_aux_f(i) = 0.0
               IF (te_aux_f(i) < 0.0) te_aux_f(i) = 0.0
            END DO
            te_aux_f(1) = te_aux_f(1) - refit_param*(te_aux_f(1) - f_val(1))
            IF (te_aux_f(1) < 0.0) te_aux_f(1) = te_aux_f(2)
            DO i = 2, ndatafmax
              IF (te_aux_s(i) >= 1.0 .or. te_aux_s(i) < 0.0) CYCLE
              IF (te_aux_f(i) == 0.0) te_aux_f(i) = te_aux_f(i-1)*0.5
            END DO
         END IF
      END IF
      IF (ALLOCATED(s_val)) DEALLOCATE(s_val)
      IF (ALLOCATED(idx)) DEALLOCATE(idx)
      IF (ALLOCATED(f_val)) DEALLOCATE(f_val)
      !----------------------------------------------------------------
      !     ION TEMPERATURE (TI)
      !----------------------------------------------------------------
      IF (n_ti > 3)  THEN
         ALLOCATE(s_val(n_ti),idx(n_ti),f_val(n_ti))
         j = 1
         DO i = 1, nprof
            IF (sigma_ti(i) >= bigno) CYCLE
            ier = 0
            f_val(j) = target_ti(i)
            CALL get_equil_s(r_ti(i),phi_ti(i),z_ti(i),s_val(j),ier)
            IF (ier < 0 .or. s_val(j) > 1.0 .or. s_val(j) < 0.0) THEN 
               s_val(j) = 1.5+0.1*j
               f_val(j) = 0.0
            END IF
            j = j + 1
         END DO
         DO i = 1, n_ti-1  ! Replace duplicates
            DO j = i+1, n_ti
               IF (s_val(i) == s_val(j)) THEN
                  s_val(j) = s_val(j) + 1.5+0.1*j
                  f_val(j) = 0.0
               END IF
            END DO
         END DO
         idx = 1
         DO i = 1, n_ti ! Sort
            dex = MINLOC(s_val,DIM=1,MASK=idx > 0)
            temp = s_val(i)
            s_val(i) = s_val(dex)
            s_val(dex) = temp
            temp = f_val(i)
            f_val(i) = f_val(dex)
            f_val(dex) = temp
            idx(i) = 0
         END DO
         dex = MINLOC(s_val, DIM=1, MASK = s_val >= 1)
         IF (dex > 2) THEN
            DO WHILE(f_val(dex-1) == 0.0)
              dex = dex - 1
            END DO
            s_val(dex) = 0.5*(1.0-s_val(dex-1))+s_val(dex-1)
            f_val(dex) = 0.5*(f_val(dex-1))
            s_val(dex) = 1.0
            f_val(dex) = 0.0
            DO i = 2, ndatafmax  ! Adjust spline knots
               IF (ti_aux_s(i) > 1.0 .or. ti_aux_s(i) < 0.0) CYCLE
               ier = 0
               s_temp = ti_aux_s(i)
               CALL eval_prof_spline(dex,s_val(1:dex),f_val(1:dex),s_temp,temp,ier)
               ti_aux_f(i) = ti_aux_f(i) - refit_param*(ti_aux_f(i) - temp)
            END DO
            ti_aux_f(1) = ti_aux_f(1) - refit_param*(ti_aux_f(1) - f_val(1))
            DO i = 2, ndatafmax
              IF (ti_aux_s(i) >= 1.0 .or. ti_aux_s(i) < 0.0) CYCLE
              IF (ti_aux_f(i) == 0.0) ti_aux_f(i) = ti_aux_f(i-1)*0.5
            END DO
         END IF
      END IF
      IF (ALLOCATED(s_val)) DEALLOCATE(s_val)
      IF (ALLOCATED(idx)) DEALLOCATE(idx)
      IF (ALLOCATED(f_val)) DEALLOCATE(f_val)
      !----------------------------------------------------------------
      !     LINE-INTEGRATED ELECTRON DENSITY (NE)
      !----------------------------------------------------------------
      IF (n_ne_line > 0) THEN
         ! Thing are about to get messy because we need to reallocate
         ! the ne_spline array
         IF (EZspline_allocated(ne_spl)) CALL EZspline_free(ne_spl,ier)
         dex = MINLOC(ne_aux_s(2:),DIM=1)
         IF (dex > 4) THEN
            CALL EZspline_init(ne_spl,dex,bcs0,ier)
            ne_spl%x1 = ne_aux_s(1:dex)
            ne_spl%isHermite = 1
            CALL EZspline_setup(ne_spl,ne_aux_f,ier)
         END IF
         ! Now refit
         ALLOCATE(s_val(n_ne_line),f_val(n_ne_line))
         j=1
         f_val = 0.0
         s_val = 0.0
         DO i =1, nprof
            IF (sigma_ne_line(i) >= bigno_ne) CYCLE
            x0(1)=r0_ne_line(i); x1(1)=r1_ne_line(i)
            x0(2)=phi0_ne_line(i); x1(2)=phi1_ne_line(i)
            x0(3)=z0_ne_line(i); x1(3)=z1_ne_line(i)
            CALL line_int(fcn_linene,x0,x1,f_val(j),LENGTH=s_val(j))
            !CALL line_int_ne(x0,x1,f_val(j),LENGTH=s_val(j))
            IF (s_val(j) > 0.0 .and. f_val(j) > 0.0) THEN
               !PRINT *,j,target_ne_line(i),f_val(j)
               f_val(j) = target_ne_line(i)/f_val(j)
               !PRINT *,j,f_val(j)
            ELSE
               f_val(j) = 0.0
            END IF
            j=j+1
         END DO
         !dex    = COUNT(f_val > 0.0)
         !ne_err = SUM(f_val)/dex - 1
         dex     = MAXLOC(target_ne_line,DIM=1)
         ne_err  = f_val(dex) - 1
         !PRINT *,dex,ne_err
         ne_aux_f(:) = ne_aux_f(:) + refit_param*ne_aux_f(:) * ne_err
      END IF
      IF (ALLOCATED(s_val)) DEALLOCATE(s_val)
      IF (ALLOCATED(f_val)) DEALLOCATE(f_val)
      !----------------------------------------------------------------
      !     Readjust the stored energy
      !----------------------------------------------------------------
      !ALLOCATE(s_val(NPRESS),f_val(NPRESS))
      !FORALL (i=1:NPRESS) s_val(i) = REAL(i-1)/REAL(NPRESS-1)
      !f_val = 0.0
      !! TE
      !dex = MINLOC(te_aux_s(2:),DIM=1)
      !IF (dex > 4) THEN
      !   DO i =1, NPRESS
      !      CALL eval_prof_spline(dex,te_aux_s(1:dex),te_aux_f(1:dex),s_val(i),temp,ier)
      !      f_val(i) = f_val(i)+temp
      !   END DO
      !END IF
      ! TI
      !dex = MINLOC(ti_aux_s(2:),DIM=1)
      !IF (dex > 4) THEN
      !   DO i =1, NPRESS
      !      CALL eval_prof_spline(dex,ti_aux_s(1:dex),ti_aux_f(1:dex),s_val(i),temp,ier)
      !      f_val(i) = f_val(i) + temp
      !   END DO
      !END IF
      ! NI
      !dex = MINLOC(ne_aux_s(2:),DIM=1)
      !IF (dex > 4) THEN
      !   DO i =1, NPRESS
      !      CALL eval_prof_spline(dex,ne_aux_s(1:dex),ne_aux_f(1:dex),s_val(i),temp,ier)
      !      f_val(i) = temp*f_val(i)
      !   END DO
      !END IF
      !new_area = SUM(f_val)
      ! Old pressure
      !DO i =1, NPRESS
      !   CALL get_equil_p(s_val(i),f_val(i),ier)
      !END DO
      !new_area = SUM(f_val)/new_area
      !PRINT *,'new_area',new_area
      !IF (ALLOCATED(s_val)) DEALLOCATE(s_val)
      !IF (ALLOCATED(f_val)) DEALLOCATE(f_val)
      
      ! Now readjust the input variables
      DO i = 1, n
         IF (var_dex(i) == ine .and. arr_dex(i,2) == norm_dex) norm_ne = x(i)
         IF (var_dex(i) == ite .and. arr_dex(i,2) == norm_dex) norm_te = x(i)
         IF (var_dex(i) == iti .and. arr_dex(i,2) == norm_dex) norm_ti = x(i)
      END DO
      DO i = 1, n
         IF (arr_dex(i,2) == norm_dex) cycle
         IF (var_dex(i) == ine_aux_f) x(i) = ne_aux_f(arr_dex(i,1))/norm_ne
         IF (var_dex(i) == ite_aux_f) x(i) = te_aux_f(arr_dex(i,1))/norm_te
         IF (var_dex(i) == iti_aux_f) x(i) = ti_aux_f(arr_dex(i,1))/norm_ti
      !   IF (var_dex(i) == ipscale)   x(i) = new_area
      END DO
      ! Now output Stuff to screen
      IF (lscreen) THEN
         WRITE(6,*)'   NE_S             NE           TE_S            TE            TI_S             TI'
         DO i = 1, ndatafmax
            IF ((te_aux_s(i) < 0) .and. (ti_aux_s(i) < 0) .and. (ne_aux_s(i) < 0)) CYCLE
            WRITE(6,'(3(2X,F6.3,2X,E20.10))') ne_aux_s(i),ne_aux_f(i)*ne_norm,te_aux_s(i),te_aux_f(i),ti_aux_s(i),ti_aux_f(i)
         END DO
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_prof_refit
