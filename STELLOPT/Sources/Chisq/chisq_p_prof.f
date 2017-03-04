!-----------------------------------------------------------------------
!     SUBROUTINE:     CHISQ_P_PROF
!
!     PURPOSE:        This subroutine calculates the chi^2 for fitting
!                     of pressure profile data.
!
!     INPUTS:         pres_opt   VMEC Pressure Profile
!                     ivar       Index for pressure fitting
!                     num        Index for pressure points
!                     nrad       Number of VMEC radial gridpoints
!                     nopt       Number of optimizations
!                     extension  File extension
!
!     OUTPUTS:        None
!
!     LIBRARIES:      lib_opt.a - kind_spec
!                               - optim_params
!                               - safe_open_mod
!                               - ajax_mod
!                               - vmec_input
!
!     WRITTEN BY:     S. Lazerson (lazerson@pppl.gov)
!                     orriginal by M. Zarnstorff
!
!     DATE:           02/03/11
!-----------------------------------------------------------------------
      subroutine chisq_p_prof (pres_opt, ivar, num, nrad, nopt,
     1                         iflag, extension)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      use stel_kinds
      use chisq_mod
      use optim_params
      use safe_open_mod
!      use AJAX_MOD
      use vmec_utils
      use read_wout_mod
      use optim, only: lajax, bigno, nfp_opt
      use EZspline
      use EZspline_obj
!-----------------------------------------------------------------------
!     Input Arguments (see above)
!-----------------------------------------------------------------------
      implicit none
      integer, intent(in) :: ivar, nrad, nopt
      integer, intent(inout) :: num, iflag
      real(rprec) :: pres_opt(*)
      character*(*) :: extension
!-----------------------------------------------------------------------
!     Local Variables
!          i,j,k        Indexing
!          iunit        File ID for output to p_prof file.
!          s            Radial surface value [0,1]
!          sj           s on half grid
!          fract        Fraction of distance on half grid
!          factor       Scaling factor for pressures
!          pi           Pi
!          ec           Charge of an electron NIST 2008 pg 34
!          angfact      Angle factor for finding per-field-period angle
!          r_cyl        Cylindrical Coordinates for AJAX
!          r_flx        Flux coordinates from AJAX
!          s_prof       s for each pressure point
!          message      output message from AJAX
!-----------------------------------------------------------------------
      integer :: i, j, k, iunit, num1, ier
      INTEGER :: bcs1(2)
      real(rprec) :: s, sj, fract, factor, pi, ec, angfact, factor_vmec,
     1               br, bphi, bz, smin, val, bigno_ne
      real(rprec), dimension(3) :: r_cyl, r_flx
      real(rprec), dimension(np_prof) :: s_prof
      real(rprec), dimension(nte_prof) :: s_te_prof
      real(rprec), dimension(nti_prof) :: s_ti_prof
      real(rprec), dimension(nne_prof) :: s_ne_prof
      character*120 :: message
      TYPE(EZspline1_r8) :: VAL_spl
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      pi = 4*atan(1._rprec)   ! Calc Pi
      ec=1.60217653e-19       ! NIST 2008 pg 34
      angfact = 2*pi/nfp_opt  ! Angle Scaling Factor
      bigno_ne = 1.0E+25      ! Sigma ne can be greater than 1E10
!     If nopt >0 calculate values otherwise do initialization
      if (nopt > 0) then
         ! Open the VMEC wout file
         CAll readw_and_open(trim(extension),iflag)
         IF (iflag .ne. 0) THEN
            iflag = -12
            RETURN
         END IF
         ! Map each point to flux space
         IF (ANY(sigma_te_prof .lt. bigno)) THEN
            DO i= 1, nte_prof
               CALL GetBcyl(r_te_prof(i),phi_te_prof(i)*pi/180,
     1                      z_te_prof(i),
     2                      br,bphi,bz,SFLX=s,INFO=ier)
               s_te_prof(i) = s
               IF (ier .eq. -1) s_te_prof(i) = 1.5
               IF (ier .eq. -3) s_te_prof(i) = 1.5
            END DO
         END IF
         IF (ANY(sigma_ne_prof .lt. bigno_ne)) THEN
            DO i= 1, nne_prof
               CALL GetBcyl(r_ne_prof(i),phi_ne_prof(i)*pi/180,
     1                      z_ne_prof(i),
     2                      br,bphi,bz,SFLX=s,INFO=ier)
               s_ne_prof(i) = s
               IF (ier .eq. -1) s_ne_prof(i) = 1.5
               IF (ier .eq. -3) s_ne_prof(i) = 1.5
            END DO
         END IF
         IF (ANY(sigma_ti_prof .lt. bigno)) THEN
            DO i= 1, nti_prof
               CALL GetBcyl(r_ti_prof(i),phi_ti_prof(i)*pi/180,
     1                      z_ti_prof(i),
     2                      br,bphi,bz,SFLX=s,INFO=ier)
               s_ti_prof(i) = s
               IF (ier .eq. -1) s_ti_prof(i) = 1.5
               IF (ier .eq. -3) s_ti_prof(i) = 1.5
            END DO
         END IF
         IF (ANY(sigma_p_prof .lt. bigno)) THEN
            DO i= 1, np_prof
               CALL GetBcyl(r_p_prof(i),phi_p_prof(i)*pi/180,
     1                      z_p_prof(i),
     2                      br,bphi,bz,SFLX=s,INFO=ier)
               s_prof(i) = s
               IF (ier .eq. -1) s_prof(i) = 1.5
               IF (ier .eq. -3) s_prof(i) = 1.5
            END DO
         END IF
         ! Open output file
         iunit = unit_outdata
         call safe_open(iunit, k, 'p_prof.'//trim(extension),
     1                  'replace', 'formatted')
         if (k .ne. 0) then
            iflag = -12
            return
         endif
         ! Now do comparrison
         IF (ANY(sigma_p_prof .lt. bigno)) THEN                         !Compare against p_spline
            write(iunit, '(a,a)', iostat=k)
     1           '    r       z      phi       s    ne-data    te-data',
     2           '    p-data    p-vmec     sigma    wgted dev.'
            DO i = 1, np_prof
               IF (sigma_p_prof(j) .ge. bigno) CYCLE
               factor_vmec = 1._rprec/MAXVAL(pres_opt(1:nrad))
!              Handle Global Indexing for optimizer
               num = num + 1
               index_array(num) = ivar
!              Global Optimizer Values
               wegt(num) = sigma_p_prof(i) * factor_p_prof
               chisq_target(num) = p_prof(i) * factor_p_prof
!              VMEC provides pressure on half mesh
               sj    = s_prof(i) * REAL(nrad-1,rprec) + 1.5
               j     = FLOOR(sj)
               fract = sj - j
!              Handle optimizing to pressure boundary
               IF (s_prof(i) > 1) THEN
                  IF (lp_prof_incl_edge) THEN
                     chisq_match(num) = 0.0
                  ELSE
                     chisq_match(num) = chisq_target(num)
                  END IF
               ELSE IF (j .eq. nrad) THEN
                  chisq_match(num) = pres_opt(nrad)
               ELSE IF (sj .lt. 2) THEN
                  chisq_match(num) = pres_opt(2)
               ELSE
                  chisq_match(num) = pres_opt(j) +
     1                            fract*(pres_opt(j+1) - pres_opt(j))
                  ! We need to make sure the pressure doesn't go negative
                  if (chisq_match(num) < 0.0) then
                     wegt(num) = wegt(num) / 100.
                  end if
               END IF
               chisq_match(num) = chisq_match(num)*factor_vmec
!              Output data to p_prof
               write(iunit, '(4f8.3,6es10.2)', iostat=k)
     1            r_p_prof(i), z_p_prof(i), phi_p_prof(i), s_prof(i),
     2            ne_prof(i),te_prof(i),p_prof(i),
     3            chisq_match(num), sigma_p_prof(i),
     4            (chisq_match(num)-chisq_target(num))
     5             /wegt(num)
            END DO
         ELSE                                                           !Compare against composite pressure
            ! Construct NE Splines
            i = minloc(ne_aux_s(2:),DIM=1)
            IF (i > 4) THEN
               WRITE(IUNIT, '(A)')'NE'
               write(iunit, '(a,a)', iostat=k)
     1           '    r       z      phi       s    ne-data    ne-fit',
     2           '    sigma    wgted dev.'
               CALL EZspline_init(VAL_spl,i,bcs1,ier)
               IF (ier /=0) stop 'ERROR: EZspline_init(NE_spl)'
               VAL_spl%isHermite = 1
               VAL_spl%x1 = ne_aux_s(1:i)
               CALL EZspline_setup(VAL_spl,ne_aux_f(1:i),ier)
               IF (ier /=0) stop 'ERROR: EZspline_setup(NE_spl)'
               DO j = 1, nne_prof
                  IF (sigma_ne_prof(j) .ge. bigno_ne) CYCLE
                  num = num + 1
                  index_array(num) = ivar_ne
                  wegt(num) = sigma_ne_prof(j)
                  chisq_target(num) = ne_prof(j)
                  val = 0.0
                  IF (s_ne_prof(j) .gt. 1.0) THEN
                     chisq_match(num) = ne_aux_f(i)
                  ELSE
                     CALL EZspline_interp(VAL_spl,s_ne_prof(j),val,ier)
                     chisq_match(num) = val
                  END IF
                  IF (val < 0.0) wegt(num) = wegt(num) / 100.
                  WRITE(iunit, '(4f8.3,4es20.10)', iostat=k)
     1            r_ne_prof(j),z_ne_prof(j),phi_ne_prof(j),s_ne_prof(j),
     2            ne_prof(j),
     3            chisq_match(num), sigma_ne_prof(j),
     4            (chisq_match(num)-chisq_target(num))
     5             /wegt(num)
               END DO
               CALL EZspline_free(VAL_spl,ier)
               IF (ier /=0) stop 'ERROR: EZspline_free(NE_spl)'
            END IF
            ! Construct TE Splines
            i = minloc(te_aux_s(2:),DIM=1)
            IF (i > 4) THEN
               WRITE(IUNIT, '(A)')'TE'
               WRITE(iunit, '(a,a)', iostat=k)
     1           '    r       z      phi       s    te-data    te-fit',
     2           '    sigma    wgted dev.'
               CALL EZspline_init(VAL_spl,i,bcs1,ier)
               IF (ier /=0) stop 'ERROR: EZspline_init(TE_spl)'
               VAL_spl%isHermite = 1
               VAL_spl%x1 = te_aux_s(1:i)
               CALL EZspline_setup(VAL_spl,te_aux_f(1:i),ier)
               IF (ier /=0) stop 'ERROR: EZspline_setup(TE_spl)'
               DO j = 1, nte_prof
                  IF (sigma_te_prof(j) .ge. bigno) CYCLE
                  num = num + 1
                  index_array(num) = ivar_te
                  wegt(num) = sigma_te_prof(j)
                  chisq_target(num) = te_prof(j)
                  val = 0.0
                  IF (s_te_prof(j) .gt. 1.0) THEN
                     chisq_match(num) = 0.0
                  ELSE
                     CALL EZspline_interp(VAL_spl,s_te_prof(j),val,ier)
                     chisq_match(num) = val
                  END IF
                  IF (val < 0.0) wegt(num) = wegt(num) / 100.
                  WRITE(iunit, '(4f8.3,4es20.10)', iostat=k)
     1            r_te_prof(j),z_te_prof(j),phi_te_prof(j),s_te_prof(j),
     2            te_prof(j),
     3            chisq_match(num), sigma_te_prof(j),
     4            (chisq_match(num)-chisq_target(num))
     5             /wegt(num)
               END DO
               CALL EZspline_free(VAL_spl,ier)
               IF (ier /=0) stop 'ERROR: EZspline_free(TE_spl)'
            END IF
            ! Construct TI Splines
            i = minloc(ti_aux_s(2:),DIM=1)
            IF (i > 4) THEN
               WRITE(IUNIT, '(A)')'TI'
               WRITE(iunit, '(a,a)', iostat=k)
     1           '    r       z      phi       s    ti-data    ti-fit',
     2           '    sigma    wgted dev.'
               CALL EZspline_init(VAL_spl,i,bcs1,ier)
               IF (ier /=0) stop 'ERROR: EZspline_init(TI_spl)'
               VAL_spl%isHermite = 1
               VAL_spl%x1 = ti_aux_s(1:i)
               CALL EZspline_setup(VAL_spl,ti_aux_f(1:i),ier)
               IF (ier /=0) stop 'ERROR: EZspline_setup(TI_spl)'
               DO j = 1, nti_prof
                  IF (sigma_ti_prof(j) .ge. bigno) CYCLE
                  num = num + 1
                  index_array(num) = ivar_ti
                  wegt(num) = sigma_ti_prof(j)
                  chisq_target(num) = ti_prof(j)
                  val = 0.0
                  IF (s_ti_prof(j) .gt. 1.0) THEN
                     chisq_match(num) = 0.0
                  ELSE
                     CALL EZspline_interp(VAL_spl,s_ti_prof(j),val,ier)
                     chisq_match(num) = val
                  END IF
                  IF (val < 0.0) wegt(num) = wegt(num) / 100.
                  WRITE(iunit, '(4f8.3,4es20.10)', iostat=k)
     1            r_ti_prof(j),z_ti_prof(j),phi_ti_prof(j),s_ti_prof(j),
     2            ti_prof(j),
     3            chisq_match(num), sigma_ti_prof(j),
     4            (chisq_match(num)-chisq_target(num))
     5             /wegt(num)
               END DO
               CALL EZspline_free(VAL_spl,ier)
               IF (ier /=0) stop 'ERROR: EZspline_free(TI_spl)'
            END IF
         END IF
         write(iunit,*)
         write(iunit,*) ' Normalization factor = ',factor_p_prof
         write(iunit,*)
         CLOSE(iunit)
         CALL read_wout_deallocate
      else
!        Add np_prof to num (optimization) and set lajax to load AJAX
         IF (np_prof .gt. 0) THEN
            DO k = 1, np_prof
               IF (sigma_p_prof(k) .ge. bigno) CYCLE
               num = num + 1
               IF (nopt .eq. -2) chisq_descript(num) = descript(ivar)
            END DO
         END IF
         IF (nne_prof .gt. 0) THEN
            DO k = 1, nne_prof
               IF (sigma_ne_prof(k) .ge. bigno_ne) CYCLE
               num = num + 1
               IF (nopt .eq. -2) chisq_descript(num) = descript(ivar_ne)
            END DO
         END IF
         IF (nte_prof .gt. 0) THEN
            DO k = 1, nte_prof
               IF (sigma_te_prof(k) .ge. bigno) CYCLE
               num = num + 1
               IF (nopt .eq. -2) chisq_descript(num) = descript(ivar_te)
            END DO
         END IF
         IF (nti_prof .gt. 0) THEN
            DO k = 1, nti_prof
               IF (sigma_ti_prof(k) .ge. bigno) CYCLE
               num = num + 1
               IF (nopt .eq. -2) chisq_descript(num) = descript(ivar_ti)
            END DO
         END IF
      end if
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      end subroutine chisq_p_prof
