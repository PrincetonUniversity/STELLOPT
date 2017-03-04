      SUBROUTINE initialize_opt(xc_opt, var_descript, nvar, nopt)
      USE optim
      USE legendre_params                                                !LEGENDRE
      USE vmec_input, ONLY : rbc, zbs, ai, am, ac, ncurr, lfreeb,
     1     curtor, extcur, phiedge, nfp, bcrit, at, ah,
     2     pmass_type,pcurr_type,piota_type,
     3     am_aux_s, am_aux_f, ai_aux_s, ai_aux_f, ac_aux_s, ac_aux_f
      USE vparams, ONLY: one, zero
      USE coilsnamin, ONLY: lmodular, lsaddle, lbcoil, lvf, lsurfv,
     1    lsadsfv, ltfc, ltfcv, lmodcur, lsadcur, bcoil_file
      USE vacfield_mod, ONLY: nstell_coils, angles, shifts, vac_extcur
      USE mpi_params                                                     !MPI
      USE stel_constants
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      include "mpif.h"
!DEC$ ENDIF
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(out) :: nvar, nopt
      REAL(rprec), DIMENSION(*) :: xc_opt
      CHARACTER(len=*), DIMENSION(*) :: var_descript
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      LOGICAL :: lneteti
      INTEGER :: nb, mb, ik, nvariables
      INTEGER :: num_coefs, start_coefs                                  !SAL
      REAL(rprec), DIMENSION(-ntord:ntord,0:mpol1d) ::
     1   rbc_in, zbs_in
      REAL(rprec), DIMENSION(1) :: fvec = (/zero/)
      REAL(rprec) :: sum1, delta
      CHARACTER(LEN=200) :: temp
      EXTERNAL convert_boundary, convert_boundary_PG,
     1         unique_boundary, unique_boundary_PG
C-----------------------------------------------
      chisq_min = 1.0E30
      min_count = 0
      lneteti   = .FALSE.
!
!     SET UP CORRESPONDENCE BETWEEN XC_OPT ARRAY (BOUNDARY COEFFICIENTS
!     USED BY OPTIMIZATION ROUTINES) AND XC ARRAY. FIRST CONVERT TO A
!     UNIQUE ANGLE REPRESENTATION AT THE BOUNDARY (BASED ON HIRSHMAN/BRESLAU)*
!
      irm0_bdy = 0;      izm0_bdy = 0;      irho_bdy = 0;      nvar = 0

      IF (rbc(0,1).eq.zero .or. zbs(0,1).eq.zero) THEN
         IF (myid .eq. master)                                           !MPI
     1       PRINT *, 'Neither RBC nor ZBS for m=0, n=1 can be ZERO',
     2                ' in STELLOPT input file!'
          STOP
      END IF

      IF (.not.lfreeb .and. .not.lvac_opt) THEN

         rbc_in = rbc
         zbs_in = zbs                                                    !copy and store m=0 components

!
!        check if conversion was made or if original bdy already in proper form
!        IF NOPT_BOUNDARY=0, USE HIRSHMAN/BRESLAU REPRESENTATION
!                        =1  USE "ADVANCED" HIRSHMAN/BRESLAU REPRESENTATION
!                        =2  USE PG (PAUL GARABEDIAN) DELTA_MN REPRESENTATION
!
         IF (nopt_boundary .le. 1) THEN
            CALL convert_boundary (rbc, zbs, rhobc, mpol1d, ntord)
            IF (.not.lniter1) THEN
               CALL unique_boundary (rbc_in, zbs_in, rhobc, mpol1d, 
     1                        ntord, mpol1_opt, ntor_opt, mrho1_opt)
               delta = SUM((rbc - rbc_in)**2)/rbc(0,1)**2
     1               + SUM((zbs - zbs_in)**2)/zbs(0,1)**2

               IF (delta.gt.1.e-8_dp .and. myid.eq.master) 
     1            WRITE (6,10) 100*(one-delta)
            END IF
 10   FORMAT(' Input boundary representation was converted!',/,
     1       ' Accuracy of conversion = ',f7.2,'%')

         DO nb = -ntor_opt, ntor_opt
            IF (rbc(nb,0).ne.zero .and. (.not.lfix_ntor(nb)) .and.
     1         (nb.ne.0 .or. ABS(sigma_rmax).le.1.e6_dp .or.
     2                       ABS(sigma_rmin).le.1.e6_dp)) THEN
               irm0_bdy = irm0_bdy + 1
               nrz0_opt(irm0_bdy) = nb
               xc_opt(irm0_bdy) = rbc(nb,0)
               WRITE (var_descript(irm0_bdy), 100) 'rbc(n=',nb,',m=0)'
            ENDIF
         END DO

 100  FORMAT(a,i4,a)

         izm0_bdy = irm0_bdy
         DO nb = -ntor_opt, ntor_opt
            IF (zbs(nb,0).ne.zero .and. (.not.lfix_ntor(nb))
     1         .and. (nb.ne.0)) THEN
               izm0_bdy = izm0_bdy + 1
               nrz0_opt(izm0_bdy) = nb
               xc_opt(izm0_bdy) = zbs(nb,0)
               WRITE (var_descript(izm0_bdy), 100) 'zbs(n=',nb,',m=0)'
            ENDIF
         END DO

         IF (izm0_bdy .gt. ntord+ntor1d) THEN
            IF (myid .eq. master)                                        !MPI
     1          PRINT *,' Only one sign of n allowed for m=0'
            STOP ' i > 2*(1+ntord) in STELLOPT initialization'
         END IF

         izm0_bdy = izm0_bdy - irm0_bdy

         DO nb = -ntor_opt, ntor_opt
            DO mb = 0, mrho1_opt
               IF (
     1              ( (rhobc(nb,mb).ne.zero .and. 
     1              ALL(lfix_rhob(-ntor_opt:ntor_opt,0:mrho1_opt)) ) 
     2             .or. .not.lfix_rhob(nb,mb) ) 
     3             .and. (mb.ne.0 .or. nb.ge.0) ) THEN
                  irho_bdy = irho_bdy + 1
                  ik = irho_bdy + irm0_bdy + izm0_bdy
                  nbrho_opt(irho_bdy) = nb
                  mbrho_opt(irho_bdy) = mb
                  xc_opt(ik) = rhobc(nb,mb)
                  WRITE (var_descript(ik), 150) 
     1                 'rhobc(n=',nb,',m=',mb,')'
               ENDIF
            END DO
         END DO

  150 FORMAT (2(a,i4),a) 

         ELSE IF (nopt_boundary.eq.2 .and. niter_opt.gt.1) THEN
            CALL convert_boundary_PG (rbc,zbs,delta_mn,mpol1d,ntord)

            DO nb = -ntord, ntord
               DO mb = -mpol1d, mpol1d
                  IF (delta_mn(nb,mb) .ne. zero .and. 
     1            .not.lfix_ntor(nb) .and. 
     2            .not.(nb .eq.0 .and. mb .eq. 0)) THEN
                     irho_bdy = irho_bdy + 1
                     ik = irho_bdy + irm0_bdy + izm0_bdy
                     nbrho_opt(irho_bdy) = nb
                     mbrho_opt(irho_bdy) = mb
                     xc_opt(ik) = delta_mn(nb,mb)
                     WRITE (var_descript(ik), 150) 
     1                 'delta_mn(n=',nb,',m=',mb,')'
                  ENDIF
               END DO
            END DO

         END IF

         nvar = nvar + irm0_bdy + izm0_bdy + irho_bdy
      

!     Allow coils to determine free-boundary shape (lfreeb = T)
      ELSE 
         IF (lvac_opt) THEN
!        VARY TILT, ANGLE OF COILS 
!        ADD EXTERNAL CURRENTS AT END (DO IN NEXT "IF" BLOCK)
            lcoil_geom = .false.
            nb = 0
            DO ik = 1, nstell_coils
               nb = nb+1
               xc_opt(nb) = angles(ik) % as_array(1)                !theta
               WRITE (var_descript(nb), 100) 'theta(group=',ik,')'
               nb = nb+1
               xc_opt(nb) = angles(ik) % as_array(2)                !phi
               WRITE (var_descript(nb), 100) 'phi(group=',ik,')'
               nb = nb+1
               xc_opt(nb) = angles(ik) % as_array(3)                !rot-angle
               WRITE (var_descript(nb), 100) 'rot_angle(group=',ik,')'
               nb = nb+1
               xc_opt(nb) = shifts(ik) % as_array(1)                !x-shift
               WRITE (var_descript(nb), 100) 'shiftx(group=',ik,')'
               nb = nb+1
               xc_opt(nb) = shifts(ik) % as_array(2)                !y-shift
               WRITE (var_descript(nb), 100) 'shifty(group=',ik,')'
               nb = nb+1
               xc_opt(nb) = shifts(ik) % as_array(3)                !z-shift
               WRITE (var_descript(nb), 100) 'shiftz(group=',ik,')'
            END DO
            nvar = nvar+nb
            nb = 0
!           FOR NOW, USE SAME lextcur ARRAY FOR VACUUM CURRENTS!
            DO ik = 1, SIZE(lextcur)
               IF (lextcur(ik)) THEN
                  nb = nb + 1
                  xc_opt(nvar+nb) = vac_extcur(ik)
                  WRITE (var_descript(nvar+nb), 100) 
     1                  'vac_extcur(i=',ik,')'
               END IF
            END DO
            nvar = nvar+nb
         END IF
         IF (.not. lcoil_geom) THEN
            nb = 0
            DO ik = 1, SIZE(lextcur)
               IF (lextcur(ik)) THEN
                  nb = nb + 1
                  xc_opt(nvar+nb) = extcur(ik)
                  WRITE (var_descript(nvar+nb), 100) 'extcur(i=',ik,')'
               END IF
            END DO
            IF (nb .ne. nextcur_opt)
     1         STOP 'Error counting EXTERNAL coils!'
            nvar = nvar + nb

         ELSE                                                            ! lcoil_geom = T

!
!     Initialize and count modular coil variables
!
            IF (lmodular) THEN
               CALL init_modular_coils (nvariables,xc_opt(nvar+1),nfp)
               DO nb = 1, nvariables
                  WRITE (var_descript(nvar+nb), 100)
     1                 'modular coil(i=',nb,')'
               END DO
               nvar = nvar + nvariables
               IF (lmodcur) THEN
                 CALL init_modular_currents (nvariables, xc_opt(nvar+1))
                 DO nb = 1, nvariables
                    WRITE (var_descript(nvar+nb), 100) 
     1                    'modular currents(i=',nb,')'
                 END DO
                 nvar = nvar + nvariables
               END IF
            END IF
!
!     Initialize and count saddle coil variables (IF lsaddle = true)
!
            IF (lsaddle) THEN
               CALL init_saddle_coils (nvariables, xc_opt(nvar+1), nfp)
               DO nb = 1, nvariables
                  WRITE (var_descript(nvar+nb), 100) 
     1                 'saddle coil(i=',nb,')'
               END DO
               nvar = nvar + nvariables
               IF (lsadcur) THEN
                  CALL init_saddle_currents (nvariables, xc_opt(nvar+1))
                  DO nb = 1, nvariables
                     WRITE (var_descript(nvar+nb), 100) 
     1                    'saddle current(i=',nb,')'
                  END DO
                  nvar = nvar + nvariables
               END IF
            END IF
!
!     Initialize and count background coil variables (IF lbcoil = true)
!
            IF (lbcoil) THEN
               IF (myid .eq. master) THEN
                  temp = copy // ".." // sep // TRIM(bcoil_file) 
     1                        // " ." // sep // TRIM(bcoil_file)
                  CALL system(temp) 
               END IF
!DEC$ IF DEFINED (MPI_OPT)
               CALL MPI_BARRIER (MPI_COMM_WORLD, ierr_mpi)
!DEC$ ENDIF
               CALL init_bg_currents (nvariables, xc_opt(nvar+1))
               DO nb = 1, nvariables
                  WRITE (var_descript(nvar+nb), 100) 
     1                 'bg current(i=',nb,')'
               END DO
               nvar = nvar + nvariables
            END IF
!
!     Initialize and count vf coil variables (IF lvf = true)
!
            IF (lvf) THEN
               CALL init_vf_currents (nvariables, xc_opt(nvar+1))
               DO nb = 1, nvariables
                  WRITE (var_descript(nvar+nb), 100) 
     1                 'vf current(i=',nb,')'
               END DO
               nvar = nvar + nvariables
            END IF

!
!     Load tf coil currents (IF ltfc = true)
!
            IF (ltfc .and. ltfcv) THEN
               CALL init_tf_coils (nvariables, xc_opt(nvar+1))
               DO nb = 1, nvariables
                  WRITE (var_descript(nvar+nb), 100) 
     1                 'tf current(i=',nb,')'
               END DO
               nvar = nvar + nvariables
            END IF
!
!     Initialize modular winding surface variables (IF lsurfv = true)
!
            IF (lsurfv) THEN
               CALL init_modular_wsurf (nvariables, xc_opt(nvar+1))
               DO nb = 1, nvariables
                  WRITE (var_descript(nvar+nb), 100) 
     1                 'modular winding surf(i=',nb,')'
               END DO
               nvar = nvar + nvariables
            END IF
!
!     Initialize saddle winding surface variables (IF lsadsfv = true)
!
            IF (lsadsfv) THEN
               CALL init_saddle_wsurf (nvariables, xc_opt(nvar+1))
               DO nb = 1, nvariables
                  WRITE (var_descript(nvar+nb), 100) 
     1                 'saddle winding surf(i=',nb,')'
               END DO
               nvar = nvar + nvariables
            END IF

         END IF               !LCOIL_GEOM

      ENDIF                   !LFREEB

!
!     COMPUTE NUMBER OF AI (NCURR=0) or AC (NCURR=1) COEFFICIENTS
!     TO VARY FOR ACHIEVING OPTIMIZED PROFILES
!     IF NCURR = 0, AI VARIES. IF NCURR = 1, AC VARIES.
!     THE target_IOTA AND/OR TARGET_CURRENT FUNCTIONS WILL CONTRIBUTE
!     TO CHISQ, DEPENDING ON THE RESPECTIVE SIGMAS.
!
      num_ai = 0; num_am = 0
     
      IF (lprof_opt) THEN                                                !Backwards compatability
         IF (ncurr .eq. 1) lcur_prof_opt = .true.
         IF (ncurr .eq. 0) liota_prof_opt = .true.
      END IF

      IF (lreconj) THEN 
         lcur_prof_opt = .true.
         liota_prof_opt = .false.
      END IF

      IF (lbootsj_opt .and. (.not.liota_prof_opt) .and.
     1    (.not.lcur_prof_opt)) THEN
         IF (myid .eq. master)                                           !MPI
     1      PRINT *,' LCUR_PROF_OPT was set = TRUE in initialize_opt'
         lcur_prof_opt = .true.
      END IF

      IF (liota_prof_opt) THEN
         ! Handle NCURR = 1 but LIOTA_PROF_OPT = T
         IF (ncurr .ne. 0) THEN
            IF (myid .eq. master) PRINT *,' NCURR IS BEING SET TO 0',    !MPI
     1      ' BECAUSE LIOTA_PROF_OPT = TRUE'
            ncurr = 0
         END IF
         IF (curtor .ne. 0.0) THEN
            IF (myid .eq. master) PRINT *,' CURTOR IS BEING SET TO 0',    !MPI
     1      ' BECAUSE LIOTA_PROF_OPT = TRUE'
            curtor = 0.0
         END IF
         lcur_prof_opt = .false.
         ! Handle piota_type
         CALL tolower(piota_type)
         SELECT CASE(TRIM(piota_type))
         CASE('sum_atan')
            STOP 'PIOTA_TYPE = sum_atan not implemented in STELLOPT'
         CASE('akima_spline','cubic_spline')
            num_coefs=minloc(ai_aux_s(2:),dim=1)                        ! Find zeros
            IF (num_coefs > 4) THEN
               DO ik = 2, num_coefs-1
                  nvar = nvar + 1
                  xc_opt(nvar) = ai_aux_f(ik)/ai_aux_f(1)
                  WRITE (var_descript(nvar), 100) 
     1               'IOTA Profile Coefficient AI_AUX_F(',ik,')'
               END DO
               ! Calculate integral Iota
               !num_coefs=minloc(ai_aux_s(2:),dim=1)                     ! Find zeros for integral
               !nvar = nvar + 1
               !xc_opt(nvar) = 0
               !DO ik = 1, num_coefs-1
               !   xc_opt(nvar) = xc_opt(nvar) 
     1         !             + ai_aux_f(ik) * 
     2         !                  ( ai_aux_s(ik+1) - ai_aux_s(ik) )
               !END DO
               !WRITE (var_descript(nvar), '(a)') 
!     1               'Iota Normalizing Factor (Integral)'
            ELSE
               STOP 'CHECK AI_AUX_S for number of values (>4)'
            END IF
         CASE('pedestal')
            STOP 'PIOTA_TYPE = pedestal not implemented in STELLOPT'
         CASE('rational')
            STOP 'PIOTA_TYPE = rational not implemented in STELLOPT'
         CASE('nice_quadratic')
            STOP 'PIOTA_TYPE=nice_quadratic not implemented in STELLOPT'
         CASE DEFAULT                                                    ! Power Series
            IF (ai(0) .eq. 0)
     1          stop 'ai(0) == 0'
            DO ik = 0,10
               IF (ai(ik) .ne. zero) num_ai = ik + 1
            END DO

            IF (l_legendre) THEN                                             !LEGENDRE
               n_leg = num_ai-1
               ALLOCATE(a_leg(0:n_leg,0:n_leg), b_leg(0:n_leg,0:n_leg),
     1         a_leg_inv(0:n_leg,0:n_leg), b_leg_inv(0:n_leg,0:n_leg),
     2         ti(0:n_leg))

               CALL build_matrices_legendre(n_leg, a_leg, b_leg,
     1         a_leg_inv, b_leg_inv)
               CALL power_to_legendre(n_leg, a_leg, b_leg, ai(0), ti)
               DO ik = 1, num_ai
                  nvar = nvar + 1                                  !LEGENDRE
                  xc_opt(nvar) = ti(ik-1)
                  WRITE (var_descript(nvar), 100) 
     1                 'iota leg coef(i=',ik-1,')'
               END DO
            ELSE
               DO ik = 2, num_ai
                  nvar = nvar + 1
                  xc_opt(nvar) = ai(ik-1)/ai(0)
                  WRITE (var_descript(nvar), 100) 
     1                 'Iota Pow Ser Coef(i=',ik-1,')'
               END DO
               ! Load Normalizing Factor (integral of AM over s)
               nvar = nvar + 1
               xc_opt(nvar) = sum(ai(0:10) /
     1                       (/ 1,2,3,4,5,6,7,8,9,10,11 /) )
               WRITE (var_descript(nvar), '(a)') 
     1               'Iota Normalizing Factor (Integral)'
            END IF
         END SELECT
      ELSE IF (lcur_prof_opt) THEN
         ! Handle NCURR = 0 but LCUR_PROF_OPT
         IF (ncurr .ne. 1) THEN
            IF (myid .eq. master) PRINT *,' NCURR IS BEING SET TO 1',    !MPI
     1      ' BECAUSE LCUR_PROF_OPT = TRUE'
            ncurr = 1
         END IF
         ! Handle pcurr_type
         CALL tolower(pcurr_type)
         SELECT CASE(TRIM(pcurr_type))
         CASE('gauss_trunc')
            STOP 'PCURR_TYPE = gauss_trunc not implemented in STELLOPT'
         CASE('two_power')
            STOP 'PCURR_TYPE = two_power not implemented in STELLOPT'
         CASE('power_series_i')
            STOP 
     1        'PCURR_TYPE = power_series_i not implemented in STELLOPT'
         CASE('akima_spline_ip', 'akima_spline_i',
     1        'cubic_spline_ip', 'cubic_spline_i')
            nvar = nvar + 1
            xc_opt(nvar) = curtor
            WRITE (var_descript(nvar), 100) 
     1                 'Total Toroidal Current (CURTOR)'
            num_coefs=minloc(ac_aux_s(2:),dim=1)                        ! Find zeros
            IF (lcur_opt_edge0) THEN
               ac_aux_f(num_coefs) = 0
               num_coefs = num_coefs -1                                 ! Don't vary edge coefficient
            END IF
            !IF (ac_aux_f(1) == 0) stop 'ERROR: AC_AUX_F(1) == 0'
            IF (num_coefs > 4) THEN
               DO ik = 1, num_coefs
                  nvar = nvar + 1
                  xc_opt(nvar) = ac_aux_f(ik)
                  WRITE (var_descript(nvar), 100) 
     1               'Current Profile Coefficient AC_AUX_F(',ik,')'
               END DO
            ELSE
               STOP 'CHECK AC_AUX_S for number of values (>4)'
            END IF
         CASE('pedestal')
            STOP 'PCURR_TYPE = pedestal not implemented in STELLOPT'
         CASE('rational')
            STOP 'PCURR_TYPE = rational not implemented in STELLOPT'
         CASE DEFAULT
            num_ai=0
            ! Get Number of elements to vary
            IF (lreconj) THEN 
               num_ai = kjj
               ac(num_ai:) = zero
            ELSE
               DO ik = 0,10                     !!Find index of last non-zero ac
                  if (ac_mask(ik) .ne. zero) num_ai = ik + 1 !SAL
               END DO
               ! Default num_ai to 11 if any ac_mask
               if (num_ai .gt. 0) num_ai = 11
               ! If ac_mask wasn't set revert to old way
               if (num_ai .eq. 0) then
                  ac_mask = 0
                  do ik = 0,10
                     if (ac(ik) .ne. zero) then
                        num_ai = ik +1
                        ac_mask(ik) = 1
                     end if
                  end do
               end if
            END IF
            ! Handle no AC values being set
            IF (num_ai .eq. 0) THEN
               num_ai = 11
               ac(0) = 1
               ac(1) =-1
               ac_mask = 0
               ac_mask(0) = 1
               ac_mask(1) = 1
            ELSE
               num_ai = 11
            END IF
!           AC(9) and AC(10) shouldn't be varried if lcur_opt_edge0
            IF (lcur_opt_edge0) THEN
               ac_mask(9)=0
               ac_mask(10)=0
            END IF
            
!           COMPUTE EDGE (INTEGRATED) CURRENT SUM1 = SUM(AC(I)/(I+1)) FOR NORMALIZATION
!           WILL MULTIPLY IN PROFIL1D ROUTINE - BY CURTOR - TO GET PHYSICAL CURRENT
!           therefore, only num_ai-1 coefficients are actually "free" to be varied
            sum1 = zero
            DO ik = 0,num_ai-1
               sum1 = sum1 + ac(ik)/(ik+1)
            END DO
            IF (sum1.ne.zero) 
     1            ac(0:num_ai-1) = (curtor/sum1) * ac(0:num_ai-1)
            IF (lj1zero) 
     1            ac(num_ai) = -SUM(ac(0:num_ai-1))                     !NOT VARIED INDEPENDENTLY
            IF (l_legendre) THEN                                             !LEGENDRE
               n_leg = num_ai-1
               IF (lj1zero) n_leg = n_leg+1
               ALLOCATE(a_leg(0:n_leg,0:n_leg), b_leg(0:n_leg,0:n_leg),
     1         a_leg_inv(0:n_leg,0:n_leg), b_leg_inv(0:n_leg,0:n_leg),
     2         tc(0:n_leg))
               CALL build_matrices_legendre(n_leg, a_leg, b_leg,
     1         a_leg_inv, b_leg_inv)
               CALL power_to_legendre(n_leg, a_leg, b_leg, ac(0), tc)
            END IF
!           make the first element be curtor, so that the feedback on the profile
!           shape is decoupled from the magnitude
            if( num_ai >= 0) then
               nvar = nvar + 1
               if(l_legendre) then
                  xc_opt(nvar) = tc(0)
                  WRITE (var_descript(nvar), 100) 
     1                 'current leg coef(i=0)'
               else
                  ac_mask(0) = 1
                  xc_opt(nvar) = curtor
                  WRITE (var_descript(nvar), 100) 
     1                 'Total Toroidal Current (CURTOR)'
               endif
            endif
            DO ik = 1,num_ai-1
               IF(l_legendre) THEN
                  nvar = nvar + 1
                  xc_opt(nvar) = tc(ik)
                  WRITE (var_descript(nvar), 100) 
     1                 'current leg coef(i=',ik,')'
               ELSE
                  if (ac_mask(ik) .ne. zero) then
                     nvar = nvar + 1
                     xc_opt(nvar) = ac(ik)
                     WRITE (var_descript(nvar), 100) 
     1                 'Current Profile Coefficient(i=',ik,')'
                  endif
               ENDIF
            END DO
         END SELECT

      ELSE IF (ncurr .eq. 0) THEN        !DO NOT match to iota target (since ai is fixed!)
         sigma_iota = bigno
         sigma_iota_max = bigno
         sigma_iota_min = bigno
      END IF                           !End lprof_opt test

!
!     ALLOW MAGNITUDE OF THE CURRENT (CURTOR) TO VARY (IF NCURR=1)
!     TO MATCH THE TARGET_MAXCURRENT VALUE. IF lcur_prof_opt IS TRUE,
!     THEN curtor WILL BE VARIED ANYHOW
!
      IF (ledge_current) ncurr = 1
      IF (ncurr.eq.1 .and. (.not.lcur_prof_opt) .and.
     1     (     ABS(sigma_maxcurrent).lt.bigno 
     2      .or. ledge_current
     3      .or. ABS(sigma_curtor) < bigno)) THEN
         nvar = nvar + 1
         xc_opt(nvar) = curtor
         WRITE (var_descript(nvar), '(a)') 
     1              'Total Toroidal Current (CURTOR)'
      END IF   
 

!     VARY PHIEDGE
      IF (lphiedge .and. phiedge_diode) THEN
!     MAKE SURE phiedge_max >= phiedge_min
         IF (phiedge_max .le. phiedge_min) THEN
            delta = phiedge_max
            phiedge_max = phiedge_min; phiedge_min = delta
         END IF
         IF (phiedge.gt.phiedge_max .or. phiedge.lt.phiedge_min)
     1      STOP 'INITIAL PHIEDGE IS OUTSIDE PRESCRIBED RANGE!'
!     SOLVE "DIODE" EQUATION FOR STARTING VALUE OF XC
         sum1 = .5_dp*(phiedge_max+phiedge_min)
         delta= .5_dp*(phiedge_max-phiedge_min)
         IF (delta .eq. zero) THEN 
            lphiedge = .false.
         ELSE
            nvar = nvar + 1
            xc_opt(nvar) = TAN(pio2*(phiedge-sum1)/delta)
            WRITE (var_descript(nvar), '(a)') 'edge toroidal flux-diode'
         END IF
      ELSE IF (lphiedge .and. .not. phiedge_diode) THEN
         nvar = nvar +1
         xc_opt(nvar) = phiedge
         WRITE(var_descript(nvar), '(a)') 'Edge Toroidal Flux (PHIEDGE)'
      END IF

!---------------------------------------------------------------------------------
!   CODE added by S. Lazerson (05/14/12) to allow variation of ER and EZ for
!   MSE reconstructions
!---------------------------------------------------------------------------------
      IF (lmse_er) THEN
         num_coefs=minloc(er_aux_s(2:),dim=1)                        ! Find zeros
         IF (num_coefs > 4) THEN
            DO ik = 1, num_coefs
               nvar = nvar + 1
               xc_opt(nvar) = er_aux_f(ik)
               WRITE (var_descript(nvar), 100) 
     1           'Radial E-Field Profile Coefficient ER_AUX_F(',ik,')'
            END DO
         ELSE
            STOP 'CHECK ER_AUX_S for number of values (>4)'
         END IF
         num_coefs=minloc(ez_aux_s(2:),dim=1)                        ! Find zeros
         IF (num_coefs > 4) THEN
            DO ik = 1, num_coefs
               nvar = nvar + 1
               xc_opt(nvar) = ez_aux_f(ik)
               WRITE (var_descript(nvar), 100) 
     1          'Vertical E-Field Profile Coefficient EZ_AUX_F(',ik,')'
            END DO
         ELSE
            STOP 'CHECK EZ_AUX_S for number of values (>4)'
         END IF
      END IF

!---------------------------------------------------------------------------------
!   CODE added by S. Lazerson (06/04/12) to allow variation of NE, TI, TE
!   MSE reconstructions
!   BIG note:  TI TE and NE use the same grid as am_aux_s
!---------------------------------------------------------------------------------
      num_coefs = 0
      num_coefs = minloc(ne_aux_s(2:), dim=1)
      IF (num_coefs > 4) THEN
!         IF (lpres_opt_edge0) THEN
!            num_coefs = num_coefs-1
!         END IF
         sum1 = SUM(ne_aux_f(1:num_coefs))
         nvar = nvar + 1
         xc_opt(nvar) = sum1
         WRITE (var_descript(nvar), 100) 
     1           'Electon Density Profile Coef. NE_AUX_F(NORM)'
         DO ik = 1, num_coefs
            nvar = nvar + 1
            xc_opt(nvar) = ne_aux_f(ik)/sum1
               WRITE (var_descript(nvar), 100) 
     1           'Electon Density Profile Coef. NE_AUX_F(',ik,')'
         END DO
         lneteti = .true.
      END IF
      num_coefs = 0
      num_coefs = minloc(te_aux_s(2:), dim=1)
      IF (num_coefs > 4) THEN
         IF (lpres_opt_edge0) THEN
            num_coefs = num_coefs-1
         END IF
         sum1 = SUM(te_aux_f(1:num_coefs))
         nvar = nvar + 1
         xc_opt(nvar) = sum1
         WRITE (var_descript(nvar), 100) 
     1           'Electon Temperature Profile Coef. TE_AUX_F(NORM)'
         DO ik = 1, num_coefs
            nvar = nvar + 1
            xc_opt(nvar) = te_aux_f(ik)/sum1
               WRITE (var_descript(nvar), 100) 
     1           'Electon Temperature Profile Coef. TE_AUX_F(',ik,')'
         END DO
         lneteti = .true.
      END IF
      num_coefs = 0
      num_coefs = minloc(ti_aux_s(2:), dim=1)
      IF (num_coefs > 4) THEN
         IF (lpres_opt_edge0) THEN
            num_coefs = num_coefs-1
         END IF
         sum1 = SUM(ti_aux_f(1:num_coefs))
         nvar = nvar + 1
         xc_opt(nvar) = sum1
         WRITE (var_descript(nvar), 100) 
     1           'Ion Temperature Profile Coef. TI_AUX_F(NORM)'
         DO ik = 1, num_coefs
            nvar = nvar + 1
            xc_opt(nvar) = ti_aux_f(ik)/sum1
               WRITE (var_descript(nvar), 100) 
     1           'Ion Temperature Profile Coef. TI_AUX_F(',ik,')'
         END DO
         lneteti = .true.
      END IF
      
      IF (lneteti) THEN
         nvar = nvar + 1
         xc_opt(nvar) = factor_p_prof
               WRITE (var_descript(nvar), '(a)') 
     1                 'Pressure Scaling (FACTOR_P_PROF)'
      END IF
         

!---------------------------------------------------------------------------------
!   CODE added by S. Lazerson (08/16/11) to allow variation of ANIMEC input
!   parameters
!
!   For now it simply assumes you want to vary the first 10 coefficients
!   Should probably use a masking array in the future.
!---------------------------------------------------------------------------------
      IF (lanimec) THEN
         IF (lani_bcrit) THEN                                             ! Critical Field Strength
            nvar = nvar + 1
            xc_opt(nvar) = bcrit
            WRITE (var_descript(nvar), '(a)') 'ANIMEC: Critical Field'
         END IF
         IF (lani_tperp) THEN                                             ! T_perp/T_|| for hot part.
            DO ik = 0, 10
               IF (at_mask(ik) .ne. 0) THEN
                  nvar = nvar + 1
                  xc_opt(nvar) = at(ik)
                  WRITE (var_descript(nvar), 100) 
     1               'ANIMEC: T_hot-perp/T_hot-|| at(',ik,')'
               END IF
            END DO
         END IF
         IF (lani_phot) THEN                                              ! P_hot/P_iso ratio
            DO ik = 0, 10
               IF (ah_mask(ik) .ne. 0) THEN
                  nvar = nvar + 1
                  xc_opt(nvar) = ah(ik)
                  WRITE (var_descript(nvar), 100) 
     1               'ANIMEC: P_hot/P_iso ah(',ik,')'
               END IF
            END DO
         END IF
      END IF
!---------------------------------------------------------------------------------
!   CODE added by R.SANCHEZ (01/19/99) to allow the mass profile to vary when
!      tailoring pressure profile for ballooning stability. Notice that when
!      tailoring is in effect (LPRES_PROF_OPT=.T.), the mass profile is varied).
!      NOTE: AM(num_am) is forced to remain equal -SUM(AM(i=0,num_am)) to guarantee p_edge=0.
!      It is ASSUMED that this condition is satisfied for the initial equilibrium!!
!
!   Generalized by Ed Lazarus (5/25/2005) to read in number of coefficients kjj to vary
!   and to decide whether to zero edge pressure or not (controlled by lp1zero)
!
!      New profile options added by Sam Lazerson (04/20/12).
!---------------------------------------------------------------------------------
      IF (lreconp) lpres_prof_opt = .true.
      IF (lpres_prof_opt) THEN
         CALL tolower(pmass_type)
         SELECT CASE(TRIM(pmass_type))
         CASE('gauss_trunc')
            STOP 'PMASS_TYPE = gauss_trunc not implemented in STELLOPT'
         CASE('two_power')
            STOP 'PMASS_TYPE = two_power not implemented in STELLOPT'
         CASE('two_lorentz')
            STOP 'PMASS_TYPE = two_lorentz not implemented in STELLOPT'
         CASE('akima_spline','cubic_spline')
            num_coefs=minloc(am_aux_s(2:),dim=1)                        ! Find zeros
            IF (lpres_opt_edge0) THEN
               am_aux_f(num_coefs) = 0                                  ! Set edge to zero
            END IF
            IF (lpres_opt_edgegr0) THEN
               am_aux_f(num_coefs-1) = am_aux_f(num_coefs)              ! Set gradient to zero
            END IF
            IF (num_coefs > 4) THEN
               ! Vary Normalized Values
               DO ik = 2, num_coefs-2
                  nvar = nvar + 1
                  xc_opt(nvar) = am_aux_f(ik)/am_aux_f(1)               ! Vary ratios to ease workload 
                  WRITE (var_descript(nvar), 100) 
     1               'Mass Profile Coefficient AM_AUX_F(',ik,')'
               END DO
               ! Vary Edge Gradient Normalize Value
               IF (.not. lpres_opt_edgegr0) THEN
                  nvar = nvar + 1
                  xc_opt(nvar) = am_aux_f(num_coefs-1)/am_aux_f(1)
                  WRITE (var_descript(nvar), 100) 
     1              'Mass Profile Coefficient AM_AUX_F(',num_coefs-1,')'
               END IF
               ! Vary Edge Normalized Value
               IF (.not. lpres_opt_edge0) THEN
                  nvar = nvar + 1
                  xc_opt(nvar) = am_aux_f(num_coefs)/am_aux_f(1)
                  WRITE (var_descript(nvar), 100) 
     1              'Mass Profile Coefficient AM_AUX_F(',num_coefs,')'
               END IF
            ELSE
               STOP 'CHECK AM_AUX_S for number of values (>4)'
            END IF
            ! Calculate integral pressure
            num_coefs=minloc(am_aux_s(2:),dim=1)                     ! Find zeros for integral
            nvar = nvar + 1
            xc_opt(nvar) = 0
            DO ik = 1, num_coefs-1
               xc_opt(nvar) = xc_opt(nvar) 
     1                      + am_aux_f(ik) * 
     2                           ( am_aux_s(ik+1) - am_aux_s(ik) )
            END DO
            WRITE (var_descript(nvar), '(a)') 
     1               'Pressure Normalizing Factor (Integral)'
            ! Match Beta
!            IF (     abs(sigma_beta)    < bigno 
!     1          .or. abs(sigma_eplasma) < bigno) THEN
!               nvar = nvar + 1
!               xc_opt(nvar) = factor_p_prof*am_aux_f(1)
!               WRITE (var_descript(nvar), '(a)') 
!     1                 'Pressure Scaling (FACTOR_P_PROF)'
!            END IF
         CASE('pedestal')
            IF( am(0) == 0) STOP "am(0) = 0, no central pressure!"
            nvar = nvar + 1
            xc_opt(nvar) = am(0)
            ik = 0
            WRITE (var_descript(nvar), 100) 
     1            'Pedestal Mass Profile Coefficient AM(',ik,')'
            DO ik = 1, min(15, pres_opt_nmax)
               nvar = nvar + 1
               xc_opt(nvar) = am(ik) / am(0)
               WRITE (var_descript(nvar), 100) 
     1               'Pedestal Mass Profile Coefficient AM(',ik,')'
            END DO
            ! AM(16) isn't utilized
            ! AM(20) is set by code so don't vary
            DO ik = 17, 19
               nvar = nvar + 1
               xc_opt(nvar) = am(ik)
               WRITE (var_descript(nvar), 100) 
     1               'Pedestal Mass Profile Coefficient AM(',ik,')'
            END DO
            !STOP 'PMASS_TYPE = pedestal not implemented in STELLOPT'
         CASE('rational')
            STOP 'PMASS_TYPE = rational not implemented in STELLOPT'
         CASE DEFAULT                                                   ! power_series
            IF(.not. l_legendre) THEN                                   
               IF( am(0) == 0) STOP "am(0) = 0, no central pressure!"
               ! Load First 8 normalized coefficients (or lower values)
               DO ik=1, min(8, pres_opt_nmax)
                  nvar = nvar + 1
                  xc_opt(nvar) = am(ik) / am(0)
                  WRITE (var_descript(nvar), 100) 
     1               'Mass Profile Coefficient AM(',ik,')'
               END DO
               ! Load AM(9)/AM(0) 
               IF( .not.(lpres_opt_edge0 .and. lpres_opt_edgegr0) .and.
     1            pres_opt_nmax >= 9) THEN
                  nvar = nvar + 1
                  xc_opt(nvar) = am(9) / am(0)
                  WRITE (var_descript(nvar), 100) 
     1               'Mass Profile Coefficient AM(',9,')'
               END IF
               ! Load AM(10)/AM(0)
               IF( .not.(lpres_opt_edge0 .or. lpres_opt_edgegr0) .and.
     1            pres_opt_nmax >= 10) THEN
                  nvar = nvar + 1
                  xc_opt(nvar) = am(10) / am(0)
                  WRITE (var_descript(nvar),  100) 
     1               'Mass Profile Coefficient AM(',10,')'
               END IF
               ! Load Normalizing Factor (integral of AM over s)
               nvar = nvar + 1
               xc_opt(nvar) = sum(am(0:10) /
     1                       (/ 1,2,3,4,5,6,7,8,9,10,11 /) )
               WRITE (var_descript(nvar), '(a)') 
     1               'Pressure Normalizing Factor (Integral)'
            
            ELSE                                                        
               n_leg = num_am-1
               IF (lp1zero) n_leg = n_leg+1
               ALLOCATE(a_leg(0:n_leg,0:n_leg), b_leg(0:n_leg,0:n_leg),
     1         a_leg_inv(0:n_leg,0:n_leg), b_leg_inv(0:n_leg,0:n_leg),
     2         tm(0:n_leg))

               CALL build_matrices_legendre(n_leg, a_leg, b_leg,
     1         a_leg_inv, b_leg_inv)
               CALL power_to_legendre(n_leg, a_leg, b_leg, am(0), tm)

               DO ik=0, num_am-1
                  nvar = nvar + 1
                  xc_opt(nvar) = tm(ik)
                  WRITE (var_descript(nvar), 100) 
     1                 'mass leg ser coef(i=',ik,')'
               ENDDO
            ENDIF
            ! Scaling Factor for beta matching
            IF (     abs(sigma_beta)    < bigno 
     1          .or. abs(sigma_eplasma) < bigno) THEN
               nvar = nvar + 1
               xc_opt(nvar) = factor_p_prof*am(0)
               WRITE (var_descript(nvar), '(a)') 
     1                 'Pressure Scaling (FACTOR_P_PROF)'
            END IF
         END SELECT
      END IF
!-------------------------------------------------------------------------------
      IF (nvar .gt. nxc) STOP 'nvar>nxc --> xc_opt out of bounds'

!------------------------------------------------------------------------------
!     COMPUTE NOPT (NO. OF OPTIMIZATION FUNCTIONALS)
!------------------------------------------------------------------------------
      nopt = -1
      CALL load_target (xc_opt, fvec, opt_ext, nopt, 0,  ik, .true.)
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr_mpi)
      IF (ierr_mpi .ne. 0) STOP 'MPI_BARRIER error in STELLOPT Osetup'
!DEC$ ENDIF
!------------------------------------------------------------------------------
!     LOAD chisq_descript ARRAY ELEMENTS (need for master only)
!------------------------------------------------------------------------------
      IF (myid .eq. master) THEN
         nopt = -2
         CALL load_target (xc_opt, fvec, opt_ext, nopt, 0,  ik, .true.)
      END IF
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr_mpi)
      IF (ierr_mpi .ne. 0) STOP 'MPI_BARRIER error in STELLOPT Osetup'
!DEC$ ENDIF

      IF (myid .eq. master) PRINT 200, nvar, nopt                        !MPI
  200 FORMAT(/' No. Independent Variables = ',i5/,
     1   ' No. Dependent Constraints = ',i5/)
      CALL flush(6)

      END SUBROUTINE initialize_opt
