!-----------------------------------------------------------------------
!     Subroutine:    stellopt_regcoil_driver
!     Authors:       J.C.Schmitt (Auburn/PPPL) jcschmitt@auburn.edu
!     Date:          2017-2018
!     Description:   This subroutine calls the coil regularization code
!                    REGCOIL in 'target <mode>', where target_otion is one of
!                    the following:
!    "max_K": Search for the λ value such that the maximum current 
!             density over the winding surface equals target value.
!    "chi2_K": Search for the λ value such that
!             χ 2 K equals target value.
!    "rms_K": Search for the λ value such that the
!             root-mean-square current density R d 2a K2 1/2
!             (where the integral is over the current winding surface)
!             equals target value
!    "max_Bnormal": Search for the λ value such that
!             the maximum B · n over the plasma surface equals target value.
!    "chi2_B": Search for the λ value such that
!             χ 2 B equals target value
!    "rms_Bnormal": Search for the λ value such that the 
!             root-mean-square value of B·n, i.e. R d 2a B2 n 1/2 
!             (where the integral is over the plasma surface)
!             equals target value
!                    
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_regcoil_driver(file_str, lscreen, iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_input_mod
      USE stellopt_vars, my_mpol => mpol_rcws, my_ntor => ntor_rcws
      USE equil_utils, equil_nfp => nfp
      USE neswrite, ONLY: coil_separation

!DEC$ IF DEFINED (REGCOIL)
      !USE regcoil_auto_regularization_solve
      !USE regcoil_build_matrices
      !USE regcoil_compute_lambda 
      !USE regcoil_init_coil_surface
      USE regcoil_init_plasma_mod
      !USE regcoil_initupdate_nescin_coil_surface
      !USE regcoil_prepare_solve
      !USE regcoil_read_bnorm
      !USE regcoil_read_nescin_spectrum
      !USE regcoil_validate_input
      USE regcoil_variables
!DEC$ ENDIF

!-----------------------------------------------------------------------
!     Subroutine Parameters
!        iflag         Error flag
!----------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(256), INTENT(inout)    :: file_str
      LOGICAL, INTENT(inout)        :: lscreen
      INTEGER, INTENT(inout) :: iflag

!-----------------------------------------------------------------------
!     Local Variables
!        iunit       File unit number
!----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Local Variables
!        istat         Error status
!        iunit         File unit number
      ! FOR REGCOIL
      INTEGER :: istat, iunit = 12, m, n, ii, imn, nummodes1, nummodes2

!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
!      IF (iflag < 0) RETURN
      !lscreen = .true.
      IF (lscreen) then
         WRITE(6,'(a)') ' -------------  REGCOIL CALCULATION  ---------'
      ENDIF
      ! WRITE(6,'(a,a)') '<---- proc_string=', proc_string
      ! WRITE(6,'(a,i4.2)') ' -------------  REGCOIL: iflag=', iflag
!DEC$ IF DEFINED (REGCOIL)
      verbose = lscreen ! Suppress REGCOIL stdout if needed

      !write(6,*) "proc_str=",proc_string," file_str=",file_str

      ! Run bnorm if required
      if (load_bnorm) then
         coil_separation = regcoil_winding_surface_separation 
         ! Run BNORM code
         call stellopt_bnorm(file_str,lscreen)
         bnorm_filename = 'bnorm.' // TRIM(file_str)
      end if

      ! IF (lscreen) WRITE(6,'(a,a)') '<---- proc_string=', proc_string
      wout_filename = 'wout_'//TRIM(proc_string)//'.nc'
      ! STELLOPT (via lmdif->stellopt_fcn or similar) will modifiy the value of
      ! regcoil_winding_surface_separation, current_density, and/or the
      ! boundary coefficients. Here, the value of the REGCOIL variables
      ! are loaded with the new values, the correct mode of operation is
      ! determiend, control variables are set, and the regcoil-functions
      ! are called
      separation = regcoil_winding_surface_separation
      target_value = regcoil_target_value
     
      ! Loop over all of the spectral components of the winding surface
      ! and update the rc_*_stellopt  
      ! write(6,'(a)') '<----Looping over m and n'

      IF ((ANY(lregcoil_rcws_rbound_s_opt)) .or. (ANY(lregcoil_rcws_rbound_c_opt)) .or. &
          (ANY(lregcoil_rcws_zbound_s_opt)) .or. (ANY(lregcoil_rcws_zbound_c_opt)) ) THEN 
         nummodes1 = 0
         DO m = -my_mpol, my_mpol
             DO n = -my_ntor, my_ntor
                IF ( (regcoil_rcws_rbound_c(m,n) .ne. 0) .or. &
                     (regcoil_rcws_rbound_s(m,n) .ne. 0) .or. &
                     (regcoil_rcws_zbound_c(m,n) .ne. 0) .or. &
                     (regcoil_rcws_zbound_s(m,n) .ne. 0) ) THEN
                   nummodes1 = nummodes1 + 1
                END IF
             END do
         END do
         CALL safe_open(iunit, istat, TRIM('regcoil_nescout.'// &
                   TRIM(proc_string)), 'replace', 'formatted')
         !write(6,'(a)'), '<----JCSwrite_output'
         write (iunit, '(a)') '------ Plasma information from VMEC ----'
         write (iunit, '(a)') 'np     iota_edge       phip_edge       curpol'
         ! write nfp and curpol information 
         write (iunit, '(I6, 3ES20.12)') nfp, 0.0, 0.0, curpol  
         write (iunit,*)
         write (iunit, '(a, 1pe20.12, a)') '------ Current Surface: Coil-Plasma separation = ', separation,' -----'
         write (iunit, '(a)') 'Number of fourier modes in table'
         write (iunit,*) nummodes1
         write (iunit, '(a)') 'Table of fourier coefficients'
         write (iunit, '(a)') 'm,n,crc2,czs2,crs2,czc2'

         DO m = -my_mpol, my_mpol
             DO n = -my_ntor, my_ntor
                if ( (regcoil_rcws_rbound_c(m,n) .ne. 0) .or. &
                     (regcoil_rcws_rbound_s(m,n) .ne. 0) .or. &
                     (regcoil_rcws_zbound_c(m,n) .ne. 0) .or. &
                     (regcoil_rcws_zbound_s(m,n) .ne. 0) ) THEN
                   !rc_rmnc_stellopt(m,n) = regcoil_rcws_rbound_c(m,n)
                   !rc_rmns_stellopt(m,n) = regcoil_rcws_rbound_s(m,n)
                   !rc_zmnc_stellopt(m,n) = regcoil_rcws_zbound_c(m,n)
                   !rc_zmns_stellopt(m,n) = regcoil_rcws_zbound_s(m,n)
                   ! These are written in the same order as in a NESCIN
                   ! file: M N RC ZS RS ZC
                   write(iunit,*) m, n, &
                        regcoil_rcws_rbound_c(m,n), regcoil_rcws_zbound_s(m,n), &
                        regcoil_rcws_rbound_s(m,n), regcoil_rcws_zbound_c(m,n)
                        !rc_rmnc_stellopt(m,n), rc_zmns_stellopt(m,n), &
                        !rc_rmns_stellopt(m,n), rc_zmnc_stellopt(m,n)
                END IF
              END DO
          END DO
          DO imn = 1, mnmax_coil
             m = xm_coil(imn)
             n = xn_coil(imn)/(-nfp)
             IF (m < -my_mpol .or. m > my_mpol .or. n < -my_ntor .or. n > my_ntor) THEN
                WRITE(6,*) "Error! (m,n) in regcoil coil surface exceeds mpol_rcws or ntor_rcws."
                STOP
             END IF
             rmnc_coil(imn) = regcoil_rcws_rbound_c(m,n)
             rmns_coil(imn) = regcoil_rcws_rbound_s(m,n)
             zmnc_coil(imn) = regcoil_rcws_zbound_c(m,n)
             zmns_coil(imn) = regcoil_rcws_zbound_s(m,n)
          END DO
          CLOSE(iunit)
      END IF

      ! regcoil will overwrite nlambda each time - need to restore it to
      ! the original value here
      nlambda = regcoil_nlambda

      ! This should *almost* be a duplicate of the main code from
      ! regcoil.f90

      ! Validate the input
      ! write(6,'(a)') '<----Validate'
      call regcoil_validate_input()
      ! write(6,'(a)') '<----Compute lambda'
      call regcoil_compute_lambda()

      ! Define the position vector and normal vector at each grid point for
      ! the surfaces:
      ! write(6,'(a)') '<----init_plasma'
      call regcoil_init_plasma()
      ! write(6,'(a)') '<----init coil surfs'
      IF ( (lregcoil_winding_surface_separation_opt) .and. &
           ((ANY(lregcoil_rcws_rbound_s_opt)) .or. (ANY(lregcoil_rcws_rbound_c_opt)) .or. &
            (ANY(lregcoil_rcws_zbound_s_opt)) .or. (ANY(lregcoil_rcws_zbound_c_opt))) ) THEN 
        write(6,'(a)') 'K====-----<<<<REGCOIL ERROR: Do not optimize both separation AND Fourier series simultaneously'
      END IF

      IF (lregcoil_winding_surface_separation_opt) then 
         ! write(6,'(a)') '<----regcoil_init_coil_surface'
         call regcoil_init_coil_surface()
      END IF

      IF ((ANY(lregcoil_rcws_rbound_s_opt)) .or. (ANY(lregcoil_rcws_rbound_c_opt)) .or. &
          (ANY(lregcoil_rcws_zbound_s_opt)) .or. (ANY(lregcoil_rcws_zbound_c_opt)) ) THEN 
         ! write(6,'(a)') '<----regcoil initupdate_nescin_coil_surface'
         !call regcoil_initupdate_nescin_coil_surface(verbose)
         CALL regcoil_evaluate_coil_surface()
      END IF

      ! Initialize some of the vectors and matrices needed:
      ! write(6,'(a)') '<----read bnorm'
      IF (load_bnorm) THEN
         call stellopt_bnorm(proc_string,lscreen)
         bnorm_filename = 'bnorm.' // TRIM(proc_string)
      ENDIF

      call regcoil_read_bnorm()
      ! write(6,'(a)') '<----build matrices'
      call regcoil_build_matrices()
      call regcoil_prepare_solve()

      ! JCS: I disabled all options except for #5 (for now)
      ! As REGCOIL development continues, future cases can 
      ! be handled with case statements here.
      ! write(6,'(a)') '<----select a case'
      select case (general_option)
      !case (1)
      !   call regcoil_solve()
      !case (2)
      !   call regcoil_compute_diagnostics_for_nescout_potential()
      !case (3)
      !   call regcoil_svd_scan()
      !case (4)
      !   call regcoil_auto_regularization_solve()
      case (5)
         ! write(6,'(a)') '<----auto_reg solve'
         call regcoil_auto_regularization_solve()
         ! After regcoil_auto_reg...solve, the value we want should be
         ! in the variable 'chi2_B_target'. Normal termination of regcoil
         ! returns the achieved chi2_B (miniumum). If there is an
         ! 'error' (too high or too low of current), the chi2_B will
         ! contain the chi2_B that was achieved with infinite
         ! regularization ! (well-spaced apart, straight-ish) coils
         ! See regcoil_auto_regularization_solve.f90 for the assignment
         ! of variables for external optimizers
         !     chi2_B_target, max_K_target, rms_K_target,
         !     max_Bnormal_target, chi2_K_target,
         !     coil_plasma_dist_min_target, Bnormal_total_target,
         !     + volume and area targets
 case default
         print *,"Invalid general_option:",general_option
         stop
      END select

      output_filename = 'regcoil_out.'// TRIM(proc_string)//'.nc'
      !call regcoil_write_output() ! For debugging it can be useful to write the regcoil_out file.

      !=================DEBUG SECTION==========================
      !  This section is for debugging purposes. If
      !  uncommented, the regcoil input files should be written to the
      !  working job directory and an output statement(s) will be displaed
      !  to the screen, regardless of the value of 'lscreen'
      !  - JCS
      ! IF ((ANY(lregcoil_rcws_rbound_s_opt)) .or. (ANY(lregcoil_rcws_rbound_c_opt)) .or. &
      !    (ANY(lregcoil_rcws_zbound_s_opt)) .or. (ANY(lregcoil_rcws_zbound_c_opt)) ) THEN 
      ! write(6,'(a)') '<----REGCOIL DEBUG safe_open'
      ! CALL safe_open(iunit, istat, TRIM('regcoil_nescout.'// &
      !           TRIM(proc_string)), 'replace', 'formatted')
      ! write(6,'(a)'), '<----JCSwrite_output'
      ! write(iunit,*), "Number of fourier modes in table"
      ! write(iunit,*), nummodes1
      ! write(iunit,*), "Table of fourier coefficients"
      ! write(iunit,*), "m,n,crc2,czs2,crs2,czc2"

      ! ii = 0
      ! nummodes2 = 0
      ! DO m = -mpol_rcws, mpol_rcws
      !     DO n = -ntor_rcws, ntor_rcws
      !        if ( (regcoil_rcws_rbound_c(m,n) .ne. 0) .or. &
      !             (regcoil_rcws_rbound_s(m,n) .ne. 0) .or. &
      !             (regcoil_rcws_zbound_c(m,n) .ne. 0) .or. &
      !             (regcoil_rcws_zbound_s(m,n) .ne. 0) ) THEN
      !          ii = ii+1
      !          rc_xm_stellopt(ii) = m
      !          rc_xn_stellopt(ii) = n
      !          rc_rmnc_stellopt(ii) = regcoil_rcws_rbound_c(m,n)
      !          rc_rmns_stellopt(ii) = regcoil_rcws_rbound_s(m,n)
      !          rc_zmnc_stellopt(ii) = regcoil_rcws_zbound_c(m,n)
      !          rc_zmns_stellopt(ii) = regcoil_rcws_zbound_s(m,n)
      !          if ( (rc_rmnc_stellopt(ii) .ne. 0) .or.  (rc_rmns_stellopt(ii) .ne. 0) .or. &
      !               (rc_zmnc_stellopt(ii) .ne. 0) .or.  (rc_zmns_stellopt(ii) .ne. 0) ) THEN
      !             nummodes2 = nummodes2 + 1
      !             write(iunit,*), m, n, &
      !                    rc_rmnc_stellopt(m,n), rc_rmns_stellopt(m,n), &
      !                    rc_zmnc_stellopt(m,n), rc_zmns_stellopt(m,n)
      !          END if
      !      END DO
      ! END DO
      ! if (nummodes1 .ne. nummodes2) then
      !   write(6,'(a)') '<----Hmmm....different number of modes???'
      ! END if
      ! END if
      ! print *, chi2_B_target
      ! print *,"REGCOIL complete. Total time=",totalTime,"sec."
      !=================END OF DEBUG SECTION==========================

      IF (lscreen) WRITE(6,'(a)') ' ---------------------------  REGCOIL CALCULATION DONE  ---------------------'
!DEC$ ENDIF
      RETURN

!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_regcoil_driver
