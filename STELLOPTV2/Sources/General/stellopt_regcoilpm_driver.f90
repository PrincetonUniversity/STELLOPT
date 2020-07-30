!-----------------------------------------------------------------------
!     Subroutine:    stellopt_regcoilpm_driver
!     Authors:       J.C.Schmitt (Auburn/PPPL) jcschmitt@auburn.edu
!     Date:          2017-2018
!     Description:   This subroutine calls the coil regularization code
!                    REGCOILPM
!                    
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_regcoilpm_driver(file_str, lscreen, iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_input_mod
      USE stellopt_vars, my_mpol => mpol_rcpmds, my_ntor => ntor_rcpmds
      USE equil_utils
      USE neswrite, ONLY: coil_separation

!DEC$ IF DEFINED (REGCOILPM)
      USE regcoil_init_plasma_mod
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
      ! FOR REGCOILPM
      INTEGER :: istat, iunit = 82, m, n, ii, imn, nummodes1, nummodes2

!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
!      IF (iflag < 0) RETURN
      !lscreen = .true.
      IF (lscreen) then
         WRITE(6,'(a)') ' -------------  REGCOILPM CALCULATION  ---------'
      ENDIF
      ! WRITE(6,'(a,a)') '<---- proc_string=', proc_string
      ! WRITE(6,'(a,i4.2)') ' -------------  REGCOILPM: iflag=', iflag
!DEC$ IF DEFINED (REGCOILPM)
      verbose = lscreen ! Suppress REGCOILPM stdout if needed

      !write(6,*) "proc_str=",proc_string," file_str=",file_str

      ! Run bnorm if required
      if (load_bnorm) then
         coil_separation = regcoilpm_dipole_surface_separation 
         ! Run BNORM code
         call stellopt_bnorm(file_str,lscreen)
         bnorm_filename = 'bnorm.' // TRIM(file_str)
      end if

      ! IF (lscreen) WRITE(6,'(a,a)') '<---- proc_string=', proc_string
      wout_filename = 'wout_'//TRIM(proc_string)//'.nc'
      ! STELLOPT (via lmdif->stellopt_fcn or similar) will modifiy the value of
      ! regcoilpm_winding_surface_separation, current_density, and/or the
      ! boundary coefficients. Here, the value of the REGCOIL variables
      ! are loaded with the new values, the correct mode of operation is
      ! determiend, control variables are set, and the regcoil-functions
      ! are called
      separation = regcoilpm_dipole_surface_separation
      target_value = regcoilpm_target_value
     
      ! Loop over all of the spectral components of the winding surface
      ! and update the rc_*_stellopt  
      ! write(6,'(a)') '<----Looping over m and n'

      IF ((ANY(lregcoilpm_ds_rbound_s_opt)) .or. (ANY(lregcoilpm_ds_rbound_c_opt)) .or. &
          (ANY(lregcoilpm_ds_zbound_s_opt)) .or. (ANY(lregcoilpm_ds_zbound_c_opt)) ) THEN 
         nummodes1 = 0
         DO m = -my_mpol, my_mpol
             DO n = -my_ntor, my_ntor
                IF ( (regcoilpm_ds_rbound_c(m,n) .ne. 0) .or. &
                     (regcoilpm_ds_rbound_s(m,n) .ne. 0) .or. &
                     (regcoilpm_ds_zbound_c(m,n) .ne. 0) .or. &
                     (regcoilpm_ds_zbound_s(m,n) .ne. 0) ) THEN
                   nummodes1 = nummodes1 + 1
                END IF
             END do
         END do
         CALL safe_open(iunit, istat, TRIM('regcoilpm_nescout.'// &
                   TRIM(proc_string)), 'replace', 'formatted')
         !write(6,'(a)'), '<----JCSwrite_output'
         write (iunit, '(a)') '------ Plasma information from VMEC ----'
         write (iunit, '(a)') 'np     iota_edge       phip_edge       curpol'
         ! write nfp and curpol information 
         write (iunit, '(I6, 3ES20.12)') nfp, 0.0, 0.0, curpol  
         write (iunit,*)
         write (iunit, '(a, 1pe20.12, a)') '------ Current Surface: Dipole-Plasma separation = ', separation,' -----'
         write (iunit, '(a)') 'Number of fourier modes in table'
         write (iunit,*) nummodes1
         write (iunit, '(a)') 'Table of fourier coefficients'
         write (iunit, '(a)') 'm,n,crc2,czs2,crs2,czc2'

         DO m = -my_mpol, my_mpol
             DO n = -my_ntor, my_ntor
                if ( (regcoilpm_ds_rbound_c(m,n) .ne. 0) .or. &
                     (regcoilpm_ds_rbound_s(m,n) .ne. 0) .or. &
                     (regcoilpm_ds_zbound_c(m,n) .ne. 0) .or. &
                     (regcoilpm_ds_zbound_s(m,n) .ne. 0) ) THEN
                   ! These are written in the same order as in a NESCIN
                   ! file: M N RC ZS RS ZC
                   write(iunit,*) m, n, &
                        regcoilpm_ds_rbound_c(m,n), regcoilpm_ds_zbound_s(m,n), &
                        regcoilpm_ds_rbound_s(m,n), regcoilpm_ds_zbound_c(m,n)
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
             rmnc_dipole(imn) = regcoilpm_ds_rbound_c(m,n)
             rmns_dipole(imn) = regcoilpm_ds_rbound_s(m,n)
             zmnc_dipole(imn) = regcoilpm_ds_zbound_c(m,n)
             zmns_dipole(imn) = regcoilpm_ds_zbound_s(m,n)
          END DO
          CLOSE(iunit)
      END IF

      ! regcoil will overwrite nlambda each time - need to restore it to
      ! the original value here
      nlambda = regcoilpm_nlambda

      ! This should *almost* be a duplicate of the main code from
      ! regcoil_pm/regcoil.f90

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
      IF ( (lregcoilpm_dipole_surface_separation_opt) .and. &
           ((ANY(lregcoilpm_ds_rbound_s_opt)) .or. (ANY(lregcoilpm_ds_rbound_c_opt)) .or. &
            (ANY(lregcoilpm_ds_zbound_s_opt)) .or. (ANY(lregcoilpm_ds_zbound_c_opt))) ) THEN 
        write(6,'(a)') 'K====-----<<<<REGCOILPM ERROR: Do not optimize both separation AND Fourier series simultaneously'
      END IF

      call regcoil_init_coil_surface()

      ! Initialize some of the vectors and matrices needed:
      call regcoil_init_ports()
      ! write(6,'(a)') '<----read bnorm'
      IF (load_bnorm) THEN
         call stellopt_bnorm(proc_string,lscreen)
         bnorm_filename = 'bnorm.' // TRIM(proc_string)
      ENDIF
      call regcoil_read_bnorm()
      call regcoil_init_basis_functions()
      call regcoil_build_matrices()
      call regcoil_prepare_solve()


!      ! calculate toroidal field contribution
!      IF (lregcoil_toroidal_field) THEN
!         CALL regcoil_toroidal_field()
!         net_poloidal_current_Amperes = 0
!      ENDIF
!
!      !Initialize sensitivity arrays
!      if (sensitivity_option > 1) then
!         call regcoil_init_sensitivity()
!      endif
!
!      ! write(6,'(a)') '<----build matrices'
!      CALL regcoil_build_matrices()
!      CALL regcoil_prepare_solve()

      select case (trim(lambda_option))
      case (lambda_option_single,lambda_option_scan)
         call regcoil_lambda_scan()
      case (lambda_option_search)
         call regcoil_auto_regularization_solve()
      case default
         print *,"Invalid lambda_option:",lambda_option
         stop
      end select

 
!      select case (general_option)
!      case (1)
!         call regcoil_lambda_scan()
!         if (sensitivity_option > 1) then
!            call regcoil_adjoint_solve()
!         end if
!         call regcoil_optimal_lambda()
!      case (2)
!         call regcoil_compute_diagnostics_for_nescout_potential()
!      case (3)
!         call regcoil_svd_scan()
!      case (4,5)
!         call regcoil_auto_regularization_solve()
!         if (sensitivity_option > 1 .and. exit_code == 0) then
!            call regcoil_adjoint_solve()
!         end if
! case default
!         print *,"Invalid general_option:",general_option
!         stop
!      END select

      output_filename = 'regcoilpm_out.'// TRIM(proc_string)//'.nc'
      !call regcoil_write_output() ! For debugging it can be useful to write the regcoil_out file.

      call regcoil_evaluate_outer_surface()
      call regcoil_write_output()

      if (write_mgrid) call regcoil_write_mgrid()

      IF (lscreen) WRITE(6,'(a)') ' ---------------------------  REGCOILPM CALCULATION DONE  ---------------------'
!DEC$ ENDIF
      RETURN

!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_regcoilpm_driver
