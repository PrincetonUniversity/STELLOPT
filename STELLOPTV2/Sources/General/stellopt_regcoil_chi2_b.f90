!-----------------------------------------------------------------------
!     Subroutine:    stellopt_regcoil_chi2_b
!     Authors:       J.C.Schmitt (Auburn/PPPL) jcschmitt@auburn.edu
!     Date:          2017
!     Description:   This subroutine calls the coil regularization code
!                    REGCOIL in 'target sqrt(<K^2>)' mode
!                    
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_regcoil_chi2_b(lscreen, iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_input_mod
      USE stellopt_vars
      USE equil_utils
!      USE neswrite, ONLY: coil_separation

!DEC$ IF DEFINED (REGCOIL)
      ! REGCOIL files
      USE regcoil_variables
      USE regcoil_input_mod
      USE validate_regcoil_input
      USE compute_regcoil_lambda
      USE init_regcoil_plasma
      USE init_regcoil_coil_surface
      USE read_regcoil_bnorm
      USE build_regcoil_matrices
      USE regcoil_auto_regularization_solve
      USE write_regcoil_output
!DEC$ ENDIF

!-----------------------------------------------------------------------
!     Subroutine Parameters
!        iflag         Error flag
!----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(inout) :: iflag
      LOGICAL, INTENT(inout)        :: lscreen

!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Local Variables
!        iverb         REGCOIL screen control
!        istat         Error status
!        iunit         File unit number
      INTEGER ::  ier, iunit_rzuv
      ! FOR REGCOIL
      ! INTEGER(4)     :: regcoiloutTEMP,regcoilScrOut
      LOGICAL :: lexists
      INTEGER :: iverb, istat, nu, nv, mf, nf, md, nd, iunit, m, n, &
                 ivmec, ispline_file
!      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: bnfou, bnfou_c
      CHARACTER(8)   :: temp_str
      CHARACTER(256) :: copt_fext
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
!      IF (iflag < 0) RETURN
      IF (lscreen) then
         WRITE(6,'(a)') ' -------------  REGCOIL CALCULATION  ---------'
      ENDIF
!DEC$ IF DEFINED (REGCOIL)

      !IF (lscreen) WRITE(6,'(a,a)') '<---- proc_string=', proc_string
      wout_filename = 'wout_'//TRIM(proc_string)//'.nc'
      separation = regcoil_winding_surface_separation
      current_density_target = regcoil_current_density
      ! regcoil will overwrite nlambda - need to restore it to the
      ! original value here
      nlambda = nlambda_regcoil
      ! write(6,'(a)') '<----safe_open'
      CALL safe_open(iunit, iflag, TRIM('regcoil_in.'// &
               TRIM(proc_string)), 'replace', 'formatted')
      ! write(6,'(a)') '<----write_regcoil_input'
      CALL write_regcoil_input(proc_string, iunit, istat)
      ! write(6,'(a)') '<----flush'
      CALL FLUSH(iunit)
      ! write(6,'(a)') '<----close'
      CLOSE(iunit)

      ! input file should be written. Now perform regcoil operation 
      ! This should be *almost* a duplicate of the main code from
      ! regcoil.f90
      ! JCS : Probably don't need to re-read the namelist, but I like 
      ! to do things in baby-steps with lots of debugging opptions and
      ! info
      !write(6,'(a)') '<----safe_open'
      !CALL safe_open(iunit, iflag, TRIM('input.'//TRIM(proc_string)), &
      !          'old', 'formatted')
      !write(6,'(a)') '<----read_regcoil_input'
      !call read_regcoil_input(iunit, iflag)
      ! write(6,'(a)') '<----Validate'
      call validate_input()
      ! write(6,'(a)') '<----Compute lambda'
      !if (allocated(lambda)) deallocate(lambda)
      call compute_lambda(lscreen)

      ! Define the position vector and normal vector at each grid point for
      ! the surfaces:
      ! write(6,'(a)') '<----init_plasma'
      call init_plasma(lscreen)
      ! write(6,'(a)') '<----init coil surfs'
      call init_coil_surface(lscreen)

      ! Initialize some of the vectors and matrices needed:
      ! write(6,'(a)') '<----read bnorm'
      call read_bnorm(lscreen)
      ! write(6,'(a)') '<----build matrices'
      call build_matrices(lscreen)

      ! JCS: I disabled all options except for #5 (for now)
      ! write(6,'(a)') '<----select a case'
      select case (general_option)
      !case (1)
      !   call solve()
      !case (2)
      !   call compute_diagnostics_for_nescout_potential()
      !case (3)
      !   call svd_scan()
      !case (4)
      !   call auto_regularization_solve()
      case (5)
      ! write(6,'(a)') '<----auto_reg solve'
         call auto_regularization_solve(lscreen)
      case default
         print *,"Invalid general_option:",general_option
         stop
      end select
    
!      write(6,'(a)') '<----safe_open'
!      CALL safe_open(iunit, iflag, TRIM('regcoil_in.'// &
!                TRIM(proc_string)), 'replace', 'formatted')
!      write(6,'(a)') '<----write_output'
!      call write_output()
!
!      write(6,'(a)') '<----flush'
!      CALL FLUSH(iunit)

      ! print *, chi2_B_target
      ! print *,"REGCOIL complete. Total time=",totalTime,"sec."


      IF (lscreen) WRITE(6,'(a)') ' ---------------------------  REGCOIL CALCULATION DONE  ---------------------'
!DEC$ ENDIF
      RETURN

!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_regcoil_chi2_b
