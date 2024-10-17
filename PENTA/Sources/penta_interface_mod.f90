!-----------------------------------------------------------------------
!     Module:        PENTA_INTERFACE_MOD
!     Authors:       S. Lazerson
!     Date:          10/16/2022
!     Description:   Module provides a subroutine based interface to
!                    the PENTA code.
!-----------------------------------------------------------------------
MODULE PENTA_INTERFACE_MOD
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
   USE penta_kind_mod

!-----------------------------------------------------------------------
!     Module Variables
!-----------------------------------------------------------------------
   IMPLICIT NONE

   INTEGER(iknd), PARAMETER :: NUM_ROOTS_MAX = 10_iknd
   INTEGER(iknd), PARAMETER :: NUM_ION_MAX = 20_iknd

   LOGICAL ::  input_is_Er, log_interp, use_quanc8, read_U2_file, &
      flux_cap, output_QoT_vs_Er, Add_Spitzer_to_D33, use_beam
   INTEGER(iknd) ::  num_Er_test, numKsteps, kord_pprof, keord, kcord, &
      numargs, js, i_append, num_species, num_ion_species, Smax, &
      iocheck, ie, ind_X, ind_A, ispec1, min_ind, iroot, num_roots
   REAL(rknd) ::  Kmin, Kmax, epsabs, epsrel, Er_min, Er_max, B_Eprl, &
      U2, vth_e, loglambda, Er_test, abs_Er, min_Er, eaEr_o_kTe, cmin, &
      cmax, emin, emax, sigma_par, sigma_par_Spitzer, J_BS
   REAL(rknd), DIMENSION(NUM_ION_MAX) ::  Z_ion_init, miomp_init
   REAL(rknd), DIMENSION(NUM_ROOTS_MAX) ::  Er_roots
   REAL(rknd), DIMENSION(:), ALLOCATABLE :: ion_mass, Z_ion, vth_i, &
      charges, dens, masses, temps, vths, dTdrs, dndrs, Xvec, Avec, &
      Flows, Gammas, QoTs, xt_c, xt_e, Gamma_e_vs_Er, QoT_e_vs_Er, &
      Er_test_vals, Jprl_ambi, sigma_par_ambi, sigma_par_Spitzer_ambi, &
      J_BS_ambi
   REAL(rknd), DIMENSION(:,:), ALLOCATABLE :: lmat, Dspl_D11, Dspl_D13, &
      Dspl_D31, Dspl_D33, Dspl_Dex, Dspl_Dua, Dspl_Drat, Dspl_Drat2, &
      Dspl_logD11, Dspl_logD33, cmesh, gamma_i_vs_er, QoT_i_vs_Er, &
      Flows_ambi, gammas_ambi, QoTs_ambi, Jprl_parts, upol, utor
   CHARACTER(LEN=10) :: Method
   CHARACTER(LEN=100) :: arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, &
      arg9, coeff_ext, run_ident, pprof_char, fpos, fstatus, str_num

!-----------------------------------------------------------------------
!     Module Namelists
!-----------------------------------------------------------------------
   NAMELIST /ion_params/ num_ion_species, Z_ion_init, miomp_init
   NAMELIST /run_params/ input_is_Er, log_interp, use_quanc8, &
      read_U2_file, Add_Spitzer_to_D33, num_Er_test, numKsteps, &
      kord_pprof, keord, kcord, Kmin, Kmax, epsabs, epsrel, Method, &
      flux_cap, output_QoT_vs_Er, use_beam

!-----------------------------------------------------------------------
!     SUBROUTINES
!-----------------------------------------------------------------------
   CONTAINS

   SUBROUTINE init_penta_input
      IMPLICIT NONE
      input_is_Er          = .TRUE.
      log_interp           = .TRUE.
      use_quanc8           = .FALSE.
      read_U2_file         = .TRUE.
      flux_cap             = .TRUE.
      output_QoT_vs_Er     = .FALSE.
      Add_Spitzer_to_D33   = .TRUE.
      use_beam             = .FALSE.
      num_Er_test          = 50_iknd
      numKsteps            = 10000_iknd
      kord_pprof           = 3_iknd
      keord                = 2_iknd
      kcord                = 2_iknd
      Kmin                 = 1.0E-5_rknd
      Kmax                 = 2.0E+1_rknd
      epsabs               = 1.0E-8_rknd
      epsrel               = 1.0E-6_rknd
      sigma_par            = 0.0E+0_rknd
      sigma_par_Spitzer    = 0.0E+0_rknd
      J_BS                 = 0.0E+0_rknd
      method               = 'DKES'
      RETURN
   END SUBROUTINE init_penta_input

   SUBROUTINE read_penta_ion_params_namelist(filename,istat)
      USE safe_open_mod
      IMPLICIT NONE
      CHARACTER(*), INTENT(in) :: filename
      INTEGER, INTENT(inout) :: istat
      LOGICAL :: lexist
      INTEGER :: iunit
      CHARACTER(LEN=1000) :: line
      istat = 0; iunit = 12
      INQUIRE(FILE=TRIM(filename),EXIST=lexist)
      IF (.not.lexist) stop 'Could not find input file'
      CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
      IF (istat /= 0) THEN
         WRITE(6,'(A)') 'ERROR opening file: ',TRIM(filename)
         CALL FLUSH(6)
         STOP
      END IF
      READ(iunit,NML=ion_params,IOSTAT=istat)
      IF (istat /= 0) THEN
         WRITE(6,'(A)') 'ERROR reading PENTA ion_params namelist from file: ',TRIM(filename)
         backspace(iunit)
         read(iunit,fmt='(A)') line
         write(6,'(A)') 'Invalid line in namelist: '//TRIM(line)
         CALL FLUSH(6)
         STOP
      END IF
      CLOSE(iunit)
      ! Electrons + ions
      num_species = num_ion_species + 1_iknd
      RETURN
   END SUBROUTINE read_penta_ion_params_namelist

   SUBROUTINE write_ion_params_nml(iunit)
      IMPLICIT NONE
      INTEGER(iknd), INTENT(inout) :: iunit
      CHARACTER(LEN=*), PARAMETER :: outint  = "(2X,A,1X,'=',1X,I0)"
      INTEGER(iknd) :: k
      WRITE(iunit,'(A)') '&ION_PARAMS'
      WRITE(iunit,'(A)') '----------------------------------------------------------------'
      WRITE(iunit,outint) 'NUM_ION_SPECIES',num_ion_species
      WRITE(iunit,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'Z_ION_INIT',(Z_ion_init(k), k=1,num_ion_species)
      WRITE(iunit,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'MIOMP_INIT',(miomp_init(k), k=1,num_ion_species)
      WRITE(iunit,'(A)') '/\n'
   END SUBROUTINE write_ion_params_nml

   SUBROUTINE write_ion_params_namelist_byfile(filename)
      USE safe_open_mod
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER(iknd) :: iunit, istat
      LOGICAL :: lexists
      
      iunit = 100
      istat = 0
      INQUIRE(FILE=TRIM(filename),exist=lexists)
      IF (lexists) THEN
         CALL safe_open(iunit,istat,TRIM(filename),'old','formatted', access_in='append')
         !OPEN(unit=iunit, file=TRIM(filename),iostat=istat, status="old", position="append")
      ELSE
         CALL safe_open(iunit,istat,TRIM(filename),'replace','formatted')
         !OPEN(unit=iunit, file=TRIM(filename),iostat=istat, status="new")
      END IF
      IF (istat .ne. 0) RETURN
      CALL write_ion_params_nml(iunit)
      CLOSE(iunit)

      RETURN
   END SUBROUTINE write_ion_params_namelist_byfile

   SUBROUTINE read_penta_run_params_namelist(filename,istat)
      USE safe_open_mod
      IMPLICIT NONE
      CHARACTER(*), INTENT(in) :: filename
      INTEGER, INTENT(inout) :: istat
      LOGICAL :: lexist
      INTEGER :: iunit
      CHARACTER(LEN=1000) :: line
      istat = 0; iunit = 12
      INQUIRE(FILE=TRIM(filename),EXIST=lexist)
      IF (.not.lexist) stop 'Could not find input file'
      CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
      IF (istat /= 0) THEN
         WRITE(6,'(A)') 'ERROR opening file: ',TRIM(filename)
         CALL FLUSH(6)
         STOP
      END IF
      READ(iunit,NML=run_params,IOSTAT=istat)
      IF (istat /= 0) THEN
         WRITE(6,'(A)') 'ERROR reading PENTA run_params namelist from file: ',TRIM(filename)
         backspace(iunit)
         read(iunit,fmt='(A)') line
         write(6,'(A)') 'Invalid line in namelist: '//TRIM(line)
         CALL FLUSH(6)
         STOP
      END IF
      CLOSE(iunit)
      RETURN
   END SUBROUTINE read_penta_run_params_namelist

   SUBROUTINE write_run_params_nml(iunit)
      IMPLICIT NONE
      INTEGER(iknd), INTENT(inout) :: iunit
      CHARACTER(LEN=*), PARAMETER :: outboo  = "(2X,A,1X,'=',1X,L1)"
      CHARACTER(LEN=*), PARAMETER :: outint  = "(2X,A,1X,'=',1X,I0)"
      CHARACTER(LEN=*), PARAMETER :: outdbl  = "(2X,A,1X,'=',1X,ES22.14)"
      CHARACTER(LEN=*), PARAMETER :: outstr  = "(2X,A,1X,'=',1X,'''',A,'''')"
      INTEGER(iknd) :: k
      WRITE(iunit,'(A)') '&RUN_PARAMS'
      WRITE(iunit,outboo) 'INPUT_IS_ER',input_is_er
      WRITE(iunit,outboo) 'LOG_INTERP',log_interp
      WRITE(iunit,outboo) 'USE_QUANC8',use_quanc8
      WRITE(iunit,outboo) 'READ_U2_FILE',read_U2_file
      WRITE(iunit,outboo) 'ADD_SPITZER_TO_D33',Add_Spitzer_to_D33
      WRITE(iunit,outboo) 'FLUX_CAP',flux_cap
      WRITE(iunit,outboo) 'OUTPUT_QOT_VS_ER',output_QoT_vs_Er
      WRITE(iunit,outboo) 'USE_BEAM',use_beam
      WRITE(iunit,outint) 'NUM_ER_TEST',num_Er_test
      WRITE(iunit,outint) 'NUMKSTEPS',numksteps
      WRITE(iunit,outint) 'KORD_PPROF',kord_pprof
      WRITE(iunit,outint) 'KEORD',keord
      WRITE(iunit,outint) 'KCORD',kcord
      WRITE(iunit,outdbl) 'KMIN',kmin
      WRITE(iunit,outdbl) 'KMAX',kmax
      WRITE(iunit,outdbl) 'EPSABS',epsabs
      WRITE(iunit,outdbl) 'EPSREL',epsrel
      WRITE(iunit,outstr) 'METHOD',method
      WRITE(iunit,'(A)') '/\n'
   END SUBROUTINE write_run_params_nml

   SUBROUTINE write_run_params_namelist_byfile(filename)
      USE safe_open_mod
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER(iknd) :: iunit, istat
      LOGICAL :: lexists
      
      iunit = 100
      istat = 0
      INQUIRE(FILE=TRIM(filename),exist=lexists)
      IF (lexists) THEN
         CALL safe_open(iunit,istat,TRIM(filename),'old','formatted', access_in='append')
         !OPEN(unit=iunit, file=TRIM(filename),iostat=istat, status="old", position="append")
      ELSE
         CALL safe_open(iunit,istat,TRIM(filename),'replace','formatted')
         !OPEN(unit=iunit, file=TRIM(filename),iostat=istat, status="new")
      END IF
      IF (istat .ne. 0) RETURN
      CALL write_run_params_nml(iunit)
      CLOSE(iunit)

      RETURN
   END SUBROUTINE write_run_params_namelist_byfile

   SUBROUTINE penta_init_commandline
      IMPLICIT NONE
      CALL GETCARG(1, arg1, numargs)
      IF (numargs .NE. 9) THEN
         Write(6,'(A)') 'ERROR Incorrect number of input arguments, see penta.f90 for details'
         CALL FLUSH(6)
         STOP
      END IF
      CALL GETCARG(2, arg2, numargs)
      CALL GETCARG(3, arg3, numargs)
      CALL GETCARG(4, arg4, numargs)
      CALL GETCARG(5, arg5, numargs)
      CALL GETCARG(6, arg6, numargs)
      CALL GETCARG(7, arg7, numargs)
      CALL GETCARG(8, arg8, numargs)
      CALL GETCARG(9, arg9, numargs)

      READ(arg2,*) Er_min
      READ(arg3,*) Er_max
      READ(arg4,*) js
      READ(arg5,*) i_append
      READ(arg8,*) B_Eprl
      READ(arg9,*) Smax
      coeff_ext  = TRIM(ADJUSTL(arg1))
      run_ident  = TRIM(ADJUSTL(arg6))
      pprof_char = TRIM(ADJUSTL(arg7))
      RETURN
   END SUBROUTINE penta_init_commandline

   SUBROUTINE penta_allocate_species
      USE pprof_pass, ONLY: ni, Ti, dnidr, dTidr
      IMPLICIT NONE
      ALLOCATE(ni(num_ion_species))
      ALLOCATE(Ti(num_ion_species))          ! Ion profile info
      ALLOCATE(dnidr(num_ion_species))
      ALLOCATE(dTidr(num_ion_species))
      ALLOCATE(Z_ion(num_ion_species))
      ALLOCATE(ion_mass(num_ion_species)) ! Ion parameters
      ALLOCATE(vth_i(num_ion_species))
      ALLOCATE(charges(num_species))                             ! Parameters for
      ALLOCATE(dens(num_species))                                !  all species
      ALLOCATE(vths(num_species)) 
      ALLOCATE(masses(num_species))
      ALLOCATE(Temps(num_species))
      ALLOCATE(dTdrs(num_species))
      ALLOCATE(dndrs(num_species))
      ALLOCATE(lmat((Smax+1)*num_species,(Smax+1)*num_species))  ! Clas. fric.coeffs
      ALLOCATE(Xvec(num_species*2+1))
      ALLOCATE(Avec(num_species*3))        ! Thermo. Force vecs.
      ALLOCATE(Flows((Smax+1)*num_species))                      ! Prl flow moments
      ALLOCATE(Gammas(num_species))                              ! Rad fluxes
      ALLOCATE(QoTs(num_species))                                ! Rad energy fluxes
      ALLOCATE(Gamma_i_vs_Er(num_Er_test,num_ion_species))       ! Ion flux vs Er
      ALLOCATE(Gamma_e_vs_Er(num_Er_test))                       ! Electron flux vs Er
      ALLOCATE(Er_test_vals(num_Er_test))                        ! Er to loop over
      IF ( output_QoT_vs_Er ) THEN
        ALLOCATE(QoT_i_vs_Er(num_Er_test,num_ion_species))       ! Ion flux vs Er
        ALLOCATE(QoT_e_vs_Er(num_Er_test))                       ! Electron flux vs Er
      ENDIF
      RETURN
   END SUBROUTINE penta_allocate_species

   SUBROUTINE penta_allocate_dkescoeff
      USE coeff_var_pass, ONLY: num_c, num_e
      IMPLICIT NONE
      ALLOCATE(xt_c(num_c + kcord))
      ALLOCATE(xt_e(num_e + keord))
      ALLOCATE(Dspl_D11(num_c,num_e))
      ALLOCATE(Dspl_D13(num_c,num_e))
      ALLOCATE(Dspl_D31(num_c,num_e))
      ALLOCATE(Dspl_D33(num_c,num_e))
      ALLOCATE(Dspl_Dex(num_c,num_e))
      ALLOCATE(Dspl_Drat(num_c,num_e)) 
      ALLOCATE(Dspl_Drat2(num_c,num_e)) 
      ALLOCATE(Dspl_DUa(num_c,num_e))  
      ALLOCATE(Dspl_logD11(num_c,num_e))
      ALLOCATE(Dspl_logD33(num_c,num_e))
      ALLOCATE(cmesh(num_c,num_e))
      RETURN
   END SUBROUTINE penta_allocate_dkescoeff

   SUBROUTINE penta_deallocate_species
      USE pprof_pass, ONLY: ni, Ti, dnidr, dTidr
      IMPLICIT NONE
      IF (ALLOCATED(ni)) DEALLOCATE(ni)
      IF (ALLOCATED(ti)) DEALLOCATE(ti)
      IF (ALLOCATED(dnidr)) DEALLOCATE(dnidr)
      IF (ALLOCATED(dTidr)) DEALLOCATE(dTidr)
      IF (ALLOCATED(Z_ion)) DEALLOCATE(Z_ion)
      IF (ALLOCATED(ion_mass)) DEALLOCATE(ion_mass)
      IF (ALLOCATED(vth_i)) DEALLOCATE(vth_i)
      IF (ALLOCATED(charges)) DEALLOCATE(charges)
      IF (ALLOCATED(dens)) DEALLOCATE(dens)
      IF (ALLOCATED(vths)) DEALLOCATE(vths)
      IF (ALLOCATED(masses)) DEALLOCATE(masses)
      IF (ALLOCATED(Temps)) DEALLOCATE(Temps)
      IF (ALLOCATED(dTdrs)) DEALLOCATE(dTdrs)
      IF (ALLOCATED(dndrs)) DEALLOCATE(dndrs)
      IF (ALLOCATED(lmat)) DEALLOCATE(lmat)
      IF (ALLOCATED(Xvec)) DEALLOCATE(Xvec)
      IF (ALLOCATED(Avec)) DEALLOCATE(Avec)
      IF (ALLOCATED(Flows)) DEALLOCATE(Flows)
      IF (ALLOCATED(Gammas)) DEALLOCATE(Gammas)
      IF (ALLOCATED(QoTs)) DEALLOCATE(QoTs)
      IF (ALLOCATED(Gamma_i_vs_Er)) DEALLOCATE(Gamma_i_vs_Er)
      IF (ALLOCATED(Gamma_e_vs_Er)) DEALLOCATE(Gamma_e_vs_Er)
      IF (ALLOCATED(Er_test_vals)) DEALLOCATE(Er_test_vals)
      IF (ALLOCATED(QoT_i_vs_Er)) DEALLOCATE(QoT_i_vs_Er)
      IF (ALLOCATED(QoT_e_vs_Er)) DEALLOCATE(QoT_e_vs_Er)
      RETURN
   END SUBROUTINE penta_deallocate_species

   SUBROUTINE penta_deallocate_dkescoeff
      IMPLICIT NONE
      IF (ALLOCATED(xt_c)) DEALLOCATE(xt_c)
      IF (ALLOCATED(xt_e)) DEALLOCATE(xt_e)
      IF (ALLOCATED(Dspl_D11)) DEALLOCATE(Dspl_D11)
      IF (ALLOCATED(Dspl_D13)) DEALLOCATE(Dspl_D13)
      IF (ALLOCATED(Dspl_D31)) DEALLOCATE(Dspl_D31)
      IF (ALLOCATED(Dspl_D33)) DEALLOCATE(Dspl_D33)
      IF (ALLOCATED(Dspl_Dex)) DEALLOCATE(Dspl_Dex)
      IF (ALLOCATED(Dspl_Drat)) DEALLOCATE(Dspl_Drat)
      IF (ALLOCATED(Dspl_Drat2)) DEALLOCATE(Dspl_Drat2)
      IF (ALLOCATED(Dspl_DUa)) DEALLOCATE(Dspl_DUa)
      IF (ALLOCATED(Dspl_logD11)) DEALLOCATE(Dspl_logD11)
      IF (ALLOCATED(Dspl_logD33)) DEALLOCATE(Dspl_logD33)
      IF (ALLOCATED(cmesh)) DEALLOCATE(cmesh)
      RETURN
   END SUBROUTINE penta_deallocate_dkescoeff

   SUBROUTINE penta_read_input_files
      USE vmec_var_pass, ONLY: roa_surf, arad, Bsq
      USE pprof_pass
      USE coeff_var_pass, ONLY: D11_mat, cmul, num_c
      USE phys_const, ONLY: p_mass, elem_charge, e_mass
      USE read_input_file_mod
      IMPLICIT NONE
      CALL read_vmec_file_2(js,run_ident)
      CALL read_pprof_file(pprof_char,num_ion_species,roa_surf,arad,kord_pprof)
      CALL read_dkes_star_files(coeff_ext,Add_Spitzer_to_D33,Bsq)
      IF (use_beam) THEN
         CALL read_beam_file(roa_surf,kord_pprof)
      ELSE
         beam_force = 0._rknd
      ENDIF
      ! Optionally read file containing <U**2> info.  Else this is 
      ! calculated from the D11* coefficient at high nu/v and Er=0.
      IF ( read_U2_file ) THEN
         CALL read_Utilde2_file(roa_surf,U2,kord_pprof)
      ELSE
         U2=1.5d0*D11_mat(num_c,1)/cmul(num_c);
      ENDIF
      ! Change Er test range to V/cm if necessary
      IF ( .not. input_is_Er) THEN
         Er_min = Er_min * Te / arad
         Er_max = Er_max * Te / arad
      ENDIF
      ! Assign ion parameters
      Z_ion    = Z_ion_init(1:num_ion_species)
      ion_mass = miomp_init(1:num_ion_species) * p_mass
      ! Calculate thermal velocities 
      vth_i = Dsqrt(2._rknd*Ti*elem_charge/ion_mass)
      vth_e = Dsqrt(2._rknd*Te*elem_charge/e_mass)
      ! Calculate Coulomb logarithm
      IF ( Te > 50._rknd ) THEN
        loglambda = 25.3_rknd - 1.15_rknd*Dlog10(ne/1.e6_rknd) + 2.3_rknd*Dlog10(Te)
      ELSE
        loglambda = 23.4_rknd - 1.15_rknd*Dlog10(ne/1.e6_rknd) + 3.45_rknd*Dlog10(Te)
      ENDIF
      ! Assign arrays of parameters for all species (charge, mass, n, T, v_th, dTdr)
      charges=elem_charge*(/-1._rknd,Z_ion/)
      dens=(/ne, ni/)
      masses=(/e_mass, ion_mass/)
      Temps=(/Te,Ti/)
      vths=(/vth_e,vth_i/)
      dTdrs=(/dTedr,dTidr/)
      dndrs=(/dnedr,dnidr/)
      RETURN
   END SUBROUTINE penta_read_input_files

   SUBROUTINE penta_fit_DXX_coef
      USE coeff_var_pass
      USE PENTA_subroutines, ONLY : fit_coeffs
      IMPLICIT NONE

      ! Calculate fitting parameters to the D##* coefficients
      Call fit_coeffs(cmul,efield,num_c,num_e,D11_mat,log_interp, &
         kcord,keord,xt_c,xt_e,Dspl_D11,cmin,cmax,emin,emax)
      Call fit_coeffs(cmul,efield,num_c,num_e,D13_mat,log_interp, &
         kcord,keord,xt_c,xt_e,Dspl_D13,cmin,cmax,emin,emax)
      Call fit_coeffs(cmul,efield,num_c,num_e,D31_mat,log_interp, &
         kcord,keord,xt_c,xt_e,Dspl_D31,cmin,cmax,emin,emax)
      Call fit_coeffs(cmul,efield,num_c,num_e,D33_mat,log_interp, &
         kcord,keord,xt_c,xt_e,Dspl_D33,cmin,cmax,emin,emax)

      ! Fit log(D*) for D11 and D33
      Call fit_coeffs(cmul,efield,num_c,num_e,LOG(D11_mat),log_interp, &
         kcord,keord,xt_c,xt_e,Dspl_logD11,cmin,cmax,emin,emax)
      Call fit_coeffs(cmul,efield,num_c,num_e,LOG(D33_mat),log_interp, &
         kcord,keord,xt_c,xt_e,Dspl_logD33,cmin,cmax,emin,emax)
      RETURN
   END SUBROUTINE penta_fit_DXX_coef

   SUBROUTINE penta_screen_info
      IMPLICIT NONE
      If ( i_append == 0 ) Then
         WRITE(6,'(A)') ""
         WRITE(6,'(A)') "Welcome to PENTA3, please note the following settings:"
         WRITE(6,'(A)')
         WRITE(6,'(A,I3)') ' Number of ion species: ',num_ion_species
         IF ( input_is_Er ) THEN
            WRITE(6,'(A)') 'Interpreting input range as Er (V/cm)'
         ELSE
            WRITE(6,'(A)') 'Interpreting input range as e<a>Er/kTe'
         ENDIF
         IF ( log_interp ) THEN
            WRITE(6,'(A)') 'Performing logarithmic interpolation in Er, cmul'
         ELSE
            WRITE(6,'(A)') 'Performing linear interpolation in Er,cmul'
         ENDIF
         IF ( use_quanc8 ) THEN
            WRITE(6,'(A,2(A,E10.4))')                           &
               ' Using quanc8 integrator with tolerances: ',     &
               'abs: ',epsabs,' rel: ', epsrel
         ELSE
            WRITE(6,'(A,i6,A)') ' Using ',numKsteps,            &
               ' point integral approximation'
         ENDIF
         WRITE(6,'(a,2(" ",e15.4))') ' K range on convolution integral: ', &
            Kmin, Kmax
         IF ( Add_Spitzer_to_D33) &
            WRITE(6,'(A)') 'Adding collisional (Spitzer) portion to D33 coefficient'
         IF ( flux_cap ) &
            WRITE(6,'(A)') 'Enforcing minimum radial diffusion coefficient = 0'
         WRITE(6,'(A,I2)') ' Number of terms in Sonine expansion: ', Smax+1
         SELECT CASE (Method)
            CASE ('T')
               WRITE(6,'(A)') 'Using Taguchi Method'
            CASE ('SN')
               WRITE(6,'(A)') 'Using Sugama-Nishimura Method'
            CASE ('MBT')
               WRITE(6,'(A)') 'Using Maassberg-Beidler-Turkin Method'
            CASE ('DKES')
               WRITE(6,'(A)') 'Using DKES Method'
            CASE DEFAULT
               Write(6,'(3A)') ' Error: ''', Trim(Adjustl(Method)), &
                  ''' is not a valid Method'
               STOP 'Error: Exiting, method select error in penta.f90 (1)'
         END SELECT  
         WRITE(6,'(A)') ""
         WRITE(6,'(A)') " <r>/<a>","   Er root(s) (V/cm)"
      ENDIF
      RETURN
   END SUBROUTINE penta_screen_info

   SUBROUTINE penta_open_output
      USE safe_open_mod
      USE io_unit_spec
      IMPLICIT NONE
      INTEGER :: istat
      istat = 0
      ! Set write status
      IF (i_append == 0) THEN
         fstatus = "unknown"
         fpos = "SEQUENTIAL"
      ELSEIF (i_append == 1) THEN
         fstatus = "old"
         fpos = "append"
      ELSE
         Write(6,'(A)') 'PENTA: Bad value for i_append (0 or 1 expected)'
         CALL FLUSH(6)
         STOP 'Error: Exiting, i_append error in penta.f90'
      END IF
      ! Open files
      !Open(unit=iu_flux_out, file="fluxes_vs_roa", position=Trim(Adjustl(fpos)),status=Trim(Adjustl(fstatus)))
      CALL safe_open(iu_flux_out, istat, "fluxes_vs_roa", &
         Trim(Adjustl(fstatus)), 'formatted',&
         access_in=Trim(Adjustl(fpos)))
      CALL safe_open(iu_pprof_out, istat, "plasma_profiles_check", &
         Trim(Adjustl(fstatus)), 'formatted',&
         access_in=Trim(Adjustl(fpos)))
      CALL safe_open(iu_fvEr_out, istat, "fluxes_vs_Er", &
         Trim(Adjustl(fstatus)), 'formatted',&
         access_in=Trim(Adjustl(fpos)))
      CALL safe_open(iu_flows_out, istat, "flows_vs_roa", &
         Trim(Adjustl(fstatus)), 'formatted',&
         access_in=Trim(Adjustl(fpos)))
      CALL safe_open(iu_flowvEr_out, istat, "flows_vs_Er", &
         Trim(Adjustl(fstatus)), 'formatted',&
         access_in=Trim(Adjustl(fpos)))
      CALL safe_open(iu_Jprl_out, istat, "Jprl_vs_roa", &
         Trim(Adjustl(fstatus)), 'formatted',&
         access_in=Trim(Adjustl(fpos)))
      CALL safe_open(iu_contraflows_out, istat, "ucontra_vs_roa", &
         Trim(Adjustl(fstatus)), 'formatted',&
         access_in=Trim(Adjustl(fpos)))
      IF (method == 'SN') &
         CALL safe_open(iu_sigmas_out, istat, "sigmas_vs_roa", &
            Trim(Adjustl(fstatus)), 'formatted',&
            access_in=Trim(Adjustl(fpos)))
      IF (output_QoT_vs_Er) THEN
         CALL safe_open(iu_QoTvEr_out, istat, "QoTs_vs_Er", &
            Trim(Adjustl(fstatus)), 'formatted',&
            access_in=Trim(Adjustl(fpos)))
         Write(iu_QoTvEr_out,'("*",/,"r/a   Er[V/cm]   Q_e/T_e [m**-2s**-1] ",&
                               "   Q_i/T_i [m**-2s**-1]")')
      END IF
      ! WRITE Headers
      IF (i_append == 0) THEN
         ! Fluxes vs r/a
         Write(iu_flux_out,'("*",/,"r/a    Er[V/cm]    e<a>Er/kTe    ",  &
            "Gamma_e [m**-2s**-1]   Q_e/T_e [m**-2s**-1]     ",         &
            "Gamma_i [m**-2s**-1]   Q_i/T_i [m**-2s**-1]")')
         ! Flows vs r/a
         Write(iu_flows_out,'("*",/,"r/a   Er[V/cm]    e<a>Er/kTe    ",  &
            " <B*u_||ke>/<B**2> [m/sT]   <B*u_||ki>/<B**2> [m/sT]")')
         ! Plasma profile check
         Write(iu_pprof_out,'("*",/,"r/a    Te [eV]   ne [m**-3]     ",  & 
            "dnedr [m**-4]   dTedr [eV/m]  Ti [eV]   ni [m**-3]     ",  &
            "dnidr [m**-4]   dTidr [eV/m]")')
         Write(iu_Jprl_out,'("*",/,"r/a    Er [V/cm]    e<a>Er/kTe    ",  &
            "Jprl_e [A/m**2]    Jprli [A/m**2]    Jprl [A/m**2]    J_BS [A/m**2]")')
         Write(iu_contraflows_out,'("*",/,"r/a    Er [V/cm]    e<a>Er/kTe    ",  &
            "<ue^pol_contra> [1/s]     <ue^tor_contra> [1/s]       ",  &
            " <ui^pol_contra> [1/s]     <ui^tor_contra> [1/s]")')
         ! Sigmas vs r/a
         IF (Method == 'SN') &
            Write(iu_sigmas_out,'("*",/,"r/a   Er[V/cm]    sigma_par [1/Ohm.m]    ",  &
               " sigma_par_Spitzer [1/Ohm.m]")')
         ! Legend for fluxes vs Er 
         Write(iu_fvEr_out,'("*",/,"r/a   Er[V/cm]   Gamma_e [m**-2s**-1] ",&
            "   Gamma_i [m**-2s**-1]")')
         ! Legend for flows vs Er
         Write(iu_flowvEr_out,'("*",/,"r/a   Er[V/cm]  ", &
            "    <B*u_||ke>/<B**2> [m/sT]  <B*u_||ki>/<B**2> [m/sT]")')
      END IF
      RETURN
   END SUBROUTINE penta_open_output

   SUBROUTINE penta_fit_rad_trans
      USE coeff_var_pass
      USE vmec_var_pass
      USE PENTA_subroutines, ONLY : fit_coeffs, define_friction_coeffs
      IMPLICIT NONE
      ! Define matrix of friction coefficients (lmat)
      Call define_friction_coeffs(masses,charges,vths,Temps,dens,loglambda, &
                            num_species,Smax,lmat)
      ! Fit radial transport coefficients specific to different methods
      SELECT CASE (Method)
         CASE ('T', 'MBT')
            cmesh = Spread(cmul,2,num_e)
            ! Calculate the D11 coefficient minus the P-S contribution  (Dex)
            ! Also, do not allow for negative coefficients
            CALL fit_coeffs(cmul,efield,num_c,num_e, &
               Max(D11_mat-(2._rknd/3._rknd)*cmesh*U2,0._rknd), &
               log_interp,kcord,keord,xt_c,xt_e,Dspl_Dex,     &
               cmin,cmax,emin,emax)
         CASE ('SN')
            ! Calculate fits to D31*/D33*  (Drat)
            CALL fit_coeffs(cmul,efield,num_c,num_e, &
               D31_mat/D33_mat, &
               log_interp,kcord,keord,xt_c,xt_e,Dspl_Drat,     &
               cmin,cmax,emin,emax)
            ! Calculate fits to (D31*)**2/D33*   (Drat2)
            CALL fit_coeffs(cmul,efield,num_c,num_e, &
               D31_mat*D31_mat/D33_mat, &
               log_interp,kcord,keord,xt_c,xt_e,Dspl_Drat2,     &
               cmin,cmax,emin,emax)
            cmesh = Spread(cmul,2,num_e)
            ! Calculate coefficient for Ua term  (DUa)
            CALL fit_coeffs(cmul,efield,num_c,num_e, &
               (2._rknd*Bsq/(3._rknd*D33_mat) - cmesh), &
               log_interp,kcord,keord,xt_c,xt_e,Dspl_DUa,     &
               cmin,cmax,emin,emax)
            ! Calculate coefficient for radial flux  (Capped term)  (Dex)
            ! Also, do not allow for negative coefficients
            CALL fit_coeffs(cmul,efield,num_c,num_e, &
               Max(D11_mat-(2._rknd/3._rknd)*cmesh*U2+D31_mat*D31_mat/D33_mat, &
               0._rknd),log_interp,kcord,keord,xt_c,xt_e,Dspl_Dex,     &
               cmin,cmax,emin,emax)
        CASE ('DKES')
        CASE DEFAULT
          WRITE(6,'(3a)') ' Error: ''', Trim(Adjustl(Method)), &
            ''' is not a valid Method'
          STOP 'Error: Exiting, method select error in penta.f90 (2)'
      ENDSELECT
      RETURN
   END SUBROUTINE penta_fit_rad_trans

   SUBROUTINE penta_run_1_init
      IMPLICIT NONE
      INTEGER :: istat
      CALL init_penta_input
      istat = 0
      CALL read_penta_ion_params_namelist('ion_params',istat)
      istat = 0
      CALL read_penta_run_params_namelist('run_params',istat)
      CALL penta_init_commandline
      CALL penta_allocate_species
      CALL penta_read_input_files
      CALL penta_screen_info
      CALL penta_allocate_dkescoeff
      CALL penta_fit_DXX_coef
      CALL penta_open_output
      CALL penta_fit_rad_trans
      RETURN
   END SUBROUTINE penta_run_1_init

   SUBROUTINE penta_run_2_efield
      USE vmec_var_pass
      USE pprof_pass
      USE phys_const
      USE io_unit_spec
      USE coeff_var_pass
      USE penta_math_routines_mod, ONLY: rlinspace
      USE penta_functions_mod
      USE PENTA_subroutines, ONLY: form_xvec
      IMPLICIT NONE
      ! Define array of Er values to test [V/m]
      Er_test_vals = rlinspace(Er_min,Er_max,num_Er_test)*100._rknd

      ! Check for Er=0, doesn't work for log interpolation
      min_Er = Minval(Dabs(Er_test_vals),DIM=1) 
      If ((log_interp .EQV. .true. ) .AND. ( Dabs(min_Er) <= elem_charge ))  Then
         min_ind = Minloc(Dabs(Er_test_vals),DIM=1)
         If ( min_ind == Num_Er_test ) Then 
            Er_test_vals(min_ind) = Er_test_vals(min_ind - 1)/2._rknd
         Else
            Er_test_vals(min_ind) = Er_test_vals(min_ind + 1)/2._rknd
         EndIf
         Write(*,'(a,i4,a,f10.3)') 'Cannot use Er=0 with log_interp, using Er(',  &
            min_ind, ') = ', Er_test_vals(min_ind)
      EndIf
      ! Loop over Er to get fluxes as a function of Er
      Do ie = 1,num_Er_test
         Er_test = Er_test_vals(ie)
         abs_Er = Abs(Er_test)

         ! Form thermodynamic force vector (Xvec)
         Call form_Xvec(Er_test,Z_ion,B_Eprl,num_ion_species,Xvec)

         ! Form alternate thermodynamic force vector (Avec)
         Do ispec1 = 1,num_species
            ind_X = (ispec1-1)*2 + 1
            ind_A = (ispec1-1)*3 + 1

            Avec(ind_A)   = -Xvec(ind_X) / (Temps(ispec1)*elem_charge) &
               - 2.5_rknd*dTdrs(ispec1)/Temps(ispec1)
            Avec(ind_A+1)   = -Xvec(ind_X+1) / (Temps(ispec1)*elem_charge)
            Avec(ind_A+2)   = Xvec(num_species*2+1)*charges(ispec1) &
               * B0/(Temps(ispec1)*elem_charge*Sqrt(Bsq)) + &
               beam_force/(Temps(ispec1)*elem_charge*dens(ispec1))
         Enddo

         ! Select the appropriate algorithm and calculate the flows and fluxes
         SELECT CASE (Method)
            Case ('T', 'MBT')
               ! Calculate array of parallel flow moments
               Flows = calc_flows_T(num_species,Smax,abs_Er,Temps,dens,vths,charges,   &
                 masses,loglambda,B0,use_quanc8,Kmin,Kmax,numKsteps,log_interp,cmin,   &
                 cmax,emin,emax,xt_c,xt_e,Dspl_D31,Dspl_logD33,num_c,num_e,kcord,      &
                 keord,Avec,Bsq,lmat,J_BS)
               Gammas = calc_fluxes_MBT(num_species,Smax,abs_Er,Temps,dens,vths,       &
                 charges,masses,dTdrs,dndrs,loglambda,use_quanc8,Kmin,Kmax,numKsteps,  &
                 log_interp,cmin,cmax,emin,emax,xt_c,xt_e,Dspl_logD11,Dspl_D31,        &
                 Dspl_Dex,num_c,num_e,kcord,keord,Avec,lmat,Flows,U2,B0,flux_cap)   
               If ( output_QoT_vs_Er .EQV. .true. ) Then
                  QoTs = calc_QoTs_MBT(num_species,Smax,abs_Er,Temps,dens,vths,charges, &
                   masses,dTdrs,dndrs,loglambda,use_quanc8,Kmin,Kmax,numKsteps,        &
                   log_interp,cmin,cmax,emin,emax,xt_c,xt_e,Dspl_logD11,Dspl_D31,      &
                   Dspl_Dex,num_c,num_e,kcord,keord,Avec,lmat,Flows,U2,B0,flux_cap)   
               Endif    
            Case ('SN')                    
               Flows = calc_flows_SN(num_species,Smax,abs_Er,Temps,dens,vths,charges,  &
                  masses,loglambda,B0,use_quanc8,Kmin,Kmax,numKsteps,log_interp,       &
                  cmin,cmax,emin,emax,xt_c,xt_e,Dspl_Drat,Dspl_DUa,num_c,num_e,kcord,  &
                  keord,Avec,lmat,sigma_par,sigma_par_Spitzer,J_BS)                                                
               Gammas = calc_fluxes_SN(num_species,Smax,abs_Er,Temps,dens,vths,charges,&
                 masses,loglambda,use_quanc8,Kmin,Kmax,numKsteps,log_interp,cmin,cmax, &
                 emin,emax,xt_c,xt_e,Dspl_Drat,Dspl_Drat2,Dspl_Dex,Dspl_logD11,        &
                 Dspl_D31,num_c,num_e,kcord,keord,Avec,Bsq,lmat,Flows,U2,dTdrs,        &
                 dndrs,flux_cap)  
               If ( output_QoT_vs_Er .EQV. .true. ) Then
                  QoTs = calc_QoTs_SN(num_species,Smax,abs_Er,Temps,dens,vths,charges,  &
                     masses,loglambda,use_quanc8,Kmin,Kmax,numKsteps,log_interp,cmin,    &
                     cmax,emin,emax,xt_c,xt_e,Dspl_Drat,Dspl_Drat2,Dspl_Dex,Dspl_logD11, &
                     Dspl_D31,num_c,num_e,kcord,keord,Avec,Bsq,lmat,Flows,U2,dTdrs,      &
                     dndrs,flux_cap)  
               Endif    
            Case ('DKES')
               Flows = calc_flows_DKES(num_species,Smax,abs_Er,Temps,dens,vths,charges,&
                  masses,loglambda,B0,use_quanc8,Kmin,Kmax,numKsteps,log_interp,cmin,  &
                  cmax,emin,emax,xt_c,xt_e,Dspl_D31,Dspl_logD33,num_c,num_e,kcord,     &
                  keord,Avec,J_BS)
               Gammas = calc_fluxes_DKES(num_species,abs_Er,Temps,dens,vths,charges,   &
                  masses,loglambda,use_quanc8,Kmin,Kmax,numKsteps,log_interp,cmin,cmax, &
                  emin,emax,xt_c,xt_e,Dspl_logD11,Dspl_D31,num_c,num_e,kcord,keord,     &
                  Avec,B0)   
               If ( output_QoT_vs_Er .EQV. .true. ) Then
                  QoTs = calc_QoTs_DKES(num_species,abs_Er,Temps,dens,vths,charges,     &
                     masses,loglambda,use_quanc8,Kmin,Kmax,numKsteps,log_interp,cmin,    &
                     cmax,emin,emax,xt_c,xt_e,Dspl_logD11,Dspl_D31,num_c,num_e,kcord,    &
                     keord,Avec,B0)  
               Endif
            Case Default
               Write(6,'(3a)') ' Error: ''', Trim(Adjustl(Method)), &
              ''' is not a valid Method'
               Stop 'Error: Exiting, method select error in penta.f90 (3)'
         END SELECT

         Gamma_e_vs_Er(ie)   = Gammas(1)
         Gamma_i_vs_Er(ie,:) = Gammas(2:num_species)

         ! Write fluxes vs Er
         Write(str_num,*) num_ion_species + 2  ! Convert num to string
         Write(iu_fvEr_out,'(f7.4,' // trim(adjustl(str_num)) // '(" ",e15.7))') &
            roa_surf,Er_test/100._rknd,Gamma_e_vs_Er(ie),Gamma_i_vs_Er(ie,:)

         If ( output_QoT_vs_Er .EQV. .true. ) Then
            QoT_e_vs_Er(ie)   = QoTs(1)
            QoT_i_vs_Er(ie,:) = QoTs(2:num_species)
            Write(iu_QoTvEr_out,'(f7.4,' // trim(adjustl(str_num)) // '(" ",e15.7))') &
               roa_surf,Er_test/100._rknd,QoT_e_vs_Er(ie),QoT_i_vs_Er(ie,:)
         Endif

         ! Write flows vs Er
         Write(str_num,*) (Smax+1)*num_species + 2  ! Convert num to string
         Write(iu_flowvEr_out,'(f7.4,' // trim(adjustl(str_num)) // '(" ",e15.7))') &
          roa_surf,Er_test/100._rknd,Flows

      Enddo !efield loop
      RETURN
   END SUBROUTINE penta_run_2_efield

   SUBROUTINE penta_run_3_ambipolar
      USE vmec_var_pass
      USE phys_const
      USE coeff_var_pass
      USE penta_functions_mod
      USE PENTA_subroutines, ONLY: form_xvec, find_Er_roots
      IMPLICIT NONE
      ! Check for only one Er test value -- this is then used to evaluate the ambipolar fluxes QQ
      !If ( num_Er_test  == 1 ) Then
      !  Er_roots = Er_test_vals

      ! Find the ambipolar root(s) from gamma_e = sum(Z*gamma_i)
      Call find_Er_roots(gamma_e_vs_Er,gamma_i_vs_Er,Er_test_vals,Z_ion, &
         num_Er_test,num_ion_species,Er_roots,num_roots)

      ! Allocate arrays according to number of ambipolar roots
      Allocate(Flows_ambi((Smax+1)*num_species,num_roots)) ! Parallel flow moments
      Allocate(Gammas_ambi(num_species,num_roots))         ! Rad. particle fluxes
      Allocate(QoTs_ambi(num_species,num_roots))           ! Rad. energy fluxes
      Allocate(Jprl_ambi(num_roots))                       ! Parallel current density
      Allocate(J_BS_ambi(num_roots))                       ! BS current density
      Allocate(sigma_par_ambi(num_roots))                  ! Parallel conductivity
      Allocate(sigma_par_Spitzer_ambi(num_roots))          ! Spitzer Parallel conductivity
      Allocate(Jprl_parts(num_species,num_roots))          ! Par. curr. dens. per spec.
      Allocate(upol(num_species,num_roots))                ! fsa contra pol flow
      Allocate(utor(num_species,num_roots))                ! fsa contra tor flow

      ! Evaluate fluxes and flows at the ambipolar Er
      Do iroot = 1_iknd, num_roots

         Er_test = Er_roots(iroot)
         abs_Er = Dabs(Er_test)

         ! Form thermodynamic force vector (Xvec)
         Call form_Xvec(Er_test,Z_ion,B_Eprl,num_ion_species,Xvec)

         ! Form alternate thermodynamic force vector (Avec)
         Do ispec1 = 1,num_species
            ind_X = (ispec1-1)*2 + 1
            ind_A = (ispec1-1)*3 + 1

            Avec(ind_A)   = -Xvec(ind_X) / (Temps(ispec1)*elem_charge) &
               - 2.5_rknd*dTdrs(ispec1)/Temps(ispec1)
            Avec(ind_A+1)   = -Xvec(ind_X+1) / (Temps(ispec1)*elem_charge)
            Avec(ind_A+2)   = Xvec(num_species*2+1)*charges(ispec1) &
               * B0/(Temps(ispec1)*elem_charge*Sqrt(Bsq))
         Enddo

         ! Select the appropriate algorithm and calculate the flows and fluxes
         SELECT CASE (Method)
            Case ('T', 'MBT')
               ! Calculate array of parallel flow moments 
                 ! Note: Flow methods are the same for T and MBT
               Flows_ambi(:,iroot) = calc_flows_T(num_species,Smax,abs_Er,Temps,dens, &
                 vths,charges,masses,loglambda,B0,use_quanc8,Kmin,Kmax,numKsteps,     &
                 log_interp,cmin,cmax,emin,emax,xt_c,xt_e,Dspl_D31,Dspl_logD33,num_c, &
                 num_e,kcord,keord,Avec,Bsq,lmat,J_BS)
               ! Calculate array of radial particle fluxes
               Gammas_ambi(:,iroot) = calc_fluxes_MBT(num_species,Smax,abs_Er,Temps,  &
                 dens,vths,charges,masses,dTdrs,dndrs,loglambda,use_quanc8,Kmin,Kmax, &
                 numKsteps,log_interp,cmin,cmax,emin,emax,xt_c,xt_e,Dspl_logD11,      &
                 Dspl_D31,Dspl_Dex,num_c,num_e,kcord,keord,Avec,lmat,                 &
                 Flows_ambi(:,iroot),U2,B0,flux_cap)   
               ! Calculate array of radial energy fluxes
               QoTs_ambi(:,iroot) = calc_QoTs_MBT(num_species,Smax,abs_Er,Temps,dens, &
                 vths,charges,masses,dTdrs,dndrs,loglambda,use_quanc8,Kmin,Kmax,      &
                 numKsteps,log_interp,cmin,cmax,emin,emax,xt_c,xt_e,Dspl_logD11,      &
                 Dspl_D31,Dspl_Dex,num_c,num_e,kcord,keord,Avec,lmat,                 &
                 Flows_ambi(:,iroot),U2,B0,flux_cap)

               J_BS_ambi(iroot) = J_BS
            Case ('SN')
               ! Calculate array of parallel flow moments
                                                       
               Flows_ambi(:,iroot) = calc_flows_SN(num_species,Smax,abs_Er,Temps,dens,&
                  vths,charges,masses,loglambda,B0,use_quanc8,Kmin,Kmax,numKsteps,    &
                  log_interp,cmin,cmax,emin,emax,xt_c,xt_e,Dspl_Drat,Dspl_DUa,num_c,  &
                  num_e,kcord,keord,Avec,lmat,sigma_par,sigma_par_Spitzer,J_BS)                                                
               ! Calculate array of radial particle fluxes
               Gammas_ambi(:,iroot) = calc_fluxes_SN(num_species,Smax,abs_Er,Temps,   &
                 dens,vths,charges,masses,loglambda,use_quanc8,Kmin,Kmax,numKsteps,   &
                 log_interp,cmin,cmax,emin,emax,xt_c,xt_e,Dspl_Drat,Dspl_Drat2,       &
                 Dspl_Dex,Dspl_logD11,Dspl_D31,num_c,num_e,kcord,keord,Avec,Bsq,      &
                 lmat,Flows_ambi(:,iroot),U2,dTdrs,dndrs,flux_cap)  
               ! Calculate array of radial energy fluxes
               QoTs_ambi(:,iroot) = calc_QoTs_SN(num_species,Smax,abs_Er,Temps,dens,  &
                 vths,charges,masses,loglambda,use_quanc8,Kmin,Kmax,numKsteps,        &
                 log_interp,cmin,cmax,emin,emax,xt_c,xt_e,Dspl_Drat,Dspl_Drat2,       &
                 Dspl_Dex,Dspl_logD11,Dspl_D31,num_c,num_e,kcord,keord,Avec,Bsq,      &
                 lmat,Flows_ambi(:,iroot),U2,dTdrs,dndrs,flux_cap)

               sigma_par_ambi(iroot) = sigma_par
               sigma_par_Spitzer_ambi(iroot) = sigma_par_Spitzer
               J_BS_ambi(iroot) = J_BS
            Case ('DKES')
               ! Calculate array of parallel flow moments 
               Flows_ambi(:,iroot) = calc_flows_DKES(num_species,Smax,abs_Er,Temps,   &
                 dens,vths,charges,masses,loglambda,B0,use_quanc8,Kmin,Kmax,numKsteps,&
                 log_interp,cmin,cmax,emin,emax,xt_c,xt_e,Dspl_D31,Dspl_logD33,num_c, &
                 num_e,kcord,keord,Avec,J_BS)
               ! Calculate array of radial particle fluxes
               Gammas_ambi(:,iroot) = calc_fluxes_DKES(num_species,abs_Er,Temps,dens, &
                 vths,charges,masses,loglambda,use_quanc8,Kmin,Kmax,numKsteps,        &
                 log_interp,cmin,cmax,emin,emax,xt_c,xt_e,Dspl_logD11,Dspl_D31,num_c, &
                 num_e,kcord,keord,Avec,B0)  
               ! Calculate array of radial energy fluxes
               QoTs_ambi(:,iroot) = calc_QoTs_DKES(num_species,abs_Er,Temps,dens,     &
                 vths,charges,masses,loglambda,use_quanc8,Kmin,Kmax,numKsteps,        &
                 log_interp,cmin,cmax,emin,emax,xt_c,xt_e,Dspl_logD11,Dspl_D31,num_c, &
                 num_e,kcord,keord,Avec,B0)
               
               J_BS_ambi(iroot) = J_BS                
            Case Default
               Write(*,'(3a)') ' Error: ''', Trim(Adjustl(Method)), &
                  ''' is not a valid Method'
               Stop 'Error: Exiting, method select error in penta.f90 (4)'
         ENDSELECT

         ! Calculate parallel current density
         Jprl_parts(:,iroot) = dens*charges*Sqrt(Bsq) *         &
            Flows_ambi(1:(num_species-1)*(Smax+1)+1:Smax+1,iroot)
         Jprl_ambi(iroot) = Sum(Jprl_parts(:,iroot))

         ! Calculate flow components
         upol(:,iroot) = chip*(4._rknd*pi*pi/vol_p)*(                 &
            Flows_ambi(1:(num_species-1)*(Smax+1)+1:Smax+1,iroot) - &
            bzeta*Xvec(1:(num_species-1)*2+1:2)/(charges*chip*Bsq))
         utor(:,iroot) = psip*(4._rknd*pi*pi/vol_p)*(                 &
            Flows_ambi(1:(num_species-1)*(Smax+1)+1:Smax+1,iroot) + &
            btheta*Xvec(1:(num_species-1)*2+1:2)/(charges*psip*Bsq))

      END DO ! Ambipolar root loop
      RETURN
   END SUBROUTINE penta_run_3_ambipolar

   SUBROUTINE penta_run_4_cleanup
      USE io_unit_spec
      USE pprof_pass
      USE vmec_var_pass
      IMPLICIT NONE
      ! First write output files
      ! Loop over ambipolar Er for writing output files
      Do iroot = 1_iknd, num_roots

         Er_test = Er_roots(iroot)
         eaEr_o_kTe = arad*Er_test/Te

         ! Write fluxes to file "fluxes_vs_roa"
         Write(str_num,*) 2*num_species + 2
         Write(iu_flux_out,'(f7.3,' // Trim(Adjustl(str_num)) // '(" ",e15.7))') &
          roa_surf,Er_test/100._rknd,eaEr_o_kTe,Gammas_ambi(1,iroot),  &
          QoTs_ambi(1,iroot),Gammas_ambi(2:num_species,iroot),  &
          QoTs_ambi(2:num_species,iroot)

         ! Write flows to file "flows_vs_roa"
         Write(str_num,*) (Smax+1)*num_species + 2
         Write(iu_flows_out,'(f7.3,' // trim(adjustl(str_num)) // '(" ",e15.7))') &
          roa_surf,Er_test/100._rknd,eaEr_o_kTe,Flows_ambi(:,iroot)

         ! Write current densities to file "Jprl_vs_roa"
         Write(str_num,*) num_species + 4 
         Write(iu_Jprl_out,'(f7.3,' // trim(adjustl(str_num)) // '(" ",e15.7))')  & 
          roa_surf,Er_test/100._rknd,eaEr_o_kTe,Jprl_parts(:,iroot),Jprl_ambi(iroot),J_BS_ambi(iroot)

         ! Write contravariant flows to file "ucontra_vs_roa"
         Write(str_num,*) 2*num_species + 2
         Write(iu_contraflows_out,'(f7.3,' // trim(adjustl(str_num))//'(" ",e15.7))') & 
          roa_surf,Er_test/100._rknd,eaEr_o_kTe,upol(1,iroot),utor(1,iroot),         &
          upol(2:num_species,iroot),utor(2:num_species,iroot)

         ! Write sigmas to file "sigmas_vs_roa"
         If( Method == 'SN') then
          Write(str_num,*) 3
          Write(iu_sigmas_out,'(f7.3,' // trim(adjustl(str_num)) // '(" ",e15.7))') &
            roa_surf,Er_test/100._rknd,sigma_par_ambi(iroot),sigma_par_Spitzer_ambi(iroot)
         Endif
      EndDo ! Ambipolar root loop

      ! Write plasma profile information to "plasma_profiles_check"
      Write(str_num,*) 4*num_species
      Write(iu_pprof_out,'(f7.3,' // trim(adjustl(str_num)) // '(" ",e15.7))') & 
        roa_surf,Te,ne,dnedr,dTedr,Ti,ni,dnidr,dTidr

      ! QQ write file with number of roots per surface!

      ! Write screen output
      write(str_num,*) num_roots
      write(*,'(f7.3,' // trim(adjustl(str_num)) // '(" ",e15.4))') & 
        roa_surf,er_roots(1:num_roots)/100._rknd
      ! DEALLOCATE
      CALL penta_deallocate_species
      CALL penta_deallocate_dkescoeff
   END SUBROUTINE penta_run_4_cleanup



END MODULE PENTA_INTERFACE_MOD