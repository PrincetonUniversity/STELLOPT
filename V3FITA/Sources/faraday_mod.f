!     SPH: INTEGER(iprec) -> INTEGER
!============================================================================
      MODULE faraday_mod
!--------------------------------------------------------------------------
! 
! FUNCTION: Module containing routines necessary to compute the line-integrated
!           Faraday rotation angle and the line-integrated polarization phase shift.
!              Note that the code assumes that there is a small but finite e- density outside the last closed
!           flux surface.  Hence, for those cases, the external B field is used to calculate the (presumably very
!           small) contribution to net line integral (ie the code does NOT automatically assume a priori that all points
!           outside the flux surface have an exactly zero contribution).
!
! created May 2005 by J. Shields
!
! Required input files:  ip_v3rfun_*.nc (a netCDF file produced by V3RFUN containing info about the number and the
!                                          positions of the microwave beams)
!
! optional output file: frcm.out  (text file with table of Faraday rotation & Cotton-Mouton values)
!                        It is possible to inhibit the output of this file by setting the module-global
!                        logical toggle OUTPUT_RESULTS_TO_FILE to false
!                       
! ***IMPORTANT NOTE:  The currently implemented expression for the Faraday Rotation integral is
!                     the usual first order *APPROXIMATION*, and NOT the exact formula.  It would be
!                      reasonably trivial to replace the first order integral with the exact expression,
!                      but it simply has not happened yet. Certainly it is at the top of the list for the
!                      next major revision to the code.   JS 12/13/05
!
!   OTHER CAVEATS:  1) ellipticity effects in polarimetry (birefringence/ Cotton-Mouton effects) are not
!                      taken into account while computing the final Faraday rotation signal size
!                   2) Information about the physical walls are not accounted for (ie no mechanism to deal
!                      with (or even NOTICE) the case when the plasma strikes the vessel wall.
!                   3) Beam paths are currently assumed to be STRAIGHT LINES.  Hence, all deviations from
!                     a straight line path due to plasma refraction are ignored.
!                   4) The B field outside the plasma is assumed to be due to EXTERNAL field coils only (i.e.
!                      effect of field sources within the plasma on points outside the last closed flux
!                      surface are assumed to be negligible
!                    5) The Cotton-Mouton logic was written very quickly, and hence is fairly simplistic.
!                       Specifically: the routine assumes that the microwave beam is RADIAL (ie perpendicular
!                       to the toroidal direction phi).  However, the logic for the Faraday rotation
!                       is NOT bound by this limitation.                          .
!
! Reference:   Rev. Sci Instrum 66 (6) June 1995 by A.J.H. Donne
!
!------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
!      USE track_mod
      USE eq_T         ! module containing VMEC output info
      USE math_utilities   ! module containing vector manipulation routines
      USE intpol_cdf       ! module containing intpol netCDF i/o routines
      USE ip_beamline  ! module containing int_pol TYPE definition
      IMPLICIT NONE

!...........module-global variables..............!
      INTEGER         :: tokamak_integrand_toggle
      INTEGER         :: iou_outdata               ! specifies file to output results
      LOGICAL                :: OUTPUT_RESULTS_TO_FILE = .true.
!      LOGICAL                :: OUTPUT_RESULTS_TO_FILE = .false.

!...........module-global variable for tokamak test case..............!
      REAL(rprec)            :: tokamak_impact_param
      REAL(rprec), PARAMETER :: TOKAMAK_PLASMA_RADIUS = 0.2_rprec


      CONTAINS

!===================================================================================
      SUBROUTINE FARADAY_ROTATION(integral_mode)
!--------------------------------------------------------------------------------
!
! Function:  Top-level "driver" subroutine for module "faraday_mod".  It  is used to compute
!            the *forward* problem (ONLY) and is used in conjunction with the specialized "my_task"
!           options ( e.g. my_task ='faraday_rotation).  It calls many of the routines developed for 
!            the reconstruction of the interferometery/polarimetry parameters.
!               Note that the V3FIT "chi-squared" routines never actually calls this driver.  The
!             interface routine for the reconstruction is ipsl_integral()
!
! INPUTS:  integral_mode       ! integer toggle to select interferometry, polarimetry, etc
!                                     F_ROTATION_ONLY       = 1
!                                     PHASE_SHIFT_ONLY      = 2
!                                     F_ROTATION_AND_PSHIFT = 3
!                                     FR_AND_COTTON_MOUTON  = 4
!
! OUTPUTS:  (none)
!
!
!
!  created Jan, 2005 by J. Shields
!
! modification history:
!
!  2/15/07 JS: Inhibited the plasma edge-finding routines since  A) TRACK() is a pain and uses differing
!              flux coord definitions and B) the utility of the routine becomes questionable in any
!              event if plasma beam deflections become non-negligible.
!
! 3/21/07  JS:  Added capability to compute Cotton-Mouton parameters w1 & w2, and added a passed
!               toggle to allow the user to select the desired integral from the V3FIT namelist 
!
!----------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
!SPH010408
      USE safe_open_mod
!      USE track_mod
      USE eq_T         ! module containing VMEC output info
      USE ip_beamline  ! module containing int_pol TYPE definition
      USE ip_global    ! module containing interfermeter/polarimeter global variables
      USE b_transform_vmec, ONLY: VMEC_B_INIT
      IMPLICIT NONE

!........passed variables........................................................!
      INTEGER, INTENT(IN) :: integral_mode    ! passed toggle for int vs. polarimetry,etc

!........local variables........................................................!

      INTEGER           :: counter = 0
      INTEGER           :: i, j, nbeams, n_edge, iflag, istat
      INTEGER           :: n_ip_units
      INTEGER, PARAMETER :: F_ROTATION_ONLY       = 1
      INTEGER, PARAMETER :: PHASE_SHIFT_ONLY      = 2
      INTEGER, PARAMETER :: F_ROTATION_AND_PSHIFT = 3
      INTEGER, PARAMETER :: FR_AND_COTTON_MOUTON  = 4

      CHARACTER(len=80)  :: cdf_infile
!      TYPE(int_pol_coll) ::  ipc_in

      CHARACTER(len=80) :: mychar1, filestr
      CHARACTER(len=80) :: mychar2, message
      CHARACTER(len=80) :: frcm_outfile
      CHARACTER(len=*), PARAMETER :: subname = 'FARADAY_ROTATION: '

      REAL(rprec), DIMENSION(3) ::   q0_vec, qf_vec, q_unit


      write(*,*) '===================================================='    
      write(*,*) ' '
      write(*,*) 'Hi.  Now in ', subname

!........**TEMPORARY KLUDGE: initialize the reading-in of wout variables DIRECTLY FROM THE WOUT FILE.....!
!.......(this is necessary because not all wout variables needed for the flux coord transformations
!...........are currently included in the "model").  JS 6/20/07
      call VMEC_B_INIT()


      if ( OUTPUT_RESULTS_TO_FILE) then
        call construct_frcm_filename(frcm_outfile)

!..........open file to output FR angle, etc to...................................!
        CALL safe_open(iou_outdata, istat, TRIM(frcm_outfile),                 &
     &                 'unknown', 'formatted')
        CALL assert_eq(0,istat,subname //                                       &
     &     ' Safe_open of output data file failed')


         write(iou_outdata, 200)
 200     FORMAT(2x, 'Z_height', 4x, 'FR_angle_(deg)', 4x, 'w1_(rad)',          &
     &          8x, 'w2_(rad)', 8x, 'w3_(rad)' ) 

!         write(iou_outdata, 201)
! 201     FORMAT('----------------------------------------------------',        &
!     &          '--------------------')
       end if




!.........compute the FR angle for the analytic test-case of a radially-symmetric tokamak.......!
!      call ANA_TOKAMAK_INTEGRAL
      
!      n_ip_units = 1
      n_ip_units = ipcoll_g%n_ip

!.........loop over all "int/pol units"...........................!
      do j = 1, n_ip_units

        write(*,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
        write(*,*) subname, 'now using intpol unit number', j

        intpol = ipcoll_g%ip(j)
        nbeams = intpol % nbeams


!!...........call dummy test routine to test new ipsl cdf stuff.  JS 6/15/07........!
!        call IPSL_CDF_TEST()


!.......allocate the global TYPE array that stores the plasma/vac edge info.......!
!......**(note: this is deallocated & reallocated for EACH loop iteration over "int/pol units")
        ALLOCATE( myedge(nbeams) )

!.........Edge-finding stuff (based on TRACK) disabled since it's probably not useful.  JS 2/15/07.....!
!!..........find all plasma/vacuum edges for the current "int/pol unit" before starting the integration..........! 
!        do i = 1, nbeams
!          write(*,*) ' '
!          write(*,*) '---------------------------------------------'
!          write(*,*) 'NOW FINDING EDGES FOR BEAM NUMBER ', i
!          q0_vec = intpol % ipbeam(i) % q0vec
!          qf_vec = intpol % ipbeam(i) % qfvec
!
!!!..........call routine to debug the pv_edge derived TYPE.   JS 12/13/06........!
!!          call PV_EDGE_TEST(i)
!          call TRACK_EDGES(i, q0_vec, qf_vec)
!!          call VMEC_EDGES(i)
!        end do


!.........call the numerical recipes trapazoid rule integration program for each beam....!

        do i = 1, nbeams
          write(*,*) ' '
          write(*,*) '------------------------------------------------'
          write(*,*) 'NOW INTEGRATING FOR BEAM NUMBER ', i

!..........set the global beamline index to the loop index.............!
          ibeam = i


!...........set global variable "ipbeam_g" instead of intpol for compatibility with..!
!............routine signal_model_compute_ipsl.  JS 6/22/07..........................!
          ipbeam_g = intpol%ipbeam(i)


!............integrate over the ith beamline..........................!
          call RUN_TRAP(integral_mode)
        end do  ! end loop over individual beams

        DEALLOCATE(myedge)
      end do    ! end loop over intpol units

      RETURN
      END SUBROUTINE FARADAY_ROTATION




!================================================================================
      Subroutine read_intpol_data_from_cdf(cdf_infile, ipc_in)
!-----------------------------------------------------------------------------------
! FUNCTION: Reads in an "int_pol_coll" derived TYPE variable (essentially an array
!           of int_pol derived type variables) from an input netCDF file.
!
!----------------------------------------------------------------------------------
      IMPLICIT NONE
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER(len=80), INTENT(IN)  :: cdf_infile   !input netCDF filename
!      TYPE(int_pol_coll), INTENT(OUT) ::  ipc_in
      TYPE(int_pol_coll), INTENT(INOUT) ::  ipc_in  !input ip collection

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER    :: len_key = 30
      INTEGER :: arr_size, ipbeam_size
      INTEGER :: nvals = 0
      INTEGER :: i, j, istat, ip_count, iou_cdf, g_size
      INTEGER, DIMENSION(3) :: DIM1, DIM2, DIM3
      CHARACTER (len=30)    :: s_name, ipc_s_name                                 
      CHARACTER (len=80)    :: l_name, ipc_l_name
      CHARACTER (len=12)    :: mychar
      CHARACTER (len=30)    :: prefix
      CHARACTER (len=*), PARAMETER  ::                                         &
     &                   subname = 'read_intpol_data_from_cdf: '

!.......begin executable code.....................!
      write(*,*) 'Now in ', subname

!.........open the newly created file and try to extract info from it. JS 2/17/05..!
      call cdf_open(iou_cdf, cdf_infile, 'r', istat)
      call intpol_cdf_read(ipc_in, iou_cdf)
      call cdf_close(iou_cdf)
 
!      write(*,*) subname, 'Completed intpol_cdf_read CALL'
!
!!........verify that the "intpol collection" variable ipc_in was read in correctly...!
!      do i = 1,1
!        write(*,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
!        write(*,*) subname, '**** READ-IN intpol element = ', i
!        write(*,*) 'READ IN nbeams = ', ipc_in%ip(i)%nbeams   
!        write(*,*) 'READ IN q0vec = ',  ipc_in%ip(i)%ipbeam(1)% q0vec   
!        write(*,*) 'READ IN qfvec = ',  ipc_in%ip(i)%ipbeam(1)% qfvec   
!
!        write(*,*) 'READ IN q0 = ',  ipc_in%ip(i)%ipbeam(1)% q0   
!        write(*,*) 'READ IN wavelength = ',                                      &
!     &             ipc_in%ip(i)%ipbeam(1)% wavelength
!!        write(*,*) 'READ IN B_ratio = ', ipc_in%ip(i)%ipbeam(1)% B_ratio
!        write(*,*) 'READ IN s_name 1= ',  ipc_in%ip(i)%ipbeam(1)% s_name
!        write(*,*) 'READ IN l_name 1= ',  ipc_in%ip(i)%ipbeam(1)% l_name
!        write(*,*) 'READ IN s_name 2= ',  ipc_in%ip(i)%ipbeam(2)% s_name
!        write(*,*) 'READ IN l_name 2= ',  ipc_in%ip(i)%ipbeam(2)% l_name
!        write(*,*) 'READ IN s_name 3= ',  ipc_in%ip(i)%ipbeam(3)% s_name
!        write(*,*) 'READ IN l_name 3= ',  ipc_in%ip(i)%ipbeam(3)% l_name
!      end do

      RETURN
      END SUBROUTINE read_intpol_data_from_cdf



!====================================================================================
      SUBROUTINE RUN_TRAP(integral_mode)
!----------------------------------------------------------------------------------
!
! FUNCTION: Interface that provides Faraday Rotation data to the Numerical Recipes
!           trapazoidal integration subroutine "qtrap" in order to obtain the 
!           line-integrated faraday rotation angle and/or the line-integrated
!           polarization phase shift.
!
! INPUTS:   compute_integral    ! integer toggle specifying which integral (FR, CM, phase shift, etc)
!                               ! the program should compute
!
!
! OUTPUTS:  (none)
!
! GLOBAL outputs:  integrand_toggle (from module ip_global)
!
! CAVEAT:  the toggle to decide whether the phase shift (inteferometry) or
!         the Faraday rotation angle (polarimetry) is desired is currently
!          HARD-WIRED at the beginning of this subroutine.
!
!  **Also note that a logical toggle SKIP_INTEGRATION_FOR_DEBUGGING is included
!    to skip the integration altogether (ie it runs FARADAY_FUNC exactly ONCE)
!     in order to make the executable finish quickly during debugging.
!
!------------------------------------------------------------------------------------
      USE stel_kinds, only : rprec, iprec, cprec
      USE ip_global    ! module containing interfermeter/polarimeter global variables
      IMPLICIT NONE

!.........passed variables.....................................................!
      INTEGER, INTENT(IN) :: integral_mode    ! passed toggle for int vs. polarimetry,etc

!.........local variables.....................................................!
      CHARACTER(len=25) ::   mychar
!      INTEGER    :: nsteps                   ! # of steps for trapazoid rule
      INTEGER    :: compute_integral         ! local toggle for int vs. polarimetry,etc
      INTEGER    :: compute_integral_default ! default setting for integral toggle
      INTEGER    :: integration_method       ! toggle for method of numerical integration

      INTEGER, PARAMETER :: F_ROTATION_ONLY       = 1
      INTEGER, PARAMETER :: PHASE_SHIFT_ONLY      = 2
      INTEGER, PARAMETER :: F_ROTATION_AND_PSHIFT = 3
      INTEGER, PARAMETER :: FR_AND_COTTON_MOUTON  = 4
      INTEGER, PARAMETER :: FARADAY_ROTATION = 1
      INTEGER, PARAMETER :: PHASE_SHIFT      = 2
      INTEGER, PARAMETER :: W1_INTEGRAND     = 3
      INTEGER, PARAMETER :: W2_INTEGRAND     = 4
      INTEGER, PARAMETER :: TRAPEZOID_RULE    = 0
      INTEGER, PARAMETER :: LINEAR_SUMMATION = 1
      INTEGER, PARAMETER :: NSTEPS = 1000      ! number of steps for linear integration routine


      REAL(KIND=rprec) :: a, b    ! lower,upper limits of definite integral

      REAL(KIND=rprec) :: faraday_integral  ! integrated faraday rotation angle (radians)
      REAL(KIND=rprec) :: faraday_degrees   ! integrated faraday rotation angle (degrees)
      REAL(KIND=rprec) :: pshift_integral   ! integrated polarization phase shift (radians)

      REAL(KIND=rprec) :: s, exact_result, diff
      REAL(KIND=rprec) :: qpath, faraday_junk
      REAL(KIND=rprec) :: z_height
      REAL(KIND=rprec), DIMENSION(3) :: w      ! stokes space matrix elements for CM & FR effects

      LOGICAL          ::  SKIP_INTEGRATION_FOR_DEBUGGING
      character(len=*), PARAMETER :: subname= 'RUN_TRAP: '


      write(*,*) subname, 'Input integral_mode = ', integral_mode

!.........set toggle to select method of integration (linear 1d grid or trapazoid rule).........!
      integration_method = TRAPEZOID_RULE
!      integration_method = LINEAR_SUMMATION

!...........set default for the desired integrated quantities ............!
      compute_integral_default = F_ROTATION_ONLY
!      compute_integral_default = PHASE_SHIFT_ONLY 
!      compute_integral_default = F_ROTATION_AND_PSHIFT
!      compute_integral_default = FR_AND_COTTON_MOUTON

      SELECT CASE(integral_mode)
      CASE(1:4) 
        compute_integral = integral_mode
      CASE DEFAULT 
         compute_integral = compute_integral_default
      END SELECT

!.....debugging option to skip trapazoidal integration & call FARADAY_FUNC *once*....!
!      SKIP_INTEGRATION_FOR_DEBUGGING = .true.
      SKIP_INTEGRATION_FOR_DEBUGGING = .false.


!......assign values to lower/upper limits of the definite integral...........!
      a = ipbeam_g % q0
      b = a + ipbeam_g % qdist
      write(*,*) subname, 'lower limit a = ', a
      write(*,*) subname, 'upper limit b = ', b


!  debugging call to test subroutine SIMPLE_TRAPEZOID
!.......calc integral using a simple loop and the trapazoidal rule......!
! ***WARNING: be careful that the input value of nsteps does not cause the large 
!              integer "it" in trapzd to exceed the machine MAX INTEGER!!
!      nsteps = 25
!      call SIMPLE_TRAPEZOID(test_func,a,b,s,nsteps)
!
!.........debugging call to test subroutine qtrap.........................!
!......Instead of providing a fixed number of steps, let the number of required
!......steps depend upon the rate at which the integral terms converge......!
!      call qtrap(test_func,a,b,s)
!      exact_result = 0.33333333 * ( b**3.0 - a**3.0 )
!      diff = ABS( exact_result - s)
!      write(*,*) 'ESTIMATED integral value = ', s
!      write(*,*) 'ACTUAL integral value= ', exact_result
!      write(*,*) 'difference btwn estimated & exact results = ', diff



!----------------------------
!&&&&& below this point are actual Faraday rotation calls (as opposed to calls
!...........for testing subroutine qtrap().....................................!
 

      if ( SKIP_INTEGRATION_FOR_DEBUGGING ) then 
!...........then allow FARADAY_FUNC to run only ONE time (for debugging)..........!
        write(*,*) 'WARNING!  DEBUG MODE IN USE.  NO INTEGRATION DONE!'
        w = 0.0_rprec
        qpath = 0.3
        integrand_toggle = FARADAY_ROTATION
        faraday_junk = FARADAY_FUNC(qpath)
        write(*,*) 'integrated Faraday Rot. angle (in rad) = ',                &
     &              faraday_junk
      else

        if ( compute_integral == F_ROTATION_ONLY) then
          w = 0.0_rprec

          integrand_toggle = FARADAY_ROTATION
          if ( integration_method == TRAPEZOID_RULE ) then
!!           call qtrap(test_func,a,b,faraday_integral)
            call qtrap(FARADAY_FUNC,a,b,faraday_integral)
          else if ( integration_method == LINEAR_SUMMATION ) then
            call stepsum(FARADAY_FUNC,a,b,NSTEPS, faraday_integral)        
          else
            write(*,*) subname, 'WARNING! ERROR IN integration_method!'
          end if

!............compute the Faraday rotation stokes vector parameter from faraday_integral.................!
          w(3) = 2.0 * faraday_integral
          write(*,*) 'integrated Faraday Rot. angle (in rad) = ',                       &
     &               faraday_integral
          write(*,*) 'integrated Faraday Rot. angle (in deg) = ',                       &
     &               (faraday_integral *(180.0/3.14159))

        else if ( compute_integral == PHASE_SHIFT_ONLY) then
          w = 0.0_rprec
          integrand_toggle = PHASE_SHIFT
          if ( integration_method == TRAPEZOID_RULE ) then
            call qtrap(FARADAY_FUNC,a,b,pshift_integral)
          else if ( integration_method == LINEAR_SUMMATION ) then
            call stepsum(FARADAY_FUNC,a,b,NSTEPS, pshift_integral)        
          else
            write(*,*) subname, 'WARNING! ERROR IN integration_method!'
          end if

          write(*,*) 'integrated Phase Shift = ', pshift_integral

        else if ( compute_integral == F_ROTATION_AND_PSHIFT) then
          w = 0.0_rprec

!.............compute faraday rotation angle first.........................!
          integrand_toggle = FARADAY_ROTATION
          if ( integration_method == TRAPEZOID_RULE ) then
            call qtrap(FARADAY_FUNC,a,b,faraday_integral)
          else if ( integration_method == LINEAR_SUMMATION ) then
            call stepsum(FARADAY_FUNC,a,b,NSTEPS, faraday_integral)        
          else
            write(*,*) subname, 'WARNING! ERROR IN integration_method!'
          end if

!............compute the Faraday rotation stokes vector parameter from faraday_integral.................!
          w(3) = 2.0 * faraday_integral

!.............now do a second integral using the phase shift integrand........!
          integrand_toggle = PHASE_SHIFT
          if ( integration_method == TRAPEZOID_RULE ) then
            call qtrap(FARADAY_FUNC,a,b,pshift_integral)
          else if ( integration_method == LINEAR_SUMMATION ) then
            call stepsum(FARADAY_FUNC,a,b,NSTEPS, pshift_integral)        
          else
            write(*,*) subname, 'WARNING! ERROR IN integration_method!'
          end if

          write(*,*) 'integrated Faraday Rot. angle (in rad) = ',                       &
     &               faraday_integral
          write(*,*) 'integrated Faraday Rot. angle (in deg) = ',                       &
     &               (faraday_integral *(180.0/3.14159))
          write(*,*) 'integrated Phase Shift = ', pshift_integral

        else if ( compute_integral == FR_AND_COTTON_MOUTON) then

!.............compute faraday rotation angle (in radians) first.........................!
          integrand_toggle = FARADAY_ROTATION
          if ( integration_method == TRAPEZOID_RULE ) then
            call qtrap(FARADAY_FUNC,a,b,faraday_integral)
          else if ( integration_method == LINEAR_SUMMATION ) then
            call stepsum(FARADAY_FUNC,a,b,NSTEPS, faraday_integral)        
          else
            write(*,*) subname, 'WARNING! ERROR IN integration_method!'
          end if

!............compute the Faraday rotation stokes vector parameter from faraday_integral.................!
          w(3) = 2.0 * faraday_integral

!.............now do a second integral to compute the Cotton-Mouton w1 parameter........!
          integrand_toggle = W1_INTEGRAND
          if ( integration_method == TRAPEZOID_RULE ) then
            call qtrap( FARADAY_FUNC,a,b, w(1) )
          else if ( integration_method == LINEAR_SUMMATION ) then
            call stepsum(FARADAY_FUNC,a,b,NSTEPS, w(1) )        
          else
            write(*,*) subname, 'WARNING! ERROR IN integration_method!'
          end if

!.............now do a second integral to compute the Cotton-Mouton w2 parameter........!
          integrand_toggle = W2_INTEGRAND
          if ( integration_method == TRAPEZOID_RULE ) then
            call qtrap( FARADAY_FUNC,a,b, w(2) )
          else if ( integration_method == LINEAR_SUMMATION ) then
            call stepsum(FARADAY_FUNC,a,b,NSTEPS, w(2) )        
          else
            write(*,*) subname, 'WARNING! ERROR IN integration_method!'
          end if

          write(*,*) 'integrated Faraday Rot. angle (in rad) = ',                       &
     &               faraday_integral

          write(*,*) 'integrated Faraday Rot. angle (in deg) = ',                       &
     &               (faraday_integral *(180.0/3.14159))

          write(*,*) 'integrated w1 (in rad) = ', w(1)
          write(*,*) 'integrated w2 (in rad) = ', w(2)

          write(*,*) '             '

          write(*,*) 'ratio of w1/w3 = ', w(1)/w(3)
          write(*,*) 'ratio of w2/w3 = ', w(2)/w(3)
          write(*,*) 'ratio of w2/w1 = ', w(2)/w(1)
          
        else
          write(*,*) 'ERROR IN RUN_TRAP INTEGRAL ASSIGNMENT!'
        end if

      end if

!..........write out results for this particular beamline to output file.................!
      if ( OUTPUT_RESULTS_TO_FILE) then

!...........determine Z height of the beam starting point............................!
        z_height = ipbeam_g % q0vec(3)

        faraday_degrees = faraday_integral * (180.0/3.14159)
        write(iou_outdata, 100) z_height, faraday_degrees, w(1),               &
     &                          w(2), w(3)
 100    format ( f9.5, 1x, 4es16.6 )
      end if

      END SUBROUTINE RUN_TRAP

!====================================================================================
      SUBROUTINE ANA_TOKAMAK_INTEGRAL
!----------------------------------------------------------------------------------
!
! FUNCTION: driver routine to compute the Faraday Rotation angle for my analytically-derived 
!           "radially-symmetric tokamak" case (via a simple, 1d numerical integration).  
!           Purpose is to provide an independent test case for the
!           general Faraday Rotation code
!              *IMPORTANT NOTE: the beam wavelength is the ONE & ONLY quantity
!           that is taken from the data stored in the ipbeam derived TYPE.
!
! INPUTS:   (none)
! OUTPUTS:   (none)
!
! OUTPUTS (to std output):  faraday_integral  !integrated faraday rotation angle
!                           pshift_integral  !integrated polarization phase shift
!
! CAVEAT:  the toggle to decide whether the phase shift (inteferometry) or
!         the Faraday rotation angle (polarimetry) is desired is currently
!          HARD-WIRED at the beginning of this subroutine.
!
!  **Also note that a logical toggle SKIP_INTEGRATION_FOR_DEBUGGING is included
!    to skip the integration altogether (ie it runs FARADAY_FUNC exactly ONCE)
!     in order to make the executable finish quickly during debugging.
!
!------------------------------------------------------------------------------------
      USE stel_kinds, only : rprec, iprec, cprec
!      USE ip_global    ! module containing interfermeter/polarimeter global variables
      IMPLICIT NONE

!............local variables............................................!
      INTEGER    :: i
      INTEGER    :: compute_integral    ! toggle for int vs. polarimetry
!      INTEGER    :: nsteps              ! # of steps for trapazoid rule
      INTEGER    :: nbeamlines          ! # of beams passing through the plasma

      INTEGER, PARAMETER :: F_ROTATION_ONLY       = 1
      INTEGER, PARAMETER :: PHASE_SHIFT_ONLY      = 2
      INTEGER, PARAMETER :: F_ROTATION_AND_PSHIFT = 3
      INTEGER, PARAMETER :: FARADAY_ROTATION = 1
      INTEGER, PARAMETER :: PHASE_SHIFT      = 2

      REAL(KIND=rprec) :: xzero   ! right edge of plasma (assuming center is at origin)
      REAL(KIND=rprec) :: a, b    ! lower,upper limits of definite integral
      REAL(KIND=rprec) :: stepsize              ! stepsize to compute impact parameter
      REAL(KIND=rprec) :: start_height          ! initial value of impact parameter
      REAL(KIND=rprec) :: qpath, faraday_junk

      LOGICAL          ::  SKIP_INTEGRATION_FOR_DEBUGGING
      character(len=*), PARAMETER :: subname= 'ANA_TOKAMAK_INTEGRAL: '

      write(*,*) 'now in ', subname

!...........select which integrated quantities are desired............!
      compute_integral = F_ROTATION_ONLY
!      compute_integral = PHASE_SHIFT_ONLY 
!      compute_integral = F_ROTATION_AND_PSHIFT


!.....debugging option to skip trapazoidal integration & call FARADAY_FUNC *once*....!
!      SKIP_INTEGRATION_FOR_DEBUGGING = .true.
      SKIP_INTEGRATION_FOR_DEBUGGING = .false.



 

      if ( SKIP_INTEGRATION_FOR_DEBUGGING ) then 
!...........then allow FARADAY_FUNC to run only ONE time (for debugging)..........!
        write(*,*) 'WARNING!  DEBUG MODE IN USE.  NO INTEGRATION DONE!'
        qpath = 0.3
        tokamak_integrand_toggle = FARADAY_ROTATION
        faraday_junk = TOKAMAK_FUNC(qpath)
        write(*,*) 'integrated Faraday Rot. angle (in rad) = ',                &
     &              faraday_junk
      else
        start_height = 0.00
!        start_height = 0.01
        stepsize = 0.02

!        nbeamlines = 1
        nbeamlines = 11

        do i = 1,nbeamlines

!..............compute the impact parameter (ie vertical height above plasma center)....!
          tokamak_impact_param = start_height + stepsize * REAL(i-1)

          write(*,*) 'tokamak_impact_param = ', tokamak_impact_param 
!..............assign values to lower/upper limits of the definite integral...........!
          xzero = SQRT( TOKAMAK_PLASMA_RADIUS**2.0 -                               &
     &                tokamak_impact_param**2.0 )

!          a = 0.0
!          b = 0.636
          a = -1.0_rprec * xzero
          b = xzero

          call TOKAMAK_RUN_TRAP(a,b, compute_integral)
        end do
      end if

      END SUBROUTINE ANA_TOKAMAK_INTEGRAL

!====================================================================================
      SUBROUTINE TOKAMAK_RUN_TRAP(a, b, compute_integral)
!----------------------------------------------------------------------------------
!
! FUNCTION: Uses Trapazoid Rule to numerically integrate the radially symmetric
!           tokamak test case.
!
! OUTPUTS:  faraday_integral  ! integrated faraday rotation angle
!           pshift_integral   ! integrated polarization phase shift
!
!  Created 2/10/07 by J. Shields, based on RUN_TRAP
!
!------------------------------------------------------------------------------------
      USE stel_kinds, only : rprec, iprec, cprec
      USE ip_global    ! module containing interfermeter/polarimeter global variables
      IMPLICIT NONE

!.........dummy arguments.........................................!
      INTEGER, INTENT(IN)   :: compute_integral ! toggle for int vs. polarimetry
      REAL(KIND=rprec), INTENT(IN) :: a                ! lower limit of definite integral
      REAL(KIND=rprec), INTENT(IN) :: b                ! upper limit of definite integral


!..........local variables.........................................!

      INTEGER    :: nsteps              ! # of steps for trapazoid rule

      INTEGER, PARAMETER :: F_ROTATION_ONLY       = 1
      INTEGER, PARAMETER :: PHASE_SHIFT_ONLY      = 2
      INTEGER, PARAMETER :: F_ROTATION_AND_PSHIFT = 3
      INTEGER, PARAMETER :: FARADAY_ROTATION = 1
      INTEGER, PARAMETER :: PHASE_SHIFT      = 2

      REAL(KIND=rprec) :: faraday_integral  ! integrated faraday rotation angle
      REAL(KIND=rprec) :: pshift_integral   ! integrated polarization phase shift

!      REAL(KIND=rprec) :: s, exact_result, diff
      REAL(KIND=rprec) :: qpath, faraday_junk

      LOGICAL          ::  SKIP_INTEGRATION_FOR_DEBUGGING
      character(len=*), PARAMETER :: subname= 'TOKAMAK_RUN_TRAP: '



!      write(*,*) subname, 'input a = ', a, 'input b = ', b

      if ( compute_integral == F_ROTATION_ONLY) then
        tokamak_integrand_toggle = FARADAY_ROTATION
        call qtrap(TOKAMAK_FUNC,a,b,faraday_integral)
!        write(*,*) 'ANALTYIC TOKAMAK Faraday Rot. angle (in rad) = ',                       &
!     &               faraday_integral
        write(*,*) 'ANALYTIC TOKAMAK Faraday Rot. angle (in deg) = ',                       &
     &               (faraday_integral *(180.0/3.14159))

      else if ( compute_integral == PHASE_SHIFT_ONLY) then
        tokamak_integrand_toggle = PHASE_SHIFT
        call qtrap(TOKAMAK_FUNC,a,b,pshift_integral)
        write(*,*) 'ANALYTIC TOKAMAK Phase Shift = ', pshift_integral

      else if ( compute_integral == F_ROTATION_AND_PSHIFT) then

!............compute faraday rotation angle first.........................!
        tokamak_integrand_toggle = FARADAY_ROTATION
        call qtrap(TOKAMAK_FUNC,a,b,faraday_integral)

!............now do a second integral using the phase shift integrand........!
        tokamak_integrand_toggle = PHASE_SHIFT 
        call qtrap(TOKAMAK_FUNC,a,b,pshift_integral)

!        write(*,*) 'ANALYTIC TOKAMAK Faraday Rot. angle (in rad) = ',                       &
!     &               faraday_integral
        write(*,*) 'ANALYTIC TOKAMAK Faraday Rot. angle (in deg) = ',                       &
     &               (faraday_integral *(180.0/3.14159))
        write(*,*) 'ANALYTIC TOKAMAK Phase Shift = ', pshift_integral

      else
        write(*,*) subname, 'ERROR IN INTEGRAL ASSIGNMENT!'
      end if


      END SUBROUTINE TOKAMAK_RUN_TRAP


!=============================================================================
      REAL(KIND=rprec) FUNCTION FARADAY_FUNC(q)
!--------------------------------------------------------------------------
!
! FUNCTION: Determines the value of the integrand in the Faraday rotation angle integral
!           (OR the polarization phase shift integral) at a given point "q" along the 
!           beam path ( where "dq" is defined to be the line element along the [straight 
!           line] path of the beam).
!               The integrand that is calculated (Faraday Rotation or Phase shift) is 
!            determined by the value of the global variable "integrand_toggle".
!
! INPUTS: q   ! step point along the beam path
!
! module-global inputs: integrand_toggle
!
! 
! LOGIC:  Faraday Rot Equation is:   FR_angle = (const) * INTEGRAL_X{ ne *( B .DOT. dx) }
!
!         Phase Shift Equation is:   delta_phi = (const) * INTEGRAL_X( ne * dx)
!
!
!  Created 5/11/05 by J. Shields
!
!   Reference: A.J.H. Donne, Rev Sci Instrum 66(6), June 1995
!
! CAVEAT: The current version of this procedure assumes a STRAIGHT LINE PATH from
!         the beam start point (q0) and its endpoint (qf).
!
! CAVEAT #2: Note that the current equations are SIMPLIFIED versions of the full
!      Appleton-Hartree Equation.  Hence, their validity is confined to certain
!       wavelengths  (see Donne 1995 above)
!
!-------------------------------------------------------------------------
      USE stel_kinds, only : rprec, iprec, cprec
      USE ip_global    ! module containing interfermeter/polarimeter global variables
      USE b_transform_vmec, ONLY: TRANS_CAR2CYL, VMEC_CAR2FLX
      IMPLICIT NONE

!..........dummy variables..............................................!

      REAL(rprec) :: q                ! distance along the beam path

!..........local variables..............................................!

      INTEGER         :: q_inout
      INTEGER         :: ierr, iflag
      INTEGER, PARAMETER     :: FARADAY_ROTATION = 1
      INTEGER, PARAMETER     :: PHASE_SHIFT      = 2
      INTEGER, PARAMETER     :: W1_INTEGRAND     = 3
      INTEGER, PARAMETER     :: W2_INTEGRAND     = 4
      INTEGER, PARAMETER     :: INSIDE_PLASMA  = 1
      INTEGER, PARAMETER     :: OUTSIDE_PLASMA = 2
      INTEGER, PARAMETER     :: LOCATION_UNKNOWN = -99
      INTEGER, PARAMETER     :: CARTESIAN = 2

      REAL(rprec) :: omega          ! angular frequency
      REAL(rprec) :: nc             ! critical density
      REAL(rprec) :: ne             ! electron density
      REAL(rprec) :: B_parallel     ! magnitude of B field parallel to beam path
      REAL(rprec) :: frot_const, pshift_const
      REAL(rprec) :: cm_const       ! Cotton-Mouton coeff.  corresponds to "C1" in Segre, 1995
      REAL(rprec) :: func_output

      REAL(rprec), DIMENSION(3) ::                                             &
     &            q0_vec            ! Cartesian position vector of beam start pt
      REAL(rprec), DIMENSION(3) ::                                             &
     &            qf_vec            ! Cartesian position vector of beam end pt
      REAL(rprec), DIMENSION(3) ::                                             &
     &            q_unit            ! Cartesian unit vector parallel to the beam path 
      REAL(rprec), DIMENSION(3) ::                                             &
     &            xcar              ! Cartesian coords of a point a dist "q" along beam
      REAL(rprec), DIMENSION(3) ::                                             &
     &            xcyl              ! cylin coords of a point a dist "q" along beam
      REAL(rprec), DIMENSION(3) ::                                             &
     &            xflux             ! Flux coords of a point a dist "q" along beam
      REAL(rprec), DIMENSION(3) ::                                             &
     &            Bcar              ! Cartesian NET B field vector at point "q"
      REAL(rprec), DIMENSION(3) ::                                             &
     &            Bcar_ext          ! Cartesian EXTERNAL B field vector at point "q"
      REAL(rprec), DIMENSION(3) ::                                             &
     &            B_cyl             ! NET B field vector at point "q" in cylin coords
      REAL(rprec), DIMENSION(3) ::                                             &
     &            B_factor          ! B dependence in respective (w1,w2,w3) integrands

      REAL(rprec), PARAMETER :: ECHRG = 1.602176e-19
      REAL(rprec), PARAMETER :: ME = 9.1093826e-31
      REAL(rprec), PARAMETER :: EPZERO = 8.854188e-12
      REAL(rprec), PARAMETER :: C = 299792458.0      ! speed of light (in meter/sec)
      REAL(rprec)            :: lambda               ! beam wavelength (in METERS)

      CHARACTER(len=120) :: message
      CHARACTER(len=*), PARAMETER :: subname = 'FARADAY_FUNC: '


!      write(*,*) subname, 'input q value = ' , q

!.......extract the beam wavelength from memory & convert it from nm to meters.......!
      lambda = 1.0e-9 * ipbeam_g % wavelength
      omega = twopi * (C/lambda)

!..............compute the coefficients of the FR and phase shift integrals........!
      nc = (omega*omega*ME*EPZERO)/(ECHRG*ECHRG)
      frot_const   = ECHRG/(2.0*C*ME* nc)
      pshift_const = omega / (2.0*C* nc)
      cm_const = (ECHRG/(ME* omega)) * frot_const


      q0_vec = ipbeam_g % q0vec
      qf_vec = ipbeam_g % qfvec

!.........convert path variable "q" into a 3D position vector in Cartesian coords......!
      call PATH_TO_CARTESIAN(q, q0_vec, qf_vec, q_unit, xcar)   

!..........convert Cartesian xcar vector into cylin vector xcyl..................!
      call TRANS_CAR2CYL(xcar, xcyl)

!.........determine if point "q" is inside (1) or outside (2) the plasma...........!
      call PLASMA_INOUT(q, xcyl, q_inout, ierr)
      if( ierr .ne. 0) write(*,*) 'error in PLASMA_INOUT!'

!      write(*,*) subname, 'q_inout = ', q_inout

!.......convert Cartesian vector xcar into flux coordinates..........................!
      if( q_inout == INSIDE_PLASMA) then 
!        call CAR2FLX(xcar, xflux, iflag, message)
        call VMEC_CAR2FLX(xcar, xflux, iflag)
        if( iflag .ne. 0) write(*,*) 'CAR2FLX error in FARADAY_FUNC!'
!        xflux_rho = xflux(1)
      else
!...........if pt outside plasma, set components of the flux coord vector to "junk" values.  JS 12/14/06....!
!        call TRACK_TRANS2FLX(xcar, xflux, iflag, message,                      &
!     &                       K_COORD = CARTESIAN)
        xflux = -9999.0_rprec
      end if


!............determine the electron density at point "q"..............................!
      ne = E_DENSITY(q, q0_vec, qf_vec, q_unit, q_inout, xflux)
 
!      write(*,*) 'ne in FARADAY_FUNC = ' , ne

!...........calculate the net Cartesian B field vector at point "q"...........!
      call B_CARTESIAN(q, xcar, Bcar, q_inout, B_CYL=B_cyl)


      call INTEGRAND_B_FACTOR( xcyl, q_unit, Bcar, B_cyl, B_factor) 


      if( integrand_toggle == FARADAY_ROTATION ) then

        B_parallel = B_factor(3)

!........calculate the integrand of the Faraday rotation angle integral............!
        func_output = frot_const * ne * B_parallel


      else if( integrand_toggle == PHASE_SHIFT ) then
        func_output = -1.0 * pshift_const * ne 

      else if( integrand_toggle == W1_INTEGRAND ) then
        func_output = cm_const * ne * B_factor(1) 

      else if( integrand_toggle == W2_INTEGRAND ) then
        func_output = 2.0_rprec * cm_const * ne * B_factor(2) 

      else
        write(*,*) ' Error in FARADAY_FUNC integrand assignment!'
      end if

      FARADAY_FUNC = func_output
      
      return
      END FUNCTION FARADAY_FUNC

!=============================================================================
      REAL(KIND=rprec) FUNCTION TOKAMAK_FUNC(x)
!--------------------------------------------------------------------------
!
! FUNCTION: Integrand for the ANALYTICALLY-DERIVED case of a tokamak distribution with
!           a perfectly circularly symmetric poloidal B field.
!             This is merely a SPECIAL TEST CASE that is used to valid the general code.  
!
!        NOTE:    J(r) = JZERO * { 1 - (r^2/R^2)^GAMMA}          ! current density
!                 N(r) = NZERO * { 1 - (r^2/R^2)^TAU}            ! electron density
!
!                    where:  r = dist from center, R = plasma outer radius, 
!                            and GAMMA and TAU are constants.
!
! LOGIC:   Determines the value of the integrand in the Faraday rotation angle integral
!           (OR the polarization phase shift integral) at a given point "x" along the 
!           beam path ( where "dx" is defined to be the line element along the [perfectly 
!           horizontal] path of the beam).
!              The integrand that is calculated (Faraday Rotation or Phase shift) is 
!           determined by the value of the global variable "integrand_toggle".
!              The horizontal distance "x" is related to the radial distance from the plasma 
!           center (r) & the impact parameter (a) via:
!
!            x^2 = r^2 - a^2
!
!         Faraday Rot Equation is:   FR_angle = (const) * INTEGRAL_X{ (ne * B_parallel)*dx }
!
!         Phase Shift Equation is:   delta_phi = (const) * INTEGRAL_X( ne * dx)
!
!            where B_parallel = -sin_theta * B_poloidal = (-a/r) * B_poloidal
!
!
! INPUTS: x   ! step point along the horizontal beam path
!
! module-global inputs: tokamak_integrand_toggle
!                       TOKAMAK_PLASMA_RADIUS
!                       tokamak_impact_param
!                       intpol % ipbeam(1) % wavelength                       
!
!
!  Created 2/9/07 by J. Shields
!
!   Reference: A.J.H. Donne, Rev Sci Instrum 66(6), June 1995
!
! CAVEAT: The current version of this procedure assumes a STRAIGHT LINE PATH from
!         the beam start point (q0) and its endpoint (qf).
!
! CAVEAT #2: Note that the current equations are SIMPLIFIED versions of the full
!      Appleton-Hartree Equation.  Hence, their validity is confined to certain
!       wavelengths  (see Donne 1995 above)
!
! CAVEAT #3: Although this routine *technically* takes the wavelength info from the
!            ip_v3rfun_*.nc file, it assumes that the wavelength for ALL paths is the same.
!            Hence, it ONLY actually reads in the wavelength of the *FIRST* beam path
!            (which is stored in TYPE variable intpol%ipbeam(1) )
!
!-------------------------------------------------------------------------
      USE stel_kinds, only : rprec, iprec, cprec
      USE ip_global  ! interfermeter/polarimeter global variables (used for wavelength)
      IMPLICIT NONE

!..........dummy variables..............................................!
      REAL(rprec) :: x                ! distance along the horizontal beam path

!..........local variables..............................................!

      INTEGER            :: ierr, iflag
      INTEGER, PARAMETER :: FARADAY_ROTATION = 1
      INTEGER, PARAMETER :: PHASE_SHIFT      = 2

      REAL(rprec) :: omega          ! angular frequency
      REAL(rprec) :: nc             ! critical density
      REAL(rprec) :: ne             ! net electron density
      REAL(rprec) :: j_current      ! net current density
      REAL(rprec) :: jzero          ! current density at r= 0
      REAL(rprec) :: B_parallel     ! magnitude of B field parallel to beam path (i.e. B_x)
      REAL(rprec) :: Bzero          ! B(r=0)
      REAL(rprec) :: frot_const, pshift_const
      REAL(rprec) :: func_output
      REAL(rprec) :: outer_radius   ! outer radius of the tokamak plasma
      REAL(rprec) :: a              ! impact parameter (ie vertical height above plasma center)
      REAL(rprec) :: arg
      REAL(rprec) :: plasma_area    ! net area of plasma

      REAL(rprec), PARAMETER :: ECHRG = 1.602176e-19
      REAL(rprec), PARAMETER :: ME = 9.1093826e-31
      REAL(rprec), PARAMETER :: EPZERO = 8.854188e-12
      REAL(rprec), PARAMETER :: MUZERO = PI*4.0e-7
      REAL(rprec), PARAMETER :: C = 299792458.0      ! speed of light (in meter/sec)

!      REAL(rprec), PARAMETER :: I_TOT = 2.25e5              ! TOTAL toroidal current
      REAL(rprec), PARAMETER :: I_TOT = 2.25e7              ! TOTAL toroidal current
!      REAL(rprec), PARAMETER ::  I_TOT = 1.0e4               ! TOTAL toroidal current
      REAL(rprec), PARAMETER :: NZERO = 4.0e18               ! electron density at plasma center
      REAL(rprec), PARAMETER :: GAMMA = 3.0_rprec/7.0_rprec ! current  density exponent
!      REAL(rprec), PARAMETER :: GAMMA = 0.30                ! current  density exponent
      REAL(rprec), PARAMETER :: TAU   = 5.0_rprec           ! electron density exponent

      REAL(rprec) :: lambda         ! beam wavelength (in METERS)

      CHARACTER(len=120) :: message
      CHARACTER(len=*), PARAMETER :: subname = 'TOKAMAK_FUNC: '


!      write(*,*) subname, 'input q value = ' , q


!.......extract the beam wavelength from memory & convert it from nm to meters.......!
!.....wavelength temporarily HARD-WIRED for debugging.  JS 2/10/07.........!
!      lambda = 1.0e-9 * intpol % ipbeam(1) % wavelength
      lambda = 2.1e-3
      omega = twopi * (C/lambda)

!...........fill local variables using module-global values.................!
      outer_radius = TOKAMAK_PLASMA_RADIUS
      a = tokamak_impact_param

!.........compute the current density coeff in terms of the NET toroidal current (I_tot)...!
      plasma_area = PI * outer_radius**2.0
      jzero = ( (GAMMA + 1.0)/GAMMA ) * (I_TOT/plasma_area)

!..............compute the coefficients of the FR and phase shift integrals........!
      nc = (omega*omega*ME*EPZERO)/(ECHRG*ECHRG)
      frot_const   = ECHRG/(2.0*C*ME* nc)
      pshift_const = omega / (2.0*C* nc)
      Bzero = -0.5_rprec * MUZERO * jzero * a

      arg = (x**2.0 + a**2.0) / (outer_radius**2.0)

!............determine the electron density at point "x"..............................!
      ne = NZERO * ( 1.0_rprec - arg**TAU )

!      write(*,*) subname, 'ne = ' , ne


      if( tokamak_integrand_toggle == FARADAY_ROTATION ) then

!............determine the parallel B field (B_x) at point "x"..............................!
        B_parallel = Bzero*( 1.0_rprec - (arg**GAMMA)/(GAMMA+1.0_rprec))

!...........calculate the integrand of the Faraday rotation angle integral............!
        func_output = frot_const * ne * B_parallel

      else if( tokamak_integrand_toggle == PHASE_SHIFT ) then
        func_output = -1.0 * pshift_const * ne 
      else
        write(*,*) subname, 'Error in integrand assignment!'
      end if

      TOKAMAK_FUNC = func_output
      
      return
      END FUNCTION TOKAMAK_FUNC


!===================================================================================
      SUBROUTINE PATH_TO_CARTESIAN(q, q0_vec, qf_vec, q_unit, xcar)   
!--------------------------------------------------------------------------------
!
! FUNCTION:  Determines the Cartesian components of the position vector "xcar"
!            that corresponds to the point a distance "q" along the beam path.
!            
! 
! LOGIC:  The beam path is defined by the starting point corresponding to q0_vec
!          and the direction given by the unit vector q_unit.  Therefore, by simple
!          vector addition, the coordinates of the point a distance "q" from q0_vec
!          in the direction of q_unit must be given by
!           xcar(i) = q0_vec(i) + { q * q_unit(i) }
!
! Created 5/24/05 by J. Shields
!
! *NOTE: if it later turns out that more sophisticated calculations are required (such
!         as the computation of Jacobian derivatives) consider using the TRACK_D()
!         subroutine in TRACK rather than re-inventing the wheel.
!
!----------------------------------------------------------------------------------
!      USE stel_kinds
!      USE stel_constants
!      USE math_utilities
      IMPLICIT NONE

      REAL(rprec), INTENT(IN) :: q    ! distance along beam path

      REAL(rprec), DIMENSION(3), INTENT(IN) ::                                 &
     &            q0_vec           ! Cartesian POSITION vector of beam starting point
      REAL(rprec), DIMENSION(3), INTENT(IN) ::                                 &
     &            qf_vec           ! Cartesian POSITION vector of beam end point
      REAL(rprec), DIMENSION(3), INTENT(OUT) ::                                 &
     &            q_unit           ! Cartesian unit vector parallel to beam path
      REAL(rprec), DIMENSION(3), INTENT(OUT) ::                                 &
     &            xcar             ! Cartesian coordinate POSITION vector


!..............local variables......................................!
      REAL(rprec), DIMENSION(3) :: qdiff
      REAL(rprec)               :: qdiff_mag
      

!........calculate the unit vector along the line btwn q0_vec and qf_vec.........!
      qdiff = qf_vec - q0_vec 
      qdiff_mag = MAGNITUDE(qdiff)
      q_unit = qdiff / qdiff_mag

!........calculate the Cartesian vector at point "q"............................!
      xcar = q0_vec + q * q_unit
 
      RETURN
      END SUBROUTINE PATH_TO_CARTESIAN


!=============================================================================
      SUBROUTINE TRACK_EDGES(i_beam, q0_vec, qf_vec)
!-------------------------------------------------------------------------------
! FUNCTION: uses the TRACK module to find the (approximate) outer edges of the plasma & saves
!           the edge information to global derived TYPE variable "myedge".
!               The idea is to have a more sophisticated way of determining whether a given point
!           is inside or outside the plasma than simply assuming that anywhere that VMEC flux coords
!           can't converge is outside the plasma.
!
!  LOGIC:  Uses TRACK to find the intersections between the beam path and flux surfaces.

!  STATUS:  STILL UNDER DEVELOPMENT.
!            As of 12/14/06, this routine has NOT been fully tested.  The outputs of TRACK are still somewhat
!           at odds with my interpretation of what its documentation implies and I have not yet gotten around
!           to reconciling the discrepancy yet.  
! 
! A note about flux coordinate conventions:
!           TRACK uses the AJAX module for flux coordinate transformations and AJAX
!            defines its flux coordinates with an extra factor of rho^m from the VMEC convention.  However,
!            since rho^m = rho when rho=1.0, it is STILL possible to use TRACK to find the coordinates
!            of the outermost flux surface.
!
! INPUTS: i_beam, q0_vec, qf_vec
! OUTPUTS: (none)
!
! global outputs: myedge
!
! created 2005 by J. Shields
!
! revision notes: 
!    JS 6/8/05 added code from EDGE_POINTS routine into this subroutine .
!
!-------------------------------------------------------------------------------
      USE SPEC_KIND_MOD
      USE TRACK_MOD
      USE AJAX_MOD
      use ip_global
      IMPLICIT NONE


!declaration of input variables
      INTEGER :: i_beam     ! index of interferometer beam being integrated over

      REAL(KIND=rspec), DIMENSION(3), INTENT(IN) ::                            &
     & q0_vec, qf_vec       ! position vectors of beam start/end points


!Declaration of local variables
      INTEGER, PARAMETER :: INSIDE_PLASMA  = 1
      INTEGER, PARAMETER :: OUTSIDE_PLASMA = 2
      INTEGER, PARAMETER :: VACUUM_TO_PLASMA = 1
      INTEGER, PARAMETER :: PLASMA_TO_VACUUM = 2
      INTEGER, PARAMETER :: ENDPT_INSIDE_PLASMA   = -11
      INTEGER, PARAMETER :: ENDPT_OUTSIDE_PLASMA  = -21

      INTEGER ::                                                               &
     & k_equil,k_seg,                                                          &
     & n_rho, n_int, n_seg, iflag    ! TRACK i/o variables

      INTEGER ::                                                               &
     & beginpoint, endpoint, currentpoint, previous_pt

      CHARACTER(len=120) ::                                                    &
     & message

      INTEGER ::                                                               &
     & i, j,                                                                   & 
     & n_guess, n_temp

      INTEGER :: n_edge      ! number of plasma/vacuum edges (plus 2 endpoints)

      INTEGER, ALLOCATABLE ::                                                  &
     & irho_int(:),                                                            &
     & edge_index(:)        ! TRACK array index corresponding to the ith edge boundary point

      INTEGER, ALLOCATABLE ::       edge_type(:)

      REAL(KIND=rspec), ALLOCATABLE ::                                         &
     & rho(:), r_seg(:,:), s_int(:)

      REAL(KIND=rspec), ALLOCATABLE ::                                         &
     & rcar_int(:,:),rcyl_int(:,:),rflx_int(:,:),sdotb_int(:),                 &
     & sdotphi_int(:),                                                         &
     & rcar_edge(:,:), rflx_edge(:,:)

      LOGICAL          :: ALL_INSIDE_PLASMA
      CHARACTER(len=*), PARAMETER :: subname = 'TRACK_EDGES: '

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Warning and error flag and message
      iflag=0
      message=''

!......define the number of line segment end pts & the number of radial grid points
      n_seg = 2
      n_rho = 31

!........initialize n_int to 6x the number of radial grid pts...........!
      n_guess = 6 * n_rho
      n_int   = n_guess

!........define the input segment to be described in CARTESIAN coords (cylin is default)...!
      k_seg = 2

      ALLOCATE( rho(n_rho), irho_int(n_int), s_int(n_int) )
      ALLOCATE( r_seg(3, n_seg) )

      r_seg(:,1) = q0_vec
      r_seg(:,2) = qf_vec


!........fill rho_temp as an evenly spaced "half-grid", following Houlberg............!
      rho(1) = 0.0
      do i = 2,n_rho
        rho(i) = SQRT( (REAL(i) - 1.5) / REAL(n_rho - 1) )
      end do

!-------------------------------------------------------------------------------
!Allocate arrays for optional variables and call TRACK
!-------------------------------------------------------------------------------
!Dimensioned to allow up to 6 times as many intersections as surfaces
      ALLOCATE(rcar_int(1:3,1:n_guess),                                        &
     &         rflx_int(1:3,1:n_guess),                                        &
     &         sdotb_int(1:n_guess)  )

!.......find intersections btwn the beam path & the magnetic flux surfaces.......!
!.....(note that SDOTB_INT isn't really an optional variable for TRACK).........!
      CALL TRACK(n_rho,rho,n_seg,r_seg,n_int,irho_int,s_int,                   &
     &          iflag, message,                                                &
     &           K_SEG=k_seg,                                                  &
     &           RCAR_INT=rcar_int,                                            &
     &           RFLX_INT=rflx_int,                                            &
     &           SDOTB_INT=sdotb_int )

!      write(*,*) 'TRACK iflag in TRACK_EDGES = ', iflag
!      write(*,*) 'rcar_int AFTER TRACK in TRACK_EDGES = ', rcar_int
      write(*,*) subname, 'n_int AFTER TRACK = ', n_int
!      write(*,*) 's_int AFTER TRACK in TRACK_EDGES = ', s_int
!      write(*,*) '   '

      do i = 1, n_int
        write(*,*) subname, 'rflx_int rho = ', rflx_int(1, i)
      end do

      !Check messages
      IF(iflag /= 0) THEN

        write(*,*) 'non-zero iflag.  message = ', message
 
        IF(iflag > 0) GOTO 9999
        iflag=0
        message=''

      ENDIF

!............allocate boundary edge arrays & initialize to zero..........................!
!.........(note that n_int is a safe size, since n_edge <= n_int ).........!
      ALLOCATE( rflx_edge(3,n_int), rcar_edge(3,n_int) )
      ALLOCATE( edge_type(n_int), edge_index(n_int) )
      edge_type  = 0
      edge_index = 0
      rflx_edge  = 0.0
      rcar_edge  = 0.0


!.......decide if entire beam line (including endpts) is wholly inside the plasma....!
      ALL_INSIDE_PLASMA = .true.
      do i = 1,n_int
        if( rflx_int(1,i) .gt. 1.0 ) then
          ALL_INSIDE_PLASMA = .false.
          goto 100
        end if
      end do
 100  continue


!.........add the start point to the edge array.............................!
      rflx_edge(:,1) = rflx_int(:,1) 
      rcar_edge(:,1) = rcar_int(:,1)
           if( rflx_edge(1,1) <= 1.0 ) then 
             edge_type(1) = ENDPT_INSIDE_PLASMA
           else if( rflx_edge(1,1) > 1.0 ) then 
             edge_type(1) = ENDPT_OUTSIDE_PLASMA
           else
             write(*,*) 'WARNING.  Error in startpoint edge assigment!'
           end if  
      edge_index(1)  = 1


      if (ALL_INSIDE_PLASMA) then
!..........then startpt and endpt are only two edge array elements............!
!...........(so just assign the endpt to be the 2nd edge array element).......!
        n_edge = 2
        rflx_edge(:,n_edge) = rflx_int(:,n_int) 
        rcar_edge(:,n_edge) = rcar_int(:,n_int) 
        edge_type(n_edge)   = ENDPT_INSIDE_PLASMA
        edge_index(n_edge)  = n_int
      else

!............determine if the 1st line segment point is inside or outside the plasma.....!
        if( rflx_int(1,1) .le. 1.0 ) then
          beginpoint = INSIDE_PLASMA
        else
          beginpoint = OUTSIDE_PLASMA
        end if

!.............find all the plasma/vacuum boundary points...................!
        previous_pt = beginpoint
        n_edge = 1
        do i = 2,n_int
          if (previous_pt == INSIDE_PLASMA ) then
            if ( rflx_int(1,i) > 1.0 ) then
               currentpoint = OUTSIDE_PLASMA
!..............assign the PREVIOUS segment point to be an edge point.......!
              if ( i .ne. 2) n_edge = n_edge + 1
              rflx_edge(:,n_edge) = rflx_int(:,i-1) 
              rcar_edge(:,n_edge) = rcar_int(:,i-1) 
              edge_type(n_edge) = PLASMA_TO_VACUUM
              edge_index(n_edge) = i -1
              previous_pt = currentpoint
            end if        
          else if (previous_pt == OUTSIDE_PLASMA ) then
            if ( rflx_int(1,i) <= 1.0 ) then
              currentpoint = INSIDE_PLASMA
!..............assign the PREVIOUS segment point to be an edge point.......!
              if( i .ne. 2) n_edge = n_edge + 1
              rflx_edge(:,n_edge) = rflx_int(:,i-1) 
              rcar_edge(:,n_edge) = rcar_int(:,i-1) 
              edge_type(n_edge) = VACUUM_TO_PLASMA
              edge_index(n_edge) = i -1
              previous_pt = currentpoint
            end if                  
          else
            write(*,*) 'Error!  Ambiguous EDGE_POINT point assignment'
          end if  
        end do


!.........if the last element in the edge arrays isn't an actual plasma edgepoint, 
!..........add the beam endpoint to the edge arrays anyway....................!
        if( edge_index(n_edge) .NE.  n_int ) then
           n_edge = n_edge + 1
           rflx_edge(:,n_edge) = rflx_int(:,n_int) 
           rcar_edge(:,n_edge) = rcar_int(:,n_int)
           edge_index(n_edge) = n_int

           if( rflx_edge(1,n_edge) <= 1.0 ) then 
             edge_type(n_edge) = ENDPT_INSIDE_PLASMA
           else if( rflx_edge(1,n_edge) > 1.0 ) then 
             edge_type(n_edge) = ENDPT_OUTSIDE_PLASMA
           end if 
        end if         

      end if         ! end ALL_INSIDE_PLASMA if statement 

!      write(*,*) 'n_edge in TRACK_EDGES = ', n_edge
!      write(*,*) 'edge_type in TRACK_EDGES = ', edge_type
!      write(*,*) 'edge_index = ', edge_index
!
!      do i = 1, (n_edge + 1)
!        write(*,*) 'rcar_edge = ', rcar_edge(:,i)
!      end do

!........construct the global edge array for this beampath...................!
        call pv_edge_construct( myedge(i_beam), n_edge )

!.............transfer edge info to global arrays............................!
      myedge(i_beam) % n_edge = n_edge 
      myedge(i_beam) % xflx   = REAL(rflx_edge(:,1:n_edge), rprec) 
      myedge(i_beam) % xcar   = rcar_edge(:,1:n_edge) 
      myedge(i_beam) % etype  = edge_type(1:n_edge) 

!.........loop over edges to fill the edist array........................!
       do i= 1,n_edge
         myedge(i_beam)%edist(i) =                                             &
     &               DIST( intpol%ipbeam(i_beam)%q0vec,                       &
     &                     myedge(i_beam)%xcar(:,i) )              
       end do


 9999 CONTINUE



!..........DEBUGGING: verify that myedge filled correctly
      write(*,*) subname, 'rflx_edge(1,1) = ',  rflx_edge(1,1) 

      n_temp = myedge(i_beam)%n_edge
      write(*,*) subname, 'number of vac/int points = ', n_temp

      do i = 1, n_temp
        write(*,*) ' rho for edge ', i, '= ', myedge(i_beam)%xflx(1,i)
      end do


!.........deallocate local arrays.......................................!
      DEALLOCATE( rho, irho_int, s_int )
      DEALLOCATE( r_seg )
      DEALLOCATE(rcar_int, rflx_int, sdotb_int )

      DEALLOCATE( rflx_edge, rcar_edge )
      DEALLOCATE( edge_type, edge_index )

      return
      END SUBROUTINE TRACK_EDGES


!===================================================================================
      SUBROUTINE VMEC_EDGES(i_beam)
!--------------------------------------------------------------------------------
!
! Function:  Uses VMEC to "refine" the approximate vacuum/plasma interfaces found by TRACK
!            and saves values to derived TYPE variable "myedge"
!
!    global inputs : myedge
!    global outputs: myedge
!
! created 12/12/06 by J. Shields
!
!----------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
!      USE eq_T         ! module containing VMEC output info
      USE ip_beamline  ! module containing int_pol TYPE definition
      USE ip_global    ! module containing interfermeter/polarimeter global variables
      IMPLICIT NONE

!..............dummy variables.............................................!
      INTEGER, INTENT(in) :: i_beam     ! index of interferometer beam being integrated over

!............local variables.................................................!
      INTEGER           :: counter = 0
      INTEGER           :: i, j, nbeams, n_edge, iflag

      CHARACTER(len=80) :: mychar1
      CHARACTER(len=*), PARAMETER :: subname = 'VMEC_EDGES: '

      REAL(rprec), DIMENSION(3) ::   q0_vec, qf_vec, q_unit


      write(*,*) '===================================================='    
      write(*,*) ' '
      write(*,*) 'Hi.  Now in ', subname

      n_edge = myedge(i_beam)%n_edge
      write(*,*) subname, 'number of vac/int points = ', n_edge

      do i = 1, n_edge
        write(*,*) ' rho for edge ', i, '= ', myedge(i_beam)%xflx(1,i)

      end do


      RETURN
      END SUBROUTINE VMEC_EDGES



!=============================================================================
      REAL(KIND=rprec) FUNCTION E_DENSITY(q, q0_vec, qf_vec, q_unit,           &
     &                                    q_inout, xflux)
!-------------------------------------------------------------------------------------
!
! FUNCTION: determines the value of the electron density at a given distance "q"
!           along the beam path.
!
!
!  LOGIC:  Density is assumed to be maximum at the magnetic axis and falls off radially as a function
!          the VMEC radial flux coord "s".  The profile is defined as:
!                 ne = ne0 * ( 1 - s^tau)^kappa + ne_ambient,    where:
!
!                        ne0:  density at s = 0 (ie the max density)
!                 ne_ambient:  density of ambient medium outside plasma ( assumed << ne0)
!                tau, kappa:   constant scaling parameters
!
!   **Also note that there is a toggle to set a completely CONSTANT density profile
!      which may make life a bit easier for any future debugging.  JS 1/26/07
!
!
!  Created 5/12/05 by J. Shields
!   modified 7/25/07 JMS:  allowed density values to be takes from the model
!
!---------------------------------------------------------------------------------------
      USE stel_kinds, only : rprec, iprec, cprec
      USE stel_constants
      USE math_utilities
      USE v3f_global, ONLY: density_max_g, density_max_default_g,              &
     &      density_ambient_default_g, density_tau_default_g,                  &
     &      density_kappa_default_g, model_a
      IMPLICIT NONE

      REAL(rprec), INTENT(IN) :: q      ! distance along beam path
      REAL(rprec), DIMENSION(3), INTENT(IN) ::                                 &
     &            q0_vec                ! Cartesian position vector at beam start pt
      REAL(rprec), DIMENSION(3), INTENT(IN)  ::                                &
     &            qf_vec               !  Cartesian position vector at beam end pt
      REAL(rprec), DIMENSION(3), INTENT(IN) ::                                 &
     &            q_unit               ! Cartesian unit vector parallel to the beam path  
      REAL(rprec), DIMENSION(3), INTENT(IN) ::                                 &
     &            xflux                ! flux coord  position vector for point "q"

      INTEGER, INTENT(IN) :: q_inout  ! specifies if q is inside/outside of plasma


!.........local variables.....................................................!
      INTEGER :: iflag, ierr
      INTEGER :: ne_profile

      INTEGER, PARAMETER     :: INSIDE_PLASMA  = 1
      INTEGER, PARAMETER     :: OUTSIDE_PLASMA = 2
      INTEGER, PARAMETER     :: DENSITY_VALUES_FROM_MODEL = 0
      INTEGER, PARAMETER     :: HARD_WIRED_DEFAULT_VALUES = 1
      INTEGER, PARAMETER     :: CONSTANT_PLASMA_DENSITY   = 2

      REAL(rprec), PARAMETER :: ECHRG = 1.602176e-19
      REAL(rprec), PARAMETER :: ME = 9.1093826e-31
      REAL(rprec), PARAMETER :: EPZERO = 8.854188e-12
      REAL(rprec), PARAMETER :: C = 299792458.0
      
      REAL(rprec) :: ne0              ! max density (assumed to be at s = 0)
      REAL(rprec) :: ne_ambient       ! ambient density of medium surrounding plasma
      REAL(rprec) :: tau              ! exponent on the radial flux coord "s"
      REAL(rprec) :: kappa            ! exponent on the (1 - s^tau) term

      REAL(rprec) :: const, omega, nc
      REAL(rprec) :: ne, Bfield
      REAL(rprec) :: rho

 
      LOGICAL :: POINT_IN_VACUUM_REGION
      CHARACTER(len=80) :: message
      CHARACTER(len=*), PARAMETER :: subname = 'E_DENSITY: '

!      write(*,*) 'now entering', subname

!      ne_profile =  CONSTANT_PLASMA_DENSITY
!      ne_profile =  HARD_WIRED_DEFAULT_VALUES
      ne_profile =  DENSITY_VALUES_FROM_MODEL

      if( ne_profile == DENSITY_VALUES_FROM_MODEL ) then
        ne0        = model_a % density % ne_max
        ne_ambient = model_a % density % ne_ambient
        tau        = model_a % density % tau
        kappa      = model_a % density % kappa
      else if( ne_profile == HARD_WIRED_DEFAULT_VALUES) then
!.............use default values defined in module v3f_global.................!
        ne0        = density_max_default_g
        ne_ambient = density_ambient_default_g
        tau        = density_tau_default_g
        kappa      = density_kappa_default_g
      else if( ne_profile == CONSTANT_PLASMA_DENSITY) then
        ne = density_max_g
      else
        write(*,*) subname, 'ERROR in setting density profile!!!'
      end if


!.......set rho equal to the first (radial) coord of the flux coord vector.....!
      rho = xflux(1)



      if ( q_inout == OUTSIDE_PLASMA ) then
        POINT_IN_VACUUM_REGION = .true.
      else if ( q_inout == INSIDE_PLASMA ) then
        POINT_IN_VACUUM_REGION = .false.
      else
        write(*,*) ' error in E_DENSITY density determination!'
      end if

!.........compute the density at this particular point along the path.............! 
      if( POINT_IN_VACUUM_REGION ) then
        ne = ne_ambient
      else 
        ne = ne0 * ( 1.0 - rho**tau )**kappa + ne_ambient
      end if

      E_DENSITY = ne

      return
      END FUNCTION E_DENSITY


!=============================================================================
      SUBROUTINE PLASMA_INOUT(q, q_cyl, q_inout, ierr)
!--------------------------------------------------------------------------
!
! FUNCTION: determines whether point "q" is inside (1) or outside (2) the plasma
!
!
!  Created 6/30/05 by J. Shields
!
! **IMPORTANT CAVEAT: logic currently assumes that the beamline is a STRAIGHT LINE
!
!-------------------------------------------------------------------------
      USE stel_kinds, only : rprec, iprec, cprec
      USE ip_global
      USE b_transform_vmec, ONLY: INVALID_RHO_VMEC
      IMPLICIT NONE

!...........input variables..................................................!
      REAL(rprec), INTENT(IN) :: q       ! step point along beam path
      REAL(rprec), INTENT(IN) :: q_cyl(3) ! Cylin coords of step point "q"
      INTEGER, INTENT(OUT)    :: q_inout ! tells whether step point is INSIDE (1)
                                         ! or OUTSIDE (2) the plasma region
      INTEGER, INTENT(OUT)    :: ierr    ! Error flag. 0 = good, 1= error, -1 = warning

!.........local variables.....................................................!
      INTEGER, PARAMETER :: INSIDE_PLASMA  = 1
      INTEGER, PARAMETER :: OUTSIDE_PLASMA = 2
      INTEGER, PARAMETER :: LOCATION_UNKNOWN = -99

      INTEGER, PARAMETER :: VMEC_CONVERGENCE = 0
      INTEGER, PARAMETER :: TRACK_EDGE_INFO  = 1

      INTEGER :: iclose, n_edge
      INTEGER :: rubric_select

      REAL(rprec) :: qside, qdist


      REAL(rprec), ALLOCATABLE, DIMENSION(:)  :: qsubtract

      CHARACTER(len=80) :: message
      CHARACTER(len=*), PARAMETER :: subname = 'PLASMA_INOUT: '

!      write(*,*) 'HI.  Now in ', subname

!............set toggle to use VMEC or TRACK to decide if point is in/outside of plasma......!
!      rubric_select = TRACK_EDGE_INFO
      rubric_select = VMEC_CONVERGENCE


      if ( rubric_select == VMEC_CONVERGENCE) then
!        write(*,*) subname, 'VMEC toggle selected!'
        if (INVALID_RHO_VMEC(q_cyl) ) then
          q_inout = OUTSIDE_PLASMA
        else
          q_inout = INSIDE_PLASMA
        end if

      else if ( rubric_select == TRACK_EDGE_INFO) then
        call PLASMA_INOUT_TRACK(q, q_inout, ierr)
      else
        write(*,*) subname, 'ERROR in rubric_select ASSIGNMENT!!!!'
      end if

      return
      END SUBROUTINE PLASMA_INOUT


!=============================================================================
      SUBROUTINE PLASMA_INOUT_TRACK(q, q_inout, ierr)
!--------------------------------------------------------------------------
!
! FUNCTION: determines whether point "q" is inside (1) or outside (2) the plasma
!
!
!  Created 6/30/05 by J. Shields
!
! **IMPORTANT CAVEATS: 
!      1)  logic currently assumes that the beamline is a STRAIGHT LINE
!
!      2) This routine is NOT currently utilized.  However, if one does wish
!          to utilize it in conjuction with signal_model_compute_idsl() [as opposed
!          to the top-level driver routine FARADAY_ROTATION() that it was originally
!          designed for], one will need to implement
!          a routine that maps the input "signal" to the beamline index "ibeam"
!          that it corresponds to.  JS 6/15/07
!
!-------------------------------------------------------------------------
      USE stel_kinds, only : rprec, iprec, cprec
      USE ip_global
      IMPLICIT NONE

!...........input variables..................................................!
      REAL(rprec), INTENT(IN) :: q       ! step point along beam path

      INTEGER, INTENT(OUT)    :: q_inout ! tells whether step point is INSIDE (1)
                                         ! or OUTSIDE (2) the plasma region
      INTEGER, INTENT(OUT)    :: ierr    ! Error flag. 0 = good, 1= error, -1 = warning

!.........local variables.....................................................!
      INTEGER, PARAMETER :: INSIDE_PLASMA  = 1
      INTEGER, PARAMETER :: OUTSIDE_PLASMA = 2
      INTEGER, PARAMETER :: LOCATION_UNKNOWN = -99

      INTEGER, PARAMETER :: LEFT_OF_EDGE  = 3
      INTEGER, PARAMETER :: RIGHT_OF_EDGE = 4

      INTEGER, PARAMETER :: VACUUM_TO_PLASMA = 1
      INTEGER, PARAMETER :: PLASMA_TO_VACUUM = 2
      INTEGER, PARAMETER :: ENDPT_INSIDE_PLASMA   = -11
      INTEGER, PARAMETER :: ENDPT_OUTSIDE_PLASMA  = -21

      INTEGER :: iclose, n_edge
      INTEGER :: etype_q          ! edge type at pt "q" (vac/plas  or plas/vac)

      REAL(rprec) :: qside, qdist
      REAL(rprec), ALLOCATABLE, DIMENSION(:)  :: qsubtract

      CHARACTER(len=80) :: message
      CHARACTER(len=*), PARAMETER :: subname = 'PLASMA_INOUT_TRACK: '

!      write(*,*) 'HI.  Now in ', subname
      ierr = 0
      n_edge = myedge(ibeam) % n_edge
      ALLOCATE( qsubtract(n_edge) )

!..........find the edge point closest to q......................................!
      qsubtract = ABS( myedge(ibeam)%edist - q )
      iclose = TRANSFER( MINLOC(qsubtract), iclose )

      etype_q = myedge(ibeam) % etype(iclose)

!........determine if q is before or after the edgepoint along the beampath.......!
      qdist =  myedge(ibeam)%edist(iclose) - q

      if ( qdist <=  0.0 ) then
         qside = LEFT_OF_EDGE
      else if ( qdist >  0.0 ) then
         qside = RIGHT_OF_EDGE
      else
         write(*,*) 'error in qside assignment!'
      end if


!............determine whether q is inside or outside the plasma..................!
      if( qside == LEFT_OF_EDGE) then
        select case(etype_q)
        case(VACUUM_TO_PLASMA)
          q_inout = OUTSIDE_PLASMA
        case(PLASMA_TO_VACUUM)
          q_inout = INSIDE_PLASMA
        case(ENDPT_INSIDE_PLASMA)
          q_inout = INSIDE_PLASMA
        case(ENDPT_OUTSIDE_PLASMA)
          q_inout = OUTSIDE_PLASMA
        case default
          ierr = 1
          q_inout = LOCATION_UNKNOWN
        end select

      else if( qside == RIGHT_OF_EDGE) then
        select case(etype_q)
        case(VACUUM_TO_PLASMA)
          q_inout = INSIDE_PLASMA
        case(PLASMA_TO_VACUUM)
          q_inout = OUTSIDE_PLASMA
        case(ENDPT_INSIDE_PLASMA)
          q_inout = INSIDE_PLASMA
        case(ENDPT_OUTSIDE_PLASMA)
          q_inout = OUTSIDE_PLASMA
        case default
          ierr = 1
          q_inout = LOCATION_UNKNOWN
        end select

      else 
        ierr = 1
        q_inout = LOCATION_UNKNOWN
        write(*,*) 'error in plasma q_inout assignment!'
      end if
  
      return
      END SUBROUTINE PLASMA_INOUT_TRACK

!===========================================================================
      SUBROUTINE B_CARTESIAN( q, xcar, Bcar, q_inout, BCAR_EXT,                &
     &                        BCAR_INT, B_CYL)
!--------------------------------------------------------------------------------
!
! Function:  Determines the Cartesian components of the net magnetic field at a point
!            specified by the Cartesian position vector "xcar".  As an option, the
!            individual contributions of the external field coils and the internal
!            B field of the plasma itself can also be returned.
!
! 
! Created 5/15/05 by J. Shields
! Updated 11/22/05 (J.Shields) to use the B/I netCDF files created by module 
!         intpol_dot in V3RFUN. 
!
!
! CAVEAT : there is currently NO safeguard against the case in which the path point
!         is outside the field vessel and hence the database contains no info about
!         the external field at that point.  This type of error will generally
!          cause the module to crash.  JS 7/1/05
!
!----------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
!      USE track_mod
      USE eq_T         ! module containing VMEC output info
      USE v3f_global   ! module containing VMEC output info (i.e. it SAVEs model_a)
      USE b_transform_vmec, ONLY: TRANS_CAR2CYL
      USE b_transform_mod, ONLY: COMPUTE_B_CAR
      IMPLICIT NONE

!.........input variables....................................................!
      INTEGER, INTENT(IN)      :: q_inout ! specifies if q is in/outside of plasma

      REAL(rprec), INTENT(IN)  :: q     ! distance along path
      REAL(rprec), DIMENSION(3), INTENT(IN) ::                                 &
     &            xcar                 ! Cartesian coordinate POSITION vector
!      REAL(rprec), DIMENSION(3), INTENT(INOUT) ::                                 &
!     &            xflux                 ! flux coordinate POSITION vector
      REAL(rprec), DIMENSION(3), INTENT(OUT) ::                                &
     &            Bcar             ! Cartesian NET magnetic field vector at point xcar
      REAL(rprec), DIMENSION(3),  INTENT(OUT), OPTIONAL ::                     &
     &            BCAR_EXT         ! Cartesian EXTERNAL B field vector at point xcar
      REAL(rprec), DIMENSION(3),  INTENT(OUT), OPTIONAL ::                     &
     &            BCAR_INT         ! Cartesian INTERNAL B field vector at point xcar
      REAL(rprec), DIMENSION(3), INTENT(OUT), OPTIONAL ::                       &
     &            B_CYL            ! NET magnetic field vector at point xcar in Cylin coords

!...........local variables....................................................!
      INTEGER :: iflag
      INTEGER, PARAMETER :: INSIDE_PLASMA  = 1
      INTEGER, PARAMETER :: OUTSIDE_PLASMA = 2
      INTEGER, PARAMETER :: LOCATION_UNKNOWN = -99

      REAL(rprec), DIMENSION(3) ::                                              &
     &            x_cyl            ! Cylinderical coordinate POSITION vector
      REAL(rprec), DIMENSION(6) ::                                              &
     &            g_cyl            ! Cylinderical coordinate derivatives
      REAL(rprec), DIMENSION(3) ::                                             &
     &            Bcar_external    ! Cartesian EXTERNAL B field  vector at point xcar
      REAL(rprec), DIMENSION(3) ::                                             &
     &            Bcar_internal      ! Cartesian INTERNAL B vector at point xcar
      REAL(rprec), DIMENSION(3) ::                                             &
     &            B_cylin      ! net B field in cylin coords

      CHARACTER(len=80) :: message
      CHARACTER(len=*), PARAMETER :: subname = 'B_CARTESIAN:'

      LOGICAL :: COMPUTE_EXTERNAL_FIELD

!      write(*,*) 'Hi.  Now in', subname


!      if( PRESENT(BCAR_EXT) .OR. PRESENT(BCAR_INT) .OR.                        &
!     &    (q_inout .NE. INSIDE_PLASMA) .OR. (xflux(1) >= 1.0)    ) then
      if( PRESENT(BCAR_EXT) .OR. PRESENT(BCAR_INT) .OR.                        &
     &    (q_inout .NE. INSIDE_PLASMA)                 ) then

        COMPUTE_EXTERNAL_FIELD = .TRUE.

!.....obtain Cartesian B field vector from the ip_beam derived TYPE info....!
         call CARTESIAN_B_EXTERNAL(q, Bcar_external)

      else
         COMPUTE_EXTERNAL_FIELD = .FALSE.
      end if


      if( q_inout == INSIDE_PLASMA ) then 

!...........transform Cartesian  xcar() vector into cylinderical coords.....!
!        call AJAX_CAR2CYL(xcar, x_cyl)
        call TRANS_CAR2CYL(xcar, x_cyl)

!!...........transform Cylinderical vector x_cyl() into flux coords.............!
!        call AJAX_CYL2FLX(x_cyl, xflux, iflag, message, G_CYL = g_cyl)
!
!        if (iflag .ne. 0 ) then
!            write(*,*) 'WARNING! CYL2FLX iflag in B_CARTESIAN = ', iflag
!            write(*,*) 'error message is ', message
!        end if

!..........compute the components of the B field (in Cartesian coords) ...............!
!         call AJAX_B(xflux, x_cyl, g_cyl, iflag, message,                         &
!     &               B_CAR = Bcar )

         call COMPUTE_B_CAR(x_cyl, Bcar, iflag, B_CYLIN=B_cylin)

        if (iflag .ne. 0 ) then
            write(*,*) 'WARNING! iflag for VMEC call = ', iflag
!            write(*,*) 'error message is ', message
        end if
      else
         Bcar = Bcar_external
      end if   ! end if for q_inout INSIDE_PLASMA


!      write(*,*) subname, 'NET B Field = ', Bcar

      if (COMPUTE_EXTERNAL_FIELD) then
!...........calc the internal plasma B field from B_net and B_external................!
        Bcar_internal = Bcar - Bcar_external

        if( PRESENT(BCAR_EXT) ) BCAR_EXT = Bcar_external
        if( PRESENT(BCAR_INT) ) BCAR_INT = Bcar_internal
        if( PRESENT(B_CYL) )    B_CYL    = B_cylin

      end if  ! end if for COMPUTE_EXTERNAL_FIELD

 
      RETURN
      END SUBROUTINE B_CARTESIAN

!===============================================================================
      SUBROUTINE CARTESIAN_B_EXTERNAL(q, B_ext)
!--------------------------------------------------------------------------------
!
! Function:  Determines the Cartesian components of the net magnetic field due to
!            the EXTERNAL FIELD COILS ONLY at a point "q" along the path..
!  
!  Created 11/22/05 by J. Shields
!
!----------------------------------------------------------------------------------
      USE stel_kinds
      USE eq_T         ! module containing VMEC output info
      USE v3f_global   ! module containing VMEC output info (i.e. it SAVEs model_a)
      USE ip_global    ! module containing interfermeter/polarimeter global variables
      IMPLICIT NONE


!.........input variables....................................................!
      REAL(rprec), INTENT(IN)  ::  q   ! distance along the beam path  
      REAL(rprec), DIMENSION(3), INTENT(OUT) ::                                &
     &            B_ext          ! Cartesian NET magnetic field vector at point xcar

!...........local variables....................................................!
      INTEGER :: i,j, k, iflag
      INTEGER :: npathpts, ncoil
      INTEGER :: ipath, ipath_hi, ipath_low
      INTEGER :: icoil
      REAL(rprec)    :: qdist, q0
      REAL(rprec)    :: q_hi, q_low
      REAL(rprec)    :: delta_q
      REAL(rprec), DIMENSION(3) :: DIM
      REAL(rprec), DIMENSION(3) :: B_hi, B_low
      REAL(rprec), DIMENSION(3) :: B_temp, B_sum
      CHARACTER(len=*), PARAMETER :: subname = 'CARTESIAN_B_EXTERNAL: '

!      write(*,*) 'Hi.  Now in ', subname
      

!......Use the size of the 3D B_ratio array to determine # of grid points.............!
      DIM(1:3) = SHAPE(ipbeam_g%B_ratio)
      npathpts = DIM(1)
      ncoil = DIM(2)

!........compute dist btwn grid points, noting that i=1 is at q0 and......!
!........i= npathpts is at qf [so that n_segments = (npathpts -1) ].........! 
      q0    = ipbeam_g%q0
      qdist = ipbeam_g%qdist
      delta_q = qdist / REAL(npathpts - 1)
      q_low = REAL(npathpts - 1)

!............Find ipath to the left & right of point "q"........!
      ipath_low = AINT( q / delta_q) + 1
      ipath_hi  = ipath_low + 1
!      write(*,*) subname, 'ipath_low = ', ipath_low
!      write(*,*) subname, 'ipath_hi = ', ipath_hi 

!.......calculate the "q" value for points at ipath_hi, ipath_low....!
      q_low = q0 + REAL(ipath_low - 1) * delta_q
      q_hi  = q0 + REAL(ipath_hi - 1) * delta_q
 
!.......Find the B field at the path grid pts just above/below point "q"........!
      call PATHPOINT_B(ipath_low, B_low)
      call PATHPOINT_B(ipath_hi, B_hi)

!.....linearly interpolate btwn B_hi & B_low to compute B field at point "q"....!
      B_ext = ( (q_hi - q)*B_low + (q - q_low)*B_hi ) / delta_q

!      write(*,*) subname, 'B_low = ', B_low
!      write(*,*) subname, 'B_ext = ', B_ext
!      write(*,*) subname, 'B_hi = ', B_hi

      RETURN
      END SUBROUTINE CARTESIAN_B_EXTERNAL



!=================================================================================
      SUBROUTINE INTEGRAND_B_FACTOR(x_cyl, k_hat_car, B_car, B_cyl,            &
     &                              B_factor)
!--------------------------------------------------------------------------------
!
! Function:  determines the B field dependence for the Faraday Rotation (w1) and the 
!            the Cotton-Mouton (w1,w2) integrands.
!
! LOGIC: The beam propagation is assumed to travel perpendicular to the toroidal (phi) direction.  This
!        allows the creation of a very simple local coordinate system that has its "z" direction
!        defined in terms of the propagation vector (k_hat) and its "y" vector defined in terms of
!        the toroidal direction vector (phi_hat).  Specifically:
!
!        x3 = k_hat,  x2 = -1.0*phi_hat,   so then:  x1 = x2 .cross. x3 = -1.0 * (phi_hat .cross. k_hat)
!
!   As a result, the B field in this local coordinate system is given by:
!
!        B3 = B .dot. x3_hat =  B .dot. k_hat 
!        B2 = B .dot. x2_hat =  B .dot. (-phi_hat) = -1.0 * B_phi_cylin             
!        B1 = B .dot. x1_hat =  B .dot. ( -phi_hat .cross. k_hat) = ( B .cross. k_hat) .dot. phi_hat
!
!
! IMPORTANT NOTE: The above assumptions imply that the B dependence for w1 and w2 [ie B_factor(1)
!            and B_factor(2) ]  for microwave beams 
!           that are aimed at OBLIQUE angles through the plasma (ie angles that are NOT perdicular to phi_hat)
!            will be computed *INCORRECTLY* !!!!!!!!
!              ( However, the Faraday rotation B dependence [B_factor(3)] *will* remain valid )
!
! Created 3/6/06 by J. Shields
!
!
!----------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
!      USE b_transform_vmec, ONLY: TRANS_CAR2CYL
      USE b_transform_mod, ONLY: B_CAR_TO_CYL
!      USE v3f_global   ! module containing VMEC output info (i.e. it SAVEs model_a)
      IMPLICIT NONE

!.........input variables....................................................!

      REAL(rprec), DIMENSION(3), INTENT(IN) ::  x_cyl       ! current position vector (GLOBAL cylin)
      REAL(rprec), DIMENSION(3), INTENT(IN) ::  k_hat_car   ! beam propagation unit vector (GLOBAL Cartesian)
      REAL(rprec), DIMENSION(3), INTENT(IN) ::  B_car       ! B field in GLOBAL Cartesian coords
      REAL(rprec), DIMENSION(3), INTENT(IN) ::  B_cyl       ! B field in GLOBAL cylin coords

      REAL(rprec), DIMENSION(3), INTENT(OUT) ::  B_factor   ! B dependence of respective (w1,w2,w3) integrands
 


!...........local variables....................................................!
      INTEGER :: i
      REAL(rprec)    :: temp
      REAL(rprec), DIMENSION(3)   :: B_local   ! B field vector, expressed in LOCAL Cartesian coords

      REAL(rprec), DIMENSION(3)   :: car_B_cross_k   ! B .cross. k, expressed in GLOBAL Cartesian coords
      REAL(rprec), DIMENSION(3)   :: cyl_B_cross_k  ! B .cross. k, expressed in GLOBAL cylin coords

      CHARACTER(len=*), PARAMETER :: subname = 'INTEGRAND_B_FACTOR:'

!      write(*,*) 'Hi.  Now in', subname

!........compute the cross product of B and k in global cartesian coords..............!

!............compute B_cross_k in global Cartesian coords.............................!
      call CROSS_PRODUCT(B_car, k_hat_car, car_B_cross_k)

!..........convert Cartesian B_cross_k vector into cylin coords..................!
!.......***we need the B_CAR_TO_CYL routine
!      call TRANS_CAR2CYL(car_b_cross_k, cyl_B_cross_k)
      call B_CAR_TO_CYL(x_cyl, car_B_cross_k, cyl_B_cross_k)

!.........set B1 to the phi component of B_cross_k ................................! 
      B_local(1) = cyl_B_cross_k(2)


!..........set B2 = -B_phi..........................................................!
      B_local(2) = -1.0_rprec * B_cyl(2)


!........set B3 to the component of the B field parallel to the beam path.............................!
      B_local(3) = DOT_PRODUCT( k_hat_car, B_car)

!      write(*,*) subname, 'B_local = ', B_local

!.............now compute the integrand factors from the components of the local B field vector......!
      B_factor(1) = B_local(1)**2.0 - B_local(2)**2.0 
      B_factor(2) = B_local(1) * B_local(2)
      B_factor(3) = B_local(3)

 
      RETURN
      END SUBROUTINE INTEGRAND_B_FACTOR


!===============================================================================
      SUBROUTINE PATHPOINT_B(ipath, B_path)
!--------------------------------------------------------------------------------
!
! Function:  Determines the Cartesian components of the magnetic field due
!            to the combined EXTERNAL FIELD COILS ONLY at a point
!            that is "ipath" steps along the path (where the step size was determined
!            during the allocation of the B_ratio array in program V3RFUN).
!  
!  Created 11/22/05 by J. Shields
!
!----------------------------------------------------------------------------------
      USE stel_kinds
      USE eq_T         ! module containing VMEC output info
      USE v3f_global   ! module containing VMEC output info (i.e. it SAVEs model_a)
      USE ip_global    ! module containing interfermeter/polarimeter global variables
      IMPLICIT NONE


!.........input variables....................................................!
      INTEGER, INTENT(IN)  ::  ipath   ! index of point along beam path
      REAL(rprec), DIMENSION(3), INTENT(OUT) ::                                &
     &            B_path          ! Cartesian magnetic field vector at point ipath

!...........local variables....................................................!
      INTEGER :: i,j, k, iflag
      INTEGER :: npath, ncoil
      INTEGER :: nextcur
      REAL(rprec) :: current, temp
      REAL(rprec), DIMENSION(3) :: DIM
      REAL(rprec), DIMENSION(3) :: B_over_I
      REAL(rprec), DIMENSION(3) :: B_sum
      REAL(rprec), DIMENSION(3) :: B_temp
      CHARACTER(len=*), PARAMETER :: subname = 'PATHPOINT_B: '
      CHARACTER(len=*), PARAMETER ::                                           &
     &    err_msg = 'WARNING: nextcur, B_ratio array size MISMATCH!!!!!'

!      write(*,*) 'Hi.  Now in ', subname
      
!......determine the size of the 3D B_ratio array.............!
      DIM(1:3) = SHAPE(ipbeam_g%B_ratio)
      npath = DIM(1)
      ncoil = DIM(2)
!      write(*,*) subname, 'B_ratio ncoil = ', ncoil


!........debug test to verify that we can access model_a correctly.  JS 12/1/05
      nextcur = model_a%eqstate%fixp%nextcur
!      write(*,*) subname, 'nextcur = ', nextcur
!      write(*,*) subname, 'extcur = ', model_a%eqstate%varp%extcur(1)

!      if ( nextcur .ne. ncoil) then
!        write(*,*) subname, err_msg
!      end if

!.......sum the B field components over all coil collections........!
      B_sum = 0.0
      do j = 1,ncoil
        B_over_I = ipbeam_g%B_ratio(ipath,j,:)
        current  = model_a%eqstate%varp%extcur(j)
        B_temp   = (B_over_I) * current
        B_sum    = B_sum + B_temp

!        write(*,*) subname, 'B over I = ', B_over_I
!        write(*,*) subname, 'extcur = ', current
!        write(*,*) subname, 'B_temp = ', B_temp
!        write(*,*) subname,' B_sum = ', B_sum
      end do

      B_path = B_sum
!      write(*,*) subname, 'Output B_path = ', B_path

      RETURN
      END SUBROUTINE PATHPOINT_B


!=============================================================================
      SUBROUTINE TRACK_FLX2CYL(xflux, x_cyl, iflag, message, X_CAR)
!-------------------------------------------------------------------------------
! FUNCTION: uses the TRACK module to convert the input flux coordinate vector 
!           "xflux" from flux coordinates to cylinderical Coordinates.  
!              There is also an option to output a Cartesian coord vector as well.
!
! NOTE: this program is a more robust version of the AJAX subroutine AJAX_FLUX2CYL(), 
!      since the the AJAX transformation routines are only valid for vectors with a
!      radial flux coordinate ("rho")  < 1.0 .
!
! **WARNING TO USER: Track uses AJAX to do the transformations, and AJAX uses a slightly
!                    *different* definition of flux coords (ie there's an extra factor
!                     of rho^m in the AJAX definitions).  Therefore, ONLY use this
!                     routine if you really, REALLY mean it!  JS 12/15/06
!
!-------------------------------------------------------------------------------
      USE SPEC_KIND_MOD
      USE TRACK_MOD
      USE AJAX_MOD
      IMPLICIT NONE


!declaration of input variables

      REAL(rprec), DIMENSION(3), INTENT(IN)::                                  &
     &            xflux                ! flux coordinate vector
      REAL(rprec), DIMENSION(3), INTENT(OUT) ::                                &
     &            x_cyl                ! Cylindrical coordinate vector
      REAL(rprec), DIMENSION(3),  INTENT(OUT), OPTIONAL ::                     &
     &            X_CAR        ! Cartesian coordinate vector 

      INTEGER, INTENT(OUT) :: iflag
      CHARACTER(len=*), INTENT(OUT) :: message

!Declaration of local variables

      INTEGER ::                                                               &
     & k_equil,k_seg,                                                          &
     & n_rho, n_int, n_seg

      INTEGER ::                                                               &
     & i, j,                                                                   & 
     & n_guess

      INTEGER, ALLOCATABLE ::                                                  &
     & irho_int(:)

      REAL(KIND=rspec), DIMENSION(3) ::                                        &
     & xflux_end
      REAL(KIND=rspec), ALLOCATABLE ::                                         &
     & rho(:), r_seg(:,:), s_int(:), q_int(:)
      REAL(KIND=rspec), ALLOCATABLE ::                                         &
     & rcar_int(:,:),rcyl_int(:,:),rflx_int(:,:),sdotb_int(:)


!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Warning and error flag and message
      iflag=0
      message=''

!      write(*,*) ' now in TRACK_FLX2CYL!'

!......define the number of line segment end pts & the number of radial grid points
      n_seg = 2
      n_rho = 10

!........initialize n_int to 3x the number of radial grid pts...........!
      n_guess = 3 * n_rho
      n_int   = n_guess

!........define the input segment to be described in FLUX coords................! 
!........(  k_seg = 1 => FLUX,  k_seg = 2  => CARTESIAN,  k_seg = other => CYLIN ) ...!
      k_seg = 1

      ALLOCATE( rho(n_rho), irho_int(n_int), s_int(n_int) )
      ALLOCATE( r_seg(3, n_seg) )


!........set the line segment endpoints very close together.........! 
      xflux_end(1)   = xflux(1) - 0.0001
      xflux_end(2:3) = xflux(2:3)

      r_seg(:,1) = xflux
      r_seg(:,2) = xflux_end


!........fill rho_temp as an evenly spaced "half-grid", following Houlberg............!
      rho(1) = 0.0
      do i = 2,n_rho
        rho(i) = SQRT( (REAL(i) - 1.5) / REAL(n_rho - 1) )
      end do


!-------------------------------------------------------------------------------
!Allocate arrays and call TRACK
!-------------------------------------------------------------------------------
!Dimensioned to allow up to 3 times as many intersections as surfaces
      ALLOCATE(rcyl_int(1:3,1:n_guess),                                         &
     &         rflx_int(1:3,1:n_guess),                                         &
     &         sdotb_int(1:n_guess) )

      if( PRESENT(X_CAR) ) ALLOCATE( rcar_int(1:3,1:n_guess) )

      CALL TRACK(n_rho,rho,n_seg,r_seg,n_int,irho_int,s_int,                   &
     &          iflag, message,                                                &
     &           K_SEG=k_seg,                                                  &
     &           RCAR_INT=rcar_int,                                            &
     &           RCYL_INT=rcyl_int,                                            &
     &           RFLX_INT=rflx_int,                                            &
     &           SDOTB_INT=sdotb_int )

!.......set x_cyl & X_CAR to the values of the first TRACK pt............!
      x_cyl = rcyl_int(:,1)
      if( PRESENT(X_CAR) ) X_CAR = rcar_int(:,1)

      !Check messages
      IF(iflag /= 0) THEN

        write(*,*) 'non-zero iflag.  message = ', message
 
        IF(iflag > 0) GOTO 9999
        iflag=0
        message=''
      ENDIF


 9999 CONTINUE

!.........deallocate local arrays.......................................!
      DEALLOCATE( rho, irho_int, s_int, r_seg) 
      DEALLOCATE( rcyl_int, rflx_int, sdotb_int )
      if( PRESENT(X_CAR) ) DEALLOCATE( rcar_int)

      return
      END SUBROUTINE TRACK_FLX2CYL


!=============================================================================
      SUBROUTINE TRACK_TRANS2FLX(x_vec, xflux, iflag, message, K_COORD)
!-------------------------------------------------------------------------------
! FUNCTION: uses the TRACK module to convert the input coordinate vector 
!           "x_vec" from its original coordinates (default = cylinderical) to 
!           flux coordinates.  
!              There is an optional "toggle" variable (K_COORD) that can be set to 
!           select the x_vec coordinates ( K_COORD = 2  => Cartesians coords,
!           K_COORD = other [or not present]  => cylin coords )
!
! NOTE: this program is used in place of subroutine AJAX_CYL2FLX(), since the
!       the AJAX transformation routines are only valid for vectors with a radial
!        flux coordinate ("rho")  < 1.0 .
!
! **WARNING TO USER: Track uses AJAX to do the transformations, and AJAX uses a slightly
!                    *different* definition of flux coords (ie there's an extra factor
!                     of rho^m in the AJAX definitions).  Therefore, ONLY use this
!                     routine if you really, REALLY mean it!  JS 12/15/06
!
!-------------------------------------------------------------------------------
      USE SPEC_KIND_MOD
      USE TRACK_MOD
      USE AJAX_MOD
      IMPLICIT NONE


!declaration of input variables

      REAL(rprec), DIMENSION(3), INTENT(IN) ::                                &
     &            x_vec               ! coordinate vector (cylin or Cartesian,
                                      !       based on k_coord)
      REAL(rprec), DIMENSION(3), INTENT(OUT)::                                  &
     &            xflux                ! flux coordinate vector

      INTEGER, INTENT(IN), OPTIONAL :: K_COORD  ! K_COORD = 2  => Cartesian x_vec
                                                ! otherwise    => cylin x_vec

      INTEGER, INTENT(OUT) :: iflag              ! error flag. iflag = 0 => ok
      CHARACTER(len=*), INTENT(OUT) :: message   ! error message

!Declaration of local variables

      INTEGER ::                                                               &
     & k_equil,k_seg,                                                          &
     & n_rho, n_int, n_seg

      INTEGER ::                                                               &
     & i, j,                                                                   & 
     & n_guess

      INTEGER, ALLOCATABLE ::                                                  &
     & irho_int(:)

      REAL(rprec), DIMENSION(3)  ::                                            &
     &            x_cyl        ! Cylinderical coordinate vector                   
      REAL(rprec), DIMENSION(3)  ::                                            &
     &            x_car        ! Cartesian    coordinate vector
      REAL(KIND=rspec), DIMENSION(3) ::                                        &
     & x_vec_end
      REAL(KIND=rspec), ALLOCATABLE ::                                         &
     & rho(:), r_seg(:,:), s_int(:), q_int(:)
      REAL(KIND=rspec), ALLOCATABLE ::                                         &
     & rcar_int(:,:),rcyl_int(:,:),rflx_int(:,:),sdotb_int(:)


!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Warning and error flag and message
      iflag=0
      message=''

!      write(*,*) ' now in TRACK_TRANS2FLX!'

!......define the number of line segment end pts & the number of radial grid points
      n_seg = 2
      n_rho = 10

!........initialize n_int to 3 the number of radial grid pts...........!
      n_guess = 3 * n_rho
      n_int   = n_guess

!........input vector x_vec should be in either CYLIN or CARTESIAN coords.........! 
!...........default for this subroutine is k_seg = 0. [cylin coords] )...........!
!........(  k_seg = 1 => FLUX,  k_seg = 2  => CARTESIAN,  k_seg = other => CYLIN ) ...!

      if( PRESENT(K_COORD) ) then
        SELECT CASE( K_COORD)
        CASE (2)
          k_seg = 2     ! x_vec in Cartesian coords
        CASE DEFAULT
          k_seg = 0     ! x_vec in Cylinderical coords
        END SELECT
      else
        k_seg = 0       ! x_vec in Cylinderical coords
      end if

      ALLOCATE( rho(n_rho), irho_int(n_int), s_int(n_int) )
      ALLOCATE( r_seg(3, n_seg) )


!........set the line segment endpoints very close together.........! 
      x_vec_end(1)   = x_vec(1) - 0.0001
      x_vec_end(2:3) = x_vec(2:3)

      r_seg(:,1) = x_vec
      r_seg(:,2) = x_vec_end


!........fill rho_temp as an evenly spaced "half-grid", following Houlberg............!
      rho(1) = 0.0
      do i = 2,n_rho
        rho(i) = SQRT( (REAL(i) - 1.5) / REAL(n_rho - 1) )
      end do


!-------------------------------------------------------------------------------
!Allocate arrays and call TRACK
!-------------------------------------------------------------------------------
!Dimensioned to allow up to 3 times as many intersections as surfaces
      ALLOCATE(rcyl_int(1:3,1:n_guess),                                         &
     &         rflx_int(1:3,1:n_guess),                                         &
     &         sdotb_int(1:n_guess) )

      ALLOCATE( rcar_int(1:3,1:n_guess) )

      CALL TRACK(n_rho,rho,n_seg,r_seg,n_int,irho_int,s_int,                   &
     &          iflag, message,                                                &
     &           K_SEG=k_seg,                                                  &
     &           RCAR_INT=rcar_int,                                            &
     &           RCYL_INT=rcyl_int,                                            &
     &           RFLX_INT=rflx_int,                                            &
     &           SDOTB_INT=sdotb_int )

!.......set x_cyl & X_CAR to the values of the first TRACK pt............!
      x_cyl = rcyl_int(:,1)
      x_car = rcar_int(:,1)
      xflux = rflx_int(:,1)

      !Check messages
      IF(iflag /= 0) THEN

        write(*,*) 'non-zero iflag.  message = ', message
 
        IF(iflag > 0) GOTO 9999
        iflag=0
        message=''
      ENDIF

 9999 CONTINUE

!.........deallocate local arrays.......................................!
      DEALLOCATE( rho, irho_int, s_int, r_seg) 
      DEALLOCATE( rcyl_int, rflx_int, sdotb_int )
      DEALLOCATE( rcar_int)

      return
      END SUBROUTINE TRACK_TRANS2FLX

!=================================================================
      SUBROUTINE SIMPLE_TRAPEZOID(func, a,b,s,nsteps)
!------------------------------------------------------------------------ 
! FUNCTION: Numerical Recipes subroutine to compute an integral via the trapazoidal
!           rule using a simplistic fixed-number-of-steps loops.  This routine is 
!           designed primarily for debugging.  The true Numerical Recipes routine
!           of interest (which totally supercedes this one), qtrap, uses a variable
!           number of iteration steps.
!
!----------------------------------------------------------------------
      USE stel_kinds, only : rprec, iprec, cprec
      IMPLICIT NONE

!Declaration of namelist input variables
      CHARACTER(len=25) ::   mychar

      INTEGER ::                                                               &
     & i, j, k, nsteps

      REAL(KIND=rprec) :: a, b, s
      REAL(KIND=rprec) :: func
      EXTERNAL :: func

      write(*,*) 'Now in SIMPLE TRAPEZOID!'

      do j = 1,(nsteps+1)
        call trapzd(func,a,b,s,j)
      end do

      END SUBROUTINE SIMPLE_TRAPEZOID


!==============================================================================
      SUBROUTINE trapzd(func, a,b, s,n)
!-----------------------------------------------------------------------------
!
! FUNCTION: Numerical Recipes subroutine to compute the integral of function
!           "func" via the trapazoidal method.  Specifically, "this routine
!           computes the nth stage of refinement of an extended trapazoidal rule.
!           Func is input as the name of the function to be integrated between
!           limits a and b, also input.  When called with n = 1, the routine
!           returns as the crudest estimate of the definite integral of f(x)
!           over the interval  x = {a,b}.  Subsquent calls with n = 2,3,4,...
!          (in that order) will improve the accuracy of the result (s) by adding
!           2^(n-2) additional interior points.  s should not be modified between
!           subsequent calls. 
!
! NOTE:  this subroutine is from Numerical Receipes
!
! ***WARNING: be careful that the input value of n does not cause the large integer
!             "it" to exceed the machine MAX INTEGER!!
! 
!-----------------------------------------------------------------------------
      USE stel_kinds, only : rprec, iprec, cprec
      IMPLICIT NONE

      INTEGER :: n
      REAL(KIND=rprec)    :: a,b,s, func
      EXTERNAL :: func

      INTEGER :: it, j
      REAL(KIND=rprec) :: del, sum, tnm, x
      
      if( n == 1) then
        s = 0.5 * (b-a) * ( func(a) * func(b) )
      else
        it = 2**(n-2)
        if (it .gt. 10**9) then
           write(*,*) 'WARNING.  Variable it may exceed MAX INTEGER ',          &
     &                 'for this program! PROGRAM TERMINATING.'
           write(*,*) 'it = ', it
           stop
        end if
!        write(*,*) 'it = ', it
        tnm = it
        del = (b-a)/tnm                ! this is the spacing of the pts to be added
        x = a + 0.5*del
        sum = 0.0
        do j = 1,it
          sum = sum + func(x)
          x = x + del
        end do

        s = 0.5 * ( s + (b-a)*sum/tnm )  ! this replaces s by its refined value
      end if

      return
      END SUBROUTINE trapzd


!========================================================================
      SUBROUTINE qtrap(func,a,b,s)
!-----------------------------------------------------------------
!
!  FUNCTION: Returns s as the integral of the function func from a to b. The
!            parameters EPS can be set to the desired fractional accuracy
!            and JMAX so that 2^(JMAX+1) is the maximum allowed number of steps.
!            Integration is performed by the trapazoidal rule.
!
! NOTES:  1) this subroutine is from Numerical Receipes
!         2) this subroutine also requires the trapzd() subroutine from NR 
!-----------------------------------------------------------------------
      USE stel_kinds, only : rprec, iprec, cprec
      IMPLICIT NONE

      INTEGER :: JMAX
      REAL(KIND=rprec) :: a,b,func, s, EPS
      EXTERNAL :: func
!      PARAMETER (EPS=1.0e-6, JMAX = 20)
!      PARAMETER (EPS=1.0e-2, JMAX = 20)
      PARAMETER (EPS=1.0e-4, JMAX = 20)
      INTEGER :: j
      REAL(KIND=rprec) :: olds

      olds = -1.0e30

      do j=1,JMAX
!        write(*,*) 'counter in qtrap = ', j
        call trapzd( func,a,b,s,j)
        if( j .gt. 5) then
          if( ABS(s-olds) .lt. EPS*ABS(olds) .OR.                              &
     &        (s .eq. 0.0 .AND. olds .eq. 0.0 ) ) return
        end if
        olds = s 
      end do

      write(*,*) 'too many steps in qtrap.  Consider relaxing tolerance'
!      pause 'too many steps in qtrap. Consider relaxing tolerance'
      END SUBROUTINE qtrap

!=============================================================================
      REAL(KIND=rprec) FUNCTION test_func(x)
      USE stel_kinds, only : rprec, iprec, cprec
      IMPLICIT NONE
      REAL(KIND=rprec) :: x

      test_func = x**2.0

      return
      END FUNCTION test_func


!===================================================================================
      SUBROUTINE PV_EDGE_TEST(i_beam)
!--------------------------------------------------------------------------------
!
! Function:  tests pv_edge derived TYPE
!
!    global inputs : myedge
!    global outputs: myedge
!
! created 12/12/06 by J. Shields
!
!  **here's another interesting wrinkle:  "ibeam" is a GLOBAL integer from ip_global
!       (which is why our local version is i_beam)
!
!----------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
!      USE eq_T         ! module containing VMEC output info
      USE ip_beamline  ! module containing int_pol TYPE definition
      USE ip_global    ! module containing interfermeter/polarimeter global variables
      IMPLICIT NONE

!..............dummy variables.............................................!
      INTEGER, INTENT(in) :: i_beam     ! index of interferometer beam being integrated over

!............local variables.................................................!
      INTEGER           :: counter = 0
      INTEGER           :: i, j, nbeams, n_edge, iflag
      REAL(rprec), DIMENSION(3) ::   q0_vec, qf_vec, q_unit

      TYPE(pv_edge) :: pvtemp
      CHARACTER(len=80) :: mychar1
      CHARACTER(len=*), PARAMETER :: subname = 'PV_EDGE_TEST: '


    

      write(*,*) '===================================================='    
      write(*,*) ' '
      write(*,*) 'Hi.  Now in ', subname

!      n_edge = myedge(ibeam)%n_edge
!      write(*,*) subname, 'number of vac/int points = ', n_edge
!
!      do i = 1, n_edge
!        write(*,*) ' rho for edge ', i, '= ', myedge(ibeam)%xflx(1,i)
!      end do


      RETURN
      END SUBROUTINE PV_EDGE_TEST
!===================================================================================
      SUBROUTINE IPSL_CDF_TEST()
!--------------------------------------------------------------------------------
!
! Function:  tests ipsl_T and ipsl_cdf derived TYPES
!
!
! created 6/15/07 by J. Shields
!
!----------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
!      USE eq_T         ! module containing VMEC output info
      USE ip_beamline  ! module containing int_pol TYPE definition
      USE ipsl_T
      USE ipsl_cdf
      USE ip_global    ! module containing interfermeter/polarimeter global variables
      USE ezcdf
      IMPLICIT NONE


!............local variables.................................................!
      INTEGER           :: iou_cdf
      INTEGER           :: i, j, nbeams, istat, iflag
      REAL(rprec), DIMENSION(3) ::   q0_vec, qf_vec, q_unit

      TYPE(ipsl_desc) :: my_ipsl, ipsl_in, temp_ipsl
      CHARACTER(len=*), PARAMETER :: cdf_outfile = 'ipsl_tempout.nc'
      CHARACTER(len=*), PARAMETER :: subname = 'IPSL_CDF_TEST: '


    

      write(*,*) '===================================================='    
      write(*,*) ' '
      write(*,*) 'Hi.  Now in ', subname

      write(*,*) 'input qdist = ', intpol%ipbeam(1)%qdist

      my_ipsl%ipbeam = intpol % ipbeam(1)
      write(*,*) 'IPSL qdist = ', my_ipsl%ipbeam%qdist

      my_ipsl%l_name = 'TESTING ipsl l_name'
      my_ipsl%s_name = 'TESTING ipsl s_name'
      my_ipsl%ip_lname = 'IP long name IPSL test output'
      my_ipsl%ip_sname = 'IP short name'
      my_ipsl%units = 'furlongs'


      write(*,*) 'my_ipsl REGULAR short name = ', my_ipsl%s_name
      write(*,*) 'my_ipsl REGULAR long name = ', my_ipsl%l_name
      write(*,*) 'my_ipsl IP long name = ', my_ipsl%ip_lname
      write(*,*) 'my_ipsl IP short name = ', my_ipsl%ip_sname
      write(*,*) 'my_ipsl units = ', my_ipsl%units


!..........open file to output FR angle, etc to...................................!
!........open the input .nc file and read in the variables..............!

      call cdf_open(iou_cdf, cdf_outfile, 'w', istat)
!      write(*,*) 'iou_cdf = ', iou_cdf, ' istat = ', istat
      CALL assert_eq(0,istat,' JS CDF_OPEN called failed')

!      call ipsl_cdf_define_desc( my_ipsl, iou_cdf)
!      call ipsl_cdf_write_desc( my_ipsl, iou_cdf)
      call ipsl_cdf_define( my_ipsl, iou_cdf)
      call ipsl_cdf_write( my_ipsl, iou_cdf)
      call cdf_close(iou_cdf)


      call cdf_open(iou_cdf, cdf_outfile, 'r', istat)
!      call ipsl_cdf_read_desc(ipsl_in, iou_cdf)
      call ipsl_cdf_read(ipsl_in, iou_cdf)
      call cdf_close(iou_cdf)
 
      write(*,*) 'read-in qdist = ', ipsl_in%ipbeam%qdist
      write(*,*) 'read-in REGULAR short name = ', TRIM(ipsl_in%s_name)
      write(*,*) 'read-in REGULAR long name = ', TRIM(ipsl_in%l_name)
      write(*,*) 'read-in IP long name = ', TRIM(ipsl_in%ip_lname)
      write(*,*) 'read-in IP short name = ', TRIM(ipsl_in%ip_sname)
      write(*,*) 'read-in units = ', ipsl_in%units


      write(*,*) 'Now exiting', subname
      RETURN
      END SUBROUTINE IPSL_CDF_TEST


!===================================================================================
      SUBROUTINE IPSL_INTEGRAL(ipsl_in, line_integral)
!--------------------------------------------------------------------------------
!
! Function:  computes the line integral (ie either the interferometery phase shift
!             or the Faraday rotation angle) for a given input "ipsl diagnostic".  It
!             is essentially an interface between modules signal_model & faraday_mod.
!
! INPUTS: ipsl_in
!
! OUPUTS:  line_integral   ! line integral for either interferometery or polarimetry
!
! global outputs: ipbeam_g, integrand_toggle
!
! 
! created 6/15/07 by J. Shields
!
! Called by signal_model_compute_idsl() in module signal_model
!
!----------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
      USE ip_beamline  ! module containing int_pol TYPE definition
      USE ipsl_T
      USE ip_global    ! module containing interfermeter/polarimeter global variables
      IMPLICIT NONE


!............dummy variables.................................................!
      TYPE(ipsl_desc), INTENT(in) :: ipsl_in
      REAL(rprec),  INTENT(out) :: line_integral

!............local variables.................................................!

      INTEGER            :: i, j, nbeams, istat, iflag
      INTEGER, PARAMETER :: FARADAY_ROTATION = 1  ! Polarimetry phase angle
      INTEGER, PARAMETER :: PHASE_SHIFT      = 2  ! Interferometry phase angle

      REAL(rprec)                ::   a, b          ! lower, upper bounds of integral
      CHARACTER(len=*), PARAMETER :: subname = 'IPSL_INTEGRAL: '

!      write(*,*) '===================================================='    
!      write(*,*) ' '
!      write(*,*) 'Hi.  Now in ', subname



!..........use ipsl_type to determine which line integral to compute...........!
      SELECT CASE (TRIM(ADJUSTL(ipsl_in % ipsl_type)))
      CASE ('interf')
        write(*,*) subname, 'INTERF case SELECTED!!!'
        integrand_toggle = PHASE_SHIFT  
      CASE ('polarim')
        write(*,*) subname, 'POLARIM case SELECTED!!!'
        integrand_toggle = FARADAY_ROTATION  
      CASE DEFAULT
         write(*,*) subname, 'ERROR!! unrecognized s_type: ',                  &
     &               TRIM(ipsl_in % ipsl_type)
      END SELECT 


!..........use global "ipbeam_g" to convey ipsl info to FARADAY_FUNC().....!
!      write(*,*) 'ipsl_in qdist = ', ipsl_in%ipbeam%qdist
      ipbeam_g = ipsl_in % ipbeam
!      write(*,*) 'ipbeam_g qdist = ', ipbeam_g%qdist

!..........set lower/upper bounds of integration..........................!
      a = ipbeam_g % q0
      b = a + ipbeam_g % qdist

!...........use the trapazoid rule to compute the line_integral..........!
      call qtrap(FARADAY_FUNC,a,b, line_integral)

      RETURN
      END SUBROUTINE IPSL_INTEGRAL

!========================================================================
      SUBROUTINE stepsum(func, a, b, nsteps, s)
!-----------------------------------------------------------------
!
!  FUNCTION: Computes the 1-D line integral of an input function "func".
!
!  LOGIC: This is the simplest (and slowest) integration method that one can imagine.  It simply
!         computes func(x) at each grid point & sums the results.  Accuracy improves with increased
!         values of "nsteps".
!
! created 6/6/07  by  J. Shields
!
!-----------------------------------------------------------------------
      USE stel_kinds, only : rprec, iprec, cprec
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nsteps    ! number of steps along the grid of integration
      REAL(rprec), INTENT(IN)    :: a,b       ! lower, upper bounds of integral
      REAL(rprec)                :: func      ! input function name
      REAL(rprec), INTENT(OUT)   :: s         ! output integral
      EXTERNAL                   :: func

!............local variables...........................
      INTEGER :: i, j
      REAL(rprec)    :: stepsize             ! infinitesmal stepsize along path
      REAL(rprec)    :: x                    ! abscissa value at each point
      REAL(rprec)    :: sum

      if ( nsteps == 0 ) then
        s = 0.0
      else
        stepsize = ABS(b - a) / REAL(nsteps, rprec)
        sum = 0.0
        do i = 1, nsteps
          x = a + REAL(i-1)*stepsize
          sum = sum + func(x)
        end do
        s = sum
      end if

      END SUBROUTINE stepsum

!===================================================================================
      SUBROUTINE construct_frcm_filename(frcm_outfile)
!--------------------------------------------------------------------------------
!
! Function:  Constructs the output file for faraday rotation & cotton-mouton
!            computations, using the global input file ipsl_cdf_infile
!
!  LOGIC:  Assumes an input filename of form: (blah)ip_v3rfun_(stuff).nc
!          and outputs a file of form: frcm_(stuff).out
!
! INPUTS: (none)
!
! OUTPUTS: frcm_outfile
!
! global inputs: ipsl_cdf_infile
! global outputs: (none)
!
! created 6/19/07 by J. Shields
!
!----------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
      USE v3f_global, ONLY: intpol_filename
      IMPLICIT NONE

!............passed argument..............................................!
      CHARACTER(len=*)           :: frcm_outfile

!............local variables.................................................!
      INTEGER     :: i, j, slen
      INTEGER    :: ipstr
      CHARACTER(len=300)           :: tempstr, filestr
      CHARACTER(len=*), PARAMETER :: subname =                                 &
     &                                 'construct_frcm_filename: '

!      write(*,*) 'Hi.  Now in ', subname


!.......left justify original filename (that was output by V3RFUN)................!
      tempstr = ADJUSTL(TRIM(intpol_filename))

!........strip off the prefix & suffix from the orginal filename.............!
      slen = LEN_TRIM(tempstr)
      ipstr = INDEX(TRIM(tempstr),'ip_v3rfun_',back=.true.)
      filestr = tempstr( (ipstr+10):(slen - 3))



!.........concatenate the final filename.....................!
      frcm_outfile = 'frcm_'//TRIM(filestr)//'.out'

      write(*,*) subname, 'frcm_outfile = ', frcm_outfile

      RETURN
      END SUBROUTINE construct_frcm_filename


      END MODULE faraday_mod
