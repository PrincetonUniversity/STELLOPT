
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! module sxr_routines
!
! These are routines/functions/subprograms to be used by the
! Soft X-Ray chord implementation on V3FIT
!
! initial coding Greg Hartwell 7/24/09
! last modified 1/05/10
!
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      MODULE sxr_routines
      
      USE stel_kinds, only : rprec
      USE stel_constants, only : pi, zero
      USE safe_open_mod
      USE v3f_global
      USE v3_utilities
      USE sxrc_T
      USE read_wout_mod
      USE vmec_utils
      
      IMPLICIT NONE
!-------------------------------------------------------------------------------
! GLOBAL DEFINITIONS will go here
!-------------------------------------------------------------------------------
   
      PRIVATE pi, rprec, zero
      PRIVATE iou_sxr_nli
      PRIVATE sxr_int_to_string
      PRIVATE are_equal

!-------------------------------------------------------------------------------
! Derived Type Declarations
!
! sxrc_desc   - structure to hold chord information
!
!
! Subroutine and Function Definitions
!
! sxr_int_to_string             - specialized function to create a 3-digit string
!                                    from a 0-999 number
! init_sxr_chords               - reads namelist input file holding camera info
!                                    and creates chords
! are_equal                     - determines if two reals are equal
!                                    to a specified precision
! emissivity_model              - model of the emissivity
!
! cyl2flx_test                  - subroutine to test the functionality of the
!                                 cyl2flx and flx2cyl routine
!                                 Not normally needed
!
!-------------------------------------------------------------------------------

!===============================================================================

      CONTAINS
!===============================================================================
! sxr_int_to_string
! This is a PRIVATE function
!
! takes a 1-3 digit integer and converts it
! to a 3 digit string
! Used by write_sxr_chord_diagnostic routine
!-------------------------------------------------------------------------------
      FUNCTION sxr_int_to_string(chord_num)
!-------------------------------------------------------------------------------
! ARGUMENT
! INPUT
! chord - a number from 0-999
!-------------------------------------------------------------------------------
         INTEGER, INTENT(IN) :: chord_num   ! input chord number
         CHARACTER*3 sxr_int_to_string      ! output string
         CHARACTER*4 temp_string        ! temporary 4 digit string
         
         temp_string=''
         write(temp_string,'(i4)')1000+chord_num
         sxr_int_to_string=temp_string(2:4)
         
      END FUNCTION sxr_int_to_string
!-------------------------------------------------------------------------------

!===============================================================================
! are_equal
! This is a PRIVATE function
! takes two reals and a precision
! determines if they are equal to within given precision
!-------------------------------------------------------------------------------
      FUNCTION are_equal(num1,num2,prec)
!-------------------------------------------------------------------------------
! ARGUMENT
! INPUT
! num1 - real number
! num2 - real number
! prec - number of decimal places for equality test
!-------------------------------------------------------------------------------
         REAL(rprec), INTENT(IN) :: num1,num2
         INTEGER, INTENT(IN)     :: prec
         LOGICAL are_equal
         
         are_equal=.FALSE.
         IF(prec.LT.0) THEN
            WRITE(*,*)'error - precision in ARE_EQUAL < 0'
         ELSE
            are_equal=FLOOR(num1*10**prec).EQ.FLOOR(num2*10**prec)
         ENDIF
         
      END FUNCTION are_equal
!-------------------------------------------------------------------------------
!===============================================================================
! emissivity_model
! This is a PRIVATE function
! the emissivity function
! f(r)=n^3/2 * P^1/2
!-------------------------------------------------------------------------------
      FUNCTION emissivity_model(s,model_in)
!-------------------------------------------------------------------------------
! ARGUMENT
! INPUT
! s         - VMEC radial-like flux surface parameter
! model_in  - from model_T
!-------------------------------------------------------------------------------
      REAL(rprec), INTENT(IN) :: s
      TYPE (model),INTENT(IN) :: model_in
      REAL(rprec)             :: emissivity_model
      
      REAL(rprec) :: P
      REAL(rprec) :: n
      INTEGER     :: j
      
      IF (s.LE.1.0) THEN
        n=1                                                                     &
     &          * (1-s**model_in % density % tau)                               &
     &          **model_in % density % kappa                                    &
     &          + model_in % density % ne_ambient
        P=0.0
        DO j=1,SIZE(model_in % eqstate % varp %am)
          P=P+(model_in % eqstate % varp % am(j))*(s**(j-1))
        ENDDO
        emissivity_model=n**1.5+SQRT(P)
      ELSE
        emissivity_model=0.0
      ENDIF
      
      END FUNCTION emissivity_model
!-------------------------------------------------------------------------------

!===============================================================================
! init_sxr_chords
!
! reads namelist/sxr_camera_input/ input from sxr_camera Input File
! creates the sxr chords
! returns the sxr_chord info into the cameras structure
! last modified 8/27/09 GJH
!------------------------------------------------------------------------------

      SUBROUTINE init_sxr_chords(cam_file,sxr_chords)
!------------------------------------------------------------------------------
! ARGUMENTS
! INPUT
! cam_file - name of the sxr_camera file
! OUTPUT
! sxr_chords - array of sxr chords
! alloc_chords - number of chords allocated (this may not be needed)
!------------------------------------------------------------------------------
         CHARACTER (len=*) ,INTENT(IN) :: cam_file
         TYPE (sxrc_desc),DIMENSION(:),                                        &
     &                         ALLOCATABLE,INTENT(INOUT):: sxr_chords
         	
!------------------------------------------------------------------------------
! LOCAL VARIABLES
! cam          - camera counter
! chord        - chord counter
! alloc_chords - number of chords allocated
! tot_chords   - number of chords created
! istat        - status of safe_open call
! ier1         - error
! sub_name     - subroutine name
! lassert      - logical test
!------------------------------------------------------------------------------

         INTEGER cam         
         INTEGER chord        			 
         INTEGER :: tot_chords = 0
         INTEGER :: istat = 0
         INTEGER :: ier1
         CHARACTER*30,PARAMETER :: sub_name='init_sxr_chords'
         LOGICAL :: lassert
         INTEGER :: alloc_chords     

!------------------------------------------------------------------------------
! NAMELIST VARIABLES
! camera specifications read from namelist file
!
! camera_type      - type of camera
! camera_name      - name of particular camera of camera_type
! position_units   - units defining camera positions
! measurement_units- units defining diagnostic units
! n_chords         - number of chords in camera
! n_cameras
! r_major_sxr      - major radius
! camera_r         - minor radius of camera
! camara_t         - theta (poloidal) angle of camera
! camera_p         - phi (toroidal) angle of camera
! diode_d1         - dimension 1 (width) of diode
! diode_d2         - dimension 2 (length) of diode
! diode_sep        - separation between diodes
! focal_len        - distance between diode array and camera slit
! slit_d1          - dimension 1 (width) of slit
! slit_d2          - dimension 2 (length) of slit
! cal_constants    - array of calibration constants for diode array
! sigmas           - error estimates for each diode in array
! enable           - .TRUE. (default) use diode in calculations
!------------------------------------------------------------------------------

        CHARACTER(LEN=cam_type_len),                                          &
     &                DIMENSION(max_cameras)::camera_type=''
        CHARACTER(LEN=cam_name_len),                                          &
     &                DIMENSION(max_cameras)::camera_name=''
        CHARACTER(LEN=unit_len),                                              &
     &                DIMENSION(max_cameras)  :: position_units=''
        CHARACTER(LEN=unit_len),                                              &
     &                DIMENSION(max_cameras)  :: measurement_units=''
        INTEGER, DIMENSION(max_cameras)      :: n_chords=0
        INTEGER                              :: n_cameras=0
        REAL (rprec), DIMENSION(max_cameras) :: camera_r = zero
        REAL (rprec), DIMENSION(max_cameras) :: camera_t = zero
        REAL (rprec), DIMENSION(max_cameras) :: camera_p = zero
        REAL (rprec), DIMENSION(max_cameras) :: diode_d1 = zero
        REAL (rprec), DIMENSION(max_cameras) :: diode_d2 = zero
        REAL (rprec), DIMENSION(max_cameras) :: diode_sep = zero
        REAL (rprec), DIMENSION(max_cameras) :: focal_len = zero
        REAL (rprec), DIMENSION(max_cameras) :: slit_d1 = zero
        REAL (rprec), DIMENSION(max_cameras) :: slit_d2 = zero
        REAL (rprec), DIMENSION(max_chords,max_cameras)                       &
     &                                       :: cal_constants = zero
        REAL (rprec), DIMENSION(max_chords,max_cameras)                       &
     &                                       :: sigmas = zero
        LOGICAL, DIMENSION(max_chords,max_cameras)                            &
     &                                       :: enable=.TRUE.
        REAL (rprec)						 :: r_major = 0.75 


         NAMELIST /sxr_camera_input/ n_cameras,camera_type,camera_name,       &
     &      n_chords,camera_r,camera_t,camera_p,diode_d1,diode_d2,            &
     &      diode_sep,focal_len,slit_d1, slit_d2, cal_constants,              &
     &      sigmas,enable,position_units,measurement_units,r_major

!------------------------------------------------------------------------------
! start of executable code
!------------------------------------------------------------------------------
         alloc_chords=0
         
         CALL safe_open(iou_sxr_nli,istat,cam_file,'old','formatted')
         IF(istat.ne.0) THEN
            WRITE(*,*) 'problem opening file in---init_sxr_chords'
         END IF
!         WRITE(*,*)"opened the sxr_camera file"
         READ(iou_sxr_nli,NML=sxr_camera_input)
!         WRITE(*,*)"there are ",n_cameras,"cameras"

         IF(n_cameras .GT. max_cameras) THEN
            WRITE(*,*)"ERROR IN SUBROUTINE init_sxr_chords-------------"
            WRITE(*,*)"NUMBER OF CAMERAS GREATER THAN MAXIMUM ALLOWED"
            WRITE(*,*)"INCREASE MAX_CAMERAS IN SXR_ROUTINES MODULE"
         END IF
         
         CLOSE(iou_sxr_nli)
         
         !calculate number of chords to allocate
         DO cam=1,n_cameras
            DO chord=1,n_chords(cam)
               IF (enable(chord,cam)) THEN
                  alloc_chords=alloc_chords+1
               END IF  
            END DO
         END DO
        
         !allocate and initialize sxr_chords
         IF (alloc_chords .GE. 1) THEN
            ALLOCATE(sxr_chords(alloc_chords),STAT=ier1)
            CALL assert_eq(0,ier1,sub_name // 'allocate chords error')
            DO chord=1,alloc_chords
              CALL sxrc_desc_destroy(sxr_chords(chord))
            END DO
         ELSE
            WRITE(*,*)"INIT_SXR_CHORDS - NO CHORDS TO ALLOCATE"
         END IF
         
         ! next create the chords
         IF (n_cameras .GE. 1) THEN
           DO cam=1,n_cameras
             DO chord=1,n_chords(cam)
               IF (enable(chord,cam)) THEN
                  tot_chords=tot_chords+1
                  CALL sxrc_desc_construct(sxr_chords(tot_chords),             &
     &                     camera_type(cam),camera_name(cam),                  &
     &                     camera_r(cam),camera_t(cam),camera_p(cam),          &
     &                     diode_d1(cam),diode_d2(cam),diode_sep(cam),         &
     &                     focal_len(cam),slit_d1(cam),slit_d2(cam),           &
     &                     cal_constants(chord,cam),sigmas(chord,cam),         &
     &                     n_chords(cam),chord,r_major,                        &
     &                     position_units(cam),measurement_units(cam))                        
          
                  sxr_chords(tot_chords) % chord_name =                        &
     &                              'sxr'//sxr_int_to_string(tot_chords)
                END IF  
             END DO
           END DO
         ELSE
            WRITE(*,*)"ERROR IN SUBROUTINE init_sxr_chords"
            WRITE(*,*)"--------NO CAMERAS DEFINED-------"
         ENDIF
         
         lassert=alloc_chords.EQ.tot_chords
         CALL assert(lassert,sub_name //'chord allocation check failed')
         
         
      END SUBROUTINE init_sxr_chords

!-------------------------------------------------------------------------------
! Line integration
!-------------------------------------------------------------------------------
!   USE vmec_utils
!       from LIBSTELL/Sources/Modules
!       add to *.dep files
!
!    CALL cyl2flx(rzl_in, r_cyl, c_flx, ns_in, ntor_in, mpol_in, 
!     1      ntmax_in, lthreed_in, lasym_in, info, nfe, fmin) 
!     1     IGNORE THE OPTIONAL: mscale, nscale, ru, zu, rv, zv 
!
!     LOCAL PARAMETERS:
!     ftol    :   nominally, set to 1.E-16. Gives a maximum (relative)
!                 error in matching R and Z of sqrt(ftol), or 1.E-8. 
!                 To increase accuracy, ftol should be lowered, but this 
!                 may require more Newton iterations (slows code).
!       
!     INPUT:
!     rzl_in  :   4D array with r,z (lambda) Fourier coefficients vs. radius
!                 rzl_in(ns_in,0:ntor_in,0:mpol_in-1,2*ntmax_in)
!                 
!     r_cyl   :   vector specifying cylindrical point to match, R = r_cyl(1), 
!                 N*phi = r_cyl(2), Z = r_cyl(3)
!                 NOTE: N*phi (N=no. field periods) is input, NOT phi!
!     ns_in   :   number of radial nodes in input array rzl_in
!     ntmax_in:   number of different types of modes (cos-cos, sin-sin, etc)
!     mpol_in :   number of poloidal modes (0:mpol_in-1)
!     ntor_in :   number of toroidal modes = ntor_in+1 (0:ntor)
!     lthreed_in :true if this is a 3D plasma
!     lasym_in:   true if this is an asymmetric plasma
!     mscale  (nscale) : 
!                 optional scaling arrays for cos, sin arrays. Used
!                 only if this routine is called from within VMEC.
!
!     OUTPUT:
!     nfe     :   number of function evaluations
!     info    :   = 0, no errors, -1, fmin > ftol (tolerance exceeded on output)
!                 = -3, s > 1 outside plasma, probably
!     fmin    :   minimum value of f = (r - Rin)**2 + (z - Zin)**2 at c_flx

!     INPUT/OUTPUT:
!     c_flx   :   array of flux coordinates (s = c_flx(1), u=theta= c_flx(2), 
!                 v = N*phi= c_flx(3))
!                 on input, initial guess (recommend magnetic axis if "cold" start)
!                 on output, s, u values corresponding to r_cyl
!
!     need something like
!     CALL signal_mc_model_compute(s_data,s_desc(isig),model_in)
! perhaps
!     CALL signal_sc_model_compute(signal,s_desc(isig),model_in)
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
       
!===============================================================================       
! sxr_signal
! computes the signal detected by a sxr diode
! last modified 1/5/10  GJH
!------------------------------------------------------------------------------- 
 
       SUBROUTINE sxr_signal(sxr_chords,model_in)       
!-------------------------------------------------------------------------------
! ARGUMENTS
! sxr_chords - array of sxr chord information
! model_in   - VMEC model with equilibrium and density information
!-------------------------------------------------------------------------------       
     
       TYPE (model),INTENT(IN) :: model_in
       TYPE (sxrc_desc),INTENT(IN),                                            &
     &                  DIMENSION(:),ALLOCATABLE :: sxr_chords
       
! cyl2flx definitions
       REAL(rprec), DIMENSION(3) :: r_cyl(1:3)=0.0
       REAL(rprec), DIMENSION(3) :: r_cyl2(1:3)=0.0
       REAL(rprec), DIMENSION(3) :: c_flux(1:3)=0.0
       INTEGER :: nfe 
       INTEGER :: info,iflag,points
       REAL(rprec) :: fmin

! other specifications
       
       CHARACTER*300 :: diag_file='diagtest2'
       INTEGER       :: iounit=6
       INTEGER       :: ichord
       INTEGER       :: ncurrents
       INTEGER       :: n_chords
       INTEGER       :: ierr,istat
       CHARACTER     :: woutfile='_.nc'
       REAL(rprec)   :: dR,dZ,dPhi,dL       ! integration step sizes
       INTEGER       :: ipoints          ! number of integration points
       INTEGER       :: istep
       REAL(rprec),DIMENSION(:), ALLOCATABLE :: s_array
       
!-------------------------------------------------------------------------------
!  START OF EXECUTABLE CODE
!-------------------------------------------------------------------------------
       ! need sxr chord information --------------------- done
       ! need the VMEC flux surface parameters -----------done with wout
       ! need a model for the emissivity------------------done, but most basic
       !        incorporate this into model_T?
       !        base on Franz/Rice papers
       ! take into account the distance from the detector (1/r^2)
       ! take into account the area/volume of the position in question
       ! integrate emissivity over the chord
       !        the cyl2flx routine returns s. If s>1 then we are outside
       !        the last closed flux surface and the contribution will be zero.

!-------------------------------------------------------------------------------
!      make sure the sxr_chord information is available
!------------------------------------------------------------------------------- 
       
        n_chords=size(sxr_chords)
        WRITE(*,*) 'number of SXR chords in sxr_signal',n_chords
!-------------------------------------------------------------------------------
!       flux surface parameters should be available from the equilibrium model
!       at the moment (1/8/10  GJH)  they are not, but can be read
!       from the wout file.  The variable and fixed parameters ARE available
!       from the equilibrium model
!-------------------------------------------------------------------------------
!       ncurrents=model_in % eqstate % fixp % nextcur
!       WRITE(*,*)'surfaces',model_in % eqstate % aux1 % ns
!       WRITE(*,*)'currents',                                                   &
!     &       model_in % eqstate % varp % extcur(1:ncurrents)
!       WRITE(*,*)'iotaf',model_in % eqstate % aux1 % iotaf
!       WRITE(*,*)'code',model_in % eqstate % code
!       WRITE(*,*)'code',model_in % eqstate % version
       WRITE(*,*)'aux1 ldef',model_in % eqstate % l_def_aux1

       WRITE(*,*)'density model -- nemax    ambient   tau    kappa'
       WRITE(*,10)model_in % density % ne_max,                                  &
     &            model_in % density % ne_ambient,                              &
     &      model_in % density % tau,model_in % density % kappa
10     FORMAT('                ',E8.2,3('  ',f6.2))
       WRITE(*,*)'START OF WOUT INPUTS'
       CALL readw_and_open(woutfile,ierr)
       IF(ierr .ne. 0 ) write(*,*) 'read_wout error'
!-------------------------------------------------------------------------------       
! This section is just testing the wout read capability.  Wout reads will be
! used until the eq_state % aux1(2) variables are accessible  
!-------------------------------------------------------------------------------            
       WRITE(*,*)'mgrid file ',TRIM(mgrid_file)
       WRITE(*,*)'ns=',ns,'mnmax=',mnmax
       WRITE(*,*)'ntor=',ntor,'mpol=',mpol,'ntmax=',ntmax
       WRITE(*,*)'lthreed=',lthreed,'lasym=',lasym
       WRITE(*,*)size(rmnc)
       WRITE(*,*)'nfp=',nfp
       WRITE(*,*)'aspect ratio=',aspect

       !create inputs for cyl2flx call
       !create rzl_in
       IF (ALLOCATED(rzl_local)) THEN 
            WRITE(*,*)'RZL already allocated'
       ELSE
            CALL LoadRZL
       ENDIF
       IF (ALLOCATED(rzl_local)) THEN 
           WRITE(*,*)'RZL now allocated'
       ELSE
           WRITE(*,*)'RZL allocation error - returning from sxr_signal'
           RETURN
       ENDIF
       
       WRITE(*,*)'--------------SIGNAL STUFF--------------------'
       ichord=10
       WRITE(*,*)'looking at chord ',ichord
       WRITE(*,*)'Ro'
       WRITE(*,20)sxr_chords(ichord) % Ro,                                     &
     &            sxr_chords(ichord) % Zo,                                     &
     &            sxr_chords(ichord) % Phio
       WRITE(*,*)'Rf '
       WRITE(*,20)sxr_chords(ichord) % Rf,                                     &
     &            sxr_chords(ichord) % Zf,                                     &
     &            sxr_chords(ichord) % Phif
20     FORMAT(3(f7.3,'  '))

!------integration testing---------------
       CALL safe_open(iounit,istat,diag_file,'replace','formatted')
       ipoints=100
       ALLOCATE(s_array(ipoints))
       dR=(sxr_chords(ichord) % Rf - sxr_chords(ichord) % Ro)/ipoints
       dZ=(sxr_chords(ichord) % Zf - sxr_chords(ichord) % Zo)/ipoints
       dL=sqrt(dR**2+dZ**2)
       WRITE(*,*)'deltas',dR,dZ
       DO istep=1,ipoints
          r_cyl(1)=sxr_chords(ichord) % Ro + istep*dR
          r_cyl(2)=sxr_chords(ichord) % Phio *pi/180.0 * nfp
          r_cyl(3)=sxr_chords(ichord) % Zo + istep*dZ
          CALL cyl2flx(rzl_local, r_cyl, c_flux, ns, ntor, mpol,               & 
     &                         ntmax, lthreed, lasym, info, nfe, fmin)
          
          s_array(istep)=emissivity_model(c_flux(1),model_in)*dL
          WRITE(iounit,*)s_array(istep)
       ENDDO
       WRITE(*,*)'integral= ',SUM(s_array)
       CLOSE(iounit)
       DEALLOCATE(s_array)
              
       
       END SUBROUTINE sxr_signal 

!===============================================================================       
! cyl2flx_test
! tests the cyl2flux routine
! last modified 1/12/10  GJH
!-------------------------------------------------------------------------------       
!         This sections tests the cly2flx and flx2cly functions
!     it works by reading the rmnc, etc. from the wout file.  Need to fix this
!     so that it gathers information from the model.
!------------------------------------------------------------------------------- 
 
       SUBROUTINE cyl2flx_test(model_in)       
!-------------------------------------------------------------------------------
! INPUT
! model_in  -the equilibrium model
!-------------------------------------------------------------------------------       
       
       USE vmec_utils
       
       TYPE (model),INTENT(IN) :: model_in
 
       INTEGER       :: ierr
       CHARACTER     :: woutfile='_.nc'
 
! cyl2flx definitions
       REAL(rprec), DIMENSION(3) :: r_cyl(1:3)=0.0
       REAL(rprec), DIMENSION(3) :: r_cyl2(1:3)=0.0
       REAL(rprec), DIMENSION(3) :: c_flux(1:3)=0.0
       INTEGER :: nfe 
       INTEGER :: info,iflag,points
       REAL(rprec) :: fmin, dr, dz, dphi, terr, perr
!  JDH 2010-07-20. rtest and ztest, changed from REAL to INTEGER
       INTEGER :: rtest, ztest
       
!-------------------------------------------------------------------------------
!  START OF EXECUTABLE CODE
!-------------------------------------------------------------------------------  
       WRITE(*,*)'START OF WOUT INPUTS'
       CALL readw_and_open(woutfile,ierr)
       IF(ierr .ne. 0 ) write(*,*) 'read_wout error'
!-------------------------------------------------------------------------------       
! This section is just testing the wout read capability.  Wout reads will be
! used until the eq_state % aux1(2) variables are accessible  
!-------------------------------------------------------------------------------            
       WRITE(*,*)'mgrid file ',mgrid_file
       WRITE(*,*)'ns=',ns,'mnmax=',mnmax
       WRITE(*,*)'ntor=',ntor,'mpol=',mpol,'ntmax=',ntmax
       WRITE(*,*)'lthreed=',lthreed,'lasym=',lasym
       WRITE(*,*)size(rmnc)
       WRITE(*,*)'nfp=',nfp
       WRITE(*,*)'aspect ratio=',aspect

       !create inputs for cyl2flx call
       !create rzl_in
       IF (ALLOCATED(rzl_local)) THEN 
            WRITE(*,*)'RZL already allocated'
       ELSE
            CALL LoadRZL
       ENDIF
       IF (ALLOCATED(rzl_local)) WRITE(*,*)'RZL now allocated'

       terr=0.0
       perr=0.0
       points=0
       r_cyl(2)=3.14159
       c_flux(3)=3.14159
       DO rtest=0,600 
            r_cyl(1)=0.45+rtest/1000.
            DO ztest=-300,300
                r_cyl(3)=ztest/1000.
                CALL cyl2flx(rzl_local, r_cyl, c_flux, ns, ntor, mpol,         & 
     &                         ntmax, lthreed, lasym, info, nfe, fmin)
                CALL flx2cyl(rzl_local, c_flux, r_cyl2, ns, ntor,              & 
     &                         mpol, ntmax, lthreed, lasym, iflag)
                dr=sqrt((r_cyl(1)-r_cyl2(1))**2)
                dz=sqrt((r_cyl(3)-r_cyl2(3))**2)
                dphi=sqrt((r_cyl(2)-r_cyl2(2))**2)/nfp
                IF(c_flux(1).LE.1.0) THEN
                    terr=terr+sqrt(dr**2+dz**2)
                    perr=perr+dphi
                    points=points+1
                ENDIF
       !WRITE(*,10)r_cyl(1:3),c_flux(1:3),info,nfe,fmin,r_cyl2(1:3)
            END DO
       END DO
10     FORMAT(3(f6.3),'  ',3(f6.3),'  ',i3,i3,'  ',E8.2,'  ',3(f6.3))   
       WRITE(*,*)points,'    error = ',terr/points,'  ',perr/points
       
       END SUBROUTINE cyl2flx_test 

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      END MODULE sxr_routines
!-------------------------------------------------------------------------------
