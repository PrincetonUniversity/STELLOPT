!*******************************************************************************
!  File sxrc_mc.f
!  Contains module sxrc_mc
!           SUBROUTINE sxrc_mc_model_compute
!           FUNCTION emissivity_model
!
!  Module is part of the V3FIT stellarator equilibrium reconstruction code.
!  The postfix _mc in the name indicates Model Compute
!  The module contains the subroutine sxrc_mc_model_compute
!  Both an sxrc_desc and a model are needed for this computation.

!*******************************************************************************
!  MODULE sxrc_mc
!    (SXRC - Soft X-Ray Model Compute Calculations)
! SECTION I.     VARIABLE DECLARATIONS
! SECTION II.    INTERFACE BLOCKS
! SECTION III.   MAIN COMPUTATIONAL SUBROUTINE
!*******************************************************************************
      MODULE sxrc_mc

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Type declarations, constants, utilities
!-------------------------------------------------------------------------------
      USE stel_kinds
      USE v3f_global
      USE vmec_utils
      USE read_wout_mod

!-------------------------------------------------------------------------------
!  SXRC Derived Type
!-------------------------------------------------------------------------------
      USE sxrc_T
      
!-------------------------------------------------------------------------------
!  Model Derived Types
!-------------------------------------------------------------------------------
      USE model_T

!------------------------------------------------------------------------------
      IMPLICIT NONE
      
!*******************************************************************************
! SECTION II.  INTERFACE BLOCKS
!*******************************************************************************

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      CONTAINS
          
!*******************************************************************************
! SECTION III.  MAIN COMPUTATIONAL SUBROUTINE
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Compute an SXRC signal
!
!    Information comes from the sxrc_desc and the model
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  Actual computation of the model signal is in this subroutine
!    s_type = diagnostic
!      d_type = sxrc
!        Soft X-Ray Chord Diagnostic
!        signal_model_compute_sxrc 
!
!-------------------------------------------------------------------------------

      SUBROUTINE sxrc_mc_model_compute(a_sxrc,a_model,mod_signal,              &
     &   mod_sigma)

!-------------------------------------------------------------------------------
! ARGUMENTS
! a_sxrc        - type sxrc_desc - holds sxr chord info
! a_model       - type model - holds eq_state and density model information
! mod_signal    - output of the generated signal
! mod_sigma     - output sigma      
!------------------------------------------------------------------------------- 
      TYPE (sxrc_desc), INTENT (inout)     :: a_sxrc
      TYPE (model), INTENT (inout), TARGET :: a_model
      REAL(rprec), POINTER, DIMENSION(:) :: mod_signal, mod_sigma

!------------------------------------------------------------------------------- 
! local variables
!------------------------------------------------------------------------------- 
      TYPE (eq_state), POINTER  :: state => null()
      TYPE (e_density), POINTER :: density => null()

      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'sxrc_mc_model_compute'
     
!-------------------------------------------------------------------------------      
! cyl2flx declarations  - used for the subroutine cyl2flx call
!------------------------------------------------------------------------------- 
       REAL(rprec), DIMENSION(3) :: r_cyl(1:3)=0.0
       REAL(rprec), DIMENSION(3) :: r_cyl2(1:3)=0.0
       REAL(rprec), DIMENSION(3) :: c_flux(1:3)=0.0
       INTEGER :: nfe 
       INTEGER :: info,iflag,points
       REAL(rprec) :: fmin
       
!------------------------------------------------------------------------------- 
! other decalarations
!
!-------------------------------------------------------------------------------        
        CHARACTER*300 :: diag_file='diagtest2'
        INTEGER       :: iounit=6

        INTEGER       :: ncurrents
        INTEGER       :: n_chords
        INTEGER       :: ierr,istat1
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
         
!         WRITE(*,*) ' *** IN SUBROUTINE ', sub_name
         
!  Point state to the model equilibrium state
         state => a_model % eqstate
         density => a_model % density

!  Allocate data and sigma
         IF (ASSOCIATED(mod_signal)) THEN
           DEALLOCATE(mod_signal,mod_sigma, stat=istat1)
           CALL assert_eq(0,istat1,sub_name //                                 &
     &                     'mod_signal, sigma dealloc')
         ENDIF
         ALLOCATE(mod_signal(1),mod_sigma(1),stat=istat1)
         CALL assert_eq(0,istat1,sub_name // 'mod_signal, sigma alloc')

!-------------------------------------------------------------------------------      
!  Check that aux1 is defined.
!****************************************************************************
!*  Normally a call to define aux2 (which will define aux1) is made         *
!*  The call to define aux2 requires kv knowledge that I don't have or need *
!*  It was decided by JDH and myself to go straight to aux1                 *
!*  GJH 2010-01-26                                                          *
!****************************************************************************
!-------------------------------------------------------------------------------      
         IF (.not. state % l_def_aux1) THEN
           CALL eq_aux_1_define(state)
         END IF
!-------------------------------------------------------------------------------
!create inputs for cyl2flx call
!create rzl_in
! this is not very efficient.  This is done for EACH signal that needs 
! to be computed
!-------------------------------------------------------------------------------
         IF (.NOT. ALLOCATED(rzl_local)) THEN
            WRITE(*,*)'Allocating RZL'
            CALL LoadRZL
         ENDIF
         IF (.NOT. ALLOCATED(rzl_local)) THEN
            WRITE(*,*)'RZL allocation error - returning from sxr_signal'
            RETURN
         ENDIF
!-------------------------------------------------------------------------------       
!         WRITE(*,*)'--------------SIGNAL STUFF--------------------'
!         WRITE(*,*)'chord ',a_sxrc % chord_name
!         WRITE(*,20)a_sxrc % Ro,a_sxrc % Zo,a_sxrc % Phio,                     &
!     &              a_sxrc % Rf,a_sxrc % Zf,a_sxrc % Phif
!20       FORMAT('Ro- ',3(f7.3,'  '),' Rf- ',3(f7.3,'  '))
!------integration testing---------------
! this is not the integration that I want to use in the end
! will change later GJH 2010-01-26


         ipoints=100
         ALLOCATE(s_array(ipoints))
         dR=(a_sxrc % Rf - a_sxrc % Ro)/ipoints
         dZ=(a_sxrc % Zf - a_sxrc % Zo)/ipoints
         dL=sqrt(dR**2+dZ**2)
 !        WRITE(*,*)'delta dR and dZ',dR,dZ
         DO istep=1,ipoints
            r_cyl(1)=a_sxrc % Ro + istep*dR
            r_cyl(2)=a_sxrc % Phio *pi/180.0 * state % fixp % nfp
            r_cyl(3)=a_sxrc % Zo + istep*dZ
            CALL cyl2flx(rzl_local, r_cyl, c_flux, ns, ntor, mpol,             & 
     &                         ntmax, lthreed, lasym, info, nfe, fmin)
            
            s_array(istep)=emissivity_model(c_flux(1),a_model)*dL
         ENDDO
         mod_signal(1)=SUM(s_array)
!         WRITE(*,*)'integral= ',mod_signal(1)
         CLOSE(iounit)
         DEALLOCATE(s_array)
! make up a sigma         
         mod_sigma(1)=0.5
!         WRITE(*,*)'-------leaving ',sub_name
         
      END SUBROUTINE sxrc_mc_model_compute

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
         TYPE (model),INTENT(IN), TARGET :: model_in
         REAL(rprec)             :: emissivity_model
      
         TYPE (eq_state), POINTER  :: state => null()
         TYPE (e_density), POINTER :: density => null()
      
         REAL(rprec) :: P
         REAL(rprec) :: n
         INTEGER     :: j
         
         state => model_in % eqstate
         density => model_in % density
      
      IF (s.LE.1.0) THEN
        n=density % ne_max /1.0E18                                              &
     &          * (1-s**density % tau)**density % kappa                         &
     &          + density % ne_ambient/1.0E18
        P=0.0
        DO j=0,SIZE(state % varp % am)-1
          P=P+(state % varp % am(j))*(s**(j))
        ENDDO
        emissivity_model=n**1.5+SQRT(P)
      ELSE
        emissivity_model=0.0
      ENDIF
      
      END FUNCTION emissivity_model
!-------------------------------------------------------------------------------


!*******************************************************************************
! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!
!  JDH 2009-03-03. First version of module. Based on coding from previous
!    module signal_model
!  GJH 2010-01-20. SXR version of coding based on module mddc_mc

      END MODULE sxrc_mc
!-------------------------------------------------------------------------------
