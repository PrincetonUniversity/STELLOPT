!     SPH: INTEGER(iprec) -> INTEGER
      MODULE b_transform_vmec
!-----------------------------------------------------------------------------------------
! module to compute (flux coord) <---> (cylinderical coord) transformations
!
!  Code is shamelessly plagerized from similar routines
!      (with similar/identical names) in Steve Hirshman's vmec_utils and read_wout_mod modules.
!      For reference, the original VMEC routines used are:
!
!
!      vmec_utils:    Get_Bcyl_Wout, flx2cyl, cyl2flx, newt2d, get_flxcoord
!      read_wout_mod: tosuvspace, LoadRZL
!
!  **Development Notes:
!     1)  This routine, by using local versions of the VMEC routines and consolidating them into a 
!      single module, is now essentially a *stand-alone* module that only needs input wout data to function.
!      Thus it is easily portable to other applications that use VMEC output.
!
!     2) For V3FIT applications, it is ultimately envisioned that a new version of this module
!        will call Steve Hirshman's *ORIGINAL* VMEC routines instead of my local "knock-offs".
!         To this end, I have endeavoured to make as few changes as possible to the original
!        VMEC routines.
!        In addition, I anticipate that future versions will take the wout info from memory 
!        (i.e. from "the model") instead of directly from the wout file as I do now.
!            However, in order for all this to happen, there are a number of wout variables that
!        will be needed to be added to the model.   (these additional variables are used to compute
!        global variables [such as array rzl()] in Hirshman's flux coord transformation routines.
!
!  
! PUBLIC Routines: VMEC_B_INIT, VMEC_B, VMEC_CYL2FLX, VMEC_FLX2CYL, INVALID_RHO_VMEC,
!                  VMEC_CAR2FLX, VMEC_FLX2CAR, TRANS_CYL2CAR, TRANS_CAR2CYL
!                   (all other routines are intended to be for internal use only)
!
! VERSION HISTORY:
!  created 11/8/06 by J. Shields.
!
!  modified 08/04/07 by JMS:  enhanced to give VMEC_B_INIT the ability to take wout info from
!                             either the wout file or a V3FIT "eq_state" variable.
!
!  JDH 2008-01-19 - Added => null() to pointer declarations

!----------------------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants, ONLY: twopi, one, zero
      USE eq_T    ! defines derived TYPE "eq_state"
      IMPLICIT NONE

!........create a derived TYPE to store wout info...........!
      TYPE eq_rzlam
!-----------------------------------------------------
! Derived TYPE to store local versions of wout variables.
!   created 12/8/06 by J. Shields
!   modified by JMS 8/4/07 to include a pointer of TYPE eq_aux_1
!----------------------------------------------------------
        INTEGER :: ns
        INTEGER :: nfp
        INTEGER :: mnmax
        INTEGER :: mnmax_nyq
        INTEGER :: mnyq
        INTEGER :: nnyq
        INTEGER :: mpol
        INTEGER :: ntmax
        INTEGER :: ntor
        REAL(rprec), DIMENSION(:,:,:,:), POINTER :: rzl_local => null()
        LOGICAL        :: lthreed
        LOGICAL        :: lasym

         REAL(rprec), DIMENSION(:,:), POINTER ::  rmnc => null()
         REAL(rprec), DIMENSION(:,:), POINTER ::  zmns => null()
         REAL(rprec), DIMENSION(:,:), POINTER ::  gmnc => null()
         REAL(rprec), DIMENSION(:,:), POINTER ::  bsubumnc => null()
         REAL(rprec), DIMENSION(:,:), POINTER ::  bsubvmnc => null()
         REAL(rprec), DIMENSION(:,:), POINTER ::  bsupumnc => null()
         REAL(rprec), DIMENSION(:,:), POINTER ::  bsupvmnc => null()

         REAL(rprec), DIMENSION(:,:), POINTER ::  rmns => null()
         REAL(rprec), DIMENSION(:,:), POINTER ::  zmnc => null()
         REAL(rprec), DIMENSION(:,:), POINTER ::  gmns => null()
         REAL(rprec), DIMENSION(:,:), POINTER ::  bsupumns => null()
         REAL(rprec), DIMENSION(:,:), POINTER ::  bsupvmns => null()

!    m and n value arrays
         REAL(rprec), DIMENSION(:), POINTER   ::  xm => null() 
         REAL(rprec), DIMENSION(:), POINTER   ::  xn => null() 
         REAL(rprec), DIMENSION(:), POINTER   ::  xm_nyq => null() 
         REAL(rprec), DIMENSION(:), POINTER   ::  xn_nyq => null()

      END TYPE eq_rzlam

      INTEGER, PRIVATE, SAVE   :: eq_data_source
      INTEGER, PRIVATE   :: ns_loc, ntmax_loc, mpol_loc, ntor_loc
      REAL(rprec), PRIVATE :: r_target, phi_target, z_target, fnorm
      REAL(rprec), ALLOCATABLE :: rzl_loc(:,:,:,:)
      REAL(rprec), POINTER, PRIVATE :: rzl_ptr(:,:,:,:) => null()
      REAL(rprec), POINTER, PRIVATE :: mscale_loc(:) => null()
      REAL(rprec), POINTER, PRIVATE :: nscale_loc(:) => null()
      LOGICAL, PRIVATE :: lthreed_loc, lasym_loc, lscale
      INTEGER, PARAMETER, PRIVATE   :: FROM_WOUT_FILE = 0
      INTEGER, PARAMETER, PRIVATE   :: FROM_V3FIT = 1


      TYPE(eq_rzlam), PRIVATE, SAVE :: eq_rzl

      PRIVATE :: newt2d, get_flxcoord
      PRIVATE :: cyl2flx_local, flx2cyl_local
      PRIVATE :: tosuvspace_local_wout, tosuvspace_local_v3fit
      PRIVATE :: Get_Bcyl_JS, Get_Cflx_JS
      PRIVATE :: LoadRZL_loc_wout
      PRIVATE :: fill_eq_rzl_wout, eq_rzlam_construct_v3fit


!=====================================================================
      CONTAINS

!==========================================================================
      SUBROUTINE VMEC_B_INIT(eq_source_in)
!----------------------------------------------------------------
!
! FUNCTION: computes B field in cylin coords using "eq_state" info
!
! INPUT:  eq_source_in (optional):   0 = default = old wout interface
!                                    1 = new V3FIT interface
!
! created 11/7/06 by J. Shields, based on Steve Hirshman's
!         Get_Bcyl_Wout routine in module vmec_utils
!
! updated 8/7/07 by JMS to interface better with V3FIT, but with the intent to remain
!                backwards-compatible with older code.
!
!
! Development Notes: Basically, for the V3FIT version, read_wout_mod is ONLY used for VMEC_B_INIT
!                 to fill the TYPE variable eq_rzl.  After that, all flux coordinate transformations
!                 receive their information directly from eq_rzl.
!
!----------------------------------------------------------------
      USE read_wout_mod, ONLY: lwout_opened
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in), OPTIONAL :: eq_source_in 

C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: fmin_acceptable = 1.E-6_dp
      INTEGER     :: nfe, info_loc
      INTEGER  :: ntmax_temp
      INTEGER, PARAMETER :: EQ_SOURCE_DEFAULT = FROM_WOUT_FILE
      REAL(rprec) :: r_cyl(3), c_flx(3), fmin
      REAL(rprec) :: Ru1, Zu1, Rv1, Zv1
      REAL(rprec) :: bsupu1, bsupv1
      LOGICAL     :: lthreed_temp
      CHARACTER(len=*), PARAMETER :: subname = 'VMEC_B_INIT: '
C-----------------------------------------------
      write(*,*) 'now in ', subname

      IF (.not.lwout_opened) THEN
         WRITE(6, '(2a,/,a)')
     1   ' This form of GetBcyl can only be called if WOUT has been',
     2   ' previously opened!',' Try GetBcyl_VMEC form instead.'
         RETURN
      END IF

!.........set the module global toggle to indicate the V3FIT/WOUT setting..........!
      if( PRESENT(eq_source_in) ) then
        SELECT CASE (eq_source_in)
        CASE (FROM_V3FIT)
          eq_data_source = FROM_V3FIT
        CASE (FROM_WOUT_FILE)
          eq_data_source = FROM_WOUT_FILE
        CASE DEFAULT
          eq_data_source = EQ_SOURCE_DEFAULT
        END SELECT
      else
        eq_data_source = EQ_SOURCE_DEFAULT
      end if

!.........fill the local version of the VMEC RZL array ...........................!
      CALL LoadRZL_loc_wout(ntmax_temp, lthreed_temp)

!..........use the read_wout_mod variables to fill the (module-global) eq_rzl TYPE variable.......!
      SELECT CASE (eq_data_source)
      CASE (FROM_V3FIT)
        call eq_rzlam_construct_v3fit(eq_rzl,ntmax_temp,lthreed_temp)
      CASE (FROM_WOUT_FILE)
        call fill_eq_rzl_wout(ntmax_temp, lthreed_temp)
      CASE DEFAULT
        write(*,*) subname, 'ERROR in eq_data_source ASSIGNMENT!!!'
      END SELECT


!........now that eq_rzl%rzl_local has been filled, deallocate the (module-global) rzl_loc array..................!
      DEALLOCATE(rzl_loc)


      END SUBROUTINE VMEC_B_INIT
     

!===========================================================================
      SUBROUTINE VMEC_B(r_cyl, B_cyl, ierr, R_FLX, B_CON)
!--------------------------------------------------------
! FUNCTION: computes the B field (in cylinderical coords) for a given
!           input cylinderical position vector (r_cyl).
!              It is essentially just an interface portal to the VMEC-derived
!           "Get_Bcyl" routines.
!
!
! **NOTE: VMEC_B_INIT *must* be called prior to this routine!!!
!-------------------------------------------------------
      IMPLICIT NONE

!.........dummy arguments.........................!
      INTEGER, INTENT(OUT)     :: ierr
      REAL(rprec), INTENT(IN)  :: r_cyl(3)    ! postion vector (cylin coords)
      REAL(rprec), INTENT(OUT) :: B_cyl(3)    ! postion vector (cylin coords)
      REAL(rprec), OPTIONAL, INTENT(OUT) :: R_FLX(3)    ! postion vector (flux coords)
      REAL(rprec), OPTIONAL, INTENT(OUT) :: B_CON(3)    ! postion vector (flux coords)

!...........local variables.......................!
      INTEGER     :: istat
      REAL(rprec) :: R, phi, Z
      REAL(rprec) :: sflx, uflx
      REAL(rprec) :: Br, Bphi, Bz
      CHARACTER(len=*), PARAMETER :: subname = 'VMEC_B: '


      R   = r_cyl(1)
      phi = r_cyl(2)
      Z   = r_cyl(3)


!      call VMEC_B_INIT

      call Get_Bcyl_JS( R, phi, Z, B_cyl(1),B_cyl(2), B_cyl(3), sflx,          &
     &                    uflx, B_CON, istat)


      if( PRESENT(R_FLX) ) then
        R_FLX(1) = sflx
        R_FLX(2) = uflx
        R_FLX(3) = phi       
      end if

      ierr = istat

      END SUBROUTINE VMEC_B

!==========================================================================
           SUBROUTINE VMEC_CYL2FLX(r_cyl, r_flx, ierr)
!---------------------------------------------------------------------------
! FUNCTION: converts input cylin coord vector (r_cyl) into flux coords
!
!  INPUT:  r_cyl    
!  OUTPUT: r_flx, ierr
!
! Note: r_flx is defined as (s,u,phi) *NOT* (s,u, N*phi)
!
!-----------------------------------------------------------------------------
      IMPLICIT NONE

!.........dummy arguments.........................!
      INTEGER, INTENT(OUT)     :: ierr
      REAL(rprec), INTENT(IN)  :: r_cyl(3)    ! CYLIN COORDS postion vector (R,phi,Z)
      REAL(rprec), INTENT(OUT) :: r_flx(3)    ! FLUX  COORDS postion vector (s,u,phi)


!...........local variables.......................!
      INTEGER     :: istat
      REAL(rprec) :: R, phi, Z
      REAL(rprec) :: c_flx(3)        ! "interim" FLUX COORDS postion vector (s,u, N*phi)
      REAL(rprec) :: Br, Bphi, Bz

      R   = r_cyl(1)
      phi = r_cyl(2)
      Z   = r_cyl(3)

!..........compute flux coord vector c_flx = (s, u, N*phi)......................!
      call Get_Cflx_JS( R, phi, Z, c_flx, istat)

!.......convert (s,u, N*phi) to (s,u,phi)
      r_flx(1) = c_flx(1)
      r_flx(2) = c_flx(2)
      r_flx(3) = phi

      ierr = istat

      END SUBROUTINE VMEC_CYL2FLX


!====================================================================================
      SUBROUTINE VMEC_FLX2CYL_JS(r_flx, r_cyl, info)
!----------------------------------------------------------------
!
! FUNCTION: converts flux coords to cylin coords using wout file info
!
!  INPUT:  r_flx     ! flux coord vector (s,u,phi)
!  OUTPUT: r_cyl     ! cylin vector (R, phi, Z)
!          info      ! error indicator
!
! module-global inputs: eq_rzl 
!
!  ***important note: Hirshman uses v = (# field periods)*phi as his toroidal flux coordinate
!                     in the transformations, while the input r_flx vector assumes zeta=phi 
!
! created 12/5/06 by J. Shields, based on Steve Hirshman's
!         Get_Bcyl_Wout routine in module vmec_utils  
!
! **NOTE:  VMEC_B_INIT *must* be called prior to the calling of this routine!
!  
!---------------------------------------------------------------------------
!      USE read_wout_mod, ONLY: lwout_opened
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, OPTIONAL, INTENT(out) :: info
      REAL(rprec), INTENT(in) :: r_flx(3)      ! input flux vector (s,u,phi)
      REAL(rprec), INTENT(out) :: r_cyl(3)     ! cylin coord vector (R,phi,Z)

C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: fmin_acceptable = 1.E-6_dp
      INTEGER     :: nfe, info_loc
      INTEGER     :: ns_w, ntor_w, nfp_w
      INTEGER     :: mpol_w, ntmax_w
      REAL(rprec) :: c_flx(3)                ! "interim" flux vector (s, u, N*phi)
      REAL(rprec) :: r_cyl_temp(3)           ! "interim" cylin vector (R, N*phi, Z)
      REAL(rprec) :: Ru1, Zu1, Rv1, Zv1      ! derivatives:  dR/du, dZ/du, dR/dv, dZ/dv
      REAL(rprec) :: bsupu1, bsupv1
      REAL(rprec) :: fmin
      LOGICAL     :: lthreed_w, lasym_w
      CHARACTER(len=*), PARAMETER :: subname = 'VMEC_FLX2CYL_JS: '
C-----------------------------------------------
!      IF (.not.lwout_opened) THEN
!         WRITE(6, '(2a,/,a)')
!     1   ' This form of GetBcyl can only be called if WOUT has been',
!     2   ' previously opened!',' Try GetBcyl_VMEC form instead.'
!         RETURN
!      END IF

!      CALL LoadRZL

!       write(*,*) 'now in ', subname

!........copy read_wout_mod global variables into local versions............!
!...........(updated 12/9/06 to use the eq_rzl derived TYPE variable)
      ns_w      = eq_rzl%ns
      nfp_w     = eq_rzl%nfp
      ntor_w    = eq_rzl%ntor
      mpol_w    = eq_rzl%mpol
      ntmax_w   = eq_rzl%ntmax
      lthreed_w = eq_rzl%lthreed
      lasym_w   = eq_rzl%lasym


!.......convert (rho,theta, phi) to (s,u,v).................................!
      c_flx(1) = r_flx(1)
      c_flx(2) = r_flx(2)
      c_flx(3) = nfp_w * r_flx(3)


!.......transform flux vector c_flx into cylinderical coords.........................!
      CALL flx2cyl_local(eq_rzl%rzl_local, c_flx, r_cyl_temp, ns_w,             & 
     &         ntor_w, mpol_w, ntmax_w, lthreed_w, lasym_w,                     &
     &         info_loc, RU=Ru1, ZU=Zu1, RV=Rv1, ZV=Zv1)

!........deal with the annoying "phi vs nfp*phi" thing...................!
      Rv1 = nfp_w*Rv1;  Zv1 = nfp_w*Zv1
      r_cyl(1) = r_cyl_temp(1)
      r_cyl(2) = r_cyl_temp(2) / nfp_w
      r_cyl(3) = r_cyl_temp(3)



      IF (PRESENT(info)) info = info_loc
      IF (info_loc .ne. 0) RETURN


      END SUBROUTINE VMEC_FLX2CYL_JS

!===================================================================================
      SUBROUTINE VMEC_CAR2FLX(xcar, xflux, iflag)
!--------------------------------------------------------------------------------
!
! FUNCTION:  Transforms Cartesian POSITION vector "xcar" into flux coordinates
!              (note that it does NOT do B field vector transforms!!)
!
! LOGIC:  uses TRANS_CAR2CYL to transform from Cartesian to cylinderical and then
!         uses VMEC_CYL2FLX to transform from cylinderical to flux coordinates.
!
!----------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: iflag

      REAL(rprec), DIMENSION(3), INTENT(IN) ::                                 &
     &            xcar                 ! Cartesian coordinate position vector
      REAL(rprec), DIMENSION(3) ::                                             &
     &            x_cyl                ! Cylindrical coordinate position vector
      REAL(rprec), DIMENSION(3), INTENT(OUT)::                                 &
     &            xflux                ! flux coordinate vector

      CHARACTER(len=80) :: message

!      write(*,*) 'Hi.  Now in CAR2FLX.'

!...........transform Cartesian  xcar() vector into cylin coord vector x_cyl() 
      call TRANS_CAR2CYL(xcar, x_cyl)

!...........transform Cylin x_cyl() vector into flux coord vector xflux()
      call VMEC_CYL2FLX(x_cyl, xflux, iflag)
      
      RETURN
      END SUBROUTINE VMEC_CAR2FLX

!===================================================================================
      SUBROUTINE VMEC_FLX2CAR(xflux, xcar,iflag)
!--------------------------------------------------------------------------------
!
! FUNCTION:  Transforms flux coordinates POSITION vector "xflux" into Cartesian coordinates
!              (note that it does NOT do B field vector transforms!!)
!
! LOGIC:  uses VMEC_FLX2CYL  to transform from flux to cylinderical and then
!         uses TRANS_CYL2CAR to transform from cylinderical to Cartesian coordinates.
!
!----------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
      IMPLICIT NONE

      INTEGER :: iflag
      CHARACTER(len=80) :: message

      REAL(rprec), DIMENSION(3), INTENT(IN)::                                  &
     &            xflux                ! flux coordinate vector
      REAL(rprec), DIMENSION(3), INTENT(OUT) ::                                &
     &            xcar                 ! Cartesian coordinate vector
      REAL(rprec), DIMENSION(3) ::                                             &
     &            x_cyl                ! Cylindrical coordinate vector


!      write(*,*) 'Hi.  Now in FLX2CAR.'

!...........transform flux vector xflux() into cylin coord vector x_cyl() 
      call VMEC_FLX2CYL_JS(xflux, x_cyl, iflag)

!...........transform Cylin x_cyl() vector into Cartesian coord vector xcar
      call TRANS_CYL2CAR(x_cyl, xcar)
      
      RETURN
      END SUBROUTINE VMEC_FLX2CAR





!==========================================================================
      SUBROUTINE Get_Bcyl_JS(R1, Phi, Z1, Br, Bphi, Bz,                        &  
     &                        sflx, uflx, B_SUPER, info)
!----------------------------------------------------------------
!
! FUNCTION: Computes cylindrical components of the magnetic field, Br, Bphi, Bz,
!     at the specified cylindrical coordinate point (R1, Phi, Z1) using the wout info.
!     Phi is the true geometric toroidal angle (NOT N*Phi)
!
!
! created 11/7/06 by J. Shields, based on Steve Hirshman's
!         Get_Bcyl_Wout routine in module vmec_utils
!
!     INPUT
!     R1, Phi, Z1  : cylindrical coordinates at which evaluation is to take place
!     
!     OUTPUT
!     Br, Bphi, Bz : computed cylindrical components of B at input point
!     sflx, uflx   : computed flux and theta angle at the cylindrical point
!
!     module global inputs:  eq_rzl
!
! IMPORTANT:  VMEC_B_INIT *must* be run before this routine is called!!  JS 11/9/06
!--------------------------------------------------------------------
!      USE read_wout_mod, ONLY: lwout_opened
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, OPTIONAL, INTENT(out) :: info
      REAL(rprec), INTENT(in)  :: R1, Z1, Phi
      REAL(rprec), INTENT(out) :: Br, Bphi, Bz
      REAL(rprec), INTENT(out), OPTIONAL :: sflx, uflx
      REAL(rprec), INTENT(out), OPTIONAL :: B_SUPER(3)

C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: fmin_acceptable = 1.E-6_dp
      INTEGER     :: nfe, info_loc
      INTEGER     :: ns_w, ntor_w, nfp_w
      INTEGER     :: mpol_w, ntmax_w
      REAL(rprec) :: r_cyl(3), c_flx(3), fmin
      REAL(rprec) :: Ru1, Zu1, Rv1, Zv1
      REAL(rprec) :: bsupu1, bsupv1
      LOGICAL     :: lthreed_w, lasym_w
      CHARACTER(len=*), PARAMETER :: subname = 'Get_Bcyl_JS: '

C-----------------------------------------------
!      write(*,*) 'Now in ', subname

!      IF (.not.lwout_opened) THEN
!         WRITE(6, '(2a,/,a)')
!     1   ' This form of GetBcyl can only be called if WOUT has been',
!     2   ' previously opened!',' Try GetBcyl_VMEC form instead.'
!         RETURN
!      END IF

!      CALL LoadRZL



!........copy read_wout_mod global variables into local versions............!
!.........(now using eq_rzl derived TYPE instead of read_wout_mod.  JS 12/9/06 )
      ns_w      = eq_rzl%ns
      nfp_w     = eq_rzl%nfp
      ntor_w    = eq_rzl%ntor
      mpol_w    = eq_rzl%mpol
      ntmax_w   = eq_rzl%ntmax
      lthreed_w = eq_rzl%lthreed
      lasym_w   = eq_rzl%lasym

!.....Convert to point in flux-coordinates: cflux = (s, u, v=N*phi) & evaluate Ru, Zu, Rv, Zv at that pt ...!
      r_cyl(1) = R1;  r_cyl(2) = nfp_w*Phi;  r_cyl(3) = Z1
      c_flx(1) = 0;   c_flx(2) = 0;        c_flx(3) = r_cyl(2)

!.........transform flux coord vector c_flx into cylinderical coords................!
      CALL cyl2flx_local(eq_rzl%rzl_local, r_cyl, c_flx, ns_w, ntor_w,         & 
     &     mpol_w, ntmax_w, lthreed_w, lasym_w, info_loc, nfe, fmin,           &
     &     RU=Ru1, ZU=Zu1, RV=Rv1, ZV=Zv1)

      Rv1 = nfp_w*Rv1;  Zv1 = nfp_w*Zv1


      IF (info_loc.eq.-1 .and. (fmin .le. fmin_acceptable)) info_loc = 0

      IF (PRESENT(info)) info = info_loc
      IF (info_loc .ne. 0) RETURN

      IF (PRESENT(sflx)) sflx = c_flx(1)  
      IF (PRESENT(uflx)) uflx = c_flx(2)

!
!....... Evaluate Bsupu, Bsupv at this point......................................!
      if ( eq_data_source == FROM_WOUT_FILE) then
        CALL tosuvspace_local_wout(c_flx(1), c_flx(2), c_flx(3),                      & 
     &                 BSUPU=bsupu1, BSUPV=bsupv1)
      else if ( eq_data_source == FROM_V3FIT) then
        CALL tosuvspace_local_v3fit(c_flx(1), c_flx(2), c_flx(3),                      & 
     &                 BSUPU=bsupu1, BSUPV=bsupv1)
      else
        write(*,*) subname, 'ERROR in eq_data_source ASSIGNMENT!!!!'
      end if
!...........Form Br, Bphi, Bz..................................................!
      Br   = Ru1*bsupu1 + Rv1*bsupv1
      Bphi = R1 *bsupv1
      Bz   = Zu1*bsupu1 + Zv1*bsupv1
      
      if( PRESENT( B_SUPER) ) then
        B_SUPER(1) = 0.0_rprec
        B_SUPER(2) = bsupu1
        B_SUPER(3) = bsupv1
      end if

      END SUBROUTINE Get_Bcyl_JS



!==========================================================================
      SUBROUTINE Get_Cflx_JS(R1, Phi, Z1, c_flx_out, info)
!----------------------------------------------------------------
!
! FUNCTION: converts cylin coords to flux coords using wout file info
!
! INPUT: R1, Phi, Z1
! OUTPUT: c_flx_out, info    ! NOTE: flux vector c_flx_out = (s, u, N*phi)   
!
! module-global input: eq_rzl
!
! created 11/7/06 by J. Shields, based on Steve Hirshman's
!         Get_Bcyl_Wout routine in module vmec_utils
!
! IMPORTANT:  VMEC_B_INIT *must* be run before this routine is called!!  JS 11/9/06
!--------------------------------------------------------------------
!      USE read_wout_mod, ONLY: lwout_opened
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, OPTIONAL, INTENT(out) :: info
      REAL(rprec), INTENT(in)  :: R1, Z1, Phi
      REAL(rprec), INTENT(out) :: c_flx_out(3)

C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: fmin_acceptable = 1.E-6_dp
      INTEGER     :: nfe, info_loc
      INTEGER     :: ns_w, ntor_w, nfp_w
      INTEGER     :: mpol_w, ntmax_w
      REAL(rprec) :: r_cyl(3), c_flx(3), fmin
      REAL(rprec) :: Ru1, Zu1, Rv1, Zv1
      REAL(rprec) :: Br, Bphi, Bz
      REAL(rprec) :: bsupu1, bsupv1
      LOGICAL     :: lthreed_w, lasym_w
      CHARACTER(len=*), PARAMETER :: subname = 'Get_Cflx_JS: '

C-----------------------------------------------
!      write(*,*) 'Now in ', subname

!      IF (.not.lwout_opened) THEN
!         WRITE(6, '(2a,/,a)')
!     1   ' This form of GetBcyl can only be called if WOUT has been',
!     2   ' previously opened!',' Try GetBcyl_VMEC form instead.'
!         RETURN
!      END IF

!      CALL LoadRZL



!........copy read_wout_mod global variables into local versions............!
      ns_w      = eq_rzl%ns
      nfp_w     = eq_rzl%nfp
      ntor_w    = eq_rzl%ntor
      mpol_w    = eq_rzl%mpol
      ntmax_w   = eq_rzl%ntmax
      lthreed_w = eq_rzl%lthreed
      lasym_w   = eq_rzl%lasym

!........Convert to point in flux-coords: cflux = (s, u, v=N*phi) & evaluate Ru, Zu, Rv, Zv at that pt...!
      r_cyl(1) = R1;  r_cyl(2) = nfp_w*Phi;  r_cyl(3) = Z1
      c_flx(1) = 0;   c_flx(2) = 0;        c_flx(3) = r_cyl(2)

      CALL cyl2flx_local(eq_rzl%rzl_local, r_cyl, c_flx, ns_w, ntor_w, 
     1     mpol_w, ntmax_w, lthreed_w, lasym_w, info_loc, nfe, fmin, 
     2     RU=Ru1, ZU=Zu1, RV=Rv1, ZV=Zv1)

      Rv1 = nfp_w*Rv1;  Zv1 = nfp_w*Zv1


      IF (info_loc.eq.-1 .and. (fmin .le. fmin_acceptable)) info_loc = 0


      IF (PRESENT(info)) info = info_loc
      IF (info_loc .ne. 0) RETURN


      c_flx_out = c_flx

      END SUBROUTINE Get_Cflx_JS




!===================================================================================
      LOGICAL FUNCTION INVALID_RHO_VMEC(x_cyl)
!--------------------------------------------------------------------------------
!
! Function: Tests to verify that the radial flux coordinate (rho) satisfies 0.0 < rho < 1.0
!
!----------------------------------------------------------------------------------
      IMPLICIT NONE

!............dummy variables......................................!
      REAL(rprec), DIMENSION(3), INTENT(IN) :: x_cyl  ! position vector in CYLIN coords (R,phi,Z)


!............local variables......................................!
      INTEGER :: i, j, iflag
      INTEGER :: ierr_vmec
      REAL(rprec), DIMENSION(3) :: dims
      REAL(rprec), DIMENSION(3) :: xflux
      REAL(rprec), DIMENSION(3) :: xflux_vmec
      LOGICAL :: invalid_rho_temp
      CHARACTER(len=80)  :: message
      CHARACTER(len=*), PARAMETER :: subname = 'INVALID_RHO_VMEC: '
     

!      write(*,*) 'Hi.  Now in ', subname

!.........initialize logical...................................!
      invalid_rho_temp = .false.


      call VMEC_CYL2FLX(x_cyl, xflux_vmec, ierr_vmec) 
      if ( ierr_vmec .NE. 0 ) then
        invalid_rho_temp = .true.
      else
        invalid_rho_temp = .false.
      end if


!....**hard-wire invalid_rho_temp to "false" for debugging
!      invalid_rho_temp = .false. 


      INVALID_RHO_VMEC = invalid_rho_temp

      RETURN
      END FUNCTION INVALID_RHO_VMEC

!=========================================================================
      SUBROUTINE TRANS_CAR2CYL(r_car, r_cyl)
!-------------------------------------------------------------------------------
!
! FUNCTION: Converts from Cartesian to cylindrical coordinates
!
! LOGIC:   R  = SQRT(x^2 + y^2)
!         phi = ARCTAN( y/x )
!          Z  = z
!
! Created 12/06 by J. Shields
!
!-------------------------------------------------------------------------------

!..........dummy variables...........................!
      REAL(rprec), DIMENSION(3), INTENT(IN)  :: r_car    !Cartesian coordinates (x,y,z)
      REAL(rprec), DIMENSION(3), INTENT(OUT) :: r_cyl    !cylindrical coordinates (R,phi,Z) [phi in radians]

!...........local variables........................!
      REAL(rprec), PARAMETER :: Z_PI = 3.141592654 

!.........initialize output array.......!
      r_cyl = 0.0_rprec


!......Do Cartesian to cylindrical conversion......................!

      r_cyl(1) = SQRT(r_car(1)**2 + r_car(2)**2) !R
      r_cyl(2) = ATAN2(r_car(2), r_car(1))       !phi
      r_cyl(3) = r_car(3)                        !Z

!..........Ensure 0 <= phi <= 2*PI ...............................!
      IF(r_cyl(2) < 0.0) r_cyl(2) = r_cyl(2) + 2.0*Z_PI

      END SUBROUTINE TRANS_CAR2CYL


!================================================================================
      SUBROUTINE TRANS_CYL2CAR(r_cyl, r_car)
!-------------------------------------------------------------------------------
!
! FUNCTION: Converts from cylindrical to Cartesian coordinates
!
! LOGIC: x = R *COS(phi)
!        y = R *SIN(phi)
!        z = Z
!
! Created 12/06 by J. Shields
!
!-------------------------------------------------------------------------------

!.............dummy variables...............................!
      REAL(rprec), DIMENSION(3), INTENT(IN)  :: r_cyl    ! cylindrical coordinates (R,phi,Z) [phi in radians]
      REAL(rprec), DIMENSION(3), INTENT(OUT) :: r_car    ! Cartesian coordinates (x,y,z)

!........Initialize output array...........................!
      r_car = 0.0_rprec

!........ Do Cylindrical to Cartesian conversion............!
      r_car(1) = r_cyl(1) * COS(r_cyl(2))   !x
      r_car(2) = r_cyl(1) * SIN(r_cyl(2))   !y
      r_car(3) = r_cyl(3)                   !z

      END SUBROUTINE TRANS_CYL2CAR


!====================================================================
      SUBROUTINE flx2cyl_local(rzl_array, c_flux, r_cyl, ns, ntor, 
     1                   mpol, ntmax, lthreed, lasym, iflag,
     2                   mscale, nscale, Ru, Rv, Zu, Zv)
!---------------------------------------------------------------
!
!  FUNCTION: COMPUTES THE CYLINDRICAL COORDINATES R11 and Z11
!            AT A FIXED FLUX COORDINATE POINT si, ui(theta), vi(N*phi)
!            WHERE phi = geometric toroidal angle (0 to 2pi), N = no. field periods
!
!     INPUT:
!     c_flux:          array of (si,ui,vi) values to convert to cylindrical coordinates
!     rzl_array:       array of (r, z, lambda) Fourier coefficients for all radial, m, n values
!     ns, mpol,ntor, ntmax:  radial, poloidal, toroidal, type (r,z,l) dimensions of rzl_array
!     mscale, nscale:  option scale factors. Use ONLY if rzl_array comes from within VMEC.
!                      If arising from WOUT file, mscale, nscale => 1 are not passed
!
!     OUTPUT:
!     r_cyl    :       R = r_cyl(1);  N*PHI = r_cyl(2);   Z = r_cyl(3)
!
!     OPTIONAL OUTPUT
!                      Ru = dR/du;    Rv = dR/dv = dR/dphi / N
!                      Zu = dZ/du;    Zv = dZ/dv = dZ/dphi / N
!
!     NOTE:            User is responsible for multiplying Rv, Zv by N to get phi derivatives
!
!  **this version created 11/24/06 by J. Shields.  Based on S. Hirshman's flx2cyl
!     routine in module vmec_utils
!         Note that, in this version of flx2cyl, the rzl_array is a dummy input & hence a purely local variable.
!
!------------------------------------------------------------    
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(out) :: iflag
      INTEGER, INTENT(in)  :: ns, ntmax, mpol, ntor
      LOGICAL :: lthreed, lasym
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol-1,3*ntmax),
     1   INTENT(in) :: rzl_array
      REAL(rprec), INTENT(in) :: c_flux(3)
      REAL(rprec), INTENT(out) :: r_cyl(3)
      REAL(rprec), INTENT(in), OPTIONAL :: 
     1                            mscale(0:mpol-1), nscale(0:ntor)
      REAL(rprec), INTENT(out), OPTIONAL :: Ru, Rv, Zu, Zv
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: rcc = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: rss, rsc, rcs, zsc, zcs, zcc, zss
      INTEGER :: istat, jslo, jshi, mpol1, m, n
      REAL(rprec), DIMENSION(0:ntor,0:mpol-1) ::
     1           rmncc, rmnss, zmncs, zmnsc,
     2           rmncs, rmnsc, zmncc, zmnss
      REAL(rprec) :: wlo, whi, wlo_odd, whi_odd, hs1, 
     1               si, ui, vi, r11, z11
      REAL(rprec) :: wmins(0:ntor,0:mpol-1),
     1               wplus(0:ntor,0:mpol-1)
      REAL(rprec) :: cosu, sinu, cosv, sinv,
     1               cosmu(0:mpol-1), sinmu(0:mpol-1),
     2               cosnv(0:ntor),  sinnv(0:ntor), 
     3               cosnvn(0:ntor), sinnvn(0:ntor)
      REAL(rprec) :: work1(0:mpol-1,12)
      LOGICAL :: lru, lrv, lzu, lzv
C-----------------------------------------------

      iflag = -1
      si = c_flux(1);  ui = c_flux(2);  vi = c_flux(3)
      r_cyl(2) = vi

!      IF (si.lt.zero .or. si.gt.one) THEN
      IF (si .lt. zero) THEN
         WRITE(6, *)' In flx2cyl_local, s(flux) must be > 0'
         RETURN
      END IF

      lru = PRESENT(ru); lrv = PRESENT(rv)
      lzu = PRESENT(zu); lzv = PRESENT(zv)
      IF (lrv .and. .not. lthreed) rv = 0
      IF (lzv .and. .not. lthreed) zv = 0

!
!     FIND INTERPOLATED s VALUE AND COMPUTE INTERPOLATION WEIGHTS wlo, whi (for even m modes)
!     AND wlo_odd, whi_odd (for odd m modes).
!     FOR si <= 1, POINT IS INSIDE PLASMA;
!     FOR si > 1, TRY EXTRAPOLATION (WITH CONTINUOUS s, u DERIVATIVES INTO "vacuum" REGION
!
      hs1 = one/(ns-1)
      IF (si .le. one) THEN
 
         jslo = 1 + si/hs1
         jshi = jslo+1
         wlo = (hs1*(jshi-1) - si)/hs1
         whi = 1 - wlo
         IF (jslo .eq. ns) jshi = jslo

!
!     USE Rmn, Zmn ~ SQRT(s) FOR ODD-m MODES, SO INTERPOLATE Xmn/SQRT(s)
! 
         whi_odd = whi*SQRT(si/(hs1*(jshi-1)))
         IF (jslo .ne. 1) THEN
            wlo_odd = wlo*SQRT(si/(hs1*(jslo-1)))
         ELSE
            wlo_odd = 0
            whi_odd = SQRT(si/(hs1*(jshi-1)))
         END IF

      ELSE

         jshi = ns
         jslo = ns-1
         wlo  = -(si - 1)/hs1;    wlo_odd = wlo
         whi  = 1 - wlo;          whi_odd = whi

      ENDIF


      mpol1 = mpol-1

      wmins(:,0:mpol1:2) = wlo
      wplus(:,0:mpol1:2) = whi
      wmins(:,1:mpol1:2) = wlo_odd
      wplus(:,1:mpol1:2) = whi_odd



      IF (.not.lasym) THEN
         IF (lthreed) THEN
            IF (ntmax .ne. 2) STOP 'ntmax != 2 in flx2cyl_local!'
            rss = 2;  zcs = 2
         ELSE
            IF (ntmax .ne. 1) STOP 'ntmax != 1 in flx2cyl_local!'
         END IF
      ELSE
         IF (lthreed) THEN
             IF (ntmax .ne. 4) STOP 'ntmax != 4 in flx2cyl_local!'
             rss = 2;  rsc = 3;  rcs = 4
             zcs = 2;  zcc = 3;  zss = 4
         ELSE
             IF (ntmax .ne. 2) STOP 'ntmax != 2 in flx2cyl_local!'
             rsc = 2;  zcc = 2
         END IF
      END IF

      zsc = 1+ntmax; zcs = zcs+ntmax; zcc = zcc+ntmax

      rmncc = wmins*rzl_array(jslo,:,:,rcc) 
     1      + wplus*rzl_array(jshi,:,:,rcc)        !!COS(mu) COS(nv)
      zmnsc = wmins*rzl_array(jslo,:,:,zsc) 
     1      + wplus*rzl_array(jshi,:,:,zsc)        !!SIN(mu) COS(nv)

      IF (lthreed) THEN
         rmnss = wmins*rzl_array(jslo,:,:,rss) 
     1         + wplus*rzl_array(jshi,:,:,rss)     !!SIN(mu) SIN(nv)
         zmncs = wmins*rzl_array(jslo,:,:,zcs)
     1         + wplus*rzl_array(jshi,:,:,zcs)     !!COS(mu) SIN(nv)
      END IF

!
!     SETUP TRIG ARRAYS
!
      cosu = COS(ui);   sinu = SIN(ui)
      cosv = COS(vi);   sinv = SIN(vi)

      cosmu(0) = 1;    sinmu(0) = 0
      cosnv(0) = 1;    sinnv(0) = 0
      DO m = 1, mpol1
         cosmu(m) = cosmu(m-1)*cosu - sinmu(m-1)*sinu
         sinmu(m) = sinmu(m-1)*cosu + cosmu(m-1)*sinu
      END DO

      IF (PRESENT(mscale)) THEN
         cosmu = cosmu*mscale;  sinmu = sinmu*mscale
      END IF

      DO n = 1, ntor
         cosnv(n) = cosnv(n-1)*cosv - sinnv(n-1)*sinv
         sinnv(n) = sinnv(n-1)*cosv + cosnv(n-1)*sinv
      END DO

      IF (PRESENT(nscale)) THEN
         cosnv = cosnv*nscale;  sinnv = sinnv*nscale
      END IF

      IF (lrv .or. lzv) THEN
         DO n = 0, ntor
            cosnvn(n) = n*cosnv(n)
            sinnvn(n) =-n*sinnv(n)
         END DO
      END IF

      iflag = 0

!
!     COMPUTE R11, Z11 IN REAL SPACE
!
!     FIRST, INVERSE TRANSFORM IN N-V SPACE, FOR FIXED M
!
      DO m = 0, mpol1
 
         work1(m,1) = SUM(rmncc(:,m)*cosnv(:))
         work1(m,2) = SUM(zmnsc(:,m)*cosnv(:))
         IF (lru) work1(m,3) =-m*work1(m,1)
         IF (lzu) work1(m,4) = m*work1(m,2)
         IF (lthreed) THEN
            IF (lrv) work1(m,5) = SUM(rmncc(:,m)*sinnvn(:))
            IF (lzv) work1(m,6) = SUM(zmnsc(:,m)*sinnvn(:))
            work1(m,7) = SUM(rmnss(:,m)*sinnv(:))
            work1(m,8) = SUM(zmncs(:,m)*sinnv(:))
            IF (lru) work1(m,9) = m*work1(m,7)
            IF (lzu) work1(m,10) =-m*work1(m,8)
            IF (lrv) work1(m,11) = SUM(rmnss(:,m)*cosnvn(:))
            IF (lzv) work1(m,12) = SUM(zmncs(:,m)*cosnvn(:))
         END IF

      END DO

!
!     NEXT, INVERSE TRANSFORM IN M-U SPACE
!
      IF (lthreed) THEN
         r11 = SUM(work1(:,1)*cosmu(:) + work1(:,7)*sinmu(:))
         z11 = SUM(work1(:,2)*sinmu(:) + work1(:,8)*cosmu(:))
         IF (lru) ru = SUM(work1(:,3)*sinmu(:) + work1(:,9)*cosmu(:))
         IF (lzu) zu = SUM(work1(:,4)*cosmu(:) + work1(:,10)*sinmu(:))
         IF (lrv) rv = SUM(work1(:,5)*cosmu(:) + work1(:,11)*sinmu(:))
         IF (lzv) zv = SUM(work1(:,6)*sinmu(:) + work1(:,12)*cosmu(:))
      ELSE          
         r11 = SUM(work1(:,1)*cosmu(:))
         z11 = SUM(work1(:,2)*sinmu(:))
         IF (lru) ru = SUM(work1(:,3)*sinmu(:))
         IF (lzu) zu = SUM(work1(:,4)*cosmu(:))
      END IF


      IF (.not.lasym) GOTO 1000

      rmnsc = wmins*rzl_array(jslo,:,:,rsc) 
     1      + wplus*rzl_array(jshi,:,:,rsc)        !!SIN(mu) COS(nv)
      zmncc = wmins*rzl_array(jslo,:,:,zcc) 
     1      + wplus*rzl_array(jshi,:,:,zcc)        !!COS(mu) COS(nv)

      IF (lthreed) THEN
         rmncs = wmins*rzl_array(jslo,:,:,rcs) 
     1         + wplus*rzl_array(jshi,:,:,rcs)     !!COS(mu) SIN(nv)
         zmnss = wmins*rzl_array(jslo,:,:,zss)
     1         + wplus*rzl_array(jshi,:,:,zss)     !!SIN(mu) SIN(nv)
      END IF


!
!     COMPUTE R11, Z11 IN REAL SPACE
!
!     FIRST, INVERSE TRANSFORM IN N-V SPACE, FOR FIXED M
!
      DO m = 0, mpol1
 
         work1(m,1) = SUM(rmnsc(:,m)*cosnv(:))
         work1(m,2) = SUM(zmncc(:,m)*cosnv(:))
         IF (lru) work1(m,3) = m*work1(m,1)
         IF (lzu) work1(m,4) =-m*work1(m,2)

         IF (lthreed) THEN
            IF (lrv) work1(m,5) = SUM(rmnsc(:,m)*sinnvn(:))
            IF (lzv) work1(m,6) = SUM(zmncc(:,m)*sinnvn(:))
            work1(m,7) = SUM(rmncs(:,m)*sinnv(:))
            work1(m,8) = SUM(zmnss(:,m)*sinnv(:))
            IF (lru) work1(m,9) =-m*work1(m,7)
            IF (lzu) work1(m,10) = m*work1(m,8)
            IF (lrv) work1(m,11) = SUM(rmncs(:,m)*cosnvn(:))
            IF (lzv) work1(m,12) = SUM(zmnss(:,m)*cosnvn(:))
         END IF

      END DO

!
!     NEXT, INVERSE TRANSFORM IN M-U SPACE
!
      IF (lthreed) THEN
         r11 = r11 + SUM(work1(:,1)*sinmu(:) + work1(:,7)*cosmu(:))
         z11 = z11 + SUM(work1(:,2)*cosmu(:) + work1(:,8)*sinmu(:))
         IF (lru) ru = SUM(work1(:,3)*cosmu(:) + work1(:,9)*sinmu(:))
         IF (lzu) zu = SUM(work1(:,4)*sinmu(:) + work1(:,10)*cosmu(:))
         IF (lrv) rv = SUM(work1(:,5)*sinmu(:) + work1(:,11)*cosmu(:))
         IF (lzv) zv = SUM(work1(:,6)*cosmu(:) + work1(:,12)*sinmu(:))
      ELSE          
         r11 = r11 + SUM(work1(:,1)*sinmu(:))
         z11 = z11 + SUM(work1(:,2)*cosmu(:))
         IF (lru) ru = ru + SUM(work1(:,3)*cosmu(:))
         IF (lzu) zu = zu + SUM(work1(:,4)*sinmu(:))
      END IF

 1000 CONTINUE

      r_cyl(1) = r11;  r_cyl(3) = z11

      END SUBROUTINE flx2cyl_local

!=====================================================================================
      SUBROUTINE cyl2flx_local(rzl_in, r_cyl, c_flx, ns_in, ntor_in, 
     1      mpol_in, ntmax_in, lthreed_in, lasym_in, info, nfe, fmin, 
     1      mscale, nscale, ru, zu, rv, zv)
!------------------------------------------------------------------------------------
!  Function: converts cylin coords to flux coords
!
!  This version created 11/6/06 by J. Shields.  It is blatantly plagerized from S. Hirshman's
!     cyl2flx() routine in module vmec_utils.
! 
!     LOCAL PARAMETERS:
!     ftol    :   nominally, set to 1.E-16. Gives a maximum (relative)
!                 error in matching R and Z of sqrt(ftol), or 1.E-8. 
!                 To increase accuracy, ftol should be lowered, but this 
!                 may require more Newton iterations (slows code).
!       
!     INPUT:
!     rzl_in  :   4D array with r,z (lambda) Fourier coefficeints vs. radius
!                 
!     r_cyl   :   array specifying cylindrical point to match, R = r_cyl(1), 
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
C-----------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(out)       :: info, nfe
      INTEGER, INTENT(in)        :: ns_in, ntor_in, mpol_in, ntmax_in
      REAL(rprec), INTENT(in)    :: r_cyl(3)
      REAL(rprec), INTENT(inout) :: c_flx(3)
      REAL(rprec), INTENT(in), TARGET :: 
     1                 rzl_in(ns_in,0:ntor_in,0:mpol_in-1,2*ntmax_in)
      REAL(rprec), TARGET, OPTIONAL :: 
     1                 mscale(0:mpol_in-1), nscale(0:ntor_in)
      REAL(rprec), INTENT(out), OPTIONAL :: ru, zu, rv, zv
      REAL(rprec), INTENT(out)   :: fmin
      LOGICAL, INTENT(in) :: lthreed_in, lasym_in
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: nvar = 2
      REAL(rprec), PARAMETER :: ftol = 1.e-16_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec) :: xc_opt(nvar), r_cyl_out(3), fmin0
      INTEGER     :: iflag, itry, nfe_out

!     Initialize global variables
      rzl_ptr => rzl_in
      lthreed_loc = lthreed_in;  lasym_loc = lasym_in
      mpol_loc = mpol_in;  ntor_loc = ntor_in
      ns_loc = ns_in; ntmax_loc = ntmax_in
      lscale = PRESENT(mscale)
      IF (lscale) THEN
         mscale_loc => mscale;  nscale_loc => nscale
      END IF
      r_target = r_cyl(1); phi_target = r_cyl(2);  z_target = r_cyl(3)

!     Initialize local variables
      xc_opt(1) = c_flx(1); xc_opt(2) = c_flx(2)

!     Avoid exact magnetic axis, which is singular point
      IF (c_flx(1) .eq. zero) xc_opt(1) = one/(ns_loc-1)   

      fnorm = r_target**2 + z_target**2
      IF (fnorm .lt. EPSILON(fnorm)) fnorm = 1
      fnorm = one/fnorm


      nfe = 0
      fmin0 = 1

      DO itry = 1, 4

         CALL newt2d(xc_opt, fmin, ftol, nfe_out, nvar, info)
         nfe = nfe + nfe_out

         IF (fmin.le.ftol .or. info.eq.-3) EXIT
!
!        JOG POINT (BY ROTATING ANGLE) TO IMPROVE CONVERGENCE
!
         IF (fmin .gt. 1.E-3*fmin0) THEN
            xc_opt(2) = xc_opt(2) + twopi/20
         ELSE 
            xc_opt(2) = xc_opt(2) + twopi/40
         END IF

         fmin0 = MIN(fmin, fmin0)
!        PRINT *,' ITRY = ', itry+1,' FMIN = ', fmin
            
      END DO
         
      c_flx(1) = xc_opt(1); c_flx(2) = xc_opt(2); c_flx(3) = phi_target
      IF (info.eq.0 .and. c_flx(1).gt.one) c_flx(1) = one

      c_flx(2) = MOD(c_flx(2), twopi)
      DO WHILE (c_flx(2) .lt. zero)
         c_flx(2) = c_flx(2) + twopi
      END DO

!
!     COMPUTE Ru, Zu, Rv, Zv IF REQUIRED
!
      IF ((PRESENT(ru) .or. PRESENT(zu) .or. 
     1     PRESENT(rv) .or. PRESENT(zv)) .and. info.eq.0)
     2    CALL flx2cyl_local(rzl_in, c_flx, r_cyl_out, ns_loc, ntor_loc, 
     3         mpol_loc, ntmax_loc, lthreed_loc, lasym_loc, 
     4         iflag, MSCALE=mscale, NSCALE=nscale, 
     5         RU=ru, ZU=zu, RV=rv, ZV=zv)

    
      END SUBROUTINE cyl2flx_local


!====================================================================================
      SUBROUTINE newt2d(xc_opt, fmin, ftol, nfe, nvar, iflag)
!-------------------------------------------------------------------------------------
!
! FUNCTION:  FINDS FLUX COORDINATES (s,u) WHICH CORRESPOND TO ZERO OF TARGET FUNCTION
!               F == (R - R_TARGET)**2 + (Z - Z_TARGET)**2
!              FOR A GIVEN CYLINDRICAL COORDINATE POINT (R_TARGET, N*PHI=PHI_TARGET, Z_TARGET)
!
!     Reference:  S.E.Attenberger, W.A.Houlberg, S.P.Hirshman, J. Comp Phys 72 (1987) 435.
!
! LOGIC: The algorithm used here modifies this slightly to improve "faltering" convergence
!        by choosing a steepest-descent path when the step size has been decreased sufficiently
!        without yielding a lower value of F.
!
! INPUT/OUTPUT:
!     xc_opt:   sflux = xc_opt(1), uflux = xc_opt(2) are the toroidal flux
!               coordinate and poloidal angle coordinate, respectively
!     iflag:    = 0, successfully find s,u point
!               =-1, did not converge
!               =-3, sflux > 1, probably
!
!     LOCAL VARIABLES:
!     tau:      d(R,Z)/d(s,u) (Jacobian)
!     isgt1:    counter for number of times s>1
!
!
!
!  This version created 11/6/06 by J. Shields.  It is blatantly plagerized from S. Hirshman's
!     newt2d() routine in module vmec_utils.  (In fact, it's probably an exact copy). 
!--------------------------------------------------------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nvar
      INTEGER, INTENT(out) :: nfe, iflag
      REAL(rprec), INTENT(inout) :: xc_opt(nvar)
      REAL(rprec), INTENT(in)    :: ftol
      REAL(rprec), INTENT(out)   :: fmin
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, PARAMETER :: niter = 50
      INTEGER     :: ieval, isgt1
      REAL(rprec) :: c_flx(3), r_cyl_out(3), fvec(nvar), sflux, 
     1               uflux, eps0, eps, epu, xc_min(2), factor
      REAL(rprec) :: x0(3), xs(3), xu(3), dels, delu, tau, fmin0,
     1               ru1, zu1, edge_value

      iflag = -1      
      eps0 = SQRT(EPSILON(eps))
      xc_min = xc_opt

      c_flx(3) = phi_target
      fmin0 = 1.E10_dp
      factor = 1
      nfe = 0
      edge_value = one + one/(ns_loc-1)
      isgt1 = 0

      DO ieval = 1, niter
         nfe = nfe + 1

         sflux = MAX(xc_opt(1), zero)
!        sflux = MIN(MAX(xc_opt(1), zero), one)
         uflux = xc_opt(2)
         c_flx(1) = sflux;  c_flx(2) = uflux

!        COMPUTE R,Z, Ru, Zu
         CALL get_flxcoord(x0, c_flx, ru=ru1, zu=zu1)
         xu(1) = ru1; xu(3) = zu1

!        MAKE SURE sflux IS LARGE ENOUGH
!        TO COMPUTE d(sqrt(s))/ds ACCURATELY NEAR ORIGIN
         IF (sflux .ge. 1000*eps0) THEN
            eps = eps0
         ELSE
            eps = eps0*sflux
         END IF

!        COMPUTE Rs, Zs NUMERICALLY
         eps = ABS(eps)
         IF (sflux .ge. 1-eps) eps = -eps
         c_flx(1) = sflux + eps
         CALL get_flxcoord(r_cyl_out, c_flx)
         xs = (r_cyl_out - x0)/eps
         c_flx(1) = sflux

         x0(1) = x0(1) - r_target
         x0(3) = x0(3) - z_target
         fmin = (x0(1)**2 + x0(3)**2)*fnorm

         IF (fmin .gt. fmin0) THEN
            factor = (2*factor)/3
            xc_opt = xc_min
!           REDIRECT ALONG STEEPEST-DESCENT PATH
            IF (6*factor .lt. one) THEN
               dels =-(x0(1)*xs(1) + x0(3)*xs(3))/(xs(1)**2 + xs(3)**2)
               delu =-(x0(1)*xu(1) + x0(3)*xu(3))/(xu(1)**2 + xu(3)**2)
            END IF
         ELSE
            fmin0 = fmin
            factor = 1
            xc_min = xc_opt

!           NEWTON STEP
            tau = xu(1)*xs(3) - xu(3)*xs(1)
            IF (ABS(tau) .le. ABS(eps)*r_target**2) THEN
               iflag = -2
               EXIT
            END IF
            dels = ( x0(1)*xu(3) - x0(3)*xu(1))/tau
            delu = (-x0(1)*xs(3) + x0(3)*xs(1))/tau
            IF (fmin .gt. 1.E-3_dp) THEN
               dels = dels/2; delu = delu/2
            END IF
 
         END IF
 
         IF (fmin .le. ftol) EXIT

         IF (ABS(dels) .gt. one)   dels = SIGN(one, dels)
!         IF (ABS(delu) .gt. twopi/2) delu = SIGN(twopi/2, delu)

         xc_opt(1) = xc_opt(1) + dels*factor
         IF (xc_opt(1) .lt. zero) THEN
            xc_opt(1) = -xc_opt(1)/2               !Prevents oscillations around origin s=0
            xc_opt(2) = xc_opt(2) + twopi/2 
            delu = -delu
         END IF
         xc_opt(2) = xc_opt(2) + delu*factor

         IF (xc_opt(1) .gt. edge_value) THEN
            isgt1 = isgt1+1
            IF (xc_opt(1) .gt. 2._dp) isgt1 = isgt1+1
            IF (isgt1 .gt. 5) EXIT
         END IF

      END DO

      IF (isgt1.gt.5) THEN
         iflag = -3
         xc_min = xc_opt
      ELSE IF (xc_min(1) .gt. edge_value) THEN
         iflag = -3
      ELSE IF (fmin0 .le. ftol) THEN
         iflag = 0
      ELSE
         iflag = -1
      END IF
      
      fmin = fmin0
      xc_opt = xc_min
      xc_opt(2) = MOD(xc_opt(2), twopi)

      END SUBROUTINE newt2d

!=============================================================================================
      SUBROUTINE get_flxcoord(x1, c_flx, ru, zu)
!--------------------------------------------------------------------------------------------
!  FUNCTION: interface for flx2cyl_local routine.
!
! **This version created 11/6/06 by J. Shields.  It is blatantly plagerized from S. Hirshman's
!     get_flxcoord() routine in module vmec_utils.  (In fact, it's probably an exact copy). 
!------------------------------------------------------------------------------------

C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), INTENT(out) :: x1(3)
      REAL(rprec), INTENT(in)  :: c_flx(3)
      REAL(rprec), INTENT(out), OPTIONAL :: ru, zu
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: iflag
C-----------------------------------------------
      IF (lscale) THEN
         CALL flx2cyl_local(rzl_ptr, c_flx, x1, ns_loc, ntor_loc, 
     1              mpol_loc, ntmax_loc, lthreed_loc, lasym_loc, iflag, 
     2              MSCALE=mscale_loc, NSCALE=nscale_loc, RU=ru, ZU=zu)
      ELSE
         CALL flx2cyl_local(rzl_ptr, c_flx, x1, ns_loc, ntor_loc, 
     1              mpol_loc, ntmax_loc, lthreed_loc, lasym_loc, iflag, 
     2              RU=ru, ZU=zu)
      END IF

      END SUBROUTINE get_flxcoord

!=================================================================================
      SUBROUTINE tosuvspace_local_v3fit(s_in, u_in, v_in, gsqrt,                     &
     &     bsupu, bsupv)
!--------------------------------------------------------------------------------
!
!  FUNCTION:  COMPUTE VARIOUS HALF-RADIAL GRID QUANTITIES AT THE INPUT POINT
!            (S, U, V) , WHERE 
!        S = normalized toroidal flux (0 - 1),
!        U = poloidal angle 
!        V = N*phi = toroidal angle * no. field periods
!
!  This version created 8/7/07 by J. Shields.  It is shamelessly plagerized from S. Hirshman's
!    tosuvspace() routine in module read_wout_mod.  Basically, the only changes involve explicitly copying
!    read_wout_mod variables into local ones.  I also made the local arrays explicitly ALLOCATABLE.
!
!
! module global inputs: eq_rzl  ! derived TYPE that includes eq_state info, etc
!
!  
!  NOTE: this version will gets its info from the (module-global) TYPE variable "eq_rzl" instead
!             of directly from the wout (as is done in tosuvspace_local_wout).  Hence, it does NOT
!             need a "USE read_wout_mod" call.  Obviously, it used 
!             tosuvspace_local_wout as a template.
!
!
!------------------------------------------------------------------------------------
      USE stel_constants, ONLY: zero, one
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), INTENT(in) :: s_in, u_in, v_in
      REAL(rprec), INTENT(out), OPTIONAL :: gsqrt, bsupu, bsupv
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: c1p5 = 1.5_dp
      INTEGER :: m, n, n1, mn, ipresent, jslo, jshi
      INTEGER :: mnmax_nyq_w, mnyq_w, nnyq_w, ns_w, nfp_w    ! local versions of read_wout_mod variables
      REAL(rprec) :: hs1, wlo, whi, wlo_odd, whi_odd
!      REAL(rprec), DIMENSION(mnmax_nyq) :: gmnc1, gmns1, bsupumnc1,
!     1   bsupumns1, bsupvmnc1, bsupvmns1, wmins, wplus
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: gmnc1, gmns1, bsupumnc1,       &
     &             bsupumns1, bsupvmnc1, bsupvmns1, wmins, wplus

      REAL(rprec) :: cosu, sinu, cosv, sinv, tcosmn, tsinmn, sgn
!      REAL(rprec) :: cosmu(0:mnyq), sinmu(0:mnyq),
!     1               cosnv(0:nnyq), sinnv(0:nnyq)
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: cosmu, sinmu,                  &
     &             cosnv, sinnv
      LOGICAL :: lgsqrt, lbsupu, lbsupv
      LOGICAL :: lasym_w                                     ! local version of read_wout lasym variable
      CHARACTER(len=*), PARAMETER :: subname ='tosuvspace_local_v3fit: '
C-----------------------------------------------

!      write(*,*) 'now in ', subname

      IF (s_in.lt.zero .or. s_in.gt.one) THEN
         WRITE(6, *)
     1   ' In tosuvspace, s(flux) must be between 0 and 1'
         RETURN
      END IF


!.........copy wout variables to their local equivalents.......!
      ns_w        = eq_rzl%ns
      nfp_w       = eq_rzl%nfp
      mnyq_w      = eq_rzl%mnyq
      nnyq_w      = eq_rzl%nnyq
      mnmax_nyq_w = eq_rzl%mnmax_nyq
      lasym_w     = eq_rzl%lasym

      
!...........allocate local arrays............................!

      ALLOCATE( gmnc1(mnmax_nyq_w), gmns1(mnmax_nyq_w) )
      ALLOCATE( bsupumnc1(mnmax_nyq_w), bsupumns1(mnmax_nyq_w) )
      ALLOCATE( bsupvmnc1(mnmax_nyq_w), bsupvmns1(mnmax_nyq_w) )
      ALLOCATE( wmins(mnmax_nyq_w), wplus(mnmax_nyq_w) )

      ALLOCATE( cosmu(0:mnyq_w), sinmu(0:mnyq_w) )
      ALLOCATE( cosnv(0:nnyq_w), sinnv(0:nnyq_w) )


!
!     FIND INTERPOLATED s VALUE AND COMPUTE INTERPOLATION WEIGHTS wlo, whi
!     RECALL THAT THESE QUANTITIES ARE ON THE HALF-RADIAL GRID...
!     s-half(j) = (j-1.5)*hs, for j = 2,...ns_w
!
      hs1 = one/(ns_w-1)
      jslo = INT(c1p5 + s_in/hs1)
      jshi = jslo+1
      wlo = (hs1*(jshi-c1p5) - s_in)/hs1
      whi = 1 - wlo
      IF (jslo .eq. ns_w) THEN
!        USE Xhalf(ns_w+1) = 2*Xhalf(ns_w) - Xhalf(ns_w-1) FOR "GHOST" POINT VALUE 1/2hs OUTSIDE EDGE
!        THEN, X = wlo*Xhalf(ns_w) + whi*Xhalf(ns_w+1) == Xhalf(ns_w) + whi*(Xhalf(ns_w) - Xhalf(ns_w-1)) 
         jshi = jslo-1
         wlo = 1+whi; whi = -whi
      ELSE IF (jslo .eq. 1) THEN
         jslo = 2
      END IF

!
!     FOR ODD-m MODES X ~ SQRT(s), SO INTERPOLATE Xmn/SQRT(s)
! 
      whi_odd = whi*SQRT(s_in/(hs1*(jshi-c1p5)))
      IF (jslo .ne. 1) THEN
         wlo_odd = wlo*SQRT(s_in/(hs1*(jslo-c1p5)))
      ELSE
         wlo_odd = 0
         whi_odd = SQRT(s_in/(hs1*(jshi-c1p5)))
      END IF

      WHERE (MOD(NINT(eq_rzl%xm_nyq(:)),2) .eq. 0)
         wmins = wlo
         wplus = whi
      ELSEWHERE
         wmins = wlo_odd
         wplus = whi_odd
      END WHERE

      ipresent = 0
      lgsqrt = PRESENT(gsqrt)
      IF (lgsqrt) THEN
         gsqrt = 0 ;  ipresent = ipresent+1
         gmnc1 = wmins*eq_rzl%gmnc(:,jslo) +                              &
     &             wplus*eq_rzl%gmnc(:,jshi)
         IF (lasym_w)
     1   gmns1 = wmins*eq_rzl%gmns(:,jslo) +                              &
     &             wplus*eq_rzl%gmns(:,jshi)
      END IF
      lbsupu = PRESENT(bsupu)
      IF (lbsupu) THEN
         bsupu = 0 ;  ipresent = ipresent+1
         bsupumnc1 = wmins*eq_rzl%bsupumnc(:,jslo) +                      &
     &                 wplus*eq_rzl%bsupumnc(:,jshi)
         IF (lasym_w)
     1   bsupumns1 = wmins*eq_rzl%bsupumns(:,jslo) +                      &
     &                 wplus*eq_rzl%bsupumns(:,jshi)
      END IF
      lbsupv = PRESENT(bsupv)
      IF (lbsupv) THEN
         bsupv = 0 ;  ipresent = ipresent+1
         bsupvmnc1 = wmins*eq_rzl%bsupvmnc(:,jslo) +                      &
     &                 wplus*eq_rzl%bsupvmnc(:,jshi)
         IF (lasym_w)
     1   bsupvmns1 = wmins*eq_rzl%bsupvmns(:,jslo) +                      &
     &                 wplus*eq_rzl%bsupvmns(:,jshi)
      END IF

      IF (ipresent .eq. 0) RETURN

!
!     SETUP TRIG ARRAYS
!
      cosu = COS(u_in);   sinu = SIN(u_in)
      cosv = COS(v_in);   sinv = SIN(v_in)

      cosmu(0) = 1;    sinmu(0) = 0
      cosnv(0) = 1;    sinnv(0) = 0
      DO m = 1, mnyq_w
         cosmu(m) = cosmu(m-1)*cosu - sinmu(m-1)*sinu
         sinmu(m) = sinmu(m-1)*cosu + cosmu(m-1)*sinu
      END DO

      DO n = 1, nnyq_w
         cosnv(n) = cosnv(n-1)*cosv - sinnv(n-1)*sinv
         sinnv(n) = sinnv(n-1)*cosv + cosnv(n-1)*sinv
      END DO

!
!     COMPUTE GSQRT, ... IN REAL SPACE
!     tcosmn = cos(mu - nv);  tsinmn = sin(mu - nv)
!
      DO mn = 1, mnmax_nyq_w
         m = NINT(eq_rzl%xm_nyq(mn))
         n = NINT(eq_rzl%xn_nyq(mn))/nfp_w
         n1 = ABS(n);   sgn = SIGN(1,n)
         tcosmn = cosmu(m)*cosnv(n1) + sgn*sinmu(m)*sinnv(n1)   
         IF (lgsqrt) gsqrt = gsqrt + gmnc1(mn)*tcosmn
         IF (lbsupu) bsupu = bsupu + bsupumnc1(mn)*tcosmn
         IF (lbsupv) bsupv = bsupv + bsupvmnc1(mn)*tcosmn
      END DO

      IF (.not.lasym_w) GOTO 1000

      DO mn = 1, mnmax_nyq_w
         m = NINT(eq_rzl%xm_nyq(mn))
         n = NINT(eq_rzl%xn_nyq(mn))/nfp_w
         n1 = ABS(n);   sgn = SIGN(1,n)
         tsinmn = sinmu(m)*cosnv(n1) - sgn*cosmu(m)*sinnv(n1)
         IF (lgsqrt) gsqrt = gsqrt + gmns1(mn)*tsinmn
         IF (lbsupu) bsupu = bsupu + bsupumns1(mn)*tsinmn
         IF (lbsupv) bsupv = bsupv + bsupvmns1(mn)*tsinmn
      END DO

 1000 CONTINUE

      END SUBROUTINE tosuvspace_local_v3fit

!=================================================================================
      SUBROUTINE tosuvspace_local_wout(s_in, u_in, v_in, gsqrt,                     &
     &     bsupu, bsupv)
!--------------------------------------------------------------------------------
!
!  FUNCTION:  COMPUTE VARIOUS HALF-RADIAL GRID QUANTITIES AT THE INPUT POINT
!            (S, U, V) , WHERE 
!        S = normalized toroidal flux (0 - 1),
!        U = poloidal angle 
!        V = N*phi = toroidal angle * no. field periods
!
!  This version created 12/7/06 by J. Shields.  It is shamelessly plagerized from S. Hirshman's
!    tosuvspace() routine in module read_wout_mod.  Basically, the only changes involve explicitly copying
!    read_wout_mod variables into local ones.  I also made the local arrays explicitly ALLOCATABLE.
!
! IMPORTANT NOTE: this routine requires the use of GLOBAL read_wout_mod variables
!                 and hence is NOT independent of that module!!!  JS 12/7/06
!                     (on the bright side, at least it doesn't call any other
!                        subroutines....)
!  
!-----------------------------------------------
      USE stel_constants, ONLY: zero, one
      USE read_wout_mod, ONLY: gmnc, gmns, bsupumnc, bsupumns,                 &
     &    bsupvmnc, bsupvmns, lwout_opened, xm_nyq, xn_nyq
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), INTENT(in) :: s_in, u_in, v_in
      REAL(rprec), INTENT(out), OPTIONAL :: gsqrt, bsupu, bsupv
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: c1p5 = 1.5_dp
      INTEGER :: m, n, n1, mn, ipresent, jslo, jshi
      INTEGER :: mnmax_nyq_w, mnyq_w, nnyq_w, ns_w, nfp_w    ! local versions of read_wout_mod variables
      REAL(rprec) :: hs1, wlo, whi, wlo_odd, whi_odd
!      REAL(rprec), DIMENSION(mnmax_nyq) :: gmnc1, gmns1, bsupumnc1,
!     1   bsupumns1, bsupvmnc1, bsupvmns1, wmins, wplus
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: gmnc1, gmns1, bsupumnc1,       &
     &             bsupumns1, bsupvmnc1, bsupvmns1, wmins, wplus

      REAL(rprec) :: cosu, sinu, cosv, sinv, tcosmn, tsinmn, sgn
!      REAL(rprec) :: cosmu(0:mnyq), sinmu(0:mnyq),
!     1               cosnv(0:nnyq), sinnv(0:nnyq)
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: cosmu, sinmu,                  &
     &             cosnv, sinnv
      LOGICAL :: lgsqrt, lbsupu, lbsupv
      LOGICAL :: lasym_w                                     ! local version of read_wout lasym variable
      CHARACTER(len=*), PARAMETER :: subname = 'tosuvspace_local_wout: '
C-----------------------------------------------

!      write(*,*) 'now in ', subname

      IF (s_in.lt.zero .or. s_in.gt.one) THEN
         WRITE(6, *)
     1   ' In tosuvspace, s(flux) must be between 0 and 1'
         RETURN
      END IF

      IF (.not.lwout_opened) THEN
         WRITE(6, *)
     1   ' tosuvspace can only be called AFTER opening wout file!'
         RETURN
      END IF

!.........copy wout variables to their local equivalents.......!
      ns_w        = eq_rzl%ns
      nfp_w       = eq_rzl%nfp
      mnyq_w      = eq_rzl%mnyq
      nnyq_w      = eq_rzl%nnyq
      mnmax_nyq_w = eq_rzl%mnmax_nyq
      lasym_w     = eq_rzl%lasym

      
!...........allocate local arrays............................!

      ALLOCATE( gmnc1(mnmax_nyq_w), gmns1(mnmax_nyq_w) )
      ALLOCATE( bsupumnc1(mnmax_nyq_w), bsupumns1(mnmax_nyq_w) )
      ALLOCATE( bsupvmnc1(mnmax_nyq_w), bsupvmns1(mnmax_nyq_w) )
      ALLOCATE( wmins(mnmax_nyq_w), wplus(mnmax_nyq_w) )

      ALLOCATE( cosmu(0:mnyq_w), sinmu(0:mnyq_w) )
      ALLOCATE( cosnv(0:nnyq_w), sinnv(0:nnyq_w) )


!
!     FIND INTERPOLATED s VALUE AND COMPUTE INTERPOLATION WEIGHTS wlo, whi
!     RECALL THAT THESE QUANTITIES ARE ON THE HALF-RADIAL GRID...
!     s-half(j) = (j-1.5)*hs, for j = 2,...ns_w
!
      hs1 = one/(ns_w-1)
      jslo = INT(c1p5 + s_in/hs1)
      jshi = jslo+1
      wlo = (hs1*(jshi-c1p5) - s_in)/hs1
      whi = 1 - wlo
      IF (jslo .eq. ns_w) THEN
!        USE Xhalf(ns_w+1) = 2*Xhalf(ns_w) - Xhalf(ns_w-1) FOR "GHOST" POINT VALUE 1/2hs OUTSIDE EDGE
!        THEN, X = wlo*Xhalf(ns_w) + whi*Xhalf(ns_w+1) == Xhalf(ns_w) + whi*(Xhalf(ns_w) - Xhalf(ns_w-1)) 
         jshi = jslo-1
         wlo = 1+whi; whi = -whi
      ELSE IF (jslo .eq. 1) THEN
         jslo = 2
      END IF

!
!     FOR ODD-m MODES X ~ SQRT(s), SO INTERPOLATE Xmn/SQRT(s)
! 
      whi_odd = whi*SQRT(s_in/(hs1*(jshi-c1p5)))
      IF (jslo .ne. 1) THEN
         wlo_odd = wlo*SQRT(s_in/(hs1*(jslo-c1p5)))
      ELSE
         wlo_odd = 0
         whi_odd = SQRT(s_in/(hs1*(jshi-c1p5)))
      END IF

      WHERE (MOD(NINT(xm_nyq(:)),2) .eq. 0)
         wmins = wlo
         wplus = whi
      ELSEWHERE
         wmins = wlo_odd
         wplus = whi_odd
      END WHERE

      ipresent = 0
      lgsqrt = PRESENT(gsqrt)
      IF (lgsqrt) THEN
         gsqrt = 0 ;  ipresent = ipresent+1
         gmnc1 = wmins*gmnc(:,jslo) + wplus*gmnc(:,jshi)
         IF (lasym_w)
     1   gmns1 = wmins*gmns(:,jslo) + wplus*gmns(:,jshi)
      END IF
      lbsupu = PRESENT(bsupu)
      IF (lbsupu) THEN
         bsupu = 0 ;  ipresent = ipresent+1
         bsupumnc1 = wmins*bsupumnc(:,jslo) + wplus*bsupumnc(:,jshi)
         IF (lasym_w)
     1   bsupumns1 = wmins*bsupumns(:,jslo) + wplus*bsupumns(:,jshi)
      END IF
      lbsupv = PRESENT(bsupv)
      IF (lbsupv) THEN
         bsupv = 0 ;  ipresent = ipresent+1
         bsupvmnc1 = wmins*bsupvmnc(:,jslo) + wplus*bsupvmnc(:,jshi)
         IF (lasym_w)
     1   bsupvmns1 = wmins*bsupvmns(:,jslo) + wplus*bsupvmns(:,jshi)
      END IF

      IF (ipresent .eq. 0) RETURN

!
!     SETUP TRIG ARRAYS
!
      cosu = COS(u_in);   sinu = SIN(u_in)
      cosv = COS(v_in);   sinv = SIN(v_in)

      cosmu(0) = 1;    sinmu(0) = 0
      cosnv(0) = 1;    sinnv(0) = 0
      DO m = 1, mnyq_w
         cosmu(m) = cosmu(m-1)*cosu - sinmu(m-1)*sinu
         sinmu(m) = sinmu(m-1)*cosu + cosmu(m-1)*sinu
      END DO

      DO n = 1, nnyq_w
         cosnv(n) = cosnv(n-1)*cosv - sinnv(n-1)*sinv
         sinnv(n) = sinnv(n-1)*cosv + cosnv(n-1)*sinv
      END DO

!
!     COMPUTE GSQRT, ... IN REAL SPACE
!     tcosmn = cos(mu - nv);  tsinmn = sin(mu - nv)
!
      DO mn = 1, mnmax_nyq_w
         m = NINT(xm_nyq(mn));  n = NINT(xn_nyq(mn))/nfp_w
         n1 = ABS(n);   sgn = SIGN(1,n)
         tcosmn = cosmu(m)*cosnv(n1) + sgn*sinmu(m)*sinnv(n1)   
         IF (lgsqrt) gsqrt = gsqrt + gmnc1(mn)*tcosmn
         IF (lbsupu) bsupu = bsupu + bsupumnc1(mn)*tcosmn
         IF (lbsupv) bsupv = bsupv + bsupvmnc1(mn)*tcosmn
      END DO

      IF (.not.lasym_w) GOTO 1000

      DO mn = 1, mnmax_nyq_w
         m = NINT(xm_nyq(mn));  n = NINT(xn_nyq(mn))/nfp_w
         n1 = ABS(n);   sgn = SIGN(1,n)
         tsinmn = sinmu(m)*cosnv(n1) - sgn*cosmu(m)*sinnv(n1)
         IF (lgsqrt) gsqrt = gsqrt + gmns1(mn)*tsinmn
         IF (lbsupu) bsupu = bsupu + bsupumns1(mn)*tsinmn
         IF (lbsupv) bsupv = bsupv + bsupvmns1(mn)*tsinmn
      END DO

 1000 CONTINUE

      END SUBROUTINE tosuvspace_local_wout


!==============================================================
      SUBROUTINE LoadRZL_loc_wout(ntmax, lthreed)
!-------------------------------------------------------------------
!
! FUNCTION: loads R,Z data from wout files into arrays for B field transforms
!
! This version created 12/7/06 by J. Shields.  It is blatantly plagerized from S. Hirshman's
!      LoadRZL routine in module read_wout_mod
!
! std inputs: (none)
! std outputs:  ntmax, lthreed
! global inputs:
! module global outputs: rzl_loc
!
!------------------------------------------------------------------------
      USE read_wout_mod, ONLY: ns, xm, xn, mnmax, mpol, lasym,                 &
     &      nfp, ntor, rmns, rmnc, zmns, zmnc

      IMPLICIT NONE

!.........dummy variables.................................!
      INTEGER, INTENT(OUT) :: ntmax
      LOGICAL,        INTENT(OUT) :: lthreed

C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER     :: rcc, rss, zsc, zcs, rsc, rcs, zcc, zss
      INTEGER     :: mpol1, mn, m, n, n1
      REAL(rprec) :: sgn
C-----------------------------------------------
!
!     Arrays must be stacked (and ns,ntor,mpol ordering imposed)
!     as coefficients of cos(mu)*cos(nv), etc
!     Only need R, Z components(not lambda, for now anyhow)
!      
      IF (ALLOCATED(rzl_loc)) RETURN

      lthreed = ANY(NINT(xn) .ne. 0)
      mpol1 = mpol-1
      rcc = 1;  zsc = 1
      IF (.not.lasym) THEN
         IF (lthreed) THEN
            ntmax = 2
            rss = 2;  zcs = 2
         ELSE
            ntmax = 1
         END IF
      ELSE
         IF (lthreed) THEN
            ntmax = 4
            rss = 2;  rsc = 3;  rcs = 4
            zcs = 2;  zcc = 3;  zss = 4
         ELSE
            ntmax = 2
            rsc = 2;  zcc = 2
         END IF
      END IF

      zsc = 1+ntmax; zcs = zcs+ntmax; zcc = zcc+ntmax
      ALLOCATE(rzl_loc(ns,0:ntor,0:mpol1,2*ntmax))
      rzl_loc = 0

      DO mn = 1, mnmax
         m = NINT(xm(mn));  n = NINT(xn(mn))/nfp; n1 = ABS(n)
         sgn = SIGN(1, n)
         rzl_loc(:,n1,m,rcc) = rzl_loc(:,n1,m,rcc) + rmnc(mn,:)
         rzl_loc(:,n1,m,zsc) = rzl_loc(:,n1,m,zsc) + zmns(mn,:)
         IF (lthreed) THEN
            rzl_loc(:,n1,m,rss) = rzl_loc(:,n1,m,rss) 
     1                            + sgn*rmnc(mn,:)
            rzl_loc(:,n1,m,zcs) = rzl_loc(:,n1,m,zcs)
     1                            - sgn*zmns(mn,:)
         END IF
         IF (lasym) THEN
            rzl_loc(:,n1,m,rsc) = rzl_loc(:,n1,m,rsc) 
     1                            + rmns(mn,:)
            rzl_loc(:,n1,m,zcc) = rzl_loc(:,n1,m,zcc) 
     1                            + zmnc(mn,:)
            IF (lthreed) THEN
                rzl_loc(:,n1,m,rcs) = rzl_loc(:,n1,m,rcs)
     1                                - sgn*rmns(mn,:)
                rzl_loc(:,n1,m,zss) = rzl_loc(:,n1,m,zss) 
     1                                + sgn*zmnc(mn,:)
            END IF
         END IF
      END DO


      END SUBROUTINE LoadRZL_loc_wout


!-------------------------------------------------------------------------------
      SUBROUTINE eq_rzlam_construct_v3fit(this, ntmax_in, lthreed_in)
!----------------------------------------------------------------
!
! FUNCTION: Uses wout variables to assign values to the (module-global) derived TYPE variable "eq_rzl".
!
!
! INPUTS:  wout variables from S. Hirshman's read_wout_mod
!          (module global) temporary array "rzl_loc"
!
! OUTPUT: module-global derived TYPE variable "eq_rzl"
!
! created 8/7/07 by J. Shields  (usingeq_aux_1_construct as a template)
!
! **called by VMEC_B_INIT 
!
!---------------------------------------------------------------------------
      USE read_wout_mod, only : read_wout_file, rmnc, zmns,                    &
     &   gmnc, lasym, rmns,                                                    &
     &   zmnc, gmns,                                                           &
     &   xm, xn, xm_nyq, xn_nyq, ns, mnmax, mnmax_nyq,                         &
     &   bsupumnc, bsupvmnc, bsupumns, bsupvmns,                               &
     &   lasym, mpol, nfp, ntor, lwout_opened

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (eq_rzlam), INTENT (inout)       :: this
      INTEGER, INTENT(IN) :: ntmax_in
      LOGICAL,        INTENT(IN) :: lthreed_in
      
!  Declare local variables
      INTEGER nwout, ierr
      INTEGER, DIMENSION(3)   :: dims
      INTEGER :: mnyq, nnyq
      LOGICAL        :: lthreed_w, lasym_w
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'eq_rzlam_construct: '

!  Start of executable code


      IF (.not.lwout_opened) THEN
         WRITE(6, '(2a,/,a)')
     1   ' This form of GetBcyl can only be called if WOUT has been',
     2   ' previously opened!',' Try GetBcyl_VMEC form instead.'
         RETURN
      END IF


!  Clear out this
      CALL eq_rzlam_destroy(this)



      mnyq = INT(MAXVAL(xm_nyq))
      nnyq = INT(MAXVAL(ABS(xn_nyq)))/nfp 
 
!........copy read_wout_mod global variables into local versions............!
      this%ns        = ns
      this%nfp       = nfp
      this%mnmax     = mnmax
      this%mnmax_nyq = mnmax_nyq
      this%mnyq      = mnyq
      this%nnyq      = nnyq
      this%mpol      = mpol
      this%ntmax     = ntmax_in
      this%ntor      = ntor
      this%lthreed   = lthreed_in
      this%lasym     = lasym

!.........allocate this%rzl_local..............................................!
      ALLOCATE( this%rzl_local(ns,0:ntor,0:mpol-1,2*ntmax_in) )

!.........set eq_rzl%rzl_local to the (module global) rzl_loc array that was filled by LoadRZL_loc_wout.......!
      this%rzl_local = rzl_loc

      
!  Copy wout arrays in to this
      dims(1:2) = shape(rmnc)
      ALLOCATE(this % rmnc(dims(1),dims(2)))
      this % rmnc = rmnc
      dims(1:2) = shape(zmns)
      ALLOCATE(this % zmns(dims(1),dims(2)))
      this % zmns = zmns
      dims(1:2) = shape(gmnc)
      ALLOCATE(this % gmnc(dims(1),dims(2)))
      this % gmnc = gmnc

      dims(1:2) = shape(bsupumnc)
      ALLOCATE(this % bsupumnc(dims(1),dims(2)))
      this % bsupumnc = bsupumnc
      dims(1:2) = shape(bsupvmnc)
      ALLOCATE(this % bsupvmnc(dims(1),dims(2)))
      this % bsupvmnc = bsupvmnc

      dims(1:1) = shape(xm)
      ALLOCATE(this % xm(dims(1)))
      this % xm = xm
      dims(1:1) = shape(xn)
      ALLOCATE(this % xn(dims(1)))
      this % xn = xn
      dims(1:1) = shape(xm_nyq)
      ALLOCATE(this % xm_nyq(dims(1)))
      this % xm_nyq = xm_nyq
      dims(1:1) = shape(xn_nyq)
      ALLOCATE(this % xn_nyq(dims(1)))
      this % xn_nyq = xn_nyq


!      IF (lasym) THEN
         dims(1:2) = shape(rmns)
         ALLOCATE(this % rmns(dims(1),dims(2)))
         this % rmns = rmns
         dims(1:2) = shape(zmnc)
         ALLOCATE(this % zmnc(dims(1),dims(2)))
         this % zmnc = zmnc
         dims(1:2) = shape(gmns)
         ALLOCATE(this % gmns(dims(1),dims(2)))
         this % gmns = gmns

         dims(1:2) = shape(bsupumns)
         ALLOCATE(this % bsupumns(dims(1),dims(2)))
         this % bsupumns = bsupumns
         dims(1:2) = shape(bsupvmns)
         ALLOCATE(this % bsupvmns(dims(1),dims(2)))
         this % bsupvmns = bsupvmns
!      END IF

      
      RETURN
      
      END SUBROUTINE eq_rzlam_construct_v3fit


!====================================================================================
      SUBROUTINE fill_eq_rzl_wout(ntmax_in, lthreed_in)
!----------------------------------------------------------------
!
! FUNCTION: Uses wout variables to assign values to the (module-global) derived TYPE variable "eq_rzl".
!            It also allocates and initializes (but does *NOT* fill) the module-global rzl_local array.
!
!                **(TEMPORARY NOTE: actually, right now we're *ONLY* using the wout for testing,
!                    so it DOES, in fact, currently fill the rzl_local array!! JS 12/8/06)
!
! INPUT:  wout variables from S. Hirshman's read_wout_mod
! OUTPUT: module-global derived TYPE variable "eq_rzl"
!
! created 12/8/06 by J. Shields
!
! **called by VMEC_B_INIT 
!
!---------------------------------------------------------------------------
      USE read_wout_mod, ONLY: ns, nfp, mnmax, mnmax_nyq,                       &
     &                         xm_nyq, xn_nyq, mpol, ntor,                      &
     &                         lasym, lwout_opened
      IMPLICIT NONE

!..........dummy variables.................!
      INTEGER, INTENT(IN) :: ntmax_in
      LOGICAL,        INTENT(IN) :: lthreed_in

!......L o c a l   V a r i a b l e s.......................!
      INTEGER :: mnyq, nnyq
      LOGICAL        :: lthreed_w, lasym_w
      CHARACTER(len=*), PARAMETER :: subname = 'fill_eq_rzl_wout: '
C-----------------------------------------------
      IF (.not.lwout_opened) THEN
         WRITE(6, '(2a,/,a)')
     1   ' This form of GetBcyl can only be called if WOUT has been',
     2   ' previously opened!',' Try GetBcyl_VMEC form instead.'
         RETURN
      END IF


      mnyq = INT(MAXVAL(xm_nyq))
      nnyq = INT(MAXVAL(ABS(xn_nyq)))/nfp 
 
!........copy read_wout_mod global variables into local versions............!
      eq_rzl%ns        = ns
      eq_rzl%nfp       = nfp
      eq_rzl%mnmax     = mnmax
      eq_rzl%mnmax_nyq = mnmax_nyq
      eq_rzl%mnyq      = mnyq
      eq_rzl%nnyq      = nnyq
      eq_rzl%mpol      = mpol
      eq_rzl%ntmax     = ntmax_in
      eq_rzl%ntor      = ntor
      eq_rzl%lthreed   = lthreed_in
      eq_rzl%lasym     = lasym

!.........allocate eq_rzl%rzl_local
      ALLOCATE( eq_rzl%rzl_local(ns,0:ntor,0:mpol-1,2*ntmax_in) )

!.........set eq_rzl%rzl_local to the (module global) rzl_loc array that was filled by LoadRZL_loc_wout.......!
      eq_rzl%rzl_local = rzl_loc



      END SUBROUTINE fill_eq_rzl_wout


!===============================================================================
      SUBROUTINE eq_rzlam_destroy(this)
!-------------------------------------------------------------------
! FUNCTION: destroys a variable of TYPE eq_rzlam
!
! created 8/7/07 by J. Shields  (usingeq_aux_1_destroy as a template)
!--------------------------------------------------------------------
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (eq_rzlam), INTENT (inout)       :: this

!  Start of executable code

      this%ns        = 0
      this%nfp       = 0
      this%mnmax     = 0
      this%mnmax_nyq = 0
      this%mnyq      = 0
      this%nnyq      = 0
      this%mpol      = 0
      this%ntmax     = 0
      this%ntor      = 0

      IF (ASSOCIATED(this % rmnc)) THEN
!         DEALLOCATE(this % rmnc, this % zmns, this % lmns,                     &
!     &      this % gmnc, this % currumnc, this % currvmnc,                     &
!     &      this % bsubumnc, this % bsubvmnc,                                  &
!     &      this % bsupumnc, this % bsupvmnc)

         DEALLOCATE(this % rmnc, this % zmns,                                  &
     &      this % gmnc,                                                       &
     &      this % bsupumnc, this % bsupvmnc)
      END IF


      IF (ASSOCIATED(this % rmns)) THEN
!         DEALLOCATE(this % rmns, this % zmnc, this % lmnc,                     &
!     &      this % gmns, this % currumns, this % currvmns,                     &
!     &      this % bsubumns, this % bsubvmns,                                  &
!     &      this % bsupumns, this % bsupvmns)

         DEALLOCATE(this % rmns, this % zmnc,                                  &
     &      this % gmns,                                                       &
     &      this % bsupumns, this % bsupvmns)
      END IF

      IF (ASSOCIATED(this % xm)) THEN
         DEALLOCATE(this % xm, this % xn, this % xm_nyq,                       &
     &      this % xn_nyq)
      END IF

!      IF (ASSOCIATED(this % phi)) THEN
!         DEALLOCATE(this % phi, this % iotaf)
!      END IF

      IF (ASSOCIATED(this % rzl_local)) THEN
         DEALLOCATE(this % rzl_local)
      END IF

      END SUBROUTINE eq_rzlam_destroy


      END MODULE b_transform_vmec
