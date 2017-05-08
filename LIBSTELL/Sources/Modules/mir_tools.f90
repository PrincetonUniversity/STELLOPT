! ------------------------------------------------------------------------------
!     Module:        mir_tools
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!                    N. Pablant (npablant@pppl.gov)
!     Date:          2017-04-26
!
!     Description:
!       This module is designed as a replacement to the  AJAX stellarator
!       utilities.  It provides routines for the inversion of stellerator 
!       equilibira.
!
!       This is rewrite of the stel_tools modules designed to have a more
!       generalized interface.
!
!     Initialization:
!       This module must be initialized before being used.  This consists
!       of loading a wout file and initialzing the module.
!
!       Optionally the splines can be initialized as well by calling:
!         initialize_splines(error_status, nu, nv)
!       or calling one or more of the individual spline initialization
!       routines. The optional inputs nu, nv control the number of
!       points in the poloidal and toroidal used to construct the splines.
!
!       If the splines are not initialized manually, the will be
!       initialized on the fly when needed, using the default values
!       of nu and nv.
!
!
!     Programming Notes:
!       Where the code takes v_val as an input it is asking for a value
!       running from 0 to 2*pi over a field period.  However,
!       where the code takes phi_val as an input it is asking for
!       the real toroidal angle, 0 to 2*pi over the device.  So
!       phi_val = v_val/nfp.
!
!       VMEC calculates some quanties (b-field, jacobian) on a half integer
!       mesh rather than on a full integer mesh.  Where possible the
!       spline in this module are constructed by directly using the
!       VMEC output, however in a few cases it is necessary to linearly
!       interpolate from half mesh -> full mesh, and calculate everything
!       using full mesh quantities instead.
!
!       When programming this I made the choise to allow nu and nv
!       to be different for each spline initialization routine.
!       Because of this, it is no longer possible to reuse the value arrays
!       calculated for one set of splines in the initialization for another
!       set of splines.  This means that the value arrays must always be
!       recalculated, potentially doubling initialization time.
!       This only affects the fsa and susceptance spline initialization.
!
!     Error Status:
!       A error_status of zero means success.
!       Postive values of error_status are informational but may indicate
!       a non-terminal error of some type.
!
!       Here are a list of module wide specific error codes:
!      
!        0  : Success
!       -1  : Error (General/non-specific)
!       -2  : Outside of flux domain
!       -3  : Wout file not loaded.
!       -4  : Module variables not loaded.
!       -5  : Calculation not valid at s=0.
!       -10 : Not enough radial surfaces in input arrays.
!       -11 : Not enough fourier modes in input arrays.
!
!       -100 to -199 : Error from read_wout_file. See read_wout_mod.f for
!                      information about specific codes.
!       -200 to -399 : Error from EZSpline.  See EZspline.f90 for information
!                      about specific codes.
!
!     Programming ToDo:
!       * Fix handling of error codes from EZSpline to be consistent with
!         error_status scheme defined above.
!       * I should go back to making error_status an inout variable to make
!         stringing together commands without error checking for each line.
!       * I sholud split up the fsa spline initialization.
!
!     Programming Notes:
!        Loads a given set of Fourier variables into memory.  The
!        B-Field Fourier Harmonics must be supplied if the user wishes
!        to evaluate the field.  The Lambda variables are only needed
!        if calling the vmec2pest routine.  Assmues Kernel of the form
!        (mu+nv).
!        K1     First radial index (usually 1)
!        K2     Last radial index
!        MNMAX  Number of Fourier modes
!        NU     Number of real space poloidal gridpoints
!        NV     Number of real space toroidal gridpoints
!        XM     Array of poloidal modes (1:MNMAX)
!        XN     Array of toroidal modes (1:MNMAX)
!        IFLAG  Error flag
!        RMNC   Array of R modes (1:MNMAX,K1:K2) (cos)
!        ZMNS   Array of Z modes (1:MNMAX,K1:K2) (sin)
!        (OPTIONAL VARS)
!        RMNS   Array of R modes (1:MNMAX,K1:K2) (sin)
!        ZMNC   Array of Z modes (1:MNMAX,K1:K2) (cos)
!        BSMNS  Array of B^s modes (1:MNMAX,K1:K2) (sin)
!        BUMNC  Array of B^u modes (1:MNMAX,K1:K2) (cos)
!        BVMNC  Array of B^v modes (1:MNMAX,K1:K2) (cos)
!        BSMNC  Array of B^s modes (1:MNMAX,K1:K2) (cos)
!        BUMNS  Array of B^u modes (1:MNMAX,K1:K2) (sin)
!        BVMNS  Array of B^v modes (1:MNMAX,K1:K2) (sin)
!        LMNS   Array of Lambda modes (1:MNMAX,K1:K2) (sin)
!        LMNC   Array of Lambda modes (1:MNMAX,K1:K2) (cos)
!
! ------------------------------------------------------------------------------
      MODULE mir_tools
      
        !-----------------------------------------------------------------------
        !     Libraries
        !-----------------------------------------------------------------------
        USE EZspline_obj
        USE read_wout_mod, ONLY: &
             ns_vmec => ns &
             ,rmnc_vmec => rmnc &
             ,rmns_vmec => rmns &
             ,zmnc_vmec => zmnc &
             ,gmns_vmec => gmns &
             ,gmnc_vmec => gmnc &
             ,lmns_vmec => lmns &
             ,lmnc_vmec => lmnc &
             ,bmns_vmec => bmns &
             ,bmnc_vmec => bmnc &
             ,zmns_vmec => zmns &
             ,mnmax_vmec => mnmax &
             ,bsupumnc_vmec => bsupumnc &
             ,bsupvmnc_vmec => bsupvmnc &
             ,bsupumns_vmec => bsupumns &
             ,bsupvmns_vmec => bsupvmns &
             ,xm_vmec => xm &
             ,xn_vmec => xn &
             ,lasym_vmec => lasym &
             ,mpol_vmec => mpol &
             ,ntor_vmec => ntor &
             ,nfp_vmec => nfp &
             ,read_wout_file &
             ,read_wout_deallocate
      
        !-----------------------------------------------------------------------
        !     Module Variables
        !     PUBLIC
        !        LINTSTEPS     Number of discrete steps for line integrals (256)
        !                      
        !     PRIVATE
        !        DOMAIN_FLAG     Used to determin if in domain durring search
        !        NFP             Field periodicity determined from XN array
        !        R(PHI/Z)_target Used durring search
        !        BCS(0/1)        Spline Boundary Conditions
        !        PI2             2*pi
        !        X_SPL           Spline Variables
        !-----------------------------------------------------------------------
        IMPLICIT NONE
        INTEGER, PARAMETER ::  lintsteps=512
        INTEGER, PARAMETER, PRIVATE ::  neval_max=10000
        INTEGER, PRIVATE      ::  domain_flag, nfp
        DOUBLE PRECISION, PRIVATE ::  R_target, PHI_target, Z_target
        INTEGER, PRIVATE     ::  bcs0(2) = (/ 0, 0/)
        INTEGER, PRIVATE     ::  bcs1(2) = (/-1,-1/)
        DOUBLE PRECISION, PARAMETER, PRIVATE      :: pi2 = 6.283185482025146D+00
        DOUBLE PRECISION, PARAMETER, PRIVATE      :: one = 1.000000000000000D+00
        DOUBLE PRECISION, PARAMETER, PRIVATE      :: search_tol = 1.000000000000000D-12
        TYPE(EZspline1_r8), PRIVATE :: Vp_spl, b_fsa_spl, grho_fsa_spl, grho2_fsa_spl, b2overgrho_fsa_spl
        TYPE(EZspline1_r8), PRIVATE :: S11_spl, S12_spl, S21_spl, S22_spl
        TYPE(EZspline3_r8), PRIVATE :: R_spl, Z_spl, G_spl
        TYPE(EZspline3_r8), PRIVATE :: Ru_spl, Zu_spl
        TYPE(EZspline3_r8), PRIVATE :: Rv_spl, Zv_spl
        TYPE(EZspline3_r8), PRIVATE :: Bs_spl, Bu_spl, Bv_spl, B_spl
        TYPE(EZspline3_r8), PRIVATE :: L_spl, Lu_spl, Lv_spl
        
        DOUBLE PRECISION, ALLOCATABLE, PRIVATE :: &
          rmnc(:,:) &
          ,zmns(:,:) &
          ,rmns(:,:) &
          ,zmnc(:,:) &
          ,bsmns_half(:,:) &
          ,bsupumnc_half(:,:) &
          ,bsupvmnc_half(:,:) &
          ,bsmnc_half(:,:) &
          ,bsupumns_half(:,:) &
          ,bsupvmns_half(:,:) &
          ,lmns_half(:,:) &
          ,lmnc_half(:,:) &
          ,bmnc_half(:,:) &
          ,bmns_half(:,:) &
          ,gmnc_half(:,:) &
          ,gmns_half(:,:) &
          ,bsmns_full(:,:) &
          ,bsupumnc_full(:,:) &
          ,bsupvmnc_full(:,:) &
          ,bsmnc_full(:,:) &
          ,bsupumns_full(:,:) &
          ,bsupvmns_full(:,:) &
          ,lmns_full(:,:) &
          ,lmnc_full(:,:) &
          ,bmnc_full(:,:) &
          ,bmns_full(:,:) &
          ,gmnc_full(:,:) &
          ,gmns_full(:,:)
        
        INTEGER, ALLOCATABLE, PRIVATE :: &
          xm(:) &
          ,xn(:)

        LOGICAL, PRIVATE :: &
          lasym
             
        INTEGER, PRIVATE :: &
          k1 &
          ,k2 &
          ,mnmax &
          ,isherm &
          ,ns_t
             
        !-----------------------------------------------------------------------
        !     Private Subroutines
        !-----------------------------------------------------------------------
        PRIVATE :: mntouv
        
        !-----------------------------------------------------------------------
        !     INTERFACE Modules
        !-----------------------------------------------------------------------
        INTERFACE guess_flx_from_cyl
           MODULE PROCEDURE guess_flx_from_cyl_dbl, guess_flx_from_cyl_sgl
        END INTERFACE guess_flx_from_cyl
      CONTAINS


        SUBROUTINE initialize (error_status, filename)
          ! Initialize the mir_tools.
          !
          ! Optionally, this will also load also a given wout file.

          IMPLICIT NONE

          INTEGER, INTENT(OUT) :: error_status          
          CHARACTER(len=*), INTENT(IN), OPTIONAL :: filename

          IF (PRESENT(filename)) THEN
             CALL read_wout_file (TRIM(filename), error_status)
             IF (error_status .NE. 0) THEN
                error_status = -100 - ABS(error_status)
                RETURN
             ENDIF
          ENDIF
          
          CALL initialize_module(error_status)

        END SUBROUTINE initialize

        
        SUBROUTINE initialize_module (error_status)
          IMPLICIT NONE
          INTEGER, INTENT(OUT) :: error_status
          INTEGER :: mn
          
          error_status = 0
          
          k1 = 1
          k2 = ns_vmec
          mnmax = mnmax_vmec
          
          ns_t = k2-k1+1
          isherm = 0
          
          ! Preform checks
          IF (.NOT. ALLOCATED(xm_vmec)) THEN
             error_status = -3
             RETURN
          ENDIF
          
          IF (ns_t < 1) error_status = -10
          IF (mnmax < 1) error_status = -11
          IF (error_status < 0) RETURN
          
          ! Allocations
          IF (ALLOCATED(xm)) CALL mirtools_deallocate(error_status)
          IF (error_status < 0) RETURN
          
          ALLOCATE( &
            xm(1:mnmax_vmec) &
            ,xn(1:mnmax_vmec) &
            ,rmnc(1:mnmax_vmec,1:ns_vmec) &
            ,zmns(1:mnmax_vmec,1:ns_vmec) &
            ,rmns(1:mnmax_vmec,1:ns_vmec) &
            ,zmnc(1:mnmax_vmec,1:ns_vmec) &
            ,bsmns_half(1:mnmax_vmec,1:ns_vmec) &
            ,bsupumnc_half(1:mnmax_vmec,1:ns_vmec) &
            ,bsupvmnc_half(1:mnmax_vmec,1:ns_vmec) &
            ,bsmnc_half(1:mnmax_vmec,1:ns_vmec) &
            ,bsupumns_half(1:mnmax_vmec,1:ns_vmec) &
            ,bsupvmns_half(1:mnmax_vmec,1:ns_vmec) &
            ,lmns_half(1:mnmax_vmec,1:ns_vmec) &
            ,lmnc_half(1:mnmax_vmec,1:ns_vmec) &
            ,bmnc_half(1:mnmax_vmec,1:ns_vmec) &
            ,bmns_half(1:mnmax_vmec,1:ns_vmec) &
            ,gmnc_half(1:mnmax_vmec,1:ns_vmec) &
            ,gmns_half(1:mnmax_vmec,1:ns_vmec) &
            ,bsmns_full(1:mnmax_vmec,1:ns_vmec) &
            ,bsupumnc_full(1:mnmax_vmec,1:ns_vmec) &
            ,bsupvmnc_full(1:mnmax_vmec,1:ns_vmec) &
            ,bsmnc_full(1:mnmax_vmec,1:ns_vmec) &
            ,bsupumns_full(1:mnmax_vmec,1:ns_vmec) &
            ,bsupvmns_full(1:mnmax_vmec,1:ns_vmec) &
            ,lmns_full(1:mnmax_vmec,1:ns_vmec) &
            ,lmnc_full(1:mnmax_vmec,1:ns_vmec) &
            ,bmnc_full(1:mnmax_vmec,1:ns_vmec) &
            ,bmns_full(1:mnmax_vmec,1:ns_vmec) &
            ,gmnc_full(1:mnmax_vmec,1:ns_vmec) &
            ,gmns_full(1:mnmax_vmec,1:ns_vmec) &
            )          

          ! Copy arrays from the the READ_WOUT_MOD module.
          lasym = lasym_vmec
          
          rmnc(:,:) = rmnc_vmec(:,:)
          zmns(:,:) = zmns_vmec(:,:)

          bsmns_half(:,:) = 0
          bsupumnc_half(:,:) = bsupumnc_vmec(:,:)
          bsupvmnc_half(:,:) = bsupvmnc_vmec(:,:)

          lmns_half(:,:) = lmns_vmec(:,:)
          bmnc_half(:,:) = bmnc_vmec(:,:)
          gmnc_half(:,:) = gmnc_vmec(:,:)

          IF (lasym_vmec) THEN
             rmns(:,:) = rmns_vmec(:,:)
             zmnc(:,:) = zmnc_vmec(:,:)
             
             bsmnc_half(:,:) = 0                         
             bsupumns_half(:,:) = bsupumns_vmec(:,:)
             bsupvmns_half(:,:) = bsupvmns_vmec(:,:)
             
             lmnc_half(:,:) = lmnc_vmec(:,:)
             bmns_half(:,:) = bmns_vmec(:,:)
             gmns_half(:,:) = gmns_vmec(:,:)
          ENDIF

          ! Convert half mesh values to full mesh values.

          ! bsmns, bsmnc                   
          bsmns_full(:,1)     = 1.5*bsmns_half(:,2) - 0.5*bsmns_half(:,3)        
          FORALL(mn = 1:mnmax_vmec) bsmns_full(mn,2:ns_vmec-1) = &
               0.5*(bsmns_half(mn,2:ns_vmec-1) + bsmns_half(mn,3:ns_vmec))
          bsmns_full(:,ns_vmec) = 1.5*bsmns_half(:,ns_vmec) - 0.5*bsmns_half(:,ns_vmec-1)
          
          IF (lasym_vmec) THEN
             bsmnc_full(:,1) = 1.5*bsmnc_half(:,2) - 0.5*bsmnc_half(:,3)
             FORALL(mn = 1:mnmax_vmec) bsmnc_full(mn,2:ns_vmec-1)     = &
                  0.5*(bsmnc_half(mn,2:ns_vmec-1) + bsmnc_half(mn,3:ns_vmec))
             bsmnc_full(:,ns_vmec) = 1.5*bsmnc_half(:,ns_vmec) - 0.5*bsmnc_half(:,ns_vmec-1)
          ENDIF
          
          ! bsupumnc, bsupumns 
          bsupumnc_full(:,1)     = 1.5*bsupumnc_half(:,2) - 0.5*bsupumnc_half(:,3)        
          FORALL(mn = 1:mnmax_vmec) bsupumnc_full(mn,2:ns_vmec-1) = &
               0.5*(bsupumnc_half(mn,2:ns_vmec-1) + bsupumnc_half(mn,3:ns_vmec))
          bsupumnc_full(:,ns_vmec) = 1.5*bsupumnc_half(:,ns_vmec) - 0.5*bsupumnc_half(:,ns_vmec-1)
          
          IF (lasym_vmec) THEN
             bsupumns_full(:,1) = 1.5*bsupumns_half(:,2) - 0.5*bsupumns_half(:,3)
             FORALL(mn = 1:mnmax_vmec) bsupumns_full(mn,2:ns_vmec-1)     = &
                  0.5*(bsupumns_half(mn,2:ns_vmec-1) + bsupumns_half(mn,3:ns_vmec))
             bsupumns_full(:,ns_vmec) = 1.5*bsupumns_half(:,ns_vmec) - 0.5*bsupumns_half(:,ns_vmec-1)
          ENDIF
          
          ! bsupvmnc, bsupvmns                    
          bsupvmnc_full(:,1)     = 1.5*bsupvmnc_half(:,2) - 0.5*bsupvmnc_half(:,3)        
          FORALL(mn = 1:mnmax_vmec) bsupvmnc_full(mn,2:ns_vmec-1) = &
               0.5*(bsupvmnc_half(mn,2:ns_vmec-1) + bsupvmnc_half(mn,3:ns_vmec))
          bsupvmnc_full(:,ns_vmec) = 1.5*bsupvmnc_half(:,ns_vmec) - 0.5*bsupvmnc_half(:,ns_vmec-1)
          
          IF (lasym_vmec) THEN
             bsupvmns_full(:,1) = 1.5*bsupvmns_half(:,2) - 0.5*bsupvmns_half(:,3)
             FORALL(mn = 1:mnmax_vmec) bsupvmns_full(mn,2:ns_vmec-1)     = &
                  0.5*(bsupvmns_half(mn,2:ns_vmec-1) + bsupvmns_half(mn,3:ns_vmec))
             bsupvmns_full(:,ns_vmec) = 1.5*bsupvmns_half(:,ns_vmec) - 0.5*bsupvmns_half(:,ns_vmec-1)
          ENDIF
          
          ! bmnc, bmns
          bmnc_full(:,1)     = 1.5*bmnc_half(:,2) - 0.5*bmnc_half(:,3)        
          FORALL(mn = 1:mnmax_vmec) bmnc_full(mn,2:ns_vmec-1) = &
               0.5*(bmnc_half(mn,2:ns_vmec-1) + bmnc_half(mn,3:ns_vmec))
          bmnc_full(:,ns_vmec) = 1.5*bmnc_half(:,ns_vmec) - 0.5*bmnc_half(:,ns_vmec-1)
          
          IF (lasym_vmec) THEN
             bmns_full(:,1) = 1.5*bmns_half(:,2) - 0.5*bmns_half(:,3)
             FORALL(mn = 1:mnmax_vmec) bmns_full(mn,2:ns_vmec-1)     = &
                  0.5*(bmns_half(mn,2:ns_vmec-1) + bmns_half(mn,3:ns_vmec))
             bmns_full(:,ns_vmec) = 1.5*bmns_half(:,ns_vmec) - 0.5*bmns_half(:,ns_vmec-1)
          ENDIF
          
          ! gmnc, gmns
          gmnc_full(:,1)     = 1.5*gmnc_half(:,2) - 0.5*gmnc_half(:,3)        
          FORALL(mn = 1:mnmax_vmec) gmnc_full(mn,2:ns_vmec-1) = &
               0.5*(gmnc_half(mn,2:ns_vmec-1) + gmnc_half(mn,3:ns_vmec))
          gmnc_full(:,ns_vmec) = 1.5*gmnc_half(:,ns_vmec) - 0.5*gmnc_half(:,ns_vmec-1)
          
          IF (lasym_vmec) THEN
             gmns_full(:,1) = 1.5*gmns_half(:,2) - 0.5*gmns_half(:,3)
             FORALL(mn = 1:mnmax_vmec) gmns_full(mn,2:ns_vmec-1)     = &
                  0.5*(gmns_half(mn,2:ns_vmec-1) + gmns_half(mn,3:ns_vmec))
             gmns_full(:,ns_vmec) = 1.5*gmns_half(:,ns_vmec) - 0.5*gmns_half(:,ns_vmec-1)
          ENDIF

          ! lmns, lmnc
          lmns_full(:,1)     = 1.5*lmns_half(:,2) - 0.5*lmns_half(:,3)        
          FORALL(mn = 1:mnmax_vmec) lmns_full(mn,2:ns_vmec-1) = &
               0.5*(lmns_half(mn,2:ns_vmec-1) + lmns_half(mn,3:ns_vmec))
          lmns_full(:,ns_vmec) = 1.5*lmns_half(:,ns_vmec) - 0.5*lmns_half(:,ns_vmec-1)
          
          IF (lasym_vmec) THEN
             lmnc_full(:,1) = 1.5*lmnc_half(:,2) - 0.5*lmnc_half(:,3)
             FORALL(mn = 1:mnmax_vmec) lmnc_full(mn,2:ns_vmec-1)     = &
                  0.5*(lmnc_half(mn,2:ns_vmec-1) + lmnc_half(mn,3:ns_vmec))
             lmnc_full(:,ns_vmec) = 1.5*lmnc_half(:,ns_vmec) - 0.5*lmnc_half(:,ns_vmec-1)
          ENDIF
          
          ! Find NFP
          !
          ! I don't know why sam did this rather than just use nfp_vmec
          !   --novi.
          nfp = 1
          nfp = MINVAL(ABS(xn_vmec), MASK=(xn_vmec>0))
          IF (nfp == 0) nfp = 1

          !set xn and xm for this module.
          xn = -1 * INT(xn_vmec) / nfp
          xm = INT(xm_vmec)
          
        END SUBROUTINE initialize_module

        
        SUBROUTINE mirtools_deallocate (error_status)
          INTEGER, INTENT(OUT) :: error_status

          error_status = 0
          
          DEALLOCATE( &
            rmnc &
            ,zmns &
            ,rmns &
            ,zmnc &
            ,bsmns_half &
            ,bsupumnc_half &
            ,bsupvmnc_half &
            ,bsmnc_half &
            ,bsupumns_half &
            ,bsupvmns_half &
            ,lmns_half &
            ,lmnc_half &
            ,bmnc_half &
            ,bmns_half &
            ,gmnc_half &
            ,gmns_half &
            ,bsmns_full &
            ,bsupumnc_full &
            ,bsupvmnc_full &
            ,bsmnc_full &
            ,bsupumns_full &
            ,bsupvmns_full &
            ,lmns_full &
            ,lmnc_full &
            ,bmnc_full &
            ,bmns_full &
            ,gmnc_full &
            ,gmns_full &
            )
          
          DEALLOCATE( &
            xm &
            ,xn)
               
        END SUBROUTINE mirtools_deallocate

        
        SUBROUTINE guess_nu_nv(nu, nv, error_status)
          ! These values were chosen emperically to produce accuracy
          ! in going from flux space to real space of about 1e-5 m.
          ! Testing was done for LHD and W7-X.
          
          IMPLICIT NONE
          INTEGER, INTENT(OUT) :: nu, nv
          INTEGER, INTENT(OUT) :: error_status
          error_status = 0
          nu = 50+1
          nv = 500/nfp+1
        END SUBROUTINE guess_nu_nv

        
        SUBROUTINE initialize_splines(error_status, nu, nv)
          ! Initialize all of the splines available in the module.
          
          USE EZspline
          IMPLICIT NONE
          INTEGER, INTENT(out)       :: error_status
          INTEGER, INTENT(in), OPTIONAL :: nu
          INTEGER, INTENT(in), OPTIONAL :: nv

          CALL initialize_splines_rz(error_status, nu, nv)
          CALL initialize_splines_rzderiv(error_status, nu, nv)
          CALL initialize_splines_b(error_status, nu, nv)
          CALL initialize_splines_modb(error_status, nu, nv)
          CALL initialize_splines_jacobian(error_status, nu, nv)
          CALL initialize_splines_fsa(error_status, nu, nv)
          CALL initialize_splines_lambda(error_status, nu, nv)
          CALL initialize_splines_susceptance(error_status, nu, nv)
          
        END SUBROUTINE initialize_splines

        
        SUBROUTINE initialize_splines_rz(error_status, nu, nv)  
          USE EZspline
          IMPLICIT NONE
          INTEGER, INTENT(out)       :: error_status
          INTEGER, INTENT(in), OPTIONAL :: nu
          INTEGER, INTENT(in), OPTIONAL :: nv
          INTEGER :: nu_local, nv_local
          INTEGER :: ez_status
          INTEGER ::  u, mn
          DOUBLE PRECISION, ALLOCATABLE :: xu(:), xv(:), rho(:)
          DOUBLE PRECISION, ALLOCATABLE :: f_temp(:,:,:)
          
          error_status = 0
          
          IF (PRESENT(nu) .AND. PRESENT(nv)) THEN
             nu_local = nu
             nv_local = nv
          ELSE
             CALL guess_nu_nv(nu_local, nv_local, error_status)
          END IF

          PRINT '("Inititalizing R, Z splines with nu=", i0, " nv=", i0)', nu_local, nv_local
          
          ! Check module initalization.
          IF (.NOT. ALLOCATED(xm)) THEN
             error_status = -4
             RETURN
          ENDIF
          
          !Allocations
          ALLOCATE(xu(nu_local),xv(nv_local),rho(k1:k2))
          ALLOCATE(f_temp(nu_local,nv_local,k1:k2))
          
          f_temp = 0
          FORALL(u=k1:k2) rho(u) = REAL(u-1)/REAL(ns_t-1)
          FORALL(u=1:nu_local) xu(u) = REAL(u-1)/REAL(nu_local-1)
          FORALL(u=1:nv_local) xv(u) = REAL(u-1)/REAL(nv_local-1)
          
          ! Free Memory
          IF (EZspline_allocated(R_spl)) CALL EZspline_free(R_spl,ez_status)
          IF (EZspline_allocated(Z_spl)) CALL EZspline_free(Z_spl,ez_status)
          
          ! Preform Init
          CALL EZspline_init(R_spl,nu_local,nv_local,ns_t,bcs1,bcs1,bcs0,ez_status)
          CALL EZspline_init(Z_spl,nu_local,nv_local,ns_t,bcs1,bcs1,bcs0,ez_status)
          
          ! R/Z
          R_spl%x1 = xu*pi2; R_spl%x2 = xv*pi2; R_spl%x3 = rho; R_spl%isHermite = isherm
          CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,rmnc,xm,xn,f_temp,0,1)
          IF (lasym) CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,rmns,xm,xn,f_temp,1,0)
          CALL EZspline_setup(R_spl,f_temp,ez_status); f_temp=0
          
          Z_spl%x1 = xu*pi2; Z_spl%x2 = xv*pi2; Z_spl%x3 = rho; Z_spl%isHermite = isherm
          CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,zmns,xm,xn,f_temp,1,0)
          IF (lasym) CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,zmnc,xm,xn,f_temp,0,0)
          CALL EZspline_setup(Z_spl,f_temp,ez_status); f_temp = 0

          IF (ez_status .NE. 0) error_status = -200 - ABS(ez_status)
          
          ! deallocations
          DEALLOCATE(xu,xv,rho)
          DEALLOCATE(f_temp)


        END SUBROUTINE initialize_splines_rz

        
        SUBROUTINE initialize_splines_rzderiv(error_status, nu, nv)  
          USE EZspline
          IMPLICIT NONE
          INTEGER, INTENT(out)       :: error_status
          INTEGER, INTENT(in), OPTIONAL :: nu
          INTEGER, INTENT(in), OPTIONAL :: nv
          INTEGER :: nu_local, nv_local
          INTEGER :: ez_status
          INTEGER ::  u, mn
          DOUBLE PRECISION, ALLOCATABLE :: xu(:), xv(:), rho(:)
          DOUBLE PRECISION, ALLOCATABLE :: f_temp(:,:,:)
          DOUBLE PRECISION, ALLOCATABLE :: fmn_temp(:,:)
      
          error_status = 0

          IF (PRESENT(nu) .AND. PRESENT(nv)) THEN
             nu_local = nu
             nv_local = nv
          ELSE
             CALL guess_nu_nv(nu_local, nv_local, error_status)
          END IF

          PRINT '("Inititalizing R, Z derivative splines with nu=", i0, " nv=", i0)', nu_local, nv_local
          
          ! Check module initalization.
          IF (.NOT. ALLOCATED(xm)) THEN
             error_status = -4
             RETURN
          ENDIF
          
          !Allocations
          ALLOCATE(xu(nu_local),xv(nv_local),rho(k1:k2))
          ALLOCATE(f_temp(nu_local,nv_local,k1:k2))
          ALLOCATE(fmn_temp(1:mnmax,k1:k2))
          
          f_temp = 0
          fmn_temp = 0
          FORALL(u=k1:k2) rho(u) = REAL(u-1)/REAL(ns_t-1)
          FORALL(u=1:nu_local) xu(u) = REAL(u-1)/REAL(nu_local-1)
          FORALL(u=1:nv_local) xv(u) = REAL(u-1)/REAL(nv_local-1)
          
          ! Free Memory
          IF (EZspline_allocated(Ru_spl)) CALL EZspline_free(Ru_spl,ez_status)
          IF (EZspline_allocated(Rv_spl)) CALL EZspline_free(Rv_spl,ez_status)
          IF (EZspline_allocated(Zu_spl)) CALL EZspline_free(Zu_spl,ez_status)
          IF (EZspline_allocated(Zv_spl)) CALL EZspline_free(Zv_spl,ez_status)
          
          ! Preform Init
          CALL EZspline_init(Ru_spl,nu_local,nv_local,ns_t,bcs1,bcs1,bcs0,ez_status)
          CALL EZspline_init(Rv_spl,nu_local,nv_local,ns_t,bcs1,bcs1,bcs0,ez_status)
          CALL EZspline_init(Zu_spl,nu_local,nv_local,ns_t,bcs1,bcs1,bcs0,ez_status)
          CALL EZspline_init(Zv_spl,nu_local,nv_local,ns_t,bcs1,bcs1,bcs0,ez_status)
          
          ! dR/Du Derivatives
          Ru_spl%x1 = xu*pi2; Ru_spl%x2 = xv*pi2; Ru_spl%x3 = rho; Ru_spl%isHermite = isherm
          FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -rmnc(mn,:)*xm(mn)
          CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,f_temp,1,0)
          IF (lasym) THEN
             FORALL(mn = 1:mnmax) fmn_temp(mn,:) = rmns(mn,:)*xm(mn)
             CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,f_temp,0,0)
          END IF
          CALL EZspline_setup(Ru_spl,f_temp,ez_status); f_temp = 0
          ! dR/Dv Derivatives
          Rv_spl%x1 = xu*pi2; Rv_spl%x2 = xv*pi2; Rv_spl%x3 = rho; Rv_spl%isHermite = isherm
          FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -rmnc(mn,:)*xn(mn)
          CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,f_temp,1,0)
          IF (lasym) THEN
             FORALL(mn = 1:mnmax) fmn_temp(mn,:) = rmns(mn,:)*xn(mn)
             CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,f_temp,0,0)
          END IF
          CALL EZspline_setup(Rv_spl,f_temp,ez_status); f_temp = 0
          ! dZ/Du Derivatives
          Zu_spl%x1 = xu*pi2; Zu_spl%x2 = xv*pi2; Zu_spl%x3 = rho; Zu_spl%isHermite = isherm
          FORALL(mn = 1:mnmax) fmn_temp(mn,:) = zmns(mn,:)*xm(mn)
          CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,f_temp,0,0)
          IF (lasym) THEN
             FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -zmnc(mn,:)*xm(mn)
             CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,f_temp,1,0)
          END IF
          CALL EZspline_setup(Zu_spl,f_temp,ez_status); f_temp = 0
          ! dZ/Dv Derivatives
          Zv_spl%x1 = xu*pi2; Zv_spl%x2 = xv*pi2; Zv_spl%x3 = rho; Zv_spl%isHermite = isherm
          FORALL(mn = 1:mnmax) fmn_temp(mn,:) = zmns(mn,:)*xn(mn)
          CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,f_temp,0,0)
          IF (lasym) THEN
             FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -zmnc(mn,:)*xn(mn)
             CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,f_temp,1,0)
          END IF
          CALL EZspline_setup(Zv_spl,f_temp,ez_status); f_temp = 0
         
          IF (ez_status .NE. 0) error_status = -200 - ABS(ez_status) 
          
          ! deallocations
          DEALLOCATE(xu,xv,rho)
          DEALLOCATE(f_temp)
          DEALLOCATE(fmn_temp)

        END SUBROUTINE initialize_splines_rzderiv

        
        SUBROUTINE initialize_splines_b(error_status, nu, nv)  
          USE EZspline
          IMPLICIT NONE
          INTEGER, INTENT(out)       :: error_status
          INTEGER, INTENT(in), OPTIONAL :: nu
          INTEGER, INTENT(in), OPTIONAL :: nv
          INTEGER :: nu_local, nv_local
          INTEGER :: ez_status
          INTEGER ::  u, mn
          DOUBLE PRECISION, ALLOCATABLE :: xu(:), xv(:), rho(:)
          DOUBLE PRECISION, ALLOCATABLE :: f_temp(:,:,:)
      
          error_status = 0

          IF (PRESENT(nu) .AND. PRESENT(nv)) THEN
             nu_local = nu
             nv_local = nv
          ELSE
             CALL guess_nu_nv(nu_local, nv_local, error_status)
          END IF

          PRINT '("Inititalizing b-field splines with nu=", i0, " nv=", i0)', nu_local, nv_local
          
          ! Check module initalization.
          IF (.NOT. ALLOCATED(xm)) THEN
             error_status = -4
             RETURN
          ENDIF
          
          !Allocations
          ALLOCATE(xu(nu_local),xv(nv_local),rho(k1:k2))
          ALLOCATE(f_temp(nu_local,nv_local,k1:k2))
          
          f_temp = 0
          rho(:) = 0.0
          FORALL(u=k1+1:k2) rho(u) = REAL(u-1)/REAL(ns_t-1) - 0.5/REAL(ns_t-1)
          FORALL(u=1:nu_local) xu(u) = REAL(u-1)/REAL(nu_local-1)
          FORALL(u=1:nv_local) xv(u) = REAL(u-1)/REAL(nv_local-1)

          ! Find the values at s=0 through linear interpolation.
          rho(k1) = 0.0
          bsmnc_half(:,1) = 1.5*bsmnc_half(:,2) - 0.5*bsmnc_half(:,3)
          bsmns_half(:,1) = 1.5*bsmns_half(:,2) - 0.5*bsmns_half(:,3)
          bsupumnc_half(:,1) = 1.5*bsupumnc_half(:,2) - 0.5*bsupumnc_half(:,3)
          bsupvmnc_half(:,1) = 1.5*bsupvmnc_half(:,2) - 0.5*bsupvmnc_half(:,3)
          bsupumns_half(:,1) = 1.5*bsupumns_half(:,2) - 0.5*bsupumns_half(:,3)
          bsupvmns_half(:,1) = 1.5*bsupvmns_half(:,2) - 0.5*bsupvmns_half(:,3)
             
          ! Free Memory
          IF (EZspline_allocated(Bs_spl)) CALL EZspline_free(Bs_spl,ez_status)
          IF (EZspline_allocated(Bu_spl)) CALL EZspline_free(Bu_spl,ez_status)
          IF (EZspline_allocated(Bv_spl)) CALL EZspline_free(Bv_spl,ez_status)
          
          ! Preform Init
          CALL EZspline_init(Bs_spl,nu_local,nv_local,ns_t,bcs1,bcs1,bcs0,ez_status)
          CALL EZspline_init(Bu_spl,nu_local,nv_local,ns_t,bcs1,bcs1,bcs0,ez_status)
          CALL EZspline_init(Bv_spl,nu_local,nv_local,ns_t,bcs1,bcs1,bcs0,ez_status)
          
          ! B^s
          ! In VMEC B^s is defined as zero, so we don't need to calculate this one.
          Bs_spl%x1 = xu*pi2; Bs_spl%x2 = xv*pi2; Bs_spl%x3 = rho; Bs_spl%isHermite = isherm
          !IF (lasym) CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,bsmns_half,xm,xn,f_temp,1,0)
          !IF (lasym) CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,bsmnc_half,xm,xn,f_temp,0,0)
          CALL EZspline_setup(Bs_spl,f_temp,ez_status); f_temp = 0
          ! B^u
          Bu_spl%x1 = xu*pi2; Bu_spl%x2 = xv*pi2; Bu_spl%x3 = rho; Bu_spl%isHermite = isherm
          CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,bsupumnc_half,xm,xn,f_temp,0,0)
          IF (lasym) CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,bsupumns_half,xm,xn,f_temp,1,0)
          CALL EZspline_setup(Bu_spl,f_temp,ez_status); f_temp = 0
          ! B^v
          Bv_spl%x1 = xu*pi2; Bv_spl%x2 = xv*pi2; Bv_spl%x3 = rho; Bv_spl%isHermite = isherm
          CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,bsupvmnc_half,xm,xn,f_temp,0,0)
          IF (lasym) CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,bsupvmns_half,xm,xn,f_temp,1,0)
          CALL EZspline_setup(Bv_spl,f_temp,ez_status); f_temp = 0

          
          IF (ez_status .NE. 0) error_status = -200 - ABS(ez_status) 
          
          ! deallocations
          DEALLOCATE(xu,xv,rho)
          DEALLOCATE(f_temp)

        END SUBROUTINE initialize_splines_b        

        
        SUBROUTINE initialize_splines_modb(error_status, nu, nv)  
          USE EZspline
          IMPLICIT NONE
          INTEGER, INTENT(out)       :: error_status
          INTEGER, INTENT(in), OPTIONAL :: nu
          INTEGER, INTENT(in), OPTIONAL :: nv
          INTEGER :: nu_local, nv_local
          INTEGER :: ez_status
          INTEGER ::  u, mn
          DOUBLE PRECISION, ALLOCATABLE :: xu(:), xv(:), rho(:)
          DOUBLE PRECISION, ALLOCATABLE :: f_temp(:,:,:)
      
          error_status = 0

          IF (PRESENT(nu) .AND. PRESENT(nv)) THEN
             nu_local = nu
             nv_local = nv
          ELSE
             CALL guess_nu_nv(nu_local, nv_local, error_status)
          END IF

          PRINT '("Inititalizing |b| splines with nu=", i0, " nv=", i0)', nu_local, nv_local
          
          ! Check module initalization.
          IF (.NOT. ALLOCATED(xm)) THEN
             error_status = -4
             RETURN
          ENDIF
          
          !Allocations
          ALLOCATE(xu(nu_local),xv(nv_local),rho(k1:k2))
          ALLOCATE(f_temp(nu_local,nv_local,k1:k2))
          
          f_temp = 0
          rho(:) = 0.0
          FORALL(u=k1+1:k2) rho(u) = REAL(u-1)/REAL(ns_t-1) - 0.5/REAL(ns_t-1)
          FORALL(u=1:nu_local) xu(u) = REAL(u-1)/REAL(nu_local-1)
          FORALL(u=1:nv_local) xv(u) = REAL(u-1)/REAL(nv_local-1)

          ! Find the coefficent values at s=0 through linear interpolation.
          rho(k1) = 0.0
          bmnc_half(:,1) = 1.5*bmnc_half(:,2) - 0.5*bmnc_half(:,3)
          bmns_half(:,1) = 1.5*bmns_half(:,2) - 0.5*bmns_half(:,3)
             
          ! Free Memory
          IF (EZspline_allocated(B_spl)) CALL EZspline_free(B_spl,ez_status)
          
          ! Preform Init
          CALL EZspline_init(B_spl,nu_local,nv_local,ns_t,bcs1,bcs1,bcs0,ez_status)

          ! ModB
          B_spl%x1 = xu*pi2; B_spl%x2 = xv*pi2; B_spl%x3 = rho; B_spl%isHermite = isherm
          CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,bmnc_half,xm,xn,f_temp,0,0)
          IF (lasym) CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,bmns_half,xm,xn,f_temp,1,0)
          CALL EZspline_setup(B_spl,f_temp,ez_status); f_temp = 0
          
          IF (ez_status .NE. 0) error_status = -200 - ABS(ez_status) 
          
          ! deallocations
          DEALLOCATE(xu,xv,rho)
          DEALLOCATE(f_temp)

        END SUBROUTINE initialize_splines_modb
        
        
        SUBROUTINE initialize_splines_lambda(error_status, nu, nv)  
          USE EZspline
          IMPLICIT NONE
          INTEGER, INTENT(out)       :: error_status
          INTEGER, INTENT(in), OPTIONAL :: nu
          INTEGER, INTENT(in), OPTIONAL :: nv
          INTEGER :: nu_local, nv_local
          INTEGER :: ez_status
          INTEGER ::  u, mn
          DOUBLE PRECISION, ALLOCATABLE :: xu(:), xv(:), rho(:)
          DOUBLE PRECISION, ALLOCATABLE :: f_temp(:,:,:)
          DOUBLE PRECISION, ALLOCATABLE :: fmn_temp(:,:)
   
          error_status = 0

          IF (PRESENT(nu) .AND. PRESENT(nv)) THEN
             nu_local = nu
             nv_local = nv
          ELSE
             CALL guess_nu_nv(nu_local, nv_local, error_status)
          END IF

          PRINT '("Inititalizing flux surface average splines with nu=", i0, " nv=", i0)', nu_local, nv_local
          
          ! Check module initalization.
          IF (.NOT. ALLOCATED(xm)) THEN
             error_status = -4
             RETURN
          ENDIF
          
          !Allocations
          ALLOCATE(xu(nu_local),xv(nv_local),rho(k1:k2))
          ALLOCATE(f_temp(nu_local,nv_local,k1:k2))
          
          f_temp = 0
          rho(:) = 0.0
          FORALL(u=k1+1:k2) rho(u) = REAL(u-1)/REAL(ns_t-1) - 0.5/REAL(ns_t-1)
          FORALL(u=1:nu_local) xu(u) = REAL(u-1)/REAL(nu_local-1)
          FORALL(u=1:nv_local) xv(u) = REAL(u-1)/REAL(nv_local-1)
          
          ! Find the coefficent values at s=0 through linear interpolation.
          rho(k1) = 0.0
          lmnc_half(:,1) = 1.5*lmnc_half(:,2) - 0.5*lmnc_half(:,3)
          lmns_half(:,1) = 1.5*lmns_half(:,2) - 0.5*lmns_half(:,3)

          ! Free Memory
          IF (EZspline_allocated(L_spl)) CALL EZspline_free(L_spl,ez_status)
          IF (EZspline_allocated(Lu_spl)) CALL EZspline_free(Lu_spl,ez_status)
          IF (EZspline_allocated(Lv_spl)) CALL EZspline_free(Lv_spl,ez_status)
             
          ! Preform Init
          CALL EZspline_init(L_spl,nu_local,nv_local,ns_t,bcs1,bcs1,bcs0,ez_status)
          CALL EZspline_init(Lu_spl,nu_local,nv_local,ns_t,bcs1,bcs1,bcs0,ez_status)
          CALL EZspline_init(Lv_spl,nu_local,nv_local,ns_t,bcs1,bcs1,bcs0,ez_status)
          
          ! Lambda
          L_spl%x1 = xu*pi2; L_spl%x2 = xv*pi2; L_spl%x3 = rho; L_spl%isHermite = isherm
          CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,lmns_half,xm,xn,f_temp,1,0)
          IF (lasym) CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,lmnc_half,xm,xn,f_temp,0,0)
          CALL EZspline_setup(L_spl,f_temp,ez_status); f_temp = 0
          ! Lambda/u
          Lu_spl%x1 = xu*pi2; Lu_spl%x2 = xv*pi2; Lu_spl%x3 = rho; Lu_spl%isHermite = isherm
          FORALL(mn = 1:mnmax) fmn_temp(mn,:) = lmns_half(mn,:)*xm(mn)
          CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,f_temp,0,0)
          IF (lasym) THEN
             FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -lmnc_half(mn,:)*xm(mn)
             CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,f_temp,1,0)
          END IF
          CALL EZspline_setup(Lu_spl,f_temp,ez_status); f_temp = 0
          ! Lambda/v
          Lv_spl%x1 = xu*pi2; Lv_spl%x2 = xv*pi2; Lv_spl%x3 = rho; Lv_spl%isHermite = isherm
          FORALL(mn = 1:mnmax) fmn_temp(mn,:) = lmns_half(mn,:)*xn(mn)
          CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,f_temp,0,0)
          IF (lasym) THEN
             FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -lmnc_half(mn,:)*xn(mn)
             CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,f_temp,1,0)
          END IF
          CALL EZspline_setup(Lv_spl,f_temp,ez_status); f_temp = 0
          
          IF (ez_status .NE. 0) error_status = -200 - ABS(ez_status)
          
          ! deallocations
          DEALLOCATE(xu,xv,rho)
          DEALLOCATE(f_temp)
          DEALLOCATE(fmn_temp)

        END SUBROUTINE initialize_splines_lambda

        
        SUBROUTINE initialize_splines_jacobian(error_status, nu, nv)
          !
          ! Programming Notes:
          !   gmnc_half, gmns_half are defined on the half mesh.
          !
          USE EZspline
          IMPLICIT NONE
          INTEGER, INTENT(out)       :: error_status
          INTEGER, INTENT(in), OPTIONAL :: nu
          INTEGER, INTENT(in), OPTIONAL :: nv
          INTEGER :: nu_local, nv_local
          INTEGER :: ez_status
          INTEGER ::  u, mn
          DOUBLE PRECISION, ALLOCATABLE :: xu(:), xv(:), rho(:)
          DOUBLE PRECISION, ALLOCATABLE :: f_temp(:,:,:)
   
          error_status = 0

          IF (PRESENT(nu) .AND. PRESENT(nv)) THEN
             nu_local = nu
             nv_local = nv
          ELSE
             CALL guess_nu_nv(nu_local, nv_local, error_status)
          END IF

          PRINT '("Inititalizing jacobian splines with nu=", i0, " nv=", i0)', nu_local, nv_local
          ! Check module initalization.
          IF (.NOT. ALLOCATED(xm)) THEN
             error_status = -4
             RETURN
          ENDIF
          
          !Allocations
          ALLOCATE(xu(nu_local),xv(nv_local),rho(k1:k2))
          ALLOCATE(f_temp(nu_local,nv_local,k1:k2))
          
          f_temp = 0
          ! Evaluate on the half mesh.
          FORALL(u=k1+1:k2) rho(u) = REAL(u-1)/REAL(ns_t-1) - 0.5/REAL(ns_t-1)
          FORALL(u=1:nu_local) xu(u) = REAL(u-1)/REAL(nu_local-1)
          FORALL(u=1:nv_local) xv(u) = REAL(u-1)/REAL(nv_local-1)

          ! Find the coefficent values at s=0 through linear interpolation.
          rho(k1) = 0.0
          gmnc_half(:,1) = 1.5*gmnc_half(:,2) - 0.5*gmnc_half(:,3)
          gmns_half(:,1) = 1.5*gmns_half(:,2) - 0.5*gmns_half(:,3)
          
          ! Free Memory
          IF (EZspline_allocated(G_spl)) CALL EZspline_free(G_spl,ez_status)
             
          ! Preform Init
          CALL EZspline_init(G_spl,nu_local,nv_local,ns_t,bcs1,bcs1,bcs0,ez_status)
             
          ! Jacobian sqrt(g)
          G_spl%x1 = xu*pi2; G_spl%x2 = xv*pi2; G_spl%x3 = rho; G_spl%isHermite = isherm
          CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,gmnc_half,xm,xn,f_temp,0,0)
          IF (lasym) CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,gmns_half,xm,xn,f_temp,1,0)
          CALL EZspline_setup(G_spl,f_temp,ez_status); f_temp = 0
          
          IF (ez_status .NE. 0) error_status = -200 - ABS(ez_status)
          
          ! deallocations
          DEALLOCATE(xu,xv,rho)
          DEALLOCATE(f_temp)

        END SUBROUTINE initialize_splines_jacobian

        
        SUBROUTINE initialize_splines_fsa(error_status, nu, nv)
          !
          ! Initialize splines for flux surface averaged quantities.
          !
          ! Programming Notes:
          !   For this calculation full mesh values are used for all
          !   quantities.
          ! 
          USE EZspline
          IMPLICIT NONE
          INTEGER, INTENT(out)       :: error_status
          INTEGER, INTENT(in), OPTIONAL :: nu
          INTEGER, INTENT(in), OPTIONAL :: nv
          INTEGER :: nu_local, nv_local
          INTEGER :: ez_status
          INTEGER ::  u, mn
          DOUBLE PRECISION, ALLOCATABLE :: xu(:), xv(:), rho(:)
          DOUBLE PRECISION, ALLOCATABLE :: &
            r_temp(:,:,:) &
            ,ru_temp(:,:,:) & 
            ,rv_temp(:,:,:) & 
            ,zu_temp(:,:,:) & 
            ,zv_temp(:,:,:) &         
            ,g_temp(:,:,:) &
            ,b_temp(:,:,:)  
          DOUBLE PRECISION, ALLOCATABLE :: fmn_temp(:,:)
          DOUBLE PRECISION, ALLOCATABLE :: vp(:), fsa_bmod(:), grho(:), grho2(:), b2overgrho(:)
          DOUBLE PRECISION, ALLOCATABLE :: gsr(:,:,:), gsp(:,:,:), gsz(:,:,:), gs(:,:,:)
          error_status = 0
          
          IF (PRESENT(nu) .AND. PRESENT(nv)) THEN
             nu_local = nu
             nv_local = nv
          ELSE
             CALL guess_nu_nv(nu_local, nv_local, error_status)
          END IF

          PRINT '("Inititalizing flux surface average (fsa) splines with nu=", i0, " nv=", i0)', nu_local, nv_local
          
          ! Check module initalization.
          IF (.NOT. ALLOCATED(xm)) THEN
             error_status = -4
             RETURN
          ENDIF

          !Allocations
          ALLOCATE(xu(nu_local),xv(nv_local),rho(k1:k2))
          ALLOCATE( &
            r_temp(nu_local,nv_local,k1:k2) &
            ,ru_temp(nu_local,nv_local,k1:k2) & 
            ,rv_temp(nu_local,nv_local,k1:k2) & 
            ,zu_temp(nu_local,nv_local,k1:k2) & 
            ,zv_temp(nu_local,nv_local,k1:k2) & 
            ,g_temp(nu_local,nv_local,k1:k2) &
            ,b_temp(nu_local,nv_local,k1:k2) &
            )
          ALLOCATE(fmn_temp(1:mnmax,k1:k2))
       
          FORALL(u=k1:k2) rho(u) = REAL(u-1)/REAL(ns_t-1)
          FORALL(u=1:nu_local) xu(u) = REAL(u-1)/REAL(nu_local-1)
          FORALL(u=1:nv_local) xv(u) = REAL(u-1)/REAL(nv_local-1)

          ! Create temporary arrays using only full mesh quantities.
          !
          ! Because we allow splines to be initialized with different numbers
          ! for nu and nv, we can't use the values of these arrays that were
          ! created when initializing other splines.
          
          ! R values
          r_temp = 0
          CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,rmnc,xm,xn,r_temp,0,1)
          IF (lasym) CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,rmns,xm,xn,r_temp,1,0)
                    
          ! dR/Du Derivatives
          ru_temp = 0
          FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -rmnc(mn,:)*xm(mn)
          CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,ru_temp,1,0)
          IF (lasym) THEN
             FORALL(mn = 1:mnmax) fmn_temp(mn,:) = rmns(mn,:)*xm(mn)
             CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,ru_temp,0,0)
          END IF
          
          ! dR/Dv Derivatives
          rv_temp = 0
          FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -rmnc(mn,:)*xn(mn)
          CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,rv_temp,1,0)
          IF (lasym) THEN
             FORALL(mn = 1:mnmax) fmn_temp(mn,:) = rmns(mn,:)*xn(mn)
             CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,rv_temp,0,0)
          END IF
          
          ! dZ/Du Derivatives
          zu_temp = 0
          FORALL(mn = 1:mnmax) fmn_temp(mn,:) = zmns(mn,:)*xm(mn)
          CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,zu_temp,0,0)
          IF (lasym) THEN
             FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -zmnc(mn,:)*xm(mn)
             CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,zu_temp,1,0)
          END IF
          
          ! dZ/Dv Derivatives
          zv_temp = 0
          FORALL(mn = 1:mnmax) fmn_temp(mn,:) = zmns(mn,:)*xn(mn)
          CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,zv_temp,0,0)
          IF (lasym) THEN
             FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -zmnc(mn,:)*xn(mn)
             CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,zv_temp,1,0)
          END IF
          
          ! Jacobian
          g_temp = 0
          CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,gmnc_full,xm,xn,g_temp,0,0)
          IF (lasym) CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,gmns_full,xm,xn,g_temp,1,0)
          
          ! |B|
          b_temp = 0
          CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,bmnc_full,xm,xn,b_temp,0,0)
          IF (lasym) CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,bmns_full,xm,xn,b_temp,1,0)
                          
          ! Free Memory
          IF (EZspline_allocated(Vp_spl)) CALL EZspline_free(Vp_spl,ez_status)
          IF (EZspline_allocated(b_fsa_spl)) CALL EZspline_free(b_fsa_spl,ez_status)
          IF (EZspline_allocated(grho_fsa_spl)) CALL EZspline_free(grho_fsa_spl,ez_status)
          IF (EZspline_allocated(grho2_fsa_spl)) CALL EZspline_free(grho2_fsa_spl,ez_status)
          IF (EZspline_allocated(b2overgrho_fsa_spl)) CALL EZspline_free(b2overgrho_fsa_spl,ez_status)
          
          ! Preform Init
          CALL EZspline_init(Vp_spl,ns_t,bcs0,ez_status)
          CALL EZspline_init(grho_fsa_spl,ns_t,bcs0,ez_status)
          CALL EZspline_init(b_fsa_spl,ns_t,bcs0,ez_status)
          CALL EZspline_init(grho2_fsa_spl,ns_t,bcs0,ez_status)
          CALL EZspline_init(b2overgrho_fsa_spl,ns_t,bcs0,ez_status)
          
          ALLOCATE( &
            Vp(k1:k2) &
            ,fsa_bmod(k1:k2) &
            ,grho(k1:k2) &
            ,grho2(k1:k2) &
            ,b2overgrho(k1:k2))
          ALLOCATE( &
            gsr(nu_local,nv_local,k1:k2) &
            ,gsp(nu_local,nv_local,k1:k2) &
            ,gsz(nu_local,nv_local,k1:k2),&
            gs(nu_local,nv_local,k1:k2))
          
          ! Calc grad(s) components dR/du X dR/dv / sqrt(g)
          !    Note component of R_spl comes from dphi/dphi and cyl coordinates
          gsr = - zu_temp*r_temp
          gsp = (zu_temp*Rv_temp - Ru_temp*Zv_temp)*nfp
          gsz = Ru_temp*R_temp
          
          gs  = (gsr*gsr+gsp*gsp+gsz*gsz)/(g_temp*g_temp)  !|grad(s)|^2
          FORALL(u=k1:k2) gs(:,:,u) = gs(:,:,u)/(4*rho(u)) !|grad(rho)|^2

          ! dV/ds
          Vp = SUM(SUM(g_temp,DIM=1),DIM=1)
          !Vp(1) = 2*Vp(2) - Vp(3)

          ! <|B|>
          fsa_bmod = SUM(SUM(b_temp*g_temp,DIM=1),DIM=1)
          fsa_bmod = fsa_bmod / Vp
          fsa_bmod(1) = 2*fsa_bmod(2) - fsa_bmod(3)
          
          ! <|grad(rho)|^2>
          grho2 = SUM(SUM(gs*g_temp,DIM=1),DIM=1)
          grho2 = grho2 / Vp
          grho2(1) = 2*grho2(2) - grho2(3)
          
          ! <|grad(rho)|>|
          gs = sqrt(gs) !|grad(rho)|
          grho = SUM(SUM(gs*g_temp,DIM=1),DIM=1)
          grho = grho / Vp
          grho(1) = 2*grho(2) - grho(3)

          ! <|B|^2/grad(rho)>
          b2overgrho = SUM(SUM(b_temp**2*g_temp/gs,DIM=1),DIM=1)
          b2overgrho = b2overgrho / Vp
          b2overgrho(1) = 2*b2overgrho(2) - b2overgrho(3)        
          
          ! Construct splines
          CALL EZspline_setup(Vp_spl,ABS(Vp*pi2*pi2/(nu_local*nv_local)),ez_status) ! ABS because of negative Jacobian
          CALL EZspline_setup(b_fsa_spl,fsa_bmod,ez_status)
          CALL EZspline_setup(grho_fsa_spl,grho,ez_status)
          CALL EZspline_setup(grho2_fsa_spl,grho2,ez_status)
          CALL EZspline_setup(b2overgrho_fsa_spl,b2overgrho,ez_status)
          
          
          IF (ez_status .NE. 0) error_status = -200 - ABS(ez_status)
          
          ! deallocations
          DEALLOCATE(gsr,gsp,gsz,gs,Vp,fsa_bmod,grho,grho2,b2overgrho)
          DEALLOCATE(xu,xv,rho)
          DEALLOCATE( &
            r_temp &
            ,ru_temp & 
            ,rv_temp & 
            ,zu_temp & 
            ,zv_temp & 
            ,g_temp &
            ,b_temp &
            )
          DEALLOCATE(fmn_temp)
        
        END SUBROUTINE initialize_splines_fsa

        
        SUBROUTINE initialize_splines_susceptance(error_status, nu, nv)
          !
          ! Initialize splines for susceptance matrix quantities.
          !
          ! Programming Notes:
          !   For this calculation full mesh values are used for all
          !   quantities.
          ! 
          USE EZspline
          IMPLICIT NONE
          INTEGER, INTENT(out)       :: error_status
          INTEGER, INTENT(in), OPTIONAL :: nu
          INTEGER, INTENT(in), OPTIONAL :: nv
          INTEGER :: nu_local, nv_local
          INTEGER :: ez_status
          INTEGER ::  u, mn
          DOUBLE PRECISION, ALLOCATABLE :: xu(:), xv(:), rho(:)
          DOUBLE PRECISION, ALLOCATABLE :: f_temp(:,:,:)
          DOUBLE PRECISION, ALLOCATABLE :: &
            r_temp(:,:,:) &
            ,ru_temp(:,:,:) & 
            ,rv_temp(:,:,:) & 
            ,zu_temp(:,:,:) & 
            ,zv_temp(:,:,:) &         
            ,g_temp(:,:,:) &
            ,lu_temp(:,:,:) & 
            ,lv_temp(:,:,:)
          DOUBLE PRECISION, ALLOCATABLE :: fmn_temp(:,:)
          DOUBLE PRECISION, ALLOCATABLE :: vp(:), grho(:), grho2(:)  
          error_status = 0

          IF (PRESENT(nu) .AND. PRESENT(nv)) THEN
             nu_local = nu
             nv_local = nv
          ELSE
             CALL guess_nu_nv(nu_local, nv_local, error_status)
          END IF

          PRINT '("Inititalizing flux surface average splines with nu=", i0, " nv=", i0)', nu_local, nv_local
          
          ! Check module initalization.
          IF (.NOT. ALLOCATED(xm)) THEN
             error_status = -4
             RETURN
          ENDIF
          
          !Allocations
          ALLOCATE(xu(nu_local),xv(nv_local),rho(k1:k2))
          ALLOCATE(f_temp(nu_local,nv_local,k1:k2))
          ALLOCATE( &
            r_temp(nu_local,nv_local,k1:k2) &
            ,ru_temp(nu_local,nv_local,k1:k2) & 
            ,rv_temp(nu_local,nv_local,k1:k2) & 
            ,zu_temp(nu_local,nv_local,k1:k2) & 
            ,zv_temp(nu_local,nv_local,k1:k2) & 
            ,g_temp(nu_local,nv_local,k1:k2) &
            ,lu_temp(nu_local,nv_local,k1:k2) & 
            ,lv_temp(nu_local,nv_local,k1:k2) & 
            )
          ALLOCATE(fmn_temp(1:mnmax,k1:k2))
          
          f_temp = 0
          FORALL(u=k1:k2) rho(u) = REAL(u-1)/REAL(ns_t-1)
          FORALL(u=1:nu_local) xu(u) = REAL(u-1)/REAL(nu_local-1)
          FORALL(u=1:nv_local) xv(u) = REAL(u-1)/REAL(nv_local-1)


          ! Create temporary arrays using only full mesh quantities.
          !
          ! Because we allow splines to be initialized with different numbers of
          ! nu and nv, we can't use the values of these arrays that were created
          ! when initializing other splines.
          
          ! R values
          r_temp = 0
          CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,rmnc,xm,xn,r_temp,0,1)
          IF (lasym) CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,rmns,xm,xn,r_temp,1,0)
                    
          ! dR/Du Derivatives
          ru_temp = 0
          FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -rmnc(mn,:)*xm(mn)
          CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,ru_temp,1,0)
          IF (lasym) THEN
             FORALL(mn = 1:mnmax) fmn_temp(mn,:) = rmns(mn,:)*xm(mn)
             CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,ru_temp,0,0)
          END IF
          
          ! dR/Dv Derivatives
          rv_temp = 0
          FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -rmnc(mn,:)*xn(mn)
          CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,rv_temp,1,0)
          IF (lasym) THEN
             FORALL(mn = 1:mnmax) fmn_temp(mn,:) = rmns(mn,:)*xn(mn)
             CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,rv_temp,0,0)
          END IF
          
          ! dZ/Du Derivatives
          zu_temp = 0
          FORALL(mn = 1:mnmax) fmn_temp(mn,:) = zmns(mn,:)*xm(mn)
          CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,zu_temp,0,0)
          IF (lasym) THEN
             FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -zmnc(mn,:)*xm(mn)
             CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,zu_temp,1,0)
          END IF
          
          ! dZ/Dv Derivatives
          zv_temp = 0
          FORALL(mn = 1:mnmax) fmn_temp(mn,:) = zmns(mn,:)*xn(mn)
          CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,zv_temp,0,0)
          IF (lasym) THEN
             FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -zmnc(mn,:)*xn(mn)
             CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,zv_temp,1,0)
          END IF
          
          ! Jacobian
          g_temp = 0
          CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,gmnc_full,xm,xn,g_temp,0,0)
          IF (lasym) CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,gmns_full,xm,xn,g_temp,1,0)

          ! Lambda/u
          lu_temp = 0
          FORALL(mn = 1:mnmax) fmn_temp(mn,:) = lmns_full(mn,:)*xm(mn)
          CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,lu_temp,0,0)
          IF (lasym) THEN
             FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -lmnc_full(mn,:)*xm(mn)
             CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,lu_temp,1,0)
          END IF
          
          ! Lambda/v
          lv_temp = 0
          FORALL(mn = 1:mnmax) fmn_temp(mn,:) = lmns_full(mn,:)*xn(mn)
          CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,lv_temp,0,0)
          IF (lasym) THEN
             FORALL(mn = 1:mnmax) fmn_temp(mn,:) = -lmnc_full(mn,:)*xn(mn)
             CALL mntouv(k1,k2,mnmax,nu_local,nv_local,xu,xv,fmn_temp,xm,xn,lv_temp,1,0)
          END IF

          
          ! Free Memory
          IF (EZspline_allocated(S11_spl)) CALL EZspline_free(S11_spl,ez_status)
          IF (EZspline_allocated(S12_spl)) CALL EZspline_free(S12_spl,ez_status)
          IF (EZspline_allocated(S21_spl)) CALL EZspline_free(S21_spl,ez_status)
          IF (EZspline_allocated(S22_spl)) CALL EZspline_free(S22_spl,ez_status)
             
          ! Preform Init
          CALL EZspline_init(S11_spl,ns_t,bcs0,ez_status)
          CALL EZspline_init(S12_spl,ns_t,bcs0,ez_status)
          CALL EZspline_init(S21_spl,ns_t,bcs0,ez_status)
          CALL EZspline_init(S22_spl,ns_t,bcs0,ez_status)
             
          ! Calc S11
          f_temp = (Ru_temp*Ru_temp+ &
               Zu_temp*Zu_temp)
          f_temp = f_temp / g_temp
          grho   = SUM(SUM(f_temp,DIM=1),DIM=1)/(nu_local*nv_local)
          grho(1) = 2*grho(2) - grho(3)
          CALL EZspline_setup(S11_spl,grho,ez_status); f_temp = 0; grho = 0
          ! Calc S21
          f_temp = (Ru_temp*Rv_temp+ &
               Zu_temp*Zv_temp)*nfp
          f_temp = f_temp / g_temp
          grho   = SUM(SUM(f_temp,DIM=1),DIM=1)/(nu_local*nv_local)
          grho(1) = 2*grho(2) - grho(3)
          CALL EZspline_setup(S21_spl,grho,ez_status); f_temp = 0; grho = 0
          ! Calc S12
          f_temp = (Ru_temp*Rv_temp+ &
               Zu_temp*Zv_temp)* &
               (one+Lu_temp)*nfp
          f_temp = f_temp - (Ru_temp*Ru_temp+ &
               Zu_temp*Zu_temp)*&
               Lv_temp
          f_temp = f_temp / g_temp
          grho   = SUM(SUM(f_temp,DIM=1),DIM=1)/(nu_local*nv_local)
          grho(1) = 2*grho(2) - grho(3)
          CALL EZspline_setup(S12_spl,grho,ez_status); f_temp = 0; grho = 0
          ! Calc S22
          f_temp = (Rv_temp*Rv_temp*nfp*nfp+ &
               Zv_temp*Zv_temp*nfp*nfp+ &
               R_temp* R_temp)* &
               (one+Lu_temp)
          f_temp = f_temp - (Ru_temp*Rv_temp+ &
               Zu_temp*Zv_temp)*&
               Lv_temp*nfp
          f_temp = f_temp / g_temp
          grho   = SUM(SUM(f_temp,DIM=1),DIM=1)/(nu_local*nv_local)
          !grho(1) = 2*grho(2) - grho(3)
          CALL EZspline_setup(S22_spl,grho,ez_status); f_temp = 0; grho = 0
          
          IF (ez_status .NE. 0) error_status = -200 - ABS(ez_status)
          
          ! deallocations
          DEALLOCATE(xu,xv,rho)
          DEALLOCATE(f_temp)
          DEALLOCATE( &
            r_temp &
            ,ru_temp & 
            ,rv_temp & 
            ,zu_temp & 
            ,zv_temp & 
            ,g_temp &
            ,lu_temp & 
            ,lv_temp & 
            )
          DEALLOCATE(fmn_temp)

        END SUBROUTINE initialize_splines_susceptance

        
        SUBROUTINE rzfunct_stel_tool(m,n,x,fvec,fjac,ldfjac,iflag)
          USE EZspline
          IMPLICIT NONE
          INTEGER :: m,n,ldfjac,iflag, ier
          DOUBLE PRECISION :: x(n),fvec(m),fjac(ldfjac,n)
          DOUBLE PRECISION :: R_temp, Z_temp
          DOUBLE PRECISION :: R_grad(3), Z_grad(3)
          IF (x(2) < 0.0) x(2) = x(2) + pi2
          x(2) = MOD(x(2),pi2)
          IF (x(1) < 0) THEN
             x(1) = ABS(x(1))
             x(2) = x(2)+pi2*0.5
             x(2) = MOD(x(2),pi2)
             !x(1) = 0
          END IF
          ier = 0
          CALL EZspline_isInDomain(R_spl,x(2),PHI_Target,x(1),ier)
          IF (ier .ne. 0) THEN
             iflag = -1
             domain_flag = -1
          END IF
          IF (iflag == 1) THEN
             CALL EZspline_interp(R_spl,x(2),PHI_Target,x(1),R_temp,iflag)
             CALL EZspline_interp(Z_spl,x(2),PHI_Target,x(1),Z_temp,iflag)
             fvec(1) = (R_temp - R_target)
             fvec(2) = (Z_temp - Z_target)
          ELSE IF (iflag == 2) THEN
             CALL EZspline_gradient(R_spl,x(2),PHI_Target,x(1),R_grad,iflag)
             CALL EZspline_gradient(Z_spl,x(2),PHI_Target,x(1),Z_grad,iflag)
             !CALL EZspline_interp(Ru_spl,x(2),PHI_Target,x(1),R_temp,iflag)
             !CALL EZspline_interp(Zu_spl,x(2),PHI_Target,x(1),Z_temp,iflag)
             fjac(1,1) = R_grad(3) !dR/ds
             fjac(1,2) = R_grad(1) !dR/du
             !fjac(1,2) = R_temp
             fjac(2,1) = Z_grad(3) !dZ/ds
             fjac(2,2) = Z_grad(1) !dZ/du
             !fjac(2,2) = Z_temp
          END IF
          RETURN
        END SUBROUTINE rzfunct_stel_tool

        
        SUBROUTINE flx_from_cyl(point_cyl, point_flx, error_status, guess)
          ! Purpose:
          !   Given a set of real space (cylindrical) coordinates return
          !   a set of flux coordinates.
          !
          ! Error codes:
          !   error_status = -2  Value is outside of flux space domain.
          !  
          USE EZspline
          IMPLICIT NONE
          ! point_cyl = (R, phi, z)
          !   Realspace cylindrical cooridantes.
          DOUBLE PRECISION, INTENT(in) ::  point_cyl(3)
          ! point_flx = (s, u, phi)
          !   VMEC flux space coordinates, but with realspace phi.
          DOUBLE PRECISION, INTENT(out) ::  point_flx(3)
          INTEGER, INTENT(out) ::  error_status
          DOUBLE PRECISION, INTENT(in), OPTIONAL ::  guess

          INTEGER, PARAMETER :: mfunct=2
          INTEGER, PARAMETER :: nvars=2
          INTEGER, PARAMETER :: ldfjac=2
          INTEGER :: ik, maxfev_local, nfev, info, njev, maxfev, nprint, mode
          INTEGER, DIMENSION(nvars) :: ipvt
          DOUBLE PRECISION :: ftol,xtol,gtol,factor
          DOUBLE PRECISION, DIMENSION(nvars)  :: xc_opt, diag, qtf, wa1, wa2, wa3
          DOUBLE PRECISION, DIMENSION(mfunct) :: fval,wa4
          DOUBLE PRECISION, DIMENSION(ldfjac,nvars) :: fjac
          DOUBLE PRECISION :: guess_flx(3)
          LOGICAL :: l_make_guess
      
          DOUBLE PRECISION :: enorm

          ! In the future we can have this flag controled by an optional
          ! parameter.
          l_make_guess = .true.

          error_status = 0
          point_flx(:) = 0

          IF (.NOT. EZspline_allocated(R_spl)) CALL initialize_splines_rz(error_status)
          IF (error_status < 0) RETURN
          
          IF (PRESENT(guess)) THEN
             PRINT *, 'flx_from_cyl: Using user guess.'
             guess_flx = guess
             
             IF (guess_flx(1) > 1 .or. guess_flx(1) < 0) guess_flx = 0
             guess_flx(2) = MOD(guess_flx(2), pi2)
             if (guess_flx(2) < 0) guess_flx(2) = guess_flx(2) + pi2
             guess_flx(3) = MOD(guess_flx(3), pi2)
             if (guess_flx(3) < 0) guess_flx(3) = guess_flx(3) + pi2
          ELSE IF (l_make_guess) THEN
             CALL guess_flx_from_cyl(point_cyl, guess_flx, error_status)
          ELSE
             guess_flx(1) = 0.5
             guess_flx(2) = pi2*0.6
             guess_flx(3) = point_cyl(2)
          END IF

          error_status = 0
          IF (EZspline_allocated(R_spl) .and. EZspline_allocated(Z_spl)) THEN
             ! These target variables are modules level variables that are used
             ! in rzfunct_stel_tool as part of the optimization loop.
             R_target = point_cyl(1)
             Z_target = point_cyl(3)
             PHI_target = point_cyl(2)
             IF (PHI_target < 0) PHI_target = PHI_target + pi2
             PHI_target = MOD(PHI_target,pi2/nfp)*nfp
             xc_opt(1) = guess_flx(1)
             xc_opt(2) = guess_flx(2)
             fval = 1.0E-30
             fjac = 0
             ftol = search_tol
             xtol = search_tol
             gtol = 1.0E-30
             maxfev_local = neval_max
             diag = (/1.0,4.0/)
             mode = 2
             factor = 0.1
             nprint = 0
             info   = -2
             nfev   = 0
             njev   = 0
             ipvt   = 0
             qtf    = 0.0
             wa1 = 0; wa2 = 0; wa3 = 0; wa4 = 0
             domain_flag = 0
             DO ik = 1, 4
                info = 0
                nfev = 0
                njev = 0
                CALL lmder_serial( &
                  rzfunct_stel_tool &
                  ,mfunct &
                  ,nvars &
                  ,xc_opt &
                  ,fval &
                  ,fjac &
                  ,ldfjac &
                  ,ftol &
                  ,xtol &
                  ,gtol &
                  ,maxfev_local &
                  ,diag &
                  ,mode &
                  ,factor &
                  ,nprint &
                  ,info &
                  ,nfev &
                  ,njev &
                  ,ipvt &
                  ,qtf &
                  ,wa1 &
                  ,wa2 &
                  ,wa3 &
                  ,wa4)
                ! Check for absolute chi-sq convergence. This test is not done within the L-M routine.
                IF (info > 0) THEN
                   IF (enorm(mfunct, fval)**2 < ftol) THEN
                      info = 1
                   ENDIF
                ENDIF
                IF (info < 4) EXIT
                ! If we have not achived success, try to bump the fitter into a better space.
                IF (ik < 4) THEN 
                   xc_opt(1) = xc_opt(1) + 1e-4 * 10**(ik-1)
                   xc_opt(2) = xc_opt(2) + 1e-4*pi2 * 10**(ik-1)
                END IF
             END DO
             error_status = info
             ! Any values of info greater than 0 may indicate success, however
             ! values greater than 3 may indicate a problem in convergence.
             IF ((info > 0) .AND. (info < 4)) error_status = 0
             IF (domain_flag .ne. 0) THEN
                error_status = -2
                point_flx(1) = 1.5
                point_flx(2) = 2*pi2
                point_flx(3) = 2*pi2
                RETURN
             END IF
             point_flx(1) = xc_opt(1)
             point_flx(2) = xc_opt(2)
             point_flx(2) = MOD(point_flx(2), pi2)
             IF (point_flx(2) < 0) point_flx(2) = point_flx(2) + pi2
             point_flx(3) = point_cyl(2)
          ELSE
             error_status = -1
          END IF
          RETURN
        END SUBROUTINE flx_from_cyl

        
        SUBROUTINE guess_flx_from_cyl_dbl(point_cyl, guess_flx, error_status)
          ! Purpose:
          !   Make a guess of flux space coordinates given a set of real space
          !   (cylindrical) coordinates.  This is only valid for toroidal geometries.
          !
          !   For certain geometrys this guess will be reasonably good, however for
          !   many stellarator geometrys this will only be a crude guess.
          IMPLICIT NONE
          ! point_cyl = (R, phi, z)
          !   Realspace cylindrical cooridantes.
          DOUBLE PRECISION, INTENT(in) ::  point_cyl(3)
          ! guess_flx = (s, u, phi)
          !   VMEC flux space coordinates, but with realspace phi.
          DOUBLE PRECISION, INTENT(out) ::  guess_flx(3)
          INTEGER, INTENT(out) ::  error_status
      
          DOUBLE PRECISION :: axis(3)
          DOUBLE PRECISION :: edge(3)
          DOUBLE PRECISION :: s_temp, u_temp, v_temp
          DOUBLE PRECISION :: temp_flx(3)
          
          error_status = 0

          guess_flx(3) = point_cyl(2)

          ! Find the axis at the given phi_value.
          temp_flx(:) = 0
          temp_flx(3) = point_cyl(2)
          
          CALL cyl_from_flx(temp_flx, axis, error_status)
          IF (error_status < 0) RETURN
      
          ! Calculate the u_val guess.
          guess_flx(2) = atan2(point_cyl(3)-axis(3), point_cyl(1)-axis(1))
          IF (guess_flx(2) < 0D0) THEN
             guess_flx(2) = guess_flx(2)+pi2
          END IF
      
          ! Find the edge at the given u and v values.
          temp_flx(1) = 1.0
          CALL cyl_from_flx(temp_flx, edge, error_status)
          IF (error_status < 0) RETURN
      
          ! Make a guess for s. This assumes that sqrt(s) is approx proportional
          ! to real space distance. 
          guess_flx(1) = ((point_cyl(1)-axis(1))**2 + (point_cyl(3)-axis(3))**2) &
                         /((edge(1)-axis(1))**2+(edge(3)-axis(3))**2)
      
          IF (guess_flx(1) > 1) THEN
             guess_flx(1) = 1.0
          ENDIF

        END SUBROUTINE guess_flx_from_cyl_dbl

        SUBROUTINE guess_flx_from_cyl_sgl(point_cyl, guess_flx, error_status)
          IMPLICIT NONE
          REAL, INTENT(in) ::  point_cyl(3)
          REAL, INTENT(out) ::  guess_flx(3)
          INTEGER, INTENT(out) ::  error_status
          
          DOUBLE PRECISION ::  point_cyl_dbl(3)
          DOUBLE PRECISION ::  guess_flx_dbl(3)
          
          point_cyl_dbl(:) = point_cyl(:)
          CALL guess_flx_from_cyl_dbl(point_cyl_dbl, guess_flx_dbl, error_status)
          guess_flx(:) = guess_flx_dbl(:)
          RETURN
        END SUBROUTINE guess_flx_from_cyl_sgl


        SUBROUTINE cyl_from_flx(point_flx, point_cyl, error_status)
          ! Purpose:
          !   Given a set of real space (cylindrical) coordinates return
          !   a set of flux coordinates.
          !
          ! Error codes:
          !   error_status = -2  Value is outside of flux space domain.
          !  
          USE EZspline
          IMPLICIT NONE
          ! point_flx = (s, u, phi)
          !   VMEC flux space coordinates, but with realspace phi.
          DOUBLE PRECISION, INTENT(in) ::  point_flx(3)
          ! point_cyl = (R, phi, z)
          !   Realspace cylindrical cooridantes.
          DOUBLE PRECISION, INTENT(out) ::  point_cyl(3)
          INTEGER, INTENT(out) ::  error_status
          DOUBLE PRECISION :: v_val

          error_status = 0
          point_cyl(:) = 0
          
          IF (.NOT. EZspline_allocated(R_spl)) CALL initialize_splines_rz(error_status)
          IF (error_status < 0) RETURN
          
          v_val = MOD(point_flx(3) , pi2/nfp)*nfp
          
          CALL EZspline_isInDomain(R_spl, point_flx(2), v_val, point_flx(1), error_status)
          IF (error_status == 0) THEN
             CALL EZspline_interp(R_spl, point_flx(2), v_val, point_flx(1), point_cyl(1), error_status)
             CALL EZspline_interp(Z_spl, point_flx(2), v_val, point_flx(1), point_cyl(3), error_status)
             point_cyl(2) = point_flx(3)
          ELSE
             error_status = -2
          END IF
          RETURN
        END SUBROUTINE cyl_from_flx

        
        SUBROUTINE modb_from_flx(point_flx, modb, error_status)
          USE EZspline
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(in) :: point_flx(3)
          DOUBLE PRECISION, INTENT(out) :: modb
          INTEGER, INTENT(out) :: error_status
         
          DOUBLE PRECISION :: v_val

          error_status = 0
          modb = 0
          
          v_val = MOD(point_flx(3), pi2/nfp)*nfp

          IF (.NOT. EZspline_allocated(B_spl)) CALL initialize_splines_modb(error_status)
          IF (error_status < 0) RETURN
          
          CALL EZSPLINE_isInDomain(B_spl, point_flx(2), v_val, point_flx(1), error_status)
          IF (error_status == 0) THEN
             CALL EZspline_interp(B_spl, point_flx(2), v_val, point_flx(1), modb, error_status)
          ELSE
             error_status = -2
          END IF

        END SUBROUTINE modb_from_flx
        
        
        SUBROUTINE b_flx_from_flx(point_flx, b_flx, error_status)
          USE EZspline
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(in) :: point_flx(3)
          DOUBLE PRECISION, INTENT(out) :: b_flx(3)
          INTEGER, INTENT(out) :: error_status
         
          DOUBLE PRECISION :: v_val

          error_status = 0
          b_flx(:) = 0
          
          v_val = MOD(point_flx(3), pi2/nfp)*nfp

          IF (.NOT. EZspline_allocated(Bs_spl)) CALL initialize_splines_b(error_status)
          IF (error_status < 0) RETURN
          
          CALL EZSPLINE_isInDomain(Bs_spl, point_flx(2), v_val, point_flx(1), error_status)
          IF (error_status == 0) THEN
             CALL EZspline_interp(Bs_spl, point_flx(2), v_val, point_flx(1), b_flx(1), error_status)
             CALL EZspline_interp(Bu_spl, point_flx(2), v_val, point_flx(1), b_flx(2), error_status)
             CALL EZspline_interp(Bv_spl, point_flx(2), v_val, point_flx(1), b_flx(3), error_status)
          ELSE
             error_status = -2
          END IF

        END SUBROUTINE b_flx_from_flx

        
        SUBROUTINE b_cyl_from_flx(point_flx, b_cyl, error_status)
          USE EZspline
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(in) :: point_flx(3)
          DOUBLE PRECISION, INTENT(out) :: b_cyl(3)
          INTEGER, INTENT(out) :: error_status

          DOUBLE PRECISION :: point_cyl(3)
          DOUBLE PRECISION :: b_flx(3)
          DOUBLE PRECISION :: R_grad(3), Z_grad(3)
         
          DOUBLE PRECISION :: v_val

          error_status = 0
          b_cyl(:) = 0
          
          IF (.NOT. EZspline_allocated(Ru_spl)) CALL initialize_splines_rzderiv(error_status)
          IF (.NOT. EZspline_allocated(Bs_spl)) CALL initialize_splines_b(error_status)
          IF (error_status < 0) RETURN
          
          v_val = MOD(point_flx(3), pi2/nfp)*nfp

          CALL b_flx_from_flx(point_flx, b_flx, error_status)
          IF (error_status < 0) RETURN

          CALL cyl_from_flx(point_flx, point_cyl, error_status)
          IF (error_status < 0) RETURN

          CALL EZSPLINE_isInDomain(R_spl, point_flx(2), v_val, point_flx(1), error_status)
          IF (error_status == 0) THEN
             R_grad = 0; Z_grad = 0
             CALL EZspline_interp(Ru_spl, point_flx(2), v_val, point_flx(1), R_grad(1), error_status)
             CALL EZspline_interp(Rv_spl, point_flx(2), v_val, point_flx(1), R_grad(2), error_status)
             CALL EZspline_interp(Zu_spl, point_flx(2), v_val, point_flx(1), Z_grad(1), error_status)
             CALL EZspline_interp(Zv_spl, point_flx(2), v_val, point_flx(1), Z_grad(2), error_status)
             b_cyl(1) = R_grad(3)*b_flx(1) + R_grad(1)*b_flx(2) + R_grad(2)*b_flx(3)*nfp
             b_cyl(2) = point_cyl(1) * b_flx(3)
             b_cyl(3) = Z_grad(3)*b_flx(1) + Z_grad(1)*b_flx(2) + Z_grad(2)*b_flx(3)*nfp
          ELSE
             error_status = -2
          END IF

        END SUBROUTINE b_cyl_from_flx

        
        SUBROUTINE b_car_from_flx(point_flx, b_car, error_status)
          USE EZspline
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(in) :: point_flx(3)
          DOUBLE PRECISION, INTENT(out) :: b_car(3)
          INTEGER, INTENT(out) :: error_status

          DOUBLE PRECISION :: b_cyl(3)

          error_status = 0
          b_car(:) = 0

          CALL b_cyl_from_flx(point_flx, b_cyl, error_status)
          IF (error_status < 0) RETURN
          
          b_car(1) = b_cyl(1) * cos(point_flx(3)) - b_cyl(2) * sin(point_flx(3))
          b_car(2) = b_cyl(1) * sin(point_flx(3)) + b_cyl(2) * cos(point_flx(3))
          b_car(3) = b_cyl(3)

        END SUBROUTINE b_car_from_flx

        
        SUBROUTINE b_cyl_from_cyl(point_cyl, b_cyl, error_status)
          USE EZspline
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(in) :: point_cyl(3)
          DOUBLE PRECISION, INTENT(out) :: b_cyl(3)
          INTEGER, INTENT(out) :: error_status
          DOUBLE PRECISION :: point_flx(3)
          error_status = 0
          b_cyl(:) = 0
          CALL flx_from_cyl(point_cyl, point_flx, error_status)
          IF (error_status < 0) RETURN
          CALL b_cyl_from_flx(point_flx, b_cyl, error_status)
        END SUBROUTINE b_cyl_from_cyl

        
        SUBROUTINE b_car_from_cyl(point_cyl, b_car, error_status)
          USE EZspline
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(in) :: point_cyl(3)
          DOUBLE PRECISION, INTENT(out) :: b_car(3)
          INTEGER, INTENT(out) :: error_status
          DOUBLE PRECISION :: point_flx(3)
          error_status = 0
          b_car(:) = 0
          CALL flx_from_cyl(point_cyl, point_flx, error_status)
          IF (error_status < 0) RETURN
          CALL b_car_from_flx(point_flx, b_car, error_status)
        END SUBROUTINE b_car_from_cyl

        
        SUBROUTINE grads_cyl_from_flx(point_flx, grads_cyl, error_status)
          USE EZspline
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(in) :: point_flx(3)
          DOUBLE PRECISION, INTENT(out) :: grads_cyl(3)
          INTEGER, INTENT(out) :: error_status

          DOUBLE PRECISION :: r_cyl, g_cyl
          DOUBLE PRECISION :: R_grad(3), Z_grad(3)
         
          DOUBLE PRECISION :: v_val

          error_status = 0
          grads_cyl(:) = 0

          IF (.NOT. EZspline_allocated(R_spl)) CALL initialize_splines_rz(error_status)
          IF (.NOT. EZspline_allocated(Ru_spl)) CALL initialize_splines_rzderiv(error_status)
          IF (.NOT. EZspline_allocated(G_spl)) CALL initialize_splines_jacobian(error_status)
          IF (error_status < 0) RETURN
                    
          v_val = MOD(point_flx(3), pi2/nfp)*nfp

          CALL EZSPLINE_isInDomain(R_spl, point_flx(2), v_val, point_flx(1), error_status)
          IF (error_status == 0) THEN
             R_grad = 0; Z_grad = 0
             CALL EZspline_interp(R_spl, point_flx(2), v_val, point_flx(1), r_cyl, error_status)
             CALL EZspline_interp(Ru_spl, point_flx(2), v_val, point_flx(1), R_grad(1), error_status)
             CALL EZspline_interp(Rv_spl, point_flx(2), v_val, point_flx(1), R_grad(2), error_status)
             CALL EZspline_interp(Zu_spl, point_flx(2), v_val, point_flx(1), Z_grad(1), error_status)
             CALL EZspline_interp(Zv_spl, point_flx(2), v_val, point_flx(1), Z_grad(2), error_status)
             CALL EZspline_interp(G_spl, point_flx(2), v_val, point_flx(1), g_cyl, error_status)
             
             ! Calc grad(s) components dR/du X dR/dv / sqrt(g)
             !    Note component of R_spl comes from dphi/dphi and cyl coordinates
             grads_cyl(1) = - Z_grad(1)*r_cyl
             grads_cyl(2) = (Z_grad(1)*R_grad(2) - R_grad(1)*Z_grad(1))*nfp
             grads_cyl(3) = R_grad(1)*r_cyl
             grads_cyl(:) = grads_cyl(:)/g_cyl
          ELSE
             error_status = -2
          END IF
          
        END SUBROUTINE grads_cyl_from_flx

        
        SUBROUTINE grads_car_from_flx(point_flx, grads_car, error_status)
          USE EZspline
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(in) :: point_flx(3)
          DOUBLE PRECISION, INTENT(out) :: grads_car(3)
          INTEGER, INTENT(out) :: error_status

          DOUBLE PRECISION :: grads_cyl(3)
          DOUBLE PRECISION :: point_cyl(3)

          error_status = 0
          grads_car(:) = 0

          IF (.NOT. EZspline_allocated(R_spl)) CALL initialize_splines_rz(error_status)

          CALL grads_cyl_from_flx(point_flx, grads_cyl, error_status)
          IF (error_status < 0 ) RETURN

          CALL cyl_from_flx(point_flx, point_cyl, error_status)

          grads_car(1) = COS(point_cyl(2))*grads_cyl(1) - SIN(point_cyl(2))*grads_cyl(2)
          grads_car(2) = SIN(point_cyl(2))*grads_cyl(1) + COS(point_cyl(2))*grads_cyl(2)
          grads_car(3) = grads_cyl(3)
          
        END SUBROUTINE grads_car_from_flx

        
        SUBROUTINE gradrho_car_from_flx(point_flx, gradrho_car, error_status)
          USE EZspline
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(in) :: point_flx(3)
          DOUBLE PRECISION, INTENT(out) :: gradrho_car(3)
          INTEGER, INTENT(out) :: error_status

          DOUBLE PRECISION :: grads_car(3)

          error_status = 0
          gradrho_car(:) = 0

          IF (point_flx(1) .EQ. 0) THEN
             error_status = -5
             RETURN
          ENDIF
          
          CALL grads_car_from_flx(point_flx, grads_car, error_status)
          IF (error_status < 0 ) RETURN

          gradrho_car = 0.5*grads_car/SQRT(point_flx(1))
          
        END SUBROUTINE gradrho_car_from_flx
        
        
        SUBROUTINE jacobian_from_flx(point_flx, jacobian, error_status)
          USE EZspline
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(in) :: point_flx(3)
          DOUBLE PRECISION, INTENT(out) :: jacobian
          INTEGER, INTENT(out) :: error_status
         
          DOUBLE PRECISION :: v_val

          error_status = 0
          jacobian = 0

          IF (.NOT. EZspline_allocated(G_spl)) CALL initialize_splines_jacobian(error_status)
          IF (.NOT. EZspline_allocated(Ru_spl)) CALL initialize_splines_rzderiv(error_status)
          IF (error_status < 0) RETURN
                    
          v_val = MOD(point_flx(3), pi2/nfp)*nfp

          CALL EZSPLINE_isInDomain(G_spl, point_flx(2), v_val, point_flx(1), error_status)
          IF (error_status == 0) THEN
             CALL EZspline_interp(G_spl, point_flx(2), v_val, point_flx(1), jacobian, error_status)
          ELSE
             error_status = -2
          END IF
          
        END SUBROUTINE jacobian_from_flx

        
        SUBROUTINE fsa_modb_from_s(s_val, fsa_modb, error_status)
          USE EZspline
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(in) :: s_val
          DOUBLE PRECISION, INTENT(out) :: fsa_modb
          INTEGER, INTENT(out) :: error_status

          error_status = 0
          fsa_modb = 0
          
          IF (.NOT. EZspline_allocated(b_fsa_spl)) CALL initialize_splines_fsa(error_status)
          IF (error_status < 0) RETURN
          
          IF (s_val >= 0 .and. s_val <= 1) THEN
             CALL EZspline_interp(b_fsa_spl, s_val, fsa_modb, error_status)
          ELSE
             error_status = -2
          END IF

        END SUBROUTINE fsa_modb_from_s

        
        SUBROUTINE fsa_gradrho_from_s(s_val, fsa_gradrho, error_status)
          USE EZspline
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(in) :: s_val
          DOUBLE PRECISION, INTENT(out) :: fsa_gradrho
          INTEGER, INTENT(out) :: error_status

          error_status = 0
          fsa_gradrho = 0
          
          IF (.NOT. EZspline_allocated(grho_fsa_spl)) CALL initialize_splines_fsa(error_status)
          IF (error_status < 0) RETURN
          
          IF (s_val >= 0 .and. s_val <= 1) THEN
             CALL EZspline_interp(grho_fsa_spl, s_val, fsa_gradrho, error_status)
          ELSE
             error_status = -2
          END IF

        END SUBROUTINE fsa_gradrho_from_s

        
        SUBROUTINE fsa_gradrho2_from_s(s_val, fsa_gradrho2, error_status)
          USE EZspline
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(in) :: s_val
          DOUBLE PRECISION, INTENT(out) :: fsa_gradrho2
          INTEGER, INTENT(out) :: error_status

          error_status = 0
          fsa_gradrho2 = 0
          
          IF (.NOT. EZspline_allocated(grho2_fsa_spl)) CALL initialize_splines_fsa(error_status)
          IF (error_status < 0) RETURN
          
          IF (s_val >= 0 .and. s_val <= 1) THEN
             CALL EZspline_interp(grho2_fsa_spl, s_val, fsa_gradrho2, error_status)
          ELSE
             error_status = -2
          END IF

        END SUBROUTINE fsa_gradrho2_from_s

        
        SUBROUTINE fsa_b2overgradrho_from_s(s_val, fsa_b2overgradrho, error_status)
          USE EZspline
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(in) :: s_val
          DOUBLE PRECISION, INTENT(out) :: fsa_b2overgradrho
          INTEGER, INTENT(out) :: error_status

          error_status = 0
          fsa_b2overgradrho = 0
          
          IF (.NOT. EZspline_allocated(b2overgrho_fsa_spl)) CALL initialize_splines_fsa(error_status)
          IF (error_status < 0) RETURN
          
          IF (s_val >= 0 .and. s_val <= 1) THEN
             CALL EZspline_interp(b2overgrho_fsa_spl, s_val, fsa_b2overgradrho, error_status)
          ELSE
             error_status = -2
          END IF

        END SUBROUTINE fsa_b2overgradrho_from_s
        
        
        SUBROUTINE cyl_from_car(point_car, point_cyl, error_status)
          ! CAR_TO_CYL converts from Cartesian to cylindrical coordinates
          !
          ! References:
          !   W.A.Houlberg, F90 free format 8/2004
          
          IMPLICIT NONE
          !Cartesian coordinates (x,y,z) [m,m,m]
          DOUBLE PRECISION, INTENT(IN) :: point_car(3)               
          !cylindrical coordinates (R,phi,Z) [m,rad,m]
          DOUBLE PRECISION, INTENT(OUT) :: point_cyl(3)               
          INTEGER, INTENT(out) ::  error_status

          error_status = 0
          point_cyl(:) = 0

          ! Cartesian to cylindrical conversion
          point_cyl(1)=SQRT(point_car(1)**2+point_car(2)**2) !R
          point_cyl(2)=ATAN2(point_car(2),point_car(1))      !phi
          point_cyl(3)=point_car(3)                      !Z
          
          !Ensure 0 <= phi <= 2*pi
          IF (point_cyl(2) < 0.0) point_cyl(2)=point_cyl(2)+pi2
          
        END SUBROUTINE cyl_from_car


        SUBROUTINE car_from_cyl(point_cyl, point_car, error_status)
          ! CYL_TO_CAR converts from cylindrical to Cartesian coordinates
          !
          ! References:
          !   W.A.Houlberg, F90 free format 8/2004
          
          IMPLICIT NONE
          !cylindrical coordinates (R,phi,Z) [m,rad,m]
          DOUBLE PRECISION, INTENT(IN) :: point_cyl(3)               
          !Cartesian coordinates (x,y,z) [m,m,m]
          DOUBLE PRECISION, INTENT(OUT) :: point_car(3)
          INTEGER, INTENT(out) ::  error_status

          error_status = 0
          point_car(:) = 0

          ! Cylindrical to Cartesian conversion
          point_car(1) = point_cyl(1)*COS(point_cyl(2))   !x
          point_car(2) = point_cyl(1)*SIN(point_cyl(2))   !y
          point_car(3) = point_cyl(3)                     !z

        END SUBROUTINE car_from_cyl

        
        !-----------------------------------------------------------------------
        SUBROUTINE mntouv(k1,k,mnmax,nu,nv,xu,xv,fmn,xm,xn,f,signs,calc_trig)
          IMPLICIT NONE
          ! INPUT VARIABLES
          INTEGER, INTENT(in) :: k1
          INTEGER, INTENT(in) :: k
          INTEGER, INTENT(in) :: mnmax
          INTEGER, INTENT(in) :: nu
          INTEGER, INTENT(in) :: nv
          DOUBLE PRECISION, INTENT(in) :: xu(1:nu)
          DOUBLE PRECISION, INTENT(in) :: xv(1:nv)           
          DOUBLE PRECISION, INTENT(in) :: fmn(1:mnmax,k1:k)
          INTEGER, INTENT(in) :: xm(1:mnmax)
          INTEGER, INTENT(in) :: xn(1:mnmax)
          DOUBLE PRECISION, INTENT(inout) :: f(1:nu,1:nv,k1:k)
          INTEGER, INTENT(in) :: signs
          INTEGER, INTENT(in) :: calc_trig
          ! LOCAL VARIABLES
          INTEGER     :: mn, i, ier, ik
          DOUBLE PRECISION :: xm_temp(1:mnmax,1)
          DOUBLE PRECISION :: xn_temp(1:mnmax,1)
          DOUBLE PRECISION :: mt(1:mnmax,1:nu)
          DOUBLE PRECISION :: nz(1:mnmax,1:nv)
          DOUBLE PRECISION :: fmn_temp(1:mnmax,1:nu)
          DOUBLE PRECISION :: xu_temp(1,1:nu)
          DOUBLE PRECISION :: xv_temp(1,1:nv)
          DOUBLE PRECISION :: fmn_help(1:mnmax)
          DOUBLE PRECISION, ALLOCATABLE, SAVE :: cosmt(:,:)
          DOUBLE PRECISION, ALLOCATABLE, SAVE :: sinmt(:,:)
          DOUBLE PRECISION, ALLOCATABLE, SAVE :: cosnz(:,:)
          DOUBLE PRECISION, ALLOCATABLE, SAVE :: sinnz(:,:)
          ! BEGIN SUBROUTINE
          IF (calc_trig == 1) THEN
             IF (ALLOCATED(cosmt)) DEALLOCATE(cosmt)
             IF (ALLOCATED(sinmt)) DEALLOCATE(sinmt)
             IF (ALLOCATED(cosnz)) DEALLOCATE(cosnz)
             IF (ALLOCATED(sinnz)) DEALLOCATE(sinnz)
             ALLOCATE(cosmt(1:mnmax,1:nu),sinmt(1:mnmax,1:nu),&
                  cosnz(1:mnmax,1:nv),sinnz(1:mnmax,1:nv),STAT=ier)
             FORALL(i=1:mnmax) xm_temp(i,1)=DBLE(xm(i))
             FORALL(i=1:mnmax) xn_temp(i,1)=DBLE(xn(i))
             FORALL(i=1:nu) xu_temp(1,i)=xu(i)
             FORALL(i=1:nv) xv_temp(1,i)=xv(i)
             mt = MATMUL(xm_temp,xu_temp)
             nz = MATMUL(xn_temp,xv_temp)
             FORALL(mn=1:mnmax,i=1:nu) cosmt(mn,i) = dcos(pi2*mt(mn,i))
             FORALL(mn=1:mnmax,i=1:nu) sinmt(mn,i) = dsin(pi2*mt(mn,i))
             FORALL(mn=1:mnmax,i=1:nv) cosnz(mn,i) = dcos(pi2*nz(mn,i))
             FORALL(mn=1:mnmax,i=1:nv) sinnz(mn,i) = dsin(pi2*nz(mn,i))
          END IF
          IF (SIGNS == 0) THEN
             DO ik = k1,k
                FORALL(mn=1:mnmax) fmn_help(mn)=fmn(mn,ik)
                fmn_temp=SPREAD(fmn_help,2,nu)
                f(1:nu,1:nv,ik) = f(1:nu,1:nv,ik)  + MATMUL(TRANSPOSE((fmn_temp*cosmt)),cosnz) &
                     - MATMUL(TRANSPOSE((fmn_temp*sinmt)),sinnz)
             END DO
          ELSE IF (SIGNS == 1) THEN
             DO ik = k1,k
                FORALL(mn=1:mnmax) fmn_help(mn)=fmn(mn,ik)
                fmn_temp=SPREAD(fmn_help,2,nu)
                f(1:nu,1:nv,ik) = f(1:nu,1:nv,ik) + MATMUL(TRANSPOSE((fmn_temp*sinmt)),cosnz) &
                     + MATMUL(TRANSPOSE((fmn_temp*cosmt)),sinnz)
             END DO
          END IF

        END SUBROUTINE mntouv

      END MODULE mir_tools
      
