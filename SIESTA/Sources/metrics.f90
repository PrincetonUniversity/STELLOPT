!>
!!  \brief Module for reading in VMEC input and computing metric elements
!!   on the SIESTA (sqrt-flux) radial grid.
!<
      MODULE metrics
      USE stel_kinds
      USE stel_constants
      USE island_params
      USE read_wout_mod, ns_vmec=>ns, mpol_vmec=>mpol, ntor_vmec=>ntor, &
          rmnc_vmec=>rmnc, zmns_vmec=>zmns, lmns_vmec=>lmns,            &
          xm_vmec=>xm, xn_vmec=>xn, iotaf_vmec=>iotaf,                  &
          phipf_vmec=>phipf, presf_vmec=>presf, nfp_vmec=>nfp,          &
          wb_vmec=>wb, wp_vmec=>wp, gamma_vmec=>gamma
      USE descriptor_mod, ONLY: iam
      USE timer_mod

      IMPLICIT NONE
!     
!     WRITTEN 06-27-06 BY S. P. HIRSHMAN AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: COMPUTES AND STORES THE REAL-SPACE METRIC ELEMENTS, JACOBIAN BASED 
!     ON A SQRT(FLUX) SPLINED VMEC - COORDINATE SYSTEM
!

!     VARIABLE DECLARATIONS (3D for now; will convert to 2D or 1D as needed)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::                       &
                   sqrtg,                                               & !sqrt(g): Jacobian on half grid
                   gss, gsu, gsv, guu, guv, gvv,                        & !symmetric elements of lower metric tensor (half mesh)
                   hss, hsu, hsv, huu, huv, hvv                           !symmetric elements of upper metric tensor (full mesh)

      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::                       &
                   rmnc_spline, zmns_spline, lmns_spline

      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE ::                     &
                   rmnc_i, zmns_i, lmns_i

      REAL(rprec) :: rmax, rmin, zmax, zmin

      REAL(rprec), DIMENSION(:,:), ALLOCATABLE:: gssf, guuf, gvvf,      & !Lower Metric elements (full)
                   gsuf, gsvf, guvf   

      INTEGER :: iflipj=1                                                 !=-1 to flip sign of vmec jac
!
!     LOCAL (PRIVATE) HELPER ROUTINES (NOT ACCESSIBLE FROM OUTSIDE THIS MODULE)
!
      PRIVATE spline_fourier_modes, loadrzl_vmec, add_ghost_points
      INTEGER, PARAMETER, PRIVATE :: iScale=0                             !=0, s ~ SQRT(PHI) - default; = 1, s ~ PHI (VMEC)

      REAL(dp) :: skston, skstoff

      CONTAINS

      SUBROUTINE init_metric_elements(ns_in, mpol_in, ntor_in, wout_file)
      USE timer_mod, ONLY: init_timers
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER*(*)            :: wout_file
      INTEGER, INTENT(in)      :: ns_in, mpol_in, ntor_in
      INTEGER                  :: istat, ntype, imesh, js
      REAL(rprec)              :: t1, t2
      INTEGER                  :: m, n, mn
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec),ALLOCATABLE  :: lmns_half(:,:)
!-----------------------------------------------
!
!     LOADS VALUES FOR MESHES TO BE USED IN VMECPP ISLAND SOLVER
!     (GENERALLY DIFFERENT FROM VMEC MESHES. USE SUBSCRIPT _i FOR ISLAND VARIABLES)
!     READS wout_file FROM VMEC TO GET FOURIER COMPONENTS ON VMEC MESH
!     SPLINES THE VMEC COMPONENTS TO RADIAL ISLAND MESH (PHI->SQRT(PHI))
!     CALLS Fourier MODULE fixarray TO LOAD TRIG ARRAYS FOR COMPUTING ISLAND METRICS
!     COMPUTES gij (sub/sup) AND JACOBIAN ON ISLAND MESHES IN REAL SPACE
!

#if defined(SKS) 
      CALL second0(skston)
#endif
      CALL init_timers
#if defined(SKS) 
      CALL second0(skstoff)
      init_timers_time=init_timers_time+(skstoff-skston)
#endif

      IF (iScale .NE. 0) THEN
      PRINT *,'RESET ISCALE=0 AFTER TESTING!'
      !PAUSE
      END IF

      ns_i = ns_in
      nsh  = ns_i-1
      mpol_i = mpol_in
      ntor_i = ntor_in

!Set number of points == number of modes for now! (May want mid-points for flux conservation)
      nu_i = mpol_i+2  
      nv_i = 2*ntor_i+2

!USE 3/2 (ORSZAG) RULE FOR ANTI-ALIASING OF EVOLUTION EQNS
!Suppresses RADIAL grid separation in pressure
      nu_i = 3*nu_i/2
      nu_i = nu_i+MOD(nu_i,2)

      nv_i = 3*nv_i/2
      nv_i = nv_i+MOD(nv_i,2)

      IF (ntor_i .eq. 0) nv_i = 1
      nuv_i = nu_i*nv_i
      mnmax_i = (mpol_i + 1)*(2*ntor_i + 1)             ! Added RS. Contains total number of modes.

!      mpol32_i = nu_i                                  ! Extended mesh dims for evolved jpmnch
!      ntor32_i = nv_i/2

!
!     READ-IN DATA FROM VMEC-PRODUCED WOUT FILE (LIBSTELL ROUTINE)
!
#if defined(SKS) 
      CALL second0(skston)
#endif
      CALL read_wout_file(wout_file, istat)
#if defined(SKS) 
      CALL second0(skstoff)
      read_wout_file_time=read_wout_file_time+(skstoff-skston)
#endif
      IF (istat .ne. 0) STOP 'Read-wout error in INIT_METRIC_ELEMENTS'

      nfp_i = nfp_vmec
      wb_i  = (4*pi*pi)*wb_vmec
      wp_i  = (4*pi*pi)*wp_vmec

      IF (wb_i .eq. 0._dp) STOP 'wb_vmec = 0!'
      gnorm_i = wb_i + wp_i/(gamma-1)
      gnorm_i = 1._dp/gnorm_i

      IF (iam .eq. 0) THEN
         WRITE (33, *) 'INITIAL VMEC PARAMETERS'
         WRITE (33, 100) wp_vmec/wb_vmec, phi(ns_vmec)
      ENDIF

 100  FORMAT(' <BETA>: ',1pe12.4,' TFLUX: ',1pe12.4,/, 21('-'))
      
!
!     Allocate space for splined arrays
!     NOTE: keep lmns_spline on half mesh->ghost points at js=1, js=ns+1
!
      ALLOCATE(rmnc_spline(mnmax,ns_i), zmns_spline(mnmax,ns_i),        &
               lmns_spline(mnmax,ns_i+1), phipf_i(ns_i),                &
               chipf_i(ns_i), presf_i(ns_i), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error 1 in INIT_METRIC_ELEMENTS'

!     EVENTUALLY, WE SHOULD UPDATE WOUT TO WRITE OUT CHIPF. FOR NOW, WE COMPUTE IT
!     AND STORE IT IN IOTAF_VMEC
!
      iotaf_vmec = iotaf_vmec*phipf_vmec
!
!     SPLINE R, Z. L FOURIER COMPONENTS in s FROM ORIGINAL VMEC MESH (s ~ phi, ns_vmec points) 
!     TO A "POLAR" MESH [s ~ sqrt(phi), ns_i POINTS] WITH BETTER AXIS RESOLUTION
!
      DO ntype = 1, 3
         IF (ntype .eq. 1) THEN
            istat = 0
#if defined(SKS) 
            CALL second0(skston)
#endif
            CALL Spline_Fourier_Modes(rmnc_vmec, rmnc_spline,ns_vmec, ns_i, istat)
#if defined(SKS) 
            CALL second0(skstoff)
            Spline_Fourier_Modes_time=Spline_Fourier_Modes_time+(skstoff-skston)
#endif
            IF (istat .ne. 0) STOP 'Error splining rmnc'
         ELSE IF (ntype .eq. 2) THEN
            istat = 0
#if defined(SKS) 
            CALL second0(skston)
#endif
            CALL Spline_Fourier_Modes(zmns_vmec, zmns_spline,           &
                                      ns_vmec, ns_i, istat)   
#if defined(SKS) 
            CALL second0(skstoff)
            Spline_Fourier_Modes_time=Spline_Fourier_Modes_time+(skstoff-skston)
#endif
            IF (istat .ne. 0) STOP 'Error splining zmns'
         ELSE 
            ALLOCATE(lmns_half(mnmax,ns_vmec+1), stat=istat)
#if defined(SKS) 
            CALL second0(skston)
#endif
            CALL Add_Ghost_Points(lmns_vmec, lmns_half)
#if defined(SKS) 
            CALL second0(skstoff)
            Add_Ghost_Points_time=Add_Ghost_Points_time+(skstoff-skston)
#endif
            istat = 0
#if defined(SKS) 
            CALL second0(skston)
#endif
            CALL Spline_Fourier_Modes(lmns_half, lmns_spline,           &
                                      ns_vmec+1, ns_i+1, istat)   
#if defined(SKS) 
            CALL second0(skstoff)
            Spline_Fourier_Modes_time=Spline_Fourier_Modes_time+(skstoff-skston)
#endif
            IF (istat .ne. 0) STOP 'Error splining lmns'
            DEALLOCATE(lmns_half)
         END IF
      END DO
!
!     Spline 1-D arrays: careful -> convert phipf VMEC and multiply
!     by ds-vmec/ds-island, since phipf_i = d(PHI)/ds-island
!
#if defined(SKS) 
      CALL second0(skston)
#endif
      CALL Spline_OneD_Array (iotaf_vmec, chipf_i, istat)
      CALL Spline_OneD_Array (phipf_vmec, phipf_i, istat)
      presf_vmec = mu0 * presf_vmec
      CALL Spline_OneD_Array (presf_vmec, presf_i, istat)
#if defined(SKS) 
      CALL second0(skstoff)
      Spline_OneD_Array_time=Spline_OneD_Array_time+(skstoff-skston)
#endif
!
!     Scale phipf_i and convert to sqrt(flux) mesh by multiplying by 2*s
!
      phipf_i = phipf_i / (2*pi)
      chipf_i = chipf_i / (2*pi)
      IF (iScale == 0) THEN
         DO js = 1, ns_i
            phipf_i(js) = 2 * hs_i*(js-1) * phipf_i(js)
            chipf_i(js) = 2 * hs_i*(js-1) * chipf_i(js)
         END DO
      END IF
!
!     CONSTRUCT R, Z, L REAL-SPACE ARRAYS ON SQRT(FLUX) - "POLAR" - MESH
!     AND COMPUTE METRIC ELEMENTS AND JACOBIAN
!
#if defined(SKS) 
      CALL second0(skston)
#endif
      CALL LoadRZL_VMEC(istat)
#if defined(SKS) 
      CALL second0(skstoff)
      LoadRZL_VMEC_time=LoadRZL_VMEC_time+(skstoff-skston)
#endif
      IF (istat .ne. 0) STOP 'LoadRZL error in INIT_METRIC_ELEMENTS'

      END SUBROUTINE init_metric_elements
      
      SUBROUTINE Spline_Fourier_Modes(ymn_vmec, ymn_spline, nsv, nsi,   &
                                      istat)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(iN)        :: nsv, nsi
      INTEGER, INTENT(inout)     :: istat
      REAL(rprec), DIMENSION(mnmax, nsv), TARGET,                       &
                   INTENT(in)  :: ymn_vmec
      REAL(rprec), DIMENSION(mnmax, nsi), TARGET,                       &
                   INTENT(out) :: ymn_spline
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), PARAMETER          :: one=1
      INTEGER                         :: js, modes, ntype, mp
      REAL(rprec), DIMENSION(nsv)     :: snodes_vmec, y2_vmec, y_save, fac1, y_vmec
      REAL(rprec), DIMENSION(nsi)     :: snodes, fac2, y_spline
      REAL(rprec)                     :: hsv, yp1, ypn, expm=0
      LOGICAL                         :: lfull, lscale
!-----------------------------------------------
!
!     CALL LIBRARY SPLINE ROUTINES TO SPLINE FROM VMEC (s~phi) TO s~sqrt(phi) MESH
!

!
!     1. Set up "knots" on initial (vmec, svmec ~ phi) mesh 
!        and factor for taking out (putting back) [sqrt(s)]**m factor for odd m
!

      IF (nsv .le. 1) THEN
         istat = 1
         RETURN
      END IF

      lfull = (nsv .eq. ns_vmec)
      lscale = (iscale .eq. 0)

      hsv = one/(ns_vmec-1)
      snodes_vmec(1) = 0;  fac1(1) = 1;
      DO js = 2, ns_vmec
         IF (lfull) THEN
            snodes_vmec(js) = hsv*(js-1)
         ELSE
            snodes_vmec(js) = hsv*(js-1.5_dp)
         END IF
         fac1(js) = one/SQRT(snodes_vmec(js))     !regularization factor: sqrt(FLUX) for odd modes
      END DO

      IF (.not.lfull) THEN
         snodes_vmec(nsv) = 1
         fac1(nsv) = 1;
      END IF

!
!     2. Set up s-nodes on final (snodes, splined) mesh [s_siesta(sfinal) ~ sfinal**2]
!        if s_siesta ~ sfinal, this is then the original vmec mesh
!
      IF (ns_i .le. 1) THEN
         istat = 2
         RETURN
      END IF

      ohs_i = ns_i-1
      hs_i = one/ohs_i
      fac2(1) = 0;  snodes(1) = 0
      DO js = 2, ns_i
         IF (lfull) THEN
            fac2(js) = hs_i*(js-1)
         ELSE
            fac2(js) = hs_i*(js-1.5_dp)
         END IF
         IF (lscale) THEN
            snodes(js) = fac2(js)*fac2(js)      !SIESTA s==fac2 ~ sqrt(s-vmec) mesh
         ELSE
            snodes(js) = fac2(js)               !SIESTA s==fac2 ~ s-vmec mesh
         END IF
         fac2(js) = SQRT(snodes(js))
      END DO

      IF (.not.lfull) THEN
         fac2(nsi) = 1
         snodes(nsi) = 1
      END IF

      DO modes = 1, mnmax
         y_vmec = ymn_vmec(modes,:)
         mp = xm_vmec(modes)

         IF (istat.eq.0 .and. mp.gt.0) THEN
            IF (MOD(mp,2) .eq. 1) THEN 
               expm = 1
            ELSE
               expm = 2
            END IF
            IF (expm .ne. 0) y_vmec = y_vmec*(fac1**expm)
            IF (mp .le. 2) y_vmec(1) = 2*y_vmec(2) - y_vmec(3)
         END IF

!
!      3. Initialize spline for each mode amplitude (factor out sqrt(s) factor for odd-m)
!
         yp1 = -1.e30_dp;  ypn = -1.e30_dp
         CALL spline (snodes_vmec, y_vmec, nsv, yp1, ypn, y2_vmec)

!
!      4. Interpolate onto snodes mesh
!
         DO js = 1, nsi
            CALL splint (snodes_vmec, y_vmec, y2_vmec, nsv,             &
                         snodes(js), y_spline(js))
         END DO

         IF (istat.eq.0 .and. mp.gt.0) THEN
            IF (expm .ne. 0) y_spline = y_spline*(fac2**expm)
         END IF

         ymn_spline(modes,:) = y_spline(:)


!IDENTITY TEST THAT INTERPOLATION WORKS
         IF (nsi == nsv .and. .not.lscale) THEN
            IF (ANY(ABS(ymn_vmec(modes,2:) - ymn_spline(modes,2:)).gt.1.E-10_dp)) THEN
               PRINT *,'MODE: ',modes,'y_vmec .ne. y_spline'
            ENDIF
         END IF
!
!        PRINT OUT FOR CHECKING
!
!         IF (xn_vmec(modes) .eq. 0) THEN
!             WRITE (36, *) mp, INT(xn_vmec(modes))
!             DO js = 1, nsv/10
!               WRITE (36, '(1p2e14.4)') snodes_vmec(js), ymn_vmec(modes,js)
!            END DO
!             WRITE (36, *)
!            DO js = 1, nsi/3
!               WRITE (36, '(1p2e14.4)') snodes(js), ymn_spline(modes,js)
!            END DO
!             WRITE (36, *)
!         END IF

      END DO

      istat = 0

      END SUBROUTINE Spline_Fourier_Modes
      

      SUBROUTINE Spline_OneD_Array (y_vmec, y_spline, istat)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(out)     :: istat
      REAL(rprec), DIMENSION(ns_vmec), INTENT(in)  :: y_vmec
      REAL(rprec), DIMENSION(ns_i)   , INTENT(out) :: y_spline
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), PARAMETER          :: one = 1
      INTEGER                         :: js, modes, ntype, mp
      REAL(rprec), DIMENSION(ns_vmec) :: snodes_vmec, y2_vmec
      REAL(rprec), DIMENSION(ns_i)    :: snodes, fac2
      REAL(rprec)                     :: hs_vmec, yp1, ypn
!-----------------------------------------------
!
!     CALL LIBRARY SPLINE ROUTINES TO SPLINE FROM VMEC (s~phi) TO s~sqrt(phi) MESH
!

!
!     1. Set up "knots" on initial (vmec, svmec ~ phi) mesh 
!
      IF (ns_vmec .le. 1) THEN
         istat = 1
         RETURN
      END IF

      hs_vmec = one/(ns_vmec-1)
      DO js = 1, ns_vmec
         snodes_vmec(js) = hs_vmec*(js-1)
      END DO

!
!     2. Set up s-nodes on final (snodes, splined) mesh [s_siesta(sfinal) ~ sfinal**2]
!        if s_siesta ~ sfinal, this is then the original vmec mesh
!
      IF (ns_i .le. 1) THEN
         istat = 2
         RETURN
      END IF

      ohs_i = ns_i-1
      hs_i = one/ohs_i
      DO js = 1, ns_i
         fac2(js) = hs_i*(js-1)
         IF (iScale .EQ. 0) THEN
            snodes(js) = fac2(js)*fac2(js)        !polar mesh, s-siesta~s-vmec**2
         ELSE
            snodes(js) = fac2(js)                 !vmec mesh   s-siesta~s-vmec
         END IF
      END DO



!     4. Initialize spline for each mode amplitude (factor out sqrt(s) factor for odd-m)
!
      yp1 = -1.e30_dp;  ypn = -1.e30_dp
      CALL spline (snodes_vmec, y_vmec, ns_vmec, yp1, ypn, y2_vmec)

!
!     5. Interpolate onto snodes mesh
!
      DO js = 1, ns_i
         CALL splint (snodes_vmec, y_vmec, y2_vmec, ns_vmec,         &
                      snodes(js), y_spline(js))
      END DO


      END SUBROUTINE Spline_OneD_Array

      
      SUBROUTINE Add_Ghost_Points(ymn_vmec, ymn_ghost)
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), INTENT(in)            :: ymn_vmec(mnmax,ns_vmec)
      REAL(rprec), INTENT(out)           :: ymn_ghost(mnmax,ns_vmec+1)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), PARAMETER             :: p5 = 0.5_dp
      INTEGER                            :: modes, js, mp
      REAL(rprec), DIMENSION(ns_vmec)    :: y_vmec
!-----------------------------------------------

      ymn_ghost(:,2:ns_vmec) = ymn_vmec(:,2:ns_vmec)
!
!     ADDS GHOST POINTS AT js=1 (s=0) AND js=ns+1 (s=1)
!
      ymn_ghost(:,ns_vmec+1) = 2*ymn_vmec(:,ns_vmec) - ymn_vmec(:,ns_vmec-1)

      DO modes = 1, mnmax
         mp = xm_vmec(modes)
         IF (mp .eq. 0) THEN
            ymn_ghost(modes,1) = 2*ymn_vmec(modes,2)-ymn_vmec(modes,3)
         ELSE
            ymn_ghost(modes,1) = 0
         END IF
      END DO

      END SUBROUTINE


      SUBROUTINE LoadRZL_VMEC(istat)
      USE fourier, ONLY: toijsp, init_fourier
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(out)     :: istat
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), DIMENSION(ns_vmec)  :: fac1, fac2
      INTEGER                          :: iparity
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE    ::                  &
                   r1_i, z1_i, ru_i, zu_i, rv_i, zv_i
!-----------------------------------------------
!
!     1. LOAD FOURIER FACTORS USING VMEC mpol,ntor AND ISLAND nu,nv
!        IDEA: USE VMEC mpol, ntor VALUES AND CALL fixarray, FLIP VMEC n -> -n in REPACK
!     2. COMPUTE METRIC ELEMENTS, JACOBIAN, LAMBDA ON ISLAND MESH
      
      CALL init_fourier

!
!     COMPUTE R, Z (Metric elements, jacobian), LAMBDA (For initializing B)
!
!     FIRST, REPACK SPLINED ARRAYS FROM VMEC-ORDERING TO ISLAND-ORDERING
!     AND SWAP THE N->-N MODES TO BE CONSISTENT WITH mu+nv ISLAND ARGUMENT
!
      CALL repack (istat)

!
!     COMPUTE AND STORE R, Z AND THEIR ANGULAR DERIVATIVES, AND LAMBDA (NEED FOR VECTOR POT)
!
      ALLOCATE (r1_i(nu_i, nv_i, ns_i), z1_i(nu_i, nv_i, ns_i),         &
                ru_i(nu_i, nv_i, ns_i), zu_i(nu_i, nv_i, ns_i),         &
                rv_i(nu_i, nv_i, ns_i), zv_i(nu_i, nv_i, ns_i),         &
                stat = istat)

      IF (istat .ne. 0) STOP 'Allocation failed in LoadRZL_VMEC'

      iparity = 0
      CALL toijsp (rmnc_i, r1_i, 0, 0, iparity, 0)
      CALL toijsp (rmnc_i, ru_i, 1, 0, iparity, 0)
      CALL toijsp (rmnc_i, rv_i, 0, 1, iparity, 0)
      
      iparity = 1
      CALL toijsp (zmns_i, z1_i, 0, 0, iparity, 0)
      CALL toijsp (zmns_i, zu_i, 1, 0, iparity, 0)
      CALL toijsp (zmns_i, zv_i, 0, 1, iparity, 0)

      rmax = MAXVAL(r1_i);      rmin = MINVAL(r1_i)
      zmax = MAXVAL(z1_i);      zmin = MINVAL(z1_i)

!
!     COMPUTE HALF-MESH LOWER/UPPER METRIC ELEMENTS AND JACOBIAN
!
      CALL half_mesh_metrics (r1_i, ru_i, rv_i, z1_i, zu_i, zv_i)

!
!     COMPUTE FULL-MESH METRIC ELEMENTS
      CALL full_mesh_metrics

!
!     CLEAN-UP EXTRA ARRAYS
!
      DEALLOCATE (r1_i, z1_i, ru_i, zu_i, rv_i, zv_i, stat=istat)

      END SUBROUTINE LoadRZL_VMEC


      SUBROUTINE repack (istat)
      USE fourier, ONLY: orthonorm
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(out)     :: istat
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER                  :: modes, js, m, n, ntype, n1
!-----------------------------------------------
!
!     The splined arrays (rmnc_spline, etc) are all in VMEC ordering
!
!         rmnc_spline(1:mnmax,ns_i)
!
!     Now pack them so they can be used by ISLAND Fourier routines
!
!         rmnc_i(ns_i, 0:mpol, -ntor:ntor)
!
!     NOTE: mpol_i == mpol_vmec, ntor_i == ntor_vmec here
!
!     KEEP lmns_i ON HALF MESH AND ADD EXTRAPOLATION POINTS AT 1 (rho=0)
!          AND NS+1 (=1)!
!
      ALLOCATE(rmnc_i(0:mpol_i,-ntor_i:ntor_i,ns_i),                    & 
               zmns_i(0:mpol_i,-ntor_i:ntor_i,ns_i),                    &
               lmns_i(0:mpol_i,-ntor_i:ntor_i,ns_i+1),                  &
               stat=istat)

      IF (istat .ne. 0) STOP 'Allocation error in REPACK'
      rmnc_i = 0;  zmns_i = 0;  lmns_i = 0
      
      IF (iam .eq. 0) THEN
         IF (MAXVAL(xm_vmec).gt.mpol_i .or.                            &
            MAXVAL(ABS(xn_vmec))/nfp_vmec.gt.ntor_i)                   & 
            PRINT *,'It is recommended to increase the number of ',    &
                    'modes to at least VMEC values!'
      ENDIF
!
!     FLIP VMEC n -> SIESTA -n SO TRIG ARG IS mu+nv IN SIESTA
!     SET iflipj -> -1 FOR POSITIVE JACOBIAN SYSTEM
!
      DO modes = 1,mnmax
         m = xm_vmec(modes)
         n = xn_vmec(modes)/nfp_vmec
         IF (m.gt.mpol_i .or. ABS(n).gt.ntor_i) CYCLE
!
!     LOAD n>=0 ONLY FOR M=0
!
         IF (m .eq. 0) THEN
            n1 = ABS(n)
            rmnc_i(m,n1,:) = rmnc_i(m,n1,:)                            &
                           + rmnc_spline(modes,:)
            zmns_i(m,n1,:) = zmns_i(m,n1,:)                            &
                           -SIGN(1,n)*zmns_spline(modes,:)
            lmns_i(m,n1,:) = lmns_i(m,n1,:)                            &  
                           -SIGN(1,n)*lmns_spline(modes,:)*iflipj
         ELSE
            rmnc_i(m,-n*iflipj,:) = rmnc_spline(modes,:)
            zmns_i(m,-n*iflipj,:) = zmns_spline(modes,:)*iflipj
            lmns_i(m,-n*iflipj,:) = lmns_spline(modes,:)            !Implied phip in here?
         ENDIF
      END DO
!
! RS (10/03/06)  Divide out "orthonorm" factor used in Island FFT routines
! 
      DO js = 1, ns_i
         rmnc_i(:,:,js) = rmnc_i(:,:,js)/orthonorm(0:mpol_i,-ntor_i:ntor_i)
         zmns_i(:,:,js) = zmns_i(:,:,js)/orthonorm(0:mpol_i,-ntor_i:ntor_i)
         lmns_i(:,:,js) = lmns_i(:,:,js)/orthonorm(0:mpol_i,-ntor_i:ntor_i)
      END DO

      DEALLOCATE (rmnc_spline, zmns_spline, lmns_spline, stat=istat)

      END SUBROUTINE repack


      SUBROUTINE half_mesh_metrics (r1_i, ru_i, rv_i, z1_i, zu_i, zv_i)
      USE island_params, ns=>ns_i, nuv=>nuv_i
      IMPLICIT NONE  
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(nuv,ns), INTENT(in) ::                    &
         r1_i, ru_i, rv_i, z1_i, zu_i, zv_i
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: p5 = 1.d0/2.d0, zero = 0
      INTEGER  :: js, lk, js1, istat
      REAL(rprec)            :: mintest, maxtest
      REAL(rprec), DIMENSION(nuv) ::                                    &
                   r12, ru12, rv12, zu12, zv12, rs12, zs12
!-----------------------------------------------
!
!     ALLOCATE METRIC ELEMENT ARRAYS
!
      ALLOCATE(sqrtg(nuv,ns),                                           &
               gss(nuv,ns), gsu(nuv,ns), gsv(nuv,ns),                   &
               guu(nuv,ns), guv(nuv,ns), gvv(nuv,ns), stat=istat)
      IF (istat .ne. 0) STOP 'ALLOCATION ERROR1 IN HALF_MESH_METRICS'

!     COMPUTE ALL ON THE HALF MESH
      sqrtg(:,1) = 0

      DO js = 2, ns
         js1 = js-1
         r12 = p5*(r1_i(:,js) + r1_i(:,js1))
         rs12 = (r1_i(:,js) - r1_i(:,js1))*ohs_i
         ru12= p5*(ru_i(:,js) + ru_i(:,js1))
         rv12= p5*(rv_i(:,js) + rv_i(:,js1))
         zs12 = (z1_i(:,js) - z1_i(:,js1))*ohs_i
         zu12= p5*(zu_i(:,js) + zu_i(:,js1))
         zv12= p5*(zv_i(:,js) + zv_i(:,js1))
         guu(:,js) = ru12*ru12 + zu12*zu12
         guv(:,js) = ru12*rv12 + zu12*zv12
         gvv(:,js) = r12*r12 + rv12*rv12 + zv12*zv12
         gsu(:,js) = rs12*ru12 + zs12*zu12
         gsv(:,js) = rs12*rv12 + zs12*zv12
         gss(:,js) = rs12*rs12 + zs12*zs12
         sqrtg(:,js) = r12*(ru12*zs12 - rs12*zu12)
      END DO

      IF (iScale .EQ. 1) THEN
         sqrtg(:,1) = sqrtg(:,2)
!!         gss(:,2) = gss(:,2)/2           
      END IF

!NEED THESE FOR CURRENT CALCULATION AT ORIGIN
      gss(:,1) = gss(:,2);  gsu(:,1) = 0;  gsv(:,1) = gsv(:,2)
      guu(:,1) = 0;         guv(:,1) = 0

      mintest = MINVAL(sqrtg(:,2:))
      maxtest = MAXVAL(sqrtg(:,2:))

      IF (mintest*maxtest .le. zero) THEN
         IF (iam .eq. 0) PRINT *,"Jacobian changed sign in half_mesh_metrics!"
         STOP 
      END IF

      END SUBROUTINE half_mesh_metrics


      SUBROUTINE full_mesh_metrics
      USE island_params, ns=>ns_i, nuv=>nuv_i, nzeta=>nv_i, ntheta=>nu_i
      IMPLICIT NONE  
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER:: istat, js
      REAL(dp)                                :: maxtest, eps, r1, s1, t1
      REAL(dp)                                :: det(nuv)
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: temps, tempu, tempv
!-----------------------------------------------
!
!     This subroutine gets the lower metric elements on the full mesh
!     Preserves positive definiteness of metric tensor
!
      ALLOCATE(gssf(nuv,ns), guuf(nuv,ns), gvvf(nuv,ns),                &
               gsuf(nuv,ns), gsvf(nuv,ns), guvf(nuv,ns),                &
               hss(nuv,ns),  huu(nuv,ns),  hvv(nuv,ns),                 & 
               hsu(nuv,ns),  hsv(nuv,ns),  huv(nuv,ns), stat=istat)

      IF (istat .ne. 0) STOP 'ALLOCATION ERROR IN FULL_MESH_METRICS'

      CALL to_full_mesh(gss,gssf,nuv)
      CALL to_full_mesh(guu,guuf,nuv)
      CALL to_full_mesh(gvv,gvvf,nuv)
      CALL to_full_mesh(gsu,gsuf,nuv)
      CALL to_full_mesh(gsv,gsvf,nuv)
      CALL to_full_mesh(guv,guvf,nuv)
      guuf(:,1) = 0;  gsuf(:,1) = 0;  guvf(:,1) = 0

!
!     Compute upper metric elements (inverse of lower matrix) on full mesh
!     NOTE: DET == 1/det{|Gij| == 1/sqrtg(full)**2
!
      hss(:,1) = 0;  hsu(:,1) = 0;  hsv(:,1) = 0
      huu(:,1) = 0;  huv(:,1) = 0;  hvv(1,:) = 0

      DO js = 2, ns
         det(:) = gssf(:,js)*(guuf(:,js)*gvvf(:,js)-guvf(:,js)*guvf(:,js)) &
                + gsuf(:,js)*(guvf(:,js)*gsvf(:,js)-gsuf(:,js)*gvvf(:,js)) &
                + gsvf(:,js)*(gsuf(:,js)*guvf(:,js)-gsvf(:,js)*guuf(:,js))
         IF (ANY(det .le. zero)) THEN
            IF (iam .eq. 0)PRINT *,'Determinant |gijf| <= 0 at js = ',js
            IF (js .eq. ns) THEN
               det(:) = sqrtg(:,ns)**2
               gssf(:,ns) = gss(:,ns)
               guuf(:,ns) = guu(:,ns)
               gvvf(:,ns) = gvv(:,ns)
               gsuf(:,ns) = gsu(:,ns)
               guvf(:,ns) = guv(:,ns)
               gsvf(:,ns) = gsv(:,ns)
            ELSE
               det(:) = ((sqrtg(:,js)+sqrtg(:,js+1))/2)**2
               STOP
            END IF
         END IF
         det = one/det
         hss(:,js) = det(:)*                                            &
                    (guuf(:,js)*gvvf(:,js) - guvf(:,js)*guvf(:,js))
         hsu(:,js) = det(:)*                                            &
                    (guvf(:,js)*gsvf(:,js) - gsuf(:,js)*gvvf(:,js))
         hsv(:,js) = det(:)*                                            &
                    (gsuf(:,js)*guvf(:,js) - gsvf(:,js)*guuf(:,js))
         huu(:,js) = det(:)*                                            &
                    (gssf(:,js)*gvvf(:,js) - gsvf(:,js)*gsvf(:,js))
         huv(:,js) = det(:)*                                            &
                    (gsvf(:,js)*gsuf(:,js) - gssf(:,js)*guvf(:,js))
         hvv(:,js) = det(:)*                                            &
                    (gssf(:,js)*guuf(:,js) - gsuf(:,js)*gsuf(:,js))

      END DO

!
!     CONFIRM ACCURACY OF INVERSE
!
      eps = 0.01_dp*SQRT(EPSILON(eps))

      ALLOCATE(temps(ntheta,nzeta,ns),                                  &
               tempu(ntheta,nzeta,ns),                                  &
               tempv(ntheta,nzeta,ns), stat=istat)

!     gijf * hj1   (first column of hij matrix)
      CALL tolowerf(hss, hsu, hsv, temps, tempu, tempv)
      r1 = MAXVAL(ABS(temps(:,:,2:)-1))
      s1 = MAXVAL(ABS(tempu(:,:,2:)))
      t1 = MAXVAL(ABS(tempv(:,:,2:)))
      maxtest = MAX(r1,s1,t1)
      IF (maxtest.GT.eps .AND. iam.EQ.0) STOP "ERROR1 IN FULL_MESH_METRICS"

!     gijf * hj2   (second column of hij matrix)
      CALL tolowerf(hsu, huu, huv, temps, tempu, tempv)
      r1 = MAXVAL(ABS(temps(:,:,2:)))
      s1 = MAXVAL(ABS(tempu(:,:,2:)-1))
      t1 = MAXVAL(ABS(tempv(:,:,2:)))
      maxtest = MAX(r1,s1,t1)
      IF (maxtest.GT.eps .AND. iam.EQ.0) STOP "ERROR2 IN FULL_MESH_METRICS"

!     gijf * hj3   (third column of hij matrix)
      CALL tolowerf(hsv, huv, hvv, temps, tempu, tempv)
      r1 = MAXVAL(ABS(temps(:,:,2:)))
      s1 = MAXVAL(ABS(tempu(:,:,2:)))
      t1 = MAXVAL(ABS(tempv(:,:,2:)-1))
      maxtest = MAX(r1,s1,t1)
      IF (maxtest.GT.eps .AND. iam.EQ.0) STOP "ERROR3 IN FULL_MESH_METRICS"
      
      DEALLOCATE(temps, tempu, tempv)

      END SUBROUTINE full_mesh_metrics


      SUBROUTINE cleanup_metric_elements
      USE fourier, ONLY: dealloc_fourier
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
	  INTEGER         :: istat
!-----------------------------------------------
	
!     Note: sqrtg is deallocated in init_bcovar and stored in jacob variable

      DEALLOCATE(gss, gsu, gsv, guu, guv, gvv,                          &
                 hss, hsu, hsv, huu, huv, hvv,                          &
                 phipf_i, chipf_i, presf_i, stat=istat)

      CALL dealloc_fourier
      DEALLOCATE(rmnc_i, zmns_i)     ! Used by output subroutines.        
	
      END SUBROUTINE cleanup_metric_elements

       
        SUBROUTINE dealloc_full_lower_metrics
        USE stel_kinds
        IMPLICIT NONE
        INTEGER :: istat
       
        DEALLOCATE(gssf, guuf, gvvf, gsuf, gsvf, guvf, stat = istat)
        IF (istat .ne. 0) STOP 'Problem dealloc. in LOWER_METRIC_FULL'
        
        END SUBROUTINE dealloc_full_lower_metrics

      END MODULE metrics
