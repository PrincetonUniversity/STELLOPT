      SUBROUTINE read_wout_booz(extension, iread, ierr)
      USE read_wout_mod
      USE booz_params, mpol_b => mpol, nfp_b => nfp, ns_b => ns,
     1   ntor_b => ntor, mnmax_b => mnmax, phip_b => phip,
     2   mpol_nyq_b => mpol_nyq, ntor_nyq_b => ntor_nyq,
     3   mnmax_nyq_b => mnmax_nyq, 
     4   pres_b => pres, beta_b => beta_vol, phi_b => phi,
     5   buco_b => buco, bvco_b => bvco,
     6   bsubumnc_b => bsubumnc, bsubvmnc_b => bsubvmnc,
     6   bsubumns_b => bsubumns, bsubvmns_b => bsubvmns
      USE booz_persistent, rmnc_b => rmnc, zmns_b => zmns,
     1   lmns_b => lmns, rmns_b => rmns, zmnc_b => zmnc,
     2   lmnc_b => lmnc, 
     3   xm_b => xm, xn_b => xn, xm_nyq_b => xm_nyq,
     4   xn_nyq_b => xn_nyq
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: ierr, iread
      CHARACTER(LEN=*) :: extension
C-----------------------------------------------
      IF (.not.lwout_opened) THEN
         CALL read_wout_file(extension, ierr)

         IF (ierr .ne. 0) THEN
            PRINT *,' ierr = ', ierr,
     1   ' error in read_wout_file called from xbooz_xform'
         END IF
         IF (ierr_vmec.ne.0 .or. ierr.ne.0) GOTO 1000
      END IF

!
!     Load booz_params module from wout file module
!
      mpol_b = mpol;  mpol_nyq_b = mnyq
      mpol1  = mpol - 1
      ntor_b = ntor;  ntor_nyq_b = nnyq
      mnmax_b= mnmax; mnmax_nyq_b = mnmax_nyq
      nfp_b  = nfp
      ns_b   = ns
      lasym_b= lasym

!
!     COMPUTE ACTUAL NO. THETA, PHI POINTS FOR INTEGRATIONS
!     NEEDED FOR DYNAMIC MEMORY ALLOCATION
!
      mboz = MAX(6*mpol, 2, mboz)             
      nboz = MAX(2*ntor-1, 0, nboz)           
      nu_boz   = 2*(2*mboz+1)                     !CHANGED THIS FROM 2*(3*mboz+1)
      nv_boz   = 2*(2*nboz+1)                     !CHANGED THIS FROM 2*(2*nboz+1)
      IF (nboz .eq. 0) nv_boz = 1
!      nu_boz   = nu_boz + MOD(nu_boz,2)          !nu_boz, nv_boz MUST be even
!      nv_boz   = nv_boz + MOD(nv_boz,2)
      nunv = nu_boz*nv_boz
      mnboz = nboz+1 + (mboz-1)*(1+2*nboz)
      nu2_b = nu_boz/2+1                          !pi

!
!     ALLOCATE ARRAYS FIRST TIME THRU AND THEN STORE INPUT DATA IN PERSISTENT ARRAYS
!
      ! CALL allocate_boozer (iread)
      ! Modified so function can be called externally
      CALL allocate_boozer (iread)

      xm_b(:mnmax) = xm(:mnmax)
      xn_b(:mnmax) = xn(:mnmax)
      xm_nyq_b(:mnmax_nyq) = xm_nyq(:mnmax_nyq)
      xn_nyq_b(:mnmax_nyq) = xn_nyq(:mnmax_nyq)

      rmnc_b(:mnmax,:ns) = rmnc(:mnmax,:ns)
      zmns_b(:mnmax,:ns) = zmns(:mnmax,:ns)
      lmns_b(:mnmax,:ns) = lmns(:mnmax,:ns)

      IF (lasym_b) THEN
         rmns_b(:mnmax,:ns) = rmns(:mnmax,:ns)
         zmnc_b(:mnmax,:ns) = zmnc(:mnmax,:ns)
         lmnc_b(:mnmax,:ns) = lmnc(:mnmax,:ns)
      ENDIF

      bsubumnc_b(:mnmax_nyq,:ns) = bsubumnc(:mnmax_nyq,:ns)
      bsubvmnc_b(:mnmax_nyq,:ns) = bsubvmnc(:mnmax_nyq,:ns)
      bmodmnc(:mnmax_nyq,:ns)  = bmnc(:mnmax_nyq,:ns)

      IF (lasym_b) THEN
         bsubumns_b(:mnmax_nyq,:ns) = bsubumns(:mnmax_nyq,:ns)
         bsubvmns_b(:mnmax_nyq,:ns) = bsubvmns(:mnmax_nyq,:ns)
         bmodmns(:mnmax_nyq,:ns)  = bmns(:mnmax_nyq,:ns)
      END IF

      hiota(2:ns) = iotas(2:ns)
      phip_b(2:ns) = phip(2:ns)
      pres_b(2:ns) = pres(2:ns)
      beta_b(2:ns) = beta_vol(2:ns)
      phi_b(2:ns) = phi(2:ns)
      buco_b(2:ns) = buco(2:ns)
      bvco_b(2:ns) = bvco(2:ns)

!
!     Deallocate memory in read_wout module
!     (Don't do this so we can call this from STELLOPT as a function
! 1000 CALL read_wout_deallocate
 1000 IF (ierr .eq. 0) THEN
         ierr = -ierr_vmec
      ELSE
         PRINT *,' ierr = ', ierr, ' writing in READ_WOUT_BOOZ'
      END IF

      END SUBROUTINE read_wout_booz
