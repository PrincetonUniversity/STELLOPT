      SUBROUTINE read_wout_opt(ictrl, extension, ierr, iopen)
      USE optim, ierr_vmec_opt=>ierr_vmec
      USE read_wout_mod, version_wout => version_
      USE vmec_input, ONLY: nfp_in => nfp, mpol_in => mpol,
     1   ntor_in => ntor, lfreeb_in => lfreeb, rbc, zbs, 
     2   raxis_cc, raxis_cs, zaxis_cs, zaxis_cc
      USE vparams, ONLY: mpold, ntord, zero, one
      USE optim_params, ONLY: lv3post
      USE boozer_params
      USE mpi_params                                                     !MPI
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(out) :: ierr, iopen
      INTEGER, INTENT(in)  :: ictrl
      CHARACTER(LEN=*), INTENT(in) :: extension
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: p03=3.e-2_dp, p5=0.5_dp, two=2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mn, js, n1, mb, nb
      REAL(rprec) :: dum
C-----------------------------------------------

!     COMPUTE ACTUAL NO. THETA, PHI POINTS FOR INTEGRATIONS
!     NEEDED FOR DYNAMIC MEMORY ALLOCATION

      IF (IAND(ictrl,1) == 0) THEN
         mboz = MAX(6*mpol_in,2,mboz_opt)        !USER MAY SPECIFY mboz_opt
         nboz = MAX(2*ntor_in - 1,1,nboz_opt)    !USER MAY SPECIFY nboz_opt
         nu_boz   = 6*mboz+2
         nv_boz   = 4*nboz+2
         nu_boz   = nu_boz + MOD(nu_boz,2)       !nu_boz, nv_boz MUST be even
         nv_boz   = nv_boz + MOD(nv_boz,2)
         mnboz = nboz + 1 + (mboz - 1)*(1 + 2*nboz)
         nu2_b = nu_boz/2 + 1
         nunv = nu2_b*nv_boz
!        nunv = nu_boz*nv_boz
         ierr = 0
         iopen = 0
         RETURN
      END IF

!
!     READ_WOUT_FILE: IN READ_WOUT_MOD.F (LIBSTELL/Modules)
!
      ierr_vmec = 46
      CALL read_wout_file (extension, ierr, iopen)
      ierr_vmec_opt = ierr_vmec
      CALL FLUSH(6)
      IF (iopen .ne. 0) RETURN

      IF (ierr_vmec.ne.0 .or. ierr.ne.0) THEN
         PRINT *,' For processor ', myid+1, ' extnsn=',TRIM(extension),
     1   ' ierr_vmec = ', ierr_vmec,' in read_wout_opt'
         PRINT *,' optimization will continue,',
     1   ' but this search direction will be deprecated!'
         IF (ierr_vmec_opt .eq. 0) ierr_vmec_opt = ierr
         GOTO 1000
      END IF

      version_opt = version_wout                                   !!COBRA
      wb_opt   = wb
      wp_opt   = wp
      rmax_opt = rmax_surf
      rmin_opt = rmin_surf
      zmax_opt = zmax_surf
      aminor_opt = aminor
      nfp_opt  = nfp
      mpol_opt = mpol
      ntor_opt = ntor
      mnmax_opt= mnmax
      mpol1_opt = mpol_opt - 1
      ntor1_opt = ntor_opt + 1
      curtor_opt = Itor
      nextcur_vmec = ABS(nextcur)

!
!     Confirm that WOUT file contains same data read in from input file
!
      IF (nfp_opt.ne.nfp_in .or. ns.ne.nrad .or. mpol_opt.ne.mpol_in
     1    .or. ntor_opt.ne.ntor_in) THEN
         IF (myid .eq. master) THEN                                      !START MPI
            PRINT *,' NFP_IN  = ',nfp_in, ' NFP_WOUT  = ',nfp_opt
            PRINT *,' NS_IN   = ',nrad,   ' NS_WOUT   = ',ns
            PRINT *,' MPOL_IN = ',mpol_in,' MPOL_WOUT = ',mpol_opt
            PRINT *,' NTOR_IN = ',ntor_in,' NTOR_WOUT = ',ntor_opt
         ENDIF                                                           !END MPI
         ierr = -1
         RETURN
      END IF
      IF (mpol_opt.gt.mpold .or. ntor_opt.gt.ntord)
     1  STOP 'mpold or ntord too small in read_wout'

      raxis_cs = 0;   zaxis_cc = 0

      DO mn = 1, mnmax
         xm_bdy(mn)   = xm(mn)
         xn_bdy(mn)   = xn(mn)
         IF (NINT(xm(mn)) .eq. 0) THEN
            n1 = ABS(NINT(xn(mn)))/nfp
            raxis_cc(n1) = rmnc(mn,1)
            zaxis_cs(n1) = zmns(mn,1)
            IF (ALLOCATED(rmns)) raxis_cs(n1) = rmns(mn,1)
            IF (ALLOCATED(zmnc)) zaxis_cc(n1) = zmnc(mn,1)
         END IF
         rmnc_bdy(mn) = rmnc(mn,nrad)
         zmns_bdy(mn) = zmns(mn,nrad)
      END DO

      rmnc_opt(:mnmax,:nrad) = rmnc(:mnmax,:nrad)
      zmns_opt(:mnmax,:nrad) = zmns(:mnmax,:nrad)
      lmns_opt(:mnmax,:nrad) = lmns(:mnmax,:nrad)
      
      IF (lasym) THEN
         rmns_opt(:mnmax,:nrad) = rmns(:mnmax,:nrad)
         zmnc_opt(:mnmax,:nrad) = zmnc(:mnmax,:nrad)
         lmnc_opt(:mnmax,:nrad) = lmnc(:mnmax,:nrad)
      END IF

!     Account for free boundary possibly changing boundary coefficients
!     need to write out correctly in write_rbzb routine. Only call from clean_up
      IF (lfreeb_in .and. IAND(ictrl,2) == 2) THEN
         rbc = 0; zbs = 0
         DO mn = 1, mnmax
            nb = NINT(xn_bdy(mn)/nfp)
            mb = NINT(xm_bdy(mn))
            rbc(nb,mb) = rmnc_bdy(mn)
            zbs(nb,mb) = zmns_bdy(mn)
         END DO
      END IF

      pres_opt(2:nrad) = pres(2:nrad)                              !!COBRA
      iota_opt(2:nrad) = iotas(2:nrad)
      phip_opt(2:nrad) = phip(2:nrad)
      buco_opt(2:nrad) = buco(2:nrad)
      vp_opt(2:nrad)   = vp(2:nrad)
      jcurv_opt(1:nrad) = jcurv(1:nrad)                            !!local <current density>, NOT integrated in s
      jdotb_opt(1:nrad) = jdotb(1:nrad)
      aspect_opt = aspect
      rbtor_opt  = rbtor

!
!     MERCIER CRITERION (NORM IT FOR USE IN OPTIMIZER)
!
      Dmerc_opt(2:nrad-1) = Dmerc(2:nrad-1)
      Dmerc_opt(1) = 0
      Dmerc_opt(nrad) = 0

      DO js = 2,nrad-1
         dum = Dshear(js) + ABS(Dwell(js)) + ABS(Dcurr(js))
         dum = 0.1_dp*MAX(dum, ABS(Dgeod(js)))
         IF (dum .gt. zero) Dmerc_opt(js) = Dmerc_opt(js)/dum
      END DO

!
!     IF POOR FORCE BALANCE, DO NOT BELIEVE MERCIER...
!
      equif(2:nrad-1) = ABS(equif(2:nrad-1))/p03
      WHERE (equif(2:nrad-1) .gt. one)
     1    Dmerc_opt(2:nrad-1) = Dmerc_opt(2:nrad-1) / equif(2:nrad-1)

!
!     PUT BUCO_OPT ON FULL MESH
!
      DO js = 2,nrad-1
        buco_opt(js) = p5*(buco_opt(js) + buco_opt(js+1))
      END DO

      buco_opt(1) = zero
      buco_opt(nrad) = two*buco_opt(nrad) - buco_opt(nrad-1)

!
!     DEALLOCATE MEMORY
!
 1000 CALL read_wout_deallocate

! If we are doing a magnetic reconstruction we need to penalize poor convergence
! May need MAX(fsqr,fsqz,fsql) .GT. ftolv) 
! but we will try MAX(fsqr,fsqz,fsql) .GT. 2 * ftolv) 
      IF (lv3post .and. (MAX(fsqr,fsqz,fsql) .gt. 2*ftolv )) ierr = -3
      IF (ierr_vmec.ne.0 .or. ierr.ne.0) THEN
         PRINT *,' For processor ', myid+1, ' extnsn=',TRIM(extension),
     1   ' MAX(fsqr,fsqz,fsql) > 2*ftol in read_wout_opt'
      ENDIF

      END SUBROUTINE read_wout_opt
