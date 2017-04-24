      SUBROUTINE get_ballooning_grate(grate)
      USE stel_kinds
      USE ballooning_data
      USE normalize_data
      USE readin_data
      USE general_dimensions
      USE fmesh_quantities
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), INTENT(OUT), DIMENSION(ns_cob) :: grate
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), DIMENSION(ns_cob) ::  radii_flux, h0, xmax, xmin
      REAL(rprec) :: fullp, semip, h, eigf, err, rerr, ohs2, ohs, pi
      INTEGER :: np, npm1, npm2, npm3, inc
      INTEGER :: j, kl, i, k
      INTEGER :: dinf, dinf1, dinfp, dinfp1, dinfpp, dinfm1
      REAL(rprec),ALLOCATABLE, DIMENSION(:) :: xf, xh, pf, ph,
     1    qf, qh, rf, rh, feigfun, eigfun, i_xf, i_xh, i_pf, i_ph,
     2    i_qf, i_qh, i_rf, i_rh, hv, eigv, a1f, a2f, a3f
      REAL(rprec) :: feigenv
      LOGICAL :: l_success
!-----------------------------------------------

      grate = 0
      DO i=1,ns_cob
!        radii(i)=SQRT(REAL(i-1,rprec)/(ns_cob-1)                          !  r/a= (radii_flux)**1/2
         radii_flux(i) = REAL(i-1,rprec)/(ns_cob-1)                        !  s = radii_flux
         radios(i) = radii_flux(i)
      ENDDO

!=================
!      BEGIN MIGRATION TO FULL MESH
!=================

!...  Store surface quantites on RADIAL full mesh

      ALLOCATE(iotaf(ns_cob), presf(ns_cob), phipf(ns_cob), stat=k)
      IF (k .ne. 0) STOP 'Allocation error 1 in get_ballooning_grate'

      iotaf = 0 ; presf = 0; phipf = 0

      iotaf(2:ns_cob-1) = 0.5_dp*(hiota(2:ns_cob-1) + hiota(3:ns_cob))
      presf(2:ns_cob-1) = 0.5_dp*(hpres(2:ns_cob-1) + hpres(3:ns_cob))
      phipf(2:ns_cob-1) = 0.5_dp*(hphip(2:ns_cob-1) + hphip(3:ns_cob))

!...  Evaluate and store surface quantities derivatives on RADIAL full mesh

      ALLOCATE(iotapf(ns_cob), prespf(ns_cob), stat=k)
      IF (k .ne. 0) STOP 'Allocation error 2 in get_ballooning_grate'

      iotapf = 0; prespf = 0

      ohs = ns_cob-1                                                       ! ds to differentiate on radial HALF mesh
      ohs2 = ohs/2                                                         ! 2 comes from differentiating on radial FULL mesh

      iotapf(2:ns_cob-1) = ohs*(hiota(3:ns_cob) - hiota(2:ns_cob-1))
      prespf(2:ns_cob-1) = ohs*(hpres(3:ns_cob) - hpres(2:ns_cob-1))


!...  EVALUATE AND STORE VMEC Fourier coefficients and
!           their derivatives on RADIAL full mesh

      ALLOCATE (lmnsf(mndim), bmncf(mndim), rmncpf(mndim),
     1          zmnspf(mndim), lmnspf(mndim), bmncpf(mndim),
     2          bsupvmncf(mndim), bsupumncf(mndim), stat=k)
      IF (k .ne. 0) STOP 'Allocation error 3 in get_ballooning_grate'
      lmnsf=0; bmncf=0; bsupvmncf=0; bsupumncf=0
      rmncpf=0; zmnspf=0; lmnspf=0; bmncpf=0

      IF (lasym_v) THEN                                                    ! 110909 RS = for ASYMMETRIC input
        ALLOCATE (lmncf(mndim), bmnsf(mndim), rmnspf(mndim),
     1            zmncpf(mndim), lmncpf(mndim), bmnspf(mndim),
     2            bsupvmnsf(mndim), bsupumnsf(mndim), stat=k)
        IF (k .ne. 0) STOP 'Allocation error 4 in get_ballooning_grate'
        lmncf=0; bmnsf=0; bsupvmnsf=0; bsupumnsf=0
        rmnspf=0; zmncpf=0; lmncpf=0; bmnspf=0
      ENDIF

      DO kl = 1, nlist
         dinf  = (list(kl)-1)*mnmax_v
         dinf1 = dinf+1
         dinfp = dinf+mnmax_v
         dinfp1= dinfp+1
         dinfpp= dinf+2*mnmax_v
         dinfm1= dinf-mnmax_v+1

!...   VMEC Fourier coefficients on RADIAL full mesh

         lmnsf(dinf1:dinfp) = 0.5_dp*(lmnsh(dinfp1:dinfpp)
     1    +lmnsh(dinf1:dinfp))
         bmncf(dinf1:dinfp) = 0.5_dp*(bmnch(dinfp1:dinfpp)
     1    +bmnch(dinf1:dinfp))
         bsupvmncf(dinf1:dinfp) = 0.5_dp*(bsupvmnch(dinfp1:dinfpp)
     1    +bsupvmnch(dinf1:dinfp))
         bsupumncf(dinf1:dinfp) = 0.5_dp*(bsupumnch(dinfp1:dinfpp)
     1    +bsupumnch(dinf1:dinfp))

!...   VMEC Fourier coefficients radial derivatives on RADIAL full mesh

         rmncpf(dinf1:dinfp) = ohs2*(rmncf(dinfp1:dinfpp)
     1    -rmncf(dinfm1:dinf))
         zmnspf(dinf1:dinfp) = ohs2*(zmnsf(dinfp1:dinfpp)
     1    -zmnsf(dinfm1:dinf))
         lmnspf(dinf1:dinfp) = ohs*(lmnsh(dinfp1:dinfpp)
     1    -lmnsh(dinf1:dinfp))
         bmncpf(dinf1:dinfp) = ohs*(bmnch(dinfp1:dinfpp)
     1    -bmnch(dinf1:dinfp))

         IF (lasym_v) THEN                                                 ! 110909 RS = for ASYMMETRIC input

!...   VMEC Asymmetric Fourier coefficients on RADIAL full mesh

           lmncf(dinf1:dinfp) = 0.5_dp*(lmnch(dinfp1:dinfpp)
     1      +lmnch(dinf1:dinfp))
           bmnsf(dinf1:dinfp) = 0.5_dp*(bmnsh(dinfp1:dinfpp)
     1      +bmnsh(dinf1:dinfp))
           bsupvmnsf(dinf1:dinfp) = 0.5_dp*(bsupvmnsh(dinfp1:dinfpp)
     1      +bsupvmnsh(dinf1:dinfp))
           bsupumnsf(dinf1:dinfp) = 0.5_dp*(bsupumnsh(dinfp1:dinfpp)
     1      +bsupumnsh(dinf1:dinfp))

!...   VMEC Asymmetric Fourier coefficients radial derivatives on RADIAL full mesh

           rmnspf(dinf1:dinfp) = ohs2*(rmnsf(dinfp1:dinfpp)
     1      -rmnsf(dinfm1:dinf))
           zmncpf(dinf1:dinfp) = ohs2*(zmncf(dinfp1:dinfpp)
     1      -zmncf(dinfm1:dinf))
           lmncpf(dinf1:dinfp) = ohs*(lmnch(dinfp1:dinfpp)
     1      -lmnch(dinf1:dinfp))
           bmnspf(dinf1:dinfp) = ohs*(bmnsh(dinfp1:dinfpp)
     1      -bmnsh(dinf1:dinfp))

         ENDIF

       ENDDO

       pi = 4*ATAN(one)
       xmax(1) = 0
       xmax(ns_cob) = 0
       xmax(2:ns_cob-1)=pi*(2*k_w-1)/(2*nfp_v*iotaf(2:ns_cob-1))           ! account for K_w potential wells!!!
       fullp = 360._dp/nfp_v
       semip = fullp/2

!...   Check for stellarator symmetry

       tsymm = 1
       IF (.NOT.lasym_v) THEN
         IF (l_geom_input .AND.
     1     ((init_zeta.eq.zero .or. init_zeta.eq.semip 
     2     .or. init_zeta.eq.fullp) .and. init_theta.eq.zero) ) 
     3     tsymm = 0
         IF (.NOT. l_geom_input .AND. .NOT.l_tokamak_input 
     1     .AND. (init_zetak_st == zero .AND. init_alpha_st == zero))
     2     tsymm = 0
       ENDIF

       IF (tsymm .eq. 0) THEN
         xmin = 0
         np0  = 6*np0_in/12+1                                              ! always odd + 3-multiple  number of points; (1/2) in symmetric CASE
         h0   = xmax/(np0-1)
         inc  = 1                                                          ! solution vector includes F(0)
       ELSE
         np0  = 6*(np0_in/6)+1                                             ! always odd + 3-multiple number of points
         xmin = -xmax
         h0   = 2*(xmax/(np0-1))
         inc  = 0                                                          ! solution vector does not include f(0) since f(0)=0.
       ENDIF

!================
!      END MIGRATION TO FULL MESH
!================
       IF (lscreen) THEN
          WRITE(*,48)
          IF (l_geom_input) THEN
            WRITE(*,49)
          ELSE
           IF (l_tokamak_input) THEN
             WRITE(*,47)
           ELSE
             WRITE(*,46)
           ENDIF
          ENDIF
          WRITE(*,48)
       END IF

!=============
!     BEGIN BALLOONING LOOP
!=============

       ballooning: DO kl=1,nlist

         ALLOCATE (hv(lmax+1), eigv(lmax), stat=k)
         IF (k .ne. 0) STOP 'Allocation 4 error in get_ballooning_grate'

         hv(1) = 1
         h = h0(list(kl))                                                 ! initialize step SIZE
         np = np0                                                         ! initialize number of points
         npm1 = np-1
         npm2 = np-2
         npm3 = np-3
         l_success =.false.

!-----------------------------
!      BEGIN RICHARDSON'S LOOP
!-----------------------------

         richardson: DO j=1,lMAX                                          ! lmax= MAX. #refinements in Richardson's scheme

           IF (j .gt. 1) THEN                                             ! reduce step ALWAYS by 2, so
             h = h/2                                                      ! new full mesh  =  previous half+full meshes
             np = np*2-1                                                  ! np =  total number of points in current iteration
             npm1 = np-1
             npm2 = np-2
             npm3 = np-3
           ENDIF

           ALLOCATE(a1f(npm2+inc), a2f(npm3+inc), a3f(npm3+inc),          !a1= diag; a2=sup-diag; a3=sub-diag
     1      xf(np), xh(npm1), pf(np), ph(npm1), qf(np), qh(npm1),         !xf=LINE full mesh; xh=LINE half mesh
     2      rf(np), rh(npm1), eigfun(np),feigfun(npm2+inc), stat=k)
           IF (k .ne. 0)
     1      STOP 'Allocation error 5 in get_ballooning_grate'

           IF(j.gt.1)then
              xf(1:np:2) = i_xf
              pf(1:np:2) = i_pf
              qf(1:np:2) = i_qf
              rf(1:np:2) = i_rf                                          ! IF refinement, USE previous....
              xf(2:npm1:2) = i_xh                                        ! ... evaluations to form new LINE full mesh
              pf(2:npm1:2) = i_ph
              qf(2:npm1:2) = i_qh
              rf(2:npm1:2) = i_rh
              DEALLOCATE(i_xf, i_pf, i_qf, i_rf, i_xh, i_ph,
     1          i_qh, i_rh, stat=k)                                      ! DEALLOCATE passing vectors
           ENDIF

           CALL getmatrix(list(kl), a1f, a2f, a3f, h, np, j, xf,
     1          xh, pf, ph, qf, qh, rf, rh, inc, xmin(list(kl)))
           IF (lfail_balloon) THEN
              grate = 100.0
              DEALLOCATE (lmnsf, bmncf, rmncpf, zmnspf,
     1                    lmnspf, bmncpf, iotapf, prespf, iotaf,
     2               presf, phipf, bsupvmncf, bsupumncf, stat=k)
              IF (lasym_v) THEN                                                 ! 110909 RS = Allow ASYMMETRIC input
                 DEALLOCATE (lmncf, bmnsf, rmnspf, zmncpf,
     1                lmncpf, bmnspf, bsupvmnsf, bsupumnsf, stat=k)
              ENDIF
              RETURN
           END IF
           CALL geteigm(a1f, a2f, a3f, npm2+inc, feigenv, feigfun)       ! get eigenvalue and eigenvector: 2-nd order
           IF (lfail_balloon) THEN
              grate = 100.0
              DEALLOCATE (lmnsf, bmncf, rmncpf, zmnspf,
     1                    lmnspf, bmncpf, iotapf, prespf, iotaf,
     2               presf, phipf, bsupvmncf, bsupumncf, stat=k)
              IF (lasym_v) THEN                                                 ! 110909 RS = Allow ASYMMETRIC input
                 DEALLOCATE (lmncf, bmnsf, rmnspf, zmncpf,
     1                lmncpf, bmnspf, bsupvmnsf, bsupumnsf, stat=k)
              ENDIF
              RETURN
           END IF

           IF (tsymm .eq. 0) THEN
             eigfun(1:np-1) = feigfun(1:np-1)
           ELSE
             eigfun(1) = 0
             eigfun(2:np-1) = feigfun(1:np-2)
           ENDIF
           eigfun(np) = 0
           CALL variat_eig_full(np, h, eigfun, pf, qf, rf, eigv(j))      !variational eigenvalue: 4th-order accurate

           DEALLOCATE (a1f, a2f, a3f, feigfun, stat=k)                   ! DEALLOCATE matrix memory space

           IF (j .ge. krich) THEN
              CALL polINT(hv(j-km), eigv(j-km), krich, zero, eigf,
     1          err)       ! carry out extrapolation
              rerr = ABS(err/eigf)
              IF(rerr .lt. tole) THEN
                l_success =.true.
                IF (eigf.lt.0) grate(list(kl)) = SQRT(-eigf)             ! IF conveged, EXIT richardson loop
                IF (eigf.ge.0) grate(list(kl)) =-SQRT(eigf)
                EXIT richardson
              ENDIF
           ENDIF

           hv(j+1) = 0.0625_dp*hv(j)                                     ! 2nd. order discretization =>  variat.eigenv. error is also 4th. order

           ALLOCATE (i_xf(np), i_xh(npm1), i_pf(np), i_ph(npm1),
     1        i_qf(np), i_qh(npm1), i_rf(np), i_rh(npm1), stat=k)        ! ALLOCATE passing vectors
           IF (k .ne. 0)
     1       STOP 'Allocation error 6 in get_ballooning_grate'

           i_xf = xf
           i_xh = xh
           i_pf = pf
           i_ph = ph                                                     ! store LINE mesh info in passing vectors
           i_qf = qf
           i_qh = qh
           i_rf = rf
           i_rh = rh

           DEALLOCATE(xf, xh, pf, ph, qf, qh, rf, rh, eigfun, stat=k)    ! DEALLOCATE old LINE mesh vectors

         END DO richardson

!--------------------------
!     END RICHARDSON'S LOOP
!--------------------------

         DEALLOCATE (hv, eigv, stat=k)                                   ! DEALLOCATE interpolation variables
         IF (ALLOCATED(i_xf))
     1      DEALLOCATE(i_xf,i_xh,i_pf,i_ph,i_qf,i_qh,i_rf,i_rh, stat=k)
         IF (ALLOCATED(xf))
     1      DEALLOCATE (xf, xh, pf, ph, qf, qh, rf, rh, eigfun, stat=k)

         IF (lscreen) THEN
           IF (l_geom_input) THEN
             WRITE(*, 50) list(kl), radii_flux(list(kl)),
     1       init_zeta, init_theta, grate(list(kl)), j, np, 
     2       xmax(list(kl)), l_success,
     3       tsymm, presf(list(kl)), mercierf(list(kl))
           ELSE
             IF (l_tokamak_input) THEN
               WRITE(*, 50) list(kl), radii_flux(list(kl)),
     1         init_thetak_tok, init_alpha_tok, grate(list(kl)), 
     2         j, np, xmax(list(kl)), l_success,
     3         tsymm, presf(list(kl)), mercierf(list(kl))
             ELSE
               WRITE(*, 50) list(kl), radii_flux(list(kl)),
     1         init_zetak_st, init_alpha_st, grate(list(kl)), 
     2         j, np, xmax(list(kl)), l_success,
     3         tsymm, presf(list(kl)), mercierf(list(kl))
             ENDIF
           ENDIF
         ENDIF


       END DO ballooning

!============
!     END BALLOONING LOOP
!============

       IF (lscreen) WRITE(*, 48)

48     FORMAT(112('-'))
49     FORMAT(3x,'NS',5x,'FLUX-s',5x,'ZT_0',8x,'TH_0',7x,
     1  'GR. RATE',4x,'IT',3x,'POINTS',6x,'XMAX',3x,'OK?',
     2  2x,'SYMM',4x,'PRES',8x,'MERC')
47     FORMAT(3x,'NS',5x,'FLUX-s',5x,'ZETAK',8x,'ALPHA_ST',7x,
     1  'GR. RATE',4x,'IT',3x,'POINTS',6x,'XMAX',3x,'OK?',
     2  2x,'SYMM',4x,'PRES',8x,'MERC')
46     FORMAT(3x,'NS',5x,'FLUX-s',5x,'THETAK',8x,'ALPHA_TOK',7x,
     1  'GR. RATE',4x,'IT',3x,'POINTS',6x,'XMAX',3x,'OK?',
     2  2x,'SYMM',4x,'PRES',8x,'MERC')
50     FORMAT(i5, 2x, 1pe10.2, 2(1pe10.2,2x), 1pe12.4, 3x, i2,
     1   2x, i6, 3x, 1pe10.2, 2x, l1, 3x, i2, 2(2x,1pe10.2))

!...  DEALLOCATE ALL variables

      DEALLOCATE (lmnsf, bmncf, rmncpf, zmnspf,
     1  lmnspf, bmncpf, iotapf, prespf, iotaf,
     2  presf, phipf, bsupvmncf, bsupumncf, stat=k)
      IF (lasym_v) THEN                                                 ! 110909 RS = Allow ASYMMETRIC input
        DEALLOCATE (lmncf, bmnsf, rmnspf, zmncpf,
     1    lmncpf, bmnspf, bsupvmnsf, bsupumnsf, stat=k)
      ENDIF

      END SUBROUTINE get_ballooning_grate
