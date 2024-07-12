      SUBROUTINE order_input(extension, nl, bsurf, ierr)
      USE stel_kinds
      USE read_wout_mod
      USE normalize_data
      USE readin_data
      USE fmesh_quantities, ONLY: rmncf, zmnsf, mercierf
     1  , rmnsf, zmncf                                                        ! 110909 RS: Asymmetric input
      USE vparams, ONLY : mu0
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER*(*), INTENT(IN) :: extension
      INTEGER, INTENT(IN) :: nl
      INTEGER, INTENT(OUT) :: ierr
      INTEGER, DIMENSION(nl), INTENT(IN):: bsurf
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: lk, i, j, k
      INTEGER, DIMENSION(:), ALLOCATABLE :: plist
      REAL(rprec) :: aspect_v, r0max_v, r0min_v, betaxis_v
      LOGICAL, ALLOCATABLE, DIMENSION(:):: lballoon
!-----------------------------------------------
!
!      Variables allocated in this routine (deallocated elsewhere)
!         hiota, hpres, hphip, xm_v, xn_v, rmncf, zmnsf, mercierf,
!         lmnsh, bmnch, list, bsupumnch, bsupvmnch
!
!      110909 RS: IF LASYM_V = TRUE, also allocate arrays for asymmetric input
!         lmnch, bmnsh, bsupumnsh, bsupvmnsh, rmnsf, zmncf
!
!      All half mesh quantities have first surface (js =1) equal to zero
!

      ! SAL - PPPL Modification to avoid bad allocations
       IF (.not. lwout_opened) THEN
          call read_wout_file('wout_'//trim(extension), ierr)
          if (ierr .ne. 0) then
             ierr = 0
             call read_wout_file('wout.'//trim(extension), ierr)
             if (ierr .ne. 0) then
                print *,' ierr = ', ierr,
     1          'error in read_wout_file called from COBRA'
             end if
          end if
       END IF
!      IF (.not. lwout_opened) THEN
!         CALL read_wout_file(extension, ierr)                                    !...   Read surface quantities and normalization data
!         IF (ierr .ne. 0) THEN
!            PRINT *,' ierr = ', ierr,
!     1              'error in read_wout_file called from COBRA'
!         END IF
!      END IF
      IF (ierr_vmec.ne.0 .or. ierr.ne.0) GOTO 1000

      nfp_v = nfp
      ns_cob = ns
      aspect_v = aspect
      r0max_v = rmax_surf
      r0min_v = rmin_surf
      betaxis_v = betaxis
      mnmax_v = mnmax
      mnmax_vnyq = mnmax_nyq
      lasym_v = lasym                                                         ! 111009: RS, Logical for Asymmetric input

!...   IOTA, PRES, PHIP are ALL on HALF-MESH!!!

       k=0
       IF (.not. ALLOCATED(hiota)) ALLOCATE(hiota(ns_cob),stat=k)
       if (k .ne. 0) stop 'Allocation error 1 in cobra order_input'
       IF (.not. ALLOCATED(hpres)) ALLOCATE(hpres(ns_cob),stat=k)
       if (k .ne. 0) stop 'Allocation error 1 in cobra order_input'
       IF (.not. ALLOCATED(hphip)) ALLOCATE(hphip(ns_cob),stat=k)
       if (k .ne. 0) stop 'Allocation error 1 in cobra order_input'
       IF (.not. ALLOCATED(mercierf)) ALLOCATE(mercierf(ns_cob),stat=k)
       if (k .ne. 0) stop 'Allocation error 1 in cobra order_input'
!       allocate (hiota(ns_cob), hpres(ns_cob), hphip(ns_cob),
!     1   mercierf(ns_cob), stat=k)
!       if (k .ne. 0) stop 'Allocation error 1 in cobra order_input'

      hiota = iotas(1:ns_cob)
      IF (version_ .gt. 6.0) THEN
          hpres = mu0*pres(1:ns_cob)                                          ! in VMEC versions > 6.00 pressure is given in pascals
      ELSE
         hpres = pres(1:ns_cob)                                               ! in VMEC versions  .LE. 6.00 pressure is p= mu_o * p(pascals)
      ENDIF
      hphip = phip(1:ns_cob)

      mercierf = Dmerc(1:ns_cob)

      r0 = (r0max_v+r0min_v)/2
      amin = r0/aspect_v
      beta0 = betaxis_v
      IF (beta0 .le. zero) THEN
          beta0 = EPSILON(beta0)
          hpres = beta0*( /(1 - (j - 1)/REAL(ns - 1,rprec), j=1,ns)/ )
      END IF
      b0_v = SQRT((2._dp/beta0)*(1.5_dp*hpres(2)-.5_dp*hpres(3)))

       k=0
       IF (.not. ALLOCATED(xm_v)) ALLOCATE(xm_v(mnmax_v),stat=k)
       if (k .ne. 0) stop 'Allocation error 2 in cobra order_input'
       IF (.not. ALLOCATED(xn_v)) ALLOCATE(xn_v(mnmax_v),stat=k)
       if (k .ne. 0) stop 'Allocation error 2 in cobra order_input'
       IF (.not. ALLOCATED(xm_vnyq)) 
     1           ALLOCATE(xm_vnyq(mnmax_vnyq),stat=k)
       if (k .ne. 0) stop 'Allocation error 2 in cobra order_input'
       IF (.not. ALLOCATED(xn_vnyq)) 
     1           ALLOCATE(xn_vnyq(mnmax_vnyq),stat=k)
       if (k .ne. 0) stop 'Allocation error 2 in cobra order_input'
       IF (.not. ALLOCATED(plist)) ALLOCATE(plist(ns_cob),stat=k)
       if (k .ne. 0) stop 'Allocation error 2 in cobra order_input'
!       allocate (xm_v(mnmax_v), xn_v(mnmax_v), plist(ns_cob), stat=k)
!       if (k .ne. 0) stop 'Allocation error 2 in cobra order_input'
      xm_v = xm(1:mnmax_v)
      xn_v = xn(1:mnmax_v)
      xm_vnyq = xm_nyq(1:mnmax_vnyq)
      xn_vnyq = xn_nyq(1:mnmax_vnyq)

      mndim = ns_cob*mnmax_v
      mndim_nyq = ns_cob*mnmax_vnyq

!...   RMN, ZMN are on FULL-MESH; LMNS, BSUPUMN (c), BSUPVMN (c),
!      and BMNC are on HALF-MESH!!!

       k=0
       IF (.not. ALLOCATED(rmncf)) ALLOCATE(rmncf(mndim),stat=k)
       if (k .ne. 0) stop 'Allocation error 3 in cobra order_input'
       IF (.not. ALLOCATED(zmnsf)) ALLOCATE(zmnsf(mndim),stat=k)
       if (k .ne. 0) stop 'Allocation error 3 in cobra order_input'
       IF (.not. ALLOCATED(lmnsh)) ALLOCATE(lmnsh(mndim),stat=k)
       if (k .ne. 0) stop 'Allocation error 3 in cobra order_input'
       IF (.not. ALLOCATED(bmnch)) ALLOCATE(bmnch(mndim_nyq),stat=k)
       if (k .ne. 0) stop 'Allocation error 3 in cobra order_input'
       IF (.not. ALLOCATED(bsupumnch)) 
     1       ALLOCATE(bsupumnch(mndim_nyq),stat=k)
       if (k .ne. 0) stop 'Allocation error 3 in cobra order_input'
       IF (.not. ALLOCATED(bsupvmnch)) 
     1       ALLOCATE(bsupvmnch(mndim_nyq),stat=k)
       if (k .ne. 0) stop 'Allocation error 3 in cobra order_input'
!       allocate (rmncf(mndim), zmnsf(mndim), lmnsh(mndim),
!     1    bmnch(mndim),  bsupumnch(mndim), bsupvmnch(mndim),
!     2    stat=k)
!       if (k .ne. 0) stop 'Allocation error 3 in cobra order_input'

      zmnsf=0; rmncf=0; lmnsh=0; bmnch=0; bsupvmnch=0; bsupumnch=0

      IF (lasym_v) THEN                                                      !  110909 RS: For asymmetric input
         k=0
         IF (.not. ALLOCATED(rmnsf)) ALLOCATE(rmnsf(mndim),stat=k)
         if (k .ne. 0) stop 'Allocation error 3 in cobra order_input'
         IF (.not. ALLOCATED(zmncf)) ALLOCATE(zmncf(mndim),stat=k)
         if (k .ne. 0) stop 'Allocation error 3 in cobra order_input'
         IF (.not. ALLOCATED(lmnch)) ALLOCATE(lmnch(mndim),stat=k)
         if (k .ne. 0) stop 'Allocation error 3 in cobra order_input'
         IF (.not. ALLOCATED(bmnsh)) ALLOCATE(bmnsh(mndim_nyq),stat=k)
         if (k .ne. 0) stop 'Allocation error 3 in cobra order_input'
         IF (.not. ALLOCATED(bsupumnsh)) ALLOCATE(bsupumnsh(mndim_nyq),
     1       stat=k)
         if (k .ne. 0) stop 'Allocation error 3 in cobra order_input'
         IF (.not. ALLOCATED(bsupvmnsh)) ALLOCATE(bsupvmnsh(mndim_nyq),
     1       stat=k)
         if (k .ne. 0) stop 'Allocation error 3 in cobra order_input'
!         ALLOCATE (rmnsf(mndim), zmncf(mndim), lmnch(mndim),
!     1    bmnsh(mndim),  bsupumnsh(mndim), bsupvmnsh(mndim), stat=k)
!         IF (k .ne. 0) STOP 'Allocation error 3 in cobra order_input'
         zmncf=0; rmnsf=0; lmnch=0; bmnsh=0; bsupvmnsh=0; bsupumnsh=0
      ENDIF

      DO k=1, ns_cob
         lk = (k-1)*mnmax_v
         rmncf(lk+1:lk+mnmax_v) = rmnc(1:mnmax_v,k)
         zmnsf(lk+1:lk+mnmax_v) = zmns(1:mnmax_v,k)
         lmnsh(lk+1:lk+mnmax_v) = lmns(1:mnmax_v,k)
         IF (lasym_v) THEN                                                   !  110909 RS: For asymmetric input
           rmnsf(lk+1:lk+mnmax_v) = rmns(1:mnmax_v,k)
           zmncf(lk+1:lk+mnmax_v) = zmnc(1:mnmax_v,k)
           lmnch(lk+1:lk+mnmax_v) = lmnc(1:mnmax_v,k)
         ENDIF

      ENDDO

      DO k=1, ns_cob
         lk = (k-1)*mnmax_vnyq
         bmnch(lk+1:lk+mnmax_vnyq) = bmnc(1:mnmax_vnyq,k)
         bsupvmnch(lk+1:lk+mnmax_vnyq) = bsupvmnc(1:mnmax_vnyq,k)
         bsupumnch(lk+1:lk+mnmax_vnyq) = bsupumnc(1:mnmax_vnyq,k)
         IF (lasym_v) THEN                                                   !  110909 RS: For asymmetric input
           bmnsh(lk+1:lk+mnmax_vnyq) = bmns(1:mnmax_vnyq,k)
           bsupvmnsh(lk+1:lk+mnmax_vnyq) = bsupvmns(1:mnmax_vnyq,k)
           bsupumnsh(lk+1:lk+mnmax_vnyq) = bsupumns(1:mnmax_vnyq,k)
         ENDIF

      ENDDO

!...   identify wanted surfaces

      ALLOCATE (lballoon(ns_cob), stat=k)
      lballoon=.false.
      DO i=1,nl
         lballoon(bsurf(i))=.true.
      ENDDO

!...   check for feasibility of ballooning calculation

      nlist=0
      DO i = 2, ns_cob - 1                                                  ! exclude Boundary (i=ns)
         IF (lballoon(i)) THEN
            nlist=nlist+1
            plist(nlist)=i                                                  ! i-th surface on full is between i-th and (i+1)-th on half
         ENDIF
      ENDDO

!...   store NLIST surfaces WHERE calculation can be done in LIST(1:NLIST) vector

      IF (.not.ALLOCATED(list)) allocate (list(nlist), stat=k)
      list(1:nlist) = plist(1:nlist)
      DEALLOCATE (plist,lballoon)
!1000  CALL read_wout_deallocate
1000  IF (ierr .eq. 0) THEN
         ierr = -ierr_vmec
      ELSE
         PRINT *,' ierr = ', ierr, ' writing in READ_WOUT_BOOZ'
      END IF

      END SUBROUTINE order_input
