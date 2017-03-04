      MODULE diagnostics_mod
!      USE stel_kinds
#if defined(SKS)
      USE fourier, ONLY: tomnsp, toijsp, tomnsp_par, toijsp_par, tomnsp_pest
#else
      USE fourier, ONLY: tomnsp, toijsp, tomnsp_pest
#endif
      USE quantities
      USE descriptor_mod, ONLY: iam
      USE timer_mod
      IMPLICIT NONE
      INTEGER     :: unit35=35, unit36=36
      REAL(rprec) :: divb_rms, avbeta, bdotj_rms, divj_rms, bdotj2_rms
      REAL(rprec) :: toroidal_flux, toroidal_flux0=0, bgradp_rms, wbgradp
      REAL(rprec) :: max_bgradp, min_bgradp
      REAL(rprec) :: tnorm_bgradp
!
!     The next three allocatable variables are kept in memory and used to
!     produce output whenever L_DUMP_DIAGNO = .TRUE. 
!      
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: bdotjmnch
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: divjmnsh, divbmnsf
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: bgradpf, jpmnch_loc
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:), TARGET :: fgradpmnch
      REAL(rprec), DIMENSION(:), POINTER     :: xc, gc
      LOGICAL :: lextend = .FALSE.                   !Set in call to EVOLVE_BGRADP
      LOGICAL :: lcurr_init = .FALSE.
!
!     CALL THIS AFTER THE B-FIELD COMPONENTS bi = Jac*B-sup-i HAVE BEEN UPDATED
!     IN UPDATE_STATE (CALL FROM UPDATE_STATE IS SAFEST WAY!) 

      CONTAINS

          SUBROUTINE divb
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER                  :: m, n, istat
        REAL(rprec)              :: ton, toff, tnorm
!-----------------------------------------------
!
!     COMPUTE divB ON FULL MESH POINTS
!
        CALL second0(ton)

        IF (.NOT. ALLOCATED(divbmnsf))  THEN                                  
          ALLOCATE(divbmnsf(0:mpol,-ntor:ntor,ns),  stat=istat)          ! Otherwise, allocate first time called
          IF (istat .NE. 0) STOP 'Allocation failed in divb'
          divbmnsf = 0
        ENDIF

        DO m = 0, mpol
          DO n = -ntor,ntor
            divbmnsf(m,n,2:nsh) =                                       &
                   -m*(jbsupumnch(m,n,3:ns)+jbsupumnch(m,n,2:nsh))/2    &
               -n*nfp*(jbsupvmnch(m,n,3:ns)+jbsupvmnch(m,n,2:nsh))/2    &
               + ohs *(jbsupsmnsh(m,n,3:ns)-jbsupsmnsh(m,n,2:nsh))
          END DO
        END DO

        tnorm = hs_i*SUM(jbsupumnch**2 + jbsupvmnch**2 + jbsupsmnsh**2)
        divbmnsf = divbmnsf/SQRT(tnorm)
        divb_rms = SQRT(SUM(divbmnsf**2))

!        IF (iam .EQ. 0) WRITE(4000,*)'DIVB: tnorm=',SQRT(tnorm),' divb_rms: ',divb_rms
        CALL second0(toff)
        time_divb = time_divb + (toff-ton)

        END SUBROUTINE divb

        SUBROUTINE divb_par
        USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR
#if defined(SKS)
        INCLUDE 'mpif.h'
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER  :: istat, m, n, nsmin, nsmax, n1, n2
        REAL(dp) :: tnorm, s1(2), r1(2), ton, toff
!-----------------------------------------------
!
!     COMPUTE divB ON FULL MESH POINTS
!
        CALL second0(ton)

        nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow,ns)
        IF (.NOT. ALLOCATED(divbmnsf)) THEN
        ALLOCATE(divbmnsf(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)   ! Otherwise, allocate first time called
        IF (istat .NE. 0) STOP 'Allocation failed in divb'
        divbmnsf = 0
        END IF

        n1 = MAX(2, nsmin)
        n2 = MIN(ns-1, nsmax)
        DO m = 0, mpol
           DO n = -ntor, ntor
           divbmnsf(m,n,n1:n2) =                                        &
                -m*(jbsupumnch(m,n,n1+1:n2+1)+jbsupumnch(m,n,n1:n2))/2  &
            -n*nfp*(jbsupvmnch(m,n,n1+1:n2+1)+jbsupvmnch(m,n,n1:n2))/2  &
            + ohs *(jbsupsmnsh(m,n,n1+1:n2+1)-jbsupsmnsh(m,n,n1:n2))
           END DO
        END DO

        tnorm = hs_i*SUM(jbsupumnch(:,:,nsmin:nsmax)**2                 &
              +     jbsupvmnch(:,:,nsmin:nsmax)**2                      &
              +     jbsupsmnsh(:,:,nsmin:nsmax)**2)
       
        divb_rms = SUM(divbmnsf**2)
        s1(1) = divb_rms; s1(2) = tnorm
        CALL MPI_ALLREDUCE(s1,r1,2,MPI_REAL8, MPI_SUM,                  &
                           MPI_COMM_WORLD,MPI_ERR)
        divb_rms=r1(1);   tnorm = r1(2)
        IF (tnorm .NE. zero) THEN
           divb_rms = SQRT(divb_rms/tnorm)
           divbmnsf = divbmnsf/SQRT(tnorm)
        END IF

        CALL second0(toff)
        time_divb = time_divb + (toff-ton)
#endif
        END SUBROUTINE divb_par

         
        SUBROUTINE WRITE_PROFILES(fsq_total)
        USE safe_open_mod
        USE metrics, ONLY: lmns_i
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
        REAL(rprec) :: fsq_total
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER                  :: istat, m, js
        REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: pijch, pmnch
        REAL(rprec)              :: tnorm
!-----------------------------------------------
        ALLOCATE(pijch(ntheta,nzeta,ns), stat=istat)
        ALLOCATE(pmnch(0:mpol,-ntor:ntor,ns), stat=istat)      
        IF (istat .ne. 0) STOP 'Allocation failed in WRITE_PROFILES' 

        CALL safe_open(unit35, istat, "siesta_profiles.txt", 'replace', &
                       'formatted')
        CALL safe_open(unit36, istat, "siesta_profiles_pest.txt",       &
                       'replace','formatted')

        IF (istat .NE. 0) STOP 'Error writing siesta profiles data file'

        pijch = 0
        CALL toijsp(jpmnch, pijch, 0, 0, 0, 1)                           
        pijch = pijch/jacobh                                            ! Need first to remove jacobian...
        pmnch = 0
        CALL tomnsp(pijch, pmnch, 0, 1) 
        CALL siesta_profiles (unit35,pmnch, fsq_total, 'pmnch(jrad);  fsq_total: ')
        CALL tomnsp_pest(pijch, pmnch, lmns_i, 0, 1)
        CALL siesta_profiles (unit36, pmnch, fsq_total, 'pmnch(jrad);  fsq_total: ')


!USE pijch, pmnch to for components of bsupX on half-mesh
!BSUPS
        CALL toijsp(jbsupsmnsh, pijch, 0, 0, 1, 1)
        pijch = pijch/jacobh
        CALL tomnsp(pijch, pmnch, 1, 1)
        CALL siesta_profiles (unit35, pmnch, fsq_total, 'B^s_mn;  fsq_total: ')

!BSUPU
        CALL toijsp(jbsupumnch, pijch, 0, 0, 0, 1)
        pijch = pijch/jacobh
        CALL tomnsp(pijch, pmnch, 0, 1)
        CALL siesta_profiles (unit35, pmnch, fsq_total, 'B^u_mn;  fsq_total: ')

!BSUPV
        CALL toijsp(jbsupvmnch, pijch, 0, 0, 0, 1)
        pijch = pijch/jacobh
        CALL tomnsp(pijch, pmnch, 0, 1)
        CALL siesta_profiles (unit35, pmnch, fsq_total, 'B^v_mn;  fsq_total: ')

!KSUPS CRCook
        CALL toijsp(ksupsmnsf, pijch, 0, 0, 1, 1)
        pijch = pijch/jacobh
        CALL tomnsp(pijch, pmnch, 1, 1)
        CALL siesta_profiles (unit35, pmnch, fsq_total, 'K^s_mn;  fsq_total: ')

!KSUPU CRCook
        CALL toijsp(ksupumncf, pijch, 0, 0, 0, 1)
        pijch = pijch/jacobh
        CALL tomnsp(pijch, pmnch, 0, 1)
        CALL siesta_profiles (unit35, pmnch, fsq_total, 'K^u_mn;  fsq_total: ')

!KSUPV CRCook
        CALL toijsp(ksupvmncf, pijch, 0, 0, 0, 1)
        pijch = pijch/jacobh
        CALL tomnsp(pijch, pmnch, 0, 1)
        CALL siesta_profiles (unit35, pmnch, fsq_total, 'K^v_mn;  fsq_total: ')

!Displacements X jac (SPH:052114)
        CALL siesta_profiles (unit35, jvsupsmncf, fsq_total, 'vel^s_mn(jrad);  fsq_total: ')
        CALL siesta_profiles (unit35, jvsupumnsf, fsq_total, 'vel^u_mn(jrad);  fsq_total: ')
        CALL siesta_profiles (unit35, jvsupvmnsf, fsq_total, 'vel^v_mn(jrad);  fsq_total: ')

        DEALLOCATE (pijch, pmnch)

        CLOSE (35)

      END SUBROUTINE WRITE_PROFILES

      SUBROUTINE SIESTA_PROFILES (iunit, value, fsq_total, label)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in)       :: iunit
      REAL(rprec), INTENT(in)   :: value(0:mpol,-ntor:ntor,ns)
      REAL(rprec), INTENT(in)   :: fsq_total
      CHARACTER*(*), INTENT(in) :: label
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: m, js 
!-----------------------------------------------

      WRITE (iunit, *)
      WRITE (iunit, *) TRIM(label), fsq_total
      DO m=0,mpol,6
         IF (m+5 .gt. mpol) EXIT
         WRITE (iunit, '("MPOL",6i10)') m, m+1, m+2, m+3, m+4, m+5
         DO js=2,ns
            WRITE (iunit,'(i4,1p6e10.2)') js, SUM(value(m,:,js)),       &
            SUM(value(m+1,:,js)), SUM(value(m+2,:,js)),                 &
            SUM(value(m+3,:,js)), SUM(value(m+4,:,js)),                 &
            SUM(value(m+5,:,js))
         END DO
      END DO

      END SUBROUTINE SIESTA_PROFILES


      SUBROUTINE bgradp
!     
!     WRITTEN 08-16-07 BY R. SANCHEZ AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: Computes B*GRAD P (normalized to B*P) on the FULL mesh
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER                  :: istat
!-----------------------------------------------
!
!     COMPUTE B DOT GRAD P AT FULL-GRID POINTS (sans endpts)
!
      ALLOCATE(bgradpf(ntheta,nzeta,ns), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error in BGRADP'

      CALL get_bgradp

      DEALLOCATE (bgradpf)

      END SUBROUTINE bgradp


      SUBROUTINE get_bgradp
!     
!     PURPOSE: Computes B*GRAD P on the FULL mesh
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER     :: istat
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)  ::                       &
          dpduijsh, dpdvijsh, pijch,                                    &
          bsupsijsh0, bsupuijch0, bsupvijch0
      REAL(dp)    :: ton, toff
!-----------------------------------------------
!
!     COMPUTE B DOT GRAD P AT FULL-GRID POINTS (sans endpts)
!
      CALL second0(ton)

      bgradpf = 0 

!     Compute pressure and angle derivatives
!
      ALLOCATE(pijch(ntheta,nzeta,ns),                                  &
               dpduijsh(ntheta,nzeta,ns),                               &
               dpdvijsh(ntheta,nzeta,ns),                               &
               bsupsijsh0(ntheta,nzeta,ns),                             &
               bsupuijch0(ntheta,nzeta,ns),                             &
               bsupvijch0(ntheta,nzeta,ns), stat=istat)                         
      IF (istat .ne. 0) STOP 'Allocation failed in GET_BGRADP'

      CALL GetPressure (pijch, dpduijsh, dpdvijsh)
! 
!     Compute contravariant components of magnetic field
!           
      IF (.NOT.ALLOCATED(jbsupsmnsh)) STOP 'jbsupXmn UNALLOCATED IN GET_BGRADP'
      CALL toijsp(jbsupsmnsh, bsupsijsh0, 0, 0, 1, 1)                           
      bsupsijsh0 = bsupsijsh0/jacobh
      CALL toijsp(jbsupumnch, bsupuijch0, 0, 0, 0, 1)       
      bsupuijch0 = bsupuijch0/jacobh
      CALL toijsp(jbsupvmnch, bsupvijch0, 0, 0, 0, 1) 
      bsupvijch0 = bsupvijch0/jacobh

      bgradpf(:,:,2:nsh) =                                              &   ! Get all values except first and
               (bsupsijsh0(:,:,3:ns) + bsupsijsh0(:,:,2:nsh))           &   ! last surface from half mesh values.
               *    (pijch(:,:,3:ns) -      pijch(:,:,2:nsh))*ohs/2     &
             + (bsupuijch0(:,:,3:ns) + bsupuijch0(:,:,2:nsh))           &
              *(  dpduijsh(:,:,3:ns) +   dpduijsh(:,:,2:nsh))/4         &
             + (bsupvijch0(:,:,3:ns) + bsupvijch0(:,:,2:nsh))           &
              *(  dpdvijsh(:,:,3:ns) +   dpdvijsh(:,:,2:nsh))/4                   
    
      tnorm_bgradp = SUM((bsupsijsh0**2 + bsupuijch0**2 + bsupvijch0**2)*pijch**2)
      IF (tnorm_bgradp .eq. zero) tnorm_bgradp = EPSILON(tnorm_bgradp)

      wbgradp = SUM(bgradpf(:,:,3:nsh)**2)
      bgradp_rms = SQRT(wbgradp/tnorm_bgradp)

      tnorm_bgradp = SQRT(tnorm_bgradp/(ns*ntheta*nzeta))

!      IF (iam .EQ. 0) WRITE (4200, *)' tnorm_bgradp: ',tnorm_bgradp,' bgradp_rms: ', bgradp_rms

      bgradpf(:,:,2:nsh) = bgradpf(:,:,2:nsh)/tnorm_bgradp
!IGNORE FIRST POINT
      max_bgradp = MAXVAL(bgradpf(:,:,3:nsh))
      min_bgradp = MINVAL(bgradpf(:,:,3:nsh))

      DEALLOCATE (pijch, dpdvijsh, dpduijsh,                            &
                  bsupsijsh0, bsupuijch0, bsupvijch0)

      CALL second0(toff)
      time_bgradp = time_bgradp+(toff-ton)

      END SUBROUTINE get_bgradp

      SUBROUTINE bgradp_par
!     
!     WRITTEN 03-21-13 BY S. HIRSHMAN AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: Computes B*GRAD P (normalized to B*P) on the FULL mesh
      USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR
#if defined(SKS)
      INCLUDE 'mpif.h'
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER     :: istat, nsmin, nsmax
      REAL(dp)    :: s1(2), r1(2), ton, toff
!-----------------------------------------------
!
!     COMPUTE B DOT GRAD P AT FULL-GRID POINTS (sans endpts)
!
      nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow,ns)
      
      ALLOCATE(bgradpf(ntheta,nzeta,nsmin:nsmax), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error in BGRADP'
      CALL get_bgradp_par

      DEALLOCATE (bgradpf, stat=istat)
#endif
      END SUBROUTINE bgradp_par


      SUBROUTINE get_bgradp_par
#if defined(SKS)
      USE nscalingtools, ONLY: startglobrow, endglobrow, mpi_err
      INCLUDE 'mpif.h'
!     
!     PURPOSE: Computes B*GRAD P on the FULL mesh
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER  :: jsh=3
      INTEGER     :: istat, nsmin, nsmax, ns_span, n1, n2
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)  ::                       &
          dpduijsh, dpdvijsh, pijch,                                    &
          bsupsijsh0, bsupuijch0, bsupvijch0
      REAL(dp)    :: s1(2), r1(2), ton, toff
!-----------------------------------------------
!
!     COMPUTE B DOT GRAD P AT FULL-GRID POINTS (sans endpts)
!
      CALL second0(ton)

      nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow+1,ns)
      ns_span = nsmax-nsmin+1

      bgradpf = 0 

!     Compute pressure and angle derivatives
!
      ALLOCATE(pijch(ntheta,nzeta,nsmin:nsmax),                         &
               dpduijsh(ntheta,nzeta,nsmin:nsmax),                      &
               dpdvijsh(ntheta,nzeta,nsmin:nsmax),                      &
               bsupsijsh0(ntheta,nzeta,nsmin:nsmax),                    &
               bsupuijch0(ntheta,nzeta,nsmin:nsmax),                    &
               bsupvijch0(ntheta,nzeta,nsmin:nsmax), stat=istat)                         
      IF (istat .NE. 0) STOP 'Allocation failed in GET_BGRADP_PAR'

      IF (.NOT.ALLOCATED(jbsupsmnsh_p)) STOP 'jbsupXmn_p UNALLOCATED IN GET_BGRADP_PAR'

      CALL GetPressure_par (pijch, dpduijsh, dpdvijsh, nsmin, nsmax)
! 
!     Compute contravariant components of magnetic field
!           
      CALL toijsp_par(jbsupsmnsh_p(:,:,nsmin:nsmax), bsupsijsh0, 0, 0, 1, 1, 1, ns_span)                           
      bsupsijsh0 = bsupsijsh0/jacobh_p(:,:,nsmin:nsmax)
      CALL toijsp_par(jbsupumnch_p(:,:,nsmin:nsmax), bsupuijch0, 0, 0, 0, 1, 1, ns_span)       
      bsupuijch0 = bsupuijch0/jacobh_p(:,:,nsmin:nsmax)
      CALL toijsp_par(jbsupvmnch_p(:,:,nsmin:nsmax), bsupvijch0, 0, 0, 0, 1, 1, ns_span) 
      bsupvijch0 = bsupvijch0/jacobh_p(:,:,nsmin:nsmax)

      nsmax=MIN(endglobrow,ns)
      n1 = MAX(2, nsmin); n2 = MIN(nsh, nsmax)
      bgradpf(:,:,n1:n2) =                                              &   ! Get all values except first and
               (bsupsijsh0(:,:,n1+1:n2+1) + bsupsijsh0(:,:,n1:n2))      &   ! last surface from half mesh values.
             *ohs/2*(pijch(:,:,n1+1:n2+1) -      pijch(:,:,n1:n2))      &
             + (bsupuijch0(:,:,n1+1:n2+1) + bsupuijch0(:,:,n1:n2))      &
              *(  dpduijsh(:,:,n1+1:n2+1) +   dpduijsh(:,:,n1:n2))/4    &
             + (bsupvijch0(:,:,n1+1:n2+1) + bsupvijch0(:,:,n1:n2))      &
              *(  dpdvijsh(:,:,n1+1:n2+1) +   dpdvijsh(:,:,n1:n2))/4                   
    
      tnorm_bgradp = SUM((bsupsijsh0(:,:,n1:nsmax)**2                   &
                   +      bsupuijch0(:,:,n1:nsmax)**2                   &
                   +      bsupvijch0(:,:,n1:nsmax)**2)                  &
                   *           pijch(:,:,n1:nsmax)**2)

      n1 = MAX(3, nsmin)
      wbgradp = SUM(bgradpf(:,:, n1:n2)**2)
      s1(1)=tnorm_bgradp;  s1(2)=wbgradp
      CALL MPI_ALLREDUCE(s1,r1,2,MPI_REAL8, MPI_SUM,                    &
                         MPI_COMM_WORLD,MPI_ERR)
      tnorm_bgradp=r1(1); wbgradp=r1(2)

      IF (tnorm_bgradp .GT. zero) THEN
         bgradp_rms = SQRT(wbgradp/tnorm_bgradp)
         tnorm_bgradp = SQRT(tnorm_bgradp/(ns*nuv))

!         IF (iam .EQ. 0) WRITE (3200, *)' tnorm_bgradp: ',tnorm_bgradp, &
!                         ' bgradp_rms: ', bgradp_rms

!IGNORE FIRST POINT
         max_bgradp = MAXVAL(bgradpf(:,:,n1:n2))/tnorm_bgradp
         min_bgradp = MINVAL(bgradpf(:,:,n1:n2))/tnorm_bgradp
         s1(1)=max_bgradp
         CALL MPI_ALLREDUCE(s1,r1,1,MPI_REAL8, MPI_MAX,                 &
                            MPI_COMM_WORLD,MPI_ERR)
         max_bgradp=r1(1)
         s1(1)=min_bgradp
         CALL MPI_ALLREDUCE(s1,r1,1,MPI_REAL8, MPI_MIN,                 &
                            MPI_COMM_WORLD,MPI_ERR)
         min_bgradp=r1(1)
      END IF

      DEALLOCATE (pijch, dpdvijsh, dpduijsh,                            &
                  bsupsijsh0, bsupuijch0, bsupvijch0, stat=istat)

      CALL second0(toff)
      time_bgradp = time_bgradp+(toff-ton)
#endif
      END SUBROUTINE get_bgradp_par
      
      SUBROUTINE GetPressure (pijch, dpduijsh, dpdvijsh)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(ntheta,nzeta,ns), INTENT(OUT)  ::          &
           dpduijsh, dpdvijsh, pijch
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  :: istat
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: pmnch
!-----------------------------------------------
!      mpol_s = mpol;  mpol = mpol32
!      ntor_s = ntor;  ntor = ntor32

      IF (.NOT. ALLOCATED(jpmnch_loc)) THEN
         ALLOCATE(jpmnch_loc(0:mpol,-ntor:ntor,ns), stat=istat)
         IF (istat .ne. 0) STOP 'Allocation error in GetPressure!'
         jpmnch_loc = 0
      END IF

      jpmnch_loc = jpmnch

      CALL toijsp(jpmnch_loc, pijch, 0, 0, 0, 1)                           
      pijch = pijch/jacobh                                              ! Remove jacobian factor
      
      ALLOCATE(pmnch(0:mpol,-ntor:ntor,ns), stat=istat)      
      CALL tomnsp(pijch, pmnch, 0, 1) 
      dpduijsh = 0; dpdvijsh = 0
      CALL toijsp(pmnch, pijch, 0, 0, 0, 1)
      CALL toijsp(pmnch, dpduijsh, 1, 0, 0, 1)                          ! Get theta and phi pressure derivatives            
      CALL toijsp(pmnch, dpdvijsh, 0, 1, 0, 1)                           
     
      DEALLOCATE (pmnch)

      END SUBROUTINE GetPressure

#if defined(SKS)
      SUBROUTINE GetPressure_par (pijch, dpduijsh, dpdvijsh, nsmin, nsmax)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)   :: nsmin, nsmax
      REAL(dp), DIMENSION(ntheta,nzeta, nsmin:nsmax),                   &
               INTENT(OUT)  :: dpduijsh, dpdvijsh, pijch
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  :: istat, ns_span
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: pmnch
!-----------------------------------------------
      ns_span = nsmax-nsmin+1

      CALL toijsp_par(jpmnch_p(:,:,nsmin:nsmax), pijch, 0, 0, 0, 1, 1, ns_span)                           
      pijch = pijch/jacobh_p(:,:,nsmin:nsmax)                           ! Remove jacobian factor
      
      ALLOCATE(pmnch(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)      
      IF (istat .NE. 0) STOP 'Allocation error in GetPressure_par'

      CALL tomnsp_par(pijch, pmnch, 0, 1, 1, ns_span) 

      CALL toijsp_par(pmnch, pijch, 0, 0, 0, 1, 1, ns_span)
      CALL toijsp_par(pmnch, dpduijsh, 1, 0, 0, 1, 1, ns_span)          ! Get theta and phi pressure derivatives
      CALL toijsp_par(pmnch, dpdvijsh, 0, 1, 0, 1, 1, ns_span)                           
     
      DEALLOCATE (pmnch)

      END SUBROUTINE GetPressure_par
#endif

      SUBROUTINE BETA

        avbeta = wp/wb
 
      END SUBROUTINE BETA
 
      
      SUBROUTINE TFLUX
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: m0=0, n0=0
!-----------------------------------------------

!Averages over all toroidal cross sections (which should be the same)
      toroidal_flux = SUM(jbsupvmnch(m0,n0,2:ns))
      
      toroidal_flux = signjac*twopi*toroidal_flux*hs_i/b_factor  
      
      IF (toroidal_flux0 .EQ. zero) toroidal_flux0 = toroidal_flux

      END SUBROUTINE TFLUX


      SUBROUTINE bdotj
!
!     WRITTEN 08-20-07 BY R. SANCHEZ AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: Computes parallel current normalized to the total current (J*B/|JxB|) on the half mesh
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER      :: istat
        REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: bdotjijch
        REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: bsupsijsh,        &
          bsupuijch, bsupvijch,                                         &
          bsubsijsh, bsubuijch, bsubvijch, ksubsijsh, ksubuijch,        &
          ksubvijch, ksupsijsh, ksupuijch, ksupvijch, work1, work2,     &
          work3, work4, work5, work6
        REAL(rprec)  :: tnorm, tnorm2, ton, toff
!-----------------------------------------------
        CALL second0(ton)

        IF (.NOT. lcurr_init) STOP 'MUST CALL init_state(.TRUE.) BEFORE bdotj'

! 
!     Compute contravariant components of magnetic field
!       
        ALLOCATE(bsupsijsh(ntheta,nzeta,ns),                            &
                 bsupuijch(ntheta,nzeta,ns),                            &
                 bsupvijch(ntheta,nzeta,ns), stat=istat)
        bsupsijsh=0; bsupuijch=0; bsupvijch=0
        IF (istat .NE. 0) STOP 'Allocation #1 failed in BDOTJ' 

        CALL toijsp(jbsupsmnsh, bsupsijsh, 0, 0, 1, 1)                           
        bsupsijsh = bsupsijsh/jacobh
        CALL toijsp(jbsupumnch, bsupuijch, 0, 0, 0, 1)       
        bsupuijch = bsupuijch/jacobh
        CALL toijsp(jbsupvmnch, bsupvijch, 0, 0, 0, 1) 
        bsupvijch = bsupvijch/jacobh
!
!   Compute covariant field. Used for normalization.
!        
        ALLOCATE(bsubsijsh(ntheta,nzeta,ns),                            &
                 bsubuijch(ntheta,nzeta,ns),                            &
                 bsubvijch(ntheta,nzeta,ns), stat=istat)
        IF (istat .NE. 0) STOP 'Allocation #2 failed in BDOTJ' 
        bsubsijsh=0; bsubvijch=0; bsubuijch=0
      
        CALL tolowerh(bsupsijsh, bsupuijch, bsupvijch,                  &
                      bsubsijsh, bsubuijch, bsubvijch)

!
!    In initial state, calculate currents (remember they are multiplied by a jacobian!)
!
!
!    Move contravariant currents (*jacob) to the half mesh
!       
        ALLOCATE(ksupsijsh(ntheta,nzeta,ns),                            &
                 ksupuijch(ntheta,nzeta,ns),                            &
                 ksupvijch(ntheta,nzeta,ns), stat=istat)
        IF (istat .NE. 0) STOP 'Allocation #3 failed in BDOTJ'

        CALL to_half_mesh(ksupsijsf0, ksupsijsh, nuv)
        CALL to_half_mesh(ksupuijcf0, ksupuijch, nuv)
        CALL to_half_mesh(ksupvijcf0, ksupvijch, nuv)
!
!   Remove jacobian 
!        
        ksupsijsh = ksupsijsh/jacobh      
        ksupuijch = ksupuijch/jacobh  
        ksupvijch = ksupvijch/jacobh
!          
!   Get covariant currents
!  
        ALLOCATE(ksubsijsh(ntheta,nzeta,ns),                            &
                 ksubuijch(ntheta,nzeta,ns),                            &
                 ksubvijch(ntheta,nzeta,ns), stat=istat)
        IF (istat .NE. 0) STOP 'Allocation #4 failed in BDOTJ'
        ksubsijsh=0; ksubuijch=0; ksubvijch=0

        CALL tolowerh(ksupsijsh, ksupuijch, ksupvijch,                  &
                      ksubsijsh, ksubuijch, ksubvijch)
!
!   Compute B*J in the half mesh.
!      
        IF (ALLOCATED(bdotjmnch) .AND. SIZE(bdotjmnch,3).NE. ns)        &
            DEALLOCATE(bdotjmnch)
        IF (.NOT. ALLOCATED(bdotjmnch))  THEN                            ! Should already be allocated in QUANTITIES       
          ALLOCATE(bdotjmnch(0:mpol,-ntor:ntor,ns), stat=istat)          ! Otherwise, allocate first time called
          IF (istat .NE. 0) STOP 'Allocation #5 failed in BDOTJ'
        ENDIF
        ALLOCATE(bdotjijch(ntheta,nzeta,ns), stat=istat)
        IF (istat .NE. 0) STOP 'Allocation #6 failed in BDOTJ'
        bdotjijch=0
        
        bdotjijch(:,:,3:nsh-1) =                                        &   ! RS: If I include 2 and nsh, it gets huge
          bsupsijsh(:,:,3:nsh-1)*ksubsijsh(:,:,3:nsh-1) +               &   ! something funny happens at the boundaries.
          bsupuijch(:,:,3:nsh-1)*ksubuijch(:,:,3:nsh-1) +               &
          bsupvijch(:,:,3:nsh-1)*ksubvijch(:,:,3:nsh-1)                     ! B*J = (B^s*J_s + B^u*J_u + B^v*J_v)
! 
!   Compute tnorm = |JXB|^2 for normalization; thus we print out J_parallel/|J_perp|.
!
        ALLOCATE(work1(ntheta,nzeta,ns), work2(ntheta,nzeta,ns),        &
                 work3(ntheta,nzeta,ns), work4(ntheta,nzeta,ns),        &
                 work5(ntheta,nzeta,ns), work6(ntheta,nzeta,ns),        &
                 stat=istat)
        IF (istat .NE. 0) STOP 'Allocation #7 failed in BDOTJ'
        work1=0; work2=0; work3=0
        work4=0; work5=0; work6=0
        
        work1 = (bsubuijch*ksubvijch - bsubvijch*ksubuijch)/jacobh   ! (JxB)^s  = (1/sqrt(g))*(K_uB_b-K_vB_u)
        work2 = (bsubvijch*ksubsijsh - bsubsijsh*ksubvijch)/jacobh   ! (JxB)^u  = .....
        work3 = (bsubsijsh*ksubuijch - bsubuijch*ksubsijsh)/jacobh   ! (JxB)^v  = .....
        CALL tolowerh(work1, work2, work3, work4, work5, work6)
        
        work1 = work1*work4 + work2*work5 + work3*work6                      ! |JxB|**2
        tnorm = SUM(work1(:,:,3:nsh-1))
        IF (iam.EQ.0 .AND. ANY(work1 .LT. zero))                        &
           PRINT *,'ERROR1 IN BDOTJ, JXB^2 < 0!'
        tnorm = SQRT(ABS(tnorm))
        DEALLOCATE(work3, work4, work5, work6)
        
        work1 = bsupsijsh*bsubsijsh + bsupuijch*bsubuijch +             &     ! |B|**2
                bsupvijch*bsubvijch
        work2 = ksupsijsh*ksubsijsh + ksupuijch*ksubuijch +             &     ! |J|**2
                ksupvijch*ksubvijch
       
        tnorm2 = SUM(work1(:,:,3:nsh-1)*work2(:,:,3:nsh-1))
        work1(:,:,1) = 0;  work2(:,:,1) = 0
        IF (iam.EQ.0 .AND. ANY(work1.LT.zero))                          &
           PRINT *, 'ERROR2 IN BDOTJ: B^2 < 0!'               
        IF (iam.EQ.0 .AND. ANY(work2.LT.zero))                          &
           PRINT *, 'ERROR3 IN BDOTJ: J^2 < 0!'               
        tnorm2 = SQRT(ABS(tnorm2))                                              ! |J|*|B|
        bdotj_rms = SUM(bdotjijch*bdotjijch)

        DEALLOCATE (bsupsijsh, bsupuijch, bsupvijch, bsubsijsh,         &
                    bsubuijch, bsubvijch, ksupsijsh, ksupuijch,         &
                    ksupvijch, ksubsijsh, ksubuijch, ksubvijch,         &      
                    work1, work2)

        tnorm = MAX(tnorm, EPSILON(tnorm))
        tnorm2= MAX(tnorm2,EPSILON(tnorm2))
        bdotj2_rms = SQRT(ABS(bdotj_rms))/tnorm2                             ! RMS of bdotj/|J|*|B|  
        bdotj_rms  = SQRT(ABS(bdotj_rms))/tnorm                              ! RMS of bdotj/|JxB|  

        bdotjijch = bdotjijch/tnorm
        
        CALL tomnsp(bdotjijch, bdotjmnch, 0, 1)                              ! Keep harmonics of BDOTJ for output   

        DEALLOCATE (bdotjijch)

        CALL second0(toff)
        time_bdotj = time_bdotj+(toff-ton)

        END SUBROUTINE bdotj

      SUBROUTINE bdotj_par
#if defined(SKS)
      USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR
      INCLUDE 'mpif.h'
!     
!     WRITTEN 08-20-07 BY R. SANCHEZ AS PART OF THE ORNL SIESTA PROJECT (c)
!     UPDATED 03-21-13 BY S. HIRSHMAN FOR PARALLEL EXECUTION
!     
!     PURPOSE: Computes parallel current normalized to the total current (J*B/|JxB|) on the half mesh
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER      :: istat, nsmin, nsmax, ns_span, n1, n2
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: bdotjijch
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: bsupsijsh,             &
          bsupuijch, bsupvijch,                                         &
          bsubsijsh, bsubuijch, bsubvijch, ksubsijsh, ksubuijch,        &
          ksubvijch, ksupsijsh, ksupuijch, ksupvijch, work1, work2,     &
          work3, work4, work5, work6
      REAL(rprec)  :: tnorm, tnorm2, s1(3), r1(3), ton, toff
!-----------------------------------------------
      CALL second0(ton)

      IF (.NOT. lcurr_init) STOP 'MUST CALL init_state_par(.TRUE.) BEFORE bdotj_par'
      nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow,ns)
      ns_span=nsmax-nsmin+1
! 
!     Compute contravariant components of magnetic field
!       
        ALLOCATE(bsupsijsh(ntheta,nzeta,nsmin:nsmax),                   &
                 bsupuijch(ntheta,nzeta,nsmin:nsmax),                   &
                 bsupvijch(ntheta,nzeta,nsmin:nsmax), stat=istat)
        IF (istat .NE. 0) STOP 'Allocation #1 failed in BDOTJ_PAR' 

        CALL toijsp_par(jbsupsmnsh_p(:,:,nsmin:nsmax), bsupsijsh, 0, 0, 1, 1, 1, ns_span)                           
        bsupsijsh = bsupsijsh/jacobh_p(:,:,nsmin:nsmax)
        CALL toijsp_par(jbsupumnch_p(:,:,nsmin:nsmax), bsupuijch, 0, 0, 0, 1, 1, ns_span)       
        bsupuijch = bsupuijch/jacobh_p(:,:,nsmin:nsmax)
        CALL toijsp_par(jbsupvmnch_p(:,:,nsmin:nsmax), bsupvijch, 0, 0, 0, 1, 1, ns_span) 
        bsupvijch = bsupvijch/jacobh_p(:,:,nsmin:nsmax)
!
!   Compute covariant field. Used for normalization.
!        
        ALLOCATE(bsubsijsh(ntheta,nzeta,nsmin:nsmax),                   &
                 bsubuijch(ntheta,nzeta,nsmin:nsmax),                   &
                 bsubvijch(ntheta,nzeta,nsmin:nsmax), stat=istat)
        IF (istat .NE. 0) STOP 'Allocation #2 failed in BDOTJ_PAR' 
      
        CALL tolowerh_par(bsupsijsh, bsupuijch, bsupvijch,              &
                          bsubsijsh, bsubuijch, bsubvijch, nsmin, nsmax)
!
!    In initial state, calculate currents (remember they are multiplied by a jacobian!)
!
!
!    Move contravariant currents (*jacob) to the half mesh
!       
        nsmin=MAX(1,startglobrow-1)
        ALLOCATE(ksupsijsh(ntheta,nzeta,nsmin:nsmax),                   &
                 ksupuijch(ntheta,nzeta,nsmin:nsmax),                   &
                 ksupvijch(ntheta,nzeta,nsmin:nsmax), stat=istat)
        IF (istat .NE. 0) STOP 'Allocation #3 failed in BDOTJ_PAR'

        CALL to_half_mesh_par(ksupsijsf0, ksupsijsh, nuv, nsmin, nsmax)
        CALL to_half_mesh_par(ksupuijcf0, ksupuijch, nuv, nsmin, nsmax)
        CALL to_half_mesh_par(ksupvijcf0, ksupvijch, nuv, nsmin, nsmax)
!
!   Remove jacobian 
!        
        ksupsijsh = ksupsijsh/jacobh_p(:,:,nsmin:nsmax)
        ksupuijch = ksupuijch/jacobh_p(:,:,nsmin:nsmax)
        ksupvijch = ksupvijch/jacobh_p(:,:,nsmin:nsmax)

!          
!   Get covariant currents
!  
        ALLOCATE(ksubsijsh(ntheta,nzeta,nsmin:nsmax),                   &
                 ksubuijch(ntheta,nzeta,nsmin:nsmax),                   &
                 ksubvijch(ntheta,nzeta,nsmin:nsmax), stat=istat)
        IF (istat .NE. 0) STOP 'Allocation #4 failed in BDOTJ_PAR'

        CALL tolowerh_par(ksupsijsh, ksupuijch, ksupvijch,              &
                          ksubsijsh, ksubuijch, ksubvijch, nsmin, nsmax)
        nsmin=MAX(1,startglobrow)
!
!   Compute B*J in the half mesh.
!      
        IF (ALLOCATED(bdotjmnch) .AND. SIZE(bdotjmnch,3).NE.ns_span)    &
            DEALLOCATE(bdotjmnch)
        IF (.NOT. ALLOCATED(bdotjmnch))  THEN                            
          ALLOCATE(bdotjmnch(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
          IF (istat .NE. 0) STOP 'Allocation #5 failed in BDOTJ_PAR'
        ENDIF
        ALLOCATE(bdotjijch(ntheta,nzeta,nsmin:nsmax), stat=istat)
        IF (istat .NE. 0) STOP 'Allocation #6 failed in BDOTJ_PAR'
        bdotjijch=0
        
        n1 = MAX(3,nsmin); n2 = MIN(nsh-1, nsmax)
        bdotjijch(:,:,n1:n2) =                                          &   ! RS: If I include 2 and nsh, it gets huge
          bsupsijsh(:,:,n1:n2)*ksubsijsh(:,:,n1:n2) +                   &   ! something funny happens at the boundaries.
          bsupuijch(:,:,n1:n2)*ksubuijch(:,:,n1:n2) +                   &
          bsupvijch(:,:,n1:n2)*ksubvijch(:,:,n1:n2)                         ! B*J = (B^s*J_s + B^u*J_u + B^v*J_v)
! 
!   Compute tnorm = |JXB|^2 for normalization; thus we print out J_parallel/|J_perp|.
!
        ALLOCATE(work1(ntheta,nzeta,nsmin:nsmax),                       &
                 work2(ntheta,nzeta,nsmin:nsmax),                       &
                 work3(ntheta,nzeta,nsmin:nsmax),                       &
                 work4(ntheta,nzeta,nsmin:nsmax),                       &
                 work5(ntheta,nzeta,nsmin:nsmax),                       &
                 work6(ntheta,nzeta,nsmin:nsmax),                       &
                 stat=istat)
        IF (istat .NE. 0) STOP 'Allocation #7 failed in BDOTJ_PAR'
        
        work1 = bsubuijch*ksubvijch(:,:,nsmin:nsmax)                    &
              - bsubvijch*ksubuijch(:,:,nsmin:nsmax)                      ! (JxB)^s  = (1/sqrt(g))*(K_uB_b-K_vB_u)
        work2 = bsubvijch*ksubsijsh(:,:,nsmin:nsmax)                    & 
              - bsubsijsh*ksubvijch(:,:,nsmin:nsmax)                      ! (JxB)^u  = .....
        work3 = bsubsijsh*ksubuijch(:,:,nsmin:nsmax)                    &
              - bsubuijch*ksubsijsh(:,:,nsmin:nsmax)                      ! (JxB)^v  = .....
        work1 = work1/jacobh_p(:,:,nsmin:nsmax)
        work2 = work2/jacobh_p(:,:,nsmin:nsmax)
        work3 = work3/jacobh_p(:,:,nsmin:nsmax)
        CALL tolowerh_par(work1, work2, work3,                          &
                          work4, work5, work6, nsmin, nsmax)
        
        work1 = work1*work4 + work2*work5 + work3*work6                     ! |JxB|**2
        tnorm = SUM(work1(:,:,n1:n2))
        IF (nsmin .EQ. 1) work1(:,:,1) = 0
        IF (iam.EQ.0 .AND. ANY(work1 .LT. zero))                        &
           PRINT *,'ERROR1 IN BDOTJ_PAR, JXB^2 < 0!'
        DEALLOCATE(work3, work4, work5, work6)
        
        work1 = bsupsijsh*bsubsijsh + bsupuijch*bsubuijch +             &   ! |B|**2
                bsupvijch*bsubvijch
        work2 = ksupsijsh(:,:,nsmin:nsmax)*ksubsijsh(:,:,nsmin:nsmax)   &
              + ksupuijch(:,:,nsmin:nsmax)*ksubuijch(:,:,nsmin:nsmax)   &
              + ksupvijch(:,:,nsmin:nsmax)*ksubvijch(:,:,nsmin:nsmax)         ! |J|**2
       
        IF (nsmin .EQ. 1) THEN
           work1(:,:,1) = 0;  work2(:,:,1) = 0
        END IF
        IF (iam.EQ.0 .AND. ANY(work1.LT.zero))                          &
           PRINT *, 'ERROR2 IN BDOTJ_PAR: B^2 < 0!'               
        IF (iam.EQ.0 .AND. ANY(work2.LT.zero))                          &
           PRINT *, 'ERROR3 IN BDOTJ_PAR: J^2 < 0!'               

        tnorm2 = SUM(work1(:,:,n1:n2)*work2(:,:,n1:n2))
        bdotj_rms = SUM(bdotjijch*bdotjijch)
        s1(1)=tnorm; s1(2)=tnorm2; s1(3)=bdotj_rms
        CALL MPI_ALLREDUCE(s1,r1,3,MPI_REAL8, MPI_SUM,                  &
                           MPI_COMM_WORLD,MPI_ERR)
        tnorm  = SQRT(ABS(r1(1)))
        tnorm2 = SQRT(ABS(r1(2)))                                             ! |J|*|B|
        bdotj_rms = r1(3) 

        DEALLOCATE (bsupsijsh, bsupuijch, bsupvijch, bsubsijsh,         &
                    bsubuijch, bsubvijch, ksupsijsh, ksupuijch,         &
                    ksupvijch, ksubsijsh, ksubuijch, ksubvijch,         &      
                    work1, work2)

!        IF (iam .eq. 0) WRITE (3700, *) 'tnorm: ',tnorm,' tnorm2: ', tnorm2, &
!                                        ' bdotj_rms: ', bdotj_rms

        tnorm = MAX(tnorm, EPSILON(tnorm))
        tnorm2= MAX(tnorm2,EPSILON(tnorm2))
        bdotj2_rms = SQRT(ABS(bdotj_rms))/tnorm2                             ! RMS of bdotj/|J|*|B|  
        bdotj_rms  = SQRT(ABS(bdotj_rms))/tnorm                              ! RMS of bdotj/|JxB|  

        bdotjijch = bdotjijch/tnorm
        
        CALL tomnsp_par(bdotjijch, bdotjmnch, 0, 1, 1, ns_span)              ! Keep harmonics of BDOTJ for output   

        DEALLOCATE (bdotjijch)

        CALL second0(toff)
        time_bdotj = time_bdotj+(toff-ton)
#endif
        END SUBROUTINE bdotj_par
      
        SUBROUTINE divj
!     
!     WRITTEN 09-10-07 BY R. SANCHEZ AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: Computes divergence of the current normalized to the total current on the half mesh
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER               :: m, n, istat
        REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: temp1, temp2
        REAL(dp)              :: tnorm, ton, toff
!-----------------------------------------------
!
!     COMPUTE divJ ON HALF MESH POINTS
!
        CALL second0(ton)

        IF (.NOT. lcurr_init) STOP 'MUST CALL init_state(.TRUE.) BEFORE divJ'

        IF (.NOT. ALLOCATED(divjmnsh))  THEN                                  
          ALLOCATE(divjmnsh(0:mpol,-ntor:ntor,ns),  stat=istat)          ! Otherwise, allocate first time called
          IF (istat .NE. 0) STOP 'Allocation failed in DIVJ'
          divjmnsh = zero
        ENDIF
!
!  Compute divergence of J. It is a SINE parity quantity.
!
        DO m = 0, mpol
           DO n = -ntor,ntor
             divjmnsh(m,n,3:nsh) =                                      &
               -m*(ksupumncf(m,n,3:nsh)+ksupumncf(m,n,2:nsh-1))/2       &
               -n*nfp*(ksupvmncf(m,n,3:nsh)+ksupvmncf(m,n,2:nsh-1))/2   &
               + ohs *(ksupsmnsf(m,n,3:nsh)-ksupsmnsf(m,n,2:nsh-1))
           END DO
        END DO

 
        ALLOCATE(temp1(ntheta,nzeta,ns),                                &
                 temp2(ntheta,nzeta,ns), stat=istat)                    
        IF (istat .NE. 0) STOP 'Allocation failed in DIVJ'
        temp1 = 0; temp2 = 0
 
        temp1 = ksupsijsf0*ksubsijsf + ksupuijcf0*ksubuijcf +           & ! Norm: |J|^2
                ksupvijcf0*ksubvijcf
        CALL to_half_mesh(temp1, temp2, nuv)
        tnorm = SUM(temp2/jacobh**2)
        IF (iam.EQ. 0 .AND. ANY(temp2.lt.zero)) PRINT *,'ERROR1 in DIVJ: |J|^2 < 0!'
        tnorm = SQRT(ABS(tnorm))
 
        DEALLOCATE(temp1, temp2, stat=istat)
 !
 !   Compute rms of divergence of J
 !     
        divjmnsh = divjmnsh/tnorm
        divj_rms = SQRT(SUM(divjmnsh**2))

!        IF (iam .EQ. 0) WRITE(4100,*)'DIVJ: tnorm=',tnorm,' divj_rms: ',divj_rms
        CALL second0(toff)
        time_divj = time_divj+(toff-ton)

        END SUBROUTINE divj

        SUBROUTINE divj_par
!     
!     WRITTEN 03-18-13 BY S. HIRSHMAN AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: Computes divergence of the current normalized to the total current on the half mesh
!
#if defined(SKS)
        USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR
        INCLUDE 'mpif.h'
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER                  :: m, n, istat, nsmin, nsmax, n1, n2
        REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE ::  temp1, temp2
        REAL(rprec)              :: tnorm, s1(2), r1(2), ton, toff
!-----------------------------------------------
!
!     COMPUTE divJ ON HALF MESH POINTS
!
        CALL second0(ton)

        IF (.NOT. lcurr_init) STOP 'MUST CALL init_state_par(.TRUE.) BEFORE divJ_par'

        nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow,ns)

        IF (.NOT. ALLOCATED(divjmnsh))  THEN                                  
          ALLOCATE(divjmnsh(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)   ! Otherwise, allocate first time called
          IF (istat .NE. 0) STOP 'Allocation failed in DIVJ'
          divjmnsh = 0
        ENDIF

!
!  Compute divergence of J. It is a SINE parity quantity.
!
        n1 = MAX(3, nsmin); n2 = MIN(ns-1, nsmax)
        DO m = 0, mpol
           DO n = -ntor,ntor
             divjmnsh(m,n,n1:n2) =                                      &
               -m*(ksupumncf(m,n,n1:n2)+ksupumncf(m,n,n1-1:n2-1))/2     &
           -n*nfp*(ksupvmncf(m,n,n1:n2)+ksupvmncf(m,n,n1-1:n2-1))/2     &
           + ohs *(ksupsmnsf(m,n,n1:n2)-ksupsmnsf(m,n,n1-1:n2-1))
           END DO
        END DO
!
!    Compute covariant components of current (times jacobian). Used for normalization. 
!  
        nsmin=MAX(1,startglobrow-1); nsmax=MIN(endglobrow+1,ns)
        ALLOCATE(temp1(ntheta,nzeta,nsmin:nsmax),                       &
                 temp2(ntheta,nzeta,nsmin:nsmax), stat=istat)                    
        IF (istat .NE. 0) STOP 'Allocation failed in DIVJ'
        temp1 = 0; temp2 = 0

        temp1 = ksupsijsf0(:,:,nsmin:nsmax)*ksubsijsf(:,:,nsmin:nsmax)  &
              + ksupuijcf0(:,:,nsmin:nsmax)*ksubuijcf(:,:,nsmin:nsmax)  & ! Norm: |J|^2
              + ksupvijcf0(:,:,nsmin:nsmax)*ksubvijcf(:,:,nsmin:nsmax)

        CALL to_half_mesh_par(temp1, temp2, nuv, nsmin, nsmax)

        temp2 = temp2/jacobh_p(:,:,nsmin:nsmax)**2

        nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow,ns)
        tnorm = SUM(temp2(:,:,nsmin:nsmax))

        DEALLOCATE(temp1, temp2, stat=istat)

        divj_rms = SUM(divjmnsh(:,:,nsmin:nsmax)**2)

        s1(1) = divj_rms; s1(2) = tnorm
        CALL MPI_ALLREDUCE(s1,r1,2,MPI_REAL8, MPI_SUM,                  &
                           MPI_COMM_WORLD,MPI_ERR)
        divj_rms=r1(1);   tnorm = r1(2)
   
 !   Compute rms of divergence of J
 !     
        IF (tnorm .GE. zero) THEN
           divjmnsh = divjmnsh/SQRT(tnorm)
           divj_rms = SQRT(divj_rms/tnorm)
        END IF

!        IF (iam .EQ. 0) WRITE(3100,*)'DIVJ_PAR: tnorm=',SQRT(tnorm),' divj_rms: ',divj_rms

        CALL second0(toff)
        time_divj = time_divj+(toff-ton)
#endif

        END SUBROUTINE divj_par

       SUBROUTINE dealloc_diagnostics

       IF (ALLOCATED(divbmnsf)) DEALLOCATE(divbmnsf)
       IF (ALLOCATED(divjmnsh)) DEALLOCATE(divjmnsh)
       IF (ALLOCATED(bdotjmnch)) DEALLOCATE(bdotjmnch)

       END SUBROUTINE dealloc_diagnostics


      END MODULE diagnostics_mod
