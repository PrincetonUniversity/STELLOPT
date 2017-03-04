      SUBROUTINE allocate_boozer (iread)
      USE booz_params
      USE booz_persistent
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, jrad, istat1=0, istat2=0, istat3=0, iread, index
      INTEGER :: ier_arr(32)
      CHARACTER(LEN=1000) :: temp
      CHARACTER(LEN=10)   :: scanset='0123456789'
c-----------------------------------------------
      IF (iread < 0) THEN
         jsize = -iread
         IF (ALLOCATED(jlist)) DEALLOCATE(jlist)
         ALLOCATE(jlist(jsize), stat=istat1)
         IF (istat1 .ne. 0) STOP 'Unable to allocate jlist/lsurf_boz'
         GOTO 327
      END IF
      IF (.not.ALLOCATED(jlist)) ALLOCATE (jlist(ns), lsurf_boz(ns), 
     1    stat=istat1)
      IF (istat1 .ne. 0) STOP 'Unable to allocate jlist/lsurf_boz'

!
!     Read in (and parse) list of radial surfaces to compute
!     Surfaces SHOULD ALL be on a single line, but IFORT doesn't write it out that way
!     so the loop below will correctly read multiple lines from the in_booz file
!
      jlist = 0
      i = 1
      READ (iread, '(a)', iostat=istat1) temp
      IF (istat1 .eq. 0) THEN
         DO WHILE (istat1 .eq. 0)
            DO jrad = i, ns
               index = SCAN(temp,scanset)
               IF (index < 1) EXIT
               temp = temp(index:)
               READ(temp, *) jlist(jrad)
               IF (jlist(jrad) < 10) THEN
                  temp = temp(3:)
               ELSE IF (jlist(jrad) < 100) THEN
                  temp = temp(4:)
               ELSE
                  temp = temp(5:)
               END IF
            END DO
            READ (iread, '(a)', iostat=istat1) temp
            i = jrad
         END DO

         lsurf_boz = .FALSE.
         DO jrad = 1, ns
            i = jlist(jrad)
            IF (i.gt.1 .and. i.le.ns) lsurf_boz(i) = .TRUE.
         END DO

      ELSE IF (istat1 .ne. 0) THEN 
         WRITE(6, '(a,/,a,i4)')
     1    ' No jlist data was found in Boozer input file.',
     1    ' Will assume that all surfaces are needed.',
     1    ' Iostat: ', istat1
         lsurf_boz = .TRUE. 
         lsurf_boz(1) = .FALSE.

      END IF
      
      jsize = COUNT(lsurf_boz(1:ns))
!
!     Recompute jlist, just in case user used unordered (or repeated) indices
!
      DEALLOCATE (jlist)
      ALLOCATE (jlist(jsize), stat=istat1)
      IF (istat1 .ne. 0) STOP 'Unable to allocate jlist'

  327 i = 1
      DO jrad = 2, ns
         IF (lsurf_boz(jrad)) THEN
            jlist(i) = jrad
            i = i+1
         END IF
      END DO
      
      !  ADDED by SAL for consistency
      ier_arr = 0
      IF (ALLOCATED(bsubumnc)) DEALLOCATE(bsubumnc)
      ALLOCATE(bsubumnc(mnmax_nyq,ns), stat=ier_arr(1))
      IF (ALLOCATED(bsubvmnc)) DEALLOCATE(bsubvmnc);
      ALLOCATE(bsubvmnc(mnmax_nyq,ns), stat=ier_arr(2))
      IF (ALLOCATED(bmodmnc)) DEALLOCATE(bmodmnc);
      ALLOCATE(bmodmnc(mnmax_nyq,ns), stat=ier_arr(3))
      IF (ALLOCATED(rmnc)) DEALLOCATE(rmnc); 
      ALLOCATE(rmnc(mnmax,ns), stat=ier_arr(4))
      IF (ALLOCATED(zmns)) DEALLOCATE(zmns); 
      ALLOCATE(zmns(mnmax,ns), stat=ier_arr(5))
      IF (ALLOCATED(lmns)) DEALLOCATE(lmns); 
      ALLOCATE(lmns(mnmax,ns), stat=ier_arr(6))
      IF (ALLOCATED(xm)) DEALLOCATE(xm); 
      ALLOCATE(xm(mnmax), stat=ier_arr(7))
      IF (ALLOCATED(xn)) DEALLOCATE(xn); 
      ALLOCATE(xn(mnmax), stat=ier_arr(8))
      IF (ALLOCATED(xm_nyq)) DEALLOCATE(xm_nyq);
      ALLOCATE(xm_nyq(mnmax_nyq), stat=ier_arr(9))
      IF (ALLOCATED(xn_nyq)) DEALLOCATE(xn_nyq);
      ALLOCATE(xn_nyq(mnmax_nyq), stat=ier_arr(10))
      IF (ALLOCATED(pmns)) DEALLOCATE(pmns);
      ALLOCATE(pmns(mnmax_nyq), stat=ier_arr(11))
      IF (ALLOCATED(hiota)) DEALLOCATE(hiota);
      ALLOCATE(hiota(ns), stat=ier_arr(12))
      IF (ALLOCATED(phip)) DEALLOCATE(phip);
      ALLOCATE(phip(ns), stat=ier_arr(13))
      IF (ALLOCATED(gpsi)) DEALLOCATE(gpsi); 
      ALLOCATE(gpsi(ns), stat=ier_arr(14))
      IF (ALLOCATED(ipsi)) DEALLOCATE(ipsi);
      ALLOCATE(ipsi(ns), stat=ier_arr(15))
      IF (ALLOCATED(beta_vol)) DEALLOCATE(beta_vol);
      ALLOCATE(beta_vol(ns), stat=ier_arr(16))
      IF (ALLOCATED(pres)) DEALLOCATE(pres);
      ALLOCATE(pres(ns), stat=ier_arr(17))
      IF (ALLOCATED(phi)) DEALLOCATE(phi);
      ALLOCATE(phi(ns), stat=ier_arr(18))
      IF (ALLOCATED(buco)) DEALLOCATE(buco);
      ALLOCATE(buco(ns), stat=ier_arr(19))
      IF (ALLOCATED(bvco)) DEALLOCATE(bvco);
      ALLOCATE(bvco(ns), stat=ier_arr(20))
      IF (ALLOCATED(rmncb)) DEALLOCATE(rmncb);
      ALLOCATE(rmncb(mnboz,jsize), stat=ier_arr(21))
      IF (ALLOCATED(zmnsb)) DEALLOCATE(zmnsb);
      ALLOCATE(zmnsb(mnboz,jsize), stat=ier_arr(22))
      IF (ALLOCATED(pmnsb)) DEALLOCATE(pmnsb);
      ALLOCATE(pmnsb(mnboz,jsize), stat=ier_arr(23))
      IF (ALLOCATED(gmncb)) DEALLOCATE(gmncb); 
      ALLOCATE(gmncb(mnboz,jsize), stat=ier_arr(24))
      IF (ALLOCATED(bmncb)) DEALLOCATE(bmncb); 
      ALLOCATE(bmncb(mnboz,jsize), stat=ier_arr(25))
      IF (ALLOCATED(bmod_b)) DEALLOCATE(bmod_b);
      ALLOCATE(bmod_b(nv_boz,nu_boz), stat=ier_arr(26))
      IF(ANY(ier_arr.ne.0)) istat1=1
      
!      IF (.not.ALLOCATED(bsubumnc)) ALLOCATE(
!     1    bsubumnc(mnmax_nyq,ns), bsubvmnc(mnmax_nyq,ns),
!     1    bmodmnc(mnmax_nyq,ns),
!     2    rmnc(mnmax,ns), zmns(mnmax,ns), lmns(mnmax,ns),
!     3    xm(mnmax), xn(mnmax),
!     3    xm_nyq(mnmax_nyq), xn_nyq(mnmax_nyq),
!     4    hiota(ns), phip(ns), gpsi(ns), ipsi(ns), pmns(mnmax_nyq),
!     5    beta_vol(ns), pres(ns), phi(ns), buco(ns), bvco(ns),
!     5    rmncb(mnboz,jsize), zmnsb(mnboz,jsize), pmnsb(mnboz,jsize), 
!     6    gmncb(mnboz,jsize), bmncb(mnboz,jsize), 
!     7    bmod_b(nv_boz,nu_boz), stat=istat1 )
   
      ier_arr=0
      IF (ALLOCATED(sfull)) DEALLOCATE(sfull);
      ALLOCATE(sfull(ns), stat=ier_arr(1))
      IF (ALLOCATED(scl)) DEALLOCATE(scl);
      ALLOCATE(scl(mnboz), stat=ier_arr(2))
      IF (ALLOCATED(xmb)) DEALLOCATE(xmb);
      ALLOCATE(xmb(mnboz), stat=ier_arr(3))
      IF (ALLOCATED(xnb)) DEALLOCATE(xnb);
      ALLOCATE(xnb(mnboz), stat=ier_arr(4))
      IF(ANY(ier_arr.ne.0)) istat2=1

!      IF (.not.ALLOCATED(sfull)) ALLOCATE(
!     1    sfull(ns), scl(mnboz), xmb(mnboz), xnb(mnboz), stat=istat2)
   
      ier_arr=0
      IF (ALLOCATED(bsubumns)) DEALLOCATE(bsubumns);
      ALLOCATE(bsubumns(mnmax_nyq,ns), stat=ier_arr(1))
      IF (ALLOCATED(bsubvmns)) DEALLOCATE(bsubvmns);
      ALLOCATE(bsubvmns(mnmax_nyq,ns), stat=ier_arr(2))
      IF (ALLOCATED(bmodmns)) DEALLOCATE(bmodmns);
      ALLOCATE(bmodmns(mnmax_nyq,ns), stat=ier_arr(3))
      IF (ALLOCATED(rmns)) DEALLOCATE(rmns);
      ALLOCATE(rmns(mnmax,ns), stat=ier_arr(4))
      IF (ALLOCATED(zmnc)) DEALLOCATE(zmnc);
      ALLOCATE(zmnc(mnmax,ns), stat=ier_arr(5))
      IF (ALLOCATED(lmnc)) DEALLOCATE(lmnc);
      ALLOCATE(lmnc(mnmax,ns), stat=ier_arr(6))
      IF (ALLOCATED(pmnc)) DEALLOCATE(pmnc);
      ALLOCATE(pmnc(mnmax_nyq), stat=ier_arr(7))
      IF (ALLOCATED(rmnsb)) DEALLOCATE(rmnsb);
      ALLOCATE(rmnsb(mnboz,jsize), stat=ier_arr(8))
      IF (ALLOCATED(zmncb)) DEALLOCATE(zmncb);
      ALLOCATE(zmncb(mnboz,jsize), stat=ier_arr(9))
      IF (ALLOCATED(pmncb)) DEALLOCATE(pmncb);
      ALLOCATE(pmncb(mnboz,jsize), stat=ier_arr(10))
      IF (ALLOCATED(gmnsb)) DEALLOCATE(gmnsb);
      ALLOCATE(gmnsb(mnboz,jsize), stat=ier_arr(11))
      IF (ALLOCATED(bmnsb)) DEALLOCATE(bmnsb);
      ALLOCATE(bmnsb(mnboz,jsize), stat=ier_arr(12))
      IF(ANY(ier_arr.ne.0)) istat3=1

!      IF (.not.ALLOCATED(bsubumns)) ALLOCATE(
!!      IF (lasym_b .AND. .not.ALLOCATED(bsubumns)) ALLOCATE(
!     1    bsubumns(mnmax_nyq,ns), bsubvmns(mnmax_nyq,ns),
!     1    bmodmns(mnmax_nyq,ns), 
!     1    rmns(mnmax,ns), zmnc(mnmax,ns), lmnc(mnmax,ns),
!     4    pmnc(mnmax_nyq),
!     5    rmnsb(mnboz,jsize), zmncb(mnboz,jsize), pmncb(mnboz,jsize), 
!     6    gmnsb(mnboz,jsize), bmnsb(mnboz,jsize), 
!     1    stat=istat3)

      IF (istat1.ne.0 .or. istat2.ne.0 .or. istat3.ne.0) THEN
          PRINT *,' problem allocating boozer memory'
          PRINT *,' istat1 = ',istat1,' istat2 = ',istat2,
     1            ' istat3 = ',istat3
          STOP
      ENDIF

      END SUBROUTINE allocate_boozer
