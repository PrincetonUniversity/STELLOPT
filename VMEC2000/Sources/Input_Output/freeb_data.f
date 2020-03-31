      SUBROUTINE freeb_data (rmnc, zmns, rmns, zmnc, bmodmn, bmodmn1)
      USE vmec_main
      USE vacmod
      USE realspace, ONLY: r1, z1
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(mnmax) :: rmnc, zmns, rmns, zmnc, 
     1                                 bmodmn, bmodmn1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: iprint, nzskip, i, l, k, lk, mn,
     1           mn0, n, nedge, nedge0 = 99, iu, iv, nl, lkr
      REAL(rprec) :: zeta, potsin, potcos
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: rb, phib, zb
C-----------------------------------------------
!
!     WRITE OUT EDGE VALUES OF FIELDS TO FORT.NEDGE0 (INCLUDE REFLECTED POINT)
!
!     NOTE: BR, BPHI, BZ WERE COMPUTED IN BSS, CALLED FROM EQFOR
!
      IF (ivac.le.0 .or. (.not.lfreeb .and. .not.ledge_dump)) RETURN

      ALLOCATE (rb(2*nznt), phib(2*nznt), zb(2*nznt), stat=l)
      IF (l .ne. 0) STOP 'allocation error in freeb_data'

      nedge = 0
      lkr = nznt
      DO iv = 1,nzeta
         zeta = (twopi*(iv-1))/(nzeta*nfp)
         DO iu = 1,ntheta3
            lk = iv + nzeta*(iu-1)
            nl = ns*lk
            nedge = nedge+1
            rb(lk)   = r1(nl,0) + r1(nl,1)
            phib(lk) = zeta
            zb(lk)   = z1(nl,0) + z1(nl,1)
!
!      INCLUDE -u,-v POINTS HERE BY STELLARATOR SYMMETRY
!
            IF (.not.lasym .and. (iu.ne.1 .and. iu.ne.ntheta2)) THEN
              lkr = lkr + 1
              nedge = nedge+1
              rb(lkr)   = rb(lk)
              phib(lkr) =-phib(lk)
              zb(lkr)   =-zb(lk)
              bredge(lkr) = -bredge(lk)
              bpedge(lkr) =  bpedge(lk)
              bzedge(lkr) =  bzedge(lk)
            ENDIF
         END DO
      END DO

      IF (ledge_dump) THEN
        WRITE(NEDGE0,*) 'INPUT FILE = ',arg1
        WRITE(NEDGE0,*) 'NEDGE = ',nedge
        WRITE(NEDGE0,*) 'RB = ',  (rb(i), i=1,nedge)
        WRITE(NEDGE0,*) 'PHIB = ',(phib(i), i=1,nedge)
        WRITE(NEDGE0,*) 'ZB = ',  (zb(i), i=1,nedge)
        WRITE(NEDGE0,*) 'BREDGE = ', (bredge(i), i=1,nedge)
        WRITE(NEDGE0,*) 'BPEDGE = ', (bpedge(i), i=1,nedge)
        WRITE(NEDGE0,*) 'BZEDGE = ', (bzedge(i), i=1,nedge)
      END IF

!
!     WRITE OUT (TO THREED1 FILE) VACUUM INFORMATION
!

      IF (.not.lfreeb) THEN
         DEALLOCATE (rb, phib, zb, stat=l)
         RETURN
      END IF

      DO iprint = 1, 2
         IF (iprint .eq. 1) WRITE (nthreed, 750)
         IF (iprint .eq. 2) WRITE (nthreed, 760)
         nzskip = 1 + nzeta/6
         DO l = 1, nzeta, nzskip
            zeta = (360.0_dp*(l - 1))/nzeta
            IF (iprint .eq. 1) THEN
               DO k = 1, ntheta2
                  lk = l + nzeta*(k - 1)
                  WRITE (nthreed, 770) zeta, rb(lk),
     1            zb(lk), (bsqsav(lk,n),n=1,3), bsqvac(lk)
               END DO
            ELSE
               DO k = 1, ntheta2
                  lk = l + nzeta*(k - 1)
                  WRITE (nthreed, 780) zeta, rb(lk), zb(lk),
     1            bredge(lk), bpedge(lk), bzedge(lk),
     2            brv(lk), bphiv(lk), bzv(lk)
               END DO
            ENDIF
         END DO
      END DO

      DEALLOCATE (rb, phib, zb, bredge, bpedge, bzedge, stat=l)

      IF (lasym) THEN
         WRITE (nthreed, 900)
         DO mn = 1, mnmax
            potsin = 0; potcos = 0
            DO mn0 = 1, mnpd
               IF ( (NINT(xnpot(mn0)).eq.NINT(xn(mn))) .and.
     1              (NINT(xmpot(mn0)).eq.NINT(xm(mn))) ) THEN
                  potsin = potvac(mn0)
                  potcos = potvac(mn0+mnpd)
                  EXIT
               END IF
            END DO
         WRITE (nthreed, 910) NINT(xn(mn)/nfp), NINT(xm(mn)), rmnc(mn),
     1      zmns(mn), rmns(mn), zmnc(mn), potsin, potcos, 
     2      bmodmn(mn), bmodmn1(mn)
         END DO

      ELSE
         WRITE (nthreed, 800)
         DO mn = 1, mnmax
            potsin = 0
            DO mn0 = 1, mnpd
               IF ( (NINT(xnpot(mn0)).eq.NINT(xn(mn))) .and.
     1              (NINT(xmpot(mn0)).eq.NINT(xm(mn))) ) THEN
                  potsin = potvac(mn0)
                  EXIT
               END IF
            END DO
            WRITE (nthreed, 810) NINT(xn(mn)/nfp), NINT(xm(mn)), 
     1      rmnc(mn), zmns(mn), potsin, bmodmn(mn), bmodmn1(mn)
         END DO
      END IF

      WRITE (nthreed, *)

  750 FORMAT(/,3x,'NF*PHI',7x,' Rb ',8x,' Zb ',6x,'BSQMHDI',5x,'BSQVACI'
     1   ,5x,'BSQMHDF',5x,'BSQVACF',/)
  760 FORMAT(/,3x,'NF*PHI',7x,' Rb ',8x,' Zb ',6x,'BR',8x,'BPHI',6x,'BZ'
     1   ,8x,'BRv',7x,'BPHIv',5x,'BZv',/)
  770 FORMAT(1p,e10.2,6e12.4)
  780 FORMAT(1p,e10.2,2e12.4,6e10.2)
  790 FORMAT(i5,/,(1p,3e12.4))
  800 FORMAT(//,3x,'nb',2x,'mb',6x,'rbc',9x,'zbs',6x,'vacpot_s',
     1          2x, '|B|_c(s=.5)',1x,'|B|_c(s=1.)'/)
  810 FORMAT(i5,i4,1p,7e12.4)
  900 FORMAT(//,3x,'nb',2x,'mb',6x,'rbc',9x,'zbs',9x,'rbs',9x,'zbc',
     1          6x,'vacpot_s',4x,'vacpot_c',2x,'|B|_c(s=.5)',
     2          1x,'|B|_c(s=1.)'/)
  910 FORMAT(i5,i4,1p,10e12.4)

      END SUBROUTINE freeb_data
