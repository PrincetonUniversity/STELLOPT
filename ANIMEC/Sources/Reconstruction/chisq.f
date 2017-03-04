      SUBROUTINE chisq(amat_i, amat_p, data, idata, isize, itotal)
      USE vmec_main
      USE vsvd
      USE vspline
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: itotal
      INTEGER, DIMENSION(*) :: idata, isize
      REAL(rprec), DIMENSION(isnodes,*) :: amat_i
      REAL(rprec), DIMENSION(ipnodes,*) :: amat_p
      REAL(rprec), DIMENSION(*) :: data
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, n, i1, ispec, is(5), ip(5), j2
      REAL(rprec) :: delsq, delp, dels
      CHARACTER(LEN=130) label
C-----------------------------------------------
!
!       COMPUTES CHI**2 FROM DIFFERENT SUBROUTINES AT VARIOUS TIME-STEPS
!       WRITTEN BY D.K. LEE (3/93)
!

         chisqerr(:jchix) = zero

         DO i = 1, itotal
            delsq = (SUM(ystark(:isnodes)*amat_i(:isnodes,i)) +
     1               SUM(ythom(:ipnodes)*amat_p(:ipnodes,i))-data(i))**2
            IF (i.ge.idata(ithom0) .and. i<idata(ithom0)+isize(ithom0))
     1         THEN
               chisqerr(ithom0) = chisqerr(ithom0) + delsq
            ELSE IF (i.ge.idata(istark0) .and. i<idata(istark0)
     1            +isize(istark0)) THEN
               i1 = i - idata(istark0) + 1
               IF (i1 .eq. islope) THEN
                  chisqerr(islope0) = delsq
               ELSE IF (i1 .eq. icurrout) THEN
                  chisqerr(icurr0) = delsq
               ELSE
                  chisqerr(istark0) = chisqerr(istark0) + delsq
               ENDIF
            ELSE IF (i.ge.idata(idiam0) .and. i<idata(idiam0)
     1            +isize(idiam0)) THEN
               chisqerr(idiam0) = chisqerr(idiam0) + delsq
            ELSE IF (i.ge.idata(iflxs0) .and. i<idata(iflxs0)
     1            +isize(iflxs0)) THEN
               chisqerr(iflxs0) = chisqerr(iflxs0) + delsq
            ELSE IF (i.ge.idata(ibrzfld) .and. i<idata(ibrzfld)
     1            +isize(ibrzfld)) THEN
               chisqerr(ibrzfld) = chisqerr(ibrzfld) + delsq
            ENDIF
         END DO

!
         errsvd = SUM(chisqerr(:jchix))
         IF (.not.lpprof) errsvd = errsvd - chisqerr(ithom0)

      IF (iequi.ne.1 .or. .not.lrecon) RETURN
         IF (.not.lpprof) THEN
            WRITE (nthreed, 15)
         ELSE
            WRITE (nthreed, 10)
         ENDIF
         IF (lpprof) THEN
            DO n = 1, nchistp
               WRITE (nthreed, 20) nchi2(n), chi2(ithom0,n),
     1            chi2(istark0,n), chi2(icurr0,n), chi2(idiam0,n),
     2            chi2(iflxs0,n), chi2(ibrzfld,n), chi2(jchix1,n)
            END DO
         ELSE
            DO n = 1, nchistp
               WRITE (nthreed, 20) nchi2(n), chi2(istark0,n),
     1            chi2(icurr0,n), chi2(idiam0,n), chi2(iflxs0,n),
     2            chi2(ibrzfld,n), chi2(jchix1,n)
            END DO
         ENDIF

!
!       PRINT OUT MATRIX ELEMENTS (5 EACH FOR PRESSURE, IOTA)
!
         WRITE (nthreed, 200)
         delp = (ipnodes - 1)/4.
         dels = (isnodes - 1)/4.
         ip(1) = 1
         is(1) = 1
         ip(2:4) = ip(1) + INT((((/(j2,j2=2,4)/)) - 1)*delp)
         is(2:4) = is(1) + INT((((/(j2,j2=2,4)/)) - 1)*dels)
         ip(5) = ipnodes
         is(5) = isnodes
         WRITE (label, 210) ip(1), ip(2), ip(3), ip(4), ip(5), is(1),
     1      is(2), is(3), is(4), is(5)
         WRITE (nthreed, 220) label
         ispec = 0
         DO i = 1, itotal
            IF (i.ge.idata(ithom0) .and.
     1         i.lt.idata(ithom0)+isize(ithom0)) THEN
               i1 = i - idata(ithom0) + 1
               CALL printmatrix (amat_p(1,i), amat_i(1,i), data(i), i,
     1            i1, ip, is, '   PRES  (')
            ELSE IF (i.ge.idata(istark0) .and. i<idata(istark0)
     1            +isize(istark0)) THEN
               i1 = i - idata(istark0) + 1
               IF (i1 .eq. islope) THEN
                  ispec = ispec - 1
                  CALL printmatrix (amat_p(1,i), amat_i(1,i), data(i),
     1               i, 1, ip, is, '  IOTA0  (')
               ELSE IF (i1 .eq. icurrout) THEN
                  ispec = ispec - 1
                  CALL printmatrix (amat_p(1,i), amat_i(1,i), data(i),
     1               i, 1, ip, is, ' CURRENT (')
               ELSE
                  CALL printmatrix (amat_p(1,i), amat_i(1,i), data(i),
     1               i, i1 + ispec, ip, is, '   MSE   (')
               ENDIF
            ELSE IF (i.ge.idata(idiam0) .and. i<idata(idiam0)
     1            +isize(idiam0)) THEN
               i1 = i - idata(idiam0) + 1
               CALL printmatrix (amat_p(1,i), amat_i(1,i), data(i), i,
     1            i1, ip, is, ' DIAMAG  (')
            ELSE IF (i.ge.idata(iflxs0) .and. i<idata(iflxs0)
     1            +isize(iflxs0)) THEN
               i1 = i - idata(iflxs0) + 1
               CALL printmatrix (amat_p(1,i), amat_i(1,i), data(i), i,
     1            i1, ip, is, ' FLUXES  (')
            ELSE IF (i.ge.idata(ibrzfld) .and. i<idata(ibrzfld)
     1            +isize(ibrzfld)) THEN
               i1 = i - idata(ibrzfld) + 1
               CALL printmatrix (amat_p(1,i), amat_i(1,i), data(i), i,
     1            i1, ip, is, '  BR-BZ  (')
            ENDIF
         END DO
   20 FORMAT(i6,1p,8e12.4)
   10 FORMAT(/,30x,'ABSOLUTE CHI-SQ ERROR BY DATA TYPE'/,30x,
     1   '(NOT NORMED BY NUMBER DATA POINTS)'/,20x,
     2   'NOTE: STARK CHISQ MAY BE EVALUATED AT REDISTRIBUTED KNOTS'/,
     3   '  ITER   Thomscat      Stark',5x,'Current',5x,'Diamag.',5x,
     4   'Saddle',6x,' B-Loops',6x,'TOTAL',/,1x,5('-'),7(2x,10('-')))
   15 FORMAT(/,30x,'ABSOLUTE CHI-SQ ERROR BY DATA TYPE'/,30x,
     1   '(NOT NORMED BY NUMBER DATA POINTS)'/,20x,
     2   'NOTE: STARK CHISQ MAY BE EVALUATED AT REDISTRIBUTED KNOTS'/,
     3   '  ITER      Stark',5x,'Current',5x,'Diamag.',5x,'Saddle',6x,
     4   ' B-Loops',6x,'TOTAL',/,1x,5('-'),7(2x,10('-')))
  200 FORMAT(//,38x,'SPLINE MATRIX ELEMENTS BY DATA TYPE'/,30x,
     1   ' AI(i,j)*iota(j) + AP(i,j)*[mu0*pres(j)] = DATA(i)'/,30x,
     2   ' NOTE: DATA(I) IS THE RAW DATA NORMED TO SIGMA(I)'//)
  210 FORMAT('   I   TYPE         DATA(I)','  AP(I,',i2,')  AP(I,',i2,
     1   ')  AP(I,',i2,')  AP(I,',i2,')','  AP(I,',i2,')  AI(I,',i2,
     2   ')  AI(I,',i2,')  AI(I,',i2,')','  AI(I,',i2,')  AI(I,',i2,')')
  220 FORMAT(a,/,3x,'-',3x,4('-'),9x,7('-'),10(2x,8('-')))

      END SUBROUTINE chisq
