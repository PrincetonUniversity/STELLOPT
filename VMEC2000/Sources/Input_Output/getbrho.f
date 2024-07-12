      SUBROUTINE getbsubs (bsubsmn, frho, bsupu, bsupv, mmax, 
     1                     nmax, info)
      USE stel_kinds
      USE vmec_input, ONLY: nfp, nzeta, lasym
      USE vmec_dim, ONLY: ntheta1, ntheta2, ntheta3
      USE vmec_persistent, ONLY: cosmu, sinmu, cosnv, sinnv
      USE vmec_main, ONLY: r0scale
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: mmax, nmax
      INTEGER, INTENT(out) :: info
      REAL(rprec), INTENT(out) :: bsubsmn(0:mmax, -nmax:nmax, 0:1)
      REAL(rprec), DIMENSION(nzeta, ntheta3), INTENT(in) ::
     1    bsupu, bsupv, frho
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: p5 = 0.5_dp, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, j, m, n, nmax1, itotal, ijtot, mntot
      REAL(rprec) :: ccmn, ssmn, csmn, scmn, dm, dn, termsc, termcs,
     `               termcc, termss, amn
      REAL(rprec), ALLOCATABLE :: amatrix(:,:), save_matrix(:,:),
     1   brhs(:)
      LOGICAL :: lpior0
      EXTERNAL solver
C-----------------------------------------------
!
!     Solves the radial force balance B dot bsubs = Fs for bsubs in real space using collocation
!     Here, Fs = frho(mn) is the Fourier transform SQRT(g)*F (part of radial force
!     balance sans the B dot bsubs term)
!
!     Storage layout: bsubsmn(0:mmax, 0:nmax,0) :: coefficient of sin(mu)cos(nv)
!                     bsubsmn(0:mmax,-1:nmax,0) :: coefficient of cos(mu)sin(nv)
!     for lasym = T
!                     bsubsmn(0:mmax, 0:nmax,1) :: coefficient of cos(mu)cos(nv)
!                     bsubsmn(0:mmax,-1:nmax,1) :: coefficient of sin(mu)sin(nv)
!
!     where 0<=m<=mmax and 0<=n<=nmax 
!
      info = -3
      IF ((mmax+1 .ne. ntheta2) .or. (nmax .ne. nzeta/2)) RETURN

      nmax1 = MAX(0,nmax-1)
      itotal = ntheta3*nzeta
      IF (.not.lasym) itotal = itotal - 2*nmax1
      ALLOCATE (amatrix(itotal, itotal),
     1      brhs(itotal), save_matrix(itotal, itotal), stat=m)
      IF (m .ne. 0) STOP 'Allocation error in getbsubs'

      amatrix = 0

!
!     bsubs = BSC(M,N)*SIN(MU)COS(NV) + BCS(M,N)*COS(MU)SIN(NV)
!          + BCC(M,N)*COS(MU)COS(NV) + BSS(M,N)*SIN(MU)SIN(NV)   (LASYM=T ONLY)
!

      ijtot = 0
      brhs = 0

      DO i = 1, ntheta3
         DO j = 1, nzeta
!           IGNORE u=0,pi POINTS FOR v > pi: REFLECTIONAL SYMMETRY
            lpior0 = ((i.eq.1 .or. i.eq.ntheta2) .and. (j.gt.nzeta/2+1))
            IF (lpior0 .and. .not. lasym) CYCLE
            ijtot = ijtot + 1
            brhs(ijtot) = frho(j,i)
            mntot = 0
            DO m = 0, mmax
               DO n = 0, nmax
                  IF (mntot .ge. itotal) EXIT
                  IF (m.eq.0 .and. n.eq.0 .and. lasym) CYCLE
                  mntot = mntot+1
                  ccmn = cosmu(i,m)*cosnv(j,n)
                  ssmn = sinmu(i,m)*sinnv(j,n)
                  dm = m * bsupu(j,i)
                  dn = n * bsupv(j,i) * nfp
                  termsc = dm*ccmn - dn*ssmn
                  termcs =-dm*ssmn + dn*ccmn
                  IF (n.eq.0 .or. n.eq.nmax) THEN
                     IF (m .gt. 0) THEN
                        amatrix(ijtot,mntot) = termsc     !!ONLY bsc != 0 for n=0, nmax1
                     ELSE IF (n .eq. 0) THEN
                        amatrix(ijtot,mntot) = bsupv(j,i) !!pedestal for m=0,n=0 mode, which should = 0
                     ELSE
                        amatrix(ijtot,mntot) = termcs     !!bcs(m=0,n=nmax)
                     END IF
                  ELSE IF (m.eq.0 .or. m.eq.mmax) THEN
                     amatrix(ijtot,mntot) = termcs        !!ONLY bcs != 0 for m=0,mmax
                  ELSE
                     amatrix(ijtot,mntot) = termsc
                     mntot = mntot+1
                     amatrix(ijtot,mntot) = termcs
                  END IF

                  IF (.not.lasym) CYCLE
                  IF (m.eq.0 .and. (n.eq.0 .or. n.eq.nmax)) CYCLE

                  IF (mntot .ge. itotal) EXIT
                  mntot = mntot+1
                  csmn = cosmu(i,m)*sinnv(j,n)
                  scmn = sinmu(i,m)*cosnv(j,n)
                  termcc =-dm*scmn - dn*csmn
                  termss = dm*csmn + dn*scmn

                  IF ((n.eq.0 .or. n.eq.nmax) .or.
     1                (m.eq.0 .or. m.eq.mmax)) THEN
                      amatrix(ijtot,mntot) = termcc    !!ONLY bcc != 0 for m=0 or mmax
                  ELSE
                     amatrix(ijtot,mntot) = termcc
                     mntot = mntot+1
                     amatrix(ijtot,mntot) = termss
                  END IF

               END DO
            END DO
         END DO
      END DO

      save_matrix = amatrix

      info = -1
      IF (ijtot .ne. itotal .or. mntot .ne. itotal) THEN
         PRINT *,' itotal = ', itotal,' ijtot = ', ijtot,
     1   ' mntot = ', mntot
         PRINT *,' ntheta3: ',ntheta3,' nzeta: ', nzeta,
     1           ' mnyq: ', mmax,' nnyq: ', nmax
         GOTO 200
      ELSE
         CALL solver (amatrix, brhs, itotal, 1, info)
         IF (info .ne. 0) GOTO 200
      END IF

!
!     CHECK SOLUTION FROM SOLVER
!
!      GOTO 100
      ijtot = 0
      DO i = 1, ntheta3
         DO j = 1, nzeta
            lpior0 = ((i.eq.1 .or. i.eq.ntheta2) .and. (j.gt.nzeta/2+1))
            IF (lpior0 .and. .not.lasym) CYCLE
            ijtot = ijtot + 1
            amn = SUM(save_matrix(ijtot,:)*brhs(:))
            IF (ABS(amn) .lt. 1.E-12_dp) CYCLE
            IF (ABS(frho(j,i) - amn) .gt. 1.e-8_dp*ABS(amn)) THEN
               PRINT 50,'In GETbsubs, i = ',i,' j = ',j,
     1         ' Original force = ', frho(j,i),' Final force = ', amn
            END IF
        END DO
      END DO
  50  FORMAT(a,i5,a,i5,a,1p,e10.3,a,1p,e10.3)
 100  CONTINUE
!
!     CONVERT BACK TO BS*SIN(MU - NV) REPRESENTATION 
!     AND (FOR lasym) BC*COS(MU - NV)
!
      mntot = 0
      bsubsmn = 0
      DO m = 0, mmax
         DO n = 0, nmax
            IF (mntot .ge. itotal) EXIT
            IF (m.eq.0 .and. n.eq.0 .and. lasym) CYCLE
            mntot = mntot+1
            IF (n.eq.0 .or. n.eq.nmax) THEN
               IF (m .gt. 0) THEN
                  bsubsmn(m,n,0) = brhs(mntot)
               ELSE IF (n .eq. 0) THEN
                  bsubsmn(m,n,0) = brhs(mntot)
               ELSE
                  bsubsmn(m,-n,0) = brhs(mntot)
               END IF
            ELSE IF (m.eq.0 .or. m.eq.mmax) THEN
               bsubsmn(m,-n,0) = brhs(mntot)
            ELSE
               bsubsmn(m,n,0) = brhs(mntot)
               mntot = mntot+1
               bsubsmn(m,-n,0) = brhs(mntot)
            END IF

            IF (.not.lasym) CYCLE
            IF (m.eq.0 .and. (n.eq.0 .or. n.eq.nmax)) CYCLE
            IF (mntot .ge. itotal) EXIT
            mntot = mntot+1

            IF ((n.eq.0 .or. n.eq.nmax) .or.
     1          (m.eq.0 .or. m.eq.mmax)) THEN
               bsubsmn(m,n,1) = brhs(mntot)
            ELSE
               bsubsmn(m,n,1) = brhs(mntot)
               mntot = mntot+1
               bsubsmn(m,-n,1)= brhs(mntot)
            END IF

         END DO
      END DO

      IF (mntot .ne. ijtot) info = -2

 200  CONTINUE

      DEALLOCATE (amatrix, save_matrix, brhs)

      END SUBROUTINE getbsubs
