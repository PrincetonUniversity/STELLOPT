      SUBROUTINE eval_saddle_rz (v0, npts, rpts, zpts)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE stel_constants
      USE boundary, ONLY:nfp
      USE saddle_coils
      USE saddle_surface
      USE Vcoilpts
      USE Vwire
      USE coils
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, k, n, npt, npoints, ierr
      INTEGER, DIMENSION(ncdim) :: npts
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: uc, us
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: vc, vs
      REAL(rprec) :: u1, u2, v0, v1, v2, ds
      REAL(rprec) :: r1, s1, d1
      REAL(rprec) :: r2, s2, d2
      REAL(rprec), DIMENSION(2,3,5) :: x1,  x2, y1,  y2, z1, z2
      REAL(rprec), DIMENSION(2) :: c
      REAL(rprec), DIMENSION(ncdim, 5*ncdim) :: rpts, zpts
      REAL(rprec) :: alf, bet
      LOGICAL :: la, lb, lc
!-----------------------------------------------
      ALLOCATE (uc(0:nsad_u), us(0:nsad_u+4),
     1          vc(0:nsad_v), vs(0:nsad_v+4))

      alf = 1
      bet = 0
      la = .true.
      lb = .false.
      lc = lctrlpt
      npoints = 2*nwire

      ds = one/npoints
!     loop over unique coils
      DO n=1, nsmid
         DO k=0, nsad_u
            uc(k) = saddle(n)%u_c(k)
            us(k) = saddle(n)%u_s(k)
         END DO
         DO k=0, nsad_v
            vc(k) = saddle(n)%v_c(k)
            vs(k) = saddle(n)%v_s(k)
         END DO
         npt = 0
!        loop over number of poloidal points
         DO i=1, npoints
            s1 = (i - 1)*ds
            CALL coil_curve( nsad_u, uc, us, nsad_v, vc, vs, la,
     1         lb, lc, lspline, nfp, numsurf_sad, m_sad, n_sad,
     2         rmn_sad, zmn_sad, nfils, s1, alf, bet, deln, delt,
     3         u1, v1, r1, c, x1, y1, z1, ierr )
            s2 = i*ds
            CALL coil_curve( nsad_u, uc, us, nsad_v, vc, vs, la,
     1         lb, lc, lspline, nfp, numsurf_sad, m_sad, n_sad,
     2         rmn_sad, zmn_sad, nfils, s2, alf, bet, deln, delt,
     3         u2, v2, r2, c, x2, y2, z2, ierr )
            d1 = v1 - v0
            d2 = v2 - v0
            IF (d1*d2 .lt. 0) THEN
               npt = npt + 1
               rpts(n,npt) = (r1 + r2)/2
               zpts(n,npt) = (z1(1,1,1) + z2(1,1,1))/2
               IF (nfils .ge. 3) THEN
                  npt = npt + 1
                  r1 = SQRT(x1(1,1,2)**2 + y1(1,1,2)**2)
                  r2 = SQRT(x2(1,1,2)**2 + y2(1,1,2)**2)
                  rpts(n,npt) = (r1 + r2)/2
                  zpts(n,npt) = (z1(1,1,2) + z2(1,1,2))/2
                  npt = npt + 1
                  r1 = SQRT(x1(1,1,3)**2 + y1(1,1,3)**2)
                  r2 = SQRT(x2(1,1,3)**2 + y2(1,1,3)**2)
                  rpts(n,npt) = (r1 + r2)/2
                  zpts(n,npt) = (z1(1,1,3) + z2(1,1,3))/2
               END IF
               IF (nfils .eq. 5) THEN
                  npt = npt + 1
                  r1 = SQRT(x1(1,1,4)**2 + y1(1,1,4)**2)
                  r2 = SQRT(x2(1,1,4)**2 + y2(1,1,4)**2)
                  rpts(n,npt) = (r1 + r2)/2
                  zpts(n,npt) = (z1(1,1,4) + z2(1,1,4))/2
                  npt = npt + 1
                  r1 = SQRT(x1(1,1,5)**2 + y1(1,1,5)**2)
                  r2 = SQRT(x2(1,1,5)**2 + y2(1,1,5)**2)
                  rpts(n,npt) = (r1 + r2)/2
                  zpts(n,npt) = (z1(1,1,5) + z2(1,1,5))/2
               END IF
            END IF
         END DO
         npts(n) = npt
      END DO

      DEALLOCATE (uc, us, vc, vs)

      END SUBROUTINE eval_saddle_rz
