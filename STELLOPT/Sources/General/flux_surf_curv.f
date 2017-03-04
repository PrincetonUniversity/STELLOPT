      SUBROUTINE flux_surf_curv(curv_kur, rmnc_a, zmns_a, xm, xn)
      USE optim
      USE safe_open_mod
      USE vparams, ONLY: zero, one
      USE mpi_params                                                     !MPI
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(4) :: curv_kur
      REAL(rprec), DIMENSION(mnmax_opt) :: rmnc_a, zmns_a, xm, xn
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: four = 4.0_dp
      LOGICAL, PARAMETER :: lprint = .false.         !!Set .true. to get debugging info
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: lzth, l, ith, i, iunit1, iunit2, iunit3, istat
      REAL(rprec), DIMENSION(100) :: r_major, z_elev
      REAL(rprec), DIMENSION(100,4) :: curv
      REAL(rprec), DIMENSION(4) ::
     1   avg_curv, curv_dev2, curv_dev4, curv_max
      REAL(rprec) :: xpi, zeta_1, theta_1, arg, s_all, sa2,
     1   sb2, s_x2, s_y2, delt_1, ratio
      CHARACTER(LEN=15) :: filename
C-----------------------------------------------
!
!     Open up several diagnostic data files.
!
      IF (lprint) THEN
         iunit1 = 34
         iunit2 = iunit1+1
         iunit3 = iunit2+1

         CALL safe_open(iunit1, istat, 'curv_dev', 'replace',
     1      'formatted')
         CALL safe_open(iunit2, istat, 'outer_surf', 'replace',
     1      'formatted')
      ENDIF

      xpi = four*ATAN(one)
      curv_kur = zero
!
!     check curvature at phi=0,90,180,270 degree planes (lzth = 1,2,3,4):
!
      DO 100 lzth = 1, 4
         avg_curv(lzth) = zero
         IF (lprint) THEN
            IF (lzth .eq. 1) filename = 'test_out_0'
            IF (lzth .eq. 2) filename = 'test_out_90'
            IF (lzth .eq. 3) filename = 'test_out_180'
            IF (lzth .eq. 4) filename = 'test_out_270'
            CALL safe_open(iunit3, istat, filename, 'replace',
     1         'formatted')
            IF (istat.ne.0 .and. myid.eq.master)                         !MPI
     1      PRINT *, ' Error opening ',TRIM(filename),
     2               ' in FLUX_SURF_CURV'
         END IF
         zeta_1 = ((lzth - 1)*xpi)/REAL((2*nfp_opt),rprec)
         IF (lzth.eq.1 .and. lprint) THEN
            DO l = 1, mnmax_opt
               WRITE (iunit2, *) xm(l), xn(l), rmnc_a(l), zmns_a(l)
            END DO
            CLOSE(iunit2)
         ENDIF
c
c        Construct REAL space R, Z coordinates of outer flux surface:
c
         DO ith = 1,100
            theta_1 = (2.0_dp*xpi*ith)/100.0_dp
            R_major(ith) = zero
            z_elev(ith) = zero
            DO l=1,mnmax_opt
               arg = xm(l)*theta_1 - xn(l)*zeta_1
               R_major(ith) = R_major(ith) + rmnc_a(l)*COS(arg)
               z_elev(ith) = z_elev(ith) + zmns_a(l)*SIN(arg)
            END DO
         END DO
c
c        Calculate local curvature around flux surface by finite
c        differencing the local unit TANgent vectors:
c
         DO i=2,99
            s_all = zero
            sa2 = SQRT((R_major(i)-R_major(i-1))**2 +
     1            (z_elev(i)-z_elev(i-1))**2)
            sb2 = SQRT((R_major(i+1)-R_major(i))**2 +
     1            (z_elev(i+1)-z_elev(i))**2)
            s_x2 =  ((R_major(i) - R_major(i-1))/sa2
     1          -   (R_major(i+1) - R_major(i))/sb2)**2
            s_y2 =  ((z_elev(i) - z_elev(i-1))/sa2
     1          - (z_elev(i+1) - z_elev(i))/sb2)**2
            s_all = SQRT(s_x2 + s_y2)
            delt_1 = .5_dp*SQRT((R_major(i-1) - R_major(i+1))**2
     1             + (z_elev(i-1) - z_elev(i+1))**2)

            curv(i,lzth) = zero
            IF (delt_1 .ne. zero) curv(i,lzth) = s_all/delt_1
            avg_curv(lzth) = avg_curv(lzth) + curv(i,lzth)/98.0_dp
            IF (lprint) WRITE (iunit3,998) R_major(i),
     1         z_elev(i),curv(i,lzth)
         END DO
         IF (lprint) CLOSE(iunit3)
 100  CONTINUE

      DO lzth = 1,4
         curv_dev2(lzth) = zero
         curv_dev4(lzth) = zero
         curv_max(lzth) = -1.e30_dp
c        Calculate average curvature, variance of curvature, and kurtosis
c        (4th moment/2nd moment) around flux surface.
         DO i=2,99
            curv_max(lzth) = MAX(curv_max(lzth),curv(i,lzth))
            curv_dev2(lzth) = curv_dev2(lzth) + ((curv(i,lzth)
     1                     - avg_curv(lzth))**2)/98.0_dp
            curv_dev4(lzth) = curv_dev4(lzth) + ((curv(i,lzth)
     1                     - avg_curv(lzth))**4)/98.0_dp
         END DO
         curv_kur(lzth) = curv_dev4(lzth)/(curv_dev2(lzth)**2)
         IF(lprint) THEN
            ratio = curv_max(lzth)/avg_curv(lzth)
            WRITE(iunit1,997) lzth, avg_curv(lzth),
     1          curv_dev2(lzth), ratio, curv_kur(lzth)
         END IF
      END DO

  998 FORMAT(3(1x,e15.6))
  997 FORMAT(1x,i2,4(1x,e15.6))
      IF (lprint) CLOSE(iunit1)

      END SUBROUTINE flux_surf_curv
