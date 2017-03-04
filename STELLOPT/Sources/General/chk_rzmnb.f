      SUBROUTINE chk_rzmnb(xcb_opt, rmnc_a, zmns_a, xm, xn, iflag)
      USE optim
      USE vmec_input, ONLY : rbc, zbs
      USE vparams, ONLY: zero
      USE mpi_params                                                     !MPI
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(out) :: iflag
      REAL(rprec), DIMENSION(*) :: xcb_opt
      REAL(rprec), DIMENSION(mnmax_opt), INTENT(in) ::
     1   rmnc_a, zmns_a, xm, xn
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: eps = 1.0e-10_dp
      CHARACTER(LEN=*), PARAMETER :: error_message =
     1   'Boundary array inconsistency in CHK_RZMNB'
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: nb, mb, ikt, mn, j1
      REAL(rprec) :: rzero, tol, compit
      LOGICAL :: lerror
C-----------------------------------------------
      iflag = 0;  lerror = .false.
      IF (irm0_bdy.eq.0 .or. lniter1) RETURN
      rzero = eps

      zbs(0,0) = zero

!
!     NOTE: THIS ASSUMES THAT m=0 R,Z COMPONENTS ARE THE FIRST ELEMENTS IN XCB_OPT ARRAY
!
      DO ikt = 1, irm0_bdy
         rbc(nrz0_opt(ikt),0) = xcb_opt(ikt)
      END DO
      DO ikt = 1, izm0_bdy
         j1 = ikt+irm0_bdy
         zbs(nrz0_opt(j1),0) = xcb_opt(j1)
      END DO
      DO ikt = 1, irho_bdy
         j1 = ikt + irm0_bdy + izm0_bdy
         rhobc(nbrho_opt(ikt),mbrho_opt(ikt)) = xcb_opt(j1)
      END DO

      CALL unique_boundary(rbc, zbs, rhobc, mpol1d, ntord, 
     1                     mpol1_opt, ntor_opt, mrho1_opt)

      tol = SQRT(rzero)
      compit = rzero + ABS(rbc(0,0)) * tol

      DO mn = 1,mnmax_opt
         nb = NINT(xn(mn))/nfp_opt
         mb = NINT(xm(mn))
         IF(ABS(rbc(nb,mb) - rmnc_a(mn)).gt.compit) THEN
            IF (myid .eq. master)                                        !MPI
     1         WRITE(*,57)nb,mb,rbc(nb,mb),rmnc_a(mn)
            lerror = .true.
         END IF
         IF(ABS(zbs(nb,mb) - zmns_a(mn)).gt.compit) THEN
           IF (myid .eq. master)                                         !MPI
     1        WRITE(*,58)nb,mb,zbs(nb,mb),zmns_a(mn)
           lerror = .true.
         END IF
      ENDDO

 57   FORMAT(' nb = ',i5,' mb = ',i5,' rmnc diff: ',2e15.5)
 58   FORMAT(' nb = ',i5,' mb = ',i5,' zmns diff: ',2e15.5)

      IF (lerror) THEN
         PRINT *, TRIM(error_message),' in process ', myid
         WRITE (iunit_opt_local,*) TRIM(error_message),
     1      ' in process ', myid
         iflag = -15
      END IF

      END SUBROUTINE chk_rzmnb
