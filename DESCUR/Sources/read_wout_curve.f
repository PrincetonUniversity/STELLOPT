      SUBROUTINE read_wout_curve (wout_file, ierr,
     1           nphi20, rbdy_3d, zbdy_3d)
      USE Vname0
      USE Vname1, ONLY : mpol_in => mpol, nfp_in => nfp,
     1    nphi2_in => nphi2
      USE read_wout_mod
      IMPLICIT NONE

      INTEGER, INTENT(out) :: ierr
      INTEGER :: mn, m, n, nphi20
      REAL(rprec) :: rbdy_3d(0:mu-1,-nphi20:nphi20),
     1               zbdy_3d(0:mu-1,-nphi20:nphi20)
      CHARACTER*(*) :: wout_file

      ierr = -1
      IF (INDEX(wout_file,'wout') .eq. 1)
     1    CALL read_wout_file(wout_file(6:), ierr)
      IF (ierr .ne. 0) RETURN

      mpol_in  = mpol
      nfp_in   = nfp
      nphi2_in = ntor

      IF (mpol .gt. mu) STOP 'mpol-input > mu'
      IF (nfp .le. 0) nfp = 1

      DO mn = 1, mnmax
         m = NINT(xm(mn))
         n = NINT(xn(mn))/nfp
         rbdy_3d(m,n) = rmnc(mn,ns)
         zbdy_3d(m,n) = zmns(mn,ns)
      END DO

      END SUBROUTINE read_wout_curve
