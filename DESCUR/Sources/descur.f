      PROGRAM descur
c       THIS IS PROGRAM DESCUR - WHICH USES A STEEPEST DESCENT
c       ALGORITHM TO FIND A LEAST SQUARES APPROXIMATION TO AN
c       ARBITRARY 3-D SPACE CURVE. ANGLE CONSTRAINTS, BASED ON A
c       MINIMUM SPECTRAL WIDTH CRITERION, ARE APPLIED.
c       THE CONSTRAINT IS SATISFIED BY TANGENTIAL VARIATIONS
c       ALONG THE CURVE.
c
c       THE MAIN PROGRAM SETS UP THE INITIAL POINTS TO BE FIT AND
c       THEN CALLS THE SUBROUTINE SCRUNCH, WHICH DOES THE ACTUAL
c       CURVE-FITTING.
c
c***********************************************************************
c       REVISION 1:  January 26, 1989
c                    Compute Fit to individual toroidal planes separately
c                    and THEN Fourier decompose in phi, using optimized
c                    angle representation
c
c       REVISION 2:  July, 1995
c                    Generalized for up-down asymmetric CASE
c
c       REVISION 3:  July, 1997
c                    Converted to Fortran-90
c                    Eliminated Lagrange multiplier constraint
c                    by using an explicit representation
c
c
c       PARAMETER VALUE INTERPRETATIONS:
c
c       MPOL:   ACTUAL NO. OF POLOIDAL MODES USED IN FIT FOR R and Z
c       MRHO:   ACTUAL NO. OF POLOIDAL MODES IN FIT FOR QUASI-POLAR RHO
c       NTHETA: ACTUAL NO. OF THETA (POLOIDAL) POINTS AUSED IN FIT
c       NPHI:   ACTUAL NO. OF TOROIDAL PLANES IN FIT
c
c***********************************************************************
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE Vname0
      USE Vname1
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: nphi20 = 1 + nv/2
      INTEGER, PARAMETER :: mntot = mu*(2*nphi20 + 1)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: niter, nstep, itype, ncount,
     1   i, k, j, m, n, mexp = 4, ns, iread
      REAL(rprec), DIMENSION(nu*nv) :: rin, zin
      REAL(rprec), DIMENSION(mntot) :: rbc, rbs, zbc, zbs
      REAL(rprec), DIMENSION(0:nv) :: rmnaxis, zmnaxis
      REAL(rprec), DIMENSION(0:mu-1,-nphi20:nphi20) ::
     1 rbdy_3d, zbdy_3d
      REAL(rprec), DIMENSION(4) :: rbdy, zbdy, rdshape, zdshape,
     1 rbean, zbean, rbelt, zbelt, rellip, zellip, rsqr, zsqr
      REAL(rprec), DIMENSION(0:3,0:2) :: rheliac, zheliac
      REAL(rprec) :: ftol, zeta, texp, dels, A2, VP,
     1   dumr, dumphi, dumz, arg, PSI, EPS, KAPPA, BDIA, R1, A1, QPSI
!     2  , XI, DELTA_R, rmin, rmax, denom
      CHARACTER :: datafile*80
      CHARACTER*1 :: ch_yn
C-----------------------------------------------

      data rdshape/3.0, 0.991, 0.136, 0./
      data zdshape/0.0, 1.409, -0.118, 0./
      data rbean/3.0, 1.042, 0.502,  - 0.0389/
      data zbean/0.0, 1.339, 0.296,  - 0.0176/
      data rbelt/3.0, 0.453, 0., 0./
      data zbelt/0.0, 0.60, 0., 0.196/
      data rellip/3.0E0, 1.000E0, 0., 0./
      data zellip/0.0, 6.00E0, 0., 0./
      data rheliac/4.115, 0.4753, 0., 0., 0.3225, -0.06208, -0.1358,
     1   -0.04146, 0.0136, 0.0806, -0.0205, 0.0445/
      data zheliac/0., -0.5045, 0., 0., 0.337, 0.0637, 0.1325, -0.04094
     1   , 0.01152, 0.06186, -0.03819,  - 0.02366/
      data rsqr/    3.0, 0.4268, 0., 0.07322/,
     1     zsqr/    0.0, 0.4268, 0., -0.07322/

c**********************************************************************
c       CONTROL DATA THAT CAN BE SPECIFIED BY USER
c**********************************************************************
      data niter/100000/
      data nstep/1000/

      data ftol/1.E-15_dp/

      data datafile/'none'/

      twopi = 8*ATAN(one)

c**********************************************************************
c       CREATE ARRAYS OF CURVE POINTS, RIN(THETA,PHI) , ZIN(THETA,PHI),
c       TO BE FIT TO AN OPTIMIZED FOURIER SERIES BY CALLING "SCRUNCH".
c       NOTE: N IS IN UNITS NORMED TO NFP (NUMBER OF FIELD PERIODS) IN 20 LOOP
c**********************************************************************
c
c     Read in data
c

      WRITE (6, '(a)') ' Enter spectral convergence parameter'
      WRITE (6, '(a)') ' =0 for polar,  =1 for equal arclength'
      WRITE (6, '(a)') ' >1 for smaller spectral width'
      WRITE (6, '(a)', advance='no')
     1      ' =4 corresponds to STELLOPT choice): '
      READ (*,*) texp
      mexp = texp

      HB_Parameter = 1

      WRITE (6, '(a,/,a)', advance='no') 
     1    ' Use (default) VMEC-compatible compression (V)',
     2    ' or Advanced Hirshman-Breslau compression (A): '
      READ (*,'(a)') ch_yn

      IF (ch_yn == 'A' .or. ch_yn == 'a') HB_Parameter = 0
      nfp = 1

      WRITE (6, '(5(a,/))')
     1   ' Select source of curve data to be fit: ',
     2   ' 0 :  Points from file',
     3   ' 1 :  Fourier coeffs from file (wout file allowed)',
     4   ' 2 :  Solove''ev Equilibrium',
     5   ' 3 :  Assorted 2-D shapes'

      READ (*,*) itype

      IF (itype .eq. 3) THEN
         WRITE (*, *)
     1    ' Shape to fit: 1=ellipse; 2=D; 3=bean; 4=belt; 5=square;',
     2    ' 6=heliac)'
         READ (*, *) itype
         IF (itype<1 .or. itype>6) itype = 1
         IF (itype .eq. 6) THEN
            WRITE (*, *)
     1      ' Enter toroidal cross section to plot in degrees: '
            READ (*, *) zeta
            zeta = zeta*twopi/360
         ENDIF
         GOTO 1001

      ELSE IF (itype .eq. 0) THEN
         WRITE(*, '(a)', advance='no')
     1      ' Enter file name containing R, PHI, Z data to fit:'
         READ(*, *) datafile
         IF (datafile .eq. 'none') STOP 'Must enter datafile name!'
         iread = 0
         GOTO 1002
      ELSE IF (itype .eq. 1) THEN
         WRITE(6, '(a)', advance='no')
     1   ' Enter full file name containing Rmn, Zmn',
     2      ' coefficients to convert: '
         READ(*, *) datafile
         IF (datafile .eq. 'none') STOP 'Must enter datafile name!'
         GOTO 1003
      ELSE IF (itype .eq. 2) THEN
         WRITE(6, '(a,/,a)')
     1 ' Normalized Soloveev boundary: ',
!     2 ' Psi = .5*(1 + xi*r**2)*z**2 + (Psi/delta_r**4)*(r**2 - 1)**2',
     2 ' Psi = eps^-2[(r^2 + bdia^2)(z/kappa)^2 + (r^2 - 1)^2/4] '
         WRITE(6, '(a)')
!     1    ' Psi (>0), Xi, Delta-r (<1): (blank or tab delimited)'
     1     ' Psi (0<Psi<1), eps (<0.5), kappa, bdia: (blank delimited)'
!         READ (*,*) PSI, XI, DELTA_R
         READ (*,*) PSI, EPS, KAPPA, BDIA
         GOTO 1004
      ELSE
         STOP 'Must enter a number between 0 and 2!'
      END IF


 1001 CONTINUE

      ntheta = nu
      nphi = 1
      nphi2 = 1+nphi/2
      mpol = mu
      mpol1 = mpol-1

      IF (itype .ne. 6) THEN
         SELECT CASE (itype)
         CASE (1)
            rbdy = rellip
            zbdy = zellip
         CASE (2)
            rbdy = rdshape
            zbdy = zdshape
         CASE (3)
            rbdy = rbean
            zbdy = zbean
         CASE (4)
            rbdy = rbelt
            zbdy = zbelt
         CASE (5)
            rbdy = rsqr
            zbdy = zsqr
         END SELECT
         rin(:ntheta) = rbdy(1)
         zin(:ntheta) = zbdy(1)
         DO m = 2, 4
            DO j = 1, ntheta
               arg = (m - 1)*twopi*(j - 1)/REAL(ntheta,rprec)
               rin(j) = rin(j) + rbdy(m)*COS(arg)
               zin(j) = zin(j) + zbdy(m)*SIN(arg)
            END DO
         END DO
  
      ELSE

         rin(:ntheta) = zero
         zin(:ntheta) = zero

         DO m = 0, 2
            DO n = 0, 3
               DO j = 1, ntheta
                  arg = twopi*(j - 1)/REAL(ntheta,rprec)
                  rin(j) = rin(j) + rheliac(n,m)*COS(m*arg - n*zeta)
                  zin(j) = zin(j) + zheliac(n,m)*SIN(m*arg - n*zeta)
               END DO
            END DO
         END DO

      ENDIF
 3852 FORMAT(4e18.11)

      GOTO 4000

 1002 CONTINUE

      iread = iread+1
      OPEN(unit=20,file=TRIM(datafile),status='old', iostat=itype)
      IF (itype .ne. 0) STOP 'error reading datafile'
      READ (20, *, iostat=itype) ntheta, nphi, nfp
      IF (itype .ne. 0) THEN
         PRINT *,'Error reading first line of ',TRIM(datafile)
         STOP
      END IF
      IF (ntheta .gt. nu) STOP 'ntheta > nu'
      IF (nphi .gt. nv) STOP 'nphi > nv'
      mpol = mu
      mpol1 = mpol-1
      nphi2 = 1+nphi/2
      i = 0
      DO k = 1, nphi
         DO j = 1, ntheta
            i = i + 1
	      IF (iread .eq. 1) THEN
               READ (20, *, iostat=itype) rin(i), dumphi, zin(i)
	      ELSE
               READ (20, *, iostat=itype) rin(i), zin(i)
            END IF
            IF (itype .ne. 0) THEN
               IF (iread .eq. 1) THEN
	            CLOSE (20)
	            GOTO 1002
	         ELSE
                  PRINT *,' Error reading datafile: iostat = ',itype
                  STOP
               END IF
	      END IF
          END DO
      END DO

      CLOSE(20)

      IF (nphi .eq. 1 
     1   .and. (ALL(zin .ge. 0.) .or. ALL(zin .le. 0.))) THEN
         ncount = ntheta
         DO i = ntheta,1,-1
            IF (zin(i) .eq. 0.) CONTINUE
            ncount = ncount+1
            rin(ncount) = rin(i)
            zin(ncount) =-zin(i)
         END DO
         ntheta = ncount
         IF (ntheta .gt. nu) STOP 'ntheta > nu'
      END IF

      GOTO 4000

 1003 CONTINUE

      ntheta = nu
      rbdy_3d = 0
      zbdy_3d = 0
!
!     try to read data from a wout file
!
      CALL read_wout_curve (datafile, i, nphi20, rbdy_3d, zbdy_3d)

      IF (i .ne. 0) THEN

         OPEN(unit=20,file=datafile,iostat=itype)
         IF (itype .ne. 0) STOP 'error reading datafile'
         READ (20, *) mpol, nphi2, nfp
         IF (mpol .gt. mu) STOP 'mpol-input > mu'

         DO WHILE (itype .eq. 0)
             READ(20, *, iostat=itype, END=200) n, m, dumr, dumz
             IF (itype .ne. 0) CYCLE
             IF (nfp .ne. 0) n = n/nfp
             rbdy_3d(m,n) = dumr
             zbdy_3d(m,n) = dumz
         END DO

 200     CONTINUE
         CLOSE(20)

      END IF

      mpol1 = mpol-1
      nphi  = 2*(nphi2+1)
      IF (nphi2 .le. 0) nphi = 1
      IF (nphi .gt. nv)  STOP 'nphi > nv'

      rin = zero
      zin = zero

      DO m = 0, mpol1
         DO n = -nphi2, nphi2
            i = 0
            DO k = 1, nphi
               zeta = twopi*(k-1)/REAL(nphi,rprec)
               DO j = 1, ntheta
                  i = i + 1
                  arg = twopi*(j - 1)/REAL(ntheta,rprec)
                  rin(i) = rin(i) + rbdy_3d(m,n)*COS(m*arg - n*zeta)
                  zin(i) = zin(i) + zbdy_3d(m,n)*SIN(m*arg - n*zeta)
               END DO
            END DO
         END DO
      END DO

      IF (HB_Parameter == 0) THEN        !!extra resolution, NOT VMEC compatible
         mpol = mpol+2             
      ELSE
         mpol = mpol+1
      END IF

      mpol1 = mpol-1

      GOTO 4000

 1004 CONTINUE

      mpol = mu
      mpol1 = mpol-1
      ntheta = nu
      nphi = 1
      nphi2 = 1+nphi/2
!      IF (DELTA_R .LE. ZERO) STOP 'Delta_r must be > 0'
!      IF (DELTA_R .GT. one) STOP 'Delta_r must be < 1'
      IF (PSI .LE. ZERO) STOP 'Psi must be > 0'
      IF (PSI .GT. ONE)  STOP 'Psi must be <= 1'
      IF (EPS .GT. 0.5)  STOP 'Eps must be <= 0.5'

      A1 = EPS*SQRT(ABS(PSI))
      DO j = 1, ntheta
         arg = twopi*(j - 1)/REAL(ntheta,rprec)
         rin(j) = SQRT(1 + 2*A1*COS(arg))
         zin(j) = A1*KAPPA*SIN(arg)/SQRT(rin(j)**2 + bdia**2)
      END DO

!COMPUTE q(PSI)/q0
      ns = 11
      dels = one/(ns-1)
      DO i = 1, ns
         A1 = (i-1)*dels
         A2 = EPS*SQRT(A1)
         QPSI = 0
         VP   = 0
         DO j = 1, ntheta
            arg = twopi*(j-1)/ntheta
            r1 = SQRT(1 + 2*A2*COS(arg))
            QPSI = QPSI + one/r1**3
            VP   = VP   + one/r1
         END DO
         
         QPSI = QPSI/ntheta
         VP   = (VP/ntheta) * (KAPPA*EPS**2)/2
         WRITE (6, '(a,1p,3e12.3)')' PSI,Q/Q0, VP: ', A1, QPSI, VP
         WRITE (33,'(a,1p,3e12.3)')' PSI,Q/Q0, VP: ', A1, QPSI, VP
      END DO
          
       

!      rmin = SQRT(one - delta_r**2)
!      rmax = SQRT(one + delta_r**2)
!      zin = zero
!      DO j = 1, ntheta/2
!         rin(j) = rmin + REAL(j - 1)*(rmax - rmin)/REAL(ntheta/2)
!         rin(j+ntheta/2) =
!     1      rmax - REAL(j - 1)*(rmax - rmin)/REAL(ntheta/2)
!      END DO
!      DO j = 1, ntheta
!         denom = one + XI*rin(j)**2
!         IF (denom .le. zero) CYCLE            !!avoid separatrix...
!         zin(j) = 2*PSI*(one - (rin(j)**2 - one)**2/delta_r**4)/denom
!         IF (zin(j) .lt. zero) zin(j) = zero
!         zin(j) = SQRT(zin(j))
!         IF (j .gt. ntheta/2) zin(j) = -zin(j)
!      END DO

!
!     Attempt to redistribute points more evenly, based on arc-length...
!
!      DO j = 1, ntheta/2
!         length(j) = (rin(j+1) - rin(j))**2 + (zin(j+1) - zin(j))**2
!         length(j) = SQRT(length(j)) + 1.E-10_dp
!      END DO

!      denom = SUM(one/length(:ntheta/2))
!      denom = (rmax - rmin)/denom

!      DO j = 1, ntheta/2
!         rin(j+1) = rin(j) + denom/length(j)
!      END DO

!      DO j = 2, ntheta/2
!         rin(ntheta - j + 2) = rin(j)
!      END DO

!      DO j = 1,ntheta
!         denom = one + XI*rin(j)**2
!         IF (denom .le. zero) CYCLE            !!avoid separatrix...
!         zin(j) = 2*PSI*(one - (rin(j)**2 - one)**2/delta_r**4)/denom
!         IF (zin(j) .lt. 0.) zin(j) = 0.
!         zin(j) = SQRT(zin(j))
!         IF (j .gt. ntheta/2) zin(j) = -zin(j)
!      END DO



 4000 CONTINUE
c**********************************************************************
c     PERFORM OPTIMIZING CURVE-FIT
c**********************************************************************
      CALL scrunch (rin, zin, rbc, zbs, rbs, zbc, rmnaxis,
     1   zmnaxis, ftol, niter, nstep, mexp, mntot)

c**********************************************************************
c     PRINT OPTIMIZED TRANSFORM COEFFICIENTS TO "OUTCURVE" FILE
c     AND STORE RIN, ZIN, RBC, ZBS DATA IN PLOTOUT FILE
c**********************************************************************
      CALL printit(rin,zin,rbc,zbs,rbs,zbc,rmnaxis,zmnaxis)
      CLOSE(10)

!DEC$ IF .NOT.DEFINED(WIN32)
c**********************************************************************
c     PLOT FIT TO DATA POINTS
c**********************************************************************
      WRITE (6, '(/,a)', advance='no')
     1      ' Do you wish to plot this data (Y/N)? '
      READ (*,*) ch_yn

      IF (ch_yn == 'Y' .or. ch_yn == 'y') CALL system ("xdes_plot")

!DEC$ ENDIF

      END PROGRAM descur
