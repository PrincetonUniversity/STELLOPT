      SUBROUTINE dkes_input_prepare (arg, numargs, dkes_input_file, 
     1                               dealloc)
!
!   This code prepares an input file for DKES based on data
!   in the boozmn.* file which is periodically written out
!   by the VMEC optimizer.  After compiling (see instructions
!   given below) and running a file called input.dkes is written.
!   DKES can then be run by typing "xdkes input.dkes".
!   this input file is based on a fixed set of parameters
!   specified below that is intended to provide
!   a rapidly evaluated (high collisionality) optimization target.
!   these can be modified to consider other parameter ranges and
!   a more complete mode spectrum as the need arises.  parameters
!   which may be of interest to vary are:
!
!   The variables in the call to dkes_input_prepare are:
!
!     booz_file_name (=arg(1))
!             = name of boozmn file containing |b| data at various surfaces
!
!     nsurf (=arg(2))
!             = flux surface index (in boozer file) where
!                   DKES is to evaluate transport coefficients
!
!     cmul (arg(3))
!             = DKES collisionality parameter: nu/v  (in meter**-1)
!
!     efield (arg(4))
!             = DKES electric field parameter: E_s/v
!
!     lscreen (arg(5))
!             = logical variable which controls screen output
!            (= .true. screen output on, = .false. screen output off)
!
!     Optional
!     filename modifer (arg(6)) - see description below
!     coupling order   (arg(7))
!     legendre_modes   (arg(8))
!
!
!     max_bmns = Number of Bmnc's which are retained for the
!                  DKES spectrum (this also influences the
!                  Fourier spectrum used for the distribution
!                  function in DKES)
!
!     filename_modifer = append a modifer to the output filename
!
!     legendre_modes = number of Legendre polynomials used in
!                  DKES to represent the pitch angle variation
!                  of the distribution function
!
!     coupling_order = parameter for mode coupling order.
!                  this is the number of iterations used to
!                  generate the optimal m,n Fourier spectrum for
!                  the DKES distribution function.  at high
!                  collisionalities (cmul > 0.1), 2 iterations seems
!                  to be adequate (here adequate means that the upper
!                  and lower bounds coming out of DKES are CLOSE to
!                  to each other for the transport coefficient of
!                  interest).  As one goes to lower collisionalities,
!                  progressively more iterations must be used to get
!                  convergence of the upper/lower bounds.  This will
!                  imply increasingly longer run times for DKES.
!
!
!  Authors:  W.I. van Rij  --->  wrote original version of SPECFIL code (on which this is based)
!            H. Maassberg  --->  Adapted SPECFIL to work within an interactive
!                                code to prepare input for DKES
!            D. Spong      --->  converted to f90, made to read new
!                                    form of Bmn file, different spectrum
!                                    generation rules tried, new output files,
!                                    interactive input turned off.  Adapted
!                                    for running with VMEC optimizer. Added
!                                    comments to the spectrum optimization.
!
!
!-----------------------------------------------------------------------
!
!     Description of the DKES_INPUT_PREP code:      (W.I. van Rij)
!
!     Purpose:
!     Estimate the Fourier spectrum for the DKES2 distribution FUNCTION
!
!     Let S denote the fourier spectrum of 1/B [except for (0,0) term]
!     Let S(zeta) denote the dominant zeta-dependent spectrum of 1/B
!     Define S(theta) = S - S(zeta)
!
!     Let nszt (<=nhi) be the number of Fourier modes in S(zeta)
!     Let nsth (<=mhi) be the number of Fourier modes in S(theta)
!              (nhi, mhi are defined in PARAMETER statement)
!
!     For i-th mode of S(zeta) load m value in mzeta(i) and n value in
!     nzeta(i) [i=1,2,....,nszt]
!     For i-th mode of S(theta) load m value in mtheta(i) and n value in
!     ntheta(i) [i=1,2,....,nsth]
!
!     n values are expressed in units of nzperiod
!
!     n values must satisfy   -nhi <= n <= nhi-1
!     m values must satisfy      0 <= m <= mhi-1
!     These bounds may be violated depending on the value of the
!     coupling order NORD (error message is printed on terminal).
!     in this case, you are prompted for a reduced value of NORD.
!
!-----------------------------------------------------------------------
!     Subroutines required:
!     FILIT, SORTI
!-----------------------------------------------------------------------
!
      USE stel_kinds
      USE date_and_computer
      USE read_boozer_mod
      USE read_wout_mod, rmnc_w=>rmnc, zmns_w=>zmns, lmns_w=>lmns,
     &    xm_w=>xm, xn_w=>xn, phip_w=>phip, mpol_w=>mpol,
     &    ntor_w=>ntor, nfp_w=>nfp, ns_w=>ns, mnmax_w=>mnmax
      USE safe_open_mod
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER       :: numargs, dealloc
      CHARACTER*(*) :: arg(numargs), dkes_input_file
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!     REAL(rprec), PARAMETER :: radial_position = 0.5_dp
!     REAL(rprec), PARAMETER :: cmul = 0.05_dp
!     REAL(rprec), PARAMETER :: efield = 0._dp
      INTEGER, PARAMETER :: max_bmns = 25
      INTEGER            :: legendre_modes = 100
      INTEGER            :: coupling_order = 4
      INTEGER, PARAMETER :: mhi = 100
      INTEGER, PARAMETER :: ibbi = 1
      INTEGER, PARAMETER :: nhi = 120
      INTEGER, PARAMETER :: nfcd = 100
      INTEGER, PARAMETER :: mtd = 15
      INTEGER, PARAMETER :: nzd = 15
      INTEGER, PARAMETER :: mtd1 = mtd + 1
      INTEGER, PARAMETER :: nzd1 = 2*nzd + 1
      REAL(rprec), PARAMETER :: aval1 = 1.e-4_dp
      REAL(rprec), PARAMETER :: aval2 = 1.e-4_dp
      REAL(rprec), PARAMETER :: one = 1, two = 2,
     1  half = 0.5_dp, zero = 0
      CHARACTER*(50), PARAMETER ::
     1   banner = 'THIS IS THE DKES input prep CODE Version 1.0'
!
!    aval1 = lower limit for the |B| Fourier coefficients
!    aval2 = lower limit for dominant zeta dependent spectrum
!         (note: either one or the other of these is active, but not
!          both.- DAS)
!    The parameters nhi and mhi are described in comments to follow.
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: bmn_local,
     1   abs_blocal, blocal_ordered
      INTEGER, DIMENSION(:), ALLOCATABLE :: m_ordered, n_ordered
      INTEGER, DIMENSION(:), ALLOCATABLE :: nsurf
      INTEGER :: iunit, js, iloc, n_norm, i, j, k, ierr, istat, maxl(1)
      REAL(rprec) btheta, bzeta, iota_local, psip, chip, twopi
      REAL(rprec) phi_edge, cmul, efield, bmn_test
      INTEGER, DIMENSION(mhi,nhi,3) :: nplus, nmins
      INTEGER, DIMENSION(2*nhi,mhi) :: mns_out = 0
      INTEGER, DIMENSION(2*nhi) :: n_out = 0
      INTEGER, DIMENSION(nfcd) :: mtheta, ntheta, mzeta, nzeta, ivec
      INTEGER, DIMENSION(nhi) :: nvalsb
      INTEGER, DIMENSION(nfcd) :: intor, impol
      INTEGER :: m, n, mp, nt, mp0, nt0, mpu, ntu, ntorb, mpolb,
     1  nszt, nsth, kk, jj, ii, ntor, mpol, ij, ik, im, in, iout,
     2  ier, nord, nmin, nmax, mmin, mmax
      REAL(rprec), DIMENSION(nfcd) :: afc
      REAL(rprec), DIMENSION(mtd1,nzd1) :: alf
      REAL(rprec) a, alfm, alfu, eps, au, aa, alfz
      REAL(rprec) hs, vnorm, vol_inner, vol_outer,
     1  rminor_i, rminor_o
      INTEGER :: index_dat, index_end
      REAL(rprec) time_begin, time_end
      CHARACTER*(120) :: extension
      CHARACTER*(120) :: extension_mod
      CHARACTER*(10) :: date0, time0, zone0
      CHARACTER*(40) :: dateloc
!     CHARACTER*30 dkes_input_file
      INTEGER :: imon, index_surf
      LOGICAL :: lscreen = .true.
!-----------------------------------------------
      twopi = 8*ATAN(one)
      extension_mod = ''

      CALL second0(time_begin)

      IF (numargs < 4) STOP 'Require at least 4 command line arguments!'

      READ (arg(2),*) js
      READ (arg(3),*) cmul
      READ (arg(4),*) efield
      IF (numargs>4 .and. (arg(5)(1:1).eq.'f' .or. arg(5)(1:1).eq.'F'))
     1   lscreen = .false.
      IF (numargs > 5) READ (arg(6), *) extension_mod
      IF (numargs > 6) READ (arg(7), *) coupling_order
      IF (numargs > 7) READ (arg(8), *) legendre_modes

      index_dat = INDEX(arg(1),'.')
      index_end = LEN_TRIM(arg(1))
      IF (index_dat == index_end) index_dat = INDEX(arg(1),'_')   !SAL BOOZER netCDF file
      extension  = arg(1)(index_dat+1:index_end)
      k=0
      IF (dealloc .ne. 0) CALL read_boozer_file (extension, k)
      IF (k .ne. 0) STOP 'Error reading boozmn file in DKES_INPUT_PREP'
      ALLOCATE(bmn_local(mnboz_b), abs_blocal(mnboz_b),
     1    blocal_ordered(mnboz_b), stat=ierr)
      ALLOCATE(m_ordered(mnboz_b), n_ordered(mnboz_b), stat=ierr)
      ALLOCATE(nsurf(ns_b), stat=ierr)

      IF (lscreen) THEN
         PRINT 48
         CALL DATE_AND_TIME(date0,time0,zone0)
         READ (date0(5:6),'(i2)') imon
         WRITE (dateloc,100) months(imon),date0(7:8),date0(1:4),
     1         time0(1:2),time0(3:4),time0(5:6)
         WRITE (*,'(1x,a,/,1x,a)') banner, dateloc
         PRINT *
         WRITE(*,'(" js = ",i4,", No. Bmn-s = ", i4)') js, 
     1         MIN(max_bmns,mnboz_b)
         WRITE(*,'(" cmul = ",f10.6,", efield = ",f10.6)')cmul,efield
         WRITE(*,'(" coupling order = ",i2,
     1    ", Legendre modes = ",i4)') coupling_order,legendre_modes
         PRINT *
      END IF
 100  FORMAT('DATE = ',a3,' ',a2,',',a4,' ',' TIME = ',2(a2,':'),a2)
 120  FORMAT(/,' TIME IN DKES input prep CODE:',1pe12.2,' SEC')
  48  FORMAT('====================================================')

!
!     Read data from the wout file and allocate storage:
!
      k=0
      IF (dealloc .ne. 0) CALL read_wout_file (extension, k)
      IF (k .ne. 0)STOP 'Error reading wout file in DKES_INPUT_PREPARE'
!
!     Find out which surfaces are not included in the boozmn file
!         (i.e., by checking if Bmn's are zero)
!
      index_surf = 0
      DO j = 1, ns_b
         bmn_test = zero
         DO i = 1,mnboz_b
            abs_blocal(i) = ABS(bmnc_b(i,j))
         END DO
         bmn_test = SUM(abs_blocal)
         IF(bmn_test .gt. 1.e-8_dp) THEN
            index_surf = index_surf + 1
            nsurf(index_surf) = j
         ENDIF
      END DO


      DO i = 1,mnboz_b
         bmn_local(i) = bmnc_b(i,js)
         abs_blocal(i) = ABS(bmn_local(i))
      END DO
      bmn_test = SUM(abs_blocal)
      IF(bmn_test .lt. 1.e-8_dp .and. lscreen) THEN
         WRITE(*,'(" This surface is not in the boozmn file")')
         WRITE(*,'("  - need to pick a different surface")')
         WRITE(*,'(" The following surfaces are available:")')
         WRITE(*,'(6(3x,i2))') (nsurf(j), j=1,index_surf-1)
 
         CALL second0(time_end)
         IF (lscreen) THEN
            WRITE (*,120) time_end - time_begin
            WRITE (*,48)
         END IF

         STOP 23
      ENDIF

      iota_local = iota_b(js)
      btheta = half*(buco(js) + buco(js-1))
      bzeta = half*(bvco(js) + bvco(js-1))
      phi_edge = ABS(phi_b(ns_b))

      hs = one/(ns_w - 1)
      vnorm = TWOPI*TWOPI*hs
      vol_inner = SUM(vp(2:js-1))
      vol_outer = vol_inner + vp(js)
      vol_inner = vol_inner*vnorm
      vol_outer = vol_outer*vnorm
      rminor_i = SQRT(two*vol_inner/(TWOPI*TWOPI*Rmajor))
      rminor_o = SQRT(two*vol_outer/(TWOPI*TWOPI*Rmajor))
      psip = -phip_w(js)*hs/(rminor_o - rminor_i)
      chip = psip*iota_local
!       WRITE(*,'(" ns_w = ",i3)') ns_w

!
!     Sort the Bmn's (along with associated m's and n's in order
!     of increasing ABS(Bmn):
!
      blocal_ordered = 0.0
      m_ordered = 0.0
      n_ordered = 0.0
      DO i=1,mnboz_b
         maxl = MAXLOC(abs_blocal)
         iloc = maxl(1)

         blocal_ordered(i) = bmn_local(iloc)
         m_ordered(i) = ixm_b(iloc)
         n_ordered(i) = ixn_b(iloc)
         abs_blocal(iloc) = zero
      END DO
!
!     Find MIN/MAX m and n out of the first max_bmns of Bmn
!     modes. Compute ntorb and mpolb for DKES.
!
      nmin = 100
      mmin = 100
      nmax = -100
      mmax = -100
      DO i = 1,MIN(max_bmns,mnboz_b)
         n_norm = n_ordered(i)/nfp_b
         IF(m_ordered(i) .lt. mmin) mmin = m_ordered(i)
         IF(m_ordered(i) .gt. mmax) mmax = m_ordered(i)
         IF(n_norm .lt. nmin) nmin = n_norm
         IF(n_norm .gt. nmax) nmax = n_norm
      END DO
      ntorb = nmax - nmin + 1
      mpolb = mmax - mmin + 1
!
!     Write out the initial part of the DKES input file.
!
      dkes_input_file = 'input_dkes.' // TRIM(extension) 
     1                                // TRIM(extension_mod)
      iunit = 15
      CALL safe_open(iunit, istat, dkes_input_file, 'replace',
     1    'formatted')
      WRITE (iunit,'(1x,"&dkes_indata")')
      WRITE (iunit,'(1x,"nzperiod= ",i2,",")') nfp_b
      WRITE (iunit,'(1x,"lalpha= ",i3,", nrun = 1,")') legendre_modes
      WRITE (iunit,'(1x,"cmul = ",e12.4,",")') cmul
!      WRITE (iunit,'(1x,"efield = ",f7.4,",")') efield
      WRITE (iunit,'(1x,"efield = ",e12.4,",")') efield
      WRITE (iunit,'(1x,"mpolb = ",i2,",",2x,"ntorb = ",
     1     i2,",",2x,"ibbi = 1,")') mpolb, ntorb
      WRITE (iunit,'(1x,"chip = ",f7.4,",","  psip = ",f7.4,",")')
     1     chip, psip
      WRITE (iunit,'(1x,"btheta = ",es20.10,","," bzeta = ",
     1                    es20.10,",")')
     1     btheta, bzeta
!
!     Write out ONLY the top max_bmns largest Bmn's:
!
      DO i=1,MIN(max_bmns,mnboz_b)
         n_norm = n_ordered(i)/nfp_b
         IF (m_ordered(i) .le. 9) THEN
            j = 1
         ELSE IF (m_ordered(i).gt.9 .and. m_ordered(i).lt.100) THEN
            j = 2
         ELSE
            j = 3
         ENDIF
         IF (n_norm .ge. 0  .and. n_norm .le. 9) THEN
            k = 1
         ELSE IF ((n_norm.lt.0 .and. n_norm.ge.-9) .or. 
     1            (n_norm.ge.10 .and. n_norm.lt.100)) THEN
            k = 2
         ELSE IF((n_norm.le.-10 .and. n_norm.ge.-99) .or.
     1           (n_norm.ge.100)) THEN
            k = 3
         ELSE 
            k = 4
         END IF

         WRITE (extension,'(a,i1,a,i1,a)') 
     1          '(" borbi(",i', k, '",",i', j, ',")= ",1x,e11.5,",")'
         WRITE (iunit,extension) n_norm,m_ordered(i),blocal_ordered(i)

      END DO
!
!     Translate m, n, and Bmn arrays into arrays used in the
!       GIDKES2 code:
!
      nord = coupling_order
      eps = MAX(zero,MIN(aval1,0.01_dp))
      alfm = zero
      alf  = 0.0
      DO i=1,MIN(max_bmns,mnboz_b)
         mp = m_ordered(i)
         nt = n_ordered(i)/nfp_b
         a = blocal_ordered(i)
         impol(i) = mp
         intor(i) = nt
         afc(i) = a
         IF (mp<=mtd .and. ABS(nt)<=nzd) alf(mp+1,nzd+1+nt) = a
         IF (mp.eq.0 .and. nt.eq.0) THEN
            mp0 = mp
            nt0 = nt
         ELSE IF (ABS(a) > alfm) THEN
            alfm = ABS(a)
            alfu = a
            mpu = mp
            ntu = nt
         ENDIF
      END DO
      alfz = aval2
!     alfz = MIN(half*alfm,MAX(0.2_dp*alfm,10*eps))
!
!     The following loop assigns modes to S(zeta) (non-zero n modes through
!     the mzeta and nzeta arrays) and S(theta) (non-zero m modes through
!     the mtheta and ntheta arrays).  there is some flexibility as how this
!     is done through the alfz parameter defined above.  if alfz is
!     relatively large (i.e., if above alfz is used) not too many modes
!     make it into the mzeta and nzeta arrays, more go into the mtheta,
!     ntheta array.  if alfz is small (comment out the above alfz line)
!     then more modes go into the mzeta and nzeta arrays and not so many
!     into the mtheta and ntheta arrays.
!
      nszt = 0
      nsth = 0
      ntorb = 0
      nvalsb = 0
      DO m = 0, mtd
         au = zero
         l12: DO n = -nzd, nzd
            a = alf(m+1,nzd+1+n)
            aa = ABS(a)
            IF (aa > eps) THEN
               IF (n.ne.0 .and. aa>alfz) THEN
                  nszt = nszt + 1
                  mzeta(nszt) = m
                  nzeta(nszt) = n
               ELSE IF (m.ne.0 .or. n.ne.0) THEN
                  nsth = nsth + 1
                  mtheta(nsth) = m
                  ntheta(nsth) = n
               ENDIF
!
!       Keep track of the number of modes (mpolb, ntorb); also set
!       up the array of n's (nvalsb) as used in the old dkes input
!       file FORMAT (DAS).
!
               mpolb = m + 1
               IF (ntorb .eq. 0) THEN
                  ntorb = 1
                  nvalsb(1) = n
               ELSE
                  DO i = 1, ntorb
                     IF (nvalsb(i) .eq. n) CYCLE  l12
                  END DO
                  ntorb = ntorb + 1
                  nvalsb(ntorb) = n
               ENDIF
            ENDIF
         END DO l12
      END DO
      
      CALL sorti (nvalsb, ntorb)
!
!     estimate the fourier spectrum for the distribution function
!     -----   original SPECFIL version   -----
!
!     the nplus and nmins arrays are used to indicate where m,n
!     modes are used for the distribution function. nplus holds the
!     n .ge. 0 modes while nmins holds the n .lt. 0 modes. these
!     arrays contain the current 3 interates of s through the third
!     index.
!
      nord = MAX(2,nord)

      nplus(:mhi,:nhi,:) = 0
      nmins(:mhi,:nhi,:) = 0

      nplus(1,1,1) = 1
!
!     Generate S(0) = S(theta) + S(zeta)
!
      ier = 0
      DO i = 1, nsth
         m = mtheta(i)
         n = ntheta(i)
         CALL filit (m, n, 1, nplus, nmins, mhi, nhi, ier)
      END DO
      IF (ier .ne. 0) STOP 'FILIT ERROR #1 IN DKES_INPUT_PREPARE'

      DO i = 1, nszt
         m = mzeta(i)
         n = nzeta(i)
         CALL filit (m, n, 1, nplus, nmins, mhi, nhi, ier)
      END DO
      IF (ier .ne. 0) STOP 'FILIT ERROR #2 IN DKES_INPUT_PREPARE'

      nplus(:mhi,:nhi,2) = nplus(:mhi,:nhi,1)
      nmins(:mhi,:nhi,2) = nmins(:mhi,:nhi,1)
!
!     Generate S(1) = S(0)*S(zeta)
!

      DO i = 1, nszt
         DO in = 1, nhi
            DO im = 1, mhi
               IF (nplus(im,in,1) .eq. 1) THEN
                  m = mzeta(i) + im - 1
                  n = nzeta(i) + in - 1
                  CALL filit (m, n, 2, nplus, nmins, mhi, nhi, ier)
                  IF (ier .ne. 0) 
     1               STOP 'FILIT ERROR #3 IN DKES_INPUT_PREPARE'
                  m = mzeta(i) - im + 1
                  n = nzeta(i) - in + 1
                  CALL filit (m, n, 2, nplus, nmins, mhi, nhi, ier)
                  IF (ier .ne. 0) 
     1               STOP 'FILIT ERROR #4 IN DKES_INPUT_PREPARE'
               ENDIF
               IF (nmins(im,in,1) .eq. 1) THEN
                  m = mzeta(i) + im - 1
                  n = nzeta(i) - in
                  CALL filit (m, n, 2, nplus, nmins, mhi, nhi, ier)
                  IF (ier .ne. 0)
     1               STOP 'FILIT ERROR #5 IN DKES_INPUT_PREPARE'
                  m = mzeta(i) - im + 1
                  n = nzeta(i) + in
                  CALL filit (m, n, 2, nplus, nmins, mhi, nhi, ier)
                  IF (ier .ne. 0)
     1               STOP 'FILIT ERROR #6 IN DKES_INPUT_PREPARE'

               ENDIF
            END DO
         END DO
      END DO

      kk = 2
      IF (nord .ne. 2) THEN
!
!     Generate S(n) = S(zeta)*S(n-1) + S(theta)*S(n-2):
!
!     Note: For the S(n) = S(zeta)*S(n-1) + S(theta)*S(n-2) recurrence
!      relation (this orders eps_tor << eps_hel) set ii = 0 and jj = 1 initially.
!      For the S(n) = S(zeta)*S(n-1) + S(theta)*S(n-1) recurrence
!      relation set ii = 1 and jj = 0 initially.
!
         ii = 1
         jj = 0
         DO k = 3, nord
            ii = ii + 1
            IF (ii .eq. 4) ii = 1
            jj = jj + 1
            IF (jj .eq. 4) jj = 1
            kk = kk + 1
            IF (kk .eq. 4) kk = 1
            nplus(:mhi,:nhi,kk) = nplus(:mhi,:nhi,jj)
            nmins(:mhi,:nhi,kk) = nmins(:mhi,:nhi,jj)
!
!        Generate the S(zeta)S(n-1) convolution:
!
            DO i = 1, nszt
               DO in = 1, nhi
                  DO im = 1, mhi
                     IF (nplus(im,in,jj) .eq. 1) THEN
                        m = mzeta(i) + im - 1
                        n = nzeta(i) + in - 1
                        CALL filit(m,n,kk,nplus,nmins,mhi,nhi,ier)
                        IF (ier .ne. 0) GO TO 72
                        m = mzeta(i) - im + 1
                        n = nzeta(i) - in + 1
                        CALL filit(m,n,kk,nplus,nmins,mhi,nhi,ier)
                        IF (ier .ne. 0) GO TO 72
                     ENDIF
                     IF (nmins(im,in,jj) .eq. 1) THEN
                        m = mzeta(i) + im - 1
                        n = nzeta(i) - in
                        CALL filit(m,n,kk,nplus,nmins,mhi,nhi,ier)
                        IF (ier .ne. 0) GO TO 72
                        m = mzeta(i) - im + 1
                        n = nzeta(i) + in
                        CALL filit(m,n,kk,nplus,nmins,mhi,nhi,ier)
                        IF (ier .ne. 0) GO TO 72
                     ENDIF
                  END DO
               END DO
            END DO
!
!        Generate the S(theta)S(n-2) convolution:
!

            DO i = 1, nsth
               DO in = 1, nhi
                  DO im = 1, mhi
                     IF (nplus(im,in,ii) .eq. 1) THEN
                        m = mtheta(i) + im - 1
                        n = ntheta(i) + in - 1
                        CALL filit(m,n,kk,nplus,nmins,mhi,nhi,ier)
                        IF (ier .ne. 0) GO TO 72
                        m = mtheta(i) - im + 1
                        n = ntheta(i) - in + 1
                        CALL filit(m,n,kk,nplus,nmins,mhi,nhi,ier)
                        IF (ier .ne. 0) GO TO 72
                     ENDIF
                     IF (nmins(im,in,ii) .eq. 1) THEN
                        m = mtheta(i) + im - 1
                        n = ntheta(i) + in
                        CALL filit(m,n,kk,nplus,nmins,mhi,nhi,ier)
                        IF (ier .ne. 0) GO TO 72
                        m = mtheta(i) - im + 1
                        n = ntheta(i) - in
                        CALL filit(m,n,kk,nplus,nmins,mhi,nhi,ier)
                        IF (ier .ne. 0) GO TO 72
                     ENDIF
                  END DO
               END DO
            END DO
         END DO
      ENDIF

      ntor = 0
      mpol = 0
      ii = 0
      DO n = 1, nhi
         jj = 0
         DO m = 1, mhi
            IF (nplus(m,n,kk) .eq. 1) THEN
               jj = jj + 1
               ivec(jj) = m - 1
            ENDIF
         END DO
         IF (jj .ne. 0) THEN
            ii = ii + jj
            ntor = ntor + 1
            mpol = MAX(mpol,jj)
!         Capture the n's and m's into arrays n_out and mns_out
!           for later printout:
            n_out(n+nhi) = n-1
            DO iout=1,jj
              mns_out(n+nhi,iout) = ivec(iout)
            END DO

            ij = MIN(13,jj)
            IF (ntor <= 9) THEN
               WRITE (iunit, 101) ntor, n - 1, jj, (ivec(i),i=1,ij)
            ELSE IF (ntor <= 99) THEN
               WRITE (iunit, 201) ntor, n - 1, jj, (ivec(i),i=1,ij)
            ELSE
               WRITE (iunit, 301) ntor, n - 1, jj, (ivec(i),i=1,ij)
            ENDIF
            IF (jj > 13) THEN
               ik = -4
   61          CONTINUE
               ik = ik + 18
               ij = MIN(jj,ik + 17)
               WRITE (iunit, 202) (ivec(i),i=ik,ij)
               IF (ij < jj) GO TO 61
            ENDIF
         ENDIF
      END DO

      DO n = 1, nhi
         jj = 0
         DO m = 1, mhi
            IF (nmins(m,n,kk) .eq. 1) THEN
               jj = jj + 1
               ivec(jj) = m - 1
            ENDIF
         END DO
         IF (jj .ne.0) THEN
            ii = ii + jj
            ntor = ntor + 1
            mpol = MAX(mpol,jj)
!         Capture the n's and m's into arrays n_out and mns_out
!           for later printout:
            n_out(nhi-n+1) = -n
            DO iout=1,jj
              mns_out(nhi-n+1,iout) = ivec(iout)
            END DO

            ij = MIN(13,jj)
            IF (ntor <= 9) THEN
               WRITE (iunit, 101) ntor, (-n), jj, (ivec(i),i=1,ij)
            ELSE IF (ntor <= 99) THEN
               WRITE (iunit, 201) ntor, (-n), jj, (ivec(i),i=1,ij)
            ELSE
               WRITE (iunit, 301) ntor, (-n), jj, (ivec(i),i=1,ij)
            ENDIF
            IF (jj > 13) THEN
               ik = -4
   64          CONTINUE
               ik = ik + 18
               ij = MIN(jj,ik + 17)
               WRITE (iunit, 202) (ivec(i),i=ik,ij)
               IF (ij < jj) GO TO 64
            ENDIF
         ENDIF
      END DO

      WRITE (iunit, 203) mpol, ntor
  101 FORMAT(' mmnn(1,',i1,')=',15(i3,','))
  201 FORMAT(' mmnn(1,',i2,')=',15(i3,','))
  301 FORMAT(' mmnn(1,',i3,')=',15(i3,','))
  202 FORMAT(1x,18(i3,','))
  203 FORMAT(' mpol=',i2,', ntor=',i3,' /')

      GO TO 73
  72  WRITE(6, *) "error in CALL to filit in dkes_input_prepare"
  73  CONTINUE
      CLOSE(unit= iunit)

      CALL second0(time_end)
      IF (lscreen) THEN
         WRITE (6,120) time_end - time_begin
         WRITE (6,48)
      END IF


      DEALLOCATE(bmn_local, abs_blocal, blocal_ordered,
     >   stat=istat)
      DEALLOCATE(m_ordered, n_ordered, stat=istat)
      DEALLOCATE(nsurf, stat=istat)

      IF (dealloc .ne. 0) CALL read_boozer_deallocate
      IF (dealloc .ne. 0) CALL read_wout_deallocate

      END SUBROUTINE dkes_input_prepare
