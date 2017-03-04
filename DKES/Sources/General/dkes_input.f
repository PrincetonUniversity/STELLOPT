      MODULE dkes_input
!     Modified indata namelist to be named dkes_input so as not to
!     conflict with VMEC input namelist.
!
!     CONTAINS INPUT PARAMETERS AND ARRAYS READ IN FROM THE DATAIN NAMELIST
!
!     mpol      : Number of Fourier harmonics used to represent the theta
!                 dependence of the distribution functions; 1 < mpol <= mpold
!     mpolb     : Number of Fourier harmonics used to represent the theta
!                 dependence of the equilibrium magnetic field; 1 < mpolb <= mpolbd
!     ntor      : Number of Fourier harmonics used to represent the zeta
!                 dependence of the distributions; 0 < ntor <= ntord
!     ntorb     : Number of Fourier harmonics used to represent the zeta
!                 dependence of the equilibrium field; 0 < ntorb <= ntorbd
!     mmnn      : Matrix containing the distribution of toroidal mode numbers
!                 (expressed in units of nzperiod) in column 1, and the
!                 number of poloidal modes associated with each toroidal
!                 mode in column 2. the poloidal mode numbers associated
!                 with each toroidal mode are stored in the same row,
!                 starting in column 3. WARNING : this spectrum must
!                 encompass the spectrum of the input magnetic field.

!     nrun      : Number of cmul = nu/v  and  efield = Erad/v values
!                 for which solutions are to be calculated; must be < krun
!     cmul      : Corresponding array of cmul values; should be non-zero.
!     efield    : Corresponding array of efield values.
!     nzperiod  : Number of toroidal field periods.
!     chip      : chip = d(chi)/d(rho), radial derivative of the poloidal flux.
!     psip      : psip = d(psi)/d(rho), radial derivative of the toroidal flux
!     btheta    : covariant poloidal component of B. In Boozer coordinates,
!                 a flux function (I, the TOROIDAL current density)
!     bzeta     : covariant toroidal component of B. In Boozer coordinates,
!                 a flux function (J, the POLOIDAL current density)
!     borbi     : Array of nonzero Fourier coefficients of |B| (or 1/|B|) defined by
!                 borbi_realspace (theta,zeta) = SUM borbi(n,m)*cos[m*theta - nzperiod*n*zeta],
!                 where m = 0 to mpolb, n = -ntorb to ntorb
!                 Note: borbi represents "|B|" or 1/|B| ("binverse"), depending on ibbi)
!                 OBSOLETE (see nvalsb array below):
!                 borbi_realspace(theta,zeta) = borbi(n,m)*cos[(m-1)*theta - nzperiod*nvalsb(n)*zeta].
!     ibbi  = 1 : borbi represents B (default, comes from boozmn file, for example)
!     ibbi  = 2 : borbi represents 1/B
!     nvalsb    : (NOW OBSOLETE) Array of ntorb toroidal mode numbers
!                 expressed in units of nzperiod. For an axisymmetric
!                 device, ntorb = 1 and nvalsb(1) = 0.

!     ipmb  = 0 : DKES computes both sine and cosine solutions (default)
!     ipmb  = 1 : DKES computes only sine            solutions.
!     ipmb  = 2 : DKES computes only cosine          solutions.

!     idisk = 0 : DKES computes all l values.
!     idisk = 1 : DKES computes only l = 0,1,2,3,4,5 (default)

!     lfout = 0 : Final distributions not written in the file DKESF. (default)
!     lfout > 0 : Final distributions     written in the file DKESF,
!                 for l+1 = 1,2,.....,lfout <= lalpha.

!     meshtz    : Determines theta and zeta meshes for Fourier transforms.
!                 See ntheta and nzeta in the MNSET subroutine. Used to
!                 ensure accuracy of spectral matrix elements.

!     lalpha    : Number of orthonormalized Legendre polynomials used to
!                 represent the dependence of the distributions on the
!                 pitch-angle variable cos(alpha); lalpha >= 6  
!                 NOTE OF IMPORTANCE: the sqrt(l+1/2) factor is built in!

      USE stel_kinds
      INTEGER, PARAMETER :: mpold = 200
      INTEGER, PARAMETER :: ntord = 200
      INTEGER, PARAMETER :: mpolbd = 100
      INTEGER, PARAMETER :: ntorbd = 100
      INTEGER, PARAMETER :: bigint = 2**16
      INTEGER, PARAMETER :: krun = 100

      INTEGER :: mpol, ntor, lalpha, ipmb, meshtz, idisk,
     1   lfout, ifscl, nrun, nzperiod, mpolb, ntorb, ibbi

      INTEGER, DIMENSION(mpold + 2,ntord) :: mmnn

      REAL(rprec), DIMENSION(krun) :: cmul, efield
      REAL(rprec) :: chip, psip, btheta, bzeta
      REAL(rprec), DIMENSION(-ntorbd:ntorbd,0:mpolbd) :: borbi
      INTEGER, DIMENSION(ntorbd) :: nvalsb                              !OBSOLETE
      CHARACTER*1, PARAMETER :: dashes(131) = (/ ('-', idisk=1,131) /)
      CHARACTER*80 :: summary_file            !record file addition
      INTEGER :: itab_out                     !record file addition
      LOGICAL :: lscreen

      NAMELIST /dkes_indata/mpol, ntor, lalpha, ipmb, mmnn, meshtz, idisk,
     1   lfout, ifscl, nrun, cmul, efield, nzperiod, chip, psip, mpolb,
     2   ntorb, nvalsb, borbi, btheta, bzeta, ibbi

      CONTAINS

      SUBROUTINE read_dkes_input
      USE safe_open_mod
      USE vimatrix, ONLY: ioout, ioout_opt
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: idata = 7
      INTEGER, PARAMETER :: iout = 20
      INTEGER, PARAMETER :: iout_opt = 14
      REAL (rprec), PARAMETER :: zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: numargs, icount, istat, index_blank, iodata
      INTEGER :: m, n, nmax, mmax
      CHARACTER :: output_file*55, opt_file*55, input_file*50
      CHARACTER*50 :: arg1(7)
      CHARACTER :: tb*1            !record file addition
      LOGICAL :: lexist
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL getcarg, dkes_input_prepare
C-----------------------------------------------
      tb = char(9)            !record file addition
      lscreen = .true.

c  Read input data from "datain" namelist
      CALL getcarg(1, arg1(1), numargs)

      DO icount = 2, MIN(7,numargs)
         CALL getcarg(icount, arg1(icount), istat)
      END DO
      IF (numargs .ge. 4) THEN
         IF (numargs .eq. 4) arg1(5) = 'T'
         CALL dkes_input_prepare (arg1, numargs, input_file, 1)
         IF (arg1(5)(1:1).eq.'f' .or. arg1(5)(1:1).eq.'F')
     1      lscreen = .false.
      ELSE IF (numargs .eq. 1) THEN
         input_file = TRIM(arg1(1))
      ELSE
         PRINT *, ' MUST ENTER INPUT FILENAME ON COMMAND LINE'
         STOP
      END IF

      index_blank = INDEX(input_file,'input.')

      IF (index_blank .eq. 0) THEN
         INQUIRE (file=TRIM(input_file), exist=lexist)
         IF (.not.lexist) THEN
            input_file = 'input.' // TRIM(input_file)
            index_blank = 1
         ELSE
            index_blank = INDEX(input_file,'.') - 5
         END IF
      END IF

      index_blank = MAX(index_blank + 5,1)
      output_file= 'dkesout' // input_file(index_blank:)
      opt_file= 'opt_dkes' // input_file(index_blank:)             !DAS 2/21/2000
      summary_file = 'results' // input_file(index_blank:)         !record file addition

!
!     OPEN INPUT AND OUTPUT FILES FOR READING (AND WRITING OUTPUT)
!
      iodata = idata
      CALL safe_open(iodata, istat, input_file, 'old', 'formatted')
      IF (istat .ne. 0) STOP 'Error reading input file in DKES'
      ioout = iout
      CALL safe_open(ioout, istat, output_file, 'replace', 'formatted')
      IF (istat .ne. 0) STOP 'Error writing output file'
      ioout_opt = iout_opt
      CALL safe_open(ioout_opt, istat, opt_file, 'replace','formatted')
      IF (istat .ne. 0) STOP 'Error writing opt_output file'

!     Read namelist (datain) input
      nvalsb(1) = -bigint-1
      idisk = 1
      lfout = 0
      ibbi = 1
      borbi = 0
      READ (iodata, nml=dkes_indata, iostat=istat)
      IF (istat .ne. 0) STOP'Error reading dkes_indata NAMELIST in DKES'
      CLOSE (iodata)

!
!     Recompute ntorb, mpolb for new style input where
!     borbi is input with actual index value, borbi(n,m)
!
      IF (nvalsb(1) <= -bigint) THEN
         nmax = -bigint;  mmax = -bigint
         DO n = -ntorbd, ntorbd
            DO m = 0, mpolbd
               IF (borbi(n,m) .ne. zero) THEN
                  nmax = max (nmax, abs(n))
                  mmax = max (mmax, abs(m))
               END IF
            END DO
         END DO

!        User MAY input smaller values if he does not want to use entire array         
         ntorb = min (abs(ntorb), nmax)
         mpolb = min (abs(mpolb-1), mmax)
      
      ELSE 
         IF (ntorb > ntorbd) THEN
            WRITE (ioout, 45) ntorb, ntorbd
            STOP ' ntorb > ntorbd in DKES input'
         END IF
         IF (mpolb < 2) THEN
            WRITE (ioout, 20) mpolb
            STOP ' mpolb < 2 in DKES input'
         ENDIF
  20     FORMAT(' mpolb = ',i5,'  is less than 2')
  45     FORMAT(' ntorb = ',i5,'  is greater than ntorbd = ',i5)

      END IF
!
!     perform error checking on input data

      IF (nzperiod .le. 0) nzperiod = 1
      IF (ibbi<1 .or. ibbi >2) THEN
         WRITE (ioout,'(a)') ' ibbi must =1 or =2'
         STOP ' ibbi <1 or ibbi >2 in DKES input'
      END IF

      IF (mpol < 1) THEN
         WRITE (ioout, 10) mpol
  10     FORMAT(' mpol = ',i5,'  is less than 1')
         STOP ' mpol < 1 in DKES input'
      ENDIF

      IF (mpol > mpold) THEN
         WRITE (ioout, 15) mpol, mpold
  15     FORMAT(' mpol = ',i5,'  is greater than mpold = ',i5)
         STOP ' mpol > mpold in DKES input'
      ENDIF


      IF (mpolb > mpolbd) THEN
         WRITE (ioout, 25) mpolb, mpolbd
  25     FORMAT(' mpolb = ',i5,'  is greater than mpolbd = ',i5)
         STOP ' mpolb > mpolbd in DKES input'
      ENDIF

      IF (ntor < 1) THEN
         WRITE (ioout, 30) ntor
  30     FORMAT(' ntor = ',i5,'  is less than 1')
         STOP ' ntor < 1 in DKES input'
      ENDIF

      IF (ntor > ntord) THEN
         WRITE (ioout, 35) ntor, ntord
  35     FORMAT(' ntor = ',i5,'  is greater than ntord = ',i5)
         STOP ' ntor > ntord in DKES input'
      ENDIF

      IF (lalpha < 6) THEN
         WRITE (ioout, 50) lalpha
  50     FORMAT(' lalpha = ',i5,'  is less than 6')
         STOP ' lalpha < 6 in DKES input'
      ENDIF

      IF (ipmb<0 .or. ipmb>2) ipmb = 0

      meshtz = MAX(0,meshtz)

      END SUBROUTINE read_dkes_input

      END MODULE dkes_input
