!-----------------------------------------------------------------------
!     Program:       xbooz_xform_ad
!     Purpose:       Command-line access to the primal, tangent-linear,
!                    adjoint and central-difference Boozer spectra.
!
!     Usage:  xbooz_xform_ad in_booz.ext MODE seed.bin out.bin [eps]
!             xbooz_xform_ad in_booz.ext chk [eps]
!        MODE = fwd  primal from boozer_coords (seed ignored)
!               tl   tangent-linear model (seed: input tangents)
!               ad   adjoint model        (seed: output cotangents)
!               fd   central difference of boozer_coords with step eps
!                    along the input tangent (default eps = 1e-6)
!               chk  self-test on every surface of jlist with pseudo-random
!                    seeds: adjoint duality <v,J u> = <J^T v,u> and tl
!                    versus fd; prints the errors and stops with a nonzero
!                    exit code if duality > 1e-10 or tl-fd > 1e-6
!                    (seed and out files are ignored)
!
!     Files are unformatted stream, native-endian.  out.bin starts with
!     the header  mnboz, jsize, lasym (int32), jlist(jsize) (int32),
!     xmb(mnboz), xnb(mnboz) (real64).  It then holds one record per
!     surface of jlist, as does seed.bin (which has no header).
!     Input record (tl/fd seed, ad result), real64:
!        bmnc(mnmax_nyq) bsubumnc(mnmax_nyq) bsubvmnc(mnmax_nyq)
!        lmns(mnmax) iota  [bmns bsubumns bsubvmns (mnmax_nyq) lmnc(mnmax)]
!     Output record (fwd/tl/fd result, ad seed), real64:
!        bmnc_b(mnboz) iota G I  [bmns_b(mnboz)]
!     Bracketed parts are present only for lasym = T.
!-----------------------------------------------------------------------
      PROGRAM booz_ad_driver
      USE stel_kinds, ONLY: rprec
      USE booz_params
      USE booz_persistent, ONLY: xmb, xnb, lmns, lmnc
      USE booz_ad
      USE safe_open_mod
      USE mpi_params, ONLY: master, myid
      IMPLICIT NONE
      INTEGER :: istat, iread, jsurf, jrad, nargs, nin, nout, iseed, iout
      INTEGER :: nq, nm, k
      CHARACTER(LEN=256) :: infile, mode, seedfile, outfile, arg
      CHARACTER(LEN=120) :: extension
      REAL(rprec) :: eps
      REAL(rprec), ALLOCATABLE :: vin(:), vout(:), vp(:), vm(:),        &
     &   vx(:), vy(:)
      REAL(rprec) :: dual, errfd, wdual, wfd
      INTEGER :: iseedrng

      myid = master
      lscreen = .false.
      nargs = COMMAND_ARGUMENT_COUNT()
      CALL GET_COMMAND_ARGUMENT(2, mode)
      IF (TRIM(mode) .eq. 'chk' .and. nargs .ge. 2) nargs = MAX(nargs, 4)
      IF (nargs .lt. 4) THEN
         PRINT *, 'Usage: xbooz_xform_ad in_booz.ext fwd|tl|ad|fd ',     &
     &            'seed.bin out.bin [eps]'
         STOP 1
      END IF
      CALL GET_COMMAND_ARGUMENT(1, infile)
      CALL GET_COMMAND_ARGUMENT(2, mode)
      CALL GET_COMMAND_ARGUMENT(3, seedfile)
      CALL GET_COMMAND_ARGUMENT(4, outfile)
      eps = 1.0E-6_rprec
      IF (TRIM(mode) .eq. 'chk') THEN
         IF (COMMAND_ARGUMENT_COUNT() .ge. 3) READ (seedfile, *) eps
      ELSE IF (nargs .ge. 5) THEN
         CALL GET_COMMAND_ARGUMENT(5, arg)
         READ (arg, *) eps
      END IF

      iread = unit_booz-1
      CALL safe_open(iread, istat, TRIM(infile), 'old', 'formatted')
      IF (istat .ne. 0) STOP 'Error opening input file'
      READ (iread, *, iostat=istat) mboz, nboz
      READ (iread, *, iostat=istat) extension
      IF (istat .ne. 0) STOP 'Error reading input file'
      CALL read_wout_booz(extension, iread, istat)
      CLOSE (unit=iread)
      IF (istat .ne. 0) STOP 'Error in read_wout_booz'
      CALL boozer_setup

      nq = mnmax_nyq; nm = mnmax
      nin = 3*nq + nm + 1
      nout = mnboz + 3
      IF (lasym_b) THEN
         nin = nin + 3*nq + nm
         nout = nout + mnboz
      END IF
      ALLOCATE(vin(nin), vout(nout), vp(nout), vm(nout), vx(nin),       &
     &   vy(nout))

      IF (TRIM(mode) .eq. 'chk') THEN
         iseedrng = 12345
         wdual = 0; wfd = 0
         DO jsurf = 1, jsize
            jrad = jlist(jsurf)
            CALL random_tangent(jrad, vin)
            CALL random_vec(vy)
            CALL tl(jrad, vin, vout)
            CALL ad(jrad, vy, vx)
            dual = ABS(DOT_PRODUCT(vy, vout) - DOT_PRODUCT(vx, vin))     &
     &           / MAX(ABS(DOT_PRODUCT(vy, vout)), TINY(dual))
            CALL perturb(jrad, vin, eps)
            CALL primal(jrad, jsurf, vp)
            CALL perturb(jrad, vin, -2*eps)
            CALL primal(jrad, jsurf, vm)
            CALL perturb(jrad, vin, eps)
            vp = (vp - vm)/(2*eps)
!           spectral blocks (bmnc_b, bmns_b) are measured in one norm,
!           iota, G and I each relative to themselves
            vm = vp - vout
            vy = vout
            vm(mnboz+1:mnboz+3) = 0
            vy(mnboz+1:mnboz+3) = 0
            errfd = SQRT(SUM(vm**2)/MAX(SUM(vy**2), TINY(errfd)))
            DO k = mnboz+1, mnboz+3
               errfd = MAX(errfd, ABS(vp(k) - vout(k))                   &
     &               / MAX(ABS(vout(k)), TINY(errfd)))
            END DO
            WRITE (6, '(a,i5,a,es10.2,a,es10.2)') ' js = ', jrad,       &
     &         '  duality rel err = ', dual, '  tl vs fd rel err = ',   &
     &         errfd
            wdual = MAX(wdual, dual); wfd = MAX(wfd, errfd)
         END DO
         IF (wdual .gt. 1.0E-10_rprec .or. wfd .gt. 1.0E-6_rprec) THEN
            WRITE (6, *) 'booz_ad check FAILED'
            STOP 1
         END IF
         WRITE (6, *) 'booz_ad check passed'
         STOP
      END IF

      OPEN (newunit=iout, file=TRIM(outfile), access='stream',          &
     &      form='unformatted', status='replace', action='write')
      WRITE (iout) mnboz, jsize, MERGE(1, 0, lasym_b), jlist(1:jsize),  &
     &             xmb(1:mnboz), xnb(1:mnboz)
      iseed = -1
      IF (TRIM(mode) .ne. 'fwd')                                        &
     &   OPEN (newunit=iseed, file=TRIM(seedfile), access='stream',     &
     &         form='unformatted', status='old', action='read')

      DO jsurf = 1, jsize
         jrad = jlist(jsurf)
         SELECT CASE (TRIM(mode))
         CASE ('fwd')
            CALL primal(jrad, jsurf, vout)
            WRITE (iout) vout
         CASE ('tl')
            READ (iseed) vin
            CALL tl(jrad, vin, vout)
            WRITE (iout) vout
         CASE ('ad')
            READ (iseed) vout
            CALL ad(jrad, vout, vin)
            WRITE (iout) vin
         CASE ('fd')
            READ (iseed) vin
            CALL perturb(jrad, vin, eps)
            CALL primal(jrad, jsurf, vp)
            CALL perturb(jrad, vin, -2*eps)
            CALL primal(jrad, jsurf, vm)
            CALL perturb(jrad, vin, eps)
            vout = (vp - vm)/(2*eps)
            WRITE (iout) vout
         CASE DEFAULT
            STOP 'Unknown mode'
         END SELECT
      END DO
      CLOSE (iout)
      IF (iseed .ne. -1) CLOSE (iseed)

      CONTAINS

      SUBROUTINE primal(jrad, jsurf, v)
      INTEGER, INTENT(in) :: jrad, jsurf
      REAL(rprec), INTENT(out) :: v(:)
      CALL boozer_coords(jrad, jsurf)
      v(1:mnboz) = bmncb(:,jsurf)
      v(mnboz+1) = hiota(jrad)
      v(mnboz+2) = gpsi(jrad)
      v(mnboz+3) = ipsi(jrad)
      IF (lasym_b) v(mnboz+4:) = bmnsb(:,jsurf)
      END SUBROUTINE primal

!     Shift the jrad columns of the differentiated inputs by eps*v
      SUBROUTINE perturb(jrad, v, eps)
      INTEGER, INTENT(in) :: jrad
      REAL(rprec), INTENT(in) :: v(:), eps
      INTEGER :: k
      k = 0
      bmodmnc(:,jrad)  = bmodmnc(:,jrad)  + eps*v(k+1:k+nq); k = k+nq
      bsubumnc(:,jrad) = bsubumnc(:,jrad) + eps*v(k+1:k+nq); k = k+nq
      bsubvmnc(:,jrad) = bsubvmnc(:,jrad) + eps*v(k+1:k+nq); k = k+nq
      lmns(:,jrad)     = lmns(:,jrad)     + eps*v(k+1:k+nm); k = k+nm
      hiota(jrad)      = hiota(jrad)      + eps*v(k+1);      k = k+1
      IF (.not.lasym_b) RETURN
      bmodmns(:,jrad)  = bmodmns(:,jrad)  + eps*v(k+1:k+nq); k = k+nq
      bsubumns(:,jrad) = bsubumns(:,jrad) + eps*v(k+1:k+nq); k = k+nq
      bsubvmns(:,jrad) = bsubvmns(:,jrad) + eps*v(k+1:k+nq); k = k+nq
      lmnc(:,jrad)     = lmnc(:,jrad)     + eps*v(k+1:k+nm)
      END SUBROUTINE perturb

!     Deterministic pseudo-random numbers in [-1,1] (Park-Miller)
      SUBROUTINE random_vec(v)
      REAL(rprec), INTENT(out) :: v(:)
      INTEGER :: k
      DO k = 1, SIZE(v)
         iseedrng = INT(MOD(16807_8*iseedrng, 2147483647_8))
         v(k) = 2*REAL(iseedrng, rprec)/2147483647.0_rprec - 1
      END DO
      END SUBROUTINE random_vec

!     Tangent proportional to each input harmonic plus 1% of the block max
      SUBROUTINE random_tangent(jrad, v)
      INTEGER, INTENT(in) :: jrad
      REAL(rprec), INTENT(out) :: v(:)
      INTEGER :: k
      CALL random_vec(v)
      k = 0
      CALL scale_block(v, k, bmodmnc(:,jrad))
      CALL scale_block(v, k, bsubumnc(:,jrad))
      CALL scale_block(v, k, bsubvmnc(:,jrad))
      CALL scale_block(v, k, lmns(:,jrad))
      v(k+1) = v(k+1)*ABS(hiota(jrad)); k = k+1
      IF (lasym_b) THEN
         CALL scale_block(v, k, bmodmns(:,jrad))
         CALL scale_block(v, k, bsubumns(:,jrad))
         CALL scale_block(v, k, bsubvmns(:,jrad))
         CALL scale_block(v, k, lmnc(:,jrad))
      END IF
      END SUBROUTINE random_tangent

      SUBROUTINE scale_block(v, k, x)
      REAL(rprec), INTENT(inout) :: v(:)
      INTEGER, INTENT(inout) :: k
      REAL(rprec), INTENT(in) :: x(:)
      v(k+1:k+SIZE(x)) = v(k+1:k+SIZE(x))*(ABS(x) + 0.01_rprec*MAXVAL(ABS(x)))
      k = k + SIZE(x)
      END SUBROUTINE scale_block

      SUBROUTINE tl(jrad, vi, vo)
      INTEGER, INTENT(in) :: jrad
      REAL(rprec), INTENT(in) :: vi(:)
      REAL(rprec), INTENT(out) :: vo(:)
      REAL(rprec) :: z(MAX(nq,nm)), bs(mnboz)
      INTEGER :: o
      z = 0
      o = 3*nq + nm + 1
      IF (lasym_b) THEN
         CALL boozer_coords_tl(jrad, vi(1:nq), vi(nq+1:2*nq),           &
     &      vi(2*nq+1:3*nq), vi(3*nq+1:3*nq+nm), vi(o),                 &
     &      vi(o+1:o+nq), vi(o+nq+1:o+2*nq), vi(o+2*nq+1:o+3*nq),       &
     &      vi(o+3*nq+1:o+3*nq+nm), vo(1:mnboz), bs, vo(mnboz+1),       &
     &      vo(mnboz+2), vo(mnboz+3))
         vo(mnboz+4:) = bs
      ELSE
         CALL boozer_coords_tl(jrad, vi(1:nq), vi(nq+1:2*nq),           &
     &      vi(2*nq+1:3*nq), vi(3*nq+1:3*nq+nm), vi(o),                 &
     &      z(1:nq), z(1:nq), z(1:nq), z(1:nm), vo(1:mnboz), bs,        &
     &      vo(mnboz+1), vo(mnboz+2), vo(mnboz+3))
      END IF
      END SUBROUTINE tl

      SUBROUTINE ad(jrad, vo, vi)
      INTEGER, INTENT(in) :: jrad
      REAL(rprec), INTENT(in) :: vo(:)
      REAL(rprec), INTENT(out) :: vi(:)
      REAL(rprec) :: bs(mnboz), a_bmns(nq), a_bsubus(nq), a_bsubvs(nq), &
     &   a_lmnc(nm), a_bmnc(nq), a_bsubu(nq), a_bsubv(nq), a_lmns(nm),  &
     &   a_iota
      INTEGER :: o
      bs = 0
      IF (lasym_b) bs = vo(mnboz+4:mnboz+3+mnboz)
      CALL boozer_coords_ad(jrad, vo(1:mnboz), bs, vo(mnboz+1),         &
     &   vo(mnboz+2), vo(mnboz+3), a_bmnc, a_bsubu, a_bsubv, a_lmns,    &
     &   a_iota, a_bmns, a_bsubus, a_bsubvs, a_lmnc)
      vi(1:nq) = a_bmnc
      vi(nq+1:2*nq) = a_bsubu
      vi(2*nq+1:3*nq) = a_bsubv
      vi(3*nq+1:3*nq+nm) = a_lmns
      o = 3*nq + nm + 1
      vi(o) = a_iota
      IF (lasym_b) THEN
         vi(o+1:o+nq) = a_bmns
         vi(o+nq+1:o+2*nq) = a_bsubus
         vi(o+2*nq+1:o+3*nq) = a_bsubvs
         vi(o+3*nq+1:o+3*nq+nm) = a_lmnc
      END IF
      END SUBROUTINE ad

      END PROGRAM booz_ad_driver
