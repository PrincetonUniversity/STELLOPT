      MODULE dkes_realspace
      USE stel_kinds
C-----------------------------------------------
C   M o d u l e    V a r i a b l e s
C-----------------------------------------------
!
!     mvalue, nvalue    :    array of m,n values to be used for distribution functions
!                            computed from mmnn array in namelist
!
      INTEGER, DIMENSION(:), ALLOCATABLE :: mvalue, nvalue
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: exbgrad, bgrad, bstrs,
     1   srces0, borbi1, auxs3m
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: auxs1, auxs3p
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: diagl, diagle
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: bmat2, bmat3,
     1   bmat4, bmat5, bmat6, matjac
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: bmat1
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: blk1, blk2, blk3, blk4,
     1   blk5, blk6, blk7
      INTEGER :: mn0, mpnt, mpntp1, mpnt2, mpnt3, mpnt4,
     1   mpnt5, mpntsq, mmax, nmax, ntheta, nzeta
     
      ! ADDED BY SAL to ease optimizer use
      INTEGER :: DKES_rad_dex
      REAL(rprec), ALLOCATABLE, DIMENSION(:) ::
     1   DKES_L11p, DKES_L33p, DKES_L31p,
     2   DKES_L11m, DKES_L33m, DKES_L31m,
     3   DKES_scal11, DKES_scal33, DKES_scal31

      CONTAINS

      SUBROUTINE set_mndim
      USE Vimatrix
      USE Vnamecl2
      USE dkes_input
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mn, mtop, ntop, nn, mbig, mm, istat, m1, n1
C-----------------------------------------------
c  This subroutine sets up fourier mode number arrays and block dimensions
c  and allocates memory for spatial block arrays

!
!     Constants for Legendre polynomial norms
!
      p2 = one/SQRT(20._dp)
      root2 = SQRT(two);          rthalf = one/root2
      rt3o2 = SQRT(3._dp/2._dp);  rt5o2 = SQRT(5._dp/2._dp)

c  stack Fourier modes [(m=0, n<0) components are discarded]
c  ONLY m < mpol modes retained. ONLY ntor toroidal modes retained.

      mn = 0

c  first compute mpnt from mmnn array in NAMELIST

      ntor = MAX(1, ABS(ntor))
      ntor = MIN(ntor, ntord)
      DO nn = 1, ntor
         mbig = 2 + ABS(mmnn(2,nn))
         DO mm = 3, mbig
            m1 = mmnn(mm,nn)
            n1 = mmnn(1,nn)
            IF ((m1 .lt. mpol) .and. (m1.ne.0 .or. n1.ge.0) .and.
     1          ALL(m1 .ne. mmnn(3:mm-1,nn))) mn = mn + 1
         END DO
      END DO

      mpnt = mn

      ALLOCATE (mvalue(mpnt), nvalue(mpnt), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error(1) in DKES2'

      mvalue = -100000;   nvalue = -100000
      mn = 0
      mtop = 0
      ntop = 0
      DO nn = 1, ntor
         mbig = 2 + ABS(mmnn(2,nn))
         DO mm = 3, mbig
            m1 = mmnn(mm,nn)
            n1 = mmnn(1,nn)
            IF ((m1.lt.mpol) .and. (m1.ne.0 .or. n1.ge.0)) THEN
               IF (ANY(m1.eq.mvalue .and. n1.eq.nvalue)) CYCLE           !!do not READ same m,n more than once
               mn = mn + 1
               mvalue(mn) = m1
               nvalue(mn) = n1                                           !!Per/period n
               mtop = MAX(mtop, ABS(m1))
               ntop = MAX(ntop, ABS(n1))
               IF (m1.eq.0 .and. n1.eq.0) mn0 = mn
            ENDIF
         END DO
      END DO

      IF (mn .ne. mpnt)
     1   STOP 'Error counting mmnn modes in DKES set_mndim'

      mpnt2 = 2*mpnt
      mpnt3 = 3*mpnt
      mpnt4 = 4*mpnt
      mpnt5 = 5*mpnt
      mpntp1 = mpnt + 1
      mpntsq = mpnt*mpnt
      mmax = 2*mtop
      nmax = 2*ntop

      ALLOCATE (exbgrad(mpnt), bgrad(mpnt), bstrs(mpnt), auxs3p(mpnt,2),
     1         auxs1(mpnt,3), auxs3m(mpnt), diagl(mpnt,mpnt,2),
     2         diagle(mpnt,mpnt, 2), srces0(16*mpnt),
     3         blk1(mpntsq), blk2(mpntsq), blk3(mpntsq),
     4         blk4(mpntsq), blk5(mpntsq), blk6(mpntsq),
     5         blk7(mpntsq), borbi1(mpnt),
     6         bmat2(mpnt,mpnt,2), bmat4(mpnt,mpnt,2), bmat1(mpnt,mpnt),
     7         bmat3(mpnt,mpnt,2), bmat5(mpnt,mpnt,2),
     8         bmat6(mpnt,mpnt,2), matjac(mpnt,mpnt,2), stat = istat)
      IF (istat .ne. 0) STOP 'allocation error(1) in set_mndim!'


      END SUBROUTINE set_mndim

      SUBROUTINE free_mndim

      DEALLOCATE (exbgrad, bgrad, bstrs, borbi1,
     1   auxs1, auxs3p, auxs3m, diagl, diagle,
     2   srces0, blk1, blk2, blk3, blk4, blk5, blk6, blk7,
     3   bmat2, bmat4, bmat1, matjac, bmat3, bmat5, bmat6)
      DEALLOCATE (mvalue, nvalue)

      END SUBROUTINE free_mndim

      SUBROUTINE ftconv
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE Vnamecl2
      USE dkes_input
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: plus = 1, minus = 2
      REAL(rprec), PARAMETER :: specac = 1.e-12_dp
      REAL(rprec), PARAMETER :: pt1 = 0.1_dp, half = 0.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, DIMENSION(:), ALLOCATABLE :: nval_bfld, mval_bfld
      INTEGER :: nznt, mn, n, m, nn, mm, k, j, l, mnp, mp, np,
     1          ms, ns, md, nd, is, imask(1)
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: trigs, trigc,
     1    rmn, bmn, qmn, q1mn, q2mn, dnorm
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: blank
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: overb1, bsqi1,
     1   trigtc, trigts, trigzc, trigzs
      REAL(rprec) :: dnrm0, dnrm, twopi, dangt, dangz, ang,
     1   bi0, bb1, bi1, bisq1, fac, bi20, facjac,
     2   vpnz, jac, bgradb, blebi, qs, 
     3   jaci, gradbgi, bgbi2gi, bigi, bgbibigi, blbi, qbi, qlbi,
     4   b1, trc, trs, facc, c2, c3, stfac, bi2gi
C-----------------------------------------------
!     This subroutine calculates magnetic field spectral arrays (Fourier
!     transforms that are needed for convolutions).
!
!     overb1        :  2D array of Fourier coefficients of 1/|B|
!                      1/|B|  = SUM [ overb1(m,n) * cos(mu - nv) ]
!
!     bsqi1         :  2D array of Fourier coefficients of 1/|B|*2
!                      1/|B|**2  = SUM [ bsqi1(m,n) * cos(mu - nv) ]
!
!     bgrad         :  m*chip - n*phip (Fourier representation of B dot grad operator)
!
!     exbgrad       : (m*bzeta + n*btheta)/|B00|**2 (multiplier of -Erad X B drift term; multiplied
!                                                   later by -efield to get correct over all sign)
!
!     srces0(1:mpnt):  (m*bzeta + n*btheta)*(1/B**2)(m,n)/3 (Vdrift source for l=0)
!     srces0(mpnt2+1: 3*mpnt)
!                   :  (1/sqrt(20)*sreces0(1:mpnt)  (Vdrift source for l=2 Legendre)
!
!     bgbi2gi       :  (jacobian * B dot grad 1/|B|)**2 / jacobian
!                      (Note mneumonics: bgbi -> jacobian * B dot grad 1/|B|, gi -> 1/jacobian...)
!
!     bgbibigi      :  (jacobian * B dot grad (1/|B|)] / (jacobian |B|)  (bi -> 1/|B|...)
!
!     bi2gi         :  1/(jacobian*|B|**2)
!
!     mval_bfld     :  array of m-vals for |B|, 0-mpolb
!
!     nval_bfld     :  array of n-vals/nfp for |B|, -ntorb:ntorb
!
!     trigtc, trigts:  cos(mu), sin(mu) 2D arrays
!
!     trigzc, trigzs:  cos(nv), sin(nv) 2D arrays

!

!     Prepare m,n arrays for stacking magnetic field coefficients

      ALLOCATE (nval_bfld(-ntorb:ntorb), mval_bfld(0:mpolb))

      IF (nvalsb(1) .gt. -bigint) THEN
         IF (lscreen) THEN
         PRINT *,' This is an old-style input file.'
         PRINT *,' Convert to new-style by eliminating nvalsb array'
         PRINT *,' and then index borbi array with actual (n,m) indices'
         PRINT *,' Note: n is in units of field period'
         END IF
         nval_bfld = -bigint
         DO l = 1, ntorb
            nval_bfld(l) = nvalsb(l)
         END DO
         mval_bfld(0) = -bigint
         DO l = 1, mpolb
            mval_bfld(l) = l - 1
         END DO
      ELSE
         DO l = -ntorb, ntorb
            nval_bfld(l) = l
         END DO
         DO l = 0, mpolb
            mval_bfld(l) = l
         END DO
      END IF


      mmax = MAX(mmax, MAXVAL(mval_bfld))
      nmax = MAX(nmax, MAXVAL(nval_bfld))
      imask = MINLOC(nval_bfld, nval_bfld > -bigint) - 1 - ntorb
      nn  = MIN(-nmax, nval_bfld(imask(1)))
      nmax = MAX(nmax,ABS(nn))

      ntheta = (3 + meshtz)*(mmax/2 + 1) - 2 - meshtz
      nzeta  = (3 + meshtz)*(nmax/2 + 1) - 2 - meshtz

      ALLOCATE (trigtc(0:mmax,ntheta), trigts(0:mmax,ntheta),
     1          trigzc(-nmax:nmax,nzeta), trigzs(-nmax:nmax,nzeta),
     2          blank(0:mmax,-nmax:nmax,7), overb1(0:mmax,-nmax:nmax),
     3          bsqi1(0:mmax,-nmax:nmax), rmn(mpnt), bmn(mpnt),
     4          trigs(mpnt), trigc(mpnt),
     5          qmn(mpnt), q1mn(mpnt), q2mn(mpnt), dnorm(mpnt), stat=j)
      IF (j .ne. 0) STOP 'Allocation error in ftconv!'


!  initialize to zero; Fourier transform normalization arrays

      nznt = ntheta*nzeta
      dnrm0 = SQRT(one/nznt)
      dnrm = root2*dnrm0

      srces0 = 0;      borbi1 = 0;      overb1 = 0;      bsqi1 = 0
      rmn    = 0;      bmn    = 0;      qmn    = 0
      q1mn   = 0;      q2mn   = 0
      s1cs10 = 0;      dnorm  = dnrm;   blank  = 0
      dnorm(mn0) = dnrm0


!     stack magnetic field Fourier coefficients into 1D array borbi1

      DO nn = -ntorb, ntorb
         mloop: DO mm = 0, mpolb
            DO mn = 1, mpnt
               IF (mval_bfld(mm).eq.mvalue(mn) .and.
     1             nval_bfld(nn).eq.nvalue(mn)) THEN
                  borbi1(mn) = borbi(nn,mm)
                  CYCLE mloop
               ENDIF
            END DO
         END DO mloop
      END DO


c  theta and zeta trig functions over one field period

      twopi = 8*ATAN(one)
      dangt = twopi/ntheta
      dangz = twopi/nzeta

      DO k = 1, ntheta
         trigtc(0,k) = 1
         trigts(0,k) = 0
         DO m = 1, mmax
            ang = dangt*MOD(m*(k-1), ntheta)
            trigtc(m,k) = COS(ang)
            trigts(m,k) = SIN(ang)
         END DO
      END DO
      DO j = 1, nzeta
         trigzc(0,j) = 1
         trigzs(0,j) = 0
         DO n = 1, nmax
            ang = dangz*MOD(n*(j-1), nzeta)
            trigzc(n,j) = COS(ang)
            trigzs(n,j) = SIN(ang)
            trigzc(-n,j) =  trigzc(n,j)
            trigzs(-n,j) = -trigzs(n,j)
         END DO
      END DO

!     form bi1 = 1/B and bisq1 = 1/B**2 cosine series in real space to compute
!     respective Fourier coefficients for 1/B (overb1) and 1/B**2 (bsqi1).
!     need for ALL m = 0, mmax, n= -nmax,nmax, not JUST mpnt space

      DO k = 1, ntheta
         DO j = 1, nzeta
            bb1 = 0
            DO n = -ntorb, ntorb
               np = nval_bfld(n)
               IF (np .le. -bigint) CYCLE
               DO m = 0, mpolb
                  mp = mval_bfld(m)
                  IF (mp .le. -bigint) CYCLE
                  trc = trigtc(mp,k)*trigzc(np,j)
     1                + trigts(mp,k)*trigzs(np,j)
                  bb1 = bb1 + trc*borbi(n,m)
                END DO
            END DO
            IF (bb1 .le. zero)
     1          STOP ' 1/|B| <= 0: Check BORBI array input'
            bi1 = bb1
            IF (ibbi .eq. 1) bi1 = one/bi1
            bisq1 = bi1*bi1
            DO n = -nmax, nmax
               DO m = 0, mmax
                  trc=trigtc(m,k)*trigzc(n,j)+trigts(m,k)*trigzs(n,j)
                  overb1(m,n) = overb1(m,n) + trc*bi1
                  bsqi1(m,n)  = bsqi1(m,n)  + trc*bisq1
               END DO
            END DO
         END DO
      END DO

      overb1(0,-nmax:-1) = 0
      bsqi1(0,-nmax:-1) = 0

      DEALLOCATE (nval_bfld, mval_bfld)

      overb1 = dnrm**2 * overb1
      bsqi1  = dnrm**2 * bsqi1
      overb1(0,0) = overb1(0,0)/2
      bsqi1(0,0)  = bsqi1(0,0)/2
      IF (overb1(0,0) .eq. zero) THEN
         STOP 'Error: |B|(0,0) = 0'
      ELSE
         b00 = one/overb1(0,0)
      END IF

      IF (ibbi .eq. 2) THEN
         fac = specac*bsqi1(0,0)
         WHERE (ABS(bsqi1) < fac) bsqi1 = 0
         fac = specac*overb1(0,0)
         WHERE (ABS(overb1) < fac) overb1 = 0
      END IF

!     compute source S1+(for l=0,2 Legendre pitch harmonics) sine series, Jacobian cosine series,
!     and parallel-streaming and electric (EXB-drift) matrix elements
!     NOTE: S1 source should be multiplied by B00*Rhoi*VTi (VTi=ion thermal speed),
!        and f1 (distribution) by B00*Rhoi to be in proper units. The implied
!        scaling of the transport matrix elements (to get into REAL units) is:
!
!        L11(real) = L11(DKES2) * (B00*Rhoi)**2 * VTi
!        L13(real) = L13(DKES2) * (B00*Rhoi) * VTi
!        L33(real) = L33(DKES2) * VTi
!

      bi0 = overb1(0,0)                                                  !m=0, n=0 1/|B| coefficient
      bi20 = bsqi1(0,0)                                                  !m=0, n=0 1/|B|**2 coefficient
      facjac = btheta*chip + bzeta*psip                                  !jacobian = facjac / |B|**2
      IF (facjac .eq. zero) STOP 'facjac = 0 in ftconv'
      bsqav = one/bi20                                                   !<B**2> in Boozer coordinates
      vp = facjac*bi20                                                   !<<jacobian>>, <<...>> u,v average
      vpnz = vp*nznt
      wtov = bi0*ABS(chip/vp)                                            !transit freq, divided by 1/v
      IF (chip .ne. zero) THEN
         bpfac = bzeta/(chip*bsqav)
      ELSE
         bpfac = bzeta/(EPSILON(chip)*bsqav)
      END IF

      DO mn = 1, mpnt
         m = mvalue(mn);   n = nvalue(mn)
         srces0(mn) = (m*bzeta + nzperiod*n*btheta)
     1               * bsqi1(m,n)/3                                      !l=0 source (sine coefficients): S1
         srces0(mn+mpnt2) = p2*srces0(mn)                                !l=2 source ~ S1
         auxs3m(mn) = overb1(m,n)
         bgrad(mn) = m*chip - nzperiod*n*psip                            !jacobian*(B dot grad)
         exbgrad(mn) = bi20*(m*bzeta + nzperiod*n*btheta)                !Eradial X B term
      END DO

      bi2gi = one/facjac

!     Compute required Fourier transforms

      THETA: DO k = 1, ntheta
         ZETA: DO j = 1, nzeta
            bi1 = 0;  jac = 0;  bgradb = 0;  blebi = 0; qs = 0
            DO n = -nmax, nmax
               DO m = 0, mmax
                  trc=trigtc(m,k)*trigzc(n,j)+trigts(m,k)*trigzs(n,j)
                  trs=trigts(m,k)*trigzc(n,j)-trigtc(m,k)*trigzs(n,j)
                  bi1 = bi1 + overb1(m,n)*trc                            !1/B(nu,nv)
                  jac = jac + bsqi1(m,n)*trc                             !jacobian(nu,nv)
                  bgradb = bgradb -
     1                   (m*chip - nzperiod*n*psip)*overb1(m,n)*trs      !jacobian B dot grad (1/|B|)
                  trs = (m*bzeta + nzperiod*n*btheta)*trs
                  blebi = blebi - overb1(m,n)*trs                        !grad of 1/B in E X B direction
                  qs  = qs - bsqi1(m,n)*trs                              !-source (l=0)
               END DO
            END DO

            jac = jac*facjac;  qs = qs/3;  blebi = blebi*bi20
            b1 = one/bi1                      !|B|(nu,nv)
            jaci = one/jac                    !1/jacobian(nu,nv)
            bigi = bi1*jaci                   !1/(jac*B)(nu,nv)

            gradbgi = bgradb*jaci             !B dot grad (1/|B|))
            bgbi2gi = bgradb*gradbgi          ![jacobian * B dot grad (1/|B|)]**2 / jacobian
            bgbibigi = bgradb*bigi            ![jacobian * B dot grad (1/|B|)] / (jacobian |B|)
            blbi = bgradb*b1                  !|B| (jac * B dot grad (1/|B|))

            qbi = qs*bigi                     !-source /(B*jac)
            qs = qs*jaci                      !-source /jac
            qlbi = qs*bgradb                  !-source * B dot grad (1/|B|)
            blebi = blebi*b1*b1


!  source terms {sigma+(i), C[-1] sigma+(j): s1cs1

            s1cs10 = s1cs10 - jac*qs**2       !-jac*(source/jac)**2

!  series for parallel stress, S3+(l=1,2), S1-(l=0,1,2,3), and
!  S3-(l=0,1,2,3)

            DO mn = 1, mpnt
               n = nvalue(mn)
               m = mvalue(mn)
               trigc(mn) = trigtc(m,k)*trigzc(n,j)   !cos(mu-nv)
     1                   + trigts(m,k)*trigzs(n,j)
               trigs(mn) = trigts(m,k)*trigzc(n,j)   !sin(mu-nv)
     1                   - trigtc(m,k)*trigzs(n,j)
            END DO

            rmn(:) = rmn(:) + trigs(:)*blbi   !rmn(mn): fourier coefficients stress tensor multiplier
            bmn(:) = bmn(:) + trigs(:)*blebi
            qmn(:) = qmn(:) + trigs(:)*qs
            q1mn(:) = q1mn(:) + trigs(:)*qbi
            q2mn(:) = q2mn(:) + trigc(:)*qlbi

!  cosine and sine integrals for matrix elements

            DO n = -nmax, nmax
               DO m = 0, mmax
                  trc=trigtc(m,k)*trigzc(n,j)+trigts(m,k)*trigzs(n,j)   !cos(mu-nv)
                  trs=trigts(m,k)*trigzc(n,j)-trigtc(m,k)*trigzs(n,j)   !sin(mu-nv)
                  blank(m,n,1) = blank(m,n,1) + trc*jac
                  blank(m,n,2) = blank(m,n,2) + trc*jaci
                  blank(m,n,3) = blank(m,n,3) + trc*bigi
                  blank(m,n,4) = blank(m,n,4) + trc*bgbi2gi
                  blank(m,n,5) = blank(m,n,5) + trs*gradbgi
                  blank(m,n,6) = blank(m,n,6) + trs*bgbibigi
                  blank(m,n,7) = blank(m,n,7) + trc*bi1
               END DO
            END DO

         END DO ZETA
      END DO THETA

      DEALLOCATE (trigs, trigc, trigtc, trigts, trigzc, trigzs)

      rmn = rmn*dnorm;   bmn = bmn*dnorm;   qmn = qmn*dnorm
      q1mn = q1mn*dnorm; q2mn = q2mn*dnorm

      dnorm = rthalf * dnorm

      l = 0
      DO mnp = 1, mpnt
         mp = mvalue(mnp)
         np = nvalue(mnp)
         DO mn = 1, mpnt
            facc = dnorm(mn)*dnorm(mnp)
            ms = mvalue(mn) + mp
            ns = nvalue(mn) + np
            md = mvalue(mn) - mp
            nd = nvalue(mn) - np
            is = 1
            IF (md < 0) is = -1
            md = is*md
            nd = is*nd

c  matrix elements <ei M ej> WHERE ei, ej are e+ = SIN, e- = COS
c  symmetric sine-sine and cosine-cosine matrix elements

            IF (mn .le. mnp) THEN
               matjac(mn,mnp,1) = facc*(blank(md,nd,1)-blank(ms,ns,1))   !<e+ M e+>, M = jacobian
               matjac(mn,mnp,2) = facc*(blank(md,nd,1)+blank(ms,ns,1))   !<e- M e->
               bmat2(mn,mnp,1)  = facc*(blank(md,nd,4)-blank(ms,ns,4))   !<e+ M e+>, M = (bgrad(1/B))**2/jacobian
               bmat2(mn,mnp,2)  = facc*(blank(md,nd,4)+blank(ms,ns,4))   !<e- M e->

               bmat3(mn,mnp,1) = facc*(blank(md,nd,2)+blank(ms,ns,2))    !<e- M e->, M = 1/jacobian
               bmat3(mn,mnp,2) = facc*(blank(md,nd,2)-blank(ms,ns,2))    !<e+ M e+>
               bmat5(mn,mnp,1) = facc*(blank(md,nd,3)+blank(ms,ns,3))    !<e- M e->, M = 1/(jacobian*B)
               bmat5(mn,mnp,2) = facc*(blank(md,nd,3)-blank(ms,ns,3))    !<e+ M e+>

               diagl(mn,mnp,1) = facc*(blank(md,nd,7)-blank(ms,ns,7))
               diagl(mn,mnp,2) = facc*(blank(md,nd,7)+blank(ms,ns,7))
            ENDIF

c  sine-cosine matrix elements

            bmat4(mn,mnp,1) = facc*(is*blank(md,nd,6) + blank(ms,ns,6))
            bmat6(mn,mnp,1) = facc*(is*blank(md,nd,5) + blank(ms,ns,5))

         END DO
      END DO

      bmat4(:,:,2) = -TRANSPOSE(bmat4(:,:,1))
      bmat6(:,:,2) = -TRANSPOSE(bmat6(:,:,1))

      DEALLOCATE (blank)

c  complete integrals, sources, and fill out symmetric matrix elements
      c2 = one/(rt3o2**2*vpnz)
      s1cs10 = (root2*p2)**2*s1cs10/vpnz

      c2 = dnrm0/rt3o2                                                   !rt3o2 comes from norm of p(l=1) source
      auxs3m(:) = auxs3m(:)*facjac*rthalf/rt3o2                           !assumes jac = facjac/B**2
      auxs3m(mn0) = auxs3m(mn0)/rthalf

      bmat1 = 0
      DO mn = 1, mpnt
         bmat1(mn,mn) = bgrad(mn)*bgrad(mn)*bi2gi
      END DO

      DO l = plus, minus
         DO mn = 1, mpnt
            diagl(mn:mpnt,mn,l) = diagl(mn,mn:mpnt,l)
            bmat2(mn:mpnt,mn,l) = bmat2(mn,mn:mpnt,l)
            bmat3(mn:mpnt,mn,l) = bmat3(mn,mn:mpnt,l)
            bmat5(mn:mpnt,mn,l) = bmat5(mn,mn:mpnt,l)
            matjac(mn:mpnt,mn,l) = matjac(mn,mn:mpnt,l)
         END DO

         DO mnp = 1, mpnt
            bmat3(:,mnp,l) = exbgrad(:)*exbgrad(mnp)*bmat3(:,mnp,l)
            bmat4(:,mnp,l) = bgrad(mnp)*bmat4(:,mnp,l)
            bmat5(:,mnp,l) = bgrad(:)*exbgrad(mnp)*bmat5(:,mnp,l)
            bmat6(:,mnp,l) = exbgrad(mnp)*bmat6(:,mnp,l)
            diagl(:,mnp,l) = bgrad(mnp) * diagl(:,mnp,l)
         END DO
      END DO


      stfac = -dnrm0/(rt5o2*vp)
      c3 = p2*root2*dnrm0                   !p2 arises from l=2 component of source(+)

      bstrs(:)    = stfac*rmn(:)            !for stell sym, sin(mu-nv) coefficients for stress tensor
!
!     note: these sources are all negative, since they get divided
!     in cescale by -1/nu(l) so the (-) signs cancel..
!
      auxs3p(:,1) = c2*bmn(:)
      auxs3p(:,2) = c2*rmn(:)
      auxs1(:,1)  = c3*(2*q1mn(:)*bgrad(:) + q2mn(:))
      auxs1(:,2)  = c3*qmn(:)*exbgrad(:)
      auxs1(:,3)  = 2*c3*(q1mn(:)*bgrad(:)-2*q2mn(:))


      DEALLOCATE (overb1, rmn, bmn, qmn, q1mn, q2mn, bsqi1, dnorm)

      END SUBROUTINE ftconv

      END MODULE dkes_realspace
