
! ---------------------------------------------------------------------
      subroutine solver_nescoil
!................................................................
c                                                             01/01/89
c     purpose:
c     Set up g and h matrices, and invert g * phi = h to get phi
c     Uses matrix and fourier routines
c     Ref. P. Merkel, Nucl Fus 27, 5, 1987, pg 867-871
c................................................................
c     Modified by PMV and SPH:
c     SVD scan, Xerr target, weighted optimization
c
c ---------------------------------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      use Vmeshes
      use NumParams
      use SvdCtrl
      use OutCtrl, ONLY: w_svd, w_xerr
      use LoopCtrl
      USE Vvacuum3
      USE Vprecal1
      USE Vprecal2
      USE Vprecal3
      USE Vsurface9, ONLY: x1, y1, z1
      USE Vsurfac13, ONLY: snx1, sny1, snz1, dsur1
      USE Vsurfac14, ONLY: bex, bey, bez, ben
      USE Vsolver1
      USE Vdiagno2, ONLY: pot
      use Vbnorm, ONLY: bn_ext
      use Vsurfcurdiag
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, k, mf1, ndim,
     1   nbdim, mfnv, l, m, nfl, n, kl, k1,
     2   mstrta, mstepa, mkeepa, mdspwa, mexp, ifail, nsvdw, icmplxmin,
     3   icurvamax, icdmaxmin, icdavemin, iemnabs, iemnrms, iemnmax,
     4   iwemnabs, iwemnrms, iwemnmax, j, info
      real(rprec) ::
     1   t1, t2, err, stat, epsilon, fla, eabsmin, ermsmin,
     2   emaxmin, weabsmin, wermsmin, wemaxmin, cmplxmin, curvamax,
     3   cdmaxmin, cdavemin, bmodav, bmodavinv, potpotm, potpotmm,
     4   potsq, complexity, errmax, errabs, errrms, werrmax, werrabs,
     5   werrrms, epsm
c................................................................
c   x1,y1,z1 hold co-ordinates of nuvh1 points on half surface
c   snx1,sny1,snz1 hold normal vector and dsur1 holds surface element
c   bex1,bey1,bez1 hold b field and ben holds its normal component
c   a(nmax) and work vectors are used in matrix and fourier
c................................................................

      ndim  = nf+mf*(2*nf+1)
      nbdim = 2*(mf+1)*(2*nf+1)

      i=0
      if (.not.allocated(a))
     1 allocate (a(nmax), work(nmax), ab(ndim,ndim+1), result(ndim),
     1          bfx(nuvh1,nbdim), bfy(nuvh1,nbdim), bfz(nuvh1,nbdim),
     2          bfn(nuvh1,nbdim), bex(nuvh1), bey(nuvh1), bez(nuvh1),
     3          bpx(nuvh1), bpy(nuvh1), bpz(nuvh1), dumv(nuvh1),
     4          ben(nuvh1), daccu(nuvh1), bmag(nuvh1),
     5        pot(0:mf,-nf:nf), svdv(ndim,ndim), svdw(ndim), stat=i)
      if (i .ne. 0) stop 'allocation error in nescoil solver'

c     Zero all allocated arrays
      a(:) = zero
      work(:) = zero
      result(:) = zero
      bex(:) = zero
      bey(:) = zero
      bez(:) = zero
      bpx(:) = zero
      bpy(:) = zero
      bpz(:) = zero
      dumv(:) = zero
      ben(:) = zero
      daccu(:) = zero
      bmag(:) = zero
      svdw(:) = zero

      bfx(:nuvh1,:nbdim) = zero
      bfy(:nuvh1,:nbdim) = zero
      bfz(:nuvh1,:nbdim) = zero
      bfn(:nuvh1,:nbdim) = zero
      ab(:ndim,:ndim+1) = zero

      pot(0:mf,-nf:nf) = zero
      svdv(:ndim,:ndim) = zero
c.........................
c     Either zero or get background b field at all points  depending on
      if (ibex .ne. 0) then
         do i = 1, nuvh1
            call bexter (x1(i), y1(i), z1(i), bex(i), bey(i), bez(i))
         end do
      endif

c.........................
      write (inesc, 10) '---- Calling MATRIX+FOURIER ----'
      call second0(t1)
c.........................
      do k = 1, nuvh1

c.....   Zero a(i) vector every time
         a(:nmax) = zero

c        Calculate g(u',v', u,v) and h(u,v) of Eqs. 7 and 8
         call matrix (a, k)

c        Integral over u',v' of {g(u v, u' v') * sin(m u'+ n v')} in eq. 6
         call fourier (3)

c        Set the info from a vector into bfx,bfy,bfz
         mf1 = 2*(mf + 1)
         mfnv = mf1*nv
         l = 0

         do m = 0, mf
            nfl = -nf
            if (m .eq. 0) nfl = 1
            do n = nfl, nf
               i = mf1*n + 2*(m + 1)
               if (n .lt. 0) i = i + mfnv
               l = l + 1
               bfx(k,l) = -a(i       )*fnv
               bfy(k,l) = -a(i+  mfnv)*fnv
               bfz(k,l) = -a(i+2*mfnv)*fnv
               bfn(k,l) = snx1(k)*bfx(k,l) + sny1(k)*bfy(k,l)
     1                   + snz1(k)*bfz(k,l)
            end do
         end do

      end do
c.....End of big loop over index k

      call second0(t2)
      write (inesc,"('Time in MATRIX+FOURIER: ',g12.3,' sec')") t2-t1
c.........................

c.....Calculate bnormal from bexternal dot normal
         ben(:nuvh1) = snx1(:nuvh1)*bex(:nuvh1) + sny1(:nuvh1)*
     1      bey(:nuvh1) + snz1(:nuvh1)*bez(:nuvh1) + bn_ext(:nuvh1)
c
c     At this point, g_m_n(u,v) and h(u,v) are known:
c     bfn(i,l) is g_m_n(u,v) (Green's function) and
c     ben(i)   is h(u,v) (External B + homogeneous part from Ip,It)
c
c     bfn(p,mn) * phi(mn) gives bnorm field at plasma point p
c         due to phi(mn) of coilcurrent, and
c     ben(p) is external+vmec normal field on plasma surface
c
c     Do not change bfn, ben if you want to use them later
c     If you do change them temporarily, make sure
c     you change them back.

c*****************************************
c
c     Begin Berror or Xerror targets modifications
c                     P. Valanju, (PMV) Mar 99
         mstrta = abs(mstrt)
         mstepa = abs(mstep)
         mkeepa = abs(mkeep)
c        Extract correct exponent of dsur from MDSPW
         if (mdspw == 0) mdspw = 3        !Default gives dsur to power 1
         mdspwa = abs(mdspw)
         mexp = mdspwa - 2                !exponent in dsur**mexp
         if (mdspwa >= 2) then
            write (inesc, '(a,i4)') 'Using fac*dsur**', mexp
         else
            write (inesc, 10) 'Using no fac or dsur**mexp'
         endif
c        Make sure LSQ branch is chosen if curr density is targeted
         if (curwt > 0) mstep = -iabs(mstep)
c.................................................................
c      NESCOIL control settings explained:
c
c 1. MSTRT = Method + svdscan start if >1
c            >=0:Berr, <0:Xerr, unless MSTEP <=0: Least square
c 2. MSTEP = Method + svdscan stepsize:
c            <=0: LeastSquare, =0: use old f04abe (now SOLVER), no svd
c 3. MKEEP = svd/scan control: =0: svdscan, else keep |nkeep| wgts
c            <0: write all weights to output
c 4. MDSPW = 2+exponent of dsur multiplying bfn,ben:
c            <0: post-calculate Xerr svdscan and write to output
c 5. CURWT = Weight for surface current minimization
c            Works ONLY in LSQ branch
c 6. TRGWT : Not implemented yet (PMV)
c           = 1 : use uniform weight,
c             else Weighted (over plasma surface) error minimization
c           = wtmn: file contains mn coeffs of weight
c           = wtuv: file contains uv values of weights
c......................................
c            MSTRT  MKEEP  MSTEP
c      LSQ     -      -      0    ->  LSQ  with f04abe (solver)
c              -      n     -s    ->  LSQ  with svd=n weights
c             m>1   +-n     -s    ->  LSQ  with svdscan m,n,s
c
c      Berr    1      n      s    ->  Berr with svd=n weights
c             m>1   +-n      s    ->  Berr with svdscan
c
c      Berr   -1      n      s    ->  Berr with svd=n weights
c            -m>1   +-n      s    ->  Berr with svdscan
c.................................................................
c      Note: LSQ cannot run with Xerr since that requires many extra
c            inverse fourier transforms inside min_xerr_phimn.
c            We can add that if needed later, but
c            Solution for phi(m,n) directly in fourier space is faster.
c......................................

c      Note about MDSPW: (controls factor(i)*dsur**mexp) :
c      We calculate all things over half+1 grid points nuvh1
c      This is ok for all nuv points due to stellrator symmetry.
c      What this means is when we want to sum anything over whole fs,
c      we MUST USE factor(i) which weighs the edge points correctly
c      in a trapezoidal sum. This is true EVEN IF WE DO NOT USE DSUR.
c
c      Due to this, the matrix in the linear (svd/Berr or Xerr) case:
c      bfn(:nuvh1,ndim)* phi(ndim) = ben(:nuvh1)
c      should be modified to
c       bfn(:nuvh1,:)*sqrt(factor(:nuvh1)) and
c       ben(:nuvh1)*sqrt(factor(:nuvh1))
c         if RMS ERROR (root of sum of square err) is minimized,
c      or should be modified to
c       bfn(:nuvh1,:)*factor(:nuvh1) and
c       ben(:nuvh1)*factor(:nuvh1)
c         if ABS ERROR (sum of |err|) is minimized.
c
c      We will always minimize RMS ERROR, as in original nescoil.
c      This is what Least-sq method does, so Linear will do same.
c      However, old accuracy reported ABS ERROR, which was not
c      eactly what was minimized.
c      We will report both.
c
c      This would make NO difference if bfn were a square matrix so
c      the error would be zero, but it DOES make a difference when
c      the there are more plasma points than fourier modes in phi
c      so that an exact fit is not posiible. This is similar to
c      a weighted least-square minimization.

c------------------------------------------------------
         if (mstep<=0 .or. mstrt==0) then
c        The least-square nescoil branch using f04abe (solver) or svd:
            call second0(t1)

c......     Calculate ab matrix:
c           use MDSPWA  to control power of dsur,
c             CURWT   to control amount of J minimization
            call ab_matrix (ab, bfn, dumv, ndim, mdspwa, curwt)

c......
            if (mstep == 0) then
c           Solve using u-v solver fo4abe as in old nescoil:
c           This is good choice since no svd or svdscan is done
               write (inesc, 10) 'Case: Standard Nescoil LSQ solver'

               result(:ndim) = ab(:ndim,ndim+1)
               call solver(ab, result, ndim, 1, info)    !Do not use old f04abe solver

               call second0(t2)
               write (inesc,
     1           "('Berr Least SQ took ',g12.3,' sec')") t2-t1
c              There will be no svdscan, since svd is not used here
               go to 123            !Only this one avoids svd or svdscan

            else                                 !MSTEP not 0
c           Solve least-sq matrix eq using svd:
c           MSTRT, MKEEP, MSTEP will be used later for svdscan
               write (inesc, 10) 'Case: Standard Nescoil LSQ with svd'
               call svd_solve (ndim, ndim, ndim, ndim, ab,
     1             ab(1,ndim+1), svdv, svdw, nsvdw)
               call second0(t2)
               write (inesc,
     1             "('Berr LSQ SVD took ',g12.3,' sec')") t2-t1
            endif

         else
c        The Linear branch with Berr or Xerr:
c................................
            if (mstrt > 0) then
c           Target Berr, use svd to solve problem:
c           MSTRT, MKEEP, MSTEP will be used later for svdscan
               call second0(t1)
               write (inesc,10) 'Case: Targeting Berror Linear with svd'

c              Modify bfn, ben by proper factors of dsur and fac
c              so they agree with what is done in lsq method
c              Temporarily use dumv for this here and later
c              when un-modifying bfn, ben
               if (mdspwa >= 2) then
                  dumv(:nuvh1)=sqrt(factor(:nuvh1)*dsur1(:nuvh1)**mexp)
                  ben(:nuvh1) =ben(:nuvh1)*dumv(:nuvh1)
                  do i = 1, ndim
                     bfn(:nuvh1,i) = bfn(:nuvh1,i)*dumv(:nuvh1)
                  end do
               endif

c              Solve the linear svd problem
               call svd_solve (nuvh1, ndim, nuvh1, nbdim, bfn, -ben,
     1            svdv, svdw, nsvdw)

c              UnModify bfn, ben by proper factors of dsur and fac
c              so they can be used again
c              Use dumv calculated earlier
               if (mdspwa >= 2) then
                  ben(:nuvh1) = ben(:nuvh1)/dumv(:nuvh1)
                  do i = 1, ndim
                     bfn(:nuvh1,i) = bfn(:nuvh1,i)/dumv(:nuvh1)
                  end do
               endif
c              Bfn and Ben are again B fields at plasma points,
c              not b*ds. So they can be used to calc errors etc.

               call second0(t2)
               write (inesc,
     1             "('Berr Linear SVD took ',g12.3,' sec')") t2-t1
            endif

c................................
            if (mstrt < 0) then
c           Target Xerr, use svd to solve problem:
c           Note: we cannot use non-uniform weights in this case,
c                 since calculating xerr itself corresponds to using
c                 non-uniform weights on plasma surface
               call second0(t1)
               write (inesc,10) 'Case: Targeting Xerror Linear with svd'
               epsilon = 0

               call min_xerr_phimn (1, epsilon, nu1, nv1, nuvh1,
     1            nbdim, ndim, bfn, ben, dsur1, svdv, svdw, nsvdw)
               call second0(t2)
               write (inesc,
     1             "('Xerr Linear SVD took ',g12.3,' sec')") t2-t1
            endif

         endif

c------------------------------------------------------
c        At this point, the solution has been found (in some way)
c        All that remains is reporting the errors
c......................................................
c        Report svd errors, single if MSTRTA=1, scan if MSTRTA > 1:

c        Limit MKEEPA, use all weights if set to zero:
         if (mkeepa == 0) mkeepa = nsvdw
         mkeepa = min(nsvdw,mkeepa)

c        All branches arrive here with array of answers in svdv
c        At this point, we can do one of 3 three things:
c        1. Scan all svd weights and calculate err, jmax, complexity ...
c        2. Just use one of the weights (MKEEPA) specified by user
c        3. Find the optimum number of weights (isvdw) to keep
c           using ?? some ?? criterion and solutions svdv

         if (mstrta > 1) then
c        SVDSCAN: Find error, jmax, complexity, and anything else
c        for i=method,MKEEPA svd wgts kept

c....................................................................
c        Write all titles for scan (include identifiers for grep):
            write (inesc, '(4(a,i4),a)')
     1          '---- Svdscan from ',mstrta,' to ',mkeepa,' at ',
     2          mstepa,' steps out of ',ndim,' total wgts ----'

c           do not call accuracy later in nescoil (main) at each wght
            noaccuracy = 1

c           Calculate total area of plasma surface
            fla = sum(dsur1(:nu1)) + 2.*sum(dsur1(nu1+1:nuvh1-nu1))
     1          + sum(dsur1(nuvh1-nu1+1:nuvh1))

c           Set starting values for min/max error calculation in svdloop
            eabsmin = 1.0e14_dp                  !for min values of various errors
            ermsmin = 1.0e14_dp                  !min of rms err
            emaxmin = 1.0e14_dp                  !min of max err
            weabsmin = 1.0e14_dp                 !min of abs err
            wermsmin = 1.0e14_dp                 !weighted errors
            wemaxmin = 1.0e14_dp
            cmplxmin = 1.0e14_dp                 !min of complexity
            curvamax = zero                      !max of current curvature radius
            cdmaxmin = 1.0e14_dp                 !min of current density max
            cdavemin = 1.0e14_dp                 !min of current density ave

c...svdscan loop........................................................
            do i = mkeepa, mstrta, -mstepa     !scan in decreasing order
c           Go over all weights from best to worst, because
c           the calculation of bmod is done only first time, and
c           the best weight will give the best estimate for it.

c......        Calculate Bmag using green's functions rather than besucu
c              This is much faster since green's functions already known
               Bmag(:nuvh1) = sqrt(
     1         (bex(:nuvh1)+Matmul(bfx(:nuvh1,:ndim),svdv(:ndim,i)))**2
     2        +(bey(:nuvh1)+Matmul(bfy(:nuvh1,:ndim),svdv(:ndim,i)))**2
     3        +(bez(:nuvh1)+Matmul(bfz(:nuvh1,:ndim),svdv(:ndim,i)))**2)
c              Calculate fs average of Bmag
               dumv(:nuvh1) = dsur1(:nuvh1) * Bmag(:nuvh1)
               bmodav = ( sum(dumv(:nu1))
     1        + 2*sum(dumv(nu1+1:nuvh1-nu1))
     2        + sum(dumv(nuvh1-nu1+1:nuvh1)) )/ fla
               bmodavinv = one/bmodav

c.......       Calculate and write the maximum current density
c              Put svdv into phi(m,n) for sucudiag
c              Also calculate 1,2 moments for coil complexity
c              Save all pot(m,n) for later use by postprocessor.
               l = 0
               m = 0
               potpotm = zero       !for moments of pot(m,n) for complexity
               potpotmm = zero

c              Set and Write pot coeffs for this svd weight
               write (inesc, '(a,i4,a)') '---- Phi(m,n) for isvd = ',
     1            i,'----'
               do n = 1, nf
                  l = l + 1
                  pot(0,n)  = -svdv(l, i)
                  pot(0,-n) =  svdv(l, i)
               end do
               do n = -nf, nf
                  write (inesc,"(2i6,g25.16)") 0, n, pot(0,n)/two
               end do
               do m = 1, mf
                  do n = -nf, nf
                     l = l + 1
                     pot(m,n) = -svdv(l,i)
                     write (inesc,"(2i6,g25.16)") m, n, pot(m,n)
                     potsq = pot(m,n) * pot(m,n) * m
                     potpotm = potpotm + potsq
                     potpotmm = potpotmm + potsq*m
                  end do
               end do
               write (inesc, 10) '---- end Phi(m,n) ----'

               complexity = potpotmm/potpotm
               if (complexity < cmplxmin) then
                  cmplxmin = complexity
                  icmplxmin = i
               endif

               call surfcur_diag

               if (curvamax < curvra) then
                  curvamax = curvra
                  icurvamax = i
               endif
               if (cdmaxmin > cudema) then
                  cdmaxmin = cudema
                  icdmaxmin = i
               endif
               if (cdavemin > avcude) then
                  cdavemin = avcude
                  icdavemin = i
               endif

c......        Calculate errors for bnorm with <bmodavinv>, i.e.:
c              <|berr|>/<Bmod>, sqrt<|berr|^2>/<Bmod>, and
c              max<|berr|>/<Bmod>, where <> is fs average
c
c              1. Calculate unweighted residual |BFN*PHI+BEN|
               daccu(:nuvh1) = abs( ben(:nuvh1) +
     1            Matmul(bfn(:nuvh1,:ndim),svdv(:ndim,i)) )
c              2. Calc fs averages and max of residual
               dumv(:nuvh1) = daccu(:nuvh1)   !temporary storage
               call calc_percent_errs (dumv, dsur1, fla, nu1, nuvh1,
     1            errmax, errabs, errrms)
c              3. Calc globally weighted errors (with <Bmag>)
               errabs = errabs*bmodavinv
               errrms = errrms*bmodavinv
               errmax = errmax*bmodavinv

c.......       Calculate errors for locally wgted bnorm/bmod:
c              <|berr|/Bmod>, sqrt<|berr/Bmod|^2>>, and
c              max<|berr/Bmod|>, where <> is fs average
c
c              1. Calculate locally weighted residual |BFN*PHI-BEN|/Bmag
               dumv(:nuvh1) = daccu(:nuvh1)/bmag(:nuvh1)
c              2. Calc fs averages and max of locally weighted errors
               call calc_percent_errs (dumv, dsur1, fla, nu1, nuvh1,
     1            werrmax, werrabs, werrrms)

c.......       Write errors and j to output
               write (inesc, 10)
     1            'Info for isvd, Berrors( abs, rms, max):'
               write (inesc, '(a,i4,1p3e10.3)')
     1             '<|g*phi+h|> ',i, errabs, errrms, errmax
               write (inesc, '(a,i4,1p3e10.3)')
     1             '<|g*phi+h|/|B|> ', i, werrabs, werrrms, werrmax
               write (inesc, 10)
     1             'isvd, Jmax, Jave, J curv rad, Complexity'
               write (inesc, '(a,i4,1p4e10.3)')
     1             'J ',i, cudema, avcude, curvra, complexity

c.......       Find svd num and value of various min errors
               if (errabs < eabsmin) then
                  eabsmin = errabs
                  iemnabs = i
               endif
               if (errrms < ermsmin) then
                  ermsmin = errrms
                  iemnrms = i
               endif
               if (errmax < emaxmin) then
                  emaxmin = errmax
                  iemnmax = i
               endif

c              Find svd num and value of various weighted min errors
               if (werrabs < weabsmin) then
                  weabsmin = werrabs
                  iwemnabs = i
               endif
               if (werrrms < wermsmin) then
                  wermsmin = werrrms
                  iwemnrms = i
               endif
               if (werrmax < wemaxmin) then
                  wemaxmin = werrmax
                  iwemnmax = i
               endif

               write (inesc, '(a,i4,a)')
     1            '---- end info for isvd = ', i,' ----'
            end do
c...end of svdscan loop.....................................

            write (inesc, 10)
     1          '---- %err_U Minimum (Bmod at local) Table ----'
            write (inesc, 20) 'Abs_U ', iemnabs, ' : ', eabsmin
            write (inesc, 20) 'Rms_U ', iemnrms, ' : ', ermsmin
            write (inesc, 20) 'Max_U ', iemnmax, ' : ', emaxmin
            write (inesc, 10)
     1          '---- %err_W Minimum (Bmod at global) Table: ----'
            write (inesc, 20) 'Abs_W ', iwemnabs, ' : ', weabsmin
            write (inesc, 20) 'Rms_W ', iwemnrms, ' : ', wermsmin
            write (inesc, 20) 'Max_W ', iwemnmax, ' : ', wemaxmin
            write (inesc, 10) '---- Current density Minimum Table: ----'
            write (inesc, 20) 'Jave ', icdavemin, ' : ', cdavemin
            write (inesc, 20) 'max ', icdmaxmin, ' : ', cdmaxmin
            write (inesc, 20) 'JRad MAX at ', icurvamax, ' : ', curvamax
            write (inesc, 20) 'Complexity Min ',icmplxmin,' : ',cmplxmin

         endif                                   !End of svdscan if

c.....   Special cases to write various things:
         if( w_svd > 1 .or. w_svd == - 2 ) then
c        write all svd weights to output
            m = abs(mkeep)
            write (inesc, '(2(a,i4),a)') '---- svd weights: ', m,
     1        ' out of ', nsvdw,' ----'
            write (inesc, *) (svdw(i),i=1,m)
            write (inesc, 10) '---- end svd weights ----'
         endif

         write (inesc, 10) '---- end SVD Scan ----'

c.....   Calc Xerr with full svd scan from given solutions svdv
c        write all xerrors to output
         if (mdspw<0 .and. mstrt>=0 .and. mstep/=0) call
     1      min_xerr_phimn (0, epsilon, nu1, nv1, nuvh1, nbdim, ndim
     2      , bfn, ben, dsur1, svdv, svdw, nsvdw)

c....................................................
c        Always fill "result" array for postprocessing with MKEEPA wgts
         result(:ndim) = svdv(1,mkeepa:ndim-1+mkeepa)

c        Throw away svd arrays now
  123    continue
c..................................................................
c        All branches end up here with "result"
c        which is then used by the standard nescoil same way as always

c*********** End of Berr and Xerr modifications PMV ***********

c        Postprocessing:
c        All braches (original, Berr, Xerr) arrive here with "result"
c        By this time the phi(m,n) is known in array result(ndim),
c        all that needs to be done is unpack it into phi(m,n) etc

c        first zero all fields b
         bpx(:nuvh1) = zero
         bpy(:nuvh1) = zero
         bpz(:nuvh1) = zero
c...................................
c        Calculate b fields at all points from phi_m_n, i.e., pot(m,n)
         do l = 1, ndim
            bpx(:nuvh1) = bpx(:nuvh1) + bfx(:nuvh1,l)*result(l)
            bpy(:nuvh1) = bpy(:nuvh1) + bfy(:nuvh1,l)*result(l)
            bpz(:nuvh1) = bpz(:nuvh1) + bfz(:nuvh1,l)*result(l)
         end do
c..........................
c        First set pot(m,n) for m=0, n=1,nf. Array index l is just n here
         l = 0
         m = 0
         pot(0,1:nf) = -result(1:nf)
         pot(0,-1:(-nf):(-1)) = result(1:nf)
         l = nf + l
c.............................
c        Next set pot(m,n) for m=1,mf and n=1,nf. Array index l keeps incrementing
         do m = 1, mf
            pot(m,(-nf):nf) = -result(l+1:l+nf*2+1)
            l = l + nf*2 + 1
         end do
C.................................
C    Write pot.coeffs to output file if this was NOT an svdscan
         if (mstrta < 1) then
            potpotm = zero       !for moments of pot(m,n) for complexity
            potpotmm = zero
            write (inesc, 10) '---- Phi(m,n) for least squares ---'
            do m = 0, mf
               epsm = one
               if (m == 0) epsm = 0.5_dp
               do n = -nf, nf
                  write (inesc,'(i3,2x,i3,2x,g25.16)') m,n,epsm*pot(m,n)
                  potsq = pot(m,n) * pot(m,n) * m
                  potpotm = potpotm + potsq
                  potpotmm = potpotmm + potsq*m
               end do
            end do
            write (inesc, 10) '---- end Phi(m,n) for least squares ----'
            complexity = -1
            if (potpotm .ne. zero) complexity = potpotmm/potpotm
            write (inesc, '(a,1pe16.8)')'Complexity = ', complexity
         endif

      if (iloop .le. 0)
     1  deallocate (a, work, svdv, svdw, ab, result, bfx, bfy, bfz,
     1   dumv, bpx, bpy, bpz, daccu, bmag, bfn, bex, bey, bez, ben)

 10   format (a)
 20   format (a,i4,a,1pe16.8)

      end subroutine solver_nescoil
