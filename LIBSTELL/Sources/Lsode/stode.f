      SUBROUTINE stode(neq, y, yh, nyh, yh1, ewt, savf, acor, wm, iwm,
     1    f, jac, pjac, slvs)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER nyh, jac
      INTEGER, DIMENSION(*) :: neq, iwm
      REAL(rprec), DIMENSION(*) :: y
      REAL(rprec), DIMENSION(nyh,*) :: yh
      REAL(rprec), DIMENSION(*) :: yh1, ewt, savf, acor, wm
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
C...  /ls0001/
      common /ls0001/ conit, crate, el(13), elco(13,12), hold, rmax,
     1   tesco(3,12), ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     2   iownd(14), ialth, ipup, lmax, meo, nqnyh, nslp, icf, ierpj,
     3   iersl, jcur, jstart, kflag, l, meth, miter, maxord, maxcor,
     4   msbp, mxncf, n, nq, nst, nfe, nje, nqu
      INTEGER   iownd, ialth, ipup, lmax, meo, nqnyh, nslp, icf, ierpj,
     1   iersl, jcur, jstart, kflag, l, meth, miter, maxord, maxcor,
     2   msbp, mxncf, n, nq, nst, nfe, nje, nqu
      REAL(rprec) ::conit, crate, el, elco, hold, rmax, tesco,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: one = 1
      INTEGER :: i1, iredo, iret, j, jb, m, ncf, newq
      REAL(rprec) :: dcon, ddn, del, delp, dsm, dup, exdn, exsm,
     1   exup, r, rh, rhdn, rhsm, rhup, told
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      REAL(rprec) , EXTERNAL :: vnorm
      EXTERNAL f, jac, pjac, slvs
C-----------------------------------------------
clll. optimize
c-----------------------------------------------------------------------
c stode performs one step of the integration of an initial value
c problem for a system of ordinary differential equations.
c note.. stode is independent of the value of the iteration method
c indicator miter, when this is .ne. 0, and hence is independent
c of the TYPE of chord method used, or the jacobian structure.
c communication with stode is done with the following variables..
c
c neq    = integer array containing problem size in neq(1), and
c          passed as the neq argument in all calls to f and jac.
c y      = an array of length .ge. n used as the y argument in
c          all calls to f and jac.
c yh     = an nyh by lmax array containing the dependent variables
c          and their approximate scaled derivatives, where
c          lmax = maxord + 1.  yh(i,j+1) contains the approximate
c          j-th derivative of y(i), scaled by h**j/factorial(j)
c          (j = 0,1,...,nq).  on entry for the first step, the first
c          two columns of yh must be set from the initial values.
c nyh    = a constant integer .ge. n, the first dimension of yh.
c yh1    = a one-dimensional array occupying the same space as yh.
c ewt    = an array of length n containing multiplicative weights
c          for local error measurements.  local errors in y(i) are
c          compared to 1.0/ewt(i) in various error tests.
c savf   = an array of working storage, of length n.
c          also used for input of yh(*,maxord+2) when jstart = -1
c          and maxord .lt. the current order nq.
c acor   = a work array of length n, used for the accumulated
c          corrections.  on a successful return, acor(i) contains
c          the estimated one-step local error in y(i).
c wm,iwm = real and integer work arrays associated with matrix
c          operations in chord iteration (miter .ne. 0).
c pjac   = name of routine to evaluate and preprocess jacobian matrix
c          and p = i - h*el0*jac, IF a chord method is being used.
c slvs   = name of routine to solve linear system in chord iteration.
c ccmax  = maximum relative change in h*el0 before pjac is called.
c h      = the step size to be attempted on the next step.
c          h is altered by the error control algorithm during the
c          problem.  h can be either positive or negative, but its
c          sign must remain constant throughout the problem.
c hmin   = the minimum absolute value of the step size h to be used.
c hmxi   = inverse of the maximum absolute value of h to be used.
c          hmxi = 0.0 is allowed and corresponds to an infinite hmax.
c          hmin and hmxi may be changed at any time, but will not
c          take effect until the next change of h is considered.
c tn     = the independent variable. tn is updated on each step taken.
c jstart = an integer used for input only, with the following
c          values and meanings..
c               0  perform the first step.
c           .gt.0  take a new step continuing from the last.
c              -1  take the next step with a new value of h, maxord,
c                    n, meth, miter, and/or matrix parameters.
c              -2  take the next step with a new value of h,
c                    but with other inputs unchanged.
c          on return, jstart is set to 1 to facilitate continuation.
c kflag  = a completion code with the following meanings..
c               0  the step was succesful.
c              -1  the requested error could not be achieved.
c              -2  corrector convergence could not be achieved.
c              -3  fatal error in pjac or slvs.
c          a return with kflag = -1 or -2 means either
c          abs(h) = hmin or 10 consecutive failures occurred.
c          on a return with kflag negative, the values of tn and
c          the yh array are as of the beginning of the last
c          step, and h is the last step size attempted.
c maxord = the maximum order of integration method to be allowed.
c maxcor = the maximum number of corrector iterations allowed.
c msbp   = maximum number of steps between pjac calls (miter .gt. 0).
c mxncf  = maximum number of convergence failures allowed.
c meth/miter = the method flags.  see description in driver.
c n      = the number of first-order differential equations.
c-----------------------------------------------------------------------
      kflag = 0
      told = tn
      ncf = 0
      ierpj = 0
      iersl = 0
      jcur = 0
      icf = 0
      delp = 0
      IF (jstart > 0) GOTO 200
      IF (jstart /= (-1)) THEN
         IF (jstart == (-2)) GOTO 160
c-----------------------------------------------------------------------
c on the first CALL, the order is set to 1, and other variables are
c initialized.  rmax is the maximum ratio by which h can be increased
c in a single step.  it is initially 1.e4 to compensate for the small
c initial h, but THEN is normally equal to 10.  IF a failure
c occurs (in corrector convergence or error test), rmax is set at 2
c for the next increase.
c-----------------------------------------------------------------------
         lmax = maxord + 1
         nq = 1
         l = 2
         ialth = 2
         rmax = 10000
         rc = 0
         el0 = 1
         crate = 0.7_dp
         hold = h
         meo = meth
         nslp = 0
         ipup = miter
         iret = 3
      ELSE
c-----------------------------------------------------------------------
c the following block handles preliminaries needed when jstart = -1.
c ipup is set to miter to force a matrix update.
c IF an order increase is about to be considered (ialth = 1),
c ialth is reset to 2 to postpone consideration one more step.
c IF the caller has changed meth, cfode is called to reset
c the coefficients of the method.
c IF the caller has changed maxord to a value less than the current
c order nq, nq is reduced to maxord, and a new h chosen accordingly.
c IF h is to be changed, yh must be rescaled.
c IF h or meth is being changed, ialth is reset to l = nq + 1
c to prevent further changes in h for that many steps.
c-----------------------------------------------------------------------
         ipup = miter
         lmax = maxord + 1
         IF (ialth == 1) ialth = 2
         IF (meth /= meo) THEN
            CALL cfode (meth, elco, tesco)
            meo = meth
            IF (nq > maxord) GOTO 120
            ialth = l
            iret = 1
            GOTO 150
         END IF
         IF (nq <= maxord) GOTO 160
  120    CONTINUE
         nq = maxord
         l = lmax
         el(:l) = elco(:l,nq)
         nqnyh = nq*nyh
         rc = rc*el(1)/el0
         el0 = el(1)
         conit = 0.5_dp/(nq + 2)
         ddn = vnorm(n,savf,ewt)/tesco(1,l)
         exdn = one/(l)
         rhdn = one/(1.3_dp*ddn**exdn + 0.0000013_dp)
         rh = MIN(rhdn,one)
         iredo = 3
         IF (h == hold) GOTO 170
         rh = MIN(rh,ABS(h/hold))
         h = hold
         GOTO 175
c-----------------------------------------------------------------------
c cfode is called to get all the integration coefficients for the
c current meth.  then the el vector and related constants are reset
c whenever the order nq is changed, or at the start of the problem.
c-----------------------------------------------------------------------
      END IF
      CALL cfode (meth, elco, tesco)
  150 CONTINUE
      el(:l) = elco(:l,nq)
      nqnyh = nq*nyh
      rc = rc*el(1)/el0
      el0 = el(1)
      conit = 0.5_dp/(nq + 2)
      GOTO (160,170,200) iret
c-----------------------------------------------------------------------
c IF h is being changed, the h ratio rh is checked against
c rmax, hmin, and hmxi, and the yh array rescaled.  ialth is set to
c l = nq + 1 to prevent a change of h for that many steps, unless
c forced by a convergence or error test failure.
c-----------------------------------------------------------------------
  160 CONTINUE
      IF (h == hold) GOTO 200
      rh = h/hold
      h = hold
      iredo = 3
      GOTO 175
  170 CONTINUE
      rh = MAX(rh,hmin/ABS(h))
  175 CONTINUE
      rh = MIN(rh,rmax)
      rh = rh/MAX(one,ABS(h)*hmxi*rh)
      r = 1
      DO j = 2, l
         r = r*rh
         yh(:n,j) = yh(:n,j)*r
      END DO
      h = h*rh
      rc = rc*rh
      ialth = l
      IF (iredo == 0) GOTO 690
c-----------------------------------------------------------------------
c this section computes the predicted values by effectively
c multiplying the yh array by the pascal triangle matrix.
c rc is the ratio of new to old values of the coefficient  h*el(1).
c when rc differs from 1 by more than ccmax, ipup is set to miter
c to force pjac to be called, IF a jacobian is involved.
c in ANY CASE, pjac is called at least every msbp steps.
c-----------------------------------------------------------------------
  200 CONTINUE
      IF (ABS(rc - 1) > ccmax) ipup = miter
      IF (nst >= nslp + msbp) ipup = miter
      tn = tn + h
      i1 = nqnyh + 1
      DO jb = 1, nq
         i1 = i1 - nyh
cdir$ ivdep
         yh1(i1:nqnyh) = yh1(i1:nqnyh) + yh1(i1+nyh:nqnyh+nyh)
      END DO
c-----------------------------------------------------------------------
c up to maxcor corrector iterations are taken.  a convergence test is
c made on the r.m.s. norm of each correction, weighted by the error
c weight vector ewt.  the SUM of the corrections is accumulated in the
c vector acor(i).  the yh array is not altered in the corrector loop.
c-----------------------------------------------------------------------
  220 CONTINUE
      m = 0
      y(:n) = yh(:n,1)
      CALL f (neq, tn, y, savf)
      nfe = nfe + 1
      IF (ipup <= 0) GOTO 250
c-----------------------------------------------------------------------
c IF indicated, the matrix p = i - h*el(1)*j is reevaluated and
c preprocessed before starting the corrector iteration.  ipup is set
c to 0 as an indicator that this has been done.
c-----------------------------------------------------------------------
      CALL pjac (neq, y, yh, nyh, ewt, acor, savf, wm, iwm, f, jac)
      ipup = 0
      rc = 1
      nslp = nst
      crate = 0.7_dp
      IF (ierpj /= 0) GOTO 430
  250 CONTINUE
      acor(:n) = 0
  270 CONTINUE
      IF (miter == 0) THEN
c-----------------------------------------------------------------------
c in the CASE of functional iteration, update y directly from
c the result of the last function evaluation.
c-----------------------------------------------------------------------
         savf(:n) = h*savf(:n) - yh(:n,2)
         y(:n) = savf(:n) - acor(:n)
         del = vnorm(n,y,ewt)
         y(:n) = yh(:n,1) + el(1)*savf(:n)
         acor(:n) = savf(:n)
      ELSE
c-----------------------------------------------------------------------
c in the case of the chord method, compute the corrector error,
c and solve the linear system with that as right-hand side and
c p as coefficient matrix.
c-----------------------------------------------------------------------
         y(:n) = h*savf(:n) - (yh(:n,2)+acor(:n))
         CALL slvs (wm, iwm, y, savf)
         IF (iersl < 0) GOTO 430
         IF (iersl > 0) GOTO 410
         del = vnorm(n,y,ewt)
         acor(:n) = acor(:n) + y(:n)
         y(:n) = yh(:n,1) + el(1)*acor(:n)
c-----------------------------------------------------------------------
c test for convergence.  IF m.gt.0, an estimate of the convergence
c rate constant is stored in crate, and this is used in the test.
c-----------------------------------------------------------------------
      END IF
      IF (m /= 0) crate = MAX(0.2_dp*crate,del/delp)
      dcon = del*MIN(one,1.5_dp*crate)/(tesco(2,nq)*conit)
      IF (dcon <= one) GOTO 450
      m = m + 1
      IF (m == maxcor) GOTO 410
      IF (m>=2 .and. del>2.0_dp*delp) GOTO 410
      delp = del
      CALL f (neq, tn, y, savf)
      nfe = nfe + 1
      GOTO 270
c-----------------------------------------------------------------------
c the corrector iteration failed to converge.
c IF miter .ne. 0 and the jacobian is out of date, pjac is called for
c the next try.  otherwise the yh array is retracted to its values
c before prediction, and h is reduced, IF possible.  IF h cannot be
c reduced or mxncf failures have occurred, EXIT with kflag = -2.
c-----------------------------------------------------------------------
  410 CONTINUE
      IF (miter==0 .or. jcur==1) GOTO 430
      icf = 1
      ipup = miter
      GOTO 220
  430 CONTINUE
      icf = 2
      ncf = ncf + 1
      rmax = 2.0_dp
      tn = told
      i1 = nqnyh + 1
      DO jb = 1, nq
         i1 = i1 - nyh
cdir$ ivdep
         yh1(i1:nqnyh) = yh1(i1:nqnyh) - yh1(i1+nyh:nqnyh+nyh)
      END DO
      IF (ierpj<0 .or. iersl<0) GOTO 680
      IF (ABS(h) <= hmin*1.00001_dp) GOTO 670
      IF (ncf == mxncf) GOTO 670
      rh = 0.25_dp
      ipup = miter
      iredo = 1
      GOTO 170
c-----------------------------------------------------------------------
c the corrector has converged.  jcur is set to 0
c to signal that the jacobian involved may need updating later.
c the local error test is made and control passes to statement 500
c IF it fails.
c-----------------------------------------------------------------------
  450 CONTINUE
      jcur = 0
      IF (m == 0) dsm = del/tesco(2,nq)
      IF (m > 0) dsm = vnorm(n,acor,ewt)/tesco(2,nq)
      IF (dsm <= one) THEN
c-----------------------------------------------------------------------
c after a successful step, update the yh array.
c consider changing h IF ialth = 1.  otherwise decrease ialth by 1.
c IF ialth is then 1 and nq .lt. maxord, then acor is saved for
c USE in a possible order increase on the next step.
c IF a change in h is considered, an increase or decrease in order
c by one is considered also.  a change in h is made ONLY IF it is by a
c factor of at least 1.1.  IF not, ialth is set to 3 to prevent
c testing for that many steps.
c-----------------------------------------------------------------------
         kflag = 0
         iredo = 0
         nst = nst + 1
         hu = h
         nqu = nq
         DO j = 1, l
            yh(:n,j) = yh(:n,j) + el(j)*acor(:n)
         END DO
         ialth = ialth - 1
         IF (ialth == 0) GOTO 520
         IF (ialth > 1) GOTO 700
         IF (l == lmax) GOTO 700
         yh(:n,lmax) = acor(:n)
         GOTO 700
c-----------------------------------------------------------------------
c the error test failed.  kflag keeps track of multiple failures.
c restore tn and the yh array to their previous values, and prepare
c to try the step again.  compute the optimum step SIZE for this or
c one lower order.  after 2 or more failures, h is forced to decrease
c by a factor of 0.2 or less.
c-----------------------------------------------------------------------
      END IF
      kflag = kflag - 1
      tn = told
      i1 = nqnyh + 1
      DO jb = 1, nq
         i1 = i1 - nyh
cdir$ ivdep
         yh1(i1:nqnyh) = yh1(i1:nqnyh) - yh1(i1+nyh:nqnyh+nyh)
      END DO
      rmax = 2.0_dp
      IF (ABS(h) <= hmin*1.00001_dp) GOTO 660
      IF (kflag <= (-3)) GOTO 640
      iredo = 2
      rhup = 0
      GOTO 540
c-----------------------------------------------------------------------
c regardless of the success or failure of the step, factors
c rhdn, rhsm, and rhup are computed, by which h could be multiplied
c at order nq - 1, order nq, or order nq + 1, respectively.
c in the CASE of failure, rhup = 0.0 to avoid an order increase.
c the largest of these is determined and the new order chosen
c accordingly.  IF the order is to be increased, we compute one
c additional scaled derivative.
c-----------------------------------------------------------------------
  520 CONTINUE
      rhup = 0
      IF (l /= lmax) THEN
         savf(:n) = acor(:n) - yh(:n,lmax)
         dup = vnorm(n,savf,ewt)/tesco(3,nq)
         exup = one/(l + 1)
         rhup = one/(1.4_dp*dup**exup + 0.0000014_dp)
      END IF
  540 CONTINUE
      exsm = one/(l)
      rhsm = one/(1.2_dp*dsm**exsm + 0.0000012_dp)
      rhdn = 0
      IF (nq /= 1) THEN
         ddn = vnorm(n,yh(1,l),ewt)/tesco(1,nq)
         exdn = one/(nq)
         rhdn = one/(1.3_dp*ddn**exdn + 0.0000013_dp)
      END IF
      IF (rhsm < rhup) THEN
         IF (rhup > rhdn) GOTO 590
      ELSE
         IF (rhsm < rhdn) GOTO 580
         newq = nq
         rh = rhsm
         GOTO 620
      END IF
  580 CONTINUE
      newq = nq - 1
      rh = rhdn
      IF (kflag<0 .and. rh>one) rh = 1
      GOTO 620
  590 CONTINUE
      newq = l
      rh = rhup
      IF (rh < 1.1_dp) GOTO 610
      r = el(l)/(l)
      yh(:n,newq+1) = acor(:n)*r
      GOTO 630
  610 CONTINUE
      ialth = 3
      GOTO 700
  620 CONTINUE
      IF (kflag==0 .and. rh<1.1_dp) GOTO 610
      IF (kflag <= (-2)) rh = MIN(rh,0.2_dp)
c-----------------------------------------------------------------------
c IF there is a change of order, reset nq, l, and the coefficients.
c in ANY CASE h is reset according to rh and the yh array is rescaled.
c THEN EXIT from 690 IF the step was ok, or reDO the step otherwise.
c-----------------------------------------------------------------------
      IF (newq == nq) GOTO 170
  630 CONTINUE
      nq = newq
      l = nq + 1
      iret = 2
      GOTO 150
c-----------------------------------------------------------------------
c control reaches this section IF 3 or more failures have occured.
c IF 10 failures have occurred, EXIT with kflag = -1.
c it is assumed that the derivatives that have accumulated in the
c yh array have errors of the wrong order.  hence the first
c derivative is recomputed, and the order is set to 1.  THEN
c h is reduced by a factor of 10, and the step is retried,
c until it succeeds or h reaches hmin.
c-----------------------------------------------------------------------
  640 CONTINUE
      IF (kflag == (-10)) GOTO 660
      rh = 0.1_dp
      rh = MAX(hmin/ABS(h),rh)
      h = h*rh
      y(:n) = yh(:n,1)
      CALL f (neq, tn, y, savf)
      nfe = nfe + 1
      yh(:n,2) = h*savf(:n)
      ipup = miter
      ialth = 5
      IF (nq == 1) GOTO 200
      nq = 1
      l = 2
      iret = 3
      GOTO 150
c-----------------------------------------------------------------------
c All returns are made through this section.  h is saved in hold
c to allow the caller to change h on the next step.
c-----------------------------------------------------------------------
  660 CONTINUE
      kflag = -1
      GOTO 720
  670 CONTINUE
      kflag = -2
      GOTO 720
  680 CONTINUE
      kflag = -3
      GOTO 720
  690 CONTINUE
      rmax = 10
  700 CONTINUE
      r = one/tesco(2,nqu)
      acor(:n) = acor(:n)*r
  720 CONTINUE
      hold = h
      jstart = 1

      END SUBROUTINE stode
