      module penta_functions_mod
	use penta_kind_mod
	implicit none
	contains


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Miscellaneous functions
c   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c----------------------------------------------------------------------
c   This function is used in conjunction with 'zeroin' to find
c    the solution to gamma_e=sum(Z*gamma_i).  Returns difference
c    in fluxes at Er_test, using spline parameters from Er_fit_pass.
c   JL 7/2009
c----------------------------------------------------------------------
c
      real(rknd) function rad_flux(Er_test)
      use penta_kind_mod
      use Er_fit_pass
      implicit none
      real(rknd) :: Er_test, Er_tmp, err
      call pol_int(Er_fit,diff_qg,size(Er_fit),Er_test,Er_tmp,err)
      rad_flux=Er_tmp
      return
      end function rad_flux
c
c-----------------------------------------------------------------------
c  Function zeroin
c  a zero of the function  f(x)  is computed in the interval ax,bx
c
c  input..
c
c  ax     left endpoint of initial interval
c  bx     right endpoint of initial interval
c  f      function subprogram which evaluates f(x) for any x in
c         the interval  ax,bx
c  tol    desired length of the interval of uncertainty of the
c         final result ( .ge. 0.0)
c
c
c  output..
c
c  zeroin abcissa approximating a zero of  f  in the interval ax,bx
c
c
c      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
c  without  a  check.  zeroin  returns a zero  x  in the given interval
c  ax,bx  to within a tolerance  4*macheps*abs(x) + tol, where macheps
c  is the relative machine precision.
c      this function subprogram is a slightly  modified  translation  of
c  the algol 60 procedure  zero  given in  richard brent, algorithms for
c  minimization without derivatives, prentice - hall, inc. (1973).
c
c   modified by JL for kind precision and iteration limit 7/09
c
c----------------------------------------------------------------------------
c
      real(rknd) function zeroin(ax,bx,f,tol)
      use penta_kind_mod
      implicit none
      real(rknd) ax, bx, f, tol
      integer(iknd),save :: jfirst
      integer(iknd) :: iter
      integer(iknd), parameter :: max_iter = 1000_iknd
      real(rknd) :: a, b, c, d, e, eps, fa, fb, fc, tol1,
     1   xm, p, q, r, s
      external :: f

      !compute eps, the relative machine precision
      iter=0_iknd
      eps = epsilon(1._rknd)
      tol1 = 1._rknd + eps

      ! initialization
      a = ax
      b = bx
      fa = f(a)
      fb = f(b)
 
      ! begin step
   20 c = a
      fc = fa
      d = b - a
      e = d
   30 if (dabs(fc) .ge. dabs(fb)) go to 40
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa

      ! convergence test
   40 tol1 = 2._rknd*eps*dabs(b) + 0.5_rknd*tol
      xm = .5_rknd*(c - b)
      iter=iter+1_iknd
      if (iter .gt. max_iter) then
        write(*,*) 'Maximum iterations exceeded in function zeroin:',
     1    max_iter
        stop
      endif
      if (dabs(xm) .le. tol1) go to 90
      if (fb .eq. 0._rknd) go to 90
 
      !is bisection necessary
      if (dabs(e) .lt. tol1) go to 70
      if (dabs(fa) .le. dabs(fb)) go to 70

      ! is quadratic interpolation possible
      if (a .ne. c) go to 50

      ! linear interpolation
      s = fb/fa
      p = 2._rknd*xm*s
      q = 1._rknd - s
      go to 60

      ! inverse quadratic interpolation
   50 q = fa/fc
      r = fb/fc
      s = fb/fa
      p = s*(2._rknd*xm*q*(q - r) - (b - a)*(r - 1._rknd))
      q = (q - 1._rknd)*(r - 1._rknd)*(s - 1._rknd)
    
      ! adjust signs
   60 if (p .gt. 0._rknd) q = -q
      p = dabs(p)

      ! is interpolation acceptable
      if ((2._rknd*p) .ge. (3._rknd*xm*q - dabs(tol1*q))) go to 70
      if (p .ge. dabs(0.5_rknd*e*q)) go to 70
      e = d
      d = p/q
      go to 80

      ! bisection
   70 d = xm
      e = d

      !complete step
   80 a = b
      fa = fb
      if (dabs(d) .gt. tol1) b = b + d
      if (dabs(d) .le. tol1) b = b + dsign(tol1, xm)
      fb = f(b)
      if ((fb*(fc/dabs(fc))) .gt. 0._rknd) go to 20
      go to 30

      ! done
   90 zeroin = b
      return
      end function zeroin
c
c-----------------------------------------------------------------------
c   This function calculates the error function for input xin
c-----------------------------------------------------------------------
c
      real(rknd) function derf(xin)
      use penta_kind_mod
      implicit none
      integer(iknd) :: k
      real(rknd) :: eps, pi, xin2, r, er, c0
      real(rknd), intent(in) :: xin

      eps=1.e-15_rknd
      pi=4._rknd*datan(1._rknd)
      xin2=xin*xin
      if (dabs(xin).lt. 3.5_rknd) then
        er=1._rknd
        r=1._rknd
        do k=1,50
            r=r*xin2/(k+0.5_rknd)
            er=er+r
            if (dabs(r).le.dabs(er)*eps) exit
        enddo
        c0=2._rknd/dsqrt(pi)*xin*dexp(-xin2)
        derf=c0*er
      else
        er=1._rknd
        r=1._rknd
        do k=1,12
            r=-r*(k-0.5_rknd)/xin2
            er=er+r
        enddo
        c0=dexp(-xin2)/(dabs(xin)*dsqrt(pi))
        derf=1.0D0-c0*er
        if (xin.lt.0.0) derf=-derf
      endif
      end function derf

c
c-----------------------------------------------------------------------
c   This function returns the integrand of the energy integral.  
c       Collision frequencies are calculated using calc_perp_coll_freq.
c     Note the value set below for requested monoenergetic coefficients
c       out of range.
c   JL 7/2009 modified from don's fun123
c-----------------------------------------------------------------------
c
      real(rknd) function intfun(Ka)
      use penta_kind_mod
      use intfun_pass
      use thermal_pass
      use bspline
      use parameter_pass
      implicit none
      !dummy variable
      real(rknd) :: Ka
      !local variables
      real(rknd) :: efield, xa, va, nu_perp_aa, xb, Kb,
     1  nu_perp_a, cmul_K, ff, enrm
      real(rknd), allocatable :: nu_perp_ab(:)
      integer(iknd) :: ispec, svtb

      !normalized velocity and velocity
      xa=dsqrt(Ka)
      va=xa*vta

      !define |Er|/v and coefficient range
      efield=abs_Er/va

      !Define energy dependent collision frequencies
      svtb=size(vtb)
      call calc_perp_coll_freq(va,xa,qa,qa,ma,na,Ka,loglam,
     1  vta,Ka,nu_perp_aa)
      allocate(nu_perp_ab(svtb))

      !loop over field species to get collision frequency
      do ispec=1,svtb
        ! Field species norm. energy and velocity
        xb=xa*(vta/vtb(ispec))
        Kb=xb**2_iknd

        !Define energy dependent collision frequency
        call calc_perp_coll_freq(va,xb,qa,qb(ispec),ma,nb(ispec),Kb,
     1    loglam,vta,Ka,nu_perp_ab(ispec))
      enddo

      !define collisionality
      nu_perp_a=nu_perp_aa + sum(nu_perp_ab)
      cmul_K = nu_perp_a/va
      deallocate(nu_perp_ab)

      !use log interp or linear
      if ( log_interp ) then
        efield = dlog(efield)
        enrm = (efield - emin_pass)
     1    /(emax_pass - emin_pass)
        cmul_K=dlog(cmul_K)
      endif

      !Check for cmul or efield values out of range.  If efield is below the smallest database
      !  efield then use the smallest database efield.
      !Then look up coefficient at efield and cmul values from fits.  If the request
      ! is out of range use the value set below.
      if(efield .lt. emin_pass) efield = emin_pass                   !check for efield too small
      if(cmul_K .ge. cmin_pass .and. cmul_K .le. cmax_pass .and.     !check for out of range
     1       efield .ge. emin_pass .and. efield .le. emax_pass) then

        !interpolate database
        ff = dbs2vl(cmul_K,enrm,kcord,keord,xt_c_pass,xt_e_pass,
     1     num_c_pass,num_e_pass,c_spl_pass)
       
        !handle log of coefficient
        if( log_coeff ) then
          intfun=dexp(ff - Ka)*Ka**2*(Ka - 2.5_rknd)  !the K^2 instead of sqrt(K) accounts for K^1.5 dependance of the coefficients
     1      **(jval-1_iknd)
        else 
          intfun=ff*dexp(-Ka)*Ka**2*(Ka - 2.5_rknd)
     1      **(jval-1_iknd)
        endif  !log coefficient
      else 
        !handle out of range requests
        intfun = 0._rknd   
      endif !range check
      end function intfun



      end module penta_functions_mod