!>    \brief Module containing Fletcher-Reeves (non-linear Conjugate Gradient)
!!    routines including linear search algorithm
      MODULE f1dim_mod 
      USE stel_kinds, ONLY: dp
      INTEGER :: ncom
      REAL(dp), DIMENSION(:), POINTER :: pcom,xicom
!      INTERFACE
!         FUNCTION func(p)
!         USE stel_kinds, ONLY: dp
!         IMPLICIT NONE
!         REAL(dp), DIMENSION(:), INTENT(IN) :: p
!         REAL(dp) :: func
!         END FUNCTION func
!      END INTERFACE
      
      CONTAINS
      
!>    \brief One-dimension function used by linmin, passed to mnbrak/brent to find minimum for line search.
      FUNCTION f1dim(x)
      IMPLICIT NONE
      REAL(dp), INTENT(IN) :: x
      REAL(dp) :: f1dim
!Used by linmin as the one-dimensional function passed to mnbrak and brent.
      REAL(dp), DIMENSION(:), ALLOCATABLE :: xt
      ALLOCATE(xt(ncom))
      xt(:)=pcom(:)+x*xicom(:)
      f1dim=wfunc(xt)
      DEALLOCATE(xt)
      END FUNCTION f1dim

!>    \brief Fletcher-Reeves-Polak-Ribiere minimization
      SUBROUTINE Fletcher_Reeves (p,ftol,iter,fret,n)
      USE stel_kinds, ONLY: dp
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: n
      INTEGER, INTENT(OUT) :: iter
      REAL(dp), INTENT(IN) :: ftol
      REAL(dp), INTENT(OUT) :: fret
      REAL(dp), DIMENSION(n), INTENT(INOUT) :: p
      
      INTEGER, PARAMETER :: itmax=50
      REAL(dp), PARAMETER :: eps=1.0e-10_dp

!Given a starting point p that is a vector of length N, Fletcher-Reeves-Polak-Ribiere minimization
!is performed on a function func, using its gradient as calculated by a routine
!dfunc. The convergence tolerance on the function value is input as ftol. Returned quantities
!are p (the location of the minimum), iter (the number of iterations that were
!performed), and fret (the minimum value of the function). The routine linmin is called
!to perform line minimizations.
!Parameters: ITMAX is the maximum allowed number of iterations; EPS is a small number
!to rectify the special case of converging to exactly zero function value.
      INTEGER :: its
      REAL(dp) :: dgg,fp,gam,gg
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: g, h, xi
!      REAL(dp), DIMENSION(SIZE(p)) :: g,h,xi
      
      ALLOCATE (g(SIZE(p)), h(SIZE(p)), xi(SIZE(p)))
      fp=wfunc(p)                             !Initializations.
      xi=dfunc(p)
      g=-xi
      h=g
      xi=h
      DO its=1,itmax                          !Loop over iterations.
         iter=its
         CALL linmin(p,xi,fret)               !Next statement is the normal return:
         IF (2*ABS(fret-fp) <= ftol*(ABS(fret)+ABS(fp)+eps)) RETURN
         fp=fret
         xi=dfunc(p)
         gg=DOT_PRODUCT(g,g)
!         dgg=DOT_PRODUCT(xi,xi)                      !This statement for Fletcher-Reeves.
         dgg=DOT_PRODUCT(xi+g,xi)             !This statement for Polak-Ribiere.
         IF (gg == 0._dp) RETURN                !Unlikely. If gradient is exactly zero then we are 
         gam=dgg/gg                           !already done.
         g=-xi
         h=g+gam*h
         xi=h
      END DO
 
      DEALLOCATE (g, h, xi)
      
      END SUBROUTINE Fletcher_Reeves

!>    \brief Function to compute the MHD energy used in line search
      FUNCTION wfunc(p)
      USE stel_kinds, ONLY: dp
      USE shared_functions, ONLY: getwmhd
!      USE diagnostics_mod, ONLY: getbgradp2
      IMPLICIT NONE
      REAL(dp), DIMENSION(:), INTENT(IN) :: p
      REAL(dp) :: wfunc

      wfunc = getwmhd(p, .FALSE.)

      END FUNCTION wfunc

!>    \brief Function to compute gradient of MHD energy
      FUNCTION dfunc(p)
      USE stel_kinds, ONLY: dp
      USE shared_functions, ONLY: funct_island
      USE shared_data, ONLY: xc, gc
!      USE diagnostics_mod, ONLY: funct_fgradp, xc, gc
      IMPLICIT NONE
      REAL(dp), DIMENSION(:), INTENT(IN) :: p
      REAL(dp), DIMENSION(SIZE(p)) :: dfunc

      xc = p
      CALL funct_island
      dfunc = gc

      END FUNCTION dfunc

!>    \brief Routine to find minimum of func along a given direction xi from point p
      SUBROUTINE linmin(p,xi,fret)
      IMPLICIT NONE
      REAL(dp), INTENT(OUT) :: fret
      REAL(dp), DIMENSION(:), TARGET, INTENT(INOUT) :: p,xi
      REAL(dp), PARAMETER :: tol=1.0e-4_dp
!Given an N-dimensional point p and an N-dimensional direction xi, both vectors of length
!N, moves and resets p to where the fixed-name function func takes on a minimum along
!the direction xi from p, and replaces xi by the actual vector displacement that p was
!moved. Also returns as fret the value of func at the returned location p. This is actually
!all accomplished by calling the routines mnbrak and brent.
!Parameter TOL: Tolerance passed to brent.
      REAL(dp) :: ax,bx,fa,fb,fx,xmin,xx
!      REAL(dp), EXTERNAL :: brent

      ncom=SIZE(p)
      pcom=>p                        !Communicate the global variables to f1dim.
      xicom=>xi
      ax=0                           !Initial guess for brackets.
      xx=1
      CALL mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
      fret=brent(ax,xx,bx,f1dim,tol,xmin)
      xi=xmin*xi                     !Construct the vector results to return.
      p=p+xi
      END SUBROUTINE linmin

!>    \brief Routine to bracket the minimum of func in the downhill direction
      SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
      USE stel_kinds, ONLY: dp
      IMPLICIT NONE
      REAL(dp), INTENT(INOUT) :: ax,bx
      REAL(dp), INTENT(OUT) :: cx,fa,fb,fc
      INTERFACE
         FUNCTION func(x)
         USE stel_kinds, ONLY: dp
         IMPLICIT NONE
         REAL(dp), INTENT(IN) :: x
         REAL(dp) :: func
         END FUNCTION func
      END INTERFACE
      REAL(dp), PARAMETER :: GOLD=1.618034_dp,GLIMIT=100.0_dp,TINY=1.0e-20_dp, zero=0
!Given a function func, and given distinct initial points ax and bx, this routine searches
!in the downhill direction (defined by the function as evaluated at the initial points) and
!returns new points ax, bx, cx that bracket a minimum of the function. Also returned are
!the function values at the three points, fa, fb, and fc.
!Parameters: GOLD is the default ratio by which successive intervals are magnified; GLIMIT
!is the maximum magnification allowed for a parabolic-fit step.
      REAL(dp) :: fu,q,r,u,ulim
      fa=func(ax)
 100      fb=func(bx)
      if (fb > fa) then                            !Switch roles of a and b so that we
                                                   !can go downhill in the direction from a to b.
!         bx = bx/4; goto 100  
         call swap(ax,bx)
         call swap(fa,fb)
      end if
      cx=bx+GOLD*(bx-ax)                           !First guess for c.
      fc=func(cx)
      do                                           !Do-while-loop: Keep returning here
         if (fb < fc) RETURN                       !until we bracket.
!Compute u by parabolic extrapolation from a, b, c. TINY is used to prevent any possible division by zero.
         r=(bx-ax)*(fb-fc)
         q=(bx-cx)*(fb-fa)
         u=bx-((bx-cx)*q-(bx-ax)*r)/(2*sign(max(abs(q-r),TINY),q-r))
         ulim=bx+GLIMIT*(cx-bx)
!We won\92t go farther than this. Test various possibilities:
         if ((bx-u)*(u-cx) > zero) then             !Parabolic u is between b and c: try it.
            fu=func(u) 
            if (fu < fc) then                      !Got a minimum between b and c.
               ax=bx
               fa=fb
               bx=u
               fb=fu
               RETURN
            else if (fu > fb) then                 !Got a minimum between a and u.
               cx=u
               fc=fu
               RETURN
            end if
            u=cx+GOLD*(cx-bx)                      !Parabolic fit was no use. Use default magnification.
            fu=func(u)                                                             
         else if ((cx-u)*(u-ulim) > zero) then      !Parabolic fit is between c and its allowed limit.
            fu=func(u) 
            if (fu < fc) then
               bx=cx
               cx=u
               u=cx+GOLD*(cx-bx)
               call shft(fb,fc,fu,func(u))
            end if
         else if ((u-ulim)*(ulim-cx) >= zero) then  !Limit parabolic u to maximum allowed value.
            u=ulim 
            fu=func(u)
         else                                      !Reject parabolic u, use default magnification.
            u=cx+GOLD*(cx-bx) 
            fu=func(u)
         end if
         call shft(ax,bx,cx,u)
         call shft(fa,fb,fc,fu)                       !Eliminate oldest point and continue.
  
      end do

      END SUBROUTINE mnbrak

!>    \brief Isolates the minimum of func 
      FUNCTION brent(ax,bx,cx,func,tol,xmin)
      USE stel_kinds, ONLY: dp
      IMPLICIT NONE
      REAL(dp), INTENT(IN) :: ax,bx,cx,tol
      REAL(dp), INTENT(OUT) :: xmin
      REAL(dp) :: brent
      REAL(dp), EXTERNAL :: func
      
      INTEGER, PARAMETER :: ITMAX=100
      REAL(dp), PARAMETER :: CGOLD=0.3819660_dp,ZEPS=1.0e-3_dp*EPSILON(ax)
!Given a function func, and given a bracketing triplet of abscissas ax, bx, cx (such that bx
!is between ax and cx, and func(bx) is less than both func(ax) and func(cx)), this
!routine isolates the minimum to a fractional precision of about tol using Brent\92s method.
!The abscissa of the minimum is returned as xmin, and the minimum function value is
!returned as brent, the returned function value.
!Parameters: ITMAX, Maximum allowed number of iterations;gol den ratio; and a small number that
!protects against trying to achieve fractional accuracy for a minimum that happens to be
!exactly zero.
      INTEGER :: iter
      REAL(dp) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
      
      a=min(ax,cx)                 !a and b must be in ascending order, though
      b=max(ax,cx)                 !the input abscissas need not be.
      v=bx                         !Initializations...
      w=v
      x=v
      e=0                          !This will be the distance moved on the step before last.
      fx=func(x) 
      fv=fx
      fw=fx
      do iter=1,ITMAX              !Main program loop.
         xm=0.5_dp*(a+b)
         tol1=tol*abs(x)+ZEPS
         tol2=2*tol1
         if (abs(x-xm) <= (tol2-0.5_dp*(b-a))) then    !Test for done here.
            xmin=x                                     !Arrive here ready to exit with best values.
            brent=fx
            RETURN
         end if
         if (abs(e) > tol1) then                       !Construct a trial parabolic fit.
            r=(x-w)*(fx-fv)
            q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r
            q=2*(q-r)
            if (q > 0.0_dp) p=-p
            q=abs(q)
            etemp=e
            e=d
            if (abs(p) >= abs(0.5_dp*q*etemp) .or. &
            p <= q*(a-x) .or. p >= q*(b-x)) then
!The above conditions determine the acceptability of the parabolic fit. Here it is
!not o.k., so we take the golden section step into the larger of the two segments.
               e=merge(a-x,b-x, x >= xm)
               d=CGOLD*e
            else                                      !Take the parabolic step.
               d=p/q
               u=x+d
               if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
            end if
         else                                         !Take the golden section step into the larger
            e=merge(a-x,b-x, x >= xm )                !of the two segments.
            d=CGOLD*e
         end if
         u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
!Arrive here with d computed either from parabolic fit, or else from golden section.
         fu=func(u)
!This is the one function evaluation per iteration.
         if (fu <= fx) then                           !Now we have to decide what to do with our
            if (u >= x) then                          !function evaluation. Housekeeping follows:
               a=x
            else
               b=x
            end if
            call shft(v,w,x,u)
            call shft(fv,fw,fx,fu)
         else
            if (u < x) then
               a=u
            else
               b=u
            end if
            if (fu <= fw .or. w == x) then
               v=w
               fv=fw
               w=u
               fw=fu
            else if (fu <= fv .or. v == x .or. v == w) then
               v=u
               fv=fu
            end if
         end if
      end do                                            !Done with housekeeping. Back for another iteration

      END FUNCTION brent

!>    \brief Swaps a,b => b,a
      SUBROUTINE swap(a,b)
      USE stel_kinds, ONLY: dp
      REAL(dp), INTENT(INOUT) :: a,b
      REAL(dp) :: dum
      
      dum=a
      a=b
      b=dum
      END SUBROUTINE swap

!>    \brief Shifts a=b, b=c, c=d
      SUBROUTINE shft(a,b,c,d)
      USE stel_kinds, ONLY: dp
      REAL(dp), INTENT(OUT) :: a
      REAL(dp), INTENT(INOUT) :: b,c
      REAL(dp), INTENT(IN) :: d
      
      a=b
      b=c
      c=d
      END SUBROUTINE shft
      
      END MODULE f1dim_mod
