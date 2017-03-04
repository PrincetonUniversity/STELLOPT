      SUBROUTINE cntourp(x,nw,y,nh,cspln,xc,yc,gsq,ipts,dx,dy,ntmx,psivl
     $                    ,nxd)
c...........................................................................
c     02/04/97 (yr)
c     add a new output array gsq to cntour and rename it to cntourp
c     the original cntour is still avaible in the lib
c     06/14/96
c     cntour generates a contour of ordered points (xc(i),yc(i)) for
c     i = 1,ipts along a contour psivl of a FUNCTION psi(x,y). The function
c     values psi(x,y) are defined from the spline coefficients in cspln.
c     The contour psivl must fully encircle a given poINT (xaxd,yaxd).
c
c     The contour is determined by defining a series of rays from
c     (xaxd,yaxd) INTersecting with the contour.
c     Each ray is defined by an equation y = a*x + b  or  x = a*y + b,
c     depending on whether the angle is in the left or right quadrants
c     or in the top or bottom quadrants respectively.
c     For a given ray, the INTersection is first bounded by a coarse sliding
c     INTerval search along the ray to find two points on the ray bounding
c     the INTersection.  The search is limited to a rectangle specified by
c     xmin,xmax,ymin,ymax.  A sliding INTerval search, rather than a
c     binary search, is USEd for the coarse grid so that non-monotonic
c     psi can be handled.  Newtons Method is THEN USEd to converge to the
c     INTersection poINT and the coordinates are stored as (xc(i),yc(i)).
c     This task is performed by the SUBROUTINE 'rayatt'.
c
c     The ray angular spacing is controlled dynamically and so the total
c     number of contour points (xc(i),yc(i)) is determined by the routine.
c     Normally, the increment is defined as the minimum of arcl/rad and
c     2pi*dang*(rad0/rad), WHERE dang is an angle fraction of 2pi and
c     arcl is an arclength spacing between successive points, both input
c     through common/cntd/, and rad and rad0 are the computed radii at the
c     respective and first (i.e. theta = 0) contour points respectively.
c     After each poINT on the contour is found, the relative change in the
c     poloidal field from the previous poINT is checked to ensure that it
c     is less than bpers.  If the change in poloidal field exceeds bpers,
c     the angular increment is successively reduced, up to ihlfmx times,
c     and the poINT is recalculated.
c     If the number of contour points exceeds the maximum ALLowed DIMENSION
c     ntmx, the mapping for this contour is restarted with bpers increased
c     and appropriate warnings are PRINTed.  The mapping can restart
c     up to jcmax tries, after which, the routine aborts.
c
c     Input parameters :
c     x,y :        Vectors containing the rectangular grid poINT values on
c                  which the contoured FUNCTION is defined.
c     nw,nh :      Respective DIMENSIONs of the vectors x and y.
c     cspln :      Spline coefficients for the contoured FUNCTION.
c     dx,dy :      Increment (in metres) for the coarse grid search of
c                  the INTersection between the contour and rays.
c                  dx and dy should normally be chosen so that, over the
c                  distance dx or dy, psi is SINgle valued.  This can
c                  usually be ensured by setting dx and dy to the x and y
c                  grid increments.
c     ntmx :       Maximum DIMENSION for the vector containing the contour
c                  points.
c     psivl :      Value of the contour to be mapped.
c
c     Input paramaters passed through common/cntd/ :
c     xaxd,yaxd :  Interior poINT for the contour.  The rays INTersecting
c                  the contour emmanate from the poINT (xaxd,yaxd).
c     xmin,xmax :  The search for INTersection points is limited
c                  horizontally to (xmin,xmax).
c     ymin,ymax :  The search for INTersection points is limited vertically
c                  to (ymin,ymax).
c     dang :       The maximum angular increment ALLowed for the generated
c                  contour points is limited to
c                       2pi*dang*(rad(theta=0)/rad(theta)).
c                  In general, dang should be set according to shape of the
c                  equilibrium.  If dang is too small, mANY unnecessary
c                  points will be generated.  If dang is too large for
c                  bperr (defined below),  the routine wastes computation
c                  time in cutting the step angle DOwn to SIZE.
c                  dang = 0.03 near psilim and dang = 0.10 near psimax
c                  is usually satisfactory.
c     arcl :       The maximum arclength ALLowed between generated contour
c                  points is limited to arcl.  For highly elongated plasmas
c                  control with dang is difficult to set for ALL contours.
c                  For such CASEs, arcl should be USEd instead to set the
c                  angular increment.  arcl can be set to a large number to
c                  overide this option so that dang is always USEd.
c     bpers :      Relative change in poloidal b field between (xc(i),yc(i))
c                  and (xc(i+1),yc(i+1)).  If dang yields a relative change
c                  in poloidal b less than bperr, THEN the the computed
c                  poINT is retained.  Otherwise, dang is successively
c                  reduced by stfrac until this condition is met.  bpers
c                  is started at bperr (set in chali).  After each failure
c                  (i.e. IF ipts gets too large), THEN the calculation is
c                  restarted with larger bpers, up to a maximum of jcmax
c                  tries.
c
c     Output parameters :
c     xc,yc :      Contour points generated at INTersections of the
c                  psivl contour and the rays emmanating from (xaxd,yaxd).
c     gsq   :      square of grad(psi)--changed to bpsquared rlm 2/11/97
c     ipts :       Number of contour points generated.  ipts must be less
c                  than or equal to ntmx.
c
c     Output parameters passed through common/cntd/ :
c     xemin,xemax :  Minimum and maximum x values of the contour.
c     yemin,yemax :  Minimum and maximum y values of the contour.
c     xymin,xymax :  x values at y=yemin and y=yemax respectively.
c     yxmin,yxmax :  y values at x=xemin and x=xemax respectively.
c
c     jcmax is the maximum number of restarts ALLowed.
c     ihlfmx is the maximum number of increment halvings ALLowed for
c     each contour poINT calculation.
c     rndoff0 is a roundoff error for checking for division by zero.
c     stfrac is the fractional reduction in step SIZE imposed when the
c     test for the poloidal field fails.
c
      USE precision
      IMPLICIT NONE
      INTEGER nw,nh,nxd,ntmx,ipts
      REAL(rprec) x(nw),y(nh),cspln(2,nxd,nh),xc(ntmx)
     $	,yc(ntmx),gsq(ntmx), dx,dy
      REAL(rprec) psivl
      REAL(rprec) pds(6)
      REAL(rprec) xaxd,yaxd,psiaxd
      REAL(rprec) xemin,xemax,yemin,yemax,yxmin,yxmax,xymin,xymax
      REAL(rprec) dang,arcl,bperr,xmin,xmax,ymin,ymax
      common/cntd/xaxd,yaxd,
     $     xemin,xemax,yemin,yemax,yxmin,yxmax,xymin,xymax,
     $     dang,arcl,bperr,xmin,xmax,ymin,ymax
      INTEGER jcmax,ihlfmx
      REAL(rprec) stfrac,delta
      PARAMETER (jcmax=5, ihlfmx=4)
      PARAMETER (stfrac=0.5,delta=0.1)
c........................................................................
c     local variables:
      REAL(rprec) pi,twopi
      REAL(rprec) thet,dthet,theta0,dthet0,dthet1,dthet2,xn,yn,rad0,rad
      REAL(rprec) bp1,bp2,bpers,bpdiff
      INTEGER IFail
      INTEGER ntmx1,ier,jcount,ihalf,ipt1,ik
      REAL(rprec) xcfn,ycfn,deli2,deln2,delfn2,delti2,deltn2
c.........................................................................
c     Initialization
c.........................................................................
      pi      = ACOS(-1._dbl)
      twopi   = 2.00*pi
      xemin   = +1.0e+10
      xemax   = -1.0e+10
      yemin   = +1.0e+10
      yemax   = -1.0e+10
      xymin   =  0.0
      xymax   =  0.0
      yxmin   =  0.0
      yxmax   =  0.0
      ntmx1   = ntmx - 1
      dthet0  = dang*twopi
      CALL bcsplINT(x,y,cspln,nw,nh,nxd,xaxd,yaxd,pds,1,ier)
      psiaxd = pds(1)
      IF(psiaxd .gt. psivl) THEN
         WRITE(6,*)'errorin cntourp-cnt1'
         STOP
      ENDIF
c.........................................................................
c     find (x,y) at thet=0 and initialize some quantities
c.........................................................................
      thet=0.
      CALL rayatt(x,nw,y,nh,cspln,nxd,dx,dy,xaxd,yaxd,psiaxd,
     $     xmin,xmax,ymin,ymax,psivl,thet,xn,yn,pds,ifail)
      IF (ifail .gt. 0)stop
      bp1   = SQRT(pds(2)*pds(2) + pds(3)*pds(3))/xn
      xemin  = min(xn,xemin)
      xemax  = max(xn,xemax)
      yemin  = min(yn,yemin)
      yemax  = max(yn,yemax)
      IF(xn .eq. xemax) yxmax  = yn
      IF(xn .eq. xemin) yxmin  = yn
      IF(yn .eq. yemax) xymax  = xn
      IF(yn .eq. yemin) xymin  = xn
      xc(1)=xn
      yc(1)=yn
      gsq(1)=(pds(2)**2+pds(3)**2)/xn**2
      rad0=SQRT((xn-xaxd)**2+(yn-yaxd)**2)
c...........................................................................
c     main loop
c...........................................................................
c     Set up the loop to increment bpers up to jcmax tries, in CASE that 
c     bpers is too samll and there are too mANY contour points found.
      jcount = 0
 5    jcount = jcount + 1    ! restart poINT    
      bpers  = bperr*jcount
      ipts   = 1
      ihalf  = 0
      dthet  = min(dthet0,arcl/rad0)
      thet   = dthet
      theta0 = thet
  10  CONTINUE               ! starting poINT of main loop
         CALL rayatt(x,nw,y,nh,cspln,nxd,dx,dy,xaxd,yaxd,psiaxd,
     $     xmin,xmax,ymin,ymax,psivl,thet,xn,yn,pds,ifail)
         IF (ifail .gt. 0)stop 
c
c        Check for sufficient accuracy in the poINT spacing for theta.
c        The accuracy test is based on a relative error in the poloidal
c        magnetic field of bpers.  If the spacing is too large, set theta
c        back to its previous value theta0, decrease dtheta, and recalculate
c        the point.  The spacing can be repeatedly reduced ihlfmx times.
c        Excessive accumulation near the x points is avoided by keeping
c        ihlfmx small; IF ihlfmx is exceeded, the mapping simply CONTINUEs
c        to the next point.
c
         bp2   = SQRT(pds(2)*pds(2) + pds(3)*pds(3))/xn
         bpdiff = ABS(bp2-bp1)/max(bp2,bp1)
         IF(bpdiff .ge. bpers) THEN
            ihalf  = ihalf + 1
            IF(ihalf .le. ihlfmx) THEN
               thet  = theta0
               dthet = stfrac*dthet
               thet  = thet + dthet
               go to 10
            ENDIF
         ENDIF
c
c        Increment ipts and check IF it has exceeded the DIMENSIONs.
c        Fill in the new poINT (xc,zc) and CONTINUE.
c
         ipts   = ipts + 1
         IF(ipts .gt. ntmx1)then
            IF (jcount .ge. jcmax)then 
               WRITE(6,*)'error in cntourp--cnt2' ! quit
               STOP
            ENDIF
            go to 5                                     ! restart 
         ENDIF
         bp1    = bp2
         xemin  = min(xn,xemin)
         xemax  = max(xn,xemax)
         yemin  = min(yn,yemin)
         yemax  = max(yn,yemax)
         IF(xn .eq. xemax) yxmax  = yn
         IF(xn .eq. xemin) yxmin  = yn
         IF(yn .eq. yemax) xymax  = xn
         IF(yn .eq. yemin) xymin  = xn
         xc(ipts) = xn
         yc(ipts) = yn
         gsq(ipts)= (pds(2)**2+pds(3)**2)/xn**2
         ihalf=0
c
c        Determine the increment of thet for the next point
c        The angle increment is limited to dang*2pi*(rad0/rad) and to an
c        arclength of arcl metres.
         rad     = SQRT((xn-xaxd)**2 + (yn-yaxd)**2)
         dthet1  = dthet0*(rad0/rad)
         dthet2  = arcl/rad
         dthet   = min(dthet1,dthet2)
         theta0 = thet
         thet   = thet + dthet
      IF(thet .lt. twopi) go to 10
c     END of main loop
c.........................................................................
c     Close the contour at theta = twopi.
c     Eliminate the last calculated poINT IF it is too CLOSE to the
c     closure point.
c.........................................................................
      IF(ipts .gt. 1) THEN
         ipt1   = ipts - 1
         xcfn   = xc(1)
         ycfn   = yc(1)
         deli2  = (xc( 2  )-xc( 1  ))*(xc( 2  )-xc( 1  ))
     $          + (yc( 2  )-yc( 1  ))*(yc( 2  )-yc( 1  ))
         deln2  = (xc(ipts)-xc(ipt1))*(xc(ipts)-xc(ipt1))
     $          + (yc(ipts)-yc(ipt1))*(yc(ipts)-yc(ipt1))
         delfn2 = (xcfn    -xc(ipts))*(xcfn    -xc(ipts))
     $          + (ycfn    -yc(ipts))*(ycfn    -yc(ipts))
         delti2 = deli2*delta*delta
         deltn2 = deln2*delta*delta
         IF((delfn2 .ge. delti2) .and. (delfn2 .ge. deltn2))
     $                            ipts   = ipts + 1
      ENDIF
      thet   = twopi
      xc(ipts) = xc(1)
      yc(ipts) = yc(1)
      gsq(ipts) = gsq(1)
c..........................................................................
c     Write error messages IF bpers needed to be increased.
c..........................................................................
      IF(jcount .gt. 1) THEN
         WRITE(3,2010) psivl,jcount,jcmax,bperr,bpers
         WRITE(3,2020) xaxd,yaxd,psiaxd,psivl,ipts
         WRITE(3,2021) ipts,(xc(ik),ik=1,ipts)
         WRITE(3,2022) ipts,(yc(ik),ik=1,ipts)
      ENDIF
 2010 FORMAT(1x,'needed to increment bpers for psivl =',e12.5,/,1x
     $     ,'jcount,jcmax,bperr,bpers ='
     $     ,1x,i3,1x,i3,2x,e12.5,2x,e12.5)
 2020 FORMAT(1x,/,'xaxd/yaxd/psiaxd',/,1x,3e16.8,/,
     $     'psivl',1x,e16.8,/,'ipts',i7,/)
 2021 FORMAT(/,1x,'(xc(i),i=1,',i5,')',/,(10(1x,e11.4)))
 2022 FORMAT(/,1x,'(yc(i),i=1,',i5,')',/,(10(1x,e11.4)))
c..........................................................................
      END SUBROUTINE cntourp

