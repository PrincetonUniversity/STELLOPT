      SUBROUTINE garc(tp,xp,zp,csx,csz,npmax,arc,npc)
      USE precision
      IMPLICIT NONE
      INTEGER :: npmax,npc
      REAL(rprec) :: tp(*),xp(*),zp(*),arc(*)
      REAL(rprec) :: csx(*),csz(*)
c      REAL(rprec) :: csx(3,npmax),csz(3,npmax)
      REAL(rprec) :: yg(4),wg(4)
      REAL(rprec) :: el,summ,a,b
      INTEGER ngaus,j,is,i
      REAL(rprec) :: ws1,ws2,ws,delta
c ----------------------------------------------------------
c gaussian quadrature constants for the INTerval from 0 to 1
c ----------------------------------------------------------
      data yg/.069431844202974,.330009478207572,.669990521792428,
     1       .930568155797027/,wg/.173927422568727,
     2       2*.326072577431273,.173927422568727/,ngaus/4/
      arc(1)=0.
      el=0.
      DO  j=2,npc
         is=j-1
c --------------------------------------------------------------------
c this routine USEs a four poINT gaussian quadrature to compute the
c INTegral of SQRT((dx/dt)**2+(dz/dt)**2) with respect to t around the
c curve. the quadrature is on poINT number tp.
c --------------------------------------------------------------------
         summ=0.
         delta=tp(j)-tp(is) ! delta is h
         DO  i=1,ngaus
            b=yg(i)
            a=1.-yg(i)
            ws1=(xp(j)-xp(is))/delta-(3.*a**2-1.)/6.*delta*csx(is)
     $           +(3.*b**2-1.)/6.*delta*csx(j)
            ws2=(zp(j)-zp(is))/delta-(3.*a**2-1.)/6.*delta*csz(is)
     $           +(3.*b**2-1.)/6.*delta*csz(j)
c            yq=yg(i)*delta      !yq is a*h, yg(i) is b, 
c            ws1=(3.*csx(3,is)*yq+2.*csx(2,is))*yq+csx(1,is)
c            ws2=(3.*csz(3,is)*yq+2.*csz(2,is))*yq+csz(1,is)
            ws=SQRT(ws1*ws1+ws2*ws2)
            summ=summ+ws*wg(i)
         END DO
         el=el+summ*delta
         arc(j)=el
      END DO
      RETURN
      END
