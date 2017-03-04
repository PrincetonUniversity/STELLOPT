      SUBROUTINE initpsi(psic,psiv,xval,chipsi,chipsipsi,dpsi
     $     ,npsi,alpsi,dchi,itype)
c------------------------------------------------------------------------
c     SUBROUTINE which determines psi coordinate variables
c     CALLed by init and by wrtcamino
c     itype=0 gives equally spaced psic
c     itype=1 gives equally spaced psiv
c     note that not ALL quantities are calculated for itype=1
c------------------------------------------------------------------------
      USE precision
      IMPLICIT NONE
      REAL(rprec) :: psic(*),psiv(*),xval(*),chipsi(*),chipsipsi(*)
      REAL(rprec) :: alpsi, dchi, dpsi
      INTEGER :: npsi, itype, i
      REAL :: fpi, xend, pi
      pi=ACOS(-1._dbl)
      IF(itype.eq.0) THEN
         dpsi=1./(npsi-1.)
         DO  i=1,npsi
            psic(i)=(i-1)*dpsi
            IF(alpsi.ge.0.) THEN
               xval(i)=(1.+alpsi)*psic(i)**2/(1.+alpsi*psic(i))
               chipsi(i)=dchi*(1.+alpsi)*(2.+alpsi*psic(i))*
     $              psic(i)/(1.+alpsi*psic(i))**2
               chipsipsi(i)=dchi*2.*(1.+alpsi)
     $              /(1.+alpsi*psic(i))**3
            ELSE
               IF(alpsi.lt.(-1.)) THEN
                  PRINT *,"error in alpsi"
                  STOP
               ENDIF
               fpi=-alpsi*pi
               xend=SIN(fpi*0.5)**2
               xval(i)=(SIN(psic(i)*fpi*0.5))**2/xend
               chipsi(i)=dchi*fpi*SIN(psic(i)*fpi)*0.5/xend
               chipsipsi(i)=dchi*fpi**2*COS(psic(i)*fpi)*0.5/xend
            END IF
            psiv(i)=xval(i)*dchi
         END DO
         xval(npsi)=1.0
         psic(npsi)=1.
         psiv(npsi)=dchi
      ELSE
         dpsi=1./(npsi-1.)
         DO  i=2,npsi
            xval(i)=(i-1)*dpsi
            IF(alpsi.ge.0.) THEN
               psic(i)=(alpsi*xval(i)+SQRT((alpsi*xval(i))**2+
     $              4.*(1+alpsi)*xval(i)))/(2.*(1.+alpsi))
            ELSE
               IF(alpsi.lt.(-1.)) THEN
                  PRINT *,"error in alpsi"
                  STOP
               ENDIF
               fpi=-alpsi*pi
               xend=SIN(fpi*0.5)**2
               psic(i)=ASIN(SQRT(xval(i)*xend))*2./fpi
            END IF
            psiv(i)=xval(i)*dchi
         END DO
         xval(1)=0.
         psic(1)=0.
         psiv(1)=0.
         xval(npsi)=1.0
         psic(npsi)=1.
         psiv(npsi)=dchi
      END IF
      RETURN
      END  SUBROUTINE initpsi
