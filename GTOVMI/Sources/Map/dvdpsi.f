      SUBROUTINE dvdpsi
      USE mapg_mod
      IMPLICIT NONE
      REAL(rprec) :: chi(nthet), mu0
      INTEGER :: ipsi, jthet, j, i
      pi=ACOS(-1._dbl)
      vprime=0
      tflx=0
      volume=0
      DO jthet=1,nthet
       chi(jthet)=2.*pi*(jthet-1)/(nthet-1)
      ENDdo
      dthe=2.*pi/(nthet-1.)
      dpsi=1./(npsi-1.)
!      dpsi=dpsi*(psilim-psiaxis)
      arcrad=0.
      DO jthet=2,nthet
        DO ipsi=2,npsi
          arcrad(ipsi,jthet)=arcrad(ipsi-1,jthet)+
     &	SQRT((xs(ipsi,jthet)-xs(ipsi-1,jthet))**2
     &              +(zs(ipsi,jthet)-zs(ipsi-1,jthet))**2)
        ENDdo
      ENDdo
      mu0=4e-7*pi
      DO  j=1,nthet
         jacob(1,j)=0
         DO  i=2,npsi
            jacob(i,j)=arcsur(i,nthet)/(2.*pi*bps(i,j))
         ENDdo
      ENDdo

      DO ipsi=1,npsi
        vprime(ipsi)=2*pi*trap(nthet,chi(1:nthet),jacob(ipsi,1:nthet))
      ENDdo
      DO ipsi=2,npsi
       volume(ipsi)=trap(ipsi,psiv(1:ipsi),vprime(1:ipsi))
      ENDdo
      vprime(1)=axisv(psiv,vprime,SIZE(psiv))
!      PRINT*,"V = ",REAL(volume(SIZE(volume))),' m^3'
      DO ipsi=2,npsi	! 2pi*<x>*<lp>*<lr>
         tflx(ipsi)=trap(ipsi,psiv(1:ipsi),qsfin(1:ipsi))
      ENDdo
      END SUBROUTINE dvdpsi
