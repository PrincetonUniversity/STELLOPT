      SUBROUTINE mapperb
      USE precision
      USE mapg_mod
c-----------------------------------------------------------------------
c     changed dpsi to dpsic  : rlm 7/3/96
c     changed dpsi2 to dpsiv : same
c-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec) :: dpsic,xval(kpsi),chipsipsi(kpsi)
      INTEGER :: itype, i,j
      REAL(rprec) :: psival2(nxd),dum2(nxd)
      REAL(rprec) :: dpsiv
      REAL(rprec) :: dchi,psitest
c     have trouble finding psilim contour on mANY eqdsks
c     percenflux ALLows edge scrapeoff--rlm 10/2/96
      dchi=(psilim-psiaxis)*percenflux
      itype=0                   ! makes equally spaced psic
      CALL initpsi(psic,psiv,xval,chipsi,chipsipsi,dpsic,npsi,
     $     alpsi,dchi,itype)
c      WRITE(6,'("dchi=",e12.4)') dchi
      DO j=2,npsi
         psiv(j)=psiv(j)+psiaxis
      ENDdo
      psiv(1)=psiaxis
      psiv(npsi)=psiaxis+dchi
      dpsiv=(psilim-psiaxis)/(nx-1)
      DO  i=1,nx
         psival2(i)=psiaxis+(i-1.)*dpsiv
      END DO
      psival2(nx)=psilim
c
c set up pprime, ffprime, f, and qsfin arrays
      CALL spline1d(pprime,psiv,npsi,spp,psival2,nx,dum2)
      CALL spline1d(ffprime,psiv,npsi,sffp,psival2,nx,dum2)
      CALL spline1d(fval,psiv,npsi,sf,psival2,nx,dum2)
      CALL spline1d(qsfin,psiv,npsi,qpsi,psival2,nx,dum2)
      CALL spline1d(press,psiv,npsi,sp,psival2,nx,dum2)
      CALL spline1d(pw,psiv,npsi,pressw,psival2,nx,dum2)
      CALL spline1d(pwp,psiv,npsi,pwprim,psival2,nx,dum2)
      CALL spline1d(rho,psiv,npsi,rho0,psival2,nx,dum2)
      CALL spline1d(rhop,psiv,npsi,rho0p,psival2,nx,dum2)
c
c     bicubic spline on psi used by cntour which is called by equalarc
      CALL bcspline(xgrid,zgrid,psixz,nx,nz,nxd,csplpsi)
c     first surface is easy
      i=1
      DO j=1,nthet
         xs(i,j)=xaxis
         zs(i,j)=zaxis
         bps(i,j)=0.
      END DO
      DO i=2,npsi
         psitest=psiv(i)
!         WRITE(*,'("i=",i5," psitest=",e12.4)') i,psitest
         psitest=psiv(i)
         CALL equalarc(psitest,i)
      END DO
      END SUBROUTINE mapperb
