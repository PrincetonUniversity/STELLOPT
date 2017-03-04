      SUBROUTINE equalarc(psivin,ipsi)
      USE mapg_mod
c-----------------------------------------------------------------------
c     stripped DOwn version of chali from mbc
c     rlm 6/22/94
c     this routine finds psi contour psivin=psiv(ipsi) and RETURN through
c     common blocks the equal arc LENgth values of xs(ipsi,j) and zs(ipsi,j)
c     arcsur(ipsi,j), the arc LENgth along the contour, is also RETURNed
c-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec) :: psivin
      INTEGER :: ipsi, icount=0

      INTEGER is,js,npc,ntw,npco,npcm1
      REAL(rprec) :: delxz
      REAL(rprec) :: xp(nxzd),zp(nxzd),bpsqp(nxzd),tp(nxzd)
c      REAL(rprec) :: csx(3,nxzd),csz(3,nxzd),csbp(3,nxzd)
      REAL(rprec) :: csx(nxzd),csz(nxzd),csbp(nxzd)
      REAL(rprec) :: arc(nxzd)
      REAL(rprec) :: dl,ds,bpsqval
      INTEGER j,nlow
c      REAL(rprec) :: sterpl,f2s,bb(4)
c      INTEGER ier
c-----------------------------------------------------------------------
c     the following common is shared with routine cntour
c-----------------------------------------------------------------------
      REAL(rprec) :: xemin,xemax          !output form cntour
      REAL(rprec) :: yemin,yemax          !output form cntour
      REAL(rprec) :: yxmin,yxmax          !output form cntour
      REAL(rprec) :: xymin,xymax          !output form cntour
      REAL(rprec) :: xaxd,yaxd            !input to cntour
      REAL(rprec) :: dang,arcl,bperr      !input to cntour
      REAL(rprec) :: xmin,xmax            !input to cntour
      REAL(rprec) :: ymin,ymax            !input to cntour
      common/cntd/xaxd,yaxd,
     $     xemin,xemax,yemin,yemax,yxmin,yxmax,xymin,xymax,
     $     dang,arcl,bperr,xmin,xmax,ymin,ymax
c-----------------------------------------------------------------------
c     attempt to find contour using furpl
c     is and js are starting values for search for flux surfaces
c-----------------------------------------------------------------------
      npc=1
      ntw=nxzd
      is=nx/2
      js=nz/2
      IF(ipsi .eq. 2222) THEN !	set to 2 for PRINTing
        PRINT*," shape(psixz) ",shape(psixz)
        PRINT*," shape(bpsq) ",shape(bpsq)
        PRINT*," shape(xgrid) ",shape(xgrid)
        PRINT*," shape(zgrid) ",shape(zgrid)
        PRINT*," shape(psivin) ",shape(psivin)
        PRINT*," shape(xp) ",shape(xp)
        PRINT*," shape(zp) ",shape(zp)
        PRINT*," shape(bpsqp) ",shape(bpsqp)
        PRINT*,"npc,nx,nz,ntw,nxd,nzd,is,js=",
     &		npc,nx,nz,ntw,nxd,nzd,is,js
      ENDIF
      CALL furplm(psixz,bpsq,xgrid,zgrid,psivin,xp,zp,bpsqp,npc,
     $     nx,nz,ntw,nxd,nzd,is,js)
      IF(ipsi .eq. 2222) THEN ! set to 2 for PRINTing
       PRINT*,npc,psivin, MAXVAL(psixz),MINVAL(psixz)
       PRINT*,(psivin > MAXVAL(psixz)).and.(psivin > MINVAL(psixz))
      ENDIF

c      npc=0
c      npfit=400
      npfit=40
      IF(npc.le.npfit) THEN
c         WRITE(6,'("used contour",2i5)') ipsi,npc
         xaxd=xaxis
         yaxd=zaxis
         npco=npc
         xmin=xgrid(1)
         xmax=xgrid(nx)
         ymin=zgrid(1)
         ymax=zgrid(nz)
         bperr=0.1
         IF(ipsi.eq.2) THEN
            arcl=max(dx,dz)*0.1
         ELSE
            arcl=arcsur(ipsi-1,nthet)/npfit
         ENDIF
         npfit=400
         dang=1./npfit
         CALL cntourp(xgrid,nx,zgrid,nz,csplpsi,xp,zp,bpsqp,
     $        npc,dx,dz,nxzd,psivin,nxd)
         CALL sortr(xp,zp,bpsqp,npc,xaxis,zaxis,xaxis,zaxis,0,1)
      ELSE
         IF(xp(1).ne.xp(npc).or.zp(1).ne.zp(npc)) THEN
            WRITE(6,'("error in equalarc finding contour ",
     $           i3," of ",i3)') ipsi,npsi
            WRITE(6,'("IF this is final contour,"
     $           ," try reducing percenflux")')
            STOP
         ENDIF
         CALL sortr(xp,zp,bpsqp,npc,xaxis,zaxis,xaxis,zaxis,0,1)
         delxz=0.20*min(xgrid(2)-xgrid(1),zgrid(2)-zgrid(1))
         CALL packk(bpsqp,xp,zp,npc,delxz)
      ENDIF
      IF(npc.le.5) THEN
         WRITE(6,'("error in equalarc--npc too small=")') npc
         STOP
      ENDIF
c-----------------------------------------------------------------------
c     set up arclengths on surface
c-----------------------------------------------------------------------
!      icount=icount+1; PRINT*,'HERE',icount
      tp(1)=0.
      npcm1=npc-1
      DO  j=2,npc
         tp(j)=tp(j-1)+SQRT((xp(j)-xp(j-1))**2+(zp(j)-zp(j-1))**2)
      END DO
c-----------------------------------------------------------------------
c     spline xp wrt t 
c-----------------------------------------------------------------------
      CALL spline(tp,xp,npc,-1.e31_dbl,-1.e31_dbl,csx)
c-----------------------------------------------------------------------
c     spline zp wrt t 
c-----------------------------------------------------------------------
      CALL spline(tp,zp,npc,-1.e31_dbl,-1.e31_dbl,csz)
c-----------------------------------------------------------------------
c     calculate arclength to each point
c-----------------------------------------------------------------------
      CALL garc(tp,xp,zp,csx,csz,nxzd,arc,npc)
c-----------------------------------------------------------------------
c     spline xp wrt arc
c-----------------------------------------------------------------------
      CALL spline(arc,xp,npc,-1.e31_dbl,-1.e31_dbl,csx)
c-----------------------------------------------------------------------
c     spline zp wrt arc
c-----------------------------------------------------------------------
      CALL spline(arc,zp,npc,-1.e31_dbl,-1.e31_dbl,csz)
c-----------------------------------------------------------------------
c     spline bpsqp wrt arc
c-----------------------------------------------------------------------
      CALL spline(arc,bpsqp,npc,-1.e31_dbl,-1.e31_dbl,csbp)
c-----------------------------------------------------------------------
c     convert to equal arc LENgth along contour
c-----------------------------------------------------------------------
      ds=arc(npc)/(nthet-1.)
      DO  j=1,nthet-1
         dl=(j-1)*ds
         nlow=0
         CALL splINT(arc,xp,csx,npc,dl,xs(ipsi,j),nlow)
         CALL splINT(arc,zp,csz,npc,dl,zs(ipsi,j),nlow)
         CALL splINT(arc,bpsqp,csbp,npc,dl,bpsqval,nlow)
         IF(bpsqval.gt.0.) THEN
            bps(ipsi,j)=SQRT(bpsqval)
         ELSE
            bps(ipsi,j)=0.
         ENDIF
         arcsur(ipsi,j)=dl
      END DO
      xs(ipsi,nthet)=xs(ipsi,1)
      zs(ipsi,nthet)=zs(ipsi,1)
      bps(ipsi,nthet)=bps(ipsi,1)
      arcsur(ipsi,nthet)=arc(npc)
      END SUBROUTINE equalarc
