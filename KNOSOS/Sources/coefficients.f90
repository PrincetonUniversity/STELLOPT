!Calculate quantities related to bounce-averages: find wells, bounce points, bounce-averages, etc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
SUBROUTINE FIND_WELLS(na,mturn,nw0,z1,t1,B1,hBpp1,vd1,&
     &                       zb,tb,Bb,hBppb,vdb,&
     &                       z2,t2,B2,hBpp2,vd2,&
     &                       nw,alphap,offset)

!----------------------------------------------------------------------------------------------- 
!Find nw wells, from nw0 to nw0+nw-1, characterized by the Boozer toroidal and poloidal angles
!z and t, and alphap, of its top points, 1 and 2, and bottom b, as well as the value of the half
!-second derivative of the magnetic field strength along the magnetic field line hBpp, the 
!value of the magnetic drift vd, and offset
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER na,mturn,nw0
  !Input/output
  REAL*8 z1(nwx),z2(nwx),zb(nwx),hBpp1(nwx),vd1(nqv,nwx)
  REAL*8 t1(nwx),t2(nwx),tb(nwx),hBppb(nwx),vdb(nqv,nwx)
  REAL*8 B1(nwx),B2(nwx),Bb(nwx),hBpp2(nwx),vd2(nqv,nwx)
  REAL*8 alphap(nwx),offset
  !Output
  INTEGER nw
  !Others
  INTEGER ialpha,iw,flag
  REAL*8 z_ini,t_ini,dtheta0,theta0(nwx),zeta0,Bp1,Bp2,dummy,vdummy(Nnmp)
!  REAL*8, SAVE :: theta_ext=-1E3
  !Time
  CHARACTER*30, PARAMETER :: routine="FIND_WELLS"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  !Follow na field-lines along nzperiods and collapses everything into
  !nalpha field-lines in the first period
!  IF(theta_ext.LT.-200) THEN
!     IF(siota.GT.0) THEN
!        theta_ext=1E2
!     ELSE
!        theta_ext=-1E2
!     END IF
!  END IF
!  IF(NTV) THEN
!     dtheta0=TWOPI/na
  !  ELSE
  IF(na.EQ.1) THEN
     theta0=0
  ELSE
     dtheta0=iota*(TWOPI/nzperiod)/na
     DO ialpha=1,na
        theta0(ialpha)=ialpha*dtheta0
     END DO
  END IF
  
  !Calculate B_0 at (0,0) and (0,pi) and find location of maximum
  CALL CALCB(ZERO       ,ZERO,0,.FALSE.,Bp1,&
       & dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,vdummy)

  !Set starting point around a maximum
!  IF(NTV) THEN
!     CALL CALCB(PI/nzperiod,ZERO,0,.USE_B0.,Bp2,&
!          & dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy)
!     IF(Bp2.GT.Bp1) offset=PI/nzperiod
!     zeta0=offset
!  ELSE
  CALL CALCB(ZERO,PI,0,.FALSE.,Bp2,&
       & dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,vdummy)
  IF(Bp2.GT.Bp1) THEN
     offset=offset+PI
  ELSE IF(iota.LT.0) THEN
     offset=offset+TWOPI
  END IF
     
  IF(NTV) THEN
     zeta0=0!PI
  ELSE
     zeta0=0
  END IF
  theta0=theta0+offset

!  IF(siota.GT.0) THEN
!     theta_ext=MIN(theta0(1),theta_ext)
!  ELSE
!     theta_ext=MAX(theta0(1),theta_ext)
!  END IF

  iw=nw0-1  !nw0 wells have already been found, look for more
  DO ialpha=1,na
     
     IF(nw0.GT.1.AND.MOD(ialpha,2).EQ.0) CYCLE
!     IF(.FALSE..AND.na.EQ.1) THEN
!        t_ini=tmax
!        z_ini=zmax
!     ELSE
     t_ini=theta0(ialpha)
     z_ini=zeta0
     iw=iw+1
     !Find maxima and minima of the magnetic field strength along line that contains (z_ini,t_ini)
     !     CALL EXTREME_POINT(z_ini,t_ini,0, &
     CALL EXTREME_POINT(z_ini,t_ini,-1, &
          & z1(iw),t1(iw),B1(iw),hBpp1(iw),vd1(:,iw),flag)
     IF(flag.EQ.1) THEN  !Depending on flag, the extreme found was a maximum or a minimum
        zb(iw)=z1(iw)    !If mimimum, go backwards
        tb(iw)=t1(iw)
        Bb(iw)=B1(iw)
        hBppb(iw)=hBpp1(iw)
        vdb(:,iw)=vd1(:,iw)
        CALL EXTREME_POINT(zb(iw),tb(iw),-1, &
             & z1(iw),t1(iw),B1(iw),hBpp1(iw),vd1(:,iw),flag)
     ELSE
        CALL EXTREME_POINT(z1(iw),t1(iw),+1, &
             & zb(iw),tb(iw),Bb(iw),hBppb(iw),vdb(:,iw), flag)
     END IF
     CALL EXTREME_POINT(zb(iw),tb(iw),+2, &
          & z2(iw),t2(iw),B2(iw),hBpp2(iw),vd2(:,iw),flag)   
     
     DO WHILE (.TRUE.) 
        !Left extreme of new well is right extreme of the previous one
        iw=iw+1
        z1(iw)   =z2(iw-1)
        t1(iw)   =t2(iw-1)
        B1(iw)   =B2(iw-1)
        hBpp1(iw)=hBpp2(iw-1)
        vd1(:,iw)=vd2(:,iw-1)
        !Find new points until a large enough region has been covered
        !        IF(ABS(t1(iw)-theta_ext).GT.(NTURN+1)*TWOPI) THEN
        IF(ABS(t1(iw)-theta0(1)).GT.(mturn+1)*TWOPI) THEN
           iw=iw-1
           EXIT
        END IF
        CALL EXTREME_POINT(z1(iw),t1(iw),+1, &
             & zb(iw),tb(iw),Bb(iw), hBppb(iw),vdb(:,iw),flag)
        CALL EXTREME_POINT(zb(iw),tb(iw),+2, &
             & z2(iw),t2(iw),B2(iw),hBpp2(iw),vd2(:,iw),flag)
     END DO
  END DO
  nw=iw

  alphap(nw0:nw)=t1(nw0:nw)-iota*z1(nw0:nw)

  IF(DEBUG) THEN
     DO iw=nw0,nw 
        WRITE(2000+myrank,'(I6,12(1pe13.5))') iw,z1(iw),t1(iw),B1(iw),zb(iw),tb(iw),Bb(iw), &
             &     z2(iw),t2(iw),B2(iw),alphap(iw),vd1(2,iw),vd1(3,iw)
     END DO
     WRITE(2900+myrank,*) 'Found   ',nw-nw0+1,'local minima'
  END IF

  nw=nw-nw0+1

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE FIND_WELLS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE EXTREME_POINT(z_in,t_in,flag_in,z_out,t_out,B_out,hBpp_out,vd,flag_out)

!-----------------------------------------------------------------------------------------------
!Calculate next extreme point, characterized by (z,t,B,hBpp,vd), starting from (z_in,t_in)
!As input:
!-flag_in==+0 if looking for left maximum
!-flag_in==-1 if looking for left maximum going backwards
!-flag_in==+1 if looking for minimum
!-flag_in==+2 if looking for rigth maximum
!As output, flag_out<0(>0) indicates that a maximum/minimum has been found
!-----------------------------------------------------------------------------------------------   

  USE GLOBAL  
  IMPLICIT NONE
  !Input
  INTEGER flag_in
  REAL*8 z_in,t_in
  !Output
  INTEGER flag_out
  REAL*8 z_out,t_out,B_out,hBpp_out,vd(nqv)
  !Others
  LOGICAL DEL
  REAL*8 cosnm(Nnm),sinnm(Nnm),cosnm_del(Nnm),sinnm_del(Nnm)!,B
  REAL*8 z_l,dz_l,t_l,dt_l,dBdz,dBdt,dB,dBold,dummy,vdummy(Nnmp),dzstept

1 IF(DELTA) THEN
     DEL=.TRUE.
  ELSE
     DEL=.FALSE.
  END IF
  flag_out=0

  dzstept=dzstep
  dz_l=dzstept
  dt_l=dz_l*iota

  !Set direction along l (forward or backward)
  IF(flag_in.GE.0) THEN
     z_l=z_in+5*PREC_EXTR
     t_l=t_in+5*PREC_EXTR*iota
  ELSE
     z_l=z_in-5*PREC_EXTR
     t_l=t_in-5*PREC_EXTR*iota
     dz_l=-dz_l
     dt_l=-dt_l
  END IF

  IF(DEL) THEN
     CALL FILL_PHASE(z_l,t_l,cosnm,sinnm)
     CALL CALCB_DEL(cosnm,sinnm,1,.FALSE.,dummy,dBdz,dBdt,dummy,dummy,dummy,dummy,dummy,dummy,dummy)
     CALL FILL_PHASE(dz_l,dt_l,cosnm_del,sinnm_del)
  ELSE
     CALL CALCB(z_l,t_l,1,.FALSE.,dummy,dBdz,dBdt,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,vdummy)
  END IF
  dB=(dBdz+dBdt*iota)/iBtpBz  !not exactly dBdl, but proportional to it and with the same sign

  !Locate extrema of B by finding (z_l,t_l) where dBdl changes sign
  DO WHILE (ABS(dz_l).GT.PREC_EXTR)
     dBold=dB
     z_l=z_l+dz_l
     t_l=t_l+dt_l 
     IF(DEL) THEN
        CALL DELTA_PHASE(cosnm,sinnm,cosnm_del,sinnm_del)
        CALL CALCB_DEL(cosnm,sinnm,1,.FALSE.,dummy,dBdz,dBdt,dummy,dummy,dummy,dummy,dummy,dummy,dummy)
     ELSE
        CALL CALCB(z_l,t_l,1,.FALSE.,dummy,dBdz,dBdt,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,vdummy)
     END IF
     dB=(dBdz+dBdt*iota)/iBtpBz
     IF((dB-dBold)/dz_l/sgnB.GT.0) THEN !check sign of the second derivative
        flag_out=+1
     ELSE
        flag_out=-1
     END IF
     IF((dB*dBold).LE.0) THEN !after first change of sign, do not use precalculated cosines anynmore;
        dz_l=-dz_l/2.         !this is to avoid problems with numerical resolution when dz_l is very small
        dt_l=-dt_l/2.
        DEL=.FALSE.
     END IF
  END DO

  z_out=z_l 
  t_out=t_l
  !Calculate values of B, hBpp and vd at the extreme point
  CALL CALC_VDBP(z_out,t_out,B_out,dummy,hBpp_out,vd)

  IF((flag_in.EQ.+1.AND.flag_out.LT.0).OR.&   !if, when looking for maximum/minimum, found the other one
    &(flag_in.EQ.+2.AND.flag_out.GT.0)) THEN
     dzstept=dzstept/2.
     GOTO 1
  END IF

  
END SUBROUTINE EXTREME_POINT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALC_VDBP(z_in,t_in,B_out,Bp_out,hBpp_out,vd)
  
!-----------------------------------------------------------------------------------------------
!Calculate, at position z_in and t_in, the values of B_out, Bp_out, hBpp_out and vd
!-----------------------------------------------------------------------------------------------

  USE GLOBAL
  IMPLICIT NONE
  !Input  
  REAL*8 z_in,t_in
  !Output
  REAL*8 B_out,Bp_out,hBpp_out,vd(nqv)
  !Others
  REAL*8 B_0,dBdz_0,dBdt_0,dBdpsi
  REAL*8 B_1,dBdz_1,dBdt_1
  REAL*8 eta,dPhdz,dPhdt,vdummy(Nnmp)
  REAL*8 denom

  CALL CALCB(z_in,t_in,3,USE_B0,&
       & B_0,dBdz_0,dBdt_0,dBdpsi,hBpp_out,B_1,dBdz_1,dBdt_1,eta,dPhdz,dPhdt,vdummy)

  B_out =B_0                !quantities used for bounce integrals: if USE_B0, B_0 is what matters
  Bp_out=iota*dBdt_0+dBdz_0

  denom=B_0*aiBtpBz
  IF(USE_B1) THEN
     vd(1)=0.5*(Btheta*dBdz_1-Bzeta*dBdt_1)/denom
  ELSE
     vd(1)=0.5*(Btheta*dBdz_0-Bzeta*dBdt_0)/denom
     IF(USE_B0pB1) vd(1)=vd(1)+0.5*(Btheta*dBdz_1-Bzeta*dBdt_1)/denom
  END IF
  vd(2)=0.5*(iBtpBz*dBdpsi+(Bzeta*dBdt_0-Btheta*dBdz_0)*diotadpsi*z_in)/denom     !tangential v_M 
  vd(3)=B_0*B_0/avB2   !incompressibility factor

END SUBROUTINE CALC_VDBP


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE FILL_PHASE(z_l,t_l,cosnm,sinnm)

!-----------------------------------------------------------------------------------------------
!Precalculate cos(m*t_l+nzperiod*n*z_l) and sin(m*t_l+nzperiod*n*z_l) for bounce integrals
!-----------------------------------------------------------------------------------------------   
  
  USE GLOBAL
  IMPLICIT NONE
  !Input
  REAL*8 z_l,t_l
  !Output
  REAL*8 cosnm(Nnm),sinnm(Nnm)
  !Others
  INTEGER nm
  REAL*8 ang

  DO nm=1,Nnm
     ang=mp(nm)*t_l+nzperiod*np(nm)*z_l
     cosnm(nm)=COS(ang)
     sinnm(nm)=SIN(ang)
  END DO

END SUBROUTINE FILL_PHASE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE DELTA_PHASE(cosnm,sinnm,cosnm_del,sinnm_del)

!-----------------------------------------------------------------------------------------------
!Calculate cos(m*(t+dt)+nzperiod*n*(z+dz)) and sin(m*(t+dt)+nzperiod*n*(z+dz)) as a function of 
!cos(m*t+nzperiod*n*z), sin(m*t+nzperiod*n*z), cos(m*dt+nzperiod*n*dz), and sin(m*dt+nzperiod*n*dz),
!cosnm, sinnm,cosnm_del and sinnm_del respectively; write them back in cosnm and sinnm
!-----------------------------------------------------------------------------------------------

  USE GLOBAL
  IMPLICIT NONE
  !Input
  REAL*8 cosnm_del(Nnm),sinnm_del(Nnm)
  !Input/output
    REAL*8 cosnm(Nnm),sinnm(Nnm)
  !Others
  REAL*8 cosnm_temp(Nnm)
  
  cosnm_temp=cosnm*cosnm_del-sinnm*sinnm_del
  sinnm=cosnm*sinnm_del+sinnm*cosnm_del
  cosnm=cosnm_temp
  
END SUBROUTINE DELTA_PHASE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALCB_DEL(cosnm,sinnm,flag,flagB1,&
     & B_0,dBdz_0,dBdt_0,dBdpsi,B_1,dBdz_1,dBdt_1,eta,dPhdz,dPhdt)

!-----------------------------------------------------------------------------------------------
!Calculate magnetic field and derivatives at angular position where cosnm and sinnm were precalculated
!-flag.EQ.0: calculate only B
!-flag.EQ.1: calculate only first derivatives of B, dBdz and dBdt
!-flag.EQ.2: calculate B and its first derivatives, B and dBdz and dBdt
!------
!-IF(.NOT.flagB1) calculate only B_0 (usually B_1=0, so B=B_0) and its derivatives dBdz_0 and dBdt_0
!-IF(flagB1) calculate also B_1 and their derivatives dBdz_1 and dBdt_1
!------
!eta, dPhdz, and dPhdt not implemented (see older versions)
!-----------------------------------------------------------------------------------------------
  
  USE GLOBAL
  IMPLICIT NONE
  !Input
  LOGICAL flagB1
  INTEGER flag
  REAL*8 cosnm(Nnm),sinnm(Nnm)
  !Output
  REAL*8 B_0,dBdz_0,dBdt_0,dBdpsi
  REAL*8 B_1,dBdz_1,dBdt_1
  REAL*8 eta,dPhdz,dPhdt
  !Others
  INTEGER nm,nm2
  REAL*8 n,m
  REAL*8 qnmsinnm,qnmcosnm

  B_0   =0
  dBdz_0=0
  dBdt_0=0
  dBdpsi  =0
  dPhdz =0
  dPhdt =0
  eta =0
  IF(flagB1) THEN
     B_1   =0
     dBdz_1=0
     dBdt_1=0
  END IF

!!$  DO nm=1,Nnm
!!$     n=np(nm)
!!$     m=mp(nm)
!!$     IF(flag.EQ.0.OR.flag.EQ.2) THEN
!!$        B_0=B_0+bnmc0(nm)*cosnm(nm)
!!$        IF(flagB1) B_1=B_1+bnmc1(nm)*cosnm(nm)
!!$        IF(TANG_VM.AND.flag.GT.1) dBdpsi=dBdpsi+dbnmcdpsi(nm)*cosnm(nm)
!!$     END IF
!!$     IF(flag.NE.0) THEN
!!$        qnmsinnm=bnmc0(nm)*sinnm(nm)
!!$        dBdz_0=dBdz_0-qnmsinnm*n*nzperiod          
!!$        dBdt_0=dBdt_0-qnmsinnm*m
!!$        IF(flagB1) THEN
!!$           qnmsinnm=bnmc1(nm)*sinnm(nm)
!!$           dBdz_1=dBdz_1-qnmsinnm*n*nzperiod          
!!$           dBdt_1=dBdt_1-qnmsinnm*m
!!$        END IF
!!$     END IF
!!$  END DO


  DO nm=1,Nnm
     n=np(nm)
     m=mp(nm)
     IF(STELL_ANTISYMMETRIC) THEN

        IF(flag.EQ.0.OR.flag.EQ.2) THEN
           B_0=B_0+bnmc0(nm)*cosnm(nm)+bnms0(nm)*sinnm(nm)
           IF(flagB1) B_1=B_1+bnmc1(nm)*cosnm(nm)+bnms1(nm)*sinnm(nm)
           IF(TANG_VM.AND.flag.GT.1) THEN
              dBdpsi=dBdpsi+dbnmcdpsi(nm)*cosnm(nm)+dbnmsdpsi(nm)*sinnm(nm)
              eta=eta+enmc(nm)*cosnm(nm)+enms(nm)*sinnm(nm)
           END IF
        END IF
        IF(flag.NE.0) THEN
           qnmsinnm=bnmc0(nm)*sinnm(nm)
           dBdz_0=dBdz_0-qnmsinnm*n*nzperiod          
           dBdt_0=dBdt_0-qnmsinnm*m
           qnmcosnm=bnms0(nm)*cosnm(nm)
           dBdz_0=dBdz_0+qnmcosnm*n*nzperiod          
           dBdt_0=dBdt_0+qnmcosnm*m
           IF(flagB1) THEN
              qnmsinnm=bnmc1(nm)*sinnm(nm)
              dBdz_1=dBdz_1-qnmsinnm*n*nzperiod          
              dBdt_1=dBdt_1-qnmsinnm*m
              qnmcosnm=bnmc1(nm)*cosnm(nm)
              dBdz_1=dBdz_1+qnmcosnm*n*nzperiod          
              dBdt_1=dBdt_1+qnmcosnm*m
           END IF
        END IF
        
     ELSE

        IF(flag.EQ.0.OR.flag.EQ.2) THEN
           B_0=B_0+bnmc0(nm)*cosnm(nm)
           IF(flagB1) B_1=B_1+bnmc1(nm)*cosnm(nm)
           IF(TANG_VM.AND.flag.GT.1) THEN
              dBdpsi=dBdpsi+dbnmcdpsi(nm)*cosnm(nm)
              eta=eta+enmc(nm)*cosnm(nm)+enms(nm)*sinnm(nm)
           END IF
        END IF
        IF(flag.NE.0) THEN
           qnmsinnm=bnmc0(nm)*sinnm(nm)
           dBdz_0=dBdz_0-qnmsinnm*n*nzperiod          
           dBdt_0=dBdt_0-qnmsinnm*m
           IF(flagB1) THEN
              qnmsinnm=bnmc1(nm)*sinnm(nm)
              dBdz_1=dBdz_1-qnmsinnm*n*nzperiod          
              dBdt_1=dBdt_1-qnmsinnm*m
           END IF
        END IF

     END IF

     IF(PHI1_READ) THEN
        nm2=Nnm+nm
        IF(flag.EQ.0.OR.flag.EQ.2) THEN
           B_0=B_0+bnmc0(nm2)*sinnm(nm)
        END IF
        IF(flag.NE.0) THEN
           qnmcosnm=bnmc0(nm2)*cosnm(nm)
           dBdz_0=dBdz_0+qnmcosnm*n*nzperiod          
           dBdt_0=dBdt_0+qnmcosnm
        END IF
     END IF

  END DO

END SUBROUTINE CALCB_DEL


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE MATCH_WELLS(kw,&
     & z1_I  ,t1_I  ,B1_I  ,hBpp1_I  ,vd1_I  , & 
     & zb_I  ,tb_I  ,Bb_I  ,hBppb_I  ,vdb_I  , &
     & z2_I  ,t2_I  ,B2_I  ,hBpp2_I  ,vd2_I  , &
     & z1_II ,t1_II ,B1_II ,hBpp1_II ,vd1_II , &
     & zb_II ,tb_II ,Bb_II ,hBppb_II ,vdb_II , &
     & z2_II ,t2_II ,B2_II ,hBpp2_II ,vd2_II , &
     & z1_III,t1_III,B1_III,hBpp1_III,vd1_III, &
     & zb_III,tb_III,Bb_III,hBppb_III,vdb_III, &
     & z2_III,t2_III,B2_III,hBpp2_III,vd2_III, &
     & flag_in,flag_out)

!----------------------------------------------------------------------------------------------- 
!Try to match well I and II into well kw/III, characterized by the values at points
!1,b,2 of z,t,B,hBpp and vdb.
!-if flag_in<0, well III contains all I and II
!-if flag_in>0, well III only defined for B>B2_I=B1_II
!flag_out!=0 if there is a match
!-----------------------------------------------------------------------------------------------   

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER flag_in,kw
  REAL*8 z1_I  ,t1_I  ,B1_I  ,hBpp1_I  ,vd1_I(nqv)
  REAL*8 zb_I  ,tb_I  ,Bb_I  ,hBppb_I  ,vdb_I(nqv)
  REAL*8 z2_I  ,t2_I  ,B2_I  ,hBpp2_I  ,vd2_I(nqv)
  REAL*8 z1_II ,t1_II ,B1_II ,hBpp1_II ,vd1_II(nqv)
  REAL*8 zb_II ,tb_II ,Bb_II ,hBppb_II ,vdb_II(nqv)
  REAL*8 z2_II ,t2_II ,B2_II ,hBpp2_II ,vd2_II(nqv)
  !Output
  INTEGER flag_out
  REAL*8 z1_III,t1_III,B1_III,hBpp1_III,vd1_III(nqv)
  REAL*8 zb_III,tb_III,Bb_III,hBppb_III,vdb_III(nqv)
  REAL*8 z2_III,t2_III,B2_III,hBpp2_III,vd2_III(nqv)
  !Others
  REAL*8 dtl,dzl,dBl

  !Estimates "distance" between point 2 of well I and point 1 of well II
  dzl=ABS(z2_I-z1_II)
  dtl=ABS(t2_I-t1_II)
  dBl=ABS(B2_I-B1_II)
  !If
  !-point 2 is the lowest maximum of well I
  !-point 1 is the lowest maximum of well II
  !-the distance if small enough
  !match wells
  IF((B2_I.LE.B1_I).AND.(B1_II.LE.B2_II).AND.(dBl/B2_I.LT.ALMOST_ZERO).AND. &
       & (dtl/t2_I.LT.ALMOST_ZERO).AND.(dzl/z2_I.LT.ALMOST_ZERO)) THEN
     !Match
     flag_out=1
     !Point 1 of well III is point 1 of well I
     z1_III   =z1_I  
     t1_III   =t1_I  
     B1_III   =B1_I  
     hBpp1_III=hBpp1_I  
     vd1_III  =vd1_I  
     !Point 2 of well III is point 2 of well II
     z2_III   =z2_II
     t2_III   =t2_II
     B2_III   =B2_II 
     hBpp2_III=hBpp2_II 
     vd2_III  =vd2_II 
    !Determine new b point
     IF(flag_in.GT.0) THEN !Could be point 2 of well I or point 1 of well II
        zb_III   =z2_I  
        tb_III   =t2_I  
        Bb_III   =B2_I  
        hBppb_III=0.5*(hBpp2_I+hBpp1_II)
        vdb_III  =0.5*(vd2_I+vd1_II)
     ELSE
        !Use b of smallest B NOT USED AT THE MOMENT
        IF(Bb_I.LT.Bb_II) THEN
           zb_III   =zb_I  
           tb_III   =tb_I  
           Bb_III   =Bb_I  
           hBppb_III=hBppb_I  
           vdb_III  =vdb_I  
        ELSE
           zb_III   =zb_II 
           tb_III   =tb_II 
           Bb_III   =Bb_II 
           hBppb_III=hBppb_II 
           vdb_III  =vdb_II 
        END IF
     END IF
     IF(DEBUG) WRITE(2000+myrank,'(I6,12(1pe13.5))') -kw,z1_III,t1_III,B1_III,zb_III,tb_III,Bb_III, &
                                            &  z2_III,t2_III,B2_III,0.0,vd1_III(3),vd2_III(3) 
  ELSE     
     !If no match
     flag_out=0
  END IF
  
END SUBROUTINE MATCH_WELLS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE BOUNCES(iw,zx,&
                       & z1x,t1x,B1x,hBpp1x,vd1x, &
                       & zbx,tbx,Bbx,hBppbx,vdbx, &
                       & z2x,t2x,B2x,hBpp2x,vd2x, &
                       & Bbounce,top,nq,Q, &
                       & z1,t1,z2,t2)

!-----------------------------------------------------------------------------------------------
!Calculate, for well iw defined by z,t,B,hBpp,vd at tops 1 and 2 and bottom, and for lambda=1/Bbounce
!-bounce points (z1,t1) and (z2,t2)
!-several bounce integrals Q
!If top, one bounce point is 1 or 2
!-----------------------------------------------------------------------------------------------
  
  USE GLOBAL
  IMPLICIT NONE
  !Input
  LOGICAL top
  INTEGER iw,nq
  REAL*8 zx
  REAL*8 zbx,tbx,Bbx,hBppbx,vdbx(nqv)
  REAL*8 z1x,t1x,B1x,hBpp1x,vd1x(nqv)
  REAL*8 z2x,t2x,B2x,hBpp2x,vd2x(nqv)
  REAL*8 Bbounce
  !Output
  REAL*8 Q(nq),z1,z2,t1,t2
  !Others
  LOGICAL topl,topr
  INTEGER it
  REAL*8 Bp1,hBpp1,vd1(nqv)
  REAL*8 Bp2,hBpp2,vd2(nqv),NAN
  !Time
  CHARACTER*30, PARAMETER :: routine="BOUNCES"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  NAN=0.0 !avoid 1/0
    
  !Find bounce points
  topl=top.AND.(B1x.LE.B2x)  !if one bounce point is very close to 1
  IF(topl) THEN              !BOUNCE_POINT there may be numerical problems
     z1=z1x                  !-> set by hand
     t1=t1x
     Bp1=0
     hBpp1=hBpp1x
     vd1=vd1x
  ELSE
     CALL BOUNCE_POINT(zbx,tbx,Bbounce,z1x,t1x,z1,t1,Bp1,hBpp1,vd1,-1)
  END IF

  topr=top.AND.(B2x.LE.B1x) !same with 2
  IF(topr) THEN
     z2=z2x
     t2=t2x
     Bp2=0
     hBpp2=hBpp2x
     vd2=vd2x
  ELSE
     CALL BOUNCE_POINT(zbx,tbx,Bbounce,z2x,t2x,z2,t2,Bp2,hBpp2,vd2,1)
  END IF

  IF(DEBUG) WRITE(2100+myrank,'(I6,11(1pe13.5))') &
       & ABS(iw),z1 ,t1 ,Bbounce,zbx,tbx,NaN, & 
       &         z2 ,t2 ,Bbounce,vd1(3),vd2(3)                             

  it=1
  DO WHILE((it.EQ.1.OR.ISNAN(Q(1))).AND.it.LE.10)
  !Calculate bounce integrals
     CALL BOUNCE_INTEGRAL(iw,zx,&
          &     z1,t1,z2,t2,1./Bbounce, &
          &             Bp1,hBpp1,vd1,  &
          &             Bp2,hBpp2,vd2,  &
          &             zbx,bbx,hBppbx,vdbx,nq,Q)
     it=it+1
     IF(ISNAN(Q(1))) WRITE(6200+myrank,*) 'NAN',z1,z2,Bbounce
  END DO
  
  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)
  
END SUBROUTINE BOUNCES


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE BOUNCE_POINT(z_in,t_in,Bbounce,z_lim,t_lim, &
     & z_out,t_out,Bp_out,hBpp_out,vd_out,flag)

!-----------------------------------------------------------------------------------------------
!Calculate bounce point, characterized by (z_out,t_out,Bp_out,hBpp_out and vd_out) by looking 
!for B=Bbounce, going forward/backward if flag>1/<1, between (z_in,t_in) and (z_lim,t_lim)
!-----------------------------------------------------------------------------------------------   

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER flag
  REAL*8 z_in,t_in,z_lim,t_lim,Bbounce
  !Output
  REAL*8 z_out,t_out,Bp_out,hBpp_out,vd_out(nqv)
  !Others
  LOGICAL DEL
  REAL*8 cosnm(Nnm),sinnm(Nnm),cosnm_del(Nnm),sinnm_del(Nnm)
  REAL*8 z_l,dz_l,t_l,dt_l,dz_l_ini
  REAL*8 B,B_old,dummy,vdummy(Nnmp)
  !Time
!  CHARACTER*30, PARAMETER :: routine="BOUNCE_POINT"
!  INTEGER, SAVE :: ntotal=0
!  REAL*8,  SAVE :: ttotal=0
!  REAL*8,  SAVE :: t0=0
!  REAL*8 tstart

!  CALL CPU_TIME(tstart)

  IF(DELTA) THEN
     DEL=.TRUE.
  ELSE
     DEL=.FALSE.
  END IF

  z_l=z_in
  t_l=t_in

  dz_l_ini=1.01*ABS(z_lim-z_in)/2
  dz_l=dz_l_ini
  dt_l=dz_l*iota
  
  !Set direction along l (forward  or backward)
  IF(flag.LT.0) THEN
     dz_l=-dz_l
     dt_l=-dt_l
   END IF

  IF(DEL) THEN
     CALL FILL_PHASE(z_l,t_l,cosnm,sinnm)
     CALL CALCB_DEL(cosnm,sinnm,0,.FALSE.,B,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy)
     CALL FILL_PHASE(dz_l,dt_l,cosnm_del,sinnm_del)
  ELSE
     CALL CALCB(z_l,t_l,0,.FALSE.,B,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,vdummy)
  END IF

  !Locate extrema of B by finding (z_l,t_l) where B-Bbounce changes sign
  DO WHILE (ABS(dz_l).GT.PREC_EXTR)
     B_old=B
     z_l=z_l+dz_l
     t_l=t_l+dt_l      !t_lim and z_lim are probably redundant
     IF ((ABS(t_l-t_in).GT.ABS(t_lim-t_in)).OR.(ABS(z_l-z_in).GT.ABS(z_lim-z_in))) THEN
        z_l=z_l-dz_l     !do not go beyond (z_lim,t_lim)
        t_l=t_l-dt_l     !(there, no particles for small Bbounce, other wells for large Bbounce)
        dz_l=z_lim-z_l
        dt_l=t_lim-t_l
        CALL FILL_PHASE(dz_l,dt_l,cosnm_del,sinnm_del)
        z_l=z_l+dz_l
        t_l=t_l+dt_l        
     END IF
     IF ((ABS(t_l-t_lim).GT.ABS(t_lim-t_in)).OR.(ABS(z_l-z_lim).GT.ABS(z_lim-z_in))) THEN
        z_l=z_l-dz_l     !do not go beyond (z_in,t_in)
        t_l=t_l-dt_l     !(that region is explored with opposite sign of flag)
        dz_l=z_in-z_l
        dt_l=t_in-t_l
        CALL FILL_PHASE(dz_l,dt_l,cosnm_del,sinnm_del)
        z_l=z_l+dz_l
        t_l=t_l+dt_l
     END IF 
     IF(DEL) THEN
        CALL DELTA_PHASE(cosnm,sinnm,cosnm_del,sinnm_del)
        CALL CALCB_DEL(cosnm,sinnm,0,.FALSE.,B,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy) 
     ELSE
        CALL CALCB(z_l,t_l,0,.FALSE.,B,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,vdummy)
     END IF
     IF((B-Bbounce)*(B_old-Bbounce).LE.0) THEN
        dz_l=-dz_l/2.
        dt_l=-dt_l/2.
!        IF(ABS(dz_l).LT.PREC_EXTR) DEL=.FALSE.
        IF(DEL) THEN
           CALL FILL_PHASE(z_l,t_l,cosnm,sinnm)
           CALL CALCB_DEL(cosnm,sinnm,0,.FALSE.,B,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy)
           CALL FILL_PHASE(dz_l,dt_l,cosnm_del,sinnm_del)  
        END IF
     END IF
  END DO
  z_out=z_l
  t_out=t_l

  !Calculate values of Bp, hBpp and vd at bounce point
  CALL CALC_VDBP(z_out,t_out,dummy,Bp_out,hBpp_out,vd_out)

!  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal) 

END SUBROUTINE BOUNCE_POINT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE BOUNCE_INTEGRAL(iw,zx,z_ini,t_ini,z_fin,t_fin,lambd, &
     &      Bp_ini,hBpp_ini,vd_ini, &
     &      Bp_fin,hBpp_fin,vd_fin, &
     & z_bot,B_bot,hBpp_bot,vd_bot,nq,Q)

!-----------------------------------------------------------------------------------------------
!Calculate bounce integrals Q for well iw and lambda lambd.
!The bounce points (z_ini,t_ini) and (z_fin,t_fin) ,characterized by Bp, hBp and vd, and 
!the bottom of the well z_bot, characterized by B, hBpp and vd, are also inputs
!-----------------------------------------------------------------------------------------------   

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER iw,nq
  REAL*8 z_ini,t_ini,Bp_ini,hBpp_ini,vd_ini(nqv)
  REAL*8 z_fin,t_fin,Bp_fin,hBpp_fin,vd_fin(nqv)
  REAL*8 z_bot,      B_bot,hBpp_bot,vd_bot(nqv)
  REAL*8 zx,lambd
  !Output
  REAL*8 Q(nq)
  !Parameters
  INTEGER, PARAMETER :: nmin=3  !calculate with at least 9 points
  INTEGER, PARAMETER :: nmax=20 !large enough (it means ~10^9 points, but PREC_EXTR reached before)
  !Others
  INTEGER nint,iq,ifrac,nfrac
  REAL*8 dzl,ddzl,tdzl,tdzlo3,z_l
  REAL*8 dtl,ddtl,tdtl,tdtlo3,t_l
  REAL*8 cosnm(Nnm),cosnm_ini(Nnm),cosnm_del(Nnm),cosnm_del2(Nnm),cosnm_del4(Nnm)
  REAL*8 sinnm(Nnm),sinnm_ini(Nnm),sinnm_del(Nnm),sinnm_del2(Nnm),sinnm_del4(Nnm)
  REAL*8 Qold(nq),Qsum(nq),Qint(nq),Qana(nq0),deltas,maxdeltas,mindeltas
  REAL*8, ALLOCATABLE :: ds(:)
  !Time
!  CHARACTER*30, PARAMETER :: routine="BOUNCE_INTEGRAL"
!  INTEGER, SAVE :: ntotal=0
!  REAL*8,  SAVE :: ttotal=0
!  REAL*8,  SAVE :: t0=0
!  REAL*8 tstart

!  CALL CPU_TIME(tstart)

  !Calculate analytically logarithmic divergence
  IF(REMOVE_DIV) THEN
     Qana=0
     CALL ANA_INTEGRAL(iw,z_fin-z_ini,lambd,ZERO ,Bp_ini,hBpp_ini,vd_ini,nq0,Qana)
     CALL ANA_INTEGRAL(iw,z_ini-z_fin,lambd,ZERO ,Bp_fin,hBpp_fin,vd_fin,nq0,Qana)
     CALL ANA_INTEGRAL(iw,z_bot-z_ini,lambd,B_bot,ZERO  ,hBpp_bot,vd_bot,nq0,Qana)
     CALL ANA_INTEGRAL(iw,z_fin-z_bot,lambd,B_bot,ZERO  ,hBpp_bot,vd_bot,nq0,Qana)
  END IF
  
  !Start with the central point
  z_l=0.5*(z_ini+z_fin)
  t_l=0.5*(t_ini+t_fin)
  CALL FILL_PHASE(z_l,t_l,cosnm,sinnm)
  !Calculate integrand
  CALL BOUNCE_INTEGRAND(iw,z_ini,z_l,t_l,cosnm,sinnm,lambd,nq,Qint)
  
  !First calculation removing the divergence
  IF(REMOVE_DIV) THEN 
     CALL BOUNCE_INTEGRAND_MINF(iw,z_l,z_ini,lambd,ZERO  ,Bp_ini,hBpp_ini,vd_ini,nq0,Qint(1:nq0)) 
     CALL BOUNCE_INTEGRAND_MINF(iw,z_l,z_fin,lambd,ZERO  ,Bp_fin,hBpp_fin,vd_fin,nq0,Qint(1:nq0)) 
     CALL BOUNCE_INTEGRAND_MINF(iw,z_l,z_bot,lambd,B_bot ,ZERO  ,hBpp_bot,vd_bot,nq0,Qint(1:nq0)) 
  END IF
  tdzl=(z_fin-z_ini)
  tdtl=(t_fin-t_ini)
  Q=qint*tdzl

  tdzlo3=tdzl/3.
  tdtlo3=tdtl/3.
  nfrac=1
  Qold=Q
  IF(DELTA) CALL FILL_PHASE(z_ini,t_ini,cosnm_ini,sinnm_ini)
  !Triple the number of points until convergence is reached
  DO nint=2,nmax 
     dzl=tdzlo3/nfrac
     dtl=tdtlo3/nfrac
     ddzl=dzl+dzl
     ddtl=dtl+dtl
     IF(DELTA) THEN
        cosnm=cosnm_ini
        sinnm=sinnm_ini
        CALL FILL_PHASE(0.5*dzl,0.5*dtl,cosnm_del,sinnm_del)
        cosnm_del2=cosnm_del*cosnm_del-sinnm_del*sinnm_del
        sinnm_del2=2*cosnm_del*sinnm_del
        cosnm_del4=cosnm_del2*cosnm_del2-sinnm_del2*sinnm_del2
        sinnm_del4=2*cosnm_del2*sinnm_del2
        CALL DELTA_PHASE(cosnm,sinnm,cosnm_del,sinnm_del)
     END IF
     z_l=z_ini+0.5*dzl
     t_l=t_ini+0.5*dtl
     Qsum=0
     ALLOCATE(ds(nfrac))
     ds=0
     DO ifrac=1,nfrac
        CALL BOUNCE_INTEGRAND(iw,z_ini,z_l,t_l,cosnm,sinnm,lambd,nq,Qint)
        IF(ISNAN(Qint(1))) EXIT
        CALL DELTA_PHASE(cosnm,sinnm,cosnm_del4,sinnm_del4)
        IF(REMOVE_DIV) THEN
           CALL BOUNCE_INTEGRAND_MINF(iw,z_l,z_ini,lambd,MONE  ,Bp_ini,hBpp_ini,vd_ini,nq0,Qint(1:nq0)) 
           CALL BOUNCE_INTEGRAND_MINF(iw,z_l,z_fin,lambd,MONE  ,Bp_fin,hBpp_fin,vd_fin,nq0,Qint(1:nq0)) 
           CALL BOUNCE_INTEGRAND_MINF(iw,z_l,z_bot,lambd,B_bot ,ZERO  ,hBpp_bot,vd_bot,nq0,Qint(1:nq0)) 
        END IF
        z_l=z_l+ddzl 
        t_l=t_l+ddtl
        Qsum=Qsum+qint
        CALL BOUNCE_INTEGRAND(iw,z_ini,z_l,t_l,cosnm,sinnm,lambd,nq,Qint)
        IF(ISNAN(Qint(1))) EXIT
        ds(ifrac)=Qint(3)
        CALL DELTA_PHASE(cosnm,sinnm,cosnm_del2,sinnm_del2)
        IF(REMOVE_DIV) THEN
           CALL BOUNCE_INTEGRAND_MINF(iw,z_l,z_ini,lambd,MONE  ,Bp_ini,hBpp_ini,vd_ini,nq0,Qint(1:nq0)) 
           CALL BOUNCE_INTEGRAND_MINF(iw,z_l,z_fin,lambd,MONE  ,Bp_fin,hBpp_fin,vd_fin,nq0,Qint(1:nq0)) 
           CALL BOUNCE_INTEGRAND_MINF(iw,z_l,z_bot,lambd,B_bot ,ZERO  ,hBpp_bot,vd_bot,nq0,Qint(1:nq0)) 
        END IF
        z_l=z_l+dzl 
        t_l=t_l+dtl
        Qsum=Qsum+qint
     END DO
     IF(ISNAN(Qint(1))) EXIT
     Q=(Qold+tdzl*Qsum/nfrac)/3.
     IF(DEBUG.AND.(iw.EQ.I0.OR.I0.EQ.0)) THEN
        WRITE(2400+myrank,'(I6,7(1pe13.5))') iw,dzl,lambd,(Q(iq),iq=1,nq0)
        WRITE(2500+myrank,'(I6,7(1pe13.5))') iw,dzl,lambd,(Q(iq)+Qana(iq),iq=1,nq0)
     END IF
     
     deltas=0
     maxdeltas=0
     mindeltas=0
     DO ifrac=1,nfrac
        deltas=deltas+ds(ifrac)*tdzl/nfrac
        IF(deltas.GT.maxdeltas) maxdeltas=deltas
        IF(deltas.LT.mindeltas) mindeltas=deltas
     END DO
     DO ifrac=nfrac,1,-1
        deltas=deltas+ds(ifrac)*tdzl/nfrac
        IF(deltas.GT.maxdeltas) maxdeltas=deltas
        IF(deltas.LT.mindeltas) mindeltas=deltas
     END DO
     DEALLOCATE(ds)
 
     !Convergence is checked looking at the total integral (Q+Qana)
!     IF(FAST_IONS) THEN
!        IF(nint.GE.nmin.AND.(dzl.LT.PREC_EXTR.OR.&
!             & (ABS(1-(Q(6)+Qana(6))/(Qold(6)+Qana(6))).LT.PREC_BINT))) EXIT
     IF(nint.GE.nmin.AND.(dzl.LT.PREC_EXTR.OR.&
           & (ABS(1-(Q(1)+Qana(1))/(Qold(1)+Qana(1))).LT.PREC_BINT.AND.&
           & (ABS(1-(Q(2)+Qana(2))/(Qold(2)+Qana(2))).LT.PREC_BINT)))) EXIT
     Qold=Q
     nfrac=nfrac*3
!!     
  END DO

  IF(ISNAN(Qint(1))) THEN
     IF((z_ini-z_l)*(z_l-zx).GE.0) THEN
        z_ini=z_l
        t_ini=t_l
     ELSE
        z_fin=z_l
        t_fin=t_l
     END IF
     Q=Qint
     RETURN
  END IF

  IF(FAST_IONS) Q(2)=maxdeltas-mindeltas

  !Warnings
  IF(dzl.LT.PREC_EXTR) WRITE(1100+myrank,*) 'Bounce integral not converged, try reducing DTMIN'
  IF(nint.GE.nmax) WRITE(1100+myrank,*) 'Bounce integral not converged, try increasing NMAX!!'

  IF(DEBUG.AND.(iw.EQ.I0.OR.I0.EQ.0)) THEN
     WRITE(2400+myrank,'(I6,7(1pe13.5))') iw,-dzl,lambd,(Q(iq),iq=1,nq0)
     WRITE(2500+myrank,'(I6,7(1pe13.5))') iw,-dzl,lambd,(Q(iq)+Qana(iq),iq=1,nq0)
  END IF

  !Sum contributions
  Q=Q+Qana
!  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE BOUNCE_INTEGRAL


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE BOUNCE_INTEGRAND(iw,z_ini,z_l,t_l,cosnm,sinnm,lambd,nq,Qint)

!-----------------------------------------------------------------------------------------------
!Calculate nq integrands Qint at point (z_l,t_l) at well iw for lambda lambd
!using (if DELTA) precalculated cosnm and sinm. z_ini is relevant for the secular term of 
!d_psi J that is proportional to the magnetic shear
!-----------------------------------------------------------------------------------------------

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER iw,nq
  REAL*8 z_ini,z_l,t_l
  REAL*8 cosnm(Nnm),sinnm(Nnm)
  REAL*8 lambd
  !Output
  REAL*8 Qint(nq)
  !Others
  CHARACTER*100 serr
  INTEGER iq
  REAL*8 B_0,dBdz_0,dBdt_0,dBdpsi,dummy,vdummy(Nnmp)
  REAL*8 B_1,dBdz_1,dBdt_1
  REAL*8 eta,dPhdz,dPhdt 
  REAL*8 lambdaB0,sqrt1mlb
  REAL*8 vds,vda,denom,factB,factnm(Nnm)
  !Time
!  CHARACTER*30, PARAMETER :: routine="BOUNCE_INTEGRAND"
!  INTEGER, SAVE :: ntotal=0
!  REAL*8,  SAVE :: ttotal=0
!  REAL*8,  SAVE :: t0=0
!  REAL*8 tstart

!  CALL CPU_TIME(tstart)

  Qint=0
  IF(DELTA) THEN
     CALL CALCB_DEL(cosnm,sinnm,2,USE_B0,B_0,dBdz_0,dBdt_0,dBdpsi,B_1,dBdz_1,dBdt_1,eta,dPhdz,dPhdt)
  ELSE
     CALL CALCB(z_l,t_l,2,USE_B0,B_0,dBdz_0,dBdt_0,dBdpsi,dummy,B_1,dBdz_1,dBdt_1,eta,dPhdz,dPhdt,vdummy)
  END IF
  lambdaB0=lambd*B_0
  sqrt1mlb=SQRT(1.-lambdaB0)
  IF(lambdaB0.GE.1) THEN
!     IF(lambdaB0-1.GT.1E-4) THEN
     serr="lambda*B>1"
!     WRITE(iout,*) 'lambda*B>1'
!        CALL END_ALL(serr,.FALSE.)
!     END IF
     RETURN
  END IF
  denom=aiBtpBz*B_0
  factB=(1.0-0.5*lambdaB0)/denom
  !  vda=factB*(iBtpBz*dBdpsi+(Bzeta*dBdt_0-Btheta*dBdz_0)*diotadpsi*(z_l-z_ini))
  vda=iBtpBz*dBdpsi
  vda=vda+(Bzeta*dBdt_0-Btheta*dBdz_0)*diotadpsi*(z_l-z_ini)
  vda=vda-(iota*dBdt_0+dBdz_0)*eta
  vda=vda+((1.0-lambdaB0)/(1.0-0.5*lambdaB0))*(iBtpBz/B_0)*dmu0Pdpsi
  vda=vda*factB
  IF(USE_B1) THEN
     vds=factB*(Btheta*dBdz_1-Bzeta*dBdt_1) 
  ELSE
     vds=factB*(Btheta*dBdz_0-Bzeta*dBdt_0) 
     IF(USE_B0pB1) vds=vds+factB*(Btheta*dBdz_1-Bzeta*dBdt_1) 
  END IF
  IF(.NOT.PHI1_READ.OR..NOT.IMP1NU) THEN
     Qint(1)=1./sqrt1mlb              ! tau
     Qint(2)=lambd*sqrt1mlb/B_0       ! vpar
     Qint(3)=vds/sqrt1mlb             ! radial magnetic drift
!     IF(USE_B1) Qint(3)=factB*(Btheta*dBdz_1-Bzeta*dBdt_1)     
     Qint(4)=vda/sqrt1mlb             ! tangential magnetic drift
     IF(INC_EXB) THEN
        Qint(5)=(B_0*B_0/avB2)/sqrt1mlb ! factor for incompressible ExB tangential drift
     ELSE
        Qint(5)=(B_0/avB)/sqrt1mlb ! ExB as in d'Herbemont et al.
     END IF
     Qint(6)=sqrt1mlb                 ! J/2v
     Qint(7)=B_0/sqrt1mlb/2/sgnB      ! -(dJ/dlambda)/2v
!     Qint(7)=B_0/sqrt1mlb/2      ! -(dJ/dlambda)/2v
     IF(SOLVE_QN) THEN                ! contribution of each mode to radial ExB drift
        factnm(1:Nnm)=(Btheta*np(1:Nnm)*nzperiod-Bzeta*mp(1:Nnm))/sqrt1mlb/aiBtpBz
        Qint(8    :7+Nnm )=-factnm(1:Nnm)*sinnm(1:Nnm)
        Qint(8+Nnm:7+Nnmp)=+factnm(1:Nnm)*cosnm(1:Nnm)
     END IF
  ELSE
     vds=((Btheta* dBdz_1-        Bzeta* dBdt_1        )*(2+lambd*(B_1-2*B_0))/lambd/B_1-&
&         (Btheta*(dBdz_0-dBdz_1)-Bzeta*(dBdt_0-dBdt_1)))/aiBtpBz
     Qint(1)=B_0/sqrt1mlb                
     Qint(2)=lambd*SQRT(lambd)*sqrt1mlb  
     Qint(3)=vds/sqrt1mlb 
     Qint(4)=0.0
  END IF

  Qint(1:nq)=Qint(1:nq)*aiBtpBz/B_0   !dl=dz*(dz/dl) 
  IF(DEBUG.AND.(iw.EQ.L0.OR.L0.EQ.0)) & 
       & WRITE(2200+myrank,'(I6,7(1pe13.5))') iw,z_l,(Qint(iq),iq=1,nq0),lambd

!  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)
 
END SUBROUTINE BOUNCE_INTEGRAND


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BOUNCE_INTEGRAND_MINF(iw,z_l,zb,lambda,Bb,Bpb,hBppb,vdb,nq,Qint)

!-----------------------------------------------------------------------------------------------
!Remove infinites from qn integrands Qint, corresponding to point (z_l,t_l) at well iw and lambda. 
!Use as input values of the magnetic field and derivatives, Bb, Bpb and hBppb at the bounce
!point zb, Bb,Bpb
!-----------------------------------------------------------------------------------------------   
  
  USE GLOBAL
  IMPLICIT NONE

  !Input
  INTEGER iw,nq
  REAL*8 z_l,zb,lambda,Bb,Bpb,hBppb,vdb(nqv)
  !Input/output
  REAL*8 Qint(nq)
  !Others
  REAL*8 dzl,denom,qdiv,dlambda


  dzl=z_l-zb
  IF(Bb.GT.ALMOST_ZERO) THEN         !if (z_l,t_l) is an extreme of B
     IF(hBppb.LE.0) THEN                 !if (z_l,t_l) is a maximum of B
        dlambda=lambda-1./Bb
        denom=SQRT(ABS(dlambda)*Bb+ABS(hBppb)/Bb*dzl*dzl)
     ELSE
        RETURN
     END IF
  ELSE
     IF(hBppb.LE.0) THEN
        denom=SQRT(-lambda*dzl*(Bpb+hBppb*dzl))
     ELSE
        RETURN
!        denom=SQRT(-lambda*dzl*Bpb)
     END IF
  END IF
  
  qdiv=aiBtpBz*lambda/denom

  Qint(1)=Qint(1)-qdiv
  Qint(3)=Qint(3)-qdiv*vdb(1)
  Qint(4)=Qint(4)-qdiv*vdb(2)
  Qint(5)=Qint(5)-qdiv*vdb(3)

  IF(DEBUG.AND.(iw.EQ.L0.OR.L0.EQ.0)) WRITE(2300+myrank,'(I6,10(1pe13.5))') &
  & iw,z_l,qdiv,qdiv-qdiv,qdiv*vdb(1),qdiv*vdb(2),qdiv*vdb(3),lambda

END SUBROUTINE BOUNCE_INTEGRAND_MINF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE ANA_INTEGRAL(iw,dzb,lambda,Bb,Bpb,hBppb,vdb,nq,Qana)

!-----------------------------------------------------------------------------------------------
!Calculate analytically the divergent contribution Qana to nq bounce integrals at lambda, well iw; 
!dzb is the size of the integration region and Bb, Bpb, hBppb and vdb the main quantities 
!at the bounce point.
!-----------------------------------------------------------------------------------------------   
  
  USE GLOBAL   
  IMPLICIT NONE
  !Input
  INTEGER iw,nq
  REAL*8 dzb,lambda,Bb,Bpb,hBppb,vdb(nqv)
  !Output
  REAL*8 Qana(nq)
  !Others
  REAL*8 Ib,dlambda,dlambdaBb2ohBppb

  IF(Bb.GT.ALMOST_ZERO) THEN         !if (z_l,t_l) is an extreme of B
     IF(hBppb.LE.0) THEN             !if (z_l,t_l) is a maximum of B
        dlambda=lambda-1./Bb
        dlambdaBb2ohBppb=Bb*Bb*ABS(dlambda/hBppb)
        Ib=SQRT(Bb/ABS(hBppb))*LOG((dzb+SQRT(dzb*dzb+dlambdaBb2ohBppb))/SQRT(dlambdaBb2ohBppb))
     ELSE
        RETURN
     END IF
  ELSE
     IF(hBppb.LE.0) THEN 
        Ib=SQRT(-1./hBppb/lambda)*LOG((2*lambda*(SQRT(hBppb*dzb*(Bpb+hBppb*dzb)))-lambda*(hBppb*2*dzb+Bpb)) &
             & /(-lambda*Bpb))
        IF(dzb.LT.0) Ib=-Ib
     ELSE
        RETURN
!        Ib=2*SQRT(-dzb/(lambda*Bpb))     
     END IF
  END IF

  Ib=Ib*aiBtpBz*lambda

  Qana(1)=Qana(1)+Ib
  Qana(3)=Qana(3)+Ib*vdb(1)
  Qana(4)=Qana(4)+Ib*vdb(2)
  Qana(5)=Qana(5)+Ib*vdb(3)

  IF(DEBUG.AND.(iw.EQ.L0.OR.L0.EQ.0)) WRITE(2600+myrank,'(I6,10(1pe13.5))') iw,lambda,Ib

END SUBROUTINE ANA_INTEGRAL


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

