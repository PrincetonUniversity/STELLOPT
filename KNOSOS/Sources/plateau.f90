!Calculate neoclassical transport in the plateau regime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALC_PLATEAU(jv,Epsi,D11,dn1nm)

!----------------------------------------------------------------------------------------------- 
!Calculate monoenergetic transport coefficient D11 in the plateau regime for
!collisionality cmul=nu(jv)/v(jv) and normalized radial electric field efied=Epsi/v(jv)
!-----------------------------------------------------------------------------------------------
  
  USE GLOBAL
  USE KNOSOS_STELLOPT_MOD
  IMPLICIT NONE
  !Input
  INTEGER jv
  REAL*8 Epsi
  !Output
  REAL*8 D11(Nnmp,Nnmp),dn1nm(Nnmp,Nnmp)
  !Others
  INTEGER m,n,nm
  REAL(rprec) phrbic(-ntorbd:ntorbd,0:mpolbd),phrbis(-ntorbd:ntorbd,0:mpolbd)
  COMPLEX*16 bmn_e(-ntorbd:ntorbd,-mpolbd:mpolbd),pmn_e(-ntorbd:ntorbd,-mpolbd:mpolbd)
  REAL*8 nBtmmBz,anpim,G11dkes,cmul,D11t(Nnmp,Nnmp)
  REAL*8, SAVE :: D11p=1

  WRITE(iout,*) 'Calculating PLATEAU'

  IF(D11p.LT.0.AND..NOT.QN) GO TO 999
  
  !Change Fourier representation
  bmn_e=0
  DO m=-mpolb,mpolb
     DO n=-ntorb,ntorb
        IF(m.LT.0) bmn_e(n,m)= borbic(-n,-m)/2.
        IF(m.GT.0) bmn_e(n,m)= borbic( n, m)/2.
        IF(m.EQ.0) bmn_e(n,0)=(borbic( n,0)+borbic(-n,0))/2.
        IF(STELL_ANTISYMMETRIC) THEN
           IF(m.LT.0) bmn_e(n,m)=bmn_e(n,m)-NUMI*borbis(-n,-m)/2.
           IF(m.GT.0) bmn_e(n,m)=bmn_e(n,m)+NUMI*borbis( n, m)/2.
           IF(m.EQ.0) bmn_e(n,0)=bmn_e(n,0)+NUMI*(borbis( n,0)-borbic(-n,0))/2.
        END IF
     END DO
  END DO
  pmn_e=0
  IF(QN) THEN
     phrbic=0
     phrbis=0
     phrbic(:,0:mpolbd)=1
     phrbis(:,0:mpolbd)=1
     DO m=-mpolb,mpolb
        DO n=-ntorb,ntorb
           IF(m.LT.0) pmn_e(n,m)= (phrbic(-n,-m)-NUMI*phrbis(-n,-m))/2.
           IF(m.GT.0) pmn_e(n,m)= (phrbic( n, m)+NUMI*phrbis( n, m))/2.
           IF(m.EQ.0) pmn_e(n,0)=((phrbic( n,0)+phrbic(-n,0))+NUMI*(phrbis( n,0)-phrbic(-n,0)))/2.
        END DO
     END DO
  END IF
  D11t=0
  dn1nm=0
  DO nm=1,2*Nnm
     n=INT(np(nm))
     m=INT(mp(nm))
     IF(n.EQ.0.AND.m.EQ.0) CYCLE
     nBtmmBz=n*nzperiod*Btheta-m*Bzeta     
     anpim=ABS(n*nzperiod+iota*m)
     IF(anpim.LT.2E-1) CYCLE
     D11t(1,1)=D11t(1,1)+nBtmmBz*nBtmmBz*REAL(bmn_e(n,m)*bmn_e(-n,-m))/anpim
     IF(nm.NE.1.AND.nm.NE.Nnm+1.AND.QN) THEN
        IF(nm.LE.Nnm) THEN
           dn1nm(nm+Nnm, 1)=-2*nBtmmBz*REAL(bmn_e(n,m))/anpim
           dn1nm(nm+Nnm,nm)=-2*nBtmmBz*REAL(pmn_e(n,m))/anpim
        ELSE
           dn1nm(nm-Nnm, 1)=2*nBtmmBz*AIMAG(bmn_e(n,m))/anpim
           dn1nm(nm-Nnm,nm)=2*nBtmmBz*AIMAG(pmn_e(n,m))/anpim
        END IF
        D11t(nm,1)=        +nBtmmBz*nBtmmBz*REAL(pmn_e(n,m)*bmn_e(-n,-m))/anpim
        D11t(1,nm)=        +nBtmmBz*nBtmmBz*REAL(bmn_e(n,m)*pmn_e(-n,-m))/anpim
        D11t(nm,nm)=       +nBtmmBz*nBtmmBz*REAL(pmn_e(n,m)*pmn_e(-n,-m))/anpim
     END IF
     IF(.NOT.QN.AND.nm.EQ.Nnm) THEN
        D11p=D11t(1,1)       
        EXIT
     END IF
  END DO
  
999 IF(.NOT.QN) D11t=D11p
  D11t(1,:)=D11t(1,:)*vdconst(jv)/borbic(0,0)/2.
  D11t(:,1)=D11t(:,1)*vdconst(jv)/borbic(0,0)/2.
  dn1nm(:,1)=dn1nm(:,1)*vdconst(jv)/borbic(0,0)/2.
        
  D11t=pi*(rad_R/aiBtpBz/aiBtpBz/v(jv))*D11t
  dn1nm=dn1nm*(-pi*rad_R/v(jv)/aiBtpbz)
  cmul=nu(jv)/v(jv)/2.
  !Connect with Pfirsch-Schlueter or 1/nu
  IF(FACT_CON.GT.0.AND.cmul_PS.GT.0) THEN
     IF(cmul.GT.cmul_PS /FACT_CON) D11t(1,1)=D11t(1,1)*(1+cmul/cmul_PS)
     IF(cmul.LT.cmul_1NU*FACT_CON) D11t(1,1)=D11t(1,1)*(1+cmul_1NU/cmul)
  END IF

  G11dkes=fdkes(jv)*D11t(1,1)
  IF(.NOT.KNOSOS_STELLOPT) WRITE(200+myrank,'(3(1pe13.5)," NaN ",2(1pe13.5),"&
       & NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN")') &
       & nu(jv)/v(jv)/2,Epsi/v(jv)*psip,vmconst(jv)/v(jv),G11dkes,G11dkes
  IF(DEBUG) THEN
     IF(cmul_PS.GT.0) THEN !to be checked
        WRITE(10000+myrank,'("2 ",6(1pe13.5))') nu(jv)/v(jv)/2,Epsi/v(jv)*psip,vmconst(jv)/v(jv),&
             & G11dkes,weight(jv)/fdkes(jv),weight(jv)/fdkes(jv)*v(jv)*v(jv)
     ELSE
        WRITE(10000+myrank,'("0 ",6(1pe13.5))') nu(jv)/v(jv)/2,Epsi/v(jv)*psip,vmconst(jv)/v(jv),&
             & G11dkes,weight(jv)/fdkes(jv),weight(jv)/fdkes(jv)*v(jv)*v(jv)
     END IF
  END IF

  D11=D11t

END SUBROUTINE CALC_PLATEAU


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!$SUBROUTINE CALC_PLATEAU_OLD(jv,Epsi,D11)
!!$
!!$!----------------------------------------------------------------------------------------------- 
!!$!Calculate monoenergetic transport coefficient D11 in the plateau regime for
!!$!collisionality cmul=nu(jv)/v(jv) and normalized radial electric field efied=Epsi/v(jv)
!!$!-----------------------------------------------------------------------------------------------
!!$  
!!$  USE GLOBAL
!!$  IMPLICIT NONE
!!$  !Input
!!$  INTEGER jv
!!$  REAL*8 Epsi
!!$  !Output
!!$  REAL*8 D11
!!$  !Others
!!$  INTEGER m,n
!!$  REAL*8 G11dkes,cmul
!!$  COMPLEX*16 bmn_e(-ntorbd:ntorbd,-mpolbd:mpolbd)
!!$  REAL*8, SAVE :: D11p=1
!!$
!!$  WRITE(iout,*) 'Calculating PLATEAU_OLD'
!!$
!!$  IF(D11p.LT.0) GO TO 999
!!$
!!$  !Change Fourier representation
!!$  bmn_e=0
!!$  DO m=-mpolb,mpolb
!!$     DO n=-ntorb,ntorb
!!$        IF(m.LT.0) bmn_e(n,m)= borbic(-n,-m)/2.
!!$        IF(m.GT.0) bmn_e(n,m)= borbic( n,+m)/2.
!!$        IF(m.EQ.0) bmn_e(n,0)=(borbic( n,0)+borbic(-n,0))/2.
!!$        IF(STELL_ANTISYMMETRIC) THEN
!!$           IF(m.LT.0) bmn_e(n,m)=bmn_e(n,m)-NUMI*borbis(-n,-m)/2.
!!$           IF(m.GT.0) bmn_e(n,m)=bmn_e(n,m)+NUMI*borbis( n,+m)/2.
!!$           IF(m.EQ.0) bmn_e(n,0)=bmn_e(n,0)+NUMI*(borbis( n,0)-borbic(-n,0))/2.
!!$        END IF
!!$     END DO
!!$  END DO
!!$  D11=0 
!!$  DO m=-mpolb,mpolb
!!$     DO n=-ntorb,ntorb
!!$        IF(n.EQ.0.AND.m.EQ.0) CYCLE
!!$        D11=D11-(m*Bzeta-n*nzperiod*Btheta)*(m*Bzeta-n*nzperiod*Btheta)&
!!$             & *REAL(bmn_e(n,m)*bmn_e(-n,-m))/ABS(iota*m+nzperiod*n) 
!!$     END DO
!!$  END DO 
!!$
!!$  D11p=D11
!!$
!!$999 D11=-0.125*pi*(vdconst(jv)*vdconst(jv)/v(jv))*(rad_R/borbic(0,0)/borbic(0,0)/aiBtpBz/aiBtpBz)*D11p 
!!$
!!$  cmul=nu(jv)/v(jv)/2.
!!$  !Connect with Pfirsch-Schlueter or 1/nu
!!$  IF(FACT_CON.GT.0.AND.cmul_PS.GT.0) THEN
!!$     IF(cmul.GT.cmul_PS /FACT_CON) D11=D11*(1+cmul/cmul_PS)
!!$     IF(cmul.LT.cmul_1NU*FACT_CON) D11=D11*(1+cmul_1NU/cmul)
!!$  END IF
!!$
!!$  G11dkes=fdkes(jv)*D11
!!$  WRITE(200+myrank,'(3(1pe13.5)," NaN ",2(1pe13.5),"  NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN")') &
!!$       & nu(jv)/v(jv)/2,Epsi/v(jv)*psip,vmconst(jv)/v(jv),G11dkes,G11dkes
!!$  IF(DEBUG) THEN
!!$     IF(cmul_PS.GT.0) THEN !to be checked
!!$        WRITE(10000+myrank,'("2 ",6(1pe13.5))') nu(jv)/v(jv)/2,Epsi/v(jv)*psip,vmconst(jv)/v(jv),&
!!$             & G11dkes,weight(jv)/fdkes(jv),weight(jv)/fdkes(jv)*v(jv)*v(jv)
!!$     ELSE
!!$        WRITE(10000+myrank,'("0 ",6(1pe13.5))') nu(jv)/v(jv)/2,Epsi/v(jv)*psip,vmconst(jv)/v(jv),&
!!$             & G11dkes,weight(jv)/fdkes(jv),weight(jv)/fdkes(jv)*v(jv)*v(jv)
!!$     END IF
!!$  END IF
!!$
!!$END SUBROUTINE CALC_PLATEAU_OLD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


