!Calculate neoclassical transport in the Pfirsch-Schlueter regime  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALC_PS(jv,Epsi,D11)

!----------------------------------------------------------------------------------------------- 
!Calculate monoenergetic transport coefficient D11 in the Pfirsch-Schlueter regime for
!collisionality cmul=nu(jv)/v(jv) and normalized radial electric field efied=Epsi/v(jv)
!----------------------------------------------------------------------------------------------- 
  
  USE GLOBAL
  USE KNOSOS_STELLOPT_MOD
  IMPLICIT NONE
  !Input
  INTEGER jv
  REAL*8 Epsi
  !Output
  REAL*8 D11
  !Others
  INTEGER m,n
  COMPLEX*16 bmn_e(-ntorbd:ntorbd,-mpolbd:mpolbd),vdmn(-ntorbd:ntorbd,-mpolbd:mpolbd)
  REAL*8 omn(-ntorbd:ntorbd,-mpolbd:mpolbd) ,umn(-ntorbd:ntorbd,-mpolbd:mpolbd)
  REAL*8 umn2,fact,G11dkes,cmul

  WRITE(iout,*) 'Calculating PS'

  nu=nu/2. !Different definition than standard
  !Change Fourier representation
  bmn_e=0
  DO m=-mpolb,mpolb
     DO n=-ntorb,ntorb
        IF(m.LT.0) bmn_e(n,m)= borbic(-n,-m)/2.
        IF(m.GT.0) bmn_e(n,m)= borbic(+n,+m)/2.
        IF(m.EQ.0) bmn_e(n,m)=(borbic(+n,+m)+borbic(-n,-m))/2.
        IF(STELL_ANTISYMMETRIC) THEN
           IF(m.LT.0) bmn_e(n,m)=bmn_e(n,m)-NUMI*borbis(-n,-m)/2.
           IF(m.GT.0) bmn_e(n,m)=bmn_e(n,m)+NUMI*borbis(+n,+m)/2.
           IF(m.EQ.0) bmn_e(n,0)=bmn_e(n,0)+NUMI*(borbis(+n,0)-borbic(-n,0))/2.
        END IF
        umn(n,m) =sgnB*v(jv)*(iota*m+n*nzperiod)/rad_R
        omn(n,m) =-Epsi*(m*Bzeta-n*nzperiod*Btheta)/aiBtpBz
       vdmn(n,m)=bmn_e(n,m)*(m*Bzeta-n*nzperiod*Btheta)
     END DO
  END DO

  D11=0
  DO m=-mpolb,mpolb 
     DO n=-ntorb,ntorb
        IF((n.EQ.0.AND.m.EQ.0).OR.ABS(bmn_e(n,m)).LT.ALMOST_ZERO.OR.&
             & (ABS(iota*m+nzperiod*n).LT.5E-2)) CYCLE
        umn2=umn(n,m)*umn(n,m)
        D11=D11-REAL(vdmn(n,m)*vdmn(-n,-m))/(umn2+9.*omn(n,m)*omn(n,m)*nu(jv)*nu(jv)/umn2)
     END DO
  END DO
  fact=vdconst(jv)/(2.*borbic(0,0)*aiBtpBz)
  D11=(16./3.)*D11*(nu(jv))*fact*fact

  nu=nu*2.
  cmul=nu(jv)/v(jv)/2.
  IF(FACT_CON.GT.0.AND.cmul_PS.GT.0.AND.cmul.LT.cmul_PS*FACT_CON) D11=D11+D11pla/fdkes(jv)

  G11dkes=(2*v(jv)/vdconst(jv)/vdconst(jv))/(psip*psip)*D11
  IF(.NOT.KNOSOS_STELLOPT) WRITE(200+myrank,'(3(1pe13.5)," NaN ",2(1pe13.5),&
       & "  NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN")') &
       & nu(jv)/v(jv)/2.,Epsi/v(jv)*psip,vmconst(jv)/v(jv),G11dkes,G11dkes
  IF(DEBUG) THEN
     IF(cmul_PS.GT.0) THEN
        WRITE(10000+myrank,'("1 ",7(1pe13.5))') nu(jv)/v(jv)/2.,Epsi/v(jv)*psip,vmconst(jv)/v(jv),&
             & G11dkes,weight(jv)/fdkes(jv),weight(jv)/fdkes(jv)*v(jv)*v(jv)
     ELSE
        WRITE(10000+myrank,'("0 ",7(1pe13.5))') nu(jv)/v(jv)/2.,Epsi/v(jv)*psip,vmconst(jv)/v(jv),&
             & G11dkes,weight(jv)/fdkes(jv),weight(jv)/fdkes(jv)*v(jv)*v(jv)
     END IF
  END IF

  !Coefficient for transport simulations
  IF(ABS(Epsi).LT.ALMOST_ZERO) etet=1.5*fdkes(jv)*D11*v(jv)/nu(jv) !check
  
END SUBROUTINE CALC_PS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
