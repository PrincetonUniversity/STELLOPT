!Calculate neoclassical transport for trace impurities in the presence of low collisionality
!bulk ions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE TRACE_IMPURITIES(jt,ib,nbb,Zb,Ab,regb,s,nb,dnbdpsi,Tb,dTbdpsi,Epsi,&
     & phi1c,Mbbnm,trMnm,nalphab,zeta,theta,phi1,Mbb,trM,Gb,Db,Vb)

!----------------------------------------------------------------------------------------------- 
!For iteration jt and species ib of nbb of regb, characterized by Zb, Ab, nb, dnbdpsi, Tb, and 
!dTbdpsi at surface s, for radial electric field Epsi and varphi1 at a na x nl grid, calculate
!impurity flux Gb and coefficients Gb, Db, Vb
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Inputs
  INTEGER jt,ib,regb(nbb),nbb,nalphab
  REAL*8 Zb(nbb),Ab(nbb),s,nb(nbb),dnbdpsi(nbb),Tb(nbb),dTbdpsi(nbb),Epsi
  REAL*8 phi1c(Nnmp),Mbbnm(Nnmp),trMnm(Nnmp)
  REAL*8 zeta(nalphab),theta(nalphab),phi1(nalphab,nalphab)
  REAL*8 Mbb(nalphab,nalphab),trM(nalphab,nalphab)
  !Output
  REAL*8 Gb,Db,Vb!,ipf,f_eta
  !Others
  INTEGER, PARAMETER :: nix=5
  INTEGER it,nit!it0
  REAL*8, SAVE :: f_s,f_c=-1
  REAL*8, SAVE, ALLOCATABLE :: absnablapsi2oB2(:,:)
  REAL*8 dBdz(nalphab,nalphab)   ,dBdt(nalphab,nalphab)          , B(nalphab,nalphab)
  REAL*8 dphi1dz(nalphab,nalphab),dphi1dt(nalphab,nalphab),expmZepoT(nalphab,nalphab)
  REAL*8 u0(nalphab,nalphab),u1(nalphab,nalphab),u2(nalphab,nalphab)
  REAL*8 f(nix,nix),Gt(nix),sEpsi,sdnbdpsi(nbb),sdTbdpsi(nbb)
  REAL*8 sMbb(nalphab,nalphab),strM(nalphab,nalphab),Gnct(3)
  REAL*8 Gc,Gnc,DEr,DTi,Dni,Gan,aleph,rhostar_i,rhostar_z,fnorm,nustar_z,nustar_i
  !Time
  CHARACTER*30, PARAMETER :: routine="TRACE_IMPURITIES"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  aleph=eps32*Zb(ib)*Zb(ib)*Ab(2)/Ab(ib)
  rhostar_i=vth(2)*Ab(2)*m_e/Zb(2)/borbic(0,0)/rad_R
  rhostar_z=vth(ib)*Ab(ib)*m_e/Zb(ib)/borbic(0,0)/rad_R
  fnorm=Zb(ib)*rhostar_z*rhostar_z*vth(ib)/eps
  nustar_z=rad_R*nuzi(ib)/aiota/vth(ib)
  nustar_i=rad_R*nuth(2)/aiota/vth(2)
  !Calculate |nablapsi|/B^2 transport for classical transport
  IF(CLASSICAL.AND..NOT.ALLOCATED(absnablapsi2oB2)) THEN
     ALLOCATE(absnablapsi2oB2(nalphab,nalphab))
     CALL CALC_ABSNABLAPSI(nalphab,zeta,theta,absnablapsi2oB2)
  END IF

  IF(CLASSICAL.OR.PLATEAU_OR_PS) expmZepoT=EXP(-Zb(ib)*phi1/Tb(ib))

  IF(PLATEAU_OR_PS) THEN
     !Calculate derivatives of B and varphi_1
     CALL CALC_DERIVATIVES(nalphab,zeta,theta,phi1,B,dphi1dz,dphi1dt,dBdz,dBdt)
     !Solve equations on the flux-surface to determine U,U1 and U2
     CALL CALCULATE_U(nalphab,zeta,theta,B,dBdz,dBdt,Zb(ib)/Tb(ib),&
       & phi1,dphi1dz,dphi1dt,expmZepoT,u0,u1,u2)
     !Calculate f_c and f_s for neoclassical transport
     IF(f_c.LT.0) CALL CALCULATE_F(nalphab,zeta,theta,B,f_c,f_s)          
  END IF

  !Iterations are done in order to calculate D, V, and contributions to V...
  f=0.0
  Gt=0
  IF(D_AND_V.EQ.2) THEN
     nit=2
     f(  1,1)=1   !Gamma_b
     f(2:5,2)=1   !Gamma_a
  ELSE IF(D_AND_V.EQ.3) THEN
     nit=3
     f(1  ,1)=1 !Gamma_b
     f(2:4,2)=1 !Gamma_a proportional to E_r, T' and n_i'
     f(5  ,3)=1 !Gamma_a proportional to n_z'
  ELSE IF(D_AND_V.EQ.5) THEN
     nit=5
     f(1,1)=1 !Gamma_b
     f(2,2)=1 !Gamma_a proportional to E_r
     f(3,3)=1 !Gamma_a proportional to T'
     f(4,4)=1 !Gamma_a proportional to n_i'
     f(5,5)=1 !Gamma_a proportional to n_z'
  END IF

  !Small so that it does not affect Gamma but can be used to calculate D and V
  IF(ABS(dnbdpsi(ib)/nb(ib)).LT.1E-6) dnbdpsi(ib)=1E-5*nb(ib)     
  sEpsi   =Epsi
  sdTbdpsi=dTbdpsi
  sdnbdpsi=dnbdpsi
  sMbb=Mbb
  strM=trM

  DO it=1,nit
     
     Mbb       =         sMbb*f(1,it)
     trM       =         strM*f(1,it)
     Epsi      =        sEpsi*f(2,it)
     dTbdpsi( 2)=sdTbdpsi(2 )*f(3,it)
     dTbdpsi(ib)=sdTbdpsi(ib)*f(3,it)
     dnbdpsi( 2)=sdnbdpsi(2 )*f(4,it)
     dnbdpsi(ib)=sdnbdpsi(ib)*f(5,it)


     !Classical transport should always be included
     IF(CLASSICAL) THEN
        CALL CALC_TRACE_CLASSICAL(it,ib,NBB,ZB,AB,s,nb,dnbdpsi,Tb,dTbdpsi,&
             & Epsi,nalphab,B,absnablapsi2oB2,expmZepoT,Gc)
        Gt(it)=Gc
     END IF

     IF (regb(ib).EQ.10) THEN
        Gnct=0
        CALL CALC_TRACE_PS_O_LOWCOLL(it,ib,NBB,ZB,AB,s,nb,dnbdpsi,Tb,dTbdpsi,Epsi,&
             & nalphab,zeta,theta,B,dBdz,dBdt,dphi1dz,dphi1dt,expmZepoT,&
             & f_c,f_s,u0,u1,u2,Mbb,trM,Gnct(1))
        CALL CALC_TRACE_PLATEAU_O_LOWCOLL(it,ib,NBB,ZB,AB,s,nb,dnbdpsi,Tb,dTbdpsi,Epsi,&
             & nalphab,B,phi1,f_c,f_s,u0,Mbb,trM,Gnct(2))
        IF(it.EQ.1) THEN
!           Gnc=Gnct(MINLOC(ABS(Gnct(1:2)),1))
           Gnc=1./(SUM(1./Gnct(1:2)))
        ELSE
!           Gnc=Gnct(MAXLOC(ABS(Gnct(1:2)),1))
           Gnc=SUM(Gnct(1:2))
        END IF
        
!        CALL CALC_TRACE1NU_O_LOWCOLL(it,ib,NBB,ZB,AB,REGB,s,nb,dnbdpsi,Tb,dTbdpsi,Epsi,Gb)
        ! Gnc=SUM(Gnct(:))
        !Gnct(MAXLOC(ABS(Gnct(:)),1)) !Maxima of the three contributions
 !       
     ELSE IF(regb(ib).EQ.11) THEN 
        CALL CALC_TRACE_PS_O_LOWCOLL(it,ib,NBB,ZB,AB,s,nb,dnbdpsi,Tb,dTbdpsi,Epsi,&
             & nalphab,zeta,theta,B,dBdz,dBdt,dphi1dz,dphi1dt,expmZepoT,&
             & f_c,f_s,u0,u1,u2,Mbb,trM,Gnc)
     ELSE IF(regb(ib).EQ.12) THEN
        CALL CALC_TRACE_PLATEAU_O_LOWCOLL(it,ib,NBB,ZB,AB,s,nb,dnbdpsi,Tb,dTbdpsi,Epsi,&
             & nalphab,B,phi1,f_c,f_s,u0,Mbb,trM,Gnc)
     ELSE IF(regb(ib).EQ.22) THEN
        CALL CALC_TRACE_PLATEAU_O_LOWCOLL_OLD(it,ib,NBB,ZB,AB,s,nb,dnbdpsi,Tb,dTbdpsi,Epsi,&
             & nalphab,B,dBdz,dBdt,dphi1dz,dphi1dt,f_c,f_s,u0,Mbb,trM,Gnc)
     ELSE IF(regb(ib).EQ.13) THEN
        CALL CALC_TRACE1NU_O_LOWCOLL(it,ib,NBB,ZB,AB,REGB,s,nb,dnbdpsi,Tb,dTbdpsi,&
             & Epsi,phi1c,Mbbnm,trMnm,Gnc)
     END IF

     Gt(it)=Gt(it)+Gnc

     Epsi   =sEpsi
     dTbdpsi=sdTbdpsi
     dnbdpsi=sdnbdpsi
     Mbb    =sMbb
     trM    =strM

     IF(DEBUG) WRITE(5800+myrank,'(2I5,6(1pe13.5))') jt,it,Zb(ib),Ab(ib),s,Epsi,Gt(it)/psip

  END DO

  Gb=SUM(Gt(1:nit))
  Gan=Gt(1)
  IF(D_AND_V.GT.2) THEN
     Vb=SUM(Gt(1:nit-1))
     Db=-Gt(nit)/(dnbdpsi(ib)/nb(ib))
     IF(D_AND_V.EQ.5) THEN
        DEr=Gt(2)/(   Epsi    /Tb(ib))
        DTi=-Gt(3)/(dTbdpsi(ib)/Tb(ib))
        Dni=-Gt(4)/(dnbdpsi( 2)/nb( 2))
     END IF
  END IF
  
  IF(ABS(sdnbdpsi(ib)/nb(ib)-1E-5).LT.1E-7) THEN
     sdnbdpsi(ib)=0
     Gb=Gb-Gt(nit)
  END IF
     
  WRITE(900+myrank,'(22(1pe13.5),I3)') s,Zb(ib),Ab(ib),Gb/psip,&
       & Vb/psip,Db/psip/psip,psip*dnbdpsi(ib)/nb(ib),&
       & DEr/psip/psip,      Epsi*psip/Tb(ib),&
       & DTi/psip/psip, psip*dTbdpsi(2)/Tb(2),&
       & Dni/psip/psip, psip*dnbdpsi(2)/nb(2),&
       & Gan/psip,nustar_z,aleph,fnorm,eps32,nustar_i/eps32,&
       & nustar_i/(rhostar_i/eps),nustar_i/rhostar_i,&
       & Zb(ib)*(MAXVAL(phi1)-MINVAL(phi1))/Tb(2)/2.,regb(ib)

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE TRACE_IMPURITIES


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALC_ABSNABLAPSI(nalphab,zeta,theta,absnablapsi2oB2)

!----------------------------------------------------------------------------------------------- 
!For (zeta,theta) square grid of size nalphab, calculate factor absnablapsi2oB2 require for
!calculating classical fluxes
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER nalphab
  REAL*8 zeta(nalphab),theta(nalphab)
  !Output
  REAL*8 absnablapsi2oB2(nalphab,nalphab)
  !Others
  INTEGER iz,it,n,m
  REAL*8 R,dRdtheta,dAzdtheta,dZdtheta,cosine,sine

  DO iz=1,nalphab
     DO it=1,nalphab

        R=0
        dRdtheta= 0
        dAzdtheta=0
        dZdtheta= 0
        
        DO m=0,mpolb
           DO n=-ntorb,ntorb 
              cosine=COS(m*theta(it)+nzperiod*n*zeta(iz))
              sine  =SIN(m*theta(it)+nzperiod*n*zeta(iz))
              R=R                  +rorbic(n,m)*cosine
              dRdtheta =dRdtheta -m*rorbic(n,m)*sine
              dAzdtheta=dAzdtheta-m*porbis(n,m)*cosine
              dZdtheta =dZdtheta +m*zorbis(n,m)*cosine
              IF(STELL_ANTISYMMETRIC) THEN
                 R=R                  +rorbis(n,m)*sine
                 dRdtheta =dRdtheta +m*rorbis(n,m)*cosine
                 dAzdtheta=dAzdtheta+m*porbic(n,m)*sine
                 dZdtheta =dZdtheta -m*zorbic(n,m)*sine
              END IF
           END DO
        END DO
        absnablapsi2oB2(iz,it)=dRdtheta*dRdtheta+&
             & R*R*dAzdtheta*dAzdtheta+&
             & dZdtheta*dZdtheta

        IF(DEBUG) WRITE(5000+myrank,'(4(1pe13.5))') zeta(iz),theta(it),absnablapsi2oB2(iz,it)
        
     END DO
  END DO
  
END SUBROUTINE CALC_ABSNABLAPSI


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALC_DERIVATIVES(nalphab,zeta,theta,phi1,B,dphi1dz,dphi1dt,dBdz,dBdt)
 
!----------------------------------------------------------------------------------------------- 
!For nalphab square grid in (zeta,theta) calculate B and its derivaves, and also derivatives of
!given varphi1
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER nalphab
  REAL*8 zeta(nalphab),theta(nalphab),phi1(nalphab,nalphab)
  !Output
  REAL*8 B(nalphab,nalphab),dBdz(nalphab,nalphab),dBdt(nalphab,nalphab)
  REAL*8 dphi1dz(nalphab,nalphab),dphi1dt(nalphab,nalphab)
  !Others
  INTEGER iz,izp1,izm1,it,itp1,itm1
  REAL*8 dzeta,dtheta,dummy,vdummy(Nnmp)

  !Recalculate magnetic field in grid
  DO iz=1,nalphab
     DO it=1,nalphab
        CALL CALCB(zeta(iz),theta(it),0,.FALSE.,B(iz,it),dummy,dummy,&
             & dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,vdummy)        
     END DO
  END DO

  !Calculate derivatives
  DO iz=1,nalphab
     IF(iz.EQ.nalphab) THEN
        izp1=1
     ELSE
        izp1=iz+1
     END IF
     IF(iz.EQ.1) THEN
        izm1=nalphab
     ELSE
        izm1=iz-1
     END IF
     DO it=1,nalphab
        IF(it.EQ.nalphab) THEN
           itp1=1
        ELSE
           itp1=it+1
        END IF
        IF(it.EQ.1) THEN
           itm1=nalphab
        ELSE
           itm1=it-1
        END IF
        dbdz(iz,it)    =   B(izp1,it)-   B(izm1,it)
        dPhi1dz(iz,it) =phi1(izp1,it)-phi1(izm1,it)
        dbdt(iz,it)   =   B(iz,itp1)-   B(iz,itm1)
        dphi1dt(iz,it)=phi1(iz,itp1)-phi1(iz,itm1)
     END DO
  END DO
  dzeta= 2*TWOPI/nalphab/nzperiod
  dtheta=2*TWOPI/nalphab  
  dbdz   =dbdz   /dzeta
  dbdt   =dbdt   /dtheta
  dPhi1dz=dPhi1dz/dzeta
  dphi1dt=dphi1dt/dtheta

END SUBROUTINE CALC_DERIVATIVES


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALC_TRACE_CLASSICAL(jt,ib,NBB,ZB,AB,s,nb,dnbdpsi,Tb,dTbdpsi,Epsi,&
     & nalphab,B,absnablapsi2oB2,expmZepoT,Gb)

!----------------------------------------------------------------------------------------------- 
!For iteration iz and species ib of nbb, characterized by Zb, Ab, nb, dnbdpsi, Tb, and dTbdpsi
!at surface s, for radial electric field Epsi, for exponential of varphi1 at a nalphab square 
!grid, calculate impurity flux Gb
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Inputs
  INTEGER jt,ib,nbb,nalphab
  REAL*8 Zb(nbb),Ab(nbb),s,nb(nbb),dnbdpsi(nbb),Tb(nbb),dTbdpsi(nbb),Epsi,expmZepoT(nalphab,nalphab)
  REAL*8 B(nalphab,nalphab),absnablapsi2oB2(nalphab,nalphab)
  !Output
  REAL*8 Gb
  !Others
  REAL*8 fsacla,coeffcla,FSA
  !Time
  CHARACTER*30, PARAMETER :: routine="CALC_TRACE_CLASSICAL"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)
  jt=jt
  s=s
  fsacla=FSA(nalphab,nalphab,absnablapsi2oB2*expmZepoT,iBtpBz/B/B,1)
  coeffcla=nuzi(ib)*Ab(ib)*Ab(2)*Tb(ib)*m_e/(Zb(ib)*Zb(2))*fsacla 
  Gb=coeffcla*(dnbdpsi(2)/nb(2)-0.5*dTbdpsi(2)/Tb(2)-(Zb(2)/Zb(ib))*dnbdpsi(ib)/nb(ib))

  IF(DEBUG) WRITE(5500+myrank,'(I5,6(1pe13.5))') jt,Zb(ib),Ab(ib),s,Epsi,Gb/psip
  
  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE CALC_TRACE_CLASSICAL

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALC_TRACE_PS_O_LOWCOLL(jt,ib,nbb,Zb,Ab,s,nb,dnbdpsi,Tb,dTbdpsi,Epsi,&
     & nalphab,zeta,theta,B,dBdz,dBdt,dphi1dz,dphi1dt,expmZepoT,&
     & f_c,f_s,u0,u1,u2,Mbb,trM,Gb)

!----------------------------------------------------------------------------------------------- 
!For iteration jt and species ib of nbb, characterized by Zb, Ab, nb, dnbdpsi, Tb, and dTbdpsi
!at surface s, for radial electric field Epsi, for B, phi1 and their derivatives and exponential,
!and Mbb and trM at a na x nl grid and factors f_c and f_s, calculate PS impurity flux Gb
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Inputs
  INTEGER jt,ib,nbb,nalphab
  REAL*8 Zb(nbb),Ab(nbb),s,nb(nbb),dnbdpsi(nbb),Tb(nbb),dTbdpsi(nbb),Epsi
  REAL*8 zeta(nalphab),theta(nalphab),B(nalphab,nalphab),dBdz(nalphab,nalphab),dBdt(nalphab,nalphab)
  REAL*8 dphi1dz(nalphab,nalphab),dphi1dt(nalphab,nalphab),expmZepoT(nalphab,nalphab)
  REAL*8 f_c,f_s,u0(nalphab,nalphab),u1(nalphab,nalphab),u2(nalphab,nalphab)
  REAL*8 Mbb(nalphab,nalphab),trM(nalphab,nalphab)
  !Output
  REAL*8 Gb
  !Others
  REAL*8 Gf,Ga
  !Time
  CHARACTER*30, PARAMETER :: routine="CALC_TRACE_PS_O_LOWCOLL"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  Gb=0
  !Friction
  IF(FRICTION) THEN
     CALL FRICTION_TRACE_PS(jt,ib,nbb,Zb,Ab,nb,dnbdpsi,Tb,dTbdpsi,Epsi,&
          & f_c,f_s,nalphab,zeta,theta,B,expmZepoT,u0,u1,u2,Gf)
     Gb=Gb+Gf
  END IF
  !Pressure anisotropy
  IF(ANISOTROPY) THEN
     CALL ANISOTROPY_TRACE_PS(jt,Zb(ib),Ab(2),Tb(ib),&
          & nalphab,B,dBdz,dBdt,dphi1dz,dphi1dt,Mbb,trM,Ga)
     Gb=Gb+Ga
  END IF

  IF(DEBUG) WRITE(5400+myrank,'(I5,6(1pe13.5))') jt,Zb(ib),Ab(ib),s,Epsi,Gf/psip,Ga/psip

!  WRITE(iout,1002) jt,Zb(ib),s,Epsi*psip,Gb/psip!,ipf,&
!       & f_eta,Db,Vb,(MAXVAL(phi1)-MINVAL(phi1))/(Tb(2)*2.),&
!       & ipf0,Db0,Vb0
!1002 FORMAT('Gz',I1,30(1pe13.5)) 

!  IF(DEBUG) WRITE(1100+myrank,*) 'GNCPS',Gb,Gf
        
  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)
  
END SUBROUTINE CALC_TRACE_PS_O_LOWCOLL

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALCULATE_F(nalphab,zeta,theta,B,f_c,f_s)

!----------------------------------------------------------------------------------------------- 
!For uniform grid nalphab x nalphab, calculate magnetic field B, factor absnablapsi2oB2 for classical 
!flux and f_c and f_s, see [Helander 2017 JPP]
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER nalphab
  REAL*8 zeta(nalphab),theta(nalphab)
  !Output
  REAL*8 B(nalphab,nalphab),f_c,f_s
  !Others 
  INTEGER, PARAMETER :: nlambda=2056!1028!256
  INTEGER, PARAMETER :: nlx=10000!5000!2560
  INTEGER ilp,jlp,ila,iz,it,frac,numint
  REAL*8 z_l(nlx,1),t_l(nlx,1),B_l(nlx,1),dzeta,dBdz,dBdt,lambda,dlambda,sqrt1mlb(nlambda,nlx),Bmax,FSA
  REAL*8 intg4(nlambda,nlx),g4(nlx,1,nlambda)!,g4zt(nalphab,nalphab,nlambda)
  REAL*8 fact,dummy,vdummy(Nnmp)
  REAL*8 numm,demm,num(nlx/nalphab),denom(nlx/nalphab),FSAg4(nlambda)
  INTEGER nlp
  !Time
  CHARACTER*30, PARAMETER :: routine="CALCULATE_F"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  !Locate maximum of B
  Bmax=0
  DO iz=1,nalphab
     DO it=1,nalphab
        IF(B(iz,it).GT.Bmax) THEN
           Bmax=B(iz,it)
           z_l(1,1)= zeta(iz)
           t_l(1,1)=theta(it)
        END IF
     END DO
  END DO

  intg4=0
  g4=0
  dlambda=1./Bmax/nlambda
  dzeta=zeta(2)-zeta(1)
  DO ilp=1,nlx
     CALL CALCB(z_l(ilp,1),t_l(ilp,1),2,.FALSE.,B_l(ilp,1),dBdz,dBdt,&
          & dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,vdummy)
     DO ila=1,nlambda
        lambda=1/Bmax-ila*dlambda
        sqrt1mlb(ila,ilp)=SQRT(1-lambda*B_l(ilp,1))
        intg4(ila,ilp)=sgnB*(Bzeta*dBdt-Btheta*dBdz)&
             & /sqrt1mlb(ila,ilp)/sqrt1mlb(ila,ilp)/sqrt1mlb(ila,ilp)
        IF(ilp.GT.2.AND.MOD(ilp,2).EQ.1) THEN
           g4(ilp,1,ila)=g4(ilp-2,1,ila)+intg4(ila,ilp-1)*2*dzeta
           g4(ilp-1,1,ila)=0.5*(g4(ilp,1,ila)+g4(ilp-2,1,ila))
        END IF
     END DO
     IF(ilp.EQ.nlx) EXIT
     z_l(ilp+1,1)=z_l(ilp,1)+dzeta
     t_l(ilp+1,1)=t_l(ilp,1)+dzeta*iota
  END DO
  
  DO ila=1,nlambda
     g4(:,1,ila)=g4(:,1,ila)*(sqrt1mlb(ila,:)/2.)*(1/Bmax-ila*dlambda)

     num=0
     denom=0
     DO ilp=1,INT(nlx/nalphab)*nalphab
        jlp=(ilp-1)/nalphab+1
        num(jlp)  =num(jlp)+g4(ilp,1,ila)*aiBtpBz/B_l(ilp,1)/B_l(ilp,1)
        denom(jlp)=denom(jlp)            +aiBtpBz/B_l(ilp,1)/B_l(ilp,1)
!        IF(MOD(ilp,nalphab).EQ.1) alpha(jlp)=MOD(t_l(ilp,1),TWOPI)
     END DO
     nlp=INT(nlx/nalphab)
     numm=0
     demm=0
     DO jlp=1,nlp
        numm=numm+num(jlp)
        demm=demm+denom(jlp)
     END DO
     FSAg4(ila)=numm/demm !this could be improved

  END DO
  
  numint=nlambda
  DO frac=1,nlambda
     f_c=0
     f_s=0
     numint=numint/2
     DO ila=nlambda-numint/2,numint/2,-numint
        lambda=1/Bmax-ila*dlambda
        fact=lambda/FSA(nalphab,nalphab,SQRT(1-lambda*B),aiBtpBz/B/B,1)
        f_c=f_c+fact
        f_s=f_s+fact*FSAg4(ila)
     END DO
     f_c=f_c*(3./4.)*avB2*numint*dlambda
     f_s=f_s*(3./4.)*avB2*numint*dlambda
     IF(DEBUG) WRITE(5100+myrank,'(10(1pe13.5))') numint*dlambda,f_c,f_s,lambda*Bmax
!     print *,'f_s',f_s,(f_c-avB2/Bmax/Bmax)*Bzeta/iota
     IF(numint.EQ.8) EXIT
  END DO
  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)
  
END SUBROUTINE CALCULATE_F


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALCULATE_U(nalphab,zeta,theta,B,dBdz,dBdt,Ze_T,phi1,dphi1dz,dphi1dt,expmZepoT,&
     & u0,u1,u2)
 
!----------------------------------------------------------------------------------------------- 
!For given Ze/T, B and varphi1 defined in a nalphab x nalphab grid, calculate functions u, U1 and
!U2 (related to friction)
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER nalphab
  REAL*8 zeta(nalphab),theta(nalphab),B(nalphab,nalphab),dbdz(nalphab,nalphab)
  REAL*8 dphi1dz(nalphab,nalphab),expmZepoT(nalphab,nalphab)
  REAL*8 Ze_T,phi1(nalphab,nalphab),dbdt(nalphab,nalphab),dphi1dt(nalphab,nalphab)
  !Output
  REAL*8 u0(nalphab,nalphab),u1(nalphab,nalphab),u2(nalphab,nalphab)
  !Others
  INTEGER, SAVE :: i1,i2=0
  INTEGER iz,it,n,m,np1,mp1
  REAL*8, SAVE :: Bmax
  COMPLEX*16 denom
  REAL*8 rhs0(nalphab,nalphab),rhs1(nalphab,nalphab),rhs2(nalphab,nalphab)  
  COMPLEX*16 rhs0nm(nalphab,nalphab),u0nm(nalphab,nalphab)
  COMPLEX*16 rhs1nm(nalphab,nalphab),u1nm(nalphab,nalphab)
  COMPLEX*16 rhs2nm(nalphab,nalphab),u2nm(nalphab,nalphab)
  !Time
  CHARACTER*30, PARAMETER :: routine="CALCULATE_U"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)


  !Construct right-hand side of equation determining U0, U1 and U2
  rhs0=(Bzeta*dbdt -Btheta*  dbdz)*2/aiBtpBz/B
  rhs1=expmZepoT*&
       & (rhs0+(Bzeta*dphi1dt-Btheta*dphi1dz)*Ze_T/aiBtpBz)
  rhs2=Ze_T*phi1*rhs1-&
       & Ze_T*expmZepoT*(Bzeta*dphi1dt-Btheta*dphi1dz)/aiBtpBz
  rhs0=rhs0/B/B
  rhs1=rhs1/B/B
  rhs2=rhs2/B/B

  !Solve by means of Fourier transformation
  CALL FFTF_KN(nalphab,rhs0,rhs0nm)
  CALL FFTF_KN(nalphab,rhs1,rhs1nm)
  CALL FFTF_KN(nalphab,rhs2,rhs2nm)
  u0nm=0
  u1nm=0
  u2nm=0
  DO np1=1,nalphab
     IF(np1.LE.nalphab/2) THEN
        n=np1-1
     ELSE
        n=-nalphab+np1-1
     END IF
     DO mp1=1,nalphab
        IF(mp1.LE.nalphab/2) THEN
           m=mp1-1
        ELSE
           m=-nalphab+mp1-1
        END IF
        denom=NUMI*(n*nzperiod+iota*m)/iBtpBz
        IF(n.EQ.0.AND.m.EQ.0) CYCLE
        u0nm(np1,mp1)=rhs0nm(np1,mp1)/denom
        u1nm(np1,mp1)=rhs1nm(np1,mp1)/denom
        u2nm(np1,mp1)=rhs2nm(np1,mp1)/denom
     END DO
  END DO
  CALL FFTB_KN(nalphab,u0nm,u0)
  CALL FFTB_KN(nalphab,u1nm,u1)
  CALL FFTB_KN(nalphab,u2nm,u2)

  !Substract value at (zeta,theta) where B=Bmax
  IF(i2.EQ.0) THEN
     Bmax=0
     DO iz=1,nalphab
        DO it=1,nalphab
           IF(B(iz,it).GT.Bmax) THEN
              Bmax=B(iz,it)
              i1=iz
              i2=it
           END IF
        END DO
     END DO
  END IF
  u0=u0-u0(i1,i2)
  u1=u1-u1(i1,i2)
  u2=u2-u2(i1,i2)

  IF(DEBUG) THEN
     DO iz=1,nalphab
        DO it=1,nalphab
           WRITE(5200+myrank,'(20(1pe13.5))') &
                & zeta(iz),theta(it),u0(iz,it),u1(iz,it),u2(iz,it),phi1(iz,it),&
                & B(iz,it),Bmax,(1./Bmax/Bmax-1./B(iz,it)/B(iz,it))*(Bzeta/iota)
        END DO
     END DO
  END IF

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE CALCULATE_U


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE FRICTION_TRACE_PS(jt,ib,nbb,Zb,Ab,nb,dnbdpsi,Tb,dTbdpsi,Epsi,&
       & f_c,f_s,nalphab,zeta,theta,B,expmZepoT,u0,u1,u2,Gb)
!----------------------------------------------------------------------------------------------- 
!For iteration it and species ib of nbb, characterized by Zb, Ab, nb, dnbdpsi, Tb, and dTbdpsi
!at surface s, for radial electric field Epsi, for B, exp(-Zephi1/T) and and U, U1 and U2 defined
!in a square nalphab grid, calculate PS impurity flux driven by friction with the main ions
!----------------------------------------------------------------------------------------------- 
  
  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER jt,ib,nbb,nalphab
  REAL*8 Zb(nbb),Ab(nbb),nb(nbb),dnbdpsi(nbb),Tb(nbb),dTbdpsi(nbb),Epsi
  REAL*8 f_c,f_s,zeta(nalphab),theta(nalphab),B(nalphab,nalphab),expmZepoT(nalphab,nalphab)
  REAL*8 u0(nalphab,nalphab),u1(nalphab,nalphab),u2(nalphab,nalphab)
  !Output
  REAL*8 Gb
  !Others
  INTEGER iz,it
  REAL*8 fact1,fact2,A1i,A2i,A1z,A2z
  REAL*8 FSA,B2exp(nalphab,nalphab),fsa1,fsa2,fsa3,fsa4,fsa5,fsa5t,fsa6
  REAL*8 func(nalphab,nalphab),Jac(nalphab,nalphab),A(nalphab,nalphab)

  jt=jt

  A1i=dnbdpsi(2 )/nb(2) +dTbdpsi(2)/Tb(2)-Zb(2)*Epsi/Tb(ib)
  A2i=dTbdpsi(2) /Tb(2)
  A1z=dnbdpsi(ib)/nb(ib)+dTbdpsi(ib)/Tb(ib)-Zb(ib)*Epsi/Tb(ib)
  A2z=dTbdpsi(ib)/Tb(ib)
  A1z=dnbdpsi(ib)/nb(ib)+dTbdpsi(ib)/Tb(ib)-Zb(ib)*Epsi/Tb(ib) 
   
  B2exp=B*B
  Jac=aiBtpBz/B/B

  B2exp=B*B/expmZepoT
  Jac=iBtpBz/B/B

  func=B2exp*(A1z*u1+A2z*u2)*U1
  fsa1=FSA(nalphab,nalphab,func,Jac,1)

  func=B2exp*(A1z*u1+A2z*u2)
  fsa2=FSA(nalphab,nalphab,func,Jac,1)

  func=B2exp*U1
  fsa3=FSA(nalphab,nalphab,func,Jac,1)

  func=B2exp
  fsa4=FSA(nalphab,nalphab,func,Jac,1)

  func=B*B*u0
  fsa5t=FSA(nalphab,nalphab,func,Jac,1)
  func=((A1i-1.5 *A2i)*(f_s+avB2*u0) &
    & +(A1i-1.17*A2i)*f_c*(f_s+fsa5t)/(1-f_c))*B*B/avB2*U1  
  fsa5=FSA(nalphab,nalphab,func,Jac,1)

  func=((A1i-1.5 *A2i)*(f_s+avB2*u0) &
    & +(A1i-1.17*A2i)*f_c*(f_s+fsa5t)/(1-f_c))*B*B/avB2   
  fsa6=FSA(nalphab,nalphab,func,Jac,1)
  

  fact1=Tb(ib)*(Ab(ib)*m_e/Zb(ib)/Zb(ib))*nuzi(ib)
  fact2=fact1*Zb(ib)/Zb(2)

  Gb=-fact1*(fsa1-fsa2*fsa3/fsa4)+fact2*(fsa5-fsa6*fsa3/fsa4)

  IF(DEBUG) THEN     
     !Calculate A for comparison with EUTERPE
     A=psip*(Tb(2)/Zb(2))*(B/avB2)*((A1i-1.5 *A2i)*(f_s+avB2*u0) &
          & +(A1i-1.17*A2i)*f_c*(f_s+fsa5t)/(1-f_c))


     func=u0*u0*B*B
     fsa1=FSA(nalphab,nalphab,func,Jac,1)
     DO iz=1,nalphab
        DO it=1,nalphab
           WRITE(5300+myrank,'(4(1pe13.5))') zeta(iz),theta(it),A(iz,it)
        END DO
     END DO
     
  END IF
!  WRITE(iout,*) 'GZ',Zb(ib),nuzi(ib)*Tb(ib)*Ab(ib)*m_e/Zb(ib)/Zb(2)*(A1i-1.5 *A2i)/avB2/psip,&
!       & FSA(nalphab,nalphab,u0*u0*B*B,Jac,1)-fsa5t*fsa5t/avB2,fsa5t*fsa5t/avB2

END SUBROUTINE FRICTION_TRACE_PS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE ANISOTROPY_TRACE_PS(jt,Zz,Ai,Tz,nalphab,B,dBdz,dBdt,dphi1dz,dphi1dt,Mbb,trM,Gb)

!----------------------------------------------------------------------------------------------- 
!For iteration jt, calculate impurity flux Gb of impurities with Zz, and Tz, in the presence
!of bulk ions of mass Ai and with B, dBdz, dBdt, dphi1dz, dphi1dt, Mbb, trM defined in 
!a square grid of size nalphab
!----------------------------------------------------------------------------------------------- 
  
  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER jt,nalphab
  REAL*8 Zz,Ai,Tz,B(nalphab,nalphab),dBdz(nalphab,nalphab),dBdt(nalphab,nalphab)
  REAL*8 dphi1dz(nalphab,nalphab),dphi1dt(nalphab,nalphab),Mbb(nalphab,nalphab),trM(nalphab,nalphab)
  !Output
  REAL*8 Gb
  !Others
  REAL*8 fact1,fact2,FSA,rhsB(nalphab,nalphab),rhsphi(nalphab,nalphab)

  jt=jt

  fact1=-Ai*m_e/Tz/borbic(0,0)
  fact2= Ai*m_e/Zz/2/borbic(0,0)/borbic(0,0)
  rhsB  =(Btheta*   dbdz-Bzeta   *dbdt)*B/aiBtpBz
  rhsphi=(Btheta*dphi1dz-Bzeta*dphi1dt)*B/aiBtpBz

  Gb=fact1*FSA(nalphab,nalphab,Mbb*rhsphi,aiBtpBz/B/B,1)+&
    &fact2*FSA(nalphab,nalphab,(trM-3*Mbb)*rhsB,aiBtpBz/B/B,1)
       
END SUBROUTINE ANISOTROPY_TRACE_PS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALC_TRACE_PLATEAU_O_LOWCOLL_OLD(jt,ib,nbb,Zb,Ab,s,nb,dnbdpsi,Tb,dTbdpsi,Epsi,&
     & nalphab,B,dBdz,dBdt,dphi1dz,dphi1dt,f_c,f_s,u0,Mbb,trM,Gb)

!----------------------------------------------------------------------------------------------- 
!For iteration jt and species ib of nbb, characterized by Zb, Ab, nb, dnbdpsi, Tb, and dTbdpsi
!at surface s, for radial electric field Epsi, for B, and their derivatives and exponential,
!and Mbb, trM at a na x nl grid, calculate plateau impurity flux Gb
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Inputs
  INTEGER jt,ib,nbb,nalphab
  REAL*8 Zb(nbb),Ab(nbb),s,nb(nbb),dnbdpsi(nbb),Tb(nbb),dTbdpsi(nbb),Epsi
  REAL*8 B(nalphab,nalphab),dBdz(nalphab,nalphab),dBdt(nalphab,nalphab)
  REAL*8 dphi1dz(nalphab,nalphab),dphi1dt(nalphab,nalphab)
  REAL*8 f_c,f_s,u0(nalphab,nalphab),Mbb(nalphab,nalphab),trM(nalphab,nalphab)
  !Output
  REAL*8 Gb
  !Others
  INTEGER n,pn,mn,m,pm,mm
  REAL*8, SAVE :: Gf1B,Gf2B,Gf3B,Gf4B,Ga1B,Ga2B=-1E10
  REAL*8, SAVE :: Gf1E,Gf2E,Gf3E,Gf4E,Ga1E,Ga2E=-1E10
  REAL*8 avB2u0,A0_T,anpim!A1i,A2i,Azt_T(nalphab,nalphab)
  REAL*8 fmu0,fmu1,fmu2,fmu3,factFM,factmu,factB,FSA,sdke1,sdke2,Gf,Ga
  REAL*8 fmu10,fmu21,fmu32
  REAL*8 nSf1(nalphab,nalphab),nSf2(nalphab,nalphab),nSf3(nalphab,nalphab)
  REAL*8 nSf4(nalphab,nalphab),nSa1(nalphab,nalphab),nSa2(nalphab,nalphab)
  COMPLEX*16 nSf1nm(nalphab,nalphab),nSf2nm(nalphab,nalphab),nSf3nm(nalphab,nalphab)
  COMPLEX*16 nSf4nm(nalphab,nalphab),nSa1nm(nalphab,nalphab),nSa2nm(nalphab,nalphab)
  !Time
  CHARACTER*30, PARAMETER :: routine="CALC_TRACE_PLATEAU_O_LOWCOLL"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  Gf=0
  Ga=0

!  A1i=dnbdpsi(2)/nb(2)+dTbdpsi(2)/Tb(2)-Zb(2)*Epsi/Tb(2)
!  A2i=dTbdpsi(2)/Tb(2)
  avB2u0=FSA(nalphab,nalphab,U0*B*B,aiBtpBz/B/B,1)
!  print *,f_s,avB2u0
!  Azt_T=(B/avB2/Zb(2))*((A1i-1.5*A2i)*(f_s+avB2*u0)+&
!       & (A1i-1.17*A2i)*(f_s+avB2u0)*f_c/(1.-f_c))
!  A0_T=FSA(nalphab,nalphab,Azt_T,aiBtpBz/B/B,1)
!  f_s=(f_c-avB2/MAXVAL(B*B))*Bzeta/iota
!  print *,'tok f_s',f_s,(f_c-avB2/MAXVAL(B*B))*Bzeta/iota
!  print *,'tok uB2',avB2u0,(avB2/MAXVAL(B*B)-1)*Bzeta/iota
!  print *,'QS  f_s',f_s,(f_c-avB2/MAXVAL(B*B))*(helM*Bzeta+helN*Btheta)/(helN+iota*helM)/sgnB
!  print *,'QS  uB2',avB2u0,(avB2/MAXVAL(B*B)-1)*(helM*Bzeta+helN*Btheta)/(helN+iota*helM)/sgnB
!  print *,'ratio',f_c,(f_s*borbic(0,0)+avB2*FSA(nalphab,nalphab,B*U0,aiBtpBz/B/B,1))*(1-f_c)/f_c/(f_s+avb2u0)/borbic(0,0)
  A0_T=(dnbdpsi(2)/nb(2)-0.17*dTbdpsi(2)/Tb(2)-Zb(2)*Epsi/Tb(2))*(1./Zb(2)/borbic(0,0))*(f_s+avB2u0)/(1.-f_c)

!  A0qs_T= (dnbdpsi(2)/nb(2)-0.17*dTbdpsi(2)/Tb(2)-Zb(2)*Epsi/Tb(2))*(1./Zb(2)/borbic(0,0))*(-helM*Bzeta/(helN+iota*helM)/sgnB)
!  A0qas_T=(dnbdpsi(2)/nb(2)-0.17*dTbdpsi(2)/Tb(2)-Zb(2)*Epsi/Tb(2))*(1./Zb(2)/borbic(0,0))*(-Bzeta/iota/sgnB)
!  print *,'tok A',jt,A0_T/A0ana_T,A0ana_T/A0qas_T
!  print *,'QS  A',jt,A0_T/A0ana_T,A0ana_T/A0qs_T
!  print *,'sum',FSA(nalphab,nalphab,(B/avB2/Zb(2))*((A1i-1.5*A2i)*(f_s+avB2*u0)),aiBtpBz/B/B,1),&
!       FSA(nalphab,nalphab,(B/avB2/Zb(2))*((A1i-1.17*A2i)*(f_s+avB2u0)*f_c/(1-f_c)),aiBtpBz/B/B,1)
!  A0_T=A0qas_T

  nSf1=m_e*(Bzeta*dBdt   -Btheta*dBdz   )/(aiBtpBz)*avB2/borbic(0,0)/borbic(0,0)    
  nSf2=    (Bzeta*dphi1dt-Btheta*dphi1dz)/(aiBtpBz)*avB2/borbic(0,0)/borbic(0,0)    
!  nSf3=A0_T*m_e*B*(   dBdz+iota   *dBdt)/iBtpBz                  
  nSf3=A0_T*m_e*borbic(0,0)*(   dBdz+iota   *dBdt)/iBtpBz                  
  nSf4=A0_T    *B*(dphi1dz+iota*dphi1dt)/iBtpBz
  nSa1=-m_e*(trM-FSA(nalphab,nalphab,trM,aiBtpBz/B/B,1)) 
  nSa2=+m_e*m_e*borbic(0,0)*((trM-Mbb)-FSA(nalphab,nalphab,trM-Mbb,aiBtpBz/B/B,1))
  CALL FFTF_KN(nalphab,nSf1,nSf1nm)
  CALL FFTF_KN(nalphab,nSf2,nSf2nm)
  CALL FFTF_KN(nalphab,nSf3,nSf3nm)
  CALL FFTF_KN(nalphab,nSf4,nSf4nm)
  CALL FFTF_KN(nalphab,nSa1,nSa1nm)
  CALL FFTF_KN(nalphab,nSa2,nSa2nm)
  
  Gf1B=0
  Gf2B=0
  Gf3B=0
  Gf4B=0
  Ga1B=0
  Ga2B=0
  Gf1E=0
  Gf2E=0
  Gf3E=0
  Gf4E=0
  Ga1E=0
  Ga2E=0
  DO n=-nalphab/2+1,nalphab/2-1

     IF(n.EQ.0) THEN
        pn=1
        mn=1
     ELSE IF(n.GT.0) THEN
        pn=1+n
        mn=nalphab+1-n
     ELSE
        pn=nalphab+1+n
        mn=1-n
     END IF
     DO m=-nalphab/2+1,nalphab/2-1
        IF(n.EQ.0.AND.m.EQ.0) CYCLE  
        IF(m.EQ.0) THEN
           pm=1
           mm=1
        ELSE IF(m.GT.0) THEN
           pm=1+m
           mm=nalphab+1-m
        ELSE
           pm=nalphab+1+m
           mm=1-m
        END IF
        anpim=ABS(n*nzperiod+iota*m)
        Gf1B=Gf1B-REAL(nSf1nm(pn,pm)*nSf1nm(mn,mm))/anpim
        Gf2B=Gf2B-REAL(nSf2nm(pn,pm)*nSf1nm(mn,mm))/anpim
        Gf3B=Gf3B-REAL(nSf3nm(pn,pm)*nSf1nm(mn,mm))/anpim
        Gf4B=Gf4B-REAL(nSf4nm(pn,pm)*nSf1nm(mn,mm))/anpim
        Ga1B=Ga1B-REAL(nSa1nm(pn,pm)*nSf1nm(mn,mm))/anpim
        Ga2B=Ga2B-REAL(nSa2nm(pn,pm)*nSf1nm(mn,mm))/anpim
        Gf1E=Gf1E-REAL(nSf1nm(pn,pm)*nSf2nm(mn,mm))/anpim
        Gf2E=Gf2E-REAL(nSf2nm(pn,pm)*nSf2nm(mn,mm))/anpim
        Gf3E=Gf3E-REAL(nSf3nm(pn,pm)*nSf2nm(mn,mm))/anpim
        Gf4E=Gf4E-REAL(nSf4nm(pn,pm)*nSf2nm(mn,mm))/anpim
        Ga1E=Ga1E-REAL(nSa1nm(pn,pm)*nSf2nm(mn,mm))/anpim
        Ga2E=Ga2E-REAL(nSa2nm(pn,pm)*nSf2nm(mn,mm))/anpim
     END DO
  END DO
  
  factB=PI*TWOPI*aiBtpBz*borbic(0,0)*borbic(0,0)/avB2
  factFM=(Ab(ib)*m_e/Tb(ib)/TWOPI)*SQRT(Ab(ib)*m_e/Tb(ib)/TWOPI)
  factmu=borbic(0,0)*Ab(ib)*m_e/Tb(ib)
  fmu0= 1.*factB*factFM/(factmu)
  fmu1= 1.*factB*factFM/(factmu*factmu)
  fmu2= 2.*factB*factFM/(factmu*factmu*factmu)
  fmu3= 6.*factB*factFM/(factmu*factmu*factmu*factmu)
  
  fmu10=1!./factmu
  fmu21=2!./factmu
  fmu32=3!./factmu

!  print *,fmu0,fmu1,fmu2,fmu3
!  stop
  sdke1=dnbdpsi(ib)/nb(ib)-Zb(ib)*Epsi/Tb(ib)-1.5*dTbdpsi(ib)/Tb(ib)
  sdke2=dTbdpsi(ib)/Tb(ib)
  
  IF(FRICTION) &
  &Gf=Gf1B*fmu2*(sdke1+fmu32*sdke2)*Ab(ib)*Ab(ib)/Zb(ib)/Zb(ib)+Gf1E*fmu1*(sdke1+fmu21*sdke2)*Ab(ib)/Zb(ib)+&
  &  +Gf2B*fmu1*(sdke1+fmu21*sdke2)*Ab(ib)/Zb(ib)              +Gf2E*fmu0*(sdke1+fmu10*sdke2)+&
  & +(Gf3B*fmu21/factmu            *Ab(ib)*Ab(ib)/Zb(ib)       +Gf3E*Ab(ib))*fmu1+&
  & +(Gf4B*fmu10/factmu            *Ab(ib)                     +Gf4E*Zb(ib))*fmu0

!  &Gf=Gf1B*(fmu2*sdke1+fmu3*sdke2)*Ab(ib)*Ab(ib)/Zb(ib)/Zb(ib)+Gf1E*(fmu1*sdke1+fmu2*sdke2)*Ab(ib)/Zb(ib)+&
!  &  +Gf2B*(fmu1*sdke1+fmu2*sdke2)*Ab(ib)/Zb(ib)              +Gf2E*(fmu0*sdke1+fmu1*sdke2)+&
!  &  +Gf3B *fmu2                  *Ab(ib)*Ab(ib)/Zb(ib)       +Gf3E*fmu1*Ab(ib)+&
!  &  +Gf4B *fmu1                  *Ab(ib)                     +Gf4E*fmu0*Zb(ib)

  IF(ANISOTROPY) Ga=(nuzi(ib)*Ab(2)/Tb(ib))*((Ga1B*fmu1*Ab(ib)/Zb(ib)+Ga1E*fmu0)+&
                                           & (Ga2B*fmu2*Ab(ib)/Zb(ib)+Ga2E*fmu1)*Ab(ib)/Tb(ib))

  IF(DEBUG) WRITE(5500+myrank,'(I5,6(1pe13.5))') jt,Zb(ib),Ab(ib),s,Epsi,Gf/psip,Ga/psip

  Gb=Gf+Ga
        
  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)
  
END SUBROUTINE CALC_TRACE_PLATEAU_O_LOWCOLL_OLD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALC_TRACE_PLATEAU_O_LOWCOLL(jt,ib,nbb,Zb,Ab,s,nb,dnbdpsi,Tb,dTbdpsi,Epsi,&
     & nalphab,B,phi1,f_c,f_s,u0,Mbb,trM,Gb)

!----------------------------------------------------------------------------------------------- 
!For iteration jt and species ib of nbb, characterized by Zb, Ab, nb, dnbdpsi, Tb, and dTbdpsi
!at surface s, for radial electric field Epsi, for B, and their derivatives and exponential,
!and Mbb, trM at a na x nl grid, calculate plateau impurity flux Gb
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Inputs
  INTEGER jt,ib,nbb,nalphab
  REAL*8 Zb(nbb),Ab(nbb),s,nb(nbb),dnbdpsi(nbb),Tb(nbb),dTbdpsi(nbb),Epsi
  REAL*8 B(nalphab,nalphab),phi1(nalphab,nalphab)
  REAL*8 f_c,f_s,u0(nalphab,nalphab),Mbb(nalphab,nalphab),trM(nalphab,nalphab)
  !Output
  REAL*8 Gb
  !Others
  INTEGER n,pn,mn,m,pm,mm
  REAL*8 Vp,factf,avB2u0,nBtmmBz,npim,anpim        
  REAL*8 Dnz,Dph,Dti,Dni,Gf,Ga,FSA
  COMPLEX*16 facta
  COMPLEX*16  B1nmoB0(nalphab,nalphab),Zephi1nmoT(nalphab,nalphab)
  COMPLEX*16 Mbbnm(nalphab,nalphab), trMnm(nalphab,nalphab)
  !Time
  CHARACTER*30, PARAMETER :: routine="CALC_TRACE_PLATEAU_O_LOWCOLL2"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  Gf=0
  Ga=0

  CALL FFTF_KN(nalphab,(B   -FSA(nalphab,nalphab,   B,aiBtpBz/B/B,1))/borbic(0,0)  ,B1nmoB0)
  CALL FFTF_KN(nalphab,(phi1-FSA(nalphab,nalphab,phi1,aiBtpBz/B/B,1))*Zb(ib)/Tb(ib),Zephi1nmoT)
  CALL FFTF_KN(nalphab,(Mbb -FSA(nalphab,nalphab, Mbb,aiBtpBz/B/B,1))*Ab(2)/Ab(ib) ,Mbbnm)
  CALL FFTF_KN(nalphab,(trM -FSA(nalphab,nalphab, trM,aiBtpBz/B/B,1))*Ab(2)/Ab(ib) ,trMnm)
  avB2u0=FSA(nalphab,nalphab,U0*B*B,aiBtpBz/B/B,1)
  Vp=TWOPI*TWOPI*aiBtpBz/FSA(nalphab,nalphab,B*B,aiBtpBz/B/B,1)
  factf=2*SQ2*PI*PI*SQPI*Tb(ib)*sqrt(Tb(ib)*Ab(ib)*m_e)/&
       & (Zb(ib)*borbic(0,0)*borbic(0,0)*borbic(0,0)*Vp)
  facta=-NUMI*SQPI*nuzi(ib)*Ab(ib)*m_e*SQRT(Ab(ib)*m_e/Tb(ib))/(SQ2*Zb(ib)*borbic(0,0))
  Dnz=0
  Dph=0
  Dti=0
  Dni=0
  DO n=-nalphab/2+1,nalphab/2-1
     IF(n.EQ.0) THEN
        pn=1
        mn=1
     ELSE IF(n.GT.0) THEN
        pn=1+n
        mn=nalphab+1-n
     ELSE
        pn=nalphab+1+n
        mn=1-n
     END IF
     DO m=-nalphab/2+1,nalphab/2-1
        IF(n.EQ.0.AND.m.EQ.0) CYCLE     
        IF(m.EQ.0) THEN
           pm=1
           mm=1
        ELSE IF(m.GT.0) THEN
           pm=1+m
           mm=nalphab+1-m
        ELSE
           pm=nalphab+1+m
           mm=1-m
        END IF
        nBtmmBz=n*nzperiod*Btheta-m*Bzeta
        npim=n*nzperiod+iota*m
        anpim=ABS(npim)

        IF(FRICTION) THEN
           Dnz=Dnz+(nBtmmBz*nBtmmBz/anpim/Zb(ib))*&
                &  REAL(2*B1nmoB0(pn,pm)*CONJG(B1nmoB0(pn,pm))+&
                &         B1nmoB0(pn,pm)*CONJG(Zephi1nmoT(pn,pm))+&
                &         Zephi1nmoT(pn,pm)*CONJG(B1nmoB0(pn,pm))+&
                &         Zephi1nmoT(pn,pm)*CONJG(Zephi1nmoT(pn,pm)))
           Dph=Dph+(nBtmmBz/anpim)*(nBtmmBz-npim*sgnB*(f_s+avB2u0)/(1-f_c))*&
                &  REAL(2*B1nmoB0(pn,pm)*CONJG(B1nmoB0(pn,pm))+&
                &         B1nmoB0(pn,pm)*CONJG(Zephi1nmoT(pn,pm))+&
                &         Zephi1nmoT(pn,pm)*CONJG(B1nmoB0(pn,pm))+&
                &         Zephi1nmoT(pn,pm)*CONJG(Zephi1nmoT(pn,pm)))
           Dti=Dti+(nBtmmBz*nBtmmBz/anpim/2./Zb(ib))*& 
                &  REAL(6*B1nmoB0(pn,pm)*CONJG(B1nmoB0(pn,pm))+&
                &         B1nmoB0(pn,pm)*CONJG(Zephi1nmoT(pn,pm))+&
                &         Zephi1nmoT(pn,pm)*CONJG(B1nmoB0(pn,pm))-&
                &         Zephi1nmoT(pn,pm)*CONJG(Zephi1nmoT(pn,pm)))+&
                & (0.17*nBtmmBz*sgnB*(f_s+avB2u0)*npim/(anpim*(1.0-f_c)*Zb(2)))*&
                &  REAL(2*B1nmoB0(pn,pm)*CONJG(B1nmoB0(pn,pm))+&
                &         B1nmoB0(pn,pm)*CONJG(Zephi1nmoT(pn,pm))+&
                &         Zephi1nmoT(pn,pm)*CONJG(B1nmoB0(pn,pm))+&
                &         Zephi1nmoT(pn,pm)*CONJG(Zephi1nmoT(pn,pm)))
           
           Dni=Dni-(nBtmmBz*npim*sgnB*(f_s+avB2u0)/(1.0-f_c)/anpim/Zb(2))*&
                &  REAL(2*B1nmoB0(pn,pm)*CONJG(B1nmoB0(pn,pm))+&
                &         B1nmoB0(pn,pm)*CONJG(Zephi1nmoT(pn,pm))+&
                &         Zephi1nmoT(pn,pm)*CONJG(B1nmoB0(pn,pm))+&
                &         Zephi1nmoT(pn,pm)*CONJG(Zephi1nmoT(pn,pm)))

        END IF

        IF(ANISOTROPY) Ga=Ga+REAL(facta*(nBtmmBz/anpim)*&
             &     ((2*B1nmoB0(mn,mm)+Zephi1nmoT(mn,mm))*(trMnm(pn,pm)-Mbbnm(pn,pm))-&
             &        (B1nmoB0(mn,mm)+Zephi1nmoT(mn,mm))* trMnm(pn,pm)))
     END DO
  END DO
  
  IF(FRICTION) Gf=-factf*(Dnz*dnbdpsi(ib)/nb(ib)-Dph*Epsi/Tb(ib)&
       &                 +Dti*dTbdpsi(ib)/Tb(ib)+Dni*dnbdpsi(2)/nb(2))

  IF(DEBUG) WRITE(5500+myrank,'(I5,6(1pe13.5))') jt,Zb(ib),Ab(ib),s,Epsi,Gf/psip,Ga/psip

  Gb=Gf+Ga
        
  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)
  
END SUBROUTINE CALC_TRACE_PLATEAU_O_LOWCOLL


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



SUBROUTINE CALC_TRACE1NU_O_LOWCOLL(jt,ib,NBB,ZB,AB,REGB,s,nb,dnbdpsi,Tb,dTbdpsi,Epsi,&
     & phi1c,Mbbnm,trMnm,Gb)

!----------------------------------------------------------------------------------------------- 
!For iteration jt and species ib of nbb, characterized by Zb, Ab, REGB, nb, dnbdpsi, Tb, and dTbdpsi
!at surface s, for radial electric field Epsi, calculate 1/nu impurity flux
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Inputs
  INTEGER jt,ib,nbb,REGB(NBB)
  REAL*8 Zb(nbb),Ab(nbb),s,nb(nbb),dnbdpsi(nbb),Tb(nbb),dTbdpsi(nbb),Epsi
  REAL*8 phi1c(Nnmp),Mbbnm(Nnmp),trMnm(Nnmp)
  !Output
  REAL*8 Gb
  !Others
  INTEGER iv,nalphab
!  REAL*8 L1b(1,1,nbx,10),L2b(1,1,nbx,10)
!  REAL*8 D11(1,1),Gbt(1,1),Qb(1,1)
!  REAL*8 zeta(nax),theta(nax),mdummy(1,1)
!here
  REAL*8 L1b(Nnmp,Nnmp,nbx,10),L2b(Nnmp,Nnmp,nbx,10)
  REAL*8 D11(Nnmp,Nnmp),Gbt(Nnmp,Nnmp),Qb(Nnmp,Nnmp)!,D31
  REAL*8 zeta(nax),theta(nax),dummy,mdummy(Nnmp,Nnmp)
!here
  !Time
  CHARACTER*30, PARAMETER :: routine="CALC_TRACE1NU_O_LOWCOLL"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  PHI1_READ=.TRUE.
  CALCULATED_INT=.FALSE.
  QN=.FALSE.
  TANG_VM=.FALSE.    !JL
  REMOVE_DIV=.FALSE. !JL
  phi1c=0 !JL

  !Calculate (v,species)-dependent constants
  CALL DKE_CONSTANTS(ib,NBB,ZB,AB,REGB,nb,dnbdpsi,Tb,dTbdpsi,Epsi,.TRUE.)      

  D11=0
  Gbt=0
  Qb=0
  L1b(:,1,ib,jt)=0
  L2b(:,1,ib,jt)=0
  !Scan in v for the calculation of dQ/dv, dGamma/dv, etc
  DO iv=1,nv
     !Perform monoenergetic calculation
     CALL CALC_MONOENERGETIC(ib,ZB(ib),AB(ib),13,.TRUE.,jt,iv,ZERO,&
          & phi1c(1:Nnmp),Mbbnm(1:Nnmp),trMnm(1:Nnmp),&
          & D11(1:Nnmp,1:Nnmp),nalphab,zeta,theta,mdummy(1:Nnmp,1:Nnmp),dummy)
     !Calculate thermal transport coefficients and radial fluxes
     CALL INTEGRATE_V(jt,ib,Ab(ib),.TRUE.,Tb(ib),iv,D11/weight(iv),mdummy(1:Nnmp,1:Nnmp),&
          & L1b(1:Nnmp,1:Nnmp,ib,jt),L2b(1:Nnmp,1:Nnmp,ib,jt),Gbt,Qb,&
          & mdummy(1:Nnmp,1:Nnmp),mdummy(1:Nnmp,1:Nnmp),mdummy(1:Nnmp,1:Nnmp),dummy,dummy,dummy)
     !Check convergence
     IF(DEBUG) WRITE(5600+myrank,'(2I4,30(1pe13.5))') ib,iv,s,Epsi*psip,D11(1,1),L1b(1,1,ib,jt)
     !Check convergence
     IF(iv.GT.iv0.AND.PREC_INTV.GT.0) THEN
        IF(D11(1,1).LT.PREC_DQDV*L1b(1,1,ib,jt)) THEN
           WRITE(iout,'(" Integral in of species #",I1,"&
                & , converged for iv=",I2,", v/v_th=",f7.4)') ib,iv,v(iv)/vth(ib)
           EXIT
        ELSE IF(iv.EQ.nv) THEN
           WRITE(1100+myrank,*) 'Not converged in v'
        END IF
     END IF

  END DO

  Gb=Gbt(1,1)

  WRITE(5700+myrank,'(I5,6(1pe13.5))') jt,Zb(ib),Ab(ib),s,Epsi,Gb/psip
 
  QN=SOLVE_QN.OR.TRIVIAL_QN
  TANG_VM=.TRUE.    !JL
  REMOVE_DIV=.TRUE. !JL



  Gb=-L1b(1,1,ib,jt)*(dnbdpsi(ib)/nb(ib)-Zb(ib)*Epsi/Tb(ib)-1.5*dTbdpsi(ib)/Tb(ib))&
     & -L2b(1,1,ib,jt)*dTbdpsi(ib)/Tb(ib) !check E/v^2

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE CALC_TRACE1NU_O_LOWCOLL
