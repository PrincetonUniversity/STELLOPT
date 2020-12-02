!DEBUGGING Calculates NTV in a large aspect-ratio circular tokamak using analytical formulas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CALC_NTV(jt,s,ib,regb,Zb,Ab,nb,dnbdpsi,Tb,dTbdpsi,Epsi,L1b,Gb)              

!----------------------------------------------------------------------------------------------- 
!Calculates elliptic integrals of first and second kind
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER jt,ib,regb
  REAL*8 s,Zb,Ab,nb,dnbdpsi,Tb,dTbdpsi,Epsi
  !Output
  REAL*8 L1b,Gb
  !Others
  INTEGER, PARAMETER :: nk2=1000
  !All
  INTEGER n,ik2,ikmin,iv,ivmin
  REAL*8 x1(nv),x2(nv),vmin,detadv(nv),eta1,eta2,k,k2,dk2,eps_s,nu_s,dummy
  REAL*8 ftheta,alphan,betan,eK,eE
!  REAL*8 kx1(nv),kx2(nv),emc(nv),emcx(nv)
  !1/nu
  REAL*8 I1nu,Gb_1nu,L1b_1nu
  !sqrtnu
  REAL*8 nustard,logn,abs,Isqrtnu,Gb_sqrtnu,L1b_sqrtnu
  !sb-p
  REAL*8 k0,k02,Isbp,Gb_sbp,L1b_sbp
  !non-resonant
  INTEGER m
  REAL*8 k1nu(nv),ksqrtnu(nv),ksbp(nv)
  INTEGER, PARAMETER :: npt=128
  INTEGER iz,it
  REAL*8 g(nv,nk2),dgdk2,ming,Epsi_e(nv)
  REAL*8 dzeta,dtheta
  REAL*8 zeta(npt,npt),theta(npt,npt),B(npt,npt)
  REAL*8 zeta_h(npt,npt),theta_h(npt,npt),B_h(npt,npt)
  REAL*8 dG,func,modec(-ntorbd:ntorbd,-mpolbd:mpolbd),modes(-ntorbd:ntorbd,-mpolbd:mpolbd)
!  REAL*8 th,An,Bn,An0,Bn0
  
  jt=jt
  regb=regb
  s=s
  iota=-iota
  
!  DO n=-ntorb,ntorb
!     DO m=0,mpolb
!        WRITE(iout,*) '1',borbic(n,m),borbis(n,m)
!     END DO
  !  END DO


!!!!!  !!!!!  !!!!!  !!!!!  !!!!!  !!!!!  !!!!!  !!!!!
!!$  dtheta=TWOPI/npt
!!$  DO it=1,npt
!!$     th=-PI+(it-0.5)*dtheta     
!!$     DO n=0,ntorb
!!$        CALL CALC_AN_AND_BN(th,n,An,Bn)
!!$        IF(n.EQ.0) THEN
!!$           An0=An
!!$           Bn0=Bn
!!$        END IF
!!$        WRITE(iout,'(1(1pe13.5),I4,1000(1pe13.5))') th,n,-An,Bn,1-An0
!!$        CALL CALC_AN_AND_BN(th,-n,An,Bn)
!!$        WRITE(iout,'(1(1pe13.5),I4,1000(1pe13.5))') th,-n,-An,Bn,1-An0
!!$     END DO
!!$  END DO
!!$  stop
  !!!!!  !!!!!  !!!!!  !!!!!  !!!!!  !!!!!  !!!!!  !!!!!
  
  dzeta =TWOPI/npt/nzperiod
  dtheta=TWOPI/npt
  modec=0
  modes=0
  DO iz=1,npt
     DO it=1,npt
        zeta(iz,it) =(iz-1)*dzeta
        theta(iz,it)=(it-1)*dtheta
        CALL SUM_BORBI(zeta(iz,it),theta(iz,it),B(iz,it))
        func=avB2/B(iz,it)/B(iz,it)
        DO n=-ntorb,ntorb
           DO m=0,mpolb
              modec(n,m)=modec(n,m)+func*COS(m*theta(iz,it)-n*nzperiod*zeta(iz,it))
              modes(n,m)=modes(n,m)+func*SIN(m*theta(iz,it)-n*nzperiod*zeta(iz,it))
           END DO
        END DO
     END DO
  END DO
  modec=modec*dzeta*dtheta/PI/PI
  modes=modes*dzeta*dtheta/PI/PI
  modec(0,0)=modec(0,0)/2.0

  zeta_h =zeta
  theta_h=theta
  DO n=-ntorb,ntorb
     DO m=0,mpolb
        IF(n.EQ.0.AND.m.EQ.0) CYCLE
        DO iz=1,npt
           DO it=1,npt
              dG=(modes(n,m)*COS(m*theta(iz,it)-n*nzperiod*zeta(iz,it))&
                   & -modec(n,m)*SIN(m*theta(iz,it)-n*nzperiod*zeta(iz,it)))/(m*iota-n*nzperiod)
              zeta_h(iz,it) = zeta_h(iz,it)+dG
              theta_h(iz,it)=theta_h(iz,it)+dG*iota
              B_h(iz,it)    =B(iz,it)
           END DO
        END DO
     END DO
  END DO
  
  CALL INTERPOLATE_2DMAP(npt,zeta_h,theta_h,B_h,B)

  borbic=0
  borbis=0
  DO iz=1,npt
     DO it=1,npt
        zeta(iz,it) =(iz-1)*dzeta
        theta(iz,it)=(it-1)*dtheta
        func=B(iz,it)
        DO n=-ntorb,ntorb
           DO m=0,mpolb
              borbic(n,m)=borbic(n,m)+func*COS(m*theta(iz,it)-n*nzperiod*zeta(iz,it))
              borbis(n,m)=borbis(n,m)+func*SIN(m*theta(iz,it)-n*nzperiod*zeta(iz,it))
           END DO
        END DO
     END DO
  END DO
  borbic=borbic*dzeta*dtheta/PI/PI
  borbis=borbis*dzeta*dtheta/PI/PI
  borbic(0,0)=borbic(0,0)/2.0

!    DO n=-ntorb,ntorb
!     DO m=0,mpolb
!        WRITE(iout,*) '2',borbic(n,m),borbis(n,m)
!     END DO
!  END DO
!!!!
  
  dzeta =TWOPI/npt/nzperiod
  dtheta=TWOPI/npt
  modec=0
  modes=0
  DO iz=1,npt
     DO it=1,npt
        zeta(iz,it) =(iz-1)*dzeta
        theta(iz,it)=(it-1)*dtheta
        CALL SUM_BORBI(zeta(iz,it),theta(iz,it),B(iz,it))
        func=B(iz,it)*B(iz,it)/avB2
        DO n=-ntorb,ntorb
           DO m=0,mpolb
              modec(n,m)=modec(n,m)+func*COS(m*theta(iz,it)-n*nzperiod*zeta(iz,it))
              modes(n,m)=modes(n,m)+func*SIN(m*theta(iz,it)-n*nzperiod*zeta(iz,it))
           END DO
        END DO
     END DO
  END DO
  modec=modec*dzeta*dtheta/PI/PI
  modes=modes*dzeta*dtheta/PI/PI
  modec(0,0)=modec(0,0)/2.0

!!$  zeta_h =zeta
!!$  theta_h=theta
!!$  DO n=-ntorb,ntorb
!!$     DO m=0,mpolb
!!$        IF(n.EQ.0.AND.m.EQ.0) CYCLE
!!$        DO iz=1,npt
!!$           DO it=1,npt
!!$              dG=(modes(n,m)*COS(m*theta(iz,it)-n*nzperiod*zeta(iz,it))&
!!$                   & -modec(n,m)*SIN(m*theta(iz,it)-n*nzperiod*zeta(iz,it)))/(m*iota-n*nzperiod)
!!$              zeta_h(iz,it) = zeta_h(iz,it)+dG
!!$              theta_h(iz,it)=theta_h(iz,it)+dG*iota
!!$              B_h(iz,it)    =B(iz,it)
!!$           END DO
!!$        END DO
!!$     END DO
!!$  END DO
!!$
!!$  CALL INTERPOLATE_2DMAP(npt,zeta_h,theta_h,B_h,B)
!!$
!!$  borbic=0
!!$  borbis=0
!!$  DO iz=1,npt
!!$     DO it=1,npt
!!$        zeta(iz,it) =(iz-1)*dzeta
!!$        theta(iz,it)=(it-1)*dtheta
!!$        func=B(iz,it)
!!$        DO n=-ntorb,ntorb
!!$           DO m=0,mpolb
!!$              borbic(n,m)=borbic(n,m)+func*COS(m*theta(iz,it)-n*nzperiod*zeta(iz,it))
!!$              borbis(n,m)=borbis(n,m)+func*SIN(m*theta(iz,it)-n*nzperiod*zeta(iz,it))
!!$           END DO
!!$        END DO
!!$     END DO
!!$  END DO
!!$  borbic=borbic*dzeta*dtheta/PI/PI
!!$  borbis=borbis*dzeta*dtheta/PI/PI
!!$  borbic(0,0)=borbic(0,0)/2.0
!!$
!!$    DO n=-ntorb,ntorb
!!$     DO m=0,mpolb
!!$        WRITE(iout,*) '3',borbic(n,m),borbis(n,m)
!!$     END DO
!!$  END DO
!!$!!!!!
!!$  stop
  
  x1=v/vth(ib)
  x2=x1*x1
  nu_s=nuth(ib)/2.
  eps_s=ABS(borbic(0,1))
  vmin=vth(ib)*SQRT(ABS(2*Epsi*rad_R*psip*Zb/m_e/Ab/vth(ib)/vth(ib)))
  DO iv=1,nv
     IF(v(iv).GT.vmin) EXIT
  END DO
  ivmin=iv
  IF(.NOT.TANG_VM) ivmin=nv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!test
!  IF(Zb.LT.0) THEN
!     kx2=x2*(1.0073/4.4858E-4)
!  ELSE
!     kx2=x2/(1.0073/4.4858E-4)
!  END IF
!  kx1=SQRT(kx2)
!  DO iv=1,nv
!     emc(iv)= ERF( x1(iv))*(1.-0.5/x2(iv))+EXP( -x2(iv))/(SQPI* x1(iv))
!     emcx(iv)=ERF(kx1(iv))*(1.-.5/kx2(iv))+EXP(-kx2(iv))/(SQPI*kx1(iv))
!  END DO
!  detadv=0.5*x2*x2*x1*(x1*x2/(emc+emcx))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !1/nu regime

  detadv=0.5*x2*x2*x1*nuth(ib)/nu
  detadv=detadv*SQPI/2./x1
  CALL INT_VMAX(detadv,weight,eta1)
  detadv=detadv*(x2-2.5)
  CALL INT_VMAX(detadv,weight,eta2)
  dk2=1./nk2  
  I1nu=0
  DO ik2=1,nk2
     k2=(ik2-0.5)*dk2
     k=SQRT(k2)
     ftheta=0
     DO n=-ntorb,ntorb
        CALL CALC_ALPHABETA(k,n,alphan,betan,dummy,dummy)
        ftheta=ftheta+n*n*nzperiod*nzperiod*(alphan*alphan+betan*betan)
     END DO
     CALL ELLIPTIC(k,eK,eE)
     I1nu=I1nu+ftheta/(eE-(1.-k2)*eK)
  END DO
  I1nu=I1nu*dk2  

  L1b_1nu=-(1/4./SQ2/PI/SQPI)*(vth(ib)*vth(ib)*vth(ib)*vth(ib)/nu_s)&
       & *(Ab*m_e/Zb/iota)*(Ab*m_e/Zb/iota)*(eps_s)*SQRT(eps_s)*&
       & I1nu*eta1
  Gb_1nu=-(1/4./SQ2/PI/SQPI)*(vth(ib)*vth(ib)*vth(ib)*vth(ib)/nu_s)&
       & *(Ab*m_e/Zb/iota)*(Ab*m_e/Zb/iota)*(eps_s)*SQRT(eps_s)*&
       & I1nu*(eta1*(dnbdpsi/nb+dTbdpsi/Tb-Zb*Epsi/Tb)+eta2*dTbdpsi/Tb)

  k1nu=I1nu*eps_s/(nu/2.)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !sqrt(nu)

  IF(ABS(Epsi).GT.ALMOST_ZERO.OR.TANG_VM) THEN
  
     detadv=0.5*x2*x2*x1*SQRT(nu/nuth(ib))
     detadv=detadv*SQPI/2./x1
     CALL INT_VMAX(detadv,weight,eta1)
     detadv=detadv*(x2-2.5)
     CALL INT_VMAX(detadv,weight,eta2)
     nustard=ABS(4*nu_s*iota/(eps_s*Epsi))
     !logn=ABS(LOG(16./SQRT(nustard)))
     logn=ABS(LOG(16./SQRT(nustard/(1+nustard))))
     dk2=SQRT(nustard/logn)
     IF(dk2.GT.1) dk2=1e-9 !JL
     Isqrtnu=0
     DO n=-ntorb,ntorb
        k=SQRT(1.-dk2)
        CALL CALC_ALPHABETA(k,n,dummy,dummy,alphan,betan)     
        Isqrtnu=Isqrtnu+SQRT(ABS(1.0*n*nzperiod))*(alphan*alphan+betan*betan)
     END DO
     CALL ELLIPTIC(k,eK,eE)
     Isqrtnu=Isqrtnu/16./eK/eK

     L1b_sqrtnu=(1./4./SQ2/PI/SQPI)*(vth(ib)*vth(ib)*vth(ib)*vth(ib))*(Ab*m_e/Zb/iota)*(Ab*m_e/Zb/iota)* &
          & ABS(iota/Epsi)*SQRT(nustard*eps_s*ABS(logn))*Isqrtnu*eta1
     Gb_sqrtnu=-(1./4./SQ2/PI/SQPI)*(vth(ib)*vth(ib)*vth(ib)*vth(ib))*(Ab*m_e/Zb/iota)*(Ab*m_e/Zb/iota)* &
          & ABS(iota/Epsi)*SQRT(nustard*eps_s*ABS(logn))*Isqrtnu* &
          & (eta1*(dnbdpsi/nb+dTbdpsi/Tb-Zb*Epsi/Tb)+eta2*dTbdpsi/Tb)
     Epsi_e=Epsi
     IF(TANG_VM) Epsi_e=Epsi+Ab*iota*m_e*v*v/(2.*Zb*psip*rad_R)     
     ksqrtnu=Isqrtnu*SQRT(ABS(4.*(nu/2.)*logn/eps_s))*&
          & SQRT(ABS(iota*iota*iota/(Epsi_e*Epsi_e*Epsi_e)))

  ELSE

     L1b_sqrtnu=0.0
     Gb_sqrtnu=0.0
     ksqrtnu=1e2*k1nu
     
  END IF
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !superbanana-plateau for Er=0

!!$  k02=0.827
!!$  k0=SQRT(k02)
!!$  Isbp=0
!!$  DO n=-ntorb,ntorb
!!$     CALL CALC_ALPHABETA(k0,n,dummy,dummy,alphan,betan)
!!$     Isbp=Isbp+ABS(n*nzperiod)*(alphan*alphan+betan*betan)
!!$  END DO
!!$  CALL ELLIPTIC(k0,eK,eE)
!!$  Isbp=Isbp/4./eK
!!$  eta1=(3./4.)*SQPI*k02*(1-k02)
!!$  Gb_sbp=-(1/SQ2/SQPI)*vth(ib)*vth(ib)*(Ab*m_e/ABS(Zb*iota))*SQRT(eps_s)*(rad_R*psip)* &
!!$       & Isbp*eta1*(dnbdpsi/nb+dTbdpsi/Tb-Zb*Epsi/Tb)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !superbanana-plateau

  IF(TANG_VM) THEN
     
     dk2=1./nk2
     DO ik2=1,nk2
        k2=(ik2-0.5)*dk2
        k=SQRT(k2)
        CALL ELLIPTIC(k,eK,eE)
        DO iv=1,nv
           g(iv,ik2)=-2*Epsi*rad_R*psip*Zb/m_e/Ab/v(iv)/v(iv)-(2.*eE/eK-1.0)
        END DO
     END DO
     DO iv=ivmin,nv
        ming=ABS(g(iv,1))
        ikmin=1
        DO ik2=2,nk2
           IF(ABS(g(iv,ik2)).LT.ming) THEN
              ikmin=ik2
              ming=ABS(g(iv,ik2))
           END IF
        END DO
        k02=(ikmin-0.5)*dk2
        k0=SQRT(k2)
        CALL ELLIPTIC(k0,eK,eE)
        IF(ikmin.EQ.1) THEN
           dgdk2=(g(iv,2)-g(iv,1))/dk2
        ELSE IF(ikmin.EQ.nk2) THEN
           dgdk2=(g(iv,nk2)-g(iv,nk2-1))/dk2
        ELSE
           dgdk2=(g(iv,ikmin+1)-g(iv,ikmin-1))/(2*dk2)
        END IF
        Isbp=0
        DO n=-ntorb,ntorb
           CALL CALC_ALPHABETA(k0,n,dummy,dummy,alphan,betan)
           Isbp=Isbp+ABS(n*nzperiod)*(alphan*alphan+betan*betan)
        END DO
        detadv(iv)=Isbp/4/dgdk2
        detadv(iv)=detadv(iv)*rad_R*psip
        ksbp(iv)=8*pi*iota*ABS(Zb)/(Ab*m_e*vth(ib)*vth(ib)*x2(iv))*detadv(iv)     
     END DO
     detadv=detadv*x2*x1
     detadv=detadv*SQPI/2./x1
     detadv(1:ivmin)=0
     CALL INT_VMAX(detadv,weight,eta1)
     detadv=detadv*(x2-2.5)
     CALL INT_VMAX(detadv,weight,eta2)

     L1b_sbp=(1/SQ2/SQPI)*vth(ib)*vth(ib)*(Ab*m_e/ABS(Zb*iota))*SQRT(eps_s)*eta1
     Gb_sbp=-(1/SQ2/SQPI)*vth(ib)*vth(ib)*(Ab*m_e/ABS(Zb*iota))*SQRT(eps_s)* &
          & (eta1*(dnbdpsi/nb+dTbdpsi/Tb-Zb*Epsi/Tb)+eta2*dTbdpsi/Tb)
     
  ELSE

     L1b_sbp=0
     Gb_sbp=0
     
  END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  WRITE(iout,*) 'k',ib,ivmin
  DO iv=1,nv
     WRITE(iout,*) 'k',k1nu(iv),ksqrtnu(iv),ksbp(iv)
  END DO

 !!$  IF(ABS(Epsi).GT.1E-9) THEN
!!$     k1nu=1e2*ksqrtnu
!!$     ksbp=1e2*ksqrtnu
!!$     ivmin=nv
!!$  ELSE
!!$     IF(TANG_VM) THEN
!!$        k1nu=1e2*ksbp
!!$        ksqrtnu=1e2*k1nu       
!!$     ELSE
!!$        ksqrtnu=1e2*k1nu
!!$        ksbp=1e2*k1nu
!!$        ivmin=nv
!!$     END IF
!!$  END IF

  !Connecting formula for 1/nu, sqrtnu and sb-p
  detadv(1:ivmin)=   0.5*   x2(1:ivmin)*x2(1:ivmin)   *x1(1:ivmin)*&
       & ksqrtnu(1:ivmin)*k1nu(1:ivmin)/(ksqrtnu(1:ivmin)+k1nu(1:ivmin))
  IF(ivmin.LT.nv) detadv(ivmin+1:nv)=0.5*x2(ivmin+1:nv)*x2(ivmin+1:nv)*x1(ivmin+1:nv)*&
       ksbp(ivmin+1:nv)*k1nu(ivmin+1:nv)/(ksbp(ivmin+1:nv)+k1nu(ivmin+1:nv))
  detadv=detadv*SQPI/2./x1
  CALL INT_VMAX(detadv,weight,eta1)
  detadv=detadv*(x2-2.5)
  CALL INT_VMAX(detadv,weight,eta2)
  L1b=(1./4./SQ2/PI/SQPI)*SQRT(eps_s)*&
       & (vth(ib)*vth(ib)*vth(ib)*vth(ib))*(Ab*m_e/Zb/iota)*(Ab*m_e/Zb/iota)*eta1
  Gb=-(1./4./SQ2/PI/SQPI)*SQRT(eps_s)*&
       & (vth(ib)*vth(ib)*vth(ib)*vth(ib))*(Ab*m_e/Zb/iota)*(Ab*m_e/Zb/iota)* &
       & (eta1*(dnbdpsi/nb+dTbdpsi/Tb-Zb*Epsi/Tb)+eta2*dTbdpsi/Tb)


  IF(SHAING_1NU) THEN
     L1b=L1b_1nu
     Gb=Gb_1nu
  ELSE IF(SHAING_SQRTNU) THEN
     L1b=L1b_sqrtnu
     Gb=Gb_sqrtnu
  ELSE IF(SHAING_SBP) THEN
     L1b=L1b_sbp
     Gb=Gb_sbp
  END IF

!!$  IF(ABS(Epsi).GT.1E-5) THEN
!!$     IF(.NOT.TANG_VM) Gb=Gb/Gb_sqrtnu*psip
!!$  ELSE
!!$     IF(TANG_VM) THEN
!!$        Gb=Gb/Gb_sbp*psip
!!$     ELSE
!!$        Gb=Gb/Gb_1nu*psip
!!$     END IF
!!$  END IF

  iota=-iota
  
END SUBROUTINE CALC_NTV


!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$
SUBROUTINE CALC_ALPHABETA(k,n,alpha_1nu,beta_1nu,alpha_sbp,beta_sbp)

!----------------------------------------------------------------------------------------------- 
!Calculates alpha_n and beta_n
!----------------------------------------------------------------------------------------------- 
  
  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER n
  REAL*8 k
  !Output
  REAL*8 alpha_1nu,beta_1nu,alpha_sbp,beta_sbp
  !Others
  INTEGER, PARAMETER :: ntheta=10000
  INTEGER itheta
  REAL*8 k2,thetamin,thetamax,theta,dtheta,An,Bn,sine,sqr

  k2=k*k
  thetamax=+2.*ASIN(k)
  thetamin=-thetamax
  alpha_1nu=0
  alpha_sbp=0
  beta_1nu=0
  beta_sbp=0
  dtheta=(thetamax-thetamin)/ntheta
  DO itheta=1,ntheta
     theta=thetamin+(itheta-0.5)*dtheta
     sine=SIN(theta/2.)
     sqr=SQRT(k2-sine*sine)
     CALL CALC_AN_AND_BN(theta,n,An,Bn)
     alpha_1nu=alpha_1nu+An*sqr
     beta_1nu =beta_1nu +Bn*sqr
     alpha_sbp=alpha_sbp+An/sqr
     beta_sbp =beta_sbp +Bn/sqr
  END DO
  alpha_1nu=alpha_1nu*dtheta
  beta_1nu =beta_1nu *dtheta
  alpha_sbp=alpha_sbp*dtheta
  beta_sbp =beta_sbp *dtheta

END SUBROUTINE CALC_ALPHABETA


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALC_AN_AND_BN(theta,n,An,Bn)

!----------------------------------------------------------------------------------------------- 
!Calculates An(theta) and Bn(theta)
!----------------------------------------------------------------------------------------------- 
  
  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER n
  REAL*8 theta
  !Output
  REAL*8 An,Bn
  !Others
  INTEGER m
  REAL*8 sine,cosine

  An=0
  Bn=0
  DO m=0,mpolb
     IF(n.EQ.0.AND.(m.EQ.0.OR.m.EQ.1)) CYCLE
     cosine=COS((m-nzperiod*n/iota)*theta)
     sine  =SIN((m-nzperiod*n/iota)*theta)
     An=An-borbic(n,m)*cosine-borbis(n,m)*sine
     Bn=Bn+borbic(n,m)*sine  -borbis(n,m)*cosine
  END DO
  An=An/borbic(0,0)
  Bn=Bn/borbic(0,0)

END SUBROUTINE CALC_AN_AND_BN


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE ELLIPTIC(k,eK,eE)

!----------------------------------------------------------------------------------------------- 
!Calculates elliptic integrals of first and second kind
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Input
  REAL*8 k,eK,eE
  !Others
  INTEGER, PARAMETER :: ntheta=10000
  INTEGER itheta
  REAL*8 k2,dtheta,theta,sine,sqroot

  eK=0.0
  eE=0.0
  k2=k*k

  dtheta=PI/2./ntheta
  DO itheta=1,ntheta
     theta=(itheta-0.5)*dtheta
     sine=SIN(theta)
     sqroot=SQRT(1.-k2*sine*sine)
     eK=eK+1./sqroot
     eE=eE+   sqroot
  END DO
  eK=eK*dtheta
  eE=eE*dtheta

END SUBROUTINE ELLIPTIC


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE INT_VMAX(Q,weightb,intQ)

!----------------------------------------------------------------------------------------------- 
!Calculates the integral in v intQ of Q, weighed by weightb (convolution with Maxwellian)
!-----------------------------------------------------------------------------------------------
  
  USE GLOBAL
  IMPLICIT NONE
  !Input
  REAL*8 Q(nv),weightb(nv)
  !Output
  REAL*8 intQ
  !Others
  INTEGER iv

  intQ=0
  DO iv=1,nv
     intQ=intQ+weightb(iv)*Q(iv)
  END DO
  
END SUBROUTINE INT_VMAX


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


