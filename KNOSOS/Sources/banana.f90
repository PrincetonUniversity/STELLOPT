!Calculate neoclassical transport in the banana regime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALC_BANANA(jv,Epsi,D11)

!-----------------------------------------------------------------------------------------------
!Calculate monoenergetic transport coefficient in the banana regime for collisionality
!cmul=nu(jv)/v(jv) and normalized radial electric field efied=Epsi/v(jv)
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
!  CHARACTER*100 serr
  INTEGER, PARAMETER :: nal=64
  INTEGER, PARAMETER :: nla=64
  INTEGER, PARAMETER :: nz_l=100
  INTEGER iz,it,ila,iz_l,flag
  REAL*8 G11dkes,cmul
  REAL*8, SAVE :: D11p=1
  REAL*8 dlambda,lambda(nla),Bbounce,Bval,dBdz,dBdt,B1,Bb,B2,Bzt(nal,nal),Jac(nal,nal)
  REAL*8 dz, zeta(nal),z1,zb,z2,z_ini,z_l,z_fin,dz_l
  REAL*8 dt,theta(nal),t1,tb,t2,t_ini,t_l,t_fin,dt_l
  REAL*8 Sov(nal,nal,nla)
  REAL*8 vds(Nnmp),vd(nqv),dummy,vdummy(Nnmp)

  IF(.NOT.USE_B0) THEN
     bnmc0(1:Nnm)=bnmc(1:Nnm)
     bnms0(1:Nnm)=bnms(1:Nnm)
     bnmc1(1:Nnm)=0
     bnms1(1:Nnm)=0
  END IF
  RETURN
!  WRITE(iout,*) 'Calculating BANANA'
!  serr="Not implemented"
!  CALL END_ALL(serr,.FALSE.)

  IF(D11p.LT.0) GO TO 999

  dz=TWOPI/nal/nzperiod
  dt=TWOPI/nal
  DO iz=1,nal
     zeta(iz)=(iz-1)*dz
     theta(iz)=(iz-1)*dt
  END DO
  IF(Bmin_av.LT.ALMOST_ZERO) Bmin_av=1E3
  DO iz=1,nal
     DO it=1,nal
        CALL FILL_BNODE(zeta(iz),theta(it),Jac(iz,it),Bzt(iz,it),vds,.FALSE.)
        IF(Bzt(iz,it).LT.Bmin_av) Bmin_av=Bzt(iz,it)
        IF(Bzt(iz,it).GT.Bmax_av) Bmax_av=Bzt(iz,it)
     END DO
  END DO
  dlambda=1./Bmin_av/nla!-1/Bmax_av
  DO ila=1,nla
     lambda(ila)=(ila-0.5)*dlambda!+1/Bmax_av
  END DO

  Sov=0
  DO iz=1,1
     DO it=1,nal
        CALL EXTREME_POINT(zeta(iz),theta(it),0,z1,t1,B1,dummy,vd,flag)
        IF(flag.EQ.1) THEN  !Depending on flag, the extreme found was a maximum or a minimum
           zb=z1           !If mimimum, go backwards
           tb=t1
           Bb=B1
           CALL EXTREME_POINT(zb,tb,-1,z1,t1,B1,dummy,vd,flag)
        ELSE
           CALL EXTREME_POINT(z1,t1,+1,zb,tb,Bb,dummy,vd,flag)
        END IF
        CALL EXTREME_POINT(zb,tb,+2,z2,t2,B2,dummy,vd,flag)
        z_ini=zeta(iz)
        t_ini=theta(it)
!        IF(B1.GT.B2) THEN
!           CALL BOUNCE_POINT(zeta(iz),theta(it),Bzt(iz,it),z1,t1,z_ini,t_ini,dummy,dummy,vds,-1)
!        ELSE
!           CALL BOUNCE_POINT(zeta(iz),theta(it),Bzt(iz,it),z2,t2,z_ini,t_ini,dummy,dummy,vds,+1)
!        END IF
        DO ila=1,nla
           Bbounce=1./lambda(ila)
           IF(Bbounce.LE.Bzt(iz,it)) CYCLE
           IF(B1.GT.B2) THEN
              IF(Bbounce.GE.B1) THEN
                 z_fin=z1
                 t_fin=t1
              ELSE IF(Bbounce.GT.Bzt(iz,it)) THEN
                 CALL BOUNCE_POINT(zeta(iz),theta(it),Bbounce,z1,t1,z_fin,t_fin,dummy,dummy,vds,-1)
!              ELSE
!                 CALL BOUNCE_POINT(z2,t2,Bbounce,z1,t1,z_fin,t_fin,dummy,dummy,vds,-1)
              END IF
           ELSE
              IF(Bbounce.GE.B2) THEN
                 z_fin=z2
                 t_fin=t2
              ELSE IF(Bbounce.GT.Bzt(iz,it)) THEN
                 CALL BOUNCE_POINT(zeta(iz),theta(it),Bbounce,z2,t2,z_fin,t_fin,dummy,dummy,vds,+1)
 !             ELSE
 !                CALL BOUNCE_POINT(z1,t1,Bbounce,z2,t2,z_fin,t_fin,dummy,dummy,vds,+1)
              END IF
           END IF
           WRITE(8000+myrank,'(10(1pe13.5),3I4)') zeta(iz),theta(it),z_ini,t_ini,z_fin,t_fin,z1,t1,z2,t2,iz,it,ila
           dz_l=(z_fin-z_ini)/nz_l
           dt_l=iota*dz_l
           DO iz_l=1,nz_l
              z_l=z_ini+(iz_l-0.5)*dz_l
              t_l=t_ini+(iz_l-0.5)*dt_l
              CALL CALCB(z_l,t_l,2,.FALSE.,Bval,dBdz,dBdt,&
                   & dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,vdummy)
              Sov(iz,it,ila)=Sov(iz,it,ila)+(dBdz*helM+dBdt*helN)*&
                   & (2-Bval/Bbounce)/(2*Bval*Bval*SQRT(1-Bval/Bbounce))
           END DO
!           WRITE(8000+myrank,'(5(1pe13.5))') zeta(iz),theta(it),lambda(ila),Sov(iz,it,ila)
!           WRITE(3300+myrank,'(I6,1(1pe13.5),2I6)') ipoint,BI3(ipoint)*vdconst(jv)/nu(jv),nal,nlambda
        END DO
     END DO
  END DO
  Sov=Sov*dz_l*iBtpBz/(helN-iota*helM)
  !Sov=Sov*Ab(ib)*m_e/Zb(ib)
  !Sov=Sov*sigma




999 cmul=nu(jv)/v(jv)/2.
  !Connect with Pfirsch-Schlueter or 1/nu
!  IF(FACT_CON.GT.0.AND.cmul_PS.GT.0) THEN
!     IF(cmul.GT.cmul_PS /FACT_CON) D11=D11*(1+cmul/cmul_PS)
!     IF(cmul.LT.cmul_1NU*FACT_CON) D11=D11*(1+cmul_1NU/cmul)
!  END IF

  G11dkes=fdkes(jv)*D11
  IF(.NOT.KNOSOS_STELLOPT) WRITE(200+myrank,'(2(1pe13.5)," NaN NaN ",2(1pe13.5)," &
       & NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN")') &
       & nu(jv)/v(jv)/2.,Epsi/v(jv)*psip,fdkes(jv)*D11,fdkes(jv)*D11
  IF(DEBUG) THEN
     IF(cmul_PS.GT.0) THEN
        WRITE(10000+myrank,'("3 ",5(1pe13.5))') nu(jv)/v(jv),Epsi/v(jv)*psip,G11dkes,&
             & weight(jv)/fdkes(jv),weight(jv)/fdkes(jv)*v(jv)*v(jv)
     ELSE
        WRITE(10000+myrank,'("0 ",5(1pe13.5))') nu(jv)/v(jv),Epsi/v(jv)*psip,G11dkes,&
             & weight(jv)/fdkes(jv),weight(jv)/fdkes(jv)*v(jv)*v(jv)
     END IF
  END IF

  IF(.NOT.USE_B0) THEN
     bnmc0(1:Nnm)=bnmc(1:Nnm)
     bnms0(1:Nnm)=bnms(1:Nnm)
     bnmc1(1:Nnm)=0
     bnms1(1:Nnm)=0
  END IF

END SUBROUTINE CALC_BANANA


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALC_B0

!-----------------------------------------------------------------------------------------------
!Calculate B_0
!-----------------------------------------------------------------------------------------------

  USE GLOBAL
  USE KNOSOS_STELLOPT_MOD
  IMPLICIT NONE
  !Others
  INTEGER, PARAMETER :: npar=1
  INTEGER, PARAMETER :: npt=64
  INTEGER, PARAMETER :: nsurf=npt/2
  LOGICAL shift_twopi
  INTEGER iz,iz0,it,ieta,iturn,ibranch,imat,ipar,n,m,np1,mp1,fint
  INTEGER izmin(npt),hel_N,hel_M
  REAL*8 epsilon,iotat,dist,distb
  REAL*8 Bmax_th(npt),Bmin_th(npt),val_B(nsurf),Bzt(npt,npt),B0(npt,npt),B0zt(npt,npt),Bg(npt)!,Btemp
  REAL*8 Bmax_var,Bmin_var
  REAL*8 Bmax_av2! used for new VBT edi.sanchez@ciemat.es
  REAL*8 val_eta(nsurf),zeta(npt),theta(npt),temp(npt),zetat(npt,npt),zetal(nsurf,npt),zeta0(npt,npt)
  !Matrix
  INTEGER ierr,lwork,rank
  INTEGER, ALLOCATABLE :: iwork(:)
  REAL*8, ALLOCATABLE :: rwork(:),work(:)
  REAL*8 rcond,s_svd(npar),rhs(nsurf*npt),mat(nsurf*npt,npar),parB(npar)
  COMPLEX*16 B0mn(npt,npt)
  REAL(rprec) , SAVE :: save_borbic0(-ntorbd:ntorbd,0:mpolbd)
  LOGICAL, SAVE :: FIRST_TIME=.TRUE.

  WRITE(iout,*) 'Calculating B0'
  borbic0=0!borbic
  borbis0=0!borbis
  dborbic0dpsi=0
  dborbis0dpsi=0
  !Determine helicity of B_0
  Bmax=0
  IF(KN_STELLOPT(6).OR.KN_STELLOPT(7).OR.KN_STELLOPT(8).OR.KN_STELLOPT(9)) THEN
     helN=nzperiod
     helM=0
  ELSE
    DO m=0,mpolbd
      DO n=-ntorbd,ntorbd
        IF(n.EQ.0.AND.m.EQ.0) CYCLE
!        IF(n.EQ.0) CYCLE
        IF(ABS(borbic(n,m)).GT.Bmax) THEN
           Bmax=ABS(borbic(n,m))
           helN=n*nzperiod
           helM=m
        END IF
      END DO
    END DO
  END IF
  iotat=iota/(helN-iota*helM)

  !Look for a QS omnigenous
  IF(QS_B0.OR.helN.EQ.0) THEN
     borbic0(0,0)=borbic(0,0)
     dborbic0dpsi(0,0)=dborbicdpsi(0,0)
     IF(.NOT.QS_B0_1HEL) borbis0(0,0)=borbis(0,0)
     DO fint=1,MIN(ntorbd,mpolbd)
!        IF(QS_B0_1HEL.AND.fint.GT.1) EXIT
        n=helN*fint!/nzperiod
        m=helM*fint
        IF(m.LT.0) THEN
           m=-m
           n=-n
        END IF
        IF(m.GT.mpolbd.OR.ABS(n).GT.ntorbd) CYCLE
        IF(QS_B0_1HEL.AND.fint.GT.1) THEN
           borbic(n,m)=0
           borbis(n,m)=0
        ELSE
           borbic0(n,m)=borbic(n,m)
           dborbic0dpsi(n,m)=dborbicdpsi(n,m)
           IF(.NOT.QS_B0_1HEL) THEN
              borbis0(n,m)=borbis(n,m)
              dborbis0dpsi(n,m)=dborbisdpsi(n,m)
           END IF
        END IF
     END DO
     dborbicdpsi=dborbic0dpsi
     dborbisdpsi=dborbis0dpsi
     RETURN
  END IF

  !Rewrite B in other coordinates
  Bmax=0
  Bmax_th=0
  Bmax_av=0
  Bmax_av2=0!edi.sanchez@ciemat.es
  Bmax_var=0
  Bmin_th=1E3
  Bmin_av=0
  Bmin_var=0
   DO iz=1,npt
     zeta(iz)=(iz-1)*TWOPI/(npt*nzperiod)
  END DO
  DO it=1,npt
     theta(it)=(it-1)*TWOPI/npt
     IF(helN*zeta(1)-helM*theta(1).LT.helN*zeta(2)-helM*theta(1)) THEN
        temp(:)=helN*zeta(:)-helM*theta(it)+10*PI
     ELSE
        temp(:)=helM*theta(it)-helN*zeta(:)+10*PI
     END IF
     DO iz=1,npt
        temp(iz)=MOD(temp(iz),TWOPI)
     END DO
     DO iz=1,npt
        zetat(iz,it)=temp(iz)
        CALL SUM_BORBI(zeta(iz),theta(it),Bzt(iz,it))

     END DO
  END DO

  DO it=1,npt
     DO iz=1,npt
        iz0=MINLOC(zetat(iz:npt,it),1)+iz-1
        temp(1)=zetat(iz,it)
        zetat(iz,it)=zetat(iz0,it)
        zetat(iz0,it)=temp(1)
        temp(1)=Bzt(iz,it)
        Bzt(iz,it)=Bzt(iz0,it)
        Bzt(iz0,it)=temp(1)
     END DO
  END DO

  !Find maximum of B, minimum of B and contours of minima B
  DO it=1,npt
     DO iz=1,npt
        IF(KN_STELLOPT(8).AND.iz.LE.npt/2) KN_WBW=KN_WBW+Bzt(iz,it)*zeta(iz)*zeta(iz)
        IF(Bzt(iz,it).GT.Bmax_th(it)) Bmax_th(it)=Bzt(iz,it)
        IF(Bzt(iz,it).LT.Bmin_th(it)) THEN
           Bmin_th(it)=Bzt(iz,it)
           izmin(it)=iz
        END IF
     END DO
     IF(Bmax_th(it).GT.Bmax) Bmax=Bmax_th(it)
!     Bmax_av=Bmax_av+Bmax_th(it)
      Bmax_av2=Bmax_av2+Bmax_th(it) ! used for new VBT edi.sanchez@ciemat.es
     Bmax_av=Bmax_av+Bzt(1,it)
     Bmin_av=Bmin_av+Bmin_th(it)
     !Bmax_var=Bmax_var+Bzt(1,it)*Bzt(1,it)! the value of B at zeta=0 is used (because for a omnigenous maxima are aligned at zeta=0) edi.sanchez@ciemat.es
     Bmax_var=Bmax_var+Bmax_th(it)*Bmax_th(it)! the maximum value at a given theta is used instead of value at zeta=0 edi.sanchez@ciemat.es
     Bmin_var=Bmin_var+Bmin_th(it)*Bmin_th(it)
  END DO
  Bmax_av=Bmax_av/npt
  Bmax_av2=Bmax_av2/npt! used for new VBT edi.sanchez@ciemat.es
  Bmin_av=Bmin_av/npt
  IF(KN_STELLOPT(8)) KN_WBW=KN_WBW/(borbic(0,0)*TWOPI*TWOPI*npt*npt/2)
  IF(KN_STELLOPT(6)) KN_VBT=(Bmax_var/npt-Bmax_av2*Bmax_av2)/borbic(0,0)/borbic(0,0)
  IF(KN_STELLOPT(7)) KN_VBB=(Bmin_var/npt-Bmin_av*Bmin_av)/borbic(0,0)/borbic(0,0)
  borbic0=borbic
  borbis0=borbis
  dborbic0dpsi=dborbicdpsi
  dborbis0dpsi=dborbisdpsi
  IF(.NOT.KN_STELLOPT(9)) RETURN

  epsilon=((Bmax_av/Bmin_av)-1.)/2.
  shift_twopi=.FALSE.
  IF(ABS(zetat(izmin(npt/2),npt/2)-PI).GT.(0.1*PI)) THEN
     shift_twopi=.TRUE.
     zetat=MOD(zetat+PI,TWOPI)
  END IF

  Bg=0
  DO iz=1,npt
     DO it=1,npt
        Bg(iz)=Bg(iz)+Bzt(iz,it)
     END DO
  END DO
  Bg=Bg/npt

  lwork=-1
  ierr=0
  rcond=-1
  ALLOCATE(rwork(1000),iwork(1000),work(1000))
  DO iturn=1,3
     DO iz=1,npt
        IF(iturn.LE.2) rhs(iz)=Bg(iz)-borbic(0,0)
        DO ipar=1,npar
           mat(iz,ipar)=COS(ipar*zetat(iz,1))
        END DO
     END DO
     IF(iturn.LE.2) THEN
        IF(iturn.EQ.2) lwork=MIN(1000,INT(work(1)))
        CALL DGELSD(npt,npar,1,mat(1:npt,:),npt,rhs(1:npt),npt,&
             & s_svd,rcond,rank,work,lwork,rwork,iwork,ierr)
     END IF
  END DO
  parB=rhs(1:npar)

  !Find target contour lines
  DO ieta=1,nsurf
     val_eta(ieta)=ieta*PI/nsurf
     val_B(ieta)=borbic(0,0)
     DO ipar=1,npar
        val_B(ieta)=val_B(ieta)+parB(ipar)*COS(ipar*val_eta(ieta))
     END DO
     DO it=1,npt
        IF(ieta.EQ.nsurf.OR.val_B(ieta).LT.Bmin_th(it)) THEN
           zetal(ieta,it)=zetat(izmin(it),it)
        ELSE IF(val_B(ieta).GT.Bmax_th(it)) THEN
           zetal(ieta,it)=zeta(1)
        ELSE
           dist=1E3
           DO iz=1,npt
              distB=ABS(Bzt(iz,it)-val_B(ieta))
              IF(distB.LT.dist.AND.zetat(iz,it).LT.zetat(izmin(it),it)) THEN
                 zetal(ieta,it)=zetat(iz,it)
                 dist=distB
              END IF
           END DO
        END IF
     END DO
  END DO
!  DO it=1,npt
!     DO ieta=nsurf,1,-4
!        WRITE(6,'(5(1pe13.5),I3)') zetal(ieta,it),theta(it),val_b(ieta)
!     END DO
!  END DO

  !Find parameters
!  dist=0 !new12
  lwork=-1
  ierr=0
  rcond=-1
  DEALLOCATE(rwork,iwork,work)
  ALLOCATE(rwork(1000),iwork(1000),work(1000))
  DO iturn=1,3
     DO ieta=1,nsurf
        DO it=1,npt
           imat=(ieta-1)*npt+it
           IF(iturn.LE.2) rhs(imat)=zetal(ieta,it)-val_eta(ieta)
           DO ipar=1,npar
              mat(imat,ipar)=val_eta(ieta)*SIN(ipar*(theta(it)+iotat*(PI-val_eta(ieta))))
           END DO
        END DO
     END DO
     IF(iturn.LE.2) THEN
        IF(iturn.EQ.2) lwork=MIN(1000,INT(work(1)))
        CALL DGELSD(nsurf*npt,npar,1,mat,nsurf*npt,rhs,nsurf*npt,&
             & s_svd,rcond,rank,work,lwork,rwork,iwork,ierr)
     ELSE
        !Write B0
!        OPEN(unit=1,file="Bomni.dat",form='formatted',action='write')
        DO ieta=1,nsurf
           DO ibranch=1,-1,-2
              val_eta(ieta)=-ibranch*(PI-ieta*PI/nsurf)+PI
              val_B(ieta)=borbic(0,0)
              DO ipar=1,npar
                 val_B(ieta)=val_B(ieta)+parB(ipar)*COS(ipar*val_eta(ieta))
              END DO
              IF(ibranch.EQ.1) THEN
                 iz=ieta+1
              ELSE
                 iz=npt-ieta+1
              END IF
              B0(iz,:)=val_B(ieta)
              DO it=1,npt
                 IF(ibranch.EQ.-1.AND.ieta.EQ.nsurf) THEN
                    iz=1
                    zeta0(iz,it)=TWOPI
                    B0(iz,it)=Bmax_av
                 ELSE
                    zeta0(iz,it)=val_eta(ieta)
                    DO ipar=1,npar
                       zeta0(iz,it)=zeta0(iz,it)+rhs(ipar)*&
                            & (val_eta(ieta)+PI*(ibranch-1))*&
                            & SIN(ibranch*ipar*(theta(it)+iotat*(PI-val_eta(ieta))))
                    END DO
                 END IF
                 IF(ibranch.EQ.1) WRITE(iout,'(3(1pe13.5),I3)') MOD(zeta0(iz,it),TWOPI),zetal(ieta,it),theta(it),iz
                 IF(shift_twopi) zeta0(iz,it)=MOD(zeta0(iz,it)+PI,TWOPI)
                 zeta0(iz,it)=MOD((zeta0(iz,it)+helM*theta(it))/helN+TWOPI,TWOPI/nzperiod)
                 IF(helM.NE.0.AND.Bzt(1,1).LT.Bzt(1,npt/2)) zeta0(iz,it)=MOD(zeta0(iz,it)+PI,TWOPI/nzperiod)
                 !dist=dist+(zeta0(iz,it)-val_eta(ieta))*(zeta0(iz,it)-val_eta(ieta)) !new2

!new1
!                 CALL SUM_BORBI(zeta0(iz,it),theta(it),Btemp)
!                 dist=dist+(Btemp-val_B(ieta))*(Btemp-val_B(ieta))
!new1
!                 WRITE(1,'(3(1pe13.5),I3)') zeta0(iz,it),theta(it),B0(iz,it),-iz
              END DO
           END DO
        END DO
!        CLOSE(1)
     END IF
  END DO

!  KN_DBO=SQRT(dist/(npt*(nsurf+1)))/TWOPI !new2
!  KN_DBO=SQRT(dist/(npt*npt))/borbic(0,0) !new1
!  RETURN !new12

  !Interpolate
  DO it=1,npt
     DO iz=1,npt
        iz0=MINLOC(zeta0(iz:npt,it),1)+iz-1
        temp(1)=zeta0(iz,it)
        zeta0(iz,it)=zeta0(iz0,it)
        zeta0(iz0,it)=temp(1)
        temp(1)=B0(iz,it)
        B0(iz,it)=B0(iz0,it)
        B0(iz0,it)=temp(1)
        temp(1)=Bzt(iz,it)
        Bzt(iz,it)=Bzt(iz0,it)
        Bzt(iz0,it)=temp(1)
     END DO
  END DO
!  OPEN(unit=1,file="Bomni2.dat",form='formatted',action='write')
  DO iz=1,npt
     DO it=1,npt
        CALL LAGRANGE(zeta0(:,it),B0(:,it),npt,zeta(iz),B0zt(iz,it),2)
        CALL SUM_BORBI(zeta(iz),theta(it),Bzt(iz,it))
!        WRITE(1,'(4(1pe13.5))') zeta(iz),theta(it),B0zt(iz,it),Bzt(iz,it)
     END DO
  END DO
!  CLOSE(1)

  !Determine Fourier modes in DKES style
  CALL FFTF_KN(npt,B0zt,B0mn)
  plan_fwd=0
  borbic0=0
  borbis0=0
  DO np1=1,npt
     IF(np1.LE.npt/2) THEN
        n=np1-1
     ELSE
        n=-npt+np1-1
     END IF
     DO mp1=1,npt
        IF(mp1.LE.npt/2) THEN
           m=mp1-1
        ELSE
           m=-npt+mp1-1
        END IF
        IF(n.EQ.0.AND.m.EQ.0) THEN
           borbic0(0,0)=REAL(B0mn(np1,mp1))
        ELSE IF(m.EQ.0) THEN
           borbic0(ABS(n),0)=borbic0(ABS(n),0)+REAL(B0mn(np1,mp1))
        ELSE IF(n.EQ.0) THEN
           borbic0(0,ABS(m))=borbic0(0,ABS(m))+REAL(B0mn(np1,mp1))
        ELSE
           borbic0(-n*m/ABS(m),ABS(m))=borbic0(-n*m/ABS(m),ABS(m))+REAL(B0mn(np1,mp1))
        END IF
     END DO
  END DO

  !Reescale B_0
  borbic0(0,0)=borbic(0,0)
  IF(FIRST_TIME) save_borbic0=borbic0
  FIRST_TIME=.FALSE.
  borbic0=save_borbic0

  !Calculate distance between B and B_0
  dist=0
  distB=0
  DO n=-ntorbd,ntorbd
     DO m=0,mpolbd
        IF(ABS(borbic0(n,m)).LT.3E-5) CYCLE
        dist=dist+(borbic(n,m)-borbic0(n,m))*(borbic(n,m)-borbic0(n,m))
        IF(ABS(borbic(n,m)-borbic0(n,m)).GT.distB) THEN
           distB=ABS(borbic(n,m)-borbic0(n,m))
           hel_N=n
           hel_M=m
        END IF
     END DO
  END DO
  dist=SQRT(dist)/borbic(0,0)
  IF(KN_STELLOPT(9)) KN_DBO=dist
  !  dist=SQRT(dist)/ABS(borbic(helN,helM))

  borbic0=borbic
  borbis0=borbis
  dborbic0dpsi=dborbicdpsi
  dborbis0dpsi=dborbisdpsi

  IF ( ALLOCATED(rwork) ) DEALLOCATE(rwork,iwork,work)

END SUBROUTINE CALC_B0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE SUM_BORBI(zeta,theta,B)

!-----------------------------------------------------------------------------------------------
!Calculate magnetic field B at (zeta,theta) using borbic and borbis
!-----------------------------------------------------------------------------------------------

  USE GLOBAL
  IMPLICIT NONE
  !Input
  REAL*8 zeta,theta
  !Output
  REAL*8 B
  !Others
  INTEGER n,m
  REAL*8 arg

  B=0
  DO n=-ntorbd,ntorbd
     DO m=0,mpolbd
        arg=m*theta-n*nzperiod*zeta
        B=B+borbic(n,m)*COS(arg)
        B=B+borbis(n,m)*SIN(arg)
     END DO
  END DO

END SUBROUTINE SUM_BORBI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

