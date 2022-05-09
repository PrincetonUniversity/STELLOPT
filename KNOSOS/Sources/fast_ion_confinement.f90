
!Calculate neoclassical transport of fast ions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALC_FAST_ION_CONFINEMENT(s,is,ns,nal,nlambda)
  
!--------------------------------------------------------------------------------------------- 
!Calculate neoclassical transport of fast ions
!The DKE is solved in a nalxnlambda grid in (alpha,lambda)
!-----------------------------------------------------------------------------------------------

  USE GLOBAL
  USE KNOSOS_STELLOPT_MOD
#ifdef IPPorNIFS
  USE petscsys
  USE petscksp
#endif
  IMPLICIT NONE
!  Input
  INTEGER is,ns,nlambda,nal
  REAL*8 s(ns)
  !Others
  INTEGER, SAVE :: nalpha,nalphab,nalphab_save,npoint
  REAL*8, SAVE :: offset,theta_save(nax)
  REAL*8, SAVE, ALLOCATABLE :: lambda(:)
  !Wells and bounce points
  LOGICAL, ALLOCATABLE :: connected(:,:),bottom(:),ltemp(:),ltemp2(:,:)
  INTEGER nw,na
  REAL*8, ALLOCATABLE :: z1(:),t1(:),B1(:),hBpp1(:),vd1(:,:)
  REAL*8, ALLOCATABLE :: zb(:),tb(:),Bb(:),hBppb(:),vdb(:,:) 
  REAL*8, ALLOCATABLE :: z2(:),t2(:),B2(:),hBpp2(:),vd2(:,:) 
  REAL*8, ALLOCATABLE :: alphap_w(:),Bt(:),Btt(:),temp(:),temp2(:,:)
  !Angular and lambda grid
  INTEGER, ALLOCATABLE :: i_w(:),itemp(:)
  INTEGER, SAVE, ALLOCATABLE :: i_l(:),i_p(:,:,:),j_al(:,:)
  REAL*8, SAVE, ALLOCATABLE :: zetap(:),thetap(:,:),zetax(:,:),thetax(:,:),B_al(:,:),vds_al(:,:,:)
  REAL*8, ALLOCATABLE :: one_o_lambda(:),alphap(:),lambdab_w(:),lambdac_w(:)
  REAL*8 zeta(nax),theta(nax),dlambdap,dalphap(nturn+1)
  REAL*8, SAVE :: lambdac
   !Alpha neighbours
  REAL*8, ALLOCATABLE, SAVE :: zlw(:),zrw(:)
  REAL*8, ALLOCATABLE       :: tlw(:),trw(:)
  INTEGER ila,jla,il,ia
  REAL*8, SAVE, ALLOCATABLE :: BI1(:),BI2(:),BI3(:),BI4(:),BI5(:),BI6(:),BI7(:),BI8(:,:)
  INTEGER, ALLOCATABLE :: ia_out(:)
  REAL*8, SAVE, ALLOCATABLE :: tau(:)
  !Time
  CHARACTER*30, PARAMETER :: routine="CALC_FAST_ION_CONFINEMENT"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  ns=ns
  s=s
  WRITE(iout,*) 'Calculating fast ion confinement'

!  xxx xxx 

 !Find and characterize wells
  ALLOCATE(connected(nwx,nwx),bottom(nwx),&
       & z1(nwx),t1(nwx),B1(nwx),hBpp1(nwx),vd1(nqv,nwx),&
       & zb(nwx),tb(nwx),Bb(nwx),hBppb(nwx),vdb(nqv,nwx),&
       & z2(nwx),t2(nwx),B2(nwx),hBpp2(nwx),vd2(nqv,nwx),& 
       & alphap_w(nwx),Bt(nwx),Btt(nwx),&
       & lambdab_w(nwx),lambdac_w(nwx))
  CALL CHARACTERIZE_WELLS(nal,na,nalpha,nw,z1,t1,B1,hBpp1,vd1, &
       & zb,tb,Bb,hBppb,vdb, &
       & z2,t2,B2,hBpp2,vd2, &
       & Bt,Btt,alphap_w,dalphap(1),bottom,connected,offset)

  !Resize arrays (nwx->nw)
  ALLOCATE(temp(nw),temp2(nqv,nw),ltemp(nw))
  temp=z1(1:nw);DEALLOCATE(z1);ALLOCATE(z1(nw));z1=temp
  temp=zb(1:nw);DEALLOCATE(zb);ALLOCATE(zb(nw));zb=temp
  temp=z2(1:nw);DEALLOCATE(z2);ALLOCATE(z2(nw));z2=temp
  temp=t1(1:nw);DEALLOCATE(t1);ALLOCATE(t1(nw));t1=temp
  temp=tb(1:nw);DEALLOCATE(tb);ALLOCATE(tb(nw));tb=temp
  temp=t2(1:nw);DEALLOCATE(t2);ALLOCATE(t2(nw));t2=temp
  temp=B1(1:nw);DEALLOCATE(B1);ALLOCATE(B1(nw));B1=temp
  temp=Bb(1:nw);DEALLOCATE(Bb);ALLOCATE(Bb(nw));Bb=temp
  temp=B2(1:nw);DEALLOCATE(B2);ALLOCATE(B2(nw));B2=temp
  temp=hBpp1(1:nw);DEALLOCATE(hBpp1);ALLOCATE(hBpp1(nw));hBpp1=temp
  temp=hBppb(1:nw);DEALLOCATE(hBppb);ALLOCATE(hBppb(nw));hBppb=temp
  temp=hBpp2(1:nw);DEALLOCATE(hBpp2);ALLOCATE(hBpp2(nw));hBpp2=temp
  temp=alphap_w(1:nw);DEALLOCATE(alphap_w);ALLOCATE(alphap_w(nw));alphap_w=temp
  temp=lambdab_w(1:nw);DEALLOCATE(lambdab_w);ALLOCATE(lambdab_w(nw));lambdab_w=temp
  temp=lambdac_w(1:nw);DEALLOCATE(lambdac_w);ALLOCATE(lambdac_w(nw));lambdac_w=temp
  temp= Bt(1:nw);DEALLOCATE(Bt) ;ALLOCATE(Bt(nw)) ;Bt=temp
  temp=Btt(1:nw);DEALLOCATE(Btt);ALLOCATE(Btt(nw));Btt=temp
  ltemp=bottom(1:nw);DEALLOCATE(bottom);ALLOCATE(bottom(nw));bottom=ltemp
  temp2=vd1(1:nqv,1:nw);DEALLOCATE(vd1);ALLOCATE(vd1(nqv,nw));vd1=temp2
  temp2=vdb(1:nqv,1:nw);DEALLOCATE(vdb);ALLOCATE(vdb(nqv,nw));vdb=temp2
  temp2=vd2(1:nqv,1:nw);DEALLOCATE(vd2);ALLOCATE(vd2(nqv,nw));vd2=temp2
  DEALLOCATE(temp,temp2,ltemp)  
  ALLOCATE(ltemp2(nw,nw))
  ltemp2=connected(1:nw,1:nw);DEALLOCATE(connected);ALLOCATE(connected(nw,nw));connected=ltemp2
  DEALLOCATE(ltemp2)

  !Create grid in alpha, then (zeta,theta)
  !Determine number of modes
  nalphab=1
  DO WHILE(nalphab.LT.nalpha*1.2)
     nalphab=nalphab*2
  END DO
  nalphab=nalphab/2
  nalphab_save=nalphab
  IF(ALLOCATED(zetap)) DEALLOCATE(zetap,thetap,zetax,thetax,B_al,vds_al,j_al)
  ALLOCATE(zetax(nalpha,nalphab),thetax(nalpha,nalphab),alphap(nalpha),&
         & zetap(nalphab),thetap(nalpha,nalphab),&
         & B_al(nalpha,nalphab),vds_al(Nnmp,nalpha,nalphab),j_al(nalpha,nalphab))
  CALL CREATE_ANGULAR_GRID(na,nalpha,nalphab,alphap,dalphap,offset,&
       & zetap,thetap,zetax,thetax,B_al,vds_al,j_al)
  zeta(1:nalphab) =zetap  !square grid
  theta(1:nalphab)=zetap*nzperiod 
  theta_save=theta

  CALL EXCLUDE_WELLS(na,nalpha,nalphab,nw,bottom,connected,&
       & alphap_w,z1,zb,z2,Bb,Bt,zetax,thetax)
  !Set global grid in lambda
  IF(ALLOCATED(lambda)) DEALLOCATE(lambda,i_p)
  ALLOCATE(lambda(nlambdax),one_o_lambda(nlambdax))
  !Set global grid in lambda
  CALL CREATE_LAMBDA_GRID(nlambda,nw,Bb,Bt,&
       & lambdab_w,lambdac_w,lambdac,dlambdap,lambda,one_o_lambda)
  !Resize arrays (nlambdax->nlambda)
  ALLOCATE(temp(nlambda))
  temp=lambda(1:nlambda);      DEALLOCATE(lambda);      ALLOCATE(lambda(nlambda));      lambda=temp
  temp=one_o_lambda(1:nlambda);DEALLOCATE(one_o_lambda);ALLOCATE(one_o_lambda(nlambda));one_o_lambda=temp
  DEALLOCATE(temp)
  ALLOCATE(i_p(nlambda,nalpha,nalphab))
  IF(LJMAP.GT.0) THEN
     jla=1
     DO ila=2,nlambda
        IF(ABS(lambda(ila)-LJMAP).LT.ABS(lambda(jla)-LJMAP)) jla=ila
     END DO
     lambda(jla)=LJMAP
  END IF
  !For each point in the (zeta,theta) grid, determine well and absolute point
  !For each absolute point, determine alpha, lambda and well number
  IF(ALLOCATED(i_l)) DEALLOCATE(i_l)
  ALLOCATE(i_l(npointx),i_w(npointx))
  CALL LABEL_GRIDPOINTS(nalpha,nalphab,nlambda,nw,bottom,connected,&
       & alphap_w,z1,z2,Bb,Bt,lambda,&
       & zetap,zetax,thetax,B_al,npoint,i_l,i_w,i_p)
  ALLOCATE(itemp(npoint))
  itemp=i_l(1:npoint);DEALLOCATE(i_l);ALLOCATE(i_l(npoint));i_l=itemp
  itemp=i_w(1:npoint);DEALLOCATE(i_w);ALLOCATE(i_w(npoint));i_w=itemp
  DEALLOCATE(itemp)

  IF(ALLOCATED(BI1)) THEN
     DEALLOCATE(BI1,BI2,BI3,BI4,BI5,BI6,BI7,BI8,zlw,zrw,tau)
  END IF
  ALLOCATE(BI1(npoint),BI2(npoint),BI3(npoint))
  ALLOCATE(BI4(npoint),BI5(npoint),BI6(npoint),BI7(npoint),BI8(npoint,Nnmp))
  ALLOCATE(zlw(npoint),zrw(npoint),tlw(npoint),trw(npoint)) 
  ALLOCATE(ia_out(npoint),tau(npoint))
  
  !Order alphas in interval [0,2*pi]
  CALL SORT_ALPHA(nalpha,nalphab,alphap,zetax,thetax,thetap,B_al,vds_al,j_al,nlambda,i_p)
  DO WHILE(MAXVAL(thetap(:,1))-MINVAL(thetap(:,1)).GT.TWOPI)
     DO ia=1,nalpha
        IF(ABS(thetap(ia,1)-thetap(1,1)).GT.TWOPI) THEN
           thetap(ia,:)=thetap(ia,:)-SIOTA*TWOPI
           alphap(ia)  =alphap(ia)  -SIOTA*TWOPI
        END IF
     END DO
     CALL SORT_ALPHA(nalpha,nalphab,alphap,zetax,thetax,thetap,B_al,vds_al,j_al,nlambda,i_p)
  END DO
  IF(DEBUG) THEN
     DO ia=1,nalpha
        DO il=1,nalphab
           WRITE(3000+myrank,'(6(1pe13.5),2I5)') zetap(il),thetap(ia,il),&
                &  zetax(ia,il),thetax(ia,il),alphap(ia),B_al(ia,il),ia,il
        END DO
     END DO
     CALL FLUSH(3000+myrank)
  END IF

  !Calculate coefficients of the drift kinetic equation
  CALL COEFFICIENTS_DKE(npoint,i_w,i_l,nw,&
                 &  z1,t1,B1,hBpp1,vd1,&
                 &  zb,tb,Bb,hBppb,vdb,&
                 &  z2,t2,B2,hBpp2,vd2,&
                 &  nlambda,lambda,zlw,tlw,zrw,trw,& 
                 &  BI1,BI2,BI3,BI4,BI5,BI6,BI7,Nnmp,BI8)

  IF(MODELFI) THEN
     CALL FAST_ION_MODELS(s,is,ns,nalpha,nalphab,nlambda,lambda,i_p,npoint,&
          & BI1,BI3,BI4,BI6,zlw,zrw,thetap,theta,B_al,vds_al,tau,ia_out)
  ELSE
!     CALL FAST_ION_JMAP(s,is,ns,nalpha,nalphab,nlambda,lambda,i_p,npoint,&
!          & BI1,BI3,BI4,BI6,zlw,zrw,thetap,theta,tau)
#ifdef MPIandPETSc
     CALL FAST_ION_ORBITS(s,is,ns,nalpha,nalphab,nlambda,lambda,i_p,npoint,&
          &   BI1,BI3,BI4,BI6,zlw,tlw,zrw,trw,tau)!
#endif
  END IF
  IF(myrank.EQ.0) CALL CALCULATE_FRACTIONS(s(is),nalpha,nalphab,nlambda,lambda,i_p,npoint,&
       & thetap,B_al,vds_al,tau,ia_out)
  
  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE CALC_FAST_ION_CONFINEMENT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE FAST_ION_MODELS(vs,is,ns,nalpha,nalphab,nlambda,lambda,i_p,npoint,&
     BI1,BI3,BI4,BI6,zlw,zrw,thetap,theta,B_al,vds_al,tau,ia_out)

!-----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  USE KNOSOS_STELLOPT_MOD  
  IMPLICIT NONE
  !Input
  INTEGER is,ns,nalpha,nalphab,nlambda,i_p(nlambda,nalpha,nalphab),npoint
  REAL*8 vs(ns),thetap(nalpha,nalphab),theta(nalphab),B_al(nalpha,nalphab),vds_al(Nnmp,nalpha,nalphab)
  REAL*8 lambda(nlambda)
  REAL*8 BI1(npoint),BI3(npoint),BI4(npoint),BI6(npoint),zlw(npoint),zrw(npoint)
  !Output
  INTEGER ia_out(npoint)
  REAL*8 tau(npoint)
  !Others
  INTEGER, PARAMETER :: nmaps=8
  INTEGER, PARAMETER :: nalphaturn=1000
  INTEGER il,ial,jal,kal,na,ialp1,ila,ilap1,jla,il0(nalpha),ipoint,jpoint,sign,ifile,ig
  REAL*8 maxg(nlambda)!,maxgt(nlambda),maxgs(nlambda)
  REAL*8 BI(nmaps,nlambda,nalpha),BIe(nmaps,3*nalpha),BIi(nmaps)
  REAL*8 g(npoint),save_g(npoint),gpl(npoint),gsl(npoint)
  REAL*8 alph,dalpha,talpha,vda,s0,s,ds,tau_s,tau_a,tau_t,dtau,D(nlambda),tau_d(nlambda)
  REAL*8 ran,NAN
  LOGICAL bif(nlambda)
  INTEGER norb,iaorb(nalpha),iorb,itrans
  INTEGER, SAVE :: tnalpha
  REAL*8, ALLOCATABLE :: thetape(:,:)
  REAL*8 f_t,Gamma_c,Gamma_cc,Gamma_delta,Gamma_alpha,Gamma_sl
  !Time
  CHARACTER*30, PARAMETER :: routine="FAST_ION_MODELS"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart
#ifdef MPIandPETSc
  !Others
!  INTEGER ierr
  INCLUDE "mpif.h"
#endif

  CALL CPU_TIME(tstart)

  NAN=SQRT(MONE)
  s0=vs(is)

  !Model -1: f_trapped
  g=0.5
  CALL INTEGRATE_G_NEW(nalpha,nalphab,nlambda,lambda,i_p,npoint,&
          & g,.TRUE.,thetap,B_al,vds_al,f_t)   !Calculate f_{trapped}
  KN_FTR=f_t
!  gth=(2./PI)*ATAN((1.-s0)/(2.0*PI))
 
  !Model 0: Nemov's Gamma_c
  g=2.0*ATAN(BI3/ABS(BI4*atorflux))/PI
  CALL INTEGRATE_G_NEW(nalpha,nalphab,nlambda,lambda,i_p,npoint,g*g,.TRUE.,&
       & thetap,B_al,vds_al,Gamma_c)
  IF(KN_STELLOPT(4)) THEN
     KN_GMC=gamma_c
     IF(.NOT.KN_STELLOPT(5)) RETURN
  END IF
!  gamma_c=gamma_c*2.*SQ2/PI
  !Model 0b: Corrected Gamma_c
  g=2.0*ATAN(BI3/ABS(BI4*atorflux))/PI
  CALL INTEGRATE_G_NEW(nalpha,nalphab,nlambda,lambda,i_p,npoint,ABS(g),.TRUE.,&
       & thetap,B_al,vds_al,Gamma_cc)
  
  !Create table of gamma_c^* and other quantities (ignoring small ripples)
  !and find maximum gamma_c^*
  maxg=-1
  BI=NAN
  il0=0
  DO ila=nlambda,2,-1
     DO ial=1,nalpha 
        DO il=1,nalphab 
           jpoint=i_p(ila,ial,il)
           IF(jpoint.LE.1) CYCLE
           IF(il0(ial).EQ.0) il0(ial)=il
           IF(g(jpoint).GT.maxg(ila)) maxg(ila)=g(jpoint)
           IF(il.EQ.il0(ial)) THEN
              BI(1,ila,ial)=BI1(jpoint)
              BI(2,ila,ial)=2.0*ATAN(BI3(jpoint)/ABS(BI4(jpoint)*atorflux))/PI
              BI(3,ila,ial)=BI3(jpoint)
              BI(4,ila,ial)=BI4(jpoint)
              BI(6,ila,ial)=2*BI6(jpoint)/rad_R
           END IF
        END DO
        ipoint=i_p(ila,ial,il0(ial))
        IF(ipoint.LE.1) CYCLE
        DO il=1,nalphab !Ignore alphas covered by barely trapped orbit
           jpoint=i_p(ila,ial,il)
           DO jal=1,nalpha 
              IF(il0(jal).EQ.0) CYCLE  
              IF(jpoint.EQ.i_p(ila,jal,il0(jal))) jpoint=1
           END DO
           IF(jpoint.LE.1) CYCLE
           WRITE(iout,*) 'Ripple at',ila,ial,il
           IF(g(ipoint).LT.gth.AND.g(jpoint).GT.gth) WRITE(iout,*) 'WARNING: ripple has sb',g(ipoint),g(jpoint)
           IF(g(ipoint).GT.gth.AND.g(ipoint).LT.g(jpoint)) WRITE(iout,*) 'WARNING: ripple largest sb',g(ipoint),g(jpoint)
        END DO
     END DO
  END DO

  !Model I: largest superbanana
  save_g=g
  g=0
  DO ila=1,nlambda
!     WRITE(6200+myrank,'(3(1pe13.5),I4)') lambda(ila),maxg(ila),eps,ila-ila
     IF(maxg(ila).LT.gth) CYCLE
     DO ial=1,nalpha
        DO il=1,nalphab 
           jpoint=i_p(ila,ial,il)
           IF(jpoint.LE.1) CYCLE
           g(jpoint)=0.5
           IF(il.EQ.il0(ial)) BI(5,ila,ial)=1
        END DO
     END DO
  END DO
  CALL INTEGRATE_G_NEW(nalpha,nalphab,nlambda,lambda,i_p,npoint,g,.TRUE.,&
       & thetap,B_al,vds_al,Gamma_delta)
!  maxgs=maxg
!  DO jla=1,nlambda/16 !smooth curves
!     maxgt=maxg
!     DO ila=jla+1,nlambda-jla
!        maxg(ila)=(maxgt(ila-1)+maxgt(ila)+maxgt(ila+1))/3.
!        WRITE(6200+myrank,'(3(1pe13.5),I4)') lambda(ila),maxg(ila),eps,jla
!     END DO
!  END DO
!  maxg=maxgs
  g=save_g

  !Model II or III
  ia_out=0
  tau=10*TENDFI
  gpl=0
  !Check prompt losses
  DO jla=1,nlambda
     DO jal=1,nalpha
        tau_t=10*TENDFI*TWOEFIoZ
        tau_s=10*TENDFI*TWOEFIoZ
        tau_a=10*TENDFI*TWOEFIoZ
        DO sign=-1,1,+2
           IF(tau_t.LT.TENDFI*TWOEFIoZ) EXIT
           s=s0
           tau_a=0
           vda=0
           ialp1=jal
           ilap1=jla
           IF(sign.EQ.-1) WRITE(6200+myrank,'(5(1pe13.5),2I8)') tau_a/TWOEFIoZ,thetap(ialp1,1),lambda(ilap1),s,BI(6,ilap1,ialp1),jla,jal
           DO kal=1,nalpha
              ial=ialp1
              ila=ilap1
!              ipoint=i_p(ila,ial,il0(ial))
!              IF(ipoint.LE.1) EXIT
              IF((vda*siota*sign.LT.0.AND.ial.NE.jal).OR.BI(2,ila,ial).LT.-gth) THEN
                 tau_t=10*TENDFI*TWOEFIoZ !out of alpha-cone
                 EXIT
              END IF
              ialp1=ial+sign
              IF(ialp1.GT.nalpha) ialp1=1
              IF(ialp1.LT.1) ialp1=nalpha
              ilap1=jla
              IF(ila.GT.nlambda/2) THEN
                 DO WHILE(i_p(ilap1,ialp1,il0(ialp1)).LT.1)
                    ilap1=ilap1-1
                 END DO
              END IF
!              jpoint=i_p(ilap1,ialp1,il0(ialp1))
              dalpha=sign*siota*(thetap(ialp1,1)-thetap(ial,1))
              IF(dalpha.LT.0) dalpha=dalpha+TWOPI
              vda=vda+(BI(4,ila,ial)/BI(1,ila,ial))*dalpha
!              IF(jpoint.LE.1) EXIT
              tau_a=tau_a+2*dalpha/ABS(BI(4,ila,ial)+BI(4,ilap1,ialp1)) !time that it takes to reach alpha_out
              s=s+dalpha*(BI(3,ila,ial)+BI(3,ilap1,ialp1))/(BI(4,ila,ial)+BI(4,ilap1,ialp1)) !flux surface
!              tau_a=tau_a+dalpha/ABS(BI(4,ilap1,ialp1)) !time that it takes to reach alpha_out
!              s=s+dalpha*(BI(3,ilap1,ialp1))/(BI(4,ilap1,ialp1)) !flux surface
              IF(vda*siota*sign.GT.0) WRITE(6200+myrank,'(5(1pe13.5),2I8)') tau_a/TWOEFIoZ,thetap(ialp1,1),lambda(ilap1),s,BI(6,ilap1,ialp1),jla,jal
              IF(s.GT.1.OR.(BI(2,ila,ial).GT.gth.AND.BI(2,ilap1,ialp1).LT.BI(2,ila,ial))) THEN !superbanana found
                 IF(s.LT.1) THEN
                    tau_s=(1-s)*atorflux/(BI(3,ila,ial)/BI(1,ila,ial))  !select small ripple? JL
                 ELSE
                    tau_s=0
                 END IF
                 tau_t=tau_s+tau_a
                 EXIT
              END IF
              IF(s.LT.0) s=-s
           END DO
           IF(kal.GE.nalpha) EXIT
        END DO
        tau_s=tau_s/TWOEFIoZ
        tau_a=tau_a/TWOEFIoZ
        tau_t=tau_t/TWOEFIoZ
        IF(tau_t.LT.2E4) BI(7,jla,jal)=tau_t
!        ipoint=i_p(jla,jal,il0(jal))
        IF(tau_t.LT.TENDFI) THEN
           DO il=1,nalphab
              jpoint=i_p(jla,jal,il)
              IF(jpoint.LE.1) CYCLE
!              IF(jpoint.GT.1.AND.jpoint.EQ.ipoint) tau(jpoint)=tau_t
              tau(jpoint)=tau_t
              gpl(jpoint)=0.5
              ia_out(jpoint)=ial
           END DO
        END IF
        WRITE(6300+myrank,'(8(1pe13.5),2(I4))') thetap(jal,1),lambda(jla),&
             & tau_t,tau_s,tau_a
     END DO
  END DO
  
  !Check stochastic losses
  IF(.FALSE.) THEN
  gsl=0
  bif=.FALSE.
  DO jla=2,nlambda
     itrans=0
     talpha=0
!     IF(jla.GT.nlambda/2) CYCLE  !confined or prompt losses
     vda=(BI(4,jla,1)/BI(1,jla,1))
     IF(vda*siota.GT.0) THEN
        sign=1
     ELSE
        sign=-1
     END IF
     ila=jla
     ial=1
     D(jla)=0
     s=s0
     ds=0
     tau_a=0
     tau_t=100*TENDFI*TWOEFIoZ
     DO jal=1,nalphaturn*nalpha
        ialp1=ial+sign
        IF(ialp1.GT.nalpha) ialp1=1
        IF(ialp1.LT.1) ialp1=nalpha
        ilap1=jla
        IF(ila.GT.nlambda/2) THEN
           DO WHILE(i_p(ilap1,ialp1,il0(ialp1)).LT.1)
              ilap1=ilap1-1
           END DO
        END IF
        dalpha=sign*siota*(thetap(ialp1,1)-thetap(ial,1))
        IF(dalpha.LT.0) dalpha=dalpha+TWOPI
!        rdJda=(BI(6,ilap1,ialp1)-BI(6,ila,ial))/(thetap(ialp1,1)-thetap(ial,1))
!        pdJda=(BI(3,ila,ial)+BI(3,ilap1,ialp1))/rad_R
!        IF(ABS(rdJda/pdJda-1).GT.2) THEN
        ipoint=i_p(ila,ial,il0(ial))
        iorb=0
        DO kal=1,nalpha
           jpoint=i_p(ila,kal,il0(kal))
           IF(ipoint.EQ.jpoint) THEN
              iorb=iorb+1
              iaorb(iorb)=kal
           END IF
        END DO
        itrans=itrans+1
!        IF(itrans.GT.1) iaorb=ial!no random        
        norb=iorb
        IF(norb.GT.1) THEN
           bif(jla)=.TRUE.
           CALL RANDOM_NUMBER(ran)
           ran=ran*norb
           DO iorb=1,norb
              IF(iorb-1.LT.ran.AND.ran.LT.iorb) THEN
                 ial=iaorb(iorb)
                 EXIT                 
              END IF
           END DO
           ialp1=ial+sign
           IF(ialp1.GT.nalpha) ialp1=1
           IF(ialp1.LT.1) ialp1=nalpha
        END IF
        talpha=talpha+dalpha
        tau_a=tau_a+2*dalpha/ABS(BI(4,ilap1,ialp1)+BI(4,ila,ial))
        s=s+dalpha*(BI(3,ilap1,ialp1)+BI(3,ila,ial))/(BI(4,ilap1,ialp1)+BI(4,ila,ial)) !flux surface
        ial=ialp1
        ila=ilap1
        IF(MOD(jal,nalpha).EQ.0) THEN
           IF(bif(jla)) THEN
              D(jla)=D(jla)+(s-s0)*(s-s0)
              ds=ds+ABS(s-s0)
              s0=s
           ELSE
              EXIT
           END IF
!           IF(WRITE(6200+myrank,'(5(1pe13.5),2I8,3(1pe13.5))') tau_a/TWOEFIoZ,thetap(ialp1,1),lambda(jla),s,BI(6,ilap1,ialp1),&
!                & -ilap1,-jal,D/(2*tau_a/TWOEFIoZ),ds/(jal/nalpha),tau_a/TWOEFIoZ/(jal/nalpha)
        END IF
!        IF(jal.LT.10*nalpha) THEN
        IF(talpha.LT.TWOPI) THEN
           ipoint=i_p(jla,ialp1,il0(ialp1))
           WRITE(6200+myrank,'(5(1pe13.5),3I8,5(1pe13.5))') &
             & tau_a/TWOEFIoZ,thetap(ialp1,1),lambda(jla),s,BI(6,ilap1,ialp1),-ilap1,-jal,&
             & itrans,talpha,zlw(ipoint),zrw(ipoint)
        END IF
     END DO
     D(jla)=D(jla)/(2*tau_a/TWOEFIoZ)
     s0=vs(is)
     tau_d(jla)=tau_a/TWOEFIoZ
  END DO

  DO jla=2,nlambda
     IF(.NOT.bif(jla).AND.bif(jla-1)) WRITE(iout,'(" (lambda_s~=",1pe23.16,") ")') 0.5*(lambda(jla-1)+lambda(jla))
!     IF(.NOT.bif(jla).AND.bif(jla-1)) THEN
!        bif(jla)=.TRUE.
!        DO kla=jla+1,nlambda
!           IF(bif(kla)) THEN
!              D(jla)=(D(kla  )*(lambda(jla)-lambda(jla-1))+ &
!                   &  D(jla-1)*(lambda(kla)-lambda(jla  )))/(lambda(kla)-lambda(jla-1))
!              EXIT
!           END IF
!        END DO
!     END IF

     tau_t=(1+s0)*(1-s0)/(2*D(jla))
     IF(tau_t.LT.5*TENDFI) THEN
        IF(RANDOMSL) THEN
           dtau=tau_d(jla)/nalphaturn
           IF(dtau.LT.tau_t/100.) dtau=tau_t/100.
           ds=SQRT(2*D(jla)*dtau)
        END IF
        DO ial=1,nalpha
           IF(RANDOMSL) THEN
              tau_t=0
              s=s0
              DO WHILE(s.LT.1)
                 CALL RANDOM_NUMBER(ran)                 
                 IF(ran.GT.0.5) THEN
                    s=s+ds
                 ELSE
                    s=s-ds
                 END IF
                 IF(s.LT.0) s=-s
                 tau_t=tau_t+dtau
              END DO
           END IF
           BI(7,jla,ial)=tau_t
           DO il=1,nalphab
              jpoint=i_p(jla,ial,il)
              IF(jpoint.LE.1) CYCLE
              IF(tau(jpoint).GE.TENDFI) THEN
                 tau(jpoint)=-tau_t
                 IF(tau_t.LT.TENDFI) gsl(jpoint)=0.5
              END IF
           END DO
           WRITE(6300+myrank,'(8(1pe13.5),2(I4))') thetap(ial,1),-lambda(jla),&
                & tau_t,tau_s,tau_d(jla)/nalphaturn,D(jla)
        END DO
     END IF
  END DO
  END IF
  
  !Model II: Gamma_alpha
  CALL INTEGRATE_G_NEW(nalpha,nalphab,nlambda,lambda,i_p,npoint,gpl,.TRUE.,&
       & thetap,B_al,vds_al,Gamma_alpha)
  KN_GMA=gamma_alpha
  
  CALL INTEGRATE_G_NEW(nalpha,nalphab,nlambda,lambda,i_p,npoint,gsl,.TRUE.,&
       & thetap,B_al,vds_al,Gamma_sl)

  tnalpha=2*nalpha
  IF(aiota/nzperiod.GE.1) tnalpha=3*nalpha
  ALLOCATE(thetape(tnalpha,nalphab))
  DO il=1,nalphab
     IF(aiota/nzperiod.LT.1) THEN
        IF(iota.GT.0) THEN
           thetape(     1:nalpha,il)=thetap(1:nalpha,il)-TWOPI
        ELSE
           thetape(     1:nalpha,il)=thetap(1:nalpha,il)+TWOPI
        END IF
        thetape(nalpha+1:tnalpha,il)=thetap(1:nalpha,il)
     ELSE
        IF(iota.GT.0) THEN
           thetape(       1:  nalpha,il)=thetap(1:nalpha,il)-2*TWOPI
           thetape(nalpha+1:2*nalpha,il)=thetap(1:nalpha,il)-  TWOPI
        ELSE
           thetape(       1:  nalpha,il)=thetap(1:nalpha,il)+2*TWOPI
           thetape(nalpha+1:2*nalpha,il)=thetap(1:nalpha,il)+  TWOPI
        END IF
        thetape(2*nalpha+1: tnalpha,il)=thetap(1:nalpha,il)
     END IF
  END DO

  DO ila=1,nlambda
     IF(aiota/nzperiod.LT.1) THEN
        DO ig=1,nmaps
           BIe(ig,        1: nalpha)=BI(ig,ila,:)
           BIe(ig, nalpha+1:tnalpha)=BI(ig,ila,:)
        END DO
     ELSE
        DO ig=1,nmaps
           BIe(ig,         1:  nalpha)=BI(ig,ila,:)
           BIe(ig,  nalpha+1:2*nalpha)=BI(ig,ila,:)
           BIe(ig,2*nalpha+1: tnalpha)=BI(ig,ila,:)
        END DO
     END IF
     IF(JMAP) THEN
        na=INT(2*nalphab*SQRT(s0))
     ELSE
        na=nalpha!INT(nalpha*1.5)
     END IF
     DO ial=1,na
        alph=theta(1)+(ial-1)*TWOPI/na
        DO ig=1,nmaps
           CALL LAGRANGE(thetape(1:tnalpha,1),BIe(ig,1:tnalpha),tnalpha,&
                & alph,BIi(ig),1)
        END DO
        IF(JMAP) THEN
           ifile=6000+myrank
        ELSE
           ifile=6100+myrank
        END IF
         WRITE(ifile,'(20(1pe13.5))') &
             & vs(is),alph,lambda(ila),&
             & 2.0*ATAN(BIi(3)/ABS(BIi(4)*atorflux))/PI,BIi(3),BIi(4),BIi(6),BIi(5),BIi(7),BIi(8)!,1/(TWOEFIoZ*BIi(3)/BIi(1)/atorflux)
      END DO
   END DO

   WRITE(6700+myrank,'(8(1pe13.5))') Gamma_c,Gamma_cc,Gamma_delta,Gamma_alpha,f_t
   
   CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE FAST_ION_MODELS



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef MPIandPETSc     
!!$
!!$SUBROUTINE FAST_ION_JMAP(vs,is,ns,nalpha,nalphab,nlambda,lambda,i_p,npoint,&
!!$     & BI1,BI3,BI4,BI6,zlw,zrw,thetap,theta,tau)
!!$
!!$!-----------------------------------------------------------------------------------------------
!!$!----------------------------------------------------------------------------------------------- 
!!$
!!$  USE GLOBAL
!!$  USE KNOSOS_STELLOPT_MOD  
!!$  IMPLICIT NONE
!!$  !Input
!!$  INTEGER is,ns,nalpha,nalphab,nlambda,i_p(nlambda,nalpha,nalphab),npoint
!!$  REAL*8 vs(ns),thetap(nalpha,nalphab),theta(nalphab)
!!$  REAL*8 lambda(nlambda)
!!$  REAL*8 BI1(npoint),BI3(npoint),BI4(npoint),BI6(npoint),zlw(npoint),zrw(npoint)
!!$  !Output
!!$  REAL*8 tau(npoint)
!!$  !Others
!!$  LOGICAL sgrid
!!$  INTEGER, PARAMETER :: nmaps=8
!!$  INTEGER, PARAMETER :: nalphaturn=1000
!!$  REAL*8, PARAMETER :: F5o12 =0.416666666666667
!!$  REAL*8, PARAMETER :: F13o12=1.083333333333333
!!$  INTEGER, ALLOCATABLE :: j_p(:,:,:)
!!$  INTEGER il,ial,jal,na,ila,ilamin,ilamax,jla,il0(nalpha),ipoint,jpoint,ifile,ig,nla
!!$  INTEGER ia_out(npoint),ia0,ia1,is0,is1
!!$  REAL*8 fa0,fa1,dta,signa,fs0,fs1,dts,signs
!!$  REAL*8, ALLOCATABLE :: BI(:,:,:)
!!$  REAL*8 BIe(nmaps,3*nalpha),BIi(nmaps),BIg(ns,nalphab,nmaps)
!!$  REAL*8 alpha,dsda,dadt,dsdt,Jsa,dlambda
!!$  REAL*8 g(npoint),gpl(npoint),gsl(npoint)
!!$  REAL*8 da,va(nalphab),talpha,s0,s,ds,tau_t
!!$!  LOGICAL bif(nlambda)
!!$!  INTEGER norb,iaorb(nalpha),iorb,itrans,it
!!$  INTEGER, SAVE :: tnalpha
!!$  REAL*8, ALLOCATABLE :: thetape(:,:)
!!$  REAL*8 MODANG2!,dummy
!!$  REAL*8 la(2,1),lambdac,lambdab
!!$  !Time
!!$  CHARACTER*30, PARAMETER :: routine="FAST_ION_JMAP"
!!$  INTEGER, SAVE :: ntotal=0
!!$  REAL*8,  SAVE :: ttotal=0
!!$  REAL*8,  SAVE :: t0=0
!!$  REAL*8 tstart
!!$!#ifdef MPIandPETSc
!!$  !Others
!!$  INTEGER ierr
!!$  INCLUDE "mpif.h"
!!$!#endif
!!$
!!$  zrw=zrw
!!$  zlw=zlw
!!$  TENDFI=1E-3
!!$  
!!$  CALL CPU_TIME(tstart)
!!$
!!$  !Create table of gamma_c^* and other quantities (ignoring small ripples)
!!$  !and find maximum gamma_c^*
!!$  s0=vs(is)
!!$  g=2.0*ATAN(BI3/ABS(BI4*atorflux))/PI
!!$  il0=0
!!$
!!$  la(2,1)=myrank
!!$  la(1,1)=lambda(1)
!!$  CALL MPI_ALLREDUCE(MPI_IN_PLACE,la,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD,ierr)
!!$  lambdac=la(1,1)
!!$  la(1,1)=lambda(nlambda)
!!$  CALL MPI_ALLREDUCE(MPI_IN_PLACE,la,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD,ierr)
!!$  lambdab=la(1,1)
!!$  dlambda=lambda(2)-lambda(1)
!!$  nla=(lambdab-lambdac+ALMOST_ZERO)/dlambda+1
!!$
!!$  ALLOCATE(j_p(nla,nalpha,nalphab))
!!$  j_p=0
!!$  DO ila=1,nla
!!$     IF(ABS(lambdac+(ila-1)*dlambda-lambda(1)).LT.0.5*dlambda) THEN
!!$        ilamin=ila
!!$        EXIT
!!$     END IF
!!$  END DO
!!$  ilamax=nlambda+ilamin-1
!!$  DO ila=1,nlambda
!!$     j_p(ila+ilamin-1,:,:)=i_p(ila,:,:)
!!$  END DO
!!$  ALLOCATE(BI(nmaps,nla,nalpha))
!!$  BI=NAN
!!$
!!$  DO ila=nla,2,-1
!!$     DO ial=1,nalpha
!!$        IF(j_p(ila,ial,nalphab/2).GT.1)  il0(ial)=nalphab/2
!!$        DO il=1,nalphab 
!!$           jpoint=j_p(ila,ial,il)
!!$           IF(jpoint.LE.1) CYCLE
!!$           IF(il0(ial).EQ.0) il0(ial)=il
!!$           IF(il.EQ.il0(ial)) THEN
!!$              BI(1,ila,ial)=BI1(jpoint)
!!$              BI(2,ila,ial)=BI3(jpoint)/ABS(BI4(jpoint)*atorflux)
!!$              BI(3,ila,ial)=BI3(jpoint)/atorflux
!!$              BI(4,ila,ial)=BI4(jpoint)
!!$              BI(6,ila,ial)=2*BI6(jpoint)/rad_R
!!$           END IF
!!$        END DO
!!$        ipoint=j_p(ila,ial,il0(ial))
!!$        IF(ipoint.LE.1) CYCLE
!!$        DO il=1,nalphab 
!!$           jpoint=j_p(ila,ial,il)
!!$           DO jal=1,nalpha
!!$              IF(jpoint.EQ.j_p(ila,jal,il0(jal))) jpoint=1
!!$           END DO
!!$           IF(jpoint.LE.1) CYCLE
!!$!           WRITE(iout,*) 'Ripple at',ila,ial,il
!!$!           IF(g(ipoint).LT.gth.AND.g(jpoint).GT.gth) WRITE(iout,*) 'WARNING: ripple has sb',g(ipoint),g(jpoint)
!!$!           IF(g(ipoint).GT.gth.AND.g(ipoint).LT.g(jpoint)) WRITE(iout,*) 'WARNING: ripple largest sb',g(ipoint),g(jpoint)
!!$        END DO
!!$     END DO
!!$        
!!$  END DO
!!$
!!$  !Figures
!!$  tnalpha=2*nalpha
!!$  IF(aiota/nzperiod.GE.1) tnalpha=3*nalpha
!!$  ALLOCATE(thetape(tnalpha,nalphab))
!!$  DO il=1,nalphab
!!$     IF(aiota/nzperiod.LT.1) THEN
!!$        IF(iota.GT.0) THEN
!!$           thetape(     1:nalpha,il)=thetap(1:nalpha,il)-TWOPI
!!$        ELSE
!!$           thetape(     1:nalpha,il)=thetap(1:nalpha,il)+TWOPI
!!$        END IF
!!$        thetape(nalpha+1:tnalpha,il)=thetap(1:nalpha,il)
!!$     ELSE
!!$        IF(iota.GT.0) THEN
!!$           thetape(       1:  nalpha,il)=thetap(1:nalpha,il)-2*TWOPI
!!$           thetape(nalpha+1:2*nalpha,il)=thetap(1:nalpha,il)-  TWOPI
!!$        ELSE
!!$           thetape(       1:  nalpha,il)=thetap(1:nalpha,il)+2*TWOPI
!!$           thetape(nalpha+1:2*nalpha,il)=thetap(1:nalpha,il)+  TWOPI
!!$        END IF
!!$        thetape(2*nalpha+1: tnalpha,il)=thetap(1:nalpha,il)
!!$     END IF
!!$  END DO
!!$  IF(JMAP) THEN        
!!$     na=INT(2*nalphab*SQRT(s0))
!!$  ELSE
!!$     na=nalphab!INT(nalpha*1.5)
!!$  END IF
!!$      
!!$  ia_out=0
!!$  tau=10*TENDFI
!!$  gpl=0
!!$  gsl=0
!!$
!!$
!!$  DO jla=1,nla
!!$     BIg=0
!!$     IF(aiota/nzperiod.LT.1) THEN
!!$        DO ig=3,nmaps
!!$           BIe(ig,        1: nalpha)=BI(ig,jla,:)
!!$           BIe(ig, nalpha+1:tnalpha)=BI(ig,jla,:)
!!$        END DO
!!$     ELSE
!!$        DO ig=3,nmaps
!!$           BIe(ig,         1:  nalpha)=BI(ig,jla,:)
!!$           BIe(ig,  nalpha+1:2*nalpha)=BI(ig,jla,:)
!!$           BIe(ig,2*nalpha+1: tnalpha)=BI(ig,jla,:)
!!$        END DO
!!$     END IF
!!$     
!!$     DO ial=1,na
!!$        va(ial)=theta(1)+(ial-1)*TWOPI/na
!!$        DO ig=3,6
!!$           CALL LAGRANGE(thetape(1:tnalpha,1),BIe(ig,1:tnalpha),tnalpha,&
!!$                & va(ial),BIi(ig),0)
!!$!           dummy=BIi(ig)
!!$!           IF(ISNAN(dummy)) BIi(ig)=0
!!$           BIg(myrank+1,ial,ig)=BIi(ig)
!!$           CALL MPI_ALLREDUCE(MPI_IN_PLACE,BIg(:,ial,ig),ns,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
!!$!           IF(ISNAN(dummy)) BIi(ig)=dummy
!!$        END DO
!!$
!!$        IF(JMAP) THEN
!!$           ifile=6000+myrank
!!$        ELSE
!!$           ifile=6100+myrank
!!$        END IF
!!$        WRITE(ifile,'(20(1pe13.5))') &
!!$             & vs(is),va(ial),lambdac+(jla-1)*dlambda,&
!!$             & 2.0*ATAN(BIi(3)/ABS(BIi(4)))/PI,BIi(3),BIi(4),BIi(6),BIi(5),BIi(7),BIi(8),BIi(3)/BIi(1)/atorflux
!!$        CALL FLUSH(ifile)
!!$     END DO
!!$
!!$     IF(jla.LT.ilamin.OR.jla.GT.ilamax) CYCLE
!!$
!!$     ds=vs(myrank+2)-vs(myrank+1)
!!$     da=TWOPI/na
!!$
!!$!     IF(jla.NE.nlambda/2.OR.myrank.NE.12) CYCLE
!!$     
!!$     DO jal=1,nalpha
!!$
!!$        IF(myrank.NE.12) CYCLE
!!$!        IF(j_p(jla,jal,il0(jal)).LE.1) CYCLE
!!$!        IF(jla.NE.16) CYCLE
!!$!        IF(jal.NE.13) CYCLE
!!$        
!!$        talpha=0
!!$        tau_t=0
!!$        s=s0
!!$        alpha=MODANG2(thetap(jal,1),TWOPI)
!!$        sgrid=.TRUE.
!!$
!!$
!!$        DO WHILE(s.LT.1.AND.tau_t.LT.TENDFI*TWOEFIoZ)
!!$           
!!$           IF(s.LT.0) s=-s
!!$           CALL FLUSH(iout)
!!$           ia0=INT(1+alpha/da)
!!$           ia1=ia0+1
!!$           IF(ia1.GT.na) ia1=1
!!$           is0=INT(ns*s+0.5)
!!$           is1=is0+1
!!$           IF(sgrid) THEN
!!$              fa1=(alpha-va(ia0))/da
!!$              fa0=1-fa1
!!$              dsdt=fa0*BIg(is0,ia0,3)+fa1*BIg(is0,ia1,3)
!!$              dadt=fa0*BIg(is0,ia0,4)+fa1*BIg(is0,ia1,4)
!!$              Jsa =fa0*BIg(is0,ia0,6)+fa1*BIg(is0,ia1,6)
!!$              IF(ISNAN(Jsa)) WRITE(iout,*) 'fa',is0,ia0,ia1,fa0,fa1,BIg(is0,ia0,6),BIg(is0,ia1,6)
!!$                 
!!$           ELSE
!!$              IF(is1.GT.ns) THEN
!!$                 fs1=(vs(ns)-s)/ds
!!$                 fs0=1-fs1
!!$                 dsdt=fs0*BIg(ns,ia0,3)+fs1*BIg(ns-1,ia0,3)
!!$                 dadt=fs0*BIg(ns,ia0,4)+fs1*BIg(ns-1,ia0,4)
!!$                 Jsa =fs0*BIg(ns,ia0,6)+fs1*BIg(ns-1,ia0,6)
!!$              ELSE
!!$                 fs1=(s-vs(is0))/ds
!!$                 fs0=1-fs1
!!$                 dsdt=fs0*BIg(is0,ia0,3)+fs1*BIg(is1,ia0,3)
!!$                 dadt=fs0*BIg(is0,ia0,4)+fs1*BIg(is1,ia0,4)
!!$                 Jsa =fs0*BIg(is0,ia0,6)+fs1*BIg(is1,ia0,6)
!!$              END IF
!!$              IF(ISNAN(Jsa)) WRITE(iout,*) 'fs',ia0,is0,is1,fs0,fs1,BIg(is0,ia0,6),BIg(is1,ia0,6)
!!$           END IF
!!$
!!$           dsda=dsdt/dadt
!!$
!!$           WRITE(6200+myrank,'(7(1pe13.5),2I8,10(1pe13.5))') &
!!$                & tau_t/TWOEFIoZ,alpha,lambdac+(jla-1)*dlambda,s,Jsa,dsda,dadt,jla,jal
!!$           CALL FLUSH(6200+myrank)
!!$                      
!!$           IF(DEBUG) WRITE(6200+myrank,'(7(1pe13.5),2I8,10(1pe13.5))') &
!!$                & tau_t/TWOEFIoZ,alpha,lambdac+(jla-1)*dlambda,s,Jsa,dsda,dadt,jla,jal
!!$           IF(DEBUG) CALL FLUSH(6200+myrank)
!!$           
!!$           dts=ds/ABS(dsdt)
!!$           dta=da/ABS(dadt)           
!!$           signa=dadt/ABS(dadt)
!!$           signs=dsdt/ABS(dsdt)
!!$           IF(tau_t.GT.ALMOST_ZERO) THEN
!!$              IF(sgrid) THEN
!!$                 IF(signa.GT.0) THEN
!!$                    IF(ia1.NE.1) THEN
!!$                       dta=(va(ia1)-alpha)/dadt
!!$                    ELSE
!!$                       dta=(TWOPI-alpha)/dadt
!!$                    END IF
!!$                 ELSE
!!$                    dta=(va(ia0)-alpha)/dadt
!!$                 END IF
!!$              ELSE
!!$                 IF(signs.GT.0) THEN
!!$                    IF(is1.GT.ns) THEN
!!$                       dts=(1-s)/dsdt
!!$                    ELSE
!!$                       dts=(vs(is1)-s)/dsdt
!!$                    END IF
!!$                 ELSE
!!$                    dts=(vs(is0)-s)/dsdt
!!$                 END IF
!!$              END IF
!!$           END IF
!!$
!!$!           IF(dts.LT.0) THEN
!!$!              WRITE(iout,*) s,is0,is1,vs(is0),vs(is1),signs,dsdt
!!$!              stop
!!$!           END IF
!!$
!!$           
!!$           IF(dts.LT.dta) THEN
!!$              IF(sgrid) THEN
!!$                 s=s+signs*ds
!!$              ELSE
!!$                 sgrid=.TRUE.
!!$                 IF(signs.GT.0) THEN
!!$                    IF(is1.GT.ns) THEN
!!$                       s=1.00001
!!$                    ELSE
!!$                       s=vs(is1)
!!$                    END IF
!!$                 ELSE
!!$                    s=vs(is0)
!!$                 END IF
!!$              END IF
!!$              alpha=MODANG2(alpha+dadt*dts,TWOPI)
!!$              tau_t=tau_t+dts
!!$           ELSE
!!$              IF(.NOT.sgrid) THEN
!!$                 alpha=MODANG2(alpha+signa*da,TWOPI)
!!$              ELSE
!!$                 sgrid=.FALSE.
!!$                 IF(signa.GT.0) THEN
!!$                    alpha=va(ia1)
!!$                 ELSE
!!$                    alpha=va(ia0)
!!$                 END IF
!!$              END IF
!!$              s=s+dsdt*dta
!!$              tau_t=tau_t+dta
!!$
!!$           END IF
!!$           
!!$        END DO
!!$
!!$        tau_t=tau_t/TWOEFIoZ
!!$!        ipoint=j_p(jla,jal,il0(jal))
!!$        IF(tau_t.LT.TENDFI) THEN
!!$           BI(7,jla,jal)=tau_t
!!$           DO il=1,nalphab
!!$              jpoint=j_p(jla,jal,il)
!!$              IF(jpoint.LE.1) CYCLE
!!$              tau(jpoint)=tau_t
!!$!              IF(jpoint.GT.1.AND.jpoint.EQ.ipoint) tau(jpoint)=tau_t
!!$!              ia_out(jpoint)=ial
!!$           END DO
!!$        END IF
!!$        WRITE(6300+myrank,'(3(1pe13.5),2(I4))') thetap(jal,1),lambdac+(jla-1)*dlambda,tau_t,jal,jla
!!$           
!!$     END DO
!!$        
!!$
!!$  END DO
!!$
!!$  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!$
!!$  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)
!!$
!!$END SUBROUTINE FAST_ION_JMAP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



SUBROUTINE FAST_ION_ORBITS(vs,is,ns,nalpha,nalphab,nlambda,lambda,i_p,npoint,&
     & BI1,BI3,BI4,BI6,zlw,tlw,zrw,trw,tau)

!-----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  USE KNOSOS_STELLOPT_MOD  
  IMPLICIT NONE
  !Input
  INTEGER is,ns,nalpha,nalphab,nlambda,i_p(nlambda,nalpha,nalphab),npoint
  REAL*8 vs(ns),lambda(nlambda),zlw(npoint),tlw(npoint),zrw(npoint),trw(npoint)
  REAL*8 BI1(npoint),BI3(npoint),BI4(npoint),BI6(npoint)
  !Output
  REAL*8 tau(npoint)
  !Others
  INTEGER irank(npoint)
  LOGICAL CONVERGED_Q
  LOGICAL computed(npoint),transition,conv_step,stochastic,step_a,split
  LOGICAL time_off,lost,prompt_loss,stochastic_loss,periodic,confined_sb,error_J,end_orbit
  INTEGER ila,ial,il,ipoint,it,iorbit,norbit,naturn
  REAL*8 da,ds,dat,dst,smin,smax,fd,dt,dalpha,dalpha_max
  REAL*8 t,s0,alpha0,mu_o_E,E_o_mu,J0!,alpha_sep
  REAL*8 s,alpha,theta,zeta,zl,tl,zr,tr,dla,J
  REAL*8 dJds,dJda,dsdt,dadt,dsda,ptran,dban
  REAL*8 s_new,alpha_new,theta_new,zeta_new,zl_new,tl_new,zr_new,tr_new,dB_new,J_new
  REAL*8 dJds_new,dJda_new,dsdt_new,dadt_new,dsda_new,ptran_new,dban_new
  REAL*8 s_test,alpha_test,theta_test,zeta_test,zl_test,tl_test,zr_test,tr_test,dB_test,J_test
  REAL*8 dJds_test,dJda_test,dJdla_test,dsdt_test,dadt_test,dsda_test,ptran_test,dban_test
  REAL*8 ran,J_trans,dsdt_old,zeta_split,J_split,dJds_split,dJda_split
  !Time
  CHARACTER*30, PARAMETER :: routine="FAST_ION_ORBITS"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart
!#ifdef MPIandPETSc
  !Others
  INTEGER ierr
  INCLUDE "mpif.h"
!#endif
  
  CALL CPU_TIME(tstart)
  
  s0=vs(is)

  da=(TWOPI/nalpha)*FIDELTA
  ds=0.05*FIDELTA
  tau=-TENDFI
  
  iorbit=0
  norbit=0
  computed=.FALSE.
  IF(myrank.EQ.0) THEN
     DO ipoint=2,npoint
        CALL RANDOM_NUMBER(ran)
        irank(ipoint)=INT(ran*numprocs)
     END DO
  END IF
  CALL MPI_BCAST(irank,npoint,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  DO ipoint=2,npoint
     IF(irank(ipoint).EQ.myrank) THEN
        norbit=norbit+1
     ELSE
        computed(ipoint)=.TRUE.
     END IF
  END DO
  CALL RANDOM_NUMBER(ran)
  WRITE(iout,'(" Process ",I8," calculates ",I8," orbits (",F6.4,")")') myrank,norbit,ran

  DO ila=1,nlambda
     
     mu_o_E=lambda(ila)
     E_o_mu=1./mu_o_E

     CALL CALC_SMIN_SMAX(E_o_mu,s0,smin,smax)
     
     DO ial=1,nalpha

        DO il=1,nalphab

           ipoint=i_p(ila,ial,il)
!           IF(ipoint.NE.3153) CYCLE
           IF(JORBIT.GT.0.AND.ipoint.NE.JORBIT) CYCLE
           IF(ipoint.LE.1.OR.computed(ipoint)) CYCLE
           iorbit=iorbit+1
           WRITE(iout,'(" Calculating orbit ",I8," out of ",I8)') iorbit,norbit
           WRITE(iout,'(" s_min= ",F6.4," s_max=",F6.4)') smin,smax
           
           IF(.NOT.RSEED) CALL INIT_RANDOMSEED(myrank+ipoint)
           computed(ipoint)=.TRUE.
           stochastic=.FALSE.
           it=0
           fd=1.0
           end_orbit=.FALSE.

           
           CALL SET_INITIAL_CONDITION(s0,zlw(ipoint),tlw(ipoint),zrw(ipoint),trw(ipoint),&
                & E_o_mu,BI1(ipoint),BI3(ipoint),BI4(ipoint),BI6(ipoint),&
                & t,s,alpha,theta,zeta,zl,tl,zr,tr,J,dJds,dJda,dsdt,dadt,dsda)

           IF(DEBUG) WRITE(6200+myrank,'(17(1pe13.5),3L,F8.4,4(I6))') t,s,alpha,theta,zeta,&
                mu_o_E,zl,tl,zr,tr,J,J0,dJda,dJds,dsdt,dadt,fd,step_a,&
                & transition,stochastic,ran,ila,ial,il,ipoint

           naturn=0
           alpha0=alpha
           dalpha_max=0
           J0    =J
           
           DO WHILE(.NOT.end_orbit)

              it=it+1
              split=.FALSE.
              conv_step=.FALSE.
              fd=1.0
              ptran_new=1.0
              J_trans=-1
              dsdt_old=0.0

              DO WHILE(.NOT.conv_step)
                 ptran=ptran_new
                 CALL FORWARD_STEP(it,s,alpha,theta,zeta,E_o_mu,zl,tl,zr,tr,dla,J,dJds,dJda,dsdt,dadt,dsda,&
                      & ds,da,dst,dat,step_a,smin,smax,fd,&
                      & s_new,alpha_new,theta_new,zeta_new,zl_new,tl_new,zr_new,tr_new,dB_new,J_new,&
                      & transition,J0,dJds_new,dJda_new,dsdt_new,dadt_new,dsda_new,ptran_new,dban_new)
                 IF((transition.OR.ISNAN(dsdt_new)).AND.fd.LT.1000) THEN
                    IF(transition) stochastic=.TRUE.
                    IF(.NOT.ISNAN(dsdt_new).AND.&
                         &(fd.GT.REAL(FDMAX).OR.&
                         & ISNAN(dsdt_old).OR.&
                         & (fd.GT.1.AND.CONVERGED_Q(ptran_new,ptran,PREC_TRANS).AND.&
                         & CONVERGED_Q(J_new,J_trans,PREC_BINT)))) THEN
!                       EXIT
                       fd=1.0
                       conv_step=.TRUE.
!                       WRITE(6200,*) dsdt_new
!                       WRITE(6200,*) fd,FDMAX
!                       WRITE(6200,*) dsdt_old
!                       WRITE(6200,*) CONVERGED_Q(ptran_new,ptran,PREC_TRANS),ptran_new,ptran                      
!                       WRITE(6200,*) CONVERGED_Q(J_new,J_trans,PREC_BINT),J_new,J_trans                       
                    ELSE
                       fd=fd*2.0
                       J_trans=J_new
                       CALL INTERPOLATE_FIELD(s,.TRUE.,E_o_mu)
                    END IF
                 ELSE
!                    EXIT
                    fd=1.0
                    conv_step=.TRUE.
                 END IF

              END DO

              J_split=0
              dJds_split=0
              dJda_split=0
              IF(transition) THEN

                 IF(ABS(zl_new-zl).LT.ABS(zr_new-zr)) THEN
                    zeta_test =0.5*(zr+zr_new)
                 ELSE
                    zeta_test =0.5*(zl+zl_new)
                 END IF                 
                 IF(J_new.GT.J) THEN
                    split=.FALSE.
                    IF(JTRANS) THEN !
                       s_test=s
                       CALL INTERPOLATE_FIELD(s,.TRUE.,E_o_mu)
                       theta_test=theta+iota*(zeta_test-zeta)
                    END IF !
                 ELSE
                    split=.TRUE.
                    s_test=s_new
                    theta_test=theta_new+iota*(zeta_test-zeta_new)
                    IF(JTRANS) zeta_split=zeta_test
                 END IF
                 IF(split.OR.JTRANS) CALL CALC_DSDA(it,theta_test,zeta_test,E_o_mu,zl_test,tl_test,zr_test,tr_test,dB_test,&
                      &  J_test,dJds_test,dJda_test,dJdla_test,dsdt_test,dadt_test,dsda_test,ptran_test,dban_test)
                 IF(split) THEN
                    IF(.NOT.FITRANS) THEN
                       ptran_test=0!J_test
                       ptran_new =1!J_new
                    END IF
                    CALL RANDOM_NUMBER(ran)
                    IF(ran.LT.ptran_test/(ptran_new+ptran_test)) THEN
                       alpha_test=alpha_new
                       IF(JTRANS) THEN
                          J_split=J_new
                          dJds_split=dJds_new
                          dJda_split=dJda_new
                       END IF
                       CALL COPY_FI(s_test,alpha_test,theta_test,zeta_test,zl_test,tl_test,zr_test,tr_test,&
                            & dB_test,J_test,dJds_test,dJda_test,dsdt_test,dadt_test,dsda_test,ptran_test,dban_test,&
                            & s_new,alpha_new,theta_new,zeta_new,zl_new,tl_new,zr_new,tr_new,&
                            & dB_new,J_new,dJds_new,dJda_new,dsdt_new,dadt_new,dsda_new,ptran_new,dban_new)
                    ELSE IF(JTRANS) THEN
                       J_split=J_test 
                       dJds_split=dJds_test
                       dJda_split=dJda_test
                    END IF
                 ELSE

                    IF(JTRANS) J0=J0+J_test
                    IF(split.OR.JTRANS) CALL INTERPOLATE_FIELD(s_new,.TRUE.,E_o_mu)

                 END IF

              END IF

              IF(JCORRECTION.AND.J_new.GT.ALMOST_ZERO.AND.s_new.GT.0.AND.s_new.LT.1) THEN

                 IF(JTRANS.OR..NOT.transition) THEN

                    CALL CORRECTION_STEP(split,it,s_new,alpha_new,theta_new,zeta_new,E_o_mu,&
                         & zl_new,tl_new,zr_new,tr_new,dB_new,&
                         & J_new,dJds_new,dJda_new,dsdt_new,dadt_new,dsda_new,&
                         & J0,ds,da,dst,dat,step_a,smin,smax,fd,&
                         & zeta_split,J_split,dJds_split,dJda_split)
!                    IF(J_new.LT.0) THEN
!                       J_new=-J_new
!                       CALL CORRECTION_STEP(split,it,s_new,alpha_new,theta_new,zeta_new,E_o_mu,&
!                            & zl_new,tl_new,zr_new,tr_new,dB_new,&
!                            & J_new,dJds_new,dJda_new,dsdt_new,dadt_new,dsda_new,&
!                            & J,ds,da,dst,dat,step_a,smin,smax,fd,&
!                            & zeta_split,J_split,dJds_split,dJda_split)
!                       IF(J_new.LT.0) J_new=-J_new
!                    END IF
                       
                 ELSE

                    J0=J_new
                    
                 END IF
                    
              END IF

              IF(step_a) THEN
                 dt=(alpha_new-alpha)/(dadt*TWOEFIoZ)
              ELSE
                 dt=(s_new-s)/(dsdt*TWOEFIoZ)
              END IF
              IF(dt.LT.0) WRITE(6200+myrank,*) 'dt',dt
              t =t+dt

              CALL COPY_FI(s_new,alpha_new,theta_new,zeta_new,zl_new,tl_new,zr_new,tr_new,&
                   & dB_new,J_new,dJds_new,dJda_new,dsdt_new,dadt_new,dsda_new,ptran_new,dban_new,&
                   & s,alpha,theta,zeta,zl,tl,zr,tr,dla,J,dJds,dJda,dsdt,dadt,dsda,ptran,dban)
              
              IF(DEBUG) WRITE(6200+myrank,'(17(1pe13.5),3L,F8.4,4(I6))') t,s,alpha,theta,zeta,&
                   mu_o_E,zl,tl,zr,tr,J,J0,dJda,dJds,dsdt,dadt,fd,step_a,transition,stochastic,ptran,ila,ial,il,ipoint
              dalpha=ABS(alpha-alpha0)

              IF(dalpha.GT.dalpha_max) dalpha_max=dalpha
!!$              IF(INT(dalpha/TWOPI).GT.naturn.AND.stochastic) THEN
!!$                 D=D+(s_new-sa0)*(s_new-sa0)
!!$                 sa0=s_new
!!$                 naturn=naturn+1
!!$              END IF
              
              time_off=t.GT.TENDFI
              lost=s.GT.1.0-dban
              prompt_loss=lost.AND..NOT.stochastic
              stochastic_loss=lost.AND.stochastic
              periodic=.NOT.stochastic.AND.dalpha.GT.2*TWOPI
              confined_sb=.NOT.stochastic.AND.t.GT.10*tsl.AND.dalpha_max.LT.2*TWOPI
              error_J=(J.LE.ALMOST_ZERO).OR.fd.GT.1000
              end_orbit=time_off.OR.lost.OR.periodic.OR.confined_sb.OR.error_J
              IF(end_orbit) THEN
                 WRITE(6300+myrank,'(I6,2(1pe13.5),8L)') ipoint,t,s,time_off,lost,&
                      & prompt_loss,stochastic_loss,periodic,confined_sb,error_J
                 EXIT
              END IF
           END DO
           IF(lost.OR.error_J) tau(ipoint)=t
           WRITE(iout,'(" Calculated ",F5.1,"% of orbits")') 100*REAL(iorbit)/norbit
        END DO
     END DO
  END DO

  DO ipoint=1,npoint
       IF(tau(ipoint).LT.0) tau(ipoint)=0
    END DO
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,tau,npoint,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  DO ipoint=1,npoint
     IF(tau(ipoint).LT.DTFI) tau(ipoint)=10*TENDFI
  END DO

  CALL INTERPOLATE_FIELD(s0,.TRUE.,ZERO)
  
  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE FAST_ION_ORBITS

#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALC_SMIN_SMAX(E_o_mu,s0,smin,smax)
  
  USE GLOBAL
  IMPLICIT NONE
  !Input
  REAL*8 E_o_mu,s0
  !Output
  REAL*8 smin,smax
  !Others
  INTEGER is,js
  REAL*8 Bold,stest,ds,Bmax_axis,Bmin_axis
  !Time
  CHARACTER*30, PARAMETER :: routine="CALC_SMIN_SMAX"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart
  
  CALL CPU_TIME(tstart)

!!$
  IF(ALLOCATED(Bmax_b)) DEALLOCATE(Bmax_b,Bmin_b)
  ALLOCATE(Bmax_b(ns_b),Bmin_b(ns_b))
  Bmax_b=0
  Bmin_b=0
  DO is=1,ns_b
     IF(js_b(is).EQ.0) CYCLE
     CALL INTERPOLATE_FIELD(s_b(js_b(is)),.TRUE.,E_o_mu)
!     CALL FIND_BMIN_BMAX(js_b(is),MAL,MAL,.TRUE.)
  END DO
!!$
  
  smin=1E-4
  smax=1.1  
  CALL INTERPOLATE_FIELD(smin,.TRUE.,E_o_mu)
  Bmin_axis=Bmin
  Bmax_axis=Bmax
  stest=smin
  !Look for B_max(s)=E_o_mu
  IF((Bmax_b(js_b(1))-E_o_mu)*(Bmax_axis-E_o_mu).LT.0) THEN
     is=1
     stest=smin
     Bold=Bmax_axis
  ELSE
     DO is=2,ns_b
        js=js_b(is)
        IF((Bmax_b(js)-E_o_mu)*(Bmax_b(js-1)-E_o_mu).LT.0) EXIT
     END DO
     IF(is.LT.ns_b+1) THEN
        stest=s_b(js-1)
        Bold=Bmax_b(js-1)
     END IF
  END IF
  IF(is.LT.ns_b+1) THEN
     ds=(s_b(js)-stest)/2
     DO WHILE(ds.GT.SMALL)
        stest=stest+ds
        CALL INTERPOLATE_FIELD(stest,.TRUE.,E_o_mu)
        IF((Bold-E_o_mu)*(Bmax-E_o_mu).LT.0) THEN
           stest=stest-ds
        ELSE
           Bold=Bmax
        END IF
        ds=ds/2
     END DO
     IF(stest.LT.s0) THEN
        smin=stest
     ELSE
        smax=stest
     END IF
  END IF

  stest=smin
  !Look for B_min(s)=E_o_mu
  IF((Bmin_b(js_b(1))-E_o_mu)*(Bmin_axis-E_o_mu).LT.0) THEN
     is=1
     stest=smin
     Bold=Bmin_axis
  ELSE
     DO is=2,ns_b
        js=js_b(is)
        IF((Bmin_b(js)-E_o_mu)*(Bmin_b(js-1)-E_o_mu).LT.0) EXIT
     END DO
     IF(is.LT.ns_b+1) THEN
        stest=s_b(js-1)
        Bold=Bmin_b(js-1)
     END IF
  END IF
  IF(is.LT.ns_b+1) THEN
     ds=(s_b(js)-stest)/2
     DO WHILE(ds.GT.SMALL)
        stest=stest+ds
        CALL INTERPOLATE_FIELD(stest,.TRUE.,E_o_mu)
        IF((Bold-E_o_mu)*(Bmin-E_o_mu).LT.0) THEN
           stest=stest-ds
        ELSE
           Bold=Bmin
        END IF
        ds=ds/2
     END DO
     IF(stest.LT.s0.AND.stest.GT.smin) THEN
        smin=stest
     ELSE IF(stest.GT.s0.AND.stest.LT.smax) THEN
        smax=stest
     END IF
  END IF  

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)
  
END SUBROUTINE CALC_SMIN_SMAX


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE SET_INITIAL_CONDITION(s0,zlw,tlw,zrw,trw,E_o_mu,BI1,BI3,BI4,BI6,&
                & t,s,alpha,theta,zeta,zl,tl,zr,tr,J,dJds,dJda,dsdt,dadt,dsda)

  USE GLOBAL
  IMPLICIT NONE
  !Input
  REAL*8 s0,zlw,tlw,zrw,trw,E_o_mu,BI1,BI3,BI4,BI6
  !Output
  REAL*8 t,s,alpha,theta,zeta,zl,tl,zr,tr,J,dJds,dJda,dsdt,dadt,dsda
  !Others
  REAL*8 taub,MODANG2
  
  t=0
  s=s0
  CALL INTERPOLATE_FIELD(s,.TRUE.,E_o_mu)
  theta =MODANG2(0.5*(tlw+trw),TWOPI)
  zeta  =MODANG2(0.5*(zlw+zrw),TWOPI/nzperiod)
  alpha=MODANG2(theta-iota*zeta,TWOPI)
  zl=zlw+ zeta-0.5*(zlw+zrw)
  tl=tlw+theta-0.5*(tlw+trw)
  zr=zrw+ zeta-0.5*(zlw+zrw)
  tr=trw+theta-0.5*(tlw+trw)

  taub= BI1
  dJda=  BI3*sgnB
  dJds= -BI4*torflux
  J   = BI6    !J/v=2*BI6(jpoint)/rad_R
  dsdt= dJda/(taub*torflux)
  dadt=-dJds/(taub*torflux)
  dsda=dsdt/dadt

END SUBROUTINE SET_INITIAL_CONDITION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE FORWARD_STEP(it,s,alpha,theta,zeta,E_o_mu,zl,tl,zr,tr,dB,J,dJds,dJda,dsdt,dadt,dsda,&
                    & ds,da,dst,dat,step_a,smin,smax,fd,&
                    & s_new,alpha_new,theta_new,zeta_new,zl_new,tl_new,zr_new,tr_new,dB_new,&
                    & J_new,transition,J0,dJds_new,dJda_new,dsdt_new,dadt_new,dsda_new,ptran_new,dban_new)
  
  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER it
  REAL*8 s,alpha,theta,zeta,E_o_mu,zl,tl,zr,tr,dB,J,dJds,dJda,dsdt,dadt,dsda,ds,da,smin,smax
  !Input/output
  LOGICAL step_a
  REAL*8 dst,dat,fd
  !Output
  LOGICAL transition,CONVERGED_Q
  REAL*8 s_new,alpha_new,theta_new,zeta_new
  REAL*8 zl_new,tl_new,zr_new,tr_new,dB_new
  REAL*8 J_new,J0,dJds_new,dJda_new
  REAL*8 dsdt_new,dadt_new,dsda_new
  REAL*8 ptran_new,dban_new
  !Others
  REAL*8 J_ext,one_o_lambda,dlambda,dJdla_new,dadla,dad1la!,dsdla,dsd1la
  REAL*8 NAN
  REAL*8, SAVE :: damax,dsmax
  REAL*8, PARAMETER :: dJtrans=0.2
  !Time
  CHARACTER*30, PARAMETER :: routine="FORWARD_STEP"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart
  
  CALL CPU_TIME(tstart)

  NAN=SQRT(MONE)
  
  IF(it.GT.0) THEN
     IF(ABS(dsda).LT.ds/da) THEN
        step_a=.TRUE.
     ELSE
        step_a=.FALSE.
     END IF
  END IF

  
  IF(it.GT.0) THEN
     IF(step_a) THEN
        dat=da*dadt/(fd*ABS(dadt))
        dst=dat*dsda
     ELSE
        dst=ds*dsdt/(fd*ABS(dsdt))
        dat=dst/dsda
     END IF
     dsmax=ABS(dst)
     damax=ABS(dat)
  ELSE
     IF(step_a) THEN
        dst=(J0-J)/dJds !JL
        dsmax=ds!/fd
!           dsmax=ABS(dsda*da/fd)/2
        IF(ABS(dst).GT.dsmax) dst=dsmax*dst/ABS(dst)
        dat=0
     ELSE
        dst=0
        dat=(J0-J)/dJda
        damax=da!/fd
!           damax=ABS(ds/(dsda*fd))/2
        IF(ABS(dat).GT.damax) dat=damax*dat/ABS(dat)
     END IF
  END IF

  s_new    =s    +dst

  IF(s_new.GT.1) THEN!s_b(ns_b)) THEN !LCFS
     RETURN
  ELSE IF(s_new.LT.smin.OR.s_new.GT.smax) THEN
     IF(it.GT.0) THEN
        IF(s_new.LT.smin) THEN
           dst=(smin-s)
           s_new=smin
        END IF
!        step_a=.TRUE.        
!        dat=2*dat
!        dst=0
!        s_new=s
     ELSE
        s_new=s
        dst=0
        dat=0
        fd=-fd
        RETURN
     END IF
  END IF

  alpha_new=alpha+dat
  theta_new=theta+dat
  IF(ABS(dst).GT.ALMOST_ZERO) THEN
     theta_new=theta_new-iota*zeta
     CALL INTERPOLATE_FIELD(s_new,.TRUE.,E_o_mu)
     theta_new=theta_new+iota*zeta
  END IF
!  zl_new=zl
!  tl_new=tl+dat
!  zr_new=zr
  !  tr_new=tr+dat

  CALL CALC_DSDA(it,theta_new,zeta,E_o_mu,zl_new,tl_new,zr_new,tr_new,dB_new,&
       & J_new,dJds_new,dJda_new,dJdla_new,dsdt_new,dadt_new,dsda_new,ptran_new,dban_new)

  IF(it.GT.0) THEN
     IF(dB_new.GE.0) THEN
        IF(ABS(dst).GT.ALMOST_ZERO) THEN
           J_ext=J+dJda*dat
        ELSE
           J_ext=J+(dJds*dsda+dJda)*dat
        END IF
     ELSE
        IF(fd.GT.4) THEN
!           WRITE(6200+myrank,*) 'E_o_mu0',E_o_mu,s_new,theta,J0
           !           WRITE(6200+myrank,*) 'E_o_mu1',E_o_mu,s_new,theta_new,J_new
           IF(dB.LT.0) THEN
              one_o_lambda=E_o_mu+dB-dB_new
           ELSE
              one_o_lambda=E_o_mu-1.5*dB_new
           END IF
           CALL CALC_DSDA(it,theta_new,zeta,one_o_lambda,zl_new,tl_new,zr_new,tr_new,dB_new,&
                & J_new,dJds_new,dJda_new,dJdla_new,dsdt_new,dadt_new,dsda_new,ptran_new,dban_new)
           dlambda=(J0-J_new)/dJdla_new !           d(1/lambda)=-dlambda/lambda^2
!           WRITE(6200+myrank,*) 'E_o_mu2',one_o_lambda,s_new,theta_new,J_new
           one_o_lambda=one_o_lambda*(1-dlambda*one_o_lambda)
           IF(one_o_lambda.LT.Bmin.OR.one_o_lambda.GT.Bmax) THEN
              dsdt_new=NAN
              RETURN
           END IF
           CALL CALC_DSDA(it,theta_new,zeta,one_o_lambda,zl_new,tl_new,zr_new,tr_new,dB_new,&
                & J_new,dJds_new,dJda_new,dJdla_new,dsdt_new,dadt_new,dsda_new,ptran_new,dban_new)

!           IF(step_a) THEN
!              dsdla=-dJdla_new/dJds_new
!              dsd1la=-dsdla/(one_o_lambda*one_o_lambda)      !d_{1/x}x=d_y(1/y)=-1/y^2=-x^2,  y=1/x
!              s_new=s_new+(E_o_mu-one_o_lambda)*dsd1la
!              IF(s_new.LT.smin.OR.s_new.GT.smax) THEN
!                 dsdt_new=NAN
!                 RETURN
!              END IF
!              theta_new=theta_new-iota*zeta
!              CALL INTERPOLATE_FIELD(s_new,.TRUE.,E_o_mu)
!              theta_new=theta_new+iota*zeta
!           ELSE
           dadla=-dJdla_new/dJda_new
           dad1la=-dadla/(one_o_lambda*one_o_lambda)      !d_{1/x}x=d_y(1/y)=-1/y^2=-x^2,  y=1/x
           theta_new=theta_new+(E_o_mu-one_o_lambda)*dad1la
!           END IF
           CALL CALC_DSDA(it,theta_new,zeta,E_o_mu,zl_new,tl_new,zr_new,tr_new,dB_new,&
                & J_new,dJds_new,dJda_new,dJdla_new,dsdt_new,dadt_new,dsda_new,ptran_new,dban_new)
           J_ext=J0
           alpha_new=alpha+(theta_new-theta)
        END IF
     END IF        
  END IF

  zeta_new =0.5*(zl_new+zr_new)
  theta_new=0.5*(tl_new+tr_new)      
  IF(it.GT.0) transition=.FALSE.
  IF(it.GT.0.AND.(zr-zl)*nzperiod/TWOPI.GT.0.05.AND. &
       &(ABS(zl_new-zl).GT.PI/nzperiod.OR.ABS(zr_new-zr).GT.PI/nzperiod.OR.&
       &((.NOT.CONVERGED_Q(zl,zl_new,dJtrans).OR..NOT.CONVERGED_Q(zr,zr_new,dJTRANS)).AND.&
          .NOT.CONVERGED_Q(J_new,J_ext,dJtrans).AND..NOT.CONVERGED_Q(J_new,J0,dJTRANS)))) &
          transition=.TRUE.
  tl=tl
  tr=tr

  IF(JORBIT.GT.0.AND.DEBUG) WRITE(6200+myrank,'(17(1pe13.5),3L,F8.4,4(I6))') REAL(it),&
       s_new,alpha_new,theta_new,zeta_new,1/E_o_mu,&
       & zl_new,tl_new,zr_new,tr_new,J_new,J0,dJda_new,dJds_new,dsdt_new,dadt_new,fd,step_a,&
       & transition,transition,ptran_new!,ila,ial,il,ipoint

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)
  
END SUBROUTINE FORWARD_STEP



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALC_DSDA(it,theta,zeta,E_o_mu,zl,tl,zr,tr,dB,J,dJds,dJda,dJdla,dsdt,dadt,dsda,ptran,dban)
!-----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER it
  REAL*8 theta,zeta,E_o_mu
  !Input/output
  REAL*8 zl,tl,zr,tr
  !Output
  REAL*8 dB,J,dJds,dJda,dJdla,dsdt,dadt,dsda,ptran,dban
  !Others
  INTEGER flag,jt
  REAL*8 zeta_ini,theta_ini,taub,Q(nq0)!,z_l,t_l
  REAL*8 z1,t1,B1,hBpp1,vd1(nqv)
  REAL*8 zb,tb,Bb,hBppb,vdb(nqv)
  REAL*8 z2,t2,B2,hBpp2,vd2(nqv)
  REAL*8 B,dBdpsi,dBdz,dBdt,hBpp,dummy,vdummy(Nnmp)
  !Time
  CHARACTER*30, PARAMETER :: routine="CALC_DSDA"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)
!  IF(it.GE.90) WRITE(6200+myrank,*) Bmin,E_o_mu,Bmax
  
!!$  DO iextr=1,2
!!$     IF(iextr.EQ.1) THEN
!!$        z_l=zl
!!$        t_l=tl
!!$     ELSE
!!$        z_l=zr
!!$        t_l=tr
!!$     END IF
!!$     dzl=0.5*(zr-zl)
!!$     DO WHILE(dzl.GT.PREC_EXTR)
!!$        dzl_old=ABS(dzl)
!!$        CALL CALCB(z_l,t_l,2,.FALSE.,B,dBdz,dBdt,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,vdummy)
!!$        lambda1=1/B
!!$        dBdzl=(dBdz+dBdt*iota)   
!!$        dzl=(E_o_mu-lambda)/dBdzl
!!$        IF(ABS(dzl).GT.dzl_old) EXIT
!!$        z_l=z_l+dzl
!!$        t_l=t_l+dzl*iota
!!$     END DO
!!$     IF(ABS(E_o_mu-lambda)/E_o_mu.GT.1E-3) EXIT
!!$     IF(iextr.EQ.1) THEN
!!$        zl=z_l
!!$        tl=t_l
!!$     ELSE
!!$        zr=z_l
!!$        tr=t_l
!!$     END IF
!!$  EnD DO
!!$  CALL BOUNCE_INTEGRAL(it,z1,t1,z2,t2,E_mu,&
!!$       &             Bp1,hBpp1,vd1,  &
!!$       &             Bp2,hBpp2,vd2,  &
!!$       &             zbx,bbx,hBppbx,vdbx,nq,Q)


  Q=0
  zl=0
  zr=0
  tl=0
  tr=0
  
  zeta_ini=zeta
  theta_ini=theta
  B1=-1
  Bb=-1
  zb=0
  jt=0
  flag=1
  DO WHILE((B1.LT.E_o_mu.OR.flag.GT.0).AND.jt.LT.MAL)
     jt=jt+1
     WRITE(6200+myrank,*) 'here',B1,E_o_mu,Bmax,flag,jt,MAL
     CALL EXTREME_POINT(zeta_ini,theta_ini,-1,z1,t1,B1,hBpp1,vd1,flag)
     zeta_ini=z1
     theta_ini=t1
!     WRITE(6200+myrank,*) 'Bb',Bb,B1,E_o_mu
     IF(Bb.LT.B1.AND.B1.LT.E_o_mu) THEN
        zb=z1
        tb=t1
        Bb=B1
        hBppb=hBpp1
        vdb=vd1
     END IF
!     IF(ABS(it).eq.16) WRITE(6200+myrank,*) 'B1',B1,Bb
  END DO
  zeta_ini=zeta
  theta_ini=theta
  B2=-1
  flag=1
  jt=0
  DO WHILE((B2.LT.E_o_mu.OR.flag.GT.0).AND.jt.LT.MAL)
     jt=jt+1
     CALL EXTREME_POINT(zeta_ini,theta_ini,3,z2,t2,B2,hBpp2,vd2,flag)
     zeta_ini=z2
     theta_ini=t2
     IF(Bb.LT.B2.AND.(B2.LT.E_o_mu.OR.flag.GT.0)) THEN
        zb=z2
        tb=t2
        Bb=B2
        hBppb=hBpp2
        vdb=vd2
     END IF
!     IF(ABS(it).eq.16) WRITE(6200+myrank,*) 'B2',B2,Bb
  END DO

!!$     IF(jt.LT.MAL.AND.(dla.GT.0.OR.hBppb.LT.0.OR.Bb.GT.1E2)) THEN
!!$        IF(dla.LT.0.AND.Bb.LT.1E2) THEN
!!$           IF(zeta.LT.zb) THEN
!!$              z2=zb
!!$              t2=tb
!!$              B2=Bb
!!$              hBpp2=hBppb
!!$              vd2=vdb
!!$           ELSE
!!$              z1=zb
!!$              t1=tb
!!$              B1=Bb
!!$              hBpp1=hBppb
!!$              vd1=vdb
!!$           END IF
!!$        ELSE IF(dla.LT.0.AND.Bb.GT.1E2) THEN
!!$           zb=0.5*(z1+z2)
!!$           tb=0.5*(t1+t2)
!!$           Bb=E_o_mu/1.001
!!$           hBppb=0
!!$           vdb=0
!!$           IF(it.EQ.258) WRITE(6200+myrank,*) 'here'
!!$        END IF
     !     IF(Bb.LE.E_o_mu.AND.jt.LT.MAL) THEN

!  IF(ABS(it).eq.16) THEN
!     WRITE(6200+myrank,*) z1,zb,z2
!     WRITE(6200+myrank,*) B1,Bb,B2
!     WRITE(6200+myrank,*) hBpp1,hBppb,hBpp2
!     WRITE(6200+myrank,*) zeta,E_o_mu
!  END IF

  dB=E_o_mu-Bb
  
  IF(Bb.GT.0.AND.dB.GE.0.AND.jt.LT.MAL) &
       CALL BOUNCES(it,z1,t1,B1,hBpp1,vd1,&
       zb,tb,Bb,hBppb,vdb,&
       z2,t2,B2,hBpp2,vd2,&
       & E_o_mu,.FALSE.,nq0,Q,&
       & zl,tl,zr,tr)
  
  taub=  Q(1)
  dban=  Q(2)*VDoV/atorflux/2
  dJda=  Q(3)*sgnB
  dJds= -Q(4)*torflux
  J=     Q(6)
  dJdla=-ABS(Q(7))
  dsdt= dJda/(taub*torflux)
  dadt=-dJds/(taub*torflux)
  dsda=dsdt/dadt
!  WRITE(iout,*) 'dban',dban
  IF(B1.LT.B2) THEN
     CALL CALCB(z1,t1,3,.FALSE.,B,dBdz,dBdt,dBdpsi,hBpp,dummy,dummy,dummy,dummy,dummy,dummy,vdummy)
  ELSE
     CALL CALCB(z2,t2,3,.FALSE.,B,dBdz,dBdt,dBdpsi,hBpp,dummy,dummy,dummy,dummy,dummy,dummy,vdummy)
  END IF
  ptran=dBdpsi*atorflux*dJda-dBdt*dJds
  
  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)
     
END SUBROUTINE CALC_DSDA



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE COPY_FI(s_ori,alpha_ori,theta_ori,zeta_ori,zl_ori,tl_ori,zr_ori,tr_ori,&
     & dB_ori,J_ori,dJds_ori,dJda_ori,dsdt_ori,dadt_ori,dsda_ori,ptran_ori,dban_ori,&
     & s_end,alpha_end,theta_end,zeta_end,zl_end,tl_end,zr_end,tr_end,&
     & dB_end,J_end,dJds_end,dJda_end,dsdt_end,dadt_end,dsda_end,ptran_end,dban_end)

  USE GLOBAL
  IMPLICIT NONE
  !Input
  REAL*8 s_ori,alpha_ori,theta_ori,zeta_ori,zl_ori,tl_ori,zr_ori,tr_ori
  REAL*8 dB_ori,J_ori,dJds_ori,dJda_ori,dsdt_ori,dadt_ori,dsda_ori,ptran_ori,dban_ori
  !Output
  REAL*8 s_end,alpha_end,theta_end,zeta_end,zl_end,tl_end,zr_end,tr_end
  REAL*8 dB_end,J_end,dJds_end,dJda_end,dsdt_end,dadt_end,dsda_end,ptran_end,dban_end

  s_end    =s_ori
  alpha_end=alpha_ori
  theta_end=theta_ori
  zeta_end =zeta_ori
  zl_end   =zl_ori
  tl_end   =tl_ori
  zr_end   =zr_ori
  tr_end   =tr_ori
  dB_end   =dB_ori
  J_end    =J_ori
  dJds_end =dJds_ori
  dJda_end =dJda_ori
  dsdt_end =dsdt_ori
  dadt_end =dadt_ori
  dsda_end =dsda_ori
  ptran_end=ptran_ori
  dban_end =dban_ori

END SUBROUTINE COPY_FI


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CORRECTION_STEP(split,it,s,alpha,theta,zeta,E_o_mu,zl,tl,zr,tr,&
                         & dB,J,dJds,dJda,dsdt,dadt,dsda,&
                         & J0,ds,da,dst,dat,step_a,smin,smax,fd,&
                         & zeta_split,J_split,dJds_split,dJda_split)
  
  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER it
  LOGICAL split,step_a
  REAL*8 J0,ds,da,dst,dat,smin,smax,fd,zeta_split,J_split,dJds_split,dJda_split
  !Input/output
  REAL*8 s,alpha,theta,zeta,E_o_mu
  REAL*8 zl,tl,zr,tr,dB,J,dJds,dJda,dsdt,dadt,dsda,ptran,dban
  !Others
  LOGICAL transition
  REAL*8 dJ,dJ_test
  REAL*8 s_test,alpha_test,theta_test,zeta_test,zl_test,tl_test,zr_test,tr_test,dB_test,J_test
  REAL*8 dJds_test,dJda_test,dsdt_test,dadt_test,dsda_test,ptran_test,dban_test
  REAL*8 theta_split,zl_split,tl_split,zr_split,tr_split,dB_split
  REAL*8 dJdla_split,dsdt_split,dadt_split,dsda_split,ptran_split,dban_split
  REAL*8 REL_DIST
  !Time
  CHARACTER*30, PARAMETER :: routine="CORRECTION_STEP"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  dJ=REL_DIST(J0,J+J_split)
  DO WHILE(dJ.GT.PREC_J.AND.s.LT.1)!s_b(ns_b))!LCFS
     CALL FORWARD_STEP(-it,s,alpha,theta,zeta,E_o_mu,&
          & zl,tl,zr,tr,dB,&
          & J+J_split,dJds+dJds_split,dJda+dJda_split,dsdt,dadt,dsda,&
          & ds,da,dst,dat,step_a,smin,smax,fd,&
          & s_test,alpha_test,theta_test,zeta_test,zl_test,tl_test,zr_test,tr_test,dB_test,J_test,&
          & transition,J0,dJds_test,dJda_test,dsdt_test,dadt_test,dsda_test,ptran_test,dban_test)

     dJ_test=REL_DIST(J0,J_test)
     IF(split.AND.dJ_test.GT.dJ) THEN
        theta_split=theta+iota*(zeta_split-zeta)
        CALL CALC_DSDA(-it,theta_split,zeta_split,E_o_mu,zl_split,tl_split,zr_split,tr_split,&
             & dB_split,J_split,dJds_split,dJda_split,dJdla_split,&
             & dsdt_split,dadt_split,dsda_split,ptran_split,dban_split)
        dJ_test=REL_DIST(J0,J_test+J_split)
     END IF
     IF(dJ_test.GT.dJ.OR.fd.LT.0) THEN
!        J=-J
        IF(fd.LT.0) fd=-fd
        IF(step_a) CALL INTERPOLATE_FIELD(s,.TRUE.,E_o_mu)
        EXIT
     END IF
     
     dJ=dJ_test
     CALL COPY_FI(s_test,alpha_test,theta_test,zeta_test,zl_test,tl_test,zr_test,tr_test,&
          & dB_test,J_test,dJds_test,dJda_test,dsdt_test,dadt_test,dsda_test,ptran_test,dban_test,&
          & s,alpha,theta,zeta,zl,tl,zr,tr,&
          & dB,J,dJds,dJda,dsdt,dadt,dsda,ptran,dban)

  END DO

  IF(split.AND.JTRANS) J0=J

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE CORRECTION_STEP


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALCULATE_FRACTIONS(vs,nalpha,nalphab,nlambda,lambda,i_p,npoint,&
     & thetap,B_al,vds_al,tau,ia_out)

  USE GLOBAL
  USE KNOSOS_STELLOPT_MOD  
  IMPLICIT NONE
  !Input
  INTEGER nalpha,nalphab,nlambda,i_p(nlambda,nalpha,nalphab),npoint
  INTEGER ia_out(npoint)
  REAL*8 vs,thetap(nalpha,nalphab),lambda(nlambda),B_al(nalpha,nalphab),vds_al(Nnmp,nalpha,nalphab),tau(npoint)
  !Others
  INTEGER, PARAMETER :: ntime=10
  REAL*8 tlosses(ntime) /1E-4,2E-4,3E-4,4E-4,5E-4,6E-4,7E-4,8E-4,9E-4,1E-3/
  INTEGER il,ial,ila,jla,ipoint,it,nt
  REAL*8 g(npoint),newg(npoint),gpl(npoint),gsl(npoint),gt(npoint)
  REAL*8 t,f_t,f_pl,f_sl,f_loss,f_loss1ms,f_tl,f_lossl(nlambda),f_losslt(nlambda)
  REAL*8 f_pll(nlambda),f_pllt(nlambda),f_sll(nlambda),f_sllt(nlambda),f_lt(ntime,nlambda)
  REAL*8 f_ta,f_pla,f_sla,f_lossa,f_enda
  !Time
  CHARACTER*30, PARAMETER :: routine="CALCULATE_FRACTIONS"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart
#ifdef MPIandPETSc
  !Others
!  INTEGER ierr
  INCLUDE "mpif.h"
#endif

  CALL CPU_TIME(tstart)

  IF(.NOT.MODELFI) THEN
     DO ipoint=2,npoint
        IF(tau(ipoint).GT.tsl.AND.tau(ipoint).LT.TENDFI) tau(ipoint)=-tau(ipoint)
     END DO
  END IF

  !Calculate f_{loss,t}
  g=0
  f_t=-1
  t=DTFI
  IF(LINEART) nt=TENDFI/DTFI  
  DO WHILE(t.LT.TENDFI)
     IF(f_t.LT.0) g=0.5
     gpl=0
     gsl=0
     DO ipoint=2,npoint
        IF(tau(ipoint).GT.0.AND.tau(ipoint).LT.t) THEN
           gpl(ipoint)=0.5
        ELSE IF(ABS(tau(ipoint)).LT.t) THEN
           gsl(ipoint)=0.5
        END IF
     END DO
     CALL INTEGRATE_G_NEW(nalpha,nalphab,nlambda,lambda,i_p,npoint,gpl,.TRUE.,&
          & thetap,B_al,vds_al,f_pl)
     CALL INTEGRATE_G_NEW(nalpha,nalphab,nlambda,lambda,i_p,npoint,gsl,.TRUE.,&
          & thetap,B_al,vds_al,f_sl)
     f_loss=f_pl+f_sl
     IF((t-1e-3)*(t*1.2-1e-3).LT.0) f_loss1ms=f_loss
     IF(f_t.LT.0) CALL INTEGRATE_G_NEW(nalpha,nalphab,nlambda,lambda,i_p,npoint,&
          & g,.TRUE.,thetap,B_al,vds_al,f_t)   !Calculate f_{trapped}
     WRITE(6400+myrank,'(8(1pe13.5))') t,f_loss,f_pl,f_sl,f_t
     IF(LINEART) THEN
        t=t+DTFI
     ELSE
        t=t*1.2
     END IF
  END DO
  
  !Calculate f_{loss,lambda}
  DO ila=1,nlambda
     gsl=0
     gpl=0
     g=0
     DO it=0,ntime
        gt=0
        DO ial=1,nalpha
           DO il=1,nalphab
              ipoint=i_p(ila,ial,il)
              IF(ipoint.GT.1) THEN
                 IF(it.EQ.0) THEN
                    g(ipoint)=1
                    IF(tau(ipoint).GT.0.AND.tau(ipoint).LT.TENDFI) THEN
                       gpl(ipoint)=1
                    ELSE IF(ABS(tau(ipoint)).LT.TENDFI) THEN
                       gsl(ipoint)=1
                    END IF
                 ELSE IF(ABS(tau(ipoint)).LT.tlosses(it)) THEN
                    gt(ipoint)=1
                 END IF
              END IF
           END DO
        END DO
        IF(it.EQ.0) THEN
           IF(ila.EQ.1) DEBUG=.TRUE.
           CALL INTEGRATE_G_NEW(nalpha,nalphab,nlambda,lambda,i_p,npoint,gpl,.TRUE.,&
                & thetap,B_al,vds_al,f_pll(ila))
           DEBUG=.FALSE.
           CALL INTEGRATE_G_NEW(nalpha,nalphab,nlambda,lambda,i_p,npoint,gsl,.TRUE.,&
                & thetap,B_al,vds_al,f_sll(ila))
           CALL INTEGRATE_G_NEW(nalpha,nalphab,nlambda,lambda,i_p,npoint,g,.TRUE.,&
                & thetap,B_al,vds_al,f_tl)
           f_pll(ila)=f_pll(ila)/f_tl
           f_sll(ila)=f_sll(ila)/f_tl
           f_lossl(ila)=f_pll(ila)+f_sll(ila)
        ELSE
           IF(it.EQ.ntime.AND.ila.EQ.1) DEBUG=.TRUE.
           CALL INTEGRATE_G_NEW(nalpha,nalphab,nlambda,lambda,i_p,npoint,gt,.TRUE.,&
                & thetap,B_al,vds_al,f_lt(it,ila))
           DEBUG=.FALSE.
           f_lt(it,ila)=f_lt(it,ila)/f_tl
        END IF
     END DO
     IF(MODELFI) THEN
        WRITE(6500+myrank,'(5(1pe13.5),I4,10(1pe13.5))') vs,lambda(ila),f_lossl(ila),&
             & f_pll(ila),f_sll(ila),jla-jla
     ELSE
        WRITE(6500+myrank,'(4(1pe13.5),I4,10(1pe13.5))') vs,lambda(ila),f_lossl(ila),&
             & f_pll(ila),f_sll(ila),jla-jla,(f_lt(it,ila),it=1,ntime)
     END IF
  END DO
  DO jla=1,nlambda/16 !smooth curves
     f_pllt=f_pll
     f_sllt=f_sll
     f_losslt=f_lossl
     DO ila=jla+1,nlambda-jla
        f_pll(ila)=(f_pllt(ila-1)+f_pllt(ila)+f_pllt(ila+1))/3.
        f_sll(ila)=(f_sllt(ila-1)+f_sllt(ila)+f_sllt(ila+1))/3.        
        f_lossl(ila)=(f_losslt(ila-1)+f_losslt(ila)+f_losslt(ila+1))/3.
        WRITE(6500+myrank,'(5(1pe13.5),I4)') vs,lambda(ila),f_lossl(ila),f_pll(ila),f_sll(ila),jla
     END DO
  END DO
  
  !Calculate f_{loss,alpha}, f_{birth,alpha} and f_{escape,alpha}
  DO ial=1,nalpha
     gpl=0
     gsl=0
     g=0
     DO ila=1,nlambda
        DO il=1,nalphab
           ipoint=i_p(ila,ial,il)
           IF(ipoint.GT.1) THEN
              g(ipoint)=1
              IF(tau(ipoint).GT.0.AND.tau(ipoint).LT.TENDFI) THEN                
                 gpl(ipoint)=1
              ELSE IF(ABS(tau(ipoint)).LT.TENDFI) THEN
                 gsl(ipoint)=1
              END IF
           END IF
        END DO
     END DO
     newg=0
     DO ipoint=2,npoint
        IF(ia_out(ipoint).EQ.ial) newg(ipoint)=1
     END DO
     CALL INTEGRATE_G_NEW(nalpha,nalphab,nlambda,lambda,i_p,npoint,gpl,.TRUE.,&
          & thetap,B_al,vds_al,f_pla)
     CALL INTEGRATE_G_NEW(nalpha,nalphab,nlambda,lambda,i_p,npoint,gsl,.TRUE.,&
          & thetap,B_al,vds_al,f_sla)
     f_lossa=f_pla+f_sla
     CALL INTEGRATE_G_NEW(nalpha,nalphab,nlambda,lambda,i_p,npoint,g,.TRUE.,&
          & thetap,B_al,vds_al,f_ta)     
     CALL INTEGRATE_G_NEW(nalpha,nalphab,nlambda,lambda,i_p,npoint,newg,.TRUE.,&
          & thetap,B_al,vds_al,f_enda)
     WRITE(6600+myrank,'(8(1pe13.5))') thetap(ial,1),f_lossa/f_ta,f_pla/f_ta,f_sla/f_ta,f_lossa/2,f_enda/2,f_ta/2
  END DO

  IF(.NOT.MODELFI) WRITE(6700+myrank,'(8(1pe13.5))') f_pl,f_loss-f_pl,f_loss

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE CALCULATE_FRACTIONS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
