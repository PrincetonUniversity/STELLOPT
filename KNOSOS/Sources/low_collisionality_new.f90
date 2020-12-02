
!Calculate neoclassical transport at low collisionalities
!Include 1/nu,sqrtnu and superbanana-plateau; do not need splitting B=B_0+delta*B_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALC_LOW_COLLISIONALITY_NEW(jv,Epsi,phi1c,Mbbnm,trMnm,D11,D31,&
     & nalphab,zeta,theta,dn1dv,dn1nm)

!--------------------------------------------------------------------------------------------- 
!Calculate monoenergetic transport coefficients D11 and D31 and contribution quasineutrality dn1dv and 
!dn1nm at nalphabxnalphab grid in (zeta,theta)  for collisionality cmul=nu(jv)/v(jv) and 
!normalized radial electric field Epsi/v(jv), and in the presence of phi1c
!-----------------------------------------------------------------------------------------------

  USE GLOBAL
  IMPLICIT NONE
!  Input
  INTEGER jv
  REAL*8 Epsi,phi1c(Nnmp),Mbbnm(Nnmp),trMnm(Nnmp)
  !Output
  INTEGER nalphab
  REAL*8 D11(Nnmp,Nnmp),D31,zeta(nax),theta(nax),dn1dv(nax,nax),dn1nm(Nnmp,Nnmp)
  !Others
  CHARACTER*100 serr
  INTEGER ial,ilambda
  INTEGER, SAVE :: nal,nlambda
  REAL*8 D11r(100,100)

  IF(.NOT.QN) CONVERGED=.FALSE.

  IF(CONVERGED.OR.(MLAMBDA.GT.0.AND.MAL.GT.0)) THEN

     IF(MLAMBDA.GT.0.AND.MAL.GT.0) THEN
        nlambda=MLAMBDA
        nal=MAL
        CONVERGED=.TRUE.
     END IF
     !Write monoenergetic transport coefficients using DKES normalization
     CALL CALC_LOW_COLLISIONALITY_NANL_NEW(nal,nlambda,jv,Epsi,phi1c,Mbbnm,trMnm,&
          & D11,D31,nalphab,zeta,theta,dn1dv,dn1nm)
     WRITE(200+myrank,'(3(1pe13.5)," NaN ",2(1pe13.5)," &
          & NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN")') &
          & nu(jv)/v(jv)/2.,Epsi/v(jv)*psip,vmconst(jv)/v(jv),fdkes(jv)*D11(1,1),fdkes(jv)*D11(1,1)
  ELSE

     D11=0
     CONVERGED=.FALSE.
     ial=1
     nal=32
     DO WHILE (nal.LE.nax)
        ilambda=1
        nlambda=32
        DO WHILE (nlambda.LE.nlambdax)
           CALL CALC_LOW_COLLISIONALITY_NANL_NEW(nal,nlambda,jv,Epsi,phi1c,Mbbnm,trMnm,&
                & D11r(ial,ilambda),D31,nalphab,zeta,theta,dn1dv,dn1nm)
           IF(ilambda.GT.1.AND.ABS(D11r(ial,ilambda)/D11r(ial,ilambda-1)-1.0).LT.PREC_DQDV) THEN
              EXIT
           END IF
           ilambda=ilambda+1
           nlambda=nlambda*2  
        END DO
        D11r(ial,ilambda+1:100)=D11r(ial,ilambda)
        IF(ial.GT.1.AND.nlambda.LE.nlambdax.AND.&
             & ABS(D11r(ial,ilambda)/D11r(ial-1,100)-1.0).LT.PREC_DQDV) THEN
           WRITE(200+myrank,'(3(1pe13.5)," NaN ",2(1pe13.5),&
                & " NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN")') &
                & nu(jv)/v(jv)/2.,Epsi/v(jv)*psip,vmconst(jv)/v(jv),fdkes(jv)*D11r(ial,ilambda),fdkes(jv)*D11r(ial,ilambda)
           EXIT
        END IF
        ial=ial+1
        nal=nal*2
     END DO
     
     IF(nal.GT.nax.OR.nlambda.GT.nlambdax) THEN
        IF(nal.GE.nax) THEN
           WRITE(1100+myrank,*) 'Increase nax'
        ELSE IF(nlambda.GE.nlambdax) THEN
           WRITE(1100+myrank,*) 'Increase nlambdax'
        END IF
        serr="DKE not converged"
        IF(.NOT.ONLY_DB) CALL END_ALL(serr,.FALSE.)
     ELSE
        CONVERGED=.TRUE.
        D11=D11r(ial,ilambda)
        WRITE(iout,'(" DKE converged with nal=",I6," and nlambda=",I6)') nal,nlambda
        CALL CALC_LOW_COLLISIONALITY_NANL_NEW(nal,nlambda,jv,Epsi,phi1c,Mbbnm,trMnm,&
             D11,D31,nalphab,zeta,theta,dn1dv,dn1nm)
     END IF

  END IF

END SUBROUTINE CALC_LOW_COLLISIONALITY_NEW


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALC_LOW_COLLISIONALITY_NANL_NEW(nal,nlambda,jv,Epsi,phi1c,Mbbnm,trMnm,&
     & D11,D31,nalphab,zeta,theta,dn1dv,dn1nm)

!--------------------------------------------------------------------------------------------- 
!Calculate monoenergetic transport coefficients D11 and D31 and contribution quasineutrality dn1dv and 
!dn1nm at nalphabxnalphab grid in (zeta,theta), for collisionality cmul=nu(jv)/v(jv) and normalized 
!radial electric field Epsi/v(jv)m and in the presence of phi1c
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
  INTEGER jv,nlambda,nal
  REAL*8 Epsi,phi1c(Nnmp),Mbbnm(Nnmp),trMnm(Nnmp),zeta(nax),theta(nax)
  !Output
  INTEGER nalphab
  REAL*8 D11(Nnmp,Nnmp),D31,dn1dv(nax,nax),dn1nm(Nnmp,Nnmp)
  !Others
  INTEGER, SAVE :: nalpha,nalphab_save,npoint
  INTEGER, SAVE, ALLOCATABLE :: nbif(:),i_p_ap1(:),i_p_am1(:),i_p_ap2(:),i_p_am2(:)
  INTEGER, SAVE, ALLOCATABLE :: i_p_ap1I(:),i_p_ap1II(:),i_p_am1I(:),  i_p_am1II(:)
  INTEGER, SAVE, ALLOCATABLE :: i_p_ap2I(:),i_p_ap2II(:),i_p_ap2III(:),i_p_ap2IV(:)  
  INTEGER, SAVE, ALLOCATABLE :: i_p_am2I(:),i_p_am2II(:),i_p_am2III(:),i_p_am2IV(:)
  REAL*8, SAVE, ALLOCATABLE :: wp1I(:),wp1II(:),wp2I(:),wp2II(:),wp2III(:),wp2IV(:)
  REAL*8, SAVE, ALLOCATABLE :: wm1I(:),wm1II(:),wm2I(:),wm2II(:),wm2III(:),wm2IV(:)
  REAL*8, SAVE :: offset,theta_save(nax)
  REAL*8, SAVE, ALLOCATABLE :: lambda(:),dalpha_am1(:),dalpha_ap1(:),dalpha_am2(:),dalpha_ap2(:)
  REAL*8, SAVE, ALLOCATABLE :: dlambda_lm1(:),dlambda_lp1(:)
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
  REAL*8 lambdac,dlambdap,dalphap(nturn+1)
  !Lambda neighbours
  INTEGER, PARAMETER :: nbifx=5  !maximum number of bifurcations allowed
  INTEGER, SAVE, ALLOCATABLE :: i_p_lm1(:),i_p_lp1(:,:)
  !Alpha neighbours
  REAL*8, ALLOCATABLE, SAVE :: zlw(:),zrw(:)
  REAL*8, ALLOCATABLE       :: tlw(:),trw(:)
  INTEGER ipoint,ila,il,ia
  !Drift-kinetic equation
  INTEGER, SAVE, ALLOCATABLE :: nnz(:)
  REAL*8, SAVE, ALLOCATABLE :: BI1(:),BI2(:),BI3(:),BI3f(:),BI3b(:),BI4(:),BI5(:),BI6(:),BI7(:),BI8(:,:),gint(:,:),factnu(:)
  REAL*8 omega,cmul
#ifdef MPIandPETSc
  !Petsc
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscpc.h>
#include <petsc/finclude/petscksp.h>
  PetscErrorCode ierr
  PetscInt ipointm1,jpointm1
  PetscInt, PARAMETER :: oneps=1
  PetscInt, SAVE :: innz(0:npointx-1)
  KSP, SAVE :: ksp
  Mat, SAVE :: matCOL,matVEAf,matVEAb,matVMAf,matVMAb
  Mat matA
  INTEGER jpoint
  PetscScalar mat_entry
#else
  REAL*8, SAVE, ALLOCATABLE :: COL(:,:),VEAf(:,:),VEAb(:,:),VMAf(:,:),VMAb(:,:)
  REAL*8, ALLOCATABLE :: mat(:,:)
#endif
  REAL*8, ALLOCATABLE :: rowCOL(:),rowVEAf(:),rowVEAb(:),rowVMAf(:),rowVMAb(:)
  !Time
  CHARACTER*40, PARAMETER :: routine="CALC_LOW_COLLISIONALITY_NANL_NEW"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart
!!$
!!$  Mbbnm=Mbbnm !To be removed
!!$  trMnm=trMnm !To be removed
!!$  D31=0
!!$
!!$  CALL CPU_TIME(tstart)
!!$
!!$  WRITE(iout,*) 'Calculating low collisionality regimes, NEW'
!!$
!!$  IF(CALCULATED_INT) GOTO 123
!!$
!!$  IF(PHI1_READ) bnmc0=bnmc0+2*borbic(0,0)*phi1c/vdconst(jv) 
!!$
!!$ !Find and characterize wells
!!$  ALLOCATE(connected(nwx,nwx),bottom(nwx),&
!!$       & z1(nwx),t1(nwx),B1(nwx),hBpp1(nwx),vd1(nqv,nwx),&
!!$       & zb(nwx),tb(nwx),Bb(nwx),hBppb(nwx),vdb(nqv,nwx),&
!!$       & z2(nwx),t2(nwx),B2(nwx),hBpp2(nwx),vd2(nqv,nwx),& 
!!$       & alphap_w(nwx),Bt(nwx),Btt(nwx),&
!!$       & lambdab_w(nwx),lambdac_w(nwx))
!!$  CALL CHARACTERIZE_WELLS(nal,na,nalpha,nw,z1,t1,B1,hBpp1,vd1, &
!!$       & zb,tb,Bb,hBppb,vdb, &
!!$       & z2,t2,B2,hBpp2,vd2, &
!!$       & Bt,Btt,alphap_w,dalphap(1),bottom,connected,offset)
!!$  !Resize arrays (nwx->nw)
!!$  ALLOCATE(temp(nw),temp2(nqv,nw),ltemp(nw))
!!$  temp=z1(1:nw);DEALLOCATE(z1);ALLOCATE(z1(nw));z1=temp
!!$  temp=zb(1:nw);DEALLOCATE(zb);ALLOCATE(zb(nw));zb=temp
!!$  temp=z2(1:nw);DEALLOCATE(z2);ALLOCATE(z2(nw));z2=temp
!!$  temp=t1(1:nw);DEALLOCATE(t1);ALLOCATE(t1(nw));t1=temp
!!$  temp=tb(1:nw);DEALLOCATE(tb);ALLOCATE(tb(nw));tb=temp
!!$  temp=t2(1:nw);DEALLOCATE(t2);ALLOCATE(t2(nw));t2=temp
!!$  temp=B1(1:nw);DEALLOCATE(B1);ALLOCATE(B1(nw));B1=temp
!!$  temp=Bb(1:nw);DEALLOCATE(Bb);ALLOCATE(Bb(nw));Bb=temp
!!$  temp=B2(1:nw);DEALLOCATE(B2);ALLOCATE(B2(nw));B2=temp
!!$  temp=hBpp1(1:nw);DEALLOCATE(hBpp1);ALLOCATE(hBpp1(nw));hBpp1=temp
!!$  temp=hBppb(1:nw);DEALLOCATE(hBppb);ALLOCATE(hBppb(nw));hBppb=temp
!!$  temp=hBpp2(1:nw);DEALLOCATE(hBpp2);ALLOCATE(hBpp2(nw));hBpp2=temp
!!$  temp=alphap_w(1:nw);DEALLOCATE(alphap_w);ALLOCATE(alphap_w(nw));alphap_w=temp
!!$  temp=lambdab_w(1:nw);DEALLOCATE(lambdab_w);ALLOCATE(lambdab_w(nw));lambdab_w=temp
!!$  temp=lambdac_w(1:nw);DEALLOCATE(lambdac_w);ALLOCATE(lambdac_w(nw));lambdac_w=temp
!!$  temp= Bt(1:nw);DEALLOCATE(Bt) ;ALLOCATE(Bt(nw)) ;Bt=temp
!!$  temp=Btt(1:nw);DEALLOCATE(Btt);ALLOCATE(Btt(nw));Btt=temp
!!$  ltemp=bottom(1:nw);DEALLOCATE(bottom);ALLOCATE(bottom(nw));bottom=ltemp
!!$  temp2=vd1(1:nqv,1:nw);DEALLOCATE(vd1);ALLOCATE(vd1(nqv,nw));vd1=temp2
!!$  temp2=vdb(1:nqv,1:nw);DEALLOCATE(vdb);ALLOCATE(vdb(nqv,nw));vdb=temp2
!!$  temp2=vd2(1:nqv,1:nw);DEALLOCATE(vd2);ALLOCATE(vd2(nqv,nw));vd2=temp2
!!$  DEALLOCATE(temp,temp2,ltemp)  
!!$  ALLOCATE(ltemp2(nw,nw))
!!$  ltemp2=connected(1:nw,1:nw);DEALLOCATE(connected);ALLOCATE(connected(nw,nw));connected=ltemp2
!!$  DEALLOCATE(ltemp2)
!!$  !Create grid in alpha, then (zeta,theta)
!!$  !Determine number of modes
!!$  nalphab=1
!!$  DO WHILE(nalphab.LT.nalpha*1.2)
!!$     nalphab=nalphab*2
!!$  END DO
!!$  nalphab=nalphab/2
!!$  nalphab_save=nalphab
!!$  IF(ALLOCATED(zetap)) DEALLOCATE(zetap,thetap,zetax,thetax,B_al,vds_al)
!!$  ALLOCATE(zetax(nalpha,nalphab),thetax(nalpha,nalphab),alphap(nalpha),&
!!$         & zetap(nalphab),thetap(nalpha,nalphab),&
!!$         & B_al(nalpha,nalphab),vds_al(Nnmp,nalpha,nalphab))
!!$  CALL CREATE_ANGULAR_GRID(na,nalpha,nalphab,alphap,dalphap,offset,&
!!$       & zetap,thetap,zetax,thetax,B_al,vds_al)
!!$  zeta(1:nalphab) =zetap  !square grid
!!$  theta(1:nalphab)=zetap*nzperiod 
!!$  theta_save=theta
!!$
!!$  CALL EXCLUDE_WELLS(na,nalpha,nalphab,nw,bottom,connected,&
!!$       & alphap_w,z1,zb,z2,Bb,Bt,zetax,thetax)
!!$  
!!$  IF(ALLOCATED(lambda)) DEALLOCATE(lambda,i_p)
!!$  ALLOCATE(lambda(nlambda),one_o_lambda(nlambda),i_p(nlambda,nalpha,nalphab))
!!$
!!$  !Set global grid in lambda
!!$  CALL CREATE_LAMBDA_GRID(nlambda,nw,Bb,Bt,&
!!$       & lambdab_w,lambdac_w,lambdac,dlambdap,lambda,one_o_lambda)  
!!$
!!$  !For each point in the (zeta,theta) grid, determine well and absolute point
!!$  !For each absolute point, determine alpha, lambda and well number
!!$  IF(ALLOCATED(i_l)) DEALLOCATE(i_l)
!!$  ALLOCATE(i_l(npointx),i_w(npointx))
!!$  CALL LABEL_GRIDPOINTS(nalpha,nalphab,nlambda,nw,bottom,connected,&
!!$       & alphap_w,z1,z2,Bb,Bt,lambda,&
!!$       & zetap,zetax,thetax,B_al,npoint,i_l,i_w,i_p)
!!$  ALLOCATE(itemp(npoint))
!!$  itemp=i_l(1:npoint);DEALLOCATE(i_l);ALLOCATE(i_l(npoint));i_l=itemp
!!$  itemp=i_w(1:npoint);DEALLOCATE(i_w);ALLOCATE(i_w(npoint));i_w=itemp
!!$  DEALLOCATE(itemp)
!!$
!!$  IF(ALLOCATED(BI1)) THEN
!!$     DEALLOCATE(nbif,i_p_ap1,i_p_am1,i_p_ap2,i_p_am2,i_p_lp1,i_p_lm1,&
!!$          & dalpha_am1,dalpha_ap1,dalpha_am2,dalpha_ap2,dlambda_lm1,dlambda_lp1,&
!!$          & BI1,BI2,BI3,BI3b,BI3f,BI4,BI5,BI6,BI7,BI8,zlw,zrw,nnz,&
!!$          & i_p_ap1I,i_p_ap1II,i_p_ap2I,i_p_ap2II,i_p_ap2III,i_p_ap2IV,&
!!$          & i_p_am1I,i_p_am1II,i_p_am2I,i_p_am2II,i_p_am2III,i_p_am2IV,&
!!$          & wm1I,wm1II,wm2I,wm2II,wm2III,wm2IV,wp1I,wp1II,wp2I,wp2II,wp2III,wp2IV)
!!$#ifdef MPIandPETSc
!!$     CALL MatDestroy(matCOL,ierr)
!!$     CALL MatDestroy(matVEAf,ierr)
!!$     CALL MatDestroy(matVEAb,ierr)
!!$     IF(TANG_VM) THEN
!!$        CALL MatDestroy(matVMAf,ierr)
!!$        CALL MatDestroy(matVMAb,ierr)
!!$     END IF
!!$#endif
!!$  END IF
!!$  ALLOCATE(nbif(npoint),i_p_ap1(npoint),i_p_am1(npoint),i_p_ap2(npoint),i_p_am2(npoint))
!!$  ALLOCATE(i_p_lm1(npoint),i_p_lp1(nbifx,npoint))
!!$  ALLOCATE(dalpha_am1(npoint),dalpha_ap1(npoint),dalpha_am2(npoint),dalpha_ap2(npoint))
!!$  ALLOCATE(dlambda_lm1(npoint),dlambda_lp1(npoint))
!!$  ALLOCATE(BI1(npoint),BI2(npoint),BI3(npoint),BI3b(npoint),BI3f(npoint),factnu(npoint))
!!$  ALLOCATE(BI4(npoint),BI5(npoint),BI6(npoint),BI7(npoint),BI8(npoint,Nnmp))
!!$  ALLOCATE(zlw(npoint),zrw(npoint),tlw(npoint),trw(npoint)) 
!!$  ALLOCATE(nnz(npoint))
!!$  ALLOCATE(i_p_ap1I(npoint),i_p_ap1II(npoint),i_p_am1I(npoint),i_p_am1II(npoint))
!!$  ALLOCATE(i_p_ap2I(npoint),i_p_ap2II(npoint),i_p_ap2III(npoint),i_p_ap2IV(npoint))
!!$  ALLOCATE(i_p_am2I(npoint),i_p_am2II(npoint),i_p_am2III(npoint),i_p_am2IV(npoint))
!!$  ALLOCATE(wm1I(npoint),wm1II(npoint),wm2I(npoint),wm2II(npoint),wm2III(npoint),wm2IV(npoint))
!!$  ALLOCATE(wp1I(npoint),wp1II(npoint),wp2I(npoint),wp2II(npoint),wp2III(npoint),wp2IV(npoint))
!!$
!!$  !Order alphas in interval [0,2*pi]
!!$  CALL SORT_ALPHA(nalpha,nalphab,alphap,zetax,thetax,thetap,B_al,vds_al,nlambda,i_p)
!!$  DO WHILE(MAXVAL(thetap(:,1))-MINVAL(thetap(:,1)).GT.TWOPI)
!!$     DO ia=1,nalpha
!!$        IF(ABS(thetap(ia,1)-thetap(1,1)).GT.TWOPI) THEN
!!$           thetap(ia,:)=thetap(ia,:)-SIOTA*TWOPI
!!$           alphap(ia)  =alphap(ia)  -SIOTA*TWOPI
!!$        END IF
!!$     END DO
!!$     CALL SORT_ALPHA(nalpha,nalphab,alphap,zetax,thetax,thetap,B_al,vds_al,nlambda,i_p)
!!$  END DO
!!$  IF(DEBUG) THEN
!!$     DO ia=1,nalpha
!!$        DO il=1,nalphab
!!$           WRITE(3000+myrank,'(6(1pe13.5),2I5)') zetap(il),thetap(ia,il),&
!!$                &  zetax(ia,il),thetax(ia,il),alphap(ia),B_al(ia,il),ia,il
!!$        END DO
!!$     END DO
!!$  END IF
!!$  !Find neighbours in lambda
!!$  CALL FIND_LAMBDA_NEIGHBOURS(npoint,nalpha,nalphab,nlambda,nw,nbifx,i_p,&
!!$       & bottom,i_w,i_p_lm1,nbif,i_p_lp1)
!!$  !Correct delta lambda at the bottom
!!$  dlambda_lm1=dlambdap 
!!$  dlambda_lp1=dlambdap
!!$  DO ipoint=1,npoint-1
!!$     IF(i_p_lp1(1,ipoint).EQ.npoint) THEN 
!!$        dlambda_lp1(ipoint)=lambdab_w(i_w(ipoint))-lambda(i_l(ipoint))
!!$     ELSE IF(i_p_lm1(ipoint).EQ.npoint) THEN                        
!!$        dlambda_lm1(ipoint)=lambda(i_l(ipoint))-lambdac
!!$     END IF
!!$  END DO
!!$
!!$  !Calculate coefficients of the drift kinetic equation
!!$  CALL COEFFICIENTS_DKE(npoint,i_w,i_l,nw,&
!!$                 &  z1,t1,B1,hBpp1,vd1,&
!!$                 &  zb,tb,Bb,hBppb,vdb,&
!!$                 &  z2,t2,B2,hBpp2,vd2,&
!!$                 &  nlambda,lambda,zlw,tlw,zrw,trw,& 
!!$                 &  BI1,BI2,BI3,BI4,BI5,BI6,BI7,Nnmp,BI8)
!!$
!!$  !Find neighbours in alpha
!!$  CALL FIND_ALPHA_NEIGHBOURS(npoint,i_p,i_w,i_l,i_p_lm1,i_p_lp1,nbifx,nbif,zlw,zrw,&
!!$  & nw,connected,BI6,nlambda,nalpha,nalphab,alphap,&
!!$  & i_p_am1,i_p_ap1,dalpha_am1,dalpha_ap1,&
!!$  & i_p_am2,i_p_ap2,dalpha_am2,dalpha_ap2,&
!!$  & i_p_am1I,i_p_am1II,i_p_am2I,i_p_am2II,i_p_am2III,i_p_am2IV,&
!!$  & i_p_ap1I,i_p_ap1II,i_p_ap2I,i_p_ap2II,i_p_ap2III,i_p_ap2IV,&
!!$  & wm1I,wm1II,wm2I,wm2II,wm2III,wm2IV,wp1I,wp1II,wp2I,wp2II,wp2III,wp2IV,lambda,&
!!$  & BI7,BI3,BI3b,BI3f,dlambda_lm1,dlambda_lp1)
!!$
!!$  DEALLOCATE(connected,bottom,z1,t1,B1,hBpp1,vd1,zb,tb,Bb,hBppb,vdb,&
!!$       &    z2,t2,B2,hBpp2,vd2,alphap_w,Bt,Btt,lambdab_w,lambdac_w) 
!!$
!!$  !Find non-zero elements of the DKE matrix and initialize
!!$  IF(ALLOCATED(gint)) DEALLOCATE(gint)
!!$  ALLOCATE(rowCOL(npoint),rowVEAf(npoint),rowVEAb(npoint),rowVMAf(npoint),rowVMAb(npoint),&
!!$       & gint(npoint,Nnmp))
!!$#ifdef MPIandPETSc
!!$  DO ipoint=1,npoint
!!$     nnz(ipoint)=1
!!$     IF(ipoint.NE.npoint.AND.i_l(ipoint).NE.0) THEN
!!$        CALL FILL_DKE_ROW(ipoint,npoint,&
!!$             & dalpha_am1(ipoint),dalpha_ap1(ipoint),dalpha_am2(ipoint),dalpha_ap2(ipoint),&
!!$             & i_p_am1(ipoint),i_p_ap1(ipoint),i_p_am2(ipoint),i_p_ap2(ipoint),&
!!$             & i_p_am1I(ipoint),i_p_am1II(ipoint),&
!!$             & i_p_am2I(ipoint),i_p_am2II(ipoint),i_p_am2III(ipoint),i_p_am2IV(ipoint),&
!!$             & i_p_ap1I(ipoint),i_p_ap1II(ipoint),&
!!$             & i_p_ap2I(ipoint),i_p_ap2II(ipoint),i_p_ap2III(ipoint),i_p_ap2IV(ipoint),&
!!$             & wm1I(ipoint),wm1II(ipoint),wm2I(ipoint),wm2II(ipoint),wm2III(ipoint),wm2IV(ipoint),&
!!$             & wp1I(ipoint),wp1II(ipoint),wp2I(ipoint),wp2II(ipoint),wp2III(ipoint),wp2IV(ipoint),&
!!$             & lambda(i_l(ipoint)),dlambda_lm1,dlambda_lp1,i_p_lm1,nbif,nbifx,i_p_lp1,&
!!$             & BI1(ipoint),BI2,BI4(ipoint),BI5(ipoint),BI3(ipoint)/BI7(ipoint),&
!!$             & rowCOL,rowVEAf,rowVEAb,rowVMAf,rowVMAb,nnz(ipoint),.TRUE.)
!!$     END IF
!!$     innz(ipoint-1)=nnz(ipoint)
!!$  END DO
!!$  innz(npoint-1)=1
!!$  CALL INIT_LINEAR_PROBLEM(npoint,nnz,matCOL,matVEAf,matVEAb,matVMAf,matVMAb)
!!$#else
!!$  IF(ALLOCATED(COL)) DEALLOCATE(COL,VEAf,VEAb,VMAf,VMAb)
!!$  ALLOCATE(COL(npoint,npoint),&
!!$       & VEAf(npoint,npoint),VEAb(npoint,npoint),VMAf(npoint,npoint),VMAb(npoint,npoint))
!!$#endif
!!$
!!$  !For each point number, fill one row of the matrix with quantities
!!$  !that depend on the configuration only
!!$  DO ipoint=1,npoint
!!$     rowCOL=0
!!$     rowVEAf=0
!!$     rowVEAb=0
!!$     rowVMAf=0
!!$     rowVMAb=0
!!$     ila =i_l(ipoint)
!!$     IF(ipoint.EQ.npoint.OR.ila.EQ.0) THEN
!!$        rowCOL(ipoint)=1./dlambdap/dlambdap
!!$     ELSE
!!$        CALL FILL_DKE_ROW(ipoint,npoint,&
!!$             & dalpha_am1(ipoint),dalpha_ap1(ipoint),dalpha_am2(ipoint),dalpha_ap2(ipoint),&
!!$             & i_p_am1(ipoint),i_p_ap1(ipoint),i_p_am2(ipoint),i_p_ap2(ipoint),&
!!$             & i_p_am1I(ipoint),i_p_am1II(ipoint),&
!!$             & i_p_am2I(ipoint),i_p_am2II(ipoint),i_p_am2III(ipoint),i_p_am2IV(ipoint),&
!!$             & i_p_ap1I(ipoint),i_p_ap1II(ipoint),&
!!$             & i_p_ap2I(ipoint),i_p_ap2II(ipoint),i_p_ap2III(ipoint),i_p_ap2IV(ipoint),&
!!$             & wm1I(ipoint),wm1II(ipoint),wm2I(ipoint),wm2II(ipoint),wm2III(ipoint),wm2IV(ipoint),&
!!$             & wp1I(ipoint),wp1II(ipoint),wp2I(ipoint),wp2II(ipoint),wp2III(ipoint),wp2IV(ipoint),&
!!$             & lambda(i_l(ipoint)),dlambda_lm1,dlambda_lp1,i_p_lm1,nbif,nbifx,i_p_lp1,&
!!$             & BI1(ipoint),BI2,BI4(ipoint),BI5(ipoint),BI3(ipoint)/BI7(ipoint),&
!!$             & rowCOL,rowVEAf,rowVEAb,rowVMAf,rowVMAb,nnz(ipoint),.FALSE.)
!!$     END IF
!!$
!!$#ifdef MPIandPETSc 
!!$     ipointm1=ipoint-1    
!!$     DO jpoint=1,npoint
!!$        IF((ABS(rowCOL(jpoint)).GT.ZERO).OR.&
!!$             & (ABS(rowVEAf(jpoint)).GT.ZERO).OR. &
!!$             & (ABS(rowVEAb(jpoint)).GT.ZERO).OR. &
!!$             & (TANG_VM.AND.ABS(rowVMAf(jpoint)).GT.ZERO).OR. &
!!$             & (TANG_VM.AND.ABS(rowVMAb(jpoint)).GT.ZERO)) THEN
!!$           jpointm1=jpoint-1
!!$           mat_entry=rowCOL(jpoint)
!!$           CALL MatSetValues(matCOL,oneps,ipointm1,oneps,jpointm1,mat_entry,INSERT_VALUES,ierr)
!!$        END IF
!!$        IF(ABS(rowVEAf(jpoint)).GT.ZERO) THEN
!!$           mat_entry=rowVEAf(jpoint)
!!$           CALL MatSetValues(matVEAf,oneps,ipointm1,oneps,jpointm1,mat_entry,INSERT_VALUES,ierr)
!!$        END IF
!!$        IF(ABS(rowVEAb(jpoint)).GT.ZERO) THEN
!!$           mat_entry=rowVEAb(jpoint)
!!$           CALL MatSetValues(matVEAb,oneps,ipointm1,oneps,jpointm1,mat_entry,INSERT_VALUES,ierr)
!!$        END IF
!!$        IF(TANG_VM.AND.ABS(rowVMAf(jpoint)).GT.ZERO) THEN
!!$           mat_entry=rowVMAf(jpoint)
!!$           CALL MatSetValues(matVMAf,oneps,ipointm1,oneps,jpointm1,mat_entry,INSERT_VALUES,ierr)
!!$        END IF
!!$        IF(TANG_VM.AND.ABS(rowVMAb(jpoint)).GT.ZERO) THEN
!!$           mat_entry=rowVMAb(jpoint)
!!$           CALL MatSetValues(matVMAb,oneps,ipointm1,oneps,jpointm1,mat_entry,INSERT_VALUES,ierr)
!!$        END IF
!!$     END DO
!!$#else
!!$     COL(ipoint,:) =rowCOL
!!$     VEAf(ipoint,:)=rowVEAf
!!$     VEAb(ipoint,:)=rowVEAb
!!$     IF(TANG_VM) THEN
!!$        VMAf(ipoint,:)=rowVMAf
!!$        VMAb(ipoint,:)=rowVMAb
!!$     END IF
!!$#endif
!!$
!!$  END DO
!!$#ifdef MPIandPETSc 
!!$  
!!$  CALL MatAssemblyBegin(matCOL,MAT_FINAL_ASSEMBLY,ierr)
!!$  CALL MatAssemblyEnd(  matCOL,MAT_FINAL_ASSEMBLY,ierr)
!!$  CALL MatAssemblyBegin(matVEAf,MAT_FINAL_ASSEMBLY,ierr)
!!$  CALL MatAssemblyEnd(  matVEAf,MAT_FINAL_ASSEMBLY,ierr)
!!$  CALL MatAssemblyBegin(matVEAb,MAT_FINAL_ASSEMBLY,ierr)
!!$  CALL MatAssemblyEnd(  matVEAb,MAT_FINAL_ASSEMBLY,ierr)
!!$  IF(TANG_VM) THEN
!!$     CALL MatAssemblyBegin(matVMAf,MAT_FINAL_ASSEMBLY,ierr)
!!$     CALL MatAssemblyEnd(  matVMAf,MAT_FINAL_ASSEMBLY,ierr)
!!$     CALL MatAssemblyBegin(matVMAb,MAT_FINAL_ASSEMBLY,ierr)
!!$     CALL MatAssemblyEnd(  matVMAb,MAT_FINAL_ASSEMBLY,ierr)
!!$  END IF
!!$  
!!$#endif
!!$
!!$123 nalphab=nalphab_save 
!!$  theta=theta_save
!!$  IF(DEBUG) THEN
!!$     DO ipoint=1,npoint
!!$        WRITE(3300+myrank,'(I6,1(1pe13.5),2I6)') ipoint,BI3(ipoint)*vdconst(jv)/nu(jv),nal,nlambda
!!$     END DO
!!$  END IF
!!$
!!$  !Use linear combinations of precalculated matrices to fill actual matrix for given values
!!$  !the collisionality, radial electric field, etc
!!$#ifdef MPIandPETSc 
!!$  CALL FILL_MATRIX_PETSC(matCOL,jv,Epsi,matVEAf,matVEAb,matVMAf,matVMAb,ksp)
!!$#else
!!$  ALLOCATE(mat(npoint,npoint))
!!$  CALL FILL_MATRIX(npoint,COL,jv,Epsi,VEAf,VEAb,VMAf,VMAb,mat)
!!$#endif
!!$
!!$  factnu=1
!!$  IF(FLUX_NU) THEN
!!$     DO ipoint=1,npoint-1
!!$        factnu(ipoint)=-lambda(i_l(ipoint))*avB*nu(jv)/(-2*Epsi*BI3(ipoint))
!!$     END DO
!!$  END IF
!!$
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$!Invert linear system
!!$#ifdef MPIandPETSc
!!$  IF(CALC_RHS) THEN
!!$     CALL INVERT_MATRIX_PETSC(nalphab,jv,npoint,matA,BI3,BI7,factnu,phi1c,ksp,gint)
!!$  ELSE
!!$     CALL INVERT_MATRIX_PETSC(nalphab,jv,npoint,matCOL,BI3,BI7,factnu,phi1c,ksp,gint)
!!$  END IF
!!$#else
!!$  IF(FLUX_NU.OR.CALC_RHS) STOP
!!$  CALL INVERT_MATRIX(nalphab,jv,npoint,BI3,BI7,phi1c,mat,gint)
!!$  DEALLOCATE(mat)
!!$#endif 
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$  
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  
!!$  IF(PHI1_READ.AND.IMP1NU) THEN
!!$     gint(:,2)=gint(:,1)*ftrace1nu(jv)*EXP(-mmuoT(jv)/lambda(i_l(:)))*mmuoT(jv)/lambda(i_l(:))
!!$     gint(:,1)=gint(:,1)*ftrace1nu(jv)*EXP(-mmuoT(jv)/lambda(i_l(:)))
!!$     WRITE(iout,*) 'Using mu and lambda'
!!$  END IF
!!$
!!$  !Calculate dn_1dv and D_{11} integrating it the velocity space and taking the flux-surface-average
!!$  CALL INTEGRATE_G(nalpha,nalphab,nlambda,lambda,i_p,npoint,gint, &
!!$       & zetap,thetap,theta(1:nalphab),B_al,vds_al,D11,dn1dv(1:nalphab,1:nalphab),dn1nm)
!!$
!!$  IF(.NOT.PHI1_READ.OR..NOT.IMP1NU) THEN
!!$     D11(1,:)=D11(1,:)*vdconst(jv)
!!$     IF(PHI1_READ) D11     =D11     *weight(jv)
!!$  END IF
!!$
!!$  omega=ABS(Epsi)*psip/v(jv)
!!$  cmul=nu(jv)/v(jv)/2.  
!!$
!!$  !Connection with plateau regime
!!$  IF(FACT_CON.GT.0.AND.cmul_1NU.GT.0.&
!!$       & .AND.cmul.GT.cmul_1NU/FACT_CON.AND.omega.LT.1E-2) D11=D11+D11pla/fdkes(jv)
!!$  
!!$  !Write monoenergetic transport coefficients using DKES normalization
!!$  IF(DEBUG) THEN!.OR.(ONLY_DB.AND..NOT.KNOSOS_STELLOPT)) THEN
!!$     CALL FLUSH(10000+myrank)
!!$     IF(cmul_1NU.GT.0) THEN
!!$        WRITE(10000+myrank,'("4 ",6(1pe13.5),3I5,1pe13.5,I5)') &
!!$             & nu(jv)/v(jv)/2.,Epsi/v(jv)*psip,vmconst(jv)/v(jv),&
!!$             & fdkes(jv)*D11(1,1),weight(jv)/fdkes(jv),weight(jv)/fdkes(jv)*v(jv)*v(jv),&
!!$             & nal,nlambda,nalphab,D11pla
!!$     ELSE
!!$        WRITE(10000+myrank,'("0 ",6(1pe13.5),3I5,1pe13.5,I5)') &
!!$             & nu(jv)/v(jv)/2.,Epsi/v(jv)*psip,vmconst(jv)/v(jv),&
!!$             & fdkes(jv)*D11(1,1),weight(jv)/fdkes(jv),weight(jv)/fdkes(jv)*v(jv)*v(jv),&
!!$             & nal,nlambda,nalphab,D11pla
!!$     END IF
!!$     CALL FLUSH(10000+myrank)
!!$  END IF
!!$  
!!$  IF(PHI1_READ) bnmc0=bnmc0-2*borbic(0,0)*phi1c/vdconst(jv)
!!$
!!$  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE CALC_LOW_COLLISIONALITY_NANL_NEW


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
