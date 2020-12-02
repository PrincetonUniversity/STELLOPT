
!Calculate neoclassical transport at low collisionalities
!Include 1/nu,sqrtnu and superbanana-plateau; do not need splitting B=B_0+delta*B_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALC_LOW_COLLISIONALITY(jv,Epsi,phi1c,Mbbnm,trMnm,D11,&
     & nalphab,zeta,theta,dn1dv,dn1nm)

!--------------------------------------------------------------------------------------------- 
!Calculate monoenergetic transport coefficient D11 and contribution quasineutrality dn1dv and 
!dn1nm at nalphabxnalphab grid in (zeta,theta)  for collisionality cmul=nu(jv)/v(jv) and 
!normalized radial electric field Epsi/v(jv), and in the presence of phi1c
!-----------------------------------------------------------------------------------------------

  USE GLOBAL
  USE KNOSOS_STELLOPT_MOD
  IMPLICIT NONE  
!  Input
  INTEGER jv
  REAL*8 Epsi,phi1c(Nnmp),Mbbnm(Nnmp),trMnm(Nnmp)
  !Output
  INTEGER nalphab
  REAL*8 D11(Nnmp,Nnmp),zeta(nax),theta(nax),dn1dv(nax,nax),dn1nm(Nnmp,Nnmp)
  !Others
  CHARACTER*100 serr
  INTEGER ial,ila
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
     CALL CALC_LOW_COLLISIONALITY_NANL(nal,nlambda,jv,Epsi,phi1c,Mbbnm,trMnm,&
          & D11,nalphab,zeta,theta,dn1dv,dn1nm)
     IF(.NOT.KNOSOS_STELLOPT) WRITE(200+myrank,'(3(1pe13.5)," NaN ",2(1pe13.5)," &
          & NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN")') &
          & nu(jv)/v(jv)/2.,Epsi/v(jv)*psip,vmconst(jv)/v(jv),fdkes(jv)*D11(1,1),fdkes(jv)*D11(1,1)
  ELSE

     D11=0
     CONVERGED=.FALSE.
     ial=1
     nal=32
     DO WHILE (nal.LE.nax)
        ila=1
        nlambda=32
        DO WHILE (nlambda.LE.nlambdax)
           CALL CALC_LOW_COLLISIONALITY_NANL(nal,nlambda,jv,Epsi,phi1c,Mbbnm,trMnm,&
                & D11r(ial,ila),nalphab,zeta,theta,dn1dv,dn1nm)
           IF(ila.GT.1.AND.ABS(D11r(ial,ila)/D11r(ial,ila-1)-1.0).LT.PREC_DQDV) THEN
              EXIT
           END IF
           ila=ila+1
           nlambda=nlambda*2  
        END DO
        D11r(ial,ila+1:100)=D11r(ial,ila)
        IF(ial.GT.1.AND.nlambda.LE.nlambdax.AND.&
             & ABS(D11r(ial,ila)/D11r(ial-1,100)-1.0).LT.PREC_DQDV) THEN
           IF(.NOT.KNOSOS_STELLOPT) WRITE(200+myrank,'(3(1pe13.5)," NaN ",2(1pe13.5),&
                & " NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN")') &
                & nu(jv)/v(jv)/2.,Epsi/v(jv)*psip,vmconst(jv)/v(jv),fdkes(jv)*D11r(ial,ila),fdkes(jv)*D11r(ial,ila)
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
        D11=D11r(ial,ila)
        WRITE(iout,'(" DKE converged with nal=",I6," and nlambda=",I6)') nal,nlambda
        CALL CALC_LOW_COLLISIONALITY_NANL(nal,nlambda,jv,Epsi,phi1c,Mbbnm,trMnm,&
             D11,nalphab,zeta,theta,dn1dv,dn1nm)
     END IF

  END IF

END SUBROUTINE CALC_LOW_COLLISIONALITY


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALC_LOW_COLLISIONALITY_NANL(nal,nlambda,jv,Epsi,phi1c,Mbbnm,trMnm,&
     & D11,nalphab,zeta,theta,dn1dv,dn1nm)

!--------------------------------------------------------------------------------------------- 
!Calculate monoenergetic transport coefficient D11 and contribution quasineutrality dn1dv and 
!dn1nm at nalphabxnalphab grid in (zeta,theta), for collisionality cmul=nu(jv)/v(jv) and normalized 
!radial electric field Epsi/v(jv)m and in the presence of phi1c
!The DKE is solved in a nalxnlambda grid in (alpha,lambda)
!-----------------------------------------------------------------------------------------------

  USE GLOBAL
  USE KNOSOS_STELLOPT_MOD
#ifdef IPPorNIFS
  USE petscsys
  USE petscksp
!  USE petscvec
!  USE petscmat
!  USE petscpc
#endif
  IMPLICIT NONE
!  Input
  INTEGER jv,nlambda,nal
  REAL*8 Epsi,phi1c(Nnmp),Mbbnm(Nnmp),trMnm(Nnmp),zeta(nax),theta(nax)
  !Output
  INTEGER nalphab
  REAL*8 D11(Nnmp,Nnmp),dn1dv(nax,nax),dn1nm(Nnmp,Nnmp)
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
  REAL*8 dlambdap,dalphap(nturn+1)
  REAL*8, SAVE :: lambdac
  !Lambda neighbours
  INTEGER, PARAMETER :: nbifx=200  !maximum number of bifurcations allowed
  INTEGER, SAVE, ALLOCATABLE :: i_p_lm1(:),i_p_lp1(:,:)
  INTEGER ibif,imax,jpoint
  REAL*8 dz,dzmax
  !Alpha neighbours
  REAL*8, ALLOCATABLE, SAVE :: zlw(:),zrw(:)
  REAL*8, ALLOCATABLE       :: tlw(:),trw(:)
  INTEGER ipoint,ila,il,ia
!!$  !Second adiabatic invariant
!!$  LOGICAL, ALLOCATABLE :: readpoint(:)
!!$  INTEGER ja,ia,il
!!$  REAL*8, ALLOCATABLE :: Jnorm(:,:)
  !Drift-kinetic equation
  INTEGER, SAVE, ALLOCATABLE :: nnz(:)
  REAL*8, SAVE, ALLOCATABLE :: BI3f(:),BI3b(:)
  REAL*8, SAVE, ALLOCATABLE :: BI1(:),BI2(:),BI3(:),BI4(:),BI5(:),BI6(:),BI7(:),BI8(:,:)
  REAL*8, SAVE, ALLOCATABLE :: factnu(:),gint(:,:)
  REAL*8 omega,cmul
#ifdef MPIandPETSc
  !Petsc
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscpc.h>
#include <petsc/finclude/petscksp.h>
!#include <petsc/finclude/petscviewer.h>
!  PetscViewer, SAVE :: viewer
  PetscErrorCode ierr
  PetscInt ipointm1,jpointm1
  PetscInt, PARAMETER :: oneps=1
  PetscInt, SAVE :: innz(0:npointx-1)
  KSP, SAVE :: ksp
  Mat, SAVE :: matCOL,matVEAf,matVEAb,matVMAf,matVMAb
  Mat, SAVE :: matDIFb,matDIFf
  Mat matA
  PetscScalar factor
  PetscScalar mat_entry
#else
  REAL*8, SAVE, ALLOCATABLE :: COL(:,:),VEAf(:,:),VEAb(:,:),VMAf(:,:),VMAb(:,:)
  REAL*8, ALLOCATABLE :: mat(:,:)
#endif
  REAL*8, ALLOCATABLE :: rowCOL(:),rowVEAf(:),rowVEAb(:),rowVMAf(:),rowVMAb(:)
  REAL*8, ALLOCATABLE :: rowDIFb(:),rowDIFf(:)
  !Time
  CHARACTER*30, PARAMETER :: routine="CALC_LOW_COLLISIONALITY"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart
!  CHARACTER*30, PARAMETER :: routine2="FILL_GRID"
!  INTEGER, SAVE :: ntotal2=0
!  REAL*8,  SAVE :: ttotal2=0
!  REAL*8,  SAVE :: t02=0
!  REAL*8 tstart2

  Mbbnm=Mbbnm !To be removed
  trMnm=trMnm

  CALL CPU_TIME(tstart)

  WRITE(iout,*) 'Calculating low collisionality regimes'

  IF(CALCULATED_INT) GOTO 123

!  CALL CPU_TIME(tstart2)

  IF(PHI1_READ) bnmc0=bnmc0+2*borbic(0,0)*phi1c/vdconst(jv) 

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
  
  IF(ALLOCATED(lambda)) DEALLOCATE(lambda,i_p)
  ALLOCATE(lambda(nlambda),one_o_lambda(nlambda),i_p(nlambda,nalpha,nalphab))

  !Set global grid in lambda
  CALL CREATE_LAMBDA_GRID(nlambda,nw,Bb,Bt,&
       & lambdab_w,lambdac_w,lambdac,dlambdap,lambda,one_o_lambda) 

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
     DEALLOCATE(nbif,i_p_ap1,i_p_am1,i_p_ap2,i_p_am2,i_p_lp1,i_p_lm1,&
          & dalpha_am1,dalpha_ap1,dalpha_am2,dalpha_ap2,dlambda_lm1,dlambda_lp1,&
          & BI1,BI2,BI3,BI3f,BI3b,BI4,BI5,BI6,BI7,BI8,factnu,zlw,zrw,nnz,&          
          & i_p_ap1I,i_p_ap1II,i_p_ap2I,i_p_ap2II,i_p_ap2III,i_p_ap2IV,&
          & i_p_am1I,i_p_am1II,i_p_am2I,i_p_am2II,i_p_am2III,i_p_am2IV,&
          & wm1I,wm1II,wm2I,wm2II,wm2III,wm2IV,wp1I,wp1II,wp2I,wp2II,wp2III,wp2IV)
#ifdef MPIandPETSc
     CALL MatDestroy(matCOL,ierr)
     CALL MatDestroy(matVEAf,ierr)
     CALL MatDestroy(matVEAb,ierr)
     IF(CALC_DIFF) THEN
        CALL MatDestroy(matDIFf,ierr)
        CALL MatDestroy(matDIFb,ierr)
     END IF
     IF(TANG_VM) THEN
        CALL MatDestroy(matVMAf,ierr)
        CALL MatDestroy(matVMAb,ierr)
     END IF
#endif
  END IF
  ALLOCATE(nbif(npoint),i_p_ap1(npoint),i_p_am1(npoint),i_p_ap2(npoint),i_p_am2(npoint))
  ALLOCATE(i_p_lm1(npoint),i_p_lp1(nbifx,npoint))
  ALLOCATE(dalpha_am1(npoint),dalpha_ap1(npoint),dalpha_am2(npoint),dalpha_ap2(npoint))
  ALLOCATE(dlambda_lm1(npoint),dlambda_lp1(npoint))
  ALLOCATE(BI1(npoint),BI2(npoint),BI3(npoint),BI3b(npoint),BI3f(npoint),factnu(npoint))
  ALLOCATE(BI4(npoint),BI5(npoint),BI6(npoint),BI7(npoint),BI8(npoint,Nnmp))
  ALLOCATE(zlw(npoint),zrw(npoint),tlw(npoint),trw(npoint)) 
  ALLOCATE(nnz(npoint))
  ALLOCATE(i_p_ap1I(npoint),i_p_ap1II(npoint),i_p_am1I(npoint),i_p_am1II(npoint))
  ALLOCATE(i_p_ap2I(npoint),i_p_ap2II(npoint),i_p_ap2III(npoint),i_p_ap2IV(npoint))
  ALLOCATE(i_p_am2I(npoint),i_p_am2II(npoint),i_p_am2III(npoint),i_p_am2IV(npoint))
  ALLOCATE(wm1I(npoint),wm1II(npoint),wm2I(npoint),wm2II(npoint),wm2III(npoint),wm2IV(npoint))
  ALLOCATE(wp1I(npoint),wp1II(npoint),wp2I(npoint),wp2II(npoint),wp2III(npoint),wp2IV(npoint))

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


!!$  DO ipoint=1,nalpha*nalphab
!!$     DO ia=1,nalpha
!!$        DO il=1,nalphab
!!$           IF(j_al(ia,il).EQ.ipoint) THEN
!!$              WRITE(iout,'(2(1pe13.5),1I10)') zetax(ia,il),thetax(ia,il),ipoint
!!$              EXIT
!!$           END IF
!!$        END DO
!!$     END DO
!!$  END DO
!!$  stop
  
  !Find neighbours in lambda
  CALL FIND_LAMBDA_NEIGHBOURS(npoint,nalpha,nalphab,nlambda,nw,nbifx,i_p,&
       & bottom,i_w,i_p_lm1,nbif,i_p_lp1)
 !Correct delta lambda at the bottom
  dlambda_lm1=dlambdap 
  dlambda_lp1=dlambdap
  DO ipoint=2,npoint
     IF(i_p_lp1(1,ipoint).EQ.0) THEN
        dlambda_lp1(ipoint)=lambdab_w(i_w(ipoint))-lambda(i_l(ipoint))
     ELSE IF(i_p_lm1(ipoint).EQ.1) THEN                        
        dlambda_lm1(ipoint)=lambda(i_l(ipoint))-lambdac
     END IF
  END DO

  !Calculate coefficients of the drift kinetic equation
  CALL COEFFICIENTS_DKE(npoint,i_w,i_l,nw,&
                 &  z1,t1,B1,hBpp1,vd1,&
                 &  zb,tb,Bb,hBppb,vdb,&
                 &  z2,t2,B2,hBpp2,vd2,&
                 &  nlambda,lambda,zlw,tlw,zrw,trw,& 
                 &  BI1,BI2,BI3,BI4,BI5,BI6,BI7,Nnmp,BI8)

  IF(KN_STELLOPT(4)) THEN
     IF(ALLOCATED(gint)) DEALLOCATE(gint)
     ALLOCATE(gint(npoint,Nnmp))
     gint=0
     gint(:,1)=2.0*ATAN(BI3/BI4/atorflux)/PI
     gint(:,1)=gint(:,1)*gint(:,1)
     CALL INTEGRATE_G(nalpha,nalphab,nlambda,lambda,i_p,npoint,gint,.TRUE., &
          & zetap,thetap,theta(1:nalphab),B_al,vds_al,D11,dn1dv(1:nalphab,1:nalphab),dn1nm)
     KN_FIC=D11(1,1)
     D11=0
     IF(.NOT.KN_STELLOPT(1)) RETURN
  END IF
!!$  ALLOCATE(Jnorm(nlambda,nalpha))
!!$  ALLOCATE(readpoint(npoint)) 
!!$  Jnorm=0
!!$  DO ila=1,nlambda
!!$     DO ia=1,na
!!$        readpoint=.FALSE.
!!$        DO ja=ia,nalpha,na
!!$           DO il=1,nalphab
!!$              ipoint=i_p(ila,ja,il)
!!$              IF(readpoint(ipoint)) CYCLE
!!$              Jnorm(ila,ia)=Jnorm(ila,ia)+BI6(ipoint)
!!$              readpoint(ipoint)=.TRUE.
!!$           END DO
!!$        END DO
!!$     END DO
!!$  END DO

    !Put longest well first
  DO ipoint=2,npoint
     imax=0
     dzmax=0
     DO ibif=1,nbif(ipoint)
        jpoint=i_p_lp1(ibif,ipoint)
        dz=zrw(jpoint)-zlw(jpoint)
        IF(dz.GT.dzmax) THEN
           imax=ibif
           dzmax=dz
        END IF
     END DO
     IF(imax.GT.1) THEN
        jpoint=i_p_lp1(1,ipoint)
        i_p_lp1(1,ipoint)=i_p_lp1(imax,ipoint)
        i_p_lp1(imax,ipoint)=jpoint
     END IF
  END DO

  !Find neighbours in alpha
  CALL FIND_ALPHA_NEIGHBOURS(npoint,i_p,i_w,i_l,i_p_lm1,i_p_lp1,nbifx,nbif,zlw,zrw,&
       & nw,connected,BI6,nlambda,nalpha,nalphab,alphap,&
       & i_p_am1,i_p_ap1,dalpha_am1,dalpha_ap1,&
       & i_p_am2,i_p_ap2,dalpha_am2,dalpha_ap2,&
       & i_p_am1I,i_p_am1II,i_p_am2I,i_p_am2II,i_p_am2III,i_p_am2IV,&
       & i_p_ap1I,i_p_ap1II,i_p_ap2I,i_p_ap2II,i_p_ap2III,i_p_ap2IV,&
       & wm1I,wm1II,wm2I,wm2II,wm2III,wm2IV,wp1I,wp1II,wp2I,wp2II,wp2III,wp2IV,lambda,&
       & BI7,BI3,BI3b,BI3f,dlambda_lm1,dlambda_lp1)
  
  DEALLOCATE(connected,bottom,z1,t1,B1,hBpp1,vd1,zb,tb,Bb,hBppb,vdb,&
       &    z2,t2,B2,hBpp2,vd2,alphap_w,Bt,Btt,lambdab_w,lambdac_w) 

  !Find non-zero elements of the DKE matrix and initialize
  IF(ALLOCATED(gint)) DEALLOCATE(gint)
  ALLOCATE(rowCOL(npoint),rowVEAf(npoint),rowVEAb(npoint),rowVMAf(npoint),rowVMAb(npoint),&
       & gint(npoint,Nnmp))
  ALLOCATE(rowDIFb(npoint),rowDIFf(npoint))
#ifdef MPIandPETSc
  DO ipoint=1,npoint
     nnz(ipoint)=1
     IF(ipoint.NE.1.AND.(i_p_lm1(ipoint).NE.1.OR.i_p_lp1(1,ipoint).NE.0)) THEN
        CALL FILL_DKE_ROW(ipoint,npoint,&
             & dalpha_am1(ipoint),dalpha_ap1(ipoint),dalpha_am2(ipoint),dalpha_ap2(ipoint),&
             & i_p_am1(ipoint),i_p_ap1(ipoint),i_p_am2(ipoint),i_p_ap2(ipoint),&
             & i_p_am1I(ipoint),i_p_am1II(ipoint),&
             & i_p_am2I(ipoint),i_p_am2II(ipoint),i_p_am2III(ipoint),i_p_am2IV(ipoint),&
             & i_p_ap1I(ipoint),i_p_ap1II(ipoint),&
             & i_p_ap2I(ipoint),i_p_ap2II(ipoint),i_p_ap2III(ipoint),i_p_ap2IV(ipoint),&
             & wm1I(ipoint),wm1II(ipoint),wm2I(ipoint),wm2II(ipoint),wm2III(ipoint),wm2IV(ipoint),&
             & wp1I(ipoint),wp1II(ipoint),wp2I(ipoint),wp2II(ipoint),wp2III(ipoint),wp2IV(ipoint),&
             & lambda(i_l(ipoint)),dlambda_lm1,dlambda_lp1,i_p_lm1,nbif,nbifx,i_p_lp1,&
             & BI1(ipoint),BI2,BI4(ipoint),BI5(ipoint),BI3(ipoint)/BI7(ipoint),&
             & rowCOL,rowVEAf,rowVEAb,rowVMAf,rowVMAb,nnz(ipoint),.TRUE.,&
             & rowDIFf,rowDIFb)
     END IF
     innz(ipoint-1)=nnz(ipoint)
  END DO
  innz(0)=1
  CALL INIT_LINEAR_PROBLEM(npoint,nnz,matCOL,matVEAf,matVEAb,matVMAf,matVMAb,matDIFf,matDIFb)
#else
  IF(ALLOCATED(COL)) DEALLOCATE(COL,VEAf,VEAb,VMAf,VMAb)
  ALLOCATE(COL(npoint,npoint),&
       & VEAf(npoint,npoint),VEAb(npoint,npoint),VMAf(npoint,npoint),VMAb(npoint,npoint))
#endif

  !For each point number, fill one row of the matrix with quantities
  !that depend on the configuration only
  DO ipoint=1,npoint
     rowCOL=0
     rowVEAf=0
     rowVEAb=0
     rowVMAf=0
     rowVMAb=0
     IF(ipoint.EQ.1.OR.(i_p_lm1(ipoint).EQ.1.AND.i_p_lp1(1,ipoint).EQ.0)) THEN
        rowCOL(ipoint)=1.
     ELSE
        CALL FILL_DKE_ROW(ipoint,npoint,&
             & dalpha_am1(ipoint),dalpha_ap1(ipoint),dalpha_am2(ipoint),dalpha_ap2(ipoint),&
             & i_p_am1(ipoint),i_p_ap1(ipoint),i_p_am2(ipoint),i_p_ap2(ipoint),&
             & i_p_am1I(ipoint),i_p_am1II(ipoint),&
             & i_p_am2I(ipoint),i_p_am2II(ipoint),i_p_am2III(ipoint),i_p_am2IV(ipoint),&
             & i_p_ap1I(ipoint),i_p_ap1II(ipoint),&
             & i_p_ap2I(ipoint),i_p_ap2II(ipoint),i_p_ap2III(ipoint),i_p_ap2IV(ipoint),&
             & wm1I(ipoint),wm1II(ipoint),wm2I(ipoint),wm2II(ipoint),wm2III(ipoint),wm2IV(ipoint),&
             & wp1I(ipoint),wp1II(ipoint),wp2I(ipoint),wp2II(ipoint),wp2III(ipoint),wp2IV(ipoint),&
             & lambda(i_l(ipoint)),dlambda_lm1,dlambda_lp1,i_p_lm1,nbif,nbifx,i_p_lp1,&
             & BI1(ipoint),BI2,BI4(ipoint),BI5(ipoint),BI3(ipoint)/BI7(ipoint),&
             & rowCOL,rowVEAf,rowVEAb,rowVMAf,rowVMAb,nnz(ipoint),.FALSE.,&
             & rowDIFf,rowDIFb)
     END IF
#ifdef MPIandPETSc 
     ipointm1=ipoint-1    
     DO jpoint=1,npoint
        IF((ABS(rowCOL(jpoint)).GT.ZERO).OR.&
             & (ABS(rowVEAf(jpoint)).GT.ZERO).OR. &
             & (ABS(rowVEAb(jpoint)).GT.ZERO).OR. &
             & (TANG_VM.AND.ABS(rowVMAf(jpoint)).GT.ZERO).OR. &
             & (TANG_VM.AND.ABS(rowVMAb(jpoint)).GT.ZERO)) THEN
           jpointm1=jpoint-1
           mat_entry=rowCOL(jpoint)
           CALL MatSetValues(matCOL,oneps,ipointm1,oneps,jpointm1,mat_entry,INSERT_VALUES,ierr)
        END IF
        IF(ABS(rowVEAf(jpoint)).GT.ZERO) THEN
           mat_entry=rowVEAf(jpoint)
           CALL MatSetValues(matVEAf,oneps,ipointm1,oneps,jpointm1,mat_entry,INSERT_VALUES,ierr)
           IF(CALC_DIFF) THEN
              mat_entry=rowDIFf(jpoint)
              CALL MatSetValues(matDIFf,oneps,ipointm1,oneps,jpointm1,mat_entry,INSERT_VALUES,ierr)
           END IF
        END IF
        IF(ABS(rowVEAb(jpoint)).GT.ZERO) THEN
           mat_entry=rowVEAb(jpoint)
           CALL MatSetValues(matVEAb,oneps,ipointm1,oneps,jpointm1,mat_entry,INSERT_VALUES,ierr)
           IF(CALC_DIFF) THEN
              mat_entry=rowDIFb(jpoint)
              CALL MatSetValues(matDIFb,oneps,ipointm1,oneps,jpointm1,mat_entry,INSERT_VALUES,ierr)
           END IF
        END IF
        IF(TANG_VM.AND.ABS(rowVMAf(jpoint)).GT.ZERO) THEN
           mat_entry=rowVMAf(jpoint)
           CALL MatSetValues(matVMAf,oneps,ipointm1,oneps,jpointm1,mat_entry,INSERT_VALUES,ierr)
        END IF
        IF(TANG_VM.AND.ABS(rowVMAb(jpoint)).GT.ZERO) THEN
           mat_entry=rowVMAb(jpoint)
           CALL MatSetValues(matVMAb,oneps,ipointm1,oneps,jpointm1,mat_entry,INSERT_VALUES,ierr)
        END IF
     END DO
#else
     COL(ipoint,:) =rowCOL
     VEAf(ipoint,:)=rowVEAf
     VEAb(ipoint,:)=rowVEAb
     IF(TANG_VM) THEN
        VMAf(ipoint,:)=rowVMAf
        VMAb(ipoint,:)=rowVMAb
     END IF
#endif

  END DO
  IF(DEBUG) CALL FLUSH(3200+myrank)

#ifdef MPIandPETSc 
  
  CALL MatAssemblyBegin(matCOL,MAT_FINAL_ASSEMBLY,ierr)
  CALL MatAssemblyEnd(  matCOL,MAT_FINAL_ASSEMBLY,ierr)
  CALL MatAssemblyBegin(matVEAf,MAT_FINAL_ASSEMBLY,ierr)
  CALL MatAssemblyEnd(  matVEAf,MAT_FINAL_ASSEMBLY,ierr)
  CALL MatAssemblyBegin(matVEAb,MAT_FINAL_ASSEMBLY,ierr)
  CALL MatAssemblyEnd(  matVEAb,MAT_FINAL_ASSEMBLY,ierr)
  IF(CALC_DIFF) THEN
     CALL MatAssemblyBegin(matDIFf,MAT_FINAL_ASSEMBLY,ierr)
     CALL MatAssemblyEnd(  matDIFf,MAT_FINAL_ASSEMBLY,ierr)
     CALL MatAssemblyBegin(matDIFb,MAT_FINAL_ASSEMBLY,ierr)
     CALL MatAssemblyEnd(  matDIFb,MAT_FINAL_ASSEMBLY,ierr)
  END IF
  IF(TANG_VM) THEN
     CALL MatAssemblyBegin(matVMAf,MAT_FINAL_ASSEMBLY,ierr)
     CALL MatAssemblyEnd(  matVMAf,MAT_FINAL_ASSEMBLY,ierr)
     CALL MatAssemblyBegin(matVMAb,MAT_FINAL_ASSEMBLY,ierr)
     CALL MatAssemblyEnd(  matVMAb,MAT_FINAL_ASSEMBLY,ierr)
  END IF
  
#endif

  IF(DEBUG) THEN
     DO ipoint=1,npoint
        WRITE(3300+myrank,'(2I6,3(1pe13.5))') &
             & ipoint,i_l(ipoint),BI3(ipoint),BI3b(ipoint),BI3f(ipoint)
     END DO
     CALL FLUSH(3300+myrank)
  END IF

  !Set source to zero for isolated particles, if they exist
  DO ipoint=2,npoint
     IF(i_p_lm1(ipoint).EQ.1.AND.i_p_lp1(1,ipoint).EQ.0)  THEN
        BI3(ipoint)=0
        BI3f(ipoint)=0
        BI3b(ipoint)=0
        BI8(ipoint,:)=0
        WRITE(iout,"(A7,I5,A9,I3,A9)") "Orbit ",ipoint, " at i_la=",i_l(ipoint)," isolated"
     END IF
  END DO

123 nalphab=nalphab_save 
  theta=theta_save

  !Set source to zero for passing particles
  BI3(1)=0
  BI8(1,:)=0
  IF(RE_SOURCE) THEN
     IF(sgnB*Epsi.LT.0) THEN
        BI3=BI3b
     ELSE
        BI3=BI3f
     END IF
  END IF
  !Use linear combinations of precalculated matrices to fill actual matrix for given values
  !the collisionality, radial electric field, etc
#ifdef MPIandPETSc 
  CALL FILL_MATRIX_PETSC(matCOL,jv,Epsi,matVEAf,matVEAb,matVMAf,matVMAb,ksp)

  IF(CALC_DG) THEN
     CALL MatDuplicate(matCOL,MAT_DO_NOT_COPY_VALUES,matA,ierr)
!     factor=0
!     CALL MatAXPY(matA,factor,matCOL,SUBSET_NONZERO_PATTERN,IERR)
     IF (CALC_COL.OR.FLUX_NU) THEN
        factor=1.0
        CALL MatAXPY(matA,factor,matCOL,SUBSET_NONZERO_PATTERN,IERR)
     ELSE IF(CALC_DIFF) THEN
        factor=sgnB*Epsi/nu(jv)
        IF(sgnB*Epsi.LT.0) THEN
           CALL MatAXPY(matA,factor,matDIFb,SUBSET_NONZERO_PATTERN,IERR)
        ELSE IF(sgnB*Epsi.GT.0) THEN
           CALL MatAXPY(matA,factor,matDIFf,SUBSET_NONZERO_PATTERN,IERR)
        END IF
     ELSE
        IF(CALC_DA) THEN
           factor=1.0
        ELSE IF(CALC_RHS) THEN
           factor=sgnB*Epsi/nu(jv)
           IF(sgnB*Epsi.LT.0) THEN
              CALL MatAXPY(matA,factor,matVEAb,SUBSET_NONZERO_PATTERN,IERR)
           ELSE IF(sgnB*Epsi.GT.0) THEN
              CALL MatAXPY(matA,factor,matVEAf,SUBSET_NONZERO_PATTERN,IERR)
           END IF
        END IF
     END IF
  END IF
#else
  ALLOCATE(mat(npoint,npoint))
  CALL FILL_MATRIX(npoint,COL,jv,Epsi,VEAf,VEAb,VMAf,VMAb,mat)
#endif
  factnu=1
  IF(FLUX_NU) THEN
     DO ipoint=2,npoint
        factnu(ipoint)=-lambda(i_l(ipoint))*avB*nu(jv)/(-2*Epsi*BI3(ipoint))
     END DO
  END IF
  !stop
!Invert linear system
#ifdef MPIandPETSc
  CALL INVERT_MATRIX_PETSC(nalphab,jv,npoint,matA,BI3,BI8,factnu,phi1c,ksp,gint)
#else
  IF(CALC_DG) STOP
  CALL INVERT_MATRIX(nalphab,jv,npoint,BI3,BI8,phi1c,mat,gint)
  DEALLOCATE(mat)
#endif  

  IF(DEBUG.AND..NOT.KNOSOS_STELLOPT) THEN
     DO ila=1,nlambda
        DO ia=1,nalpha
           DO il=1,nalphab
           ipoint=i_p(ila,ia,il)
           IF(zetap(il).GT.5.6451E-01.OR.zetap(il).LT.5.6449E-01) CYCLE
           IF(thetap(ia,il).LT.-1.357E+00.OR.thetap(ia,il).GT.-1.355E+00) CYCLE
            IF(ipoint.EQ.0) CYCLE
           IF(.NOT.CALC_DG) THEN
              WRITE(iout,'("size ",3(I5),10(1pe23.15))') ipoint,ia,ila,nu(jv)/v(jv)/2.,Epsi/v(jv)*psip,&
                   & lambda(ila)-lambdac,gint(ipoint,1),(lambda(ila)-lambdac)*avB*0.5*vdconst(jv)/Epsi,zetap(il),&
                   & thetap(ia,il),gint(i_p_lm1(ipoint),1),gint(i_p_lp1(1,ipoint),1)

           ELSE IF(CALC_DA) THEN
              WRITE(iout,'("size ",3(I5),10(1pe13.5))') ipoint,ia,ila,nu(jv)/v(jv)/2.,Epsi/v(jv)*psip,&
                   & lambda(ila)-lambdac,gint(ipoint,1)*nu(jv)/vdconst(jv)/BI5(ipoint),BI3(ipoint),BI6(ipoint)
           ELSE
              WRITE(iout,'("size ",3(I5),10(1pe13.5))') ipoint,ia,ila,nu(jv)/v(jv)/2.,Epsi/v(jv)*psip,&
                   & lambda(ila)-lambdac,gint(ipoint,1)*nu(jv)/vdconst(jv),BI3(ipoint),BI6(ipoint)
           END IF
           END DO
        END DO
     END DO
  END IF
!!$  

  IF(PHI1_READ.AND.IMP1NU) THEN
     gint(:,2)=gint(:,1)*ftrace1nu(jv)*EXP(-mmuoT(jv)/lambda(i_l(:)))*mmuoT(jv)/lambda(i_l(:))
     gint(:,1)=gint(:,1)*ftrace1nu(jv)*EXP(-mmuoT(jv)/lambda(i_l(:)))
     WRITE(iout,*) 'Using mu and lambda'
  END IF

  !Calculate dn_1dv and D_{11} integrating it the velocity space and taking the flux-surface-average
  IF(GEN_FLAG(4)) THEN
     CALL INTEGRATE_G(nalpha,nalphab,nlambda,lambda,i_p,npoint,gint,.FALSE., &
          & zetap,thetap,theta(1:nalphab),B_al,vds_al,D11,dn1dv(1:nalphab,1:nalphab),dn1nm)
  ELSE
     CALL INTEGRATE_G_NEW(nalpha,nalphab,nlambda,lambda,i_p,npoint,gint(:,1),.FALSE., &
          & thetap,B_al,vds_al(1,:,:),D11(1,1))
  END IF
     
  IF(.NOT.PHI1_READ.OR..NOT.IMP1NU) THEN
     D11(1,:)=D11(1,:)*vdconst(jv)
     IF(PHI1_READ) D11     =D11     *weight(jv)
  END IF

!  ELSE
!     gint=gint*weight(jv)
!  END IF
  omega=ABS(Epsi)*psip/v(jv)
  cmul=nu(jv)/v(jv)/2.  

  !Connection with plateau regime
  IF(FACT_CON.GT.0.AND.cmul_1NU.GT.0.&
       & .AND.cmul.GT.cmul_1NU/FACT_CON.AND.omega.LT.1E-2) D11=D11+D11pla/fdkes(jv)
  
  !Write monoenergetic transport coefficients using DKES normalization
  IF(DEBUG) THEN!.OR.(ONLY_DB.AND..NOT.KNOSOS_STELLOPT)) THEN
     CALL FLUSH(10000+myrank)
     IF(cmul_1NU.GT.0) THEN
        WRITE(10000+myrank,'("4 ",6(1pe13.5),3I5,1pe13.5,I5)') &
             & nu(jv)/v(jv)/2.,Epsi/v(jv)*psip,vmconst(jv)/v(jv),&
             & fdkes(jv)*D11(1,1),weight(jv)/fdkes(jv),weight(jv)/fdkes(jv)*v(jv)*v(jv),&
             & nal,nlambda,nalphab,D11pla
     ELSE
        WRITE(10000+myrank,'("0 ",6(1pe13.5),3I5,1pe13.5,I5)') &
             & nu(jv)/v(jv)/2.,Epsi/v(jv)*psip,vmconst(jv)/v(jv),&
             & fdkes(jv)*D11(1,1),weight(jv)/fdkes(jv),weight(jv)/fdkes(jv)*v(jv)*v(jv),&
             & nal,nlambda,nalphab,D11pla
     END IF
     CALL FLUSH(10000+myrank)
  END IF
  
  IF(PHI1_READ) bnmc0=bnmc0-2*borbic(0,0)*phi1c/vdconst(jv)

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE CALC_LOW_COLLISIONALITY_NANL


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CHARACTERIZE_WELLS(nal,na,nalpha,nw,&
                     & z1,t1,B1,hBpp1,vd1, &
                     & zb,tb,Bb,hBppb,vdb, &
                     & z2,t2,B2,hBpp2,vd2, &
                     & Bt,Btt,alphap_w,dalphap,bottom,connected,offset)

!-----------------------------------------------------------------------------------------------
!Find and characterize nw wells (including matched regions) characterized by the Boozer
!toroidal and poloidal angles z and t, and alpha, of its top points, 1 and 2, and bottom
!b, as well as the value of the half-second derivative of the magnetic field strength along
!the magnetic field line hBpp, the value of the magnetic drift vd, alphap_w, their top according
!to two definitions Bt and Btt, whether they are the the bottom.
!Determine a map of connections, connected, and offset TODO
!It exits when nalpha.GT.nal
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL  
  IMPLICIT NONE
  !Input
  INTEGER nal
  !Output
  LOGICAL matched(nwx),connected(nwx,nwx),bottom(nwx)
  INTEGER na,nalpha,nw
  REAL*8 z1(nwx),t1(nwx),B1(nwx),hBpp1(nwx),vd1(nqv,nwx)
  REAL*8 zb(nwx),tb(nwx),Bb(nwx),hBppb(nwx),vdb(nqv,nwx) 
  REAL*8 z2(nwx),t2(nwx),B2(nwx),hBpp2(nwx),vd2(nqv,nwx) 
  REAL*8 Bt(nwx),Btt(nwx),alphap_w(nwx),dalphap,offset
  !Others
  CHARACTER*100 serr
  INTEGER namin,namax,iw,jw,kw,lw,nwperiods,iturn,flag,fath(nwx),nfound,nfoundt,nw0,nw1,nw2,nwmax,nalphas,mturn
  REAL*8 maxBt,maxBb,da,dat
  REAL*8, PARAMETER :: diota=1E-3
  !Time
  CHARACTER*30, PARAMETER :: routine="CHARACTERIZE_WELLS"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)
  
  IF(ONE_ALPHA) THEN
     namax=INT(nal*1.5)
     DO na=1,namax
        dat=aiota*na-INT(aiota*na)
        IF(dat.GT.0.5) dat=1.0-dat
        IF(dat.LT.diota) THEN
!           namax=na
           EXIT
        END IF
     END DO
     da=1E3
     namin=1
     IF(namax.GT.nal) namin=nal
     DO na=namin,namax        
        dat=MOD(TWOPI*na*aiota/nzperiod,TWOPI)
        IF(dat.GT.PI) dat=TWOPI-dat
        IF(dat.LT.da) THEN
           da=dat
           nalphas=na
        END IF
     END DO
     mturn=INT(nalphas*aiota/nzperiod+1)
     na=1
  ELSE
     nal=nal/NTURN
     na=3  !Follow na field lines along several periods and collapses them into one period;
     !need na>=3 when using i+/-1 and i+/-2 for the derivatives at i;
     IF(nalpha.GT.nal) THEN
        serr="increase nal"
        !     CALL END_ALL(serr,.FALSE.)
     END IF
     mturn=NTURN
  END IF
  nalpha=0

  !Initialize variables
  connected=.FALSE.
  matched  =.FALSE.
  bottom   =.FALSE.
  iturn=0
  nw0=1
  nw2=0
  nwmax=0
  nfoundt=0
  fath=0
  alphap_w=-1000.0 
  !Maximum possible number of wells determined by size of arrays
  !  DO WHILE((nalpha/NTURN.LT.nal.AND.(nwmax.LE.nwx.OR.(NTV.AND.nwmax.LE.nwx))).OR.iturn.NE.0)
  DO WHILE((nalpha/NTURN.LT.nal.AND.nwmax.LE.nwx).OR.iturn.NE.0)
     iturn=iturn+1
     !Find wells along the field lines
     IF(iturn.EQ.1) THEN        
        nw0=nw2+1
        !Characterize wells with extreme values
        IF(na.GT.1) THEN
           offset=-iota*PI/nzperiod
        ELSE
           dat=MOD(TWOPI*nalphas*aiota/nzperiod,TWOPI)
           IF(dat.LT.PI) THEN
              offset=-siota*dat/2.
           ELSE
              dat=TWOPI-dat
              offset=siota*dat/2.
           END IF
        END IF
        CALL FIND_WELLS(na,mturn,nw0,z1,t1,B1,hBpp1,vd1, &
             &                 zb,tb,Bb,hBppb,vdb, &
             &                 z2,t2,B2,hBpp2,vd2, &
             &           nfound,alphap_w,offset)
        nw1=nw0
        nw2=nw0+nfound-1
        DO iw=nw1,nw2
           Bt(iw) =MIN(B1(iw),B2(iw))!smallest maximum
           Btt(iw)=MAX(B1(iw),B2(iw))!largest  maximum
           bottom(iw)=.TRUE.
        END DO
        nfoundt=nfoundt+nfound !to nfoundt local mimima, by matching
        nwmax=4*nfoundt        !we may obtain up to 2*nfoundt regions
                               !exit loop if nwmax.GT.NWX
     !Match wells  
     ELSE
        kw=nw2+1
        nw1=kw 
        !Try to match well iw...
        DO iw=nw0,nw2 
!           IF(NTV) CYCLE
           IF(matched(iw)) CYCLE 
           !...with well jw
           DO jw=nw0,nw2
              IF(jw.EQ.iw.OR.matched(jw)) CYCLE 
              CALL MATCH_WELLS(kw, & 
                   & z1(iw),t1(iw),B1(iw),hBpp1(iw),vd1(:,iw), &
                   & zb(iw),tb(iw),Bb(iw),hBppb(iw),vdb(:,iw), &
                   & z2(iw),t2(iw),B2(iw),hBpp2(iw),vd2(:,iw), &
                   & z1(jw),t1(jw),B1(jw),hBpp1(jw),vd1(:,jw), &
                   & zb(jw),tb(jw),Bb(jw),hBppb(jw),vdb(:,jw), &
                   & z2(jw),t2(jw),B2(jw),hBpp2(jw),vd2(:,jw), &
                   & z1(kw),t1(kw),B1(kw),hBpp1(kw),vd1(:,kw), &
                   & zb(kw),tb(kw),Bb(kw),hBppb(kw),vdb(:,kw), &
                   & z2(kw),t2(kw),B2(kw),hBpp2(kw),vd2(:,kw), &
                   & 1,flag)
              !If they match, set limits, matrix of relations....
              IF(flag.EQ.1) THEN 
                 Bt(kw) =MIN(B1(kw),B2(kw))
                 Btt(kw)=MAX(B1(kw),B2(kw))
                 nwperiods=INT(ABS(z2(kw)-z1(kw))*nzperiod/TWOPI/1.5)+1
                 IF(DEBUG) WRITE(2900+myrank,'("    Region ",I6," (from ",I6," and" ,I6,   &
                      & "contained in",I6," periods")') kw,iw,jw,nwperiods
                 matched(iw)=.TRUE. 
                 matched(jw)=.TRUE. 
                 fath(iw)=kw
                 fath(jw)=kw
                 alphap_w(kw)=alphap_w(iw)
                 connected(kw,iw)=.TRUE. !needed for determing that well iw
                 connected(kw,jw)=.TRUE. !cover alphas that jw and kw covered
                 DO lw=1,nw            
                    IF(connected(iw,lw).OR.connected(jw,lw)) connected(kw,lw)=.TRUE.
                 END DO
                 kw=kw+1
              END IF
           END DO
        END DO
        nw2=kw-1 
     END IF
     nw=nw2
     
     !If there were matches, continue matching for the same alphas
!     IF(.NOT.NTV.AND.nw2.GE.nw1) THEN
     IF(nw2.GE.nw1) THEN

        IF(DEBUG) THEN
           maxBt=-1E3
           DO iw=1,nw
              IF(.NOT.matched(iw).AND.Bt(iw).GT.maxBt) maxBt=Bt(iw)
           END DO
           WRITE(2900+myrank,'(" Labelled regions ",I6," to ",I6," out of a maximum of ", &           
                & I6,":",1pe13.5," < 1/lambda <",1pe13.5)')                     &
                & nw1,nw2,nwmax,MINVAL(Bb),maxBt
        END IF
     !If there were no matches, turn to new alphas
     ELSE
        WRITE(iout,'(" Calculating a total of",I6," regions found following",I6," lines")') nw,na

        IF(DEBUG) THEN
           WRITE(2900+myrank,'(" Calculating a total of",I6," regions out of a maximum of",I6)') &
                & nw,nwmax
           maxBb=MAX(MAXVAL(B1),MAXVAL(B2))
           WRITE(2900+myrank,'(" Covering most of region   ",1pe13.5," < 1/lambda <",1pe13.5)') &
                & MINVAL(Bb),maxBt
           WRITE(2900+myrank,'(" Neglecting at least region",1pe13.5," < 1/lambda <",1pe13.5,   &
                & " (",1pe9.2," % of the lambda space)")')                            &
                & maxBb,maxBt,(1./maxBt-1./maxBb)/(1./MINVAL(Bb(nw0:nw))-1./maxBb)*100
        END IF
        IF(nw.GT.nwx) THEN
           serr="nw<nwx"
           CALL END_ALL(serr,.FALSE.)
        END IF
        !Double the number of alphas and start over
        iturn=0  
!        IF(NTV) THEN
!           dalphap=siota*TWOPI/na
!           nalpha=na
!        ELSE
        dalphap=iota*TWOPI/nzperiod/na  !dalphap has sign
        na=na*2
        IF(na.EQ.2) THEN
           nalpha=nalphas
           EXIT
        ELSE
           nalpha=INT(TWOPI/ABS(dalphap)+1)*NTURN
        END IF
     ENDIF
  END DO
  na=na/2
  IF(EXTRA_ALPHA) nalpha=na*INT(nzperiod/aiota+1)*NTURN

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE CHARACTERIZE_WELLS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CREATE_ANGULAR_GRID(na,nalpha,nalphab,alphap,dalphap,offset,&
     & zetap,thetap,zetax,thetax,B_al,vds_al,j_al)

!-----------------------------------------------------------------------------------------------
!From na independent field lines, create:
!- a nalphaxnalpha (zetap,thetap) grid aligned with the magnetic field lines 
! and labelled with alphap (dalphap is the spacing);
!- a nax(nalphaxnalpha/na) grid (zetax,thetap) extended over several field periods.
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL  
  IMPLICIT NONE
  !Input
  INTEGER na,nalpha,nalphab
  REAL*8 alphap(nalpha),dalphap(nturn+1),offset
  !Output
  INTEGER j_al(nalpha,nalphab)
  REAL*8 zetap(nalphab),thetap(nalpha,nalphab),zetax(nalpha,nalphab),thetax(nalpha,nalphab)
  REAL*8 B_al(nalpha,nalphab),vds_al(Nnmp,nalpha,nalphab)
  !Others
  INTEGER ja,ia,il,iturn
  REAL*8 dzeta,dummy
!  REAL*8 zetal(nalphab*nalpha),thetal(nalphab*nalpha)
!  INTEGER imin,np1,jp1,intnumn(nax),nmj
!  REAL*8 offsetb,doff,zetat(nalpha),thetat(nalphab),temp(nalphab,nalphab)
  !Time
  CHARACTER*30, PARAMETER :: routine="CREATE_ANGULAR_GRID"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  WRITE(iout,'(" Angular grid has size ",I4," in alpha and ",I4," in l. &
       & Maximum helicity in B is ",I4)') nalpha,nalphab,INT(MAX(MAXVAL(np),MAXVAL(mp)))

!  IF(NTV) THEN
!     offtheta=0
!     offzeta  =PI/nzperiod
!     dzeta=TWOPI/nalphab
!     DO ia=1,nalpha
!        DO il=1,nalphab
!           zetap(il) =offzeta            +(il-1)*dzeta
!           IF(iota.GT.0) THEN
!              thetap(ia,il)=offtheta+ia*dalphap+(il-1)*dzeta*iota
!           ELSE
!              thetap(ia,il)=TWOPI+offtheta+ia*dalphap+(il-1)*dzeta*iota
!           END IF
!           zetax(ia,il) =zetap(il)
!           CALL FILL_BNODE(zetap(il),thetap(ia,il),dummy,B_al(ia,il),vds_al(ia,il),dummy,.FALSE.)
!           IF(DEBUG) WRITE(3000+myrank,'(5(1pe13.5),2I5)') zetax(ia,il),&
!                &  zetap(il),thetap(ia,il),B_al(ia,il),ia,il
!        END DO
!     END DO

  !  ELSE

  dzeta=TWOPI/nzperiod/nalphab

  DO il=1,nalphab
     zetap(il)=(il-1)*dzeta

!     IF(NTV) zetap(il)=zetap(il)+PI/nzperiod
  END DO

!!$  DO il=1,nalphab
!!$!     IF(NTV) zetap(il)=zetap(il)+PI/nzperiod
!!$     DO ia=1,nalpha !note that dalphap can be negative, if iota is, but...
!!$        thetap(ia,il)=offset+ia*dalphap+(il-1)*iota*dzeta
!!$!        IF(NTV) thetap(ia,il)=thetap(ia,il)+PI 
!!$        zetax(ia,il) =zetap(il)+INT((ia-1)/na)*TWOPI/nzperiod
!!$        IF(NTV) thetap(ia,il)=thetap(ia,il)+iota*PI/nzperiod
!!$        CALL FILL_BNODE(zetap(il),thetap(ia,il),dummy,B_al(ia,il),vds_al(:,ia,il),.FALSE.)
!!$        IF(USE_B0) CALL FILL_BNODE(zetap(il),thetap(ia,il),dummy,dummy,vds_al(:,ia,il),.TRUE.)
!!$        IF(DEBUG) WRITE(3000+myrank,'(5(1pe13.5),2I5)') zetap(il),thetap(ia,il),&
!!$             &  zetax(ia,il),thetap(ia,il),B_al(ia,il),ia,il
!!$     END DO
!!$  END DO
!!$!  END IF

  DO iturn=1,NTURN
     DO ia=1,nalpha/NTURN !note that dalphap can be negative, if iota is, but...
        ja=ia+(iturn-1)*nalpha/NTURN
        DO il=1,nalphab
           IF(na.EQ.1) THEN
              thetax(ja,il)=offset+(ja-1)*dalphap(1)+(il-1)*iota*dzeta
           ELSE
              thetax(ja,il)=offset+ja*dalphap(1)+(il-1)*iota*dzeta
           END IF
           thetap(ja,il)=thetax(ja,il)-SIOTA*TWOPI*(iturn-1)
!           IF(NTV) thetax(ja,il)=thetap(ja,il)+iota*PI/nzperiod
           zetax(ja,il) =zetap(il)+INT((ja-1)/na)*TWOPI/nzperiod
           j_al(ja,il)=(ja-1)*nalphab+il
!           zetal(jl(ja,il))=zetax(ja,il)
!           thetal(jl(ja,il))=thetax(ja,il)
           CALL FILL_BNODE(zetap(il),thetap(ja,il),dummy,B_al(ja,il),vds_al(:,ja,il),.FALSE.)
           IF(USE_B0) CALL FILL_BNODE(zetap(il),thetap(ja,il),dummy,dummy,vds_al(:,ja,il),.TRUE.)
        END DO
     END DO
  END DO

  DO ia=1,nalpha
     alphap(ia)=thetap(ia,1)-iota*zetap(1)  !...alpha has the correct sign
  END DO

  IF(NTURN.EQ.2) THEN
     dalphap(2)=thetap(nalpha/NTURN+1,1)-thetap(1,1)
     dalphap(3)=thetap(2,1)-thetap(nalpha/NTURN+1,1)
  END IF

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE CREATE_ANGULAR_GRID


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE EXCLUDE_WELLS(na,nalpha,nalphab,nw,bottom,connected,&
     & alphap_w,z1,zb,z2,Bb,Bt,zetax,thetax)

!-----------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------- 

  USE GLOBAL  
  IMPLICIT NONE
  !Input
  LOGICAL bottom(nw),connected(nw,nw)
  INTEGER na,nalpha,nalphab,nw
  REAL*8 alphap_w(nw),z1(nw),zb(nw),z2(nw),Bb(nw),Bt(nw)
  REAL*8 zetax(nalpha,nalphab),thetax(nalpha,nalphab)
  !Others
  LOGICAL left(nalpha),right(nalpha)
  INTEGER ja,ia,iw,jw
  REAL*8 alpp,tempr(nalpha),templ(nalpha)
  !Time
  CHARACTER*30, PARAMETER :: routine="EXCLUDE_WELLS"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  left=.FALSE.
  right=.FALSE.
  templ=thetax(:,1)
  tempr=thetax(:,nalphab)
  DO ia=1,na
     IF(siota.GT.0) THEN
        ja=MINLOC(templ,1)
        templ(ja)=1E10
        left(ja)=.TRUE.
        ja=MAXLOC(tempr,1)
        tempr(ja)=-1E10
        right(ja)=.TRUE.
     ELSE
        ja=MAXLOC(templ,1)
        templ(ja)=-1E10
        left(ja)=.TRUE.
        ja=MINLOC(tempr,1)
        tempr(ja)=+1E10
        right(ja)=.TRUE.
     END IF
  END DO
     
  !Ignore wells that are outside region of interest
  DO iw=1,nw
     IF(.NOT.bottom(iw)) CYCLE
     DO ia=1,nalpha
        IF(.NOT.left(ia)) CYCLE
        alpp=(thetax(ia,1)-iota*zetax(ia,1))
        IF((ABS(alphap_w(iw)-alpp).LT.PREC_EXTR).AND.&
         & (z2(iw).LT.zetax(ia,1).OR.(zb(iw).LT.zetax(ia,1).AND.z2(iw)-z1(iw).GT.PI/nzperiod))) THEN
           Bb(iw)=1E10
           Bt(iw)=1E-9
           DO jw=iw+1,nw 
              IF(connected(jw,iw)) THEN
                 Bb(jw)=1E10
                 Bt(jw)=1E-9
              END IF
           END DO
        END IF
     END DO
     DO ia=1,nalpha
        IF(.NOT.right(ia)) CYCLE
        alpp=(thetax(ia,nalphab)-iota*zetax(ia,nalphab))
!        IF(ABS(alphap_w(iw)-alpp).LT.PREC_EXTR.AND.z1(iw).GT.zetax(ia,nalphab)) THEN
        IF((ABS(alphap_w(iw)-alpp).LT.PREC_EXTR).AND.&
             & (z1(iw).GT.zetax(ia,nalphab).OR.&
             & (zb(iw).GT.zetax(ia,nalphab).AND.z2(iw)-z1(iw).GT.PI/nzperiod))) THEN
           Bb(iw)=1E10
           Bt(iw)=1E-9
           DO jw=iw+1,nw 
              IF(connected(jw,iw)) THEN
                 Bb(jw)=1E10
                 Bt(jw)=1E-9
              END IF
           END DO
        END IF
     END DO
  END DO

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)
  
END SUBROUTINE EXCLUDE_WELLS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CREATE_LAMBDA_GRID(nlambda,nw,Bb,Bt,&
     & lambdab_w,lambdac_w,lambdac,dlambdap,lambda,one_o_lambda)

!-----------------------------------------------------------------------------------------------
!Create a lambda and one_lambda grids with nlambda elements and dlambdap spacing
!from the values of the bottoms Bb and tops Bt of nw wells
!(the extremes in lambda, lambdac_w and lambdab_w for each well are calculated)
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL  
  IMPLICIT NONE
  !Input
  INTEGER nlambda,nw
  REAL*8 Bb(nw),Bt(nw)
  !Output
  REAL*8 lambdab_w(nw),lambdac_w(nw),lambdac,dlambdap,lambda(nlambda),one_o_lambda(nlambda)
  !Others
  CHARACTER*100 serr
  LOGICAL passing
  INTEGER ila,iw
  REAL*8 lambdab,B1,B2,dummy,vdummy(Nnmp),dlambda(nlambda),fact
  !Time
  CHARACTER*30, PARAMETER :: routine="CREATE_LAMBDA_GRID"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  !Create uniform grid
  lambdab_w=1./Bb
  lambdab=MAXVAL(lambdab_w)
  lambdac_w=1./Bt
  lambdac=MINVAL(lambdac_w)+PREC_EXTR
  dlambdap=(lambdab-lambdac)/nlambda
  DO ila=1,nlambda
     lambda(ila)=lambdac+(ila-1)*dlambdap
  END DO
  !Change to grid that is thinner closer to lambda
  IF(ILAGRID) THEN
     lambda(1) =lambdac  
     dlambda(1)=PREC_B*lambdac*lambdac
     DO ila=2,nlambda
        lambda(ila)=lambda(ila-1)+dlambda(ila-1)
        dlambda(ila)=MIN(lambda(ila)-lambdac,dlambdap)
     END DO
     fact=(lambda(nlambda)-lambdac)/(lambdab-lambdac-dlambdap)
     lambda=lambdac+(lambda-lambdac)/fact
  END IF
  one_o_lambda=1./lambda
  DO ila=1,nlambda
     passing=.TRUE.
     DO iw=1,nw
        IF(lambda(ila).GE.1/Bt(iw)) passing=.FALSE.
     END DO
     IF(passing) THEN
        serr="Passing trajectories"
        CALL END_ALL(serr,.FALSE.)
!        WRITE(iout,*) 'Passing trajectories at ila=',ila
     END IF
  END DO

  CALL CALCB(ZERO       ,ZERO,0,.FALSE.,B1,&
       & dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,vdummy)
  CALL CALCB(ZERO,PI,0,.FALSE.,B2,&
       & dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,vdummy)
  lambdac=1./MAX(B1,B2)

  WRITE(iout,'(" Velocity grid has size ",I4," in lambda")') nlambda
  WRITE(iout,'(" (lambda_c=",1pe23.16,", lambda_b=",1pe13.6,") ")') lambdac,lambdab
  
  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)
  
END SUBROUTINE CREATE_LAMBDA_GRID


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE LABEL_GRIDPOINTS(nalpha,nalphab,nlambda,nw,bottom,connected,&
     & alphap_w,z1,z2,Bb,Bt,lambda,zetap,zetax,thetax,B_al,npoint,i_l,i_w,i_p)

!-----------------------------------------------------------------------------------------------
!For a nalpha x nlambda grid, and nw wells characterized by the matrix of connections 
!connected, alphap_w, z1, z2, Bt Bb, determine number of points npoint, alpha, lambda and well
!of each of them, i_l and i_w, and determine point number. i_p for each point of the extended
!grid for each value of lambda
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL  
  IMPLICIT NONE
  !Input
  LOGICAL bottom(nw),connected(nw,nw)
  INTEGER nalpha,nalphab,nlambda,nw
  REAL*8 alphap_w(nw),z1(nw),z2(nw),Bb(nw),Bt(nw),lambda(nlambda)
  REAL*8 zetap(nalphab),zetax(nalpha,nalphab)
  REAL*8 thetax(nalpha,nalphab),B_al(nalpha,nalphab)
  !Output
  INTEGER npoint,i_l(npointx),i_w(npointx),i_p(nlambda,nalpha,nalphab)
  !Others
  CHARACTER*100 serr
  LOGICAL passing(nw)
  INTEGER ja,ia,il,jl,ila,iw,jw,assigned(nw),iwell(nalpha,nalphab,nlambda)
  REAL*8 one_o_lambda(nlambda),alpp,d1,d2,dist,distt
  !Time
  CHARACTER*30, PARAMETER :: routine="LABEL_GRIDPOINTS"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  !For each point in the (zeta,theta) extended grid, determine point
  
  !Start by determine well label for each point
  one_o_lambda=1/lambda
  iwell=0
  DO ia=1,nalpha
     DO il=1,nalphab
        alpp=(thetax(ia,il)-iota*zetax(ia,il))
        DO iw=1,nw
           IF(bottom(iw).AND.ABS(alphap_w(iw)-alpp)&
                & .LT.PREC_EXTR.AND.zetax(ia,il).GT.z1(iw).AND.zetax(ia,il).LE.z2(iw)) EXIT
        END DO
        DO jw=iw,nw
           IF(iw.EQ.jw.OR.connected(jw,iw)) THEN
              DO ila=1,nlambda
                 IF(MAX(bb(jw),B_al(ia,il)).LT.one_o_lambda(ila)&
                      & .AND.one_o_lambda(ila).LE.bt(jw)) iwell(ia,il,ila)=jw
              END DO
           END IF
        END DO
     END DO
  END DO

  !Trapped trajectories closest to 1/Bmax (also longest) are considered passing
  npoint=0
  DO WHILE(npoint.LE.1)
     passing=.FALSE.
     i_p=0
     DO ia=1,nalpha
        DO il=1,nalphab
           iw=iwell(ia,il,1)
           IF(iw.EQ.0) CYCLE
           IF(npoint.EQ.0) passing(iw)=.TRUE.
           i_p(1,ia,il)=1
        END DO
     END DO
     i_l(1)=1
     dist=0
     DO iw=1,nw
        IF(passing(iw).AND.(z2(iw)-z1(iw).GT.dist)) THEN
           dist=z2(iw)-z1(iw)
           i_w(1)=iw
        END IF
     END DO
     i_l(2:npointx)=0
     i_w(2:npointx)=0
     i_p(2:nlambda,:,:)=0
     npoint=1
     DO ila=2,nlambda
        assigned=0
        DO ia=1,nalpha
           DO il=1,nalphab
              iw=iwell(ia,il,ila)
              IF(iw.EQ.0) CYCLE
              IF(passing(iw)) THEN
                 i_p(ila,ia,il)=1
                 CYCLE
              END IF
              IF(assigned(iw).GT.0) THEN
                 i_p(ila,ia,il)=assigned(iw)
              ELSE
                 npoint=npoint+1
                 i_p(ila,ia,il)=npoint
                 i_l(npoint)   =ila
                 i_w(npoint)   =iw
                 assigned(iw)=npoint
              END IF
           END DO
        END DO
     END DO
  END DO
  
  IF(npoint.GT.npointx) THEN
     serr="npoint>npointx"
     CALL END_ALL(serr,.FALSE.)
  END IF

  !Contour conditions: for the points ignored in the first loop of the routine,
  !use periodicity to find equivalent points that have been labelled 
  DO ia=1,nalpha
     DO il=1,nalphab
        alpp=(thetax(ia,il)-iota*zetax(ia,il))
        !alpha=~0,=~0 
        DO iw=1,nw
           IF(bottom(iw).AND.ABS(alphap_w(iw)-alpp).LT.PREC_EXTR.AND.&
                & ((zetax(ia,il).GT.z1(iw)))) EXIT
        END DO
        IF(iw.EQ.nw+1) THEN           
           dist=1E5
           DO ja=1,nalpha
              DO jl=1,nalphab
                 IF(MOD(zetax(ja,jl),TWOPI/nzperiod)-zetap(nalphab-1).LT.ALMOST_ZERO) CYCLE
                 d1= zetax(ia,il)- zetax(ja,jl)+TWOPI/nzperiod
                 d2=thetax(ia,il)-thetax(ja,jl)+TWOPI*siota
                 distt=d1*d1+d2*d2
                 IF(distt.LT.dist) THEN
                    DO ila=1,nlambda
                       IF(one_o_lambda(ila).GT.B_al(ia,il)) THEN
                          dist=distt
                          i_p(  ila,ia,il)=i_p(  ila,ja,jl)
                          iwell(ia,il,ila)=iwell(ja,jl,ila)                          
                       END IF
                    END DO
                 END IF
              END DO
           END DO
        END IF
        !alpha=~2pi,zeta=~2pi/nzperiod  
        DO iw=1,nw
           IF(bottom(iw).AND.ABS(alphap_w(iw)-alpp).LT.PREC_EXTR.AND.&
                & ((zetax(ia,il).LT.z2(iw)))) EXIT
        END DO
        IF(iw.EQ.nw+1) THEN           
           dist=1E5
           DO ja=1,nalpha
              DO jl=1,nalphab
                 IF(zetap(2)-zetax(ja,jl).LT.ALMOST_ZERO) CYCLE
                 d1= zetax(ia,il)- zetax(ja,jl)-TWOPI/nzperiod
                 d2=thetax(ia,il)-thetax(ja,jl)-TWOPI*siota
                 distt=d1*d1+d2*d2
                 IF(distt.LT.dist) THEN
                    DO ila=1,nlambda
                       IF(one_o_lambda(ila).GT.B_al(ia,il)) THEN
                          dist=distt
                          i_p(  ila,ia,il)=i_p(  ila,ja,jl)
                          iwell(ia,il,ila)=iwell(ja,jl,ila)
                       END IF
                    END DO
                 END IF
              END DO
           END DO
        END IF
        
     END DO
  END DO
  
  WRITE(iout,'(" Global grid ",I6," points")') npoint

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE LABEL_GRIDPOINTS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE SORT_ALPHA(nalpha,nalphab,alphap,zetax,thetax,thetap,B_al,vds_al,j_al,nlambda,i_p)

  USE GLOBAL  
  IMPLICIT NONE
  !Input
  INTEGER nalpha,nalphab,nlambda
  !Input/output
  INTEGER i_p(nlambda,nalpha,nalphab),j_al(nalpha,nalphab)
  REAL*8 zetax(nalpha,nalphab),thetax(nalpha,nalphab),thetap(nalpha,nalphab)
  REAL*8 alphap(nalpha),B_al(nalpha,nalphab),vds_al(Nnmp,nalpha,nalphab)
  !Others
  INTEGER inm,isave,ialpha,ila,itemp(nalphab)
  REAL*8 temp(nalphab)
  
  DO ialpha=1,nalpha
     IF(siota.GT.0) THEN
        isave=MINLOC(alphap(ialpha:nalpha),1)+ialpha-1
     ELSE
        isave=MAXLOC(alphap(ialpha:nalpha),1)+ialpha-1
     END IF
     temp(1)        =alphap(ialpha)
     alphap(ialpha) =alphap(isave)
     alphap(isave)  =temp(1)
     temp           =zetax(ialpha,:)
     zetax(ialpha,:)=zetax(isave,:)
     zetax(isave,:) =temp
     temp            =thetax(ialpha,:)
     thetax(ialpha,:)=thetax(isave,:)
     thetax(isave,:) =temp
     temp            =thetap(ialpha,:)
     thetap(ialpha,:)=thetap(isave,:)
     thetap(isave,:) =temp
     temp          =B_al(ialpha,:)
     B_al(ialpha,:)=B_al(isave,:)
     B_al(isave,:) =temp
     itemp         =j_al(ialpha,:)
     j_al(ialpha,:)=j_al(isave,:)
     j_al(isave,:) =itemp
     DO inm=1,Nnmp
        temp                =vds_al(inm,ialpha,:)
        vds_al(inm,ialpha,:)=vds_al(inm,isave,:)
        vds_al(inm,isave,:) =temp
     END DO
     DO ila=1,nlambda
        itemp                =i_p(ila,ialpha,:)
        i_p(ila,ialpha,:)=i_p(ila,isave,:)
        i_p(ila,isave,:) =itemp
     END DO
  END DO

END SUBROUTINE SORT_ALPHA


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 SUBROUTINE FIND_LAMBDA_NEIGHBOURS(npoint,nalpha,nalphab,nlambda,nw,nbifx,i_p,bottom,i_w,&
        & i_p_lm1,nbif,i_p_lp1)

!----------------------------------------------------------------------------------------------- 
!For each point i_p of npoint, find in nalpha x nlambda grid (with nbifx as maximum number of 
!inmediate neighbours at larger lambda), the neighbours in lambda i_p_lm1, nbif, i_p_lp1 are found
!using that there are nw wells characterized by i_w and bottom
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL  
  IMPLICIT NONE
  !Input
  LOGICAL bottom(nw)
  INTEGER npoint,nalpha,nalphab,nlambda,nw,nbifx,i_p(nlambda,nalpha,nalphab),i_w(npoint)
  !Output
  INTEGER i_p_lm1(npoint),i_p_lp1(nbifx,npoint),nbif(npoint)
  !Others
  CHARACTER*100 serr
  INTEGER ia,il,ila,ipoint,jpoint
  !Time
  CHARACTER*30, PARAMETER :: routine="FIND_LAMBDA_NEIGHBOURS"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  !Default values are always npoint
  i_p_lm1(1)=0    
  i_p_lm1(2:npoint)=1
  !Decreasing lambda: one neighbour, since no bifurcations
  DO ia=1,nalpha
     DO il=1,nalphab
        DO ila=2,nlambda
           ipoint=i_p(ila,ia,il)
           IF(ipoint.GT.1) i_p_lm1(ipoint)=MAX(i_p(ila-1,ia,il),i_p_lm1(ipoint))
        END DO 
     END DO
  END DO
  
  !Increasing lambda: several neighours are possible in a bifurcation
  nbif=0
  i_p_lp1=0
  DO ipoint=1,npoint
     DO jpoint=2,npoint
        IF (i_p_lm1(jpoint).EQ.ipoint) THEN
           nbif(ipoint)=nbif(ipoint)+1
           IF(nbif(ipoint).GT.nbifx) THEN
              WRITE(serr,"(A12,I5,A4,I3,A6)") "nbif(ipoint=",ipoint,") = ",nbif(ipoint)," > nbifx"
              CALL END_ALL(serr,.FALSE.)
           END IF
           i_p_lp1(nbif(ipoint),ipoint)=jpoint
        END IF
     END DO
     IF(i_p_lp1(1,ipoint).EQ.0.AND..NOT.bottom(i_w(ipoint))) THEN
        bottom(i_w(ipoint))=.TRUE.

     END IF
  END DO

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE FIND_LAMBDA_NEIGHBOURS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE FIND_ALPHA_NEIGHBOURS(npoint,i_p,i_w,i_l,i_p_lm1,i_p_lp1,nbifx,nbif,zlw,zrw,&
     & nw,connected,BI6,nlambda,nalpha,nalphab,alphap,&
     & i_p_am1,i_p_ap1,dalpha_am1,dalpha_ap1,&
     & i_p_am2,i_p_ap2,dalpha_am2,dalpha_ap2,&
     & i_p_am1I,i_p_am1II,i_p_am2I,i_p_am2II,i_p_am2III,i_p_am2IV,&
     & i_p_ap1I,i_p_ap1II,i_p_ap2I,i_p_ap2II,i_p_ap2III,i_p_ap2IV,&
     & wm1I,wm1II,wm2I,wm2II,wm2III,wm2IV,wp1I,wp1II,wp2I,wp2II,wp2III,wp2IV,lambda,&
     & BI7,BI3,BI3b,BI3f,dlambda_lm1,dlambda_lp1)
  
!----------------------------------------------------------------------------------------------- 
!For each point of npoint, located in nalpha x nlambda (alphap,one_o_lambda) grid,
!i_p_am1, i_p_am2, i_p_ap1, i_p_ap2 are found
!using that there are nw wells characterized by i_w, i_l, zlw, zrw, connected
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL  
  IMPLICIT NONE
  !Input
  INTEGER i_p(nlambda,nalpha,nalphab),i_w(npoint),i_l(npoint),i_p_lm1(npoint),i_p_lp1(nbifx,npoint)
  INTEGER nbifx,nbif(npoint),npoint,nw,nalpha,nalphab,nlambda
  LOGICAL connected(nw,nw)
  REAL*8 zlw(npoint),zrw(npoint)
  REAL*8 BI6(npoint),alphap(nalpha),lambda(nlambda),BI7(npoint),BI3(npoint)
  REAL*8 dlambda_lm1(npoint),dlambda_lp1(npoint)
  !Output
  INTEGER i_p_am1(npoint),i_p_ap1(npoint),i_p_am2(npoint),i_p_ap2(npoint)
  REAL*8 dalpha_am1(npoint),dalpha_ap1(npoint),dalpha_am2(npoint),dalpha_ap2(npoint)
  INTEGER i_p_am1I(npoint),i_p_am1II(npoint),i_p_ap1I(npoint),i_p_ap1II(npoint)
  INTEGER i_p_am2I(npoint),i_p_am2II(npoint),i_p_am2III(npoint),i_p_am2IV(npoint)
  INTEGER i_p_ap2I(npoint),i_p_ap2II(npoint),i_p_ap2III(npoint),i_p_ap2IV(npoint)
  REAL*8 wm1I(npoint),wm1II(npoint),wm2I(npoint),wm2II(npoint),wm2III(npoint),wm2IV(npoint)
  REAL*8 wp1I(npoint),wp1II(npoint),wp2I(npoint),wp2II(npoint),wp2III(npoint),wp2IV(npoint)
  REAL*8 BI3b(npoint),BI3f(npoint)
  !Others
  INTEGER, PARAMETER :: nrep=5
  LOGICAL abottom,gotit
  INTEGER ibif,ipoint,jpoint,ila,jla,ial,jal,il,id,ntot,nvec,ia_min
  REAL*8 da,dat,dz(-1:1,npoint)
  INTEGER kpoint,lpoint,opoint
  INTEGER jwp,iwp,wp(npoint,4),nwp(4),imax(4)
  REAL*8 Jlambda,Jlambda2,Jnew,Jold,Jmax,Jmax2
  REAL*8 mdlBI1opsitf(npoint),mdlBI1opsitb(npoint)
  !Time
  CHARACTER*30, PARAMETER :: routine="FIND_ALPHA_NEIGHBOURS"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  !Determine where the alpha grid is thinner
  ia_min=-100
  IF(MOD(nalpha,2).EQ.0.AND.CENTERED_ALPHA) THEN
     da=1E3
     DO ial=1,nalpha
        jal=ial+1
        IF(jal.EQ.nalpha+1) jal=1
        dat=ABS(alphap(ial)-alphap(jal))
        IF(dat.GT.PI) dat=TWOPI-dat
        IF(dat.LT.da) THEN
           ia_min=ial
           da=dat
        END IF
     END DO
  END IF

  !For each point in the (zeta,theta) grid, determine two neighbours in alpha
  i_p_ap1=0
  i_p_am1=0
  dalpha_ap1=0
  dalpha_am1=0
  dz=0
  DO id=-1,1,2
     DO ila=1,nlambda
        DO ial=1,nalpha
           DO il=1,nalphab
              ipoint=i_p(ila,ial,il)
              IF(ipoint.LE.1) CYCLE !JL
              IF(id*siota.LT.0) THEN
                 jal=ial-1
                 IF(jal.EQ.0) jal=nalpha
              ELSE
                 jal=ial+1
                 IF(ial.EQ.ia_min.OR.ial.EQ.ia_min+1) jal=jal+1
                 IF(jal.GT.nalpha) jal=jal-nalpha
              END IF
              da=ABS(alphap(ial)-alphap(jal))
              IF(da.GT.PI) da=TWOPI-da
              IF(id.LT.0) THEN
                 dalpha_am1(ipoint)=da
              ELSE
                 dalpha_ap1(ipoint)=da
              END IF
              jpoint=i_p(ila,jal,il)
              IF(jpoint.EQ.ipoint.OR.jpoint.EQ.0) CYCLE
              IF(zrw(jpoint)-zlw(jpoint).GT.dz(id,ipoint)) THEN
                 dz(id,ipoint)=zrw(jpoint)-zlw(jpoint)
                 IF(id.LT.0) THEN
                    i_p_am1(ipoint)=jpoint
                 ELSE
                    i_p_ap1(ipoint)=jpoint
                 END IF
              END IF
           END DO
        END DO
     END DO
  END DO

  !Make sure each neighbours are consistent
  DO ipoint=2,npoint
     IF(i_p_ap1(ipoint).EQ.0) THEN
        DO jpoint=2,npoint
           IF(i_p_am1(jpoint).EQ.ipoint) THEN
              i_p_ap1(ipoint)=jpoint
              dalpha_ap1(ipoint)=dalpha_am1(jpoint)
           END IF
        END DO
     END IF
     IF(i_p_am1(ipoint).EQ.0) THEN
        DO jpoint=2,npoint
           IF(i_p_ap1(jpoint).EQ.ipoint) THEN
              i_p_am1(ipoint)=jpoint
              dalpha_am1(ipoint)=dalpha_ap1(jpoint)
           END IF
        END DO
     END IF
  END DO

  !When trajectories do not close in alpha, there may be problems, because no neighbour exist
  !at same value of lambda. Here, two different models are implemented  
  IF(.NOT.CLOSEST_LAMBDA) THEN !Find largest lambda with closed trajectories in alpha
     jla=0
     DO ila=1,nlambda
        ntot=0
        nvec=0
        DO ipoint=2,npoint
           IF(i_l(ipoint).NE.ila) CYCLE
           ntot=ntot+1
           IF(i_p_ap1(ipoint).NE.0) nvec=nvec+1
        END DO
        IF(nvec.EQ.ntot.AND.ila.GT.jla) jla=ila
     END DO
  END IF  
  DO ipoint=2,npoint
     IF(i_p_ap1(ipoint).EQ.0) THEN
        jpoint=ipoint        
        DO jla=i_l(ipoint),1,-1
           jpoint=i_p_lm1(jpoint)
           IF(jpoint.EQ.1) EXIT
           IF((CLOSEST_LAMBDA.AND.i_p_ap1(jpoint).NE.0).OR.&
                (.NOT.CLOSEST_LAMBDA.AND.i_l(jpoint).EQ.ila)) THEN
              i_p_ap1(ipoint)   =i_p_ap1(jpoint)
              dalpha_ap1(ipoint)=dalpha_ap1(jpoint)
              EXIT
           END IF
        END DO
     END IF          
     IF(i_p_am1(ipoint).EQ.0) THEN
        jpoint=ipoint
        DO jla=i_l(ipoint),1,-1
           jpoint=i_p_lm1(jpoint)
           IF(jpoint.EQ.1) EXIT
           IF((CLOSEST_LAMBDA.AND.i_p_am1(jpoint).NE.0).OR.&
                (.NOT.CLOSEST_LAMBDA.AND.i_l(jpoint).EQ.ila)) THEN
              i_p_am1(ipoint)   =i_p_am1(jpoint)
              dalpha_am1(ipoint)=dalpha_am1(jpoint)
              EXIT
           END IF
        END DO
     END IF
  END DO

  !Trajectories very close to the top
  DO ipoint=npoint,2,-1
     IF(nbif(ipoint).EQ.0) CYCLE
     IF(i_p_am1(ipoint).EQ.0) i_p_am1(ipoint)=i_p_am1(i_p_lp1(1,ipoint))
     IF(i_p_ap1(ipoint).EQ.0) i_p_ap1(ipoint)=i_p_ap1(i_p_lp1(1,ipoint))
  END DO

  !Initial values
  i_p_ap1I  =i_p_ap1
  i_p_ap1II =0
  i_p_ap2I  =0
  i_p_ap2II =0
  i_p_ap2III=0
  i_p_ap2IV =0
  i_p_am1I  =i_p_am1
  i_p_am1II =0
  i_p_am2I  =0
  i_p_am2II =0
  i_p_am2III=0
  i_p_am2IV =0
  wp1I  =1.0
  wp1II =0.0
  wp2I  =0.0
  wp2II =0.0
  wp2III=0.0
  wp2IV =0.0
  wm1I  =1.0
  wm1II =0.0
  wm2I  =0.0
  wm2II =0.0
  wm2III=0.0
  wm2IV =0.0

  IF(MODEL_ALPHAD) GOTO 123

  !Find neighbours at constant J
  DO ipoint=2,npoint
     !id=-1, alpha_{i-1}; id=1, alpha_{i+1}     
     DO id=-1,1,2
        nwp=0
        wp=0
        imax=0
        Jlambda=0
        Jlambda2=0
        Jmax=0
        Jmax2=0
        IF(id.LT.0) THEN
           jpoint=i_p_am1(ipoint)
        ELSE
           jpoint=i_p_ap1(ipoint)
        END IF
        !Sum over W in old alpha, if needed
        DO kpoint=2,npoint
           IF(id.LT.0) THEN
              lpoint=i_p_am1(kpoint)
           ELSE
              lpoint=i_p_ap1(kpoint)
           END IF
           IF(lpoint.EQ.jpoint) THEN
              Jlambda=Jlambda+BI6(kpoint)
              nwp(1)=nwp(1)+1
              wp(nwp(1),1)=kpoint
              IF(BI6(kpoint).GT.Jmax) THEN
                 Jmax=BI6(kpoint)
                 imax(1)=kpoint
              END IF
           END IF
        END DO
        !Sum over W in old alpha, lambda+1, if needed
        DO iwp=1,nwp(1) 
           kpoint=wp(iwp,1)
           DO opoint=2,npoint
              IF((i_w(kpoint).EQ.i_w(opoint).OR.connected(i_w(kpoint),i_w(opoint)))&
                   & .AND.i_l(opoint).EQ.i_l(kpoint)+1) THEN
                 Jlambda2=Jlambda2+BI6(opoint)
                 nwp(2)=nwp(2)+1
                 wp(nwp(2),2)=opoint
                 IF(BI6(opoint).GT.Jmax2) THEN
                    Jmax2=BI6(opoint)
                    imax(2)=opoint
                 END IF
              END IF
           END DO
        END DO
        !Sum over W in new alpha, if needed        
        DO kpoint=2,npoint
           IF(id.LT.0) THEN
              lpoint=i_p_ap1(kpoint)
           ELSE
              lpoint=i_p_am1(kpoint)
           END IF
           IF(lpoint.EQ.ipoint) THEN
              nwp(3)=nwp(3)+1
              wp(nwp(3),3)=kpoint
           END IF
        END DO
        IF(nwp(3).EQ.0) THEN
           nwp(3)=1
           IF(id.LT.0) THEN
              wp(1,3)=i_p_am1(ipoint)
           ELSE
              wp(1,3)=i_p_ap1(ipoint)
           END IF           
        END IF

        !Scan in lambda towards lambda_c at new alpha
        DO ila=1,nlambda
           IF(i_p_lm1(wp(1,3)).LE.1) EXIT !JL
           DO iwp=1,nwp(3)
              kpoint=wp(iwp,3)
              gotit=.TRUE.
              IF(nbif(i_p_lm1(kpoint)).GT.1) THEN !JL
                 DO ibif=1,nbif(i_p_lm1(kpoint))
                    gotit=.FALSE.
                    DO jwp=1,nwp(3)
                       gotit=gotit.OR.(wp(jwp,3).EQ.i_p_lp1(ibif,i_p_lm1(kpoint)))
                    END DO
                    IF(.NOT.gotit) EXIT
                 END DO
              END IF
              IF(.NOT.gotit) EXIT
           END DO
           IF(.NOT.gotit) EXIT
           DO kpoint=2,npoint
              DO iwp=1,nwp(3)
                 IF(i_p_lm1(wp(iwp,3)).EQ.kpoint) THEN
                    nwp(4)=nwp(4)+1
                    wp(nwp(4),4)=kpoint
                    EXIT
                 END IF
              END DO
           END DO
           IF(nwp(3).GT.nwp(4)) wp(nwp(4)+1:nwp(3),3)=1
           wp(1:nwp(4),3)=wp(1:nwp(4),4)
           wp(1:nwp(4),4)=1
           nwp(3)=nwp(4)
           nwp(4)=0
           Jnew=0
           DO iwp=1,nwp(3)
              Jnew=Jnew+BI6(wp(iwp,3))
           END DO
           IF(Jnew.GT.Jlambda) EXIT           
        END DO
        !Calculate J in new alpha
        Jnew=0
        Jmax=0
        abottom=.FALSE.
        DO iwp=1,nwp(3)
           kpoint=wp(iwp,3)
           Jnew=Jnew+BI6(kpoint)
           IF(BI6(kpoint).GT.Jmax) THEN
              Jmax=BI6(kpoint)
              imax(3)=kpoint
           END IF
           abottom=abottom.OR.(i_p_lp1(1,kpoint).NE.0)
        END DO
        !Scan in lambda towards 1/Bmin at new alpha
        DO ila=1,nlambda
           Jold=Jnew
           IF(ila.GT.1) THEN
              IF(nwp(3).GT.nwp(4)) wp(nwp(4)+1:nwp(3),3)=0
              wp(1:nwp(4),3)=wp(1:nwp(4),4)
              wp(1:nwp(4),4)=0
              nwp(3)=nwp(4)
              nwp(4)=0
              imax(3)=imax(4)
              imax(4)=0
           END IF
           Jnew=0
           Jmax=0
           IF(.NOT.abottom) THEN
              Jnew=0
              nwp(4)=0
              wp(1,4)=0
              imax(4)=0
           ELSE
              DO iwp=1,nwp(3)
                 jpoint=wp(iwp,3)
                 DO ibif=1,nbif(jpoint)
                    kpoint=i_p_lp1(ibif,jpoint)
                    Jnew=Jnew+BI6(kpoint)
                    nwp(4)=nwp(4)+1
                    wp(nwp(4),4)=kpoint
                    IF(BI6(kpoint).GT.Jmax) THEN
                       Jmax=BI6(kpoint)
                       imax(4)=kpoint
                    END IF
                 END DO
              END DO
           END IF
           !Select neighbours and weights
           IF(id.LT.0) THEN
              IF(((Jnew.LT.Jlambda.AND.Jlambda.LE.Jold).OR.&
                   (Jnew.GT.Jlambda.AND.i_p_lp1(1,imax(4)).EQ.0).OR.&
                   (Jold.LT.Jlambda.AND.ila.EQ.1))) THEN
                 DO iwp=1,nwp(1)
                    jpoint=wp(iwp,1)
                    IF(Jnew.GT.Jlambda.AND.i_p_lp1(1,imax(4)).EQ.0& 
                         & .AND.wm1I(jpoint).LT.1) CYCLE
                    IF(Jold.LT.Jlambda.AND.ila.EQ.1 &
                         & .AND.wm1I(jpoint).LT.1) CYCLE
                    IF(Jnew.GT.PREC_EXTR) THEN
                       i_p_am1I(jpoint) =i_p_lm1(imax(4))
                       i_p_am1II(jpoint)=imax(4)
                       wm1I(jpoint)=(Jlambda-Jnew)/(Jold-Jnew)
                    ELSE
                       i_p_am1I(jpoint) =i_p_lm1(imax(3))
                       i_p_am1II(jpoint)=imax(3)
                       wm1I(jpoint)=(Jlambda-Jold)/(Jold-BI6(i_p_lm1(imax(3))))
                    END IF
                 END DO
              END IF
              IF(nwp(2).NE.0) THEN
                 IF(Jlambda.GT.Jnew.AND.Jnew.GT.Jlambda2) THEN
                    DO iwp=1,nwp(4)
                       jpoint=wp(iwp,4)
                       i_p_ap1I(jpoint) =i_p_lm1(imax(2))
                       i_p_ap1II(jpoint) =imax(2)
                       wp1I(jpoint)=(Jnew-Jlambda2)/(Jlambda-Jlambda2)
                    END DO
                 END IF
                 IF(Jlambda.GT.Jold.AND.Jold.GT.Jlambda2) THEN
                    DO iwp=1,nwp(3)
                       jpoint=wp(iwp,3)
                       i_p_ap1I(jpoint) =i_p_lm1(imax(2))
                       i_p_ap1II(jpoint) =imax(2)
                       wp1I(jpoint)=(Jold-Jlambda2)/(Jlambda-Jlambda2)
                    END DO
                 END IF
              END IF
           ELSE
              IF(((Jnew.LT.Jlambda.AND.Jlambda.LE.Jold).OR.&
                   (Jnew.GT.Jlambda.AND.i_p_lp1(1,imax(4)).EQ.0).OR.&
                   (Jold.LT.Jlambda.AND.ila.EQ.1))) THEN
                 DO iwp=1,nwp(1)
                    jpoint=wp(iwp,1)
                    IF(Jnew.GT.Jlambda.AND.i_p_lp1(1,imax(4)).EQ.0& 
                         & .AND.wp1I(jpoint).LT.1) CYCLE
                    IF(Jold.LT.Jlambda.AND.ila.EQ.1 &
                         & .AND.wp1I(jpoint).LT.1) CYCLE
                    IF(Jnew.GT.PREC_EXTR) THEN
                       i_p_ap1I(jpoint) =i_p_lm1(imax(4))
                       i_p_ap1II(jpoint)=imax(4)
                       wp1I(jpoint)=(Jlambda-Jnew)/(Jold-Jnew)
                    ELSE
                       i_p_ap1I(jpoint) =i_p_lm1(imax(3))
                       i_p_ap1II(jpoint)=imax(3)
                       wp1I(jpoint)=(Jlambda-Jold)/(Jold-BI6(i_p_lm1(imax(3))))
                    END IF
                 END DO
              END IF
              IF(nwp(2).NE.0) THEN
                 IF(Jlambda.GT.Jnew.AND.Jnew.GT.Jlambda2) THEN
                    DO iwp=1,nwp(4)
                       jpoint=wp(iwp,4)
                       i_p_am1I(jpoint) =i_p_lm1(imax(2))
                       i_p_am1II(jpoint) =imax(2)
                       wm1I(jpoint)=(Jnew-Jlambda2)/(Jlambda-Jlambda2)
                    END DO
                 END IF
                 IF(Jlambda.GT.Jold.AND.Jold.GT.Jlambda2) THEN
                    DO iwp=1,nwp(3)
                       jpoint=wp(iwp,3)
                       i_p_am1I(jpoint) =i_p_lm1(imax(2))
                       i_p_am1II(jpoint) =imax(2)
                       wm1I(jpoint)=(Jold-Jlambda2)/(Jlambda-Jlambda2)
                    END DO
                 END IF
              END IF
           END IF
           IF(Jold.LT.Jlambda2) EXIT
           IF(imax(4).EQ.0.OR..NOT.abottom) EXIT
        END DO
     END DO
  END DO

  !Transition from trapped to passing
  DO ipoint=2,npoint
     IF(i_p_lm1(ipoint).NE.1) CYCLE
     Jlambda=BI6(ipoint)
     !id=-1
     Jnew=0
     Jmax=0
     imax=0     
     DO jpoint=2,npoint
        IF(i_p_ap1(jpoint).EQ.ipoint) THEN
           Jnew=Jnew+BI6(jpoint)
           IF(BI6(jpoint).GT.Jmax) THEN
              imax(1)=jpoint
              Jmax=BI6(jpoint)
           END IF
        END IF
     END DO
     IF(imax(1).EQ.0) THEN
        imax(1)=i_p_am1(ipoint)
        Jnew=BI6(imax(1))
     END IF
     IF(Jlambda.GE.Jnew) THEN
        i_p_am1I(ipoint) =1
        i_p_am1II(ipoint)=imax(1)
        wm1I(ipoint)=(Jlambda-Jnew)/(BI6(1)-Jnew)
     END IF
     !id=+1
     Jnew=0
     Jmax=0
     imax=0
     DO jpoint=2,npoint
        IF(i_p_am1(jpoint).EQ.ipoint) THEN
           Jnew=Jnew+BI6(jpoint)
           IF(BI6(jpoint).GT.Jmax) THEN
              imax(1)=jpoint
              Jmax=BI6(jpoint)
           END IF
        END IF
     END DO
     IF(imax(1).EQ.0) THEN
        imax(1)=i_p_ap1(ipoint)
        Jnew=BI6(imax(1))
     END IF
     IF(Jlambda.GE.Jnew) THEN
        i_p_ap1I(ipoint) =1
        i_p_ap1II(ipoint)=imax(1)
        wp1I(ipoint)=(Jlambda-Jnew)/(BI6(1)-Jnew)
     END IF
  END DO

  !Separate passing from trapped
!  DO ipoint=1,npoint
!     IF(i_p_am1I(ipoint).EQ.1) wm1I(ipoint)=0.0
!     IF(i_p_ap1I(ipoint).EQ.1) wp1I(ipoint)=0.0
!  END DO
  
  !Check JL
  DO ipoint=1,npoint
     IF(i_p_am1II(ipoint).EQ.0) i_p_am1II(ipoint)=1
     IF(i_p_ap1II(ipoint).EQ.0) i_p_ap1II(ipoint)=1
  END DO
  
  !Look for value of dalpha at other values of lambda
123  DO ipoint=2,npoint
     IF(dalpha_am1(ipoint).LT.PREC_EXTR) THEN
        jpoint=ipoint
        DO ila=i_l(ipoint),nlambda
           jpoint=i_p_lp1(1,jpoint)
           IF(jpoint.EQ.0) EXIT
           IF(dalpha_am1(jpoint).GT.PREC_EXTR) THEN
              dalpha_am1(ipoint)=dalpha_am1(jpoint)
              EXIT
           END IF
        END DO
        jpoint=ipoint
        DO ila=1,i_l(ipoint)
           IF(dalpha_am1(ipoint).GT.PREC_EXTR) EXIT
           jpoint=i_p_lm1(jpoint)
           IF(dalpha_am1(jpoint).GT.PREC_EXTR) THEN
              dalpha_am1(ipoint)=dalpha_am1(jpoint)
              EXIT
           END IF
        END DO
     END IF
     IF(dalpha_ap1(ipoint).LT.PREC_EXTR) THEN
        jpoint=ipoint
        DO ila=i_l(ipoint),nlambda
           jpoint=i_p_lp1(1,jpoint)
           IF(jpoint.EQ.0) EXIT
           IF(dalpha_ap1(jpoint).GT.PREC_EXTR) THEN
              dalpha_ap1(ipoint)=dalpha_ap1(jpoint)
              EXIT
           END IF
        END DO
        jpoint=ipoint
        DO ila=1,i_l(ipoint)
           IF(dalpha_ap1(ipoint).GT.PREC_EXTR) EXIT
           jpoint=i_p_lm1(jpoint)
           IF(dalpha_ap1(jpoint).GT.PREC_EXTR) THEN
              dalpha_ap1(ipoint)=dalpha_ap1(jpoint)
              EXIT
           END IF
        END DO
     END IF
  END DO

  

  IF(NEW_DALPHA) THEN
     DO ipoint=2,npoint
        DO id=-1,1,2
           IF(id.EQ.-1) THEN
              IF(i_p_am1(ipoint).EQ.0) CYCLE
              IF(i_l(ipoint).NE.i_l(i_p_am1(ipoint))) CYCLE !JLNEW
              IF(BI3(ipoint)/BI7(ipoint).GT.0) THEN
                 IF(i_p_lm1(ipoint).EQ.0) CYCLE
                 jpoint=i_p_lm1(ipoint)
                 i_p_am1II(ipoint)=jpoint
                 mdlBI1opsitb(ipoint)=-BI6(jpoint)
                 mdlBI1opsitb(ipoint)=mdlBI1opsitb(ipoint)+BI6(ipoint)
!                 DO ibif=1,nbif(jpoint)                    
!                    mdlBI1opsitb(ipoint)=mdlBI1opsitb(ipoint)+BI6(i_p_lp1(ibif,jpoint))
!                 END DO
                 mdlBI1opsitb(ipoint)=mdlBI1opsitb(ipoint)/dlambda_lm1(ipoint)/(-sgnb)
              ELSE
                 IF(i_p_lp1(1,ipoint).EQ.0) CYCLE
                 jpoint=i_p_lp1(1,ipoint)
                 i_p_am1II(ipoint)=jpoint
                 mdlBI1opsitb(ipoint)=-BI6(ipoint)
                 DO ibif=1,1!nbif(ipoint)                    
                    mdlBI1opsitb(ipoint)=mdlBI1opsitb(ipoint)+BI6(i_p_lp1(ibif,ipoint))
                 END DO
                 mdlBI1opsitb(ipoint)=mdlBI1opsitb(ipoint)/dlambda_lp1(ipoint)/(-sgnb)
              END IF
              i_p_am1I(ipoint)=i_p_am1(ipoint)
              wm1I(ipoint)= -1E3
           ELSE
              IF(i_p_ap1(ipoint).EQ.0) CYCLE
              IF(i_l(ipoint).NE.i_l(i_p_ap1(ipoint))) CYCLE
              IF(BI3(ipoint)/BI7(ipoint).LT.0) THEN
                 IF(i_p_lm1(ipoint).EQ.0) CYCLE
                 jpoint=i_p_lm1(ipoint)
                 i_p_ap1II(ipoint)=jpoint
                 mdlBI1opsitf(ipoint)=-BI6(jpoint)
                 mdlBI1opsitf(ipoint)=mdlBI1opsitf(ipoint)+BI6(ipoint)
!                 DO ibif=1,nbif(jpoint)                    
!                    mdlBI1opsitf(ipoint)=mdlBI1opsitf(ipoint)+BI6(i_p_lp1(ibif,jpoint))
!                 END DO
                 mdlBI1opsitf(ipoint)=mdlBI1opsitf(ipoint)/dlambda_lm1(ipoint)/(-sgnb)
              ELSE
                 IF(i_p_lp1(1,ipoint).EQ.0) CYCLE
                 jpoint=i_p_lp1(1,ipoint)
                 i_p_ap1II(ipoint)=jpoint
                 mdlBI1opsitf(ipoint)=-BI6(ipoint)
                 DO ibif=1,1!nbif(ipoint)                    
                    mdlBI1opsitf(ipoint)=mdlBI1opsitf(ipoint)+BI6(i_p_lp1(ibif,ipoint))
                 END DO
                 mdlBI1opsitf(ipoint)=mdlBI1opsitf(ipoint)/dlambda_lp1(ipoint)/(-sgnb)
              END IF
              i_p_ap1I(ipoint)=i_p_ap1(ipoint)
              wp1I(ipoint)= -1E3
           END IF
        END DO
     END DO
  END IF

  !Print
  IF(DEBUG) THEN
!     IF(ABS(wm1I(ipoint)).LT.ALMOST_ZERO) WRITE(iout,*) ipoint,i_p_lp1(1,ipoint),i_p_lm1(ipoint)
!     IF(ABS(wp1I(ipoint)).LT.ALMOST_ZERO) WRITE(iout,*) -ipoint,i_p_lp1(1,ipoint),i_p_lm1(ipoint)
     DO ipoint=2,npoint
        WRITE(3400+myrank,'(2I6,1pe23.15,2I6,1pe23.15,2I6,2(1pe23.15))') &
             & -ipoint          ,i_l(ipoint)          ,BI6(ipoint),&
             & i_p_am1I(ipoint) ,i_l(i_p_am1I(ipoint)) ,BI6(i_p_am1I(ipoint)),&
             & i_p_am1II(ipoint),i_l(i_p_am1II(ipoint)),BI6(i_p_am1II(ipoint)),wm1I(ipoint)
        WRITE(3400+myrank,'(2I6,1pe23.15,2I6,1pe23.15,2I6,2(1pe23.15))') &
             & ipoint           ,i_l(ipoint)          ,BI6(ipoint),&
             & i_p_ap1I(ipoint) ,i_l(i_p_ap1I(ipoint)) ,BI6(i_p_ap1I(ipoint)),&
             & i_p_ap1II(ipoint),i_l(i_p_ap1II(ipoint)),BI6(i_p_ap1II(ipoint)),wp1I(ipoint)
     END DO
     CALL FLUSH(3400+myrank)
  END IF

  !Other weights and second neighbours
  wm1II=1.0-wm1I
  wp1II=1.0-wp1I
  DO ipoint=2,npoint
     i_p_am2I(  ipoint)=i_p_am1I( i_p_am1I(ipoint))
     i_p_am2II( ipoint)=i_p_am1II(i_p_am1I(ipoint))
     i_p_am2III(ipoint)=i_p_am1I( i_p_am1II(ipoint))
     i_p_am2IV( ipoint)=i_p_am1II(i_p_am1II(ipoint))
     wm2I(  ipoint)=wm1I( i_p_am1I(ipoint ))*wm1I(ipoint)
     wm2II( ipoint)=wm1II(i_p_am1I(ipoint ))*wm1I(ipoint)
     wm2III(ipoint)=wm1I( i_p_am1II(ipoint))*wm1II(ipoint)
     wm2IV( ipoint)=wm1II(i_p_am1II(ipoint))*wm1II(ipoint)
     dalpha_am2(ipoint)=dalpha_am1(i_p_am1I(ipoint))
     i_p_ap2I(  ipoint)=i_p_ap1I( i_p_ap1I(ipoint))
     i_p_ap2II( ipoint)=i_p_ap1II(i_p_ap1I(ipoint))
     i_p_ap2III(ipoint)=i_p_ap1I( i_p_ap1II(ipoint))
     i_p_ap2IV( ipoint)=i_p_ap1II(i_p_ap1II(ipoint))
     wp2I(  ipoint)=wp1I( i_p_ap1I(ipoint ))*wp1I(ipoint)
     wp2II( ipoint)=wp1II(i_p_ap1I(ipoint ))*wp1I(ipoint)
     wp2III(ipoint)=wp1I( i_p_ap1II(ipoint))*wp1II(ipoint)
     wp2IV( ipoint)=wp1II(i_p_ap1II(ipoint))*wp1II(ipoint)
     dalpha_ap2(ipoint)=dalpha_ap1(i_p_ap1I(ipoint))
  END DO
  i_p_am2=i_p_am2I
  i_p_ap2=i_p_ap2I

  !Recalculate BI3 to be consistent with d_alpha J
  BI3f(1)=0
  BI3b(1)=0
  BI3f(2:npoint)=BI3(2:npoint)
  BI3b(2:npoint)=BI3(2:npoint)
  DO ipoint=2,npoint
!     WRITE(iout,*) ipoint,i_l(ipoint),BI7(ipoint)/dalpha_ap1(ipoint)
!     CALL FLUSH(iout)
     IF(i_p_ap1I(ipoint).NE.0.AND.i_p_ap1I(ipoint).NE.0.AND.wp1I(ipoint).GT.-100) &
          & BI3f(ipoint)=(wp1I(ipoint)*lambda(i_l((i_p_ap1I(ipoint))))&
                       & +wp1II(ipoint)*lambda(i_l((i_p_ap1II(ipoint))))&
                       &               -lambda(i_l((ipoint))))&
                       & *BI7(ipoint)/dalpha_ap1(ipoint)
!     IF(i_p_ap1I(ipoint).NE.0.AND.i_p_ap1I(ipoint).NE.0.AND.wp1I(ipoint).LT.-100) &
!   & BI3f(ipoint)=BI3(ipoint)/(mdlBI1opsitf(ipoint)/BI7(ipoint))
        
     IF(i_p_am1I(ipoint).NE.0.AND.i_p_am1I(ipoint).NE.0.AND.wm1I(ipoint).GT.-100) &
          &  BI3b(ipoint)=(     lambda(i_l((ipoint)))&
          &       -wm1I(ipoint)*lambda(i_l((i_p_am1I(ipoint))))&
          &      -wm1II(ipoint)*lambda(i_l((i_p_am1II(ipoint)))))&
          &  *BI7(ipoint)/dalpha_am1(ipoint)
!     IF(i_p_am1I(ipoint).NE.0.AND.i_p_am1I(ipoint).NE.0.AND.wm1I(ipoint).LT.-100) &
!          &  BI3b(ipoint)=BI3(ipoint)/(mdlBI1opsitb(ipoint)/BI7(ipoint))
  END DO 

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE FIND_ALPHA_NEIGHBOURS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE COEFFICIENTS_DKE(npoint,i_w,i_l,nw,z1,t1,B1,hBpp1,vd1,&
                         & zb,tb,Bb,hBppb,vdb,z2,t2,B2,hBpp2,vd2,&
                         & nlambda,lambda,zlw,tlw,zrw,trw, &
                         & BI1,BI2,BI3,BI4,BI5,BI6,BI7,nmodes,BI8)
!----------------------------------------------------------------------------------------------- 
!For npoints characterized by i_w,i_l, with bounce points zlw, tlw, zrw and trw, located in nw 
!wells characterized by z,t,B,hBpp and vd at extremes 1, 2 and b, and in a lambda grid of size
!nlambda, calculate bounce integrals BI1, BI2, BI3, BI4, BI5, BI6, BI7 and BI8, the latter with
!size nmodes
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL  
  IMPLICIT NONE
  !Input
  INTEGER npoint,i_w(npoint),i_l(npoint),nw,nlambda,nmodes
  REAL*8 z1(nw),t1(nw),B1(nw),hBpp1(nw),vd1(nqv,nw)
  REAL*8 zb(nw),tb(nw),Bb(nw),hBppb(nw),vdb(nqv,nw) 
  REAL*8 z2(nw),t2(nw),B2(nw),hBpp2(nw),vd2(nqv,nw) 
  REAL*8 lambda(nlambda),zlw(npoint),tlw(npoint),zrw(npoint),trw(npoint)
  !Output
  REAL*8 BI1(npoint),BI2(npoint),BI3(npoint),BI4(npoint),BI5(npoint),BI6(npoint),BI7(npoint)
  REAL*8 BI8(npoint,nmodes)
  !Others
  INTEGER ipoint,iw,ila,nq
  REAL*8, ALLOCATABLE :: Q(:)
  !Time
  CHARACTER*30, PARAMETER :: routine="COEFFICIENTS_DKE"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  nq=nq0
  IF(SOLVE_QN) nq=nq+nmodes
  ALLOCATE(Q(nq))

  DO ipoint=1,npoint
     iw =i_w(ipoint)
     ila=i_l(ipoint)
     CALL BOUNCES(iw,z1(iw),t1(iw),B1(iw),hBpp1(iw),vd1(:,iw), &
          &          zb(iw),tb(iw),Bb(iw),hBppb(iw),vdb(:,iw), &
          &          z2(iw),t2(iw),B2(iw),hBpp2(iw),vd2(:,iw), &
          & 1./lambda(ila),ipoint.EQ.1,nq,Q,&
          & zlw(ipoint),tlw(ipoint),zrw(ipoint),trw(ipoint))
     BI1(ipoint)=Q(1)
     BI2(ipoint)=Q(2)
     BI3(ipoint)=Q(3)
     BI4(ipoint)=Q(4)
     BI5(ipoint)=Q(5)  
     BI6(ipoint)=Q(6)
     BI7(ipoint)=Q(7)
     IF(SOLVE_QN) BI8(ipoint,1:nmodes)=Q(8:nq)
     IF(DEBUG) WRITE(3100+myrank,'(3I6,8(1pe13.5),I5,4(1pe13.5))') &
          &  ipoint,ila,iw,BI1(ipoint),BI2(ipoint),BI3(ipoint),&
          &  BI4(ipoint),BI5(ipoint),BI6(ipoint),BI7(ipoint),BI8(ipoint,1),&
          & npoint,zlw(ipoint),tlw(ipoint),zrw(ipoint),trw(ipoint)
  END DO
  IF(DEBUG) CALL FLUSH(3100+myrank)

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE COEFFICIENTS_DKE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef MPIandPETSc

SUBROUTINE INIT_LINEAR_PROBLEM(npoint,nnz,matCOL,matVEAf,matVEAb,matVMAf,matVMAb,matDIFf,matDIFb)
 
!----------------------------------------------------------------------------------------------- 
!Initialize linear problem of size npoint and nnz with PETSc: ksp and matrices matCOL, matVEAf
!matVEAb, matVMAf, matVMAb
!-----------------------------------------------------------------------------------------------    
                    
  USE GLOBAL  
  USE petscsys
  USE petscksp
  IMPLICIT NONE
!#include <petsc/finclude/petscksp.h>
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!#include <petsc/finclude/petscpc.h>
  !Input
  INTEGER npoint,nnz(npoint)
  !Output
!  KSP ksp
  Mat matCOL,matVEAf,matVEAb,matVMAf,matVMAb,matDIFf,matDIFb!,matA
  !Others
  PetscErrorCode ierr
  PetscInt iz,innz(0:npoint-1),npoint_ps,npalpha
!  INTEGER, PARAMETER :: MAXITS= 1000 !these choices could go in the parameters namelist
!  REAL, PARAMETER :: ATOL=1E-2
!  REAL, PARAMETER :: TOL =1E-1
!  PetscReal, PARAMETER :: factor=10.
!  PetscInt       Istart,Iend
!  PC  pc
  !Time
  CHARACTER*30, PARAMETER :: routine="INIT_LINEAR_PROBLEM"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  npoint_ps=npoint
  
  !Number of non-zero elements per row
  DO iz=1,npoint
     innz(iz-1)=nnz(iz)
  END DO
  !Collision operator
  CALL MatCreate(PETSC_COMM_WORLD,matCOL,ierr)
  CALL MatSetSizes(matCOL,PETSC_DECIDE,PETSC_DECIDE,npoint_ps,npoint_ps,ierr)
  CALL MatSetType(matCOL,MATAIJ,ierr)
!  CALL MatSeqAIJSetPreallocation(matCOL,PETSC_NULL_INTEGER,innz,ierr)
  CALL MatSeqAIJSetPreallocation(matCOL,0,innz,ierr)
  CALL MatSetup(matCOL,ierr)

  IF(CENTERED_ALPHA) THEN
     npalpha=4
  ELSE IF(SECOND_ORDER_ALPHA) THEN
     npalpha=7
  ELSE
     npalpha=3
  END IF

  !Tangential ExB drift, forward derivative
  CALL MatCreate(PETSC_COMM_WORLD,matVEAf,ierr)
  CALL MatSetSizes(matVEAf,PETSC_DECIDE,PETSC_DECIDE,npoint_ps,npoint_ps,ierr)
  CALL MatSetType(matVEAf,MATAIJ,ierr)
  CALL MatSeqAIJSetPreallocation(matVEAf,npalpha,PETSC_NULL_INTEGER,ierr)
  CALL MatSetup(matVEAf,ierr)
  !Tangential ExB drift, backward derivative
  CALL MatCreate(PETSC_COMM_WORLD,matVEAb,ierr)
  CALL MatSetSizes(matVEAb,PETSC_DECIDE,PETSC_DECIDE,npoint_ps,npoint_ps,ierr)
  CALL MatSetType(matVEAb,MATAIJ,ierr)
  CALL MatSeqAIJSetPreallocation(matVEAb,npalpha,PETSC_NULL_INTEGER,ierr)
  CALL MatSetup(matVEAb,ierr)

  IF(TANG_VM) THEN
     !Tangential magnetic drift, forward derivative
     CALL MatCreate(PETSC_COMM_WORLD,matVMAf,ierr)
     CALL MatSetSizes(matVMAf,PETSC_DECIDE,PETSC_DECIDE,npoint_ps,npoint_ps,ierr)
     CALL MatSetType(matVMAf,MATAIJ,ierr)
     CALL MatSeqAIJSetPreallocation(matVMAf,npalpha,PETSC_NULL_INTEGER,ierr)
     CALL MatSetup(matVMAf,ierr)
     !Tangential magnetic drift, backward derivative
     CALL MatCreate(PETSC_COMM_WORLD,matVMAb,ierr)
     CALL MatSetSizes(matVMAb,PETSC_DECIDE,PETSC_DECIDE,npoint_ps,npoint_ps,ierr)
     CALL MatSetType(matVMAb,MATAIJ,ierr)
     CALL MatSeqAIJSetPreallocation(matVMAb,npalpha,PETSC_NULL_INTEGER,ierr)
     CALL MatSetup(matVMAb,ierr)
  END IF

  IF(CALC_DIFF) THEN
     npalpha=3
     CALL MatCreate(PETSC_COMM_WORLD,matDIFf,ierr)
     CALL MatSetSizes(matDIFf,PETSC_DECIDE,PETSC_DECIDE,npoint_ps,npoint_ps,ierr)
     CALL MatSetType(matDIFf,MATAIJ,ierr)
     CALL MatSeqAIJSetPreallocation(matDIFf,npalpha,PETSC_NULL_INTEGER,ierr)
     CALL MatSetup(matDIFf,ierr)
     CALL MatCreate(PETSC_COMM_WORLD,matDIFb,ierr)
     CALL MatSetSizes(matDIFb,PETSC_DECIDE,PETSC_DECIDE,npoint_ps,npoint_ps,ierr)
     CALL MatSetType(matDIFb,MATAIJ,ierr)
     CALL MatSeqAIJSetPreallocation(matDIFb,npalpha,PETSC_NULL_INTEGER,ierr)
     CALL MatSetup(matDIFb,ierr)
  END IF
  
!!$  CALL KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
!!$  CALL KSPSetTolerances(ksp,TOL,ATOL,PETSC_DEFAULT_REAL,MAXITS,ierr)
!!$!  CALL KSPSetType(ksp,KSPPREONLY,ierr)
!!$!  CALL KSPGetPC(ksp,pc,ierr)
!!$!  CALL PCSetType(pc,PCILU,ierr)
!!$!  CALL PCFactorSetLevels(pc,10,ierr)
!!$!  CALL PCFactorSetMatOrderingType(pc,"natural",ierr)
!!$
!!$  IF(DIRECT_SOLVER) THEN
!!$     CALL KSPSetType(ksp,KSPPREONLY,ierr)
!!$     CALL KSPGetPC(ksp,pc,ierr)
!!$     CALL PCSetType(pc,PCLU,ierr)
!!$     CALL PCFactorSetMatOrderingType(pc,"natural",ierr)
!!$!     CALL PCFactorSetMatOrderingType(pc,"nd",ierr)
!!$!     CALL PCFactorSetReuseOrdering(pc,PETSC_TRUE,ierr)
!!$!     CALL PCFactorSetReuseFill(pc,PETSC_TRUE,ierr)
!!$!#ifdef PETSC_HAVE_MUMPS
!!$!     call PCFactorSetMatSolverType(pc,MATSOLVERMUMPS,ierr)
!!$!     call PCFactorSetUpMatSolverType(pc,ierr)
!!$!#endif
!!$  ELSE
!!$     CALL KSPSetType(ksp,KSPGMRES,ierr)
!!$     CALL KSPGetPC(ksp,pc,ierr)     
!!$!     CALL PCSetType(pc,PCJACOBI,ierr)
!!$     CALL PCSetType(pc,PCILU,ierr)
!!$     CALL PCFactorSetLevels(pc,2,ierr)
!!$     CALL PCFactorSetMatOrderingType(pc,"natural",ierr)
!!$  END IF

!From more to less stable
!  CALL PCFactorSetMatOrderingType(pc,"rcm",ierr)
!  CALL PCFactorSetMatOrderingType(pc,"1wd",ierr)
!  CALL PCFactorSetMatOrderingType(pc,"nd",ierr) !default
!  CALL PCFactorSetMatOrderingType(pc,"qmd",ierr)
!  CALL PCFactorSetMatOrderingType(pc,"rowlength",ierr)
!  CALL PCFactorSetMatOrderingType(pc,"wbm",ierr)
!  CALL PCFactorSetMatOrderingType(pc,"spectral",ierr)
!  CALL PCFactorSetMatOrderingType(pc,"amd",ierr)
!  CALL PCFactorSetReuseOrdering(pc,.TRUE.,ierr)
!Other options
!  CALL PCFactorSetPivotInBlocks(pc,.TRUE.,ierr)
!  CALL PCFactorSetMatSolverType(pc,MATSOLVERMUMPS,ierr)
!  CALL PCFactorSetUpMatSolverType(pc,ierr)
!  CALL PCFactorSetFill(pc,factor,ierr)
  
  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE INIT_LINEAR_PROBLEM

#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE FILL_DKE_ROW(ipoint,npoint,&
     & dalpha_am1,dalpha_ap1,dalpha_am2,dalpha_ap2,i_p_am1,i_p_ap1,i_p_am2,i_p_ap2,&
     & i_p_am1I,i_p_am1II,i_p_am2I,i_p_am2II,i_p_am2III,i_p_am2IV,&
     & i_p_ap1I,i_p_ap1II,i_p_ap2I,i_p_ap2II,i_p_ap2III,i_p_ap2IV,&
     & wm1I,wm1II,wm2I,wm2II,wm2III,wm2IV,&
     & wp1I,wp1II,wp2I,wp2II,wp2III,wp2IV,&
     & lambda,dlambda_lm1,dlambda_lp1,i_p_lm1,nbif,nbifx,i_p_lp1,&
     & BI1,BI2,BI4,BI5,BI3oBI7,matCOL,matVEAf,matVEAb,matVMAf,matVMAb,nnz,flag_nnz,matDIFf,matDIFb)

!----------------------------------------------------------------------------------------------- 
!Fill row ipoint of npoint of several matrixes (matCOL, matVEAf, matVEAb, matVMAf, matVMAb), 
!with nnz non-zero elements, using spacing dlambda_lm1, dlambda_lp1, dalpha, factors
!factDKErhs, lambda, Epsi/vdconst, neighbours i_p_am1, i_p_ap1, i_p_lm1, nbif neighbours
!i_p_lp1(:) and coefficients BI[1-5]
!If flag_nnz, use this routine to calculate nnz only
!-----------------------------------------------------------------------------------------------  
                      
  USE GLOBAL  
  IMPLICIT NONE
  !Input
  LOGICAL flag_nnz
  INTEGER ipoint,npoint,nbifx,nbif(npoint),i_p_lp1(nbifx,npoint)
  INTEGER i_p_am1,i_p_ap1,i_p_am2,i_p_ap2,i_p_lm1(npoint)
  INTEGER i_p_am1I,i_p_am1II,i_p_am2I,i_p_am2II,i_p_am2III,i_p_am2IV
  INTEGER i_p_ap1I,i_p_ap1II,i_p_ap2I,i_p_ap2II,i_p_ap2III,i_p_ap2IV
  REAL*8 wm1I,wm1II,wm2I,wm2II,wm2III,wm2IV,wp1I,wp1II,wp2I,wp2II,wp2III,wp2IV
  REAL*8 lambda,dlambda_lm1(npoint),dlambda_lp1(npoint),dalpha_am1,dalpha_ap1,dalpha_am2,dalpha_ap2
  REAL*8 BI1,BI2(npoint),BI4,BI5,BI3oBI7
  !Output
  INTEGER nnz
  REAL*8 matCOL(npoint),matVEAf(npoint),matVEAb(npoint),matVMAf(npoint),matVMAb(npoint)
  REAL*8 matDIFf(npoint),matDIFb(npoint)
  !Others
  CHARACTER*300 serr
  INTEGER ibif,jbif,j_p_lm1,j_p_lm2,j_p_lp1,j_p_lp2,jpoint
  REAL*8 denom,a,b,bm,bp,suma,sumb,d
  REAL*8 dlambda,twodlambda,dlambda2,fourdlambda2
 !Time
!  CHARACTER*30, PARAMETER :: routine="FILL_DKE_ROW"
!  INTEGER, SAVE :: ntotal=0
!  REAL*8,  SAVE :: ttotal=0
!  REAL*8,  SAVE :: t0=0
!  REAL*8 tstart

!  CALL CPU_TIME(tstart)

  j_p_lm1=i_p_lm1(ipoint)
  IF(j_p_lm1.EQ.1.AND.i_p_lp1(1,ipoint).EQ.0) THEN
     WRITE(serr,"(A7,I5,A9)") "Orbit ",ipoint, " isolated"
     CALL END_ALL(serr,.FALSE.)
  END IF

  IF(nbif(ipoint).EQ.1.AND.nbif(j_p_lm1).EQ.1) THEN
     nnz=3
  ELSE IF(j_p_lm1.EQ.1.AND.nbif(1).GT.1.AND..NOT.PREC_TOP) THEN
     nnz=2
  ELSE
     nnz=2
     DO ibif=1,nbif(j_p_lm1)
        jpoint=i_p_lp1(ibif,j_p_lm1)
        DO jbif=1,nbif(jpoint)
           j_p_lp1=i_p_lp1(jbif,jpoint)
           j_p_lp2=i_p_lp1(1  ,j_p_lp1)
           IF(j_p_lp2.EQ.0) CYCLE 
           nnz=nnz+1
           IF(jpoint.NE.ipoint) nnz=nnz+1
        END DO
     END DO
  END IF
  
  IF(nnz.NE.1) THEN
     IF(CENTERED_ALPHA) THEN
        nnz=nnz+4
     ELSE IF(SECOND_ORDER_ALPHA) THEN
        nnz=nnz+12
     ELSE
        nnz=nnz+4
     END IF
  END IF
  IF(flag_nnz) RETURN

  !Collision operator
  twodlambda=dlambda_lm1(ipoint)+dlambda_lp1(ipoint)
  dlambda2  =dlambda_lm1(ipoint)*dlambda_lp1(ipoint)

  IF(.NOT.PHI1_READ) THEN
     a=(BI2(ipoint)/lambda-BI1*lambda/2.)
     b=BI2(ipoint)
  ELSE
     a=0.5*(3*BI2(ipoint)/lambda-BI1*lambda*SQRT(lambda))
     b=BI2(ipoint)        
  END IF

!  IF(KROOK_OP) THEN
!     matCOL(ipoint)=-BI1
!  ELSE IF(j_p_lm1.EQ.1.AND.i_p_lp1(1,ipoint).EQ.0) THEN
  !     matCOL(ipoint)=1


!!$  IF((nbif(ipoint).EQ.1.AND.nbif(j_p_lm1).EQ.1).OR.&
!!$  IF((nbif(ipoint).EQ.1).OR.&
!!$       (j_p_lm1.EQ.1.AND..NOT.PREC_TOP)) THEN
!!$     j_p_lp1=i_p_lp1(1,ipoint)
!!$     denom=twodlambda*(dlambda_lm1(ipoint)*dlambda_lm1(ipoint)+dlambda_lp1(ipoint)*dlambda_lp1(ipoint))/4
!!$     matCOL(j_p_lp1)= b*dlambda_lp1(ipoint)/denom+ a/twodlambda
!!$     matCOL(ipoint) =-b*twodlambda         /denom
!!$     matCOL(j_p_lm1)= b*dlambda_lm1(ipoint)/denom- a/twodlambda
!!$  ELSE
!!$     DO ibif=1,nbif(j_p_lm1)
!!$        jpoint=i_p_lp1(ibif,j_p_lm1)
!!$!        jpoint=ipoint
!!$        IF(j_p_lm1.EQ.1) THEN
!!$           bm=BI2(jpoint)
!!$           matCOL(i_p_lp1(1,jpoint))=matCOL(i_p_lp1(1,jpoint))-bm/(dlambda_lp1(jpoint)*twodlambda)
!!$        ELSE
!!$           j_p_lm2=i_p_lm1(j_p_lm1)
!!$           bm=BI2(j_p_lm1)
!!$           fourdlambda2=(dlambda_lp1(j_p_lm1)+dlambda_lm1(j_p_lm1))*twodlambda
!!$           matCOL(j_p_lm2)=matCOL(j_p_lm2)+bm/fourdlambda2
!!$           matCOL(jpoint )=matCOL(jpoint) -bm/fourdlambda2
!!$        END IF
!!$
!!$        DO jbif=1,nbif(jpoint)
!!$           j_p_lp1=i_p_lp1(jbif,jpoint)
!!$           j_p_lp2=i_p_lp1(1  ,j_p_lp1)
!!$           IF(j_p_lp2.EQ.0) CYCLE
!!$           bp=BI2(j_p_lp1)
!!$           IF(j_p_lm1.EQ.1) THEN
!!$              fourdlambda2=dlambda_lp1(jpoint)*(dlambda_lp1(j_p_lp1)+dlambda_lm1(j_p_lp1))
!!$           ELSE
!!$              fourdlambda2=twodlambda*(dlambda_lp1(j_p_lp1)+dlambda_lm1(j_p_lp1)) 
!!$           END IF
!!$           matCOL(j_p_lp2)=matCOL(j_p_lp2)+bp/fourdlambda2
!!$           matCOL(jpoint) =matCOL(jpoint) -bp/fourdlambda2
!!$        END DO
!!$     END DO
!!$
!!$     END IF

  IF(nbif(ipoint).EQ.1.AND.nbif(j_p_lm1).EQ.1) THEN
     j_p_lp1=i_p_lp1(1,ipoint)
     suma=a/twodlambda
     sumb=b/dlambda2
     matCOL(j_p_lp1)=sumb+suma
     matCOL(ipoint) =-2*sumb
     matCOL(j_p_lm1)=sumb-suma
  ELSE IF(j_p_lm1.EQ.1.AND.nbif(1).GT.1.AND..NOT.PREC_TOP) THEN
     j_p_lp1=i_p_lp1(1,ipoint)
     denom=twodlambda*(dlambda_lm1(ipoint)*dlambda_lm1(ipoint)+dlambda_lp1(ipoint)*dlambda_lp1(ipoint))/4
     matCOL(j_p_lp1)= b*dlambda_lp1(ipoint)/denom+ a/twodlambda
     matCOL(ipoint) =-b*twodlambda         /denom
  ELSE
     j_p_lp1=i_p_lp1(1,ipoint)
     IF(PREC_TOP.AND.j_p_lm1.EQ.1) THEN
        bm=BI2(ipoint)
        matCOL(j_p_lp1)=-bm/(dlambda_lp1(ipoint)*twodlambda)
     ELSE
        j_p_lm2=i_p_lm1(j_p_lm1)
        bm=BI2(j_p_lm1)
        fourdlambda2=(dlambda_lp1(j_p_lm1)+dlambda_lm1(j_p_lm1))*twodlambda 
        matCOL(j_p_lm2)=bm/fourdlambda2  
        matCOL(ipoint)=-bm/fourdlambda2
     END IF
     DO ibif=1,nbif(j_p_lm1)
        jpoint=i_p_lp1(ibif,j_p_lm1)
        IF(PREC_TOP.AND.j_p_lm1.EQ.1) jpoint=ipoint
        DO jbif=1,nbif(jpoint)
           j_p_lp1=i_p_lp1(jbif,jpoint)
           j_p_lp2=i_p_lp1(1  ,j_p_lp1)
           IF(j_p_lp2.EQ.0) CYCLE
           bp=BI2(j_p_lp1)
           IF(PREC_TOP.AND.j_p_lm1.EQ.1) THEN
              fourdlambda2=dlambda_lp1(jpoint)*(dlambda_lp1(j_p_lp1)+dlambda_lm1(j_p_lp1))
           ELSE
              fourdlambda2=twodlambda*(dlambda_lp1(j_p_lp1)+dlambda_lm1(j_p_lp1)) 
           END IF
           matCOL(j_p_lp2)=bp/fourdlambda2
           matCOL(jpoint)=matCOL(jpoint)-bp/fourdlambda2
        END DO
        IF(PREC_TOP.AND.j_p_lm1.EQ.1) EXIT
     END DO
  END IF
  !d_psi J d_alpha g
  IF(nnz.NE.1) THEN
     d=BI5
     IF(CENTERED_ALPHA) THEN
        denom=dalpha_ap1+dalpha_am1
        matVEAf(i_p_ap1I)  =wp1I  *( d/denom)
        matVEAf(i_p_ap1II) =wp1II *( d/denom)
        matVEAf(i_p_am1I)  =wm1I  *(-d/denom)
        matVEAf(i_p_am1II) =wm1II *(-d/denom)
        matVEAb(i_p_ap1I)  =wp1I  *( d/denom)
        matVEAb(i_p_ap1II) =wp1II *( d/denom)
        matVEAb(i_p_am1I)  =wm1I  *(-d/denom)
        matVEAb(i_p_am1II) =wm1II *(-d/denom)
        
!!$        denom=7*(dalpha_ap1+dalpha_am1)-dalpha_ap2-dalpha_am2
!!$        matVEAf(i_p_ap2I)  =wp2I  *(  -d/denom)
!!$        matVEAf(i_p_ap2II) =wp2II *(  -d/denom)
!!$        matVEAf(i_p_ap2III)=wp2III*(  -d/denom)
!!$        matVEAf(i_p_ap2IV) =wp2IV *(  -d/denom)
!!$        matVEAf(i_p_ap1I)  =wp1I  *(+8*d/denom)
!!$        matVEAf(i_p_ap1II) =wp1II *(+8*d/denom)
!!$        matVEAf(i_p_am1I)  =wm1I  *(-8*d/denom)
!!$        matVEAf(i_p_am1II) =wm1II *(-8*d/denom)
!!$        matVEAf(i_p_am2I)  =wm2I  *(  +d/denom)
!!$        matVEAf(i_p_am2II) =wm2II *(  +d/denom)
!!$        matVEAf(i_p_am2III)=wm2III*(  +d/denom)
!!$        matVEAf(i_p_am2IV) =wm2IV *(  +d/denom)
!!$        matVEAb(i_p_ap2I)  =wp2I  *(  -d/denom)
!!$        matVEAb(i_p_ap2II) =wp2II *(  -d/denom)
!!$        matVEAb(i_p_ap2III)=wp2III*(  -d/denom)
!!$        matVEAb(i_p_ap2IV) =wp2IV *(  -d/denom)
!!$        matVEAb(i_p_ap1I)  =wp1I  *(+8*d/denom)
!!$        matVEAb(i_p_ap1II) =wp1II *(+8*d/denom)
!!$        matVEAb(i_p_am1I)  =wm1I  *(-8*d/denom)
!!$        matVEAb(i_p_am1II) =wm1II *(-8*d/denom)
!!$        matVEAb(i_p_am2I)  =wm2I  *(  +d/denom)
!!$        matVEAb(i_p_am2II) =wm2II *(  +d/denom)
!!$        matVEAb(i_p_am2III)=wm2III*(  +d/denom)
!!$        matVEAb(i_p_am2IV) =wm2IV *(  +d/denom)
     ELSE
        IF(SECOND_ORDER_ALPHA.AND.i_p_ap2I.NE.0) THEN
           denom=dalpha_ap1*dalpha_ap2*(dalpha_ap1+dalpha_ap2)
           matVEAf(i_p_ap2I  )=wp2I  *(-d*dalpha_ap1*dalpha_ap1/denom)
           matVEAf(i_p_ap2II )=wp2II *(-d*dalpha_ap1*dalpha_ap1/denom)
           matVEAf(i_p_ap2III)=wp2III*(-d*dalpha_ap1*dalpha_ap1/denom)
           matVEAf(i_p_ap2IV )=wp2IV *(-d*dalpha_ap1*dalpha_ap1/denom)
           matVEAf(i_p_ap1I )=wp1I* (+d*(dalpha_ap1+dalpha_ap2)*&
                & (dalpha_ap1+dalpha_ap2)/denom)
           matVEAf(i_p_ap1II)=wp1II*(+d*(dalpha_ap1+dalpha_ap2)*&
                & (dalpha_ap1+dalpha_ap2)/denom)                      
           matVEAf(ipoint )=+d*&
                & (dalpha_ap1*dalpha_ap1-(dalpha_ap1+dalpha_ap2)*(dalpha_ap1+dalpha_ap2))&
                & /denom
        ELSE
           IF(NEW_DALPHA.AND.wp1I.LT.-100) THEN
              IF(BI3oBI7.GT.0) THEN
                 dlambda=dlambda_lp1(ipoint)
              ELSE
                 dlambda=dlambda_lm1(ipoint)
              END IF
              matVEAf(i_p_ap1I )= d/dalpha_ap1
              matVEAf(i_p_ap1II)= d*ABS(BI3oBI7)/dlambda
              matVEAf(ipoint )  =-d*(1./dalpha_ap1+ABS(BI3oBI7)/dlambda)
              IF(CALC_DIFF.AND.nbif(ipoint).GT.0) THEN
                 sumb=d*ABS(BI3oBI7)/dlambda2*dlambda
                 j_p_lm1=i_p_lm1(ipoint)
                 j_p_lp1=i_p_lp1(1,ipoint)
                 matDIFf(j_p_lp1)=sumb
                 matDIFf(ipoint) =-2*sumb
                 matDIFf(j_p_lm1)=sumb
              END IF
           ELSE
              matVEAf(i_p_ap1I )=wp1I *(+d/dalpha_ap1)
              matVEAf(i_p_ap1II)=matVEAf(i_p_ap1II)+wp1II*(+d/dalpha_ap1)
              matVEAf(ipoint )         =-d/dalpha_ap1
           END IF
        END IF
        IF(SECOND_ORDER_ALPHA.AND.i_p_am2I.NE.0) THEN
           denom=dalpha_am1*dalpha_am2*(dalpha_am1+dalpha_am2)
           matVEAb(i_p_am2I  )=wm2I  *(+d*dalpha_am1*dalpha_am1/denom)
           matVEAb(i_p_am2II )=wm2II *(+d*dalpha_am1*dalpha_am1/denom)
           matVEAb(i_p_am2III)=wm2III*(+d*dalpha_am1*dalpha_am1/denom)
           matVEAb(i_p_am2IV )=wm2IV *(+d*dalpha_am1*dalpha_am1/denom)
           matVEAb(i_p_am1I  )=wm1I  *(-d*(dalpha_am1+dalpha_am2)*(dalpha_am1+dalpha_am2)/denom)
           matVEAb(i_p_am1II )=wm1II *(-d*(dalpha_am1+dalpha_am2)*(dalpha_am1+dalpha_am2)/denom)
           matVEAb(ipoint )=-d*(dalpha_am1*dalpha_am1-(dalpha_am1+dalpha_am2)*(dalpha_am1+dalpha_am2))&
                & /denom
        ELSE
           IF(NEW_DALPHA.AND.wm1I.LT.-100) THEN
              IF(BI3oBI7.GT.0) THEN
                 dlambda=dlambda_lm1(ipoint)
              ELSE
                 dlambda=dlambda_lp1(ipoint)
              END IF
              matVEAb(i_p_am1I )=-d/dalpha_am1
              matVEAb(i_p_am1II)=-d*ABS(BI3oBI7)/dlambda
              matVEAb(ipoint )  = d*(1/dalpha_am1+ABS(BI3oBI7)/dlambda)
              IF(CALC_DIFF.AND.nbif(ipoint).GT.0) THEN
                 sumb=d*ABS(BI3oBI7)/dlambda2*dlambda
                 j_p_lm1=i_p_lm1(ipoint)
                 j_p_lp1=i_p_lp1(1,ipoint)
                 matDIFb(j_p_lp1)=sumb
                 matDIFb(ipoint) =-2*sumb
                 matDIFb(j_p_lm1)=sumb
              END IF
              ELSE
              matVEAb(i_p_am1I )=wm1I *(-d/dalpha_am1)
              matVEAb(i_p_am1II)=matVEAb(i_p_am1II)+wm1II*(-d/dalpha_am1)
              matVEAb(ipoint )         =+d/dalpha_am1
           END IF
        END IF
     END IF

     IF(TANG_VM) THEN
        d=-BI4
        IF(CENTERED_ALPHA) THEN
           denom=dalpha_ap1+dalpha_am1
           matVMAf(i_p_ap1I)  =wp1I  *( d/denom)
           matVMAf(i_p_ap1II) =wp1II *( d/denom)
           matVMAf(i_p_am1I)  =wm1I  *(-d/denom)
           matVMAf(i_p_am1II) =wm1II *(-d/denom)
           matVMAb(i_p_ap1I)  =wp1I  *( d/denom)
           matVMAb(i_p_ap1II) =wp1II *( d/denom)
           matVMAb(i_p_am1I)  =wm1I  *(-d/denom)
           matVMAb(i_p_am1II) =wm1II *(-d/denom)
!!$           denom=7*(dalpha_ap1+dalpha_am1)-dalpha_ap2-dalpha_am2
!!$           matVMAf(i_p_ap2I)  =wp2I  *(  -d/denom)
!!$           matVMAf(i_p_ap2II) =wp2II *(  -d/denom)
!!$           matVMAf(i_p_ap2III)=wp2III*(  -d/denom)
!!$           matVMAf(i_p_ap2IV) =wp2IV *(  -d/denom)
!!$           matVMAf(i_p_ap1I)  =wp1I  *(+8*d/denom)
!!$           matVMAf(i_p_ap1II) =wp1II *(+8*d/denom)
!!$           matVMAf(i_p_am1I)  =wm1I  *(-8*d/denom)
!!$           matVMAf(i_p_am1II) =wm1II *(-8*d/denom)
!!$           matVMAf(i_p_am2I)  =wm2I  *(  +d/denom)
!!$           matVMAf(i_p_am2II) =wm2II *(  +d/denom)
!!$           matVMAf(i_p_am2III)=wm2III*(  +d/denom)
!!$           matVMAf(i_p_am2IV) =wm2IV *(  +d/denom)
!!$           matVMAb(i_p_ap2I)  =wp2I  *(  -d/denom)
!!$           matVMAb(i_p_ap2II) =wp2II *(  -d/denom)
!!$           matVMAb(i_p_ap2III)=wp2III*(  -d/denom)
!!$           matVMAb(i_p_ap2IV) =wp2IV *(  -d/denom)
!!$           matVMAb(i_p_ap1I)  =wp1I  *(+8*d/denom)
!!$           matVMAb(i_p_ap1II) =wp1II *(+8*d/denom)
!!$           matVMAb(i_p_am1I)  =wm1I  *(-8*d/denom)
!!$           matVMAb(i_p_am1II) =wm1II *(-8*d/denom)
!!$           matVMAb(i_p_am2I)  =wm2I  *(  +d/denom)
!!$           matVMAb(i_p_am2II) =wm2II *(  +d/denom)  
!!$           matVMAb(i_p_am2III)=wm2III*(  +d/denom)
!!$           matVMAb(i_p_am2IV) =wm2IV *(  +d/denom)  
        ELSE
           IF(d.LT.0) THEN
              IF(SECOND_ORDER_ALPHA.AND.i_p_ap2I.NE.0) THEN
                 denom=dalpha_ap1*dalpha_ap2*(dalpha_ap1+dalpha_ap2)
                 matVMAf(i_p_ap2I  )=wp2I  *(-d*dalpha_ap1*dalpha_ap1/denom)
                 matVMAf(i_p_ap2II )=wp2II *(-d*dalpha_ap1*dalpha_ap1/denom)
                 matVMAf(i_p_ap2III)=wp2III*(-d*dalpha_ap1*dalpha_ap1/denom)
                 matVMAf(i_p_ap2IV )=wp2IV *(-d*dalpha_ap1*dalpha_ap1/denom)
                 matVMAf(i_p_ap1I )=wp1I* (+d*(dalpha_ap1+dalpha_ap2)*(dalpha_ap1+dalpha_ap2)/denom)
                 matVMAf(i_p_ap1II)=wp1II*(+d*(dalpha_ap1+dalpha_ap2)*(dalpha_ap1+dalpha_ap2)/denom)
                 matVMAf(ipoint )=+d*&
                      & (dalpha_ap1*dalpha_ap1-(dalpha_ap1+dalpha_ap2)*(dalpha_ap1+dalpha_ap2))/denom
              ELSE
                 IF(NEW_DALPHA.AND.wp1I.LT.-100) THEN
                    IF(BI3oBI7.GT.0) THEN
                       dlambda=dlambda_lp1(ipoint)
                    ELSE
                       dlambda=dlambda_lm1(ipoint)
                    END IF
                    matVMAf(i_p_ap1I )= d/dalpha_ap1
                    matVMAf(i_p_ap1II)= d*ABS(BI3oBI7)/dlambda
                    matVMAf(ipoint )  =-d*(1./dalpha_ap1+ABS(BI3oBI7)/dlambda)
                    IF(CALC_DIFF.AND.nbif(ipoint).GT.0) THEN
                       sumb=d*ABS(BI3oBI7)/dlambda2*dlambda
                       j_p_lm1=i_p_lm1(ipoint)
                       j_p_lp1=i_p_lp1(1,ipoint)
                    END IF
                 ELSE
                    matVMAf(i_p_ap1I )=wp1I *(+d/dalpha_ap1)
                    matVMAf(i_p_ap1II)=matVMAf(i_p_ap1II)+wp1II*(+d/dalpha_ap1)
                    matVMAf(ipoint )         =-d/dalpha_ap1
                 END IF
!                 matVMAf(i_p_ap1I )=wp1I *(+d/dalpha_ap1)
!                 matVMAf(i_p_ap1II)=wp1II*(+d/dalpha_ap1)             
!                 matVMAf(ipoint )         =-d/dalpha_ap1
              END IF
              IF(SECOND_ORDER_ALPHA.AND.i_p_am2I.NE.0) THEN
                 denom=dalpha_am1*dalpha_am2*(dalpha_am1+dalpha_am2) 
                 matVMAb(i_p_am2I  )=wm2I  *(+d*dalpha_am1*dalpha_am1/denom)
                 matVMAb(i_p_am2II )=wm2II *(+d*dalpha_am1*dalpha_am1/denom)
                 matVMAb(i_p_am2III)=wm2III*(+d*dalpha_am1*dalpha_am1/denom)
                 matVMAb(i_p_am2IV )=wm2IV *(+d*dalpha_am1*dalpha_am1/denom)
                 matVMAb(i_p_am1I  )=wm1I  *(-d*(dalpha_am1+dalpha_am2)*(dalpha_am1+dalpha_am2)/denom)
                 matVMAb(i_p_am1II )=wm1II *(-d*(dalpha_am1+dalpha_am2)*(dalpha_am1+dalpha_am2)/denom)
                 matVMAb(ipoint )=-d*&
                      & (dalpha_am1*dalpha_am1-(dalpha_am1+dalpha_am2)*(dalpha_am1+dalpha_am2))/denom
              ELSE
                 IF(NEW_DALPHA.AND.wm1I.LT.-100) THEN
                    IF(BI3oBI7.GT.0) THEN
                       dlambda=dlambda_lm1(ipoint)
                    ELSE
                       dlambda=dlambda_lp1(ipoint)
                    END IF
                    matVMAb(i_p_am1I )=-d/dalpha_am1
                    matVMAb(i_p_am1II)=-d*ABS(BI3oBI7)/dlambda
                    matVMAb(ipoint )  = d*(1/dalpha_am1+ABS(BI3oBI7)/dlambda)
                    IF(CALC_DIFF.AND.nbif(ipoint).GT.0) THEN
                       sumb=d*ABS(BI3oBI7)/dlambda2*dlambda
                       j_p_lm1=i_p_lm1(ipoint)
                       j_p_lp1=i_p_lp1(1,ipoint)
                       matDIFb(j_p_lp1)=sumb
                       matDIFb(ipoint) =-2*sumb
                       matDIFb(j_p_lm1)=sumb
                    END IF
                 ELSE
                    matVMAb(i_p_am1I )=wm1I *(-d/dalpha_am1)
                    matVMAb(i_p_am1II)=matVMAb(i_p_am1II)+wm1II*(-d/dalpha_am1)
                    matVMAb(ipoint )         =+d/dalpha_am1
                 END IF
              END IF
!              matVMAb(i_p_am1I )=wm1I *(-d/dalpha_am1)
!              matVMAb(i_p_am1II)=wm1II*(-d/dalpha_am1)
!              matVMAb(ipoint )         =+d/dalpha_am1
           ELSE IF(d.GT.0) THEN
              IF(SECOND_ORDER_ALPHA.AND.i_p_am2I.NE.0) THEN
                 denom=dalpha_am1*dalpha_am2*(dalpha_am1+dalpha_am2)
                 matVMAf(i_p_am2I  )=wm2I  *(+d*dalpha_am1*dalpha_am1/denom)
                 matVMAf(i_p_am2II )=wm2II *(+d*dalpha_am1*dalpha_am1/denom)
                 matVMAf(i_p_am2III)=wm2III*(+d*dalpha_am1*dalpha_am1/denom)
                 matVMAf(i_p_am2IV )=wm2IV *(+d*dalpha_am1*dalpha_am1/denom)
                 matVMAf(i_p_am1I  )=wm1I  *(-d*(dalpha_am1+dalpha_am2)*(dalpha_am1+dalpha_am2)/denom)
                 matVMAf(i_p_am1II )=wm1II *(-d*(dalpha_am1+dalpha_am2)*(dalpha_am1+dalpha_am2)/denom)
                 matVMAf(ipoint )=-d*&
                      & (dalpha_am1*dalpha_am1-(dalpha_am1+dalpha_am2)*(dalpha_am1+dalpha_am2))/denom
              ELSE
                 IF(NEW_DALPHA.AND.wm1I.LT.-100) THEN
                    IF(BI3oBI7.GT.0) THEN
                       dlambda=dlambda_lm1(ipoint)
                    ELSE
                       dlambda=dlambda_lp1(ipoint)
                    END IF
                    matVMAf(i_p_am1I )=-d/dalpha_am1
                    matVMAf(i_p_am1II)=-d*ABS(BI3oBI7)/dlambda
                    matVMAf(ipoint )  = d*(1/dalpha_am1+ABS(BI3oBI7)/dlambda)
                    IF(CALC_DIFF.AND.nbif(ipoint).GT.0) THEN
                       sumb=d*ABS(BI3oBI7)/dlambda2*dlambda
                       j_p_lm1=i_p_lm1(ipoint)
                       j_p_lp1=i_p_lp1(1,ipoint)
                    END IF
                 ELSE
                    matVMAf(i_p_am1I )=wm1I *(-d/dalpha_am1)
                    matVMAf(i_p_am1II)=matVMAf(i_p_am1II)+wm1II*(-d/dalpha_am1)
                    matVMAf(ipoint )         =+d/dalpha_am1
                 END IF                    
!                 matVMAf(i_p_am1I )=wm1I *(-d/dalpha_am1)
!                 matVMAf(i_p_am1II)=wm1II*(-d/dalpha_am1)
!                 matVMAf(ipoint )         =+d/dalpha_am1
              END IF
              IF(SECOND_ORDER_ALPHA.AND.i_p_ap2I.NE.0) THEN
                 denom=dalpha_ap1*dalpha_ap2*(dalpha_ap1+dalpha_ap2)
                 matVMAb(i_p_ap2I  )=wp2I  *(-d*dalpha_ap1*dalpha_ap1/denom)
                 matVMAb(i_p_ap2II )=wp2II *(-d*dalpha_ap1*dalpha_ap1/denom)
                 matVMAb(i_p_ap2III)=wp2III*(-d*dalpha_ap1*dalpha_ap1/denom)
                 matVMAb(i_p_ap2IV )=wp2IV *(-d*dalpha_ap1*dalpha_ap1/denom)
                 matVMAb(i_p_ap1I )=wp1I* (+d*(dalpha_ap1+dalpha_ap2)*(dalpha_ap1+dalpha_ap2)/denom)
                 matVMAb(i_p_ap1II)=wp1II*(+d*(dalpha_ap1+dalpha_ap2)*(dalpha_ap1+dalpha_ap2)/denom)
                 matVMAb(ipoint )=+d*&
                      & (dalpha_ap1*dalpha_ap1-(dalpha_ap1+dalpha_ap2)*(dalpha_ap1+dalpha_ap2))/denom
              ELSE
                 IF(NEW_DALPHA.AND.wm1I.LT.-100) THEN
                    IF(BI3oBI7.GT.0) THEN
                       dlambda=dlambda_lm1(ipoint)
                    ELSE
                       dlambda=dlambda_lp1(ipoint)
                    END IF
                    matVMAb(i_p_am1I )=-d/dalpha_am1
                    matVMAb(i_p_am1II)=-d*ABS(BI3oBI7)/dlambda
                    matVMAb(ipoint )  = d*(1/dalpha_am1+ABS(BI3oBI7)/dlambda)
                    IF(CALC_DIFF.AND.nbif(ipoint).GT.0) THEN
                       sumb=d*ABS(BI3oBI7)/dlambda2*dlambda
                       j_p_lm1=i_p_lm1(ipoint)
                       j_p_lp1=i_p_lp1(1,ipoint)
                    END IF
                 ELSE
                    matVMAb(i_p_am1I )=wm1I *(-d/dalpha_am1)
                    matVMAb(i_p_am1II)=matVMAb(i_p_am1II)+wm1II*(-d/dalpha_am1)
                    matVMAb(ipoint )         =+d/dalpha_am1
                 END IF   
!                 matVMAb(i_p_ap1I )=wp1I *(+d/dalpha_ap1)
!                 matVMAb(i_p_ap1II)=wp1II*(+d/dalpha_ap1)             
!                 matVMAb(ipoint )         =-d/dalpha_ap1
              END IF

           END IF
        END IF
     END IF
  END IF

  IF(DEBUG) THEN
     WRITE(3150+myrank,'(10I6,5(1pe13.5),I4)') ipoint,i_p_am2,i_p_am1,i_p_ap1,i_p_ap2
     IF(ipoint.EQ.I0.OR.I0.EQ.0) THEN
        IF(ABS(matCOL(ipoint)).GT.ZERO) THEN
           DO jpoint=1,npoint
              WRITE(3200+myrank,'(2I6,5(1pe13.5),I4)') ipoint,jpoint,matCOL(jpoint),&
                   & matVEAb(jpoint),matVEAf(jpoint),matVMAb(jpoint),matVMAf(jpoint)
           END DO
        END IF
     ELSE
        WRITE(3200+myrank,'(2I6)') ipoint,nnz
     END IF
  END IF

END SUBROUTINE FILL_DKE_ROW


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#ifdef MPIandPETSc

SUBROUTINE FILL_MATRIX_PETSC(matCOL,jv,Epsi,matVEAf,matVEAb,matVMAf,matVMAb,ksp)

!-----------------------------------------------------------------------------------------------
!Use precalculated matrices matCOL, matVEAf, matVEAb, matVMAf, matVMAb and factors Epsi, nu(jv)
!and vdconst(jv) to create linear system ksp  
!-----------------------------------------------------------------------------------------------  

  USE GLOBAL
  USE petscsys
  USE petscksp
  IMPLICIT NONE
!#include <petsc/finclude/petscksp.h>
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!#include <petsc/finclude/petscpc.h>
  !Input
  INTEGER jv
  REAL*8 Epsi
  Mat matCOL,matVEAf,matVEAb,matVMAf,matVMAb
  !Output
  KSP ksp
  !Others
  PetscErrorCode ierr
  PetscScalar factor
  PetscReal, PARAMETER :: dtcol=0.5
  PetscReal, PARAMETER :: damping=1E-10
  Mat matA!,matA2
  PC pc
  INTEGER, PARAMETER :: MAXITS= 1000000 
  REAL*8, PARAMETER :: ATOL=1E-17
  REAL*8, PARAMETER :: RTOL =1E-12
  !Time
  CHARACTER*30, PARAMETER :: routine="FILL_MATRIX_PETSC"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)
!  CALL MatCreate(PETSC_COMM_WORLD,matA,ierr)
!  CALL MatSetSizes(matA,PETSC_DECIDE,PETSC_DECIDE,npoint,npoint,ierr)
!  CALL MatSetType(matA,MATAIJ,ierr)
!  CALL MatSeqAIJSetPreallocation(matA,PETSC_NULL_INTEGER,innz(0:npoint-1),ierr)
!  CALL MatSetup(matA,ierr)

  CALL MatDuplicate(matCOL,MAT_COPY_VALUES,matA,ierr)
!  CALL MatDuplicate(matCOL,MAT_DO_NOT_COPY_VALUES,matA,ierr) 
!  factor=nu(jv)
!  CALL MatAXPY(matA,factor,matCOL,SUBSET_NONZERO_PATTERN,IERR)
  factor=sgnB*Epsi/nu(jv)
  IF(sgnB*Epsi.LT.0) THEN
     CALL MatAXPY(matA,factor,matVEAb,SUBSET_NONZERO_PATTERN,IERR)
  ELSE IF(sgnB*Epsi.GT.0) THEN
     CALL MatAXPY(matA,factor,matVEAf,SUBSET_NONZERO_PATTERN,IERR)
  END IF
  IF(TANG_VM) THEN
     factor=vmconst(jv)/nu(jv)
     IF(vmconst(jv).LT.0) THEN
        CALL MatAXPY(matA,factor,matVMAf,SUBSET_NONZERO_PATTERN,ierr)
     ELSE
        CALL MatAXPY(matA,factor,matVMAb,SUBSET_NONZERO_PATTERN,ierr)
     END IF
  END IF

  CALL KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
  CALL KSPSetTolerances(ksp,RTOL,ATOL,PETSC_DEFAULT_REAL,MAXITS,ierr)

!  CALL KSPSetType(ksp,KSPGMRES,ierr)
!  CALL KSPGetPC(ksp,pc,ierr)
!  CALL PCSetType(pc,PCJACOBI,ierr)

!  CALL KSPSetType(ksp,KSPPREONLY,ierr)
!  CALL KSPGetPC(ksp,pc,ierr)
!  CALL PCSetType(pc,PCILU,ierr)
!  CALL PCFactorSetLevels(pc,10,ierr)
!  CALL PCFactorSetMatOrderingType(pc,"natural",ierr)

  CALL KSPSetType(ksp,KSPPREONLY,ierr)
  CALL KSPGetPC(ksp,pc,ierr)
  CALL PCSetType(pc,PCLU,ierr)
  CALL PCFactorSetMatOrderingType(pc,"nd",ierr)
#if defined(PETSC_HAVE_SUPERLU)
  CALL PCFactorSetMatSolverType(pc,MATSOLVERSUPERLU,ierr)
#endif
!  CALL PCFactorReorderForNonzeroDiagonal(pc,dtcol,ierr)
!  CALL PCFactorSetShiftAmount(pc,dtcol,ierr)
!  CALL PCFactorSetShiftNonzero(pc,damping,ierr)
!  CALL PCFactorSetReuseOrdering(pc,PETSC_TRUE,ierr)
!  CALL PCFactorSetReuseFill(pc,PETSC_TRUE,ierr)
!  CALL PCFactorSetColumnPivot(pc,dtcol,ierr)
!!#ifdef PETSC_HAVE_MUMPS
!!     call PCFactorSetMatSolverType(pc,MATSOLVERMUMPS,ierr)
!!     call PCFactorSetUpMatSolverType(pc,ierr)
!!#endif
  CALL KSPSetOperators(ksp,matA,matA,ierr)
  CALL KSPSetUp(ksp,ierr)
  CALL MatDestroy(matA,ierr)

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE FILL_MATRIX_PETSC


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#else 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE FILL_MATRIX(npoint,COL,jv,Epsi,VEAf,VEAb,VMAf,VMAb,mat)

!-----------------------------------------------------------------------------------------------
!Use precalculated matrices COL, VEAf, VEAb, VMAf, VMAb and factors Epsi, nu(jv) and vmconst(jv)
!to create matrix mat
!-----------------------------------------------------------------------------------------------  

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER jv,npoint
  REAL*8 Epsi,COL(npoint,npoint)
  REAL*8 VEAf(npoint,npoint),VEAb(npoint,npoint),VMAf(npoint,npoint),VMAb(npoint,npoint)
  !Output
  REAL*8 mat(npoint,npoint)
  !Time
  CHARACTER*30, PARAMETER :: routine="FILL_MATRIX"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)
  
  mat=COL
  IF(sgnB*Epsi.LT.0) THEN
     mat=mat+VEAb*sgnB*Epsi/nu(jv)
  ELSE IF(sgnB*Epsi.GT.0) THEN
     mat=mat+VEAf*sgnB*Epsi/nu(jv)
  END IF
  IF(TANG_VM) THEN
     IF(vmconst(jv).LT.0) THEN
        mat=mat+VMAf*vmconst(jv)/nu(jv)
     ELSE
        mat=mat+VMAb*vmconst(jv)/nu(jv)
     END IF
  END IF
  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE FILL_MATRIX

#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#ifdef MPIandPETSc

SUBROUTINE INVERT_MATRIX_PETSC(nalphab,jv,npoint,mata,BI3,BI8,factnu,phi1c,ksp,gint)

!-----------------------------------------------------------------------------------------------
!Fill rhs with some linear combination (depending on phi1c,jv...) of arrays BI3 and BI8 of size
!npoint and invert linear system ksp to obtain g. Use nalphab to skip some helicities in phi1c
!-----------------------------------------------------------------------------------------------  

  USE GLOBAL
  use petscsys
  use petscksp
  IMPLICIT NONE
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!#include <petsc/finclude/petscpc.h>
!#include <petsc/finclude/petscksp.h>
!#include <petsc/finclude/petscviewer.h>
  !Input
  INTEGER nalphab,jv,npoint
  REAL*8 BI3(npoint),BI8(npoint,Nnmp),factnu(npoint),phi1c(Nnmp)
  Mat matA
  KSP ksp
  !Output
  REAL*8 gint(npoint,Nnmp)
  !Others
  CHARACTER*100 serr
  INTEGER ii,irhs!,jrhs
  PetscScalar c(npoint),g(npoint),g2(npoint)
  PetscErrorCode ierr
  PetscInt indx(npoint),npoint_ps
  Vec vecb,vecx
  Vec vecb2
  !Time
  CHARACTER*30, PARAMETER :: routine="INVERT_MATRIX"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart
  !Time
  CHARACTER*30, PARAMETER :: routine2="INVERT_MATRIX2"
  INTEGER, SAVE :: ntotal2=0
  REAL*8,  SAVE :: ttotal2=0
  REAL*8,  SAVE :: t02=0
  REAL*8 tstart2

  CALL CPU_TIME(tstart)

  npoint_ps=npoint
  phi1c=phi1c !To be removed

  DO ii=1,npoint
     indx(ii)=ii-1
  END DO
  CALL VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,npoint_ps,vecb2,ierr)
  CALL VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,npoint_ps,vecb,ierr)
  CALL VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,npoint_ps,vecx,ierr)
  
  !Scan in possible rhs of the DKE·
  !(corresponding to different contributions the the radial ExB from varphi1)
  DO irhs=1,Nnmp  !Skip helicities too large for the alpha precision
     IF(irhs.GT.1.AND.(ABS(ext_np(irhs)).GT.nalphab/nsamp.OR.ABS(ext_mp(irhs)).GT.nalphab/nsamp))&
          & CYCLE
     CALL VecZeroEntries(vecb2,ierr)
     CALL VecZeroEntries(vecb,ierr)
     CALL VecZeroEntries(vecx,ierr)
     !Fill rhs and invert
     IF(irhs.EQ.1) THEN
        c=BI3  !radial magnetic drift
!        IF(PHI1_READ) THEN
!           DO jrhs=2,Nnmp !total radial ExB drift added to the magnetic drift
!              c=c+phi1c(jrhs)*BI8(:,jrhs)*2*borbic(0,0)/vdconst(jv)
!           END DO
!        END IF
     ELSE !radial ExB drift
        c=BI8(:,irhs)/vdconst(jv)
     END IF
     CALL VecSetValues(vecb,npoint_ps,indx,c,INSERT_VALUES,ierr)      
     CALL VecAssemblyBegin(vecb,ierr)
     CALL VecAssemblyEnd(vecb,ierr)
     !Solve
     CALL CPU_TIME(tstart2)
     CALL KSPSolve(ksp,vecb,vecx,ierr)
     CALL CALCULATE_TIME(routine2,ntotal2,t02,tstart2,ttotal2)
     g=0
     IF(.NOT.CALC_DG) THEN
        CALL VecGetValues(vecx,npoint_ps,indx,g,ierr)
     ELSE IF(FLUX_NU) THEN
        CALL MatMult(matA,vecx,vecb2,ierr)
        CALL VecGetValues(vecb2,npoint_ps,indx,g,ierr)
     ELSE IF(CALC_RHS.OR.CALC_DA.OR.CALC_DIFF.OR.CALC_COL) THEN
        CALL MatMult(matA,vecx,vecb2,ierr)
        CALL VecGetValues(vecb2,npoint_ps,indx,g2,ierr)
        IF(CALC_DA) THEN
           CALL VecGetValues(vecx,npoint_ps,indx,g,ierr)
           g=g2/g
        ELSE
           g=g2
        END IF
     END IF
     
!     CALL KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)
!     CALL PetscMemoryGetCurrentUsage(ierr)
!     CALL PetscViewerASCIIOpen(MPI_COMM_WORLD,'filename.xml',viewer,ierr)
!     CALL PetscViewerPushFormat(viewer,PETSC_VIEWER_DEFAULT,ierr)
!     CALL PetscLogView(viewer,ierr)
     !Distribution function at each point
     gint(:,irhs)=factnu*g*vdconst(jv)/nu(jv)     
!     gint(:,irhs+1)=gint(:,irhs+1)+g*Sdke(jv)*weight(jv)
     IF(ierr.NE.0) THEN
        serr="Error when inverting the DKE"
        CALL END_ALL(serr,.FALSE.)
     END IF
     
     IF(.NOT.SOLVE_QN) EXIT
  END DO
  IF(ESCOTO.OR.ABS(BI3(npoint)).GT.0) CALL KSPDestroy(ksp,ierr)

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)
  
END SUBROUTINE INVERT_MATRIX_PETSC


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#else 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE INVERT_MATRIX(nalphab,jv,npoint,BI3,BI8,phi1c,mat,gint)

!-----------------------------------------------------------------------------------------------
!Fill rhs with some linear combination (depending on phi1c,jv...) of arrays BI3 and BI8 of size
!npoint and invert linear system with matrix mat to obtain g.
!Use nalphab to skip some helicities in phi1c
!-----------------------------------------------------------------------------------------------  

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER nalphab,jv,npoint
  REAL*8 BI3(npoint),BI8(npoint,Nnmp),phi1c(Nnmp),mat(npoint,npoint)
  !Output
  REAL*8 gint(npoint,Nnmp)
  !Others
  LOGICAL ierr
  CHARACTER*100 serr
  INTEGER irhs!,jrhs
  REAL*8 rhs(npoint),ipivot(npoint,npoint),g(npoint)
  !Time
  CHARACTER*30, PARAMETER :: routine="INVERT_MATRIX"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  phi1c=phi1c !To be removed
    
  DO irhs=1,Nnmp  !Skip helicities too large for the alpha precision
     IF(irhs.GT.1.AND.(ABS(ext_np(irhs)).GT.nalphab/nsamp.OR.ABS(ext_mp(irhs)).GT.nalphab/nsamp))&
          & CYCLE
     IF(irhs.EQ.1) THEN
        rhs=BI3 !radial magnetic drift
!        IF(PHI1_READ) THEN
!           DO jrhs=2,Nnmp !total radial ExB drift added to the magnetic drift
!              rhs=rhs+phi1c(jrhs)*BI8(:,jrhs)*2*borbic(0,0)/vdconst(jv)
!           END DO
!        END IF
     ELSE !radial ExB drift
        rhs=BI8(:,irhs)/vdconst(jv)
     END IF
     !Solve
     CALL DGETRF(npoint,npoint,mat,npoint,ipivot,ierr) 
     CALL DGETRS('No transpose',npoint,1,mat,npoint,ipivot,rhs,npoint,ierr)
     g=rhs
     !Distribution function at each point
     gint(:,irhs)=g*vdconst(jv)/nu(jv)     
!     gint(:,irhs+1)=gint(:,irhs+1)+g*Sdke(jv)*weight(jv)
     IF(ierr) THEN
        serr="Error when inverting the DKE"
        CALL END_ALL(serr,.FALSE.)
     END IF
     IF(.NOT.SOLVE_QN) EXIT
  END DO

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE INVERT_MATRIX

#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE INTEGRATE_G(nalpha,nalphab,nlambda,lambda,i_p,npoint,g,notg,&
     & zetap,thetap,theta,B_al,vds_al,D11,dn1,dn1nm)

!-----------------------------------------------------------------------------------------------
!Calculate contribution D11, dn1 and dn1nm to the flux and to quasineutrality,
! by integrating in lambda grid of nlambda size the distribution function g known at npoints, i_p.
! The integral is calculated in (zetap,thetap) grid of size nalphabxnalpha using precalculated values
! of B_al and vds_al, and then interpolated to a (zetap,thetap) square grid, where flux-surface 
!average is done for the flux
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  USE KNOSOS_STELLOPT_MOD  
  IMPLICIT NONE
  !Input
  LOGICAL notg
  INTEGER nalpha,nalphab,nlambda,i_p(nlambda,nalpha,nalphab),npoint
  REAL*8 zetap(nalphab),thetap(nalpha,nalphab),theta(nalphab)
  REAL*8 lambda(nlambda),g(npoint,Nnmp),B_al(nalpha,nalphab),vds_al(Nnmp,nalpha,nalphab)
  !Output
  REAL*8 D11(Nnmp,Nnmp),dn1(nalphab,nalphab),dn1nm(Nnmp,Nnmp)
  !Others
  REAL*8, PARAMETER :: F5o12 =0.416666666666667
  REAL*8, PARAMETER :: F13o12=1.083333333333333
  INTEGER ia,il,ial,ila,jla,kla,ipoint,jpoint,kpoint,ind,nm
  REAL*8 fdlambda,fhdlambda2,lambdaB,d3vdlambdadK,dlambda,fint
  REAL*8 D11_alp(nalpha*nalphab),D11_ale(3*nalpha),D11_zt(nalphab,nalphab)
  REAL*8 dn1_alp(nalpha*nalphab),dn1_ale(3*nalpha)
  REAL*8 one_o_modB,modB,lambda1,FSA,sqrt1mlB,dummy
  REAL*8 Jac(nalphab,nalphab),dn1_zt(nalphab,nalphab),dn1c_zt(nalphab,nalphab)
  COMPLEX*16 dn1c_nm(nalphab,nalphab)
  INTEGER, SAVE :: tnalpha
  REAL*8, SAVE, ALLOCATABLE :: D11p(:,:),dn1nmp(:,:)
  REAL*8, ALLOCATABLE :: thetape(:,:),vds_zt(:,:,:)
  !Time
  CHARACTER*30, PARAMETER :: routine="INTEGRATE_G"
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
     
  PRE_INTV=(SOLVE_AMB.OR.(SOLVE_QN.AND..NOT.PHI1_READ.AND..NOT.ZERO_PHI1).OR.NERR.GT.1).AND.CONVERGED
  D11=0
  dn1=0
  dn1nm=0
  IF(DEBUG) THEN
     DO ipoint=1,npoint  
        WRITE(3500+myrank,'(I6,1(1pe13.5),2I6)') ipoint,g(ipoint,1),nalpha,nlambda
     END DO
     CALL FLUSH(3500+myrank)
  END IF

  IF(CALCULATED_INT.AND.PRE_INTV) GOTO 123
#ifndef NAG
  IF(plan_fwd.NE.0) THEN
     CALL DFFTW_DESTROY_PLAN(plan_fwd)
     plan_fwd=0
  END IF
#endif

  IF(ALLOCATED(D11p)) THEN
     DEALLOCATE(D11p)
     IF(QN) DEALLOCATE(dn1nmp)
  END IF
  tnalpha=2*nalpha
  IF(aiota/nzperiod.GE.1) tnalpha=3*nalpha
  ALLOCATE(thetape(tnalpha,nalphab))
  IF(PRE_INTV) THEN
     ALLOCATE(D11p(Nnmp,npoint))
     IF(QN) ALLOCATE(dn1nmp(Nnmp,npoint))
  ELSE
     ALLOCATE(D11p(Nnmp,1))
     IF(QN) ALLOCATE(dn1nmp(Nnmp,1))
  END IF

  D11p=0
  IF(QN) dn1nmp=0
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
  
  !Precalculate quantities
  lambda1=lambda(1)
  dlambda=lambda(2)-lambda1
     ALLOCATE(vds_zt(Nnmp,nalphab,nalphab))
     DO ia=1,nalphab
        DO il=1,nalphab
           CALL FILL_BNODE(zetap(il),theta(ia),Jac(ia,il),dummy,vds_zt(:,ia,il),.FALSE.)
           IF(USE_B0) CALL FILL_BNODE(zetap(il),theta(ia),dummy,dummy,vds_zt(:,ia,il),.TRUE.)
        END DO
     END DO

  DO ipoint=1,npoint

     D11_alp=0
     IF(QN) dn1_alp=0
     !Scan in the flux surface
     DO ia=1,nalpha 
        IF(.NOT.PRE_INTV.AND.ipoint.GT.1) EXIT
        !Scan in the flux surface
        DO il=1,nalphab 
           ial=(il-1)*nalpha+ia
           kla=-1
           modB=B_al(ia,il)
           one_o_modB=1./modB
           jla=INT((one_o_modB-lambda1)/dlambda)        
           IF(jla.LE.1) jla=2
           IF(jla.NE.0) fdlambda=(one_o_modB-lambda(jla))/dlambda
           kpoint=1
           !Integral in lambda
           DO ila=jla,2,-1
              jpoint=i_p(ila,ia,il)
              fint=1
              IF(jpoint.LE.1) THEN 
                 jpoint=kpoint
                 IF(jpoint.LE.1) CYCLE
                 fint=(ila-1.)/(kla-1.)
              ELSE
                 kpoint=jpoint
                 kla=ila
              END IF
              IF(PRE_INTV.AND.jpoint.NE.ipoint) CYCLE
!              IF(DEBUG.AND.ipoint.EQ.1) WRITE(3453+myrank,'(5(1pe13.5),10I6)') &
!                   & zetap(il),thetap(ia,il),lambda(1),g(1,1)-g(1,1)              
              fhdlambda2=0.5*fdlambda*fdlambda
              IF(ila.EQ.jla) THEN
                 fint=fint*(F5o12+fdlambda+fhdlambda2)
              ELSE IF(ila.EQ.jla-1) THEN
                 fint=fint*(F13o12-fhdlambda2)
              ELSE IF (ila.EQ.2) THEN
                 fint=fint*F13o12
              END IF
              lambdaB=lambda(ila)*modB
              IF(lambdaB.GT.1) CYCLE
              sqrt1mlB=SQRT(1.-lambdaB)
              d3vdlambdadK=fint*modB/sqrt1mlB*dlambda
              IF(PRE_INTV) THEN   !Calculate the linear contribution of g(ipoint)
                 IF(notg) THEN
                    D11_alp(ial)=D11_alp(ial)+d3vdlambdadK
                 ELSE
                    D11_alp(ial)=D11_alp(ial)-d3vdlambdadK*vds_al(1,ia,il)*(1.0-0.5*lambdaB) /2.
                    IF(QN) dn1_alp(ial)=dn1_alp(ial)+d3vdlambdadK
                 END IF
              ELSE                !Accumulate the total contribution g for each point jpoint
                 IF(notg) THEN
                    D11_alp(ial)=D11_alp(ial)+g(jpoint,1)*d3vdlambdadK
                 ELSE
                    IF(.NOT.PHI1_READ.OR..NOT.IMP1NU) THEN
                       D11_alp(ial)=D11_alp(ial)-&
                            & g(jpoint,1)*d3vdlambdadK*vds_al(1,ia,il)*(1.0-0.5*lambdaB)/2.
                    ELSE
                       D11_alp(ial)=D11_alp(ial)+g(jpoint,1)*&
                            & modB/(lambda(ila)*lambda(ila)*lambda(ila)*sqrt1mlb)*vds_al(1,ia,il)
                    END IF
                    IF(QN) dn1_alp(ial)=dn1_alp(ial)+g(jpoint,1)*d3vdlambdadK        
                 END IF
              END IF
           END DO

!           IF(DEBUG.AND.ipoint.EQ.1) WRITE(3451+myrank,'(5(1pe13.5),10I6)') &
!                & il,ia,D11_alp(ial),vds_al(1,ia,il)
 !           IF(DEBUG.AND.ipoint.EQ.1) WRITE(3451+myrank,'(5(1pe13.5),10I6)') &
 !                & zetap(il),thetap(ia,il),dn1_alp(ial),D11_alp(ial),vds_al(1,ia,il)
        END DO
     END DO

     !Copy values to extended grid and interpolate to square grid
     dn1_zt=0
     DO il=1,nalphab
        ial=(il-1)*nalpha
        IF(aiota/nzperiod.LT.1) THEN
           D11_ale(        1: nalpha)=D11_alp(ial+1:ial+nalpha)
           D11_ale( nalpha+1:tnalpha)=D11_alp(ial+1:ial+nalpha)
           IF(QN) THEN
              dn1_ale(        1: nalpha)=dn1_alp(ial+1:ial+nalpha)
              dn1_ale( nalpha+1:tnalpha)=dn1_alp(ial+1:ial+nalpha)
           END IF
        ELSE
           D11_ale(         1:  nalpha)=D11_alp(ial+1:ial+nalpha)
           D11_ale(  nalpha+1:2*nalpha)=D11_alp(ial+1:ial+nalpha)
           D11_ale(2*nalpha+1: tnalpha)=D11_alp(ial+1:ial+nalpha)
           IF(QN) THEN
              dn1_ale(         1:  nalpha)=dn1_alp(ial+1:ial+nalpha)
              dn1_ale( nalpha+1: 2*nalpha)=dn1_alp(ial+1:ial+nalpha)
              dn1_ale(2*nalpha+1: tnalpha)=dn1_alp(ial+1:ial+nalpha)
           END IF
        END IF
        
        !Interpolation to square grid
        DO ia=1,nalphab
           ind=(ia-1)*nalphab+il
           D11_zt(ia,il)=0
           CALL LAGRANGE(thetape(1:tnalpha,il),D11_ale(1:tnalpha),tnalpha,&
                & theta(ia),D11_zt(ia,il),2)
           IF(QN) THEN!.OR.NTV)THEN
              dn1_zt(ia,il)=0
              CALL LAGRANGE(thetape(1:tnalpha,il),dn1_ale(1:tnalpha),tnalpha,&
                 & theta(ia),dn1_zt(ia,il),2)
           END IF
 !           IF(DEBUG.AND.ipoint.EQ.1) WRITE(3450+myrank,'(10(1pe13.5),10I6)') &
 !                & zetap(il),theta(ia),dn1_zt(ia,il),D11_zt(ia,il)
        END DO
     END DO
     !Flux surface average
     D11p(1,ipoint)=FSA(nalphab,nalphab,D11_zt,Jac,1)
     IF(QN) THEN
        DO nm=2,Nnmp
           IF(ONLY_PHI1) EXIT 
           D11p(nm,ipoint)=FSA(nalphab,nalphab,dn1_zt*2*vds_zt(nm,:,:),Jac,1)
        END DO
        DO ia=1,nalphab
           DO il=1,nalphab
              dn1c_zt(il,ia)=dn1_zt(ia,il)              
           END DO
        END DO
        CALL FFTF_KN(nalphab,dn1c_zt,dn1c_nm)
        CALL FILL_ORBI(nalphab,nalphab,dn1c_nm,Nnmp,dn1nmp(:,ipoint))
     END IF
     IF(.NOT.PRE_INTV.AND.ipoint.EQ.1) EXIT !JL
  END DO

123 IF(PRE_INTV) THEN

     DO nm=1,Nnmp
        IF(QN) CALL DGEMV('n',Nnmp,npoint,ONE,&
             & dn1nmp(:,:),Nnmp,g(:,nm),1,ZERO,dn1nm(:,nm),1)
        IF(.NOT.ONLY_PHI1) THEN
           CALL DGEMV('n',Nnmp,npoint,ONE,&
                & D11p(:,:),Nnmp,g(:,nm),1,ZERO,D11(:,nm),1)
        ELSE IF(nm.EQ.1) THEN
           CALL DGEMV('n',1,npoint,ONE,&
             & D11p(1,:),1,g(:,1),1,ZERO,D11(1,1),1)
        END IF
     END DO
!     CALL DGEMM('n','n',Nnmp,Nnmp,npoint-1,ONE,&
!          & D11p(:,1:npoint-1)  ,Nnmp,g(1:npoint-1,:),npoint-1,ZERO,D11,Nnmp)
!     IF(QN) CALL DGEMM('n','n',Nnmp,Nnmp,npoint-1,ONE,&
!          & dn1nmp(:,1:npoint-1),Nnmp,g(1:npoint-1,:),npoint-1,ZERO,dn1nm,Nnmp) 
  ELSE
     D11(1,1)=D11p(1,1)
     IF(QN) THEN
        dn1=dn1_zt
        dn1nm(:,1)=dn1nmp(:,1)
     END IF
 END IF

 CALCULATED_INT=.TRUE.

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE INTEGRATE_G


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE INTEGRATE_G_NEW(nalpha,nalphab,nlambda,lambda,i_p,npoint,g,notg,&
     & thetap,B_al,vds_al,D11)

!-----------------------------------------------------------------------------------------------
!Calculate contribution D11 to the flux
! by integrating in lambda grid of nlambda size the distribution function g known at npoints, i_p.
! The integral is calculated in (zetap,thetap) grid of size nalphabxnalpha using precalculated values
! of B_al and vds_al, and then interpolated to a (zetap,thetap) square grid, where flux-surface 
!average is done for the flux
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  USE KNOSOS_STELLOPT_MOD  
  IMPLICIT NONE
  !Input
  LOGICAL notg
  INTEGER nalpha,nalphab,nlambda,i_p(nlambda,nalpha,nalphab),npoint
  REAL*8 lambda(nlambda),g(npoint),thetap(nalpha,nalphab),B_al(nalpha,nalphab),vds_al(nalpha,nalphab)
  !Output
  REAL*8 D11
  !Others
  REAL*8, PARAMETER :: F5o12 =0.416666666666667
  REAL*8, PARAMETER :: F7o12 =0.583333333333333
  REAL*8, PARAMETER :: F13o12=1.083333333333333
  REAL*8, PARAMETER :: F23o12=1.91666666666667
  INTEGER ia,il,ila,ipoint,nla
  REAL*8 lambdaB,d3vdlambdadK,dlambda,fdlambda,one_o_modB,modB,lambda1,sqrt1mlB
  REAL*8 g1oB,D11_alp(nalpha,nalphab),dD11,doffset,offset,FSA2
  !Time
  CHARACTER*30, PARAMETER :: routine="INTEGRATE_G_NEW"
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
     
  D11=0
  IF(DEBUG) THEN
     DO ipoint=1,npoint  
        WRITE(3500+myrank,'(I6,1(1pe13.5),2I6)') ipoint,g(ipoint),nalpha,nlambda
     END DO
     CALL FLUSH(3500+myrank)
  END IF
  
  !Precalculate quantities
  lambda1=lambda(1)
  dlambda=lambda(2)-lambda1
  
  D11_alp=0
  !Scan in the flux surface
  DO ia=1,nalpha 
     DO il=1,nalphab 
        modB=B_al(ia,il)
        one_o_modB=1./modB
        DO ila=1,nlambda
           IF(lambda(ila).GT.one_o_modB) THEN
              nla=ila-1
              fdlambda=(one_o_modB-lambda(nla))/dlambda
              EXIT
           END IF
        END DO
        IF(nla.EQ.0) CYCLE
        g1oB=g(i_p(nla,ia,il))*(1+fdlambda)-g(i_p(nla-1,ia,il))*fdlambda
        IF(notg) THEN
           offset=-2*g1oB*SQRT(1-lambda(1)*modB)
        ELSE
           offset=-2*g1oB*vds_al(ia,il)*0.25*SQRT(1-lambda(1)*modB)
        END IF
        DO ila=1,nla
           lambdaB=lambda(ila)*modB
           sqrt1mlB=SQRT(1.-lambdaB)
           d3vdlambdadK=modB/sqrt1mlB
           ipoint=i_p(ila,ia,il)
           IF(notg) THEN
              dD11=-g(ipoint)*d3vdlambdadK
              doffset=-g1oB*d3vdlambdadK
           ELSE
              dD11=-g(ipoint)*(vds_al(ia,il)*(1.0-0.5*lambdaB)/2.)*d3vdlambdadK
              doffset=-g1oB*vds_al(ia,il)*0.25*d3vdlambdadK
           END IF
           IF(ila.EQ.1) THEN
              dD11   =dD11/2.
              doffset=doffset/2.
           ELSE IF(ila.EQ.nla) THEN
              dD11   =dD11   *(1.+fdlambda)/2.
              doffset=doffset*(1.+fdlambda)/2.
           END IF
           D11_alp(ia,il)=D11_alp(ia,il)+dD11-doffset
        END DO
        D11_alp(ia,il)=D11_alp(ia,il)*dlambda+offset
        
        IF(DEBUG) WRITE(3451+myrank,'(2I6,5(1pe13.5),10I6)') &
             & il,ia,D11_alp(ia,il),vds_al(ia,il)
     END DO
  END DO
  !Flux surface average
  D11=FSA2(nalpha,nalphab,thetap(:,1),D11_alp,aiBtpBz/B_al/B_al,1)
  CALCULATED_INT=.TRUE.

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE INTEGRATE_G_NEW


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE FILL_ORBI(numz,numt,qnmc_in,num,qnm_out)

!-------------------------------------------------------------------------------------------------
!Read qnmc_in(numz,numt) from FFT and write num Fourier modes in 'ddkes2.data' format
!-------------------------------------------------------------------------------------------------

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER numz,numt,num
  COMPLEX*16 qnmc_in(numz,numt)
  !Output
  REAL*8 qnm_out(num)
!  REAL qnm_out(num)
  !Others
  INTEGER n,m,nm,ind1,ind2,ind3,ind4
  REAL*8 qorbic(-ntorb:ntorb,0:mpolb),qorbis(-ntorb:ntorb,0:mpolb)
  REAL*8 qorbic_temp(-ntorb:ntorb,-mpolb:mpolb),qorbis_temp(-ntorb:ntorb,-mpolb:mpolb)
  
  qorbic_temp=0
  qorbis_temp=0
  qorbic    =0
  qorbis    =0
  
  DO m=-mpolb,mpolb
     DO n=-ntorb,ntorb
        IF(n.EQ.0.OR.m.EQ.0) THEN 
           IF(m.EQ.0.AND.n.EQ.0) THEN
              qorbic_temp(0,0)=REAL(qnmc_in(1,1))
              qorbis_temp(0,0)=AIMAG(qnmc_in(1,1))
           ELSE IF(n.EQ.0) THEN
              IF(m.GT.0) THEN
                 ind1=m+1
                 ind2=numt-m+1
                 IF(ind1.GT.numt.OR.ind2.LT.1) CYCLE
                 qorbic_temp(0,m)= REAL( qnmc_in(1,ind1)+qnmc_in(1,ind2))
                 qorbis_temp(0,m)=AIMAG(-qnmc_in(1,ind1)+qnmc_in(1,ind2))
              ELSE
                 ind1=1-m
                 ind2=numt+m+1
                 IF(ind1.GT.numt.OR.ind2.LT.1) CYCLE
                 qorbic_temp(0,-m)= REAL( qnmc_in(1,ind1)+qnmc_in(1,ind2))
                 qorbis_temp(0,-m)=AIMAG(-qnmc_in(1,ind1)+qnmc_in(1,ind2))
              END IF
           ELSE IF(m.EQ.0) THEN
              IF(n.GT.0) THEN
                 ind1=n+1
                 ind2=numz-n+1
                 IF(ind1.GT.numz.OR.ind2.LT.1) CYCLE
                 qorbic_temp(n,0)=qorbic_temp(n,0) +REAL( qnmc_in(ind1,1)-qnmc_in(ind2,1))
                 qorbis_temp(n,0)=qorbis_temp(n,0)+AIMAG(-qnmc_in(ind1,1)+qnmc_in(ind2,1))/2
              ELSE
                 ind1=1-n
                 ind2=numz+n+1
                 IF(ind1.GT.numz.OR.ind2.LT.1) CYCLE
                 qorbic_temp(-n,0)=qorbic_temp(-n,0)+ REAL( qnmc_in(ind1,1)+qnmc_in(ind2,1))
                 qorbis_temp(-n,0)=qorbis_temp(-n,0)+AIMAG(-qnmc_in(ind1,1)+qnmc_in(ind2,1))/2
              END IF
           END IF
        ELSE IF(n.GT.0) THEN
           IF(m.GT.0) THEN
              ind1=n+1
              ind2=m+1
              ind3=numz-n+1
              ind4=numt-m+1
              IF(ind1.GT.numz.OR.ind2.GT.numt.OR.ind3.LT.1.OR.ind4.LT.1) CYCLE
              qorbic_temp(n,m)= REAL(+qnmc_in(ind1,ind2)+qnmc_in(ind3,ind4))/2
              qorbis_temp(n,m)=AIMAG(-qnmc_in(ind1,ind2)+qnmc_in(ind3,ind4))
           ELSE
              ind1=numz-n+1
              ind2=1-m
              ind3=n+1
              ind4=numt+m+1
              IF(ind1.LT.1.OR.ind2.GT.numt.OR.ind3.GT.numz.OR.ind4.LT.1) CYCLE
              qorbic_temp(n,m)= REAL(+qnmc_in(ind1,ind2)+qnmc_in(ind3,ind4))/2
              qorbis_temp(n,m)=AIMAG(+qnmc_in(ind1,ind2)-qnmc_in(ind3,ind4))
           END IF
        ELSE IF(n.LT.0) THEN
           IF(m.GT.0) THEN
              ind1=1-n
              ind2=numt-m+1
              ind3=numz+n+1
              ind4=m+1
              IF(ind1.GT.numz.OR.ind2.LT.1.OR.ind3.LT.1.OR.ind4.GT.numt) CYCLE
              qorbic_temp(n,m)= REAL(+qnmc_in(ind1,ind2)+qnmc_in(ind3,ind4))/2
              qorbis_temp(n,m)=AIMAG(+qnmc_in(ind1,ind2)-qnmc_in(ind3,ind4))
           ELSE
              ind1=numz+n+1
              ind2=numt+m+1
              ind3=1-n
              ind4=1-m
              IF(ind1.LT.1.OR.ind2.LT.1.OR.ind3.GT.numz.OR.ind4.GT.numt) CYCLE
              qorbic_temp(n,m)= REAL(+qnmc_in(ind1,ind2)+qnmc_in(ind3,ind4))/2
              qorbis_temp(n,m)=AIMAG(+qnmc_in(ind1,ind2)-qnmc_in(ind3,ind4))
           END IF
        END IF
     END DO
  END DO

  DO m=-mpolb,mpolb
     DO n=-ntorb,ntorb
        IF(m.GE.0) THEN
           qorbic(n,m)=qorbic(n,m)+qorbic_temp(n,m)
           qorbis(n,m)=qorbis(n,m)+qorbis_temp(n,m)
        ELSE
           qorbic(-n,-m)=qorbic(-n,-m)+qorbic_temp(n,m)
!           qorbis(-n,-m)=qorbis(-n,-m)-qorbis_temp(n,m)
        END IF
     END DO
  END DO

  nm=0
  DO m=0,mpolb
     DO n=-ntorb,ntorb
        IF(m.EQ.0.AND.n.LT.0) CYCLE
        nm=nm+1
        qnm_out(nm)=qorbic(n,m)
     END DO
  END DO
  DO m=0,mpolb
     DO n=-ntorb,ntorb
        IF(m.EQ.0.AND.n.LT.0) CYCLE
        nm=nm+1
        qnm_out(nm)=qorbis(n,m)
     END DO
  END DO

END SUBROUTINE FILL_ORBI


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



