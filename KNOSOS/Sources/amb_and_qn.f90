

!Solve the drift kinetic equation, together with quasineutrality and ambipolarity


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE SOLVE_DKE_QN_AMB(it,NBB,ZB,AB,REGB,S,nb,dnbdpsi,Tb,dTbdpsi,Epsi,Gb,Qb)

!----------------------------------------------------------------------------------------------- 
!For a plasma of NBB species, with charge ZB, mass AB and in regime REGB,
!with kinetic profiles characterized at radial position s by nb,dnbdpsi,Tb,dTbdpsi and Epsi,
!calculates radial flux of particles and energy (and Epsi if ambipolarity is imposed).
!For different values of it, different calculations are made, see below.
!-----------------------------------------------------------------------------------------------
!?it
  
  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER it,NBB,REGB(NBB)
  REAL*8 ZB(NBB),AB(NBB),s,nb(NBB),dnbdpsi(NBB),Tb(NBB),dTbdpsi(NBB)
  !Input/output
  REAL*8 Epsi
  !Output
  REAL*8 Gb(NBB),Qb(NBB),L1b(NBB),L2b(NBB),L3b(NBB)
  !Others
  INTEGER, PARAMETER :: nrootx=5
  INTEGER ib,iEpsi,iroot,nroot
  REAL*8 ephi1oTsize,Epsiacc,Ebx,dEpsi,Jr(NER+1),Jr_old,q
  REAL*8 Epsimin,Epsimax,Epsi1(nrootx),Epsi2(nrootx)
  !Time
  CHARACTER*30, PARAMETER :: routine="SOLVE_DKE_QN_AMB"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)
  !Calculate Epsi by solving ambipolarity of the neoclassical fluxes using bisection
  IF(SOLVE_AMB.OR.SCAN_ER) THEN

     WRITE(iout,*) 'Scanning Er'
     Epsiacc=1E3*ERACC/psip
     !Find changes of sign of the radial neoclasical current (nroot, may be more than 1)
     IF(ERMAX-ERMIN.GT.-1E-3) THEN !between Epsimin and Epsimax
        Epsimin=1E3*ERMIN/psip         
        Epsimax=1E3*ERMAX/psip
     ELSE 
        !1.5 times ion root solution with ions in the 1/nu regime
        Epsimin= 1.5*(3.5*dTbdpsi(2)+Tb(2)*dnbdpsi(2)/nb(2))/ZB(2) 
        IF(Epsimin.GT.-1E+4/psip) Epsimin=-1E+4/psip
        !1.5 times electron root solution with electrons in the 1/nu regime
        Epsimax=-1.5*(3.5*dTbdpsi(1)+Tb(1)*dnbdpsi(1)/nb(1))
        IF(Epsimax.LT.+1E+4/psip) Epsimax=+1E+4/psip

        Epsi=Epsimin
        dEpsi=(Epsimax-Epsimin)/(NER-1)
        DO iEpsi=2,NER
           Epsi=Epsi+dEpsi
           IF(Epsi.GT.0) THEN
              IF(Epsi.LT.dEpsi/2) THEN
                 Epsimin=Epsimin-Epsi
                 Epsimax=Epsimax-Epsi
              ELSE
                 Epsimin=Epsimin+(dEpsi-Epsi)
                 Epsimax=Epsimax+(dEpsi-Epsi)
              END IF
              EXIT
           END IF
        END DO
     END IF
       
     Epsi=Epsimin
     WRITE(iout,'(" Calculating for Er=",1pe13.5,", kV/m")') Epsi*psip/1E3
     CALL CALC_FLUXES(it,NBB,ZB,AB,REGB,s,nb,dnbdpsi,Tb,dTbdpsi,Epsi,&
          & Gb,Qb,L1b,L2b,L3b,ephi1oTsize)
     Jr(1)=SUM(ZB(1:NBB)*nb(1:NBB)*Gb(1:NBB))
     dEpsi=(Epsimax-Epsimin)/(NER-1)
     nroot=0
     DO iEpsi=2,NER
        WRITE(iout,'(" Calculating for Er=",1pe13.5,", kV/m")') Epsi*psip/1E3
        Jr_old=Jr(iEpsi-1)
        Epsi=Epsi+dEpsi
        WRITE(iout,*) 'dEpsi',Epsi,dEpsi,psip,Epsimax,Epsimin
        CALL CALC_FLUXES(it,NBB,ZB,AB,REGB,s,nb,dnbdpsi,Tb,dTbdpsi,Epsi,&
             & Gb,Qb,L1b,L2b,L3b,ephi1oTsize)
        Jr(iEpsi)=SUM(ZB(1:NBB)*nb(1:NBB)*Gb(1:NBB))
        IF(Jr(iEpsi)*Jr_old.LT.0) THEN
           IF((Jr(iEpsi)-Jr_old)/dEpsi.GT.0) THEN !only estable roots
              nroot=nroot+1
              Epsi1(nroot)=Epsi-dEpsi  !root lies in Epsi1<Epsi<Epsi2 or Epsi1>Epsi>Epsi2 
              Epsi2(nroot)=Epsi
           END IF
        END IF
     END DO

     IF(SCAN_ER) THEN 
        CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)
        RETURN
     END IF

     WRITE(iout,*) 'Calculating Er'
     !Find nroot solutions of the ambipolar equation using bisection
     DO iroot=1,nroot        
        CALL CALC_FLUXES(it,NBB,ZB,AB,REGB,s,nb,dnbdpsi,Tb,dTbdpsi,Epsi1(iroot),&
             & Gb,Qb,L1b,L2b,L3b,ephi1oTsize)
        Jr_old=SUM(ZB(1:NBB)*nb(1:NBB)*Gb(1:NBB))
        CALL CALC_FLUXES(it,NBB,ZB,AB,REGB,s,nb,dnbdpsi,Tb,dTbdpsi,Epsi2(iroot),&
             & Gb,Qb,L1b,L2b,L3b,ephi1oTsize)
        Jr(NER+1)=SUM(ZB(1:NBB)*nb(1:NBB)*Gb(1:NBB))
        IF(Jr_old.LT.0.) THEN
           Ebx=Epsi1(iroot)
           dEpsi=Epsi2(iroot)-Epsi1(iroot) 
        ELSE
           Ebx=Epsi2(iroot)
           dEpsi=Epsi1(iroot)-Epsi2(iroot) 
        ENDIF
        DO iEpsi=1,NER
           dEpsi=dEpsi*.5 
           Epsi=Ebx+dEpsi           
           WRITE(iout,'(" Calculating for Er=",1pe13.5,", V/m")') Epsi*psip
           CALL CALC_FLUXES(it,NBB,ZB,AB,REGB,s,nb,dnbdpsi,Tb,dTbdpsi,Epsi,&
                & Gb,Qb,L1b,L2b,L3b,ephi1oTsize)
           Jr(NER+1)=SUM(ZB(1:NBB)*nb(1:NBB)*Gb(1:NBB))
           IF(Jr(NER+1).LE.0.) Ebx=Epsi
           !When the root has been found with precision dEpsi, exit
           IF(ABS(dEpsi).LT.Epsiacc) THEN
              Epsi1(iroot)=Epsi
              Epsi2(iroot)=Epsi
              EXIT
           END IF
        END DO
        
        WRITE(600+myrank,'(1000(1pe13.5))') s,Epsi*psip,&
             & (Gb(ib)/psip,Qb(ib)/psip,L2b(ib)/psip/psip,L3b(ib)/psip/psip,&
             & nb(ib),dnbdpsi(ib)/nb(ib)*psip,Tb(ib),dTbdpsi(ib)/Tb(ib)*psip,&
             & zb(ib),ib=1,MAX(2,NBB)),ephi1oTsize,iota,spol,psip        
        IF(TASK3D) THEN
           WRITE(5600+myrank,'(1000(1pe13.5))') SQRT(s),&
             & nb(1)*1E19,nb(2)*1E19,ZERO,Tb(1),Tb(2),ZERO,&
             & Epsi*psip/1E3,nb(1)*1E19*Gb(1)/psip,&
             & 1.602*nb(1)*Tb(1)*Qb(1)/psip,1.602*nb(2)*Tb(2)*Qb(2)/psip,&
             & Qb(1)/psip/psip/dTbdpsi(1)*Tb(1),Qb(2)/psip/psip/dTbdpsi(2)*Tb(2)
        ELSE IF(TANGO) THEN
           WRITE(5600+myrank,'(1000(1pe13.5))') SQRT(s),Epsi*psip,& 
                & 1.602*nb(2)*Tb(2)*Qb(2)/psip,1.602*nb(1)*Tb(1)*Qb(1)/psip
        END IF
     END DO

     !Find most stable root
     IF(nroot.GT.1) THEN
        IF(ER_ROOT.EQ.-1) THEN
           Epsi=Epsi1(1)
           CALL CALC_FLUXES(it,NBB,ZB,AB,REGB,s,nb,dnbdpsi,Tb,dTbdpsi,Epsi,&
                & Gb,Qb,L1b,L2b,L3b,ephi1oTsize)   
        ELSE IF(ER_ROOT.EQ.1) THEN
           Epsi=Epsi1(2)
        ELSE
           q=0
           dEpsi=(Epsimax-Epsimin)/(NER-1)
           DO iEpsi=1,NER
              Epsi=Epsimin+(iEpsi-1)*dEpsi
              IF(Epsi.LT.Epsi1(1)) CYCLE
              IF(ABS(Epsi).LT.SMALL.AND.(.NOT.TANG_VM)) CYCLE
              IF(Epsi.GT.Epsi2(2)) EXIT
              q=q+Jr(iEpsi)
           END DO
           IF(q.GE.0) THEN
              Epsi=Epsi1(1)
              CALL CALC_FLUXES(it,NBB,ZB,AB,REGB,s,nb,dnbdpsi,Tb,dTbdpsi,Epsi,&
                   & Gb,Qb,L1b,L2b,L3b,ephi1oTsize)   
           ELSE
              Epsi=Epsi1(2)
           END IF
        END IF
     ELSE IF(nroot.EQ.0) THEN
        WRITE(iout,*) 'No root found'
        Gb=0
        Qb=0
        Epsi=0
     END IF

  !Use pre-calculated Epsi
  ELSE

     WRITE(iout,'(" Calculating for Er=",1pe13.5,", kV/m")') Epsi*psip/1E3
     CALL CALC_FLUXES(it,NBB,ZB,AB,REGB,s,nb,dnbdpsi,Tb,dTbdpsi,Epsi,&
          & Gb,Qb,L1b,L2b,L3b,ephi1oTsize)
     WRITE(600+myrank,'(1000(1pe13.5))') s,Epsi*psip,&
          & (Gb(ib)/psip,Qb(ib)/psip,L2b(ib)/psip/psip,L3b(ib)/psip/psip, &
          & nb(ib),dnbdpsi(ib)/nb(ib)*psip,Tb(ib),dTbdpsi(ib)/Tb(ib)*psip,&
          & Zb(ib),ib=1,MIN(2,NBB)),ephi1oTsize,iota,spol,psip,(-rad_R*Epsi/vth(1)/iota),(-rad_R*Epsi/vth(2)/iota),rad_R
!          & PI*vth(1)*vth(1)*vth(1)*m_e*m_e*Ab(1)*Ab(1)/(16.*aiota*rad_R*borbic(0,0)*borbic(0,0)*Zb(1)*Zb(1))
     
  END IF
  
  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE SOLVE_DKE_QN_AMB


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALC_FLUXES(it,NBB,ZB,AB,REGB,s,nb,dnbdpsi,Tb,dTbdpsi,Epsi,Gb,Qb,L1b,L2b,L3b,ephi1oTsize)

!----------------------------------------------------------------------------------------------- 
!For a plasma of NBB species, with charge ZB, mass AB and in regime REGB,
!with kinetic profiles characterized at radial position s by nb,dnbdpsi,Tb,dTbdpsi and Epsi,
!calculate radial flux of particles  Gb and energy Qb and transport coefficients L1b, L2b and L3b
!(and size of varphi1 if quasineutrality is imposed).
!-----------------------------------------------------------------------------------------------
  
  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER it,NBB,REGB(NBB)
  REAL*8 ZB(NBB),AB(NBB),s,nb(NBB),dnbdpsi(NBB),Tb(NBB),dTbdpsi(NBB),Epsi
  !Output
  REAL*8 Gb(NBB),Qb(NBB),L1b(NBB),L2b(NBB),L3b(NBB),ephi1oTsize
  !Others
  INTEGER, SAVE :: nalphab=0
  LOGICAL regp(nbb)
  INTEGER iv,kv,ib,kb,jb(NBB),ib_kin,jt,jt0,nm,Nnmpm1,iz
  INTEGER jv(nv) /1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, & !start from high collisionality
	& 17,18,19,20,21,22,23,24,25,26,27,28/
  REAL*8, ALLOCATABLE :: Grhs(:,:,:),Qrhs(:,:,:),L1rhs(:,:,:),L2rhs(:,:,:),L3rhs(:,:,:)
  REAL*8, ALLOCATABLE :: Gphi(:,:)  ,Qphi(:,:)  ,L1phi(:,:),  L2phi(:,:),  L3phi(:,:)
  REAL*8, ALLOCATABLE :: phi1(:,:),Mbb(:,:),trM(:,:)! n1(:,:,:),
  REAL*8, ALLOCATABLE :: trig(:,:),dtrigdz(:,:),dtrigdt(:,:)
  REAL*8, ALLOCATABLE :: n1nmrhs(:,:,:),Mbbnmrhs(:,:),trMnmrhs(:,:)
  REAL*8 sumzot,x2(nv)
  REAL*8 D11(Nnmp,Nnmp),dL2dv
  REAL*8 zeta(nax),theta(nax)
  REAL*8 dn1nmdv(Nnmp,Nnmp),phi1nm(Nnmp,Nnmp),phi1anm(Nnmp,Nnmp)
  REAL*8 Mbbnm(Nnmp),trMnm(Nnmp)
  REAL*8 phi1c(100,Nnmp)
  REAL*8 Db(NBB),Vb(NBB),n1nmb(Nnmp,NBB)
  REAL*8 D31,L31(NBB),L32(NBB)
  !Time
  CHARACTER*30, PARAMETER :: routine="CALC_FLUXES"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart
#ifdef MPIandPETSc
  !Others
!  INTEGER ierr
  INCLUDE "mpif.h"
#endif
  WRITE(iout,*) 'here'
!!$  IF(USE_SHAING) THEN
!!$     DO ib=1,MIN(2,nbb)
!!$        CALL DKE_CONSTANTS(ib,NBB,ZB,AB,REGB,nb,dnbdpsi,Tb,dTbdpsi,Epsi,.TRUE.)      
!!$        CALL CALC_NTV(jt,s,ib,regb(ib),Zb(ib),Ab(ib),nb(ib),dnbdpsi(ib),Tb(ib),dTbdpsi(ib),Epsi,L1b(ib),Gb(ib))
!!$     END DO
!!$     WRITE(600+myrank,'(1000(1pe13.5))') s,Epsi*psip,&
!!$          & (Gb(ib)/psip,Qb(ib)/psip,L1b(ib)/psip/psip,L2b(ib)/psip/psip,&
!!$          & nb(ib),dnbdpsi(ib)/nb(ib)*psip,Tb(ib),dTbdpsi(ib)/Tb(ib)*psip,&
!!$          & zb(ib),ib=1,MIN(2,NBB)),ephi1oTsize,iota,spol,psip
!!$     RETURN
!!$  END IF

  CALL CPU_TIME(tstart)
  WRITE(iout,*) 'there',NBB,Nnmp,Nbulk

  ALLOCATE(Grhs(NBB,Nnmp,Nnmp),Qrhs(NBB,Nnmp,Nnmp),L1rhs(NBB,Nnmp,Nnmp),L2rhs(NBB,Nnmp,Nnmp),L3rhs(NBB,Nnmp,Nnmp))
  WRITE(iout,*) 'there',NBB,Nnmp,Nbulk
  Grhs=0
  Qrhs=0
  L1rhs=0
  L2rhs=0
  L3rhs=0
  L31=0
  L32=0
  WRITE(iout,*) 'there',NBB,Nnmp,Nbulk
  ALLOCATE(Gphi(NBB,Nnmp),Qphi(NBB,Nnmp),L1phi(NBB,Nnmp),L2phi(NBB,Nnmp),L3phi(NBB,Nnmp))
  WRITE(iout,*) 'there',NBB,Nnmp,Nbulk
  ALLOCATE(n1nmrhs(Nnmp,Nnmp,NBULK),Mbbnmrhs(Nnmp,Nnmp),trMnmrhs(Nnmp,Nnmp))
  WRITE(iout,*) 'ewhere',jt

  Nnmpm1=Nnmp-1
  it=it
  !Calculate prefactor in quasineutrality equation
  sumzot=0
  DO ib=1,2  !valid for trace impurities 
     sumzot=sumzot+ABS(ZB(ib))/Tb(ib)
  END DO

!  jv(1)=iv0   !start with v~=v_th (this can be useful to set precission with thermal particles)
!  jv(iv0)=1 
!  DO iv=1,nv  !start from low collisionalities 
!     jv(iv)=nv-iv+1
!  END DO 
  jb(1)=2     !start with bulk ions
  jb(2)=1     !then, electrons
  DO ib=3,NBB !finally, impurities
     jb(ib)=ib
  END DO
  DO ib=1,NBB !determine if contribution of species ib is considered in quasineutrality
     regp(ib)=((regb(ib).EQ.0).OR.(regb(ib).EQ.2).OR.(regb(ib).EQ.3))
  END DO
  !No QS
  IF(.NOT.SOLVE_QN.AND..NOT.TANG_VM                          ) jt0=01 !valid in large aspect ratio limit
  IF(.NOT.SOLVE_QN.AND.     TANG_VM                          ) jt0=02 !inconsistent !MANUAL
  !Solve QN, necessary close to omnigeneity
  IF(     SOLVE_QN.AND.     TANG_VM.AND.(regp(1) .OR.regp(2))) jt0=12 !one adiabatic species in QN 
  IF(     SOLVE_QN.AND.     TANG_VM.AND.(regp(1).AND.regp(2))) jt0=22 !both kinetic species in QN
  IF(    SOLVE_QN.AND..NOT.TANG_VM.AND.(regp(1) .OR.regp(2))) jt0=11 !inconsistent
  IF(    SOLVE_QN.AND..NOT.TANG_VM.AND.(regp(1).AND.regp(2))) jt0=21 !inconsistent
  WRITE(iout,*) 'ewhere',jt
  IF(REGB(1).LT.0 .AND.REGB(2).LT.0) jt0=00 !use DKES

  IF(regp(1)) ib_kin=1
  IF(regp(2)) ib_kin=2
  IF(regp(1).AND.regp(2)) THEN
     ib_kin=2
     IF((AB(2).GT.AB(1).AND.Tb(1)/Tb(2).GT.3).OR. &
      & (AB(1).GT.AB(2).AND.Tb(2)/Tb(1).LT.3)) ib_kin=1
  END IF

     
  phi1c=0
  DO jt=MOD(jt0,10),jt0
     WRITE(iout,*) 'ewhere',kb

!     IF(jt.EQ.0.AND.(.NOT.DKES_READ.OR..NOT.(REGB(1).EQ.-1.AND.REGB(2).EQ.-1))) CYCLE
!     IF(jt.GT.0.AND.SUM(regb(1:nbb)).EQ.0) CYCLE
     IF((MOD(jt0,10).NE.MOD(jt,10)).AND.jt.NE.0) CYCLE

     IF(jt.EQ.00) WRITE(iout,&
          & '(" Using DKES transport coefficients if available")')
     IF(jt.EQ.01) WRITE(iout,&
          & '(" Without tangential magnetic drift, not solving quasineutrality")')
     IF(jt.EQ.02) WRITE(iout,&
          & '(" With tangential magnetic drift, not solving quasineutrality")') 
!     IF(jt.EQ.11) WRITE(iout,&
!          & '(" Without tangential magnetic drift, solving quasineutrality")')
     IF(jt.EQ.12) WRITE(iout,&
          &'(" With tangential magnetic drift, solving quasineutrality")') 
!     IF(jt.EQ.21) WRITE(iout,&
!          & '(" Without tangential magnetic drift, solving quasineutrality, all kinetic")')
     IF(jt.EQ.22) WRITE(iout,&
          & '(" With tangential magnetic drift, solving quasineutrality, all kinetic")') 
     !Scan in species
     DO kb=1,NBB
        ib=jb(kb)
        IF(ib.EQ.1) WRITE(iout,'(" Electrons")')
        IF(ib.EQ.2) WRITE(iout,'(" Bulk ions, Z=",f3.0,", A=",f6.3)') ZB(ib),AB(ib)
        IF(ib.GT.2) WRITE(iout,'(" Impurities #",I2,", Z=",f4.0,", A=",f8.4)') ib,ZB(ib),AB(ib)
        !Calculate (v,species)-dependent constants
        CALL DKE_CONSTANTS(ib,NBB,ZB,AB,REGB,nb,dnbdpsi,Tb,dTbdpsi,Epsi,.TRUE.)      
        D11=0 
        IF(jt.LT.10.AND.REGB(ib).LT.10) THEN!.AND.ib.LE.NBULK) THEN
           x2=v*v/vth(ib)/vth(ib)
           !Scan in v for the calculation of dQ/dv, dGamma/dv, etc
           DO kv=1,nv
              !Perform monoenergetic calculation
              iv=jv(kv)
              CALL CALC_MONOENERGETIC(ib,ZB(ib),AB(ib),REGB(ib),regp(ib),jt,iv,Epsi,&
                   & phi1c(jt,:),Mbbnm,trMnm,D11,nalphab,zeta,theta,dn1nmdv,D31)
              IF(jt.GT.0.AND.(QN.OR.TRACE_IMP).AND.kb.EQ.1.AND.kv.EQ.1) THEN !.AND.REGB(ib).EQ.3
                 phi1nm=0
                 phi1anm=0
                 n1nmrhs=0
                 IF(TRACE_IMP) THEN
                    Mbbnmrhs=0
                    trMnmrhs=0
                 END IF
              END IF
              !Calculate thermal transport coefficients and radial fluxes
              CALL INTEGRATE_V(jt,ib,Ab(ib),regp(ib),Tb(ib),iv,D11,dn1nmdv,&
                   & L1rhs(ib,:,:),L2rhs(ib,:,:),L3rhs(ib,:,:),Grhs(ib,:,:),Qrhs(ib,:,:),&
                   & n1nmrhs(:,:,ib),Mbbnmrhs,trMnmrhs,D31,L31(ib),L32(ib))

              !Check convergence
              dL2dv=D11(1,1)*x2(iv)*weight(iv)
              IF(DEBUG) WRITE(4400+myrank,'(2I4,30(1pe13.5))') ib,iv,s,Epsi*psip,x2(iv),L2rhs(ib,1,1),dL2dv
              !Check convergence
              IF(regb(ib).GT.-2.AND.iv.GT.iv0.AND.jv(2).GT.jv(1).AND.PREC_INTV.GT.0) THEN
                  IF (dL2dV.LT.PREC_DQDV*L2rhs(ib,1,1)) THEN
                     WRITE(iout,'(" Integral in of species #",I1,"&
                          & , converged for iv=",I2,", v/v_th=",f7.4)') ib,iv,SQRT(x2(iv))
                     EXIT
                  ELSE IF(iv.EQ.nv) THEN
                     WRITE(1100+myrank,*) 'Not converged in v'
                  END IF
               END IF
            END DO
            WRITE(iout,'(" L2/L1 for species #",I1," is ",f4.1,&
                 & " (should be 3 for plateau, 5 for 1/nu, 2.75 for sqrt(nu) and 2.5 for sb-p)")') &
                 & ib,L2rhs(ib,1,1)/L1rhs(ib,1,1)
           !Sums over species
            IF(QN.AND.regp(ib)) THEN
               n1nmb(:,ib)=n1nmrhs(:,1,ib)
               IF(ib.EQ.ib_kin) phi1anm=phi1anm+(Zb(ib)/ABS(Zb(ib)))*n1nmrhs(:,:,ib)/sumzot 
               IF(jt0.GT.20.OR.TRIVIAL_QN) phi1nm =phi1nm +(Zb(ib)/ABS(Zb(ib)))*n1nmrhs(:,:,ib)/sumzot
            END IF

            !Give values for the first jt iteration
            Gb(ib) = Grhs(ib,1,1)
            Qb(ib) = Qrhs(ib,1,1)
            L1b(ib)=L1rhs(ib,1,1)
            L2b(ib)=L2rhs(ib,1,1)
            L3b(ib)=L3rhs(ib,1,1)    
            IF(ANISOTROPY) THEN
               Mbbnm=Mbbnmrhs(:,1)
               trMnm=trMnmrhs(:,1)
            END IF

            !After bulk species have been calculated, solve QN
            IF(QN.AND.kb.EQ.2) CALL CALC_QN(jt,jt0,phi1anm,phi1nm,phi1c)

         ELSE

            IF(regp(ib).AND.ib.LE.NBULK) THEN
               !Calculate fluxes etc using phi1c coefficients from QN
!               CALL WRITE_PHI1(nalphab,Tb(2))
               CALL READ_PHI1(nalphab,phi1c(jt,:))
               CALL CALCULATE_WITH_VARPHI1(ib,regb(ib),phi1c(jt,:),&
                    & Grhs(ib,:,:),Qrhs(ib,:,:),L1rhs(ib,:,:),L2rhs(ib,:,:),L3rhs(ib,:,:),&
                    & Gphi(ib,:)  ,Qphi(ib,:)  ,L1phi(ib,:)  ,L2phi(ib,:),  L3phi(ib,:),&
                    & Gb(ib)      ,Qb(ib)      ,L1b(ib)      ,L2b(ib),      L3b(ib),&
                    & n1nmrhs(:,:,ib),Mbbnmrhs,trMnmrhs,&
                    & n1nmb(:,ib)    ,Mbbnm   ,trMnm)
            ELSE
               Gphi(ib,:)=0
               Qphi(ib,:)=0  
               L1phi(ib,:)=0
               L2phi(ib,:)=0
               L3phi(ib,:)=0
               Gb(ib)=0
               Qb(ib)=0
               L1b(ib)=0
               L2b(ib)=0
               L3b(ib)=0
               n1nmb(:,ib)=0
            END IF

         END IF

         IF(ZERO_PHI1) phi1c=0
         !Plot modes of phi1
         IF(kb.EQ.2.AND.QN.AND.(jt.EQ.jt0.OR.COMPARE_MODELS)) THEN
            DO nm=2,Nnmp
               IF(COMPARE_MODELS) WRITE(4100+myrank,'(2I5,4(1pe13.5))') &
                    & jt,INT((nm-1)/Nnm),ext_np(nm),ext_mp(nm),phi1c(jt,nm)
               IF(jt.EQ.jt0) WRITE(400+myrank,'(1pe13.5,I5,4(1pe13.5))') &
                   & s,INT((nm-1)/Nnm),ext_np(nm),ext_mp(nm),phi1c(jt,nm)
            END DO
         END IF

         IF((jt.EQ.jt0.OR.COMPARE_MODELS).AND.(QN.OR.TRACE_IMP)) THEN
            IF(kb.EQ.2) THEN
               IF(.NOT.ALLOCATED(trig)) THEN    
                  !Once precision in alpha is decided, allocate variables and 
                  !calculate trigonometric quantities
                  IF(nalphab.EQ.0) THEN
                     nalphab=64
                     DO iz=1,nalphab
                        zeta(iz)=(iz-1)*TWOPI/nzperiod/nalphab
                     END DO
                     theta=zeta*nzperiod
                  END IF
                  ALLOCATE(trig(Nnmp,nalphab*nalphab),&
                       & dtrigdz(Nnmp,nalphab*nalphab),dtrigdt(Nnmp,nalphab*nalphab),&
                       & phi1(nalphab,nalphab),Mbb(nalphab,nalphab),trM(nalphab,nalphab))
               END IF               
               !Calculate (zeta,theta) map of varphi1, Mbb and trM
               CALL PRECALC_TRIG(nalphab,zeta(1:nalphab),theta(1:nalphab),trig,dtrigdz,dtrigdt)
               CALL PREPARE_IMP_CALC(jt,jt0,nbb,nalphab,trig,dtrigdz,dtrigdt,&
                    & n1nmb,Mbbnm,trMnm,phi1c(jt,:),Ab(2),Tb(2),Epsi,s,zeta(1:nalphab),theta(1:nalphab),&
                    & phi1,Mbb,trM)
               ephi1oTsize=(MAXVAL(phi1)-MINVAL(phi1))/Tb(2)/2.
               !Check if varphi1 and M exist from previous calculations
               IF(.NOT.SOLVE_AMB) THEN
                  CALL READ_BULKSPECIES(nalphab,"ph1",phi1,ONE)
                  CALL READ_BULKSPECIES(nalphab,"Mbb",Mbb ,ONE)
                  CALL READ_BULKSPECIES(nalphab,"trM",trM ,ONE)
               END IF
               IF(QN.AND.jt.EQ.jt0) THEN
                  CALL WRITE_BULKSPECIES(s,nalphab,phi1,"ph1")
                  CALL WRITE_BULKSPECIES(s,nalphab,Mbb ,"Mbb")
                  CALL WRITE_BULKSPECIES(s,nalphab,trM ,"trM")
               END IF

            !Calculate trace impurities
            ELSE IF(REGB(ib).GT.3.AND.TRACE_IMP) THEN
               CALL TRACE_IMPURITIES(jt,ib,NBB,ZB,AB,REGB,s,nb,dnbdpsi,Tb,dTbdpsi,Epsi,&
                       & phi1c(jt,:),Mbbnm,trMnm,&
                       & nalphab,zeta(1:nalphab),theta(1:nalphab),phi1,Mbb,trM,&
                       & Gb(ib),Db(ib),Vb(ib))
            END IF
         END IF

      END DO

      !Plots results
      CALL PLOT_FLUX(jt,jt0,nbb,s,Epsi,Gb,Qb,L1b,L2b,L3b,Zb,nb,dnbdpsi,Tb,dTbdpsi,ephi1oTsize,&
           & Grhs(:,1,1),Qrhs(:,1,1),L1rhs(:,1,1),L2rhs(:,1,1),L3rhs(:,1,1),Gphi,Qphi,L1phi,L2phi,L3phi)
      
   END DO
   DEALLOCATE(Grhs,Qrhs,L1rhs,L2rhs,L3rhs,Gphi,Qphi,L1phi,L2phi,L3phi)
   DEALLOCATE(n1nmrhs,Mbbnmrhs,trMnmrhs)
   IF(QN.OR.TRACE_IMP) DEALLOCATE(trig,dtrigdz,dtrigdt,phi1,Mbb,trM)

   CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)
   
 END SUBROUTINE CALC_FLUXES


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALC_MONOENERGETIC(ib,Zb,Ab,regb,regp,jt,iv,Epsi,phi1c,Mbbnm,trMnm,&
     & D11,nalphab,zeta,theta,dn1nmdv,D31)
 
!----------------------------------------------------------------------------------------------- 
!Calculate monoenergetic transport coefficient D11 and contribution quasineutrality dn1nmdv
!at nalphabxnalphab grid in (zeta,theta) for collisionality cmul=nu(iv)/v(iv) and 
!normalized radial electric field Epsi/v(iv) and in the presence of phi1c.
!Species quantities ib, Zb, Ab and regb are used to determine which monoenergetic DKE is solved
!-----------------------------------------------------------------------------------------------
  
  USE GLOBAL
  IMPLICIT NONE
  !Input
  LOGICAL regp
  INTEGER ib,regb,jt,iv
  REAL*8 Zb,Ab,Epsi,phi1c(Nnmp),Mbbnm(Nnmp),trMnm(Nnmp)
  !Output
  INTEGER nalphab
  REAL*8 D11(Nnmp,Nnmp),zeta(nax),theta(nax),dn1nmdv(Nnmp,Nnmp),D31
  !Others
  INTEGER ia,il
  REAL*8 efield,cmul,nustar,rhostar,dn1dv(nax,nax),dD11(1,1)
  !Time
  CHARACTER*30, PARAMETER :: routine="CALC_MONOENERGETIC"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)
  
  dn1dv=0
  dn1nmdv=0
  !Parameters of the monoenergetic drift-kinetic equation
  efield =Epsi*psip/v(iv)
  cmul   =nu(iv)/v(iv)/2.
  nustar =rad_R*nu(iv)/v(iv)/aiota
  rhostar=v(iv)*Ab*m_e/Zb/borbic(0,0)/rad_R
  IF(iv.EQ.iv0.OR.DEBUG) THEN
     WRITE(iout,'(" NU  ",1pe13.5)') nu(iv)
     IF(CALC_DB) THEN
        WRITE(iout,'(" EFIELD",1pe13.5)') efield
        WRITE(iout,'(" CMUL  ",1pe13.5)') cmul
     ELSE
        WRITE(iout,'(" EPS^3/2        ",1pe13.5)') eps32
        WRITE(iout,'(" RHOSTAR/EPS    ",1pe13.5)') rhostar/eps
        WRITE(iout,'(" RHOSTAR/EPS^1/2",1pe13.5)') rhostar/SQRT(eps)
        WRITE(iout,'(" RHOSTAR        ",1pe13.5)') rhostar
        WRITE(iout,'(" NUSTAR         ",1pe13.5)') nustar
        WRITE(iout,'(" NUSTAR(1/nu)   ",1pe13.5)') rad_R*cmul_1nu/aiota
                
     END IF
  END IF
  !Determine regime in which species ib is:
  IF(regb.LE.0.AND.regb.GT.-2) THEN  !depends on collisionality
     IF(cmul.GT.cmul_1NU) THEN
        IF(DKES_READ) THEN
           CALL INTERP_DATABASE(jt-1,iv,Epsi,D11(1,1),D31,.FALSE.)
        ELSE
           IF(cmul.GT.cmul_PS) THEN
              CALL CALC_PS(iv,Epsi,D11(1,1))    
           ELSE IF(cmul.GT.cmul_1NU) THEN   
!              CALL CALC_PLATEAU_OLD(iv,Epsi,D11(1,1))
              CALL CALC_PLATEAU(iv,Epsi,D11,dn1nmdv)
           END IF
        END IF
     ELSE
        IF(regb.EQ.-1.AND.DKES_READ) THEN
           CALL INTERP_DATABASE(jt-1,iv,Epsi,D11(1,1),D31,.FALSE.)
           !Correction of DKES with the effect of the tangential magnetic drift
           IF(TANG_VM) THEN
              CALL CALC_LOW_COLLISIONALITY(iv,Epsi,phi1c,Mbbnm,trMnm,&
                   & dD11,nalphab,zeta,theta,dn1dv,dn1nmdv)
!              D11(1,1)=D11(1,1)*dD11(1,1)
              D11(1,1)=D11(1,1)+dD11(1,1)
              TANG_VM=.FALSE.  
              INC_EXB=.TRUE.
              CALL CALC_LOW_COLLISIONALITY(iv,Epsi,phi1c,Mbbnm,trMnm,&
                   & dD11,nalphab,zeta,theta,dn1dv,dn1nmdv)
!              D11(1,1)=D11(1,1)/dD11(1,1)
              D11(1,1)=D11(1,1)-dD11(1,1)
              TANG_VM=.TRUE.
              INC_EXB=.FALSE.
           END IF
        ELSE IF(regb.EQ.0.OR..NOT.DKES_READ) THEN
           IF(FAST_AMB) THEN
              CALL INTERP_DATABASE(jt-1,iv,Epsi,D11(1,1),D31,.TRUE.) 
           ELSE
              CALL CALC_LOW_COLLISIONALITY(iv,Epsi,phi1c,Mbbnm,trMnm,&
                   & D11,nalphab,zeta,theta,dn1dv,dn1nmdv)
           END IF
        END IF
     END IF !regime set from input
  ELSE IF(regb.EQ.1) THEN
     CALL CALC_PS(iv,Epsi,D11(1,1))
  ELSE IF(regb.EQ.2) THEN
!     CALL CALC_PLATEAU_OLD(iv,Epsi,D11(1,1))
     CALL CALC_PLATEAU(iv,Epsi,D11,dn1nmdv)
  ELSE IF(MOD(regb,10).EQ.3) THEN
     CALL CALC_LOW_COLLISIONALITY(iv,Epsi,phi1c,Mbbnm,trMnm,&
          & D11,nalphab,zeta,theta,dn1dv,dn1nmdv)
  END IF

  IF(DEBUG.AND.(iv.EQ.iv0.OR.TRIVIAL_QN).AND.regp) THEN
     DO ia=1,nalphab  !contribution of thermal species to n1
        DO il=1,nalphab
           WRITE(4000+myrank,'(3I4,3(1pe13.5))') jt,ib,iv,&
                & zeta(il),theta(ia),dn1dv(ia,il)*weight(iv)*Sdke(iv)
        END DO
     END DO
  END IF

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)
  
END SUBROUTINE CALC_MONOENERGETIC

                     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE PRECALC_TRIG(nalphab,zeta,theta,trig,dtrigdz,dtrigdt)                    

!----------------------------------------------------------------------------------------------- 
!Precalculate trigonometric functions trig, dtrigdz and dtrigdt
!at nalphabxalphab (zeta,theta) array
!-----------------------------------------------------------------------------------------------
  
  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER nalphab
  REAL*8 zeta(nalphab),theta(nalphab)
  !Output
  REAL*8 trig(Nnmp,nalphab*nalphab),dtrigdz(Nnmp,nalphab*nalphab),dtrigdt(Nnmp,nalphab*nalphab)
  !Others
  INTEGER nm,ia,il
  REAL*8 arg,cosine,sine
  !Time
  CHARACTER*30, PARAMETER :: routine="PRECALC_TRIG"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  DO nm=1,Nnm
     DO ia=1,nalphab
        DO il=1,nalphab
           arg=mp(nm)*theta(ia)+np(nm)*nzperiod*zeta(il)
           cosine=COS(arg)
           sine  =SIN(arg)
           trig(    nm,(ia-1)*nalphab+il)=cosine
           trig(Nnm+nm,(ia-1)*nalphab+il)=sine
           IF(QN) THEN
              dtrigdz(    nm,(ia-1)*nalphab+il)= -sine*np(nm)*nzperiod
              dtrigdz(Nnm+nm,(ia-1)*nalphab+il)=cosine*np(nm)*nzperiod
              dtrigdt(    nm,(ia-1)*nalphab+il)= -sine*mp(nm)
              dtrigdt(Nnm+nm,(ia-1)*nalphab+il)=cosine*mp(nm)
           END IF
        END DO
     END DO
  END DO

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE PRECALC_TRIG


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE INTEGRATE_V(jt,ib,Ab,regp,Tb,iv,D11,dn1nmdv,L1rhs,L2rhs,L3rhs,Grhs,Qrhs,&
     & n1nmrhs,Mbbnmrhs,trMnmrhs,D31,L31,L32)

!----------------------------------------------------------------------------------------------- 
!At instant jt and for species ib of mass Ab and temperature Tb in regime regp
!add contribution from monoenergetic iv to L1rhs, L2rhs, L3rhs, Grhs, Qrhs, n1nmrhs, 
!Mbbnmrhs, and trMnmrhs
!-----------------------------------------------------------------------------------------------
  
  USE GLOBAL
  IMPLICIT NONE
  !Input
  LOGICAL regp
  INTEGER jt,ib,iv
  REAL*8 Ab,Tb,D11(Nnmp,Nnmp),dn1nmdv(Nnmp,Nnmp),D31
  !Output
  REAL*8 L1rhs(Nnmp,Nnmp),L2rhs(Nnmp,Nnmp),L3rhs(Nnmp,Nnmp),Grhs(Nnmp,Nnmp),Qrhs(Nnmp,Nnmp)
  REAL*8 n1nmrhs(Nnmp,Nnmp)
  REAL*8 Mbbnmrhs(Nnmp,Nnmp),trMnmrhs(Nnmp,Nnmp)
  REAL*8 L31,L32
  !Others
  REAL*8 v2,x2,wD11(Nnmp,Nnmp),wn1nmS(Nnmp,Nnmp)
  !Time
  CHARACTER*30, PARAMETER :: routine="INTEGRATE_V"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  v2=v(iv)*v(iv)
  x2=v2/vth(ib)/vth(ib)

  wD11=weight(iv)*D11

  L1rhs=L1rhs+wD11
  L2rhs=L2rhs+wD11*x2
  L3rhs=L3rhs+wD11*x2*x2
  Grhs = Grhs-wD11*Sdke(iv)
  Qrhs = Qrhs-wD11*Sdke(iv)*x2

  wn1nmS=weight(iv)*(dn1nmdv/2.)*Sdke(iv)
  IF((QN.OR.TRACE_IMP).AND.regp.AND.jt.GT.0) THEN
     n1nmrhs=n1nmrhs+wn1nmS
     IF(TRACE_IMP.AND.ib.EQ.2) THEN
        Mbbnmrhs=Mbbnmrhs+  (wn1nmS/v(ib))*(-Ab*m_e/Tb+2/v2)
        trMnmrhs=trMnmrhs-2*(wn1nmS/v(iv))*(+Ab*m_e/Tb)
     END IF    
  END IF

  L31=L31+weight(iv)*D31
  L32=L32+weight(iv)*D31*x2
  
  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE INTEGRATE_V


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



SUBROUTINE CALC_QN(jt,jt0,phi1anm,phi1nm,phi1c)

!----------------------------------------------------------------------------------------------- 
!Using matrixes ph1anm and phi1nm calculated at instant jt, fill phi1c at times jt and jt0 
!by solving the quasineutrality equation 
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER jt,jt0
  REAL*8 phi1anm(Nnmp,Nnmp),phi1nm(Nnmp,Nnmp)
  !Input/output
  REAL*8 phi1c(100,Nnmp)
  !Others
  INTEGER ikin,jtt,irow,Nnmpm1
  REAL*8 mat(Nnmp-1,Nnmp-1),rhs(Nnmp-1)
  !Time
  CHARACTER*30, PARAMETER :: routine="CALC_QN"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  Nnmpm1=Nnmp-1
  DO ikin=1,2
     jtt=100
     mat=0
     DO irow=1,Nnmp
        mat(irow,irow)=1
     END DO
     IF(ikin.EQ.1) THEN
        rhs(1:Nnmpm1)=phi1anm(2:Nnmp,1)
        DO irow=1,Nnmpm1
           mat(1:Nnmpm1,irow)=mat(1:Nnmpm1,irow)-phi1anm(2:Nnmp,irow+1)
        END DO
        jtt=jt+10
     ELSE IF(ikin.EQ.2.AND.(jt0.GT.20.OR.TRIVIAL_QN)) THEN
        rhs(1:Nnmpm1)=phi1nm(2:Nnmp,1)
        DO irow=1,Nnmpm1
           mat(1:Nnmpm1,irow)=mat(1:Nnmpm1,irow)-phi1nm(2:Nnmp,irow+1)
        END DO
        jtt=jt+20
     ELSE
        CYCLE
     END IF
     phi1c(jt,2:Nnmp)=rhs
     IF(SOLVE_QN) THEN
        CALL INVERT_QN(Nnmpm1,mat,rhs)           
        phi1c(jtt,2:Nnmp)=rhs
     END IF
  END DO

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)
  
END SUBROUTINE CALC_QN


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               
SUBROUTINE INVERT_QN(nrow,mat,rhs)

!----------------------------------------------------------------------------------------------- 
!Solve mat*x=rhs, and write solution x in rhs
!mat is a nrowxnrow matrix, and rhs has 1 column
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER nrow
  REAL*8 mat(nrow,nrow)
  !Input/output
  REAL*8 rhs(nrow)
  !Others
  LOGICAL, PARAMETER :: USE_SDV=.FALSE.
  CHARACTER*100 serr
  INTEGER ierr,rank,lwork,ipivot(nrow)
  INTEGER, ALLOCATABLE :: iwork(:)
  REAL*8 s_sdv(nrow),rcond
  REAL*8, ALLOCATABLE :: rwork(:)
  REAL*8, ALLOCATABLE :: work(:)

  IF(USE_SDV) THEN
     !Solve equation using SVD decomposition
     lwork=-1
     rcond=-1
     ALLOCATE(rwork(100000),iwork(100000),work(100000))
     
     CALL DGELSD(nrow,nrow,1,mat,nrow,rhs,nrow,&
          & s_sdv,rcond,rank,work,lwork,rwork,iwork,ierr)
     lwork=MIN(100000,INT(work(1)))
     CALL DGELSD(nrow,nrow,1,mat,nrow,rhs,nrow,&
          & s_sdv,rcond,rank,work,lwork,rwork,iwork,ierr)
  ELSE
     !Solves equation using LU decomposition        
     CALL DGESV(nrow,1,mat,nrow,ipivot,rhs,nrow,ierr)     
  END IF
  
  IF(ierr.NE.0) THEN
     serr="Error solving quasineutrality"
     CALL END_ALL(serr,.FALSE.)
  END IF

  IF (ALLOCATED(rwork) ) DEALLOCATE(rwork,iwork,work)
  
END SUBROUTINE INVERT_QN


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALCULATE_WITH_VARPHI1(ib,regb,phi1c,Grhs,Qrhs,L1rhs,L2rhs,L3rhs,&
     & Gphi,Qphi,L1phi,L2phi,L3phi,Gb,Qb,L1b,L2b,L3b,n1nmrhs,Mbbnmrhs,trMnmrhs,&
     & n1nmb,Mbbnm,trMnm)

!----------------------------------------------------------------------------------------------- 
!For species ib in regime regb, for a particular phi1c
!and matrixes Grhs, Qrhs, L1rhs, L2rhs, L3rhs, n1nmrhs, Mbbnmrhs and trMnmrhs
!calculate contributions of each varphi1, Gphi, Qphi, L1phi, L2phi, L3phi and total values 
!Gb, Qb, L1b, L2b, L3b, n1nmb, Mbbnm and trMnm
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER ib,regb
  REAL*8 phi1c(Nnmp),Grhs(Nnmp,Nnmp),Qrhs(Nnmp,Nnmp),L1rhs(Nnmp,Nnmp),L2rhs(Nnmp,Nnmp),L3rhs(Nnmp,Nnmp)
  REAL*8 n1nmrhs(Nnmp,Nnmp),Mbbnmrhs(Nnmp,Nnmp),trMnmrhs(Nnmp,Nnmp)
  !Output
  REAL*8 Gphi(Nnmp),Qphi(Nnmp),L1phi(Nnmp),L2phi(Nnmp),L3phi(Nnmp)
  REAL*8 Gb,Qb,L1b,L2b,L3b,n1nmb(Nnmp),Mbbnm(Nnmp),trMnm(Nnmp)
  !Others
  INTEGER nm,nm2
  !Time
  CHARACTER*30, PARAMETER :: routine="CALCULATE_WITH_VARPHI1"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  phi1c(1)=1
  Gb=0
  Qb=0
  L1b=0
  L2b=0
  L3b=0
  Gphi=0
  Qphi=0
  L1phi=0
  L2phi=0
  L3phi=0
  DO nm=1,Nnmp
     DO nm2=1,Nnmp
        Gphi(nm) =Gphi(nm)  +Grhs(nm,nm2)*phi1c(nm2)
        Qphi(nm) =Qphi(nm)  +Qrhs(nm,nm2)*phi1c(nm2)
        L1phi(nm)=L1phi(nm)+L1rhs(nm,nm2)*phi1c(nm2)
        L2phi(nm)=L2phi(nm)+L2rhs(nm,nm2)*phi1c(nm2)
        L3phi(nm)=L3phi(nm)+L3rhs(nm,nm2)*phi1c(nm2)
     END DO
     Gb =Gb  +phi1c(nm) *Gphi(nm)
     Qb =Qb  +phi1c(nm) *Qphi(nm)
     L1b=L1b+phi1c(nm)*L1phi(nm)
     L2b=L2b+phi1c(nm)*L2phi(nm)
     L3b=L3b+phi1c(nm)*L3phi(nm)
  END DO
  DO nm=1,Nnmp
     Gphi(nm) = Grhs(nm,nm)*phi1c(nm)*phi1c(nm2)
     Qphi(nm) = Qrhs(nm,nm)*phi1c(nm)*phi1c(nm2)
     L1phi(nm)=L1rhs(nm,nm)*phi1c(nm)*phi1c(nm2)
     L2phi(nm)=L2rhs(nm,nm)*phi1c(nm)*phi1c(nm2)
     L3phi(nm)=L3rhs(nm,nm)*phi1c(nm)*phi1c(nm2)
     DO nm2=1,Nnmp
        IF(nm2.EQ.nm) CYCLE
        Gphi(nm) =Gphi(nm) +(Grhs(nm,nm2)+Grhs(nm2,nm))*phi1c(nm)*phi1c(nm2)
        Qphi(nm) =Qphi(nm) +(Grhs(nm,nm2)+Grhs(nm2,nm))*phi1c(nm)*phi1c(nm2)
        L1phi(nm)=L1phi(nm)+(Grhs(nm,nm2)+Grhs(nm2,nm))*phi1c(nm)*phi1c(nm2)
        L2phi(nm)=L2phi(nm)+(Grhs(nm,nm2)+Grhs(nm2,nm))*phi1c(nm)*phi1c(nm2)
        L3phi(nm)=L3phi(nm)+(Grhs(nm,nm2)+Grhs(nm2,nm))*phi1c(nm)*phi1c(nm2)
     END DO
  END DO
  
  IF(regb.EQ.3) n1nmb=n1nmrhs(:,1)
  IF(TRACE_IMP.AND.ib.EQ.2) THEN
     Mbbnm=Mbbnmrhs(:,1)
     trMnm=trMnmrhs(:,1)
  END IF
  DO nm=2,Nnmp
     IF(regb.EQ.3) n1nmb=n1nmb+n1nmrhs(:,nm)*phi1c(nm)   
     IF(QN.AND.TRACE_IMP.AND.ib.EQ.2) THEN
        Mbbnm=Mbbnm+Mbbnmrhs(:,nm)*phi1c(nm)         
        trMnm=trMnm+trMnmrhs(:,nm)*phi1c(nm)         
     END IF
  END DO
  

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)  

END SUBROUTINE CALCULATE_WITH_VARPHI1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE PREPARE_IMP_CALC(jt,jt0,nbb,nalphab,trig,dtrigdz,dtrigdt,n1nmb,Mbbnm,trMnm,&
     & phi1c,Ai,Ti,Epsi,s,zeta,theta,phi1,Mbb,trM)

!----------------------------------------------------------------------------------------------- 
!For instant jt and nbb species, using precalculated quantities trig, dtrigdz and dtrigdt, use
!matrixes n1nmb, Mbbnm and trMnm and actual phi1c to calculate nalphabxnalphab
!(Zeta,theta) map of phi1, Mbb and trM at surface s. Ti is used for normalization
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER jt,jt0,nbb,nalphab
  REAL*8 trig(Nnmp,nalphab*nalphab),dtrigdz(Nnmp,nalphab*nalphab),dtrigdt(Nnmp,nalphab*nalphab)
  REAL*8 n1nmb(Nnmp,nbb),Mbbnm(Nnmp),trMnm(Nnmp)
  REAL*8 phi1c(Nnmp),Ai,Ti,Epsi,s,zeta(nalphab),theta(nalphab)
  !Output
  REAL*8 phi1(nalphab,nalphab),Mbb(nalphab,nalphab),trM(nalphab,nalphab)
  !Others
!  INTEGER, PARAMETER :: npow=10
  INTEGER, PARAMETER :: npow=4
  INTEGER, PARAMETER :: npoints=10
  INTEGER ia,il,nm,index
  REAL*8 n1(nalphab,nalphab,nbb),dphi1dz,dphi1dt,ver,pfact,fact
#ifdef MPIandPETSc
  INTEGER iarray,iparray
  REAL*8 dz,dt,zeta2(nalphab),theta2(nalphab),phi1DR,absnablarDR,dvarphi1dsDR,xDR,yDR,zDR
#endif
! Plateau test
!  REAL*8 phi1plateau(nalphab,nalphab),factphi
  !Time
  CHARACTER*30, PARAMETER :: routine="PREPARE_IMP_CALC"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart
#ifdef MPIandPETSc
  INTEGER ierr,ipow,is0,is1,is,js,ns,ms,lwork,rank
  INTEGER, ALLOCATABLE :: iwork(:)
  REAL*8 varphi1fit,dvarphi1dsfit(nalphab,nalphab),phi1coeff(numprocs),vars(numprocs),varphi1(numprocs)
  REAL*8 s_svd(npow),rcond,mats(numprocs,npow)
  REAL*8, ALLOCATABLE :: rwork(:),work(:)
  INCLUDE "mpif.h"
#endif
  CALL CPU_TIME(tstart)

  Epsi=Epsi
  n1=0
  phi1=0
  Mbb=0
  trM=0

! This is a test for the plateau. Will work for n'=T'=0, B=B_{00}+borbic(-1,1)*cos(\theta-nzperiod*\zeta)
!  factphi=(PI*PI)*(Ti/2)*SQRT(Ti/m_e/PI/PI/PI/2)*(-Epsi/Ti)*(m_e/borbic(0,0)/borbic(0,0))*(-Bzeta*borbic(-1,1))/(abs(nzperiod-iota))

  DO ia=1,nalphab
     DO il=1,nalphab
        dphi1dz=0
        dphi1dt=0
!        phi1plateau(il,ia)=factphi*sin(theta(ia)-nzperiod*zeta(il))
        DO nm=1,Nnmp
!           IF(ABS(mp(MOD(nm-1,Nnm)+1)).LT.0.1.AND.ABS(np(MOD(nm-1,Nnm)+1)).LT.0.1) CYCLE
           IF(ABS(ext_mp(nm)).LT.0.1.AND.ABS(ext_np(nm)).LT.0.1) CYCLE
           index=(ia-1)*nalphab+il
           IF(QN) THEN
              n1(il,ia,:)=n1(il,ia,:)+n1nmb(nm,:)*trig(nm,index)
              phi1(il,ia)=phi1(il,ia)+phi1c(nm)*trig(nm,index)
              dphi1dz=dphi1dz+phi1c(nm)*dtrigdz(nm,index)
              dphi1dt=dphi1dt+phi1c(nm)*dtrigdt(nm,index)
           END IF
           IF(TRACE_IMP) THEN 
              Mbb(il,ia)=Mbb(il,ia)+Mbbnm(nm)*trig(nm,index)
              trM(il,ia)=trM(il,ia)+trMnm(nm)*trig(nm,index)
           END IF
        END DO
        ver=(Btheta*dphi1dz-Bzeta*dphi1dt)/aiBtpBz
        IF(COMPARE_MODELS) WRITE(4200+myrank,'(I4,7(1pe13.5))') jt,zeta(il),theta(ia),&
             phi1(il,ia),phi1(il,ia)/Ti,ver/Ti,n1(il,ia,1),n1(il,ia,2)!,phi1plateau(il,ia)
        IF(jt.EQ.jt0) WRITE(500+myrank,'(8(1pe13.5))') &
             & s,zeta(il),theta(ia),phi1(il,ia),phi1(il,ia)/Ti,ver/Ti,&
             & n1(il,ia,1),n1(il,ia,2)!,phi1plateau(il,ia)
     END DO
  END DO
  
  pfact=SQRT(Ti/(Ai*m_e)) 
  fact=-1.5*(SQPI/SQ2)*pfact*pfact*pfact*pfact*pfact
  Mbb=Mbb*fact
  trM=trM*fact
  
#ifdef MPIandPETSc
  ns=numprocs
  is=myrank+1
  vars=0
  vars(is)=s
  IF(DR_READ) THEN
     
     dz=dzeta0_DR/nalphab
     dt=dtheta0_DR/nalphab
     DO il=1,nalphab
        zeta2(il)=zeta0_DR+(il-nalphab/2.)*dz
     END DO
     DO ia=1,nalphab
        theta2(ia)=TWOPI-(theta0_DR+(ia-nalphab/2.)*dt)
     END DO
     CALL REAL_ALLREDUCE(vars,ns)
     DO ia=1,nalphab
        DO il=1,nalphab
           !Calculate derivative
!           IF(array(il,ia).EQ.0) CYCLE
           varphi1=0
           varphi1(is)=phi1(il,ia)
!           CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
           CALL REAL_ALLREDUCE(varphi1,ns)
!           CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

           IF(npow.LT.5) THEN
              is0=is-npoints
              is1=is+npoints
              IF(is0.LT.1) THEN
                 is1=is1+(1-is0)
                 is0=1
              ELSE IF(is1.GT.ns) THEN
                 is0=is0+(ns-is1)
                 is1=ns
              END IF
              ms=2*npoints+1
           ELSE
              is0=1
              is1=ns
              ms=ns
           END IF
           mats(1:ms,1)=1
           DO ipow=2,npow
              mats(1:ms,ipow)=mats(1:ms,ipow-1)*vars(is0:is1)
           END DO
           phi1coeff(1:ms)=varphi1(is0:is1)
           IF(.NOT.ALLOCATED(rwork)) THEN
              lwork=-1
              ierr=0
              lwork=-1
              rcond=-1
              ALLOCATE(rwork(100000),iwork(100000),work(100000))  
!              CALL DGELSD(ns,npow,1,mats,ns,phi1coeff,ns,s_svd,rcond,rank,work,lwork,rwork,iwork,ierr)  
              CALL DGELSD(ms,npow,1,mats(1:ms,:),ms,phi1coeff(1:ms),ms,s_svd,rcond,rank,work,lwork,rwork,iwork,ierr)  
              lwork=MIN(100000,INT(work(1)))
              mats(1:ms,1)=1
              DO ipow=2,npow
                 mats(1:ms,ipow)=mats(1:ms,ipow-1)*vars(is0:is1)
              END DO
              phi1coeff(1:ms)=varphi1(is0:is1)
           END IF
!           CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!           CALL DGELSD(ns,npow,1,mats,ns,phi1coeff,ns,s_svd,rcond,rank,work,lwork,rwork,iwork,ierr)  
           CALL DGELSD(ms,npow,1,mats(1:ms,:),ms,phi1coeff(1:ms),ms,s_svd,rcond,rank,work,lwork,rwork,iwork,ierr)  
           DO js=1,ms
              varphi1fit=0 
              dvarphi1dsfit(il,ia)=0
              mats(js,1)=1
              DO ipow=1,npow
                 IF(ipow.GT.1) THEN
                    mats(js,ipow)=mats(js,ipow-1)*vars(is0+js-1)
                    dvarphi1dsfit(il,ia)=dvarphi1dsfit(il,ia)+phi1coeff(ipow)*mats(js,ipow-1)*(ipow-1)
                 END IF                 
                 varphi1fit=varphi1fit+phi1coeff(ipow)*mats(js,ipow)
              END DO
!              WRITE(iout,'(I4,8(1pe13.5),I4)') jt,vars(is0+js-1),varphi1(is0+js-1),varphi1fit
           END DO
           mats(1,1)=1
           varphi1fit=0
           dvarphi1dsfit(il,ia)=0
           DO ipow=1,npow
              IF(ipow.GT.1) THEN
                 mats(1,ipow)=mats(1,ipow-1)*s
                 dvarphi1dsfit(il,ia)=dvarphi1dsfit(il,ia)+phi1coeff(ipow)*mats(1,ipow-1)*(ipow-1)
              END IF
              varphi1fit=varphi1fit+phi1coeff(ipow)*mats(1,ipow)
           END DO

           IF(COMPARE_MODELS) WRITE(4600+myrank,'(I4,8(1pe13.5),3(1pe13.5))') jt,zeta(il),theta(ia),phi1(il,ia),&
                & varphi1fit,Ti,dvarphi1dsfit(il,ia),Epsi*psip,Epsi*psip*absnablar(il,ia),&
                & (Epsi*absnablar(il,ia)-dvarphi1dsfit(il,ia)/atorflux)*psip,&!array(il,ia)
                & posx(il,ia),posy(il,ia),posz(il,ia)
           IF(jt.EQ.jt0) WRITE(4500+myrank,'(10(1pe13.5),3(1pe13.5))') s,zeta(il),theta(ia),phi1(il,ia),&
                & varphi1fit,Ti,dvarphi1dsfit(il,ia),Epsi*psip,Epsi*psip*absnablar(il,ia),&
                & (Epsi-dvarphi1dsfit(il,ia)/atorflux)*absnablar(il,ia)*psip,&!array(il,ia),&
                & posx(il,ia),posy(il,ia),posz(il,ia)
        END DO
     END DO
     DO iarray=1,narray
        DO iparray=1,nparray
           IF(zetaDR(iarray,iparray).LT.-0.5) CYCLE
           CALL BILAGRANGE(zeta, theta,         phi1,nalphab,nalphab,zetaDR(iarray,iparray),thetaDR(iarray,iparray),      phi1DR,2)
           CALL BILAGRANGE(zeta, theta,dvarphi1dsfit,nalphab,nalphab,zetaDR(iarray,iparray),thetaDR(iarray,iparray),dvarphi1dsDR,2)
           CALL BILAGRANGE(zeta2,theta2,       zoomx,nalphab,nalphab,zetaDR(iarray,iparray),thetaDR(iarray,iparray),         xDR,2)
           CALL BILAGRANGE(zeta2,theta2,       zoomy,nalphab,nalphab,zetaDR(iarray,iparray),thetaDR(iarray,iparray),         yDR,2)
           CALL BILAGRANGE(zeta2,theta2,       zoomz,nalphab,nalphab,zetaDR(iarray,iparray),thetaDR(iarray,iparray),         zDR,2)
           CALL BILAGRANGE(zeta2,theta2,      zoomdr,nalphab,nalphab,zetaDR(iarray,iparray),thetaDR(iarray,iparray), absnablarDR,2)
           IF(jt.EQ.jt0) WRITE(4500+myrank,'(10(1pe13.5),3(1pe13.5))') &  
                & s,zetaDR(iarray,iparray),thetaDR(iarray,iparray),phi1DR,&
                & iarray*Ti/Ti,Ti,dvarphi1dsDR,Epsi*psip,Epsi*psip*absnablarDR,&
                & (Epsi-dvarphi1dsDR/atorflux)*absnablarDR*psip,&
                & xDR,yDR,zDR
        END DO
     END DO
  END IF
#endif
  
#ifdef MPIandPETSc
  IF(ALLOCATED(rwork)) THEN
     DEALLOCATE(rwork,iwork,work)
  END IF
#endif

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE PREPARE_IMP_CALC


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE READ_PHI1(nalphab,phi1c)

!----------------------------------------------------------------------------------------------- 
!Read varphi1(nakphab,nalphab) and write Fourier coefficients in phi1c
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER nalphab
  !Output
  REAL*8 phi1c(Nnmp)
  !Others
  CHARACTER*12 filename
  INTEGER il,ia,iostat
  REAL*8 dummy,fact,phi1(nalphab,nalphab)
  COMPLEX*16 phi1c_nm(nalphab,nalphab)
  
  !Read original file
  filename="ph1_2D.in"
  OPEN(unit=999+myrank,file=filename,action='read',iostat=iostat)
  IF(iostat.NE.0) RETURN
  WRITE(iout,*) 'File ',filename,' read'
  READ(999+myrank,*) fact
  DO il=1,nalphab
     DO ia=1,nalphab
        READ(999+myrank,*) dummy,dummy,phi1(il,ia)
     END DO
  END DO
  CLOSE(999+myrank)
  phi1=phi1*fact
  CALL FFTF_KN(nalphab,phi1,phi1c_nm)
  CALL FILL_ORBI(nalphab,nalphab,phi1c_nm,Nnmp,phi1c)

END SUBROUTINE READ_PHI1

                
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE WRITE_PHI1(nalphab,Ti)

!----------------------------------------------------------------------------------------------- 
!Read varphi1(nakphab,nalphab) and write Fourier coefficients in phi1c
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER nalphab
  REAL*8 Ti
  !Others
  CHARACTER*12 filename
  INTEGER iz,it
  REAL*8, PARAMETER :: fact=1
  REAL*8 zeta,theta,phi1

  filename="ph1_2D.in"
  OPEN(unit=999+myrank,file=filename,action='write')
  WRITE(999+myrank,*) fact
  DO iz=1,nalphab
     DO it=1,nalphab
        zeta=((iz-1.)*TWOPI/nzperiod)/nalphab
        theta=(it-1.)*TWOPI/nalphab
        phi1=0.5*Ti*cos((theta-nzperiod*zeta))
!        phi1=0.0
        WRITE(999+myrank,*) zeta,theta,phi1
        
     END DO
  END DO
  CLOSE(999+myrank)

END SUBROUTINE WRITE_PHI1

                
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE READ_BULKSPECIES(nalphab,filename,Q,fact)

!----------------------------------------------------------------------------------------------- 
!Read quantity Q (varphi1, Mbb, trM) from file filename and write it in a 
!nalphab x nalphab regular grid. Write in file divided by fact
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER nalphab
  CHARACTER*3 filename
  REAL*8 fact
  !Output
  REAL*8 Q(nalphab,nalphab)
  !Others
  CHARACTER*11 name
  CHARACTER*100 serr
  INTEGER iz,it,jz,jt,nz,nt,nzp,iostat
  REAL, ALLOCATABLE :: Qr(:,:)
  REAL*8 s,dummy

  !Read quantity in EUTERPE format
  name=filename//"_2d.in"
  OPEN(unit=1,file=name,action='read',iostat=iostat)
  IF(iostat.EQ.0) THEN
     WRITE(iout,*) 'File ',name,' read'
     name=filename//"_2d.read"
     OPEN(unit=100+myrank,file=name,form='formatted',action='write',iostat=iostat)
     READ(1,*) serr
     READ(1,*) serr
     READ(1,*) serr
     READ(1,*) serr
     READ(1,*) s
     READ(1,*) nt
     READ(1,*) nz
     READ(1,*) nzp
     nt=nt-1
     nz=nz-1
     ALLOCATE(Qr(nz,nt))
     Qr=0
     IF(nz.LT.nalphab.OR.nt.LT.nalphab) THEN
        IF(nz.LT.0) THEN
           WRITE(iout,*) 'varphi_1 set to zero'
           RETURN
        ELSE
           serr="Error in varphi_1 file"
           CALL END_ALL(serr,.FALSE.)
        END IF
     END IF
     READ(1,*) (dummy,iz=1,nz+1)
     DO it=nt,1,-1
        READ(1,*) (Qr(iz,it),iz=1,nz,1)
     END DO
     CLOSE(1)
     Q=0
     DO iz=1,nz
        IF(MOD(iz-1,nz/nalphab).NE.0) CYCLE
        jz=(iz-1)/(nz/nalphab)+1
        DO it=1,nt
           IF(MOD(it-1,nt/nalphab).NE.0) CYCLE
           jt=(it-1)/(nt/nalphab)+1
           Q(jz,jt)=Qr(iz,it)
           WRITE(100+myrank,'(4(1pe13.5))') s,jz*(TWOPI/nzp)/nz,jt*TWOPI/nt,Q(jz,jt)/fact
        END DO
     END DO
     CLOSE(100+myrank)
  END IF

  IF(ALLOCATED(Qr)) DEALLOCATE(Qr)

END SUBROUTINE READ_BULKSPECIES

                 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE WRITE_BULKSPECIES(s,nalphab,Q,filename)

!----------------------------------------------------------------------------------------------- 
!Write quantity Q (varphi1, Mbb, trM) in s and a nalphab x nalphab regular grid into 
!file filename
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER nalphab
  REAL*8 s,Q(nalphab,nalphab)
  CHARACTER*3 filename
  !Others
  CHARACTER*12 name
  INTEGER ia,il

  name=filename//"_2d.out"
  IF(numprocs.GT.1) WRITE(name,'(A9,".",I2.2)') filename,myrank

  OPEN(unit=999+myrank,file=name,form='formatted',action='write')
  WRITE(999+myrank,*) '##################################################'
  WRITE(999+myrank,*) '#WARNING: these are left-handed Boozer coordinates'
  WRITE(999+myrank,*) '#(similar, but not equal, to PEST)'
  WRITE(999+myrank,*) '##################################################'
  WRITE(999+myrank,*) s
  WRITE(999+myrank,*) nalphab+1
  WRITE(999+myrank,*) nalphab+1
  WRITE(999+myrank,*) nzperiod
  WRITE(999,1001) (Q(il,1),il=1,nalphab),Q(1,1)
  DO ia=nalphab,1,-1
     WRITE(999,1001) (Q(il,ia),il=1,nalphab),Q(1,ia)
  END DO
1001 FORMAT(257(1pe13.5))
  CLOSE(999)

END SUBROUTINE WRITE_BULKSPECIES


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE PLOT_FLUX(jt,jt0,nbb,s,Epsi,Gb,Qb,L1b,L2b,L3b,Zb,nb,dnbdpsi,Tb,dTbdpsi,ephi1oTsize,&
           & Grhs,Qrhs,L1rhs,L2rhs,L3rhs,Gphi,Qphi,L1phi,L2phi,L3phi)

!----------------------------------------------------------------------------------------------- 
!
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER jt,jt0,nbb
  REAL*8 s,Epsi,Gb(nbb),Qb(nbb),L1b(nbb),L2b(nbb),L3b(nbb),Zb(nbb),nb(nbb),dnbdpsi(nbb),Tb(nbb),dTbdpsi(nbb)
  REAL*8 ephi1oTsize,Grhs(nbb),Qrhs(nbb),L1rhs(nbb),L2rhs(nbb),L3rhs(nbb)
  REAL*8 Gphi(nbb,Nnmp),Qphi(nbb,Nnmp),L1phi(nbb,Nnmp),L2phi(nbb,Nnmp),L3phi(nbb,Nnmp)
  !Others
  INTEGER nm,ib
  
  IF(SOLVE_AMB.AND.jt.EQ.jt0) WRITE(300+myrank,'(30(1pe13.5))') s,Epsi*psip,&  
       & (Gb(ib)/psip,Qb(ib)/psip,L2b(ib)/psip/psip,L3b(ib)/psip/psip, &
       & nb(ib),dnbdpsi(ib)/nb(ib)*psip,Tb(ib),dTbdpsi(ib)/Tb(ib)*psip,&
       & Zb(ib),ib=1,MIN(2,NBB)),ephi1oTsize,psip*sgnB*iota*Zb(ib),&
       & PI*vth(ib)*vth(ib)*vth(ib)*m_e*m_e/(16.*aiota*rad_R*borbic(0,0)*borbic(0,0)*Zb(ib)*Zb(ib))
  IF(COMPARE_MODELS) WRITE(4300+myrank,'(I4,30(1pe13.5))') jt,Epsi*psip,&
       & (Gb(ib)/psip,Qb(ib)/psip,L2b(ib)/psip/psip,L3b(ib)/psip/psip, &
       & nb(ib),dnbdpsi(ib)/nb(ib)*psip,Tb(ib),dTbdpsi(ib)/Tb(ib)*psip,&
       & Zb(ib),ib=1,MIN(2,NBB)),ephi1oTsize,psip*sgnB*iota*Zb(ib),&
       & PI*vth(ib)*vth(ib)*vth(ib)*m_e*m_e/(16.*aiota*rad_R*borbic(0,0)*borbic(0,0)*Zb(ib)*Zb(ib))
  IF(TASK3D) THEN
     WRITE(5300+myrank,'(1000(1pe13.5))') SQRT(s),Epsi*psip/1E3,& 
          & nb(1)*1E19*Gb(1)/psip,nb(2)*1E19*Gb(2)/psip,ZERO,nb(2)*1E19*Gb(2)/psip+ZERO,&
          & 1.602*nb(1)*Tb(1)*Qb(1)/psip,1.602*nb(2)*Tb(2)*Qb(2)/psip,ZERO,1.602*nb(2)*Tb(2)*Qb(2)/psip
  END IF
  IF(jt.EQ.jt0.AND.SOLVE_QN) THEN
     nm=1
     WRITE(700+myrank,'(2(1pe13.5),I5,1000(1pe13.5))') &
          & s,Epsi*psip,INT((nm-1)/Nnm),ext_np(nm),ext_mp(nm),&
          & (Grhs(ib)/psip,Qrhs(ib)/psip,&
          & L2rhs(ib)/psip/psip,L3rhs(ib)/psip/psip,ib=1,MIN(2,NBB))
     DO nm=2,Nnmp
        WRITE(700+myrank,'(2(1pe13.5),I5,1000(1pe13.5))') &
             & s,Epsi*psip,INT((nm-1)/Nnm),ext_np(nm),ext_mp(nm),&
             & (Gphi(ib,nm)/psip,Qphi(ib,nm)/psip,&
             & L2phi(ib,nm)/psip/psip,L3phi(ib,nm)/psip/psip,ib=1,MIN(2,NBB))
     END DO
  END IF
  
END SUBROUTINE PLOT_FLUX


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
