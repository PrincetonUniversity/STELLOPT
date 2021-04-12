
!Solve transport

SUBROUTINE TRANSPORT(nbb,ns,dt,s,Zb,Ab,regb,nb,dnbdpsi,Gb,Sb,Tb,dTbdpsi,Qb,Pb,Epsi)

!----------------------------------------------------------------------------------------------- 
!For densities and temperatures and electrostatic potential given by nb, dnbdpsi, Tb, dTbdpsi,
!and Er for nbb species of charge Zb and mass Ab at ns surfaces s, calculate evolution 
!after time step dt
!----------------------------------------------------------------------------------------------- 
  USE GLOBAL
  USE KNOSOS_STELLOPT_MOD
  IMPLICIT NONE
  !Input
  INTEGER nbb,ns,regb(nbb)
  REAL*8 dt,s(ns),Zb(nbb),Ab(nbb),Gb(nbb,ns),Sb(nbb,ns),Qb(nbb,ns),Pb(nbb,ns)
  !Input/output
  REAL*8 nb(nbb,ns),dnbdpsi(nbb,ns),Tb(nbb,ns),dTbdpsi(nbb,ns),Epsi(ns)
  !Others
  REAL*8, PARAMETER  ::  prefact_Epsi=1213173.45142083 !e/(8pi^2m)
  REAL*8, PARAMETER :: e4omemieps02=5516.42011 !units are m^6/s^4
  REAL*8, PARAMETER :: ne_edge=0.1
  REAL*8, PARAMETER :: Ti_edge=50.0
  REAL*8, PARAMETER :: Te_edge=50.0
  REAL*8, PARAMETER :: cgb=10
  REAL*8, PARAMETER :: fgb=0.000288993908
  INTEGER ib,is!,iostat
  REAL*8 dsdV(ns),dVdpsi(ns),dsdr(ns),psi(ns)
  REAL*8 ds,Dgb(ns),Dnc(ns),dlnnbds(ns),expo(ns),nI(ns)
!!  REAL*8 extnb(nbb,ns+1),extTb(nbb,ns+1)
!  REAL*8 dGbdpsi(nbb,ns),dQbdpsi(nbb,ns),dEpsids(ns),dErdr(ns)
!  REAL*8 Pin(ns),ohm(nbb,ns),coulomb(nbb,ns),nue(ns),loglambda(ns)
!  REAL*8 dnbdt(nbb,ns),dTbdt(nbb,ns),dnbTbdt(nbb,ns),dEpsidt(ns)
!  REAL*8 fact_Epsi(ns)
#ifdef MPIandPETSc
  INCLUDE "mpif.h"

  IF(numprocs.GT.1) THEN
     DO ib=1,nbb
        CALL REAL_ALLREDUCE(nb(ib,:),ns)
        CALL REAL_ALLREDUCE(Gb(ib,:),ns)
        CALL REAL_ALLREDUCE(Sb(ib,:),ns)
        CALL REAL_ALLREDUCE(Tb(ib,:),ns)
        CALL REAL_ALLREDUCE(Qb(ib,:),ns)
        CALL REAL_ALLREDUCE(Pb(ib,:),ns)
     END DO
     CALL REAL_ALLREDUCE(Epsi,ns)
  END IF

#endif

  dt=dt
  dnbdpsi=dnbdpsi
  
  IF(.NOT.SS_IMP) RETURN
  
  dsdV=1./(TWOPI*PI*rad_R*rad_a*rad_a)
  dVdpsi=1./(dsdV*atorflux)
  dsdr=2*SQRT(s)/rad_a
  psi=atorflux*s


  DO ib=nbulk+1,nbb
     IF(regb(ib).EQ.2) THEN
        delta=1.5
     ELSE
        delta=3.5
     END IF    
     Dgb(:)=fgb*cgb*((Tb(ib,:)/Ab(ib))**1.5/(Zb(ib)*borbic(0,0)/Ab(ib))**2)*rad_a/rad_R/rad_R
     Dnc(:)=Gb(ib,:)/(Zb(ib)*Epsi(:)/Tb(ib,:))
     dlnnbds(:)=(Zb(ib)*Epsi(:)/Tb(ib,:)-delta*dTbdpsi(ib,:)/Tb(ib,:))*Dnc(:)/(Dnc(:)+Dgb(:))
     !Integral
     expo=0
     expo(1)=dlnnbds(1)*s(1)/2.
     DO is=2,ns
        ds=s(is)-s(is-1)
        expo(is)=expo(is)+0.5*(dlnnbds(is)+dlnnbds(is-1))*ds
     END DO
     nI(:)=EXP(expo(:))/EXP(expo(ns))
     DO is=1,ns
        WRITE(7000+myrank,'(1000(1pe13.5))') s,Epsi*psip,Tb(ib,is),dlnnbds(is),nI(is)
     END DO
  END DO
 
!!$  IF(myrank.EQ.0.AND.KN_STELLOPT(x)) THEN
!!$     Pin=0
!!$     DO ib=1,nbb
!!$        Pin=Pin+Qb(ib,:)*nb(ib,:)*Tb(ib,:)*1.60218*dVdpsi
!!$     END DO
!!$     CALL DERIVE(s,Epsi,ns,4,dEpsids)
!!$     dErdr=dEpsids*dsdr*dsdr*atorflux
!!$     OPEN(unit=7000+myrank,file="flux.opt",form='formatted',action='write',iostat=iostat)
!!$     WRITE(7000+myrank,'("s P_in[W] MAX(P_in)[W] E_r[V/m] dE_r/dr[V/m^2] MAX|dE_r/dr|[V/m^2]")')
!!$     DO is=1,ns
!!$        WRITE(7000+myrank,'(30(1pe13.5))') s(is),Pin(is),MAXVAL(Pin),&
!!$             & Epsi(is)*atorflux*dsdr(is),dErdr(is),MAXVAL(ABS(dErdr))
!!$     END DO
!!$  END IF
!!$
!!$  WRITE(iout,*) 'SUBROUTINE TRANSPORT not ready'
!!$  RETURN
!!$
!!$  DO ib=1,nbb
!!$     CALL DERIVE(psi,Gb(ib,:),ns,4,dGbdpsi(ib,:))
!!$     CALL DERIVE(psi,Qb(ib,:),ns,4,dQbdpsi(ib,:))
!!$  END DO
!!$
!!$  !Particle balance
!!$  DO ib=1,2
!!$     Sb(ib,:)=-dGbdpsi(ib,:) !No particle transport
!!$  END DO
!!$  DO ib=3,nbb
!!$     DO is=1,ns
!!$        CALL READ_ADAS(Zb(ib),Ab(ib),nb(1,is),Tb(1,is),Sb(ib,is))
!!$     END DO
!!$  END DO
!!$  dnbdt  = -dGbdpsi+Sb
!!$
!!$  !Energy balance
!!$  !Ohmic term
!!$  ohm(1,:)=-Epsi*Gb(1,:)*nb(1,:)*dVdpsi
!!$  ohm(2,:)=-ohm(1,:)
!!$  !Collisional transfer
!!$  logLambda=24.0-LOG(SQRT(nb(1,:)*1.0E13)/Tb(1,:))
!!$  nue=e4omemieps02*nb(1,:)*1E19*logLambda/(3.*PI*SQPI*vth(1)*vth(1)*vth(1))
!!$  coulomb(1,:)=3*nb(2,:)*(Tb(2,:)-Tb(1,:))*1.602*nue
!!$  coulomb(2,:)=-coulomb(1,:)
!!$  !
!!$  dnbTbdt=(-dQbdpsi+Pb+ohm+coulomb)*2./3
!!$  dTbdt  =(dnbTbdt-Tb*dnbdt)/nb
!!$  
!!$  !Radial electric field evlution
!!$  dEpsidt=0
!!$  DO ib=1,nbb
!!$     dEpsidt=dEpsidt+Zb(ib)*nb(ib,:)*Gb(ib,:)
!!$  END DO
!!$  fact_Epsi=1.602*torflux*dsdr/(4*PI*PI*Ab(2)*nb(1,:)*etet)
!!$  dEpsidt=dEpsidt*fact_Epsi
!!$
!!$  !Update profiles
!!$  nb=nb+dt*dnbdt
!!$  Tb=Tb+dt*dTbdt
!!$  nb(:,ns+1)    =ne_edge
!!$  Tb(1,ns+1)    =Te_edge
!!$  Tb(2:nbb,ns+1)=Ti_edge
!!$!  extpsi(1)     =0.0
!!$!  extpsi(2:ns+1)=psi(1:ns)
!!$!  extpsi(ns+2)  =atorflux
!!$!  DO ib=1,nbb
!!$!     extnb(ib,1:ns)=nb(ib,:)
!!$!     extTb(ib,1:ns)=Tb(ib,:)
!!$!     CALL DERIVE(extpsi,extnb(ib,:),ns,4,dnbdpsi(ib,:))
!!$!     CALL DERIVE(extpsi,extTb(ib,:),ns,4,dTbdpsi(ib,:))
!!$!  END DO
!!$  Epsi=Epsi+dt*dEpsidt    

END SUBROUTINE TRANSPORT


!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$
!!$SUBROUTINE DERIVE(s,p,ns,order,dpdpsi)
!!$ 
!!$  IMPLICIT NONE
!!$  
!!$  INTEGER order,ns
!!$  REAL*8 s(ns),p(ns),dpdpsi(ns)
!!$
!!$  INTEGER is
!!$  REAL*8 dp(ns),dpsi(ns)
!!$  REAL*8 c21,c41,c42,c61,c62,c63,c81,c82,c83,c84
!!$  REAl f20,f21,f22,f40,f41,f42,f43,f44,f60,f61,f62,f63,f64,f65,f66
!!$
!!$  c21=+1/2.
!!$  c41=+2/3.
!!$  c42=-1/12.
!!$  c61=+3/4.
!!$  c62=-3/20.
!!$  c63=+1/60.
!!$  c81=+4/5.
!!$  c82=-1/5.
!!$  c83=+4/105.
!!$  c84=-1/280.
!!$
!!$  f20=-3/2.
!!$  f21=+2.
!!$  f22=-1/2.
!!$  f40=-25/12.
!!$  f41=+4.
!!$  f42=-3.
!!$  f43=+4/3.
!!$  f44=-1/4.
!!$  f60=-49/20.
!!$  f61=+6.
!!$  f62=-15/2.
!!$  f63=+20/3.
!!$  f64=-15/4.
!!$  f65=+6/5.
!!$  f66=-1/6.
!!$
!!$  dpsi=s(ns)-s(ns-1)
!!$  dp(1) =p(2)-p(1)
!!$  dp(ns)=p(ns)-p(ns-1)
!!$  dpsi(1)=s(2)-s(1)
!!$  dpsi(ns)=s(ns)-s(ns-1)
!!$  DO is=1,ns
!!$     IF(is.NE.ns) dpsi(is)=(s(is+1)-s(is-1))/2.
!!$     IF(order.EQ.6) THEN
!!$        IF((is.GE.4).AND.(is.LE.(ns-3))) THEN
!!$           dp(is)=-c63*p(is-3)-c62*p(is-2)-c61*p(is-1)+c61*p(is+1)+c62*p(is+2)+c63*p(is+3)
!!$        ELSE IF (is.LT.4) THEN
!!$           dp(is)=+f60*p(is)+f61*p(is+1)+f62*p(is+2)+f63*p(is+3)+f64*p(is+4)+f65*p(is+5)+f66*p(is+6)
!!$        ELSE
!!$           dp(is)=-f60*p(is)-f61*p(is-1)-f62*p(is-2)-f63*p(is-3)-f64*p(is-4)-f65*p(is-5)-f66*p(is-6)
!!$        END IF
!!$     ELSE IF(order.EQ.4) THEN
!!$        IF((is.GE.3).AND.(is.LE.(ns-2))) THEN
!!$           dp(is)=-c42*p(is-2)-c41*p(is-1)+c41*p(is+1)+c42*p(is+2)
!!$        ELSE IF (is.LT.3) THEN
!!$           dp(is)=+f40*p(is)+f41*p(is+1)+f42*p(is+2)+f43*p(is+3)+f44*p(is+4)
!!$        ELSE
!!$           dp(is)=-f40*p(is)-f41*p(is-1)-f42*p(is-2)-f43*p(is-3)-f44*p(is-4)
!!$        END IF
!!$     ELSE IF(ABS(order).EQ.2) THEN
!!$        IF((is.GE.2).AND.(is.LE.(ns-1))) THEN
!!$           dp(is)=-c21*p(is-1)+c21*p(is+1)
!!$        ELSE IF(is.EQ.1) THEN
!!$           IF(order.LT.0) THEN
!!$              dp(is)=p(is+1)-p(is)
!!$           ELSE
!!$              dp(is)=+f20*p(is)+f21*p(is+1)+f22*p(is+2)
!!$           END IF
!!$        ELSE
!!$           dp(is)=-f20*p(is)-f21*p(is-1)-f22*p(is-2)
!!$        END IF
!!$     ELSE IF((order.EQ.0).AND.(is.LT.ns)) THEN
!!$        dp(is)=p(is+1)-p(is)
!!$     END IF
!!$  END DO
!!$  dpdpsi=dp/dpsi
!!$  
!!$END SUBROUTINE DERIVE
!!$
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$SUBROUTINE READ_ADAS(Zb,Ab,nb,Tb,Sb)
!!$
!!$!----------------------------------------------------------------------------------------------- 
!!$!
!!$!----------------------------------------------------------------------------------------------- 
!!$  USE GLOBAL
!!$  IMPLICIT NONE
!!$  !Input
!!$  INTEGER ns
!!$  REAL*8 Zb,Ab,nb,Tb
!!$  !Output
!!$  REAL*8 Sb
!!$
!!$  ns=ns
!!$  Zb=Zb
!!$  Ab=Ab
!!$  nb=nb
!!$  Tb=Tb
!!$  Sb=Sb
!!$
!!$END SUBROUTINE READ_ADAS
