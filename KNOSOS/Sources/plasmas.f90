
!Read the kinetic profiles of every species (and electrostatic potential) and
!prepare the integrals in v

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE READ_PLASMAS(nbb,fracb,s0,ZB,AB,nb,dnbdpsi,Tb,dTbdpsi,Epsi)

!-------------------------------------------------------------------------------------------------
!For nbb species of charge ZB(nbb), read density, temperature profiles:
!nb, Tb and radial derivatives, dnbdpsi, dTbdpsi and Epsi
!-------------------------------------------------------------------------------------------------

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER nbb
  REAL*8 fracb(nbb),s0,ZB(nbb),AB(nbb)
  !Output
  REAL*8 nb(nbb),dnbdpsi(nbb),Tb(nbb),dTbdpsi(nbb),Epsi
  !Others
  CHARACTER*100 serr
  INTEGER ib,iostat
  REAL*8 varphi0,ne,dnedpsi,Te,dTedpsi,Ti,dTidpsi,dPhidpsi,arg,expo,tanhh,RGAUSS

  !Read electron density ( and set n_i=n_e, corrections due to ZEFF and Z_i!=1 come later)
  CALL READ_PROFILE(s0,"ne",ne,dnedpsi,nbb)
  nb(1)       =ne     *(1+RGAUSS(iperr))
  dnbdpsi(1:2)=dnedpsi*(1+RGAUSS(iperr))
  nb(2)     =nb(1)        
  dnbdpsi(2)=dnbdpsi(1)
  !If there are more than two species, look for the effective charge 
  IF(nbb.GE.3) THEN
     DO ib=3,nbb   !Incorrect if nbb>3, especially if not trace
        IF(ib.EQ.3) THEN
           nb(ib)=ne*fracb(ib)
           dnbdpsi(ib)=dnedpsi*fracb(ib)
!           nb(ib)     =     ne*(ZEFF-1.)/ZB(ib)/(ZB(ib)-1)
!           dnbdpsi(ib)=dnedpsi*(ZEFF-1.)/ZB(ib)/(ZB(ib)-1)
        ELSE
           nb(ib)=ne/1000.
           dnbdpsi(ib)=dnedpsi/1000.
        END IF
        IF(ib.GT.nbulk.AND.SS_IMP) dnbdpsi(ib)=0        
     END DO
  END IF
  !Impose quasineutrality
  DO ib=3,nbb
     nb(2)     =     nb(2)-ZB(ib)*nb(ib)
     dnbdpsi(2)=dnbdpsi(2)-ZB(ib)*dnbdpsi(ib)
  END DO

  nb(2)     =     nb(2)/ZB(2)
  dnbdpsi(2)=dnbdpsi(2)/ZB(2)

  !Check if impurities are trace or not
!  DO ib=3,nbb
!     IF(SQRT(Ab(ib)/Ab(2))*(nb(ib)*ZB(ib)*ZB(ib))/(nb(2)*ZB(2)*ZB(2)).GT.5E-2) THEN
!        WRITE(iout,'(" ZEFF=",F6.3," (n_z/n_i)(Z_z/Z_i)^2(m_z/m_i)^{1/2}=",F6.3)') &
!             ZEFF,SQRT(Ab(ib)/Ab(2))*(nb(ib)*ZB(ib)*ZB(ib))/(nb(2)*ZB(2)*ZB(2))
!       serr="Not trace impurity"
!        IF(REGB(2).GT.0) THEN
!           CALL END_ALL(serr,.FALSE.)
!        ELSE
!           WRITE(iout,*) serr
!        END IF
!     END IF
!  END DO
!  IF(ne.LT.0) THEN
!     IF(ONLY_DB) THEN
!        ne=-ne
!     ELSE
!        serr="Input plasma profiles not read"
!        CALL END_ALL(serr,.FALSE.)
!     END IF
!  END IF

  !Read temperature profiles
  CALL READ_PROFILE(s0,"te",Te,dTedpsi,nbb)
  CALL READ_PROFILE(s0,"ti",Ti,dTidpsi,nbb)
  IF(Te.LT.0) THEN
     WRITE(iout,*) 'No electron temperature available, setting T_e=T_i'
     Te=Ti
     dTedpsi=dTidpsi
  END IF 
  Tb(1)     =Te     *(1.0+RGAUSS(iperr))
  dTbdpsi(1)=dTedpsi*(1.0+RGAUSS(iperr))
  !Impurities have the same temperature than the bulk ions
  Tb(2:nbb)     =Ti     *(1.0+RGAUSS(iperr))
  dTbdpsi(2:nbb)=dTidpsi*(1.0+RGAUSS(iperr))
  IF(nbb.GT.nbulk.AND.SS_IMP) dTbdpsi(nbulk+1:nbb)=0        
  !Read electrostatic potential (varphi0 is ignored)
  IF(.NOT.SOLVE_AMB) THEN
     IF(TRIVIAL_AMB) THEN
        !Ion root solution with ions in the 1/nu regime
        WRITE(iout,*) 'Ambipolarity is not solved:  &
             & ion root solution with ions in the 1/nu regime'
        dPhidpsi=-(3.5*dTbdpsi(2)+Tb(2)*dnbdpsi(2)/nb(2))/ZB(2) 
     ELSE
        CALL READ_PROFILE(s0,"ph",varphi0,dPhidpsi,nbb)  
     END IF
     Epsi=-dPhidpsi*(1.0+RGAUSS(iperr))
  END IF
  !The following lines correspond to specific physics studies:
  !------------------------------------------------------------------------------------------- 
  !For [Calvo 2018 JPP], use B0 from [Landreman 2011 PoP] and scan in aspect ratio
  !------------------------------------------------------------------------------------------- 
  IF(JPP) THEN

     OPEN(unit=1,file='ne_in.d',action='read',iostat=iostat)
     READ(1,*,iostat=iostat) nb(1)
     nb(2)=nb(1)
     dnbdpsi(1:2)=-nb(1:2)/(2.*SQRT(s0)*atorflux)
     CLOSE(1)
     OPEN(unit=1,file='ti_in.d',action='read',iostat=iostat)
     READ(1,*,iostat=iostat) Tb(2)
     Tb(1)=Tb(2)
     dTbdpsi(1:2)=0
     CLOSE(1)
     IF(ABS(Epsi).GT.1E-8) Epsi=-2*Tb(2)/(2.*SQRT(s0)*atorflux)

  !------------------------------------------------------------------------------------------- 
  !For benchmarking with FORTEC-3D, tokamak perturbed with a gaussian profile
  !------------------------------------------------------------------------------------------- 
  ELSE IF(SATAKE) THEN

     tanhh=TANH((sqrt(s0)-0.5)/0.2)
     arg=-0.2*6*rad_a/rad_R*tanhh
     expo=EXP(arg)
     nb(1:2)=5.045644*expo
     Tb(2)=1877.628*expo
     dnbdpsi(1:2)=-nb(1:2)*6*rad_a/rad_R*(1-tanhh*tanhh)
     dTbdpsi(2)  =-Tb(2)  *6*rad_a/rad_R*(1-tanhh*tanhh)

     tanhh=TANH((sqrt(s0)-0.5)/0.4)
     arg=-0.4*6*rad_a/rad_R*tanhh
     expo=EXP(arg)
     Tb(1)=1949.723*expo
     dTbdpsi(1)  =-Tb(1)  *6*rad_a/rad_R*(1-tanhh*tanhh)
 
     Epsi=0
     OPEN(unit=1,file='ph_in.d',action='read',iostat=iostat)
     IF(iostat.EQ.0) THEN
        READ(1,*) varphi0
        READ(1,*) varphi0,Epsi
        CLOSE(1)
     END IF
     Epsi=Epsi*1E3/psip
     dnbdpsi=dnbdpsi/rad_a/psip  
     dTbdpsi=dTbdpsi/rad_a/psip  

  END IF

  nb        =nb        *FNE
  dnbdpsi   =dnbdpsi   *FNE*FDLNE
  Tb(1)     =Tb(1)     *FTE
  dTbdpsi(1)=dTbdpsi(1)*FTE*FDLTE
  Tb(2)     =Tb(2)     *FTI
  dTbdpsi(2)=dTbdpsi(2)*FTI*FDLTI
  Epsi      =Epsi      *FER

END SUBROUTINE READ_PLASMAS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE READ_PROFILE(s0,filename,q,dqdpsi,nbb)

!-------------------------------------------------------------------------------------------------
!Look for file filename and read value of q and dqdpsi at s0
!(if input in EUTERPE style, using the number of species nbb)  
!-------------------------------------------------------------------------------------------------

  USE GLOBAL
  IMPLICIT NONE
  !Input
  CHARACTER*2 filename
  INTEGER nbb
  REAL*8 s0
  !Output
  REAL*8 q,dqdpsi
  !Others  
  INTEGER, PARAMETER :: ns0=1000
  INTEGER, PARAMETER :: nc=11
  REAL*8,  PARAMETER :: n0   = -10
  REAL*8,  PARAMETER :: dndpsi0=-1
  REAL*8,  PARAMETER :: T0   = 1
  REAL*8,  PARAMETER :: dTdpsi0=0  
  INTEGER iostat,is,ns,iprof,nprof,ic
  CHARACTER*100 serr
  CHARACTER*10 nameeut
  CHARACTER*12 namegraz
  CHARACTER*13 namegraz2
  CHARACTER*11 nametg
  CHARACTER*13 nametask3d
  CHARACTER*14 nametask3d2
  CHARACTER*15 namedr
  CHARACTER*15 nametg2
  CHARACTER*22 nameneo
  CHARACTER*10 nameprof
  CHARACTER*7 namepoint
  CHARACTER*200 line
  REAL*8 r(ns0),rho(ns0),s_pol(ns0),dummy,q_p(ns0),dqdx_p(ns0),dlnqdx_p(ns0),s(ns0),dqdpsi_p(ns0)
  REAL*8 sqs,coeff(0:nc-1)

  q=-1.0
  IF(filename.EQ."ph") dqdpsi=0

  !TASK3D profile format
  nametask3d2="input-prof.txt"
  OPEN(unit=1,file=nametask3d2,action='read',iostat=iostat)
  IF(iostat.NE.0) OPEN(unit=1,file="../"//nametask3d2,action='read',iostat=iostat)
  IF(iostat.EQ.0.AND.filename.NE."ph") THEN
     WRITE(iout,'(" File ",A22," read")') nametask3d2
     IF(filename.EQ."ti") THEN
        nprof=2
     ELSE IF(filename.EQ."te") THEN
        nprof=3
     ELSE IF(filename.EQ."ne") THEN
        nprof=1
     END IF
     DO iprof=1,nprof
        READ(1,*) line
        READ(1,*,iostat=iostat) (coeff(ic),ic=0,nc-1)
     END DO
     CLOSE(1)
     sqs=sqrt(s0)
     q=coeff(0)+sqs*(coeff(1)+sqs*(coeff(2)+sqs*(coeff(3)+sqs*(coeff(4)+sqs*(coeff(5)+&
      & sqs*(coeff(6)+sqs*(coeff(7)+sqs*(coeff(8)+sqs*(coeff(9)+sqs*(coeff(10)))))))))))
     dqdpsi=coeff(1)+sqs*(2*coeff(2)+sqs*(3*coeff(3)+sqs*(4*coeff(4)+sqs*(5*coeff(5)+&
      & sqs*(6*coeff(6)+sqs*(7*coeff(7)+sqs*(8*coeff(8)+sqs*(9*coeff(9)+sqs*(10*coeff(10))))))))))&
      & /(2.*SQRT(s0)*atorflux)  !x=rho=r/a         
     IF(filename.NE."ne") THEN
        q     =q     *1.0E3  !Temperatures and radial electric field are read in keV and kV/m
        dqdpsi=dqdpsi*1.0E3
     END IF
!     GOTO 1000    
  END IF
  
  !TASK3D profile format
  nametask3d="prof-file.txt"
  OPEN(unit=1,file=nametask3d,action='read',iostat=iostat)
  IF(iostat.NE.0) OPEN(unit=1,file="../"//nametask3d,action='read',iostat=iostat)
  IF(iostat.EQ.0.AND.filename.NE."ph") THEN
     WRITE(iout,'(" File ",A22," read")') nametask3d
     READ(1,*) line
     READ(1,*) line
     IF(filename.EQ."ti") THEN
        nprof=3
     ELSE IF(filename.EQ."te") THEN
        nprof=2
     ELSE IF(filename.EQ."ne") THEN
        nprof=1!NBB+2
     END IF
     DO is=1,ns0
        READ(1,*,iostat=iostat) s(is),(q_p(is),iprof=1,nprof)
        IF(iostat.NE.0) EXIT
     END DO
     ns=is-1
     CLOSE(1)
     dqdx_p(1)=(-3*q_p(1)+4*q_p(2)-q_p(3))/(s(3)-s(1))
     DO is=2,ns-1
        dqdx_p(is)=(q_p(is+1)-q_p(is-1))/(s(is+1)-s(is-1))
     END DO
     dqdx_p(ns)=(3*q_p(ns)-4*q_p(ns-1)+q_p(ns-2))/(s(3)-s(1))
     dqdpsi_p=dqdx_p/(2.*SQRT(s0)*atorflux)  !x=rho=r/a
     !Interpolate at s0
     CALL LAGRANGE(s,     q_p,ns,s0,     q,MIN(ns,3))
     CALL LAGRANGE(s,dqdpsi_p,ns,s0,dqdpsi,MIN(ns,3))        
     IF(filename.EQ."ne") THEN
        q     =q     /1.0E19 
        dqdpsi=dqdpsi/1.0E19
     ELSE
        q     =q     *1.0E3  !Temperatures and radial electric field are read in keV and kV/m
        dqdpsi=dqdpsi*1.0E3
     END IF
!     GOTO 1000    
  END IF

  !GRAZ profile format
  namegraz="profiles.txt"
  OPEN(unit=1,file=namegraz,action='read',iostat=iostat)
  IF(iostat.NE.0) OPEN(unit=1,file="../"//namegraz,action='read',iostat=iostat)
  IF(iostat.EQ.0) THEN
     WRITE(iout,'(" File ",A22," read")') namegraz
     READ(1,*) line
     IF(filename.EQ."ti") THEN
        nprof=3
     ELSE IF(filename.EQ."te") THEN
        nprof=4
     ELSE IF(filename.EQ."ne") THEN
        nprof=2
     ELSE IF(filename.EQ."ph") THEN
        nprof=5
     END IF
     DO is=1,ns0
        READ(1,*,iostat=iostat) s_pol(is),(q_p(is),iprof=1,nprof)
        IF(iostat.NE.0) EXIT
     END DO
     ns=is-1
     CLOSE(1)
     dqdx_p(1)=(-3*q_p(1)+4*q_p(2)-q_p(3))/(s_pol(3)-s_pol(1))
     DO is=2,ns-1
        dqdx_p(is)=(q_p(is+1)-q_p(is-1))/(s_pol(is+1)-s_pol(is-1))
     END DO
     dqdx_p(ns)=(3*q_p(ns)-4*q_p(ns-1)+q_p(ns-2))/(s_pol(3)-s_pol(1))
     dqdpsi_p=dqdx_p*dspolds/atorflux   !x=spol
     !Interpolate at spol(s0)
     CALL LAGRANGE(s_pol,     q_p,ns,spol,     q,MIN(ns,3))
     CALL LAGRANGE(s_pol,dqdpsi_p,ns,spol,dqdpsi,MIN(ns,3))        
     IF(filename.EQ."ne") THEN
        q     =q     /1.0E19 
        dqdpsi=dqdpsi/1.0E19
     ELSE IF(filename.EQ."ph") THEN
        iota=-iota
        IF(NTV) iota=sgnhel/iota
        dqdpsi=-q*iota   !dPhi/dpsi_T=-iota*omega_E
        iota=-iota
        IF(NTV) iota=sgnhel/iota
        q=0.0
     END IF
!     GOTO 1000    
  END IF
  
  !GRAZ profile format
  namegraz2="profiles2.txt"
  OPEN(unit=1,file=namegraz2,action='read',iostat=iostat)
  IF(iostat.NE.0) OPEN(unit=1,file="../"//namegraz2,action='read',iostat=iostat)
  IF(iostat.EQ.0) THEN
     WRITE(iout,'(" File ",A22," read")') namegraz2
     READ(1,*) line
     IF(filename.EQ."ti") THEN
        nprof=6
     ELSE IF(filename.EQ."te") THEN
        nprof=5
     ELSE IF(filename.EQ."ne") THEN
        nprof=3
     ELSE IF(filename.EQ."ph") THEN
        nprof=11
     END IF
     DO is=1,ns0
        READ(1,*,iostat=iostat) s(is),(q_p(is),iprof=1,nprof)
        IF(iostat.NE.0) EXIT
     END DO
     ns=is-1
     CLOSE(1)
     dqdx_p(1)=(-3*q_p(1)+4*q_p(2)-q_p(3))/(s(3)-s(1))
     DO is=2,ns-1
        dqdx_p(is)=(q_p(is+1)-q_p(is-1))/(s(is+1)-s(is-1))
     END DO
     dqdx_p(ns)=(3*q_p(ns)-4*q_p(ns-1)+q_p(ns-2))/(s(3)-s(1))
     dqdpsi_p=dqdx_p/atorflux      !x=s=rho^2=(r/a)^2
     !Interpolate at spol(s0)
     CALL LAGRANGE(s,     q_p,ns,s0,     q,MIN(ns,3))
     CALL LAGRANGE(s,dqdpsi_p,ns,s0,dqdpsi,MIN(ns,3))        
     IF(filename.EQ."ne") THEN
        q     =q     /1.0E13
        dqdpsi=dqdpsi/1.0E13
     ELSE IF(filename.EQ."ph") THEN
        dqdpsi=-q*29979.19999934/psip  !E_r read in stV/cm
        q=0
     END IF
!     GOTO 1000    
  END IF
  
  !TG profile format
  nametg="profiles.TG"
  OPEN(unit=1,file=nametg,action='read',iostat=iostat)
  IF(iostat.NE.0) OPEN(unit=1,file="../"//nametg,action='read',iostat=iostat)
  IF(iostat.EQ.0) THEN
     WRITE(iout,'(" File ",A22," read")') nametg
     IF(filename.EQ."ti") THEN
        nprof=568
     ELSE IF(filename.EQ."te") THEN
        nprof=568+305
     ELSE IF(filename.EQ."ne") THEN
        nprof=568+305*2
     ELSE IF(filename.EQ."ph") THEN
        nprof=3102
     END IF
     DO is=1,nprof
        READ(1,*) line
     END DO
     IF(filename.EQ."ph") THEN
        ns=50
        DO is=1,ns
           READ(1,*,iostat=iostat) r(is),dummy,dqdx_p(is)
        END DO
        s=r*r/rad_a/rad_a
        q_p=0
        dqdpsi_p=-dqdx_p/psip        
     ELSE  
        ns=300
        DO is=1,ns
           READ(1,*,iostat=iostat) r(is),dummy,q_p(is)
        END DO
        s=r*r/rad_a/rad_a
        dqdx_p(1)=(-3*q_p(1)+4*q_p(2)-q_p(3))/(s(3)-s(1))
        DO is=2,ns-1
           dqdx_p(is)=(q_p(is+1)-q_p(is-1))/(s(is+1)-s(is-1))
        END DO
        dqdx_p(ns)=(3*q_p(ns)-4*q_p(ns-1)+q_p(ns-2))/(s(3)-s(1))

     END IF
     CLOSE(1)
     !Interpolate at s0
     CALL LAGRANGE(s,     q_p,ns,s0,     q,MIN(ns,3))
     CALL LAGRANGE(s,dqdpsi_p,ns,s0,dqdpsi,MIN(ns,3)) 
     IF(filename.NE."ne") THEN
        q     =q     *1.0E3      !Temperatures and radial electric field are read in keV
        dqdpsi=dqdpsi*1.0E3
     END IF
!     GOTO 1000    
  END IF

  !DR@W7-X profile format
  namedr="profiles_DR.txt"
  OPEN(unit=1,file=namedr,action='read',iostat=iostat)
  IF(iostat.NE.0) OPEN(unit=1,file="../"//namedr,action='read',iostat=iostat)
  IF(iostat.EQ.0.AND.filename.NE."ph") THEN
     WRITE(iout,'(" File ",A22," read")') namedr
     IF(filename.EQ."ti") THEN
        nprof=3
     ELSE IF(filename.EQ."te") THEN
        nprof=2
     ELSE IF(filename.EQ."ne") THEN
        nprof=1
     END IF
     READ(1,*) line
     DO is=1,ns0
        READ(1,*,iostat=iostat) rho(is),(q_p(is),iprof=1,nprof)
        IF(iostat.NE.0) EXIT
     END DO
     ns=is-1
     CLOSE(1)
     IF(rho(ns).LT.0.9) rho=rho/rad_a
     s=rho*rho
     dqdx_p(1)=(-3*q_p(1)+4*q_p(2)-q_p(3))/(rho(3)-rho(1))
     DO is=2,ns-1
        dqdx_p(is)=(q_p(is+1)-q_p(is-1))/(rho(is+1)-rho(is-1))
     END DO
     dqdx_p(ns)=(3*q_p(ns)-4*q_p(ns-1)+q_p(ns-2))/(rho(3)-rho(1))
     dqdpsi_p=dqdx_p/(2.*SQRT(s0)*atorflux)  !x=rho=r/a
     !Interpolate at s0
     CALL LAGRANGE(s,     q_p,ns,s0,     q,MIN(ns,3))
     CALL LAGRANGE(s,dqdpsi_p,ns,s0,dqdpsi,MIN(ns,3))    
     IF(filename.NE."ne") THEN
        q     =q     *1.0E3      !Temperatures and radial electric field are read in keV and kV/m
        dqdpsi=dqdpsi*1.0E3
     END IF
!     GOTO 1000    
  END IF
  
  !TG profile format
  nametg2="profiles_TG.txt"
  OPEN(unit=1,file=nametg2,action='read',iostat=iostat)
  IF(iostat.NE.0) OPEN(unit=1,file="../"//nametg2,action='read',iostat=iostat)
  IF(iostat.EQ.0.AND.filename.NE."ph") THEN
     WRITE(iout,'(" File ",A22," read")') nametg2
     IF(filename.EQ."ti") THEN
        nprof=4
     ELSE IF(filename.EQ."te") THEN
        nprof=3
     ELSE IF(filename.EQ."ne") THEN
        nprof=2
     END IF
     DO is=1,15
        READ(1,*) line
     END DO
     DO is=1,ns0
        READ(1,*,iostat=iostat) rho(is),(q_p(is),iprof=1,nprof)
        IF(iostat.NE.0) EXIT
     END DO
     ns=is-1
     CLOSE(1)
     s=rho*rho
     dqdx_p(1)=(-3*q_p(1)+4*q_p(2)-q_p(3))/(rho(3)-rho(1))
     DO is=2,ns-1
        dqdx_p(is)=(q_p(is+1)-q_p(is-1))/(rho(is+1)-rho(is-1))
     END DO
     dqdx_p(ns)=(3*q_p(ns)-4*q_p(ns-1)+q_p(ns-2))/(rho(3)-rho(1))
     dqdpsi_p=dqdx_p/(2.*SQRT(s0)*atorflux)  !x=rho=r/a
     !Interpolate at s0
     CALL LAGRANGE(s,     q_p,ns,s0,     q,MIN(ns,3))
     CALL LAGRANGE(s,dqdpsi_p,ns,s0,dqdpsi,MIN(ns,3))    
     IF(filename.NE."ne") THEN
        q     =q     *1.0E3      !Temperatures and radial electric field are read in keV and kV/m
        dqdpsi=dqdpsi*1.0E3
     END IF
!     GOTO 1000    
  END IF

  !neotransp profile format
  nameneo="profiles_neotransp.dat"
  OPEN(unit=1,file=nameneo,action='read',iostat=iostat)
  IF(iostat.NE.0) OPEN(unit=1,file="../"//nameneo,action='read',iostat=iostat)
  IF(iostat.EQ.0.AND.filename.NE."ph") THEN
     WRITE(iout,'(" File ",A22," read")') nameneo
     IF(filename.EQ."ti") THEN
        nprof=3
     ELSE IF(filename.EQ."te") THEN
        nprof=2
     ELSE IF(filename.EQ."ne") THEN
        nprof=1
     END IF
     READ(1,*) line
     DO is=1,ns0
        READ(1,*,iostat=iostat) r(is),(q_p(is),iprof=1,nprof)
        IF(iostat.NE.0) EXIT
     END DO
     ns=is-1
     CLOSE(1)
     rho=r/rad_a
     s=r*r/rad_a/rad_a
     dqdx_p(1)=(-3*q_p(1)+4*q_p(2)-q_p(3))/(rho(3)-rho(1))
     DO is=2,ns-1
        dqdx_p(is)=(q_p(is+1)-q_p(is-1))/(rho(is+1)-rho(is-1))
     END DO
     dqdx_p(ns)=(3*q_p(ns)-4*q_p(ns-1)+q_p(ns-2))/(rho(3)-rho(1))
     dqdpsi_p=dqdx_p/(2.*SQRT(s0)*atorflux)  !x=rho=r/a
     !Interpolate at s0
     CALL LAGRANGE(s,     q_p,ns,s0,     q,MIN(ns,3))
     CALL LAGRANGE(s,dqdpsi_p,ns,s0,dqdpsi,MIN(ns,3))    
     IF(filename.EQ."ne".AND.(INDEX(line, "n_20m3").EQ.0)) THEN
        q     =q     *10.      !Density read in 10{^20}m^{-3}
        dqdpsi=dqdpsi*10.
     ELSE IF(filename.EQ."te".AND.(INDEX(line, "keV_e").EQ.0).OR.&
          &  filename.EQ."ti".AND.(INDEX(line, "keV_H").EQ.0).OR.&
          &  filename.EQ."ti".AND.(INDEX(line, "keV_D").EQ.0)) THEN
        q     =q     *1.0E3    !Temperatures read in keV
        dqdpsi=dqdpsi*1.0E3
     END IF
!     GOTO 1000    
  END IF

  !EUTERPE profile format
  nameeut="profiles.d"
  OPEN(unit=1,file=nameeut,action='read',iostat=iostat)
  IF(iostat.NE.0) OPEN(unit=1,file="../"//nameeut,action='read',iostat=iostat)
  IF(iostat.EQ.0.AND.filename.NE."ph") THEN
     WRITE(iout,'(" File ",A22," read")') nameeut
     IF(filename.EQ."ti") THEN
        nprof=1
     ELSE IF(filename.EQ."te") THEN
        nprof=2
     ELSE IF(filename.EQ."ne") THEN
        IF(NBB.EQ.3) THEN
           nprof=NBB+2
        ELSE
           nprof=3
        END IF
     END IF
     DO is=1,ns0
        READ(1,*,iostat=iostat) s(is),(dlnqdx_p(is),q_p(is),iprof=1,nprof)
        IF(iostat.NE.0) EXIT
     END DO
     ns=is-1
     CLOSE(1)
     dqdpsi_p=dlnqdx_p*q_p/atorflux      !x=s=rho^2=(r/a)^2
     !Interpolate at s0
     CALL LAGRANGE(s,     q_p,ns,s0,     q,MIN(ns,3))
     CALL LAGRANGE(s,dqdpsi_p,ns,s0,dqdpsi,MIN(ns,3))        
     IF(filename.EQ."ne") THEN
        q     =q     /1.0E19 
        dqdpsi=dqdpsi/1.0E19
     END IF
!     GOTO 1000    
  ELSE IF (filename.EQ."ph") THEN
     nameeut="er_s.dat"
     OPEN(unit=1,file=nameeut,action='read',iostat=iostat)
     IF(iostat.NE.0) OPEN(unit=1,file="../"//nameeut,action='read',iostat=iostat)
     IF(iostat.EQ.0) THEN
        WRITE(iout,'(" File ",A22," read")') nameeut
        READ(1,*) line
        READ(1,*) line
        READ(1,*) line
        READ(1,*) line
        DO is=1,ns0
           READ(1,*,iostat=iostat) s(is),dqdx_p(is)
           IF(iostat.NE.0) EXIT
        END DO
        ns=is-1
        CLOSE(1)
        dqdpsi_p=-dqdx_p/psip !Er=-d_r varphi_0 was read, instead of d_psi varphi_0
        !Interpolate at s0
        CALL LAGRANGE(s,dqdpsi_p,ns,s0,dqdpsi,MIN(ns,3))        
!        GOTO 1000    
     END IF
  END IF

  !neprof profile format
  nameprof=filename//"prof_r.d"
  OPEN(unit=1,file=nameprof,action='read',iostat=iostat)
  IF(iostat.NE.0) OPEN(unit=1,file="../"//nameprof,action='read',iostat=iostat)
  IF(iostat.EQ.0) THEN
     WRITE(iout,'(" File ",A22," read")') nameprof
     READ(1,*) line
     DO is=1,ns0
        READ(1,*,iostat=iostat) rho(is),dummy,q_p(is),dummy,dummy,dqdx_p(is),dummy
        IF(iostat.NE.0) EXIT
     END DO
     ns=is-1
     CLOSE(1)
     s=rho*rho
     IF(filename.EQ."ph") THEN      !Er=-d_r varphi_0 was read, instead of d_psi varphi_0
        dqdpsi_p=-dqdx_p/psip       !x=r=rho*a
        q_p=0
     ELSE
        dqdpsi_p=dqdx_p/(2.*SQRT(s0)*atorflux)  !x=rho=r/a
     END IF
     !Interpolate at s0
     CALL LAGRANGE(s,     q_p,ns,s0,    q,MIN(ns,3))
     CALL LAGRANGE(s,dqdpsi_p,ns,s0,dqdpsi,MIN(ns,3))
     IF(filename.NE."ne") THEN
        q     =q     *1.0E3      !Temperatures and radial electric field are read in keV and kV/m
        dqdpsi=dqdpsi*1.0E3
     END IF
     GOTO 1000
  END IF

  !Local values
  namepoint=filename//"_in.d"
  OPEN(unit=1,file=namepoint,action='read',iostat=iostat)
  IF(iostat.NE.0) OPEN(unit=1,file="../"//namepoint,action='read',iostat=iostat)
  IF(iostat.EQ.0) THEN
     WRITE(iout,'(" File ",A22," read")') namepoint
     READ(1,*) line
     READ(1,*,iostat=iostat) q,dlnqdx_p(1)  !x=r=rho*a
     dqdpsi=q*dlnqdx_p(1)/psip        
     CLOSE(1)
     IF(filename.NE."ne") THEN
        q     =q     *1.0E3      !Temperatures and radial electric field are read in keV and kV/m
        dqdpsi=dqdpsi*1.0E3
     END IF
  END IF

  IF(filename.NE."ph".AND.filename.NE."te".AND.q.LT.ALMOST_ZERO.AND..NOT.SATAKE) THEN
!     WRITE(1000,*) filename,q
     serr="Missing profile"
     CALL END_ALL(serr,.FALSE.)
  END IF
  
1000 IF(filename.EQ."ph") dqdpsi=dqdpsi/FE  !for scan in the aspect ratio         

      
END SUBROUTINE READ_PROFILE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE READ_SOURCES(nbb,s0,Sb,Pb)

!-------------------------------------------------------------------------------------------------
!Set particle and energy source, Sb and Pb, for nbb species at radial position s0
!-------------------------------------------------------------------------------------------------

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER nbb
  REAL*8 s0
  !Output
  REAL*8 Sb(nbb),Pb(nbb)

  REAL*8, PARAMETER :: ws     = 0.1
  REAL*8, PARAMETER :: P_ECH  = 1E6       !Wm^-3
  REAL*8, PARAMETER :: J_to_eV= 1.602e-19 

  WRITE(iout,*) 'SUBROUTINE READ_SOURCES not ready'
  RETURN
  
  Sb=0
  Pb=0
  Pb(1)=P_ECH*J_to_eV*EXP(-s0/ws)
  
END SUBROUTINE READ_SOURCES


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   
SUBROUTINE DKE_CONSTANTS(ib,nbb,ZB,AB,REGB,nb,dnbdpsi,Tb,dTbdpsi,Epsi,flag)

!-------------------------------------------------------------------------------------------------
!For species ib (of nbb) of charge ZB, mass AB, regime REGB, density nb, temperature Tb,
!and radial derivatives dnbdpsi, dTbdpsi and Epsi, calculate velocity-dependent quantities
!source of the DKE S_dke, convolution weight, normalizations, fnorm, fdkes, vdconst=mv^2/qt
!and collision frequency nu (global variables). Plot quantities if flag
!-------------------------------------------------------------------------------------------------

  USE GLOBAL
  IMPLICIT NONE
  !Input
  LOGICAL flag
  INTEGER ib,nbb,regb(nbb)
  REAL*8 ZB(nbb),AB(nbb),nb(nbb),dnbdpsi(nbb),Tb(nbb),dTbdpsi(nbb),Epsi
  !Others
  LOGICAL, SAVE :: REGIMES_CHECKED=.FALSE.
  INTEGER jb,iv
  REAL*8 x2M(nv) /2.24159e-02,1.18123e-01,2.90366e-01,5.39286e-01,8.65037e-01,    &
       & 1.26781e+00,1.74786e+00,2.30546e+00,2.94096e+00,3.65475e+00,4.44727e+00, &
       & 5.31900e+00,6.27050e+00,7.30237e+00,8.41528e+00,9.60994e+00,1.08872e+01, &
       & 1.22478e+01,1.36927e+01,1.52230e+01,1.68397e+01,1.85439e+01,2.03370e+01, & 
       & 2.22202e+01,2.41951e+01,2.62631e+01,2.84260e+01,3.06855e+01/
  REAL*8 wM(nv) /5.62529e-02,1.19024e-01,1.57496e-01,1.67547e-01,1.53353e-01,     &
       & 1.24221e-01,9.03423e-02,5.94778e-02,3.56276e-02,1.94804e-02,9.74360e-03, &
       & 4.46431e-03,1.87536e-03,7.22647e-04,2.55488e-04,8.28714e-05,2.46569e-05, &
       & 6.72671e-06,1.68179e-06,3.85081e-07,8.06883e-08,1.54573e-08,2.70450e-09, &
       & 4.31679e-10,6.27784e-11,8.30642e-12,9.98408e-13,1.08837e-13/
  REAL*8 nu0,coul_log,x(nv),x2(nv),x3(nv),kx(nv),kx2(nv),expke(nv),erfv(nv),nuzistar,aleph
  REAL*8 cmul,nustar,erstar,factPS,rhostar
  
  !Calculate thermal velocity
  vth(ib)=1.389165e4*SQRT(Tb(ib)/AB(ib)) 
  IF(flag) THEN
     WRITE(iout,*) 'Calculating source',ZB(ib)
     WRITE(iout,*) 'Temperature (keV) =',Tb(ib)/1.0E3
     WRITE(iout,*) 'Density (10^{19}m^{-3}) =',nb(ib)
     WRITE(iout,'(" Reference v_{th}=",1pe13.5)') vth(ib)
  END IF
  
  !Calculate prefactor of the collision frequency
  nu0=FN*2.0*0.242594*1E19*ZB(ib)*ZB(ib)/(AB(ib)*AB(ib)*vth(ib)*vth(ib)*vth(ib))

  x2=x2M
  x=SQRT(x2)
  x3=x2*x
  nu=0
  DO jb=1,nbb

     !Calculate Coulomb logarithm
     IF(ZB(ib).LT.0) THEN
        IF(ZB(jb).LT.0) THEN    !ee
           coul_log=23.5-LOG(SQRT(nb(ib)*1.0E13)*(Tb(ib)**(-5.0/4.0)))-&
                & SQRT( 1.0E-5 + ((LOG(Tb(ib))-2.0)**2) / 16.0)
        ELSE                    !ei and iz
           coul_log=24.0-LOG(SQRT(nb(ib)*1.0E13)/Tb(ib))
        END IF
     ELSE
        IF(ZB(jb).LT.0) THEN    !ie and ze
           coul_log=24.0-LOG(SQRT(nb(jb)*1.0E13)/Tb(jb))
        ELSE                    !ii, iz, zi, and zz
           coul_log=23.0-LOG( (ZB(ib)*ZB(jb)*(AB(ib)+AB(jb))/(AB(ib)*Tb(jb)+AB(jb)*Tb(ib)))* &
                &   SQRT(nb(ib)*1.0E13*(ZB(ib)**2)/Tb(ib) + nb(jb)*1.0E13*(ZB(jb)**2)/Tb(jb)) )
        END IF
     END IF
     IF(ONLY_DB) coul_log=20     

     !Factors for interespecies terms
     kx2=x2*(AB(jb)/AB(ib))*(Tb(ib)/Tb(jb))
     kx=SQRT(kx2)
     expke=EXP(-kx2)
     DO iv=1,nv
        erfv(iv)=ERF(kx(iv))
     END DO
     
     !Determine the collisionality regime of each ion species
     IF(ib.EQ.2.AND.jb.GT.2) THEN
        IF(jb.EQ.3) THEN
           WRITE(iout,*) "------------------------------------------------------"
           WRITE(iout,*) "Mixed collisionality?"
        END IF
        nuzi(jb)=(2./3./SQPI)*nu0*nb(ib)*ZB(jb)*ZB(jb)*coul_log*AB(ib)/AB(jb)
        nuzistar=rad_R*nuzi(jb)/(1.389165e4*SQRT(Tb(jb)/AB(jb)))           
        aleph=eps32*Zb(jb)*Zb(jb)*Ab(2)/Ab(jb)
        WRITE(iout,'(" Impurity #",I1," (Z=",f3.0,",A=",f8.4,"): (nu_zi)^*=",1pe13.5,&
             & ",  aleph=",1pe13.5,",  (nu_zi)^*/aleph^{1/2}=",1pe13.5)') jb,ZB(jb),AB(jb),&
             & nuzistar,aleph,nuzistar/SQRT(aleph)
     END IF
     nu=nu+nu0*nb(jb)*ZB(jb)*ZB(jb)*coul_log*(ERFv*(1.-.5/kx2)+expke/(sqpi*kx))/x3
     IF(ib.EQ.2.AND.jb.GT.2.AND.jb.EQ.nbb) THEN
        WRITE(iout,'(" Bulk ions: (nu_i)^*/EPS^{3/2}=",1pe13.5)') &
             & rad_R*nu(iv0)/vth(ib)/eps32*x(iv0)*x2M(iv0)*SQRT(x2M(iv0))
        WRITE(iout,*) "------------------------------------------------------"
!        IF(rad_R*nu(iv0)/vth(ib)/eps32*x(iv0)*x2M(iv0)*SQRT(x2M(iv0)).GT.10) THEN
!           TRACE_IMP=.FALSE.
!           WRITE(1100+myrank,'(" Bulk ions: (nu_i)^*/EPS^{3/2}=",1pe13.5)') &
!                & rad_R*nu(iv0)/vth(ib)/eps32*x(iv0)*x2M(iv0)*SQRT(x2M(iv0))
!        END IF
     END IF

     IF(ib.LE.2) nuth(ib)=nu0*nb(MIN(2,nbb))*ZB(MIN(2,nbb))*ZB(MIN(2,nbb))*coul_log !collision frequency of thermal species

     nustar=rad_R*nu(iv0)/vth(ib)/eps32*x(iv0)*x2M(iv0)*SQRT(x2M(iv0))
     rhostar=vth(ib)*Ab(ib)*m_e/Zb(ib)/borbic(0,0)/rad_R
  END DO

  !Quantities evaluated for all nv velocities
  v=vth(ib)*x
  vdconst=v*v*AB(ib)*m_e/ZB(ib)
  vmconst=vdconst
  Sdke=dnbdpsi(ib)/nb(ib)+dTbdpsi(ib)/Tb(ib)*(x2-1.5)-ZB(ib)*Epsi/Tb(ib)
  ! This can be used to calculate numerically the integral
  !   4\pi\int_0^\inf dv v^2 F_M f(v)=
  ! =(2/\sqrt(\pi))*\int_0^\inf dK K^{1/2}exp{-K}g(K),
  ! with f(v)=g(K),  K=mv^2/(2T) and x=K^{1/2}
  weight=(2./sqpi)*x*wM 
  fdkes=(2*v/vdconst/vdconst)/psip/psip
  fdkes2=2/vdconst/psip
  mu=v*v/2./borbic(0,0)
  ftrace1nu=-Ab(ib)*m_e*Ab(ib)*m_e*SQRT(Ab(ib)*m_e/Tb(ib))/Tb(ib)/Tb(ib)/nuzi(ib)*&
       & mu*SQRT(mu)/2./SQPI*vdconst*vdconst/4/borbic(0,0)/borbic(0,0)
  mmuoT=mu*Ab(ib)*m_e/Tb(ib)

  DO iv=1,nv
     
     IF(ONLY_DB.OR.REGIMES_CHECKED.OR.REGB(ib).EQ.0) EXIT

     IF(ib.LE.2) cmul=nuth(ib)/vth(ib)
     IF(ib.GT.2) cmul=nuzi(ib)/vth(ib)
     nustar=rad_R*cmul/iota
     erstar=Epsi*psip/vth(ib)/borbic(0,0)
     factPS=(3*nustar*Erstar/iota/eps)

     IF(factPS.GT.0.5.AND.regb(ib).GE.10) WRITE(1100+myrank,*) &
          & 'Collisionality is so high that ExB cannot be neglected',myrank,ib,regb(ib),iv
     IF(cmul.GT.cmul_PS.AND.(MOD(regb(ib),10).EQ.2.OR.MOD(regb(ib),10).EQ.3)) &
          & WRITE(1100+myrank,*) 'PS should be solved',myrank,ib,regb(ib),iv
     IF(cmul.GT.cmul_1NU.AND.cmul.LT.cmul_PS.AND.(MOD(regb(ib),10).EQ.1.OR.MOD(regb(ib),10).EQ.3)) &
          & WRITE(1100+myrank,*) 'Plateau should be solved',myrank,ib,regb(ib),iv
     IF(cmul.LT.cmul_1NU.AND.(MOD(regb(ib),10).EQ.1.OR.MOD(regb(ib),10).EQ.2)) &
          & WRITE(1100+myrank,*) 'Low collisionality should be solved',myrank,ib,regb(ib),iv

     REGIMES_CHECKED=.TRUE.
     
  END DO

  
END SUBROUTINE DKE_CONSTANTS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
