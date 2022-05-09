!Read input files

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE READ_INPUT(dt,ns,s,nbb,Zb,Ab,regb,fracb)

!-------------------------------------------------------------------------------------------------
!Read simulation input. This includes the number ns of surfaces s, the number of species nnb and
!their charge number Zb, mass number Ab and NC regime regb, and the ratio of ion species
!The rest are global variables.
!All of them are described in the user manual.  
!-------------------------------------------------------------------------------------------------

  USE GLOBAL
  USE KNOSOS_STELLOPT_MOD
  IMPLICIT NONE

  !Namelist 'transport' contains parameters of the
  !solution of theenergy transport equation
  NAMELIST /transport/ dt
  !Namelist 'parameters' contains simulation parameters, including
  !decisions on how to solve some equations and what to plot
  NAMELIST /parameters/  GEN_FLAG,DEBUG,TIME,I0,L0,PLOTG,&
       & USE_B1,USE_B0pB1,QS_B0,QS_B0_1HEL,USE_SHAING,SHAING_1NU,SHAING_SQRTNU,SHAING_SBP,&
       & REMOVE_DIV,DELTA,&
       & NTURN,MLAMBDA,MAL,ILAGRID,TRUNCATE_B,PREC_B,PREC_EXTR,PREC_BINT,PREC_DQDV,PREC_INTV,&
       & NEFIELD,EFIELD,NCMUL,CMUL,NVMAG,VMAG,&
       & NER,ERMIN,ERMAX,ERACC,&
       & NERR,IPERR,RSEED, &
       & CENTERED_ALPHA,SECOND_ORDER_ALPHA
  !Namelist 'model' contains physics parameters
  NAMELIST /model/ ONLY_B0,CALC_DB,NO_PLATEAU,ONLY_DB,INC_EXB,TANG_VM,CLASSICAL,ANISOTROPY,FRICTION,&
       & SCAN_ER,SOLVE_AMB,ER_ROOT,TRIVIAL_AMB,FAST_AMB,SOLVE_QN,TRIVIAL_QN,ZERO_PHI1,ONLY_PHI1,&
       & FACT_CON,D_AND_V,COMPARE_MODELS,&
       & FN,FI,FS,FP,FE,FR,FB,FNE,FTI,FTE,FER,&
       & TANGO,NEOTRANSP,PENTA,TASK3D,TASK3Dlike,NTV,SATAKE,ANA_NTV,JPP,ESCOTO,NEQ2,KN_STELLOPT
  !Namelist 'surfaces' contains the list of flux-surfaces calculated
  NAMELIST /surfaces/ NS,S,SMIN,SMAX,DIRDB,SDKES,DIRS
  !Namelist 'species' contains the list of species calculated
  NAMELIST /species/ NBB,ZB,AB,REGB,FRACB,NBULK,SS_IMP
  !Namelist 'fastions' contains parameters of the fast ions
  NAMELIST /fastions/ GTH,TENDFI,DTFI,EFI,ZFI,LINEART,JMAP,LJMAP,GLOBALFI,&
       & MODELFI,RANDOMSL,FIDELTA,FITRANS,JCORRECTION,JTRANS,PREC_J,PREC_TRANS,PREC_S,JORBIT,FDMAX
  !Namelist 'others' contains other variables
  NAMELIST /others/ FAST_IONS,QN,TRACE_IMP,PLATEAU_OR_PS,USE_B0,&
       & IMP1NU,CLOSEST_LAMBDA,MODEL_ALPHAD,ONE_ALPHA,EXTRA_ALPHA,RE_SOURCE,&
       & CALC_RHS,CALC_DA,CALC_DIFf,CALC_COL,CALC_DG,FLUX_NU,PREC_TOP,NEW_DALPHA,INT_G_NEW!,KROOK_OP
  !Output
  INTEGER NS,ib,NBB,REGB(nbx)
  REAL*8 dt,S(nsx),ZB(nbx),AB(nbx),FRACB(nbx)
  !Other
  CHARACTER*100 serr,line,filename
  INTEGER iline,ioutt,iostat,nefield,ncmul,nvmag,iz,ia,is
  REAL*8 SMAX,SMIN
  !DKES database
  INTEGER, PARAMETER :: ncmuld=28!20
  REAL*8 cmuld(ncmuld) /3E+2,1E+2,3E+1,1E+1,3E+0,1E+0,3E-1,1E-1,3E-2,1E-2,&
       & 3E-3,1E-3,3E-4,1E-4,3E-5,1E-5,3E-6,1E-6,3E-7,1E-7,&
       & 3E-8,1E-8,3E-9,1E-9,3E-10,1E-10,3E-11,1E-11/
!                      & 3E-3,1E-3,3E-4,1E-4,3E-5,1E-5,3E-6,1E-6,3E-7,1E-7/
  INTEGER, PARAMETER :: nefieldd=9
  REAL*8 efieldd(nefieldd) /0E-0,1E-5,3E-5,1E-4,3E-4,1E-3,3E-3,1E-2,3E-2/
  INTEGER, PARAMETER :: nvmagd=1
  REAL*8 vmagd(nvmagd) /0E+0/
  !3D database
  INTEGER, PARAMETER :: ncmulx=18
  REAL*8 cmulx(ncmulx) /3E+2,1E+2,3E+1,1E+1,3E+0,1E+0,3E-1,1E-1,3E-2,1E-2,&
                      & 3E-3,1E-3,3E-4,1E-4,3E-5,1E-5,3E-6,1E-6/
!  INTEGER, PARAMETER :: nefieldx=9
!  REAL*8 efieldx(nefieldx) /0E-0,1E-5,3E-5,1E-4,3E-4,1E-3,3E-3,1E-2,3E-2/
!  INTEGER, PARAMETER :: nvmagx=13
!  REAL*8 vmagx(nvmagx) /-3E-2,-1E-2,-3E-3,-1E-3,-3E-4,-1E-4,0E-0,&
!                      & +1E-4,+3E-4,+1E-3,+3E-3,+1E-2,+3E-2/
!  3D finer database 
!  INTEGER, PARAMETER :: ncmulx=52
!  REAL*8 cmulx(ncmulx) /5E+2,2E+2,1E+2,5E+1,2E+1,1E+1,5E+0,2E+0,1.0E+0,&
!                      & 5E-1,2E-1,1E-1,5E-2,2E-2,1E-2,5E-3,2E-3,1.5E-3,1.0E-3,&
!                      & 9E-4,8E-4,7E-4,6E-4,5E-4,4E-4,3E-4,2E-4,1.5E-4,1.2E-4,1E-4,&
!                      & 9E-5,8E-5,7E-5,6E-5,5E-5,4E-5,3E-5,2E-5,1.5E-5,1.2E-5,1E-5,&
!                      & 9E-6,8E-6,7E-6,6E-6,5E-6,4E-6,3E-6,2E-6,1.5E-6,1.2E-6,1E-6/
  INTEGER, PARAMETER :: nefieldx=50 
  REAL*8 efieldx(nefieldx) /0E-0,1E-5,1.2E-5,1.5E-5,2E-5,3E-5,4E-5,5E-5,6E-5,7E-5,8E-5,9E-5,&
                              &  1E-4,1.2E-4,1.5E-4,2E-4,3E-4,4E-4,5E-4,6E-4,7E-4,8E-4,9E-4,&
                              &  1E-3,1.2E-3,1.5E-3,2E-3,3E-3,4E-3,5E-3,6E-3,7E-3,8E-3,9E-3,&
                              &  1E-2,1.2E-2,1.5E-2,2E-2,3E-2,4E-2,5E-2,6E-2,7E-2,8E-2,9E-2,&
                              &  1E-1,1.2E-1,1.5E-1,2E-1,5E-5/
  INTEGER, PARAMETER :: nvmagx=66
  REAL*8 vmagx(nvmagx)/-9E-2,-8E-2,-7E-2,-8E-2,-5E-2,-4E-2,-3E-2,-2E-2,-1.5E-2,-1.2E-2,-1E-2,&
                     & -9E-3,-8E-3,-7E-3,-8E-3,-5E-3,-4E-3,-3E-3,-2E-3,-1.5E-3,-1.2E-3,-1E-3,&
                     & -9E-4,-8E-4,-7E-4,-8E-4,-5E-4,-4E-4,-3E-4,-2E-4,-1.5E-4,-1.2E-4,-1E-4,&
                     & +1E-4,+1.2E-4,+1.5E-4,+2E-4,+3E-4,+4E-4,+5E-4,+6E-4,+7E-4,+8E-4,+9E-4,&
                     & +1E-3,+1.2E-3,+1.5E-3,+2E-3,+3E-3,+4E-3,+5E-3,+6E-3,+7E-3,+8E-3,+9E-3,&
                     & +1E-2,+1.2E-2,+1.5E-2,+2E-2,+3E-2,+4E-2,+5E-2,+6E-2,+7E-2,+8E-2,+9E-2/
  !Database for NEOTRANSP
  INTEGER, PARAMETER :: ncmuln=17
  REAL*8 cmuln(ncmuln) /1E+1,3E+0,1E+0,3E-1,1E-1,3E-2,1E-2,&
                      & 3E-3,1E-3,3E-4,1E-4,3E-5,1E-5,3E-6,1E-6,3E-7,1E-7/
  INTEGER, PARAMETER :: nefieldn=27
  REAL*8 efieldn(nefieldn) /0.0,3E-7,1E-6,3E-6,1E-5,3E-5,1E-4,3E-4,1E-3,2E-3,5E-3,&
       & 1E-2,2E-2,3E-2,5E-2,1E-1,2E-1,3E-1,5E-1,7E-1,8E-1,1.0,1.2,1.5,2.0,3.0,5.0/
  INTEGER, PARAMETER :: nvmagn=1
  REAL*8 vmagn(nvmagn) /0.0/
!  INTEGER, PARAMETER :: nvmagn=66
!  REAL*8 vmagn(nvmagn)/-9E-2,-8E-2,-7E-2,-8E-2,-5E-2,-4E-2,-3E-2,-2E-2,-1.5E-2,-1.2E-2,-1E-2,&
!                     & -9E-3,-8E-3,-7E-3,-8E-3,-5E-3,-4E-3,-3E-3,-2E-3,-1.5E-3,-1.2E-3,-1E-3,&
!                     & -9E-4,-8E-4,-7E-4,-8E-4,-5E-4,-4E-4,-3E-4,-2E-4,-1.5E-4,-1.2E-4,-1E-4,&
!                     & +1E-4,+1.2E-4,+1.5E-4,+2E-4,+3E-4,+4E-4,+5E-4,+6E-4,+7E-4,+8E-4,+9E-4,&
!                     & +1E-3,+1.2E-3,+1.5E-3,+2E-3,+3E-3,+4E-3,+5E-3,+6E-3,+7E-3,+8E-3,+9E-3,&
!                     & +1E-2,+1.2E-2,+1.5E-2,+2E-2,+3E-2,+4E-2,+5E-2,+6E-2,+7E-2,+8E-2,+9E-2/
!  !Database for PENTA
  INTEGER, PARAMETER :: ncmulp=16
!  REAL*8 cmulp(ncmulp) /1E+1,3E+0,1E+0,3E-1,1E-1,3E-2,1E-2,&
!                      & 3E-3,1E-3,3E-4,1E-4,3E-5,1E-5,3E-6,1E-6,3E-7/
  INTEGER, PARAMETER :: nefieldp=27
!  REAL*8 efieldp(nefieldp) /0.0,3E-7,1E-6,3E-6,1E-5,3E-5,1E-4,3E-4,1E-3,2E-3,5E-3,&
!       & 1E-2,2E-2,3E-2,5E-2,1E-1,2E-1,3E-1,5E-1,7E-1,8E-1,1.0,1.2,1.5,2.0,3.0,5.0/
  INTEGER, PARAMETER :: nvmagp=1
!  REAL*8 vmagp(nvmagp) /0.0/
!  INTEGER, PARAMETER :: nvmagp=66
!  REAL*8 vmagp(nvmagp)/-9E-2,-8E-2,-7E-2,-8E-2,-5E-2,-4E-2,-3E-2,-2E-2,-1.5E-2,-1.2E-2,-1E-2,&
!                     & -9E-3,-8E-3,-7E-3,-8E-3,-5E-3,-4E-3,-3E-3,-2E-3,-1.5E-3,-1.2E-3,-1E-3,&
!                     & -9E-4,-8E-4,-7E-4,-8E-4,-5E-4,-4E-4,-3E-4,-2E-4,-1.5E-4,-1.2E-4,-1E-4,&
!                     & +1E-4,+1.2E-4,+1.5E-4,+2E-4,+3E-4,+4E-4,+5E-4,+6E-4,+7E-4,+8E-4,+9E-4,&
!                     & +1E-3,+1.2E-3,+1.5E-3,+2E-3,+3E-3,+4E-3,+5E-3,+6E-3,+7E-3,+8E-3,+9E-3,&
!                     & +1E-2,+1.2E-2,+1.5E-2,+2E-2,+3E-2,+4E-2,+5E-2,+6E-2,+7E-2,+8E-2,+9E-2/
  REAL*8   cmul(MAX(  ncmulx,  ncmuln,  ncmulp, ncmuld))
  REAL*8 efield(MAX(nefieldx,nefieldn,nefieldp, nefieldd))
  REAL*8   vmag(MAX(  nvmagx,  nvmagn,  nvmagp, nvmagd))
  REAL*8 dummy

  FAST_IONS= .FALSE.

  !-------------------------------------------------------------------------------------------
  !Set values for namelist 'model'
  !-------------------------------------------------------------------------------------------

  !Default values
  !Ambipolarity and quasineutrality
  FAST_AMB=   .FALSE.  
  SCAN_ER=    .FALSE.
  SOLVE_AMB=  .TRUE.
  ER_ROOT  = 0
  TRIVIAL_AMB=.FALSE.
  SOLVE_QN=   .TRUE.  
  TRIVIAL_QN= .FALSE. 
  ZERO_PHI1 = .FALSE.
  ONLY_PHI1 = .FALSE.
  D_AND_V   =5
  COMPARE_MODELS=.FALSE.
  !Details of the drift-kinetic equation
  INC_EXB=   .FALSE. 
  TANG_VM=   .TRUE.  
  CLASSICAL= .TRUE.  
  ANISOTROPY=.TRUE. 
  FRICTION=  .TRUE. 
  FACT_CON=  -3      
  !Equations to be solved
  DIRDB  ="./dummy/"          
  DIRS="./"
  SDKES=-1
  CALC_DB=.FALSE.
  NO_PLATEAU=.FALSE.
  ONLY_DB=.FALSE.   
  !Scan in parameters
  FN=1.0
  FI=1.0
  FS=1.0
  FP=1.0
  FE=1.0
  FR=1.0
  FB=1.0
  FNE=1.0
  FTE=1.0
  FTI=1.0
  FDLNE=1.0
  FDLTE=1.0
  FDLTI=1.0
  FER=1.0
  !Particular problems
  TANGO    = .FALSE.
  NEOTRANSP= .FALSE.
  TASK3D=    .FALSE. 
  TASK3Dlike=.FALSE. 
  PENTA=     .FALSE. 
  NTV=       .FALSE.
  SATAKE=    .FALSE. 
  ANA_NTV=   .FALSE. 
  JPP=       .FALSE.
  ESCOTO=    .FALSE.
  NEQ2=      .FALSE.

  !Read namelist 'model'
  OPEN(unit=1,file="input.model",form='formatted',action='read',iostat=iostat) 
  IF(iostat.NE.0) OPEN(unit=1,file="../input.model",form='formatted',action='read',iostat=iostat) 
  IF (iostat.EQ.0) THEN
!     IF(myrank.EQ.0) WRITE(ioutt,*) 'File "input.model" found'
     READ (1,nml=model)
     CLOSE(1)
  END IF
  
  KNOSOS_STELLOPT=KN_STELLOPT(1).OR.KN_STELLOPT(2).OR.KN_STELLOPT(3).OR.KN_STELLOPT(4).OR.KN_STELLOPT(5).OR.&
                & KN_STELLOPT(6).OR.KN_STELLOPT(7).OR.KN_STELLOPT(8).OR.KN_STELLOPT(9).OR.KN_STELLOPT(10)
  ioutt=iout

  IF(KNOSOS_STELLOPT.AND.LEN(TRIM(KN_EXT)).NE.0) THEN
     filename=TRIM(KN_EXT)
     OPEN(unit=iout,file=filename,form='formatted',action='write',iostat=iostat,&
          access='append')
  ELSE
     IF(numprocs.EQ.1) THEN
        filename="STDOUT"
        OPEN(unit=iout,file=filename,form='formatted',action='write',iostat=iostat)!,&
           !  access='append',status='old')
        filename="STDERR"
        OPEN(unit=1100+myrank,file=filename,form='formatted',action='write',iostat=iostat)
     ELSE
        WRITE(filename,'("STDOUT.",I2.2)') myrank
        OPEN(unit=iout,file=filename,form='formatted',action='write',iostat=iostat)
        WRITE(filename,'("STDERR.",I2.2)') myrank
        OPEN(unit=1100+myrank,file=filename,form='formatted',action='write',iostat=iostat)
     END IF
  END IF
!!$  IF(KNOSOS_STELLOPT.AND.LEN(TRIM(KN_EXT)).NE.0) THEN
!!$     filename=TRIM(KN_EXT)
!!$     OPEN(unit=ioutt,file=filename,form='formatted',action='write',iostat=iostat)
!!$  ELSE
!!$     OPEN(unit=ioutt,file="STDOUT",form='formatted',action='write',iostat=iostat)
!!$  END IF

  !-------------------------------------------------------------------------------------------
  !Set values for namelist 'fastions'
  !-------------------------------------------------------------------------------------------
  FAST_IONS=.FALSE.
  JMAP=.FALSE.
  GTH=0.2
  EFI=5.0E4 !(eV)
  ZFI=1.0
  AFI=1.0
  TENDFI=1.0E-1 !(s)
  DTFI=  1.0E-5 !(s)
  LINEART=.FALSE.
  MODELFI=.FALSE.
  GLOBALFI=.FALSE.
  RANDOMSL=.TRUE.
  JCORRECTION=.TRUE.
  JTRANS=.TRUE.
  FITRANS=.TRUE.
  PREC_J=2E-4
  IF(KNOSOS_STELLOPT) PREC_J=1E-4
  PREC_TRANS=1E-1
  PREC_S=1E-4
  FIDELTA=0.5
  IF(KNOSOS_STELLOPT) FIDELTA=1.0
  FDMAX=8
  JORBIT=-1
  LJMAP=-0.38   
  !Read namelist 'fastions'
  OPEN(unit=1,file="input.fastions",form='formatted',action='read',iostat=iostat) 
  IF(iostat.NE.0) OPEN(unit=1,file="../input.fastions",form='formatted',action='read',iostat=iostat) 
  IF (iostat.EQ.0) THEN
     FAST_IONS=.TRUE.
!     IF(myrank.EQ.0) WRITE(ioutt,*) 'File "input.fastions" found'
     READ (1,nml=fastions)
     CLOSE(1)
  END IF
  TWOEFI=2*EFI
  TWOEFIoZ=2*EFI/ZFI
  VDoV=1.389165e4*SQRT(EFI)*AFI*m_e/ZFI
  IF(.NOT.FAST_IONS.AND.(KN_STELLOPT(4).OR.KN_STELLOPT(5))) MODELFI=.TRUE.
  
  !-------------------------------------------------------------------------------------------
  !Set values for namelist 'parameters'
  !-------------------------------------------------------------------------------------------

  !To be discontinued
  IMP1NU        =.FALSE.
  CENTERED_ALPHA=.FALSE.
  SECOND_ORDER_ALPHA=.FALSE.
  CLOSEST_LAMBDA=.TRUE.
  MODEL_ALPHAD=.FALSE.
  INT_G_NEW=.FALSE.
  !Default values
  !Debug
  GEN_FLAG=.FALSE. 
  DEBUG=   .FALSE.  
  TIME=    .TRUE.   
  I0=      -1
  L0=      -1        
  PLOTG=   .FALSE.  
  !How to describe magnetic field structure
  USE_B1=    .FALSE.  
  USE_B0pB1 =.FALSE.
  QS_B0     =.FALSE.
  QS_B0_1HEL=.FALSE.
  USE_SHAING=.FALSE.
  SHAING_1NU=.FALSE.
  SHAING_SQRTNU=.FALSE.
  SHAING_SBP=.FALSE.
  !Determine algorithms to be used
  REMOVE_DIV=    .FALSE.  
  DELTA=         .TRUE.   
  !Set numerical resolution
  NTURN=1
  ILAGRID=.FALSE.
  IF(JMAP) THEN
     MLAMBDA= 32
     MAL    = 128
  ELSE IF(TANGO) THEN
     MLAMBDA=128
     MAL    =128
  ELSE
     MLAMBDA= 64
     MAL    = 64
  END IF
  TRUNCATE_B= -200     
  PREC_B=     1E-5     
  PREC_EXTR=  1E-6   
  IF(KNOSOS_STELLOPT) PREC_EXTR=  1E-7
  PREC_BINT=  1E-2     
  PREC_DQDV=  5E-2     
  PREC_INTV=  1E-2     
  !Set details of the monoenergetic database
  NCMUL  =ncmuld       
  NEFIELD=nefieldd     
  NVMAG  =nvmagd
  CMUL=0
  EFIELD=0
  VMAG=0
  CMUL(1:ncmuld)    =cmuld
  EFIELD(1:nefieldd)=efieldd
  VMAG(1:nvmagd)    =vmagd        
  !Set details of ambipolarity
  NER  = 21     
  ERMIN=+20.0  
  ERMAX=-20.0  
  ERACC= 0.1    
  !Estimate error bars
  NERR =1       
  IPERR=-1.0   
  RSEED=.FALSE.!TRUE.
  
  !Read namelist 'parameters'
  OPEN(unit=1,file="input.parameters",form='formatted',action='read',iostat=iostat) 
  IF(iostat.NE.0) OPEN(unit=1,file="../input.parameters",form='formatted',action='read',iostat=iostat) 
  IF(iostat.EQ.0) THEN
!     IF(myrank.EQ.0) WRITE(ioutt,*) 'File "input.parameters" found'
     READ (1,nml=parameters)
     CLOSE(1)
  END IF
  !-------------------------------------------------------------------------------------------
  !Set default values for namelists 'surfaces' and 'species'
  !-------------------------------------------------------------------------------------------

  ns=1
  s(1)=1
  s(2:nsx)=0
  DIRS=' '
  SMIN=0
  SMAX=1
  !Default values: low collisionality hydrogen + adiabatic electrons, no impurities
  nbb    = 2
  nbulk  = 2
  Zb(1)  =-1.            
  Ab(1)  = 5.4858E-4
  regb(1)=+3           
  Zb(2)  =+1.            
  Ab(2)  = 1.0073
  regb(2)=+3
  Zb(3:nbx)=0
  Ab(3:nbx)=0
  fracb(1:2)=1
  fracb(3:nbx)=0
  SS_IMP=.FALSE.

  IF(TASK3D) THEN
     !Could be overwritten by file 'input.surfaces'
     ns=40
     DO is=1,ns
        s(is)=(is-0.5)*(is-0.5)/REAL(ns*ns)
     END DO
     !Could be overwritten by file 'input.species'
     REGB(1:2)=0
     ZB(1)=-1.
     AB(1)=5.48579909E-4
     OPEN(unit=1,file="input-prof.txt",action='read',iostat=iostat)
     IF(iostat.NE.0) OPEN(unit=1,file="../input-prof.txt",action='read',iostat=iostat)     
     IF(iostat.EQ.0) THEN
!        IF(myrank.EQ.0) WRITE(ioutt,*) 'File "input-prof.txt" found'
        DO iline=1,7
           READ(1,*) line
        END DO
        READ(1,*) dummy,dummy,dummy,iz,ia
        CLOSE(1)
        Zb(2)=REAL(iz)
        Ab(2)=ia*1.00727647
     ELSE
        serr="Wrong data in input-prof.txt"
        CALL END_ALL(serr,.FALSE.)
     END IF
     NBB=6
     NBULK=NBB
     Zb(3)  =+3.         
     Ab(3)  =6.941
     Zb(4)  =+19.         
     Ab(4)  =47.867
     Zb(5)  =+20.            
     Ab(5)  =50.9415
     Zb(6)  =+23.            
     Ab(6)  =55.847
     regb(3:nbb)=10
     fracb(1:2)=1
     fracb(3:nbb)=0.0001
     D_AND_V=5
  ELSE IF(NEOTRANSP) THEN
     ns=7
     DO is=1,ns
        s(is)=(is-0.5)*(is-0.5)/REAL(ns*ns)
     END DO
     s(1)=s(1)*3
 ! ELSE IF(PENTA) THEN
 !    ns=7
 !    DO is=1,ns
 !       s(is)=(is-0.5)*(is-0.5)/REAL(ns*ns)
 !    END DO
     INT_G_NEW=.TRUE.
  ELSE IF(TANGO) THEN
     REGB(1)=-1
     REGB(2)=-1
     INT_G_NEW=.TRUE.
  END IF

  !-------------------------------------------------------------------------------------------
  !Set values for namelist 'surfaces'
  !-------------------------------------------------------------------------------------------

  !Read namelist 'surfaces'
  OPEN(unit=1,file="input.surfaces",form='formatted',action='read',iostat=iostat) 
  IF(iostat.NE.0) OPEN(unit=1,file="../input.surfaces",form='formatted',action='read',iostat=iostat) 
  IF (iostat.EQ.0) THEN
!     IF(myrank.EQ.0) WRITE(ioutt,*) 'File "input.surfaces" found'
     READ (1,nml=surfaces)
     CLOSE(1)
  END IF
  IF(ns.GT.1.AND.ABS(s(1)-1).LT.ALMOST_ZERO) THEN
     DO is=1,ns
        !        s(is)=SMIN+(SMAX-SMIN)*(is-0.5)*(is-0.5)/REAL(ns*ns)
        s(is)=SMIN+(SMAX-SMIN)*(is-0.5)/REAL(ns)
     END DO
  END IF

  !-------------------------------------------------------------------------------------------
  !Set values for namelist 'species'
  !-------------------------------------------------------------------------------------------

  !Read namelist 'species'
  OPEN(unit=1,file="input.species",form='formatted',action='read',iostat=iostat) 
  IF(iostat.NE.0) OPEN(unit=1,file="../input.species",form='formatted',action='read',iostat=iostat) 
  IF (iostat.EQ.0) THEN
!     IF(myrank.EQ.0) WRITE(ioutt,*) 'File "input.species" found'
     READ (1,nml=species)
     CLOSE(1)
  END IF


  !-------------------------------------------------------------------------------------------
  !Set values for namelist 'transport'
  !-------------------------------------------------------------------------------------------

  dt=-1.0
  !Read namelist 'transport'
  OPEN(unit=1,file="input.transport",form='formatted',action='read',iostat=iostat) 
  IF(iostat.NE.0) OPEN(unit=1,file="../input.transport",form='formatted',action='read',iostat=iostat) 
  IF (iostat.EQ.0) THEN
     READ (1,nml=transport)
     CLOSE(1)
  END IF

  !-------------------------------------------------------------------------------------------

  !These variables are determined at some point of the run
  CALCULATED_INT=.FALSE.
  PHI1_READ=     .FALSE.
  DKES_READ=     .FALSE.
  CONVERGED=     .FALSE.
  TRACE_IMP=     .FALSE.
  PLATEAU_OR_PS= .FALSE.
  DR_READ=       .FALSE.
  
  !Some of the choices above may be incompatible, have to be changed accordingly.
  DO ib=1,nbb !in case you don't remember the electron mass, just write a small number
     IF(Ab(ib).LT.1) Ab(ib)=5.48579909E-4 
     IF((ib.EQ.1.AND.Zb(ib).GT.0).OR.(ib.GT.1.AND.Zb(ib).LT.0)) THEN
        serr="Wrong data in input.species"
        CALL END_ALL(serr,.FALSE.)
     END IF
  END DO

  IF(FAST_AMB) THEN
     CALC_DB=.TRUE.
     NCMUL  =ncmulx  
     NEFIELD=nefieldx     
     NVMAG  =nvmagx
     CMUL(1:ncmulx)    =cmulx       
     EFIELD(1:nefieldx)=efieldx     
     VMAG(1:nvmagx)    =vmagx        
  END IF
  IF(TASK3D) THEN
     !Could be overwritten by 'input.model'
     SOLVE_QN  =.FALSE.
     TRIVIAL_QN=.TRUE.
  ELSE IF(TANGO) THEN
     SOLVE_QN=.FALSE.
     TANG_VM=.TRUE.
     NER=41
  ELSE IF(NEOTRANSP.OR.PENTA.OR.KNOSOS_STELLOPT) THEN
     CALC_DB=.TRUE.
     IF(NEOTRANSP.OR.PENTA) TANG_VM=.FALSE.
     IF(NEOTRANSP.OR.PENTA.OR.KN_STELLOPT(4)) THEN !.OR..NOT.(KN_STELLOPT(5))) THEN
        SOLVE_AMB=.FALSE.
        SOLVE_QN=.FALSE.   
        ONLY_DB=.TRUE.
     END IF
     IF(NEOTRANSP.OR.PENTA) TANG_VM=.FALSE.
     IF(NEOTRANSP) THEN
        ncmul  =ncmuln
        nefield=nefieldn
        nvmag  =nvmagn
        cmul(1:ncmuln)    =cmuln
        efield(1:nefieldn)=efieldn
        vmag(1:nvmagn)    =vmagn
     ELSE IF(PENTA) THEN
        ncmul=1
        nefield=1
        nvmag=1
!        ncmul  =ncmulp
!        nefield=nefieldp
!        nvmag  =nvmagp
!        cmul(1:ncmulp)    =cmulp
!        efield(1:nefieldp)=efieldp
!        vmag(1:nvmagp)    =vmagp
     ELSE IF(KNOSOS_STELLOPT) THEN
        ncmul  =1
        nefield=3
        nvmag  =1
        cmul(1)=1E-5
        efield(1)=3E-4
        efield(2)=0E+0
        efield(3)=0E+0
        vmag=0
     END IF
  END IF
  IF(.NOT.TANG_VM) THEN
     FS=0
     nvmag=1
     vmag=0.0
  END IF

  MODEL_ALPHAD=.FALSE.
  CLOSEST_LAMBDA=.TRUE.
  EXTRA_ALPHA=.FALSE.
  ONE_ALPHA=.TRUE.
!  KROOK_OP=.FALSE.
  FLUX_NU=.FALSE.
  CALC_RHS=.FALSE.
  CALC_DA=.FALSE.
  CALC_DIFF=.FALSE.
  CALC_COL=.FALSE.
  PREC_TOP=.FALSE.
  NEW_DALPHA=.TRUE.
  RE_SOURCE=.TRUE.
  IF(GEN_FLAG(3)) THEN
     MODEL_ALPHAD=.TRUE.
     NEW_DALPHA=.FALSE.
     RE_SOURCE=.FALSE.
  END IF
  IF(GEN_FLAG(4)) INT_G_NEW=.TRUE.
!  IF(TANGO) RE_SOURCE=.FALSE.

  IF(.NOT.MODEL_ALPHAD) CLOSEST_LAMBDA=.TRUE.
  CALC_DG=CALC_RHS.OR.CALC_DA.OR.CALC_COL.OR.CALC_DIFF.OR.FLUX_NU
  
  IF(KNOSOS_STELLOPT) ONE_ALPHA=.TRUE.
  USE_B0=USE_B1.OR.USE_B0pB1
  IF(QS_B0_1HEL) QS_B0=.TRUE.
  DO ib=3,nbb
     IF(REGB(ib).GE.10) THEN
        TRACE_IMP    =.TRUE.
        IF(REGB(ib).LT.13.OR.REGB(ib).EQ.22) PLATEAU_OR_PS=.TRUE.
     END IF
  END DO
  IF(TRACE_IMP) THEN
     IF(REGB(2).LE.-1) REGB(2)=+3
!     IF(.NOT.SOLVE_QN) TRIVIAL_QN=.TRUE.
     IF(ZERO_PHI1) TRIVIAL_QN=.TRUE.
  ELSE
     ANISOTROPY=.FALSE.
  END IF
!  IF(ZERO_PHI1) TRIVIAL_QN=.TRUE.
  IF(TRIVIAL_QN)  SOLVE_QN =.FALSE.
  IF(TRIVIAL_AMB) SOLVE_AMB=.FALSE.
  QN=SOLVE_QN.OR.TRIVIAL_QN
  IF(SOLVE_QN) THEN
     DELTA=.TRUE.
!     CLOSEST_LAMBDA=.FALSE.
  END IF
  IF(.NOT.SOLVE_QN) COMPARE_MODELS=.FALSE.
  IF(SATAKE) TRUNCATE_B=20
  IF(PLOT_XYZ) THEN
     PREC_B=-1E-9
     TRUNCATE_B=-10
  END IF
  IF(ALLOCATED(absnablar)) DEALLOCATE(absnablar)
  IF(MAL.GT.0) ALLOCATE(absnablar(MAL,MAL))
  
  IF(TASK3Dlike) TASK3D=.TRUE.

  ncmult  =ncmul
  nefieldt=nefield
  nvmagt  =nvmag
  IF(ALLOCATED(cmult)) THEN
     DEALLOCATE( cmult, efieldt)
     DEALLOCATE(lcmult,lefieldt)
     DEALLOCATE(lD11dkes1,lD11dkes2,lD11tab,D31dkes)
     IF(ALLOCATED(vmagt)) DEALLOCATE(vmagt,lvmagt)
  END IF
  ALLOCATE( cmult(ncmult), efieldt(nefieldt), vmagt(nvmagt))
  ALLOCATE(lcmult(ncmult),lefieldt(nefieldt),lvmagt(nvmagt))
  ALLOCATE(lD11dkes1(ncmult,nefieldt),lD11dkes2(ncmult,nefieldt),lD11tab(ncmult,nefieldt,nvmagt))
  ALLOCATE(D31dkes(ncmult,nefield))
  cmult  =cmul
  efieldt=efield
  vmagt  =vmag
  lcmult  =LOG(cmult)   
  lefieldt=LOG(efieldt)
  lvmagt  =LOG(ABS(lvmagt))
  lefieldt(1)=-1000 !dummy value to avoid log(0.0)
  
  !Write input parameters
  IF(myrank.EQ.0) THEN
     WRITE(ioutt,nml=model)
     WRITE(ioutt,nml=parameters)
     WRITE(ioutt,nml=surfaces)
     WRITE(ioutt,nml=species)
     WRITE(ioutt,nml=fastions)
     WRITE(ioutt,nml=others)
  END IF


END SUBROUTINE READ_INPUT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE INIT_FILES()

!-------------------------------------------------------------------------------------------------
!Initialize output files
!-------------------------------------------------------------------------------------------
  USE GLOBAL
  USE KNOSOS_STELLOPT_MOD
  IMPLICIT NONE
  !Others
  CHARACTER*100 filename
  INTEGER iostat
  
  IF(DEBUG) THEN
     OPEN(unit=1400+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=1500+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=2000+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=2100+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=2900+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=3000+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=3100+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=3200+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=3300+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=3400+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=4400+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=5000+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=5100+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=5200+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=5300+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=5400+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=5500+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=5600+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=7000+myrank,form='formatted',action='write',iostat=iostat)
  END IF

  IF(KNOSOS_STELLOPT.AND.LEN(TRIM(KN_EXT)).NE.0) THEN
     !     filename='kn_log.'//TRIM(KN_EXT)
     filename=TRIM(KN_EXT)
     OPEN(unit=iout,file=filename,form='formatted',action='write',iostat=iostat,&
          access='append',status='old')
  ELSE
     IF(numprocs.EQ.1) THEN
        filename="STDOUT"
        OPEN(unit=iout,file=filename,form='formatted',action='write',iostat=iostat,&
             access='append',status='old')
        filename="STDERR"
        OPEN(unit=1100+myrank,file=filename,form='formatted',action='write',iostat=iostat)
     ELSE
        WRITE(filename,'("STDOUT.",I2.2)') myrank
        OPEN(unit=iout,file=filename,form='formatted',action='write',iostat=iostat)
        WRITE(filename,'("STDERR.",I2.2)') myrank
        OPEN(unit=1100+myrank,file=filename,form='formatted',action='write',iostat=iostat)
     END IF
  END IF
  IF(.NOT.(KNOSOS_STELLOPT)) THEN
     IF(numprocs.EQ.1) filename="B.map"
     IF(numprocs.GT.1) WRITE(filename,'("B.map.",I2.2)') myrank
     OPEN(unit=1200+myrank,file=filename,form='formatted',action='write',iostat=iostat)
     WRITE(1200+myrank,&
          & '("s \zeta_{Boozer}  \theta_{Boozer}(right-handed)  B[T]  (v_B.\nabla\psi)[A.U.]")')
     
     IF(numprocs.EQ.1) filename="B0.map"
     IF(numprocs.GT.1) WRITE(filename,'("B0.map.",I2.2)') myrank
     OPEN(unit=1300+myrank,file=filename,form='formatted',action='write',iostat=iostat)
     WRITE(1300+myrank,&
          & '("s \zeta_{Boozer}  \theta_{Boozer}(right-handed)  B[T]  (v_B.\nabla\psi)[A.U.]")')
  END IF

!  IF(numprocs.EQ.1) filename="results.knosos"
!  IF(numprocs.GT.1) WRITE(filename,'("results.knosos.",I2.2)') myrank
!  OPEN(unit=200+myrank,file=filename,form='formatted',action='write',iostat=iostat)
!  WRITE(200+myrank,'("cmul efield weov wtov L11m L11p L31m L31p L33m L33p scal11&
!       & scal13 scal33 max\_residual chip psip btheta bzeta vp vmag")')

  IF(SOLVE_AMB.OR.SOLVE_QN) THEN
     IF(numprocs.EQ.1) filename="flux.amb"
     IF(numprocs.GT.1) WRITE(filename,'("flux.amb.",I2.2)') myrank
     OPEN(unit=300+myrank,file=filename,form='formatted',action='write',iostat=iostat)
     WRITE(300+myrank,'("s E_r[V/m] (\Gamma_b/n_b[m/s]  Q_b/n_b/T_b[m/s] L_1^b[m^2/s] L_2^b[m^2/s]&
          &  Z_b n_b[10^{19}m^{-3}] dlnn_b/dr[m^{-1}] T_b[eV] dlnT_b/dr[m^{-1}], b=1,NBB),&
          & size(e\varphi_1/T_i) ")')

     IF(TASK3D) THEN
        IF(numprocs.EQ.1) filename="knososTASK3D.flux"
        IF(numprocs.GT.1) WRITE(filename,'("knososTASK3D.flux.",I2.2)') myrank
        OPEN(unit=5600+myrank,file=filename,form='formatted',action='write',iostat=iostat)
        WRITE(5600+myrank,'("rho   E_r[kV/m]   Gamma_e[m^-2s^-1]   Gamma_i1[m^-2s^-1]  &
             & Gamma_i2[m^-2s^-1]   Gamma_i[m^-2s^-1]  &
             & Q_e[W/m^2]   Q_i1[W/m^2]   Q_i2[W/m^2]   Q_i[W/m^2]")')
     ELSE IF(TANGO) THEN
        WRITE(filename,'("knosos_for_tango.txt.",I2.2)') myrank
        OPEN(unit=5600+myrank,file=filename,form='formatted',action='write',iostat=iostat)
        WRITE(5600+myrank,'("x/a=sqrt(s)   E_r[V/m]   &
             & <Q_i\dot\nabla x>[W/m^2] <Q_e\dot\nabla x>[W/m^2]")')
     END IF

     IF(COMPARE_MODELS) THEN
        IF(numprocs.EQ.1) filename="flux.amb.comp"
        IF(numprocs.GT.1) WRITE(filename,'("flux.amb.comp.",I2.2)') myrank
        OPEN(unit=4300+myrank,file=filename,form='formatted',action='write',iostat=iostat)
        WRITE(4300+myrank,'("s E_r[V/m] (\Gamma_b/n_b[m/s]  Q_b/n_b/T_b[m/s] &
          &  L_1^b[m^2/s] L_2^b[m^2/s] &
          &  n_b[10^{19}m^{-3}] dlnn_b/dr T_b[eV] dlnT_b/dr Z_b, b=1,NBB), size(e\varphi_1/T_i) ")')
     END IF
     
     IF(SOLVE_QN) THEN
        IF(numprocs.EQ.1) filename="flux.modes"
        IF(numprocs.GT.1) WRITE(filename,'("flux.modes.",I2.2)') myrank
        OPEN(unit=700+myrank,file=filename,form='formatted',action='write',iostat=iostat)
        WRITE(700+myrank,'("s E_r[V/m] cosine(0)/sine(1) n  m &
	& (\Gamma_b/n_b[m/s]  Q_b/n_b/T_b[m/s] L_1^b[m^2/s] L_2^b[m^2/s], b=1,NBB)")')
     END IF
  END IF
  
  IF(QN) THEN
     IF(numprocs.EQ.1) filename="varphi1.modes"
     IF(numprocs.GT.1) WRITE(filename,'("varphi1.modes.",I2.2)') myrank
     OPEN(unit=400+myrank,file=filename,form='formatted',action='write',iostat=iostat)
     WRITE(400+myrank,&
     &'("s  cosine(0)/sine(1) n  m  \varphi_{nm} (Boozer angles are right-handed)")')

     IF(numprocs.EQ.1) filename="varphi1.map"
     IF(numprocs.GT.1) WRITE(filename,'("varphi1.map.",I2.2)') myrank
     OPEN(unit=500+myrank,file=filename,form='formatted',action='write',iostat=iostat)
     WRITE(500+myrank,'("s  \zeta_{Boozer}  \theta_{Boozer}(right-handed) &
          & \varphi_1[V]  e\varphi_1/T_i  (v_E.\nabla\psi)[A.U.]  n1e/n0e  n1i/n0i ")')

     IF(COMPARE_MODELS) THEN

        IF(numprocs.EQ.1) filename="varphi1.modes.comp"
        IF(numprocs.GT.1) WRITE(filename,'("varphi1.modes.comp.",I2.2)') myrank
        OPEN(unit=4100+myrank,file=filename,form='formatted',action='write',iostat=iostat)
        WRITE(4100+myrank,&
        &'("s cosine(0)/sine(1) n  m  \varphi_{nm} (Boozer angles are right-handed)")')
        
        IF(numprocs.EQ.1) filename="varphi1.map.comp"
        IF(numprocs.GT.1) WRITE(filename,'("varphi1.map.comp.",I2.2)') myrank
        OPEN(unit=4200+myrank,file=filename,form='formatted',action='write',iostat=iostat)
        WRITE(4200+myrank,'("s \zeta_{Boozer}  \theta_{Boozer}(right-handed) &
             & \varphi_1[V]  e\varphi_1/T_i  (v_E.\nabla\psi)[A.U.]  n1e/n0e  n1i/n0i ")')
        
     END IF

   END IF

   IF(KNOSOS_STELLOPT) THEN
      IF(numprocs.EQ.1) filename="stellopt.knosos"
      IF(numprocs.GT.1) WRITE(filename,'("stellopt.knosos.",I2.2)') myrank
      OPEN(unit=600+myrank,file=filename,form='formatted',action='write',iostat=iostat)
      WRITE(600+myrank,'("s 1NU SNU SBP GMC GMA QER VB0 VBB WBW DBO VBM FTR")')
      IF(numprocs.EQ.1) filename="gammacs.map"
      IF(numprocs.GT.1) WRITE(filename,'("gammacs.map.",I2.2)') myrank
      OPEN(unit=6100+myrank,file=filename,form='formatted',action='write',iostat=iostat)
      WRITE(6100+myrank,'("s alpha lambda[1/T] gamma_C^* vd_s[A.U.] vd_alpha[A.U.] J[A.U.] bottom_index exit_time[s]")')
      IF(numprocs.EQ.1) filename="prompt.lambda"
      IF(numprocs.GT.1) WRITE(filename,'("prompt.lambda",I2.2)') myrank
      OPEN(unit=6500+myrank,file=filename,form='formatted',action='write',iostat=iostat)
      WRITE(6500+myrank,'("s lambda[1/T] fraction fraction promt_fraction protracted_fraction")')
   ELSE
      IF(numprocs.EQ.1) filename="flux.knosos"
      IF(numprocs.GT.1) WRITE(filename,'("flux.knosos.",I2.2)') myrank
      OPEN(unit=600+myrank,file=filename,form='formatted',action='write',iostat=iostat)
      WRITE(600+myrank,'("s E_r[V/m] (\Gamma_b/n_b[m/s]  Qb/n_b/T_b[m/s] L_1^b[m^2/s] L_2^b[m^2/s]&
           &  n_b[10^{19}m^{-3}] dlnn_b/dr T_b[eV] dlnT_b/dr Z_b, b=1,NBB), size(e\varphi_1/T_i) ")')
   END IF
   
   IF(TASK3D) THEN
      IF(numprocs.EQ.1) filename="knososTASK3D.ambEr"
      IF(numprocs.GT.1) WRITE(filename,'("knososTASK3D.ambEr.",I2.2)') myrank
      OPEN(unit=5300+myrank,file=filename,form='formatted',action='write',iostat=iostat)
      WRITE(5300+myrank,'("rho   n_e[m^-3]   n_i1[m^-3   n_i2[m^-3]   T_e[eV]   T_i1[eV]   T_i2[eV]  &
           & E_r[kV/m]   Gamma[m^-2s^-1]   Q_e[W/m^2]   Q_i[W/m^2]   Chi_e[m^2/s]   Chi_i[m^2/s]")')
   END IF
   
   IF(nerr.GT.1) THEN
!      IF(numprocs.EQ.1) THEN
         filename="flux.av"
         OPEN(unit=800+myrank,file=filename,form='formatted',action='write',iostat=iostat)
         WRITE(800+myrank,'("ohos E_r[V/m] err(E_r)[V/m] (\Gamma_b/n_b[m/s] err(\Ganna_b/n_b)[m/s] &
              & Qb/n_b/T_b[m/s], err(Qb/n_b/T_b)[m/s], b=1,NBB) ")')
!      ELSE
!         WRITE(filename,'("flux.av.",I2.2)') myrank
!         OPEN(unit=800+myrank,file=filename,form='formatted',action='write',iostat=iostat)
!         WRITE(800+myrank,'("s E_r[V/m] err(E_r)[V/m] (\Gamma_b/n_b[m/s] err(\Ganna_b/n_b)[m/s] &
!              & Qb/n_b/T_b[m/s], err(Qb/n_b/T_b)[m/s], b=1,NBB) ")')
!      END IF
   END IF
   
   IF(TRACE_IMP) THEN
      IF(numprocs.EQ.1) filename="imp.knosos"
      IF(numprocs.GT.1) WRITE(filename,'("imp.knosos.",I2.2)') myrank
      OPEN(unit=900+myrank,file=filename,form='formatted',action='write',iostat=iostat)
      WRITE(900+myrank,'("s  Z_z  A_z \Gamma_z/n_z[m/s]  V_z[m/s]  D_z[m^2/s] dlnn_z/dr[m^-1]  &
        & D_{E_r}[m^2/s] eE_r/T_z[m^-1]  D_T[m^2/s] dlnT_z/dr[m^-1]  D_n[m^2/s] dlnn_i/dr[m^-1]&
        & \Gamma_{anisotrp}/n_z[m/s]")')
   END IF
  
END SUBROUTINE INIT_FILES


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

