!Read and proces results from DKES

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CALC_DATABASE(s,is,ns)

!----------------------------------------------------------------------------------------------- 
!Solve the monoenergetic drift kinetic equation at s(is) for several values of cmul and efield,
!either hard-coded or read from namelist 'parameters'
!Generate a database of transport coefficients that can be:
!-compared with DKES;
!-used in a transport simulation.
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  USE KNOSOS_STELLOPT_MOD
  IMPLICIT NONE
  !Input 
  INTEGER is,ns
  REAL*8 s(ns)
  !Others
  CHARACTER*100 filename
  INTEGER nalphab
!  INTEGER nlambda,nal
  INTEGER icmul,iefield,ivmag,idummy(1),iostat,iturn
  REAL*8 ZB(1),AB(1),nb(1),Tb(1),Epsi,dummy(1)
  REAL*8 new_nb(ncmult),cmul0,D11tab(ncmult,nefieldt,nvmagt),D11(Nnmp,Nnmp),D31
  REAL*8 zeta(nax),theta(nax),dn1(nax,nax),phi1c(Nnmp),Mbbnm(Nnmp),trMnm(Nnmp),dn1nm(Nnmp,Nnmp)
  REAL*8 g11_ft,eps_eff,rhostar
  !Time
  CHARACTER*30, PARAMETER :: routine="CALC_DATABASE"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  IF(.NOT.KNOSOS_STELLOPT) THEN
     IF(numprocs.EQ.1) filename="results.knosos"
     IF(numprocs.GT.1) WRITE(filename,'("results.knosos.",I2.2)') myrank
     OPEN(unit=200+myrank,file=filename,form='formatted',action='write',iostat=iostat)
  END IF
  !Read a DKES database of monoenergetic transport coefficients
!  IF(.NOT.(PENTA.OR.NEOTRANSP.OR.PENTA.OR.FAST_AMB.OR.&
  IF(.NOT.(PENTA.OR.NEOTRANSP.OR.PENTA.OR.KNOSOS_STELLOPT)) CALL READ_DKES_TABLE(is)

  CALCULATED_INT=.FALSE.
  !Calculate dummy sources, etc
  ZB=+1.
  AB=+1.
  nb=1E+0
  Tb=5E+3
  cmul_1NU=1E+20  
  cmul_PS= 1E+20
  CALL DKE_CONSTANTS(1,1,ZB,AB,IDUMMY,nb,dummy,Tb,dummy,dummy(1),.FALSE.)
  phi1c=0
  trMnm=0
  Mbbnm=0
  cmul0=nu(1)/v(1)/2.
  !Determine collisionalities that limit different collisionality regimes
  !if not low-collisionality problems
  cmul_PS=1E10
  cmul_1NU=1E10

  IF(KN_STELLOPT(4)) THEN
     CALL CALC_FAST_ION_CONFINEMENT(S,IS,NS,MAL,MLAMBDA)
     IF(.NOT.KN_STELLOPT(1)) RETURN
  END IF
  
  IF(.NOT.JPP.AND..NOT.NO_PLATEAU) THEN
     cmul_PS =-1.0   !this may be necessary if there is 
     cmul_1NU=-1.0   !some connection imposed between regimes
     CALL CALC_PS(1,ALMOST_ZERO,D11(1,1))
     D11onu=D11(1,1)/cmul0
     IF(DKES_READ) THEN
        D11pla=D11pla/fdkes(1)
     ELSE
        CALL CALC_PLATEAU(1,ALMOST_ZERO,D11(1,1),dn1nm)
        D11pla=D11(1,1)
     END IF
     cmul_PS =ABS(D11pla/D11onu)  !CMUL such that plateau and PS transport are equal
     vmconst=0
     !     CALL CALC_BANANA(1,Epsi,D11(1,1))
     IF(ESCOTO) THEN
        CALL CALC_LOW_COLLISIONALITY_NEW(1,ZERO,phi1c,Mbbnm,trMnm,&
             & D11,D31,nalphab,zeta,theta,dn1,dn1nm)
     ELSE
        CALL CALC_LOW_COLLISIONALITY(1,ZERO,phi1c,Mbbnm,trMnm,&
             & D11,nalphab,zeta,theta,dn1,dn1nm)
     END IF
     D11tab(1,1,1)=D11(1,1)
     CALCULATED_INT=.FALSE.
     D11nu=D11tab(1,1,1)*cmul0
     cmul_1NU=ABS(D11nu/D11pla)  !CMUL such that plateau and 1/nu transport are equal 
!     cmul_1NU=ABS(aiota/rad_R*eps32)!Use formula, because between plateau and 1/nu there
     D11pla=D11pla*fdkes(1)           !might exist banana regime

     IF(KN_STELLOPT(2)) THEN
        cmul0=nu(iv0)/v(iv0)/2
        new_nb=nb(1)*cmult/cmul0
        CALL DKE_CONSTANTS(1,1,ZB,AB,IDUMMY,new_nb(1),dummy,Tb,dummy,dummy(1),.FALSE.)
        Epsi=efieldt(1)*v(iv0)/psip
        IF(ESCOTO) THEN
           CALL CALC_LOW_COLLISIONALITY_NEW(1,ZERO,phi1c,Mbbnm,trMnm,&
                & D11,D31,nalphab,zeta,theta,dn1,dn1nm)
        ELSE
           CALL CALC_LOW_COLLISIONALITY(1,ZERO,phi1c,Mbbnm,trMnm,&
                & D11,nalphab,zeta,theta,dn1,dn1nm)
        END IF
        KN_SNU=D11(1,1)
     END IF
  END IF

!  IF(.NOT.USE_B0) THEN
!     bnmc0(1:Nnm)=bnmc(1:Nnm)
!     bnms0(1:Nnm)=bnms(1:Nnm)
!     bnmc1(1:Nnm)=0
!     bnms1(1:Nnm)=0
!  END IF

  IF(NEOTRANSP.OR.KNOSOS_STELLOPT) THEN
     g11_ft=d11nu*fdkes(1)
     IF(KN_STELLOPT(1)) THEN
        eps_eff=5.0*g11_ft*rad_R*rad_R*borbic(0,0)*borbic(0,0)
        KN_1NU=eps_eff
     END IF
     IF(KNOSOS_STELLOPT) RETURN
!        IF(.NOT.KNOSOS_STELLOPT) THEN
!           IF(numprocs.EQ.1) filename="mono.opt"
!           IF(numprocs.GT.1) WRITE(filename,'("mono.opt.",I2.2)') myrank
!           OPEN(unit=6000+myrank,file=filename,form='formatted',action='write',iostat=iostat,position='append')
!           WRITE(6000+myrank,'("s Gamma_1nu Gamma_sqrtnu Gamma_sbp")')
!           WRITE(6000+myrank,'(1(1pe13.5))') eps_eff
!           CLOSE(6000+myrank)
        !        END IF
        
  END IF
     
  IF (NEOTRANSP) THEN
     eps_eff=(5.0*g11_ft*rad_R*rad_R*borbic(0,0)*borbic(0,0))**0.66666
     IF(numprocs.EQ.1) filename="knosos.dk"
     IF(numprocs.GT.1) WRITE(filename,'("knosos.dk.",I2.2)') myrank
     OPEN(unit=6000+myrank,file=filename,form='formatted',action='write',iostat=iostat)
     IF(is.EQ.1) THEN
        WRITE(6000+myrank,'("cc")')  
        WRITE(6000+myrank,'("cc")')  
        WRITE(6000+myrank,'("cc")')  
        WRITE(6000+myrank,'("cc")')  
        WRITE(6000+myrank,'("cc")')  
     END IF
     WRITE(6000+myrank,'(5(1pe13.5),"  NaN",1pe13.5,"  r,R,B,io,xkn,ft,<b^2>")') &
          & SQRT(s(is))*rad_a,rad_R,borbic(0,0),ABS(iota),ABS(borbic(0,1))/eps,avb2
     IF(is.EQ.1) WRITE(6000+myrank,'("c         eps_eff     g11_ft      efield_u    g11_er    ex_er")')  
     WRITE(6000+myrank,'("cfit ",2(1pe13.5),"      NaN        NaN      NaN")') eps_eff,g11_ft

!  ELSE IF(PENTA) THEN
!!     IF(numprocs.EQ.1) filename="knosos.penta"
!!     IF(numprocs.GT.1) WRITE(filename,'("knosos.penta.",I2.2)') myrank
!     OPEN(unit=6000+myrank,file=filename,form='formatted',action='write',iostat=iostat)
!     IF(is.EQ.1) THEN
!!        WRITE(6000+myrank,'("cc")')  
!!        WRITE(6000+myrank,'("cc")')  
!!        WRITE(6000+myrank,'("cc")')  
!!        WRITE(6000+myrank,'("cc")')  
!!        WRITE(6000+myrank,'("cc")')  
!     END IF
!!     WRITE(6000+myrank,'(5(1pe13.5),"  NaN",1pe13.5,"  r,R,B,io,xkn,ft,<b^2>")') &
!!          & SQRT(s0)*rad_a,rad_R,borbic(0,0),ABS(iota),ABS(borbic(0,1))/eps,avb2
!!     IF(is.EQ.1) WRITE(6000+myrank,'("c         eps_eff     g11_ft      efield_u    g11_er    ex_er")')  
!!     WRITE(6000+myrank,'("cfit          NaN        NaN           NaN       NaN      NaN")')
!  ELSE IF(KNOSOS_STELLOPT) THEN
!     IF(numprocs.EQ.1) filename="mono.opt"
!     IF(numprocs.GT.1) WRITE(filename,'("mono.opt.",I2.2)') myrank
!!     OPEN(unit=6000+myrank,file=filename,form='formatted',action='write',iostat=iostat)
!     WRITE(6000+myrank,'("s Gamma_1nu Gamma_sqrtnu Gamma_sbp")')
!     WRITE(6000+myrank,'(1(1pe13.5))') eps_eff
!     WRITE(6,*) 'JL'
  END IF

  !Only continue if the goal is to compare with DKES or perform a monoenergetic calculation
  IF(.NOT.CALC_DB) RETURN

  cmul0=nu(iv0)/v(iv0)/2
  new_nb=nb(1)*cmult/cmul0
  rhostar=v(1)*Ab(1)*m_e/Zb(1)/borbic(0,0)/rad_R

!  nal=16
!  DO WHILE (nal.LE.nax)
!  nlambda=64
!  DO WHILE (nlambda.LE.nlambdax)
!  nal=64
!  DO WHILE (nal.LE.64)
!  nlambda=256
!  DO WHILE (nlambda.LE.256)
  !Scan in EFIELD
  CALCULATED_INT=.FALSE.

  IF(NEOTRANSP) efieldt=efieldt*aiota*eps

  IF(.NOT.KNOSOS_STELLOPT) THEN
     CLOSE(200+myrank)
     IF(numprocs.EQ.1) filename="results.knosos"
     IF(numprocs.GT.1) WRITE(filename,'("results.knosos.",I2.2)') myrank
     OPEN(unit=200+myrank,file=filename,form='formatted',action='write',iostat=iostat)
     WRITE(200+myrank,'("cmul efield weov wtov L11m L11p L31m L31p L33m L33p scal11&
          & scal13 scal33 max\_residual chip psip btheta bzeta vp vmag")')
  END IF
  
  DO iturn=1,2
     DO iefield=1,nefieldt
!        IF((PENTA.OR.KN_STELLOPT(1)).AND.iturn.GT.1)    EXIT
!        IF(KN_STELLOPT(1).AND.nefieldt.EQ.3) TANG_VM=.FALSE.
!        IF(KN_STELLOPT(1)) TANG_VM=.FALSE.
        DO icmul=1,ncmult
           DO ivmag=1,nvmagt
              !Calculate (v,species)-dependent constants
              CALL DKE_CONSTANTS(1,1,ZB,AB,IDUMMY,new_nb(icmul),dummy,Tb,dummy,dummy(1),.FALSE.)
!              IF(.NOT.(KN_STELLOPT(1).OR.KN_STELLOPT(2))) vmconst=vmagt(ivmag)*v
              Epsi=efieldt(iefield)*v(iv0)/psip
              IF(iturn.EQ.1) THEN
                 WRITE(iout,'(" CMUL  ",1pe13.5)') cmult(icmul)
                 WRITE(iout,'(" EFIELD",1pe13.5)') efieldt(iefield)
                 WRITE(iout,'(" VMAG  ",1pe13.5)') vmagt(ivmag)
                 WRITE(iout,'(" RHOSTAR    ",1pe13.5)') rhostar
                 IF(cmult(icmul).GT.cmul_1NU.AND.ivmag.GT.1) THEN
                    D11tab(icmul,iefield,ivmag)=D11tab(icmul,iefield,1) 
                    CYCLE
                 END IF
                 !Use different subroutines in different regimes
                 IF(cmult(icmul).GT.cmul_PS) THEN
                    CALL CALC_PS(iv0,Epsi,D11tab(icmul,iefield,ivmag))   
                 ELSE IF(cmult(icmul).GT.cmul_1NU) THEN
!                    CALL CALC_PLATEAU_OLD(iv0,Epsi,D11tab(icmul,iefield,ivmag))
                    CALL CALC_PLATEAU(iv0,Epsi,D11tab(icmul,iefield,ivmag),dn1nm)
                 ELSE
                    CONVERGED=.FALSE.
                    IF(ABS(efieldt(iefield)).LT.aiota*borbic(0,0)*eps) THEN
                       IF(ESCOTO) THEN
                          CALL CALC_LOW_COLLISIONALITY_NEW(iv0,Epsi,phi1c,Mbbnm,trMnm,&
                               & D11tab(icmul,iefield,ivmag),D31,nalphab,zeta,theta,dn1,dn1nm)
                       ELSE
                          CALL CALC_LOW_COLLISIONALITY(iv0,Epsi,phi1c,Mbbnm,trMnm,&
                               & D11tab(icmul,iefield,ivmag),nalphab,zeta,theta,dn1,dn1nm)
                       END IF
                    ELSE
                       D11tab(icmul,iefield,ivmag)=D11tab(icmul,iefield-1,ivmag)
                    END IF
                 END IF
              ELSE! IF(NEOTRANSP.OR.PENTA) THEN
                 IF(ABS(efieldt(iefield)).GT.0.2*aiota*borbic(0,0)*eps) THEN
                    D11tab(icmul,iefield,ivmag)=D11tab(icmul,iefield-1,ivmag)
                ! ELSE IF(cmult(icmul-1).LT.cmul_1NU.AND.D11tab(icmul-1,iefield,ivmag).LT.D11pla.AND.&
                !  ((D11tab(icmul,iefield,ivmag).LT.0).OR.&
                !  ((D11tab(icmul,iefield,ivmag).GT.D11tab(icmul-1,iefield,ivmag))))) THEN
                !    D11tab(icmul,iefield,ivmag)=0.5*(&
                !         D11tab(icmul-1,iefield,ivmag)*SQRT(cmult(icmul)/cmult(icmul-1))+&
                !        &D11tab(icmul-1,iefield,ivmag)*efieldt(iefield-1)*&
 	!		&SQRT(efieldt(iefield-1)/efieldt(iefield))/efieldt(iefield))
                 END IF
                 IF(NEOTRANSP.OR.PENTA) THEN
                    IF(nvmagt.EQ.1) THEN
                       WRITE(6000+myrank,'(3(1pe13.5)," NaN NaN")') &
                            & nu(iv0)/v(iv0)/2.,Epsi/v(iv0)*psip,-D11tab(icmul,iefield,ivmag)
                       WRITE(6000+myrank,'(">3                  NaN NaN 0.00000E+00 NaN NaN")')
                    ELSE
                       WRITE(6000+myrank,'(4(1pe13.5)," NaN")') &
                            & nu(iv0)/v(iv0)/2.,Epsi/v(iv0)*psip,vmagt(ivmag),-D11tab(icmul,iefield,ivmag)
                    END IF
                    IF(iefield.EQ.nefieldt.AND.icmul.EQ.ncmult.AND.ivmag.EQ.nvmagt) WRITE(6000+myrank,'("e")')
                 END IF
              END IF
           END DO
        END DO
     END DO
     !Save table of monenergetic transport coefficients (with DKES normalization)
     D11tab=D11tab*fdkes(iv0)
     lD11tab=LOG(D11tab)
!     IF(KNOSOS_STELLOPT) EXIT
  END DO

  IF(NEOTRANSP) efieldt=efieldt/(aiota*eps)

  IF(ONLY_DB) RETURN

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE CALC_DATABASE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE READ_DKES_TABLE(is)

!----------------------------------------------------------------------------------------------- 
!Read monoenergetic transport coefficients from DKES for surface s(is)
!-----------------------------------------------------------------------------------------------   

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER is
  !Others
  CHARACTER*200 line,file,dir_efield(nefieldt),dir_cmul(ncmult)
  INTEGER iefield,icmul,iostat,nefield,ncmul
  REAL*8 D11p,D11m,D31p,D31m,dummy
  !Database used by DKES
  INTEGER, PARAMETER :: ncmuld=18
  INTEGER, PARAMETER :: nefieldd=9

  IF(nefieldt.LT.9.OR.ncmult.LT.18) RETURN
  !Folders for different values of EFIELD
  dir_efield(1) ="omega_0e-0/"
  dir_efield(2) ="omega_1e-5/"
  dir_efield(3) ="omega_3e-5/"
  dir_efield(4) ="omega_1e-4/"
  dir_efield(5) ="omega_3e-4/"
  dir_efield(6) ="omega_1e-3/"
  dir_efield(7) ="omega_3e-3/"   !these lines may produce warnings depending on the compiler
  dir_efield(8) ="omega_1e-2/"
  dir_efield(9) ="omega_3e-2/"
!  dir_efield(10)="omega_1e-1/"  
  !Subfolders of omega_?e-? for different values of CMUL
  dir_cmul(1) ="cl_3e+2/"
  dir_cmul(2) ="cl_1e+2/"
  dir_cmul(3) ="cl_3e+1/"
  dir_cmul(4) ="cl_1e+1/"
  dir_cmul(5) ="cl_3e-0/"
  dir_cmul(6) ="cl_1e-0/"
  dir_cmul(7) ="cl_3e-1/"
  dir_cmul(8) ="cl_1e-1/"
  dir_cmul(9) ="cl_3e-2/"
  dir_cmul(10)="cl_1e-2/"
  dir_cmul(11)="cl_3e-3/"
  dir_cmul(12)="cl_1e-3/"
  dir_cmul(13)="cl_3e-4/"
  dir_cmul(14)="cl_1e-4/"
  dir_cmul(15)="cl_3e-5/"
  dir_cmul(16)="cl_1e-5/"
  dir_cmul(17)="cl_3e-6/"
  dir_cmul(18)="cl_1e-6/"

  file=TRIM(DIRDB)//TRIM(DIRS(is))
  WRITE(iout,*) 'Reading DKES output in folder ',TRIM(file) 
  !Read data
  nefield=0  
  ncmul=0
  D11pla=1e10
  DO iefield=1,nefieldd
     DO icmul=1,ncmuld
        file=TRIM(DIRDB)//TRIM(DIRS(is))//TRIM(dir_efield(iefield))//TRIM(dir_cmul(icmul))//"results.data"
        OPEN(unit=1,file=TRIM(file),action='read',iostat=iostat) 
        IF (iostat.EQ.0) THEN 
           WRITE(iout,*) 'Reading file ',file
           IF(icmul.EQ.1) nefield=nefield+1
           IF(iefield.EQ.1) ncmul=ncmul+1
           READ(1,*) line
           READ(1,*) line
           READ(1,*) dummy,dummy,dummy,dummy,D11m,D11p,&
                & D31m,D31p,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy
           CLOSE(1)
           !Logarithms are stored
           IF(D11p/D11m.GT.1E4.OR.D11p.LT.0) D11p=EXP(lD11dkes1(icmul-1,iefield))*cmult(icmul-1)/cmult(icmul)
           IF(iefield.EQ.1.AND.0.5*(D11m+D11p).LT.D11pla) D11pla=0.5*(D11m+D11p)
           WRITE(10000+myrank,'("-100 ",3(1pe13.5))') &
                & cmult(icmul),efieldt(iefield),0.5*(D11m+D11p)
           lD11dkes1(icmul,iefield)=LOG(0.5*(D11p+D11m)) !both values are the same, which means
           lD11dkes2(icmul,iefield)=LOG(0.5*(D11p+D11m)) !ignoring error bars in D11
           D31dkes(icmul,iefield)=0.5*(D31p+D31m)
        END IF
     END DO
  END DO
  IF(nefield.EQ.0) THEN
     WRITE(iout,*) 'No DKES output found'
  ELSE
     nefieldt=nefield
     ncmult  =ncmul
     DKES_READ=.TRUE.
  END IF

END SUBROUTINE READ_DKES_TABLE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE INTERP_DATABASE(it,jv,Epsi,D11,D31,knososdb)

!----------------------------------------------------------------------------------------------- 
!Calculate monoenergetic transport coefficient D11 from database, interpolating
!at cmul=nu(jv)/v(jv) and efield=Epsi*psip/v(jv). The database comes from KNOSOS if
!knososdb.EQ..TRUE. and from DKES otherwise
!-----------------------------------------------------------------------------------------------

  USE GLOBAL
  USE KNOSOS_STELLOPT_MOD
  IMPLICIT NONE
  !Input
  LOGICAL knososdb
  INTEGER it,jv
  REAL*8 Epsi
  !Output
  REAL*8 D11,D31
  !Others
  INTEGER nvmagtt
  REAL*8 efield,cmul,vmag,lD11,lD11dkes(ncmult,nefieldt)
  REAL*8 lvmagtt((nvmagt-1)/2),lD11tabt(ncmult,nefieldt,(nvmagt-1)/2)
  !Time
  CHARACTER*30, PARAMETER :: routine="INTERP_DATABASE"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  IF(knososDB) THEN
     WRITE(iout,*) 'Interpolating KNOSOS coefficients'
     vmag=vmconst(jv)/v(jv)
     IF(Epsi.LT.0) vmag=-vmag        
     nvmagtt=(nvmagt-1)/2
     IF(vmag.LT.0) THEN
        lvmagtt=LOG(ABS(vmagt(1:nvmagtt)))
        lD11tabt=lD11tab(1:ncmult,1:nefieldt,1:nvmagtt)
     ELSE
        lvmagtt=LOG(ABS(vmagt(nvmagt-nvmagtt+1:nvmagt)))
        lD11tabt=lD11tab(1:ncmult,1:nefieldt,nvmagt-nvmagtt+1:nvmagt)
     END IF     
  ELSE
     WRITE(iout,*) 'Interpolating DKES coefficients'
  END IF
  efield=ABS(Epsi*psip/v(jv))
  cmul=nu(jv)/v(jv)/2 !careful with the definition of collision frequency

  IF(it.EQ.-1) THEN !irrelevant, as lD11dkes1=lD11dkes2, but could be used for errorbars
     lD11dkes(1:ncmult,1:nefieldt)=lD11dkes1(1:ncmult,1:nefieldt)
  ELSE
     lD11dkes(1:ncmult,1:nefieldt)=lD11dkes2(1:ncmult,1:nefieldt)
  END IF

  IF(knososDB) THEN     !Interpolate from a database calculated with KNOSOS
     IF(nvmagt.EQ.1) THEN
       IF(efield.GT.0.1*efieldt(2)) THEN!  IF(efield.GT.1E-10) THEN
           CALL BILAGRANGE(lcmult(1:ncmult),lefieldt(2:nefieldt),lD11tab(1:ncmult,2:nefieldt,1),&
                & ncmult,nefieldt-1,LOG(cmul),LOG(efield),lD11,2)
        ELSE
           CALL LAGRANGE(lcmult(1:ncmult),lD11tab(1:ncmult,1,1),&
                & ncmult,LOG(cmul),lD11,2)
        END IF
     ELSE
      IF(efield.GT.0.1*efieldt(2)) THEN!  IF(efield.GT.1E-10) THEN
         
           CALL TRILAGRANGE(lcmult(1:ncmult),lefieldt(2:nefieldt),lvmagtt(1:nvmagtt),lD11tabt(1:ncmult,2:nefieldt,1:nvmagtt),&
                & ncmult,nefieldt-1,nvmagtt,LOG(cmul),LOG(efield),LOG(ABS(vmag)),lD11,1)
        ELSE
           CALL BILAGRANGE(lcmult(1:ncmult),lvmagtt(1:nvmagtt),lD11tabt(1:ncmult,1,1:nvmagtt),&
             & ncmult,nvmagtt,LOG(cmul),LOG(ABS(vmag)),lD11,1)
        END IF
     END IF
  ELSE                  !Interpolate from DKES
     IF(efield.GT.0.1*efieldt(2)) THEN
        CALL BILAGRANGE(lcmult(1:ncmult),lefieldt(2:nefieldt),lD11dkes(1:ncmult,2:nefieldt),&
             & ncmult,nefieldt-1,LOG(cmul),LOG(efield),lD11,1)
        CALL BILAGRANGE(lcmult(1:ncmult),lefieldt(2:nefieldt),D31dkes(1:ncmult,2:nefieldt),&
             & ncmult,nefieldt-1,LOG(cmul),LOG(efield),D31,1)
     ELSE
        !If the radial electric field is close to zero, use only coefficients of efield=0
        CALL LAGRANGE(lcmult(1:ncmult),lD11dkes(1:ncmult,1),&
             & ncmult,LOG(cmul),lD11,1)
        CALL LAGRANGE(lcmult(1:ncmult),D31dkes(1:ncmult,1),&
             & ncmult,LOG(cmul),D31,1)     
     END IF
    
  END IF

  !After logarithmic interpolation, take exponential
  D11=EXP(lD11)

  IF(.NOT.KNOSOS_STELLOPT) WRITE(200+myrank,'(3(1pe13.5)," NaN ",2(1pe13.5)," &
          & NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN")') &
          & cmul,Epsi*psip/v(jv),vmconst(jv)/v(jv),D11,D11
  WRITE(10000+myrank,'(I3,7(1pe13.5))') it,cmul,Epsi*psip/v(jv),vmconst(jv)/v(jv),D11,D31,&
       & weight(jv)/fdkes(jv),weight(jv)/fdkes(jv)*v(jv)*v(jv)

  !Normalize
  D11=D11/fdkes(jv)
  D31=D31/fdkes2(jv)

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)
  
END SUBROUTINE INTERP_DATABASE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
