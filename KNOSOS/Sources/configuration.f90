!Read and process the magnetic configuration

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE READ_BFIELD(s0)
  
!-------------------------------------------------------------------------------------------------
!Read and process the quantities at surface s=s0 (and s0+/-ds) of the magnetic configuration
!-------------------------------------------------------------------------------------------------
  USE GLOBAL
  USE KNOSOS_STELLOPT_MOD
  IMPLICIT NONE
  !Input
  REAL*8 s0 !s0 labels the flux surface s=psi/psi_{LCMS}, where psi is the toroidal flux
  !Variables from "ddkes2.data"
  NAMELIST /datain/mpol,ntor,lalpha,ipmb,mmnn,meshtz,idisk,lfout,ifscl,nrun,cmul,&
       & efield,nzperiod,chip,psip,mpolb,ntorb,nvalsb,borbi,btheta,bzeta,ibbi
  CHARACTER*100 serr
  INTEGER n,m,iostat,is
  INTEGER mmnn(mpold+2,ntord),nvalsb(ntorbd),mpol,ntor,lalpha,ipmb,meshtz,idisk,lfout,ifscl,nrun,ibbi
  REAL*8 cmul,efield
  LOGICAL addkes_read
  !Variables from "boozmn.data"
  LOGICAL boozmndata_read
  INTEGER imn,jmn
!  REAL*8 torbic(-ntorbd:ntorbd,0:mpolbd)
!  REAL*8, ALLOCATABLE :: bmnc_js(:),bmns_js(:)
  !Variables from "boozer.txt"
  LOGICAL booztxt_read
  !Find main helicity
  INTEGER hel_Nmax,hel_Mmax!,ihel,imax,is,js
  !Time
  CHARACTER*30, PARAMETER :: routine="READ_BFIELD"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart
#ifdef MPIandPETSc
  !Others
  INTEGER ierr
  INCLUDE "mpif.h"
#endif
  
  CALL CPU_TIME(tstart)

  boozmndata_read=.FALSE.
  booztxt_read=.FALSE.
  STELL_ANTISYMMETRIC=.FALSE.
  
  CALL READ_BOOZMNDATA(s0,boozmndata_read)

  IF(.NOT.boozmndata_read) CALL READ_BOOZMNNC(s0,boozmndata_read)

  IF(.NOT.boozmndata_read) CALL READ_BOOZTXT(s0,booztxt_read)

  IF(.NOT.booztxt_read.AND..NOT.boozmndata_read) CALL READ_DDKESDATA(s0)

  CALL READ_ADDDDKESDATA(addkes_read)

  DO imn=1,mnboz_b
     IF(ixm_b(imn).EQ.0.AND.ixn_b(imn).LT.0) THEN
        DO jmn=1,mnboz_b      !For m=0, keep only n>0, possibly redundant
           IF(ixm_b(jmn).EQ.0.AND.ixn_b(jmn).EQ.-ixn_b(imn)) THEN
              bmnc_b(jmn,:)=bmnc_b(jmn,:)+bmnc_b(imn,:)
              rmnc_b(jmn,:)=rmnc_b(jmn,:)+rmnc_b(imn,:)
              zmns_b(jmn,:)=zmns_b(jmn,:)-zmns_b(imn,:)
              pmns_b(jmn,:)=pmns_b(jmn,:)-pmns_b(imn,:)
              bmnc_b(imn,:)=0
              rmnc_b(imn,:)=0
              zmns_b(imn,:)=0
              pmns_b(imn,:)=0
              IF(STELL_ANTISYMMETRIC) THEN
                 bmns_b(jmn,:)=bmns_b(jmn,:)-bmns_b(imn,:)
                 rmns_b(jmn,:)=rmns_b(jmn,:)-rmns_b(imn,:)
                 zmnc_b(jmn,:)=zmnc_b(jmn,:)+zmnc_b(imn,:)
                 pmnc_b(jmn,:)=pmnc_b(jmn,:)+pmnc_b(imn,:)
                 bmns_b(imn,:)=0
                 rmns_b(imn,:)=0
                 zmnc_b(imn,:)=0
                 pmnc_b(imn,:)=0
              END IF
           END IF
        END DO
     END IF
  END DO
  
  CALL RESCALE()

  IF(bvco_b(jsn).LT.0) WRITE(iout,*) 'Configuration with dB/dl<0'
  !Change from left-handed to a right-handed coordinate system
  buco_b=-buco_b
  iota_b=-iota_b
  bmns_b=-bmns_b
  pmns_b=-pmns_b
  rmns_b=-rmns_b
  zmns_b=-zmns_b
  !If torflux(a)<torflux(0), change direction of radial coordinate
  IF(torflux.LT.0) THEN
!     torflux=-torflux
     psip=-psip
  END IF
  IF(NEOTRANSP) THEN
     torflux=torflux*rad_a*sqrt(s0)/psip
     psip=rad_a*sqrt(s0)
  END IF  
  dpsidr=psip

  !Keep the largest TRUNCATE_B values of Bmn (if positive)
!!$  IF(TRUNCATE_B>0.AND.TRUNCATE_B.LT.mnboz_b) THEN !JL Check if this works fine
!!$     STOP
!!$     ALLOCATE(bmnc_js(mnboz_b),bmns_js(mnboz_b))
!!$     bmnc_js=0
!!$     WRITE(iout,'(" Truncating spectrum of magnetic field strength from",&
!!$          & I5," to ",I3," largest cosines")') mnboz_b,TRUNCATE_B
!!$     DO js=js0,js1,js1-js0
!!$        bmnc_js(1:mnboz_b)  =bmnc_b(1:mnboz_b,js)
!!$        bmnc_b(1:mnboz_b,js)=0
!!$        IF(STELL_ANTISYMMETRIC) THEN
!!$           bmns_js(1:mnboz_b)  =bmns_b(1:mnboz_b,js)
!!$           bmns_b(1:mnboz_b,js)=0
!!$        END IF
!!$        DO ihel=1,MIN(TRUNCATE_B,mnboz_b)
!!$           imax=MAXLOC(ABS(bmnc_js(:)),1)
!!$           bmnc_b(ihel,js)=bmnc_js(imax)
!!$           bmnc_js(imax)=0
!!$           IF(STELL_ANTISYMMETRIC) THEN
!!$              bmns_b(ihel,js)=bmns_js(imax)
!!$              bmns_js(imax)=0
!!$           END IF
!!$        END DO
!!$     END DO
!!$  END IF

  IF(boozmndata_read.OR.booztxt_read) CALL JACOBIAN_AND_3DGRID(s0)

  CALL INTERPOLATE_FIELD(s0,boozmndata_read.OR.booztxt_read,ZERO)

  !Find maximum helicity and truncate spectra accordingly
  hel_Nmax=0
  hel_Mmax=0
  DO m=0,mpolb
     DO n=-ntorb,ntorb
        IF(ABS(borbic(n,m)).LT.PREC_B.AND.ABS(borbis(n,m)).LT.PREC_B) CYCLE
        IF(ABS(n).GT.hel_Nmax) hel_Nmax=ABS(n)
        IF(m.GT.hel_Mmax) hel_Mmax=ABS(m)
     END DO
  END DO
 IF(MAL.GT.0.AND.QN.AND..NOT.PLOT_XYZ) THEN
    hel_Nmax=MIN(hel_NMAX,MAL/NSAMP)
    hel_Mmax=MIN(hel_MMAX,MAL/NSAMP)
 END IF
 IF(FAST_IONS) THEN
    hel_Nmax=ntorb
    hel_Mmax=mpolb
 END IF
  IF(hel_NMAX.LT.ntorb.OR.hel_MMAX.LT.mpolb) WRITE(iout,&
     & '(" Automatically truncating spectrum of magnetic field strength from &
     & (nmax,mmax)=(",I3,",",I3,") to (",I2,",",I2,")")') ntorb,mpolb,hel_Nmax,hel_Mmax
  ntorb=MIN(hel_Nmax,ntorb)
  mpolb=MIN(hel_Mmax,mpolb)
  !Set, according to the required precision, steps along the Boozer angles, etc
  dzstep=TWOPI/(nzperiod*hel_NMAX*100)
  
  !Copy some information to arrays (Bmn, etc) with only one index
  borbic0=borbic
  borbis0=borbis
  CALL FILL_NM()
  !If USE_B0 is .TRUE. , use [Calvo 2017 PPCF], but B_0 MUST be omnigenous!
  !If USE_B0 is .FALSE., calculation using the total B=B_0+B_1
  IF(USE_B0.OR.KN_STELLOPT(6).OR.KN_STELLOPT(7).OR.KN_STELLOPT(8).OR.KN_STELLOPT(9)) CALL CALC_B0()
  CALL FILL_NM()
  
  !If rescaled field, write in file
  IF(ABS(FI*FE*FB*FR-1).GT.ALMOST_ZERO) THEN
     borbi=borbic
     OPEN(unit=1,file='ddkes3.data',form='formatted',action='write',iostat=iostat)
     WRITE(1,nml=datain)
     CLOSE(1)
  END IF  

  CALL INTERPOLATE_FIELD(s0,boozmndata_read.OR.booztxt_read,ZERO)
  
  !Calculate quantities on the flux surface, such as <B^2>, <B> or the location of the maximum of B
  CALL FILL_BGRID(MAL,MAL,s0,.FALSE.)
  IF(addkes_read.OR.USE_B0) CALL FILL_BGRID(MAL,MAL,s0,.FALSE.)
  CALL FILL_BGRID(MAL,MAL,s0,.TRUE.)
  borbic0=borbic
  borbis0=borbis
  dborbic0dpsi=dborbicdpsi
  dborbis0dpsi=dborbisdpsi
  CALL FILL_NM()

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)
  IF(ONLY_B0) THEN 
     serr="B0 calculated"
     CALL END_ALL(serr,.FALSE.)
  END IF

  IF(FAST_IONS) THEN
     IF(ALLOCATED(phi_b)) DEALLOCATE(phi_b,phip_b)
     ALLOCATE(phi_b(ns_b),phip_b(ns_b))
     phi_b=0
     phip_b=0
     IF(s_b(1).LT.0) s_b(1)=-s_b(1)
     DO is=1,ns_b
        CALL INTERPOLATE_FIELD(s_b(is),boozmndata_read.OR.booztxt_read,ZERO)
        CALL READ_PROFILE(s_b(is),"ph",phi_b(is),phip_b(is),1)
     END DO
     IF(ABS(phi_b(1)).LT.SMALL) THEN
        DO is=ns_b-1,1,-1
           phi_b(is)=phi_b(is+1)-0.5*(phip_b(is+1)+phip_b(is))*(s_b(is+1)-s_b(is))*atorflux
           WRITE(iout,*) 'ph',is,ns_b,phi_b(is),phip_b(is+1),phip_b(is),s_b(is+1),s_b(is)
        END DO
     END IF
  END IF
  CALL INTERPOLATE_FIELD(s0,boozmndata_read.OR.booztxt_read,ZERO)
  IF(FAST_IONS) THEN
     phi_b=phi_b-phi0
     phi0=0
  END IF
  
END SUBROUTINE READ_BFIELD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE READ_BOOZMNDATA(s0,boozmndata_read)
!-------------------------------------------------------------------------------------------
!Read output from BOOZER_XFORM, binary file
!-------------------------------------------------------------------------------------------

  USE GLOBAL
  USE KNOSOS_STELLOPT_MOD
  IMPLICIT NONE
  !Input
  REAL*8 s0
  !Output
  LOGICAL boozmndata_read
  !Others
  CHARACTER*100 serr
  CHARACTER*38 version
  INTEGER iostat,nsval,jsize,js,is0,imn
  REAL*8 aspect_b,rmax_b,rmin_b,betaxis_b
  
  OPEN(unit=1,file='boozmn.data',form='unformatted',iostat=iostat,action='read')
  IF(iostat.NE.0) OPEN(unit=1,file='../boozmn.data',form='unformatted',iostat=iostat,action='read')
  IF(iostat.NE.0) RETURN
  
  boozmndata_read=.TRUE.
  WRITE(iout,*) 'File "boozmn.data" found'
  !Read scalar values
  READ(1) nfp_b,ns_b,aspect_b,rmax_b,rmin_b,betaxis_b
  IF(ALLOCATED(s_b)) DEALLOCATE(s_b,js_b,iota_b,pres_b,beta_b,psip_b,psi_b,bvco_b,buco_b,&
       & bmnc_b,rmnc_b,zmns_b,pmns_b,rmns_b,zmnc_b,pmnc_b,bmns_b,gmnc_b,ixn_b,ixm_b)
  ALLOCATE(s_b(ns_b),js_b(ns_b),iota_b(ns_b),pres_b(ns_b),beta_b(ns_b),&
       &psip_b(ns_b),psi_b(ns_b),bvco_b(ns_b),buco_b(ns_b))
  !Give values to quantities at the magnetic axis
  s_b(1)=0
  iota_b(1)=0
  pres_b(1)=0
  beta_b(1)=0
  psip_b(1)=0
  psi_b(1)=0
  bvco_b(1)=0
  buco_b(1)=0
  !Read profiles
  DO nsval=2,ns_b
     READ(1) iota_b(nsval),pres_b(nsval),beta_b(nsval),&
          & psip_b(nsval),psi_b(nsval),bvco_b(nsval),buco_b(nsval)
  END DO
  READ(1) mboz_b,nboz_b,mnboz_b,jsize
  READ(1) version
  ALLOCATE(bmnc_b(mnboz_b,ns_b), rmnc_b(mnboz_b,ns_b),zmns_b(mnboz_b,ns_b), &
       &pmns_b(mnboz_b,ns_b),gmnc_b(mnboz_b,ns_b),ixm_b(mnboz_b),ixn_b(mnboz_b))
  ALLOCATE(rmns_b(mnboz_b,ns_b),zmnc_b(mnboz_b,ns_b),pmnc_b(mnboz_b,ns_b),bmns_b(mnboz_b,ns_b))
  ixn_b=0
  ixm_b=0
  rmnc_b=0
  zmns_b=0
  pmns_b=0
  bmnc_b=0
  gmnc_b=0     
  READ(1) ixn_b(1:mnboz_b),ixm_b(1:mnboz_b)
  s_b(1)=0
  is0=0
  DO js=1,jsize*2
     READ(1,iostat=iostat) nsval
     IF(iostat.NE.0) EXIT
     IF((nsval.gt.ns_b).or.(nsval.le.0)) CYCLE
     s_b(nsval)=(nsval-1.)/(ns_b-1)
     is0=is0+1
     js_b(is0)=nsval
     !Read the Boozer spectra at every flux-surface
     READ(1,iostat=iostat) bmnc_b(1:mnboz_b,nsval),rmnc_b(1:mnboz_b,nsval), &
          &zmns_b(1:mnboz_b,nsval),pmns_b(1:mnboz_b,nsval),gmnc_b(1:mnboz_b,nsval)
  END DO
  jsn=is0
  CLOSE(1)
  WRITE(iout,*) 'File "boozmn.data" read'
  DO imn=1,mnboz_b
     IF(ixn_b(imn).EQ.0.AND.ixm_b(imn).EQ.0) THEN
        rad_R=rmnc_b(imn,2)
        IF(aspect_b.LT.1e-10) THEN
           serr="aspect_b=0"
           CALL END_ALL(serr,.FALSE.)
        END IF
        rad_a=rad_R/aspect_b
        EXIT
     END IF
  END DO

  pmns_b=-pmns_b
  ixn_b=ixn_b/nfp_b
  torflux=psi_b(ns_b)/TWOPI
  psip=2*torflux*sqrt(s0)/rad_a

END SUBROUTINE READ_BOOZMNDATA


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE READ_BOOZMNNC(s0,boozmndata_read)
!-------------------------------------------------------------------------------------------
!Read output from BOOZER_XFORM, netcdf file, by Satake
!-------------------------------------------------------------------------------------------
  
  USE GLOBAL
  USE KNOSOS_STELLOPT_MOD
  IMPLICIT NONE
  !Input
  REAL*8 s0
  !Output
  LOGICAL boozmndata_read
  !Others
  CHARACTER*100 serr
  CHARACTER*38 version
  INTEGER nsval,jsize,is0,imn
  REAL*8 aspect_b,rmax_b,rmin_b,betaxis_b
  !netcdf by Satake
  INCLUDE "netcdf.inc"
  CHARACTER*60 filename
  CHARACTER*40 varname
  INTEGER  ncid,status_nc,rhid,idimid,ilist
  INTEGER, ALLOCATABLE :: jlist(:),idx_b(:)
  REAL*8, ALLOCATABLE :: packed2d(:,:,:)

  IF(KNOSOS_STELLOPT) THEN
     filename=TRIM(KN_DKESFILE)
  ELSE
     filename='boozmn.nc'
  END IF
  varname='***'
  status_nc=nf_open(filename, nf_nowrite, ncid)
  IF(status_nc.NE.nf_noerr) THEN!.AND..NOT.KNOSOS_STELLOPT) THEN
     filename='../boozmn.nc'
     status_nc=nf_open(filename,nf_nowrite,ncid)
  END IF
  IF(status_nc.NE.nf_noerr) RETURN
  boozmndata_read=.TRUE.
  WRITE(iout,*) 'File "',TRIM(filename),'" found'   
  !Read scalar values
  varname='nfp_b'
  status_nc=nf_inq_varid(ncid,varname,rhid)
  IF(status_nc.NE.nf_noerr) CALL QUIT_READ(filename,varname,status_nc) 
  status_nc=nf_get_var_int(ncid,rhid,nfp_b)
  varname='ns_b'
  status_nc=nf_inq_varid(ncid,varname,rhid)
  IF(status_nc.NE.nf_noerr) CALL QUIT_READ(filename,varname,status_nc) 
  status_nc=nf_get_var_int(ncid,rhid,ns_b)
  varname='aspect_b'
  status_nc=nf_inq_varid(ncid,varname,rhid)
  IF(status_nc.NE.nf_noerr) CALL QUIT_READ(filename,varname,status_nc) 
  status_nc=nf_get_var_double(ncid,rhid,aspect_b)
  varname='rmax_b'
  status_nc=nf_inq_varid(ncid,varname,rhid)
  IF(status_nc.NE.nf_noerr) CALL QUIT_READ(filename,varname,status_nc) 
  status_nc=nf_get_var_double(ncid,rhid,rmax_b)        
  varname='rmin_b'
  status_nc=nf_inq_varid(ncid,varname,rhid)
  IF(status_nc.NE.nf_noerr) CALL QUIT_READ(filename,varname,status_nc) 
  status_nc=nf_get_var_double(ncid,rhid,rmin_b)
  varname='betaxis_b'
  status_nc=nf_inq_varid(ncid,varname,rhid)
  IF(status_nc.NE.nf_noerr) CALL QUIT_READ(filename,varname,status_nc) 
  status_nc=nf_get_var_double(ncid,rhid,betaxis_b)
  varname='mboz_b'
  status_nc=nf_inq_varid(ncid,varname,rhid)
  IF(status_nc.NE.nf_noerr) CALL QUIT_READ(filename,varname,status_nc) 
  status_nc=nf_get_var_int(ncid,rhid,mboz_b)
  varname='nboz_b'
  status_nc=nf_inq_varid(ncid,varname,rhid)
  IF(status_nc.NE.nf_noerr) CALL QUIT_READ(filename,varname,status_nc) 
  status_nc=nf_get_var_int(ncid,rhid,nboz_b)
  varname='mnboz_b'
  status_nc=nf_inq_varid(ncid,varname,rhid)
  IF(status_nc.NE.nf_noerr) CALL QUIT_READ(filename,varname,status_nc) 
  status_nc=nf_get_var_int(ncid,rhid,mnboz_b)
  varname='version'
  status_nc=nf_inq_varid(ncid,varname,rhid)
  IF(status_nc.NE.nf_noerr) CALL QUIT_READ(filename,varname,status_nc) 
  status_nc=nf_get_var_text(ncid,rhid,version)
  varname='pack_rad'
  status_nc=nf_inq_dimid(ncid,varname,idimid)
  IF(status_nc.NE.nf_noerr) CALL QUIT_READ(filename,varname,status_nc) 
  status_nc=nf_inq_dim(ncid,idimid,varname,jsize)
  IF(ALLOCATED(s_b)) DEALLOCATE(s_b,js_b,iota_b,pres_b,beta_b,psip_b,psi_b,bvco_b,buco_b,&
       & bmnc_b,rmnc_b,zmns_b,pmns_b,rmns_b,zmnc_b,pmnc_b,bmns_b,gmnc_b,ixn_b,ixm_b)!,jlist,idx_b,packed2d)
  ! note : jsize is the number of flux surfaces contained in booz_xform (/=ns_b)
  ALLOCATE(s_b(ns_b),js_b(ns_b),iota_b(ns_b),pres_b(ns_b),beta_b(ns_b),&
       &psip_b(ns_b),psi_b(ns_b),bvco_b(ns_b),buco_b(ns_b))
  ALLOCATE(ixm_b(mnboz_b),ixn_b(mnboz_b),gmnc_b(mnboz_b,ns_b))
  ALLOCATE(bmnc_b(mnboz_b,ns_b),rmnc_b(mnboz_b,ns_b),zmns_b(mnboz_b,ns_b),pmns_b(mnboz_b,ns_b))
  ALLOCATE(rmns_b(mnboz_b,ns_b),zmnc_b(mnboz_b,ns_b),pmnc_b(mnboz_b,ns_b),bmns_b(mnboz_b,ns_b))
  ALLOCATE(jlist(jsize),packed2d(mnboz_b,jsize,5),idx_b(ns_b))
  !Give values to quantities at the magnetic axis
  s_b=-EPSILON(0d0)
  iota_b(1)=0
  pres_b(1)=0
  beta_b(1)=0
  psip_b(1)=0
  psi_b(1)=0
  bvco_b(1)=0
  buco_b(1)=0
  !Set to zero
  ixn_b=0
  ixm_b=0
  rmnc_b=0
  zmns_b=0
  pmns_b=0
  bmnc_b=0
  gmnc_b=0 
  !Read profiles
  varname='iota_b'
  status_nc=nf_inq_varid(ncid,varname,rhid)
  IF(status_nc.NE.nf_noerr) CALL QUIT_READ(filename,varname,status_nc) 
  status_nc=nf_get_var_double(ncid,rhid,iota_b)
  varname='pres_b'
  status_nc=nf_inq_varid(ncid,varname,rhid)
  IF(status_nc.NE.nf_noerr) CALL QUIT_READ(filename,varname,status_nc) 
  status_nc=nf_get_var_double(ncid,rhid,pres_b)
  varname='beta_b'
  status_nc=nf_inq_varid(ncid,varname,rhid)
  IF(status_nc.NE.nf_noerr) CALL QUIT_READ(filename,varname,status_nc) 
  status_nc=nf_get_var_double(ncid,rhid,beta_b)
  varname='phip_b'
  status_nc=nf_inq_varid(ncid,varname,rhid)
  IF(status_nc.NE.nf_noerr) CALL QUIT_READ(filename,varname,status_nc) 
  status_nc=nf_get_var_double(ncid,rhid,psip_b)  !name difference 
  varname='phi_b'
  status_nc=nf_inq_varid(ncid,varname,rhid)
  IF(status_nc.NE.nf_noerr) CALL QUIT_READ(filename,varname,status_nc) 
  status_nc=nf_get_var_double(ncid,rhid,psi_b)   !name difference 
  varname='bvco_b'
  status_nc=nf_inq_varid(ncid,varname,rhid)
  IF(status_nc.NE.nf_noerr) CALL QUIT_READ(filename,varname,status_nc) 
  status_nc=nf_get_var_double(ncid,rhid,bvco_b)
  varname='buco_b'
  status_nc=nf_inq_varid(ncid,varname,rhid)
  IF(status_nc.NE.nf_noerr) CALL QUIT_READ(filename,varname,status_nc) 
  status_nc=nf_get_var_double(ncid,rhid,buco_b)
  varname='ixm_b'
  status_nc=nf_inq_varid(ncid,varname,rhid)
  IF(status_nc.NE.nf_noerr) CALL QUIT_READ(filename,varname,status_nc) 
  status_nc=nf_get_var_int(ncid,rhid,ixm_b)
  varname='ixn_b'
  status_nc=nf_inq_varid(ncid,varname,rhid)
  IF(status_nc.NE.nf_noerr) CALL QUIT_READ(filename,varname,status_nc) 
  status_nc=nf_get_var_int(ncid,rhid,ixn_b)
  varname='jlist'
  status_nc=nf_inq_varid(ncid,varname,rhid)
  IF(status_nc.NE.nf_noerr) CALL QUIT_READ(filename,varname,status_nc) 
  status_nc=nf_get_var_int(ncid,rhid,jlist)
  idx_b=0
  DO ilist=1,jsize
     nsval=jlist(ilist)
     idx_b(nsval)=ilist
  END DO
  !2d arrays (only jlist-ed radial nodes store in file)
  varname='bmnc_b'
  status_nc=nf_inq_varid(ncid,varname,rhid)
  IF(status_nc.NE.nf_noerr) THEN
     varname='bmn_b'  !some versions has different name 
     status_nc=nf_inq_varid(ncid,varname,rhid)
     IF(status_nc.NE.nf_noerr) THEN
        filename ='bmnc_b or bmn_b'
        CALL QUIT_READ(filename,varname,status_nc)
     END IF
  END IF
  status_nc=nf_get_var_double(ncid, rhid, packed2d(:,:,1))
  varname='rmnc_b'
  status_nc=nf_inq_varid(ncid,varname,rhid)
  IF(status_nc.NE.nf_noerr) CALL QUIT_READ(filename,varname,status_nc) 
  status_nc=nf_get_var_double(ncid, rhid, packed2d(:,:,2))
  varname='zmns_b'
  status_nc=nf_inq_varid(ncid,varname,rhid)
  IF(status_nc.NE.nf_noerr) CALL QUIT_READ(filename,varname,status_nc) 
  status_nc=nf_get_var_double(ncid, rhid, packed2d(:,:,3))
  varname='pmns_b'
  status_nc=nf_inq_varid(ncid,varname,rhid)
  IF(status_nc.NE.nf_noerr) CALL QUIT_READ(filename,varname,status_nc) 
  status_nc=nf_get_var_double(ncid, rhid, packed2d(:,:,4))
  varname='gmn_b'
  status_nc=nf_inq_varid(ncid,varname,rhid)
  IF(status_nc.NE.nf_noerr) CALL QUIT_READ(filename,varname,status_nc) 
  status_nc=nf_get_var_double(ncid, rhid, packed2d(:,:,5))
  is0=0
  DO nsval=2,ns_b
     ilist=idx_b(nsval)
     IF(ilist/=0) THEN
        bmnc_b(:,nsval)=packed2d(:,ilist,1)
        rmnc_b(:,nsval)=packed2d(:,ilist,2)
        zmns_b(:,nsval)=packed2d(:,ilist,3)
        pmns_b(:,nsval)=packed2d(:,ilist,4)
        gmnc_b(:,nsval)=packed2d(:,ilist,5)
        s_b(nsval)=(nsval-1.)/(ns_b-1)
        is0=is0+1
        js_b(is0)=nsval
        ! half / full-mesh problem in s_b and psi_b is not resolved here
        ! This should be : s_b(nsval)=(nsval-1.5)/(ns_b-1) 
     END IF
  END DO
  jsn=is0
  WRITE(iout,*) 'File "boozmn.nc" read'
  DO imn=1,mnboz_b
     IF(ixn_b(imn).EQ.0.AND.ixm_b(imn).EQ.0) THEN
        rad_R=rmnc_b(imn,2)
        IF(aspect_b.LT.1e-10) THEN
           serr="aspect_b=0"
           CALL END_ALL(serr,.FALSE.)
        END IF
        rad_a=rad_R/aspect_b
        EXIT
     END IF
  END DO
  DEALLOCATE(jlist,packed2d,idx_b)
  !Close cdf File
  status_nc=nf_close(ncid)   

  pmns_b=-pmns_b
  ixn_b=ixn_b/nfp_b
  torflux=psi_b(ns_b)/TWOPI
  psip=2*torflux*sqrt(s0)/rad_a

END SUBROUTINE READ_BOOZMNNC


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE READ_BOOZTXT(s0,booztxt_read)
!------------------------------------------------------------------------------------------- 
!Reads output from Boozer txt file
!-------------------------------------------------------------------------------------------
  USE GLOBAL
  USE KNOSOS_STELLOPT_MOD
  IMPLICIT NONE
  !Input
  REAL*8 s0
  !Output
  LOGICAL booztxt_read
  !Others
  CHARACTER*38 version
  INTEGER iostat,is,imn,strumberger
  REAL*8 int_iota
  REAL*8, PARAMETER   :: mu0_o_2pi=2.000000000e-7
  REAL*8, ALLOCATABLE :: pprime(:),sqrtg00(:)

  OPEN(unit=1,file='boozer.txt',iostat=iostat,action='read')
  IF(iostat.NE.0) OPEN(unit=1,file='../boozer.txt',iostat=iostat,action='read')
  IF(iostat.NE.0) RETURN
  booztxt_read=.TRUE.
  WRITE(iout,*) 'File "boozer.txt" found'
  Strumberger=0
  DO is=1,5
     READ(1,'(A)') version
     Strumberger=MAX(Strumberger,INDEX(version,'Strumberger'))
  END DO
  DO is=1,15
     READ(1,*,iostat=iostat) mboz_b,nboz_b,ns_b,nfp_b,torflux,rad_a,rad_R!,fdummy,fdummy,fdummy
     IF(iostat.EQ.0) EXIT
  END DO
  jsn=ns_b
  mnboz_b=MAX(mboz_b*(2*nboz_b+1)+nboz_b+1,(mboz_b+1)*(2*nboz_b+1))
  IF(ALLOCATED(s_b)) DEALLOCATE(s_b,js_b,iota_b,psip_b,psi_b,bvco_b,buco_b)
  IF(ALLOCATED(spol_b)) DEALLOCATE(spol_b)
  IF(ALLOCATED(pprime)) DEALLOCATE(pprime,sqrtg00)
  IF(ALLOCATED(bmnc_b)) DEALLOCATE(bmnc_b,rmnc_b,zmns_b)
  IF(ALLOCATED(pmns_b)) DEALLOCATE(pmns_b,rmns_b,zmnc_b,pmnc_b,bmns_b,ixm_b,ixn_b)
!    IF(ALLOCATED(spol_b)) DEALLOCATE(spol_b,pprime,sqrtg00,&
!       & bmnc_b,rmnc_b,zmns_b,pmns_b,rmns_b,zmnc_b,pmnc_b,bmns_b,ixm_b,ixn_b)
  ALLOCATE(s_b(ns_b),js_b(ns_b),iota_b(ns_b),psip_b(ns_b),psi_b(ns_b),bvco_b(ns_b),buco_b(ns_b))
  ALLOCATE(spol_b(ns_b),pprime(ns_b),sqrtg00(ns_b))
  ALLOCATE(rmnc_b(mnboz_b,ns_b),zmns_b(mnboz_b,ns_b),pmns_b(mnboz_b,ns_b),bmnc_b(mnboz_b,ns_b))
  ALLOCATE(rmns_b(mnboz_b,ns_b),zmnc_b(mnboz_b,ns_b),pmnc_b(mnboz_b,ns_b),bmns_b(mnboz_b,ns_b))
  ALLOCATE(ixm_b(mnboz_b),ixn_b(mnboz_b))
  ixn_b=0
  ixm_b=0
  rmnc_b=0
  zmns_b=0
  pmns_b=0
  bmnc_b=0
  rmns_b=0
  zmnc_b=0
  pmnc_b=0
  bmns_b=0
  !Reads the profiles Boozer spectra at every flux-surface
  DO is=1,ns_b
     js_b(is)=is
     IF(is.NE.1.AND.imn.LT.mnboz_b) THEN
        mnboz_b=imn-1
     ELSE
        READ(1,*) version
     END IF
     READ(1,*) version
     READ(1,*) s_b(is),iota_b(is),bvco_b(is),buco_b(is),pprime(is),sqrtg00(is)
     IF(is.EQ.1) THEN
        spol_b(is)=iota_b(1)*s_b(1)
        int_iota=iota_b(1)*s_b(1)
     ELSE
        int_iota=int_iota+(iota_b(is)+iota_b(is-1))*(s_b(is)-s_b(is-1))/2.
        spol_b(is)=spol_b(is-1)+(iota_b(is)+iota_b(is-1))*(s_b(is)-s_b(is-1))/2.
     END IF
     READ(1,'(A)') version
     STELL_ANTISYMMETRIC=INDEX(version,"rmns",.TRUE.).GT.0
     DO imn=1,mnboz_b
        IF(STELL_ANTISYMMETRIC) THEN
           IF(is.EQ.1.AND.imn.EQ.1) WRITE(iout,*) 'Non stellarator-symmetric equilibrium'
           READ(1,*,iostat=iostat) ixm_b(imn),ixn_b(imn),&
                & rmnc_b(imn,is),rmns_b(imn,is),zmnc_b(imn,is),zmns_b(imn,is),&
                & pmnc_b(imn,is),pmns_b(imn,is),bmnc_b(imn,is),bmns_b(imn,is)
        ELSE
           IF(is.EQ.1.AND.imn.EQ.1) WRITE(iout,*) 'Stellarator-symmetric equilibrium'
           READ(1,*,iostat=iostat) ixm_b(imn),ixn_b(imn),&
                & rmnc_b(imn,is),zmns_b(imn,is),pmns_b(imn,is),bmnc_b(imn,is)
        END IF
        IF(iostat.NE.0) EXIT
     END DO
     IF(Strumberger.GT.0) ixn_b=-ixn_b           
  END DO
  CLOSE(1)
  WRITE(iout,*) 'File "boozer.txt" read'        

  spol_b=spol_b/int_iota
  buco_b=buco_b*nfp_b*mu0_o_2pi
  bvco_b=bvco_b*nfp_b*mu0_o_2pi
  pmns_b=pmns_b*TWOPI/nfp_b
  torflux=torflux/TWOPI
  psip=2*torflux*sqrt(s0)/rad_a
                
END SUBROUTINE READ_BOOZTXT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE READ_DDKESDATA(s0)
!------------------------------------------------------------------------------------------- 
!Read output from DKES_INPUT_PREPARE
!This is necesary if there is no equilibrium file and/or also for comparison with DKES
!------------------------------------------------------------------------------------------- 
  USE GLOBAL
  USE KNOSOS_STELLOPT_MOD
  IMPLICIT NONE
  !Input/output
  REAL*8 s0
  !Others
  CHARACTER*100 serr
  INTEGER iostat,imn,n,m
  INTEGER mmnn(mpold+2,ntord),nvalsb(ntorbd),mpol,ntor,lalpha,ipmb,meshtz,idisk,lfout,ifscl,nrun,ibbi
  REAL*8 cmul,efield
  NAMELIST /datain/mpol,ntor,lalpha,ipmb,mmnn,meshtz,idisk,lfout,ifscl,nrun,cmul,&
       & efield,nzperiod,chip,psip,mpolb,ntorb,nvalsb,borbi,btheta,bzeta,ibbi
  NAMELIST /radius/rad_R,rad_a

  OPEN(unit=1,file='ddkes2.data',form='formatted',action='read',iostat=iostat)
  IF(iostat.NE.0) OPEN(unit=1,file='../ddkes2.data',form='formatted',action='read',iostat=iostat)
  IF(iostat.EQ.0) THEN
     borbic=0
     WRITE(iout,*) 'File "ddkes2.data" found'     
     READ(1,nml=datain,iostat=iostat)
     CLOSE(1)
     borbic=borbi
     iota=chip/psip
     !Use large aspect ratio expression to estimate major and minor radius
     rad_R=ABS(Bzeta/borbic(0,0))
     rad_a=ABS(psip)/borbic(0,0)
     torflux=borbic(0,0)*PI*rad_a*rad_a
     !Read major radius, minor (needed for radial derivatives)
     OPEN(unit=1,file='input.radius',action='read',iostat=iostat)
     IF(iostat.NE.0) OPEN(unit=1,file='../input.radius',action='read',iostat=iostat)
     IF(iostat.EQ.0) THEN
        WRITE(iout,*) 'File "input.radius" found'
        READ(1,nml=radius)
        CLOSE(1)  
        IF(s0.LT.0) s0=psip*psip/borbic(0,0)/borbic(0,0)/rad_a/rad_a
     END IF
  ELSE
     serr="Configuration file not found"
     CALL END_ALL(serr,.FALSE.)
  END IF

  nfp_b=nzperiod
  ns_b=1
  IF(ALLOCATED(s_b)) DEALLOCATE(s_b,js_b,iota_b,psip_b,psi_b,bvco_b,buco_b,bmnc_b,bmns_b,ixn_b,ixm_b)
  ALLOCATE(s_b(ns_b),js_b(ns_b),iota_b(ns_b),psip_b(ns_b),psi_b(ns_b),bvco_b(ns_b),buco_b(ns_b))
  s_b(1)=s0
  iota_b(1)=iota
  bvco_b(1)=bzeta
  buco_b(1)=btheta
  nboz_b=ntorb
  mboz_b=mpolb
  mnboz_b=(2*ntorb+1)*mpolb+ntorb+1
  ALLOCATE(bmnc_b(mnboz_b,1),bmns_b(mnboz_b,1),ixn_b(mnboz_b),ixm_b(mnboz_b))
  imn=0
  DO m=0,mpolb
     DO n=-ntorb,ntorb
        IF(m.EQ.0.AND.n.LT.0) CYCLE
        imn=imn+1
        ixn_b(imn)=n
        ixm_b(imn)=m
        bmnc_b(imn,1)=borbic(n,m)
     END DO
  END DO
  bmns_b=0
  
END SUBROUTINE READ_DDKESDATA


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE READ_ADDDDKESDATA(addkes_read)
!------------------------------------------------------------------------------------------- 
!Add an extra magnetic field B_1 to previously read B_0
!------------------------------------------------------------------------------------------- 
  USE GLOBAL
  !Output
  LOGICAL addkes_read
  INTEGER iostat,n,m,imn
  REAL*8 borbic_add(-ntorbd:ntorbd,0:mpolbd),borbis_add(-ntorbd:ntorbd,0:mpolbd)
  NAMELIST /dataadd/borbic_add,borbis_add

  addkes_read=.FALSE.
  OPEN(unit=1,file='add_ddkes2.data',form='formatted',iostat=iostat,action='read')
  IF(iostat.EQ.0) THEN
     WRITE(iout,*) 'File "add_ddkes2.data" found'
     addkes_read=.TRUE.
     borbic_add=0
     borbis_add=0
     READ (1,nml=dataadd)
     CLOSE(1)
     IF(SUM(borbis_add).GT.PREC_B) STELL_ANTISYMMETRIC=.TRUE.
     DO m=0,mpolbd
        DO n=-ntorbd,ntorbd
           IF(m.EQ.0.AND.n.LT.0) CYCLE
           IF(ABS(borbic_add(n,m)).LT.PREC_B.AND.ABS(borbis_add(n,m)).LT.PREC_B) CYCLE
           DO imn=1,mnboz_b
              IF(ixn_b(imn).EQ.n.AND.ixm_b(imn).EQ.m) THEN
                 bmnc_b(imn,:)=bmnc_b(imn,:)+borbic_add(n,m)
                 bmns_b(imn,:)=bmns_b(imn,:)+borbis_add(n,m)
              END IF
           END DO
        END DO
     END DO
  END IF
  
END SUBROUTINE READ_ADDDDKESDATA


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE RESCALE()
!-------------------------------------------------------------------------------------------------
!Update quantities if there is scan in parameters
!-------------------------------------------------------------------------------------------------
  USE GLOBAL
  IMPLICIT NONE
  !Others
  INTEGER imn
  
  DO imn=1,mnboz_b
     IF(ixn_b(imn).EQ.0.AND.ixm_b(imn).EQ.0) THEN
        bmnc_b(imn,:)=bmnc_b(imn,:)/FE
        IF(NEOTRANSP) FB=1./bmnc_b(imn,1)
     END IF
  END DO
  nfp_b=INT(nfp_b*FP+0.1)
  iota_b=iota_b*FI
  bmnc_b=bmnc_b  *FB*FE

  bvco_b=bvco_b  *FB      *FR
  buco_b=buco_b  *FB      *FR
  psip   =psip   *FB*FE   *FR
  chip   =chip*FI*FB*FE   *FR
  torflux=torflux*FB*FE*FE*FR*FR
  rad_a  =rad_a     *FE   *FR
  rad_R  =rad_R           *FR

END SUBROUTINE RESCALE
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE INTERPOLATE_FIELD(s0,booz_read,E_o_mu)  
!-------------------------------------------------------------------------------------------
!Interpolate
!-------------------------------------------------------------------------------------------
  
  USE GLOBAL
  USE KNOSOS_STELLOPT_MOD
  IMPLICIT NONE
  !Input
  LOGICAL booz_read
  REAL*8 s0,E_o_mu
  !Others
  INTEGER imn,n,m,is,is0,is1
  REAL*8 fs0,fs1,f_ext(mboz_b)
  !Time
  CHARACTER*30, PARAMETER :: routine="INTERPOLATE_FIELD"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart
  CALL CPU_TIME(tstart)

  f_ext=0
  IF(booz_read) THEN
     !Find flux-surface and interpolate
     IF(s0.LE.s_b(js_b(1))) THEN
        is0=js_b(1)
        is1=js_b(2)
        f_ext(1)=SQRT(s0/s_b(js_b(1)))
        DO imn=2,mboz_b
           f_ext(imn)=f_ext(imn)*f_ext(imn-1)
        END DO
     ELSE IF(s0.GT.s_b(js_b(jsn))) THEN
        is0=jsn-1
        is1=jsn
     ELSE
        DO is=2,jsn
           IF((s0-s_b(js_b(is-1)))*(s0-s_b(js_b(is))).LE.0) THEN
              is0=js_b(is-1)
              is1=js_b(is)
              EXIT
           END IF
        END DO
     END IF

     fs0=(s_b(is1)-s0)/(s_b(is1)-s_b(is0))
     fs1=(s0-s_b(is0))/(s_b(is1)-s_b(is0))
   ELSE
     is0=1
     is1=1
     fs0=1
     fs1=0
  END IF

  mpolb=mboz_b
  ntorb=nboz_b
  nzperiod=nfp_b

  bzeta =bvco_b(is0)*fs0+bvco_b(is1)*fs1
  btheta=buco_b(is0)*fs0+buco_b(is1)*fs1
  iota=  iota_b(is0)*fs0+iota_b(is1)*fs1
  chip=iota*psip
  IF(ALLOCATED(spol_b)) THEN
     spol  =spol_b(is0)*fs0+spol_b(is1)*fs1
     dspolds=(spol_b(is1)-spol_b(is0))/(s_b(is1)-s_b(is0))
  END IF
  IF(ALLOCATED(phi_b)) phi0=phi_b(is0)*fs0+phi_b(is1)*fs1
     
  !Calculate some radial derivatives
  !The radial coordinate is the toroidal flux over 2pi
  atorflux =ABS(torflux)
  IF(is1.NE.is0) THEN
     dpsi=(s_b(is1)-s_b(is0))*atorflux !the radial coordinate must grow from core to edge
  ELSE IF(NEOTRANSP) THEN
     dpsi=dpsi*rad_a*sqrt(s0)/psip
  ELSE
     dpsi=1
  END IF
  diotadpsi=FS*(iota_b(is1)-iota_b(is0))/dpsi
  dBzdpsi =    (bvco_b(is1)-bvco_b(is0))/dpsi
  dBtdpsi =    (buco_b(is1)-buco_b(is0))/dpsi
  
  !Calculate several other flux-surface quantities
  psip=2*atorflux*sqrt(s0)/rad_a
  aiota=ABS(iota)
  siota=iota/aiota
  iota2=iota*iota
  iBtpBz=(Btheta*iota+Bzeta)
  aiBtpBz=ABS(iBtpBz)
  sgnB=iBtpBz/aiBtpBz
  eps=rad_a*SQRT(s0)/rad_R
  eps32=eps*SQRT(eps)

  borbic=0
  porbis=0
  rorbic=0
  zorbis=0 
  borbis=0
  porbic=0
  rorbis=0
  zorbic=0
  dborbicdpsi=0
  dborbisdpsi=0
  DO imn=1,mnboz_b
     n=ixn_b(imn)
     m=ixm_b(imn)
     IF(s0.GE.s_b(js_b(1)).OR.m.EQ.0) THEN
        borbic(n,m)=bmnc_b(imn,is0)*fs0+bmnc_b(imn,is1)*fs1
         IF(STELL_ANTISYMMETRIC) borbis(n,m)=bmns_b(imn,is0)*fs0+bmns_b(imn,is1)*fs1
        IF(booz_read) THEN
           porbis(n,m)=pmns_b(imn,is0)*fs0+pmns_b(imn,is1)*fs1
           rorbic(n,m)=rmnc_b(imn,is0)*fs0+rmnc_b(imn,is1)*fs1
           zorbis(n,m)=zmns_b(imn,is0)*fs0+zmns_b(imn,is1)*fs1
           IF(STELL_ANTISYMMETRIC) THEN
              porbic(n,m)=pmnc_b(imn,is0)*fs0+pmnc_b(imn,is1)*fs1
              rorbis(n,m)=rmns_b(imn,is0)*fs0+rmns_b(imn,is1)*fs1
              zorbic(n,m)=zmnc_b(imn,is0)*fs0+zmnc_b(imn,is1)*fs1
           END IF
!           IF(is0.LT.2.OR.is1.GT.ns_b-1) THEN
           dborbicdpsi(n,m)=(bmnc_b(imn,is1)-bmnc_b(imn,is0))/dpsi
           dborbisdpsi(n,m)=(bmns_b(imn,is1)-bmns_b(imn,is0))/dpsi
!           ELSE
!              dborbicdpsi(n,m)=(fs0*(bmnc_b(imn,is1)-bmnc_b(imn,is0-1))+fs1*(bmnc_b(imn,is1+1)-bmnc_b(imn,is0)))/(4*dpsi)
!              dborbisdpsi(n,m)=(fs0*(bmns_b(imn,is1)-bmns_b(imn,is0-1))+fs1*(bmns_b(imn,is1+1)-bmns_b(imn,is0)))/(4*dpsi)
!           END IF
        END IF
     ELSE
        borbic(n,m)=bmnc_b(imn,is0)*f_ext(m)
!        IF(ABS(borbic(n,m)).GT.10) WRITE(iout,*) n,m,borbic(n,m)
        IF(STELL_ANTISYMMETRIC) borbis(n,m)=bmns_b(imn,is0)*f_ext(m)
        IF(booz_read) THEN
           porbis(n,m)=pmns_b(imn,is0)*f_ext(m)
           rorbic(n,m)=rmnc_b(imn,is0)*f_ext(m)
           zorbis(n,m)=zmns_b(imn,is0)*f_ext(m)
           IF(STELL_ANTISYMMETRIC) THEN
              porbic(n,m)=pmnc_b(imn,is0)*f_ext(m)
              rorbis(n,m)=rmns_b(imn,is0)*f_ext(m)
              zorbic(n,m)=zmnc_b(imn,is0)*f_ext(m)
           END IF
           dborbicdpsi(n,m)=bmnc_b(imn,is0)*f_ext(m)*(m/2.)/s0
           dborbisdpsi(n,m)=bmns_b(imn,is0)*f_ext(m)*(m/2.)/s0
        END IF
     END IF
  END DO

  IF(ALLOCATED(phi_b)) THEN
     borbic(0,0)=borbic(0,0)+E_o_mu*ZFI*phi0/EFI
     dborbicdpsi(0,0)=dborbicdpsi(0,0)+E_o_mu*ZFI*(phi_b(is1)-phi_b(is0))/(EFI*dpsi)
  END IF

  IF(.NOT.KN_STELLOPT(9)) THEN
     borbic0=borbic
     borbis0=borbis
     CALL FILL_NM()
  END IF

!  IF(s0.GE.s_b(js_b(1)).AND.ALLOCATED(Bmax_b)) THEN
  IF(ALLOCATED(Bmax_b)) THEN
     IF(Bmax_b(is1).GT.ALMOST_ZERO) THEN
        Bmax  =Bmax_b(is0)*fs0+Bmax_b(is1)*fs1
        Bmin  =Bmin_b(is0)*fs0+Bmin_b(is1)*fs1
     ELSE
        CALL FIND_BMIN_BMAX(-1,MAL,MAL,.TRUE.)
     END IF
  ELSE
     CALL FIND_BMIN_BMAX(-1,MAL,MAL,.TRUE.)
  END IF

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE INTERPOLATE_FIELD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE JACOBIAN_AND_3DGRID(s0)

!-------------------------------------------------------------------------------------------
!Interpolate
!-------------------------------------------------------------------------------------------
  
  USE GLOBAL
  USE KNOSOS_STELLOPT_MOD
  IMPLICIT NONE
  !Input
  REAL*8 s0
  !Others
  LOGICAL LeftHanded
  CHARACTER*100 serr
  INTEGER it,jt,is,is0,is1,imn
  REAL*8 fs0,fs1,x1(3,MAL,MAL),x2(3,MAL,MAL),x3(3,MAL,MAL),Bzt(3,MAL,MAL)

  IF(s0.LE.s_b(js_b(1))) THEN
     is0=js_b(1)
     is1=js_b(2)
  ELSE IF(s0.GT.s_b(js_b(jsn))) THEN
     is0=jsn-1
     is1=jsn
  ELSE
     DO is=2,jsn
        IF((s0-s_b(js_b(is-1)))*(s0-s_b(js_b(is))).LE.0) THEN
           is0=js_b(is-1)
           is1=js_b(is)
           EXIT
        END IF
     END DO
  END IF
  fs0=(s_b(is1)-s0)/(s_b(is1)-s_b(is0))
  fs1=(s0-s_b(is0))/(s_b(is1)-s_b(is0))
  !Fill borbic,dborbicdps,rorbic,etc
  borbic=0
  porbis=0
  rorbic=0
  zorbis=0  
  borbis=0
  porbic=0
  rorbis=0
  zorbic=0 
  DO imn=1,mnboz_b
     borbic(ixn_b(imn),ixm_b(imn))=bmnc_b(imn,is0)
     porbis(ixn_b(imn),ixm_b(imn))=pmns_b(imn,is0)
     rorbic(ixn_b(imn),ixm_b(imn))=rmnc_b(imn,is0)
     zorbis(ixn_b(imn),ixm_b(imn))=zmns_b(imn,is0)
     IF(STELL_ANTISYMMETRIC) THEN
        borbis(ixn_b(imn),ixm_b(imn))=bmns_b(imn,is0)
        porbic(ixn_b(imn),ixm_b(imn))=pmnc_b(imn,is0)
        rorbis(ixn_b(imn),ixm_b(imn))=rmns_b(imn,is0)
        zorbic(ixn_b(imn),ixm_b(imn))=zmnc_b(imn,is0)
     END IF
  END DO
  CALL FILL_3DGRID(MAL,MAL,s_b(is0),x1(1,:,:),x2(1,:,:),x3(1,:,:),Bzt(1,:,:),.FALSE.)
  borbic=0
  porbis=0
  rorbic=0
  zorbis=0  
  borbis=0
  porbic=0
  rorbis=0
  zorbic=0  
  DO imn=1,mnboz_b
     borbic(ixn_b(imn),ixm_b(imn))=bmnc_b(imn,is1)
     porbis(ixn_b(imn),ixm_b(imn))=pmns_b(imn,is1)
     rorbic(ixn_b(imn),ixm_b(imn))=rmnc_b(imn,is1)
     zorbis(ixn_b(imn),ixm_b(imn))=zmns_b(imn,is1)
     IF(STELL_ANTISYMMETRIC) THEN
        borbis(ixn_b(imn),ixm_b(imn))=bmns_b(imn,is1)
        porbic(ixn_b(imn),ixm_b(imn))=pmnc_b(imn,is1)
        rorbis(ixn_b(imn),ixm_b(imn))=rmns_b(imn,is1)
        zorbic(ixn_b(imn),ixm_b(imn))=zmnc_b(imn,is1)
     END IF
  END DO
  CALL FILL_3DGRID(MAL,MAL,s_b(is1),x1(3,:,:),x2(3,:,:),x3(3,:,:),Bzt(3,:,:),.FALSE.)
  borbic=0
  porbis=0
  rorbic=0
  zorbis=0 
  borbis=0
  porbic=0
  rorbis=0
  zorbic=0
  DO imn=1,mnboz_b
     borbic(ixn_b(imn),ixm_b(imn))=bmnc_b(imn,is0)*fs0+bmnc_b(imn,is1)*fs1
     porbis(ixn_b(imn),ixm_b(imn))=pmns_b(imn,is0)*fs0+pmns_b(imn,is1)*fs1
     rorbic(ixn_b(imn),ixm_b(imn))=rmnc_b(imn,is0)*fs0+rmnc_b(imn,is1)*fs1
     zorbis(ixn_b(imn),ixm_b(imn))=zmns_b(imn,is0)*fs0+zmns_b(imn,is1)*fs1
     IF(STELL_ANTISYMMETRIC) THEN
        borbis(ixn_b(imn),ixm_b(imn))=bmns_b(imn,is0)*fs0+bmns_b(imn,is1)*fs1
        porbic(ixn_b(imn),ixm_b(imn))=pmnc_b(imn,is0)*fs0+pmnc_b(imn,is1)*fs1
        rorbis(ixn_b(imn),ixm_b(imn))=rmns_b(imn,is0)*fs0+rmns_b(imn,is1)*fs1
        zorbic(ixn_b(imn),ixm_b(imn))=zmnc_b(imn,is0)*fs0+zmnc_b(imn,is1)*fs1
     END IF
  END DO
  CALL FILL_3DGRID(MAL,MAL,s0,x1(2,:,:),x2(2,:,:),x3(2,:,:),Bzt(2,:,:),.TRUE.)
  IF(.NOT.ALLOCATED(posx)) ALLOCATE( posx(MAL,MAL), posy(MAL,MAL), posz(MAL,MAL))
  DO it=1,MAL  !from lhs to rhs
     jt=MAL-it+1
     IF(it.EQ.1) jt=1
     posx(:,it)=x1(2,:,jt)
     posy(:,it)=x2(2,:,jt)
     posz(:,it)=x3(2,:,jt)
  END DO
  !Check if the coordinate system is left-handed, as expected
  CALL CHECK_JACSIGN(MAL,MAL,x1,x2,x3,Bzt(2,:,:),LeftHanded) 
  IF(LeftHanded.AND..NOT.(SATAKE.OR.KNOSOS_STELLOPT.OR.dpsi/atorflux.GT.0.05)) THEN !JL
     serr="The magnetic field is now in left-handed coordinates" 
     CALL END_ALL(serr,.FALSE.)
  END IF
  IF(.NOT.KNOSOS_STELLOPT) THEN
     CALL FIND_3DPOINTS(MAL,MAL,s0,x1(1,:,:),x2(1,:,:),x3(1,:,:))
     IF(.NOT.ALLOCATED(zoomx)) ALLOCATE(zoomx(MAL,MAL),zoomy(MAL,MAL),zoomz(MAL,MAL))
     zoomx=x1(1,:,:)
     zoomy=x2(1,:,:)
     zoomz=x3(1,:,:)
  END IF

END SUBROUTINE JACOBIAN_AND_3DGRID


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE FILL_NM()

!-------------------------------------------------------------------------------------------------
!Fill arrays (Bmn, etc) with only one index
!Keep only (n,m) such that some quantity f_{mn} is large enough
!-------------------------------------------------------------------------------------------------

  USE GLOBAL
  IMPLICIT NONE
  !Others
  INTEGER n,m,nm

  nm=0
  bnmc=0
  bnmc0=0
  dbnmcdpsi=0
  bnms=0
  bnms0=0
  dbnmsdpsi=0
!  phnmc=0
!  phnms=0
!  pnm=0

  !Loop in modes; keeps only modes larger than PREC_B
  DO m=0,mpolb
     DO n=-ntorb,ntorb
        IF(m.EQ.0.AND.n.LT.0) CYCLE
        IF(.NOT.QN.AND.(ABS(borbic0(n,m)).LT.PREC_B).AND.(ABS(borbic(n,m)).LT.PREC_B)&
             &    .AND.(ABS(borbis0(n,m)).LT.PREC_B).AND.(ABS(borbis(n,m)).LT.PREC_B)) CYCLE
!        IF(.NOT.QN.AND.(ABS(borbic0(n,m)).LT.PREC_B).AND.(ABS(borbic(n,m)).LT.PREC_B)) CYCLE
!             & .AND.(ABS(phorbicc(n,m)).LT.PREC_B).AND.(ABS(phorbics(n,m)).LT.PREC_B)) CYCLE
        IF(NEQ2.AND.n.NE.1.AND.n.NE.0) CYCLE
        nm=nm+1
        bnmc(nm)     =borbic(n,m)     
        bnmc0(nm)    =borbic0(n,m)
!        WRITE(iout,*) 'bmn1',bnmc(nm)        -bnmc0(nm)        
        dbnmcdpsi(nm)=dborbicdpsi(n,m)
        IF(STELL_ANTISYMMETRIC) THEN
           bnms(nm)     =borbis(n,m)     
           bnms0(nm)    =borbis0(n,m)     
           dbnmsdpsi(nm)=dborbisdpsi(n,m)
        END IF
!        phnmc(nm)=phorbicc(n,m)
!        phnms(nm)=phorbics(n,m)
!        pnm(nm)=porbis(n,m)
!        IF(NTV) THEN
!           np(nm)=m
!           mp(nm)=sgnhel*n*nzperiod           
!        ELSE
           np(nm)=n
           mp(nm)=m
           !        END IF
!           WRITE(iout,*) 'nnm',nnm,ntorb,mpolb
     END DO
  END DO
  
  Nnm=nm                          !Nnm is the total number of modes
  bnmc1(1:Nnm)=bnmc(1:Nnm)-bnmc0(1:Nnm) !B1=B-B0, remember that B1<<B0
  bnms1(1:Nnm)=bnms(1:Nnm)-bnms0(1:Nnm)
  IF(QN.OR.TRACE_IMP) THEN
     Nnmp=2*Nnm
     DO nm=1,Nnmp
        ext_np(nm)=np(MOD(nm-1,Nnm)+1)
        ext_mp(nm)=mp(MOD(nm-1,Nnm)+1)
     END DO
     WRITE(iout,'(" Electrostatic potential characterized by ",&
          & I5," harmonics")') Nnmp
  ELSE
     Nnmp=1
  END IF

  IF(NTV) nzperiod=1
!  IF(SATAKE) NTV=.TRUE.

END SUBROUTINE FILL_NM


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALCB(z,t,flag,flagB1,&
     & B_0,dBdz_0,dBdt_0,dBdpsi,hBpp,&
     & B_1,dBdz_1,dBdt_1,Phi_1,dPhdz,dPhdt,vde)

!-------------------------------------------------------------------------------------------------
!Calculate magnetic field and derivatives at angular position (z,t)
!-flag.EQ.0: calculate only B
!-flag.EQ.1: calculate only first derivatives of B, dBdz and dBdt
!-flag.EQ.2: calculate B and its first derivatives, B and dBdz and dBdt
!-flag.EQ.3: calculate B, its first (including radial) and its second derivatives hBpp
!-flag.EQ.4: calculate the angle-dependent part of the ExB radial drift vde
!------
!-IF(.NOT.flagB1) calculate only B_0 (usually B_1=0, so B=B_0) and its derivatives dBdz_0 and dBdt_0
!-IF(flagB1) calculate also B_1 and their derivatives dBdz_1 and dBdt_1
!------
!Phi_1, dPhdz, and dPhdt not implemented (see older versions)
!-------------------------------------------------------------------------------------------------

  USE GLOBAL
  IMPLICIT NONE
  !Input
  LOGICAL flagB1 !if true, calculate quantities related to B1
  INTEGER flag
  REAL*8 z,t
  !Output
  REAL*8 B_0,dBdz_0,dBdt_0,dBdpsi,hBpp
  REAL*8 B_1,dBdz_1,dBdt_1
  REAL*8 Phi_1,dPhdz,dPhdt,vde(Nnmp)
  !Others
  INTEGER nm,nm2
  REAL*8 cosinex,sinex,d2Bdz2,d2Bdt2,d2Bdzt

  !Set to zero
  B_0=0
  dBdz_0=0
  d2Bdz2=0
  dBdt_0=0
  dBdpsi=0
  d2Bdt2=0
  d2Bdzt=0
  B_1=0
  dBdz_1=0
  dBdt_1=0
  dPhdz=0
  dPhdt=0
  Phi_1=0

!!$  IF(flag.NE.0) THEN
!!$     DO nm=1,Nnm
!!$        sinex=SIN(mp(nm)*t+nzperiod*np(nm)*z)
!!$        dBdz_0=dBdz_0-bnmc0(nm)*sinex*nzperiod*np(nm)
!!$        dBdt_0=dBdt_0-bnmc0(nm)*sinex*mp(nm)
!!$        IF(flagB1) THEN
!!$           dBdz_1=dBdz_1-bnmc1(nm)*sinex*nzperiod*np(nm)
!!$           dBdt_1=dBdt_1-bnmc1(nm)*sinex*mp(nm)
!!$        END IF       
!!$        IF(flag.EQ.4) vde(nm)=-sinex*(Btheta*np(nm)*nzperiod-Bzeta*mp(nm))/aiBtpBz
!!$     END DO
!!$  END IF
!!$  IF(flag.EQ.1) RETURN
!!$
!!$  DO nm=1,Nnm
!!$     cosinex=COS(mp(nm)*t+nzperiod*np(nm)*z)
!!$     B_0=B_0+bnmc0(nm)*cosinex
!!$     IF(flagB1) B_1=B_1+bnmc1(nm)*cosinex
!!$     IF(TANG_VM.AND.flag.GT.0) dBdpsi=dBdpsi+dbnmcdpsi(nm)*cosinex
!!$     IF(flag.EQ.3) THEN
!!$        d2Bdz2=d2Bdz2-bnmc0(nm)*cosinex*nzperiod*np(nm)*nzperiod*np(nm)
!!$        d2Bdzt=d2Bdzt-bnmc0(nm)*cosinex*nzperiod*np(nm)*mp(nm)
!!$        d2Bdt2=d2Bdt2-bnmc0(nm)*cosinex*mp(nm)*mp(nm)
!!$     END IF
!!$     IF(flag.EQ.4.AND..NOT.ONLY_PHI1) vde(nm+Nnm)=cosinex*(Btheta*np(nm)*nzperiod-Bzeta*mp(nm))/aiBtpBz
!!$  END DO
!!$  IF(flag.EQ.3) hBpp=d2Bdz2/2.


  DO nm=1,Nnm
     sinex=  -100
     cosinex=-100

     IF(STELL_ANTISYMMETRIC) THEN

        sinex  =SIN(mp(nm)*t+nzperiod*np(nm)*z)
        cosinex=COS(mp(nm)*t+nzperiod*np(nm)*z)
        IF(flag.NE.0) THEN
           dBdz_0=dBdz_0-bnmc0(nm)*sinex*nzperiod*np(nm)+bnms0(nm)*cosinex*nzperiod*np(nm)
           dBdt_0=dBdt_0-bnmc0(nm)*sinex*mp(nm)         +bnms0(nm)*cosinex*mp(nm)
           IF(flagB1) THEN
              dBdz_1=dBdz_1-bnmc1(nm)*sinex*nzperiod*np(nm)+bnms1(nm)*cosinex*nzperiod*np(nm)
              dBdt_1=dBdt_1-bnmc1(nm)*sinex*mp(nm)         +bnms1(nm)*cosinex*mp(nm)
           END IF
           IF(flag.EQ.4.AND..NOT.ONLY_PHI1) &
                & vde(nm)=-sinex*(Btheta*np(nm)*nzperiod-Bzeta*mp(nm))/aiBtpBz
        END IF
        IF(flag.NE.1) THEN
           B_0=B_0+bnmc0(nm)*cosinex+bnms0(nm)*sinex
           IF(flagB1) B_1=B_1+bnmc1(nm)*cosinex+bnms1(nm)*sinex
           IF(TANG_VM.AND.flag.GT.0) dBdpsi=dBdpsi+dbnmcdpsi(nm)*cosinex+dbnmsdpsi(nm)*sinex
           IF(flag.EQ.3) THEN
              d2Bdz2=d2Bdz2-bnmc0(nm)*cosinex*nzperiod*np(nm)*nzperiod*np(nm)&
                         & -bnms0(nm)*  sinex*nzperiod*np(nm)*nzperiod*np(nm)
              d2Bdzt=d2Bdzt-bnmc0(nm)*cosinex*nzperiod*np(nm)*mp(nm)&
                         & -bnms0(nm)*  sinex*nzperiod*np(nm)*mp(nm)
              d2Bdt2=d2Bdt2-bnmc0(nm)*cosinex*mp(nm)*mp(nm)-bnms0(nm)*sinex*mp(nm)*mp(nm)
           END IF
           IF(flag.EQ.4.AND..NOT.ONLY_PHI1) &
                & vde(nm+Nnm)=cosinex*(Btheta*np(nm)*nzperiod-Bzeta*mp(nm))/aiBtpBz           
        END IF

     ELSE

        IF(flag.NE.0) THEN
           sinex=SIN(mp(nm)*t+nzperiod*np(nm)*z)
           dBdz_0=dBdz_0-bnmc0(nm)*sinex*nzperiod*np(nm)
           dBdt_0=dBdt_0-bnmc0(nm)*sinex*mp(nm)
           IF(flagB1) THEN
              dBdz_1=dBdz_1-bnmc1(nm)*sinex*nzperiod*np(nm)
              dBdt_1=dBdt_1-bnmc1(nm)*sinex*mp(nm)
           END IF
           IF(flag.EQ.4.AND..NOT.ONLY_PHI1) &
                & vde(nm)=-sinex*(Btheta*np(nm)*nzperiod-Bzeta*mp(nm))/aiBtpBz
        END IF

        IF(flag.NE.1) THEN
           cosinex=COS(mp(nm)*t+nzperiod*np(nm)*z)
           B_0=B_0+bnmc0(nm)*cosinex
           IF(flagB1) B_1=B_1+bnmc1(nm)*cosinex
           IF(TANG_VM.AND.flag.GT.0) dBdpsi=dBdpsi+dbnmcdpsi(nm)*cosinex
           IF(flag.EQ.3) THEN
              d2Bdz2=d2Bdz2-bnmc0(nm)*cosinex*nzperiod*np(nm)*nzperiod*np(nm)
              d2Bdzt=d2Bdzt-bnmc0(nm)*cosinex*nzperiod*np(nm)*mp(nm)
              d2Bdt2=d2Bdt2-bnmc0(nm)*cosinex*mp(nm)*mp(nm)
           END IF
           IF(flag.EQ.4.AND..NOT.ONLY_PHI1) &
                & vde(nm+Nnm)=cosinex*(Btheta*np(nm)*nzperiod-Bzeta*mp(nm))/aiBtpBz           
        END IF

     END IF

     IF(PHI1_READ) THEN
        nm2=Nnm+nm
        IF(sinex.LT.-10)     sinex=SIN(mp(nm)*t+nzperiod*np(nm)*z)
        IF(cosinex.LT.-10) cosinex=COS(mp(nm)*t+nzperiod*np(nm)*z)
        IF(flag.NE.0) THEN
           dBdz_0=dBdz_0+bnmc0(nm2)*cosinex*nzperiod*np(nm)
           dBdt_0=dBdt_0+bnmc0(nm2)*cosinex*mp(nm)
        END IF
        IF(flag.NE.1) THEN
           B_0=B_0+bnmc0(nm2)*sinex
           IF(flag.EQ.3) THEN
              d2Bdz2=d2Bdz2-bnmc0(nm2)*sinex*nzperiod*np(nm)*nzperiod*np(nm)
              d2Bdzt=d2Bdzt-bnmc0(nm2)*sinex*nzperiod*np(nm)*mp(nm)
              d2Bdt2=d2Bdt2-bnmc0(nm2)*sinex*mp(nm)*mp(nm)
           END IF
        END IF
     END IF
     
  END DO

  IF(flag.EQ.3) hBpp=d2Bdz2/2.

END SUBROUTINE CALCB


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CHECK_JACSIGN(nz,nt,x1,x2,x3,Bzt,LeftHanded)

!-------------------------------------------------------------------------------------------------
!With x1, x2 and x3 evaluated in a 3 x nz x nt grid in boozer coordinates, evaluates sign of the
!jacobian
!-------------------------------------------------------------------------------------------------

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER nz,nt
  REAL*8 x1(3,nz,nt),x2(3,nz,nt),x3(3,nz,nt),Bzt(nz,nt)
  !Output
  LOGICAL LeftHanded
  !Others
  INTEGER iz,izp1,izm1,it,itp1,itm1
  REAL*8 dx1dpsi,dx1dt,dx1dz,dx2dpsi,dx2dt,dx2dz,dx3dpsi,dx3dt,dx3dz,sqrtg,k,k2,dz,dt,denom

  k =0
  k2=0
  denom=0 
  etet=0
  !Calculate derivatives
  DO iz=1,nz
     IF(iz.GT.2.AND..NOT.DEBUG.AND..NOT.ONLY_PHI1) CYCLE
     dz=1.0
     IF(iz.EQ.nz) THEN
        izp1=nz
        dz=0.5
     ELSE
        izp1=iz+1
     END IF
     IF(iz.EQ.1) THEN
        izm1=1
        dz=0.5
     ELSE
        izm1=iz-1
     END IF
     DO it=1,nt
        IF(it.GT.2.AND..NOT.DEBUG.AND..NOT.ONLY_PHI1) CYCLE
        dt=1.0
        IF(it.EQ.nt) THEN
           itp1=nt
           dt=0.5
        ELSE
           itp1=it+1
        END IF
        IF(it.EQ.1) THEN
           itm1=1
           dt=0.5
        ELSE
           itm1=it-1
        END IF
        dx1dpsi=(x1(3,iz,it)  -x1(1,iz,it))/dpsi
        dx2dpsi=(x2(3,iz,it)  -x2(1,iz,it))/dpsi
        dx3dpsi=(x3(3,iz,it)  -x3(1,iz,it))/dpsi
        dx1dz=(x1(2,izp1,it)-x1(2,izm1,it))/dz
        dx2dz=(x2(2,izp1,it)-x2(2,izm1,it))/dz
        dx3dz=(x3(2,izp1,it)-x3(2,izm1,it))/dz
        dx1dt=(x1(2,iz,itp1)-x1(2,iz,itm1))/dt
        dx2dt=(x2(2,iz,itp1)-x2(2,iz,itm1))/dt
        dx3dt=(x3(2,iz,itp1)-x3(2,iz,itm1))/dt
        sqrtg=dx1dpsi*(dx2dt*dx3dz-dx3dt*dx2dz) &
           & -dx2dpsi*(dx1dt*dx3dz-dx3dt*dx1dz) &
           & +dx3dpsi*(dx1dt*dx2dz-dx2dt*dx1dz)
        k =k +sqrtg*Bzt(iz,it)*Bzt(iz,it)
        k2=k2+sqrtg*sqrtg*Bzt(iz,it)*Bzt(iz,it)*Bzt(iz,it)*Bzt(iz,it)
        LeftHanded=(sqrtg.LT.0)
!        IF(DEBUG) 
        absnablar(iz,it)=-SQRT((dx2dt*dx3dz-dx3dt*dx2dz)*(dx2dt*dx3dz-dx3dt*dx2dz)+&
             &                (dx1dt*dx3dz-dx3dt*dx1dz)*(dx1dt*dx3dz-dx3dt*dx1dz)+&
             &                (dx1dt*dx2dz-dx2dt*dx1dz)*(dx1dt*dx2dz-dx2dt*dx1dz))/ABS(psip)/sqrtg
        IF(DEBUG) WRITE(1400+myrank,'(6(1pe13.5),L)') x1(2,iz,it),x2(2,iz,it),x3(2,iz,it),&
             & sqrtg,Bzt(iz,it),absnablar(iz,it),LeftHanded
!        IF(.NOT.LeftHanded) RETURN
        IF(iz.EQ.2.AND.it.EQ.2.AND..NOT.DEBUG.AND..NOT.ONLY_PHI1) RETURN
        etet=etet+(dx1dt*dx1dt+dx2dt*dx2dt+dx3dt*dx3dt)/(Bzt(iz,it)*Bzt(iz,it))
        denom=denom+1/(Bzt(iz,it)*Bzt(iz,it))
     END DO
  END DO
  !Check if there is too much deviation from sqrt(g)*B^2=constant
  k=k/(nz*nt)
  k2=SQRT((k2-k*k*(nz*nt))/(nz*nt-1.0))  
!  IF(k2.GT.0.5*ABS(k)) LeftHanded=.FALSE. 
  IF(k2.LT.0) LeftHanded=.FALSE. 
  IF(DEBUG) WRITE(1400+myrank,'(5(1pe13.5),L)') k,k2
  etet=etet/denom

END SUBROUTINE CHECK_JACSIGN
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALC_XYZ(s,z,t,x1,x2,x3,Bzt,flag_plot)

!-------------------------------------------------------------------------------------------------
!Calculate cartesian coordinates (x1,x2,x3) of point labelled by (s,t,z) in flux-coordinates
!-------------------------------------------------------------------------------------------------

  USE GLOBAL
  IMPLICIT NONE
  !Input
  LOGICAL flag_plot
  REAL*8 s,z,t
  !Output
  REAL*8 x1,x2,x3,Bzt
  !Others
  INTEGER n,m
  REAL*8 arg,cosine,sine,R,p

  x3=0
  R =0
  p =z
  Bzt=0
  DO m=0,mpolb
     DO n=-ntorb,ntorb 
        arg=m*t-n*nzperiod*z
        cosine=COS(arg)
        sine  =SIN(arg)
        Bzt=Bzt+borbic(n,m)*cosine
        x3=x3+zorbis(n,m)*  sine
        R =R +rorbic(n,m)*cosine
        p =p-(porbis(n,m)*  sine)
        IF(STELL_ANTISYMMETRIC) THEN
           Bzt=Bzt+borbis(n,m)*sine
           x3=x3+zorbic(n,m)*cosine
           R =R +rorbis(n,m)*  sine
           p =p-(porbic(n,m)*cosine)
        END IF
     END DO
  END DO
  x1=R*COS(p)
  x2=R*SIN(p)
  
  IF(DEBUG.AND.flag_plot) WRITE(1600+myrank,'(9(1pe13.5))') x1,x2,x3,R,p,Bzt,s,z,t

END SUBROUTINE CALC_XYZ


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


REAL*8 FUNCTION MODANG(ang_in,max)

!-------------------------------------------------------------------------------------------------
!Calculate MOD(ang_in,max) making sure that the result lies in [0,max)
!-------------------------------------------------------------------------------------------------

  USE GLOBAL
  IMPLICIT NONE

  REAL*8 ang_in
  REAL*8 max
  
  modang=MOD(ang_in-ALMOST_ZERO,max)
  IF(modang.LT.0) modang=modang+max
  modang=modang+ALMOST_ZERO

  RETURN 

END FUNCTION MODANG


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


REAL*8 FUNCTION MODANG2(ang_in,max)

!-------------------------------------------------------------------------------------------------
!Calculate MOD(ang_in,max) making sure that the result lies in [0,max)
!-------------------------------------------------------------------------------------------------

  USE GLOBAL
  IMPLICIT NONE

  REAL*8 ang_in
  REAL*8 max
  
  modang2=MOD(ang_in,max)
  IF(modang2.LT.0) modang2=modang2+max

  RETURN 

END FUNCTION MODANG2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


REAL*8 FUNCTION MODANGHEL(z_in,offset,max)
!-------------------------------------------------------------------------------------------------
!xxxxx TO BE REMOVED
!-------------------------------------------------------------------------------------------------
  
  USE GLOBAL
  IMPLICIT NONE
  !Input
  REAL*8 z_in,offset,max
  !Others
  INTEGER, PARAMETER :: nmax=100
  REAL*8 ang_hel,angtemp

  IF(nzperiod.EQ.1) THEN
     ang_hel=offset
  ELSE
     ang_hel=0
  END IF
  angtemp=z_in+nmax*max
  DO WHILE (.TRUE.)
     IF(angtemp.LT.ang_hel) EXIT
     angtemp=angtemp-max
  END DO
 modanghel=angtemp+max

 RETURN

END FUNCTION MODANGHEL



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE FILL_3DGRID(nz,nt,s,x1,x2,x3,Bzt,flag)
  
!-------------------------------------------------------------------------------------------------
!For nz x nt points uniformly distributed in the Boozer angles at flux-surface s, calculate points
!(x1,x2,x3) and plot if flag
!-------------------------------------------------------------------------------------------------

  USE GLOBAL
  IMPLICIT NONE
  !Input
  LOGICAL flag
  INTEGER nz,nt
  REAL*8 s
  !Output
  REAL*8 x1(nz,nt),x2(nz,nt),x3(nz,nt),Bzt(nz,nt)
  !Others
  INTEGER iz,it
  REAL*8 dz,dt,zeta(nz),theta(nt)
  !Time
  CHARACTER*30, PARAMETER :: routine="FILL_3DGRID"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  dz=TWOPI/nz/nzperiod
  dt=TWOPI/nt
  DO iz=1,nz
     zeta(iz)=(iz-1.)*dz
  END DO
  DO it=1,nt
     theta(it)=(it-1.)*dt
  END DO
  DO iz=1,nz
     IF(iz.GT.3.AND..NOT.DEBUG.AND..NOT.ONLY_PHI1) CYCLE
     DO it=1,nt
        IF(it.GT.3.AND..NOT.DEBUG.AND..NOT.ONLY_PHI1) CYCLE
        CALL CALC_XYZ(s,zeta(iz),theta(it),x1(iz,it),x2(iz,it),x3(iz,it),Bzt(iz,it),flag)
     END DO
  END DO

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE FILL_3DGRID


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE FIND_3DPOINTS(nz,nt,s,x1,x2,x3)
  
!-------------------------------------------------------------------------------------------------
!For nz x nt points uniformly distributed in the Boozer angles at flux-surface s, find points
!in 3D
!-------------------------------------------------------------------------------------------------

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER nz,nt
  REAL*8 s
  !Output
  REAL*8 x1(nz,nt),x2(nz,nt),x3(nz,nt)
  !Others
  CHARACTER*100 filename
  INTEGER iarray,jarray,iparray,jparray,iz,it,iz0,it0,iostat
  REAL*8 dz,dt,zeta(nz),theta(nt),dist,newdist,Bzt(nz,nt)
  REAL*8 x10(narray,nparray),x20(narray,nparray),x30(narray,nparray),rho0,s0(narray,nparray),dummy
  REAL*8 absnablapsioB2(nz,nt)
  !Time
  CHARACTER*30, PARAMETER :: routine="FILL_3DPOINTS"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  dz=dzeta0_DR/nz
  dt=dtheta0_DR/nt
  DO iz=1,nz
     zeta(iz)=zeta0_DR+(iz-nz/2.)*dz
  END DO
  DO it=1,nt
     theta(it)=theta0_DR+(it-nt/2.)*dt
  END DO
  DO iz=1,nz
     DO it=1,nt
        CALL CALC_XYZ(s,zeta(iz),theta(it),x1(iz,it),x2(iz,it),x3(iz,it),Bzt(iz,it),DEBUG)
     END DO
  END DO

!  array=0
  zetaDR=-1.0
  thetaDR=-1.0
  DO iarray=1,narray
     iz0=0
     it0=0
     jarray=0
     jparray=0
     dist=1e10
     IF(iarray.EQ.1) filename="exp1.dat"
     IF(iarray.EQ.2) filename="exp2.dat"
     IF(iarray.GT.2) filename="expn.dat"
     OPEN(unit=1,file=filename,iostat=iostat,action='read')
     IF(iostat.EQ.0) THEN
        DR_READ=.TRUE.
        DO iparray=1,nparray 
           READ(1,*,iostat=iostat) dummy,&
                & x10(iarray,iparray),x20(iarray,iparray),x30(iarray,iparray),rho0,dummy           
           IF(iostat.NE.0) EXIT
           s0(iarray,iparray)=rho0*rho0
           x10(iarray,iparray)=x10(iarray,iparray)/1e2
           x20(iarray,iparray)=x20(iarray,iparray)/1e2
           x30(iarray,iparray)=x30(iarray,iparray)/1e2
!           IF(ABS(s0(iarray,iparray)-s).GT.0.0035) CYCLE
           DO iz=1,nz
              DO it=1,nt
                 newdist=(x1(iz,it)-x10(iarray,iparray))*(x1(iz,it)-x10(iarray,iparray))&
                      & +(x2(iz,it)-x20(iarray,iparray))*(x2(iz,it)-x20(iarray,iparray))&
                      & +(x3(iz,it)-x30(iarray,iparray))*(x3(iz,it)-x30(iarray,iparray))
                 IF(newdist.LT.dist) THEN
                    iz0=iz
                    it0=it
                    dist=newdist
                    jarray =iarray
                    jparray=iparray
                 END IF
              END DO
           END DO
        END DO
        CLOSE(1)
        IF(iz0.NE.0) THEN
           WRITE(1500+myrank,'(9(1pe13.5),4I4)') s0(jarray,jparray),&
                & x10(jarray,jparray),x20(jarray,jparray),x30(jarray,jparray),& 
                & s,x1(iz0,it0),x2(iz0,it0),x3(iz0,it0),dist,jarray,jparray,iz0,it0
           !        array(iz0,it0)=jarray
           zetaDR(jarray,jparray)=zeta(iz0)
           thetaDR(jarray,jparray)=TWOPI-theta(it0) !From lhs to rhs
        END IF
     END IF
  END DO

  IF(DR_READ) THEN
     IF(numprocs.EQ.1) filename="Er.map"
     IF(numprocs.GT.1) WRITE(filename,'("Er.map.",I2.2)') myrank
     OPEN(unit=4500+myrank,file=filename,form='formatted',action='write',iostat=iostat)
     WRITE(4500+myrank,'("s  \zeta_{Boozer}  \theta_{Boozer}(right-handed) &
          & varphi_1[V] varphi_1(fit[V] T_i [eV] dvarphi_1/ds [V] &
          & -dvarphi_0/dr[V/m]  -dvarphi_0/dr*|absnablar| [V/m] -d(varphi_0+varphi_1)/dr*|absnablar| [V/m]")')

     CALL CALC_ABSNABLAPSI(MAL,zeta,TWOPI-theta,absnablapsioB2)
     IF(.NOT.ALLOCATED(zoomdr)) ALLOCATE(zoomdr(MAL,MAL))     
     zoomdr=SQRT(absnablapsioB2*Bzt*Bzt)/psip
  END IF


  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)
     
END SUBROUTINE FIND_3DPOINTS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!$
!!$SUBROUTINE FILL_3DGRID_OLD(nz,nt,s,x1,x2,x3,Bzt,flag_plot)
!!$
!!$!-------------------------------------------------------------------------------------------------
!!$!For nz x nt points uniformly distributed in the Boozer angles at flux-surface s, calculate points
!!$!(x1,x2,x3) and plot if flag_plot
!!$!-------------------------------------------------------------------------------------------------
!!$
!!$  USE GLOBAL
!!$  IMPLICIT NONE
!!$  !Inpubamt
!!$  LOGICAL flag_plot
!!$  INTEGER nz,nt
!!$  REAL*8 s
!!$  !Output
!!$  REAL*8 x1(nz,nt),x2(nz,nt),x3(nz,nt),Bzt(nz,nt)
!!$  !Others
!!$  INTEGER iz,it
!!$  REAL*8 dz,dt,zeta(nz),theta(nt)
!!$  !Others Others
!!$  INTEGER iz0,it0
!!$  REAL*8 dist,newdist,x10,x20,x30,s0
!!$
!!$  dz=TWOPI/nz/nzperiod
!!$  dt=TWOPI/nt
!!$  DO iz=1,nz
!!$     zeta(iz)=(iz-1.)*dz
!!$  END DO
!!$  DO it=1,nt
!!$     theta(it)=(it-1.)*dt
!!$  END DO
!!$
!!$  x10=5.243070999999999771e-01 
!!$  x20=1.319277999999999906e+00 
!!$  x30=2.785726000000000013e-01
!!$  s0=5.964831933927392527e-01*5.964831933927392527e-01
!!$  dist=1e10
!!$  DO iz=1,nz
!!$     DO it=1,nt
!!$        CALL CALC_XYZ(s,zeta(iz),theta(it),x1(iz,it),x2(iz,it),x3(iz,it),Bzt(iz,it),flag_plot)
!!$        newdist=(x1(iz,it)-x10)*(x1(iz,it)-x10)+(x2(iz,it)-x20)*(x2(iz,it)-x20)+(x3(iz,it)-x30)*(x3(iz,it)-x30)
!!$        IF(newdist.LT.dist) THEN
!!$           iz0=iz
!!$           it0=it
!!$           dist=newdist
!!$        END IF
!!$     END DO
!!$  END DO
!!$  IF(flag_plot) WRITE(1500+myrank,'(9(1pe13.5))') s0,x10,x20,x30,s,x1(iz0,it0),x2(iz0,it0),x3(iz0,it0),dist
!!$  
!!$END SUBROUTINE FILL_3DGRID_OLD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE FILL_BGRID(nz,nt,s,flagB1)

!-------------------------------------------------------------------------------------------------
!Create a grid of size nz x nt in Boozer angles at flux-surface s
!-if flagB1 is true, B_0+B_1 is used;
!-else, B_0 is used (note that B=B_0 if USE_B0 is false)
!Find location of maximum magnetic field
!-------------------------------------------------------------------------------------------------

  USE GLOBAL
  USE KNOSOS_STELLOPT_MOD
  IMPLICIT NONE
  !Input
  LOGICAL flagB1
  INTEGER nz,nt
  REAL*8 s
  !Others
  CHARACTER*64 filename
  INTEGER iz,it,ifile,iostat
  REAL*8 dz,dt,zeta(nz),theta(nt),Bzt(nz,nt),Jac(nz,nt),vds_Bzt(Nnmp,nz,nt)
  REAL*8 FSA,ftr
  
  Bmax=0
!  IF(KN_STELLOPT(9)) THEN
!     ifile=iout
  IF(flagB1) THEN
     ifile=1200+myrank
     IF(KNOSOS_STELLOPT.AND.LEN(TRIM(KN_EXT)).NE.0) THEN 
        filename='B.map.'//TRIM(KN_EXT)
     ELSE
        IF(numprocs.EQ.1) THEN
           filename="B.map"
        ELSE
           WRITE(filename,'("B.map.",I2.2)') myrank
        END IF
     END IF
  ELSE 
     ifile=1300+myrank
     IF(KNOSOS_STELLOPT.AND.LEN(TRIM(KN_EXT)).NE.0) THEN 
        filename='B0.map.'//TRIM(KN_EXT)
     ELSE
        IF(numprocs.EQ.1) THEN
           filename="B0.map"
        ELSE
           WRITE(filename,'("B0.map.",I2.2)') myrank
        END IF
     END IF
  END IF
  OPEN(unit=ifile,file=filename,form='formatted',action='write',iostat=iostat)
  WRITE(ifile,&
       & '("s \zeta_{Boozer}  \theta_{Boozer}(right-handed)  B[T]  (v_B.\nabla\psi)[A.U.]")')

  dz=TWOPI/nz/nzperiod
  dt=TWOPI/nt
  DO iz=1,nz
     zeta(iz)=(iz-1)*dz
  END DO
  DO it=1,nt
     theta(it)=(it-1)*dt
  END DO

  DO iz=1,nz
     DO it=1,nt
        CALL FILL_BNODE(zeta(iz),theta(it),Jac(iz,it),Bzt(iz,it),vds_Bzt(:,iz,it),flagB1)
        IF(Bzt(iz,it).GT.Bmax) THEN !Find location of maximum magnetic field
           Bmax=Bzt(iz,it)
           zmax=zeta(iz)
           tmax=theta(it)
        END IF
     END DO
  END DO
  avB =FSA(nz,nt,Bzt    ,Jac,1)
  avB2=FSA(nz,nt,Bzt*Bzt,Jac,1)
!  ftr=FSA(nz,nt,sqrt(1.-Bzt/Bmax),Jac,1)
!  WRITE(iout,*) 'ftr',s,ftr
!  IF(.NOT.KNOSOS_STELLOPT.OR.KN_STELLOPT(9).OR.DEBUG) THEN
!     IF(.NOT.KNOSOS_STELLOPT) OPEN(unit=ifile,form='formatted',action='write')
     DO iz=1,nz
        DO it=1,nt
           WRITE(ifile,'(5(1pe13.5),L)') s,zeta(iz),theta(it),Bzt(iz,it),vds_Bzt(1,iz,it),flagB1
        END DO
     END DO
!  END IF

     CLOSE(ifile)
     
  CALL CALC_ABSNABLAPSI(MAL,zeta,theta,absnablar)
  absnablar=SQRT(absnablar*Bzt*Bzt)/psip
  
END SUBROUTINE FILL_BGRID


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE FILL_BNODE(zeta,theta,Jac,Bzt,vds_Bzt,flagB1)

!-------------------------------------------------------------------------------------------------
!Calculate the magnetic field Bzt, the Jacobian Jac, and a representative radial drift vms_Bzt
!at angular position (zeta,theta) in Boozer angles:
!-if flagB1 is true, B_0+B_1 is used;
!-else, B_0 is used (note that B=B_0 if USE_B0 is false)
!-------------------------------------------------------------------------------------------------

  USE GLOBAL
  IMPLICIT NONE
  !Input
  LOGICAL flagB1
  REAL*8 zeta,theta
  !Output
  REAL*8 Jac,Bzt,vds_Bzt(Nnmp)
  !Others
  REAL*8 B_0,dBdz_0,dBdt_0
  REAL*8 B_1,dBdz_1,dBdt_1
  REAL*8 Bzt2,dBdz,dBdt,fdummy,dBdpsi

!  CALL CALCB(zeta,theta,2,flagB1,B_0,dBdz_0,dBdt_0,fdummy,fdummy,B_1,dBdz_1,dBdt_1,&
!       & fdummy,fdummy,fdummy,vds_Bzt(:))
  CALL CALCB(zeta,theta,3,flagB1,B_0,dBdz_0,dBdt_0,dBdpsi,fdummy,B_1,dBdz_1,dBdt_1,&
       & fdummy,fdummy,fdummy,vds_Bzt(:))
  IF(flagB1) THEN
     Bzt=B_0+B_1
     dBdz=dBdz_0+dBdz_1
     dBdt=dBdt_0+dBdt_1
  ELSE
     Bzt=B_0
     dBdz=dBdz_0
     dBdt=dBdt_0
  END IF
  Bzt2=Bzt*Bzt
  Jac=aiBtpBz/Bzt2
  vds_Bzt(1)=(Btheta*dBdz-Bzeta*dBdt)/(aiBtpBz*Bzt)     

END SUBROUTINE FILL_BNODE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


REAL*8 FUNCTION FSA(nz,nt,func,Jac,fdegr)
  
!----------------------------------------------------------------------------------------------- 
!Calculate the flux-surface average of func defined in a grid of size (nz/fdegr) x (nt/fdegr) 
!in the Boozer angles (zeta,theta) using the jacobian Jac
!-----------------------------------------------------------------------------------------------
  
  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER fdegr,nz,nt
  REAL*8 func(nz,nt),Jac(nz,nt)
  !Others
  INTEGER iz,it
  REAL*8 denom

  fsa=0
  denom=0
  DO iz=fdegr,nz,fdegr
     DO it=fdegr,nt,fdegr
        fsa=fsa+func(iz,it)*Jac(iz,it)
        denom=denom+Jac(iz,it)
     END DO
  END DO
  fsa=fsa/denom

  RETURN 

END FUNCTION FSA


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


REAL*8 FUNCTION FSA2(na,nz,thetap,func,Jac,fdegr)
  
!----------------------------------------------------------------------------------------------- 
!Calculate the flux-surface average of func defined in a grid of size (nz/fdegr) x (nt/fdegr) 
!in the Boozer angles (zeta,theta) using the jacobian Jac
!-----------------------------------------------------------------------------------------------
  
  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER fdegr,na,nz
  REAL*8 thetap(na),func(na,nz),Jac(na,nz)
  !Others
  INTEGER ia,iz
  REAL*8 dthetap,denom

  fsa2=0
  denom=0
  DO ia=fdegr,na,fdegr
     IF(ia.GT.1.AND.ia.LT.na) THEN
        dthetap=thetap(ia+1)-thetap(ia-1)
     ELSE IF(ia.EQ.1) THEN
        dthetap=thetap(2)-thetap(na)+siota*TWOPI
     ELSE IF(ia.EQ.na) THEN
        dthetap=thetap(1)-thetap(na-1)+siota*TWOPI
     END IF
     DO iz=fdegr,nz,fdegr
        fsa2=fsa2+func(ia,iz)*Jac(ia,iz)*dthetap
        denom=denom+Jac(ia,iz)*dthetap
     END DO
  END DO
  fsa2=fsa2/denom
  RETURN 

END FUNCTION FSA2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE FIND_BMIN_BMAX(is,nz,nt,flagb1)
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------

  USE GLOBAL
  IMPLICIT NONE
  !Input
  LOGICAL flagB1
  INTEGER is,nz,nt
  !Others
  INTEGER iz,it
  REAL*8 dz,dt,zeta,theta,B,dBdz,dBdt,dummy,vdummy(Nnmp)
  REAL*8 dzeta,dtheta,dBdz_old,dBdt_old,zmin,tmin

  Bmax=borbic(0,0)
  Bmin=borbic(0,0)
  dz=TWOPI/nz/nzperiod
  dt=TWOPI/nt
  DO iz=1,nz
     zeta=(iz-1)*dz
     DO it=1,nt
        theta=(it-1)*dt
        CALL CALCB(zeta,theta,0,flagb1,B,dBdz,dBdt,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,vdummy)
        IF(B.GT.Bmax) THEN
           Bmax=B
           zmax=zeta
           tmax=theta
        ELSE IF (B.LT.BMIN) THEN
           Bmin=B
           zmin=zeta
           tmin=theta
        END IF
     END DO
  END DO
  DO it=1,2
     IF(it.EQ.1) THEN
        zeta =zmax
        theta=tmax
     ELSE
        zeta =zmin
        theta=tmin
     END IF
     CALL CALCB(zeta,theta,2,flagb1,B,dBdz,dBdt,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,vdummy)
     dzeta =dz/2  
     dtheta=dzeta*dBdt/dBdz
     dBdz_old=dBdz
     dBdt_old=dBdt
     DO iz=1,20
        zeta =zeta +dzeta
        theta=theta+dtheta
        CALL CALCB(zeta,theta,2,flagb1,B,dBdz,dBdt,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,vdummy)
        IF((it.EQ.1.AND.B.GT.Bmax).OR.(it.EQ.2.AND.B.LT.Bmin)) THEN
           IF(it.EQ.1) THEN
              Bmax=B
           ELSE
              Bmin=B
           END IF
           dzeta =dzeta *dBdz/dBdz_old
           dtheta=dtheta*dBdt/dBdt_old
           dBdz_old=dBdz
           dBdt_old=dBdt
        ELSE
           zeta =zeta -dzeta
           theta=theta-dtheta
           dzeta=dzeta  /2
           dtheta=dtheta/2
        END IF
     END DO
  END DO
  
  Bmax=Bmax/(1.0+(Bmax/Bmin-1.0)/(2*MLAMBDA))   !move the trapped/passing boundary

  IF(is.GT.0) THEN
     Bmax_b(is)=Bmax
     Bmin_b(is)=Bmin
  END IF
  
END SUBROUTINE FIND_BMIN_BMAX


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

SUBROUTINE INTERPOLATE_2DMAP(npt,zeta_i,theta_i,func_i,func_o)

!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER npt
  REAL*8 zeta_i(npt,npt),theta_i(npt,npt),func_i(npt,npt)
  !Output
  REAL*8 func_o(npt,npt)
  !Others
  INTEGER zer_nptpiz,one_nptpiz,two_nptpiz,zer_nptpit,one_nptpit,two_nptpit,thr_npt
  INTEGER iz,it,info,ipivot(10),ipoint,jpoint
  REAL*8 dist,zeta_t,theta_t,dist_t
  REAL*8 zeta(npt),theta(npt)
  REAL*8 zeta_e(9*npt*npt),theta_e(9*npt*npt),func_e(9*npt*npt)
  REAL*8 ztl,ttl,ftl,ztr,ttr,ftr,zbr,tbr,fbr,zbl,tbl,fbl
  REAL*8 mat(10,10),rhs(10)
!  INTEGER IWK(max(31,27+4)*9*npt*npt+npt*npt)
!  REAL*8 WK(5*9*npt*npt)
  !Time
  CHARACTER*30, PARAMETER :: routine="FILL_3DPOINTS"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart
  INTEGER, PARAMETER :: md=1
  INTEGER, PARAMETER :: ncp=4

  CALL CPU_TIME(tstart)

  DO iz=1,npt
     DO it=1,npt

        zer_nptpiz=iz
        one_nptpiz=npt+iz
        two_nptpiz=2*npt+iz
        zer_nptpit=it
        one_nptpit=npt+it
        two_nptpit=2*npt+it
        thr_npt=3*npt
        
         zeta_e((one_nptpiz-1)*thr_npt+one_nptpit)= zeta_i(iz,it)
        theta_e((one_nptpiz-1)*thr_npt+one_nptpit)=theta_i(iz,it)
         func_e((one_nptpiz-1)*thr_npt+one_nptpit)= func_i(iz,it)

         zeta_e((zer_nptpiz-1)*thr_npt+zer_nptpit)= zeta_i(iz,it)-TWOPI/nzperiod
        theta_e((zer_nptpiz-1)*thr_npt+zer_nptpit)=theta_i(iz,it)-TWOPI
         func_e((zer_nptpiz-1)*thr_npt+zer_nptpit)= func_i(iz,it)
        
         zeta_e((zer_nptpiz-1)*thr_npt+one_nptpit)= zeta_i(iz,it)-TWOPI/nzperiod
        theta_e((zer_nptpiz-1)*thr_npt+one_nptpit)=theta_i(iz,it)
         func_e((zer_nptpiz-1)*thr_npt+one_nptpit)= func_i(iz,it)

         zeta_e((zer_nptpiz-1)*thr_npt+two_nptpit)= zeta_i(iz,it)-TWOPI/nzperiod
        theta_e((zer_nptpiz-1)*thr_npt+two_nptpit)=theta_i(iz,it)+TWOPI
         func_e((zer_nptpiz-1)*thr_npt+two_nptpit)= func_i(iz,it)

         zeta_e((one_nptpiz-1)*thr_npt+zer_nptpit)= zeta_i(iz,it)         
        theta_e((one_nptpiz-1)*thr_npt+zer_nptpit)=theta_i(iz,it)-TWOPI
         func_e((one_nptpiz-1)*thr_npt+zer_nptpit)= func_i(iz,it)

         zeta_e((one_nptpiz-1)*thr_npt+two_nptpit)= zeta_i(iz,it)
        theta_e((one_nptpiz-1)*thr_npt+two_nptpit)=theta_i(iz,it)+TWOPI
         func_e((one_nptpiz-1)*thr_npt+two_nptpit)= func_i(iz,it)

         zeta_e((two_nptpiz-1)*thr_npt+zer_nptpit)= zeta_i(iz,it)+TWOPI/nzperiod
        theta_e((two_nptpiz-1)*thr_npt+zer_nptpit)=theta_i(iz,it)-TWOPI
         func_e((two_nptpiz-1)*thr_npt+zer_nptpit)= func_i(iz,it)

         zeta_e((two_nptpiz-1)*thr_npt+one_nptpit)= zeta_i(iz,it)+TWOPI/nzperiod
        theta_e((two_nptpiz-1)*thr_npt+one_nptpit)=theta_i(iz,it)
         func_e((two_nptpiz-1)*thr_npt+one_nptpit)= func_i(iz,it)

         zeta_e((two_nptpiz-1)*thr_npt+two_nptpit)= zeta_i(iz,it)+TWOPI/nzperiod         
        theta_e((two_nptpiz-1)*thr_npt+two_nptpit)=theta_i(iz,it)+TWOPI
        func_e((two_nptpiz-1)*thr_npt+two_nptpit)= func_i(iz,it)

     END DO
  END DO

  DO iz=1,3*npt
     DO it=1,3*npt
        WRITE(iout,'(I2,8(1pe13.5))') iz-iz,zeta_e((iz-1)*thr_npt+it),theta_e((iz-1)*thr_npt+it),func_e((iz-1)*thr_npt+it)
     END DO
  END DO

  DO iz=1,npt
     zeta(iz) =(iz-1)*(TWOPI/nzperiod)/npt
     theta(iz)=(iz-1)* TWOPI/npt
  END DO

  DO iz=1,npt
     DO it=1,npt

        DO jpoint=1,4
           dist=1e9
           DO ipoint=1,9*npt*npt
              zeta_t =zeta_e(ipoint)- zeta(iz)
              theta_t=theta_e(ipoint)-theta(it)
              IF(jpoint.EQ.1.AND.(zeta_t.GT.0.OR.theta_t.GT.0)) CYCLE
              IF(jpoint.EQ.2.AND.(zeta_t.LT.0.OR.theta_t.GT.0)) CYCLE
              IF(jpoint.EQ.3.AND.(zeta_t.GT.0.OR.theta_t.LT.0)) CYCLE
              IF(jpoint.EQ.4.AND.(zeta_t.LT.0.OR.theta_t.LT.0)) CYCLE
              dist_t=zeta_t*zeta_t+theta_t*theta_t
              IF(dist_t.LT.dist) THEN
                 dist=dist_t
                 IF(jpoint.EQ.1) THEN
                    ztl=zeta_t
                    ttl=theta_t
                    ftl=func_e(ipoint)
                 ELSE IF(jpoint.EQ.2) THEN
                    ztr=zeta_t
                    ttr=theta_t
                    ftr=func_e(ipoint)
                 ELSE IF(jpoint.EQ.3) THEN
                    zbl=zeta_t
                    tbl=theta_t
                    fbl=func_e(ipoint)
                 ELSE IF(jpoint.EQ.4) THEN
                    zbr=zeta_t
                    tbr=theta_t
                    fbr=func_e(ipoint)
                 END IF
              END IF
           END DO
        END DO
        
!!$        ztr= zeta_e((npt+iz+1)*thr_npt+npt+it+1) -zeta(iz)
!!$        ttr=theta_e((npt+iz+1)*thr_npt+npt+it+1)-theta(it)
!!$        ftr= func_e((npt+iz+1)*thr_npt+npt+it+1)
!!$        ztl= zeta_e((npt+iz-1)*thr_npt+npt+it+1) -zeta(iz)
!!$        ttl=theta_e((npt+iz-1)*thr_npt+npt+it+1)-theta(it)
!!$        ftl= func_e((npt+iz-1)*thr_npt+npt+it+1)
!!$        zbr= zeta_e((npt+iz+1)*thr_npt+npt+it-1) -zeta(iz)
!!$        tbr=theta_e((npt+iz+1)*thr_npt+npt+it-1)-theta(it)
!!$        fbr= func_e((npt+iz+1)*thr_npt+npt+it-1)
!!$        zbl= zeta_e((npt+iz-1)*thr_npt+npt+it-1) -zeta(iz)
!!$        tbl=theta_e((npt+iz-1)*thr_npt+npt+it-1)-theta(it)
!!$        fbl= func_e((npt+iz-1)*thr_npt+npt+it-1)
        mat=0
        rhs=0
        mat(1,1)=1
        mat(2,2)=1
        mat(3,3)=1
        mat(1,7)=ztl*ztl
        mat(2,7)=ztl*ttl
        mat(3,7)=ttl*ttl
        mat(4,7)=ztl
        mat(5,7)=ttl
        mat(6,7)=1
        mat(1,8)=ztr*ztr
        mat(2,8)=ztr*ttr
        mat(3,8)=ttr*ttr
        mat(4,8)=ztr
        mat(5,8)=ttr
        mat(6,8)=1
        mat(1,9)=zbr*zbr
        mat(2,9)=zbr*tbr
        mat(3,9)=tbr*tbr
        mat(4,9)=zbr
        mat(5,9)=tbr
        mat(6,9)=1
        mat(1,10)=zbl*zbl
        mat(2,10)=zbl*tbl
        mat(3,10)=tbl*tbl
        mat(4,10)=zbl
        mat(5,10)=tbl
        mat(6,10)=1
        mat(7,1)=ztl*ztl
        mat(7,2)=ztl*ttl
        mat(7,3)=ttl*ttl
        mat(7,4)=ztl
        mat(7,5)=ttl
        mat(7,6)=1
        mat(8,1)=ztr*ztr
        mat(8,2)=ztr*ttr
        mat(8,3)=ttr*ttr
        mat(8,4)=ztr
        mat(8,5)=ttr
        mat(8,6)=1
        mat(9,1)=zbr*zbr
        mat(9,2)=zbr*tbr
        mat(9,3)=tbr*tbr
        mat(9,4)=zbr
        mat(9,5)=tbr
        mat(9,6)=1
        mat(10,1)=zbl*zbl
        mat(10,2)=zbl*tbl
        mat(10,3)=tbl*tbl
        mat(10,4)=zbl
        mat(10,5)=tbl
        mat(10,6)=1

        rhs(7)=ftl
        rhs(8)=ftr
        rhs(9)=fbr
        rhs(10)=fbl

        CALL DGESV(10,1,mat,10,ipivot,rhs,10,info)
        func_o(iz,it)=rhs(6)
        !        IF(DEBUG.OR.myrank.EQ.0)
!        CALL IDBVIP(2,4,9*npt*npt,zeta_e,theta_e,func_e,1,zeta(iz),theta(it),func_o(iz,it),iwk,wk)        
        WRITE(iout,'(I2,7(1pe13.5),I3)') iz-iz-1,zeta(iz),theta(it),func_o(iz,it),ftl,ftr,fbl,fbr
     END DO
  END DO

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)
  
END SUBROUTINE INTERPOLATE_2DMAP


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE FFTF_KN(nalphab,q,qnm)

!-------------------------------------------------------------------------------------------------
!Performs FFT of real q(nalphab,nalphab) and produces COMPLEX qnm(nalphab,nalphab)
!-------------------------------------------------------------------------------------------------

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER nalphab
  REAL*8 q(nalphab,nalphab)
  !Output
  COMPLEX*16 qnm(nalphab,nalphab)
  !Others
#ifdef NAG
  INTEGER ierr
  REAL*8 TRIGN(2*nalphab),TRIGM(2*nalphab),WORK(2*nalphab*nalphab)
  REAL*8 qr(nalphab,nalphab),qi(nalphab,nalphab)
#else
  INTEGER, PARAMETER :: FFTW_ESTIMATE =64
  INTEGER, PARAMETER :: FFTW_FORWARD=-1
  COMPLEX*16 qc(nalphab,nalphab)
#endif  

  qnm=0

#ifdef NAG

  qr=REAL(q)
  qi=0
  CALL C06FUF(nalphab,nalphab,qr,qi,'I',&
       & TRIGM(1:2*nalphab),TRIGN(1:2*nalphab),work(1:2*nalphab*nalphab),ierr)
  qnm=CMPLX(qr,qi)/nalphab

#else

  qc=0
  qc=q
  IF(plan_fwd.EQ.0) CALL DFFTW_PLAN_DFT_2d(plan_fwd,nalphab,nalphab,qc,qnm,&
       & FFTW_FORWARD,FFTW_ESTIMATE)
  CALL DFFTW_EXECUTE_DFT(plan_fwd,qc,qnm)
  qnm=qnm/(nalphab*nalphab)

#endif

END SUBROUTINE FFTF_KN



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE FFTB_KN(nalphab,qnm,q)

!-------------------------------------------------------------------------------------------------
!Performs FFT of real q(nalphab,nalphab) and produces COMPLEX qnm(nalphab,nalphab)
!-------------------------------------------------------------------------------------------------

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER nalphab
  COMPLEX*16 qnm(nalphab,nalphab)
  !Output
  REAL*8 q(nalphab,nalphab)
  !Others
#ifdef NAG
  INTEGER ierr
  REAL*8 TRIGN(2*nalphab),TRIGM(2*nalphab),WORK(2*nalphab*nalphab)
  REAL*8 qr(nalphab,nalphab),qi(nalphab,nalphab)
#else
  INTEGER, PARAMETER :: FFTW_ESTIMATE =64
  INTEGER, PARAMETER :: FFTW_BACKWARD=+1
  COMPLEX*16 qc(nalphab,nalphab)
#endif  

#ifdef NAG

  qr=REAL(qnm)
  qi=AIMAG(qnm)
  CALL C06GCF(qi,nalphab*nalphab,ierr)
  CALL C06FUF(nalphab,nalphab,qr,qi,'I',&
       & TRIGM(1:2*nalphab),TRIGN(1:2*nalphab),work(1:2*nalphab*nalphab),ierr)
!  CALL C06GCF(qi,nalphab*nalphab,ierr)
  q=qr

#else

  IF(plan_bwd.EQ.0) CALL DFFTW_PLAN_DFT_2d(plan_bwd,nalphab,nalphab,qnm,qc,&
       & FFTW_BACKWARD,FFTW_ESTIMATE)
  CALL DFFTW_EXECUTE_DFT(plan_bwd,qnm,qc)
  q=REAL(qc)

#endif

END SUBROUTINE FFTB_KN


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE QUIT_READ(filename,varname,status_nc)  !2020/06 Satake for NETCDF

  USE GLOBAL, ONLY : myrank,iout
  IMPLICIT none
  !Input
  CHARACTER*60 filename
  CHARACTER*40 varname
  INTEGER status_nc
  !Others
  CHARACTER*100 serr

  IF(myrank==0) THEN 
     WRITE(iout,*) "error read NETCDF file!"
     WRITE(iout,*) filename,varname,status_nc
  END IF
  serr="netcdf"
  CALL END_ALL(serr,.FALSE.)
  
END SUBROUTINE QUIT_READ


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!$  !The following lines correspond to specific physics studies:
!!$  !------------------------------------------------------------------------------------------- 
!!$  !For [Calvo 2018 JPP], uses B0 from [Landreman 2011 PoP] and scan in aspect ratio
!!$  !------------------------------------------------------------------------------------------- 
!!$  IF(JPP) THEN
!!$     !Undo change of coordinates (it was already right-handed)...
!!$     Btheta   =-Btheta
!!$     chip     =-chip
!!$     iota     =-iota
!!$     porbis   =-porbis
!!$     zorbis   =-zorbis
!!$     iotap=iota
!!$     iotam=iota
!!$     Bzetap=Bzeta
!!$     Bzetam=Bzeta
!!$     Bthetap=Btheta
!!$     Bthetam=Btheta
!!$     !..except for Bmn
!!$     torbic=borbic
!!$     borbic=0
!!$     DO n=-ntorb,ntorb
!!$        borbic(n,:)=torbic(-n,:)
!!$     END DO
!!$     torbic=borbic0
!!$     borbic0=0
!!$     DO n=-ntorb,ntorb
!!$        borbic0(n,:)=torbic(-n,:)
!!$     END DO
!!$     torbic=dborbicdpsi
!!$     dborbicdpsi=0
!!$     DO n=-ntorb,ntorb
!!$        dborbicdpsi(n,:)=torbic(-n,:)
!!$     END DO
!!$     !Multiply all modes except B00 times FE (first substracts 0.072 so that B00=1)
!!$!     borbi0 =borbi0*FE
!!$!     borbi  =borbi  *FE
!!$!     borbi0(0,0) =borbi0(0,0) /FE
!!$!     borbi(0,0)  =borbi(0,0)  /FE
!!$!     PREC_B=PREC_B*FE
!!$!     IF THEN
!!$!        borbi  =borbi *5
!!$!        borbi0 =borbi0*5
!!$!        PREC_B=PREC_B *5
!!$        !Change radial derivatives according to value of FE
!!$!        rad_a=0.072*rad_R*FE
!!$!     ELSE
!!$!        !All modes multiplied in order to change omegastar
!!$        borbic =borbic *10000
!!$        borbic0=borbic0*10000
!!$        PREC_B=PREC_B*1000
!!$!        !Radial derivatives changed according to value of FE
!!$        rad_a=0.072*rad_R*FE
!!$!     END IF
!!$     torflux=borbic0(0,0)*rad_a*rad_a/2
!!$     psip=borbic0(0,0)*rad_a*SQRT(s0)
!!$     chip=iota*psip
!!$     dpsi=dpsi*10000
!!$     !The radial derivative is made up
!!$     dborbicdpsi=0
!!$     dborbicdpsi=borbic0/rad_a/SQRT(s0)/psip!*1E16x
!!$     dborbicdpsi(0,0)=0
!!$     dborbicdpsi(0,2)=0
!!$
!!$  !------------------------------------------------------------------------------------------- 
!!$  !For benchmarking with FORTEC-3D, tokamak perturbed with a gaussian profile
!!$  !------------------------------------------------------------------------------------------- 
!!$  ELSE IF(SATAKE) THEN  
!!$     ntorb=MAX(ntorb,3)
!!$     mpolb=MAX(ABS(mpolb),6)
!!$     borbic(3,6) =1.90*0.002*EXP(-(SQRT(s0)-0.5)*(SQRT(s0)-0.5)/0.01)
!!$     borbic0=borbic
!!$  END IF

  !------------------------------------------------------------------------------------------- 
  !End of specific physics studies
  !------------------------------------------------------------------------------------------- 


!Change the coordinates if the main helicity is N=0 (tokamak or quasi-axissymmetric stellarator, QAS)
!!$  IF(NTV) THEN
!!$!     serr="NOT implemented, look for SATAKE"
!!$!     CALL END_ALL(serr,.FALSE.)
!!$!     siota=iota/ABS(iota)
!!$!     Bzeta=Btheta
!!$!     Btheta=-siota*fdummy
!!$!     psip=psip*ABS(iota)
!!$!     dborbicdpsi=dborbicdpsi/ABS(iota)
!!$!     diotadpsi=diotadpsi/iota/iota/iota
!!$!     iota=-siota/iota
!!$!     chip=psip*iota
!!$     sgnhel=-1
!!$     fdummy=Bzeta
!!$     Bzeta =Btheta
!!$     Btheta=sgnhel*fdummy
!!$     psip=psip*ABS(iota)
!!$     dpsi=dpsi*ABS(iota)
!!$     dborbicdpsi=dborbicdpsi/ABS(iota)
!!$     dborbisdpsi=dborbisdpsi/ABS(iota)
!!$     dBzdpsi=dBzdpsi/ABS(iota)
!!$     dBtdpsi=dBtdpsi/ABS(iota)
!!$     torflux=torflux*ABS(iota)
!!$     atorflux=atorflux*ABS(iota)
!!$     diotadpsi=diotadpsi/iota/iota/iota
!!$     iota=sgnhel/iota
!!$     chip=psip*iota
!!$  END IF
