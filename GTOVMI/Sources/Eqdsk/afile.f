
      module sparmd
      INTEGER, PARAMETER :: magpri67=29,magpri322=31,magprirdp=8
      INTEGER, PARAMETER :: magpri=magpri67+magpri322+magprirdp
      INTEGER, PARAMETER :: nfcoil=18,nsilop=41,nrogow=1,ntime=1
      INTEGER, PARAMETER :: npca=70,nparm=1,ncoef=41
      INTEGER, PARAMETER :: necoil=122,mpress=10,nvesel=24
      INTEGER, PARAMETER :: nffcur=5,nppcur=5,npcurn=nffcur+nppcur
     .     ,mfnpcr=nfcoil+npcurn,npcur2=npcurn*2
     .     ,nrsmat=nsilop+magpri+nrogow+nffcur+1+npcurn+mpress+nfcoil
      INTEGER, PARAMETER :: npoint=300
      INTEGER, PARAMETER :: nw=65,nh=65,nwnh=nw*nh
      INTEGER, PARAMETER :: nh2=2*nh,nwrk=2*(nw+1*nh)
      INTEGER, PARAMETER :: nesum=6
      INTEGER, PARAMETER :: ncurrt=nvesel+nesum+nfcoil
c     INTEGER, PARAMETER :: mbdry=MAX(nsilop+magpri+nrogow,60)
      INTEGER, PARAMETER :: mbdry=nsilop+magpri+nrogow
c     INTEGER, PARAMETER :: nbwork=MAX(nsilop,60-magpri-nrogow)
      INTEGER, PARAMETER :: nbwork=nsilop
      INTEGER, PARAMETER :: kxiter=100
      INTEGER, PARAMETER :: nlimit=120,nlimbd=6
      INTEGER, PARAMETER :: msbdry=mbdry+nsilop+nfcoil+1,
     .	msbdr2=2*msbdry
      INTEGER, PARAMETER :: nrsma2=2*nrsmat
c     INTEGER, PARAMETER :: nwf=MAX(nw,nfcoil)
      INTEGER, PARAMETER :: nwf=nw
      INTEGER, PARAMETER :: nxtram=10,nxtrap=npoint
      INTEGER, PARAMETER :: nxtlim=9,nco2v=3,nco2r=2
      INTEGER, PARAMETER :: nslit=4,nangle=16
c
      INTEGER, PARAMETER :: mxmagpri=100,mxnfcoil=50,mxnsilop=100,
     .	mxnesum=20
      end module sparmd
      module mcomd1
      use sparmd

      common/inaver/iavem,iaved,iavev
      real*8 :: rotam(npca,npca),coef(ncoef,nparm),pmean(nparm),
     .     psigma(nparm),qmean(npca),qsigma(npca),lcoef,nlcoef,
     .     eigen(npca),cfwtsi(mxnsilop),cfwtmp2(mxmagpri),mslref
      real*8 :: rotamd(npca,npca),coefd(ncoef,nparm),pmeand(nparm),
     .     psigmad(nparm),qmeand(npca),qsigmad(npca)
      integer :: lcoefd,nlcoefd
     .     eigend(npca),cfwtsid(mxnsilop),cfwtmp2d(mxmagpri)
      integer :: mslrefd
      real*8 :: ivesel,volecs(mxnesum),volecc(mxnesum)
     .  ,rsisvs(nvesel),efreq,sumvs0,volfcs(mxnfcoil)
     .  ,volfcc(mxnfcoil),rsisec(mxnesum),powvs,pvscur,pscurn,ppscur
      real*8 :: s1(ntime),s2(ntime),s3(ntime),bpolav(ntime)
      real*8 :: bpol(npoint),plengt(npoint),bpolz(npoint),siar,siaz
      real*8 :: eout(ntime),rout(ntime),zout(ntime),doutu(ntime)
     .  ,doutl(ntime),aout(ntime),vout(ntime),betat(ntime),otop(ntime)
     .  ,betap(ntime),ali(ntime),oleft(ntime),oright(ntime),qsta(ntime)
     .  ,rcurrt(ntime),zcurrt(ntime),qout(ntime),olefs(ntime)
     .  ,orighs(ntime),otops(ntime),sibdry(ntime),areao(ntime)
     .  ,wplasm(ntime),elongm(ntime),qqmagx(ntime),terror(ntime)
     .  ,rmagx(ntime),zmagx(ntime),obott(ntime),obots(ntime)
     .  ,alpha(ntime),rttt(ntime),dbpli(ntime),delbp(ntime)
     .  ,rseps(2,ntime),zseps(2,ntime),sepexp(ntime),shearb(ntime)
     .  ,xtch(ntime),ytch(ntime),qpsib(ntime),vertn(ntime),aaq1(ntime)
     .  ,aaq2(ntime),aaq3(ntime),btaxp(ntime),btaxv(ntime)
     .  ,simagx(ntime),jerror(ntime),seplim(ntime)
     .  ,wbpol(ntime),taumhd(ntime),betapd(ntime),betatd(ntime)
     .  ,alid(ntime),wplasmd(ntime),taudia(ntime),wbpold(ntime)
     .  ,qmerci(ntime),slantu(ntime),slantl(ntime),zeff(ntime)
     .  ,zeffr(ntime),tave(ntime),civs(ntime),civst(ntime),
     .   cipmp2(ntime),frvst(ntime),fzvst(ntime),vsurfa(ntime),
     .   wpdot(ntime),wbdot(ntime),rvsin(ntime),zvsin(ntime),
     .   rvsout(ntime),zvsout(ntime),cjor95(ntime),pp95(ntime),
     .   ssep(ntime),yyy2(ntime),xnnc(ntime),cprof(ntime),qqmin(ntime)
     .  ,chigamt(ntime),cjor0(ntime),fexpan(ntime),ssi01(ntime),
     .   fexpvs(ntime),sepnose(ntime),ssi95(ntime),rqqmin(ntime)
     .  ,cjor99(ntime),cj1ave(ntime)
     .  ,rmidin(ntime),rmidout(ntime),psurfa(ntime)
      real*8 :: zuperts(ntime),zlowerts,rmajts
      real*8 :: xout(npoint),yout(npoint),dpsi,rymin,rymax,
     .  zxmin,zxmax,xmin,xmax,ymin,ymax,rmaxis,zmaxis,nfound,emaxis
      real*8 :: psi(nwnh),psibry,simag,sidif,xpsi(nwnh),eouter
      real*8 :: csilop(mxnsilop,ntime),crogow(nrogow,ntime),
     .  cmpr2(mxmagpri,ntime),cpasma(ntime),xndnt(ntime)
     . ,cbetap,cli,cqqxis,cbetat,ci0
     . ,ccbrsp(mxnfcoil,ntime)
      real*8 :: pi,tmu,twopi,ibunmn,tmu2,errcut
      integer :: kinput,jwake
      real*8 :: cstabz
      integer :: ivacum
      real*8 :: rfila(npcurn),zfila(npcurn)
      integer :: idfila
      integer :: nin,nout,ntty,nrsppc,nrspfc,nttyo,neqdsk,nffile,nsave
      real*8 :: brsp(nrsmat),tsaisq(ntime),cond
      integer :: ishot,itime,nparam,nfnpcr
     .     ,kpcurn,kersil
      real*8 :: bcentr(ntime),rcentr
      real*8 :: serror,fwtsi(mxnsilop),fwtmp2(mxmagpri),
     .     fitwgt(nrogow),fwtcur,elomin,fwtbp,fwtdlc,fwtfc(mxnfcoil),
     .     rwtsi(mxnsilop),rwtmp2(mxmagpri)
      real*8 :: silopt(ntime,mxnsilop),expmpi(ntime,mxmagpri),
     .     prexpt(ntime,nrogow),pasmat(ntime),xtime(ntime),
     .     ierpsi(mxnsilop),ierpr(nrogow),iermpi(mxmagpri),ierpla
     .     ,psibit(mxnsilop),prbit(nrogow),bitmpi(mxmagpri),bitip
     .     ,fccurt(ntime,mxnfcoil),vloopt(ntime),psiref(ntime)
     .     ,denvt(ntime,nco2v),denrt(ntime,nco2r),eccurt(ntime,mxnesum)
     .     ,bitfc(mxnfcoil),ierfc(mxnfcoil),diamag(ntime)
     .     ,sigdia(ntime),ierdia(3),pbinj(ntime)
      integer :: ico2
      real*8 :: saisil(mxnsilop),saimpi(mxmagpri),saipr(nrogow)
      real*8 :: dfluxc(ntime),sigdlc,rspdlc(nffcur),cdflux(ntime)
      real*8 :: mxiter,idone,error,errorm,nitera,nxiter,ixnn,errmin
      real*8 :: rco2r(nco2r,ntime),rco2v(nco2v,ntime),chordv(nco2v)
     .              ,chordr(nco2r),zcentr,dco2r(ntime,nco2r)
     .              ,dco2v(ntime,nco2v)
      real*8 :: sbpp
      integer :: nvernum,lflag
      real*8 :: kflag(30),ktimeo
      real*8 :: darea,drgrid,dzgrid,qmaxis,cratio,dfsqe
      real*8 :: rxray(nslit),zxray(nslit),xangle(nangle,nslit)
      integer ::  icondn,itek,kdata,itrace,ierchk,iconvr,ixray
     .     ,ichisq,modep,ibound,ibatch,idite,ilaser,islant
     .     ,chimin,cutip,lookfw,iplotvs,iplotmp,iplotsi
     .     ,mplotmp(3),mplotsi(3),ifitvs,idiamag,ipca,iheader,keqdsk
     .     ,negcur,itest_nc,kbound,kplot10
      real*8 :: nqaxis,isumip,sumip,fbetap,fli,fbetat
      real*8 :: rsi(mxnsilop),zsi(mxnsilop),wsi(mxnsilop)
     .     ,hsi(mxnsilop)
      real*8 :: as(mxnsilop),as2(mxnsilop)
      integer :: nslref
      real*8 :: rf(mxnfcoil),zf(mxnfcoil),wf(mxnfcoil),hf(mxnfcoil),
     .     af(mxnfcoil),af2(mxnfcoil),rsisfc(mxnfcoil)
      real*8 :: xmp2(mxmagpri),ymp2(mxmagpri),
     .     amp2(mxmagpri),smp2(mxmagpri),nsmp2,patmp2(mxmagpri)
      real*8 :: scrape,nextra,ixstrt,iextra,iprobe,ifcoil,iecoil
     .     ,iexcal,iconsi
      real*8 :: limitr,xlim(nlimit),ylim(nlimit),iplim,limfag
     .  ,limitr_180,xlim_180(nlimit),ylim_180(nlimit)
      real*8 :: xlmin,xlmax,ylmin,ylmax,limid,limup,limbot
      real*8 :: ipsi(ntime),irogw(ntime),imag2(ntime),iplasm(ntime)
     .     ,idlopc(ntime)
      real*8 :: iopen,ifread,errbry,xncoil,itcoil,ifcurr
      real*8 :: nbbbs,rbbbs(mbdry),zbbbs(mbdry)
      real*8 :: cfcoil,fcsum(mxnfcoil),fczero(mxnfcoil),ifref
      real*8 :: zero(nwnh)
      real*8 :: xlimbd(nlimbd),ylimbd(nlimbd)
      integer :: nsilop0,magpri0,nfcoil0,nesum0,turnfc(nfcoil),turnbc
      dimension pcurrt(npcurn)
      equivalence (pcurrt(1),brsp(mxnfcoil+1))
      dimension :: mfvers(2)
      character mfvers*5,terminal*100
      dimension :: vsname(nvesel),mpnam2(mxmagpri),lpname(mxnsilop)
      character*10 vsname,mpnam2,lpname
      character*15 ifname(ntime)
      character*4 cvax
      character*3 limloc(ntime)
      character*25 filimt
      character cshot*6,mfitpop*12
      integer erflag(ntime,30)
      data nhist/37/,nslref/1/
      data mfvers(1)/'07/19'/,mfvers(2)/'/99  '/
      data  turnfc/5*58,2*55,58,55,5*58,2*55,58,55/,turnbc/144/
      end module mcomd1
