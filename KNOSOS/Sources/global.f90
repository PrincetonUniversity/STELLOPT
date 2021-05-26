!Global variables

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE GLOBAL

  IMPLICIT NONE

  !Maximum number of points in (alpha,lambda), wells, points
  INTEGER, PARAMETER :: nlambdax  =1024
  INTEGER, PARAMETER :: nax       =512
  INTEGER, PARAMETER :: npointx   =400000
  INTEGER, PARAMETER :: nwx       =5000
  INTEGER, PARAMETER :: nsamp     =8

  !Maximum number of species and radial positions
  INTEGER, PARAMETER :: nsx       =100
  INTEGER, PARAMETER :: nbx       =11

  !Bounce integrals
  INTEGER, PARAMETER :: nq0=7
  INTEGER, PARAMETER :: nqv=3

  !Algebraic constants
  REAL*8, PARAMETER :: PI          = 3.14159265358979
  REAL*8, PARAMETER :: TWOPI       = 6.28318530717959
  REAL*8, PARAMETER :: SQPI        = 1.77245385090552
  REAL*8, PARAMETER :: SQ2         = 1.4142135623731
  REAL*8, PARAMETER :: ZERO        = 0.00000000000000
  REAL*8, PARAMETER :: ONE         = 1.00000000000000
  REAL*8, PARAMETER :: MONE        =-1.00000000000000
  REAL*8, PARAMETER :: ALMOST_ZERO = 1.00000000000E-9
  REAL*8, PARAMETER :: SMALL       = 1.00000000000E-7
  COMPLEX*16, PARAMETER :: NUMI    = CMPLX(0.0,1.0)

  !Physics constants
  REAL*8, PARAMETER :: m_e =1.04396844E-08 !ratio between unitary mass and charge
  REAL*8, PARAMETER ::   e =1.60217662e-19 !electron charge
  REAL*8, PARAMETER :: tsl=1E-3
  
  !MPI constants
  INTEGER iout,myrank,numprocs

  !Input parameters
  !Namelist /param/
  LOGICAL GEN_FLAG(20),DEBUG,TIME,PLOTG,USE_SHAING,SHAING_1NU,SHAING_SQRTNU,SHAING_SBP
  LOGICAL USE_B1,USE_B0pB1,QS_B0,QS_B0_1HEL,REMOVE_DIV,DELTA,FAST_AMB,ILAGRID,RSEED
  INTEGER I0,L0,TRUNCATE_B,NER,NERR,MLAMBDA,MAL,NTURN
  REAL*8 PREC_B,PREC_EXTR,PREC_BINT,PREC_DQDV,PREC_INTV,IPERR,ERMIN,ERMAX,ERACC
  !Namelist /model/
  CHARACTER*100 DIRDB
  CHARACTER*7 DIRS(nsx)
  LOGICAL ONLY_B0,CALC_DB,NO_PLATEAU,ONLY_DB,INC_EXB,TANG_VM,CLASSICAL,FRICTION,ANISOTROPY
  LOGICAL SCAN_ER,SOLVE_AMB,TRIVIAL_AMB,SOLVE_QN,TRIVIAL_QN,ZERO_PHI1,ONLY_PHI1,COMPARE_MODELS
  LOGICAL NEOTRANSP,TASK3D,TASK3Dlike,PENTA,NTV,SATAKE,ANA_NTV,JPP,ESCOTO,NEQ2
  INTEGER D_AND_V,ER_ROOT,NBULK
  REAL*8 FN,FI,FS,FP,FE,FR,FB,FACT_CON,FNE,FTE,FTI,FDLNE,FDLTE,FDLTI,FER
  !Namelist /fastions/
  LOGICAL JMAP,LINEART,GLOBALFI,MODELFI,RANDOMSL,JCORRECTION,JTRANS,FITRANS
  INTEGER JORBIT,FDMAX
  REAL*8 GTH,TENDFI,DTFI,EFI,TWOEFI,LJMAP,FIDELTA,PREC_J,PREC_S,PREC_TRANS
  !Namelist /species/
  LOGICAL SS_IMP
  !Global
  LOGICAL FAST_IONS,DKES_READ,PHI1_READ,CALCULATED_INT,PRE_INTV,QN,USE_B0,CONVERGED,TRACE_IMP,PLATEAU_OR_PS
  LOGICAL IMP1NU,MODEL_ALPHAD,CLOSEST_LAMBDA,EXTRA_ALPHA,ONE_ALPHA,CENTERED_ALPHA,SECOND_ORDER_ALPHA,RE_SOURCE
  LOGICAL CALC_RHS,CALC_DA,CALC_DIFF,CALC_COL,CALC_DG
  LOGICAL FLUX_NU,PREC_TOP,NEW_DALPHA,INT_G_NEW,STELL_ANTISYMMETRIC,DR_READ!,KROOK_OP
  LOGICAL, PARAMETER :: PLOT_XYZ=.FALSE.

  !Global configuration constants
  INTEGER helN,helM,sgnhel
  REAL*8 rad_a,rad_R,atorflux,torflux,sgnB
  INTEGER, PARAMETER :: mpolbd=128
  INTEGER, PARAMETER :: ntorbd=128
  INTEGER, PARAMETER :: mpold = 1100
  INTEGER, PARAMETER :: ntord = 1100
  INTEGER, PARAMETER :: nmd=2*mpolbd*ntorbd
  INTEGER mpolb,ntorb,nzperiod,Nnm,Nnmp
  INTEGER, PARAMETER :: rprec   = SELECTED_REAL_KIND(12,100)
  REAL*8 mp(nmd),np(nmd),ext_mp(nmd),ext_np(nmd)

  !Flux-surface constants
  REAL*8 eps,eps32,avB,avB2,etet
  REAL*8 psip,chip,iota,siota,aiota,iota2,diotadpsi,Bzeta,Btheta
  REAL*8 iBtpBz,aiBtpBz,dBzdpsi,dBtdpsi,dB0dpsi
  REAL*8 Bmax,Bmax_av,Bmin,Bmin_av
  REAL*8 spol,dspolds
  !Magnetic field map
  REAL(rprec) borbi(-ntorbd:ntorbd,0:mpolbd)
  REAL(rprec) borbic(-ntorbd:ntorbd,0:mpolbd)  ,borbis(-ntorbd:ntorbd,0:mpolbd)
  REAL(rprec) borbic0(-ntorbd:ntorbd,0:mpolbd),borbis0(-ntorbd:ntorbd,0:mpolbd)
  REAL(rprec) dborbicdpsi(-ntorbd:ntorbd,0:mpolbd),dborbisdpsi(-ntorbd:ntorbd,0:mpolbd)
  REAL(rprec) dborbic0dpsi(-ntorbd:ntorbd,0:mpolbd),dborbis0dpsi(-ntorbd:ntorbd,0:mpolbd)   
  REAL*8 rorbic(-ntorbd:ntorbd,0:mpolbd),rorbis(-ntorbd:ntorbd,0:mpolbd)
  REAL*8 zorbic(-ntorbd:ntorbd,0:mpolbd),zorbis(-ntorbd:ntorbd,0:mpolbd)
  REAL*8 porbic(-ntorbd:ntorbd,0:mpolbd),porbis(-ntorbd:ntorbd,0:mpolbd)
  REAL*8 bnmc(nmd),bnmc0(nmd),bnmc1(nmd),dbnmcdpsi(nmd)
  REAL*8 bnms(nmd),bnms0(nmd),bnms1(nmd),dbnmsdpsi(nmd)
  REAL*8, ALLOCATABLE :: absnablar(:,:)
  REAL*8, ALLOCATABLE ::  posx(:,:), posy(:,:), posz(:,:)
  REAL*8, ALLOCATABLE :: zoomx(:,:),zoomy(:,:),zoomz(:,:),zoomdr(:,:)
  !New
  INTEGER nfp_b,ns_b,mnboz_b,mboz_b,nboz_b,jsn
  REAL*8 dpsi
  INTEGER, ALLOCATABLE :: ixm_b(:),ixn_b(:),js_b(:)
  REAL*8,  ALLOCATABLE :: Bmax_b(:),Bmin_b(:)
  REAL*8,  ALLOCATABLE :: s_b(:),spol_b(:),iota_b(:),pres_b(:),beta_b(:),psip_b(:),psi_b(:),bvco_b(:),buco_b(:)
  REAL*8,  ALLOCATABLE :: bmnc_b(:,:),rmnc_b(:,:),pmns_b(:,:),zmns_b(:,:),gmnc_b(:,:)
  REAL*8,  ALLOCATABLE :: bmns_b(:,:),rmns_b(:,:),pmnc_b(:,:),zmnc_b(:,:)
  

  !Grid, resolution and experimetal points
  REAL*8 zmax,tmax,dzstep
  INTEGER, PARAMETER :: narray=10
  INTEGER, PARAMETER :: nparray=100
  REAL, PARAMETER ::   zeta0_DR=1.14
  REAL, PARAMETER ::  dzeta0_DR=0.04
  REAL, PARAMETER ::  theta0_DR=1.88
  REAL, PARAMETER :: dtheta0_DR=2.00
  REAL*8 zetaDR(narray,nparray),thetaDR(narray,nparray)
!  INTEGER array(nax,nax)

  !FFTs
  INTEGER*8, SAVE :: plan_fwd=0
  INTEGER*8, SAVE :: plan_bwd=0

  !Velocity integral
  INTEGER, PARAMETER :: nv=28
  INTEGER, PARAMETER :: iv0=6   !Representative velocity (v~=v_th) 
  REAL*8 v(nv),weight(nv),Sdke(nv),vdconst(nv),vmconst(nv),fdkes(nv),fdkes2(nv),nu(nv),mu(nv),ftrace1nu(nv),mmuoT(nv)
  REAL*8 nuth(nbx),vth(nbx),nuzi(nbx)

  !Collisionality and normalized radial electric field
  REAL*8 cmul_PS,cmul_plateauPS,cmul_1NU
  REAL*8 D11onu,D11pla,D11nu

  !DKES/neotransp/PENTA-related variables
  INTEGER ncmult,nefieldt,nvmagt
  REAL*8, ALLOCATABLE ::  cmult(:), efieldt(:), vmagt(:)
  REAL*8, ALLOCATABLE :: lcmult(:),lefieldt(:),lvmagt(:)
  REAL*8, ALLOCATABLE :: lD11dkes1(:,:),lD11dkes2(:,:),lD11tab(:,:,:),D31dkes(:,:)
  
  INTEGER, ALLOCATABLE :: seed(:)

END MODULE GLOBAL
