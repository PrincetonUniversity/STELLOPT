
  
MODULE KNOSOS_STELLOPT_MOD

  LOGICAL KN_STELLOPT(10),KNOSOS_STELLOPT
  CHARACTER*64 KN_DKESFILE,KN_EXT
  INTEGER KNOSOS_rad_dex
  REAL*8 dpsidr
  REAL*8, ALLOCATABLE :: KNOSOS_1NU(:),KNOSOS_SNU(:),KNOSOS_SBP(:),KNOSOS_GMC(:),KNOSOS_GMA(:)
  REAL*8, ALLOCATABLE :: KNOSOS_QER(:),KNOSOS_VBM(:),KNOSOS_VB0(:),KNOSOS_VBB(:),KNOSOS_WBW(:),KNOSOS_DBO(:)
!  REAL*8, ALLOCATABLE :: KNONOS_ERSHEAR(:),KNOSOS_DJRDER(:)
  REAL*8 KN_1NU,KN_SNU,KN_SBP,KN_GMC,KN_GMA,KN_QER,KN_VBM,KN_VB0,KN_VBB,KN_WBW,KN_DBO,KN_FTR
  
END MODULE KNOSOS_STELLOPT_MOD