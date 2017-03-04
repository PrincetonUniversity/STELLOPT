      MODULE descriptor_mod
      USE stel_kinds
      INTEGER, PARAMETER ::  BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,  &
                         CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,     &
                         RSRC_ = 7, CSRC_ = 8, LLD_ = 9
      INTEGER, PARAMETER :: nrhs1 = 1
      INTEGER :: iam, nprocs

      REAL(rprec), TARGET,ALLOCATABLE, DIMENSION(:,:) :: ublkp, dblkp, lblkp
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: tempp
      INTEGER :: icontxt_global, icontxt, icontxt_px1, icontxt_1xp
      INTEGER, DIMENSION(DLEN_) :: descA, descX, descR, descA_1xp, descR_all
      INTEGER :: nprow,npcol,myrow,mycol
      INTEGER :: ineed,ineedR,mblk_size2,Locp,Locq,LocpR,LocqR,lld,info,ierr
      INTEGER :: mb,nb, rsrc,csrc,ir,jr,mm,mm0,nrhs0
      LOGICAL :: isroot
      LOGICAL :: LSCALAPACK

      LOGICAL :: DIAGONALDONE, INHESSIAN
      LOGICAL :: ININITSTATE=.FALSE.
      LOGICAL :: INUPDATEFORCE=.FALSE.
      INTEGER :: GMRESPASS, WRAPPASS, MATVECPASS
      INTEGER :: in_hess_nfunct=0, out_hess_nfunct=0
      END MODULE descriptor_mod
