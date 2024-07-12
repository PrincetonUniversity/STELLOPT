
      module parambs
      use vmec0
      use bootsj_input
      real(rprec), parameter :: pi = 3.141592653589793238462643_dp
      real(rprec), parameter :: dmu0 = 4.0e-7_dp*pi
      integer, parameter :: nlambda = 101       !this now only governs the trapped calcualtion
      integer :: nthetah, nzetah                !poloidal, toridal grid parameters
      integer, parameter :: nthetahm = 32       !poloidal, toroidal  grid refinement
      integer, parameter :: nzetahm = 16        !for bmax
      integer, parameter :: iotasign = 1
      real(rprec), parameter :: zetasign = -1
      integer :: irdim, irup
      integer, dimension(:), allocatable :: jlist, jlist_idx
      real(rprec)  delta_rho_1
      real(rprec), dimension(:), allocatable :: flux,
     1    aiogar, aipsi, gpsi, pres1,
     2    betar, dense, densi, tempe1, tempi1
      integer, dimension(:), allocatable :: idx
      real(rprec), dimension(:,:), allocatable ::
     1  bfield, gsqrt_b,  b2obm, omb32, bfieldm,    !  arrays to store quantities on the
     2  sinmi, sinnj, cosmi, cosnj,                 !  theta phi grid
     2  sinmim, sinnjm, cosmim, cosnjm              !  quantities that are on the theata

      real(rprec), dimension(:), allocatable ::
     1    dibs, aibs, dibst, aibst, phip,
     2    bsdense, bsdensi, bstempe, bstempi, bsdenste, bsdensti,
     3    bstempte, bstempti, qsafety, capr, caps, h2,
     4    ftrapped, fpassing, epsttok, fttok, b2avg,
     5    gbsnorm, aiterm1, other1, ajBbs,
     6    rhoar, bsnorm, fptok, amain, d_rho,
     7    bmax1, thetamax, zetahmax
      real(rprec)   alphae, alphai , psimax
      real(rprec)   temperho1, tempirho1, densrho1
      real(rprec), dimension(:,:,:), allocatable :: amnfit
      real(rprec), dimension(:,:,:), allocatable :: amnfit2
      real(rprec)   periods
      real(rprec), dimension(:), allocatable :: theta, zetah
      real(rprec), dimension(nthetahm) :: thetam
      real(rprec), dimension(nzetahm) :: zetahm
      complex(rprec), dimension(:,:), allocatable ::
     1    dmn, fmn, alpha1mn
      real(rprec), dimension(:,:), allocatable :: rfmn
      real(rprec)   avgbobm2, sum_gsqrt_b
      real(rprec), dimension(136) :: w1, alamb_h
CCCCfix needed--put nlambda variables into the trapped fraction calculation
      real(rprec)   dlambda, drho
      real(rprec)  sign_jacobian
      logical, dimension(:), allocatable :: lsurf
      logical l_boot_all, lscreen, lasym_bootsj
      end module parambs
