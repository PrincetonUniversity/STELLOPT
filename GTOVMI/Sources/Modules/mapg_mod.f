      MODULE mapg_mod
      USE precision
      IMPLICIT NONE
c =====  DIMENSIONal parameters =====
      INTEGER, PARAMETER :: npsi=513,  nthet=257
!npsi              !SIZE of mapped psi mesh
!nthet             !SIZE of mapped theta mesh
! these are not related to the gfile DIMENSIONs which are ALLOCATED
c =====  DIMENSIONal parameters =====
      INTEGER, PARAMETER :: kpsi=npsi, kthet=nthet, kv=kthet+2
c
      INTEGER :: nxd
      INTEGER :: nzd
      INTEGER :: nxzd
      INTEGER :: nh2
      INTEGER :: nwk
      INTEGER :: kbnd
      INTEGER :: klim
      parameter(kbnd=300,klim=300)
c =====  i/o channels =====
      INTEGER :: ninbal            !io_unit for inbal
      INTEGER :: nterm             !io_unit for terminal
      INTEGER :: noutbal           !io_unit for outbal
      INTEGER :: ndskbal           !io_unit for dskbal
      INTEGER :: neqdsk=0            !io_unit for g-eqdsk
      CHARACTER*40 filename     !name of eqdsk or dskbal
c =====  input parameters for equilibrium data =====
      LOGICAL READeqdsk         !flag for toq equilibrium
      INTEGER :: npfit             !threshhold pts for using furpl
      REAL(rprec) alpsi                !flux prop. psic^alpsi
      INTEGER :: rotate            !=1 for rotation to READ eqdsk correctly
      REAL(rprec) percenflux           !outermost flux surface is 
                                !psiv(npsi)=psiaxis+(psilim-psiaxis)*percenflux
      CHARACTER*6 contour
      REAL(rprec) epsarc
      REAL(rprec) pi                   !3.14159...

c =====  geqdsk variabls READ in for mapper =====
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: xgrid           !x array
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: zgrid           !z array
      REAL(rprec) :: xaxis                !x of mag axis
      REAL(rprec) :: zaxis                !z of mag axis
      REAL(rprec) :: psiaxis              !psi at mag axis
      REAL(rprec) :: psilim               !psi at limiter
      REAL(rprec) :: dx                   !x(2)-x(1)
      REAL(rprec) :: dz                   !z(2)-z(1)
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: psixz         !complete psi array
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: bpsq        !calculated in READeqdsk
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: sp              !pressure on equally spaced psi
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: spp             !pprime on equally spaced psi
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: sf              !f on equally spaced psi
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: sffp            !ffprime on equally spaced psi
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: qpsi            !qpsi on equally spaced psi
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: xbndry         !x position of boundary points
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: zbndry         !z position of boundary points
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: xlim           !x position of limiter points
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: zlim           !z position of limiter points
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: pressw          !see Lang's defn for press with rotation
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: pwprim          !see Lang's defn for press with rotation
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: rho0            !see Lang's defn for press with rotation
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: rho0p           !d(rho0)/d(psi)
      REAL(rprec) rvtor                !reference r for rotation
      INTEGER nbndry            !number of boundary points
      INTEGER nlim              !number of limiter points
      INTEGER kvtor             !>0 means rotation
      INTEGER nmass             !>0 means rho0 specified
      INTEGER nx                !SIZE of x mesh
      INTEGER nz                !SIZE of z mesh
c ===== variables READ in from dskbal or obtained from mapper=====
      REAL(rprec) psic(kpsi)           !the psi coordinates
      REAL(rprec) psiv(kpsi)           !the REAL poloidal flux array (a.k.a. chi)
      REAL(rprec) pprime(kpsi)         !pprime
      REAL(rprec) fval(kpsi)           !this is f as in Btor=f/R
      REAL(rprec) ffprime(kpsi)        !this is ff'
      REAL(rprec) chipsi(kpsi)         !d(psiv)/d(psic)
      REAL(rprec) qsfin(kpsi)          !the safety factor from input
      REAL(rprec) xs(kpsi,kthet)       !R-coordinate of pts. on psiv contour
      REAL(rprec) zs(kpsi,kthet)       !Z-coordinate
      REAL(rprec) bps(kpsi,kthet)    !bp in flux coordinates
      REAL(rprec) arcsur(kpsi,kthet)   !mapped arc LENgths
      REAL(rprec) arcrad(kpsi,kthet)
      REAL(rprec) press(kpsi)
      REAL(rprec) pw(kpsi)
      REAL(rprec) pwp(kpsi)
      REAL(rprec) rho(kpsi)
      REAL(rprec) rhop(kpsi)
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: csplpsi
      REAL(rprec) dthe                 !spacing between theta points
      REAL(rprec) dpsi                 !spacing between psic surfaces
      REAL(rprec) jtor(kpsi,kthet)
      REAL(rprec) eikon(kpsi,kthet)
      REAL(rprec) qsf(kpsi)
      REAL(rprec) thec(kthet)
      REAL(rprec) vprime(kpsi)
      REAL(rprec) volume(kpsi), wbav(kpsi), tflx(kpsi)
      REAL(rprec) bsqav(kpsi) 
      REAL(rprec) itor(kpsi)
      REAL(rprec) jacob(kpsi,kthet)
      REAL(rprec) jovr(kpsi), iprime(kpsi)
      REAL(rprec) kappapav(kpsi)
      REAL(rprec) fracg(kpsi)
      INTEGER iphi
      CONTAINS
        REAL FUNCTION trap(n,x,y)
          INTEGER, INTENT(IN) :: n
          REAL*8, INTENT(IN) :: x(n), y(n)
          REAL*8 :: summ=0., half=0.5
          INTEGER :: j
          summ=0.
          DO j=2,n
            summ=summ+half*(y(j)+y(j-1))*(x(j)-x(j-1))
          ENDdo
          trap=summ
        END FUNCTION trap
        REAL FUNCTION axisv(x,y,n)
          IMPLICIT NONE
          INTEGER, INTENT(in) :: n
          REAL(rprec), INTENT(in) :: x(n), y(n)
          axisv=y(2)+(y(3)-y(2))*(x(2)-x(1))/(x(3)-x(2))
        END FUNCTION axisv
      END MODULE mapg_mod
