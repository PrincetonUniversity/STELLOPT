      SUBROUTINE allocate_ns (linterp, neqs2_old)
      USE vmec_main
      USE vmec_params, ONLY: ntmax
      USE realspace
      USE vsvd
      USE vforces
      USE xstuff
      USE csplinx
      USE mgrid_mod
      USE fbal
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   V a r i a b l e s
C-----------------------------------------------
      INTEGER, INTENT(in) :: neqs2_old
      LOGICAL, INTENT(inout) :: linterp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ndim, nsp1, istat1
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: xc_old, scalxc_old
      REAL(rprec) delr_mse
C-----------------------------------------------
!
!     FIRST STORE COARSE-MESH XC FOR INTERPOLATION
!
      ndim  = 1 + nrzt
      nsp1  = 1 + ns
      delr_mse = zero

!
!     Save old xc, scalxc for possible interpolation or IF iterations restarted on same mesh...
!
      IF (neqs2_old .gt. 0 .and. ALLOCATED(scalxc) .and. linterp) THEN
         ALLOCATE(xc_old(neqs2_old), scalxc_old(neqs2_old), stat=istat1)
         IF (istat1.ne.0) STOP 'allocation error #1 in allocate_ns'
         xc_old(:neqs2_old) = xc(:neqs2_old)
         scalxc_old(:neqs2_old) = scalxc(:neqs2_old)
         IF (lrecon) delr_mse = xc(neqs2_old)
      END IF

!
!     ALLOCATES MEMORY FOR NS-DEPENDENT ARRAYS
!     FIRST BE SURE TO FREE MEMORY PREVIOUSLY ALLOCATED
!
      CALL free_mem_ns (.true.)

      ALLOCATE (phip(ndim), chip(ndim), shalf(ndim), sqrts(ndim), 
     1          wint(ndim), stat=istat1)
      IF (istat1.ne.0) STOP 'allocation error #2 in allocate_ns'
      phip=0; chip=0; shalf=0; sqrts=0; wint=0

      ALLOCATE( ireflect(ns*nzeta), indexr(2*ns) ,imid(2*ns),
     1  stat=istat1)
      IF (istat1.ne.0) STOP 'allocation error #3 in allocate_ns'

      ALLOCATE( current(ns),rm2(ns),vrm2(ns),
     1  ovrm2(ns), ochip(ns), presph(ns), presint(ns),
     2  w_ia(ns), w1_ia(ns), u_ia(ns), u1_ia(ns),
     3  w_pa(ns), w1_pa(ns), u_pa(ns), u1_pa(ns),
     4  w_ib(ns), w1_ib(ns), u_ib(ns), u1_ib(ns),
     5  w_pb(ns), w1_pb(ns), u_pb(ns), u1_pb(ns),
     6  rmid(2*ns),datamse(2*ns),qmid(2*ns),
     7  shear(2*ns),presmid(2*ns),alfa(2*ns),curmid(2*ns),
     8  curint(2*ns),psimid(2*ns),ageo(2*ns),volpsi(2*ns),
     9  isplinef(ns),isplineh(ns),psplinef(ns),psplineh(ns),
     A  phimid(2*ns),pm(ns,0:nobser+nobd),
     B  im(ns,0:nobser+nobd),stat=istat1)
      IF (istat1.ne.0) STOP 'allocation error #4 in allocate_ns'

      IF (nbsets.gt.0)
     1  ALLOCATE( pmb(ns,0:nbcoil_max,nbsets,2),
     2            imb(ns,0:nbcoil_max,nbsets,2),stat=istat1)
      IF (istat1.ne.0) STOP 'allocation error #5 in allocate_ns'

      ALLOCATE( ard(nsp1,2),arm(nsp1,2),brd(nsp1,2),brm(nsp1,2),
     1          azd(nsp1,2),azm(nsp1,2),bzd(nsp1,2), bzm(nsp1,2),
     2          sm(ns), sp(0:ns), bmin(ntheta2,ns), bmax(ntheta2,ns),
     3          stat=istat1)
      IF (istat1.ne.0) STOP 'allocation error #6 in allocate_ns'

      ALLOCATE( iotaf(nsp1), crd(nsp1), mass(ns), phi(ns), presf(ns),
     1          jcuru(ns), jcurv(ns), jdotb(ns), buco(ns), bvco(ns),
#ifdef _ANIMEC
!WAC ANISTROPIC VARIABLES
     1          phot(ns), pmap(ns), pppr(ns), papr(ns), tpotb(ns), 
     1          pd(ns),
#endif
     1          bucof(ns), bvcof(ns), chi(ns),
     2          bdotgradv(ns), equif(ns), specw(ns), tcon(ns), 
     3          psi(ns),yellip(ns),yinden(ns), ytrian(ns),yshift(ns),
     4          ygeo(ns),overr(ns), faclam(ns,0:ntor,0:mpol1,ntmax),
     5          iotas(nsp1), phips(nsp1), chips(nsp1), pres(nsp1), 
     6          beta_vol(ns), jperp2(ns), jpar2(ns), bdotb(ns),
     7          phipf(ns), chipf(ns), blam(nsp1), clam(nsp1), 
     8          dlam(nsp1), rru_fac(ns), rzu_fac(ns), frcc_fac(ns),
     9          fzsc_fac(ns), icurv(ns+1), vpphi(ns), bdamp(ns),
     9          presgrad(ns), vp(nsp1), r01(ns), z01(ns), stat=istat1)

      frcc_fac = 0; fzsc_fac = 0
#ifdef _ANIMEC
      phot=0; tpotb=0
#endif
      IF (istat1.ne.0) STOP 'allocation error #7 in allocate_ns'

      iotaf(nsp1) = 0

      ALLOCATE( rmidx(2*ns), hmidx(2*ns), wmidx(2*ns), qmidx(2*ns),
     1          tenmidx(2*ns), ymidx(2*ns), y2midx(2*ns), stat=istat1)
      IF (istat1.ne.0) STOP 'allocation error #8 in allocate_ns'

      ALLOCATE (gc(neqs2), xcdot(neqs2), xsave(neqs2), 
     1          xstore(neqs2), stat=istat1)
      xstore = zero
      IF (istat1.ne.0) STOP 'allocation error #9 in allocate_ns'

      IF (.not.ALLOCATED(xc)) THEN
         ALLOCATE (xc(neqs2), scalxc(neqs2), stat=istat1)
      IF (istat1.ne.0) STOP 'allocation error #10 in allocate_ns'
         xc(:neqs2) = zero
      END IF

      IF (ALLOCATED(xc_old)) THEN
         xstore(1:neqs2_old) = xc_old(1:neqs2_old)
         scalxc(1:neqs2_old) = scalxc_old(1:neqs2_old)
         DEALLOCATE (xc_old, scalxc_old)
      END IF


      xc(neqs2) = delr_mse

!
!     Allocate nrzt-dependent arrays (persistent) for funct3d
!
      CALL allocate_funct3d

      END SUBROUTINE allocate_ns
