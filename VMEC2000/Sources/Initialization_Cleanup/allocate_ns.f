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
      USE parallel_include_module
      USE vmec_input, ONLY: nzeta
      USE vmec_dim, ONLY: ns, ntheta3

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
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: pxc_old, pscalxc_old
      REAL(rprec) :: delr_mse
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
      IF (PARVMEC) THEN 
         IF (neqs2_old.GT.0 .AND. ALLOCATED(pscalxc) .AND. linterp) THEN
            ALLOCATE(pxc_old(neqs2_old),pscalxc_old(neqs2_old),
     &               stat=istat1)
            IF (istat1.NE.0) THEN
               STOP 'allocation error #1 in allocate_ns'
            ENDIF
            pxc_old(:neqs2_old) = pxc(:neqs2_old)
            pscalxc_old(:neqs2_old) = pscalxc(:neqs2_old)
            IF (lrecon) THEN
               delr_mse = pxc(neqs2_old)
            END IF
         END IF
      END IF

      IF (neqs2_old .GT. 0 .AND. ALLOCATED(scalxc) .AND. linterp) THEN
         ALLOCATE(xc_old(neqs2_old), scalxc_old(neqs2_old), stat=istat1)
         IF (istat1.ne.0) THEN
            STOP 'allocation error #1 in allocate_ns'
         END IF
         xc_old(:neqs2_old) = xc(:neqs2_old)
         scalxc_old(:neqs2_old) = scalxc(:neqs2_old)
         IF (lrecon) THEN
            delr_mse = xc(neqs2_old)
         END IF
      END IF

!
!     ALLOCATES MEMORY FOR NS-DEPENDENT ARRAYS
!     FIRST BE SURE TO FREE MEMORY PREVIOUSLY ALLOCATED
!
      IF (PARVMEC) THEN 
        CALL free_mem_ns_par (.true.)
      END IF
      CALL free_mem_ns (.true.)

      ALLOCATE (phip(ndim), chip(ndim), shalf(ndim), sqrts(ndim), 
     1          wint(ndim), stat=istat1)
      IF (istat1.ne.0) THEN
         STOP 'allocation error #2 in allocate_ns'
      END IF
      phip=0; chip=0; shalf=0; sqrts=0; wint=0

      IF(PARVMEC) THEN
         ALLOCATE(pshalf(nznt,ns),stat=istat1)
         ALLOCATE(pwint(nznt,ns),stat=istat1)
         ALLOCATE(pwint_ns(nznt),stat=istat1)
         ALLOCATE(ireflect_par(nzeta),stat=istat1)
         ALLOCATE(pchip(nznt,ns),stat=istat1)
         ALLOCATE(pphip(nznt,ns),stat=istat1)
         ALLOCATE(psqrts(nznt,ns),stat=istat1)
         ALLOCATE(pfaclam(0:ntor,0:mpol1,1:ns,ntmax),stat=istat1)
      END IF

      ALLOCATE(ireflect(ns*nzeta), indexr(2*ns) ,imid(2*ns),
     &         stat=istat1)
      IF (istat1.ne.0) THEN
         STOP 'allocation error #3 in allocate_ns'
      END IF

      ALLOCATE(current(ns),rm2(ns),vrm2(ns),
     &         ovrm2(ns), ochip(ns), presph(ns), presint(ns),
     &         w_ia(ns), w1_ia(ns), u_ia(ns), u1_ia(ns),
     &         w_pa(ns), w1_pa(ns), u_pa(ns), u1_pa(ns),
     &         w_ib(ns), w1_ib(ns), u_ib(ns), u1_ib(ns),
     &         w_pb(ns), w1_pb(ns), u_pb(ns), u1_pb(ns),
     &         rmid(2*ns),datamse(2*ns),qmid(2*ns),
     &         shear(2*ns),presmid(2*ns),alfa(2*ns),curmid(2*ns),
     &         curint(2*ns),psimid(2*ns),ageo(2*ns),volpsi(2*ns),
     &         isplinef(ns),isplineh(ns),psplinef(ns),psplineh(ns),
     &         phimid(2*ns),pm(ns,0:nobser+nobd),
     &         im(ns,0:nobser+nobd),stat=istat1)
      IF (istat1.ne.0) THEN
         STOP 'allocation error #4 in allocate_ns'
      END IF

      IF (nbsets.gt.0) THEN
         ALLOCATE(pmb(ns,0:nbcoil_max,nbsets,2),
     &            imb(ns,0:nbcoil_max,nbsets,2),stat=istat1)
      END IF
      IF (istat1.ne.0) THEN
         STOP 'allocation error #5 in allocate_ns'
      END IF

      ALLOCATE(ard(nsp1,2),arm(nsp1,2),brd(nsp1,2),brm(nsp1,2),
     &         azd(nsp1,2),azm(nsp1,2),bzd(nsp1,2), bzm(nsp1,2),
     &         sm(ns), sp(0:ns), bmin(ntheta2,ns), bmax(ntheta2,ns),
     &         stat=istat1)
      IF (istat1.ne.0) THEN
         STOP 'allocation error #6 in allocate_ns'
      END IF

      ALLOCATE(iotaf(nsp1), crd(nsp1), mass(ns), phi(ns), presf(ns),
     &         jcuru(ns), jcurv(ns), jdotb(ns), buco(ns), bvco(ns),
#ifdef _ANIMEC
!WAC ANISTROPIC VARIABLES
     &         phot(ns), pmap(ns), pppr(ns), papr(ns), tpotb(ns),
     &         pd(ns),
#endif
     &         bucof(ns), bvcof(ns), chi(ns),
     &         bdotgradv(ns), equif(ns), specw(ns), tcon(ns),
     &         psi(ns),yellip(ns),yinden(ns), ytrian(ns),yshift(ns),
     &         ygeo(ns),overr(ns), faclam(ns,0:ntor,0:mpol1,ntmax),
     &         iotas(nsp1), phips(nsp1), chips(nsp1), pres(nsp1),
     &         beta_vol(ns), jperp2(ns), jpar2(ns), bdotb(ns),
     &         phipf(ns), chipf(ns), blam(nsp1), clam(nsp1),
     &         dlam(nsp1), rru_fac(ns), rzu_fac(ns), frcc_fac(ns),
     &         fzsc_fac(ns), icurv(ns+1), vpphi(ns), bdamp(ns),
     &         presgrad(ns), vp(nsp1), r01(ns), z01(ns), stat=istat1)

      frcc_fac = 0; fzsc_fac = 0
#ifdef _ANIMEC
      phot=0; tpotb=0
#endif
      IF (istat1 .NE. 0) THEN
         STOP 'allocation error #7 in allocate_ns'
      END IF

      iotaf(nsp1) = 0

      ALLOCATE(rmidx(2*ns), hmidx(2*ns), wmidx(2*ns), qmidx(2*ns),
     1         tenmidx(2*ns), ymidx(2*ns), y2midx(2*ns), stat=istat1)
      IF (istat1 .NE. 0) THEN
         STOP 'allocation error #8 in allocate_ns'
      END IF

      IF(PARVMEC) THEN
         ALLOCATE(pgc(neqs2), pxcdot(neqs2), pxsave(neqs2),
     &            pxstore(neqs2), pcol_scale(neqs2), stat=istat1)
         pxstore = zero
         IF (istat1 .NE. 0) THEN
            STOP 'allocation error #9 in allocate_ns'
         END IF

         IF (.not.ALLOCATED(pxc)) THEN
            ALLOCATE (pxc(neqs2), pscalxc(neqs2), stat=istat1)
            IF (istat1 .NE. 0) THEN
               STOP 'allocation error #10 in allocate_ns'
            END IF
            pxc(:neqs2) = zero
         END IF

         IF (ALLOCATED(pxc_old)) THEN
            pxstore(1:neqs2_old) = pxc_old(1:neqs2_old)
            pscalxc(1:neqs2_old) = pscalxc_old(1:neqs2_old)
            DEALLOCATE (pxc_old, pscalxc_old)
         END IF
         pxc(neqs2) = delr_mse
      END IF

      ALLOCATE(gc(neqs2), xcdot(neqs2), xsave(neqs2),
     &         xstore(neqs2), col_scale(neqs2), stat=istat1)
      xstore = zero
      IF (istat1 .NE. 0) THEN
         STOP 'allocation error #9 in allocate_ns'
      END IF

      IF (.NOT.ALLOCATED(xc)) THEN
         ALLOCATE (xc(neqs2), scalxc(neqs2), stat=istat1)
         IF (istat1 .NE. 0) THEN
            STOP 'allocation error #10 in allocate_ns'
         END IF
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
      IF (PARVMEC) THEN
         CALL allocate_funct3d_par
      ELSE
         CALL allocate_funct3d
      END IF

      END SUBROUTINE allocate_ns
