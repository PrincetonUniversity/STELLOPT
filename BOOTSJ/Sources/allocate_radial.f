

      subroutine allocate_radial
      use parambs
      use vmec0
      use read_boozer_mod

      implicit none
      integer :: istat
      INTEGER :: ier_arr(64)

      ! ADDDED BY SAL - 08/06/12
      ier_arr = 0
      IF (ALLOCATED(flux)) DEALLOCATE(flux)
      ALLOCATE(flux(irdim), stat=ier_arr(1))
      IF (ALLOCATED(qsafety)) DEALLOCATE(qsafety)
      ALLOCATE(qsafety(irdim), stat=ier_arr(2))
      IF (ALLOCATED(aiogar)) DEALLOCATE(aiogar)
      ALLOCATE(aiogar(irdim), stat=ier_arr(3))
      IF (ALLOCATED(idx)) DEALLOCATE(idx)
      ALLOCATE(idx(irup), stat=ier_arr(4))
      IF (ALLOCATED(aipsi)) DEALLOCATE(aipsi)
      ALLOCATE(aipsi(irdim), stat=ier_arr(5))
      IF (ALLOCATED(gpsi)) DEALLOCATE(gpsi)
      ALLOCATE(gpsi(irdim), stat=ier_arr(6))
      IF (ALLOCATED(pres1)) DEALLOCATE(pres1)
      ALLOCATE(pres1(irdim), stat=ier_arr(7))
      IF (ALLOCATED(betar)) DEALLOCATE(betar)
      ALLOCATE(betar(irdim), stat=ier_arr(8))
      IF (ALLOCATED(dense)) DEALLOCATE(dense)
      ALLOCATE(dense(irdim), stat=ier_arr(9))
      IF (ALLOCATED(densi)) DEALLOCATE(densi)
      ALLOCATE(densi(irdim), stat=ier_arr(10))
      IF (ALLOCATED(tempe1)) DEALLOCATE(tempe1)
      ALLOCATE(tempe1(irdim), stat=ier_arr(11))
      IF (ALLOCATED(tempi1)) DEALLOCATE(tempi1)
      ALLOCATE(tempi1(irdim), stat=ier_arr(12))
      IF (ALLOCATED(lsurf)) DEALLOCATE(lsurf)
      ALLOCATE(lsurf(irdim), stat=ier_arr(13))
      IF (ALLOCATED(jlist)) DEALLOCATE(jlist)
      ALLOCATE(jlist(irdim), stat=ier_arr(14))
      IF (ALLOCATED(jlist_idx)) DEALLOCATE(jlist_idx)
      ALLOCATE(jlist_idx(irdim), stat=ier_arr(15))
      IF(ANY(ier_arr.ne.0)) STOP 'allocation error1 in allocate_radial'
!      allocate (flux(irdim), qsafety(irdim), aiogar(irdim),
!     1    idx(irup), aipsi(irdim), gpsi(irdim), pres1(irdim),
!     2    betar(irdim), dense(irdim), densi(irdim),
!     3    tempe1(irdim), tempi1(irdim), lsurf(irdim),
!     4    jlist(irdim), jlist_idx(irdim), stat=istat)
!      if (istat .ne. 0) return

      IF (ALLOCATED(amnfit)) DEALLOCATE(amnfit)
      allocate (amnfit(irup,-mboz_b:mboz_b,0:nboz_b), stat=istat)
      if (istat .ne. 0) STOP 'allocation error2 in allocate_radial'
      IF (lasym_b) THEN
         IF (ALLOCATED(amnfit2)) DEALLOCATE(amnfit2)
         allocate (amnfit2(irup,-mboz_b:mboz_b,0:nboz_b), stat=istat)
         if (istat .ne. 0) STOP 'allocation error2a in allocate_radial'
      END IF

      ier_arr = 0
      IF (ALLOCATED(dibs)) DEALLOCATE(dibs)
      ALLOCATE(dibs(irup), stat=ier_arr(1))
      IF (ALLOCATED(aibs)) DEALLOCATE(aibs)
      ALLOCATE(aibs(irup), stat=ier_arr(2))
      IF (ALLOCATED(dibst)) DEALLOCATE(dibst)
      ALLOCATE(dibst(irup), stat=ier_arr(3))
      IF (ALLOCATED(aibst)) DEALLOCATE(aibst)
      ALLOCATE(aibst(irup), stat=ier_arr(4))
      IF (ALLOCATED(bsdense)) DEALLOCATE(bsdense)
      ALLOCATE(bsdense(irup), stat=ier_arr(5))
      IF (ALLOCATED(bsdensi)) DEALLOCATE(bsdensi)
      ALLOCATE(bsdensi(irup), stat=ier_arr(6))
      IF (ALLOCATED(bstempe)) DEALLOCATE(bstempe)
      ALLOCATE(bstempe(irup), stat=ier_arr(7))
      IF (ALLOCATED(bstempi)) DEALLOCATE(bstempi)
      ALLOCATE(bstempi(irup), stat=ier_arr(8))
      IF (ALLOCATED(bsdenste)) DEALLOCATE(bsdenste)
      ALLOCATE(bsdenste(irup), stat=ier_arr(9))
      IF (ALLOCATED(bsdensti)) DEALLOCATE(bsdensti)
      ALLOCATE(bsdensti(irup), stat=ier_arr(10))
      IF (ALLOCATED(bstempte)) DEALLOCATE(bstempte)
      ALLOCATE(bstempte(irup), stat=ier_arr(11))
      IF (ALLOCATED(bstempti)) DEALLOCATE(bstempti)
      ALLOCATE(bstempti(irup), stat=ier_arr(12))
      IF (ALLOCATED(capr)) DEALLOCATE(capr)
      ALLOCATE(capr(irup), stat=ier_arr(13))
      IF (ALLOCATED(caps)) DEALLOCATE(caps)
      ALLOCATE(caps(irup), stat=ier_arr(14))
      IF (ALLOCATED(h2)) DEALLOCATE(h2)
      ALLOCATE(h2(irup), stat=ier_arr(15))
      IF (ALLOCATED(ftrapped)) DEALLOCATE(ftrapped)
      ALLOCATE(ftrapped(irup), stat=ier_arr(16))
      IF (ALLOCATED(fpassing)) DEALLOCATE(fpassing)
      ALLOCATE(fpassing(irup), stat=ier_arr(17))
      IF (ALLOCATED(epsttok)) DEALLOCATE(epsttok)
      ALLOCATE(epsttok(irup), stat=ier_arr(18))
      IF (ALLOCATED(fttok)) DEALLOCATE(fttok)
      ALLOCATE(fttok(irup), stat=ier_arr(19))
      IF (ALLOCATED(gbsnorm)) DEALLOCATE(gbsnorm)
      ALLOCATE(gbsnorm(irup), stat=ier_arr(20))
      IF (ALLOCATED(aiterm1)) DEALLOCATE(aiterm1)
      ALLOCATE(aiterm1(irup), stat=ier_arr(21))
      IF (ALLOCATED(other1)) DEALLOCATE(other1)
      ALLOCATE(other1(irup), stat=ier_arr(22))
      IF (ALLOCATED(rhoar)) DEALLOCATE(rhoar)
      ALLOCATE(rhoar(irup), stat=ier_arr(23))
      IF (ALLOCATED(bsnorm)) DEALLOCATE(bsnorm)
      ALLOCATE(bsnorm(irup), stat=ier_arr(24))
      IF (ALLOCATED(fptok)) DEALLOCATE(fptok)
      ALLOCATE(fptok(irup), stat=ier_arr(25))
      IF (ALLOCATED(amain)) DEALLOCATE(amain)
      ALLOCATE(amain(irup), stat=ier_arr(26))
      IF (ALLOCATED(bmax1)) DEALLOCATE(bmax1)
      ALLOCATE(bmax1(irup), stat=ier_arr(27))
      IF (ALLOCATED(thetamax)) DEALLOCATE(thetamax)
      ALLOCATE(thetamax(irup), stat=ier_arr(28))
      IF (ALLOCATED(zetahmax)) DEALLOCATE(zetahmax)
      ALLOCATE(zetahmax(irup), stat=ier_arr(29))
      IF (ALLOCATED(ajBbs)) DEALLOCATE(ajBbs)
      ALLOCATE(ajBbs(irup), stat=ier_arr(30))
      IF (ALLOCATED(phip)) DEALLOCATE(phip)
      ALLOCATE(phip(irup), stat=ier_arr(31))
      IF (ALLOCATED(d_rho)) DEALLOCATE(d_rho)
      ALLOCATE(d_rho(irdim), stat=ier_arr(32))
      IF (ALLOCATED(b2avg)) DEALLOCATE(b2avg)
      ALLOCATE(b2avg(irup), stat=ier_arr(33))
      IF (ANY(ier_arr.ne.0)) stop 'allocation error3 in allocate_radial'
!      allocate (dibs(irup), aibs(irup), dibst(irup), aibst(irup),
!     1    bsdense(irup), bsdensi(irup), bstempe(irup), bstempi(irup),
!     2    bsdenste(irup), bsdensti(irup), bstempte(irup),
!     3    bstempti(irup), capr(irup),
!     4    caps(irup), h2(irup), ftrapped(irup), fpassing(irup),
!     5    epsttok(irup), fttok(irup), gbsnorm(irup), aiterm1(irup),
!     6    other1(irup),
!     7    rhoar(irup), bsnorm(irup), fptok(irup), amain(irup),
!     8    bmax1(irup), thetamax(irup), zetahmax(irup),
!     9    ajBbs(irup), phip(irup), d_rho(irdim),
!     A    b2avg(irup), stat = istat)

!      if (istat .ne. 0) stop 'allocation error in allocate_radial'

      end subroutine allocate_radial
