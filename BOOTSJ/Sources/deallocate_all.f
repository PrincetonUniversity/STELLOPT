

      subroutine deallocate_all
      use parambs
      use vmec0
      implicit none

      if (allocated(amnfit)) deallocate (amnfit)
      if (allocated(dibs)) deallocate (dibs, aibs, dibst, aibst,
     1    bsdense, bsdensi, bstempe, bstempi, bsdenste, bsdensti,
     2    bstempte, bstempti, qsafety, capr, caps, h2,
     3    ftrapped, fpassing, epsttok, fttok, gbsnorm, aiterm1,
     4    other1, rhoar, bsnorm, fptok, amain,
     5    bmax1, thetamax, zetahmax, d_rho, b2avg,
     6    ajBbs, phip)
      IF (ALLOCATED(theta)) DEALLOCATE(theta)
      IF (ALLOCATED(zetah)) DEALLOCATE(zetah)
      if (allocated(flux)) deallocate (flux, aiogar, aipsi,
     1    gpsi, pres1, betar, dense, densi, tempe1, tempi1, lsurf,
     2    jlist, jlist_idx)

      IF (ALLOCATED(idx)) deallocate(idx)
      IF (ALLOCATED(bfield)) DEALLOCATE(bfield)
      IF (ALLOCATED(b2obm)) DEALLOCATE(b2obm)
      IF (ALLOCATED(gsqrt_b)) DEALLOCATE(gsqrt_b)
      IF (ALLOCATED(sinmi)) DEALLOCATE(sinmi)
      IF (ALLOCATED(sinnj)) DEALLOCATE(sinnj)
      IF (ALLOCATED(cosmi)) DEALLOCATE(cosmi)
      IF (ALLOCATED(cosnj)) DEALLOCATE(cosnj)
      IF (ALLOCATED(bfieldm)) DEALLOCATE(bfieldm)
      IF (ALLOCATED(sinmim)) DEALLOCATE(sinmim)
      IF (ALLOCATED(sinnjm)) DEALLOCATE(sinnjm)
      IF (ALLOCATED(cosmim)) DEALLOCATE(cosmim)
      IF (ALLOCATED(cosnjm)) DEALLOCATE(cosnjm)
!      deallocate(b2obm, gsqrt_b, sinmi, sinnj, cosmi, cosnj,
!     1  bfieldm, sinmim, sinnjm, cosmim, cosnjm)

      end subroutine deallocate_all
