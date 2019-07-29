      SUBROUTINE free_mem_ns_par(lreset)
      USE vmec_main
      USE realspace
      USE vforces
      USE vsvd
      USE xstuff
      USE csplinx
      USE fbal

      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      LOGICAL, INTENT(in) :: lreset

C-----------------------------------------------
      IF (ALLOCATED(pshalf))  DEALLOCATE(pshalf)
      IF (ALLOCATED(psqrts))  DEALLOCATE(psqrts)
      IF (ALLOCATED(pwint))   DEALLOCATE(pwint)
      IF (ALLOCATED(pfaclam)) DEALLOCATE(pfaclam)
      IF (ALLOCATED(pchip))   DEALLOCATE(pchip)
      IF (ALLOCATED(pphip))   DEALLOCATE(pphip)
      IF (ALLOCATED(pfaclam)) DEALLOCATE(pfaclam)

      IF (ALLOCATED(pxc)     .and. lreset) DEALLOCATE(pxc)
      IF (ALLOCATED(pscalxc) .and. lreset) DEALLOCATE (pscalxc)

      IF (ALLOCATED(pxstore))    DEALLOCATE (pxstore)
      IF (ALLOCATED(pxcdot))     DEALLOCATE (pxcdot)
      IF (ALLOCATED(pxsave))     DEALLOCATE (pxsave)
      IF (ALLOCATED(pgc))        DEALLOCATE (pgc)
      IF (ALLOCATED(pcol_scale)) DEALLOCATE(pcol_scale)

      END SUBROUTINE free_mem_ns_par

      SUBROUTINE free_mem_ns(lreset)
      USE vmec_main
      USE realspace
      USE vforces
      USE vsvd
      USE xstuff
      USE csplinx
      USE fbal

      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      LOGICAL, INTENT(in) :: lreset
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: istat1 = 0, istat2 = 0, istat3 = 0, istat4 = 0,
     1           istat5 = 0, istat6 = 0, istat7 = 0, istat8 = 0,
     2           istat9 = 0, istat10 = 0
C-----------------------------------------------

      IF (ALLOCATED(phip)) THEN
         DEALLOCATE(phip, chip, shalf, sqrts, wint, stat=istat1)
      END IF

      IF (ALLOCATED(ireflect)) THEN
         DEALLOCATE(ireflect, indexr, imid, stat=istat2)
      END IF

      IF (ALLOCATED(current)) THEN
         DEALLOCATE(current, rm2, vrm2, ovrm2, ochip, presph, presint,
     &              w_ia, w1_ia, u_ia, u1_ia, w_pa, w1_pa, u_pa, u1_pa,
     &              w_ib, w1_ib, u_ib, u1_ib, w_pb, w1_pb, u_pb, u1_pb,
     &              rmid, datamse, qmid, shear, presmid, alfa, curmid,
     &              curint, psimid, ageo, volpsi, isplinef, isplineh,
     &              psplinef, psplineh, phimid, pm, im, stat=istat3)
      END IF

      IF (ALLOCATED(pmb)) DEALLOCATE(pmb, imb, stat=istat4)

      IF (ALLOCATED(ard)) THEN
         DEALLOCATE(ard, arm, brd, brm, crd, azd, azm, bzd, bzm, sm,
     &              sp, bmin, bmax, stat=istat5)
      END IF

      IF (ALLOCATED(iotaf)) THEN
         DEALLOCATE(iotaf, mass, phi, presf, jcuru, jcurv, jdotb, buco,
     &              bvco, bucof, bvcof, chi,
#ifdef _ANIMEC
     &              phot, pmap, pppr, papr, tpotb, pd,
#endif
     &              bdotgradv, equif, specw, tcon, psi, yellip, yinden,
     &              ytrian, yshift, ygeo, overr, faclam, iotas, phips,
     &              chips, pres, vp, beta_vol, jperp2, jpar2, bdotb,
     &              clam, blam, dlam, phipf, chipf, rru_fac, rzu_fac,
     &              frcc_fac, fzsc_fac, icurv, vpphi, presgrad, r01,
     &              z01, bdamp, stat=istat6)
      END IF

      IF (ALLOCATED(rmidx)) THEN
         DEALLOCATE(rmidx, hmidx, wmidx, qmidx, tenmidx, ymidx, y2midx,
     &              stat=istat7)
      END IF

      IF (ALLOCATED(gc)) THEN
         DEALLOCATE(gc, xsave, xstore, xcdot, col_scale, stat=istat8)
      END IF
      IF (ALLOCATED(xc) .AND. lreset) THEN
         DEALLOCATE(xc, scalxc)
      END IF

      IF ((istat1 + istat2 + istat3 + istat4 + istat5 + istat6 +
     &     istat7 + istat8) .ne. 0) THEN
          PRINT *,' deallocation problem in free_mem_ns'
          PRINT *,' istat1 = ',istat1,' istat2 = ',istat2
          PRINT *,' istat3 = ',istat3,' istat4 = ',istat4
          PRINT *,' istat5 = ',istat5,' istat6 = ',istat6
          PRINT *,' istat7 = ',istat7,' istat8 = ',istat8
       END IF

      END SUBROUTINE free_mem_ns
