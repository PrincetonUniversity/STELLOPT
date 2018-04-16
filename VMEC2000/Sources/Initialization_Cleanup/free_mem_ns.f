      SUBROUTINE free_mem_ns_par(lreset)
      USE vmec_main
      USE realspace
      USE vforces
      USE vsvd
      USE xstuff
      USE csplinx
      USE fbal
#if defined(SKS)
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
      IF (ALLOCATED(pshalf)) DEALLOCATE(pshalf)
      IF (ALLOCATED(psqrts)) DEALLOCATE(psqrts)
      IF (ALLOCATED(pwint)) DEALLOCATE(pwint)
      IF (ALLOCATED(pfaclam)) DEALLOCATE(pfaclam)
      IF (ALLOCATED(pchip)) DEALLOCATE (pchip)
      IF (ALLOCATED(pphip)) DEALLOCATE (pphip)
      IF (ALLOCATED(pfaclam)) DEALLOCATE (pfaclam)

      IF (ALLOCATED(pxc) .and. lreset) DEALLOCATE (pxc)
      IF (ALLOCATED(pscalxc) .and. lreset) DEALLOCATE (pscalxc)
      IF (ALLOCATED(pxstore)) DEALLOCATE (pxstore)
      IF (ALLOCATED(pxcdot)) DEALLOCATE (pxcdot)
      IF (ALLOCATED(pxsave)) DEALLOCATE (pxsave)
      IF (ALLOCATED(pgc)) DEALLOCATE (pgc)
      IF (ALLOCATED(pcol_scale)) DEALLOCATE(pcol_scale)
#endif
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

      IF (ALLOCATED(phip))
     1  DEALLOCATE (phip, chip, shalf, sqrts, wint, stat=istat3)

      IF (ALLOCATED(ireflect))
     1  DEALLOCATE (ireflect,indexr,imid, stat=istat4)

      IF (ALLOCATED(current))
     1  DEALLOCATE (current, rm2, vrm2, ovrm2, ochip, presph, presint,
     2      w_ia, w1_ia, u_ia, u1_ia, w_pa, w1_pa, u_pa, u1_pa,
     3      w_ib, w1_ib, u_ib, u1_ib, w_pb, w1_pb, u_pb, u1_pb,
     4      rmid, datamse, qmid, shear, presmid, alfa, curmid,
     5      curint, psimid, ageo, volpsi, isplinef, isplineh,
     6      psplinef, psplineh, phimid, pm, im, stat=istat5)


      IF (ALLOCATED(pmb)) DEALLOCATE (pmb,imb,stat=istat6)

      IF (ALLOCATED(ard))
     1  DEALLOCATE (ard,arm,brd,brm,crd,azd,azm,bzd,bzm, sm,sp,
     2        bmin, bmax,stat=istat7)

      IF (ALLOCATED(iotaf))
     1  DEALLOCATE (iotaf,mass,phi,presf,jcuru,jcurv,jdotb,buco,bvco,
     1     bucof, bvcof, chi,
#ifdef _ANIMEC
     1     phot, pmap, pppr, papr, tpotb, pd,
#endif
     2     bdotgradv,equif,specw,tcon,psi,yellip,yinden,
     3     ytrian,yshift,ygeo,overr,faclam,iotas,phips,chips,pres,vp,
     4     beta_vol, jperp2, jpar2, bdotb, clam, blam, dlam, phipf,
     5     chipf, rru_fac, rzu_fac, frcc_fac, fzsc_fac, icurv, vpphi, 
     6     presgrad, r01, z01, bdamp, stat=istat8)

      IF (ALLOCATED(rmidx))
     1  DEALLOCATE (rmidx,hmidx,wmidx,qmidx,tenmidx,ymidx,y2midx,
     2     stat=istat9)

      IF (ALLOCATED(gc))
     1  DEALLOCATE (gc, xsave, xstore, xcdot, col_scale, stat=istat10)
      IF (ALLOCATED(xc) .AND. lreset) DEALLOCATE (xc, scalxc)

      IF (istat1.ne.0 .or. istat2.ne.0 .or. istat3.ne.0 .or.
     1      istat4.ne.0 .or. istat5.ne.0 .or. istat6.ne.0 .or.
     2      istat7.ne.0 .or. istat8.ne.0 .or. istat9.ne.0 .or.
     3      istat10.ne.0) THEN
          PRINT *,' deallocation problem in free_mem_ns'
          PRINT *,' istat1 = ',istat1,' istat2 = ',istat2
          PRINT *,' istat3 = ',istat3,' istat4 = ',istat4
          PRINT *,' istat5 = ',istat5,' istat6 = ',istat6
          PRINT *,' istat7 = ',istat7,' istat8 = ',istat8
          PRINT *,' istat9 = ',istat9,' istat10= ',istat10
       ENDIF

      END SUBROUTINE free_mem_ns
