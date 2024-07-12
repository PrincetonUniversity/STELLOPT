      SUBROUTINE free_mem_boozer
      USE booz_params
      USE booz_persistent
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: istat1=0, istat2=0, istat3=0
C-----------------------------------------------
      IF (ALLOCATED(jlist)) DEALLOCATE(jlist)
      IF (ALLOCATED(lsurf_boz)) DEALLOCATE(lsurf_boz)
      !IF (ALLOCATED(jlist)) DEALLOCATE (jlist, lsurf_boz)

      IF (ALLOCATED(bsubumnc)) DEALLOCATE(bsubumnc)
      IF (ALLOCATED(bsubvmnc)) DEALLOCATE(bsubvmnc)
      IF (ALLOCATED(bmodmnc)) DEALLOCATE(bmodmnc)
      IF (ALLOCATED(rmnc)) DEALLOCATE(rmnc)
      IF (ALLOCATED(zmns)) DEALLOCATE(zmns)
      IF (ALLOCATED(lmns)) DEALLOCATE(lmns)
      IF (ALLOCATED(xm)) DEALLOCATE(xm)
      IF (ALLOCATED(xn)) DEALLOCATE(xn)
      IF (ALLOCATED(xm_nyq)) DEALLOCATE(xm_nyq)
      IF (ALLOCATED(xn_nyq)) DEALLOCATE(xn_nyq)
      IF (ALLOCATED(hiota)) DEALLOCATE(hiota)
      IF (ALLOCATED(phip)) DEALLOCATE(phip)
      IF (ALLOCATED(gpsi)) DEALLOCATE(gpsi)
      IF (ALLOCATED(ipsi)) DEALLOCATE(ipsi)
      IF (ALLOCATED(pmns)) DEALLOCATE(pmns)
      IF (ALLOCATED(beta_vol)) DEALLOCATE(beta_vol)
      IF (ALLOCATED(pres)) DEALLOCATE(pres)
      IF (ALLOCATED(phi)) DEALLOCATE(phi)
      IF (ALLOCATED(buco)) DEALLOCATE(buco)
      IF (ALLOCATED(bvco)) DEALLOCATE(bvco)
      IF (ALLOCATED(rmncb)) DEALLOCATE(rmncb)
      IF (ALLOCATED(zmnsb)) DEALLOCATE(zmnsb)
      IF (ALLOCATED(pmnsb)) DEALLOCATE(pmnsb)
      IF (ALLOCATED(gmncb)) DEALLOCATE(gmncb)
      IF (ALLOCATED(bmncb)) DEALLOCATE(bmncb)
      IF (ALLOCATED(bmod_b)) DEALLOCATE(bmod_b)
!      IF (ALLOCATED(bsubumnc)) DEALLOCATE(
!     1    bsubumnc, bsubvmnc, bmodmnc, rmnc, zmns, lmns,
!     2    xm, xn, xm_nyq, xn_nyq, hiota, phip, gpsi, 
!     3    ipsi, pmns, beta_vol, pres, phi, buco, bvco,
!     3    rmncb, zmnsb, pmnsb, gmncb, bmncb, bmod_b, stat=istat1 )

      IF (ALLOCATED(cosm_b)) DEALLOCATE(cosm_b)
      IF (ALLOCATED(sinm_b)) DEALLOCATE(sinm_b)
      IF (ALLOCATED(cosn_b)) DEALLOCATE(cosn_b)
      IF (ALLOCATED(sinn_b)) DEALLOCATE(sinn_b)
      IF (ALLOCATED(cosm_nyq)) DEALLOCATE(cosm_nyq)
      IF (ALLOCATED(sinm_nyq)) DEALLOCATE(sinm_nyq)
      IF (ALLOCATED(cosn_nyq)) DEALLOCATE(cosn_nyq)
      IF (ALLOCATED(sinn_nyq)) DEALLOCATE(sinn_nyq)
      IF (ALLOCATED(sfull)) DEALLOCATE(sfull)
      IF (ALLOCATED(scl)) DEALLOCATE(scl)
      IF (ALLOCATED(xmb)) DEALLOCATE(xmb)
      IF (ALLOCATED(xnb)) DEALLOCATE(xnb)
      IF (ALLOCATED(thgrd)) DEALLOCATE(thgrd)
      IF (ALLOCATED(ztgrd)) DEALLOCATE(ztgrd)
     
!      IF (ALLOCATED(cosm_b)) DEALLOCATE(
!     1    cosm_b, sinm_b, cosn_b, sinn_b,
!     2    cosm_nyq, sinm_nyq, cosn_nyq, sinn_nyq,
!     2    sfull, scl, xmb, xnb, thgrd, ztgrd, stat=istat2)

      IF (ALLOCATED(bsubumns)) DEALLOCATE(bsubumns)
      IF (ALLOCATED(bsubvmns)) DEALLOCATE(bsubvmns)
      IF (ALLOCATED(bmodmns)) DEALLOCATE(bmodmns)
      IF (ALLOCATED(rmns)) DEALLOCATE(rmns)
      IF (ALLOCATED(zmnc)) DEALLOCATE(zmnc)
      IF (ALLOCATED(lmnc)) DEALLOCATE(lmnc)
      IF (ALLOCATED(pmnc)) DEALLOCATE(pmnc)
      IF (ALLOCATED(rmnsb)) DEALLOCATE(rmnsb)
      IF (ALLOCATED(zmncb)) DEALLOCATE(zmncb)
      IF (ALLOCATED(pmncb)) DEALLOCATE(pmncb)
      IF (ALLOCATED(gmnsb)) DEALLOCATE(gmnsb)
      IF (ALLOCATED(bmnsb)) DEALLOCATE(bmnsb)
!      IF (ALLOCATED(bsubumns)) DEALLOCATE(
!     1    bsubumns, bsubvmns, bmodmns, rmns, zmnc, lmnc, pmnc,
!     2    rmnsb, zmncb, pmncb, gmnsb, bmnsb, stat=istat3)

      IF (istat1 .ne.0 .or. istat2 .ne. 0)
     1  PRINT *,' Deallocation error in Free_mem_boozer'

      END SUBROUTINE free_mem_boozer
