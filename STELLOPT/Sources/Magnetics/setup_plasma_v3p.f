!----------------------------------------------------------------------
!----------------------------------------------------------------------
      SUBROUTINE setup_plasma_v3p(ib)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  COMPUTE MAGNETIC SIGNALS                      **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          setup_plasma_v3p sets grid and computes current densty      **
!**                       components                                 **
!**                                                                  **
      USE stel_kinds
      USE stel_constants
      USE v3post_rfun
      USE read_wout_mod
      IMPLICIT NONE
      INTEGER, INTENT(in) :: ib
      INTEGER(iprec) :: i, j, k, js
      REAL(rprec)    :: fperiod
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: cosz, sinz
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: cosu, sinu,
     1                                            cosv, sinv
!
      IF (eqtype.eq.'none') RETURN
!----------------------------------------------------------------------
!-- get plasma contribution                                          --
!----------------------------------------------------------------------
        ju = ns
!  JDH 06.28.2003. Getting rid of boxdims
!    kp comes through v3post_rfun, and is defined in get_plasma_v3p.
!        kv = boxdims(3)
        kv = kp
        fperiod = twopi / n_field_periods
        dels = one / (ns - 1)
        delu = twopi / ju
        delv = fperiod / kv
!----------------------------------------------------------------------
!-- allocate arrays                                                  --
!----------------------------------------------------------------------
        IF (ALLOCATED(sgrid)) DEALLOCATE (sgrid)
        IF (ALLOCATED(ugrid)) DEALLOCATE (ugrid)
        IF (ALLOCATED(vgrid)) DEALLOCATE (vgrid)
        IF (ALLOCATED(rsuv)) DEALLOCATE (rsuv)
        IF (ALLOCATED(gsuv)) DEALLOCATE (gsuv)
        IF (ALLOCATED(zsuv)) DEALLOCATE (zsuv)
!        IF (ALLOCATED(psuv)) DEALLOCATE (psuv)
        IF (ALLOCATED(currusuv)) DEALLOCATE (currusuv)
        IF (ALLOCATED(currvsuv)) DEALLOCATE (currvsuv)
        IF (ALLOCATED(rusuv)) DEALLOCATE (rusuv)
        IF (ALLOCATED(zusuv)) DEALLOCATE (zusuv)
        IF (ALLOCATED(rvsuv)) DEALLOCATE (rvsuv)
        IF (ALLOCATED(zvsuv)) DEALLOCATE (zvsuv)
        ALLOCATE (sgrid(ns),ugrid(ju),vgrid(kv),stat=i)
        IF (i .ne. 0) GOTO 1000
        ALLOCATE (rsuv(ns,ju,kv),zsuv(ns,ju,kv),stat=i)
        IF (i .ne. 0) GOTO 1000
        ALLOCATE (gsuv(ns,ju,kv),stat=i)
        IF (i .ne. 0) GOTO 1000
!        ALLOCATE (psuv(ns,ju,kv))
        ALLOCATE (currusuv(ns,ju,kv),currvsuv(ns,ju,kv),stat=i)
        IF (i .ne. 0) GOTO 1000
        ALLOCATE (rusuv(ns,ju,kv),zusuv(ns,ju,kv),stat=i)
        IF (i .ne. 0) GOTO 1000
        ALLOCATE (rvsuv(ns,ju,kv),zvsuv(ns,ju,kv),stat=i)
        IF (i .ne. 0) GOTO 1000
        ALLOCATE (cosz(mnmax), sinz(mnmax),stat=i)
        IF (i .ne. 0) GOTO 1000
        ALLOCATE (cosu(mnmax, ju), sinu(mnmax, ju),
     &            cosv(mnmax, kv), sinv(mnmax, kv),stat=i)
        IF (i .ne. 0) GOTO 1000
!----------------------------------------------------------------------
!-- Allocate arrays for surface formulation                          --
!----------------------------------------------------------------------
        IF (ALLOCATED(bsubu)) DEALLOCATE (bsubu)
        IF (ALLOCATED(bsubv)) DEALLOCATE (bsubv)
        ALLOCATE (bsubu(ns,ju,kv),bsubv(ns,ju,kv),stat=i)
        IF (i .ne. 0) GOTO 1000
!----------------------------------------------------------------------
!-- set up (s,u,v) VMEC flux grids                                   --
!----------------------------------------------------------------------
        DO i = 1, ns
          sgrid(i) = (i - 1) * dels
        ENDDO
        DO i = 1, ju
          ugrid(i) = (i - 1) * delu
          cosu(:,i) = COS(xm*ugrid(i))
          sinu(:,i) = SIN(xm*ugrid(i))
        ENDDO
        DO i = 1, kv
          vgrid(i) = (i - 1) * delv
          cosv(:,i) = COS(xn*vgrid(i))
          sinv(:,i) = SIN(xn*vgrid(i))
        ENDDO

        IF (lsurf) THEN
           bsubu(1:ns-2,:,:) = 0;    bsubv(1:ns-2,:,:) = 0
           rsuv(1:ns-1,:,:) = 0;     zsuv(1:ns-1,:,:) = 0
           rusuv(1:ns-1,:,:) = 0;    zusuv(1:ns-1,:,:) = 0
           rvsuv(1:ns-1,:,:) = 0;    zvsuv(1:ns-1,:,:) = 0
           gsuv(1:ns-1,:,:) = 0
        ENDIF

!
!       DO FOURIER MODE SUM IN INNERMOST LOOP
!       FIRST DO mnmax-sized ARRAYS (r, z)
!
        DO j = 1, ju
        DO k = 1, kv
!          zetamn = xm*ugrid(j) - xn*vgrid(k)
!          cosz = COS(zetamn)
!          sinz = SIN(zetamn)
          cosz = cosu(:,j)*cosv(:,k) + sinu(:,j)*sinv(:,k)
          sinz = sinu(:,j)*cosv(:,k) - cosu(:,j)*sinv(:,k)
          DO js = ib, ns
            rsuv(js,j,k) = SUM(rmnc(:,js)*cosz)
            zsuv(js,j,k) = SUM(zmns(:,js)*sinz)
            rusuv(js,j,k) = -SUM(xm*rmnc(:,js)*sinz)
            zusuv(js,j,k) = SUM(xm*zmns(:,js)*cosz)
            rvsuv(js,j,k) = SUM(xn*rmnc(:,js)*sinz)
            zvsuv(js,j,k) = -SUM(xn*zmns(:,js)*cosz)
          ENDDO

!----------------------------------------------------------------------
!-- stellarator asymmetric terms                                     --
!----------------------------------------------------------------------
          IF (lasym) THEN
            DO js = ib, ns
              rsuv(js,j,k) = rsuv(js,j,k) + SUM(rmns(:,js)*sinz)
              zsuv(js,j,k) = zsuv(js,j,k) + SUM(zmnc(:,js)*cosz)
              rusuv(js,j,k) = rusuv(js,j,k) + SUM(xm*rmns(:,js)*cosz)
              zusuv(js,j,k) = zusuv(js,j,k) - SUM(xm*zmnc(:,js)*sinz)
              rvsuv(js,j,k) = rvsuv(js,j,k) - SUM(xn*rmns(:,js)*cosz)
              zvsuv(js,j,k) = zvsuv(js,j,k) + SUM(xn*zmnc(:,js)*sinz)
            ENDDO
          ENDIF
        ENDDO
        ENDDO

        DEALLOCATE (cosz, sinz, cosu, cosv, sinu, sinv)

!
!       NEXT DO mnmax_nyq-sized ARRAYS (g, curru,v, bsubu,v)
!
        ALLOCATE (cosz(mnmax_nyq), sinz(mnmax_nyq),stat=i)
        IF (i .ne. 0) GOTO 1000
        ALLOCATE (cosu(mnmax_nyq, ju), sinu(mnmax_nyq, ju),
     &            cosv(mnmax_nyq, kv), sinv(mnmax_nyq, kv),stat=i)
        IF (i .ne. 0) GOTO 1000

        DO i = 1, ju
          cosu(:,i) = COS(xm_nyq(:)*ugrid(i))
          sinu(:,i) = SIN(xm_nyq(:)*ugrid(i))
        ENDDO
        DO i = 1, kv
          vgrid(i) = (i - 1) * delv
          cosv(:,i) = COS(xn_nyq(:)*vgrid(i))
          sinv(:,i) = SIN(xn_nyq(:)*vgrid(i))
        ENDDO

        DO j = 1, ju
        DO k = 1, kv
          cosz = cosu(:,j)*cosv(:,k) + sinu(:,j)*sinv(:,k)
          DO js = ib, ns
            gsuv(js,j,k) = SUM(gmnc(:,js)*cosz)
            currusuv(js,j,k) = SUM(currumnc(:,js)*cosz)
            currvsuv(js,j,k) = SUM(currvmnc(:,js)*cosz)
          ENDDO

          IF (lsurf) THEN
            DO js = ns-1, ns
               bsubu(js,j,k) = SUM(bsubumnc(:,js)*cosz)
               bsubv(js,j,k) = SUM(bsubvmnc(:,js)*cosz)
            ENDDO
          ENDIF

!----------------------------------------------------------------------
!-- stellarator asymmetric terms                                     --
!----------------------------------------------------------------------
          IF (lasym) THEN
            sinz = sinu(:,j)*cosv(:,k) - cosu(:,j)*sinv(:,k)
            DO js = ib, ns
              gsuv(js,j,k) = gsuv(js,j,k) + SUM(gmns(:,js)*sinz)
              currusuv(js,j,k) = currusuv(js,j,k)
     &                                    + SUM(currumns(:,js)*sinz)
              currvsuv(js,j,k) = currvsuv(js,j,k)
     &                                    + SUM(currvmns(:,js)*sinz)
            ENDDO
            IF (lsurf) THEN
              DO js = ns-1, ns
                 bsubu(js,j,k) = bsubu(js,j,k)+SUM(bsubumns(:,js)*sinz)
                 bsubv(js,j,k) = bsubv(js,j,k)+SUM(bsubvmns(:,js)*sinz)
              ENDDO
            ENDIF
          ENDIF
        ENDDO
        ENDDO

        DEALLOCATE (cosz, sinz, cosu, cosv, sinu, sinv)

!       CALL flush(6)

        RETURN

 1000 CONTINUE
        STOP 'Memory allocation error in setup_plasma_v3p!'
      END SUBROUTINE setup_plasma_v3p
