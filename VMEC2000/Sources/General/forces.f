      SUBROUTINE forces_par
      USE vmec_main, p5 => cp5
      USE realspace
      USE vforces, ru12 => pazmn_e, zu12 => parmn_e, 
     1             pazmn_e => pazmn_e, parmn_e => parmn_e,
     2             lv_e => pcrmn_e, lu_e => pczmn_e, lu_o => pczmn_o,
     3             pcrmn_e => pcrmn_e, pczmn_e => pczmn_e, 
     4             pczmn_o => pczmn_o
      USE parallel_include_module
      USE timer_sub
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(dp), PARAMETER :: p25 = p5*p5, dshalfds=p25
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: l, ndim
      INTEGER :: i, j, k, nsmin, nsmax
      REAL(dp), DIMENSION(:,:), POINTER :: 
     1    bsqr, gvvs, guvs, guus
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: bcastbuf
C-----------------------------------------------
      IF (.NOT.lactive .AND. .NOT.lfreeb) RETURN

      CALL second0 (tforon)

      ndim = 1+nrzt

      nsmin=tlglob; nsmax=t1rglob

!     POINTER ALIASES
      bsqr => pextra1(:,:,1);  gvvs => pextra2(:,:,1)
      guvs => pextra3(:,:,1);  guus => pextra4(:,:,1)

      lu_e(:,1) = 0; lv_e(:,1) = 0
      pguu(:,1)  = 0; pguv(:,1)  = 0; pgvv(:,1) = 0

      DO l = nsmin, nsmax      
      guus(:,l) = pguu(:,l)*pshalf(:,l)    
      guvs(:,l) = pguv(:,l)*pshalf(:,l)    
      gvvs(:,l) = pgvv(:,l)* pshalf(:,l)    

      parmn_e(:,l) = ohs*zu12(:,l)*lu_e(:,l)
      pazmn_e(:,l) =-ohs*ru12(:,l)*lu_e(:,l)
      pbrmn_e(:,l) = pbrmn_e(:,l)*lu_e(:,l)
      pbzmn_e(:,l) =-pbzmn_e(:,l)*lu_e(:,l)
      bsqr(:,l)    = dshalfds*lu_e(:,l)/pshalf(:,l)

      parmn_o(:,l) = parmn_e(:,l)*pshalf(:,l)
      pazmn_o(:,l) = pazmn_e(:,l)*pshalf(:,l)
      pbrmn_o(:,l) = pbrmn_e(:,l)*pshalf(:,l)
      pbzmn_o(:,l) = pbzmn_e(:,l)*pshalf(:,l)
      END DO

!
!     CONSTRUCT CYLINDRICAL FORCE KERNELS
!     NOTE: presg(ns+1) == 0, AND WILL BE "FILLED IN" AT EDGE
!     FOR FREE-BOUNDARY BY RBSQ
!
!DIR$ IVDEP
      nsmin=tlglob; nsmax=MIN(ns-1,trglob)
      DO l = nsmin, nsmax
      pguu(:,l) = p5*(pguu(:,l) + pguu(:,l+1))
      pgvv(:,l) = p5*(pgvv(:,l) + pgvv(:,l+1))
      bsqr(:,l) = bsqr(:,l) +     bsqr(:,l+1)
      guus(:,l) = p5*(guus(:,l) + guus(:,l+1))
      gvvs(:,l) = p5*(gvvs(:,l) + gvvs(:,l+1))
      END DO

      IF (trglob .ge. ns) THEN
      pguu(:,ns) = p5*pguu(:,ns)
      pgvv(:,ns) = p5*pgvv(:,ns)
      guus(:,ns) = p5*guus(:,ns)
      gvvs(:,ns) = p5*gvvs(:,ns)
      END IF

!DIR$ IVDEP
      nsmin=tlglob; nsmax=MIN(ns-1,trglob)
      DO l = nsmin, nsmax
      parmn_e(:,l) = parmn_e(:,l+1) - parmn_e(:,l) 
     1             + p5*(lv_e(:,l) + lv_e(:,l+1))
      pazmn_e(:,l) = pazmn_e(:,l+1) - pazmn_e(:,l)
      pbrmn_e(:,l) = p5*(pbrmn_e(:,l) + pbrmn_e(:,l+1))
      pbzmn_e(:,l) = p5*(pbzmn_e(:,l) + pbzmn_e(:,l+1))
      END DO

      parmn_e(:,ns) = - parmn_e(:,ns) + p5*lv_e(:,ns)  
      pazmn_e(:,ns) = - pazmn_e(:,ns)
      pbrmn_e(:,ns) = p5*pbrmn_e(:,ns)
      pbzmn_e(:,ns) = p5*pbzmn_e(:,ns) 

      nsmin=tlglob; nsmax=t1rglob
      DO l = nsmin, nsmax
      parmn_e(:,l) = parmn_e(:,l) 
     1             - (gvvs(:,l)*pr1(:,l,1) + pgvv(:,l)*pr1(:,l,0))
      pbrmn_e(:,l) = pbrmn_e(:,l) + bsqr(:,l)*pz1(:,l,1)
     1              - (guus(:,l)*pru(:,l,1) + pguu(:,l)*pru(:,l,0))
      pbzmn_e(:,l) = pbzmn_e(:,l) - (bsqr(:,l)*pr1(:,l,1)
     2              +  guus(:,l)*pzu(:,l,1) + pguu(:,l)*pzu(:,l,0))
      lv_e(:,l) = lv_e(:,l)*pshalf(:,l)
      lu_o(:,l) = dshalfds*lu_e(:,l)
      END DO

      nsmin=tlglob; nsmax=MIN(ns-1,trglob)
!DIR$ IVDEP
      DO l = nsmin, nsmax
      parmn_o(:,l) = parmn_o(:,l+1) - parmn_o(:,l) 
     1             - pzu(:,l,0)*bsqr(:,l) + p5*(lv_e(:,l)+lv_e(:,l+1))
      pazmn_o(:,l) = pazmn_o(:,l+1) - pazmn_o(:,l) 
     1             + pru(:,l,0)*bsqr(:,l)
      pbrmn_o(:,l) = p5*(pbrmn_o(:,l) +  pbrmn_o(:,l+1))
      pbzmn_o(:,l) = p5*(pbzmn_o(:,l) +  pbzmn_o(:,l+1))
      lu_o(:,l)   = lu_o(:,l) + lu_o(:,l+1)
      END DO

      parmn_o(:,ns) = - parmn_o(:,ns) - pzu(:,ns,0)*bsqr(:,ns)
     1                + p5*lv_e(:,ns)
      pazmn_o(:,ns) = - pazmn_o(:,ns) + pru(:,ns,0)*bsqr(:,ns)
      pbrmn_o(:,ns) = p5*pbrmn_o(:,ns)
      pbzmn_o(:,ns) = p5*pbzmn_o(:,ns) 
      lu_o(:,ns)   = lu_o(:,ns)

      nsmin=tlglob; nsmax=trglob
      DO l = nsmin, nsmax
      pguu(:,l) = pguu(:,l) * psqrts(:,l)**2
      bsqr(:,l) = pgvv(:,l) * psqrts(:,l)**2
      END DO

      DO l = nsmin, nsmax
      parmn_o(:,l) = parmn_o(:,l) - (pzu(:,l,1)*lu_o(:,l)
     1             + bsqr(:,l)*pr1(:,l,1) + gvvs(:,l)*pr1(:,l,0))
      pazmn_o(:,l) = pazmn_o(:,l) +  pru(:,l,1)*lu_o(:,l)
      pbrmn_o(:,l) = pbrmn_o(:,l) +  pz1(:,l,1)*lu_o(:,l)
     1             -(pguu(:,l)*pru(:,l,1) + guus(:,l)*pru(:,l,0))
      pbzmn_o(:,l) = pbzmn_o(:,l) - (pr1(:,l,1)*lu_o(:,l)
     1             + pguu(:,l)*pzu(:,l,1) + guus(:,l)*pzu(:,l,0))
      END DO

      IF (lthreed) THEN
!DIR$ IVDEP
        nsmin=tlglob; nsmax=MIN(ns-1,trglob)
        DO l = nsmin, nsmax
          pguv(:,l)  = p5*(pguv(:,l) + pguv(:,l+1))
          guvs(:,l) = p5*(guvs(:,l) + guvs(:,l+1))
        END DO
        pguv(:,ns) = p5*pguv(:,ns)
        guvs(:,ns) = p5*guvs(:,ns)

        nsmin=tlglob; nsmax=trglob
        DO l = nsmin, nsmax
          pbrmn_e(:,l) = pbrmn_e(:,l) 
     1                 - (pguv(:,l)*prv(:,l,0) + guvs(:,l)*prv(:,l,1))
          pbzmn_e(:,l) = pbzmn_e(:,l) 
     1          - (pguv(:,l)*pzv(:,l,0) + guvs(:,l)*pzv(:,l,1))
          pcrmn_e(:,l) = pguv(:,l)*pru(:,l,0) + pgvv(:,l)*prv(:,l,0)
     1                 + gvvs(:,l)*prv(:,l,1) + guvs(:,l)*pru(:,l,1)
          pczmn_e(:,l) = pguv(:,l)*pzu(:,l,0) + pgvv(:,l)*pzv(:,l,0)
     1                 + gvvs(:,l)*pzv(:,l,1) + guvs(:,l)*pzu(:,l,1)
          pguv(:,l) = pguv(:,l)*psqrts(:,l)*psqrts(:,l)
          pbrmn_o(:,l) = pbrmn_o(:,l) 
     1          - (guvs(:,l)*prv(:,l,0) + pguv(:,l)*prv(:,l,1))
          pbzmn_o(:,l) = pbzmn_o(:,l) 
     1          - (guvs(:,l)*pzv(:,l,0) + pguv(:,l)*pzv(:,l,1))
          pcrmn_o(:,l) = guvs(:,l)*pru(:,l,0) + gvvs(:,l)*prv(:,l,0)
     1          + bsqr(:,l)*prv(:,l,1) + pguv(:,l)*pru(:,l,1)
          pczmn_o(:,l) = guvs(:,l)*pzu(:,l,0) + gvvs(:,l)*pzv(:,l,0)
     1          + bsqr(:,l)*pzv(:,l,1) + pguv(:,l)*pzu(:,l,1)
        END DO
      ENDIF

!
!     ASSIGN EDGE FORCES (JS = NS) FOR FREE BOUNDARY CALCULATION
!
      IF (ivac .GE. 1) THEN

        DO k = 1, ntheta3
          DO j = 1, nzeta
            l = (k-1)*nzeta + j
            parmn_e(l,ns) = parmn_e(l,ns) + pzu0(l,ns)*rbsq(l)
            parmn_o(l,ns) = parmn_o(l,ns) + pzu0(l,ns)*rbsq(l)
            pazmn_e(l,ns) = pazmn_e(l,ns) - pru0(l,ns)*rbsq(l)
            pazmn_o(l,ns) = pazmn_o(l,ns) - pru0(l,ns)*rbsq(l)
          END DO
        END DO

      ENDIF

 100  CONTINUE

!
!     COMPUTE CONSTRAINT FORCE KERNELS
!
#ifndef _HBANGLE
      DO l = nsmin, nsmax
        prcon(:,l,0) = (prcon(:,l,0)-prcon0(:,l))*pgcon(:,l)
        pzcon(:,l,0) = (pzcon(:,l,0)-pzcon0(:,l))*pgcon(:,l)
        pbrmn_e(:,l) = pbrmn_e(:,l) + prcon(:,l,0)
        pbzmn_e(:,l) = pbzmn_e(:,l) + pzcon(:,l,0)
        pbrmn_o(:,l) = pbrmn_o(:,l)+ prcon(:,l,0)*psqrts(:,l)
        pbzmn_o(:,l) = pbzmn_o(:,l)+ pzcon(:,l,0)*psqrts(:,l)
        prcon(:,l,0) = pru0(:,l) * pgcon(:,l)
        pzcon(:,l,0) = pzu0(:,l) * pgcon(:,l)
        prcon(:,l,1) = prcon(:,l,0) * psqrts(:,l)
        pzcon(:,l,1) = pzcon(:,l,0) * psqrts(:,l)
      END DO
#endif
      CALL second0 (tforoff)
      timer(tfor) = timer(tfor) + (tforoff - tforon)
      forces_time = timer(tfor)

      END SUBROUTINE forces_par

      SUBROUTINE forces
      USE vmec_main, p5 => cp5
      USE realspace
      USE vforces, ru12 => azmn_e, zu12 => armn_e, 
     1             azmn_e => azmn_e, armn_e => armn_e,
     2             lv_e => crmn_e, lu_e => czmn_e, lu_o => czmn_o,
     3             crmn_e => crmn_e, czmn_e => czmn_e, czmn_o => czmn_o
      USE timer_sub

      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(dp), PARAMETER :: p25 = p5*p5, dshalfds=p25
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: l, ndim
      REAL(dp), DIMENSION(:), POINTER :: 
     1    bsqr, gvvs, guvs, guus

C-----------------------------------------------
      ndim = 1+nrzt
      CALL second0 (tforon)

!     POINTER ALIASES
      bsqr => extra1(:,1);  gvvs => extra2(:,1)
      guvs => extra3(:,1);  guus => extra4(:,1)

!
!     ON ENTRY, ARMN=ZU,BRMN=ZS,AZMN=RU,BZMN=RS,LU=R*BSQ,LV = BSQ*SQRT(G)/R12
!     HERE, XS (X=Z,R) DO NOT INCLUDE DERIVATIVE OF EXPLICIT SQRT(S)
!     BSQ = |B|**2/2 + p
!     GIJ = (BsupI * BsupJ) * SQRT(G)  (I,J = U,V)
!     IT IS ESSENTIAL THAT LU,LV AT j=1 ARE ZERO INITIALLY
!
!     SOME OF THE BIGGER LOOPS WERE SPLIT TO FACILITATE CACHE
!     HITS, PIPELINING ON RISCS
!
!     FOR OPTIMIZATION ON CRAY, MUST USE COMPILER DIRECTIVES TO
!     GET VECTORIZATION OF LOOPS INVOLVING POINTERS!
!
!
!     ORIGIN OF VARIOUS TERMS
!
!     LU :  VARIATION OF DOMINANT RADIAL-DERIVATIVE TERMS IN JACOBIAN
!
!     LV :  VARIATION OF R-TERM IN JACOBIAN
!
!     GVV:  VARIATION OF R**2-TERM AND Rv**2,Zv**2 IN gvv
!
!     GUU, GUV: VARIATION OF Ru, Rv, Zu, Zv IN guu, guv
!

      lu_e(1:ndim:ns) = 0; lv_e(1:ndim:ns) = 0
      guu(1:ndim:ns)  = 0; guv(1:ndim:ns)  = 0; gvv(1:ndim:ns) = 0
      guus = guu*shalf;    guvs = guv*shalf;    gvvs = gvv*shalf

      armn_e  = ohs*zu12 * lu_e
      azmn_e  =-ohs*ru12 * lu_e
      brmn_e  = brmn_e * lu_e
      bzmn_e  =-bzmn_e * lu_e
      bsqr    = dshalfds*lu_e/shalf

      armn_o(1:ndim)  = armn_e(1:ndim) *shalf
      azmn_o(1:ndim)  = azmn_e(1:ndim) *shalf
      brmn_o(1:ndim)  = brmn_e(1:ndim) *shalf
      bzmn_o(1:ndim)  = bzmn_e(1:ndim) *shalf
      
!
!     CONSTRUCT CYLINDRICAL FORCE KERNELS
!     NOTE: presg(ns+1) == 0, AND WILL BE "FILLED IN" AT EDGE
!     FOR FREE-BOUNDARY BY RBSQ
!
!DIR$ IVDEP
      DO l = 1, nrzt
         guu(l) = p5*(guu(l) + guu(l+1))
         gvv(l) = p5*(gvv(l) + gvv(l+1))
         bsqr(l) = bsqr(l) + bsqr(l+1)
         guus(l) = p5*(guus(l) + guus(l+1))
         gvvs(l) = p5*(gvvs(l) + gvvs(l+1))
      END DO

!DIR$ IVDEP
      DO l = 1, nrzt
         armn_e(l) = armn_e(l+1) - armn_e(l) + p5*(lv_e(l) + lv_e(l+1))
         azmn_e(l) = azmn_e(l+1) - azmn_e(l)
         brmn_e(l) = p5*(brmn_e(l) + brmn_e(l+1))
         bzmn_e(l) = p5*(bzmn_e(l) + bzmn_e(l+1))
      END DO

      armn_e(:nrzt) = armn_e(:nrzt) - (gvvs(:nrzt)*r1(:nrzt,1)
     1              + gvv(:nrzt)*r1(:nrzt,0))
      brmn_e(:nrzt) = brmn_e(:nrzt) + bsqr(:nrzt)*z1(:nrzt,1)
     1              -(guus(:nrzt)*ru(:nrzt,1) + guu(:nrzt)*ru(:nrzt,0))
      bzmn_e(:nrzt) = bzmn_e(:nrzt) - (bsqr(:nrzt)*r1(:nrzt,1)
     1              + guus(:nrzt)*zu(:nrzt,1) + guu(:nrzt)*zu(:nrzt,0))
      lv_e(1:ndim) = lv_e(1:ndim)*shalf(1:ndim)
      lu_o(1:ndim) = dshalfds*lu_e(1:ndim)

!DIR$ IVDEP
      DO l = 1, nrzt
         armn_o(l) = armn_o(l+1) - armn_o(l) - zu(l,0)*bsqr(l)
     1             + p5*(lv_e(l) + lv_e(l+1))
         azmn_o(l) = azmn_o(l+1) - azmn_o(l) + ru(l,0)*bsqr(l)
         brmn_o(l) = p5*(brmn_o(l) + brmn_o(l+1))
         bzmn_o(l) = p5*(bzmn_o(l) + bzmn_o(l+1))
         lu_o(l)   = lu_o(l) + lu_o(l+1)
      END DO

      guu(1:nrzt)  = guu(1:nrzt) * sqrts(1:nrzt)**2
      bsqr(1:nrzt) = gvv(1:nrzt) * sqrts(1:nrzt)**2

      armn_o(:nrzt) = armn_o(:nrzt) - (zu(:nrzt,1)*lu_o(:nrzt)
     1             + bsqr(:nrzt)*r1(:nrzt,1) + gvvs(:nrzt)*r1(:nrzt,0))
      azmn_o(:nrzt) = azmn_o(:nrzt) + ru(:nrzt,1)*lu_o(:nrzt)
      brmn_o(:nrzt) = brmn_o(:nrzt) + z1(:nrzt,1)*lu_o(:nrzt)
     1             -(guu(:nrzt)*ru(:nrzt,1) + guus(:nrzt)*ru(:nrzt,0))
      bzmn_o(:nrzt) = bzmn_o(:nrzt) - (r1(:nrzt,1)*lu_o(:nrzt)
     1             + guu(:nrzt)*zu(:nrzt,1) + guus(:nrzt)*zu(:nrzt,0))

      IF (lthreed) THEN
!DIR$ IVDEP
         DO l = 1, nrzt
            guv(l)  = p5*(guv(l) + guv(l+1))
            guvs(l) = p5*(guvs(l) + guvs(l+1))
         END DO

         brmn_e(:nrzt) = brmn_e(:nrzt) 
     1          - (guv(:nrzt)*rv(:nrzt,0) + guvs(:nrzt)*rv(:nrzt,1))
         bzmn_e(:nrzt) = bzmn_e(:nrzt) 
     1          - (guv(:nrzt)*zv(:nrzt,0) + guvs(:nrzt)*zv(:nrzt,1))
         crmn_e(:nrzt) = guv(:nrzt) *ru(:nrzt,0) 
     1                 + gvv(:nrzt) *rv(:nrzt,0)
     2          + gvvs(:nrzt)*rv(:nrzt,1) + guvs(:nrzt)*ru(:nrzt,1)
         czmn_e(:nrzt) = guv(:nrzt) *zu(:nrzt,0) 
     1                 + gvv(:nrzt) *zv(:nrzt,0)
     2          + gvvs(:nrzt)*zv(:nrzt,1) + guvs(:nrzt)*zu(:nrzt,1)
         guv(:nrzt) = guv(:nrzt) *sqrts(:nrzt)*sqrts(:nrzt)
         brmn_o(:nrzt) = brmn_o(:nrzt) 
     1          - (guvs(:nrzt)*rv(:nrzt,0) + guv(:nrzt)*rv(:nrzt,1))
         bzmn_o(:nrzt) = bzmn_o(:nrzt) 
     1          - (guvs(:nrzt)*zv(:nrzt,0) + guv(:nrzt)*zv(:nrzt,1))
         crmn_o(:nrzt) = guvs(:nrzt)*ru(:nrzt,0) 
     1                 + gvvs(:nrzt)*rv(:nrzt,0)
     2          + bsqr(:nrzt)*rv(:nrzt,1) + guv(:nrzt) *ru(:nrzt,1)
         czmn_o(:nrzt) = guvs(:nrzt)*zu(:nrzt,0) 
     1                 + gvvs(:nrzt)*zv(:nrzt,0)
     2          + bsqr(:nrzt)*zv(:nrzt,1) + guv(:nrzt) *zu(:nrzt,1)
      ENDIF

!
!     ASSIGN EDGE FORCES (JS = NS) FOR FREE BOUNDARY CALCULATION
!
      IF (ivac .ge. 1) THEN
         armn_e(ns:nrzt:ns) = armn_e(ns:nrzt:ns) 
     1                      + zu0(ns:nrzt:ns)*rbsq(1:nznt)
         armn_o(ns:nrzt:ns) = armn_o(ns:nrzt:ns) 
     1                      + zu0(ns:nrzt:ns)*rbsq(1:nznt)
         azmn_e(ns:nrzt:ns) = azmn_e(ns:nrzt:ns) 
     1                      - ru0(ns:nrzt:ns)*rbsq(1:nznt)
         azmn_o(ns:nrzt:ns) = azmn_o(ns:nrzt:ns) 
     1                      - ru0(ns:nrzt:ns)*rbsq(1:nznt)
!         fz00_edge = SUM(wint(ns:nrzt:ns)*ru0(ns:nrzt:ns)*rbsq(1:nznt))
      ENDIF

 100  CONTINUE
!      CALL Print1DArrayMNSP (gcon,tlglob,trglob,500,.TRUE.,"prcon")

!
!     COMPUTE CONSTRAINT FORCE KERNELS
!
#ifndef _HBANGLE
      rcon(:nrzt,0) = (rcon(:nrzt,0) - rcon0(:nrzt)) * gcon(:nrzt)
      zcon(:nrzt,0) = (zcon(:nrzt,0) - zcon0(:nrzt)) * gcon(:nrzt)
      brmn_e(:nrzt) = brmn_e(:nrzt) + rcon(:nrzt,0)
      bzmn_e(:nrzt) = bzmn_e(:nrzt) + zcon(:nrzt,0)
      brmn_o(:nrzt) = brmn_o(:nrzt)+ rcon(:nrzt,0)*sqrts(:nrzt)
      bzmn_o(:nrzt) = bzmn_o(:nrzt)+ zcon(:nrzt,0)*sqrts(:nrzt)
      rcon(:nrzt,0) =  ru0(:nrzt) * gcon(:nrzt)
      zcon(:nrzt,0) =  zu0(:nrzt) * gcon(:nrzt)
      rcon(:nrzt,1) = rcon(:nrzt,0) * sqrts(:nrzt)
      zcon(:nrzt,1) = zcon(:nrzt,0) * sqrts(:nrzt)
#endif
      
      CALL second0 (tforoff)
      timer(tfor) = timer(tfor) + (tforoff - tforon)

      END SUBROUTINE forces
