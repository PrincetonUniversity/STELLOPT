      SUBROUTINE forces
      USE vmec_main, p5 => cp5
      USE realspace
      USE vforces, ru12 => azmn_e, zu12 => armn_e, 
     1             azmn_e => azmn_e, armn_e => armn_e,
     2             lv_e => crmn_e, lu_e => czmn_e, lu_o => czmn_o,
     3             crmn_e => crmn_e, czmn_e => czmn_e, czmn_o => czmn_o
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: p25 = p5*p5, dshalfds=p25
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: l, ndim
      REAL(rprec), DIMENSION(:), POINTER :: 
     1    bsqr, gvvs, guvs, guus
C-----------------------------------------------
      ndim = 1+nrzt
!
!     POINTER ALIASES
!
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
      DO l = 1, nrzt
         guu(l) = p5*(guu(l) + guu(l+1))
         gvv(l) = p5*(gvv(l) + gvv(l+1))
         bsqr(l) = bsqr(l) + bsqr(l+1)
         guus(l) = p5*(guus(l) + guus(l+1))
         gvvs(l) = p5*(gvvs(l) + gvvs(l+1))
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
      END SUBROUTINE forces
