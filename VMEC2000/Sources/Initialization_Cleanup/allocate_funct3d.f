      SUBROUTINE allocate_funct3d
      USE vmec_main
      USE realspace
      USE vforces
#ifdef _VACUUM2
      USE vac2_vacmod
#else
      USE vacmod
#endif
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: istat1, ndim, ndim2
C-----------------------------------------------
      ndim  = 1+nrzt
      ndim2 = 2*ndim

      CALL free_mem_funct3d

      ALLOCATE (armn(ndim2), azmn(ndim2), brmn(ndim2), bzmn(ndim2),
     1   crmn(ndim2), czmn(ndim2), blmn(ndim2), clmn(ndim2),
     2   r1(nrzt,0:1), ru(nrzt,0:1), rv(nrzt,0:1),
     3   z1(nrzt,0:1), zu(nrzt,0:1), zv(nrzt,0:1),
     4   rcon(nrzt,0:1), zcon(nrzt,0:1), ru0(ndim), zu0(ndim),
     5   rcon0(ndim), zcon0(ndim), guu(ndim), guv(ndim), gvv(ndim), 
     6   gcon(ndim), sigma_an(nrzt), stat=istat1)
      IF (istat1.ne.0) STOP 'allocation error #1 in allocate_funct3d'
      armn=0; azmn=0; brmn=0; bzmn=0; crmn=0; czmn=0; blmn=0; clmn=0
      r1=0; ru=0; rv=0; z1=0; zu=0; zv=0; rcon=0; zcon=0
      ru0=0; zu0=0; guu=0; guv=0; gvv=0
      sigma_an=1

#ifdef _ANIMEC
      ALLOCATE(pperp(nrzt), ppar(nrzt), onembc(nrzt),
     1   pp1(nrzt), pp2(nrzt), pp3(nrzt), stat=istat1)
      IF (istat1.ne.0) STOP 'allocation error #1A in allocate_funct3d'
      pperp=0; ppar=0; onembc=0; pp1=0; pp2=0; pp3=0
#endif

      IF (lfreeb) THEN
         ALLOCATE (brv(nznt), bphiv(nznt), bzv(nznt), bsqvac(nznt),
     1             stat=istat1)
         IF (istat1.ne.0) STOP 'allocation error #2 in allocate_funct3d'
         brv=0; bphiv=0; bzv=0; bsqvac=0
      END IF

      ALLOCATE (extra1(ndim,0:1), stat=istat1)
      IF (istat1.ne.0) STOP 'allocation error #3 in allocate_funct3d'
      extra1=0
      
      IF (lasym) THEN
         ALLOCATE (extra2(ndim,0:1), extra3(ndim,0:1), 
     1             extra4(ndim,0:1),stat=istat1)
      ELSE
         ALLOCATE (extra2(ndim,1), extra3(ndim,1), extra4(ndim,1),
     1             stat=istat1)
      END IF
      IF (istat1.ne.0) STOP 'allocation error #3 in allocate_funct3d'
      extra2=0; extra3=0; extra4=0
!
!     Pointer alias assignments
!     NOTE: In FORCES, X_e(nrzt+1) overlaps X_o(1), which should never be used...
!
      armn_e => armn(:ndim)
      armn_o => armn(ndim:)
      armn(:ndim2) = zero
      brmn_e => brmn(:ndim)
      brmn_o => brmn(ndim:)
      brmn(:ndim2) = zero
      azmn_e => azmn(:ndim)
      azmn_o => azmn(ndim:)
      azmn(:ndim2) = zero
      bzmn_e => bzmn(:ndim)
      bzmn_o => bzmn(ndim:)
      bzmn(:ndim2) = zero
      crmn_e => crmn(:ndim)
      crmn_o => crmn(ndim:)
      crmn(:ndim2) = zero
      czmn_e => czmn(:ndim)
      czmn_o => czmn(ndim:)
      czmn(:ndim2) = zero
      blmn_e => blmn(:ndim)
      blmn_o => blmn(ndim:)
      blmn(:ndim2) = zero
      clmn_e => clmn(:ndim)
      clmn_o => clmn(ndim:)
      clmn(:ndim2) = zero
      rcon0(:ndim) = zero
      zcon0(:ndim) = zero

      END SUBROUTINE allocate_funct3d
