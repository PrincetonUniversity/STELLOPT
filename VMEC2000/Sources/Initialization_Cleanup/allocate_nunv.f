      SUBROUTINE allocate_nunv
      USE vmec_main
      USE vmec_params, ONLY: ntmax
      USE vacmod
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: istat1
C-----------------------------------------------
      CALL free_mem_nunv

      ALLOCATE (bsubu0(nznt), rbsq(nznt), dbsq(nznt), stat=istat1)
      IF (istat1.ne.0) STOP 'allocation error #1 in allocate_nunv'

#ifdef _ANIMEC
      ALLOCATE (pperp_ns(nznt), stat=istat1)
#endif
      ALLOCATE (rmn_bdy(0:ntor,0:mpol1,ntmax),
     1          zmn_bdy(0:ntor,0:mpol1,ntmax), stat=istat1)
      IF (istat1.ne.0) STOP 'allocation error #2 in allocate_nunv'

!     PERSISTENT ARRAYS (DURATION OF PROGRAM)
      IF (lfreeb)
     1   ALLOCATE (amatsav(mnpd2*mnpd2),bvecsav(mnpd2),
     2          bsqsav(nznt,3), potvac(2*mnpd), raxis_nestor(nv),
     3          zaxis_nestor(nv), stat=istat1)
         IF (istat1.ne.0) STOP 'allocation error #3 in allocate_nunv'

      END SUBROUTINE allocate_nunv
