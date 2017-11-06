      SUBROUTINE findphi(reven, rodd, rmeas, dse, dso, rmid, ismeas,
     1   iumeas, indexr, npts)
      USE vmec_main
      USE realspace
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER npts
      INTEGER, DIMENSION(npts) :: ismeas, iumeas
      INTEGER, DIMENSION(2*ns) :: indexr
      REAL(rprec), DIMENSION(ns,nzeta,*) :: reven, rodd
      REAL(rprec), DIMENSION(npts) :: rmeas, dse, dso
      REAL(rprec), DIMENSION(2*ns) :: rmid
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, k, itemp, jtemp, i1, j1, ns2
C-----------------------------------------------


!
!       THIS ROUTINE FINDS THE LOWER FLUX INDEX [=INDEXS] CORRESPONDING
!       TO THE MEASURED VALUE [=RMEAS] OF (R,Z=0) ALONG THE MIDPLANE
!       (THETA=0 OR PI).
!       THE QUANTITIES RMEAS(K=1,NPTS) ARE INTERPOLATED AS FOLLOWS:
!
!       RMEAS(K) = R[ISMEAS(K),IUMEAS(K)]*[1-DSO(K)] +
!                  R[ISMEAS(K)+1,IUMEAS(K)]*DSO(K)
!
!       BECAUSE OF THE SQRT(S) BEHAVIOUR OF R IN THE FIRST RADIAL ZONE,
!       THE ACTUAL S-INTERPOLAND IN THE FIRST ZONE IS DSO(K)**2 = DSE(K).
!       IN ALL OTHER ZONES, DSE(K) = DSO(K).
!

      IF (npts .le. 0) RETURN
      ns2 = 2*ns
!
!     COMPUTE THE GRID VALUES (S-COORDINATE) OF R ALONG THE MIDPLANE,
!     STARTING AT THETA=PI (I=NTHETHA2) AND ENDING AT THETA=0 (I=1)
!

      rmid(:ns) = reven(indexr(:ns),1,ntheta2) + sqrts(indexr(:ns))
     1   *rodd(indexr(:ns),1,ntheta2)
      rmid(ns+1:ns2) = reven(indexr(ns+1:ns2),1,1) +
     1   sqrts(indexr(ns+1:ns2))*rodd(indexr(ns+1:ns2),1,1)

!
!     FIND THE RADIAL ZONE INDEX [=ITEMP], WHICH BRACKETS THE MEASURED R-VALUE
!
!     RMID(ITEMP-1) .le. RMEAS .le. RMID(ITEMP)
!

      DO k = 1, npts
         itemp = 0
         DO i = 1, ns2 - 1
            IF (rmeas(k) .lt. rmid(i)) THEN
               itemp = i
               GOTO 100
            ENDIF
         END DO
         itemp = ns2
!
!         FIND FLUX-COORDINATE S-INDEX [=ISMEAS], POLOIDAL ANGLE
!         INDEX [=IUMEAS], AND INTERPOLAND [=DSO]
!

  100    CONTINUE
         IF (itemp.gt.1 .and. itemp.lt.ns2) THEN
            i1 = itemp - 1
            jtemp = indexr(itemp)
            j1 = indexr(i1)
            dso(k) = (rmeas(k)-rmid(i1))/(rmid(itemp)-rmid(i1))
            IF (j1 .lt. jtemp) THEN                 !THETA = 0
               ismeas(k) = j1
               iumeas(k) = 1
            ELSE
               ismeas(k) = jtemp
               dso(k) = 1.0 - dso(k)
               iumeas(k) = ntheta2
            ENDIF
         ELSE
            dso(k) = 1.0
            ismeas(k) = indexr(itemp) - 1
!           IGNORE MEASURED POINTS OUTSIDE GRID
            IF (itemp.eq.1 .or. rmeas(k).gt.rmid(ns2-1)) dso(k) = -1.0
            IF (itemp .eq. ns2) iumeas(k) = 1
         ENDIF
!
!        ACCOUNT FOR SQRT(S) SINGULARITY IN 1st ZONE
!        DSE IS THE S-FRACTIONAL INTERPOLAND
!
         dse(k) = dso(k)
         IF (ismeas(k) .eq. 1) dse(k) = dso(k)*dso(k)
      END DO

      END SUBROUTINE findphi
