      SUBROUTINE update_style (lmask, surf_mask, nrad)
      USE stel_kinds, ONLY: rprec
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: nrad
      REAL(rprec) :: surf_mask(nrad)
      LOGICAL :: lmask(nrad)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: k, js
      LOGICAL :: l1
C-----------------------------------------------
      l1 = ALL(.not.lmask)               !!Check if user is using old-style input
      IF (l1) THEN                       !!New style, find lsurf_mask array
         DO k = 1, nrad
            DO js = 2, nrad
            IF (js .eq. (1+ NINT(surf_mask(k)*(nrad-1))))
     1            lmask(js) = .true.
            END DO
         END DO
      ELSE                               !!Old style, convert lsurf to nsurf
         k = 1
         DO js = 2, nrad
            IF (lmask(js)) THEN
              surf_mask(k) = REAL(js - 1,rprec)/(nrad-1)
              k = k+1
            END IF
         END DO
      END IF
      lmask(1) = .false.           !!MUST ignore axis point in boozer code

      END SUBROUTINE update_style
