      SUBROUTINE becoil (rad, zee, br, bp, bz, brvac, bpvac, bzvac, &
     &                   lscreen)
      USE vparams, ONLY: nthreed
      USE vacmod
      USE parallel_include_module
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), DIMENSION(nuv3), INTENT(in) :: rad, zee
      REAL(dp), DIMENSION(nuv3), INTENT(out) :: br, bp, bz
      REAL(dp), DIMENSION(nr0b,nz0b,np0b), INTENT(in) ::
     1   brvac, bpvac, bzvac
      LOGICAL :: lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      CHARACTER(LEN=50), PARAMETER :: warning =
     1   'Plasma Boundary exceeded Vacuum Grid Size'
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, SAVE :: icount = 0
      INTEGER :: i, kv, ir, jz, ir1, jz1
      REAL(dp) :: rad0, zee0, ri, zj,
     1   pr, qz, w22, w21, w12, w11, tbecon, tbecoff
C-----------------------------------------------
!
!     DETERMINE THE CYLINDRICAL COMPONENTS OF THE EXTERNAL
!     MAGNETIC FIELD (BR, BP, BZ) AT A FIXED PHI PLANE BY 
!     USING 2-D INTERPOLATION BASED ON THE FOUR POINT FORMULA
!     IN ABRAMOWITZ AND STEGUN, EQ. 25.2.66 
!
!     BRVAC, BPVAC, BZVAC:       CYLINDRICAL COMPONENTS OF VACUUM B-FIELD 
!                                STORED ON R, Z, PHI GRID
!
!     RAD, ZEE:                  ARRAY OF R, Z VALUES (ON FLUX SURFACE)
!                                AT WHICH B-FIELD IS TO BE DETERMINED
!
C-----------------------------------------------
      CALL second0(tbecon)

      icount = icount + 1

      DO i = nuv3min, nuv3max
!
!       CHECK THAT BOUNDARY POINTS ARE INSIDE VACUUM GRID.  IF NOT,
!       SET THEM EQUAL TO LARGEST (OR SMALLEST) VACUUM GRID POINTS
!
         rad0 = MIN(rad(i), rmaxb)
         rad0 = MAX(rad0,   rminb)
         zee0 = MIN(zee(i), zmaxb)
         zee0 = MAX(zee0,   zminb)
!
!       DETERMINE PHI-PLANE, KV (MUST LIE IN FIRST FIELD PERIOD)
!
         kv = 1 + MOD(i - 1,nv)
         kv = MIN (kv, np0b)                       !!Axi-symmetric special case
!
!
!       DETERMINE INTEGER INDICES (IR,JZ) FOR LOWER LEFT R, Z CORNER GRID POINT
!
         ir = INT((rad0 - rminb)/delrb) + 1
         jz = INT((zee0 - zminb)/delzb) + 1
         ir1 = MIN(nr0b,ir + 1)
         jz1 = MIN(nz0b,jz + 1)
!
!       COMPUTE RI, ZJ AND PR , QZ AT GRID POINT (IR , JZ)
!       ALSO, COMPUTE WEIGHTS WIJ FOR 4 CORNER GRID POINTS
!
         ri = rminb + (ir - 1)*delrb
         zj = zminb + (jz - 1)*delzb
         pr = (rad0 - ri)/delrb
         qz = (zee0 - zj)/delzb
         w22 = pr*qz                   !p*q               
         w21 = pr - w22                !p*(1-q)
         w12 = qz - w22                !q*(1-p)
         w11 = 1 + w22 - (pr + qz)     !(1-p)*(1-q)
!
!       COMPUTE B FIELD AT R, PHI, Z BY INTERPOLATION
!
         br(i) = w11*brvac(ir,jz,kv) + w22*brvac(ir1,jz1,kv) +
     1      w21*brvac(ir1,jz,kv) + w12*brvac(ir,jz1,kv)
         bz(i) = w11*bzvac(ir,jz,kv) + w22*bzvac(ir1,jz1,kv) +
     1      w21*bzvac(ir1,jz,kv) + w12*bzvac(ir,jz1,kv)
         bp(i) = w11*bpvac(ir,jz,kv) + w22*bpvac(ir1,jz1,kv) +
     1      w21*bpvac(ir1,jz,kv) + w12*bpvac(ir,jz1,kv)

      END DO

!
!     PRINT INFO IF R, Z OUT OF BOX
!
      IF (MOD(icount,25).EQ.0 .AND. rank.EQ.0) THEN
         i = 0
         rad0 = MAXVAL(rad)
         zee0 = MAXVAL(zee)
         IF (rad0 .gt. rmaxb) i = 1
         IF (zee0 .gt. zmaxb) i = i + 2
         ri   = MINVAL(rad)
         zj   = MINVAL(zee)
         IF (ri .lt. rminb) i = i + 4
         IF (zj .lt. zminb) i = i + 8
         IF (i .ne. 0 .and. lscreen) THEN
            PRINT *, warning
            WRITE (nthreed, *) warning
            IF (i/8 .ne. 0) PRINT *,' zmin = ', zj
            i = MOD(i,8)
            IF (i/4 .ne. 0) PRINT *,' rmin = ', ri
            i = MOD(i,4)
            IF (i/2 .ne. 0) PRINT *,' zmax = ', zee0
            i = MOD(i,2)
            IF (i .ne. 0) PRINT *,' rmax = ', rad0
         END IF
      ENDIF

      CALL second0(tbecoff)
      becoil_time = becoil_time + (tbecoff - tbecon)

      END SUBROUTINE becoil
