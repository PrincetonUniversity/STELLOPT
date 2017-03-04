      SUBROUTINE init_scale(r_scale, b_scale)
      USE vmec_input, ONLY: rprec, curtor, am, raxis_cc, zaxis_cs,
     1     raxis_cs, zaxis_cc, phiedge, extcur, rbc, zbs, dp, lfreeb
      USE optim_params, ONLY: Target_MaxCurrent, lcoil_geom
      USE coilsnamin, ONLY: cursad, curmod, bcoil_cur, I_pol, cc_vf
      USE mpi_params                                                     !MPI
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), INTENT(inout) :: r_scale, b_scale
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1
C-----------------------------------------------
      IF (rbc(0,0) .le. zero) RETURN
      IF (r_scale .le. zero) r_scale = 1

      IF (r_scale.ne.one .and. myid.eq.master) THEN
         PRINT *,'  R  being scaled by ', r_scale
         IF (lfreeb) PRINT *,
     1   ' ***CAUTION*** Free-boundary mode: dimensions in mgrid file',
     2   ' may no longer be invalid!'
      END IF

      IF (b_scale.ne.one .and. myid.eq.master)                           !MPI
     1   PRINT *,' |B| being scaled by ', b_scale

      phiedge = r_scale*r_scale*b_scale * phiedge                        !Scale |B| by b_scale

      rbc = r_scale*rbc
      zbs = r_scale*zbs

      raxis_cc = r_scale*raxis_cc; raxis_cs = r_scale*raxis_cs
      zaxis_cs = r_scale*zaxis_cs; zaxis_cc = r_scale*zaxis_cc

      Target_MaxCurrent = r_scale*b_scale*Target_MaxCurrent
      curtor = r_scale*b_scale*curtor
      extcur = r_scale*b_scale*extcur
      am = b_scale*b_scale * am

!
!     Scale currents in COILSIN namelist
!
      IF (lcoil_geom) THEN
         cursad = r_scale*b_scale*cursad
         curmod = r_scale*b_scale*curmod
         bcoil_cur = r_scale*b_scale*bcoil_cur
         I_pol = r_scale*b_scale*I_pol
         cc_vf = r_scale*b_scale*cc_vf
      END IF
!
!     Reset scale factors to unity after scaling performed
!
      r_scale = 1
      b_scale = 1

      END SUBROUTINE init_scale
