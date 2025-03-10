      SUBROUTINE loadparams (xc, nvar)
      USE modular_coils
      USE saddle_coils
      USE bcoils_mod
      USE saddle_surface
      USE tf_coils
      USE vf_coils
      USE Vcoilpts
      USE coils
      USE mpi_params                                         !mpi stuff
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: nvar
      REAL(rprec), TARGET :: xc(nvar)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: n, nvariables
      REAL(rprec), POINTER :: xvariables(:)
!-----------------------------------------------

!     Load modular coefficients

      n = 0
      IF (lmodular) THEN
         xvariables => xc(n+1:)
         CALL load_modular_coils (nvariables, xvariables)
         n = n + nvariables


         IF (lmodcur) xvariables => xc(n+1:)
!        Need to CALL load_modular_currents even if lmodcur = F, in case
!        TF current has changed
         CALL load_modular_currents (nvariables, xvariables)
         n = n + nvariables

         CALL load_modular_structures
      END IF

!     Load saddle coil coefficients

      IF (lsaddle) THEN
         xvariables => xc(n+1:)
         CALL load_saddle_coils (nvariables, xvariables)
         n = n + nvariables

!        For now, can not vary TF current when lsaddle = T, unless
!        saddle currents also vary (need to fix)
         IF (lsadcur) THEN
            xvariables => xc(n+1:)
            CALL load_saddle_currents (nvariables, xvariables)
            n = n + nvariables
         END IF

         CALL load_saddle_structures
      END IF

!     Load background coil currents (IF lbcoil = true)

      IF (lbcoil) THEN
         xvariables => xc(n+1:)
         CALL load_bg_currents (nvariables, xvariables)
         n = n + nvariables
      END IF

!     Load vf coil currents

      IF (lvf) THEN
         IF (lvfvar) THEN
            xvariables => xc(n+1:)
            CALL load_vf_coils (nvariables, xvariables)
            n = n + nvariables
         END IF

         xvariables => xc(n+1:)
         CALL load_vf_currents (nvariables, xvariables)
         n = n + nvariables
!
         CALL load_vf_structures
      END IF

!     Load tf coil currents

      IF (ltfcv) THEN
         xvariables => xc(n+1:)
         CALL load_tf_coils (nvariables, xvariables)
         n = n + nvariables
      END IF

!     Load modular winding surface coefficients

      IF (lsurfv) THEN
         xvariables => xc(n+1:)
         CALL load_modular_wsurf (nvariables, xvariables)
         n = n + nvariables
      END IF

!     Load saddle winding surface coefficients

      IF (lsadsfv) THEN
         xvariables => xc(n+1:)
         CALL load_saddle_wsurf (nvariables, xvariables)
         n = n + nvariables
      END IF

      IF (myid .EQ. master) THEN
        IF (n .NE. nvar) STOP 'loadparams: n .ne. nvar'
      END IF

      END SUBROUTINE loadparams
