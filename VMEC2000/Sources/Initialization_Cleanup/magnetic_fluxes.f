      FUNCTION torflux_deriv (x)
      USE stel_kinds
      USE vmec_main, ONLY: zero
      USE vmec_input, ONLY: lRFP, tf => aphi
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), INTENT(IN) :: x
      REAL(rprec) :: torflux_deriv
      REAL(rprec), EXTERNAL :: polflux_deriv, piota
      INTEGER     :: i
C-----------------------------------------------
!     x: radial flux variable (=TOROIDAL FLUX ONLY IF APHI=1)

      IF (lRFP) THEN
!        RFP/TOKAMAK
         IF (piota(x) .eq. zero) STOP 'piota(x) = 0!'
         torflux_deriv = polflux_deriv(x)/piota(x)

      ELSE
!        TOKAMAK/STELLARATOR (default is tf(1) = 1)
         torflux_deriv = 0
         DO i = UBOUND(tf,1), LBOUND(tf,1), -1
            torflux_deriv = x*torflux_deriv + i*tf(i)
         END DO
!        torflux_deriv = 1
      END IF

      END FUNCTION torflux_deriv

      FUNCTION polflux_deriv (x)
      USE stel_kinds
      USE vmec_input, ONLY: lRFP
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), INTENT(IN) :: x
      REAL(rprec) :: tf, polflux_deriv
      REAL(rprec), EXTERNAL :: torflux, torflux_deriv, piota
C-----------------------------------------------
!     x: radial flux variable (=TOROIDAL FLUX ONLY IF APHI=1)
!     polflux_deriv == d(chi)/dx = iota(TF(x)) * torflux_deriv(x)

      IF (lRFP) THEN
!        RFP/TOKAMAK
         polflux_deriv = 1

      ELSE
!        TOKAMAK/STELLARATOR: dchi/ds = iota * dphi/ds
!        piota is assumed to be a function of the TF(x) on input
         tf = torflux(x)
         tf = MIN(tf, 1.0_dp)
         polflux_deriv = piota(tf)*torflux_deriv(x)
!        polflux_deriv = piota(x)*torflux_deriv(x)
      END IF

      END FUNCTION polflux_deriv

      FUNCTION torflux (x)
      USE stel_kinds
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), INTENT(IN)  :: x
      REAL(rprec) :: torflux, h, xi
      REAL(rprec), EXTERNAL :: torflux_deriv
      INTEGER     :: i
C-----------------------------------------------
!     x: radial flux variable (=TOROIDAL FLUX ONLY IF APHI=1)
      h = 1.E-2_dp*x
      torflux = 0
      DO i=1,101
         xi = (i-1)*h
         torflux = torflux + torflux_deriv(xi)
      END DO
      torflux = torflux-0.5_dp*(torflux_deriv(0._dp)+torflux_deriv(x))
      torflux = h*torflux

      END FUNCTION torflux

      FUNCTION polflux (x)
      USE stel_kinds
!     USE vmec_input, ONLY: af => achi
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), INTENT(IN) :: x
      REAL(rprec) :: polflux, h, xi
      REAL(rprec), EXTERNAL :: polflux_deriv
      INTEGER     :: i
C-----------------------------------------------
      h = 1.E-2_dp*x
      polflux = 0
      DO i=1,101
         xi = (i-1)*h
         polflux = polflux + polflux_deriv(xi)
      END DO
      polflux = polflux-0.5_dp*(polflux_deriv(0._dp)+polflux_deriv(x))
      polflux = h*polflux

      END FUNCTION polflux

!  function piota moved to a separate file, piota.f. J Hanson, 2010-03-16
