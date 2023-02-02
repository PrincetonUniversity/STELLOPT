!-----------------------------------------------------------------------
!     Subroutine:    thrift_run_eccd
!     Authors:       S. Lazerson
!     Date:          11/24/2022
!     Description:   This subroutine calls the relevant ECCD code.
!-----------------------------------------------------------------------
      SUBROUTINE thrift_run_eccd
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE thrift_runtime
      USE thrift_equil
      USE thrift_vars
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        Rc, w       Model profile coefficients
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: i, ier
      REAL(rprec) :: Rc, w, Ieccd, Inorm, vp, dPhidrho
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      ! Check to make sure we're not zero beta
      IF (eq_beta == 0) THEN
         THRIFT_JECCD(:,mytimestep) = 0
         RETURN
      END IF

      SELECT CASE(TRIM(eccd_type))
         CASE ('model','offaxis','test')
            ! See Turkin, Yu., Maassberg, H., Beidler, C. D., 
            !        Geiger, J. & Marushchenko, N. B. Current Control 
            !        by ECCD for W7-X. Fusion Science and Technology 50,
            !        387–394 (2006).
            Rc = 0.15
            w  = 0.1
            Ieccd = 43E3
            ! From Wolfram
            Inorm = 0.5*w*sqrt(pi)*(ERF((1-Rc)/w)+ERF(rc/w))
            DO i = 1, nrho
               ! dI/ds
               THRIFT_JECCD(i,mytimestep) = Ieccd*EXP(-(THRIFT_RHO(i)-Rc)**2/w**2)/Inorm
               ! But we need j = dI/ds*Aminor/dVds
               ! But we need j = dI/ds * ds/dA
               !               = dI/ds / (dA/ds)
               !               = dI/ds / (dV/ds / (2*pi*Rmajor)) 
               !               = dI/ds * 2 * pi * Rmajor / dV/ds
               CALL EZspline_interp(vp_spl,THRIFT_RHO(i),vp,ier)
               CALL EZspline_interp(phip_spl,THRIFT_RHO(i),dPhidrho,ier) ! dPhi/drho
               vp = vp*dPhidrho/(2*THRIFT_RHO(i)) ! dV/ds = dV/dPhi * dPhi/drho / ds/drho
               THRIFT_JECCD(i,mytimestep) = THRIFT_JECCD(i,mytimestep)*pi2*eq_Rmajor/vp
            END DO
            IF (lscreen_subcodes) THEN
               WRITE(6,*) '-------------------  ANALYTIC ECCD MODEL  ---------------------'
               WRITE(6,'(A)') '    RHO     J_ECCD'
               DO i = 1, nrho
                  WRITE(6,'(2X,F5.3,1(1X,ES10.2))') THRIFT_RHO(i),THRIFT_JECCD(i,mytimestep)
               END DO
               WRITE(6,*) '-------------------  ANALYTIC ECCD MODEL  ---------------------'
            END IF
         CASE ('travis')
            lscreen_subcodes=.TRUE.
            CALL thrift_paraexe('travis',proc_string,lscreen_subcodes)
      END SELECT

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_run_eccd

