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
      USE thrift_funcs
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        Rc, w       Model profile coefficients
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: i, ier
      INTEGER :: bcs0(2)
      REAL(rprec) :: Rc, w, Ieccd, Inorm, vp, dPhidrho, temp, &
                     s_val, rho_val
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: dIds_temp, j_temp
      TYPE(EZspline1_r8) :: dIds_spl
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

            ! Calculate dI/ds
            ALLOCATE(dIds_temp(nrho+2))
            dIds_temp(1) = 0.0
            DO i = 2, nrho+1
               dIds_temp(i) = Ieccd*EXP(-(THRIFT_RHO(i-1)-Rc)**2/w**2)/Inorm
            END DO
            dIds_temp(nrho+2) = 0.0
            ! Setup dIds spline
            bcs0=(/ 0, 0/)
            CALL EZspline_init(dIds_spl,nrho+2,bcs0,ier)
            dIds_spl%x1 = THRIFT_RHOFULL
            dIds_spl%isHermite = 1
            CALL EZspline_setup(dIds_spl,dIds_temp,ier,EXACT_DIM=.true.)

            ! Calculate J in s space = dI/ds * 1/(pi*a^2)
            ALLOCATE(j_temp(nsj))
            DO i = 1, nsj
               s_val = THRIFT_S(i)
               rho_val = SQRT(s_val)
               CALL EZspline_interp(dIds_spl,rho_val,temp,ier)
               THRIFT_JECCD(i,mytimestep) = temp/(pi*eq_Aminor**2)
            END DO
            CALL EZspline_free(dIds_spl,ier)

            DEALLOCATE(j_temp)

            IF (lscreen_subcodes) THEN
               WRITE(6,*) '-------------------  ANALYTIC ECCD MODEL  ---------------------'
               WRITE(6,'(A)') '    RHO     J_ECCD'
               DO i = 1, nrho
                  WRITE(6,'(2X,F5.3,1(1X,ES10.2))') THRIFT_RHO(i),THRIFT_JECCD(i,mytimestep)
               END DO
               WRITE(6,*) '-------------------  ANALYTIC ECCD MODEL  ---------------------'
            END IF
         CASE ('travis')
            CALL thrift_paraexe('travis',proc_string,lscreen_subcodes)
      END SELECT

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_run_eccd

