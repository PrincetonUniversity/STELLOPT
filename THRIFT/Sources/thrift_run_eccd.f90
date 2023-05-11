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
      INTEGER :: i, ier, i1,i2
      INTEGER :: bcs0(2)
      REAL(rprec) :: Rc, w, Ieccd, Inorm, vp, dPhidrho, temp, &
                     s_val, rho_val
      REAL(rprec), DIMENSION(:), ALLOCATABLE ::  j_temp
      TYPE(EZspline1_r8) :: j_spl
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      ! Check to make sure we're not zero beta
      IF (eq_beta == 0) THEN
         THRIFT_JECCD(:,mytimestep) = 0
         RETURN
      END IF

      ! Get Power at timestep
      i1 = COUNT(PECRH_AUX_T < THRIFT_T(mytimestep))
      i2 = i1+1
      IF (PECRH_AUX_F(i1)<=0 .or. PECRH_AUX_F(i2)<=0) THEN
         THRIFT_JECCD(:,mytimestep) = 0
         RETURN
      END IF

      ! Set Ieccd here so that in the future we can use this
      ! to control TRAVIS ECCD by the same code.
      Ieccd =    ( PECRH_AUX_F(i2)      - PECRH_AUX_F(i1) ) &
               * ( THRIFT_T(mytimestep) - PECRH_AUX_T(i1) ) &
               / ( PECRH_AUX_T(i2)      - PECRH_AUX_T(i1) ) &
               + PECRH_AUX_F(i1)


      SELECT CASE(TRIM(eccd_type))
         CASE ('model','offaxis','test','simple')
            ! See Turkin, Yu., Maassberg, H., Beidler, C. D., 
            !        Geiger, J. & Marushchenko, N. B. Current Control 
            !        by ECCD for W7-X. Fusion Science and Technology 50,
            !        387–394 (2006).
            Rc = 0.15
            w  = 0.1

            ! From Wolfram
            Inorm = 0.5*w*( SQRT(pi)*Rc*( ERF((1-Rc)/w) + ERF(Rc/w) )+w*( EXP(-Rc**2/w**2) - EXP(-(Rc-1)**2/w**2) ))
            Inorm = Ieccd/(2*pi*eq_Aminor**2*Inorm)

            ! Calculate J (rho-space)
            ALLOCATE(j_temp(nrho+2))
            DO i = 1, nrho+2
               j_temp(i) = Inorm*EXP( -(THRIFT_RHOFULL(i)-Rc)**2/w**2 )
            END DO

            ! Setup J spline
            bcs0=(/ 0, 0/)
            CALL EZspline_init(j_spl,nrho+2,bcs0,ier)
            j_spl%x1 = THRIFT_RHOFULL
            j_spl%isHermite = 1
            CALL EZspline_setup(j_spl,j_temp,ier,EXACT_DIM=.true.)
            DEALLOCATE(j_temp)

            ! Interpolate J in s-space
            DO i = 1, nsj
               s_val = THRIFT_S(i)
               rho_val = SQRT(s_val)
               CALL EZspline_interp(j_spl,rho_val,THRIFT_JECCD(i,mytimestep),ier)
            END DO
            CALL EZspline_free(j_spl,ier)


            IF (lscreen_subcodes) THEN
               WRITE(6,*) '-------------------  ANALYTIC ECCD MODEL  ---------------------'
               WRITE(6,'(A)') '    S     J_ECCD'
               DO i = 1, nsj
                  WRITE(6,'(2X,F5.3,1(1X,ES10.2))') THRIFT_S(i),THRIFT_JECCD(i,mytimestep)
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

