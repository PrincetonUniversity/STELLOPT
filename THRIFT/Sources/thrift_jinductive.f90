!-----------------------------------------------------------------------
!     Subroutine:    thrift_jinductive
!     Authors:       L. van Ham
!     Date:          01/23/2023
!     Description:   This subroutine updates the plasma response to
!                    source currents.  
!-----------------------------------------------------------------------
      SUBROUTINE thrift_jinductive
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE thrift_runtime
      USE thrift_vars
      USE thrift_equil
      USE thrift_profiles_mod
      USE stel_tools
      USE EZspline_obj
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: i, ier
      REAL(rprec) :: rho,s,s11,s12,Bav,Bsqav,iota,phip,etapara
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: rho_temp, I_temp, j_temp, f2
      TYPE(EZspline1_r8) :: f_spl
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
!     The general idea here is to use equation (22) in form to update
!     THRIFT_JPLASMA.  Where THRIFT_JSOURCE is the J in <J_s.B>.

!         0   |---|---|---|---|---|---|---|---|   1       rho grid
!         |---|---|---|---|---|---|---|---|---|---|       grid for f
!     |---|---|---|---|---|---|---|---|---|---|---|---|   grid for I
!     dI                                        mu0 I
!     -- = 0                                    ----- = S11 iota + S12
!     ds                                        phip
!
      ! Check to make sure we're not zero beta
      IF (eq_beta == 0) THEN
         IF (mytimestep) == 1 THEN 
            THRIFT_JPLASMA(:,mytimestep) = 0
         ELSE 
            THRIFT_JPLASMA(:,mytimestep) = THRIFT_JPLASMA(:,mytimestep-1)
         END IF
         RETURN
      END IF

      ! Allocations
      ALLOCATE(rho_temp(nrho+4),I_temp(nrho+4),f2(nrho+2),j_temp(nrho+2))
      rho_temp=0; I_temp=0, j_temp=0; f2=0;
   
      ! Temp rho grid extended beyond boundaries (nrho+4)
      ! To get on nrho+2 grid, shift index by +1
      rho_temp(1) = -THRIFT_RHO(1)
      rho_temp(2) = 0
      rho_temp(3:nrho+2) = THRIFT_RHO
      rho_temp(nrho+3) = 1
      rho_temp(nrho+4) = 2*rho_temp(nrho+3) - rho_temp(nrho+2)

      ! Setup grid with mu0 I/phip (nrho+4)
      I_temp(1) = 0 ! rho = -rho(1)
      I_temp(2) = 0 ! rho = 0
      DO i = 1, nrho+1
            rho = THRIFT_RHO(i)
            s   = rho*rho
            ier = 0
            CALL get_equil_sus(s,s11,s12,Bav,Bsqav,ier) 
            CALL EZspline_interp(iota_spl,rho,iota,ier)   
            CALL EZspline_interp(phip_spl,rho,phip,ier)
            I_temp(i+2) = (s11*iota+s12) ! mu0 I/phip|(rho)
      END DO
      I_temp(nrho+4)=I_temp(nrho+3) ! no current outside rho=1

      ! Setup grid with du/drho derivative (nrho+2)
      f2(1) = 0 ! rho = 0
      DO i = 2, nrho+2
         f2(i) = (I_temp(i+2)-I_temp(i))/(rho_temp(i+3)-rho_temp(i+1))/2 
      END DO

      ! Temporary Jsource grid with boundaries (nrho+2)
      j_temp(1)          = THRIFT_JSOURCE(1,mytimestep)
      j_temp(2:nrho+1)   = THRIFT_JSOURCE(:,mytimestep)
      j_temp(nrho+2)     = 0.0

      ! Construct f (nrho+2)
      DO i = 1, nrho+2
         rho = rho_temp(i+1)
         s   = rho*rho
         ier = 0
         CALL get_prof_etapara(rho,THRIFT_T(mytimestep),etapara) !eta = eta
         CALL EZspline_interp(vp_spl,rho,phip,ier) ! phip = V'
         CALL get_prof_pprime(rho,THRIFT_T(mytimestep),s12) !s12 = p'
         CALL get_equil_Bav(s,Bav,Bsqav,ier) ! Bav =<B>;Bsqav=<B^2>
         ! I on nrho+4 grid, f on nrho+2, so I_temp is shifted by +1
         f2(i) = etapara*phip*(s12*I_temp(i+1)     &      ! eta*V'(p'u  
                        + Bavsq/mu0*f2(i)  &          !     + <B^2>/mu0*du/dx 
                        j_temp(i)*Bav)                !         - <J.B>)
      END DO

      ! I_temp (nrho+4), j_temp(nrho+2) freed up
      I_temp = 0; j_temp = 0;
      
      ! Construct du/dt (nrho)
      DO i = 1, nrho
         rho = THRIFT_RHO(i)
         s   = rho*rho
         ier = 0
         CALL get_equil_sus(s,s11,s12,Bav,Bsqav,ier) 
         j_temp(i) = s11/(eq_phiedge**2) * &
            (f2(i+2)-f2(i))/(rho_temp(i+3)-rho_temp(i+1))/2 ! s11/phi_a^2*df/dx = du/dt
         j_temp(i) = phip/mu0*j_temp(i) ! j_temp= dI/dt, assuming d/dt(phip)=0
      END DO
      ! Note that 

      ! Need dj/dt though, not dI/dt
      ! dj/dt = pi R/(rho V') d/drho(dI/dt)
      ! Once again, need to evaluate d/drho using CD scheme (nrho+2)
      ! BCs: dI/dt(0) = 0, dI/dt(rho) = dI/dt(prev rho)
      I_temp(1) = 0 
      I_temp(2:nrho) = j_temp
      I_temp(nrho+1) = I_temp(nrho)
      ! Note that I_temp(nrho+3)=I_temp(nrho+4)=0 are never accessed
      
      ! j_temp freed up again
      j_temp = 0

      ! Calculate dj/dt (nrho)
      DO i = 1, nrho
         rho = THRIFT_RHO(i)
         CALL EZspline_interp(vp_spl,rho,phip,ier)
         j_temp(i) = (I_temp(i+2)-I_temp(i))/(rho_temp(i+3)-rho_temp(i+1))/2 ! d/rho(dI/dt)
         j_temp(i) = j_temp(i)*pi*eq_Rmajor/(rho*phip) ! dj/dt
      END DO     

      ! update plasma current
      ! s11 = dt
      IF (mytimestep==1) THEN 
         s11 = THRIFT_T(2)-THRIFT_T(1)
         THRIFT_JPLASMA(:,mytimestep) = j_temp*s11
      ELSE 
         s11 = THRIFT_T(mytimestep)-THRIFT_T(mytimestep-1)
         THRIFT_JPLASMA(:,mytimestep) = THRIFT_JPLASMA(:,mytimestep-1)+j_temp*s11
      END IF

      ! Deallocate
      DEALLOCATE(f1,f2,f3,f4)

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_jinductive

