!-----------------------------------------------------------------------
!     Subroutine:    stelltran_bal_update
!     Authors:       J. Mittelstaedt (jmittelstaedt@uchicago.edu)
!     Date:          06/21/2016
!     Description:   This subroutine updates ne, te, ti using the particle
!                    and power balance equations.
!-----------------------------------------------------------------------
      SUBROUTINE stelltran_bal_update
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
!       USE stel_kinds, ONLY: rprec
!       USE stelltran_runtime
!       USE EZspline_obj
!       USE EZspline
!       USE stelltran_equilutils
!       USE stelltran_vars
!-----------------------------------------------------------------------
!     Local Variables
!     RHS               Right hand side of the differential equation
!     init_nete         Initial value of n_e*T_e for electron power balance update
!     init_niti         Initial value of n_i*T_i for ion power balance update
!     tau_e             Collisional time for electrons for ion collisional heating
!     lambda            Coulomb logarithm
!     dr                Spatial resolution of our grid.
!-----------------------------------------------------------------------
      IMPLICIT NONE
!       REAL(rprec), DIMENSION(prof_length) :: RHS_parte, RHS_parti, RHS_powi, RHS_powe, init_nete
!       REAL(rprec), DIMENSION(prof_length) :: init_niti, tau_e, lambda, rho_coord, gradrho, gradrhosq
!       REAL(rprec) :: dr, mi, dt ! ********** dt IS ONLY TEMPORARY
!       INTEGER :: i, ier
! !-----------------------------------------------------------------------
! !     Begin Subroutine
! !-----------------------------------------------------------------------
!       print*, '********* in the balance updating thing ***************************'
!       dt = 0.001 ! Setting this to 1ms at first, should change when this is properly incorporated
!       mi = 2000.*me ! tentative until we have mi as an input parameter
!       ni(:,2) = ne(:,2) / zeff(:,2) ! Assuming quasi-neutrality for now
! 
!       ! Find needed geometric factors.
!       DO i = 1, prof_length
!          ! Given S this function returns rho,Vp,<|grad(rho)|>,<|grad(rho)|^2>
!          CALL get_equil_rho(Vp(i,1),rho_coord(i),Vp(i,2),gradrho(i),gradrhosq(i),ier)
!       END DO
! 
! 
!       dr = Vp(2,1) - Vp(1,1) ! find what the spatial resolution is
! 
!       ! NOTE: Should this be included in the ECRH subroutine?
!       ! Determine appropriate lambda, NRL Formulary p. 34
!       DO i=1,prof_length
!             IF (te(i,2) < 10*zeff(i,2)**2.) THEN
!                   lambda(i) = 23. - LOG(ne(i,2)**0.5*zeff(i,2)/te(i,2)**(1.5))
!             ELSE
!                   lambda(i) = 24. - LOG(ne(i,2)**0.5/te(i,2))
!             END IF
!       END DO
!       ! NRL Formulary p. 37
!       tau_e = (3.*me**0.5*te(:,2)**1.5)/(4*pi2**0.5*ne(:,2)*lambda*ec**2.5) !divide by ec^2.5 to convert eV to J
!       pi_col(:,2) = 3.0*me*ne(:,2)*(te(:,2)-ti(:,2))/(mi*tau_e)
! 
!       ! *************** BAD this is derivative wrt s not rho ************************************
!       ! Would probably just be easiest to take derivative of the spline since easy conversion will
!       ! leave us with the function defined at uneven points, which is difficult to do finite difference
!       ! method on. Same applies to s_bal_update
! 
!       ! Left endpoint second order finite difference derivative
!       RHS_parte(1) = S_pe(1,2) - (4.*Vp(2,2)*Ge(2,2) - 3.*Vp(1,2)*Ge(1,2) - Vp(3,2)*Ge(3,2))/(2.*dr*Vp(1,2))
!       RHS_parti(1) = S_pi(1,2) - (4.*Vp(2,2)*Gi(2,2) - 3.*Vp(1,2)*Gi(1,2) - Vp(3,2)*Gi(3,2))/(2.*dr*Vp(1,2))
!       RHS_powe(1) = pe_rf(1,2) - pi_col(1,2) - (4.*Vp(2,2)*Qe(2,2) - 3.*Vp(1,2)*Qe(1,2) - Vp(3,2)*Qe(3,2))/(2.*dr*Vp(1,2)) - Ge(1,2)*Er(1,2)
!       RHS_powi(1) = pi_col(1,2) - (4.*Vp(2,2)*Qi(2,2) - 3.*Vp(1,2)*Qi(1,2) - Vp(3,2)*Qi(3,2))/(2.*dr*Vp(1,2)) + zeff(1,2)*Gi(1,2)*Er(1,2)
! 
!       ! Right endpoint second order finite difference derivative
!       RHS_parte(prof_length) = S_pe(prof_length,2) - (3.*Vp(prof_length,2)*Ge(prof_length,2) - 4.*Vp(prof_length-1,2)*Ge(prof_length-1,2) + Vp(prof_length-2,2)*Ge(prof_length-2,2))/(2.*dr*Vp(prof_length,2))
!       RHS_parti(prof_length) = S_pi(prof_length,2) - (3.*Vp(prof_length,2)*Gi(prof_length,2) - 4.*Vp(prof_length-1,2)*Gi(prof_length-1,2) + Vp(prof_length-2,2)*Gi(prof_length-2,2))/(2.*dr*Vp(prof_length,2))
!       RHS_powe(prof_length) = pe_rf(prof_length,2) - pi_col(prof_length,2) - (3.*Vp(prof_length,2)*Qe(prof_length,2) - 4.*Vp(prof_length-1,2)*Qe(prof_length-1,2) + Vp(prof_length-2,2)*Qe(prof_length-2,2))/(2.*dr*Vp(prof_length,2)) - Ge(prof_length,2)*Er(prof_length,2)
!       RHS_powi(prof_length) = pi_col(prof_length,2) - (3.*Vp(prof_length,2)*Qi(prof_length,2) - 4.*Vp(prof_length-1,2)*Qi(prof_length-1,2) + Vp(prof_length-2,2)*Qi(prof_length-2,2))/(2.*dr*Vp(prof_length,2)) + zeff(prof_length,2)*Gi(prof_length,2)*Er(prof_length,2)
! 
! 
!       ! Middle second order finite difference derivative
!       DO i=2,prof_length-1
!             RHS_parte(i) = S_pe(i,2) - (Vp(i+1,2)*Ge(i+1,2) - Vp(i-1,2)*Ge(i-1,2))/(2.*dr*Vp(i,2))
!             RHS_parti(i) = S_pi(i,2) - (Vp(i+1,2)*Gi(i+1,2) - Vp(i-1,2)*Gi(i-1,2))/(2.*dr*Vp(i,2))
!             RHS_powe(i) = pe_rf(i,2) - pi_col(i,2) - (Vp(i+1,2)*Qe(i+1,2) - Vp(i-1,2)*Qe(i-1,2))/(2.*dr*Vp(i,2)) - Ge(i,2)*Er(i,2)
!             RHS_powi(i) = pi_col(i,2) - (Vp(i+1,2)*Qi(i+1,2) - Vp(i-1,2)*Qi(i-1,2))/(2.*dr*Vp(i,2)) + zeff(i,2)*Gi(i,2)*Er(i,2)
!       END DO
! 
!       init_nete = ne(:,2)*te(:,2)
!       init_niti = ni(:,2)*ti(:,2)
!       ! updating the density in the simplest way.
!       ne(:,2) = ne(:,2) + RHS_parte*dt
!       ni(:,2) = ni(:,2) + RHS_parti*dt
!       Te(:,2) = (init_nete + RHS_powe*dt)/ne(:,2)
!       Ti(:,2) = (init_niti + RHS_powi*dt)/ni(:,2)

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE stelltran_bal_update