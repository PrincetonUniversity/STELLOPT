!-----------------------------------------------------------------------
!     Subroutine:    stelltran_bal_update
!     Authors:       J. Mittelstaedt (jmittelstaedt@uchicago.edu)
!     Date:          07/14/2016
!     Description:   This subroutine updates ne, using the particle
!                    and power balance equations.
!-----------------------------------------------------------------------
      SUBROUTINE stelltran_bal_update(itime)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE stelltran_runtime
      USE EZspline_obj
      USE EZspline
      USE stelltran_equilutils
      USE stelltran_vars
      USE stelltran_data
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
      INTEGER, INTENT(IN) :: itime
!       REAL(rprec), DIMENSION(prof_length) :: rho_coord, gradrho, Vps
!       REAL(rprec), DIMENSION(nne) :: rho_ne
!       REAL(rprec), DIMENSION(nte) :: rho_te
!       REAL(rprec), DIMENSION(nti) :: rho_ti
!       REAL(rprec) :: dt_update, s ! ********** dt IS ONLY TEMPORARY
!       INTEGER :: j, ier
! 
!       ! LSODE Variables
!       INTEGER :: neq, itol, itask, istate, iopt, lrw, liw, mf
!       REAL(rprec) :: t, tout, rtol, atol
!       REAL(rprec) :: init_sol(prof_length*3)
!       REAL(rprec) :: rwork(22+27*prof_length+9*prof_length**2), iwork(20+2*prof_length)
!       REAL(rprec), EXTERNAL :: bal_update_rhs, bal_update_jac
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      print*, '********* in the balance updating thing ***************************'
      ! NOTE: Need to change how it gets power density, or leave this until we figure it out. Still uses
      !       pe_rf (now dPdV), should either switch to Ptot (or derived pow_density) or just switch back
      !       after we figure out what was going on with pe_rf in the first place. If we use pow_density,
      !       should probably have it calculated once (maybe in ECRH) instead of calculated in s_balance and
      !       s_bal_update.
!       dt_update = 0.300 ! Setting this to 1ms at first, should change when this is properly incorporated
!       DO j=1,prof_length
!             s = real(j-1)/real(prof_length-1) ! sqrt to go from s to rho
!             S_pe(j,2) = 1.D26*EXP(-s**2/1.0D-2)
!       END DO
! 
!       ! Find needed geometric factors.
!       DO j = 1, prof_length
!          ! Given S this function returns rho,Vp,<|grad(rho)|>,<|grad(rho)|^2>
!          CALL get_equil_rho(te(j,1),rho_coord(j),Vps(j),gradrho(j),gradrhosq(j,2),ier)
!       END DO
!       DO j=1,nne
!             CALL eval_prof_spline(prof_length,ne(:,1),rho_coord(:),ne_s(j),rho_ne(j),ier)
!       END DO
!       DO j=1,nte
!             CALL eval_prof_spline(prof_length,te(:,1),rho_coord(:),te_s(j),rho_te(j),ier)
!       END DO
!       DO j=1,nti
!             CALL eval_prof_spline(prof_length,ti(:,1),rho_coord(:),ti_s(j),rho_ti(j),ier)
!       END DO
!       print*, 'here a'
!       ! Convert dV/dPhi to dV/drho
!       Vp(:,2) = 2.0*rho_coord(:)*Vps(:)
!       Vp(1,2) = Vp(2,2) ! TENTATIVE while we decide what to do about division by zero stuff
! 
!       !Set up for LSODE
!       init_sol(1:prof_length) = ne(:,2)
!       init_sol(1+prof_length:2*prof_length) = te(:,2)
!       init_sol(1+2*prof_length:3*prof_length) = ti(:,2)
!       neq = 3*prof_length ! # of equations
!       itol = 1 ! Whether we want one abs tolerace for all or to specify it for each eqn
!       rtol = 1.D-3 ! relative error tolerance. Use 0 if only absolute
!       atol = 0.0 ! absolute error tolerance. Use 0 if only relative
!       itask = 1 ! Don't change this. Tells it to do normal output.
!       istate = 1 ! ~ error flag. 1 is first run, 2 is subsequent, negative is an error.
!       iopt = 0 ! If we want to use optional arguments
!       t = 0.
!       tout = dt_update
!       lrw = 22 + 9*neq + neq**2 ! length of real work array
!       liw = 20+neq ! length of integer work array
!       mf = 21 ! Method of solution
! 
!       print*,'here e new'
! 
!       ! Save the initial profiles
!       CALL save_ST_param2(13,'./te_bali.txt',prof_length,rho_coord(:),te(:,2))
!       CALL save_ST_param2(13,'./ti_bali.txt',prof_length,rho_coord(:),ti(:,2))
!       CALL save_ST_param2(13,'./ne_bali.txt',prof_length,rho_coord(:),ne(:,2))
! 
!       print*, 'here f'
! 
!       ! Ready to solve it (I hope...)
!       CALL DLSODE(bal_update_rhs,neq,init_sol,t,tout,itol,rtol,atol,itask,istate,iopt,rwork,lrw,iwork,liw,bal_update_jac,mf)
! 
!       print*, 'here g'
! 
!       DO j=1,nne
!             CALL eval_prof_spline(prof_length,rho_coord(:),init_sol(1:prof_length),rho_ne(j),ne_f(j,itime+1),ier)
!       END DO
!       DO j=1,nte
!             CALL eval_prof_spline(prof_length,rho_coord(:),&
!                   init_sol(prof_length+1:2*prof_length),rho_te(j),te_f(j,itime+1),ier)
!       END DO
!       DO j=1,nti
!             CALL eval_prof_spline(prof_length,rho_coord(:),&
!                   init_sol(2*prof_length+1:3*prof_length),rho_ti(j),ti_f(j,itime+1),ier)
!       END DO
!       ! Should need to update zeff somehow but idk how. Would need to come from maybe some particle source things we haven't really
!       ! included yet??
! 
!       print*,'here h'
! 
!       ! Save the final profiles and the coefficients
!       CALL save_ST_param2(12,'./Xe_bal.txt',prof_length,Xe(:,1),Xe(:,2))
!       CALL save_ST_param2(12,'./Xi_bal.txt',prof_length,Xi(:,1),Xi(:,2))
!       CALL save_ST_param2(12,'./De_bal.txt',prof_length,De(:,1),De(:,2))
!       CALL save_ST_param2(13,'./te_balf.txt',prof_length,rho_coord(:),init_sol(1:prof_length))
!       CALL save_ST_param2(13,'./ti_balf.txt',prof_length,rho_coord(:),init_sol(prof_length+1:2*prof_length))
!       CALL save_ST_param2(13,'./ne_balf.txt',prof_length,rho_coord(:),init_sol(2*prof_length+1:3*prof_length))
! 
!       print*,' here at end of predictive'

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE stelltran_bal_update

!       SUBROUTINE bal_update_rhs(neq,t,y,ydot)
!             ! SECOND ORDER OPTION
!             USE stel_kinds, ONLY: rprec
!             USE stelltran_equilutils
!             USE stelltran_vars
!             IMPLICIT NONE
!             INTEGER, INTENT(IN) :: neq
!             REAL(rprec), INTENT(IN) :: t, y(neq)
!             REAL(rprec), INTENT(OUT) :: ydot(neq)
!             INTEGER :: i
!             REAL(rprec) :: dr, temp(neq/3), temp2(neq/3), temp3(neq/3), picol(neq/3), mi
!             dr = 1./(prof_length-1) ! assuming all equally spaced
!             mi = 2000.*me
! 
!             ! n_e system
!             CALL arr_fderiv2(prof_length,Vp(:,2)*gradrhosq(:,2)*De(:,2),dr,temp)
!             CALL arr_fderiv2(prof_length,y(1:prof_length),dr,temp2)
!             CALL arr_2fderiv2(prof_length,y(1:prof_length),dr,temp3)
!             ydot(1:prof_length) = S_pe(:,2) - (temp(:)*temp2(:) + Vp(:,2)*gradrhosq(:,2)*De(:,2)*temp3(:))/Vp(:,2)
! 
!             CALL calc_pi_coll2(prof_length,y(1:prof_length),y(1+prof_length:2*prof_length),y(1+2*prof_length:3*prof_length),zeff(:,2), mi, picol)
! 
!             ! t_e system
!             CALL arr_fderiv2(prof_length,Vp(:,2)*gradrhosq(:,2)*y(1:prof_length)*Xe(:,2),dr,temp)
!             CALL arr_fderiv2(prof_length,y(1+prof_length:2*prof_length),dr,temp2)
!             CALL arr_2fderiv2(prof_length,y(1+prof_length:2*prof_length),dr,temp3)
!             ydot(1+prof_length:2*prof_length) = 2.*(pe_rf(:,2) - picol(:) - (temp(:)*temp2(:) & 
!                   + Vp(:,2)*gradrhosq(:,2)*y(1:prof_length)*Xe(:,2)*temp3)/Vp(:,2))/(3.*y(1:prof_length)) &
!                   - y(1+prof_length:2*prof_length)*ydot(1:prof_length)/y(1:prof_length)
! 
!             ! t_i system, assuming zeff is constant in time
!             CALL arr_fderiv2(prof_length,Vp(:,2)*gradrhosq(:,2)*y(1:prof_length)*Xi(:,2)/zeff(:,2),dr,temp)
!             CALL arr_fderiv2(prof_length,y(1+2*prof_length:3*prof_length),dr,temp2)
!             CALL arr_2fderiv2(prof_length,y(1+2*prof_length:3*prof_length),dr,temp3)
!             ydot(1+2*prof_length:3*prof_length) = 2.*zeff(:,2)*(picol(:) - (temp(:)*temp2(:) & 
!                   + Vp(:,2)*gradrhosq(:,2)*y(1:prof_length)*Xi(:,2)*temp3/zeff(:,2))/Vp(:,2))/(3.*y(1:prof_length)) &
!                   - y(1+2*prof_length:3*prof_length)*ydot(1:prof_length)/y(1:prof_length)
!       END SUBROUTINE bal_update_rhs
! 
!       SUBROUTINE bal_update_jac(neq, t, y, ml, mu, pd, nrowpd)
!             USE stel_kinds, ONLY: rprec
!             USE stelltran_equilutils
!             USE stelltran_vars
!             IMPLICIT NONE
!             INTEGER, INTENT(IN) :: neq, ml, mu, nrowpd
!             REAL(rprec), INTENT(IN) :: t, y(neq)
!             REAL(rprec), INTENT(OUT) :: pd(nrowpd,neq)
!             REAL(rprec), DIMENSION(neq/3) :: Ge_jac, Qe_mn, Qi_mn, DpcDti, DpcDte, DpcDn, fn, temp, temp2, temp3
!             REAL(rprec), DIMENSION(neq/3) :: picol, temp4, temp5
!             REAL(rprec) :: col_pow_const, mi, temp1, dr
!             INTEGER :: j
!             ! Just need to load pd with w/e we need
!             dr = 1./(prof_length-1) ! assuming all equally spaced
!             mi = 2000.*me ! Making up mi while we don't have a correct value
!             Ge_jac = Vp(:,2)*gradrhosq(:,2)*De(:,2)
!             Qe_mn = Vp(:,2)*gradrhosq(:,2)*Xe(:,2)
!             Qi_mn = Vp(:,2)*gradrhosq(:,2)*Xi(:,2)
!             CALL arr_fderiv2(prof_length,Vp(:,2)*gradrhosq(:,2)*De(:,2),dr,temp)
!             CALL arr_fderiv2(prof_length,y(1:prof_length),dr,temp2)
!             CALL arr_2fderiv2(prof_length,y(1:prof_length),dr,temp3)
!             CALL calc_pi_coll2(prof_length,y(1:prof_length),y(1+prof_length:2*prof_length),y(1+2*prof_length:3*prof_length),zeff(:,2), mi, picol)
!             fn = S_pe(:,2) - (temp(:)*temp2(:) + Vp(:,2)*gradrhosq(:,2)*De(:,2)*temp3(:))/Vp(:,2)
!             col_pow_const = 3.*me*ec/(mi*3.44D11)
!             DpcDti = col_pow_const*y(1:prof_length)**2.*(LOG(y(1:prof_length)**0.5/y(1+prof_length:2*prof_length)) - 24.&
!                               - LOG(1.D3))/y(1+prof_length:2*prof_length)**1.5
!             DpcDn = 2.*y(1:prof_length)*col_pow_const*(y(1+prof_length:2*prof_length)-y(1+2*prof_length:3*prof_length))*(23.25 + LOG(1.D3) &
!                               - LOG(y(1:prof_length)**0.5/y(1+prof_length:2*prof_length)))/y(1+prof_length:2*prof_length)**1.5
!             DpcDte = col_pow_const*y(1:prof_length)**2.*(2.*(y(1+prof_length:2*prof_length)-y(1+2*prof_length:3*prof_length))&
!                               + (24. + LOG(1.D3) -  LOG(y(1:prof_length)**0.5/y(1+prof_length:2*prof_length)))&
!                               *(3.*y(1+2*prof_length:3*prof_length)-y(1+prof_length:2*prof_length)))/(2.*y(1+prof_length:2*prof_length)**2.5)
! 
!             ! left endpoints
!             CALL fin_deriv2(Ge_jac(1:3),dr,'l',temp1)
!             pd(1,1) = -1.*(-3.*temp1/(2.*dr) + 2.*Ge_jac(1)/(dr**2.))/Vp(1,2)
!             pd(1,2) = -1.*(2.*temp1/(dr) - 5.*Ge_jac(1)/(dr**2.))/Vp(1,2)
!             pd(1,3) = -1.*(-1.*temp1/(2.*dr) + 4.*Ge_jac(1)/(dr**2.))/Vp(1,2)
!             pd(1,4) = -1.*(-1.*Ge_jac(1)/(dr**2.))/Vp(1,2)
!             CALL fin_deriv2(Qe_mn(1:3)*y(1:3),dr,'l',temp1)
!             pd(prof_length+1,prof_length+1) = -2.*(DpcDte(1) + (-3.*temp1/(2.*dr) + 2.*Qe_mn(1)*y(1)/(dr**2.))/Vp(1,2))/(3.*y(1)) - fn(1)/y(1)
!             pd(prof_length+1,prof_length+2) = -2.*((2.*temp1/(dr) - 5.*Qe_mn(1)*y(1)/(dr**2.))/Vp(1,2))/(3.*y(1))
!             pd(prof_length+1,prof_length+3) = -2.*((-1.*temp1/(2.*dr) + 4.*Qe_mn(1)*y(1)/(dr**2.))/Vp(1,2))/(3.*y(1))
!             pd(prof_length+1,prof_length+4) = -2.*((-1.*Qe_mn(1)*y(1)/(dr**2.))/Vp(1,2))/(3.*y(1))
!             CALL fin_deriv2(Qi_mn(1:3)*y(1:3)/zeff(1:3,2),dr,'l',temp1)
!             pd(prof_length+1,prof_length+1) = 2.*zeff(1,2)*(DpcDti(1) - (-3.*temp1/(2.*dr) + 2.*Qi_mn(1)*y(1)/(zeff(1,2)*dr**2))/Vp(1,2))/(3.*y(1)) - fn(1)/y(1)
!             pd(prof_length+1,prof_length+2) = 2.*zeff(1,2)*((2.*temp1/(dr) - 5.*Qi_mn(1)*y(1)/(zeff(1,2)*dr**2.))/Vp(1,2))/(3.*y(1))
!             pd(prof_length+1,prof_length+3) = 2.*zeff(1,2)*((-1.*temp1/(2.*dr) + 4.*Qi_mn(1)*y(1)/(zeff(1,2)*dr**2.))/Vp(1,2))/(3.*y(1))
!             pd(prof_length+1,prof_length+4) = 2.*zeff(1,2)*((-1.*Qi_mn(1)*y(1)/(zeff(1,2)*dr**2.))/Vp(1,2))/(3.*y(1))
!             ! Diagonal entries
!             DO j=2,prof_length-1
!                   ! ne
!                   pd(j,j) = 2.*Ge_jac(j)/(dr**2.*Vp(j,2))
!                   ! Te
!                   pd(j+prof_length,j+prof_length) = -2.*(DpcDte(j) - 2.*Qe_mn(j)*y(j)/(dr**2.*Vp(j,2)))/(3.*y(j)) - fn(j)/y(j)
!                   ! Ti
!                   pd(j+2*prof_length,j+2*prof_length) = 2.*zeff(j,2)*(DpcDti(j) +2.*Qi_mn(j)*y(j)/(dr**2.*zeff(j,2)*Vp(j,2)))/(3.*y(j)) - fn(j)/y(j)
!             END DO
!             ! right endpoints
!             CALL fin_deriv2(Ge_jac(prof_length-2:prof_length),dr,'r',temp1)
!             pd(prof_length,prof_length) = -1.*(3.*temp1/(2.*dr) + 2.*Ge_jac(prof_length)/(dr**2.))/Vp(prof_length,2)
!             pd(prof_length,prof_length-1) = -1.*(-2.*temp1/(dr) - 5.*Ge_jac(prof_length)/(dr**2.))/Vp(prof_length,2)
!             pd(prof_length,prof_length-2) = -1.*(1.*temp1/(2.*dr) + 4.*Ge_jac(prof_length)/(dr**2.))/Vp(prof_length,2)
!             pd(prof_length,prof_length-3) = -1.*(-1.*Ge_jac(prof_length)/(dr**2.))/Vp(prof_length,2)
!             CALL fin_deriv2(Qe_mn(prof_length-2:prof_length)*y(prof_length-2:prof_length),dr,'r',temp1)
!             pd(2*prof_length,2*prof_length) = -2.*(DpcDte(prof_length) + (3.*temp1/(2.*dr)&
!                                     + 2.*Qe_mn(prof_length)*y(prof_length)/(dr**2.))/Vp(prof_length,2))/(3.*y(prof_length))&
!                                     - fn(prof_length)/y(prof_length)
!             pd(2*prof_length,2*prof_length-1) = -2.*(-2.*temp1/(dr) - 5.*Qe_mn(prof_length)*y(prof_length)/(dr**2.))/(3.*Vp(prof_length,2)*y(prof_length))
!             pd(2*prof_length,2*prof_length-2) = -2.*(1.*temp1/(2.*dr) + 4.*Qe_mn(prof_length)*y(prof_length)/(dr**2.))/(3.*Vp(prof_length,2)*y(prof_length))
!             pd(2*prof_length,2*prof_length-3) = -2.*(-1.*Qe_mn(prof_length)*y(prof_length)/(dr**2.))/(3.*Vp(prof_length,2)*y(prof_length))
!             CALL fin_deriv2(Qi_mn(prof_length-2:prof_length)*y(prof_length-2:prof_length)/zeff(prof_length-2:prof_length,2),dr,'l',temp1)
!             pd(3*prof_length,3*prof_length) = 2.*zeff(prof_length,2)*(DpcDti(prof_length) - (3.*temp1/(2.*dr)&
!                                     + 2.*Qi_mn(prof_length)*y(prof_length)/(zeff(prof_length,2)*dr**2.))&
!                                     /Vp(prof_length,2))/(3.*y(prof_length)) - fn(prof_length)/y(prof_length)
!             pd(3*prof_length,3*prof_length-1) = 2.*zeff(prof_length,2)*(-2.*temp1/dr&
!                                     - 5.*Qi_mn(prof_length)*y(prof_length)/(zeff(prof_length,2)*dr**2.))/(3.*Vp(prof_length,2)*y(prof_length))
!             pd(3*prof_length,3*prof_length-2) = 2.*zeff(prof_length,2)*(temp1/(2.*dr)& 
!                                     + 4.*Qi_mn(prof_length)*y(prof_length)/(zeff(prof_length,2)*dr**2.))&
!                                     /(3.*y(prof_length)*Vp(prof_length,2))
!             pd(3*prof_length,3*prof_length-3) = -2.*Qi_mn(prof_length)/(3.*Vp(prof_length,2)*dr**2.)
! 
!             ! Super/sub diagonal Entries
!             DO j=2,prof_length-1
!                   CALL fin_deriv2(Ge_jac(j-1:j+1),dr,'m',temp1)
!                   ! ne
!                   pd(j,j+1) = -1.*(temp1/(2.*dr) + Ge_jac(j)/dr**2.)/Vp(j,2)
!                   pd(j,j-1) = (temp1/(2.*dr) - Ge_jac(j)/dr**2.)/Vp(j,2)
!                   ! te
!                   CALL fin_deriv2(Qe_mn(j-1:j+1)*y(j-1:j+1),dr,'m',temp1)
!                   pd(j+prof_length,j+prof_length+1) = -2.*(temp1/(2.*dr) + Qe_mn(j)*y(j)/dr**2.)/(3.*Vp(j,2)*y(j))
!                   pd(j+prof_length,j+prof_length-1) = -2.*(-1.*temp1/(2.*dr) + Qe_mn(j)*y(j)/dr**2.)/(3.*Vp(j,2)*y(j))
!                   ! ti
!                   CALL fin_deriv2(Qi_mn(j-1:j+1)*y(j-1:j+1)/zeff(j-1:j+1,2),dr,'m',temp1)
!                   pd(j+2*prof_length,j+2*prof_length+1) = -2.*zeff(j,2)*(temp1/(2.*dr) + Qe_mn(j)*y(j)/(zeff(j,2)*dr**2.))/(3.*Vp(j,2)*y(j))
!                   pd(j+2*prof_length,j+2*prof_length-1) = -2.*zeff(j,2)*(-1.*temp1/(2.*dr)+ Qe_mn(j)*y(j)/(zeff(j,2)*dr**2.))/(3.*Vp(j,2)*y(j))
!             END DO
! 
!             ! cross terms
!             CALL arr_fderiv2(prof_length,y(prof_length+1:prof_length*2),dr,temp4)
!             CALL arr_fderiv2(prof_length,y(prof_length*2+1:prof_length*3),dr,temp5)
!             DO j=1,prof_length
!                   ! te/n
!                   pd(j+prof_length,j) = 2.*((picol(j) - pe_rf(j,2))/y(j) - DpcDn(j)&
!                                     + Qe_mn(j)*temp2(j)*temp4(j)/(Vp(j,2)*y(j)))/(3.*y(j))&
!                                     + y(j+prof_length)*(fn(j)/y(j) - pd(j,j))/y(j)
!                   ! ti/n
!                   pd(j+2*prof_length,j) = 2.*zeff(j,2)*(DpcDn(j) - picol(j)/y(j)&
!                                     + Qi_mn(j)*temp2(j)*temp5(j)/(Vp(j,2)*y(j)*zeff(j,2)))/(3.*y(j))&
!                                     + y(j+2*prof_length)*(fn(j)/y(j) - pd(j,j))/y(j)
!                   ! te/ti
!                   pd(j+prof_length,j+2*prof_length) = -2.*DpcDti(j)/(3.*y(j))
!                   ! ti/te
!                   pd(j+2*prof_length,j+prof_length) = 2.*zeff(j,2)*DpcDte(j)/(3.*y(j))
!             END DO
!       END SUBROUTINE bal_update_jac
! 
!       SUBROUTINE calc_pi_coll2(numpts,n_e,t_e,t_i,z_eff,m_i,p_i)
!             USE stel_kinds, ONLY: rprec
!             USE stelltran_vars, ONLY: ec, me
!             IMPLICIT NONE
!             INTEGER, INTENT(IN) :: numpts
!             REAL(rprec), INTENT(IN) :: m_i
!             REAL(rprec), DIMENSION(numpts), INTENT(IN) :: n_e, t_e, t_i, z_eff
!             REAL(rprec), DIMENSION(numpts), INTENT(OUT) :: p_i
!             REAL(rprec), DIMENSION(numpts) :: lambda, tau_e
!             INTEGER :: j
!             ! NOTE: Formulary eqns are in cgs so need to do converting
!             ! Determine appropriate Coulomb Logarithm, NRL Formulary p. 34
!             DO j=1,numpts
!                   IF (t_e(j) < 10.*z_eff(j)**2.) THEN
!                         lambda(j) = 23. - LOG(SQRT(n_e(j)/1.D6)*z_eff(j)/t_e(j)**(1.5))
!                   ELSE
!                         lambda(j) = 24. - LOG(SQRT(n_e(j)/1.D6)/t_e(j))
!                   END IF
!             END DO
!             ! Find electron -> ion energy from collisions. NRL Formulary p. 37
!             tau_e = 3.44D5*t_e(:)**1.5/(lambda(:)*(n_e(:)/1.D6))
!             p_i = 3.0*me*n_e(:)*ec*(t_e(:)-t_i(:))/(m_i*tau_e(:)) ! Should be good to not use cgs here since tau is in sec, ec is for eV -> J conversion
!       END SUBROUTINE calc_pi_coll2
! 
!       SUBROUTINE save_ST_param2(str_len,fname_saving,npts,indep,dep)
!             INTEGER, INTENT(IN) :: npts, str_len
!             CHARACTER(str_len), INTENT(IN) :: fname_saving
!             REAL*8, DIMENSION(npts), INTENT(IN) :: indep, dep
!             INTEGER :: j
!             IF (access(fname_saving,' ') == 0) THEN
!                   OPEN(UNIT=12, FILE=fname_saving, ACTION="write", STATUS="old",POSITION="append")
!                   do j=1,npts
!                         write(12,*) dep(j), indep(j)
!                   end do
!                   CLOSE(UNIT=12)
!             ELSE
!                   OPEN(UNIT=12, FILE=fname_saving, ACTION="write", STATUS="new")
!                   write(12,'(I2.2)') npts
!                   do j=1,npts
!                         write(12,*) dep(j), indep(j)
!                   end do
!                   CLOSE(UNIT=12)
!             END IF
!       END SUBROUTINE save_ST_param2















