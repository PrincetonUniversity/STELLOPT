!-----------------------------------------------------------------------
!     Subroutine:    stelltran_bal_update
!     Authors:       J. Mittelstaedt (jmittelstaedt@uchicago.edu)
!     Date:          07/14/2016
!     Description:   This subroutine updates ne, using the particle
!                    and power balance equations.
!-----------------------------------------------------------------------
      SUBROUTINE stelltran_bal_update
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE stelltran_runtime
      USE EZspline_obj
      USE EZspline
      USE stelltran_equilutils
      USE stelltran_vars
      USE stelltran_rho_vars
!-----------------------------------------------------------------------
!     Local Variables
!     RHS               Right hand side of the differential equation
!     init_nete         Initial value of n_e*T_e for electron power balance update
!     init_niti         Initial value of n_i*T_i for ion power balance update
!     tau_e             Collisional time for electrons for ion collisional heating
!     lambda            Coulomb logarithm
!     dr                Spatial resolution of our grid.
!     num_spat_sol      Number of rho coordinates to include when solving our ode. Will end up with 3*num_spat_sol equations
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), DIMENSION(prof_length) :: Vps
      REAL(rprec) :: dt ! ********** dt IS ONLY TEMPORARY
      INTEGER :: j, ier, num_spat_sol
      INTEGER :: neq, itol, itask, istate, iopt, lrw, liw, mf
      INTEGER, DIMENSION(:), ALLOCATABLE :: iwork
      REAL(rprec) :: t, tout, rtol, atol, s
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: rwork, spat_sol_rho, init_sol
      REAL(rprec), EXTERNAL :: bal_update_rhs, bal_update_jac
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      print*, '********* in the balance updating thing ***************************'
      num_spat_sol = 20 ! Would need to be able to set this elsewhere but will put it here for now.
      dt = 0.300 ! Setting this to 1ms at first, should change when this is properly incorporated
      DO j=1,prof_length
            s = real(j-1)/real(prof_length-1) ! sqrt to go from s to rho
            S_pe(j,2) = 1.D15*EXP(-s**2/1.0D-1)
      END DO


      ! Find needed geometric factors.
      DO j = 1, prof_length
         ! Given S this function returns rho,Vp,<|grad(rho)|>,<|grad(rho)|^2>
         CALL get_equil_rho(Vp(j,1),rho_coord(j,2),Vps(j),gradrho(j,2),gradrhosq(j,2),ier)
      END DO

      ! Convert dV/dPhi to dV/drho
      Vp(:,2) = 2*rho_coord(:,2)*Vps(:)
      Vp(1,2) = Vp(2,2) ! TENTATIVE while we decide what to do about division by zero stuff

      ! Put variables onto uniform array in rho of specified size
      ALLOCATE(spat_sol_rho(1:num_spat_sol))
      DO j=1,num_spat_sol
            spat_sol_rho(j) = (j-1)/(num_spat_sol-1)
      END DO
      
      ALLOCATE(ne_rho(1:num_spat_sol))
      ALLOCATE(te_rho(1:num_spat_sol))
      ALLOCATE(ti_rho(1:num_spat_sol))
      ALLOCATE(De_rho(1:num_spat_sol))
      ALLOCATE(Xe_rho(1:num_spat_sol))
      ALLOCATE(Xi_rho(1:num_spat_sol))
      ALLOCATE(Vp_rho(1:num_spat_sol))
      ALLOCATE(Spe_rho(1:num_spat_sol))
      ALLOCATE(zeff_rho(1:num_spat_sol))
      ALLOCATE(perf_rho(1:num_spat_sol))
      ALLOCATE(gradrhosq_rho(1:num_spat_sol))
      ! NOTE: Pretty sure from implementation of eval_prof_spline it is more efficient to have all in their own loops
      DO j=1,num_spat_sol
            CALL eval_prof_spline(prof_length,rho_coord(:,2),ne(:,2),spat_sol_rho(j),ne_rho(j),ier)
      END DO
      DO j=1,num_spat_sol
            CALL eval_prof_spline(prof_length,rho_coord(:,2),te(:,2),spat_sol_rho(j),te_rho(j),ier)
      END DO
      DO j=1,num_spat_sol
            CALL eval_prof_spline(prof_length,rho_coord(:,2),ti(:,2),spat_sol_rho(j),ti_rho(j),ier)
      END DO
      DO j=1,num_spat_sol
            CALL eval_prof_spline(prof_length,rho_coord(:,2),De(:,2),spat_sol_rho(j),De_rho(j),ier)
      END DO
      DO j=1,num_spat_sol
            CALL eval_prof_spline(prof_length,rho_coord(:,2),Xe(:,2),spat_sol_rho(j),Xe_rho(j),ier)
      END DO
      DO j=1,num_spat_sol
            CALL eval_prof_spline(prof_length,rho_coord(:,2),Xi(:,2),spat_sol_rho(j),Xi_rho(j),ier)
      END DO
      DO j=1,num_spat_sol
            CALL eval_prof_spline(prof_length,rho_coord(:,2),Vp(:,2),spat_sol_rho(j),Vp_rho(j),ier)
      END DO
      DO j=1,num_spat_sol
            CALL eval_prof_spline(prof_length,rho_coord(:,2),S_pe(:,2),spat_sol_rho(j),Spe_rho(j),ier)
      END DO
      DO j=1,num_spat_sol
            CALL eval_prof_spline(prof_length,rho_coord(:,2),Zeff(:,2),spat_sol_rho(j),zeff_rho(j),ier)
      END DO
      DO j=1,num_spat_sol
            CALL eval_prof_spline(prof_length,rho_coord(:,2),pe_rf(:,2),spat_sol_rho(j),perf_rho(j),ier)
      END DO
      DO j=1,num_spat_sol
            CALL eval_prof_spline(prof_length,rho_coord(:,2),gradrhosq(:,2),spat_sol_rho(j),gradrhosq_rho(j),ier)
      END DO

      !Set up for LSODEA
      ALLOCATE(init_sol(1:num_spat_sol*3))
      init_sol(1:num_spat_sol) = ne_rho(:)
      init_sol(1+num_spat_sol:2*num_spat_sol) = te_rho(:)
      init_sol(1+2*num_spat_sol:3*num_spat_sol) = ti_rho(:)
      neq = 3*num_spat_sol ! # of equations
      itol = 1 ! Whether we want one abs tolerace for all or to specify it for each eqn
      rtol = 1.D-3 ! relative error tolerance. Use 0 if only absolute
      atol = 0.0 ! absolute error tolerance. Use 0 if only relative
      itask = 1 ! Don't change this. Tells it to do normal output.
      istate = 1 ! ~ error flag. 1 is first run, 2 is subsequent, negative is an error.
      iopt = 0 ! If we want to use optional arguments
      t = 0.
      tout = dt
      lrw = 22 + 9*neq + neq**2 ! length of real work array
      liw = 20+neq ! length of integer work array
      mf = 21 ! Method of solution
      ALLOCATE(iwork(1:liw))
      ALLOCATE(rwork(1:lrw))

      ! Ready to solve it (I hope...)
      CALL DLSODE(bal_update_rhs,neq,init_sol,t,tout,itol,rtol,atol,itask,istate,iopt,rwork,lrw,iwork,liw,bal_update_jac,mf)

      ! Save the initial profiles
      CALL save_ST_param('./te_bali.txt',prof_length,te(:,1),te(:,2))
      CALL save_ST_param('./ti_bali.txt',prof_length,ti(:,1),ti(:,2))
      CALL save_ST_param('./ne_bali.txt',prof_length,ne(:,1),ne(:,2))

      ! Update the parameters in the correct coordinate system
      DO j=1,prof_length
            CALL eval_prof_spline(num_spat_sol,spat_sol_rho(:),init_sol(1:num_spat_sol),rho_coord(j,2),ne(j,2),ier)
      END DO
      DO j=1,prof_length
            CALL eval_prof_spline(num_spat_sol,spat_sol_rho(:),init_sol(num_spat_sol+1:2*num_spat_sol),rho_coord(j,2),te(j,2),ier)
      END DO
      DO j=1,prof_length
            CALL eval_prof_spline(num_spat_sol,spat_sol_rho(:),init_sol(2*num_spat_sol+1:3*num_spat_sol),rho_coord(j,2),ti(j,2),ier)
      END DO

      ! Save the final profiles and the coefficients
      CALL save_ST_param('./Xe_bal.txt',prof_length,Xe(:,1),Xe(:,2))
      CALL save_ST_param('./Xi_bal.txt',prof_length,Xi(:,1),Xi(:,2))
      CALL save_ST_param('./De_bal.txt',prof_length,De(:,1),De(:,2))
      CALL save_ST_param('./te_balf.txt',prof_length,te(:,1),te(:,2))
      CALL save_ST_param('./ti_balf.txt',prof_length,ti(:,1),ti(:,2))
      CALL save_ST_param('./ne_balf.txt',prof_length,ne(:,1),ne(:,2))

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE stelltran_bal_update

      SUBROUTINE bal_update_rhs(neq,t,y,ydot)
            ! SECOND ORDER OPTION
            USE stel_kinds, ONLY: rprec
            USE stelltran_equilutils
            USE stelltran_rho_vars
            USE stelltran_vars, ONLY: me ! eventually mi if we have it
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: neq
            REAL(rprec), INTENT(IN) :: t, y(neq)
            REAL(rprec), INTENT(OUT) :: ydot(neq)
            INTEGER :: wpflen, i
            REAL(rprec) :: dr, temp(neq/3), temp2(neq/3), temp3(neq/3), picol(neq/3), mi
            wpflen = neq/3 ! working number of spatial points
            dr = 1./(wpflen-1) ! assuming all equally spaced
            mi = 2000.*me
      
            ! n_e system
            CALL arr_fderiv2(wpflen,Vp_rho(:)*gradrhosq_rho(:)*De_rho(:),dr,temp)
            CALL arr_fderiv2(wpflen,y(1:wpflen),dr,temp2)
            CALL arr_2fderiv2(wpflen,y(1:wpflen),dr,temp3)
            ydot(1:wpflen) = Spe_rho(:) - (temp(:)*temp2(:) + Vp_rho(:)*gradrhosq_rho(:)*De_rho(:)*temp3(:))/Vp_rho(:)
      
            CALL calc_pi_coll2(wpflen,y(1:wpflen),y(1+wpflen:2*wpflen),y(1+2*wpflen:3*wpflen),zeff_rho(:), mi, picol)
      
            ! t_e system
            CALL arr_fderiv2(wpflen,Vp_rho(:)*gradrhosq_rho(:)*y(1:wpflen)*Xe_rho(:),dr,temp)
            CALL arr_fderiv2(wpflen,y(1+wpflen:2*wpflen),dr,temp2)
            CALL arr_2fderiv2(wpflen,y(1+wpflen:2*wpflen),dr,temp3)
            ydot(1+wpflen:2*wpflen) = 2.*(perf_rho(:) - picol(:) - (temp(:)*temp2(:) & 
                  + Vp_rho(:)*gradrhosq_rho(:)*y(1:wpflen)*Xe_rho(:)*temp3)/Vp_rho(:))/(3.*y(1:wpflen)) &
                  - y(1+wpflen:2*wpflen)*ydot(1:wpflen)/y(1:wpflen)
      
            ! t_i system, assuming zeff is constant in time
            CALL arr_fderiv2(wpflen,Vp_rho(:)*gradrhosq_rho(:)*y(1:wpflen)*Xi_rho(:)/zeff_rho(:),dr,temp)
            CALL arr_fderiv2(wpflen,y(1+2*wpflen:3*wpflen),dr,temp2)
            CALL arr_2fderiv2(wpflen,y(1+2*wpflen:3*wpflen),dr,temp3)
            ydot(1+2*wpflen:3*wpflen) = 2.*zeff_rho(:)*(picol(:) - (temp(:)*temp2(:) & 
                  + Vp_rho(:)*gradrhosq_rho(:)*y(1:wpflen)*Xi_rho(:)*temp3/zeff_rho(:))/Vp_rho(:))/(3.*y(1:wpflen)) &
                  - y(1+2*wpflen:3*wpflen)*ydot(1:wpflen)/y(1:wpflen)
      END SUBROUTINE bal_update_rhs

      SUBROUTINE bal_update_jac(neq, t, y, ml, mu, pd, nrowpd)
            USE stel_kinds, ONLY: rprec
            USE stelltran_equilutils
            USE stelltran_rho_vars
            USE stelltran_vars, ONLY: ec, me ! Eventually mi if we have it
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: neq, ml, mu, nrowpd
            REAL(rprec), INTENT(IN) :: t, y(neq)
            REAL(rprec), INTENT(OUT) :: pd(nrowpd,neq)
            REAL(rprec), DIMENSION(neq/3) :: Ge_jac, Qe_mn, Qi_mn, DpcDti, DpcDte, DpcDn, fn, temp, temp2, temp3
            REAL(rprec), DIMENSION(neq/3) :: picol, temp4, temp5
            REAL(rprec) :: col_pow_const, mi, temp1, dr
            INTEGER :: wpflen, j
            ! Just need to load pd with w/e we need
            wpflen = neq/3
            dr = 1./(wpflen-1) ! assuming all equally spaced
            mi = 2000.*me ! Making up mi while we don't have a correct value
            Ge_jac = Vp_rho(:)*gradrhosq_rho(:)*De_rho(:)
            Qe_mn = Vp_rho(:)*gradrhosq_rho(:)*Xe_rho(:)
            Qi_mn = Vp_rho(:)*gradrhosq_rho(:)*Xi_rho(:)
            CALL arr_fderiv2(wpflen,Vp_rho(:)*gradrhosq_rho(:)*De_rho(:),dr,temp)
            CALL arr_fderiv2(wpflen,y(1:wpflen),dr,temp2)
            CALL arr_2fderiv2(wpflen,y(1:wpflen),dr,temp3)
            CALL calc_pi_coll2(wpflen,y(1:wpflen),y(1+wpflen:2*wpflen),y(1+2*wpflen:3*wpflen),zeff_rho(:), mi, picol)
            fn = Spe_rho(:) - (temp(:)*temp2(:) + Vp_rho(:)*gradrhosq_rho(:)*De_rho(:)*temp3(:))/Vp_rho(:)
            col_pow_const = 3.*me*ec/(mi*3.44D11)
            DpcDti = col_pow_const*y(1:wpflen)**2.*(LOG(y(1:wpflen)**0.5/y(1+wpflen:2*wpflen)) - 24.&
                              - LOG(1.D3))/y(1+wpflen:2*wpflen)**1.5
            DpcDn = 2.*y(1:wpflen)*col_pow_const*(y(1+wpflen:2*wpflen)-y(1+2*wpflen:3*wpflen))*(23.25 + LOG(1.D3) &
                              - LOG(y(1:wpflen)**0.5/y(1+wpflen:2*wpflen)))/y(1+wpflen:2*wpflen)**1.5
            DpcDte = col_pow_const*y(1:wpflen)**2.*(2.*(y(1+wpflen:2*wpflen)-y(1+2*wpflen:3*wpflen))&
                              + (24. + LOG(1.D3) -  LOG(y(1:wpflen)**0.5/y(1+wpflen:2*wpflen)))&
                              *(3.*y(1+2*wpflen:3*wpflen)-y(1+wpflen:2*wpflen)))/(2.*y(1+wpflen:2*wpflen)**2.5)

            ! left endpoints
            CALL fin_deriv2(Ge_jac(1:3),dr,'l',temp1)
            pd(1,1) = -1.*(-3.*temp1/(2.*dr) + 2.*Ge_jac(1)/(dr**2.))/Vp_rho(1)
            pd(1,2) = -1.*(2.*temp1/(dr) - 5.*Ge_jac(1)/(dr**2.))/Vp_rho(1)
            pd(1,3) = -1.*(-1.*temp1/(2.*dr) + 4.*Ge_jac(1)/(dr**2.))/Vp_rho(1)
            pd(1,4) = -1.*(-1.*Ge_jac(1)/(dr**2.))/Vp_rho(1)
            CALL fin_deriv2(Qe_mn(1:3)*y(1:3),dr,'l',temp1)
            pd(wpflen+1,wpflen+1) = -2.*(DpcDte(1) + (-3.*temp1/(2.*dr) + 2.*Qe_mn(1)*y(1)/(dr**2.))/Vp_rho(1))/(3.*y(1)) - fn(1)/y(1)
            pd(wpflen+1,wpflen+2) = -2.*((2.*temp1/(dr) - 5.*Qe_mn(1)*y(1)/(dr**2.))/Vp_rho(1))/(3.*y(1))
            pd(wpflen+1,wpflen+3) = -2.*((-1.*temp1/(2.*dr) + 4.*Qe_mn(1)*y(1)/(dr**2.))/Vp_rho(1))/(3.*y(1))
            pd(wpflen+1,wpflen+4) = -2.*((-1.*Qe_mn(1)*y(1)/(dr**2.))/Vp_rho(1))/(3.*y(1))
            CALL fin_deriv2(Qi_mn(1:3)*y(1:3)/zeff_rho(1:3),dr,'l',temp1)
            pd(wpflen+1,wpflen+1) = 2.*zeff_rho(1)*(DpcDti(1) - (-3.*temp1/(2.*dr) + 2.*Qi_mn(1)*y(1)/(zeff_rho(1)*dr**2))/Vp_rho(1))/(3.*y(1)) - fn(1)/y(1)
            pd(wpflen+1,wpflen+2) = 2.*zeff_rho(1)*((2.*temp1/(dr) - 5.*Qi_mn(1)*y(1)/(zeff_rho(1)*dr**2.))/Vp_rho(1))/(3.*y(1))
            pd(wpflen+1,wpflen+3) = 2.*zeff_rho(1)*((-1.*temp1/(2.*dr) + 4.*Qi_mn(1)*y(1)/(zeff_rho(1)*dr**2.))/Vp_rho(1))/(3.*y(1))
            pd(wpflen+1,wpflen+4) = 2.*zeff_rho(1)*((-1.*Qi_mn(1)*y(1)/(zeff_rho(1)*dr**2.))/Vp_rho(1))/(3.*y(1))
            ! Diagonal entries
            DO j=2,wpflen-1
                  ! ne
                  pd(j,j) = 2.*Ge_jac(j)/(dr**2.*Vp_rho(j))
                  ! Te
                  pd(j+wpflen,j+wpflen) = -2.*(DpcDte(j) - 2.*Qe_mn(j)*y(j)/(dr**2.*Vp_rho(j)))/(3.*y(j)) - fn(j)/y(j)
                  ! Ti
                  pd(j+2*wpflen,j+2*wpflen) = 2.*zeff_rho(j)*(DpcDti(j) +2.*Qi_mn(j)*y(j)/(dr**2.*zeff_rho(j)*Vp_rho(j)))/(3.*y(j)) - fn(j)/y(j)
            END DO
            ! right endpoints
            CALL fin_deriv2(Ge_jac(wpflen-2:wpflen),dr,'r',temp1)
            pd(wpflen,wpflen) = -1.*(3.*temp1/(2.*dr) + 2.*Ge_jac(wpflen)/(dr**2.))/Vp_rho(wpflen)
            pd(wpflen,wpflen-1) = -1.*(-2.*temp1/(dr) - 5.*Ge_jac(wpflen)/(dr**2.))/Vp_rho(wpflen)
            pd(wpflen,wpflen-2) = -1.*(1.*temp1/(2.*dr) + 4.*Ge_jac(wpflen)/(dr**2.))/Vp_rho(wpflen)
            pd(wpflen,wpflen-3) = -1.*(-1.*Ge_jac(wpflen)/(dr**2.))/Vp_rho(wpflen)
            CALL fin_deriv2(Qe_mn(wpflen-2:wpflen)*y(wpflen-2:wpflen),dr,'r',temp1)
            pd(2*wpflen,2*wpflen) = -2.*(DpcDte(wpflen) + (3.*temp1/(2.*dr)&
                                    + 2.*Qe_mn(wpflen)*y(wpflen)/(dr**2.))/Vp_rho(wpflen))/(3.*y(wpflen))&
                                    - fn(wpflen)/y(wpflen)
            pd(2*wpflen,2*wpflen-1) = -2.*(-2.*temp1/(dr) - 5.*Qe_mn(wpflen)*y(wpflen)/(dr**2.))/(3.*Vp_rho(wpflen)*y(wpflen))
            pd(2*wpflen,2*wpflen-2) = -2.*(1.*temp1/(2.*dr) + 4.*Qe_mn(wpflen)*y(wpflen)/(dr**2.))/(3.*Vp_rho(wpflen)*y(wpflen))
            pd(2*wpflen,2*wpflen-3) = -2.*(-1.*Qe_mn(wpflen)*y(wpflen)/(dr**2.))/(3.*Vp_rho(wpflen)*y(wpflen))
            CALL fin_deriv2(Qi_mn(wpflen-2:wpflen)*y(wpflen-2:wpflen)/zeff_rho(wpflen-2:wpflen),dr,'l',temp1)
            pd(3*wpflen,3*wpflen) = 2.*zeff_rho(wpflen)*(DpcDti(wpflen) - (3.*temp1/(2.*dr)&
                                    + 2.*Qi_mn(wpflen)*y(wpflen)/(zeff_rho(wpflen)*dr**2.))&
                                    /Vp_rho(wpflen))/(3.*y(wpflen)) - fn(wpflen)/y(wpflen)
            pd(3*wpflen,3*wpflen-1) = 2.*zeff_rho(wpflen)*(-2.*temp1/dr&
                                    - 5.*Qi_mn(wpflen)*y(wpflen)/(zeff_rho(wpflen)*dr**2.))/(3.*Vp_rho(wpflen)*y(wpflen))
            pd(3*wpflen,3*wpflen-2) = 2.*zeff_rho(wpflen)*(temp1/(2.*dr)& 
                                    + 4.*Qi_mn(wpflen)*y(wpflen)/(zeff_rho(wpflen)*dr**2.))&
                                    /(3.*y(wpflen)*Vp_rho(wpflen))
            pd(3*wpflen,3*wpflen-3) = -2.*Qi_mn(wpflen)/(3.*Vp_rho(wpflen)*dr**2.)

            ! Super/sub diagonal Entries
            DO j=2,wpflen-1
                  CALL fin_deriv2(Ge_jac(j-1:j+1),dr,'m',temp1)
                  ! ne
                  pd(j,j+1) = -1.*(temp1/(2.*dr) + Ge_jac(j)/dr**2.)/Vp_rho(j)
                  pd(j,j-1) = (temp1/(2.*dr) - Ge_jac(j)/dr**2.)/Vp_rho(j)
                  ! te
                  CALL fin_deriv2(Qe_mn(j-1:j+1)*y(j-1:j+1),dr,'m',temp1)
                  pd(j+wpflen,j+wpflen+1) = -2.*(temp1/(2.*dr) + Qe_mn(j)*y(j)/dr**2.)/(3.*Vp_rho(j)*y(j))
                  pd(j+wpflen,j+wpflen-1) = -2.*(-1.*temp1/(2.*dr) + Qe_mn(j)*y(j)/dr**2.)/(3.*Vp_rho(j)*y(j))
                  ! ti
                  CALL fin_deriv2(Qi_mn(j-1:j+1)*y(j-1:j+1)/zeff_rho(j-1:j+1),dr,'m',temp1)
                  pd(j+2*wpflen,j+2*wpflen+1) = -2.*zeff_rho(j)*(temp1/(2.*dr) + Qe_mn(j)*y(j)/(zeff_rho(j)*dr**2.))/(3.*Vp_rho(j)*y(j))
                  pd(j+2*wpflen,j+2*wpflen-1) = -2.*zeff_rho(j)*(-1.*temp1/(2.*dr)+ Qe_mn(j)*y(j)/(zeff_rho(j)*dr**2.))/(3.*Vp_rho(j)*y(j))
            END DO

            ! cross terms
            CALL arr_fderiv2(wpflen,y(wpflen+1:wpflen*2),dr,temp4)
            CALL arr_fderiv2(wpflen,y(wpflen*2+1:wpflen*3),dr,temp5)
            DO j=1,wpflen
                  ! te/n
                  pd(j+wpflen,j) = 2.*((picol(j) - perf_rho(j))/y(j) - DpcDn(j)&
                                    + Qe_mn(j)*temp2(j)*temp4(j)/(Vp_rho(j)*y(j)))/(3.*y(j))&
                                    + y(j+wpflen)*(fn(j)/y(j) - pd(j,j))/y(j)
                  ! ti/n
                  pd(j+2*wpflen,j) = 2.*zeff_rho(j)*(DpcDn(j) - picol(j)/y(j)&
                                    + Qi_mn(j)*temp2(j)*temp5(j)/(Vp_rho(j)*y(j)*zeff_rho(j)))/(3.*y(j))&
                                    + y(j+2*wpflen)*(fn(j)/y(j) - pd(j,j))/y(j)
                  ! te/ti
                  pd(j+wpflen,j+2*wpflen) = -2.*DpcDti(j)/(3.*y(j))
                  ! ti/te
                  pd(j+2*wpflen,j+wpflen) = 2.*zeff_rho(j)*DpcDte(j)/(3.*y(j))
            END DO
      END SUBROUTINE bal_update_jac

      SUBROUTINE calc_pi_coll2(numpts,n_e,t_e,t_i,z_eff,m_i,p_i)
            USE stel_kinds, ONLY: rprec
            USE stelltran_vars, ONLY: ec, me
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: numpts
            REAL(rprec), INTENT(IN) :: m_i
            REAL(rprec), DIMENSION(numpts), INTENT(IN) :: n_e, t_e, t_i, z_eff
            REAL(rprec), DIMENSION(numpts), INTENT(OUT) :: p_i
            REAL(rprec), DIMENSION(numpts) :: lambda, tau_e
            INTEGER :: j
            ! NOTE: Formulary eqns are in cgs so need to do converting
            ! Determine appropriate Coulomb Logarithm, NRL Formulary p. 34
            DO j=1,numpts
                  IF (t_e(j) < 10.*z_eff(j)**2.) THEN
                        lambda(j) = 23. - LOG(SQRT(n_e(j)/1.D6)*z_eff(j)/t_e(j)**(1.5))
                  ELSE
                        lambda(j) = 24. - LOG(SQRT(n_e(j)/1.D6)/t_e(j))
                  END IF
            END DO
            ! Find electron -> ion energy from collisions. NRL Formulary p. 37
            tau_e = 3.44D5*t_e(:)**1.5/(lambda(:)*(n_e(:)/1.D6))
            p_i = 3.0*me*n_e(:)*ec*(t_e(:)-t_i(:))/(m_i*tau_e(:)) ! Should be good to not use cgs here since tau is in sec, ec is for eV -> J conversion
      END SUBROUTINE calc_pi_coll2

      SUBROUTINE save_ST_param(fname,npts,indep,dep)
            CHARACTER(200), INTENT(IN) :: fname
            INTEGER, INTENT(IN) :: npts
            REAL, DIMENSION(npts), INTENT(IN) :: indep, dep
            INTEGER :: j
            IF (access(fname,' ') == 0) THEN
                  OPEN(UNIT=12, FILE=fname, ACTION="write", STATUS="old",POSITION="append")
                  do j=1,npts
                        write(12,*) dep(j), indep(j)
                  end do
                  CLOSE(UNIT=12)
            ELSE
                  OPEN(UNIT=12, FILE=fname, ACTION="write", STATUS="new")
                  write(12,'(I2.2)') prof_length
                  do j=1,npts
                        write(12,*) dep(j), indep(j)
                  end do
                  CLOSE(UNIT=12)
            END IF
      END SUBROUTINE save_ST_param















