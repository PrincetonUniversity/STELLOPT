!-----------------------------------------------------------------------
!     Module:        beams3d_physics_mod
!     Authors:       M. McMillan (matthew.mcmillan@my.wheaton.edu)
!     Date:          07/20/2012
!     Description:   This module contains the physics subroutines for
!                    manipulating particles and finding when to do so.
!-----------------------------------------------------------------------
MODULE beams3d_physics_mod
    !-----------------------------------------------------------------------
    !     Libraries
    !-----------------------------------------------------------------------
    USE stel_kinds, ONLY: rprec
    USE beams3d_runtime, ONLY: lneut, pi, pi2, dt, lverb, ADAS_ERR, dt_save

    USE beams3d_lines, ONLY: R_lines, Z_lines, PHI_lines, &
                             myline, mytdex, moment, ltherm, &
                             nsteps, nparticles, vll_lines, moment_lines, mybeam, mycharge, myZ, &
                             mymass, myv_neut, B_temp, rand_prob, cum_prob, tau, PE_lines, PI_lines

    USE beams3d_grid, ONLY: BR_spl, BZ_spl, delta_t, BPHI_spl, MODB_spl, &
                            phimax, TE_spl, NE_spl, TI_spl

    USE EZspline_obj
    USE EZspline
!DEC$ IF DEFINED (NTCC)
    USE adas_mod_simpl
    USE fpreact_calls
    USE periodic_table_mod
!DEC$ ENDIF  
      USE mpi_params 
    !-----------------------------------------------------------------------
    !     Module Variables
    !          lverb         Logical to control screen output
    !          vessel_string Limiting Surface Filename
    !          extcur        External currents for MGRID calculation
    !----------------------------------------------------------------------


CONTAINS
    !     Subroutines: _beams3d_physics: makes decisions about when to do what
    !                  _ionize, _neutralize: changes particle position, velocities,
    !                    and lneut logical to change from following a neutral path
    !                    to an ion path.
    !                  _t_flight: calculate the expectation value of the flight time
    !                    for the particle, in whatever state it is in, before switching
    !                    states (e.g. ion to neutral).

!-----------------------------------------------------------------------
!     Function:      beams3d_physics
!     Authors:       M. McMillan (matthew.mcmillan@my.wheaton.edu)
!     Date:          07/13/2012
!     Description:   Decide when to contribute various physical effects
!-----------------------------------------------------------------------
SUBROUTINE beams3d_physics(t, q)
    !-----------------------------------------------------------------------
    !     Input Parameters
    !          t          Location along fieldline in t
    !          q            (q(1),q(2),q(3),q(4)) = (R,phi,Z,vll)
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(inout) :: t
    DOUBLE PRECISION, INTENT(inout) :: q(4)
    !-----------------------------------------------------------------------
    !     Local Variables
    !     jint      Index along phi
    !     B_temp_cyl      B in cylindrical coords; B_temp is cartesian.
    !-----------------------------------------------------------------------
    INTEGER        :: ier
    DOUBLE PRECISION    :: r_temp, phi_temp, z_temp, vll, te_temp, ne_temp, ti_temp, speed, newspeed, &
                          zeta, sigma, zeta_mean, zeta_o, v_s, tau_inv, tau_spit_inv, &
                          reduction, dve,dvi, tau_spit, v_crit, coulomb_log, te_cube, &
                          inv_mymass, speed_cube, vcrit_cube, vfrac, modb
    DOUBLE PRECISION :: Ebench  ! for ASCOT Benchmark
    REAL*8         :: energy(1), temp(1), sigvii(1), sigvcx(1), sigvcxn(1), tempA(1), sigvei(1)

    ! Helpers
    DOUBLE PRECISION, PARAMETER :: half = 0.5D0
    DOUBLE PRECISION, PARAMETER :: electron_mass = 9.10938356D-31
    DOUBLE PRECISION, PARAMETER :: e_charge      = 1.60217662E-19
!    DOUBLE PRECISION, PARAMETER :: e_charge_4    = 6.5893345788E-76
!    DOUBLE PRECISION, PARAMETER :: sqrt_emass    = 9.5443090688E-16
    DOUBLE PRECISION, PARAMETER :: sqrt_pi       = 1.7724538509
!    DOUBLE PRECISION, PARAMETER :: inv_emass     = 1.0977691228E+30
      ! For splines
      INTEGER :: i,j,k
      REAL*8 :: xparam, yparam, zparam, hx, hy, hz, hxi, hyi, hzi
      INTEGER, parameter :: ict(8)=(/1,0,0,0,0,0,0,0/)
      REAL*8 :: fval(1)
    !-----------------------------------------------------------------------
    !     Begin Function
    !-----------------------------------------------------------------------

    ! Some intitial declarations. First 5 are used, remainder are (usually) reassigned below.

    ier      = 0
    r_temp   = q(1)
    phi_temp = MODULO(q(2), phimax)
    IF (phi_temp < 0) phi_temp = phi_temp + phimax
    z_temp   = q(3)
    vll      = q(4)

    te_temp  = 0.0
    ne_temp  = 0.0
    ti_temp  = 0.0  ! KeV (target temperature/nucleon)
    speed = 0.0
    reduction = 0.0

    tau_spit = 1.0E20  ! The Spitzer characteristic time. This and v_c reset below.
    tau_spit_inv = 0.0 ! So we can set dve dvi to zero
    v_crit   = 0.0     ! Critical velocity associated with the critical energy; transition to ion slowing.
    coulomb_log = 15   ! Usually reset below depending on Te etc.

    tau_inv = 10.0

    ier = 0       ! Set B_temp from splines.
    CALL EZspline_isInDomain(BR_spl,r_temp,phi_temp,z_temp,ier)
    IF (ier == 0) THEN
       ! Get the gridpoint info (this is possible since all grids are the same)
       CALL R8HERM3xyz(r_temp,phi_temp,z_temp,&
                            MODB_spl%x1(1),MODB_spl%n1,&
                            MODB_spl%x2(1),MODB_spl%n2,&
                            MODB_spl%x3(1),MODB_spl%n3,&
                            MODB_spl%ilin1,MODB_spl%ilin2,MODB_spl%ilin3,&
                            i,j,k,xparam,yparam,zparam,&
                            hx,hxi,hy,hyi,hz,hzi,ier)
       ! Evaluate the Splines
       CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                       hx,hxi,hy,hyi,hz,hzi,&
                       MODB_spl%fspl(1,1,1,1),MODB_spl%n1,MODB_spl%n2,MODB_spl%n3)
       modb = fval(1)
       CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                       hx,hxi,hy,hyi,hz,hzi,&
                       TE_spl%fspl(1,1,1,1),TE_spl%n1,TE_spl%n2,TE_spl%n3)
       te_temp = fval(1)
       CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                       hx,hxi,hy,hyi,hz,hzi,&
                       NE_spl%fspl(1,1,1,1),NE_spl%n1,NE_spl%n2,NE_spl%n3)
       ne_temp = fval(1)
       CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                       hx,hxi,hy,hyi,hz,hzi,&
                       TI_spl%fspl(1,1,1,1),TI_spl%n1,TI_spl%n2,TI_spl%n3)
       ti_temp = fval(1)
       !IF (myworkid == master) PRINT *,modb,te_temp,ne_temp,ti_temp

       ! Adjust to eV from keV
       !te_temp = te_temp !* 1000
       !ti_temp = ti_temp * 1000 ! ti needed in keV for Adas Ionization
       !ne_temp = ne_temp *1E-6 !ne in cm^-3

       ! Helpers
       te_cube = te_temp * te_temp *te_temp
       inv_mymass = 1/mymass

       !---------------------------------------------------------------
       !  Viscouse Velocity Reduction
       !---------------------------------------------------------------
       IF ((te_temp > 0).and.(ne_temp > 0)) THEN
          ! Electron-ion Coulomb Logaritm (NRL pg. 35), ignore cold electrons
          ! te in eV and ne in cm^-3
          IF (te_temp < 10*myZ*myZ) THEN
             coulomb_log = 23 - log( myZ*sqrt(ne_temp*1E-6/(te_cube) )   )
          ELSE
             coulomb_log = 24 - log( sqrt(ne_temp*1E-6)/(te_temp) )
          END IF
          IF (coulomb_log .le. 1) coulomb_log = 1
          v_crit = (( 0.75*sqrt_pi*electron_mass*inv_mymass )**0.33333333333 )*sqrt(te_temp)*5.93096892024D5
          vcrit_cube = v_crit*v_crit*v_crit
          tau_spit = 3.777183D41*mymass*SQRT(te_cube)/(ne_temp*myZ*myZ*coulomb_log)  ! note ne should be in m^-3 here
          tau_spit_inv = (1.0D0)/tau_spit
       ELSE
          coulomb_log = 15
          vcrit_cube = 0
          tau_spit_inv = 0
       END IF

!       IF (.not. lneut) THEN  
          v_s = 1.5*sqrt(e_charge*ti_temp*inv_mymass)
          speed = sqrt( vll*vll + 2*moment*modb*inv_mymass )
          !dve   = speed/tau_spit
          !dvi   = vcrit_cube/(tau_spit*speed*speed)
          dve   = speed*tau_spit_inv
          dvi   = vcrit_cube*tau_spit_inv/(speed*speed)
          reduction = dve + dvi
          newspeed = speed - reduction*dt
          !Ebench = half*mymass*newspeed*newspeed/e_charge
          !IF ((Ebench <= 1000.) .or. (Ebench <= 1.5*te_temp)) THEN !Benchmark Thermalize
          IF (newspeed < v_s) THEN  ! Thermalize
             dve       = dve/reduction
             dvi       = dvi/reduction
             reduction = (speed - v_s)/dt  ! Thermalize
             dve       = dve*reduction
             dvi       = dvi*reduction
             newspeed = speed - reduction*dt
             ltherm = .true.
             PE_lines(mytdex,myline) = PE_lines(mytdex,myline)+half*mymass*dve*dve*dt
             PI_lines(mytdex,myline) = PI_lines(mytdex,myline)+half*mymass*dve*dve*dt
             vll = (newspeed/speed)*vll
             moment = newspeed*newspeed*moment/(speed*speed)
             q(4) = vll
             RETURN
          END IF
          PE_lines(mytdex,myline) = PE_lines(mytdex,myline)+half*mymass*dve*dve*dt
          PI_lines(mytdex,myline) = PI_lines(mytdex,myline)+half*mymass*dve*dve*dt
          vll = (newspeed/speed)*vll
          moment = newspeed*newspeed*moment/(speed*speed)
          speed = newspeed

       !---------------------------------------------------------------
       !  Pitch Angle Scattering
       !---------------------------------------------------------------
          !v_s = half*vcrit_cube/tau_spit
          v_s = half*vcrit_cube*tau_spit_inv
          speed_cube = speed*speed*speed
          zeta_o = vll/speed   ! Record the current pitch.
          CALL gauss_rand(1,zeta)  ! A random from a standard normal (1,1)
          sigma = sqrt( 2*v_s*(1.0D0-zeta_o*zeta_o)*dt/speed_cube ) ! The standard deviation.
          zeta_mean = zeta_o !*(1.0 - (2*v_s*dt)/(speed*speed*speed) )  ! The new mean in the distribution.
          !zeta_mean = zeta_o *(1.0D0 - (2*v_s*dt)/(speed_cube) )  ! The new mean in the distribution.
          zeta = zeta*sigma + zeta_mean  ! The new pitch angle.
!         The pitch angle MUST NOT go outside [-1,1] nor be NaN; but could happen accidentally with the distribution.
          IF (ABS(zeta) >  0.999999999D+00) zeta =  SIGN(0.999999999D+00,zeta)
!         Comment out the next 4 lines to turn off pitch angle scattering.
          vll = zeta*speed
          moment = half*mymass*(speed*speed - vll*vll)/modb
!       END IF

       ! Now update parallel velocity
       q(4) = vll

    END IF

    RETURN
    !-----------------------------------------------------------------------
    !     End Function
    !-----------------------------------------------------------------------
END SUBROUTINE beams3d_physics

SUBROUTINE beams3d_follow_neut(t, q)
    USE beams3d_grid, ONLY: TI_spl, NE_spl, S_spl, rmin, rmax, phimax
    USE beams3d_lines, ONLY: myline,xlast,ylast,zlast
    USE beams3d_runtime, ONLY: t_end, lvessel,to3
    USE wall_mod, ONLY: collide
    !-----------------------------------------------------------------------
    !     Input Parameters
    !          t          Location along fieldline in t
    !          q            (q(1),q(2),q(3),q(4)) = (R,phi,Z,vll)
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(inout) :: t
    DOUBLE PRECISION, INTENT(inout) :: q(4)
    LOGICAL          :: ltest
    INTEGER          :: ier
    DOUBLE PRECISION :: rinv, phi_temp, dt_local, ti_temp, ne_temp, tau_inv, &
                        s_temp, x0, y0, z0, xw, yw, zw, te_temp
    DOUBLE PRECISION :: qf(3)
    REAL*8           :: energy(1), temp(1), sigvii(1), sigvcx(1), sigvcxn(1), tempA(1), sigvei(1)
    DOUBLE PRECISION, PARAMETER :: dl = 1D-3
    DOUBLE PRECISION, PARAMETER :: one = 1
    DOUBLE PRECISION, PARAMETER :: half = 0.5D0
    ! For splines
    INTEGER :: i,j,k
    REAL*8 :: xparam, yparam, zparam, hx, hy, hz, hxi, hyi, hzi
    INTEGER, parameter :: ict(8)=(/1,0,0,0,0,0,0,0/)
    REAL*8 :: fval(1)
!

    dt_local = 10*dl/q(4)  ! Timestep (10 times larger than in plasma)
    energy(1) = half*mymass*q(4)*q(4)/1000  ! Vll = V_neut (doesn't change durring integration) needed in keV
    qf(1) = q(1)*cos(q(2))
    qf(2) = q(1)*sin(q(2))
    qf(3) = q(3)
    ! First follow into plasma
    DO
       qf = qf + myv_neut*dt_local
       q(1) = sqrt(qf(1)*qf(1)+qf(2)*qf(2))
       q(2) = ATAN2(qf(2),qf(1))
       q(3) = qf(3)
       t = t + dt_local
       phi_temp = MODULO(q(2), phimax)
       IF (phi_temp < 0) phi_temp = phi_temp + phimax
       CALL EZspline_isInDomain(S_spl,q(1),phi_temp,q(3),ier)
       IF (ier==0) THEN
          s_temp =1.5
          CALL EZspline_interp(S_spl,q(1),phi_temp,q(3),s_temp,ier)
          IF (s_temp < one) EXIT
       END IF
       IF ((q(1) > 5*rmax)  .or. (q(1) < rmin)) THEN; t = t_end(myline)+dt_local; RETURN; END IF  ! We're outside the grid
    END DO
    ! Take a step back
    qf = qf - myv_neut*dt_local
    t  =  t - dt_local
    ! Setup for high resolution stepping
    dt_local = dl/q(4)  ! Set high res timestep
    CALL RANDOM_NUMBER(rand_prob)
    cum_prob = one
    ltest = .false.
    ! Follow into plasma
    DO
       qf = qf + myv_neut*dt_local
       CALL FLUSH(6)
       q(1) = sqrt(qf(1)*qf(1)+qf(2)*qf(2))
       q(2) = ATAN2(qf(2),qf(1))
       q(3) = qf(3)
       t = t + dt_local
       phi_temp = MODULO(q(2), phimax)
       IF (phi_temp < 0) phi_temp = phi_temp + phimax
       IF ((q(1) > 5*rmax)  .or. (q(1) < rmin)) THEN; t = t_end(myline)+dt_local; RETURN; END IF  ! We're outside the grid
       CALL R8HERM3xyz(q(1),phi_temp,q(3),&
                            TI_spl%x1(1),TI_spl%n1,&
                            TI_spl%x2(1),TI_spl%n2,&
                            TI_spl%x3(1),TI_spl%n3,&
                            TI_spl%ilin1,TI_spl%ilin2,TI_spl%ilin3,&
                            i,j,k,xparam,yparam,zparam,&
                            hx,hxi,hy,hyi,hz,hzi,ier)
       CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                       hx,hxi,hy,hyi,hz,hzi,&
                       TI_spl%fspl(1,1,1,1),TI_spl%n1,TI_spl%n2,TI_spl%n3)
       ti_temp = fval(1)/1000
       CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                       hx,hxi,hy,hyi,hz,hzi,&
                       TE_spl%fspl(1,1,1,1),TE_spl%n1,TE_spl%n2,TE_spl%n3)
       te_temp = fval(1)/1000
       CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                       hx,hxi,hy,hyi,hz,hzi,&
                       NE_spl%fspl(1,1,1,1),NE_spl%n1,NE_spl%n2,NE_spl%n3)
       ne_temp = fval(1)
       CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                       hx,hxi,hy,hyi,hz,hzi,&
                       S_spl%fspl(1,1,1,1),S_spl%n1,S_spl%n2,S_spl%n3)
       s_temp = fval(1)
       IF ((s_temp > one) .and. ltest) EXIT
       IF ((ti_temp>0).and.(ne_temp>0)) THEN
          IF (s_temp <= one) ltest = .true.
          temp(1) = ti_temp! KeV
          ! This formula comes from the ADAS description of how to use the functions. (M. Gorelenkova)
          tempA(1) = te_temp + to3*0.000544602984424355*energy(1)
          ! Arguments to btsigv(irtype,beamchrg,eabeam,tatarg,n,izbeam,iztarg,btsigv, istat )hatom_btsigv.f90
          ! irtype 1:CX, 2:II
          ! beamchrg 1:neutral, 2:ion
          ! eabeam beam energy in keV/nucleon
          ! tatarg target energy in keV/nucleon
          ! n array size
          ! izbeam  Beam Z
          ! iztarg  Target Z
          ! btsigv  cross section
!DEC$ IF DEFINED (NTCC)
          CALL adas_btsigv(2,1,energy,temp,1,1,1,sigvii,ier)  ! Ion Impact ionization cross-section term.
          CALL adas_btsigv(1,1,energy,temp,1,1,1,sigvcx,ier)  ! Charge Exchange ionization cross-section term.
!DEC$ ENDIF
          ! Arguments to sigvte(zneut,tevec,n1,sigv_adas,istat)
          ! zneut charge (=1)
          ! tevec electron temperature [keV]
          ! n1 array size
!DEC$ IF DEFINED (NTCC)
          CALL adas_sigvte_ioniz(1,tempA,1,sigvei,ier)        ! Electron Impact ionization cross-section term.
!DEC$ ENDIF
          tau_inv = ((sigvii(1) + sigvcx(1) + sigvei(1))*ne_temp) ! Delete a term if desired. (save a comment)
          cum_prob = cum_prob*exp(-dt_local*tau_inv)
          IF (cum_prob < rand_prob) THEN
             lneut = .false.
             RETURN
          END IF
       END IF
    END DO
    ! Follow to wall
    qf = qf - myv_neut*dt_local
    t  =  t - dt_local
    x0 = qf(1); y0 = qf(2); z0 = qf(3)
    dt_local = 100*dl/q(4)  ! Timestep (100 times larger than in plasma)
    ltest = .FALSE.
    DO
       qf = qf + myv_neut*dt_local
       t = t + dt_local
       IF (lvessel) CALL collide(x0,y0,z0,qf(1),qf(2),qf(3),xw,yw,zw,ltest)
       !WRITE(328,*) x0,y0,z0,qf(1),qf(2),qf(3),xw,yw,zw,ltest
       IF (ltest) THEN
          q(1) = sqrt(xw*xw+yw*yw)
          q(2) = atan2(yw,xw)
          q(3) = zw
          ! Next lines are so that out_beams3d_nag detects the wall.
          xlast = x0; ylast=y0; zlast=z0
          q(1) = sqrt(qf(1)*qf(1)+qf(2)*qf(2))
          q(2) = ATAN2(qf(2),qf(1))
          q(3) = qf(3)
          RETURN
       END IF
       xlast = x0; ylast=y0; zlast=z0
       x0 = qf(1); y0 = qf(2); z0 = qf(3)
       q(1) = sqrt(qf(1)*qf(1)+qf(2)*qf(2))
       q(2) = ATAN2(qf(2),qf(1))
       q(3) = qf(3)
       IF ((q(1) > 5*rmax)  .or. (q(1) < rmin)) THEN; t = t_end(myline)+dt_local; RETURN; END IF  ! We're outside the grid
    END DO

    RETURN
END SUBROUTINE beams3d_follow_neut


SUBROUTINE beams3d_ionize(t, q)
    !-----------------------------------------------------------------------
    !     Input Parameters
    !          t          Location along trajectory in t
    !          q            (q(1),q(2),q(3),q(4)) = (R,phi,Z,vll)
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(inout) :: t
    DOUBLE PRECISION, INTENT(inout) :: q(4)
    !-----------------------------------------------------------------------
    !     Local Variables
    !     binv      1/mod(B_temp)
    !-----------------------------------------------------------------------
    INTEGER          :: npi2theta, ier
    DOUBLE PRECISION :: binv, rho(3), x, y, phi_temp, br_temp, bp_temp, bz_temp,&
                        bx_temp, by_temp, modb_temp, r_temp, z_temp
    REAL*8, PARAMETER :: one = 1

    !Begin subroutine
    npi2theta = floor(q(2)/pi2)

    ! Set some stuff before we do anything
    lneut = .false.
    CALL RANDOM_NUMBER(rand_prob)

    ! Evaluate the Splines
    ier = 0
    phi_temp = MOD(q(2),phimax)
    IF (phi_temp < 0) phi_temp = phi_temp + phimax
    r_temp = q(1)
    z_temp = q(3)
    CALL EZspline_interp(BR_spl,r_temp,phi_temp,z_temp,br_temp,ier)
    CALL EZspline_interp(BPHI_spl,r_temp,phi_temp,z_temp,bp_temp,ier)
    CALL EZspline_interp(BZ_spl,r_temp,phi_temp,z_temp,bz_temp,ier)
    CALL EZspline_interp(MODB_spl,r_temp,phi_temp,z_temp,modb_temp,ier)
    bx_temp = br_temp*cos(q(2))-bp_temp*sin(q(2))
    by_temp = br_temp*sin(q(2))+bp_temp*cos(q(2))
    binv = one/modb_temp

    ! Calculate the parallel velocity
    !q(4) = binv*SUM(myv_neut(1:3)*B_temp(1:3),DIM=1)
    q(4) = binv*( myv_neut(1)*bx_temp + myv_neut(2)*by_temp + myv_neut(3)*bz_temp )
    moment = 0.5*binv*mymass*(myv_neut(1)*myv_neut(1) + myv_neut(2)*myv_neut(2) + myv_neut(3)*myv_neut(3) - q(4)*q(4) )

    ! Check to see we didn't inject perfectly parallel (negative moment possible)
    IF (moment < 0) THEN
       moment = 1000*TINY(moment)
       RETURN
    END IF

    ! Calculate Gyroradius
    rho(1) = myv_neut(2)*bz_temp - myv_neut(3)*by_temp
    rho(2) = myv_neut(3)*bx_temp - myv_neut(1)*bz_temp
    rho(3) = myv_neut(1)*by_temp - myv_neut(2)*bx_temp
    rho    = (mymass*binv*binv/mycharge)*rho

    ! Move to Gyrocenter
    x = q(1)*cos(q(2)) + rho(1)
    y = q(1)*sin(q(2)) + rho(2)
    q(1) = sqrt(x*x + y*y)
    IF (atan2(y,x) <= 0.0) THEN
        q(2) = atan2(y,x)+(npi2theta+1)*pi2
    ELSE
        q(2) = atan2(y,x)+npi2theta*pi2
    END IF
    q(3) = q(3) + rho(3)

    RETURN

END SUBROUTINE beams3d_ionize


SUBROUTINE beams3d_neutralize(t, q)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(inout) :: t
    DOUBLE PRECISION, INTENT(inout) :: q(4)
    !-----------------------------------------------------------------------
    !     Local Variables
    !     binv      1/mod(B_temp)
    !-----------------------------------------------------------------------
    INTEGER          :: npi2theta
    DOUBLE PRECISION :: binv, bnz, e1(3), e2(3), theta, rho(3), x, y, vperp
    !Begin subroutine

    npi2theta = floor(q(2)/pi2)

    binv   = 1/B_temp(4)
    bnz    = sqrt(B_temp(4)*B_temp(4) - B_temp(3)*B_temp(3))
    vperp  = sqrt( 2*B_temp(4)*moment/mymass )

    IF (bnz == 0) THEN
        e1    = (/1,0,0/)
        e2    = (/0,1,0/)
    ELSE
        e1 = (/ -B_temp(3)*B_temp(1)*binv/bnz, -B_temp(3)*B_temp(2)*binv/bnz, binv*bnz /)
        e2 = binv*(/ B_temp(2)*e1(3)-B_temp(3)*e1(2), B_temp(3)*e1(1)-B_temp(1)*e1(3), B_temp(1)*e1(2)-B_temp(2)*e1(1) /)
    END IF

    CALL RANDOM_NUMBER(theta)
    theta = theta*pi2
    rho         = ( sin(theta)*e1 + cos(theta)*e2 )

    myv_neut(1) = q(4)*binv*B_temp(1) - vperp*binv*( B_temp(2)*rho(3)-B_temp(3)*rho(2) )
    myv_neut(2) = q(4)*binv*B_temp(2) - vperp*binv*( B_temp(3)*rho(1)-B_temp(1)*rho(3) )
    myv_neut(3) = q(4)*binv*B_temp(3) - vperp*binv*( B_temp(1)*rho(2)-B_temp(2)*rho(1) )

    rho         = ( mymass*vperp*binv/mycharge )*rho
    x           = q(1)*cos(q(2)) + rho(1)
    y           = q(1)*sin(q(2)) + rho(2)

    q(1)        = sqrt(x*x + y*y)
    IF (atan2(y,x) <= 0.0) THEN
        q(2) = atan2(y,x)+(npi2theta+1)*pi2
    ELSE
        q(2) = atan2(y,x)+npi2theta*pi2
    END IF
    q(3)        = q(3) + rho(3)

    moment = 0.5*binv*mymass*( myv_neut(1)*myv_neut(1) + myv_neut(2)*myv_neut(2) + myv_neut(3)*myv_neut(3) - q(4)*q(4) )

    lneut = .true.

    CALL RANDOM_NUMBER(rand_prob)

END SUBROUTINE beams3d_neutralize

END MODULE beams3d_physics_mod
