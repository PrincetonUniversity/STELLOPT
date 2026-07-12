!-----------------------------------------------------------------------
!     Function:      out_beams3d_part
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          10/21/2021
!     Description:   Save output from full orbit pushing.
!-----------------------------------------------------------------------
SUBROUTINE out_beams3d_part(t, q)
    !-----------------------------------------------------------------------
    !     Libraries
    !-----------------------------------------------------------------------
    USE stel_kinds, ONLY: rprec
    USE beams3d_runtime, ONLY: dt, lverb, pi2, lneut, t_end, lvessel, &
                               lhitonly, npoinc, lcollision, ldepo, &
                               weight, invpi2, ndt, ndt_max, lfidasim, lfidasim_cyl
    USE beams3d_lines, ONLY: R_lines, Z_lines, PHI_lines, myline, moment, &
                             nparticles, moment_lines, myend, &
                             vr_lines, vphi_lines, vz_lines, &
                             vll_lines, neut_lines, mytdex, next_t,&
                             xlast, ylast, zlast, dense_prof, &
                             ltherm, S_lines, U_lines, B_lines, &
                             ndot_prof, partvmax, &
                             ns_prof1, ns_prof2, ns_prof3, ns_prof4, &
                             ns_prof5, mymass, mycharge, mybeam, end_state, &
                             dist5d_prof, dist5d_fida, win_dist5d, nsh_prof4, &
                             h2_prof, h3_prof, h4_prof, h5_prof, my_end, &
                             r_h, p_h, z_h, e_h, pi_h, E_by_v, h1_prof, &
                             previous_time, wall_hit_valid, wall_hit_field_valid, wall_hit_model, wall_hit_face, &
                             wall_hit_fraction, wall_hit_time, wall_hit_r, &
                             wall_hit_phi, wall_hit_z, wall_hit_vll, &
                             wall_hit_moment, wall_hit_b, wall_hit_s, &
                             wall_hit_u, wall_hit_vr, wall_hit_vphi, &
                             wall_hit_vz, wall_hit_energy, time_lines, &
                             previous_vll, previous_moment, previous_vr, &
                             previous_vphi, previous_vz
    USE beams3d_grid
    USE beams3d_physics_mod, ONLY: beams3d_physics_fo, beams3d_event_diagnostics
    USE wall_mod, ONLY: collide, get_wall_ik, get_wall_area
    USE mpi_params
    USE mpi_inc
    !-----------------------------------------------------------------------
    !     Input Parameters
    !          t          Location along fieldline in t
    !          q          (R,phi,Z,vR,vphi,vz)
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(inout) :: t
    DOUBLE PRECISION, INTENT(inout) :: q(6)
    !-----------------------------------------------------------------------
    !     Local Variables
    !     jint      Index along phi
    !-----------------------------------------------------------------------
    LOGICAL             :: lhit, lfield, loutside_grid
    INTEGER             :: ier, d1, d2, d3, d4, d5
    DOUBLE PRECISION         :: x0,y0,z0,x1,y1,z1,xw,yw,zw,vperp, &
                                br_temp, bphi_temp, bz_temp, &
                                v_total, binv, vll_temp, &
                                rho0, callback_time, hit_fraction, &
                                event_br, event_bphi, event_bz
    DOUBLE PRECISION    :: q2(6),qdot(6), q4(4)
    ! For splines
    INTEGER :: i,j,k,l
    REAL*8 :: xparam, yparam, zparam !, hx, hy, hz, hxi, hyi, hzi
    REAL*8 :: fval(1), fval2(1)
    INTEGER, parameter :: ict(8)=(/1,0,0,0,0,0,0,0/)
    REAL*8, PARAMETER :: one = 1
    !-----------------------------------------------------------------------
    !     Begin Function
    !-----------------------------------------------------------------------
    callback_time = t
    lhit = .false.
    time_lines(mytdex, myline)   = callback_time
    R_lines(mytdex, myline)      = q(1)
    PHI_lines(mytdex, myline)    = q(2)
    Z_lines(mytdex, myline)      = q(3)
    vr_lines(mytdex, myline)     = q(4)
    vphi_lines(mytdex, myline)   = q(5)
    vz_lines(mytdex, myline)     = q(6)
    S_lines(mytdex, myline)      = 1.5
    U_lines(mytdex, myline)      = 0.0
    B_lines(mytdex, myline)      = -1.0
    neut_lines(mytdex,myline)    = lneut
    x0 = MOD(q(2), phimax)
    IF (x0 < 0) x0 = x0 + phimax
    loutside_grid = .not.((q(1) >= rmin-eps1) .and. (q(1) <= rmax+eps1) .and. &
                          (x0 >= phimin-eps2) .and. (x0 <= phimax+eps2) .and. &
                          (q(3) >= zmin-eps3) .and. (q(3) <= zmax+eps3))
    rho_help = 0
    IF (.not.loutside_grid) THEN
       i = MIN(MAX(COUNT(raxis < q(1)),1),nr-1)
       j = MIN(MAX(COUNT(phiaxis < x0),1),nphi-1)
       k = MIN(MAX(COUNT(zaxis < q(3)),1),nz-1)
       xparam = (q(1) - raxis(i)) * hri(i)
       yparam = (x0 - phiaxis(j)) * hpi(j)
       zparam = (q(3) - zaxis(k)) * hzi(k)
       CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                       hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                       XRHO4D(1,1,1,1),nr,nphi,nz)
       
       CALL R8HERM3FCN(ict,1,1,fval2,i,j,k,xparam,yparam,zparam,&
                       hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                       YRHO4D(1,1,1,1),nr,nphi,nz)
       y0 = fval(1)*fval(1)+fval2(1)*fval2(1) ! s = rho*rho
       rho_help = SQRT(y0)
       z0 = ATAN2(fval2(1),fval(1))
       S_lines(mytdex, myline) = y0 
       U_lines(mytdex, myline) = z0
       ! Now We get Br, Bphi, Bz to calc B (and vll)
       CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                       hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                       BR4D(1,1,1,1),nr,nphi,nz)
       br_temp = fval(1)
       CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                       hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                       BPHI4D(1,1,1,1),nr,nphi,nz)
       bphi_temp = fval(1)
       CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                       hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                       BZ4D(1,1,1,1),nr,nphi,nz)
       bz_temp = fval(1)

       fval(1) = SQRT( br_temp*br_temp + bphi_temp*bphi_temp + bz_temp*bz_temp )
       B_lines(mytdex, myline) = fval(1)
       binv = one/fval(1)

       ! Calc Vll
       vll_temp = ( br_temp*q(4) +  bphi_temp*q(5) + bz_temp*q(6) ) * binv
       vll_lines(mytdex,myline) = vll_temp


       ! Calculate the moment
       v_total = SUM(q(4:6)*q(4:6)) !Vtotal^2
       vperp   = v_total - vll_temp*vll_temp ! Vperp^2
       moment  = 0.5*mymass*vperp*binv
       moment_lines(mytdex,myline) = moment

       ! Calc dist func bins
       x0    = MOD(q(2),pi2)
       IF (x0 < 0) x0 = x0 + pi2
       d1 = MAX(MIN(CEILING(rho_help*h1_prof     ), ns_prof1), 1) ! Rho Bin
       d2 = MAX(MIN(CEILING( z0*h2_prof           ), ns_prof2), 1) ! U Bin
       d3 = MAX(MIN(CEILING( x0*h3_prof           ), ns_prof3), 1) ! V Bin
       d4 = MAX(MIN(1+nsh_prof4+FLOOR(h4_prof*vll_temp), ns_prof4), 1) ! vll
       d5 = MAX(MIN(CEILING(SQRT(vperp)*h5_prof         ), ns_prof5), 1) ! Vperp
       xw = weight(myline)*dt
       !CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,myworkid,0,win_dist5d,ier)
       dist5d_prof(mybeam,d1,d2,d3,d4,d5) = dist5d_prof(mybeam,d1,d2,d3,d4,d5) + xw
       !CALL MPI_WIN_UNLOCK(myworkid,win_dist5d,ier)
       IF (lfidasim_cyl) THEN
         d1 = MIN(MAX(CEILING((q(1)-rmin_fida)*r_h),1),nr_fida)
         d2 = MIN(MAX(CEILING((x0-phimin_fida)*p_h),1),nphi_fida)
         d3 = MIN(MAX(CEILING((q(3)-zmin_fida)*z_h),1),nz_fida)
         y0 = (q(4)**2+vperp**2)
         d4 = MIN(MAX(CEILING((y0*E_by_v-emin_fida)*e_h),1),nenergy_fida)
         d5 = MIN(MAX(CEILING((q(4)/SQRT(y0)-pimin_fida)*pi_h),1),npitch_fida)
         dist5d_fida(d1,d3,d2,d4,d5) = dist5d_fida(d1,d3,d2,d4,d5) + xw
       END IF
       IF (lcollision) CALL beams3d_physics_fo(t,q)
       IF (ltherm) THEN
          ndot_prof(mybeam,d1)   =   ndot_prof(mybeam,d1) + weight(myline)
          end_state(myline) = 1
          t = my_end
       END IF
    ELSE
       IF (lneut) end_state(myline)=3
       v_total = SUM(q(4:6)*q(4:6))
       vll_temp = SQRT(v_total) ! Vll=Vtotal
    END IF
    IF (lvessel .and. mytdex > 0 .and. wall_hit_valid(myline) == 0 .and. &
        (rho_help > 0.5 .or. loutside_grid)) THEN
       x0    = xlast
       y0    = ylast
       z0    = zlast
       x1    = q(1)*cos(q(2))
       y1    = q(1)*sin(q(2))
       z1    = q(3)
       CALL collide(x0,y0,z0,x1,y1,z1,xw,yw,zw,lhit,hit_fraction)
       IF (lhit) THEN
          q2(1) = SQRT(xw*xw+yw*yw)
          q2(2) = atan2(yw,xw)
          q2(3) = zw
          R_lines(mytdex,myline)       = q2(1)
          PHI_lines(mytdex,myline)     = q2(2)
          Z_lines(mytdex,myline)       = zw
          t = my_end+dt
          l = get_wall_ik()
          IF (wall_hit_valid(myline) == 0) THEN
             wall_hit_valid(myline) = 1
             wall_hit_model(myline) = 2
             wall_hit_face(myline) = l
             wall_hit_fraction(myline) = hit_fraction
             wall_hit_time(myline) = previous_time + &
                                     hit_fraction*(callback_time-previous_time)
             wall_hit_r(myline) = q2(1)
             wall_hit_phi(myline) = q2(2)
             wall_hit_z(myline) = q2(3)
             wall_hit_vr(myline) = previous_vr + hit_fraction* &
                                     (vr_lines(mytdex,myline)-previous_vr)
             wall_hit_vphi(myline) = previous_vphi + hit_fraction* &
                                       (vphi_lines(mytdex,myline)-previous_vphi)
             wall_hit_vz(myline) = previous_vz + hit_fraction* &
                                     (vz_lines(mytdex,myline)-previous_vz)
             CALL beams3d_event_diagnostics(q2(1),q2(2),q2(3), &
                    wall_hit_s(myline),wall_hit_u(myline),wall_hit_b(myline), &
                    event_br,event_bphi,event_bz,lfield)
             IF (lfield) THEN
                wall_hit_field_valid(myline) = 1
                wall_hit_vll(myline) = (event_br*wall_hit_vr(myline) + &
                                         event_bphi*wall_hit_vphi(myline) + &
                                         event_bz*wall_hit_vz(myline))/wall_hit_b(myline)
                v_total = wall_hit_vr(myline)**2 + wall_hit_vphi(myline)**2 + &
                          wall_hit_vz(myline)**2
                vperp = MAX(v_total-wall_hit_vll(myline)**2,0.0D0)
                wall_hit_moment(myline) = 0.5*mymass*vperp/wall_hit_b(myline)
             END IF
             wall_hit_energy(myline) = 0.5*mymass*(wall_hit_vr(myline)**2 + &
                                        wall_hit_vphi(myline)**2 + wall_hit_vz(myline)**2)
             time_lines(mytdex,myline) = wall_hit_time(myline)
             IF (lfield) THEN
                vll_lines(mytdex,myline) = wall_hit_vll(myline)
                moment_lines(mytdex,myline) = wall_hit_moment(myline)
                B_lines(mytdex,myline) = wall_hit_b(myline)
                S_lines(mytdex,myline) = wall_hit_s(myline)
                U_lines(mytdex,myline) = wall_hit_u(myline)
             END IF
             vr_lines(mytdex,myline) = wall_hit_vr(myline)
             vphi_lines(mytdex,myline) = wall_hit_vphi(myline)
             vz_lines(mytdex,myline) = wall_hit_vz(myline)
          END IF
          IF (lneut) THEN
             wall_shine(mybeam,l) = wall_shine(mybeam,l) + &
                                    weight(myline)*wall_hit_energy(myline)/get_wall_area(l)
          ELSE
             end_state(myline) = 2
             wall_load(mybeam,l) = wall_load(mybeam,l) + &
                                    weight(myline)*wall_hit_energy(myline)/get_wall_area(l)
          END IF
          IF (lhitonly) THEN
             R_lines(0,myline)      = SQRT(xlast*xlast+ylast*ylast)
             PHI_lines(0,myline)    = ATAN2(ylast,xlast)
             Z_lines(0,myline)      = zlast
             vll_lines(0,myline)    = previous_vll
             moment_lines(0,myline) = previous_moment
             vr_lines(0,myline)     = previous_vr
             vphi_lines(0,myline)   = previous_vphi
             vz_lines(0,myline)     = previous_vz
             time_lines(0,myline)   = previous_time
             neut_lines(0,myline)   = lneut
             R_lines(2,myline)      = q(1)
             PHI_lines(2,myline)    = q(2)
             Z_lines(2,myline)      = q(3)
             vll_lines(2,myline)    = vll_temp
             moment_lines(2,myline) = moment
             vr_lines(2,myline)     = q(4)
             vphi_lines(2,myline)   = q(5)
             vz_lines(2,myline)     = q(6)
             time_lines(2,myline)   = callback_time
             neut_lines(2,myline)   = lneut
          END IF
       ELSE
          xlast = x1
          ylast = y1
          zlast = z1
       END IF
    ELSE
       xlast = q(1)*cos(q(2))
       ylast = q(1)*sin(q(2))
       zlast = q(3)
    END IF
    IF (.not. lhit) THEN
       previous_time = callback_time
       IF (B_lines(mytdex,myline) > 0.0) THEN
          previous_vll = (br_temp*q(4) + bphi_temp*q(5) + bz_temp*q(6))*binv
          v_total = SUM(q(4:6)*q(4:6))
          vperp = MAX(v_total-previous_vll*previous_vll,0.0D0)
          previous_moment = 0.5*mymass*vperp*binv
       ELSE
          previous_vll = vll_lines(mytdex,myline)
          previous_moment = moment_lines(mytdex,myline)
       END IF
       previous_vr = q(4)
       previous_vphi = q(5)
       previous_vz = q(6)
    END IF
    ndt = ndt + 1
    IF (ndt .ge. ndt_max) THEN ! ge needed if npoinc = ndt
       mytdex = mytdex + 1
       ndt = 1
    END IF
    IF (lhitonly) mytdex = 1
    IF (mytdex > npoinc) t = my_end
    t = t + dt

    RETURN
    !-----------------------------------------------------------------------
    !     End Function
    !-----------------------------------------------------------------------
END SUBROUTINE out_beams3d_part
