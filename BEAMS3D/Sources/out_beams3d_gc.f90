!-----------------------------------------------------------------------
!     Function:      out_beams3d_gc
!     Authors:       S. Lazerson (lazerson@pppl.gov), M. McMillan (matthew.mcmillan@my.wheaton.edu)
!     Date:          06/20/2012
!     Description:   Save output from field line following while running
!                    and updates progress bar.
!-----------------------------------------------------------------------
SUBROUTINE out_beams3d_gc(t, q)
    !-----------------------------------------------------------------------
    !     Libraries
    !-----------------------------------------------------------------------
    USE stel_kinds, ONLY: rprec
    USE beams3d_runtime, ONLY: dt, lverb, pi2, lneut, t_end, lvessel, &
                               lhitonly, npoinc, lcollision, ldepo, &
                               weight, invpi2, ndt, ndt_max, lfidasim, lfidasim_cyl
    USE beams3d_lines, ONLY: R_lines, Z_lines, PHI_lines, myline, moment, &
                             nparticles, moment_lines, myend, &
                             vll_lines, neut_lines, mytdex, next_t,&
                             xlast, ylast, zlast, dense_prof, &
                             ltherm, S_lines, U_lines, B_lines, &
                             ndot_prof, partvmax, &
                             ns_prof1, ns_prof2, ns_prof3, ns_prof4, &
                             ns_prof5, mymass, mycharge, mybeam, end_state, &
                             dist5d_prof, dist5d_fida, win_dist5d, nsh_prof4, &
                             h2_prof, h3_prof, h4_prof, h5_prof, my_end, &
                             r_h, p_h, z_h, e_h, pi_h, E_by_v, h1_prof
    USE beams3d_grid
    USE beams3d_physics_mod, ONLY: beams3d_physics_gc
    USE wall_mod, ONLY: collide, get_wall_ik, get_wall_area
    USE mpi_params
    USE mpi_inc
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
    !-----------------------------------------------------------------------
    LOGICAL             :: lhit
    INTEGER             :: ier, d1, d2, d3, d4, d5, d1f
    DOUBLE PRECISION         :: x0,y0,z0,x1,y1,z1,xw,yw,zw, vperp
    DOUBLE PRECISION    :: q2(4),qdot(4)
    ! For splines
    INTEGER :: i,j,k,l
    REAL*8 :: xparam, yparam, zparam !, hx, hy, hz, hxi, hyi, hzi
    REAL*8 :: fval(1), fval2(1)
    INTEGER, parameter :: ict(8)=(/1,0,0,0,0,0,0,0/)
    REAL*8, PARAMETER :: one = 1
    !-----------------------------------------------------------------------
    !     Begin Function
    !-----------------------------------------------------------------------
    R_lines(mytdex, myline)      = q(1)
    PHI_lines(mytdex, myline)    = q(2)
    Z_lines(mytdex, myline)      = q(3)
    vll_lines(mytdex, myline)    = q(4)
    moment_lines(mytdex, myline) = moment
    neut_lines(mytdex,myline)     = lneut
    x0 = MOD(q(2), phimax)
    IF (x0 < 0) x0 = x0 + phimax
    !CALL EZspline_isInDomain(S_spl,q(1),x0,q(3),ier)
    rho_help = 0  ! If we're out of domain then don't worry about collisions
    !IF (ier==0) THEN
    IF ((q(1) >= rmin-eps1) .and. (q(1) <= rmax+eps1) .and. &
        (x0 >= phimin-eps2) .and. (x0 <= phimax+eps2) .and. &
        (q(3) >= zmin-eps3) .and. (q(3) <= zmax+eps3)) THEN


       i = MIN(MAX(COUNT(raxis < q(1)),1),nr-1)
       j = MIN(MAX(COUNT(phiaxis < x0),1),nphi-1)
       k = MIN(MAX(COUNT(zaxis < q(3)),1),nz-1)
       xparam = (q(1) - raxis(i)) * hri(i)
       yparam = (x0 - phiaxis(j)) * hpi(j)
       zparam = (q(3) - zaxis(k)) * hzi(k)
       CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                       hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                       X4D(1,1,1,1),nr,nphi,nz)
       CALL R8HERM3FCN(ict,1,1,fval2,i,j,k,xparam,yparam,zparam,&
                       hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                       Y4D(1,1,1,1),nr,nphi,nz)
       y0 = SQRT(fval(1)*fval(1) + fval2(1) * fval2(1))
       rho_help = sqrt(y0)
       z0 = ATAN2(fval2(1),fval(1))
       S_lines(mytdex, myline) = y0 
       U_lines(mytdex, myline) = z0
       CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                       hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                       MODB4D(1,1,1,1),nr,nphi,nz)
       B_lines(mytdex, myline) = fval(1)
       IF (.not. lneut) THEN
       ! Calc dist func bins
          x0    = MOD(q(2),pi2)
          IF (x0 < 0) x0 = x0 + pi2
          IF (z0 < 0) z0 = z0 + pi2
          vperp = SQRT(2*moment*fval(1)/mymass)
          d1 = MAX(MIN(CEILING(rho_help*h1_prof      ), ns_prof1), 1) ! Rho Bin
          d2 = MAX(MIN(CEILING( z0*h2_prof           ), ns_prof2), 1) ! U Bin
          d3 = MAX(MIN(CEILING( x0*h3_prof           ), ns_prof3), 1) ! V Bin
          d4 = MAX(MIN(1+nsh_prof4+FLOOR(h4_prof*q(4)), ns_prof4), 1) ! vll
          d5 = MAX(MIN(CEILING(vperp*h5_prof         ), ns_prof5), 1) ! Vperp
          xw = weight(myline)*dt
          !CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,myworkid,0,win_dist5d,ier)
          dist5d_prof(mybeam,d1,d2,d3,d4,d5) = dist5d_prof(mybeam,d1,d2,d3,d4,d5) + xw
          !CALL MPI_WIN_UNLOCK(myworkid,win_dist5d,ier)
          IF (lfidasim_cyl) THEN
             !x0 = MOD(q(2), phimax_fida)
             !IF (x0 < 0) x0 = x0 + phimax_fida
            d1f = MIN(MAX(CEILING((q(1)-rmin_fida)*r_h),1),nr_fida)
            d2 = MIN(MAX(CEILING((x0-phimin_fida)*p_h),1),nphi_fida)
            d3 = MIN(MAX(CEILING((q(3)-zmin_fida)*z_h),1),nz_fida)
            y0 = MAX((q(4)**2+vperp**2),1.0)
            d4 = MIN(MAX(CEILING((y0*E_by_v-emin_fida)*e_h),1),nenergy_fida)
            d5 = MIN(MAX(CEILING((q(4)/SQRT(y0)-pimin_fida)*pi_h),1),npitch_fida)
            dist5d_fida(d1f,d3,d2,d4,d5) = dist5d_fida(d1f,d3,d2,d4,d5) + xw
          END IF
          IF (lcollision) CALL beams3d_physics_gc(t,q)
          IF (ltherm) THEN
             ndot_prof(mybeam,d1)   =   ndot_prof(mybeam,d1) + weight(myline)
             end_state(myline) = 1
             t = my_end
          END IF
       END IF
    ELSE
       IF (lneut) end_state(myline)=3
    END IF
    IF (lvessel .and. mytdex > 0 .and. rho_help > 0.5) THEN
       lhit = .false.
       x0    = xlast
       y0    = ylast
       z0    = zlast
       x1    = q(1)*cos(q(2))
       y1    = q(1)*sin(q(2))
       z1    = q(3)
       CALL collide(x0,y0,z0,x1,y1,z1,xw,yw,zw,lhit)
       IF (lhit) THEN
          q2(1) = SQRT(xw*xw+yw*yw)
          q2(2) = atan2(yw,xw)
          q2(3) = zw
          R_lines(mytdex,myline)       = q2(1)
          PHI_lines(mytdex,myline)     = q2(2)
          Z_lines(mytdex,myline)       = zw
          t = my_end+dt
          l = get_wall_ik()
          IF (lneut) THEN
             wall_shine(mybeam,l) = wall_shine(mybeam,l) + weight(myline)*0.5*mymass*q(4)*q(4)/get_wall_area(l)
          ELSE
             end_state(myline) = 2
             CALL fgc_eom(t,q2,qdot)
             qdot(4)=0
             wall_load(mybeam,l) = wall_load(mybeam,l) + weight(myline)*0.5*mymass*(SUM(qdot*qdot)+vperp*vperp)/get_wall_area(l)
          END IF
          IF (lhitonly) THEN
             R_lines(0,myline)      = SQRT(xlast*xlast+ylast*ylast)
             PHI_lines(0,myline)    = ATAN2(ylast,xlast)
             Z_lines(0,myline)      = zlast
             vll_lines(0,myline)    = q(4)
             moment_lines(0,myline) = moment
             neut_lines(0,myline)   = lneut
             R_lines(2,myline)      = q(1)
             PHI_lines(2,myline)    = q(2)
             Z_lines(2,myline)      = q(3)
             vll_lines(2,myline)    = q(4)
             moment_lines(2,myline) = moment
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
END SUBROUTINE out_beams3d_gc
