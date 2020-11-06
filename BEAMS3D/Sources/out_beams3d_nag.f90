!-----------------------------------------------------------------------
!     Function:      out_beams3d_nag
!     Authors:       S. Lazerson (lazerson@pppl.gov), M. McMillan (matthew.mcmillan@my.wheaton.edu)
!     Date:          06/20/2012
!     Description:   Save output from field line following while running
!                    and updates progress bar.
!-----------------------------------------------------------------------
SUBROUTINE out_beams3d_nag(t, q)
    !-----------------------------------------------------------------------
    !     Libraries
    !-----------------------------------------------------------------------
    USE stel_kinds, ONLY: rprec
    USE beams3d_runtime, ONLY: dt, lverb, pi2, lneut, t_end, lvessel, &
                               lhitonly, npoinc, lcollision, ldepo, &
                               weight, invpi2, ndt, ndt_max
    USE beams3d_lines, ONLY: R_lines, Z_lines, PHI_lines, myline, moment, &
                             nsteps, nparticles, moment_lines, myend, &
                             vll_lines, neut_lines, mytdex, next_t,&
                             xlast, ylast, zlast, dense_prof, &
                             ltherm, S_lines, U_lines, B_lines, &
                             j_prof, ndot_prof, partvmax, &
                             ns_prof1, ns_prof2, ns_prof3, ns_prof4, &
                             ns_prof5, mymass, mycharge, mybeam, end_state, &
                             dist5d_prof, win_dist5d, nsh_prof4, &
                             h2_prof, h3_prof, h4_prof, h5_prof
    USE beams3d_grid
    USE beams3d_physics_mod, ONLY: beams3d_physics
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
    INTEGER             :: ier, d1, d2, d3, d4, d5
    DOUBLE PRECISION         :: x0,y0,z0,x1,y1,z1,xw,yw,zw, vperp
    DOUBLE PRECISION    :: q2(4),qdot(4)
    ! For splines
    INTEGER :: i,j,k,l
    REAL*8 :: xparam, yparam, zparam !, hx, hy, hz, hxi, hyi, hzi
    REAL*8 :: fval(1)
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
    y0 = 0  ! If we're out of domain then don't worry about collisions
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
                       S4D(1,1,1,1),nr,nphi,nz)
       y0 = fval(1)
       S_lines(mytdex, myline) = y0 
       CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                       hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                       U4D(1,1,1,1),nr,nphi,nz)
       z0 = fval(1)
       U_lines(mytdex, myline) = z0
       CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                       hr(i),hri(i),hp(j),hpi(j),hz(k),hzi(k),&
                       MODB4D(1,1,1,1),nr,nphi,nz)
       B_lines(mytdex, myline) = fval(1)
       ! Calc dist func bins
       x0    = MOD(q(2),pi2)
       IF (x0 < 0) x0 = x0 + pi2
       vperp = SQRT(2*moment*fval(1)/mymass)
       d1 = MAX(MIN(CEILING(SQRT(y0)*ns_prof1     ), ns_prof1), 1) ! Rho Bin
       d2 = MAX(MIN(CEILING( z0*h2_prof           ), ns_prof2), 1) ! U Bin
       d3 = MAX(MIN(CEILING( x0*h3_prof           ), ns_prof3), 1) ! V Bin
       d4 = MAX(MIN(1+nsh_prof4+FLOOR(h4_prof*q(4)), ns_prof4), 1) ! vll
       d5 = MAX(MIN(CEILING(vperp*h5_prof         ), ns_prof5), 1) ! Vperp
       xw = weight(myline)*dt
       !j_prof(mybeam,d1)      =      j_prof(mybeam,d1) + mycharge*q(4)*xw
       !CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,myworkid,0,win_dist5d,ier)
       dist5d_prof(mybeam,d1,d2,d3,d4,d5) = dist5d_prof(mybeam,d1,d2,d3,d4,d5) + xw
       !CALL MPI_WIN_UNLOCK(myworkid,win_dist5d,ier)
       IF (lcollision) CALL beams3d_physics(t,q)
       IF (ltherm) THEN
          ndot_prof(mybeam,d1)   =   ndot_prof(mybeam,d1) + weight(myline)
          end_state(myline) = 1
          t = t_end(myline)
       END IF
    ELSE
       IF (lneut) end_state(myline)=3
    END IF
    IF (lvessel .and. mytdex > 0 .and. y0 > 0.5) THEN
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
          t = t_end(myline)+dt
          l = get_wall_ik()
          IF (lneut) THEN
             wall_shine(mybeam,l) = wall_shine(mybeam,l) + weight(myline)*0.5*mymass*q(4)*q(4)/get_wall_area(l)
          ELSE
             end_state(myline) = 2
             CALL fpart_nag(t,q2,qdot)
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
    IF (mytdex > npoinc) t = t_end(myline)
    t = t + dt

    RETURN
    !-----------------------------------------------------------------------
    !     End Function
    !-----------------------------------------------------------------------
END SUBROUTINE out_beams3d_nag

