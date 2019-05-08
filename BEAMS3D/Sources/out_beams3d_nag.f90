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
                               lhitonly, npoinc, lcollision, ldepo
    USE beams3d_lines, ONLY: R_lines, Z_lines, PHI_lines, myline, moment, &
                             nsteps, nparticles, moment_lines, myend, &
                             vll_lines, neut_lines, mytdex, next_t,&
                             lost_lines, dt_out, xlast, ylast, zlast,&
                             ltherm, S_lines, U_lines, B_lines
    USE beams3d_grid
!DEC$ IF DEFINED (NTCC)
    USE beams3d_physics_mod, ONLY: beams3d_physics
!DEC$ ENDIF
    USE wall_mod, ONLY: collide
    USE EZspline_obj
    USE EZspline
    USE mpi_params
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
    !REAL(rprec)         :: x0,y0,z0,x1,y1,z1,xw,yw,zw,delta,dl
    INTEGER             :: ier
    DOUBLE PRECISION         :: x0,y0,z0,x1,y1,z1,xw,yw,zw,delta,dl
    LOGICAL             :: lhit
    ! For splines
    INTEGER :: i,j,k
    REAL*8 :: xparam, yparam, zparam, hx, hy, hz, hxi, hyi, hzi
    REAL*8 :: fval(4)
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
       i = COUNT(raxis < q(1))
       j = COUNT(phiaxis < x0)
       k = COUNT(zaxis < q(3))
       IF (i==0) i=1
       IF (j==0) j=1
       IF (k==0) k=1
       hx     = raxis(i+1) - raxis(i)
       hy     = phiaxis(j+1) - phiaxis(j)
       hz     = zaxis(k+1) - zaxis(k)
       hxi    = one / hx
       hyi    = one / hy
       hzi    = one / hz
       xparam = (raxis(i+1) - q(1)) * hxi
       yparam = (phiaxis(j+1) - x0) * hyi
       zparam = (zaxis(k+1) - q(3)) * hzi
       CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                       hx,hxi,hy,hyi,hz,hzi,&
                       S4D(1,1,1,1),nr,nphi,nz)
       S_lines(mytdex, myline) = fval(1)
       CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                       hx,hxi,hy,hyi,hz,hzi,&
                       U4D(1,1,1,1),nr,nphi,nz)
       U_lines(mytdex, myline) = fval(1)
       CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                       hx,hxi,hy,hyi,hz,hzi,&
                       MODB4D(1,1,1,1),nr,nphi,nz)
       B_lines(mytdex, myline) = fval(1)
       !CALL EZspline_interp(S_spl,q(1),x0,q(3),y0,ier)
       !CALL EZspline_interp(U_spl,q(1),x0,q(3),z0,ier)
       !CALL EZspline_interp(MODB_spl,q(1),x0,q(3),x1,ier)
       !IF (myworkid == 0) PRINT *,'--',y0,z0,x1
    END IF
!DEC$ IF DEFINED (NTCC)
    IF (lcollision) CALL beams3d_physics(t,q)
    IF (ltherm) t = t_end(myline)
!DEC$ ENDIF
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
          R_lines(mytdex,myline)       = SQRT(xw*xw+yw*yw)
          PHI_lines(mytdex,myline)     = atan2(yw,xw)
          Z_lines(mytdex,myline)       = zw
          vll_lines(mytdex,myline)     = 0.5*(vll_lines(mytdex-1,myline)+vll_lines(mytdex,myline))
          moment_lines(mytdex,myline)     = 0.5*(moment_lines(mytdex-1,myline)+moment_lines(mytdex,myline))
          lost_lines(myline) = .TRUE.
          t = t_end(myline)+dt
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
    IF (lhitonly) mytdex = 0
    IF ((t+dt) .ge. mytdex*dt_out) mytdex = mytdex + 1
    t = t + dt

    RETURN
    !-----------------------------------------------------------------------
    !     End Function
    !-----------------------------------------------------------------------
END SUBROUTINE out_beams3d_nag

