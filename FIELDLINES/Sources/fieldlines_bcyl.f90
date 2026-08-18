!-----------------------------------------------------------------
!     Function:      fieldlines_BCYL
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          01/07/2026
!     Description:   Returns Br, Bphi, Bz
!-----------------------------------------------------------------
SUBROUTINE fieldlines_BCYL(r,phi,z,Br,Bphi,Bz)
 !--------------------------------------------------------------
 !     Input Parameters
 !         r, phi, z     Cylindrical coordiantes
 !     Output Parameters
 !         Br, Bphi, Bz  Magnetic field components
 !--------------------------------------------------------------
 USE fieldlines_globals, ONLY: nr, nphi, nz, rmin, rmax, zmin, &
            zmax, phimin, phimax
 USE fieldlines_grid, ONLY : eps1, eps2, eps3, raxis, phiaxis, &
            zaxis, BR4D, BPHI4D, BZ4D
 IMPLICIT NONE
 DOUBLE PRECISION, INTENT(IN) :: r, phi, z
 DOUBLE PRECISION, INTENT(OUT) :: Br, Bphi, Bz

 !--------------------------------------------------------------
 !     Local Variables
 !        phi_temp     Helpers (r,phi,z)
 !        i,j,k      Spline Grid indicies
 !        xparam     Spline subgrid factor [0,1] (yparam,zparam)
 !        ict        Spline output control
 !        fval       Spline output array
 !--------------------------------------------------------------
 DOUBLE PRECISION :: phi_temp
 ! For splines
 INTEGER :: i,j,k
 REAL*8 :: xparam, yparam, zparam, hx, hy, hz, hxi, hyi, hzi
 INTEGER, parameter :: ict(8)=(/1,0,0,0,0,0,0,0/)
 REAL*8 :: fval(1)
 REAL*8, PARAMETER :: one = 1

 !--------------------------------------------------------------
 !     Begin Subroutine
 !--------------------------------------------------------------

 ! Setup position in a vll arrays
 phi_temp = MODULO(phi, phimax)

 ! Initialize values
 Br = 0
 Bphi = 0
 Bz = 0

 ! Check that we're inside the domain then proceed
 IF ((r >= rmin-eps1) .and. (r <= rmax+eps1) .and. &
     (phi_temp >= phimin-eps2) .and. (phi_temp <= phimax+eps2) .and. &
     (z >= zmin-eps3) .and. (z <= zmax+eps3)) THEN
    i = MIN(MAX(COUNT(raxis < r),1),nr-1)
    j = MIN(MAX(COUNT(phiaxis < phi_temp),1),nphi-1)
    k = MIN(MAX(COUNT(zaxis < z),1),nz-1)
    hx     = raxis(i+1) - raxis(i)
    hy     = phiaxis(j+1) - phiaxis(j)
    hz     = zaxis(k+1) - zaxis(k)
    hxi    = one / hx
    hyi    = one / hy
    hzi    = one / hz
    xparam = (r - raxis(i)) * hxi
    yparam = (phi_temp - phiaxis(j)) * hyi
    zparam = (z - zaxis(k)) * hzi
    ! Evaluate the Splines
    CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hx,hxi,hy,hyi,hz,hzi,&
                         BR4D(1,1,1,1),nr,nphi,nz)
    Br = fval(1)
    CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hx,hxi,hy,hyi,hz,hzi,&
                         BPHI4D(1,1,1,1),nr,nphi,nz)
    Bphi = fval(1)
    CALL R8HERM3FCN(ict,1,1,fval,i,j,k,xparam,yparam,zparam,&
                         hx,hxi,hy,hyi,hz,hzi,&
                         BZ4D(1,1,1,1),nr,nphi,nz)
    Bz = fval(1)
 ELSE
    RETURN
 END IF

 RETURN

END SUBROUTINE fieldlines_BCYL
