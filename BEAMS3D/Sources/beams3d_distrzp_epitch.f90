!-----------------------------------------------------------------------
!     Module:        beams3d_distrzp_epitch
!     Authors:       D. Kulla (david.kulla@ipp.mpg.de)
!     Date:          24/01/2022
!     Description:   This subroutine converts the distribution function
!                    to cylindrical and Energy/pitch coordinates for FIDASIM
!                    !! (dummy)It is presently not called from the code !!
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_distrzp_epitch
        !-----------------------------------------------------------------------
        !     Libraries
        !-----------------------------------------------------------------------
              USE stel_kinds, ONLY: rprec
              USE beams3d_runtime, ONLY: nbeams, pi2, mass_beams
              USE beams3d_grid, ONLY: raxis, phiaxis, zaxis, S4D, U4D, &
                                      nr, nphi, nz, hr, hp, hz, hri, hpi, hzi
              USE beams3d_lines, ONLY: ns_prof1, ns_prof2, ns_prof3, ns_prof4, &
                                       ns_prof4, ns_prof5, dist5d_prof, &
                                       partvmax, dist5D_fida, &
                                       h2_prof, h3_prof, h4_prof, h5_prof, my_end
              USE beams3d_physics_mod, ONLY: beams3d_suv2rzp
              !USE beams3d_write_fidasim, ONLY: nr_fida, nphi_fida, nz_fida,rmin_fida, &
              !rmax_fida, zmin_fida, zmax_fida, phimin_fida, &
              !phimax_fida
              IMPLICIT NONE
              
        !-----------------------------------------------------------------------
        !     Input Variables
        !        NONE
        !-----------------------------------------------------------------------
              
        !-----------------------------------------------------------------------
        !     Local Variables
        !        s,i,j,k     Helper index
        !        nvol        Number of volume voxel vertices
        !        s1,s2       Radial Helper
        !        u1,u2       Poloidal Helper
        !        p1,p2       Toroidal Helper
        !        ds,du,dp    Delta coordiantes
        !        xt,yt,zt    Helper for xyz coordiantes of voxel
        !-----------------------------------------------------------------------
              INTEGER ::  s, i, j, k, nvol, l, m, n,b, i3, j3, k3
              INTEGER, DIMENSION(2) :: minln
              !REAL(rprec) :: s1, s2, u1, u2, p1, ds, du, dp, area, dvol, a, b, c, d, f1, f2, f3
              REAL(rprec), DIMENSION(4) :: rt,zt,pt
              REAL*8 :: xparam, yparam, zparam !, hx, hy, hz, hxi, hyi, hzi
              REAL*8 :: fval(1)
              INTEGER, parameter :: ict(8)=(/1,0,0,0,0,0,0,0/)
              DOUBLE PRECISION         :: x0,y0,z0, jac, v_parr, v_perp, pitch
        
        !-----------------------------------------------------------------------
        !     Begin Subroutine
        !-----------------------------------------------------------------------
            ! Do volume normalization
            CALL beams3d_distnorm !TODO: check if this has already been done   

         !Allocate with Radial-like dimensions for clean transfer and to avoid explicitly looping over every element
        ALLOCATE(dist5d_fida(nbeams,ns_prof1, ns_prof2, ns_prof3, nenergy_fida, npitch_fida))
        DO b=1,nbeams
           DO d1 = 1, nenergy_fida
              DO d2 = 1, npitch_fida
                 v = SQRT(2 * energy_fida(d1) *1000.0 * e_charge / mass_beams(b))
                 IF (v .gt. partvmax) THEN !) .or. (v .eq. 0.0))
                    dist5d_fida(b,:,:,:,d1,d2) = 0
                 ELSE
                    pitch = pitch_fida(d2)
                    v_parr = pitch * v
                    v_perp = SQRT(1- pitch * pitch) * v
                    !determine beams3d-grid indices (velocity space)
                    i3 = MAX(MIN(1+nsh_prof4+FLOOR(h4_prof*v_parr), ns_prof4), 1) ! vll
                    j3 = MAX(MIN(CEILING(v_perp*h5_prof         ), ns_prof5), 1) ! Vperp
                    jac = 1 / (mass_beams(b) * SQRT(1-pitch * pitch))
                    !jac = MIN(MAX(jac, 0.0),1.0E)
                    jac = jac / pi2 / REAL(1000000) * e_charge / 1000 !convert to TRANSP convention? and 1/cm^3/keV
                    dist5d_fida(b,:,:,:,d1,d2) = dist5d_prof(b,:,:,:,i3,j3) * jac !* pi2 * v /  mass_beams(b) !problematic, circular reference
                 END IF
              END DO
           END DO
        END DO
        !DEALLOCATE(dist5d_prof)
        !ALLOCATE(dist5d_prof(nbeams,nr, nphi, nz, nenergy_fida, npitch_fida))
        dist5d_prof = dist5d_fida
        DEALLOCATE(dist5d_fida)

        !Now allocate with correct dimensions
        ALLOCATE(dist5d_fida(nbeams,nenergy_fida,npitch_fida,nr_fida,nz_fida,nphi_fida))
        DO b=1,nbeams
           DO i=1,nr_fida
              DO k = 1, nz_fida
                 DO j=1,nphi_fida
                    !convert i,j,k to distribution function indices l,m,n
                    !determine beams3d-grid indices
                    i3 = MIN(MAX(COUNT(raxis < raxis_fida(i)),1),nr_fida-1)
                    j3 = MIN(MAX(COUNT(phiaxis < phiaxis_fida(j)),1),nphi_fida-1)
                    k3 = MIN(MAX(COUNT(zaxis < zaxis_fida(k)),1),nz_fida-1)
                    !setup interpolation
                    xparam = (raxis_fida(i) - raxis(i3)) * hri(i3)
                    yparam = (phiaxis_fida(j) - phiaxis(j3)) * hpi(j3)
                    zparam = (zaxis_fida(k) - zaxis(k3)) * hzi(k3)
                    CALL R8HERM3FCN(ict,1,1,fval,i3,j3,k3,xparam,yparam,zparam,&!maybe switch to x/y interpolation?
                       hr(i3),hri(i3),hp(j3),hpi(j3),hz(k3),hzi(k3),&
                       S4D(1,1,1,1),nr,nphi,nz)
                    y0 = fval(1)
                    CALL R8HERM3FCN(ict,1,1,fval2,i3,j3,k3,xparam,yparam,zparam,&
                       hr(i3),hri(i3),hp(j3),hpi(j3),hz(k3),hzi(k3),&
                       U4D(1,1,1,1),nr,nphi,nz)
                    z0 = fval2(1)

                    IF (z0 < 0) z0 = z0 + pi2
                    IF (x0 < 0) x0 = x0 + pi2
                    IF (y0 < 0) y0 = -y0
                    ! Calc dist func bins
                    x0    = phiaxis_fida(j)
                    l = MAX(MIN(CEILING(SQRT(y0)*ns_prof1     ), ns_prof1), 1) ! Rho Bin
                    m = MAX(MIN(CEILING( z0*h2_prof           ), ns_prof2), 1) ! U Bin
                    n = MAX(MIN(CEILING( x0*h3_prof           ), ns_prof3), 1) ! V Bin
                    IF (y0 .GT. 1.0) THEN !might introduce a small deviation here
                       dist5d_fida(b,:,:,i,k,j) = 0 !distribution is 0 outside plasma
                    ELSE
                       dist5d_fida(b,:,:,i,k,j) = dist5d_prof(b,l,m,n,:,:) !output in r-z-phi
                    END IF

                 END DO
              END DO
           END DO
        END DO

        !Still would have to update dist5d_prof if that is wanted
        RETURN
        !-----------------------------------------------------------------------
        !     END SUBROUTINE
        !-----------------------------------------------------------------------
        END SUBROUTINE beams3d_distrzp_epitch