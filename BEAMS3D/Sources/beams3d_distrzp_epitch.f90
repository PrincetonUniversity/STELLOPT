!-----------------------------------------------------------------------
!     Module:        beams3d_distrzp_epitch
!     Authors:       D. Kulla (david.kulla@ipp.mpg.de)
!     Date:          24/01/2022
!     Description:   This subroutine converts the distribution function
!                    to cylindrical and Energy/pitch coordinates for FIDASIM
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

            !Apply jacobian for transformation from v_parr/v_perp to E,p
            i3 = ns_prof4/2
            DO b = 1, nbeams
               DO i = 1, ns_prof4
                  DO j = 1, ns_prof5
                     ! Half grid
                      v_parr = partvmax * REAL(i - i3 -0.5) / i3
                      v_perp = partvmax * (j - 0.5) / ns_prof5 
                      pitch = v_parr / SQRT(v_parr * v_parr + v_perp * v_perp)
                      jac = 1 / (mass_beams(k) * SQRT(1-pitch * pitch))
                      dist5d_prof(b,:,:,:,i,j) = dist5d_prof(b,:,:,:,i,j) * jac !probably should inicate somewhere that this has been done, or put it in another variable
                  END DO
               END DO
            END DO

           !Todo (nearest neighbor) interpolation to background grid
           ALLOCATE(dist5d_fida(nbeams,nr,nz,nphi,ns_prof4,ns_prof5))
           DO b=1,nbeams
              DO i=1,nr
                 DO j=1,nz
                    DO k = 1, nphi
                    !convert i,j,k to distribution function indices l,m,n

                    !determine beams3d-grid indices
                    i3 = MIN(MAX(COUNT(raxis < raxis_fida(i)),1),nr-1) 
                    k3 = MIN(MAX(COUNT(phiaxis < phiaxis_fida(k)),1),nphi-1)
                    j3 = MIN(MAX(COUNT(zaxis < zaxis_fida(j)),1),nz-1)
                    !setup interpolation
                    xparam = (raxis_fida(i) - raxis(i3)) * hri(i3)
                    yparam = (phiaxis_fida(k) - phiaxis(j3)) * hpi(j3)
                    zparam = (zaxis_fida(j) - zaxis(k3)) * hzi(k3)
                    CALL R8HERM3FCN(ict,1,1,fval,i3,j3,k3,xparam,yparam,zparam,&
                                             hr(i3),hri(i3),hp(j3),hpi(j3),hz(k3),hzi(k3),&
                                             S4D(1,1,1,1),nr,nphi,nz)
                    y0 = fval(1)
                    CALL R8HERM3FCN(ict,1,1,fval,i3,j3,k3,xparam,yparam,zparam,&
                                             hr(i3),hri(i3),hp(j3),hpi(j3),hz(k3),hzi(k3),&
                                             U4D(1,1,1,1),nr,nphi,nz)
                    z0 = fval(1)
                    ! Calc dist func bins
                    x0    = phiaxis_fida(k)
                    IF (x0 < 0) x0 = x0 + pi2
                    k = MAX(MIN(CEILING(SQRT(y0)*ns_prof1     ), ns_prof1), 1) ! Rho Bin
                    l = MAX(MIN(CEILING( z0*h2_prof           ), ns_prof2), 1) ! U Bin
                    m = MAX(MIN(CEILING( x0*h3_prof           ), ns_prof3), 1) ! V Bin
                    dist5d_fida(b,i,j,k,:,:) = dist5d_prof(b,k,l,m,:,:)
                 END DO
              END DO 
           END DO
        END DO        
        
        PRINT *,'GOT HERE'
        RETURN
        !-----------------------------------------------------------------------
        !     END SUBROUTINE
        !-----------------------------------------------------------------------
        END SUBROUTINE beams3d_distrzp_epitch