!-----------------------------------------------------------------------
!     Module:        beams3d_distnorm
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          03/07/2021
!     Description:   This subroutine applies the volume normalization
!                    to the distribution function.
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_distnorm
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE beams3d_runtime, ONLY: lfidasim, lfidasim2, nbeams, pi2, &
                              EZSPLINE_ERR
      USE beams3d_grid, ONLY: raxis, phiaxis, zaxis, X4D, Y4D, S4D, &
                              nr, nphi, nz
      USE beams3d_lines, ONLY: ns_prof1, ns_prof2, ns_prof3, ns_prof4, &
                               ns_prof4, ns_prof5, dist5d_prof, &
                               partvmax, dist5d_fida
      USE beams3d_physics_mod, ONLY: beams3d_suv2rzp
      USE EZspline_obj
      USE EZspline
      USE mpi_params !Used for call to write_fidasim
      USE mpi_inc
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
!        drho,du,dp    Delta coordiantes
!        xt,yt,zt    Helper for xyz coordiantes of voxel
!-----------------------------------------------------------------------
      INTEGER ::  s, i, j, k, nvol, m
      INTEGER, DIMENSION(2) :: minln
      REAL(rprec) :: rho1, rho2, s1, s2, u1, u2, p1, ds, du, dp, area, dvol
      REAL(rprec), DIMENSION(4) :: rt,zt,pt

      INTEGER :: bcs1(2), bcs2(2), bcs3(2)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------

      ! Do physical volume elements
      nvol = (ns_prof1)*(ns_prof2)*(ns_prof3)
      ds   = 1.0/REAL(ns_prof1) !distribution defined on centers (half grid)
      du   = pi2/REAL(ns_prof2)
      dp   = pi2/REAL(ns_prof3)

      !Initial condition near axis
      rt = -1 ! will then use 0.75 point of R-grid
      zt = 0
      pt = 0

      ! Iterate over volumes
      DO s = 1, nvol
         ! Distribution grid index
         i = MOD(s-1,(ns_prof1))+1
         j = MOD(s-1,(ns_prof1)*(ns_prof2))
         j = FLOOR(REAL(j) / REAL(ns_prof1))+1
         k = CEILING(REAL(s) / REAL(ns_prof1*ns_prof2))

         ! Bounding values of rho and u, p is centered in voxel
         rho1 = (i-1)*ds
         rho2 = rho1+ds
         u1 = (j-1)*du
         u2 = u1+du
         p1 = MOD((k-1)*dp*0.5,phiaxis(nphi))

         ! rho to s
         s1 = rho1*rho1
         s2 = rho2*rho2

         ! Get first gridpoint
         CALL beams3d_suv2rzp(s1,u1,p1,rt(1),zt(1),pt(1))

         ! All points are nearby so use as guess
         rt(2:4) = rt(1)
         zt(2:4) = zt(1)

         ! If not the axis point then we evaluate the second s1 point
         IF (i .ne. 1) CALL beams3d_suv2rzp(s1,u2,p1,rt(2),zt(2),pt(2))

         ! Other two points
         CALL beams3d_suv2rzp(s2,u1,p1,rt(3),zt(3),pt(3))
         CALL beams3d_suv2rzp(s2,u2,p1,rt(4),zt(4),pt(4))

         ! Assume parallel piped
            area = 0.5 * abs( (rt(2) - rt(1)) *(zt(3) - zt(1)) - (rt(3) - rt(1)) *(zt(2) - zt(1))) + &
                  0.5 * abs( (rt(3) - rt(2)) *(zt(4) - zt(2)) - (rt(4) - rt(2)) *(zt(3) - zt(2)))

!         m = MIN(MAX(COUNT(phiaxis < p1),1),nphi-1)
!         WRITE(327,*) m,s1,s2,u1,u2,rt,zt,area
!         CALL FLUSH(327)
         dvol = area*sum(rt)*dp*0.25
         dist5d_prof(:,i,j,k,:,:) = dist5d_prof(:,i,j,k,:,:)/dvol
!        WRITE(328,*) i,j,k,dvol,rt,zt
!         CALL FLUSH(328)
      END DO

      IF (lfidasim) THEN
            CALL beams3d_write_fidasim('DENF') !write density before velocity space normalization
      END IF

      ! Do phase space volume elements 
      ds = 2.0*partvmax/ns_prof4
      du =   partvmax/ns_prof5
      dvol = pi2*ds*du
      nvol = ns_prof4*ns_prof5
      DO j = 1, ns_prof5
        u1 = REAL(j-0.5)*du
        dist5d_prof(:,:,:,:,:,j) = dist5d_prof(:,:,:,:,:,j)/(dvol*u1)
        IF (lfidasim2) THEN
            dist5d_fida(:,:,:,:,:,j) = dist5d_fida(:,:,:,:,:,j)/(dvol*u1)
        END IF
!        WRITE(329,*) j,dvol*u1
!        CALL FLUSH(329)
      END DO

      


      RETURN
!-----------------------------------------------------------------------
!     END SUBROUTINE
!-----------------------------------------------------------------------
      END SUBROUTINE beams3d_distnorm