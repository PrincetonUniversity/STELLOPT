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
      USE beams3d_runtime, ONLY: lfidasim, nbeams, pi2, &
                              EZSPLINE_ERR
      USE beams3d_grid, ONLY: raxis, phiaxis, zaxis, X4D, Y4D, S4D, &
                              nr, nphi, nz
      USE beams3d_lines, ONLY: ns_prof1, ns_prof2, ns_prof3, ns_prof4, &
                               ns_prof4, ns_prof5, dist5d_prof, &
                               partvmax
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
      INTEGER ::  s, i, j, k, nvol, l, m, n, ier
      INTEGER, DIMENSION(2) :: minln
      REAL(rprec) :: s1, s2, u1, u2, p1, drho, du, dp, area, dvol, a, b, c, d, f1, f2, f3
      REAL(rprec), DIMENSION(4) :: rt,zt,pt
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: targ
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: RHO3D

      INTEGER :: bcs1(2), bcs2(2), bcs3(2)
      !TYPE(EZspline3_r8) :: X_spl, Y_spl
      !REAL(rprec), POINTER, DIMENSION(:,:,:) :: X_ARR, Y_ARR
      !REAL(rprec), POINTER, DIMENSION(:,:,:,:) :: X4D, Y4D
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------

         !Set up Splines on the X/Y grid to fix singular interpolation and gradients

      ! ALLOCATE(X4D(8, nr, nphi, nz), Y4D(8, nr, nphi, nz))
      ! ALLOCATE(X_ARR(nr,nphi,nz),Y_ARR(nr,nphi,nz))
      ! X_ARR = S4D(1,:,:,:) * COS(U4D(1,:,:,:))
      ! Y_ARR = S4D(1,:,:,:) * SIN(U4D(1,:,:,:))
      ! bcs1=(/ 0, 0/)
      ! bcs2=(/-1,-1/)
      ! bcs3=(/ 0, 0/)
      ! ier = 0
      ! CALL EZspline_init(X_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
      ! IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_fix_poloidal:X_spl',ier)
      ! CALL EZspline_init(Y_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
      ! IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_fix_poloidal:Y_spl',ier)
      ! X_spl%isHermite   = 1
      ! X_spl%x1   = raxis
      ! X_spl%x2   = phiaxis
      ! X_spl%x3   = zaxis
      ! Y_spl%isHermite   = 1
      ! Y_spl%x1   = raxis
      ! Y_spl%x2   = phiaxis
      ! Y_spl%x3   = zaxis
      ! CALL EZspline_setup(X_spl,X_ARR,ier,EXACT_DIM=.true.)
      ! IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_fix_poloidal:X_spl',ier)
      ! CALL EZspline_setup(Y_spl,Y_ARR,ier,EXACT_DIM=.true.)
      ! IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_fix_poloidal:Y_spl',ier)
      ! X4D = X_spl%fspl
      ! Y4D = Y_spl%fspl
      ! CALL EZspline_free(X_spl,ier)
      ! CALL EZspline_free(Y_spl,ier)
      ! DEALLOCATE(X_ARR,Y_ARR)

      ! Do physical volume elements
      nvol = (ns_prof1)*(ns_prof2)*(ns_prof3)
      drho   = 1.0/REAL(ns_prof1) !distribution defined on centers (half grid)
      du   = pi2/REAL(ns_prof2)
      dp   = pi2/REAL(ns_prof3)
      !ALLOCATE(targ(nr,nz))
      ! ALLOCATE(RHO3D(nr,nphi,nz))
      ! RHO3D = S4D(1,:,:,:)
      ! WHERE (RHO3D>1) RHO3D = 2;
      ! RHO3D = sqrt(RHO3D) ! in RHO
      DO s = 1, nvol
         i = MOD(s-1,(ns_prof1))+1
         j = MOD(s-1,(ns_prof1)*(ns_prof2))
         j = FLOOR(REAL(j) / REAL(ns_prof1))+1
         k = CEILING(REAL(s) / REAL(ns_prof1*ns_prof2))
         s1 = ((i-1)*drho)**2 !Distribution function has radial rho coordinate
         s2 = ((i-1)*drho+drho)**2
         u1 = (j-1)*du
         u2 = u1+du
         p1 = MOD((k-1)*dp*0.5,phiaxis(nphi))
         ! Find helper value of grids
         m = MIN(MAX(COUNT(phiaxis < p1),1),nphi-1)
         rt = 0!5.25+COS((u1+u2)/2)*.5 !dirty estimation of R
         zt = 0!SIN((u1+u2)/2.)*.6 
         pt = 0; !Initial conditions for successful lookup?

         ! Assume parallel piped
         CALL beams3d_suv2rzp(s1,u1,p1,rt(1),zt(1),pt(1))
         CALL beams3d_suv2rzp(s1,u2,p1,rt(2),zt(2),pt(2))
         CALL beams3d_suv2rzp(s2,u1,p1,rt(3),zt(3),pt(3))
         CALL beams3d_suv2rzp(s2,u2,p1,rt(4),zt(4),pt(4))
            area = 0.5 * abs( (rt(2) - rt(1)) *(zt(3) - zt(1)) - (rt(3) - rt(1)) *(zt(2) - zt(1))) + &
                  0.5 * abs( (rt(3) - rt(2)) *(zt(4) - zt(2)) - (rt(4) - rt(2)) *(zt(3) - zt(2)))

         WRITE(327,*) m,s1,s2,u1,u2,rt,zt,area
!         CALL FLUSH(327)
         dvol = area*sum(rt)*dp/4
         dist5d_prof(:,i,j,k,:,:) = dist5d_prof(:,i,j,k,:,:)/dvol
!        WRITE(328,*) i,j,k,dvol,rt,zt
!         CALL FLUSH(328)
      END DO

      IF (lfidasim) THEN
            CALL beams3d_write_fidasim('DENF') !write density before velocity space normalization
      END IF

      ! Do phase space volume elements -- correct?
      drho = 2*partvmax/ns_prof4 !vll
      du = partvmax/ns_prof5 !v_perp
      dvol = pi2*drho*du
      nvol = ns_prof4*ns_prof5
      DO j = 1, ns_prof5
        u1 = (j-0.5)*du
        dist5d_prof(:,:,:,:,:,j) = dist5d_prof(:,:,:,:,:,j)/(dvol*u1)
!        WRITE(329,*) j,dvol*u1
!        CALL FLUSH(329)
      END DO

      
      !DEALLOCATE(targ)
      !DEALLOCATE(RHO3D)


      RETURN
!-----------------------------------------------------------------------
!     END SUBROUTINE
!-----------------------------------------------------------------------
      END SUBROUTINE beams3d_distnorm