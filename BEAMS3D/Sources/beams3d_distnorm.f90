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
      USE beams3d_runtime, ONLY: lfidasim, lfidasim_cyl, nbeams, pi2, &
                              EZSPLINE_ERR, MPI_BARRIER_ERR
      USE beams3d_grid, ONLY: raxis, phiaxis, zaxis, X4D, Y4D, S4D, &
                              nr, nphi, nz,VOL_ARR, win_VOL_ARR
      USE beams3d_lines, ONLY: ns_prof1, ns_prof2, ns_prof3, ns_prof4, &
                               ns_prof4, ns_prof5, dist5d_prof, &
                               partvmax, dist5d_fida, h1_prof, h2_prof,&
                               h3_prof
      USE beams3d_physics_mod, ONLY: beams3d_suv2rzp
      USE fidasim_input_mod, ONLY: beams3d_write_fidasim
      USE EZspline_obj
      USE EZspline
      USE mpi_params !Used for call to write_fidasim
      USE mpi_inc
      USE mpi_sharmem      
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
      INTEGER ::  s, i, j, k, nvol, mystart, myend
      REAL(rprec) :: rho1, rho2, s1, s2, u1, u2, p1, ds, du, dp, area, dvol
      REAL(rprec), DIMENSION(4) :: rt,zt,pt
      INTEGER :: numprocs_local, mylocalid, mylocalmaster
      INTEGER :: MPI_COMM_LOCAL
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      
      ! Setup Local communicator
      ! Divide up Work
      mylocalid = myworkid
      numprocs_local = 1
#if defined(MPI_OPT)
      CALL MPI_COMM_DUP( MPI_COMM_SHARMEM, MPI_COMM_LOCAL, ierr_mpi)
      CALL MPI_COMM_RANK( MPI_COMM_LOCAL, mylocalid, ierr_mpi )              ! MPI
      CALL MPI_COMM_SIZE( MPI_COMM_LOCAL, numprocs_local, ierr_mpi )         ! MPI
      CALL MPI_BARRIER(MPI_COMM_LOCAL, ierr_mpi)
#endif
      ! Do physical volume elements
      nvol = (ns_prof1)*(ns_prof2)*(ns_prof3)
      ds   = 1.0/h1_prof
      du   = 1.0/h2_prof
      dp   = 1.0/h3_prof
!      ds   = 1.0/REAL(ns_prof1) !distribution defined on centers (half grid)
!      du   = pi2/REAL(ns_prof2)
!      dp   = pi2/REAL(ns_prof3)

      !Initial condition near axis
      rt = -1 ! will then use 0.75 point of R-grid
      zt = 0
      pt = 0


      ! Allocate VOL_ARR
      IF (.NOT. ASSOCIATED(VOL_ARR)) THEN 
#if defined(MPI_OPT)         
      CALL mpialloc(VOL_ARR, ns_prof1, ns_prof2, ns_prof3, myid_sharmem, 0, MPI_COMM_LOCAL, win_VOL_ARR)
#else
      ALLOCATE(VOL_ARR(ns_prof1, ns_prof2, ns_prof3))
#endif
         VOL_ARR=0.0
      END IF

      ! Divide up work
      CALL MPI_CALC_MYRANGE(MPI_COMM_LOCAL, 1, nvol, mystart, myend)

      ! Iterate over volumes
      DO s = mystart, myend
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
         p1 = (k-0.5)*dp

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
            area = abs( (rt(2) - rt(1)) *(zt(3) - zt(1)) - (rt(3) - rt(1)) *(zt(2) - zt(1))) + &
                   abs( (rt(3) - rt(2)) *(zt(4) - zt(2)) - (rt(4) - rt(2)) *(zt(3) - zt(2)))

         dvol = area*sum(rt)
         VOL_ARR(i, j, k) = dvol/8.0
         dist5d_prof(:,i,j,k,:,:) = dist5d_prof(:,i,j,k,:,:)/dvol
      END DO

      IF (myworkid == master) THEN
         ! Constants moved here
         ! Factor 4 from average of rt
         ! Factor 2 from parallel piped
         dist5d_prof = dist5d_prof * 8 / dp
         ! Write Density before velocity space normalization
         IF (lfidasim) CALL beams3d_write_fidasim('DENF')
      END IF

      ! Do phase space volume elements
      IF (myworkid == master) THEN
         ds = 2.0*partvmax/ns_prof4
         du =   partvmax/ns_prof5
         dvol = pi2*ds*du
         DO j = 1, ns_prof5
            u1 = REAL(j-0.5)*du
            dist5d_prof(:,:,:,:,:,j) = dist5d_prof(:,:,:,:,:,j)/(dvol*u1)
         END DO
      END IF

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, VOL_ARR, ns_prof1*ns_prof2*ns_prof3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_LOCAL, ierr_mpi)
      CALL MPI_COMM_FREE(MPI_COMM_LOCAL,ierr_mpi)
      CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'beams3d_distnorm',ierr_mpi)
#endif

      RETURN
!-----------------------------------------------------------------------
!     END SUBROUTINE
!-----------------------------------------------------------------------
      END SUBROUTINE beams3d_distnorm