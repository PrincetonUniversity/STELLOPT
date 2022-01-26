!-----------------------------------------------------------------------
!     Module:        beams3d_distnorm
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          03/07/2021
!     Description:   This subroutine applies the volume normalization
!                    to the distribution fucntion.
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_distnorm
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE beams3d_runtime, ONLY: nbeams, pi2
      USE beams3d_grid, ONLY: raxis, phiaxis, zaxis, S4D, U4D, &
                              nr, nphi, nz
      USE beams3d_lines, ONLY: ns_prof1, ns_prof2, ns_prof3, ns_prof4, &
                               ns_prof4, ns_prof5, dist5d_prof, &
                               partvmax
      USE beams3d_physics_mod, ONLY: beams3d_suv2rzp
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
      INTEGER ::  s, i, j, k, nvol, l, m, n
      INTEGER, DIMENSION(2) :: minln
      REAL(rprec) :: s1, s2, u1, u2, p1, ds, du, dp, area, dvol, a, b, c, d, f1, f2, f3
      REAL(rprec), DIMENSION(4) :: rt,zt,pt
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: targ
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: RHO3D

!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------

      ! Do physical volume elements
      nvol = (ns_prof1)*(ns_prof2)*(ns_prof3)
      ds   = 1.0/REAL(ns_prof1) !distribution defined on centers (half grid)
      du   = pi2/REAL(ns_prof2)
      dp   = pi2/REAL(ns_prof3)
      ALLOCATE(targ(nr,nz))
      ALLOCATE(RHO3D(nr,nphi,nz))
      RHO3D = S4D(1,:,:,:)
      WHERE (RHO3D>1) RHO3D = 2;
      RHO3D = sqrt(RHO3D) ! in RHO
      DO s = 1, nvol
         i = MOD(s-1,(ns_prof1))+1
         j = MOD(s-1,(ns_prof1)*(ns_prof2))
         j = FLOOR(REAL(j) / REAL(ns_prof1))+1
         k = CEILING(REAL(s) / REAL(ns_prof1*ns_prof2))
         s1 = (i-1)*ds
         s2 = s1+ds
         u1 = (j-1)*du
         u2 = u1+du
         p1 = MOD((k-1)*dp*0.5,phiaxis(nphi))
         ! Find helper value of grids
         m = MIN(MAX(COUNT(phiaxis < p1),1),nphi-1)
         rt = 5.25+COS((u1+u2)/2)*.5 !dirty estimation of R
         zt = SIN((u1+u2)/2.)*.6 
         pt = 0; !Initial conditions for successful lookup?
         ! Assume a parallel piped
         ! s1 u1
         targ = (0.5*sqrt(RHO3D(:,m,:) + RHO3D(:,m+1,:))-s1)**2 &
              + (0.5*(U4D(1,:,m,:) + U4D(1,:,m+1,:))-u1)**2
         minln = MINLOC(targ)
         l = MIN(MAX(minln(1),2),nr-1); n = MIN(MAX(minln(2),2),nz-1)
         ! Linear interpolation
         f1 = sqrt(targ(l,n)); f2 = sqrt(targ(l+1,n)); f3 = sqrt(targ(l,n+1))
         b = f1; a = f2-f1 
         c = f1; d = f3-f1
         a = -b/a
         b = -d/c
         !PRINT *,l,n,a,b
         !PRINT *,l,n,targ(l,n),targ(l+1,n),targ(l,n+1),targ(l+1,n+1)
         !rt(1) = raxis(l)*(1-a)+raxis(l+1)*a; zt(1) = zaxis(n)*(1-b)+zaxis(n+1)*b
         ! s1 u2
         !targ = (0.5*sqrt(RHO3D(:,m,:) + RHO3D(:,m+1,:))-s1)**2 &
         !     + (0.5*(U4D(1,:,m,:) + U4D(1,:,m+1,:))-u2)**2
         !minln = MINLOC(targ)
         !l = MIN(MAX(minln(1),1),nr-1); n = MIN(MAX(minln(2),1),nz-1)
         !a = sqrt(targ(l,n)); b = sqrt(targ(l,n+1))
         !rt(2) = raxis(l)*(1-a)+raxis(l+1)*a; zt(2) = zaxis(n)*(1-b)+zaxis(n+1)*b
         ! s2 u1
         !targ = (0.5*sqrt(RHO3D(:,m,:) + RHO3D(:,m+1,:))-s2)**2 &
         !     + (0.5*(U4D(1,:,m,:) + U4D(1,:,m+1,:))-u1)**2
         !minln = MINLOC(targ)
         !l = MIN(MAX(minln(1),1),nr-1); n = MIN(MAX(minln(2),1),nz-1)
         !a = sqrt(targ(l,n)); b = sqrt(targ(l,n+1))
         !rt(3) = raxis(l)*(1-a)+raxis(l+1)*a; zt(3) = zaxis(n)*(1-b)+zaxis(n+1)*b
         ! s2 u2
         !targ = (0.5*sqrt(RHO3D(:,m,:) + RHO3D(:,m+1,:))-s2)**2 &
         !     + (0.5*(U4D(1,:,m,:) + U4D(1,:,m+1,:))-u2)**2
         !minln = MINLOC(targ)
         !l = MIN(MAX(minln(1),1),nr-1); n = MIN(MAX(minln(2),1),nz-1)
         !a = sqrt(targ(l,n)); b = sqrt(targ(l,n+1))
         !rt(4) = raxis(l)*(1-a)+raxis(l+1)*a; zt(4) = zaxis(n)*(1-b)+zaxis(n+1)*b
         ! Assume parallel piped
         CALL beams3d_suv2rzp(s1,u1,p1,rt(1),zt(1),pt(1))
         CALL beams3d_suv2rzp(s1,u2,p1,rt(2),zt(2),pt(2))
         CALL beams3d_suv2rzp(s2,u1,p1,rt(3),zt(3),pt(3))
         CALL beams3d_suv2rzp(s2,u2,p1,rt(4),zt(4),pt(4))
         area = 0.5*abs(rt(1)*zt(2)-zt(1)*rt(2) + &
                      rt(2)*zt(3)-zt(2)*rt(3) + &
                      rt(3)*zt(4)-zt(3)*zt(4) + &
                      rt(4)*zt(1)-zt(4)*rt(1))
         !area = 0.5*(s1+s2)*du*ds ! Missing aminor
         WRITE(327,*) m,s1,s2,u1,u2,rt,zt,area
         CALL FLUSH(327)
         dvol = area*sum(rt)*dp/4
         dist5d_prof(:,i,j,k,:,:) = dist5d_prof(:,i,j,k,:,:)/dvol
         WRITE(328,*) i,j,k,dvol,rt,zt
         CALL FLUSH(328)
      END DO

      ! Do phase space volume elements
      ds = 2*partvmax/ns_prof4
      du = partvmax/ns_prof5
      dvol = pi2*ds*du
      nvol = ns_prof4*ns_prof5
      DO j = 1, ns_prof5
         u1 = (j-0.5)*du
         dist5d_prof(:,:,:,:,:,j) = dist5d_prof(:,:,:,:,:,j)/(dvol*u1)
         WRITE(329,*) j,dvol*u1
         CALL FLUSH(329)
      END DO

      DEALLOCATE(targ)
      DEALLOCATE(RHO3D)


      RETURN
!-----------------------------------------------------------------------
!     END SUBROUTINE
!-----------------------------------------------------------------------
      END SUBROUTINE beams3d_distnorm