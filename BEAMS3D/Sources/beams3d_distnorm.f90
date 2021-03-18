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
      USE beams3d_lines, ONLY: ns_prof1, ns_prof2, ns_prof3, ns_prof4, &
                               ns_prof4, ns_prof5, dist5d_prof, &
                               partvmax
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
      INTEGER ::  s, i, j, k, nvol
      REAL(rprec) :: s1, s2, u1, u2, p1, ds, du, dp, area, dvol
      REAL(rprec), DIMENSION(4) :: rt,zt,pt

!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------

      ! Do physical volume elements
      nvol = (ns_prof1)*(ns_prof2)*(ns_prof3)
      ds   = 1.0/REAL(ns_prof1)
      du   = pi2/REAL(ns_prof2)
      dp   = pi2/REAL(ns_prof3)
      DO s = 1, nvol
         i = MOD(s-1,(ns_prof1))+1
         j = MOD(s-1,(ns_prof1)*(ns_prof2))
         j = FLOOR(REAL(j) / REAL(ns_prof1))+1
         k = CEILING(REAL(s) / REAL(ns_prof1*ns_prof2))
         s1 = (i-1)*ds
         s2 = s1+ds
         u1 = (j-1)*du
         u2 = u1+du
         p1 = (k-1)*dp*0.5
         rt = 0; zt = 0; pt = 0
         ! Assume parallel piped
         !CALL beams3d_suv2rzp(s1,u1,p1,rt(1),zt(1),pt(1))
         !CALL beams3d_suv2rzp(s1,u2,p1,rt(2),zt(2),pt(2))
         !CALL beams3d_suv2rzp(s2,u1,p1,rt(3),zt(3),pt(3))
         !CALL beams3d_suv2rzp(s2,u2,p1,rt(4),zt(4),pt(4))
         area = 0.5*(rt(1)*zt(2)-zt(1)*rt(2) + &
                      rt(2)*zt(3)-zt(2)*rt(3) + &
                      rt(3)*zt(4)-zt(3)*zt(4) + &
                      rt(4)*zt(1)-zt(4)*rt(1))
         !area = ds*0.5.*(s1+s2)*du
         !dvol = area*sum(rt)*dp/4
         area = ds*0.5*(s1+s2)*du
         dvol = area*(s1+s2)*0.5*dp
         dist5d_prof(:,i,j,k,:,:) = dist5d_prof(:,i,j,k,:,:)/dvol
         WRITE(328,*) i,j,k,dvol
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


      RETURN
!-----------------------------------------------------------------------
!     END SUBROUTINE
!-----------------------------------------------------------------------
      END SUBROUTINE beams3d_distnorm