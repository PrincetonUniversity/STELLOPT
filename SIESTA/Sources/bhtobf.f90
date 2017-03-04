#if defined(SKS)
      SUBROUTINE bhalftobfull_par (bsupsijsh, bsupuijch, bsupvijch,     &
                                   bsupsijsf, bsupuijcf, bsupvijcf,     &
                                   nsmin, nsmax)
!     
!     WRITTEN 01-04-10 BY S. HIRSHMAN AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: Store BsupX on full mesh
!     NOTE:    On Entry, bsups,u_h at js=1 contain the "evolved" values of 
!              bsups,u_f at s=0 (js=1), as computed in UPDATE_BFIELD
!
       USE stel_kinds
       USE quantities, ONLY: ns, ntheta, nzeta, jacobf_p, jacobh_p

       USE descriptor_mod,ONLY: iam, nprocs
       USE nscalingtools, ONLY: startglobrow, endglobrow
       USE timer_mod, ONLY: bhtobf_time
       IMPLICIT NONE
       INCLUDE 'mpif.h'
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       REAL(dp), DIMENSION(ntheta,nzeta,nsmin:nsmax), INTENT(IN)    ::  &
           bsupsijsh, bsupvijch
       REAL(dp), DIMENSION(ntheta,nzeta,nsmin:nsmax), INTENT(INOUT) ::  &
           bsupuijch 
       REAL(dp), DIMENSION(ntheta,nzeta,nsmin:nsmax), INTENT(OUT)   ::  &
           bsupsijsf, bsupuijcf, bsupvijcf
       INTEGER, INTENT(in) :: nsmin, nsmax
!-----------------------------------------------
       REAL(dp) :: skston, skstoff
       INTEGER  :: ns_span, nblock
!-----------------------------------------------
       CALL second0(skston)
       ns_span=nsmax-nsmin+1
       nblock = SIZE(bsupsijsh,1)*SIZE(bsupsijsh,2)

!Interpolate b^X finite at the origin
!On exit, bsupXF = [nsmin, nsmax-1]
!Radial (s)
       CALL to_full_mesh_par(bsupsijsh, bsupsijsf, nblock, nsmin, nsmax)

!Poloidal(u): on input, bsupuijch = bsupuijch*jacobh
       CALL to_full_mesh_par(bsupuijch, bsupuijcf, nblock, nsmin, nsmax)

!Toroidal (v)
       CALL to_full_mesh_par(bsupvijch, bsupvijcf, nblock, nsmin, nsmax)

       IF (nsmax .EQ. ns) bsupsijsf(:,:,ns) = 0

       bsupuijcf = bsupuijcf/jacobf_p(:,:,nsmin:nsmax)
       bsupuijch = bsupuijch/jacobh_p(:,:,nsmin:nsmax)

       CALL second0(skstoff)
       bhtobf_time=bhtobf_time+(skstoff-skston)

      END SUBROUTINE bhalftobfull_par
#endif

      SUBROUTINE bhalftobfull (bsupsijsh, bsupuijch, bsupvijch,         &
                               bsupsijsf, bsupuijcf, bsupvijcf)
!     
!     WRITTEN 01-04-10 BY S. HIRSHMAN AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: Store BsupX on full mesh
!     NOTE:    On Entry, bsups,u_h at js=1 contain the "evolved" values of 
!              bsups,u_f at s=0 (js=1), as computed in UPDATE_BFIELD
!
#if defined(SKS)
       USE timer_mod, ONLY: bhtobf_time
#endif
       USE stel_kinds
       USE quantities, ONLY: ns, ntheta, nzeta, jacobf, jacobh
       USE descriptor_mod,ONLY: DIAGONALDONE, INHESSIAN 
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       REAL(dp), DIMENSION(ntheta,nzeta,ns), INTENT(in)    ::           &
           bsupsijsh, bsupvijch
       REAL(dp), DIMENSION(ntheta,nzeta,ns), INTENT(inout) ::           &
           bsupuijch
       REAL(dp), DIMENSION(ntheta,nzeta,ns), INTENT(out)   ::           &
           bsupsijsf, bsupuijcf, bsupvijcf
       REAL(dp) :: skston, skstoff
       INTEGER  :: nblock
!-----------------------------------------------
       CALL second0(skston)
       nblock = SIZE(bsupsijsh,1)*SIZE(bsupsijsh,2)

!Interpolate b^X which is finite at the origin
!Radial (s)
       CALL to_full_mesh(bsupsijsh, bsupsijsf, nblock)
       bsupsijsf(:,:,1) = bsupsijsh(:,:,1)
       bsupsijsf(:,:,ns) = 0

!Poloidal(u): on input, bsupuijch = bsupuijch*jacobh
       CALL to_full_mesh(bsupuijch, bsupuijcf, nblock)
       bsupuijcf(:,:,1) = bsupuijch(:,:,1)
       bsupuijcf(:,:,ns) = bsupuijch(:,:,ns)        !1pt "extrap" => tridi
       bsupuijcf = bsupuijcf/jacobf
       bsupuijch = bsupuijch/jacobh

!Toroidal (v)
       CALL to_full_mesh(bsupvijch, bsupvijcf, nblock)
       bsupvijcf(:,:,1) = bsupvijch(:,:,1)
       bsupvijcf(:,:,ns) = bsupvijch(:,:,ns)        !1pt "extrap" => tridi
#if defined(SKS)
       CALL second0(skstoff)
       IF(DIAGONALDONE.AND.INHESSIAN) bhtobf_time=bhtobf_time+(skstoff-skston)
#endif
       END SUBROUTINE bhalftobfull
