!>  \brief Module contained subroutines for computing the real-space contravariant (full mesh) "displacements"
!!  jacobian*velocity*delta_t from the contravariant velocity harmonics obtained from solving
!!  the Fourier components of the MHD force balance as part of the SIESTA project.
!!  \author S. P. Hirshman and R. Sanchez
!!  \date Aug 18, 2006
      MODULE siesta_displacement
      USE stel_kinds
      USE stel_constants        
      USE fourier 
      USE quantities
      USE shared_data, ONLY: xc, l_push_s, l_push_u, l_push_v,        &
                             l_push_edge, l_natural
      USE timer_mod, ONLY: time_update_upperv
      USE descriptor_mod, ONLY: iam  

      CONTAINS
#if defined(SKS)
!>    \brief Parallel subroutine for initializing updated displacement components in real space
      SUBROUTINE update_upperv_par
      USE nscalingtools, ONLY: startglobrow, endglobrow
      USE siesta_state, ONLY: clear_field_perts
      IMPLICIT NONE
      INCLUDE 'mpif.h'
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: m0=0,m1=1,m2=2,ifull=0,pcos=0,psin=1
      INTEGER ::  js, nsmin, nsmax, ns_span, istat
      REAL(dp) :: ton, toff, r0
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                       &
                gvsupsmncf, gvsupumnsf, gvsupvmnsf
!-----------------------------------------------
      CALL second0(ton)
      nsmin=MAX(1,startglobrow-2); nsmax=MIN(endglobrow+2,ns)
      ns_span = nsmax-nsmin+1

      CALL Clear_Field_Perts

!NOTE: although sqrt(g) = 0 at origin, jvsups and jvsupv are defined with
!      jacobf which is finite (=jacobh(2)) at the origin to allow asymptotically
!      correct variations in the s and v components
       ALLOCATE (gvsupsmncf(0:mpol,-ntor:ntor,nsmin:nsmax),             &
                 gvsupumnsf(0:mpol,-ntor:ntor,nsmin:nsmax),             &
                 gvsupvmnsf(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
       IF (istat .ne. 0) STOP 'Allocation error # 1in UPDATE_UPPERV'

!Apply scaling in real space (keep xc components unscaled)
       gvsupsmncf = jvsupsmncf(:,:,nsmin:nsmax)
       gvsupumnsf = jvsupumnsf(:,:,nsmin:nsmax)
       gvsupvmnsf = jvsupvmnsf(:,:,nsmin:nsmax)

!Origin boundary conditions (evolve m=1 F_u,F_s, m=0 F_v)
       IF (nsmin .EQ. 1) THEN
       gvsupsmncf(m0,:,1) = 0;  gvsupsmncf(m2:,:,1) = 0
       IF (.NOT.l_push_s) gvsupsmncf(m1,:,1) = 0

       gvsupumnsf(m0,:,1) = 0;  gvsupumnsf(m2:,:,1) = 0
       IF (.NOT.l_push_u) gvsupumnsf(m1,:,1) = 0

       gvsupvmnsf(m1:,:,1) = 0
       IF (.NOT.l_push_v) gvsupvmnsf(m0,:,1) = 0

       IF (.NOT. l_natural) THEN
          r0 = hs_i/2
          gvsupsmncf(m1,:,1) = -r0*gvsupumnsf(m1,:,1)
       END IF
       END IF   !nsmin == 1

!m=0, n<0 (need so Hessian catch-all loop gets this...)
       gvsupsmncf(m0,-ntor:-1,nsmin:nsmax) = 0       
       gvsupumnsf(m0,-ntor:0,nsmin:nsmax) = 0       
       gvsupvmnsf(m0,-ntor:0,nsmin:nsmax) = 0       

!Edge boundary conditions
       IF (nsmax .EQ. ns) THEN
       gvsupsmncf(:,:,ns) = 0
       IF (.NOT.l_push_edge) THEN
          gvsupumnsf(:,:,ns) = 0
          gvsupvmnsf(:,:,ns) = 0
       END IF
       END IF  !!nsmax == ns

       CALL toijsp_par(gvsupsmncf, jvsupsijcf(:,:,nsmin:nsmax), 0, 0, pcos, ifull, 1, ns_span)    ! Calculate covariant velocities in real space (full mesh)
       CALL toijsp_par(gvsupumnsf, jvsupuijsf(:,:,nsmin:nsmax), 0, 0, psin, ifull, 1, ns_span)    ! NOTE: s-component is COS, remaining are SIN
       CALL toijsp_par(gvsupvmnsf, jvsupvijsf(:,:,nsmin:nsmax), 0, 0, psin, ifull, 1, ns_span)

       DEALLOCATE (gvsupsmncf, gvsupumnsf, gvsupvmnsf,stat=istat)

       CALL second0(toff)
       time_update_upperv = time_update_upperv + (toff-ton)

      END SUBROUTINE update_upperv_par
#endif

!>    \brief Serial subroutine for initializing updated displacement components in real space
      SUBROUTINE update_upperv
       IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER, PARAMETER :: m0=0,m1=1,m2=2,ifull=0,pcos=0,psin=1
       REAL(dp) :: ton, toff, r0
       REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                       &
                           gvsupsmncf, gvsupumnsf, gvsupvmnsf
       INTEGER :: i, j, k, nsmin, nsmax, istat
!-----------------------------------------------
!NOTE: although sqrt(g) = 0 at origin, jvsups and jvsupv are defined with
!      jacobf which is finite (=jacobh(2)) at the origin to allow asymptotically
!      correct variations in the s and v components
       CALL second0(ton)

       ALLOCATE (gvsupsmncf(0:mpol,-ntor:ntor,ns),                      &
                 gvsupumnsf(0:mpol,-ntor:ntor,ns),                      &
                 gvsupvmnsf(0:mpol,-ntor:ntor,ns), stat=istat)
       IF (istat .ne. 0) STOP 'Allocation error # 1in UPDATE_UPPERV'

!Apply scaling in real space (keep xc components unscaled)
       gvsupsmncf = jvsupsmncf
       gvsupumnsf = jvsupumnsf
       gvsupvmnsf = jvsupvmnsf

!Origin boundary conditions (evolve m=1 F_u,F_s, m=0 F_v)
       gvsupsmncf(m0,:,1) = 0;  gvsupsmncf(m2:,:,1) = 0
       IF (.NOT.l_push_s) gvsupsmncf(m1,:,1) = 0

       gvsupumnsf(m0,:,1) = 0;  gvsupumnsf(m2:,:,1) = 0
       IF (.NOT.l_push_u) gvsupumnsf(m1,:,1) = 0

       gvsupvmnsf(m1:,:,1) = 0
       IF (.NOT.l_push_v) gvsupvmnsf(m0,:,1) = 0

       IF (.NOT. l_natural) THEN
          r0 = hs_i/2
          gvsupsmncf(m1,:,1) = -r0*gvsupumnsf(m1,:,1)
       END IF

!m=0, n<0 (need so Hessian catch-all loop gets this...)
       gvsupsmncf(m0,-ntor:-1,:) = 0       
       gvsupumnsf(m0,-ntor:0,:) = 0       
       gvsupvmnsf(m0,-ntor:0,:) = 0       

!Edge boundary conditions
       gvsupsmncf(:,:,ns) = 0
       IF (.NOT.l_push_edge) THEN
          gvsupumnsf(:,:,ns) = 0
          gvsupvmnsf(:,:,ns) = 0
       END IF

       CALL toijsp(gvsupsmncf, jvsupsijcf, 0, 0, pcos, ifull)    ! Calculate upper velocities in real space (full mesh)
       CALL toijsp(gvsupumnsf, jvsupuijsf, 0, 0, psin, ifull)    ! NOTE: s-component is COS, remaining are SIN
       CALL toijsp(gvsupvmnsf, jvsupvijsf, 0, 0, psin, ifull)

       DEALLOCATE (gvsupsmncf, gvsupumnsf, gvsupvmnsf)

       CALL second0(toff)
       time_update_upperv = time_update_upperv + (toff-ton)

      END SUBROUTINE update_upperv

      END MODULE siesta_displacement
