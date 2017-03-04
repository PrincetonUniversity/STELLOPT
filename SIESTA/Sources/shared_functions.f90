
!>  \brief Module contained subroutines and functions updating MHD forces and Wmhd
!!  \author S. P. Hirshman
!!  \date July, 2014
      MODULE shared_functions
      USE shared_data
      USE descriptor_mod, ONLY: INHESSIAN, nprocs,                       &
                                in_hess_nfunct, out_hess_nfunct
      USE island_params, ns=>ns_i, hs=>hs_i
      USE hessian, ONLY: apply_precond, l_Compute_Hessian
      USE quantities, ONLY: fsubsmncf, fsubumnsf, fsubvmnsf, fsupsmncf,  &
                            fsupumnsf, fsupvmnsf, wb, wp
      USE timer_mod, ONLY: time_init_state, time_funci

      CONTAINS

!>  \brief Parallel routine to compute forces for perturbed state.
!!  \author S. P. Hirshman, R. Sanchez, S. K. Seal
      SUBROUTINE funct_island_par
!-----------------------------------------------
!     PARALLEL (in N rows) VERSION
!-----------------------------------------------
#if defined(SKS)      
      USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR
      USE siesta_bfield, ONLY: update_bfield_par
      USE siesta_pressure, ONLY: update_pres_par
      USE siesta_force, ONLY: update_force_par
      USE siesta_displacement, ONLY: update_upperv_par
      IMPLICIT NONE
      INCLUDE 'mpif.h'
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp) :: ton, toff, skston, skstoff
      INTEGER  :: nsmin, nsmax, js
      REAL(dp) :: stmp(3), rtmp(3)
      LOGICAL  :: ltype
!-----------------------------------------------
      CALL second0(ton)
      nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow,ns)

      IF (INHESSIAN) THEN
         in_hess_nfunct=in_hess_nfunct+(nsmax-nsmin+1)
      ELSE
         out_hess_nfunct=out_hess_nfunct+(nsmax-nsmin+1)*nprocs
      END IF

      IF (l_init_state) THEN
        CALL second0(skston)
        CALL init_state_par(.FALSE.)
        CALL second0(skstoff)
        time_init_state=time_init_state+(skstoff-skston)
        l_init_state = .FALSE.
      END IF

      CALL update_upperv_par
      CALL update_bfield_par (.FALSE.)
      CALL update_pres_par (delta_t)
      CALL update_force_par

      IF (gamma .eq. 1._dp) STOP 'SIESTA REQUIRES gamma != 1'
      wtotal = wb + wp/(gamma-1)

      IF (wtotal0 == -1) wtotal0 = wtotal
      gc = gnorm_i*gc

!     CALLED FROM GMRES: ADD BACK INITIAL FORCE gc0 or STORE UNPRECONDITIONED FORCES
      IF (ALLOCATED(gc0) .AND. l_getfsq .AND. l_linearize .AND.         &
         .NOT.l_Compute_Hessian) THEN
         gc = gc+gc0
      END IF

      IF (l_getfsq .AND. INHESSIAN) STOP 'l_getfsq must be set to FALSE in Hessian'

!     COMPUTE UN-PRECONDITIONED FORCE RESIDUES
      IF (l_getfsq) THEN
         stmp(1) = SUM(fsubsmncf(:,:,nsmin:nsmax)*fsubsmncf(:,:,nsmin:nsmax))
         stmp(2) = SUM(fsubumnsf(:,:,nsmin:nsmax)*fsubumnsf(:,:,nsmin:nsmax))
         stmp(3) = SUM(fsubvmnsf(:,:,nsmin:nsmax)*fsubvmnsf(:,:,nsmin:nsmax))
         CALL MPI_ALLREDUCE(stmp,rtmp,3,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,MPI_ERR)
         fsqvs=hs*rtmp(1); fsqvu=hs*rtmp(2); fsqvv=hs*rtmp(3)
         CALL toupper_forces_par
!     VOLUME-AVERAGE |F|**2
         stmp(1) = 0
         DO js = nsmin,nsmax
            stmp(1) = stmp(1) + vp_f(js)*SUM(                           &
                          fsupsmncf(:,:,js)*fsubsmncf(:,:,js)           & 
                        + fsupumnsf(:,:,js)*fsubumnsf(:,:,js)           &
                        + fsupvmnsf(:,:,js)*fsubvmnsf(:,:,js))
         END DO
         stmp(2) = SUM(vp_f(nsmin:nsmax))
         CALL MPI_ALLREDUCE(stmp,rtmp,2,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,MPI_ERR)
         fsq_total = rtmp(1)/rtmp(2);  CALL truncate(fsq_total,10)
         fsq_total1 = fsq_total
      END IF

!      IF (l_Gmres .AND. l_ApplyPrecon) STOP 'l_ApplyPrecon must be FALSE in funct_island_par'
!
!     GMRES handles preconditioning itself: do NOT apply
!     Do not precondition when Hessian is being computed
      IF (l_ApplyPrecon) CALL apply_precond(gc)

      CALL second0(toff)
      time_funci = time_funci+(toff-ton)
#endif
      END SUBROUTINE funct_island_par

!>  \brief Computes forces for perturbed state
!!  \author S. P. Hirshman, R. Sanchez
      SUBROUTINE funct_island
      USE siesta_bfield, ONLY: update_bfield
      USE siesta_pressure, ONLY: update_pres
      USE siesta_force, ONLY: update_force
      USE siesta_displacement, ONLY: update_upperv
!-----------------------------------------------
!     SERIAL (in N rows) VERSION
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  :: js
      REAL(dp) :: skston, skstoff, ton, toff
!-----------------------------------------------
      CALL second0(ton)
      skston=ton

      IF (INHESSIAN) THEN
         in_hess_nfunct=in_hess_nfunct+ns
      ELSE
         out_hess_nfunct=out_hess_nfunct+ns
      ENDIF

      IF (l_init_state) THEN
        CALL second0(skston)
        CALL init_state(.FALSE.)
        CALL second0(skstoff)
        time_init_state=time_init_state+(skstoff-skston)
        l_init_state = .FALSE.
      END IF

      CALL update_upperv
      CALL update_bfield (.FALSE.)
      CALL update_pres (delta_t)
      CALL update_force

      IF (gamma .eq. 1._dp) STOP 'SIESTA REQUIRES gamma != 1'
      wtotal = wb + wp/(gamma-1)
 
      IF (wtotal0 == -1) wtotal0 = wtotal
      gc = gnorm_i*gc

!     CALLED FROM GMRES: ADD BACK INITIAL FORCE gc0 or STORE UNPRECONDITIONED FORCES
      IF (ALLOCATED(gc0) .AND. l_getfsq) THEN
         IF (l_linearize .AND. .not.l_Compute_Hessian) THEN
            gc = gc+gc0
         END IF
      END IF

!     COMPUTE UN-PRECONDITIONED FORCE RESIDUES
      IF (l_getfsq) THEN
         fsqvs = hs*SUM(fsubsmncf*fsubsmncf)
         fsqvu = hs*SUM(fsubumnsf*fsubumnsf)
         fsqvv = hs*SUM(fsubvmnsf*fsubvmnsf)
         CALL toupper_forces
         DO js = 1, ns
            fsq_total = SUM((fsupsmncf(:,:,js)*fsubsmncf(:,:,js)        &
                      +      fsupumnsf(:,:,js)*fsubumnsf(:,:,js)        &
                      +      fsupvmnsf(:,:,js)*fsubvmnsf(:,:,js))*vp_f(js))
         END DO
         fsq_total = fsq_total/SUM(vp_f)
         CALL truncate(fsq_total, 10)
         fsq_total1 = fsq_total
      END IF
!
!     GMRES handles preconditioning itself: do NOT apply
!     Do not precondition when Hessian is being computed
      IF (l_ApplyPrecon) CALL apply_precond(gc)

      CALL second0(toff)
      time_funci = time_funci + (toff-ton)

      END SUBROUTINE funct_island


!>   \brief Returns perturbed MHD+KINETIC energy
!!   \author S. P. Hirshman, R. Sanchez
      FUNCTION getwmhd (p, lpar)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(in) :: p(neqs)           !<displacement vector
      LOGICAL,  INTENT(in) :: lpar              !<=T, use parallel funct_island
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp) :: getwmhd
!-----------------------------------------------
        
      xc = p
      l_init_state = .TRUE.
      l_linearize  = .FALSE.
      l_getwmhd    = .TRUE.
      l_getfsq     = .FALSE.
#if defined(SKS)
      IF (lpar) THEN
         CALL funct_island_par
      ELSE
#endif
         CALL funct_island
#if defined(SKS)
      END IF
#endif
      l_getwmhd    = .FALSE.

      getwmhd = wtotal

      RETURN

      END FUNCTION getwmhd 

      END MODULE shared_functions
