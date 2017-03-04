!>  \brief Module contained subroutines for updating (from t to t+delta_t) the magnetic field and pressure as part of the SIESTA project.
!!  Stores updated values of PMNCF, BSUBMNF, BSUPMNF, and BSUPIJF in Quantities Module
!!  \author S. P. Hirshman and R. Sanchez
!!  \date Sep 15, 2006
      MODULE siesta_state

      USE stel_kinds
      USE quantities
      USE diagnostics_mod
      USE shared_data, ONLY: ste, bs0, bu0, bsbu_ratio
      USE timer_mod, ONLY: time_update_state

      CONTAINS

!>  \brief Parallel subroutine for updating the SIESTA state
      SUBROUTINE update_state_par (lprint, fsq_total, ftol)
       USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR,         &
                             mnspcounts,mnspdisps
       IMPLICIT NONE
#if defined(MPI_OPT)
       INCLUDE 'mpif.h'
#endif
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       REAL(rprec), INTENT(in)  :: fsq_total, ftol
       LOGICAL, INTENT(in)      :: lprint
#if defined(SKS)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       LOGICAL, PARAMETER :: lpar=.FALSE.
       REAL(rprec) :: ton, toff
       REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: tmps, tmpu, tmpv, tmpp
       INTEGER :: nsmin, nsmax, nloc, js, istat
!-----------------------------------------------
      CALL second0(ton)
      nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow,ns)

      ALLOCATE (tmps(0:mpol,-ntor:ntor,ns), tmpu(0:mpol,-ntor:ntor,ns), &
                tmpv(0:mpol,-ntor:ntor,ns), tmpp(0:mpol,-ntor:ntor,ns), &
                stat=js)
      IF (js .NE. 0) STOP 'Allocation error in update_state_par'

      nloc = (nsmax-nsmin+1)*SIZE(tmps,1)*SIZE(tmps,2)
      
      CALL MPI_ALLGATHERV(djbsupsmnsh(:,:,nsmin:nsmax),nloc,MPI_REAL8,  &
                         tmps,mnspcounts,mnspdisps,MPI_REAL8,MPI_COMM_WORLD,MPI_ERR)
      CALL MPI_ALLGATHERV(djbsupumnch(:,:,nsmin:nsmax),nloc,MPI_REAL8,  &
                         tmpu,mnspcounts,mnspdisps,MPI_REAL8,MPI_COMM_WORLD,MPI_ERR)
      CALL MPI_ALLGATHERV(djbsupvmnch(:,:,nsmin:nsmax),nloc,MPI_REAL8,  &
                         tmpv,mnspcounts,mnspdisps,MPI_REAL8,MPI_COMM_WORLD,MPI_ERR)
      CALL MPI_ALLGATHERV(djpmnch(:,:,nsmin:nsmax),nloc,MPI_REAL8,      &
                         tmpp,mnspcounts,mnspdisps,MPI_REAL8,MPI_COMM_WORLD,MPI_ERR)

      CALL Init_Allocate_Arrays (lpar)

      djbsupsmnsh = tmps
      djbsupumnch = tmpu
      djbsupvmnch = tmpv
      djpmnch     = tmpp

      DEALLOCATE (tmps, tmpu, tmpv, tmpp)
      CALL second0(toff)
      time_update_state = time_update_state + (toff-ton)

      CALL update_state(.FALSE., fsq_total, ftol)


      IF (.NOT.lprint) THEN
         CALL Init_Allocate_Arrays (.TRUE.)          !Clears parallel perts
         RETURN
      END IF

      CALL second0(ton)

      CALL init_state_par (.TRUE.)                          
!   Compute co-variant and contravariant components of current (times jacobian) - need for divJ, BdotJ
      lcurr_init=.TRUE.
      CALL divb_par
      CALL divj_par
      CALL bgradp_par
      CALL tflux
      CALL bdotj_par      
      lcurr_init=.FALSE.                             

      toroidal_flux = toroidal_flux - toroidal_flux0
      IF (iam .EQ. 0) THEN
      DO ISTAT = 6, 33, 27
         WRITE (ISTAT , 110) ste(1), ste(2), ste(3), ste(4),            &
                             divb_rms, toroidal_flux, wp/wb,            &
                             bgradp_rms, max_bgradp, min_bgradp,        &
                             bdotj_rms, bdotj2_rms, divj_rms
      END DO
      ENDIF

 110   FORMAT(' SPECTRAL TRUNC ERROR - p: ',1pe11.3,' B_s: ',1pe11.3,   &
             ' B_u: ',1pe11.3,' B_v: ',1pe11.3,/,                       &
             ' DIV-B (rms): ',1pe11.3, ' DEL_TFLUX: ',1pe11.3,/,        &
             ' <BETA>: ', 1pe11.3,' B.GRAD-P (rms): ', 1pe11.3,         &
             ' B.GRAD-P (max): ', 1pe11.3,' B.GRAD-P (min): ',1pe11.3,  & 
           /,' (J*B)/|JxB| (rms): ', 1pe11.3,                           &
             ' (J_par)/|J_tot| (rms): ', 1pe11.3,                       &
             '   DIV-J (rms): ', 1pe11.3)

      CALL second0(toff)
      time_update_state = time_update_state + (toff-ton)
#endif       
      END SUBROUTINE update_state_par

!>  \brief Serial subroutine for updating the SIESTA state
      SUBROUTINE update_state (lprint, fsq_total, ftol)
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       REAL(rprec), INTENT(in)  :: fsq_total, ftol
       LOGICAL, INTENT(in)      :: lprint
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER, PARAMETER :: m1=1
       INTEGER     :: istat, js
       REAL(rprec) :: ton, toff, r0, xFSQ
       REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: work4
!-----------------------------------------------
       CALL second0(ton)

!      Update BSUP?MN?H*jacobh, pressure-h*jacobh harmonics based on vsup?mn?f perturbations
!  
       jbsupsmnsh = jbsupsmnsh + djbsupsmnsh                ! update B^s harmonics  (half mesh)     
       jbsupumnch = jbsupumnch + djbsupumnch                ! update B^u harmonics  (half mesh)  
       jbsupvmnch = jbsupvmnch + djbsupvmnch                ! update B^v harmonics  (half mesh) 
       jpmnch     = jpmnch     + djpmnch                    ! update pressure harmonics  (half mesh) 

      
!      Reset perturbations
       CALL Clear_Field_Perts

       CALL second0(toff)
       time_update_state = time_update_state + (toff-ton)

!      fsq_total == 0 when updating from add_resistivity
       IF (fsq_total.NE.zero .AND. fsq_total.LT.10*ftol) THEN
          CALL write_profiles(fsq_total)                    ! SPH: write pmn, jbsupsmn profiles
       END IF

      r0 = hs_i/2

      bs0(1) = SQRT(SUM((jbsupsmnsh(m1,:,2) - r0*jbsupumnch(m1,:,2))**2))/r0
      bu0(1) = SQRT(SUM((jbsupsmnsh(m1,:,2) + r0*jbsupumnch(m1,:,2))**2))/r0
      IF (bu0(1) .GT. 1.E-10_dp) THEN
         bsbu_ratio = bs0(1)/bu0(1)
         DO js = 2, MIN(ns,6)
            bs0(js) = SQRT(SUM(jbsupsmnsh(m1,:,js)**2))/r0
            bu0(js) = SQRT(SUM(jbsupumnch(m1,:,js)**2))
         END DO
      ELSE 
         bs0=0; bu0=0; bsbu_ratio = 1
      END IF

      IF (.NOT.lprint) RETURN

      CALL second0(ton)

      lcurr_init=.TRUE.                             
      CALL init_state(.TRUE.)                              ! Compute co/contra variant currents for divJ, bdotj

      CALL divb
      CALL divj                                            ! (RS) Added, 09/10/07  = compute rms of divergence of J
      CALL bgradp
      CALL tflux
      CALL bdotj                                           ! (RS) Added, 08/30/07  = compute rms of parallel current
      CALL beta
      lcurr_init=.FALSE.                             

      toroidal_flux = toroidal_flux - toroidal_flux0
      IF (iam .EQ. 0) THEN
      DO ISTAT = 6, 33, 27
         WRITE (ISTAT , 110) ste(1), ste(2), ste(3), ste(4),            &
                             divb_rms, toroidal_flux, wp/wb,            &
                             bgradp_rms, max_bgradp, min_bgradp,        &
                             bdotj_rms, bdotj2_rms, divj_rms
      END DO
      ENDIF

 110   FORMAT(' SPECTRAL TRUNC ERROR - p: ',1pe11.3,' B_s: ',1pe11.3,   &
             ' B_u: ',1pe11.3,' B_v: ',1pe11.3,/,                       &
             ' DIV-B (rms): ',1pe11.3, ' DEL_TFLUX: ',1pe11.3,/,        &
             ' <BETA>: ', 1pe11.3,' B.GRAD-P (rms): ', 1pe11.3,         &
             ' B.GRAD-P (max): ', 1pe11.3,' B.GRAD-P (min): ',1pe11.3,  & 
           /,' (J*B)/|JxB| (rms): ', 1pe11.3,                           &
             ' (J_par)/|J_tot| (rms): ', 1pe11.3,                       &
             '   DIV-J (rms): ', 1pe11.3)
      CALL second0(toff)
      time_update_state = time_update_state + (toff-ton)
       
      END SUBROUTINE update_state


!>  \brief Subroutine for resetting the perturbations of the SIESTA state
      SUBROUTINE Clear_Field_Perts
      IMPLICIT NONE
!-----------------------------------------------
      djbsupsmnsh = 0
      djbsupumnch = 0
      djbsupvmnch = 0
      djpmnch     = 0

      END SUBROUTINE Clear_Field_Perts

      END MODULE siesta_state
