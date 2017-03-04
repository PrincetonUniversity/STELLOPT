!>
!!  \mainpage SIESTA
!!  \brief Program for computing 3D MHD equilibria including magnetic islands
!!  \author S. P. Hirshman, R. Sanchez, and C. R. Cook (Physics of Plasmas 18, 062504, 2011)
!!  \author S. K. Seal, K. S. Perumalla and S. P. Hirshman (Concurrency and Computation: Practice and Experience, 25(15), p 2207-2223, 2013)
!!  \version 2.3 (June 2014)
!!  \bug Please report any bugs to hirshmansp@ornl.gov
      PROGRAM SIESTA
!
!     READS siesta.jcf INPUT FILE
!
      USE evolution, ONLY: bsbu_ratio, jsju_ratio, bs0, bu0, niter,     &
                           init_evolution, converge_diagonal, converge_blocks
      USE shared_data, ONLY: nprecon, ngmres_type, xc, gc, fsq_total1
      USE metrics, ONLY: init_metric_elements, cleanup_metric_elements, &
                         dealloc_full_lower_metrics
      USE quantities, ONLY: init_quantities, dealloc_quantities,        &
                            djpmnch, fbdy
      USE hessian, ONLY: levmarq_param, mupar, asym_index, removedafile,&
                         dealloc_hessian
      USE perturbation, ONLY: nsin, mpolin, ntorin, wout_file, lrestart, &
                         ftol, lresistive, init_data, read_restart_file
      USE diagnostics_mod, ONLY: dealloc_diagnostics, write_profiles
      USE dumping_mod, ONLY: write_output
      USE prof_mod
      USE descriptor_mod
      USE timer_mod
      USE island_params, mpol=>mpol_i, ntor=>ntor_i, ns=>ns_i
#if defined(SKS)
      USE nscalingtools, ONLY: MyEnvVariables, SetOutputFile, PARSOLVER,   &
                            SKSDBG, TOFU, PARGMRES, PARFUNCTISL,        &
                            finalizeRemap, bcyclicMapping, initRemap,   &
                            OUTPUT_TIMINGS, GetTimes, MPI_ERR
#else
      USE nscalingtools, ONLY: MyEnvVariables, SetOutputFile, PARSOLVER,   &
                            SKSDBG, TOFU, PARGMRES, PARFUNCTISL
#endif
      IMPLICIT NONE
#if defined(MPI_OPT)
      INCLUDE 'mpif.h'
#endif
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER       :: i,j                       !<counter
      INTEGER       :: istat                     !<status variable
      REAL(rprec)   :: ton, skston, t1           !<timer counter (on)
      REAL(rprec)   :: toff, skstoff, t2         !<timer counter (off)
      REAL(rprec), ALLOCATABLE :: work2(:,:,:)   !<work array
      INTEGER       :: mblk_size
#if defined(SKS)
      REAL(rprec)   :: starttime, endtime, usedtime, rvar
      CHARACTER*10  :: envval
#endif      

!-----------------------------------------------
      CALL second0(ton)
      skston=ton
#if defined(MPI_OPT)
      CALL profinit()
!     -----------
!     setup blacs
!     -----------
      CALL blacs_pinfo(iam,nprocs)
      CALL blacs_setup(iam,nprocs)
      LSCALAPACK=.TRUE.
      INHESSIAN=.FALSE.
#if defined(SKS)
      CALL MyEnvVariables(iam)
      IF (PARSOLVER) LSCALAPACK=.FALSE.
      CALL sks_timers
      IF (SKSDBG) CALL SetOutputFile(iam,nprocs)
#endif
      IF (LSCALAPACK) THEN
      do nprow=int(sqrt(dble(nprocs)))+1,1,-1
         npcol = int( nprocs/nprow )
         if (nprow*npcol.eq.nprocs) exit
      enddo

      call blacs_get(0,0,icontxt)
      call blacs_gridinit(icontxt,'C',nprow,npcol)

      call blacs_get(0,0,icontxt_1xp)
      call blacs_gridinit(icontxt_1xp,'R',1,nprow*npcol)

      call blacs_get(0,0,icontxt_px1)
      call blacs_gridinit(icontxt_px1,'C',nprow*npcol,1)

      icontxt_global = icontxt_1xp

      call blacs_gridinfo(icontxt, nprow,npcol,myrow,mycol)
      isroot = (myrow.eq.0).and.(mycol.eq.0)
      if (isroot) then
         write(*,*) ' nprow,npcol ',nprow,npcol
         write(*,*) 'icontxt_global,icontxt,icontxt_1xp,icontxt_px1 ',  &
                     icontxt_global,icontxt,icontxt_1xp,icontxt_px1
      endif
      END IF   !LSCALAPACK=T
#else
      LSCALAPACK=.FALSE.
      iam = 0
      nprocs = 1
#endif

#if defined(SKS)
      CALL second0(t1)
#endif
      CALL init_data (nprecon)
#if defined(SKS)
      CALL second0(t2)
      init_data_time=init_data_time+(t2-t1)
#endif
!
!     INITIALIZE METRICS, FOURIER, AND B-FIELD
!
#if defined(SKS)
      CALL second0(t1)
#endif
      CALL init_metric_elements (nsin,mpolin,ntorin,wout_file)
#if defined(SKS)
      CALL second0(t2)
      init_metric_elements_time=init_metric_elements_time+(t2-t1)
#endif
      
#if defined(SKS)
      CALL second0(t1)
#endif
      CALL init_quantities
#if defined(SKS)
      CALL second0(t2)
      init_quantities_time=init_quantities_time+(t2-t1)
#endif
      
      IF (lrestart) CALL read_restart_file (nprecon)

!      ALLOCATE (work2(ns,ntheta,nzeta))
!      CALL ZeroHalfPoint(djpmnch, work2, 2, 2, 0)
!      DEALLOCATE(work2)

      CALL second0(toff)
      time_init = toff-ton

!     Initializes the EVOLUTION module variables and
!     converges the force residual using a diagonal preconditioner
!     Applies the external island perturbation if possible
#if defined(SKS)
      CALL second0(t1)
#endif
      CALL init_evolution
#if defined(SKS)
      CALL second0(t2)
      init_evolution_time=init_evolution_time+(t2-t1)
      CALL initRemap(mpolin,ntorin,nsin,nprocs,iam)
#endif
      mblk_size=3*(mpol+1)*(2*ntor+1)
!
!     Converge initial residues with diagonal preconditioner
      DIAGONALDONE=.FALSE.
#if defined(SKS)
      CALL second0(t1)
#endif
      CALL converge_diagonal (wout_file, ftol)
#if defined(SKS)
      CALL second0(t2)
      converge_diagonal_time=converge_diagonal_time+(t2-t1)
#endif
      DIAGONALDONE=.TRUE.
!
!     Converge using block preconditioner
#if defined(SKS)
      CALL second0(t1)
#endif
      CALL converge_blocks (wout_file, ftol)
#if defined(SKS)
      CALL second0(t2)
      converge_blocks_time=converge_blocks_time+(t2-t1)
#endif

      CALL second0(toff)
      skstoff=toff
      time_total = toff-ton
#if defined(SKS)
      total_time = skstoff-skston
#endif
      CALL RemoveDAFile

#if defined(SKS)
      IF (OUTPUT_TIMINGS) CALL GetTimes
#endif

      IF (iam .EQ. 0) THEN
      DO I = 6, 33, 27
         WRITE (I, 95) nprecon, levmarq_param, mupar, asym_index
#if defined(SKS)
         WRITE (I, '(a,i5)') ' Number processors: ', nprocs
#endif
         IF (I .NE. 6) WRITE (I,120) bsbu_ratio, jsju_ratio,            &
                       (bs0(j), j=2,6), (bu0(j), j=2,6)
         WRITE (I, 100)
         WRITE (I, 105)
         WRITE (I, 100)
         WRITE (I, 110) time_total, fbdy(1), time_init, fbdy(2),        &
                        time_diag_prec, fbdy(3), time_block_prec,       & 
                        fbdy(4), time_factor_blocks, fbdy(5),           &
                        time_toijsp, fbdy(6), time_tomnsp, gmres_time,  &
                        conj_grad_time, time_update_pres, fbdy(7),      &
                        time_update_bfield, fbdy(8), time_current,      &
                        fbdy(9), get_force_harmonics_time, fbdy(10),    &
                        time_update_force, fbdy(11), time_update_upperv,&
                        fbdy(12), time_update_state, fbdy(13),          &
                        time_funci, time_apply_precon,                  &
                        (diag_add_pert_time+block_add_pert_time),       &
                        time_init_state
         WRITE (I, 115)
!----------------------------------------------
         WRITE (I, *)
         IF (PARSOLVER) THEN
           WRITE (I, *)             'PARSOLVER      : NSCALED'
         ELSE
           WRITE (I, *)             'PARSOLVER      : ScaLAPACK'
         END IF
         WRITE (I, *)               'PARFUNCTISL    :', PARFUNCTISL
         WRITE (I, '(1x,a,L2,a,I2)')'PARGMRES       :', PARGMRES , ', GMRES_TYPE:', ngmres_type
#if defined(SKS)
         WRITE (I, *)               'OUTPUT_TIMINGS :', OUTPUT_TIMINGS
#endif
         WRITE (I, *)
         WRITE (I, '(a,1p,e12.4)')' TIME DIVB  : ',time_divb
         WRITE (I, '(a,1p,e12.4)')' TIME DIVJ  : ',time_divj
         WRITE (I, '(a,1p,e12.4)')' TIME BGRADP: ',time_bgradp
         WRITE (I, '(a,1p,e12.4)')' TIME BDOTJ : ',time_bdotj
         WRITE (I, *)
         WRITE(I, *) 'M (block size)         :',mblk_size
         WRITE(I, *) 'N (block rows)         :',nsin
         WRITE(I, *) 'P (processors)         :',nprocs
#if defined(SKS)
!         WRITE (I, '(1x,a,1p,2e12.4)')'Total time (max/min)   :',total_time_max, total_time_min
!         WRITE (I, *)
#endif
      END DO
      ENDIF

 95   FORMAT (/,' nprecon: ',i3, ' LM parameter: ',1pe9.2,              &
     &        ' mu||: ', 1pe9.2, ' Asym Index: ',1pe9.2)

 100  FORMAT(/'==============================              =======================')
 105  FORMAT(/,' TIMING INFORMATION                          RMS BOUNDARY FORCES')
 110  FORMAT(' Total runtime  : ', f12.3,'               fs(1,m=1)  :', 1pe10.2/,     &
             ' Initialization : ',0pf12.3,'               fs(2,m=1)  :', 1pe10.2/,    &
             ' Diagonal prec  : ',0pf12.3,'               fs(2,m!=1) :', 1pe10.2/,    &
             ' Compute blocks : ',0pf12.3,'               fu(1,m=1)  :', 1pe10.2/,    &
             ' Factor blocks  : ',0pf12.3,'               fu(2,m=1)  :', 1pe10.2/,    &
             ' Toijsp         : ',0pf12.3,'               fu(2,m!=1) :', 1pe10.2/,    &
             ' Tomnsp         : ',0pf12.3,                                      /,    &
             ' GMRES          : ',0pf12.3,                                      /,    &
             ' Conj Gradient  : ',0pf12.3,//,  &
             ' SUBROUTINES     ',/,            &
             ' Update Pressure: ',0pf12.3,'               fv(1,m=0)  :', 1pe10.2/,    &
             ' Update Bfield  : ',0pf12.3,'               fv(2,m=0)  :', 1pe10.2/,    &
             ' CV Currents    : ',0pf12.3,'               fv(2,m!=0) :', 1pe10.2/,    &
             ' Force Harmonics: ',0pf12.3,'               fu(ns)     :', 1pe10.2/,    &
             ' Update Force   : ',0pf12.3,'               fu(ns-1)   :', 1pe10.2/,    &
             ' Update UpperV  : ',0pf12.3,'               fv(ns)     :', 1pe10.2/,    &
             ' Update State   : ',0pf12.3,'               fv(ns-1)   :', 1pe10.2/,    &
             ' Funct Island   : ',0pf12.3/,    &
             ' Apply Precon   : ',0pf12.3/,    &
             ' Add Perturb    : ',0pf12.3/,    &
             ' Init State     : ',0pf12.3      &
     )
 115  FORMAT('==============================              =======================')
 120  FORMAT(' |B^s-r0*B^u|/|B^s+r0*B^u| (m=1,r->0): ',1pe10.3,/,                     &
             ' |J^s-r0*J^u|/|J^s+r0*J^u| (m=1,r->0): ',1pe10.3,/,                     &
             ' JBSUPSH(JS=2-6,M=1)/R0 : ',1p5e10.3,/,                                 &
             ' JBSUPUH(JS=2-6,M=1)    : ',1p5e10.3,/)

      CALL write_output (wout_file, niter)
      CALL write_profiles(fsq_total1)                 ! SPH: write jpmn, jbsupXmn, kbsubXmn, jvsupXmn profiles
      
      CALL cleanup_metric_elements

      CALL dealloc_quantities
      CALL dealloc_hessian
      CALL dealloc_diagnostics
      IF (lresistive) CALL dealloc_full_lower_metrics      
      IF (ALLOCATED(xc)) DEALLOCATE(xc, gc)

#if defined(MPI_OPT)
      IF (LSCALAPACK) THEN
        CALL blacs_barrier(icontxt,'All')
        CALL blacs_gridexit(icontxt)
      END IF
      CALL blacs_exit(0)

      IF ((myrow.eq.0).and.(mycol.eq.0)) CALL profstat()

#if defined(SKS)
      IF(SKSDBG) WRITE(TOFU,*) 'Called finalizeRemap' 
      IF(SKSDBG) FLUSH(TOFU)
      CALL finalizeRemap
#endif

      IF (iam .EQ. 0) THEN
      PRINT *,' Writing output to "siesta_profiles.txt" is finished!'
      CLOSE (unit=33)
      ENDIF
#endif

      END PROGRAM SIESTA
