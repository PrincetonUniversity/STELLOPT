!-----------------------------------------------------------------------
!     Subroutine:    stellopt_terpsichore
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          06/06/2012
!     Description:   This subroutine handles calculating kink stability.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_kink(lscreen,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_input_mod
      USE stellopt_vars
      USE equil_vals
      USE safe_open_mod, ONLY: safe_open
      USE read_wout_mod, ONLY: xm_vmec => xm, xn_vmec => xn, mnmax_vmec => mnmax, &
                               wb_vmec => wb, gamma_vmec => gamma, mpol_vmec => mpol, &
                               ntor_vmec => ntor, itfsq_vmec => itfsq, niter_vmec => niter, &
                               iota_vmec => iotas, gmnc_vmec => gmnc, &
                               rmnc_vmec => rmnc, zmns_vmec => zmns, &
                               mass_vmec => mass, pres_vmec => pres, &
                               phip_vmec => phip, vp_vmec => vp
      USE vmec_input, ONLY: lfreeb_vmec => lfreeb
      USE mpi_params
      USE mpi_inc

!DEC$ IF DEFINED (TERPSICHORE)
      USE tpr_param
!DEC$ ENDIF
      
!-----------------------------------------------------------------------
!     Subroutine Parameters
!        iflag         Error flag
!----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL, INTENT(in)    :: lscreen
      INTEGER, INTENT(inout) :: iflag
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      INTEGER ::  ier, ik, mn, mystart, myend, mypace, nworkers, nmodes
      INTEGER, ALLOCATABLE :: mnum(:)
      REAL(rprec) :: sigc, tauc, pbpc, pppc, fourpi, rmu0
      CHARACTER(256) :: num_str
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (ALLOCATED(wp_kink)) DEALLOCATE(wp_kink,wk_kink,growth_kink,omega_kink)
      ALLOCATE(wp_kink(nsys),wk_kink(nsys),growth_kink(nsys),omega_kink(nsys))
      wp_kink = 0; wk_kink = 0; growth_kink = 0; omega_kink = 0;
      IF (lscreen) WRITE(6,'(a)') ' ---------------------------  KINK STABILITY  -------------------------'
      SELECT CASE(TRIM(equil_type))
         CASE('vmec2000','parvmec','paravmec','vboot','vmec2000_oneeq')
!DEC$ IF DEFINED (TERPSICHORE)
            ! Setup TERPSICHORE Variables
!DEC$ IF DEFINED (MPI_OPT)
            CALL MPI_BCAST(nrad,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_kink: BCAST nrad',ierr_mpi)
            CALL MPI_BCAST(mnmax_vmec,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_kink: BCAST mnmax_vmec',ierr_mpi)
!DEC$ ENDIF
            NI    = nrad - 1
            MMAXDF = mmaxdf_kink
            NMAXDF = nmaxdf_kink
            MLMNV = mnmax_vmec
            MLMNB = mlmnb_kink
            IVAC=IVAC_KINK
            NVI   = NI+IVAC
            ND=NVI
            ND1=ND+1
            iunit_in = 25

            ! Create TERPSICHORE equilibrium file
            IF (myworkid == master) THEN
               fourpi = 16 * ATAN(1.0)
               rmu0 = fourpi * 1.0e-7
               CALL safe_open(iunit_eq,ier,'terpsichore_eq.'//TRIM(proc_string),'unknown','formatted')
               WRITE(iunit_eq,840)wb_vmec,gamma_vmec,1./REAL(NFP),0.,mnmax_vmec,nrad,mpol_vmec,&
                               ntor_vmec,1,1,itfsq_vmec,niter_vmec,0
               pbpc=0; pppc=0
               DO ik = 1, nrad
                  DO mn = 1, mnmax_vmec
                     IF (xm_vmec(mn) == 0 .and. xn_vmec(mn) == 0) THEN
                        sigc = 1;  tauc = 1
                     ELSE
                        sigc = 0;  tauc = 0
                     END IF
                     WRITE(iunit_eq,850) xm_vmec(mn), xn_vmec(mn), rmnc_vmec(mn,ik), zmns_vmec(mn,ik), gmnc_vmec(mn,ik)
                     WRITE(iunit_eq,850) sigc ,tauc, pbpc, pppc
                  END DO !mn
               END DO !ik
               WRITE(iunit_eq,850)(iota_vmec(ik),mass_vmec(ik),pres_vmec(ik)*rmu0,-phip_vmec(ik),vp_vmec(ik),ik=1,nrad)
               CLOSE(iunit_eq)
            END IF !MPI root

            ! This section divides up the work for parallel execution
            mystart = 1; myend = nsys; mypace = 1
!DEC$ IF DEFINED (MPI_OPT)
            CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi)
            IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_kink: BARRIER1',ierr_mpi)
            nmodes = COUNT(sigma_kink<bigno)
            CALL MPI_COMM_SIZE( MPI_COMM_MYWORLD, nworkers, ierr_mpi )
            IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_kink: MPI_COMM_SIZE',ierr_mpi)
           
            IF (ALLOCATED(mnum)) DEALLOCATE(mnum)
            ALLOCATE(mnum(nworkers))
            mnum = 0
            DO
               DO ik = 1, nworkers
                  mnum(ik) = mnum(ik) + 1
                  IF (SUM(mnum) == nmodes) EXIT
               END DO
               IF (SUM(mnum) == nmodes) EXIT
            END DO
            mypace  = mnum(myworkid+1)
            IF (myworkid == 0) THEN
               mystart = 1
            ELSE
               mystart = SUM(mnum(1:myworkid))+1
            END IF
            myend   = mystart+mypace-1
            DEALLOCATE(mnum)
            IF (myend < mystart) myend = mystart
!DEC$ ENDIF

            ! Now we run TERPSICHORE
            CALL safe_open(iunit_eq,ier,'terpsichore_eq.'//TRIM(proc_string),'unknown','formatted')
            DO ik = mystart, myend
               IF (ik == 0) EXIT
               IF (sigma_kink(ik) >= bigno) CYCLE
               lrun_error = .false.
               NJ    = nj_kink(ik)
               NK    = nk_kink(ik)
               NJK   = NJ*NK
               MLMNS = mlmns_kink(ik)
               LSSL  = lssl_kink(ik)
               LSSD  = lssd_kink(ik)
               MD=MLMNS
               MDY=MLMNS
               NA=(MD+MDY)*ND+MD
               lscreen_trp = .FALSE.
               IF (lscreen .and. ik == 1) lscreen_trp = lscreen
               REWIND(iunit_eq)
               WRITE(num_str,'(1(A,I2.2))') '_',ik-1
               CALL safe_open(iunit_in,ier,'terpsichore_input'//TRIM(num_str),'unknown','formatted')
               CALL safe_open(iunit_16,ier,'terpsichore_16.'//TRIM(proc_string)//TRIM(num_str),'unknown','formatted')
               CALL safe_open(iunit_17,ier,'terpsichore_17.'//TRIM(proc_string)//TRIM(num_str),'unknown','formatted')
               CALL safe_open(iunit_19,ier,'terpsichore_19.'//TRIM(proc_string)//TRIM(num_str),'unknown','formatted')
               CALL safe_open(iunit_22,ier,'terpsichore_22.'//TRIM(proc_string)//TRIM(num_str),'unknown','formatted')
               CALL safe_open(iunit_23,ier,'terpsichore_23.'//TRIM(proc_string)//TRIM(num_str),'unknown','unformatted')
               ! SIMULATE TERPSICHORE
               call TPRALLOC
               CALL EQINVM
               CALL VEQREC
               CALL MTASKB
               IF (lrun_error) THEN
                  wp_kink(ik)     = bigno
                  wk_kink(ik)     = bigno
                  omega_kink(ik)  = bigno
                  growth_kink(ik) = bigno
                  CALL TPRDEALLOC
                  CLOSE(iunit_in)
                  CLOSE(iunit_16)
                  CLOSE(iunit_17)
                  CLOSE(iunit_19)
                  CLOSE(iunit_22)
                  CLOSE(iunit_23)
                  CYCLE
               END IF
               CALL STABIN
               CALL MTASKA
               CALL MTASKS_AP
               CALL TPRDEALLOC
               CLOSE(iunit_in)
               CLOSE(iunit_16)
               CLOSE(iunit_17)
               CLOSE(iunit_19)
               CLOSE(iunit_22)
               CLOSE(iunit_23)
               wp_kink(ik) = wp_global
               wk_kink(ik) = wk_global
               omega_kink(ik) = omega_global
               growth_kink(ik) = growth_global
            END DO
            CLOSE(iunit_eq)

            ! Get all the data to master
!DEC$ IF DEFINED (MPI_OPT)
            CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi)
            IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_kink: BARRIER2',ierr_mpi)
            !PRINT *,myworkid,wp_global,wp_kink(mystart)
            IF (myworkid == master) THEN
               CALL MPI_REDUCE(MPI_IN_PLACE,wp_kink,nsys,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,wk_kink,nsys,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,omega_kink,nsys,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,growth_kink,nsys,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            ELSE
               CALL MPI_REDUCE(wp_kink,wp_kink,nsys,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(wk_kink,wk_kink,nsys,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(omega_kink,omega_kink,nsys,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(growth_kink,growth_kink,nsys,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               DEALLOCATE(wp_kink, wk_kink, omega_kink,growth_kink)
            END IF
!DEC$ ENDIF

            ! PRINT SOME STUFF TO SCREEN
            IF (lscreen) THEN
               WRITE(6,'(a)') '    RUN       WP           WK         OMEGA     GROWTH_RATE'
               DO ik = 1, nsys
                  IF (sigma_kink(ik) >= bigno) CYCLE
                  WRITE(6,'(5X,I2.2,4(1X,E12.5))') ik,wp_kink(ik),wk_kink(ik),omega_kink(ik),growth_kink(ik)
               END DO
            END IF
!DEC$ ENDIF
         CASE('spec')
      END SELECT
      IF (lscreen) WRITE(6,'(a)') ' -------------------------  KINK STABILITY DONE  ----------------------'
      RETURN
 840  FORMAT(1x,1pe22.12,1x,1p3e15.5,8i6,i3)
 850  FORMAT(1x,1p5e22.14)
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_kink
