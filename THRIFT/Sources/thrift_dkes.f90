!-----------------------------------------------------------------------
!     Subroutine:    thrift_dkes
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          08/22/2024
!     Description:   This subroutine calculates the DKES coefficients.
!                    see stellopt_dkes.f90 for more information.
!-----------------------------------------------------------------------
      SUBROUTINE thrift_dkes(lscreen,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE thrift_runtime
      USE thrift_input_mod
      USE thrift_vars, nrho_thrift => nrho
      USE thrift_profiles_mod
      USE thrift_equil
      USE thrift_funcs
      USE mpi_params
      USE mpi_inc
      USE mpi_sharmem
      USE read_wout_mod, ONLY: read_wout_file, read_wout_deallocate
      USE read_boozer_mod, ONLY: bcast_boozer_vars
      ! DKES LIBRARIES
#if defined(DKES_OPT)
      USE Vimatrix
      USE Vnamecl2
      USE Vcmain
      USE dkes_input, lscreen_dkes => lscreen
      USE dkes_realspace
#endif
      
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
      INTEGER :: i, j, k, l, istat, neqs, ier_phi, mystart, myend
      INTEGER, DIMENSION(:), ALLOCATABLE :: ik_dkes
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: Earr_dkes, nuarr_dkes
      REAL(rprec), DIMENSION(:), POINTER :: f0p1, f0p2, f0m1, f0m2
      REAL(rprec) :: tcpu0, tcpu1, tcpui, tcput, tcpu, tcpua, &
                     phi_temp, stime, etime
      CHARACTER :: tb*1           
      CHARACTER*50 :: arg1(6)       
      INTEGER :: numargs, iodata, iout_opt, idata, iout
      INTEGER :: m, n ! nmax, mmax (revmoed due to conflict with dkes_realspace
      CHARACTER :: output_file*64, opt_file*64, dkes_input_file*64, temp_str*64
      LOGICAL :: lfirst_pass  

!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (lscreen) WRITE(6,'(a)') ' --------------------  NEOCLASSICAL CALCULATION USING DKES  -------------------'
      IF (lvmec) THEN
         ierr_mpi = 0
         ! First we deallocate the old arrays and allocate new ones
#if defined(MPI_OPT)
         CALL MPI_BCAST(nruns_dkes,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
#endif
         IF (ALLOCATED(DKES_L11p)) DEALLOCATE(DKES_L11p)
         IF (ALLOCATED(DKES_L33p)) DEALLOCATE(DKES_L33p)
         IF (ALLOCATED(DKES_L31p)) DEALLOCATE(DKES_L31p)
         IF (ALLOCATED(DKES_L11m)) DEALLOCATE(DKES_L11m)
         IF (ALLOCATED(DKES_L33m)) DEALLOCATE(DKES_L33m)
         IF (ALLOCATED(DKES_L31m)) DEALLOCATE(DKES_L31m)
         IF (ALLOCATED(DKES_scal11)) DEALLOCATE(DKES_scal11)
         IF (ALLOCATED(DKES_scal33)) DEALLOCATE(DKES_scal33)
         IF (ALLOCATED(DKES_scal31)) DEALLOCATE(DKES_scal31)
         ALLOCATE(DKES_L11p(nruns_dkes),DKES_L33p(nruns_dkes),DKES_L31p(nruns_dkes),&
               DKES_L11m(nruns_dkes),DKES_L33m(nruns_dkes),DKES_L31m(nruns_dkes),&
               DKES_scal11(nruns_dkes),DKES_scal33(nruns_dkes),DKES_scal31(nruns_dkes))
         DKES_L11p=0.0; DKES_L33p=0.0; DKES_L31p=0.0
         DKES_L11m=0.0; DKES_L33m=0.0; DKES_L31m=0.0
         DKES_scal11=0.0; DKES_scal33=0.0; DKES_scal31=0.0
         l = 1
         ! Allocate Helper arrays to loop over dkes runs
         DKES_NK = COUNT(DKES_K>0)
         DKES_NC = COUNT(DKES_Nustar<1E10)
         DKES_NE = COUNT(DKES_Erstar<1E10)
         ALLOCATE(ik_dkes(nruns_dkes),Earr_dkes(nruns_dkes),nuarr_dkes(nruns_dkes))
         IF (myworkid == master) THEN
            DO k = 1, DKES_NS_MAX
               DO i = 1, DKES_NSTAR_MAX
                  DO j = 1, DKES_NSTAR_MAX
                     IF ((DKES_K(k) > 0) .AND. (DKES_Erstar(i)<1E10) .AND. (DKES_Nustar(j) < 1E10)) THEN
                        ik_dkes(l) = DKES_K(k)
                        Earr_dkes(l) = DKES_Erstar(i)
                        nuarr_dkes(l) = DKES_Nustar(j)
                        l = l + 1
                     END IF
                  END DO
               END DO
            END DO
         ENDIF
#if defined(MPI_OPT)
         ierr_mpi = 0
         CALL MPI_BCAST(ik_dkes,nruns_dkes,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(nuarr_dkes,nruns_dkes,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(Earr_dkes,nruns_dkes,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
#endif
         ! Read the wout file
         CALL read_wout_file(proc_string, ier)
         ! Now bcast the boozer parameters
         CALL bcast_boozer_vars(master, MPI_COMM_MYWORLD, ierr_mpi)
         ! Breakup work
         CALL MPI_CALC_MYRANGE(MPI_COMM_MYWORLD,1,nruns_dkes,mystart,myend)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!! Parallel Work block
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Note there are multile file write statements which can
         ! probably be removed so we only write what we need for PENTA
         CALL second0(stime)
         DO k = mystart,myend
            ! Setup input
            arg1(1) = TRIM(proc_string)
            WRITE(arg1(2),'(i3)') ik_dkes(k)
            WRITE(arg1(3),'(e20.10)') nuarr_dkes(k)
            WRITE(arg1(4),'(e20.10)') Earr_dkes(k) !dkes_efield
            arg1(5) = 'F'
            IF (lscreen .and. lfirst_pass) arg1(5) = 'T'
            WRITE(temp_str,'(i3.3)') k
            arg1(6) = '_k' // TRIM(temp_str)
            ier_phi = 0 ! We don't read the boozmn or wout file we've done that already
            CALL dkes_input_prepare_old(arg1,6,dkes_input_file,ier_phi)
            output_file= 'dkesout.' // TRIM(proc_string) // '_k' // TRIM(temp_str)
            opt_file= 'opt_dkes.' // TRIM(proc_string) // '_k' // TRIM(temp_str)       !DAS 2/21/2000  !Probably won't need
            summary_file = 'results.' // TRIM(proc_string) //'_k' // TRIM(temp_str) !record file addition
            !  OPEN INPUT AND OUTPUT FILES FOR READING (AND WRITING OUTPUT)
            idata    = 7
            iout     = 30
            iout_opt = 14
            iodata = idata
            CALL safe_open(iodata, istat, dkes_input_file, 'old', 'formatted')
            IF (istat .ne. 0) STOP 'Error reading input file in DKES'
            ioout = iout
            CALL safe_open(ioout, istat, output_file, 'replace', 'formatted')
            IF (istat .ne. 0) STOP 'Error writing output file'
            ioout_opt = iout_opt
            CALL safe_open(ioout_opt, istat, opt_file, 'replace','formatted')
            IF (istat .ne. 0) STOP 'Error writing opt_output file'
            ! Read namelist (datain) input
            lscreen_dkes = lscreen
            nvalsb(1) = -bigint-1
            idisk = 1
            lfout = 0
            ibbi = 1
            borbi = 0
            READ (iodata, nml=dkes_indata, iostat=istat)
            IF (istat .ne. 0) STOP 'Error reading dkes_indata NAMELIST in DKES'
            CLOSE (iodata)
            ! Recompute ntorb, mpolb for new style input where
            ! borbi is input with actual index value, borbi(n,m)
            IF (nvalsb(1) <= -bigint) THEN
               nmax = -bigint;  mmax = -bigint
               DO n = -ntorbd, ntorbd
                  DO m = 0, mpolbd
                     IF (borbi(n,m) .ne. zero) THEN
                        nmax = max (nmax, abs(n))
                        mmax = max (mmax, abs(m))
                     END IF
                  END DO
               END DO
               ! User MAY input smaller values if he does not want to use entire array         
               ntorb = min (abs(ntorb), nmax)
               mpolb = min (abs(mpolb-1), mmax)
            ELSE 
               IF (ntorb > ntorbd) THEN
                  WRITE (ioout, 45) ntorb, ntorbd
                  STOP ' ntorb > ntorbd in DKES input'
               END IF
               IF (mpolb < 2) THEN
                  WRITE (ioout, 20) mpolb
                  STOP ' mpolb < 2 in DKES input'
               ENDIF
            END IF
            IF (nzperiod .le. 0) nzperiod = 1
            IF (ipmb<0 .or. ipmb>2) ipmb = 0
            meshtz = MAX(0,meshtz)
            !!!! END read_dkes_input
            CALL set_mndim
            neqs = mpnt*(lalpha + 1)
            ALLOCATE (fzerop(2*neqs), fzerom(2*neqs), stat=istat)
            !IF (istat .ne. 0) STOP 'allocation error(1) in dkes2!'
            f0p1 => fzerop(1:neqs);   f0p2 => fzerop(neqs+1:)                  !+ functions
            f0m1 => fzerom(1:neqs);   f0m2 => fzerom(neqs+1:)                  !- functions
            CALL ftconv
            CALL lcalc
            CALL second0 (tcpu1); tcpui = tcpu1 - tcpu0
            CALL header
            CALL second0 (tcpu0); tcput = zero
            ! Here things get a bit screwy
            ! DKES allows for an array of nu/v and E/v to be evaluated
            ! but dkes_input_prepare does not.  So nrun will always be 1 for stellopt
            ! which is probably fine since we can parallelize over runs easier.
            ! But in the future this logic could be improved.
            nrun = MAX(nrun, 1)
            nrun = MIN(nrun, krun)
            DO i = 1, nrun
               efield1 = efield(i)
               cmul1 = cmul(i)
               if(i .eq. 1) then
                  call safe_open(itab_out, istat, summary_file,'unknown', &
                 'formatted')
                  write(itab_out,'("*",/,"cmul",a1,"efield",a1,"weov",a1,"wtov", &
                                     & a1,"L11m",a1,"L11p",a1,"L31m",a1,"L31p",a1,"L33m",a1,"L33p", &
                                     & a1,"scal11",a1,"scal13",a1,"scal33",a1,"max_residual", &
                                     & a1,"chip",a1,"psip",a1,"btheta",a1,"bzeta",a1,"vp")') &
                        tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb
               else if(i .gt. 1) then
                  open(itab_out,file=summary_file,status='unknown',position='append',form='formatted')
               endif
               CALL cescale (srces0)
               WRITE (ioout, 950) dashes, cmul1, efield1, weov, wtov, wcyclo, vthermi
               IF (ipmb < 2) THEN
                  iswpm = 1
                  CALL blk5d (blk1, blk2, blk3, blk4, blk5, blk6, blk7, f0p1, f0p2, srces0)
                  IF (ier .ne. 0) THEN
                     WRITE (ioout, 1000) ier
                     STOP
                  ENDIF
                  CALL residue_dkes (blk1, blk2, blk3, blk4, blk5, blk6,&
                                f0p1, f0p2, srces0, rsd1p, rsd3p, g11p, g33p,&
                                g31p, g13p, crs1p, crs3p)
               ENDIF
               IF (ipmb .ne. 1) THEN
                  iswpm = 2
                  CALL blk5d (blk1, blk2, blk3, blk4, blk5, blk6, blk7, f0m1, f0m2, srces0)
                  IF (ier .ne. 0) THEN
                     WRITE (ioout, 1050) ier
                     STOP
                  ENDIF
                  CALL residue_dkes (blk1, blk2, blk3, blk4, blk5, blk6,&
                       f0m1, f0m2, srces0, rsd1m, rsd3m, g11m, g33m,&
                       g31m, g13m, crs1m, crs3m)
               ENDIF
               ! This is a trick to get the arrays corretly sorted
               DKES_rad_dex = k !i ---> I think w/ 'i' is completely wrong!!
               IF (.not. lfirst_pass) lscreen_dkes = .FALSE.
               CALL dkes_printout (f0p1, f0m1, f0p2, f0m2, srces0)
               DKES_rad_dex = ik_dkes(k) !ik_dkes(i)  -- same here !!!
               ! End trickMPI_O
               CALL second0 (tcpu1); tcpu = tcpu1 - tcpu0; tcpu0 = tcpu1; tcput = tcput + tcpu; tcpua = tcput/i
               WRITE (ioout, 1100) tcpu
               CLOSE(unit=itab_out)
            END DO
            !IF (lfout .ne. 0) CALL wrout (f0p1, f0m1, f0p2, f0m2, srces0)  ! We don't need to do this
            CALL free_mndim
            DEALLOCATE (cols, al1, al2, al3, al4, bl1, bl2, bl3, bl4, cl1,&
               cl2, cl3, cl4, cols0, omgl, al01, al02, al03, al04, bl01,&
               bl02, bl03, bl04, cl01, cl02, cl03, cl04, fzerop, fzerom)
            CLOSE(unit=ioout)
            CLOSE(unit=ioout_opt)
            lfirst_pass = .FALSE.
            CALL second0(etime)
         END DO
         PRINT *, 'I took ', etime-stime,'s.'
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!! Parallel Work block Done
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(MPI_OPT)
         CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi)
         IF (myworkid == master) THEN
            CALL MPI_REDUCE(MPI_IN_PLACE, DKES_L11p, nruns_dkes, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE, DKES_L11m, nruns_dkes, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE, DKES_L33p, nruns_dkes, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE, DKES_L33m, nruns_dkes, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE, DKES_L31p, nruns_dkes, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE, DKES_L31m, nruns_dkes, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE, DKES_scal11, nruns_dkes, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE, DKES_scal33, nruns_dkes, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE, DKES_scal31, nruns_dkes, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
         ELSE
            CALL MPI_REDUCE(DKES_L11p,    DKES_L11p, nruns_dkes, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
            CALL MPI_REDUCE(DKES_L11m,    DKES_L11m, nruns_dkes, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
            CALL MPI_REDUCE(DKES_L33p,    DKES_L33p, nruns_dkes, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
            CALL MPI_REDUCE(DKES_L33m,    DKES_L33m, nruns_dkes, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
            CALL MPI_REDUCE(DKES_L31p,    DKES_L31p, nruns_dkes, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
            CALL MPI_REDUCE(DKES_L31m,    DKES_L31m, nruns_dkes, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
            CALL MPI_REDUCE(DKES_scal11,  DKES_scal11, nruns_dkes, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
            CALL MPI_REDUCE(DKES_scal33,  DKES_scal33, nruns_dkes, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
            CALL MPI_REDUCE(DKES_scal31,  DKES_scal31, nruns_dkes, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
         END IF
#endif
         !IF (lscreen) WRITE(6,'(a)') ' -----------------  DKES CALCULATION (DONE) ----------------'
         IF (myworkid .ne. master) CALL read_wout_deallocate

         ! Sort out the DKES runs for later use
         IF (ASSOCIATED(DKES_D11)) CALL mpidealloc(DKES_D11,win_dkes_d11)
         IF (ASSOCIATED(DKES_D31)) CALL mpidealloc(DKES_D31,win_dkes_d31)
         IF (ASSOCIATED(DKES_D33)) CALL mpidealloc(DKES_D33,win_dkes_d33)
         CALL mpialloc(DKES_D11, DKES_NK, DKES_NC, DKES_NE, myid_sharmem, 0, MPI_COMM_MYWORLD, win_dkes_d11)
         CALL mpialloc(DKES_D31, DKES_NK, DKES_NC, DKES_NE, myid_sharmem, 0, MPI_COMM_MYWORLD, win_dkes_d31)
         CALL mpialloc(DKES_D33, DKES_NK, DKES_NC, DKES_NE, myid_sharmem, 0, MPI_COMM_MYWORLD, win_dkes_d33)

         IF (myworkid .eq. master) THEN
            DKES_D11 = RESHAPE( (DKES_L11p + DKES_L11m)*0.5, shape=(/DKES_NK, DKES_NC, DKES_NE/), order=(/2,3,1/) )
            DKES_D31 = RESHAPE( (DKES_L31p + DKES_L31m)*0.5, shape=(/DKES_NK, DKES_NC, DKES_NE/), order=(/2,3,1/) )
            DKES_D33 = RESHAPE( (DKES_L33p + DKES_L33m)*0.5, shape=(/DKES_NK, DKES_NC, DKES_NE/), order=(/2,3,1/) )
         END IF
      END IF
      IF (lscreen) WRITE(6,'(a)') ' -------------------  DKES NEOCLASSICAL CALCULATION DONE  ---------------------'
      RETURN
   90 format(5e16.8)
  950 FORMAT(/9x,'CMUL',7x,'EFIELD',4x,'OMEGA-E/v',4x,'OMEGA-T/v',&
               5x,'OMEGA-ci',3x,'VI-THERMAL'/&
               7x,'(nu/v)',7x,'(Es/v)',2x,'(ExB drift)',4x,'(transit)',&
               3x, '(H+, B=B00)',2x,'(H+, Ti=1keV)'/,131a,/&
               1p,6e13.4/)
   20 FORMAT(' mpolb = ',i5,'  is less than 2')
   45 FORMAT(' ntorb = ',i5,'  is greater than ntorbd = ',i5)
 1000 FORMAT(/' blk5d error in block = ',i5)
 1050 FORMAT(/' blk5d error : ierm = ',i5)
 1100 FORMAT(/' time used:    tcpu =',1p,e10.2,'  sec'/)
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_dkes
