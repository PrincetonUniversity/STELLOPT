!-----------------------------------------------------------------------
!     Subroutine:    stellopt_dkes
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          06/06/2012
!     Description:   This subroutine calculates upper and lower bounds
!                    for the diffusion coefficient.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_dkes(lscreen,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime, ONLY:  proc_string, bigno, rprec
      USE equil_utils, ONLY: get_equil_phi, nrad, shat, phi_type
      USE stellopt_targets, ONLY: nu_dkes, sigma_dkes, lbooz, nsd, E_dkes, nprof, nruns_dkes
      ! NEO LIBRARIES
!DEC$ IF DEFINED (DKES_OPT)
      USE Vimatrix
      USE Vnamecl2
      USE Vcmain
      USE dkes_input, lscreen_dkes => lscreen
      USE dkes_realspace
      USE read_wout_mod, ONLY: read_wout_file, read_wout_deallocate
      USE read_boozer_mod, ONLY: bcast_boozer_vars
!DEC$ ENDIF
      USE safe_open_mod
      USE mpi_params
      USE mpi_inc
      
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
      INTEGER :: ir, irun, istat, neqs, ik, ij, ier_phi, mystart, myend
      INTEGER, DIMENSION(:), ALLOCATABLE :: ik_dkes
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: Earr_dkes, nuarr_dkes
      REAL(rprec), DIMENSION(:), POINTER :: f0p1, f0p2, f0m1, f0m2
      REAL(rprec) :: tcpu0, tcpu1, tcpui, tcput, tcpu, tcpua, dkes_efield,&
                     phi_temp
      CHARACTER :: tb*1           
      CHARACTER*50 :: arg1(6)       
      INTEGER :: numargs, icount, index_blank, iodata, iout_opt,&
                 idata, iout
      INTEGER :: m, n ! nmax, mmax (revmoed due to conflict with dkes_realspace
      CHARACTER :: output_file*64, opt_file*64, dkes_input_file*64, temp_str*64
      LOGICAL :: lexist, lfirst_pass  
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
!DEC$ IF DEFINED (DKES_OPT)
      lscreen_dkes = lscreen
      lfirst_pass = .TRUE.
      IF (lscreen) WRITE(6,'(a)') ' ---------------------------    DKES CALCULATION     -------------------------'
!DEC$ IF DEFINED (MPI_OPT)
      ierr_mpi = 0
      CALL MPI_BCAST(nruns_dkes,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
!DEC$ ENDIF
      ! Enter the main loop
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
      ! Setup the helper arrays
      IF (ALLOCATED(ik_dkes)) DEALLOCATE(ik_dkes)
      IF (ALLOCATED(nuarr_dkes)) DEALLOCATE(nuarr_dkes)
      IF (ALLOCATED(Earr_dkes)) DEALLOCATE(Earr_dkes)
      ALLOCATE(ik_dkes(nruns_dkes),nuarr_dkes(nruns_dkes),Earr_dkes(nruns_dkes))
      IF (myworkid == master) THEN
         ik = 0
         DO ir = 1, nsd
            IF (sigma_dkes(ir) >= bigno) CYCLE
            DO ij = 1, nprof
               IF (E_dkes(ij) <= -bigno .or. nu_dkes(ij) <= -bigno) CYCLE
               ik = ik + 1
               ik_dkes(ik)    = ir
               nuarr_dkes(ik) = nu_dkes(ij)
               Earr_dkes(ik)  = E_dkes(ij)
            END DO
         ENDDO
      END IF
      ! Now read the wout file
      CALL read_wout_file(proc_string, ier)
      ! Now bcast the boozer parameters
      CALL bcast_boozer_vars(master, MPI_COMM_MYWORLD, ierr_mpi)
!DEC$ IF DEFINED (MPI_OPT)
      ierr_mpi = 0
      CALL MPI_BCAST(ik_dkes,nruns_dkes,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_BCAST(nuarr_dkes,nruns_dkes,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_BCAST(Earr_dkes,nruns_dkes,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
!DEC$ ENDIF
      CALL MPI_CALC_MYRANGE(MPI_COMM_MYWORLD,1,nruns_dkes,mystart,myend)
      ! Loop over radial surfaces
      DO ik = mystart,myend
         ipmb = 0
         DKES_rad_dex = ik_dkes(ik)
         IF(.not. lbooz(ik_dkes(ik))) CYCLE
         tb = char(9)                                   !record file addition
         CALL second0 (tcpu0)
         dkes_efield = Earr_dkes(ik)
         arg1(1) = TRIM(proc_string)
         !!!! CALL read_dkes_input
         WRITE(arg1(2),'(i3)') ik_dkes(ik)
         WRITE(arg1(3),'(e20.10)') nuarr_dkes(ik)
         WRITE(arg1(4),'(e20.10)') dkes_efield
         arg1(5) = 'F'
         IF (lscreen .and. lfirst_pass) arg1(5) = 'T'
         WRITE(temp_str,'(i3.3)') ik
         arg1(6) = '_s' // TRIM(temp_str)
         ier_phi = 0 ! We don't read the boozmn or wout file we've done that already
         CALL dkes_input_prepare(arg1,6,dkes_input_file,ier_phi)
         output_file= 'dkesout.' // TRIM(proc_string) // '_s' // TRIM(temp_str)
         opt_file= 'opt_dkes.' // TRIM(proc_string) // '_s' // TRIM(temp_str)       !DAS 2/21/2000  !Probably won't need
         summary_file = 'results.' // TRIM(proc_string) //'_s' // TRIM(temp_str) !record file addition
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
         DO ir = 1, nrun
            irun = ir
            efield1 = efield(irun)
            cmul1 = cmul(irun)
            if(ir .eq. 1) then
               call safe_open(itab_out, istat, summary_file,'unknown', &
              'formatted')
               write(itab_out,'("*",/,"cmul",a1,"efield",a1,"weov",a1,"wtov", &
                                  & a1,"L11m",a1,"L11p",a1,"L31m",a1,"L31p",a1,"L33m",a1,"L33p", &
                                  & a1,"scal11",a1,"scal13",a1,"scal33",a1,"max_residual", &
                                  & a1,"chip",a1,"psip",a1,"btheta",a1,"bzeta",a1,"vp")') &
                     tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb
            else if(ir .gt. 1) then
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
            DKES_rad_dex = ik
            IF (.not. lfirst_pass) lscreen_dkes = .FALSE.
            CALL dkes_printout (f0p1, f0m1, f0p2, f0m2, srces0)
            DKES_rad_dex = ik_dkes(ik)
            ! End trick
            CALL second0 (tcpu1); tcpu = tcpu1 - tcpu0; tcpu0 = tcpu1; tcput = tcput + tcpu; tcpua = tcput/irun
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
      END DO
!DEC$ IF DEFINED (MPI_OPT)
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
!DEC$ ENDIF
      IF (lscreen) WRITE(6,'(a)') ' -----------------  DKES CALCULATION (DONE) ----------------'
      IF (myworkid .ne. master) CALL read_wout_deallocate
!DEC$ ELSE
      IF (lscreen) WRITE(6,'(a)') ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      IF (lscreen) WRITE(6,'(a)') ' !! DKES based optimization not supported on your machine !!'
      IF (lscreen) WRITE(6,'(a)') ' !! Check for DKES directory and DKES_OPT in your makefile!!'
      IF (lscreen) WRITE(6,'(a)') ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!DEC$ ENDIF
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
      END SUBROUTINE stellopt_dkes
