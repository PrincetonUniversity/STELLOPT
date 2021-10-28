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
      USE stellopt_targets, ONLY: nu_dkes, sigma_dkes, lbooz, nsd
      ! NEO LIBRARIES
!DEC$ IF DEFINED (DKES_OPT)
      USE Vimatrix
      USE Vnamecl2
      USE Vcmain
      USE dkes_input, lscreen_dkes => lscreen
      USE dkes_realspace
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
      INTEGER :: ir, irun, istat, neqs, jk, ik, ier_phi !JLVG
      REAL(rprec), DIMENSION(:), POINTER :: f0p1, f0p2, f0m1, f0m2
      REAL(rprec) :: tcpu0, tcpu1, tcpui, tcput, tcpu, tcpua, dkes_efield,&
                     phi_temp
      CHARACTER :: tb*1           
      CHARACTER*50 :: arg1(7)       
      INTEGER :: numargs, icount, index_blank, iodata, iout_opt,&
                 idata, iout
      INTEGER :: m, n ! nmax, mmax (revmoed due to conflict with dkes_realspace
      CHARACTER :: output_file*64, opt_file*64, dkes_input_file*64, temp_str*64, temp_str2*64 !JLVG
      LOGICAL :: lexist  
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
!DEC$ IF DEFINED (DKES_OPT)
      lscreen_dkes = lscreen
      IF (lscreen) WRITE(6,'(a)') ' ---------------------------    DKES CALCULATION     -------------------------'
      ! Now we need to mimic some behavior to avoid reading from command line
      arg1(5) = 'F'
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
      ALLOCATE(DKES_L11p(nsd),DKES_L33p(nsd),DKES_L31p(nsd),&
               DKES_L11m(nsd),DKES_L33m(nsd),DKES_L31m(nsd),&
               DKES_scal11(nsd),DKES_scal33(nsd),DKES_scal31(nsd))
      DKES_L11p=0.0; DKES_L33p=0.0; DKES_L31p=0.0
      DKES_L11m=0.0; DKES_L33m=0.0; DKES_L31m=0.0
      DKES_scal11=0.0; DKES_scal33=0.0; DKES_scal31=0.0
      jk=0 !JLVG
      DO ik = 2, nsd
         ipmb = 0
         DKES_rad_dex = ik
         IF(.not. lbooz(ik)) CYCLE
         IF(sigma_dkes(ik) >= bigno) CYCLE
         tb = char(9)                                   !record file addition
         CALL second0 (tcpu0)
         ! Get ER based on phi specification
         ier_phi = 0
         ier_phi=-1 !JLVG
!         CALL get_equil_phi(rho(ik),phi_type,phi_temp,ier_phi,dkes_efield) !JLVG
         IF (ier_phi < 0) dkes_efield = 0.0
         !!!! CALL read_dkes_input
         WRITE(temp_str,'(i3.3)') ik
         arg1(1) = TRIM(proc_string) // '_s' // TRIM(temp_str)
         WRITE(arg1(2),'(i3)') ik
         WRITE(arg1(3),'(e20.10)') nu_dkes(ik)
         dkes_efield = 0.0
         WRITE(arg1(4),'(e20.10)') dkes_efield
         IF (lscreen) arg1(5) = 'T'
         ier_phi = 0
         !JLVG--------------------------------------------------------------------------------
         IF(myworkid == 0) CALL dkes_input_prepare(arg1,5,dkes_input_file,ier_phi)
         CALL MPI_BCAST(dkes_input_file,64,MPI_CHARACTER,0,MPI_COMM_MYWORLD,ierr_mpi)
         jk=jk+1
         IF(myworkid+1 == jk) temp_str2=dkes_input_file
      END DO
      dkes_input_file=temp_str2 !JLVG
      jk=0 
      DO ik = 2, nsd
         ipmb = 0
         DKES_rad_dex = ik
         WRITE(temp_str,'(i3.3)') ik
         IF(.not. lbooz(ik)) CYCLE
         IF(sigma_dkes(ik) >= bigno) CYCLE
         jk=jk+1 
         IF(myworkid+1 /= jk ) CYCLE 
         !JLVG--------------------------------------------------------------------------------
         output_file= 'dkesout.' // TRIM(proc_string) // '_s' // TRIM(temp_str)
         opt_file= 'opt_dkes.' // TRIM(proc_string) // '_s' // TRIM(temp_str)       !DAS 2/21/2000  !Probably won't need
         summary_file = 'results.' // TRIM(proc_string) //'_s' // TRIM(temp_str) !record file addition
         !  OPEN INPUT AND OUTPUT FILES FOR READING (AND WRITING OUTPUT)
         idata    = 7000+myworkid
         iout     = 30000+myworkid
         iout_opt = 14000+myworkid
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
         nrun = MAX(nrun, 1)
         nrun = MIN(nrun, krun)
         DO ir = 1, nrun
            irun = ir
            efield1 = efield(irun)
            cmul1 = cmul(irun)
            if(ir .eq. 1) then
               itab_out=90000+myworkid
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
               CALL residue (blk1, blk2, blk3, blk4, blk5, blk6,&
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
               CALL residue (blk1, blk2, blk3, blk4, blk5, blk6,&
                    f0m1, f0m2, srces0, rsd1m, rsd3m, g11m, g33m,&
                    g31m, g13m, crs1m, crs3m)
            ENDIF
            CALL dkes_printout (f0p1, f0m1, f0p2, f0m2, srces0)
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
      END DO
      !JLVG--------------------------------------------------------------------------------
      DO ik = 2, nsd
         CALL MPI_BCAST(DKES_L11p(ik),1,MPI_REAL8,ik,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(DKES_L11m(ik),1,MPI_REAL8,ik,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(DKES_L33p(ik),1,MPI_REAL8,ik,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(DKES_L33m(ik),1,MPI_REAL8,ik,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(DKES_L31p(ik),1,MPI_REAL8,ik,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(DKES_L31m(ik),1,MPI_REAL8,ik,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(DKES_scal11(ik),1,MPI_REAL8,ik,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(DKES_scal33(ik),1,MPI_REAL8,ik,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(DKES_scal31(ik),1,MPI_REAL8,ik,MPI_COMM_MYWORLD,ierr_mpi)
      END DO
      !JLVG--------------------------------------------------------------------------------
      IF (lscreen) WRITE(6,'(a)') ' -----------------  DKES CALCULATION (DONE) ----------------'
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
