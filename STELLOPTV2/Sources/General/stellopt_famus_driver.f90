!-----------------------------------------------------------------------
!     Subroutine:    stellopt_famus_driver
!     Authors:       J.C.Schmitt (Auburn/PPPL) jcschmitt@auburn.edu
!     Date:          2020
!     Description:   This subroutine calls the coil/permanent magnet
!                    optimization code, FAMUS.
!                    This driver code does not attempt to re-assign
!                    the communicators groups for the focus runs. The
!                    communicator MPI_COMM_MYWORLD is used for this 
!                    iteration of FAMUS run.
!                    
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_famus_driver(file_str, lscreen, iflag)


!DEC$ IF DEFINED (FAMUS)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_input_mod
      USE stellopt_vars, !  my_mpol => mpol_fds, my_ntor => ntor_fds
      USE equil_utils
      USE neswrite, ONLY: coil_separation
      USE read_wout_mod, ONLY: rmnc_vmec => rmnc, rmns_vmec => rmns, &
                               zmnc_vmec => zmnc, zmns_vmec => zmns, &
                               ierr_vmec, read_wout_file, &
                               read_wout_deallocate, mnmax_vmec => mnmax, &
                               mnmax_nyq_vmec => mnmax_nyq, nfp_vmec => nfp, &
                               xm_vmec => xm, xn_vmec => xn, lasym_vmec => lasym

      USE famus_globals, ONLY: nfp_famus => Nfp, IsQuiet_famus => IsQuiet, &
           famus_id => myid, famus_ext => ext, &
           famus_Nzeta => Nzeta, &
           famus_Nteta => Nteta, famus_ncpu => ncpu, famus_xdof => xdof, &
                        famus_MPI_COMM_MASTERS => MPI_COMM_MASTERS, &
                        famus_MPI_COMM_MYWORLD => MPI_COMM_MYWORLD, &
                        famus_MPI_COMM_WORKERS => MPI_COMM_WORKERS, &
                        famus_MPI_COMM_FAMUS => MPI_COMM_FAMUS, &
                        famus_call_from_ext_opt => call_from_ext_opt, &
                        famus_init_from_ext_opt => init_from_ext_opt, &
                        famus_case_surface => case_surface, &
                        famus_case_coils => case_coils, &
                        famus_input_coils => input_coils, &
                        famus_case_postproc => case_postproc, &
                        famus_case_optimize => case_optimize, &
                        famus_case_bnormal => case_bnormal, &
                        famus_input_surf => input_surf, &
                        famus_surf => surf
!                                input_surf, input_coils, fixed_coils, input_harm
!                               famusin_filename_from_stellopt, famus_called_from_stellopt => called_from_stellopt


      USE mpi_params
      USE mpi_inc
             
!DEC$ ENDIF

!-----------------------------------------------------------------------
!     Subroutine Parameters
!        iflag         Error flag
!----------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(256), INTENT(inout)    :: file_str
      LOGICAL, INTENT(inout)        :: lscreen
      INTEGER, INTENT(inout) :: iflag
!DEC$ IF DEFINED (FAMUS)

!-----------------------------------------------------------------------
!     Local Variables
!        istat         Error status
!        iunit         File unit number
      ! FOR FAMUS
      INTEGER :: istat, iunit = 112, iunit2=113,  m, n, ii, imn, num_bn_modes
      INTEGER :: ierr, num_vmec_modes, num_vmec_surfs, tmm, tnn
      REAL    :: bnc, trmnc, trmns, tzmnc, tzmns
      INTEGER :: mm, nn, famus_rankWorld, famus_color, famus_key
      CHARACTER(400) :: famusin_out_filename, famus_bnorm_filename,  &
                        famus_plasmaboundary_filename
      INTEGER, dimension(2) :: vmec_sizes
      integer :: nprocs_myworld, myrank_myworld


!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      CALL MPI_COMM_SIZE( MPI_COMM_MYWORLD, nprocs_myworld, ierr_mpi )
      print *, '<----famus_driver entry famus_id = ',famus_id,' myid=',myid, ' of ',nprocs_myworld,' in myworld', &
               ' proc_string=',trim(proc_string), ' file_str=',trim(file_str), &
               ' iflag=', iflag, ' myworkid = ', myworkid,' famus_dc_ox(1)=',famus_dc_ox(1)
      IF (iflag < 0) RETURN
      CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi)
      print *,'<---checkup post entry mpi_barrier 1: famus_dc_ox(1)=',famus_dc_ox(1)
      lscreen = .true.
      IF (lscreen) then
         WRITE(6,'(a)') 'lscreen=true : -------------  FAMUS CALCULATION 2 ---------'
      ENDIF

      ! Suppress FAMUS stdout if needed
      IsQuiet_famus = 1
      if (lscreen) IsQuiet_famus = -1

      ! Setup MPI stuff
      !myworkid = master
      famus_color = 0
      famus_key = myworkid
!      CALL MPI_COMM_FREE(FAMUS_MPI_COMM_FAMUS, ierr_mpi)
      !FAMUS_MPI_COMM_MYWORLD = MPI_COMM_MYWORLD
      FAMUS_MPI_COMM_FAMUS = MPI_COMM_MYWORLD
!      CALL MPI_COMM_SPLIT(MPI_COMM_MYWORLD, famus_color, famus_key, FAMUS_MPI_COMM_FAMUS, ierr_mpi)
      CALL MPI_COMM_SIZE(FAMUS_MPI_COMM_FAMUS, famus_ncpu, ierr_mpi)
      CALL MPI_COMM_RANK(FAMUS_MPI_COMM_FAMUS, myrank_myworld, ierr_mpi)
      famus_id = myrank_myworld
      print *, '<----famus_driver post comm redifine famus_id = ',famus_id,' of ',famus_ncpu,' myid=',myid, &
               ' proc_string=',trim(proc_string), ' file_str=',trim(file_str), &
               ' iflag=', iflag, ' myworkid = ', myworkid,' famus_dc_ox(1)=',famus_dc_ox(1)

!     famus_ncpu = famus_ncpu
      CALL MPI_BARRIER(FAMUS_MPI_COMM_FAMUS,ierr_mpi)

! STEP 1.
      print *, '<----famus_driver post mpi split myid=', myid, &
               ' proc_string=',trim(proc_string), ' file_str=',trim(file_str), &
               ' iflag=', iflag, ' myworkid = ', myworkid,' famus_id=',famus_id, &
               ' famus_ncpu=', famus_ncpu, '--->'
      ! Run bnorm if required ( JCS - is this handled correctly in the MPI
      ! environment? Do all procs help, or does only master do this?)
      !if (load_bnorm_famus) then
      if (famus_id == master) then
      print *,'<---checkup 3: famus_dc_ox(1)=',famus_dc_ox(1)
         print *, '<----stellopt_famus_driver: famus_id = ',famus_id,' (master) is calling stellopt_bnorm(', &
                    trim(file_str), ', ', lscreen,')'
         coil_separation = 0.1 ! JCS - Does this matter?
         ! Run BNORM code
         call stellopt_bnorm(file_str,lscreen)
      print *,'<---checkup 4: famus_dc_ox(1)=',famus_dc_ox(1)
         famus_bnorm_filename = 'bnorm.' // TRIM(file_str)

         IF (lscreen) WRITE(6,'(a,a)') '<----stellopt_famus_driver: proc_string=', trim(proc_string)
         !wout_filename = 'wout_'//TRIM(proc_string)//'.nc'
   
         ! STELLOPT (via lmdif->stellopt_fcn or similar) will modifiy the value of
         ! the boundary coefficients. Here, the value of the FAMUS variables
         ! are loaded with the new values, the correct mode of operation is
         ! determiend, control variables are set, and the regcoil-functions
         ! are called

!  things to update and/or write out:
! yes    input_surf = 'plasma.boundary'
! no?   fixed_coils = 'only_Bt_positive_symmetry.focus'
! yes   input_coils = 'ellipse_reg128_P64_T64.focus' 


! STEP 2. Write the dipole coil information to file so that FAMUS can read it
          write(6,*) '<----famus_driver: Step 2a: writing dipole coil file'
     
         ! The famus variables are written to file, first (read in by FAMUS later).
         ! During 'cleanup', the new minimum will be copied/moved
         ! to have the iteration number as a suffix
!      print *,'<---checkup 5: famus_dc_ox(1)=',famus_dc_ox(1)
         famus_dipolecoil_filename = trim(proc_string) // '.focus'
         CALL safe_open(iunit, istat, famus_dipolecoil_filename, 'replace', 'formatted')
         write (iunit, '(a)') ' # Total number of coils,  momentq'
         write (iunit, '(I6, (a), I6)') famus_ncoils, ',', 1
         write (iunit, '(a)') '#coiltype, symmetry,  coilname,  ox,  oy,  oz,  Ic,  M_0,  pho,  Lc,  mp,  mt'

         DO m = 1, famus_ncoils
           !print *,'<--m=',m,'type=',famus_dc_itype(m)
           write(iunit,'(2(I4, ", ")A13,", ",3(ES15.8,", "),2(I2,", ",ES15.8,", ",ES15.8,", "))') &
             famus_dc_itype(m), famus_dc_symmetry(m), trim(famus_dc_name(m)), &
             famus_dc_ox(m),famus_dc_oy(m),  famus_dc_oz(m), &
             famus_dc_Ic(m),famus_dc_M_0(m),  famus_dc_pho(m), &
             famus_dc_Lc(m),famus_dc_mp(m),  famus_dc_mt(m)
          END DO
          write(iunit, '(A)') ''
          CLOSE(iunit)

       ! 
          write(6,*) '<----famus_driver: finished writing dipole coil file: ', &
                   trim(famus_dipolecoil_filename)


          write(6,*) '<----famus_driver: Step 2b: writing plasma boundary file'

          ! Need to retrieve the plasma boundary and the BN harmonics on the
          ! boundary. The VMEC data can be loaded. The BN data file will be parsed
          ! twice: once to determine the # of modes, and a second time to copy the
          ! contents over to the plasma boundary file
   
          ! Retrieve plasma boundary information
          ! Read the VMEC output
          CALL read_wout_deallocate
          CALL read_wout_file(TRIM(proc_string),ierr)
          IF (ierr .ne. 0 .or. ierr_vmec .ne. 0) THEN
               iflag = -1      
          RETURN          
          END IF       
          ! VMEC data should now be available
          print *,'<----famus_driver: vmec mnmax=',mnmax_vmec, &
                  'mnmax_nyq=',mnmax_nyq_vmec
          vmec_sizes = shape(rmnc_vmec)
          num_vmec_modes = vmec_sizes(1)
          num_vmec_surfs = vmec_sizes(2)
          print *,'<----famus_driver: # vmec modes = ', num_vmec_modes
          print *,'<----famus_driver: # vmec surfs = ', num_vmec_surfs

          ! Retrieve BN infomration
          if (.True.) print *,"Loading B_normal on the plasma surface due to plasma current from file ",trim(famus_bnorm_filename)
   
          call safe_open(iunit2, ierr, trim(famus_bnorm_filename), 'old', 'formatted')
          if (ierr .ne. 0 ) then
             stop 'Unable to open bnorm in famus_driver.'
          end if
          ! Step 1 : Count the number of entries
          num_bn_modes = 0
          do
             read(iunit2,*,iostat = ierr) mm, nn, bnc
             if (ierr .ne. 0) exit
             num_bn_modes = num_bn_modes+1
          end do
          close(iunit2)
          print *,'<----famus_driver found this many bnc modes: ', num_bn_modes
   
          ! Ready to create plasma boundary file
          famus_plasmaboundary_filename = trim(proc_string) // '.boundary'
          CALL safe_open(iunit, istat, famus_plasmaboundary_filename, 'replace', 'formatted')
          write(iunit,'(A)' ) "#Nfou Nfp  Nbnf"
          write(iunit,'(3I6)' ) num_vmec_modes, nfp_vmec, num_bn_modes
          write(iunit,'(A)' ) "#------- plasma boundary------"
          write(iunit,'(A)' ) "#  n   m   Rbc   Rbs    Zbc   Zbs"
          do imn = 1, num_vmec_modes
             tmm = INT(xm_vmec(imn))
             tnn = INT(xn_vmec(imn) / nfp_vmec)
             trmnc = rmnc_vmec(imn, num_vmec_surfs)
             tzmns = zmns_vmec(imn, num_vmec_surfs)
             if (lasym_vmec) then
                trmns = rmns_vmec(imn, num_vmec_surfs)
                tzmnc = zmnc_vmec(imn, num_vmec_surfs)
             else
                trmns = 0
                tzmnc = 0
             end if
             write(iunit,'(2I5, 4ES15.6)') tnn, tmm, trmnc, trmns, tzmnc, tzmns
          enddo
     
          write(iunit,'(A)' ) "#-------Bn harmonics----------"
          write(iunit,'(A)' ) "#  n  m  bnc   bns"
          if (num_bn_modes .gt. 0) then
             ! Open the bn file and transfer the data (NOTE: THE ORDER OF THE
             ! COLUMNS IS DIFFERENT BETWEEN BNORM AND FAMUS.)
             call safe_open(iunit2, ierr, trim(famus_bnorm_filename), 'old', 'formatted')
             if (ierr .ne. 0 ) then
                stop 'Unable to open bnorm in famus_driver.'
             end if
             do
                read(iunit2,*,iostat = ierr) mm, nn, bnc
                if (ierr .ne. 0) exit
                write(iunit,'(2I6, 2ES15.6)') nn, mm, bnc, 0.0 
             end do
             write(iunit, '(A)') ''
             close(iunit2)
          else
             write(iunit,'(2I6, 2ES15.6)') 0, 0, 0.0, 0.0
          endif
     
          close(iunit)


          write(6,*) '<----famus_driver: finished writing plasma boundary file: ', &
                      trim(famus_plasmaboundary_filename)

      end if

      print *, '<----famus_driver: famus_id=',famus_id,' waiting at mpi_barrier #1'
      CALL MPI_BARRIER(FAMUS_MPI_COMM_FAMUS,ierr_mpi)
      print *, '<----famus_driver: famus_id=',famus_id,' moving on from mpi_barrier #1'

! STEP 3: Initialize FAMUS
       write(6,*) '<----famus_driver: Step 3'
!        call sleep(3) 

       famus_call_from_ext_opt = .True.
!       famus_color = 0
!       CALL MPI_COMM_RANK(MPI_COMM_MYWORLD, famus_rankWorld, ierr_mpi)
!       CALL MPI_COMM_SPLIT(MPI_COMM_MYWORLD, famus_color, famus_rankWorld, famus_MPI_COMM_FAMUS, ierr_mpi)
       
       famus_ext = trim(proc_string)
       famus_input_surf = trim(famus_ext) // '.boundary'
       famus_input_coils = trim(famus_ext) // '.focus'

! This was performed above
!       call famus_initialize
!      call MPI_COMM_RANK( famus_MPI_COMM_FAMUS, famus_id, ierr )
!      call MPI_COMM_SIZE( famus_MPI_COMM_FAMUS, famus_ncpu, ierr )

! STEP 4: check_input
       write(6,*) '<----famus_driver: Step 4'
       
! what is 'myid' in famus?
       call check_input
      !CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_BARRIER(FAMUS_MPI_COMM_FAMUS,ierr_mpi)

       write(6,*) '<----famus_driver: Begin execution with famus_ncpu =', famus_ncpu

! STEP 5: handle FAMUS surface
       ! JCS: This probably only needs to be done if the MHD equilibrium has changed. 
       !  Need to check for *_oneiter types of equilibrium and handle appropriately.
       !  For now, this will be called every time. 

       ! fousurf: read plasma surface from filename stored in famus variable 'input_surf'
       print *, '<----famus_driver: myid=',myid,',famus_id=',famus_id,' is handling case_surface=',famus_case_surface,' famus_input_surf = ',trim(famus_input_surf)

       !call fousurf2
       call fousurf
      !CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_BARRIER(FAMUS_MPI_COMM_FAMUS,ierr_mpi)

!       select case( famus_case_surface )
!
!           case( 0 ) ; call fousurf   ! general format (VMEC-like) plasma boundary;
!           case( 1 ) ; call rdknot    ! knototran-like plasma boundary;
!          !case( 2 ) ; call readwout  ! read vmec output for plasma boundary and Boozer
!          !coordinates; for future;
!
!       end select

       write(6,*) '<----famus_driver: myid=',myid,',famus_id=',famus_id,' is handling case_coils'
       select case( famus_case_coils )

       ! JCS: For now, only enabling the settings specified in the FAMUS/examples/ellipse folder
       ! JCS: This can be replaced by assigning the correct variables. For now, this
       ! will read the files genereated above.
       ! JCS: Updating the .famus file to read from.
       ! For now, the 'fixed' coils file will not be modified.
          !case( 0 )   ; call coilpwl ! piece-wise linear; for future;
           case( 1 )   !;
                    write(6,*) '<----famus_driver: myid=',myid,' is calling rdcoils ', famus_input_coils

                     call rdcoils
                    ! CALL famus_rdcoils_mpi(famus_dipolecoil_filename, famus_ncoils, famus_momentq)
      !                  famus_input_coils = famusin_out_filename

       end select
      !CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_BARRIER(FAMUS_MPI_COMM_FAMUS,ierr_mpi)

       write(6,*) '<----famus_driver: myid=',myid,',famus_id=',famus_id,' is handling packeddof'
       call packdof(famus_xdof)  ! packdof in xdof array;
      !CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_BARRIER(FAMUS_MPI_COMM_FAMUS,ierr_mpi)

!       write(6,*) '<----famus_driver: myid=',myid,' is calling mpi_barrier (myworld)'
!       call MPI_BARRIER( MPI_COMM_MYWORLD, ierr )

!JCS  ! if (lsqmin) call leastsq   ! least-square minimization
       famus_case_optimize = 0
       famus_case_bnormal = 0
       print *,'<----famus_case_optimize=',famus_case_optimize
       if (famus_case_optimize /= 0) call solvers       ! call different solvers;
      !CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_BARRIER(FAMUS_MPI_COMM_FAMUS,ierr_mpi)

!       write(6,*) '<----famus_driver: myid=',myid,' is calling mpi_barrier (myworld)'
      ! call MPI_BARRIER( famus_MPI_COMM_MYWORLD, ierr )
!      CALL MPI_BARRIER(FAMUS_MPI_COMM_FAMUS, ierr_mpi)


  write(6,*) '<----famus_driver: myid=',myid,',famus_id=',famus_id,' is calling unpacking'
  call unpacking(famus_xdof)  ! unpack the optimized xdof array;
      !CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_BARRIER(FAMUS_MPI_COMM_FAMUS,ierr_mpi)
 !JCS should this be myid or famus_id?
  !    write(6, *) "-----------POST-PROCESSING myid=0-----------------------------------"
  !    if (famus_id == 0) write(6, *) "-----------POST-PROCESSING famus_id=0-----------------------------------"
      print *,  '<----famus_driver: post unpacking'
      print *,  '<----famus_driver: famus_id=',famus_id,' myid=',myid,' beginning diagnos'
       call diagnos
      CALL MPI_BARRIER(FAMUS_MPI_COMM_FAMUS,ierr_mpi)
      print *,  '<----famus_driver: famus_id=',famus_id,' myid=',myid,' beginning prepare_ind'
       call prepare_inductance()
      CALL MPI_BARRIER(FAMUS_MPI_COMM_FAMUS,ierr_mpi)
      print *,  '<----famus_driver: famus_id=',famus_id,' myid=',myid,' beginning famus_bnormal'
       call famus_bnormal(0)
      CALL MPI_BARRIER(FAMUS_MPI_COMM_FAMUS,ierr_mpi)
!      select case( famus_case_postproc )
!      case( 0 ) 
!      case( 1 ) ; call diagnos ; 
!      case( 2 ) ; call diagnos ; call specinp !; call saving 
!     !case( 2 ) ; call saving  ; call diagnos ; call wtmgrid  ! write mgrid file;
!      case( 3 ) ; call diagnos ; call poinplot ! Poincare plots; for future; 
!     ! case( 3 ) ;  call poinplot ! Poincare plots; for future; 
!      case( 4 ) ; call diagnos ; call boozmn ; call poinplot ! Last closed surface
!      case( 5 ) ; call diagnos ; call wtmgrid  ! write mgrid file
!      case(6)
!         call poinplot ; call wtmgrid 
!     !case( 4 ) ; call saving  ; call diagnos ; call resonant ! resonant harmonics analysis; for future; 
!      end select

! After famus is called, the value we want should be
! in the variable 'famus_bn_target'. 'Saving' can be skipped, unless this
! information is useful for debugging. 
! To do: Ensure that the output file has the correct extension/name.
!      call saving ! save all the outputs

      write(6,*) '<----famus_driver: myid=',myid,' is almost done!'
!      call MPI_BARRIER( MPI_COMM_MYWORLD, ierr )
!      call MPI_COMM_FREE( FAMUS_MPI_COMM_FAMUS, ierr_mpi)
!JCS  if(myid == 0) write(ounit, *) "-------------------------------------------------------------"
      CALL MPI_BARRIER(FAMUS_MPI_COMM_FAMUS,ierr_mpi)

!JCS  call MPI_FINALIZE( ierr )

  !TO DO: Check to see how cleanup handles any intermediate or output files.
  

      IF (lscreen) WRITE(6,'(a)') ' -------------------  FAMUS CALCULATION FINISHED  ---------------------'
!DEC$ ELSE
      STOP "Error! stellopt_famus_driver was called, but STELLOPT was compiled without FAMUS."
!DEC$ ENDIF
      RETURN
 

!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_famus_driver
