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
      USE stellopt_vars, my_mpol => mpol_fds, my_ntor => ntor_fds
      USE equil_utils
      USE neswrite, ONLY: coil_separation

      USE famus_globals, ONLY: nfp_famus => Nfp, IsQuiet_famus => IsQuiet, &
           famus_id => myid, &
           famus_Nzeta => Nzeta, &
           famus_Nteta => Nteta, famus_ncpu => ncpu, famus_xdof => xdof, &
                        famus_MPI_COMM_MASTERS => MPI_COMM_MASTERS, &
                        famus_MPI_COMM_MYWORLD => MPI_COMM_MYWORLD, &
                        famus_MPI_COMM_WORKERS => MPI_COMM_WORKERS, &
                        famus_MPI_COMM_FAMUS => MPI_COMM_FAMUS, &
                        famus_called_from_ext_opt => called_from_ext_opt, &
                        famus_case_surface => case_surface, &
                        famus_case_coils => case_coils, &
                        famus_input_coils => input_coils, &
                        famus_case_postproc => case_postproc, &
                        famus_case_optimize => case_optimize
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
      INTEGER :: istat, iunit = 112, m, n, ii, imn, nummodes1, nummodes2
      INTEGER :: secs, mins, hrs, ierr
!JCS      REAL    :: tstart, tfinish ! local variables
      CHARACTER(400) :: famusin_out_filename, famus_bnorm_filename


!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
         WRITE(6,'(a)') ' -------------  FAMUS CALCULATION 1 ---------'
      IF (iflag < 0) RETURN
      lscreen = .true.
      IF (lscreen) then
         WRITE(6,'(a)') ' -------------  FAMUS CALCULATION 2 ---------'
      ENDIF
       WRITE(6,'(a,a)') '<----stellopt_famus_driver: proc_string=', trim(proc_string), &
                      ' file_str=', trim(file_str)
       WRITE(6,'(a,i4.2)') '<-------------stellopt_famus_driver: iflag=', iflag

      ! Suppress FAMUS stdout if needed
      IsQuiet_famus = 1
      if (lscreen) IsQuiet_famus = -1

      write(6,*) "<----proc_str=",trim(proc_string)," file_str=",file_str

! STEP 1.
       write(6,*) '<----famus_driver: Step 1'

      ! Run bnorm if required ( JCS - is this handled correctly in the MPI
      ! environment?)
      !if (load_bnorm_famus) then
      if (.True.) then
         write(6,*) '<----stellopt_famus_driver: calling stellopt_bnorm(', &
                    file_str, ', ', lscreen,')'
         coil_separation = 0.1 ! JCS - Does this matter?
         ! Run BNORM code
         call stellopt_bnorm(file_str,lscreen)
         famus_bnorm_filename = 'bnorm.' // TRIM(file_str)
      end if

      IF (lscreen) WRITE(6,'(a,a)') '<----stellopt_famus_driver: proc_string=', trim(proc_string)
      !wout_filename = 'wout_'//TRIM(proc_string)//'.nc'

      ! STELLOPT (via lmdif->stellopt_fcn or similar) will modifiy the value of
      ! the boundary coefficients. Here, the value of the FAMUS variables
      ! are loaded with the new values, the correct mode of operation is
      ! determiend, control variables are set, and the regcoil-functions
      ! are called

! STEP 2. Write out the winding surface coefficients so that FAMUS can read them
       write(6,*) '<----famus_driver: Step 2'
     
      ! Loop over all of the spectral components of the winding surface
      ! and update the foucs_ds_(rz)bound_(cs)

      IF ((ANY(lfamus_ds_rbound_s_opt)) .or. (ANY(lfamus_ds_rbound_c_opt)) .or. &
          (ANY(lfamus_ds_zbound_s_opt)) .or. (ANY(lfamus_ds_zbound_c_opt)) ) THEN 
         nummodes1 = 0
         DO m = -my_mpol, my_mpol
             DO n = -my_ntor, my_ntor
                IF ( (famus_ds_rbound_c(m,n) .ne. 0) .or. &
                     (famus_ds_rbound_s(m,n) .ne. 0) .or. &
                     (famus_ds_zbound_c(m,n) .ne. 0) .or. &
                     (famus_ds_zbound_s(m,n) .ne. 0) ) THEN
                   nummodes1 = nummodes1 + 1
                END IF
             END do
         END do


         ! The famus variables are written to file, first (read in by FAMUS later).
         ! During 'cleanup', the new minimum will be copied/moved
         ! to have the iteration number as a suffix
         famusin_out_filename = trim(proc_string) // '.famus'
         CALL safe_open(iunit, istat, famusin_out_filename, 'replace', 'formatted')
         write (iunit, '(a)') ' # Total number of coils,  momentq'
         write (iunit, '(I6, (a), I6)') famus_Nzeta*famus_Nteta, ',', 1
         write (iunit, '(a)') '#coiltype, symmetry,  coilname,  ox,  oy,  oz,  Ic,  M_0,  pho,  Lc, mp,  mt'

         !JCS: This is the incorrect range for m and n.
         DO m = -my_mpol, my_mpol
             DO n = -my_ntor, my_ntor
                if ( (famus_ds_rbound_c(m,n) .ne. 0) .or. &
                     (famus_ds_rbound_s(m,n) .ne. 0) .or. &
                     (famus_ds_zbound_c(m,n) .ne. 0) .or. &
                     (famus_ds_zbound_s(m,n) .ne. 0) ) THEN
                   ! These are written in the same order as in FAMUS  (see blah.f90)
                   ! file: M N RC ZS RS ZC
                   write(iunit,*) m, n, &
                        famus_ds_rbound_c(m,n), famus_ds_zbound_s(m,n), &
                        famus_ds_rbound_s(m,n), famus_ds_zbound_c(m,n)
                        !rc_rmnc_stellopt(m,n), rc_zmns_stellopt(m,n), &
                        !rc_rmns_stellopt(m,n), rc_zmnc_stellopt(m,n)
                END IF
              END DO
          END DO
          CLOSE(iunit)

       ! 
       write(6,*) '<----famus_driver: Step 2.5'

! STEP 2.5: Assign the FAMUS variables, or use the FAMUS functionality to the
! same?

          ! JCS: To do: Assign the famus variables that are read in during 'rdcoils'
          !DO imn = 1, mnmax_coil
          !   m = xm_coil(imn)
          !   n = xn_coil(imn)/(-nfp_regcoil)
          !   IF (m < -my_mpol .or. m > my_mpol .or. n < -my_ntor .or. n > my_ntor) THEN
          !      WRITE(6,*) "Error! (m,n) in regcoil coil surface exceeds mpol_rcws or ntor_rcws."
          !      STOP
          !   END IF
          !   rmnc_coil(imn) = regcoil_rcws_rbound_c(m,n)
          !   rmns_coil(imn) = regcoil_rcws_rbound_s(m,n)
          !   zmnc_coil(imn) = regcoil_rcws_zbound_c(m,n)
          !   zmns_coil(imn) = regcoil_rcws_zbound_s(m,n)
          !END DO
      END IF

      ! This should *almost* be a duplicate of the main code from
      ! focus.f90


! STEP 3: Initialize FAMUS
       write(6,*) '<----famus_driver: Step 3'

       famus_called_from_ext_opt = .True.
       call famus_initialize

! STEP 4: check_input
       write(6,*) '<----famus_driver: Step 4'
       call check_input

       write(6,*) '<----famus_driver: Begin execution with famus_ncpu =', famus_ncpu

! STEP 5: handle FAMUS surface
       ! JCS: This probably only needs to be done if the MHD equilibrium has changed. 
       !  Need to check for *_oneiter types of equilibrium and handle appropriately.
       !  For now, this will be called every time. 



       ! fousurf: read plasma surface from filename stored in famus variable 'input_surf'
       write(6,*) '<----famus_driver: myid=',myid,' is handling case_surface'

       select case( famus_case_surface )

           case( 0 ) ; call fousurf   ! general format (VMEC-like) plasma boundary;
           case( 1 ) ; call rdknot    ! knototran-like plasma boundary;
          !case( 2 ) ; call readwout  ! read vmec output for plasma boundary and Boozer
          !coordinates; for future;

       end select

       write(6,*) '<----famus_driver: myid=',myid,' is handling case_coils'
       select case( famus_case_coils )

       ! JCS: For now, only enabling the settings specified in the FAMUS/examples/ellipse folder
       ! JCS: This can be replaced by assigning the correct variables. For now, this
       ! will read the files genereated above.
       ! JCS: Updating the .famus file to read from.
          !case( 0 )   ; call coilpwl ! piece-wise linear; for future;
           case( 1 )   ; call rdcoils
                        famus_input_coils = famusin_out_filename
                        write(6,*) '<----famus_driver: myid=',myid,' is calling rdcoils'

       end select

       write(6,*) '<----famus_driver: myid=',myid,' is calling packdof'
       call packdof(famus_xdof)  ! packdof in xdof array;

       write(6,*) '<----famus_driver: myid=',myid,' is calling mpi_barrier (myworld)'
       call MPI_BARRIER( MPI_COMM_MYWORLD, ierr )

!JCS  ! if (lsqmin) call leastsq   ! least-square minimization

       if (famus_case_optimize /= 0) call solvers       ! call different solvers;

       write(6,*) '<----famus_driver: myid=',myid,' is calling mpi_barrier (myworld)'
       call MPI_BARRIER( MPI_COMM_MYWORLD, ierr )


  write(6,*) '<----famus_driver: myid=',myid,' is calling unpacking'
  call unpacking(famus_xdof)  ! unpack the optimized xdof array;

      !if (myid == 0) write(ounit, *) "-----------POST-PROCESSING-----------------------------------"
      write(6, *) "<----famus_driver: myid=",myid," beginning POST-PROCESSING"

      select case( famus_case_postproc )
      case( 0 ) 
      case( 1 ) ; call diagnos ; 
      case( 2 ) ; call diagnos ; call specinp !; call saving 
     !case( 2 ) ; call saving  ; call diagnos ; call wtmgrid  ! write mgrid file;
      case( 3 ) ; call diagnos ; call poinplot ! Poincare plots; for future; 
     ! case( 3 ) ;  call poinplot ! Poincare plots; for future; 
      case( 4 ) ; call diagnos ; call boozmn ; call poinplot ! Last closed surface
      case( 5 ) ; call diagnos ; call wtmgrid  ! write mgrid file
      case(6)
         call poinplot ; call wtmgrid 
     !case( 4 ) ; call saving  ; call diagnos ; call resonant ! resonant harmonics analysis; for future; 
      end select

! After famus is called, the value we want should be
! in the variable 'famus_bn_target'. 'Saving' can be skipped, unless this
! information is useful for debugging. 
! To do: Ensure that the output file has the correct extension/name.
      call saving ! save all the outputs

      write(6,*) '<----famus_driver: myid=',myid,' is calling mpi_barrier'
      call MPI_BARRIER( MPI_COMM_MYWORLD, ierr )

!JCS  if(myid == 0) write(ounit, *) "-------------------------------------------------------------"

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
