!-----------------------------------------------------------------------
!     Subroutine:    stellopt_focus_driver
!     Authors:       J.C.Schmitt (Auburn/PPPL) jcschmitt@auburn.edu
!     Date:          2020
!     Description:   This subroutine calls the coil/permanent magnet
!                    optimization code, FOCUS.
!                    
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_focus_driver(file_str, lscreen, iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_input_mod
      USE stellopt_vars, my_mpol => mpol_fds, my_ntor => ntor_fds
      USE equil_utils
      USE neswrite, ONLY: coil_separation

!DEC$ IF DEFINED (FOCUS)
      USE focus_globals, ONLY: nfp_focus => Nfp, verbose_focus => verbose, &
           load_bnorm_focus => load_bnorm, focus_id => myid, &
           focus_bnorm_filename => bnorm_filename, focus_Nzeta => Nzeta, &
           focus_Nteta => Nteta, focus_ncpu => ncpu, focus_xdof => xdof, &
           input_coils

             
!DEC$ ENDIF
      USE mpi_params
      USE mpi_inc

!-----------------------------------------------------------------------
!     Subroutine Parameters
!        iflag         Error flag
!----------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(256), INTENT(inout)    :: file_str
      LOGICAL, INTENT(inout)        :: lscreen
      INTEGER, INTENT(inout) :: iflag

!-----------------------------------------------------------------------
!     Local Variables
!        istat         Error status
!        iunit         File unit number
      ! FOR FOCUS
      INTEGER :: istat, iunit = 112, m, n, ii, imn, nummodes1, nummodes2
      INTEGER :: secs, mins, hrs
!JCS      REAL    :: tstart, tfinish ! local variables
      CHARACTER(400) :: focusin_out_filename


!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
!      IF (iflag < 0) RETURN
      lscreen = .true.
      IF (lscreen) then
         WRITE(6,'(a)') ' -------------  FOCUS CALCULATION  ---------'
      ENDIF
       WRITE(6,'(a,a)') '<----stellopt_focus_driver: proc_string=', trim(proc_string)
       WRITE(6,'(a,i4.2)') '<-------------stellopt_focus_driver: iflag=', iflag
!DEC$ IF DEFINED (FOCUS)
      verbose_focus = lscreen ! Suppress FOCUS stdout if needed

      write(6,*) "<----proc_str=",trim(proc_string)," file_str=",file_str

      ! Run bnorm if required
      if (load_bnorm_focus) then
         write(6,*) '<----stellopt_focus_driver: calling stellopt_bnorm(', &
                    file_str, ', ', lscreen,')'
         coil_separation = 0.1 ! JCS - Does this matter?
         ! Run BNORM code
         call stellopt_bnorm(file_str,lscreen)
         focus_bnorm_filename = 'bnorm.' // TRIM(file_str)
      end if

      IF (lscreen) WRITE(6,'(a,a)') '<----stellopt_focus_driver: proc_string=', trim(proc_string)
!      wout_filename = 'wout_'//TRIM(proc_string)//'.nc'
      ! STELLOPT (via lmdif->stellopt_fcn or similar) will modifiy the value of
      ! the boundary coefficients. Here, the value of the FOUCS variables
      ! are loaded with the new values, the correct mode of operation is
      ! determiend, control variables are set, and the regcoil-functions
      ! are called
     
      ! Loop over all of the spectral components of the winding surface
      ! and update the foucs_ds_(rz)bound_(cs)

      IF ((ANY(lfocus_ds_rbound_s_opt)) .or. (ANY(lfocus_ds_rbound_c_opt)) .or. &
          (ANY(lfocus_ds_zbound_s_opt)) .or. (ANY(lfocus_ds_zbound_c_opt)) ) THEN 
         nummodes1 = 0
         DO m = -my_mpol, my_mpol
             DO n = -my_ntor, my_ntor
                IF ( (focus_ds_rbound_c(m,n) .ne. 0) .or. &
                     (focus_ds_rbound_s(m,n) .ne. 0) .or. &
                     (focus_ds_zbound_c(m,n) .ne. 0) .or. &
                     (focus_ds_zbound_s(m,n) .ne. 0) ) THEN
                   nummodes1 = nummodes1 + 1
                END IF
             END do
         END do


         ! The focus variables are written to file, first (not acutally used by
         ! FOCUS). But, during 'cleanup', the new minimum will be copied/moved
         ! to have the iteration number as a suffix
         focusin_out_filename = trim(proc_string) // '.focus'
         CALL safe_open(iunit, istat, focusin_out_filename, 'replace', 'formatted')
         write (iunit, '(a)') ' # Total number of coils,  momentq'
         write (iunit, '(I6, (a), I6)') focus_Nzeta*focus_Nteta, ',', 1
         write (iunit, '(a)') '#coiltype, symmetry,  coilname,  ox,  oy,  oz,  Ic,  M_0,  pho,  Lc, mp,  mt'

         !JCS: This is the incorrect range for m and n.
         DO m = -my_mpol, my_mpol
             DO n = -my_ntor, my_ntor
                if ( (focus_ds_rbound_c(m,n) .ne. 0) .or. &
                     (focus_ds_rbound_s(m,n) .ne. 0) .or. &
                     (focus_ds_zbound_c(m,n) .ne. 0) .or. &
                     (focus_ds_zbound_s(m,n) .ne. 0) ) THEN
                   ! These are written in the same order as in FOCUS  (see blah.f90)
                   ! file: M N RC ZS RS ZC
                   write(iunit,*) m, n, &
                        focus_ds_rbound_c(m,n), focus_ds_zbound_s(m,n), &
                        focus_ds_rbound_s(m,n), focus_ds_zbound_c(m,n)
                        !rc_rmnc_stellopt(m,n), rc_zmns_stellopt(m,n), &
                        !rc_rmns_stellopt(m,n), rc_zmnc_stellopt(m,n)
                END IF
              END DO
          END DO
          CLOSE(iunit)
 
          ! JCS: To do: Assign the focus variables that are read in during 'rdcoils'
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

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!JCS  myid = 0 ; ncpu = 1

!JCS  ! MPI initialize
!JCS  call MPI_INIT( ierr )
!JCS  call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
!JCS  call MPI_COMM_SIZE( MPI_COMM_WORLD, ncpu, ierr )

!JCS  tstart =  MPI_WTIME()
!JCS  if(myid == 0) write(ounit, *) "---------------------  FOCUS ", version, "------------------------------"
!JCS  if(myid == 0) write(ounit,'("focus   : Begin execution with ncpu =",i5)') ncpu
  write(6,*) '<----focus_driver: Begin execution with focus_ncpu =', focus_ncpu
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! inital was called in stellopt_input_mod
!  call initial ! read input namelist and broadcast;

! JCS: For now, only enabling the settings specified in the FOCUS/examples/ellipse folder
!JCS  select case( case_surface )
!JCS  case( 0 ) ; call fousurf   ! general format (VMEC-like) plasma boundary;
!JCS  case( 1 ) ; call rdknot    ! knototran-like plasma boundary;
!JCS !case( 2 ) ; call readwout  ! read vmec output for plasma boundary and Boozer coordinates; for future;
!JCS  end select

  ! JCS: This probably only needs to be done if the MHD equilibrium has changed. 
  !  Need to check for *_oneiter types of equilibrium and handle appropriately.
  !  For now, this will be called every time. 
  write(6,*) '<----focus_driver: myid=',myid,' is calling fousurf'
  call fousurf

! JCS: For now, only enabling the settings specified in the FOCUS/examples/ellipse folder
!JCS  select case( case_coils )
!JCS !case( 0 )   ; call coilpwl ! piece-wise linear; for future;
!JCS  case( 1 )   ; call rdcoils
!JCS  end select

  ! JCS: This can be replaced by assigning the correct variables. For now, this
  ! will read the files genereated above.
  write(6,*) '<----focus_driver: myid=',myid,' is calling rdcoils'
  ! JCS: Updating the .focus file to read from.
  input_coils = focusin_out_filename
  call rdcoils

  write(6,*) '<----focus_driver: myid=',myid,' is calling packdof'
  call packdof(focus_xdof)  ! packdof in xdof array;

!JCS  write(6,*) '<----focus_driver: myid=',myid,' is calling mpi_barrier'
!JCS  call MPI_BARRIER( MPI_COMM_WORLD, ierr )

!JCS  tfinish = MPI_Wtime()
!JCS  time_initialize = tfinish - tstart
!JCS  if( myid  ==  0 ) write(ounit, '(A, ES12.5, A3)') "focus   : Initialization took ", time_initialize," S."
!JCS  write(6, '(A, I5, A, ES12.5, A3)') "<----focus_driver: myid=", myid, " Initialization took ", time_initialize," S."

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!JCS  ! if (lsqmin) call leastsq   ! least-square minimization

!JCS  if (case_optimize /= 0) call solvers       ! call different solvers;

!JCS  call MPI_BARRIER( MPI_COMM_WORLD, ierr )

!JCS  tstart = MPI_Wtime()
!JCS  time_optimize = tstart - tfinish
!JCS  if( myid  ==  0 ) then
!JCS     secs = int(time_optimize)
!JCS     hrs = secs/(60*60)
!JCS     mins = (secs-hrs*60*60)/60
!JCS    secs = secs-hrs*60*60-mins*60
!JCS    if(hrs>0)then
!JCS        write(ounit, '(A, 3(I6, A3))') "focus   : Optimization took ",hrs," H ", mins," M ",secs," S."
!JCS    elseif(mins>0)then
!JCS        write(ounit, '(A, 2(I6, A3))') "focus   : Optimization took ", mins," M ",secs," S."
!JCS    else
!JCS        write(ounit, '(A, ES12.5, A3)') "focus   : Optimization took ", time_optimize," S."
!JCS    endif
!JCS endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  write(6,*) '<----focus_driver: myid=',myid,' is calling unpacking'
  call unpacking(focus_xdof)  ! unpack the optimized xdof array;

!JCS  if (myid == 0) write(ounit, *) "-----------POST-PROCESSING-----------------------------------"
!JCS  write(6, *) "<----focus_driver: myid=",myid," beginning POST-PROCESSING"

!JCS: For now, skipping post processing. 
!JCS  select case( case_postproc )
!JCS  case( 0 ) 
!JCS  case( 1 ) ; call diagnos ; 
!JCS  case( 2 ) ; call diagnos ; call specinp !; call saving 
!JCS !case( 2 ) ; call saving  ; call diagnos ; call wtmgrid  ! write mgrid file;
!JCS  case( 3 ) ; call diagnos ; call poinplot ! Poincare plots; for future; 
!JCS ! case( 3 ) ;  call poinplot ! Poincare plots; for future; 
!JCS  case( 4 ) ; call diagnos ; call boozmn ; call poinplot ! Last closed surface
!JCS  case( 5 ) ; call diagnos ; call wtmgrid  ! write mgrid file
!JCS  case(6)
!JCS     call poinplot ; call wtmgrid 
!JCS !case( 4 ) ; call saving  ; call diagnos ; call resonant ! resonant harmonics analysis; for future; 
!JCS  end select

! After focus is called, the value we want should be
! in the variable 'focus_bn_target'. 'Saving' can be skipped, unless this
! information is useful for debugging. 
! To do: Ensure that the output file has the correct extension/name.
!JCS  call saving ! save all the outputs

!JCS  write(6,*) '<----focus_driver: myid=',myid,' is calling mpi_barrier'
!JCS  call MPI_BARRIER( MPI_COMM_WORLD, ierr )

!JCS  tfinish = MPI_Wtime()
!JCS  time_postproc = tfinish - tstart
!JCS if( myid  ==  0 )then
!JCS    secs = int(time_postproc)
!JCS    hrs = secs/(60*60)
!JCS    mins = (secs-hrs*60*60)/60
!JCS    secs = secs-hrs*60*60-mins*60
!JCS    if(hrs>0)then
!JCS        write(ounit, '(A, 3(I6, A3))') "focus   : Post-processing took ",hrs," H ", mins," M ",secs," S."
!JCS    elseif(mins>0)then
!JCS        write(ounit, '(A, 2(I6, A3))') "focus   : Post-processing took ", mins," M ",secs," S."
!JCS    else
!JCS        write(ounit, '(A, ES12.5, A3)') "focus   : Post-processing took ", time_postproc," S."
!JCS    endif
!JCS endif

!JCS  if(myid == 0) write(ounit, *) "-------------------------------------------------------------"

!JCS  call MPI_FINALIZE( ierr )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

      !=================DEBUG SECTION==========================
      !  output_filename = 'focus_out.'// TRIM(proc_string)//'.nc'
      ! For debugging it can be useful to write the focus output file.
      !=================END OF DEBUG SECTION==========================

  !TO DO: Check to see how cleanup handles any intermediate or output files.
  

      IF (lscreen) WRITE(6,'(a)') ' -------------------  FOCUS CALCULATION FINISHED  ---------------------'
!DEC$ ENDIF
      RETURN
 

!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_focus_driver
