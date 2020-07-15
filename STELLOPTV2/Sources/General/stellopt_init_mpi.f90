!-----------------------------------------------------------------------
!     Subroutine:    stellopt_init_mpi
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          06/04/2019
!     Description:   This subroutine initializes the vmec MPI 
!                    communicators.  The general idea is that first we
!                    figure out the shared memory groups.  From this 
!                    we define our MPI_COMM_MYWORLD and then from the 
!                    masters of MPI_COMM_MYWORLD we redefine
!                    MPI_COMM_STEL.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_init_mpi
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_input_mod, ONLY: noptimizers
      USE mpi_params
      USE mpi_inc
      USE safe_open_mod, ONLY: safe_open
!----------------------------------------------------------------------
!     Local Variables
!        nprocs_total :: total number of processors         
!----------------------------------------------------------------------
      INTEGER :: nprocs_total, vmajor, vminor, color, key, nshar, &
                 ngshar, comm_share, k
      INTEGER, EXTERNAL :: common_factor
      INTEGER :: iunit, ierr, tag, myid_world, myid_share, optimizer_color
      INTEGER, PARAMETER :: buffer_length = 1000
      CHARACTER(len=buffer_length) :: proc_assignments_string
      INTEGER :: mpi_status_local(MPI_STATUS_SIZE)
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
#if defined(MPI_OPT)
      ierr_mpi = 0
      ! Get total number of processors
      CALL MPI_COMM_SIZE( MPI_COMM_STEL, nprocs_total, ierr_mpi )
      CALL MPI_COMM_RANK( MPI_COMM_STEL, myid, ierr_mpi)
      myid_world = myid

      IF (myid == master) THEN
         iunit = 80
         ierr = 0
         CALL safe_open(iunit, ierr, "stellopt_mpi", "unknown", "formatted")

         CALL MPI_GET_VERSION(vmajor,vminor,ier)
         WRITE(6,*) '-----  MPI Params.   -----'
         WRITE(6,'(4X,A,I2,A,I2.2)') 'MPI Version: ',vmajor,'.',vminor
         WRITE(6,*)            '   Optimizers requested: ',noptimizers
         WRITE(iunit,"(A,I5)") '   Optimizers requested: ',noptimizers
         WRITE(6,*)            '   Number of Processors: ',nprocs_total
         WRITE(iunit,"(A,I5)") '   Number of Processors: ',nprocs_total
         CALL FLUSH(6)
      END IF

      ! See how many shared memory groups we have
      CALL MPI_COMM_SPLIT_TYPE( MPI_COMM_STEL, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, comm_share, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_init_mpi1',ierr_mpi)
      CALL MPI_COMM_RANK(comm_share, myworkid, ierr_mpi)
      CALL MPI_COMM_SIZE(comm_share, nshar, ierr_mpi )
      !CALL MPI_COMM_FREE(MPI_COMM_MYWORLD, ierr_mpi)
      myid_share = myworkid

      !JCS print *,'<----myid=',myid,'of ', nprocs_total,' total. myworkid=', &
      !JCS         myworkid, ' of size ', nshar
      ! Count number of shared memory groups
      ngshar = 0
      ! JCS Isn't there only 1 master, so only one cpu has ngshar  = 1?
      IF (myworkid == master) ngshar = 1
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, ngshar, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_STEL, ierr_mpi)
      IF (myid == master) THEN
         WRITE(6,*)            '   Shared memory groups: ',ngshar
         WRITE(iunit,"(A,I5)") '   Shared memory groups: ',ngshar
         WRITE(6,*)            '   Processors per group: ',nshar
         WRITE(iunit,"(A,I5)") '   Processors per group: ',nshar
         CALL FLUSH(6)
      END IF

      ! Catch default behavoir
      IF (noptimizers <= 0) noptimizers = nprocs_total
      IF (noptimizers > nprocs_total) THEN
         noptimizers = nprocs_total
         IF (myid == master) WRITE(iunit,"(A)") "Lowering noptimizers to nprocs_total"
      END IF

      ! Logic here
      ! NOPTIMIZERS <= NSHARED_GROUPS (we spread over groups nodes)
      ! NOPTIMIZERS > NSHARED_GROUPS (subdivide groups)
      IF (noptimizers <= ngshar) THEN
         IF (myid == master) WRITE(iunit,"(A)") "noptimizers <= ngshar."

         IF (MOD(ngshar,noptimizers)/=0) THEN ! we need to redefine NOPTIMIZERS
            noptimizers = common_factor(ngshar, noptimizers, 1)
         END IF 
         ! Destroy MPI_COMM_STELa
         !JCS print *,'<----Freeing mpi_comm_stell in init_mpi'
         CALL MPI_COMM_FREE(MPI_COMM_STEL, ierr_mpi)
         ! Make MPI_COMM_STEL out of just shared comm masters
         key = myid
         color = MPI_UNDEFINED
         IF (myworkid == master) color = 0
         CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, color, key, MPI_COMM_STEL, ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_init_mpi2',ierr_mpi)
         myid = -1
         IF (myworkid == master)THEN
            CALL MPI_COMM_RANK( MPI_COMM_STEL, myid, ierr_mpi )              ! MPI
            IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_init_mpi4',ierr_mpi)
            CALL MPI_COMM_SIZE( MPI_COMM_STEL, numprocs, ierr_mpi )          ! MPI
            IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_init_mpi5',ierr_mpi)
         END IF

         ! Now MPI_COMM_STEL is composed of the head nodes of each shared memory group
         ! We now create MPI_COMM_MYWORLD from bunches of these groups
         color = MPI_UNDEFINED
         IF (myworkid == master) color = MOD(myid,noptimizers)
         CALL MPI_BCAST(color, 1, MPI_INTEGER, master, comm_share, ierr_mpi)
         optimizer_color = color
         CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, color, key, MPI_COMM_MYWORLD, ierr_mpi)
         CALL MPI_COMM_RANK( MPI_COMM_MYWORLD, myworkid, ierr_mpi )              ! MPI
         CALL MPI_COMM_SIZE( MPI_COMM_MYWORLD, nshar, ierr_mpi )          ! MPI
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_init_mpi4',ierr_mpi)

         ! Now CREATE MPI_COMM_STEL from MPI_COMM_MYWORLD
         CALL MPI_COMM_FREE(MPI_COMM_STEL, ierr_mpi)
         color = MPI_UNDEFINED
         IF (myworkid == master) color=0
         CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, color, key, MPI_COMM_STEL, ierr_mpi)
         IF (myworkid == master)THEN
            CALL MPI_COMM_RANK( MPI_COMM_STEL, myid, ierr_mpi )              ! MPI
            IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_init_mpi4',ierr_mpi)
            CALL MPI_COMM_SIZE( MPI_COMM_STEL, numprocs, ierr_mpi )          ! MPI
            IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_init_mpi5',ierr_mpi)
         END IF

      ELSE ! Now handle the case noptimizers > ngshar:
         IF (myid == master) WRITE(iunit,"(A)") "noptimizers > ngshar."

         ! We subdivide the shared memory communicators
         k = MAX(noptimizers/ngshar,2) ! at least we need to divide each group by 2
         IF (myid == master) WRITE(iunit,"(A,I5,A,I5)") "Initial k:",k,", MOD(nshar,k):",MOD(nshar,k)
         IF (MOD(nshar,k)/=0) THEN ! we need to redefine k
            k = common_factor(nshar, k, 1)
         END IF 
         IF (myid == master) WRITE(iunit,"(A,I5)") "Final k:",k

         ! Subdivide the shared memory communicator according to k
         key = myid
         !color = MOD(myworkid,k)
         color = myworkid / (nshar / k) ! Note integer division here, in which FORTRAN will round down. No rounding should be needed for (nshar/k), but the preceding division may involve rounding down.
         optimizer_color = color
         CALL MPI_COMM_SPLIT(comm_share, color, key, MPI_COMM_MYWORLD, ierr_mpi)
         CALL MPI_COMM_RANK( MPI_COMM_MYWORLD, myworkid, ierr_mpi )
         CALL MPI_COMM_SIZE( MPI_COMM_MYWORLD, nshar, ierr_mpi )

         ! Now CREATE MPI_COMM_STEL from MPI_COMM_MYWORLD
         CALL MPI_COMM_FREE(MPI_COMM_STEL, ierr_mpi)
         color = MPI_UNDEFINED
         IF (myworkid == master) color=0
         CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, color, key, MPI_COMM_STEL, ierr_mpi)
         !CALL MPI_COMM_SPLIT(MPI_COMM_MYWORLD, color, key, MPI_COMM_STEL, ierr_mpi)
         IF (myworkid == master)THEN
            CALL MPI_COMM_RANK( MPI_COMM_STEL, myid, ierr_mpi )              ! MPI
            IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_init_mpi4',ierr_mpi)
            CALL MPI_COMM_SIZE( MPI_COMM_STEL, numprocs, ierr_mpi )          ! MPI
            IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_init_mpi5',ierr_mpi)
         END IF

      END IF

      ! Free the shared memory communicator
      CALL MPI_COMM_FREE(comm_share, ierr_mpi)


      IF (myid == master) THEN
         WRITE(6,*)            '  Workers per optimizer: ',nshar
         WRITE(iunit,"(A,I5)") '  Workers per optimizer: ',nshar
         WRITE(6,*)            '    Optimizers provided: ',numprocs
         WRITE(iunit,"(A,I5)") '    Optimizers provided: ',numprocs
         CALL FLUSH(6)
      END IF

      ! Record the MPI information for each processor.
      WRITE (proc_assignments_string, "(100(I5,A))") myid_world,",",myid_share,",",optimizer_color,",",myworkid
      IF (myid_world == master) THEN
         WRITE (iunit,*)
         WRITE (iunit, "(A)") "rank in MPI_COMM_WORLD, rank in comm_share, optimizer color, rank in MPI_COMM_MYWORLD"
         WRITE (iunit, "(A)") TRIM(proc_assignments_string)
         DO tag = 1,nprocs_total - 1
            CALL MPI_RECV(proc_assignments_string,buffer_length,MPI_CHAR,tag,tag,MPI_COMM_WORLD,mpi_status_local,ierr_mpi)
            WRITE(iunit,"(A)") TRIM(proc_assignments_string)
         END DO
         CLOSE(iunit)
      ELSE
         ! Every processor other than the world master executes this block
         tag = myid_world
         CALL MPI_SEND(proc_assignments_string,buffer_length,MPI_CHAR,0,tag,MPI_COMM_WORLD,ierr_mpi)
      END IF

#endif
     
      RETURN
      
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_init_mpi

