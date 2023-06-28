!-----------------------------------------------------------------------
!     Subroutine:    thrift_nclass
!     Authors:       S. Lazerson
!     Date:          11/15/2022
!     Description:   This subroutine handles running NCLASS
!-----------------------------------------------------------------------
      SUBROUTINE thrift_nclass(lscreen,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE thrift_runtime
      USE thrift_vars
      USE thrift_profiles_mod
      USE read_wout_mod
      USE nclass_mod, ONLY: NCLASS
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
      INTEGER ::  k_conf, k_order, k_potato,m_i,m_z
      REAL(rprec) :: dense_min, cpot_b, cpot_l
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: b2av, b2av_inv, fhat, &
                  ftrap, gradrhobm2, Er, gradEr, ngradtheta, Aion
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: fm, gradTi, Ti
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: ni, gradpi
      REAL(rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: fexi

      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (lscreen) WRITE(6,'(a)') ' ---------------------------  NCLASS -------------------------'
      
      ! First setup the grid of quantities based on equilibrium
      IF (lvmec .and. myworkid==master) THEN
         ! Integer
         k_conf = 1
         k_order = 2
         k_potato = 0
         m_i      = nion_prof
         m_z      = MAXVAL(Zatom_prof)
         ! Float
         dense_min = 1.0E10
         c_potb = ABS(el0*Bt0*iotaf(1)*iotaf(1)/2) !elongation thing
         c_potl = Rmajor/iotaf(1)
         ! Float array
         ALLOCATE(b2av(nrho), b2av_inv(nrho), fhat(nrho), ftrap(nrho),&
            gradrhobm2(nrho), Er(nrho), gradEr(nrho), ngradtheta(nrho),&
            Aion(nion_prof))
         Er     = 1.0
         gradEr = 0.0
         


      END IF

#if defined(MPI_OPT)
      ! Now broadcast the inputs
      CALL MPI_BARRIER(MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(k_conf,    1, MPI_INTEGER, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(k_order,   1, MPI_INTEGER, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(k_potato,  1, MPI_INTEGER, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(m_i,       1, MPI_INTEGER, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(m_z,       1, MPI_INTEGER, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(dense_min, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(c_potb,    1, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(c_potl,    1, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(b2av,       nrho, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(b2av_inv,   nrho, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(fhat,       nrho, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(ftrap,      nrho, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(gradrhobm2, nrho, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(Er,         nrho, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(gradEr,     nrho, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(ngradtheta, nrho, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(Aion,  nion_prof, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
#endif

      ! Now loop over surfaces
      CALL MPI_CALC_MYRANGE(MPI_COMM_MYWORLD,1,nrho,mystart,myend)
      DO i = mystart, myend
         CALL NCLASS(k_conf, k_order, k_potato, m_i, m_z, &
                     dense_min, c_potb, c_potl, &
                     b2av(i), b2av_inv(i), fhat(i)
      END DO

      ! Free helpers
      DEALLOCATE(b2av, b2ab_inv, fhat, ftrap, gradrhobm2, Er, gradEr, ngradtheta, Aion)

      ! Now Reduce to master so they have the full arrays

      IF (lscreen) WRITE(6,'(a)') ' -------------------------  NCLASS DONE  ----------------------'
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_nclass
