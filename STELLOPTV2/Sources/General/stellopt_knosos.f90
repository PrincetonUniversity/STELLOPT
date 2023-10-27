!-----------------------------------------------------------------------
!     Subroutine:    stellopt_knosos
!     Authors:       J.L. Velasco (joseluis.velasco@ciemat.es),
!                    building on stellopt_dkes by S. Lazerson
!     Date:          last modified, 30/06/2020
!     Description:   This subroutine calculates the diffusion
!                    coefficient.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_knosos(lscreen,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
        USE stellopt_runtime, ONLY:  proc_string, bigno
        USE equil_utils, ONLY: rho
        USE stellopt_targets, ONLY: sigma_knosos_1nu, sigma_knosos_snu, sigma_knosos_sbp, sigma_knosos_gmc, sigma_knosos_gma, &
             & sigma_knosos_qer, sigma_knosos_vb0, sigma_knosos_vbm, sigma_knosos_vbb, sigma_knosos_wbw, sigma_knosos_dbo, lbooz, nsd
!DEC$ IF DEFINED (KNOSOS_OPT)
        USE knosos_stellopt_mod
        !USE GLOBAL, ONLY: nerr
!DEC$ ENDIF
        USE safe_open_mod
        USE mpi_params
        USE mpi_inc

!-----------------------------------------------------------------------
!     Subroutine Parameters
!----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL, INTENT(in)    :: lscreen
      INTEGER, INTENT(inout) :: iflag
!-----------------------------------------------------------------------
!     Local Parameters
!----------------------------------------------------------------------
      INTEGER, PARAMETER :: nbx=11
      INTEGER, PARAMETER :: nsx=1000
      INTEGER nbb,ns,REGB(nbx), irho(nsd)
      REAL*8 dt,Ab(nbx),Zb(nbx),s(nsx),fracb(nbx)

!-----------------------------------------------------------------------
!     Local Variables
!----------------------------------------------------------------------
      INTEGER ik, jk, istat, w_u3, ierr, numprocs_local, jerr, mystart, myend
      CHARACTER(120) :: out_file
      CHARACTER :: dkes_input_file*64, temp_str*64
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: rank_local
      REAL*8, ALLOCATABLE :: s_kn(:)
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      IF (iflag < 0) RETURN
!DEC$ IF DEFINED (KNOSOS_OPT)
      IF (lscreen) WRITE(6,'(a)') ' ---------------------------    KNOSOS CALCULATION     -------------------------'

      ! Setup MPI, note this is really for PETSC we should probably have a way to create a PETSC subcommunictor
      CALL SET_KNOSOS_COMM(MPI_COMM_MYWORLD)
      
      ! Deallocate and allocate the saved data arrays and initialize
      IF (ALLOCATED(KNOSOS_1NU)) DEALLOCATE(KNOSOS_1NU)
      IF (ALLOCATED(KNOSOS_SNU)) DEALLOCATE(KNOSOS_SNU)
      IF (ALLOCATED(KNOSOS_SBP)) DEALLOCATE(KNOSOS_SBP)
      IF (ALLOCATED(KNOSOS_GMC)) DEALLOCATE(KNOSOS_GMC)
      IF (ALLOCATED(KNOSOS_GMA)) DEALLOCATE(KNOSOS_GMA)
      IF (ALLOCATED(KNOSOS_QER)) DEALLOCATE(KNOSOS_QER)
      IF (ALLOCATED(KNOSOS_VBM)) DEALLOCATE(KNOSOS_VBM)
      IF (ALLOCATED(KNOSOS_VB0)) DEALLOCATE(KNOSOS_VB0)
      IF (ALLOCATED(KNOSOS_VBB)) DEALLOCATE(KNOSOS_VBB)
      IF (ALLOCATED(KNOSOS_WBW)) DEALLOCATE(KNOSOS_WBW)
      IF (ALLOCATED(KNOSOS_DBO)) DEALLOCATE(KNOSOS_DBO)
      ALLOCATE(KNOSOS_1NU(nsd),KNOSOS_SNU(nsd),KNOSOS_SBP(nsd),KNOSOS_GMC(nsd),KNOSOS_GMA(nsd),&
           KNOSOS_QER(nsd),KNOSOS_VBM(nsd),KNOSOS_VB0(nsd),KNOSOS_VBB(nsd),KNOSOS_WBW(nsd),KNOSOS_DBO(nsd))
      KNOSOS_1NU=0.0; KNOSOS_SNU=0.0; KNOSOS_SBP=0.0; KNOSOS_GMC=0.0; KNOSOS_GMA=0.0;
      KNOSOS_QER=0.0; KNOSOS_VBM=0.0; KNOSOS_VB0=0.0; KNOSOS_VBB=0.0; KNOSOS_WBW=0.0; KNOSOS_DBO=0.0;
      KN_1NU=0.0; KN_SNU=0.0; KN_SBP=0.0; KN_GMC=0.0; KN_GMA=0.0; KN_QER=0.0; KN_VBM=0.0; KN_VB0=0.0; KN_VBB=0.0; KN_WBW=0.0; KN_DBO=0.0

      ! Now we count the number of surface we work on
      irho = -1
      ns=-1 ! not sure why
      DO ik=2,nsd
         IF(sigma_knosos_1nu(ik) >= bigno .and. sigma_knosos_snu(ik) >= bigno .and. sigma_knosos_sbp(ik) >= bigno .and. sigma_knosos_vb0(ik) >= bigno .and. &
            sigma_knosos_gmc(ik) >= bigno .and. sigma_knosos_gma(ik) >= bigno .and. sigma_knosos_qer(ik) >= bigno .and. sigma_knosos_vbm(ik) >= bigno .and. &
            sigma_knosos_vbb(ik) >= bigno .and. sigma_knosos_wbw(ik) >= bigno .and. sigma_knosos_dbo(ik) >= bigno) CYCLE
         ns=ns+1
         irho(ns+1) = ik
         IF(myworkid == ns ) WRITE(temp_str,'(A,I3.3)') '_s',ik
      END DO

      ! Read the input information
      KN_EXT = TRIM(proc_string)
      IF(ANY(sigma_knosos_1nu < bigno)) KN_STELLOPT(1)  = .TRUE.
      IF(ANY(sigma_knosos_snu < bigno)) KN_STELLOPT(2)  = .TRUE.
      IF(ANY(sigma_knosos_gmc < bigno)) KN_STELLOPT(4)  = .TRUE.
      IF(ANY(sigma_knosos_gma < bigno)) KN_STELLOPT(5)  = .TRUE.
      IF(ANY(sigma_knosos_dbo < bigno)) KN_STELLOPT(6)  = .TRUE.
      IF(ANY(sigma_knosos_vbm < bigno)) KN_STELLOPT(7)  = .TRUE.
      IF(ANY(sigma_knosos_vbb < bigno)) KN_STELLOPT(8)  = .TRUE.
      IF(ANY(sigma_knosos_wbw < bigno)) KN_STELLOPT(9)  = .TRUE.
      IF(ANY(sigma_knosos_vb0 < bigno)) KN_STELLOPT(10) = .TRUE.
      CALL READ_INPUT(dt,ns,s,nbb,Zb,Ab,regb,fracb)

      !Make an output file
      out_file ='kn_out.'//TRIM(proc_string)
      w_u3 =222
      istat=0
      IF (myworkid == master) call safe_open(w_u3,istat,out_file,'replace','formatted')
      CALL MPI_BCAST(istat,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
      IF (istat .ne. 0) THEN
         IF (myworkid == master) PRINT *,istat,out_file
         IF (myworkid == master) WRITE(6,*) 'Error opening KNOSOS output file:',TRIM(out_file),istat
         iflag=-1
         RETURN
      END IF

      ! For now we default nerr = 1
      !nerr = 1

      ! Here's where we handle PETSC
      !    This probably needs to be a bit more sophisticated
      !    For now one processor per run
      !IF (ALLOCATED(rank_local)) DEALLOCATE(rank_local)
      !ALLOCATE(rank_local(ns,nerr))
      !CALL MPI_INIT_PETSC(ns,nerr,rank,MPI_COMM_KNOSOS)

      ! Random number generator
      CALL INIT_RANDOMSEED(0)

      ! Intialize files
      CALL INIT_FILES()
      KN_DKESFILE='boozmn_'//TRIM(proc_string)//'.nc'
      KN_IFILE='input.'//TRIM(proc_string)

      ! OK now we handle dividing up the work
      ik = COUNT(irho > 0)
      CALL MPI_CALC_MYRANGE(MPI_COMM_MYWORLD, 1, ik, mystart, myend)

      ! Create s_kn
      IF (ALLOCATED(s_kn)) DEALLOCATE(s_kn)
      ALLOCATE(s_kn(nsd))
      IF(myworkid == master) s_kn=rho*rho
      CALL MPI_BCAST(s_kn,nsd,MPI_REAL8,master,MPI_COMM_MYWORLD,ierr_mpi)

      IF (lscreen) WRITE(6,'(A)') ' -- Intialized'

      ! Loop over problem
      DO ik = mystart,myend
            jk = irho(ik)
            KN_STELLOPT = .FALSE.
            IF(sigma_knosos_1nu(jk) < bigno) KN_STELLOPT(1)=.TRUE.
            IF(sigma_knosos_snu(jk) < bigno) KN_STELLOPT(2)=.TRUE.
            IF(sigma_knosos_gmc(jk) < bigno) KN_STELLOPT(4)=.TRUE.
            IF(sigma_knosos_gma(jk) < bigno) KN_STELLOPT(5)=.TRUE.
            IF(sigma_knosos_vbm(jk) < bigno) KN_STELLOPT(6)=.TRUE.
            IF(sigma_knosos_vbb(jk) < bigno) KN_STELLOPT(7)=.TRUE.
            IF(sigma_knosos_wbw(jk) < bigno) KN_STELLOPT(8)=.TRUE.
            IF(sigma_knosos_dbo(jk) < bigno) KN_STELLOPT(9)=.TRUE.
            IF(sigma_knosos_vb0(jk) < bigno) KN_STELLOPT(10)=.TRUE.
      IF (lscreen) WRITE(6,'(A)') ' -- B_FIELD'
            CALL READ_BFIELD(s_kn(jk))
            KNOSOS_VBM(jk)=KN_VBM
            KNOSOS_VB0(jk)=KN_VB0
            KNOSOS_VBB(jk)=KN_VBB
            KNOSOS_WBW(jk)=KN_WBW
            KNOSOS_DBO(jk)=KN_DBO
      IF (lscreen) WRITE(6,'(A)') ' -- CALC_DATABASE'
            CALL CALC_DATABASE(s_kn,jk,ns)
            KNOSOS_1NU(jk)=KN_1NU
            KNOSOS_SNU(jk)=KN_SNU
            KNOSOS_SBP(jk)=KN_SBP
            KNOSOS_GMC(jk)=KN_GMC
            KNOSOS_GMA(jk)=KN_GMA
            !nb=0.0;dnbdpsi=0.0/dpsidr;Tb=0.0;dTbdpsi=0.0/dpsidr:  !JLVG: to be done
            !CALL SOLVE_DKE_QN_AMB(itime,nbb,Zb,Ab,regb,s_kn(jk),nb,dnbdpsi,Tb,dTbdpsi,Epsi,Gb,Qb)
            KNOSOS_QER(jk)=KN_QER
            !KNOSOS_DJRDER(jk))=KN_DJRDER
            !KNOSOS_ERSHEAR(jk)=KN_ERSHEAR
      END DO
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,KNOSOS_1NU,nsd,MPI_REAL8,MPI_SUM,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,KNOSOS_SNU,nsd,MPI_REAL8,MPI_SUM,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,KNOSOS_SBP,nsd,MPI_REAL8,MPI_SUM,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,KNOSOS_GMC,nsd,MPI_REAL8,MPI_SUM,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,KNOSOS_GMA,nsd,MPI_REAL8,MPI_SUM,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,KNOSOS_QER,nsd,MPI_REAL8,MPI_SUM,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,KNOSOS_VBM,nsd,MPI_REAL8,MPI_SUM,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,KNOSOS_VB0,nsd,MPI_REAL8,MPI_SUM,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,KNOSOS_VBB,nsd,MPI_REAL8,MPI_SUM,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,KNOSOS_WBW,nsd,MPI_REAL8,MPI_SUM,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,KNOSOS_DBO,nsd,MPI_REAL8,MPI_SUM,MPI_COMM_MYWORLD,ierr_mpi)
!DEC$ ENDIF
      IF (myworkid == master) THEN
         DO ik=2,nsd
            IF(.not. lbooz(ik)) CYCLE
            IF(sigma_knosos_1nu(ik) >= bigno .and. sigma_knosos_snu(ik) >= bigno .and. sigma_knosos_sbp(ik) >= bigno .and. sigma_knosos_vb0(ik) >= bigno .and. &
               sigma_knosos_gmc(ik) >= bigno .and. sigma_knosos_gma(ik) >= bigno .and. sigma_knosos_qer(ik) >= bigno .and. sigma_knosos_vbm(ik) >= bigno .and. &
               sigma_knosos_vbb(ik) >= bigno .and. sigma_knosos_wbw(ik) >= bigno .and. sigma_knosos_dbo(ik) >= bigno) CYCLE
            WRITE(w_u3,'(1(1x,i8),20(1x,e17.10))') ik,KNOSOS_1NU(ik),&
                 & KNOSOS_SNU(ik),KNOSOS_SBP(ik),KNOSOS_GMC(ik),KNOSOS_GMA(ik),KNOSOS_QER(ik),KNOSOS_VB0(ik),KNOSOS_VBB(ik),KNOSOS_WBW(ik),KNOSOS_DBO(ik),KNOSOS_VBM(ik)
!            eff_ripple(ik)=KNOSOS_1NU(ik)
         END DO
         IF (lscreen) WRITE(6,'(2X,I8,1(2X,E17.10))') ik,KNOSOS_1NU(ik),&
              KNOSOS_SNU(ik),KNOSOS_SBP(ik),KNOSOS_GMC(ik),KNOSOS_GMA(ik),KNOSOS_QER(ik),KNOSOS_VB0(ik),KNOSOS_VBB(ik),KNOSOS_WBW(ik),KNOSOS_DBO(ik),KNOSOS_VBM(ik)
         CALL FLUSH(6)
         CLOSE(w_u3)
      END IF

      IF (ALLOCATED(s_kn)) DEALLOCATE(s_kn); !edi.sanchez@ciemat.es
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!KNOSOS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (lscreen) WRITE(6,'(a)') ' -----------------  KNOSOS CALCULATION (DONE) ----------------'
!DEC$ ELSE
      IF (lscreen) WRITE(6,'(a)') ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      IF (lscreen) WRITE(6,'(a)') ' !! KNOSOS based optimization not supported on your machine !!'
      IF (lscreen) WRITE(6,'(a)') ' !! Check for KNOSOS directory and KNOSOS_OPT in your makefile!!'
      IF (lscreen) WRITE(6,'(a)') ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!DEC$ ENDIF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_knosos
