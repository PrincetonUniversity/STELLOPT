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
             & sigma_knosos_qer, sigma_knosos_vbt, sigma_knosos_vbb, sigma_knosos_wbw, sigma_knosos_dbo, lbooz, nsd
!DEC$ IF DEFINED (KNOSOS_OPT)
        USE knosos_stellopt_mod
!DEC$ ENDIF
        USE safe_open_mod
        USE equil_vals, ONLY: eff_ripple
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
      INTEGER nbb,ns,REGB(nbx)
      REAL*8 Ab(nbx),Zb(nbx),s(nsx),Zeff

!-----------------------------------------------------------------------
!     Local Variables
!----------------------------------------------------------------------
      INTEGER ik, jk, istat, w_u3, ierr
      CHARACTER(120) :: out_file
      CHARACTER :: dkes_input_file*64, temp_str*64
      REAL*8, ALLOCATABLE :: s_kn(:)
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      IF (lscreen) WRITE(6,*) 'KNOSOS',myworkid,iflag
      IF (iflag < 0) RETURN
!DEC$ IF DEFINED (KNOSOS_OPT)
      IF (lscreen) WRITE(6,'(a)') ' ---------------------------    KNOSOS CALCULATION     -------------------------'
      ! Enter the main loop
      IF (ALLOCATED(KNOSOS_1NU)) DEALLOCATE(KNOSOS_1NU)
      IF (ALLOCATED(KNOSOS_SNU)) DEALLOCATE(KNOSOS_SNU)
      IF (ALLOCATED(KNOSOS_SBP)) DEALLOCATE(KNOSOS_SBP)
      IF (ALLOCATED(KNOSOS_GMC)) DEALLOCATE(KNOSOS_GMC)
      IF (ALLOCATED(KNOSOS_GMA)) DEALLOCATE(KNOSOS_GMA)
      IF (ALLOCATED(KNOSOS_QER)) DEALLOCATE(KNOSOS_QER)
      IF (ALLOCATED(KNOSOS_VBT)) DEALLOCATE(KNOSOS_VBT)
      IF (ALLOCATED(KNOSOS_VBB)) DEALLOCATE(KNOSOS_VBB)
      IF (ALLOCATED(KNOSOS_WBW)) DEALLOCATE(KNOSOS_WBW)
      IF (ALLOCATED(KNOSOS_DBO)) DEALLOCATE(KNOSOS_DBO)
      ALLOCATE(KNOSOS_1NU(nsd),KNOSOS_SNU(nsd),KNOSOS_SBP(nsd),KNOSOS_GMC(nsd),KNOSOS_GMA(nsd),&
           KNOSOS_QER(nsd),KNOSOS_VBT(nsd),KNOSOS_VBB(nsd),KNOSOS_WBW(nsd),KNOSOS_DBO(nsd))

      KNOSOS_1NU=0.0; KNOSOS_SNU=0.0; KNOSOS_SBP=0.0; KNOSOS_GMC=0.0; KNOSOS_GMA=0.0;
      KNOSOS_QER=0.0; KNOSOS_VBT=0.0; KNOSOS_VBB=0.0; KNOSOS_WBW=0.0; KNOSOS_DBO=0.0;
      KN_1NU=0.0; KN_SNU=0.0; KN_SBP=0.0; KN_GMC=0.0; KN_GMA=0.0; KN_QER=0.0; KN_VBT=0.0; KN_VBB=0.0; KN_WBW=0.0; KN_DBO=0.0

      ns=-1
      DO ik=2,nsd
         IF(.not. lbooz(ik)) CYCLE
         IF(sigma_knosos_1nu(ik) >= bigno .and. sigma_knosos_snu(ik) >= bigno .and. sigma_knosos_sbp(ik) >= bigno .and. &
            sigma_knosos_gmc(ik) >= bigno .and. sigma_knosos_gma(ik) >= bigno .and. sigma_knosos_qer(ik) >= bigno .and. sigma_knosos_vbt(ik) >= bigno .and. &
            sigma_knosos_vbb(ik) >= bigno .and. sigma_knosos_wbw(ik) >= bigno .and. sigma_knosos_dbo(ik) >= bigno) CYCLE
         ns=ns+1
         IF(myworkid == ns ) WRITE(temp_str,'(A,I3.3)') '_s',ik
      END DO
      KN_EXT='kn_log.'//TRIM(proc_string)//TRIM(ADJUSTL(temp_str))
!      IF (myworkid == ns) THEN
!         KN_EXT=
!      ELSE
!         KN_EXT='/dev/null'
!      END IF
      !Read input files (simulation parameters, models, flux-surfaces, species...)
      IF(ANY(sigma_knosos_1nu < bigno)) KN_STELLOPT(1)=.TRUE.
      IF(ANY(sigma_knosos_snu < bigno)) KN_STELLOPT(2)=.TRUE.
      IF(ANY(sigma_knosos_gmc < bigno)) KN_STELLOPT(4)=.TRUE.
      IF(ANY(sigma_knosos_gma < bigno)) KN_STELLOPT(5)=.TRUE.
      IF(ANY(sigma_knosos_dbo < bigno)) KN_STELLOPT(6)=.TRUE.
      IF(ANY(sigma_knosos_vbt < bigno)) KN_STELLOPT(7)=.TRUE.
      IF(ANY(sigma_knosos_vbb < bigno)) KN_STELLOPT(8)=.TRUE.
      IF(ANY(sigma_knosos_wbw < bigno)) KN_STELLOPT(9)=.TRUE.
      CALL READ_INPUT(ns,s,nbb,Zb,Ab,regb,Zeff)
!!$      !Allocate some transport-related quantities
!!$      ALLOCATE(nb(nbb,ns,nerr),dnbdpsi(nbb,ns,nerr),Tb(nbb,ns,nerr),dTbdpsi(nbb,ns,nerr),&
!!$           & Epsi(ns,nerr),Gb(nbb,ns,nerr),Qb(nbb,ns,nerr),Sb(nbb,ns,nerr),Pb(nbb,ns,nerr),rank(ns,nerr))
!!$      nb=0
!!$      dnbdpsi=0
!!$      Tb=0
!!$      dTbdpsi=0
!!$      Epsi=0
!!$      Gb=0
!!$      Qb=0
!!$      Sb=0
!!$      Pb=0
      !Open files
      out_file ='kn_out.'//TRIM(proc_string)
      w_u3 =222
!      w_u6 =223
      istat=0
      IF (myworkid == master) call safe_open(w_u3,istat,out_file,'replace','formatted')
      IF (istat .ne. 0) THEN
         IF (myworkid == master) PRINT *,istat,out_file
         IF (myworkid == master) WRITE(6,*) 'Error opening KNOSOS output file:',TRIM(out_file),istat
         iflag=-1
         RETURN
      END IF
      ! Loop over magnetic surfaces
      IF (ALLOCATED(eff_ripple)) DEALLOCATE(eff_ripple) ; ALLOCATE(eff_ripple(nsd)) ; eff_ripple=0
      IF (ALLOCATED(s_kn))       DEALLOCATE(s_kn);        ALLOCATE(s_kn(nsd))
      IF(myworkid == master) s_kn=rho ; CALL MPI_BCAST(s_kn,nsd,MPI_REAL8,master,MPI_COMM_MYWORLD,ierr_mpi)

      jk=0
      DO ik=2,nsd
         IF(.not. lbooz(ik)) CYCLE
         IF(sigma_knosos_1nu(ik) >= bigno .and. sigma_knosos_snu(ik) >= bigno .and. sigma_knosos_sbp(ik) >= bigno .and. &
            sigma_knosos_gmc(ik) >= bigno .and. sigma_knosos_gma(ik) >= bigno .and. sigma_knosos_qer(ik) >= bigno .and. sigma_knosos_vbt(ik) >= bigno .and. &
            sigma_knosos_vbb(ik) >= bigno .and. sigma_knosos_wbw(ik) >= bigno .and. sigma_knosos_dbo(ik) >= bigno) CYCLE
         jk=jk+1
         IF(myworkid+1 /= jk ) CYCLE
         KN_STELLOPT=.FALSE.
         IF(sigma_knosos_1nu(ik) < bigno) KN_STELLOPT(1)=.TRUE.
         IF(sigma_knosos_snu(ik) < bigno) KN_STELLOPT(2)=.TRUE.
         IF(sigma_knosos_gmc(ik) < bigno) KN_STELLOPT(4)=.TRUE.
         IF(sigma_knosos_gma(ik) < bigno) KN_STELLOPT(5)=.TRUE.
         IF(sigma_knosos_dbo(ik) < bigno) KN_STELLOPT(9)=.TRUE.
         IF(sigma_knosos_vbt(ik) < bigno) KN_STELLOPT(6)=.TRUE.
         IF(sigma_knosos_vbb(ik) < bigno) KN_STELLOPT(7)=.TRUE.
         IF(sigma_knosos_wbw(ik) < bigno) KN_STELLOPT(8)=.TRUE.
         KN_DKESFILE='boozmn_'//TRIM(proc_string)//'.nc'
         CALL READ_BFIELD(s_kn(ik))
         KNOSOS_VBT(ik)=KN_VBT
         KNOSOS_VBB(ik)=KN_VBB
         KNOSOS_WBW(ik)=KN_WBW
         KNOSOS_DBO(ik)=KN_DBO
         IF(sigma_knosos_1nu(ik) >= bigno .and. sigma_knosos_snu(ik) >= bigno .and. sigma_knosos_sbp(ik) >= bigno .and. &
            sigma_knosos_gmc(ik) >= bigno .and. sigma_knosos_gma(ik) >= bigno .and. sigma_knosos_qer(ik) >= bigno .and. sigma_knosos_vbt(ik) >= bigno .and. &
            sigma_knosos_vbb(ik) >= bigno .and. sigma_knosos_wbw(ik) >= bigno .and. sigma_knosos_dbo(ik) >= bigno) CYCLE
         CALL CALC_DATABASE(s_kn,ik,ns)
         KNOSOS_1NU(ik)=KN_1NU
         KNOSOS_SNU(ik)=KN_SNU
         KNOSOS_SBP(ik)=KN_SBP
         KNOSOS_GMC(ik)=KN_GMC
         KNOSOS_GMA(ik)=KN_GMA
         !         nb=0.0;dnbdpsi=0.0/dpsidr;Tb=0.0;dTbdpsi=0.0/dpsidr:  !JLVG: to be done
!         CALL SOLVE_DKE_QN_AMB(itime,nbb,Zb,Ab,regb,s_kn(ik),nb,dnbdpsi,Tb,dTbdpsi,Epsi,Gb,Qb)
         KNOSOS_QER(ik)=KN_QER
!         KNOSOS_DJRDER(ik))=KN_DJRDER
!         KNOSOS_ERSHEAR(ik)=KN_ERSHEAR
      END DO
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,KNOSOS_1NU,nsd,MPI_REAL8,MPI_SUM,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,KNOSOS_SNU,nsd,MPI_REAL8,MPI_SUM,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,KNOSOS_SBP,nsd,MPI_REAL8,MPI_SUM,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,KNOSOS_GMC,nsd,MPI_REAL8,MPI_SUM,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,KNOSOS_GMA,nsd,MPI_REAL8,MPI_SUM,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,KNOSOS_QER,nsd,MPI_REAL8,MPI_SUM,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,KNOSOS_VBT,nsd,MPI_REAL8,MPI_SUM,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,KNOSOS_VBB,nsd,MPI_REAL8,MPI_SUM,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,KNOSOS_WBW,nsd,MPI_REAL8,MPI_SUM,MPI_COMM_MYWORLD,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,KNOSOS_DBO,nsd,MPI_REAL8,MPI_SUM,MPI_COMM_MYWORLD,ierr_mpi)
!DEC$ ENDIF
      IF (myworkid == master) THEN
         DO ik=2,nsd
            IF(.not. lbooz(ik)) CYCLE
            IF(sigma_knosos_1nu(ik) >= bigno .and. sigma_knosos_snu(ik) >= bigno .and. sigma_knosos_sbp(ik) >= bigno .and. &
               sigma_knosos_gmc(ik) >= bigno .and. sigma_knosos_gma(ik) >= bigno .and. sigma_knosos_qer(ik) >= bigno .and. sigma_knosos_vbt(ik) >= bigno .and. &
               sigma_knosos_vbb(ik) >= bigno .and. sigma_knosos_wbw(ik) >= bigno .and. sigma_knosos_dbo(ik) >= bigno) CYCLE
            WRITE(w_u3,'(1(1x,i8),20(1x,e17.10))') ik,KNOSOS_1NU(ik),&
                 & KNOSOS_SNU(ik),KNOSOS_SBP(ik),KNOSOS_GMC(ik),KNOSOS_GMA(ik),KNOSOS_QER(ik),KNOSOS_VBT(ik),KNOSOS_VBB(ik),KNOSOS_WBW(ik),KNOSOS_DBO(ik)
            eff_ripple(ik)=KNOSOS_1NU(ik)
         END DO
         IF (lscreen) WRITE(6,'(2X,I8,1(2X,E17.10))') ik,KNOSOS_1NU(ik),&
              KNOSOS_SNU(ik),KNOSOS_SBP(ik),KNOSOS_GMC(ik),KNOSOS_GMA(ik),KNOSOS_QER(ik),KNOSOS_VBT(ik),KNOSOS_VBB(ik),KNOSOS_WBW(ik),KNOSOS_DBO(ik)
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
