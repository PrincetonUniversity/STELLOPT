      MODULE bootsj_input
      USE stel_kinds
      INTEGER :: nrho, mbuse, nbuse, isymm0
      REAL(rprec) :: tempres, ate(0:11), ati(0:11)
      REAL(rprec) :: zeff1, dens0, teti, damp, damp_bs

      NAMELIST /bootin/ nrho, mbuse, nbuse, zeff1, dens0, teti, tempres,
     1   damp, damp_bs, isymm0, ate, ati

      CONTAINS

      SUBROUTINE read_boot_namelist (iunit, istat)
      INTEGER :: iunit, istat

      !nrho = 30; mbuse = 0; nbuse = 0; zeff1 = 1
      !dens0 = 0.3; teti = 2; tempres = -1;
      !damp = -0.01; damp_bs = -0.01; isymm0 = 0; ate = -1; ati = -1

      nrho = 30                !number of rho values to use
      mbuse = 6                !number of m (poloidal) terms in B field.
      nbuse = 4                !number of nzetah (toroidal) terms in B field.
      zeff1 = 1.0_dp           !effective ion charge
      dens0 = 0.3_dp           !central electron density in 10**20 m-3
      teti = 2.0_dp            !ratio of Te/Ti for electron to ion
                               !temperature profiles
      tempres = -one           !tempe1(i)=pres(ir)**tempres
                               !if(tempres.eq.-1.0_dp) then
                               !tempe1(i)=sqrt(pres(i))
      damp = -0.01_dp          !superceded by damp_bs
      damp_bs = -0.01_dp       !damping factor for resonances
      isymm0 = 0               !if ne.0 then force a symmetric-device calculation,

      ate    = 0
      ati    = 0
      ate(0) = -1
      ati(0) = -1
      READ (iunit, nml=bootin, iostat=istat)

      END SUBROUTINE read_boot_namelist
      
      SUBROUTINE write_bootsj_input(iunit,istat)
      INTEGER, INTENT(in)    :: iunit
      INTEGER, INTENT(inout) :: istat
      INTEGER ::  n
      WRITE(iunit,'(A)') '&BOOTIN'
      WRITE(iunit,"(2X,A,1X,'=',1X,I0)") 'MBUSE',mbuse
      WRITE(iunit,"(2X,A,1X,'=',1X,I0)") 'NBUSE',nbuse
      WRITE(iunit,"(2X,A,1X,'=',1X,I0)") 'ISYMM0',isymm0
      WRITE(iunit,"(2X,A,1X,'=',1X,E22.14)") 'ZEFF1',zeff1
      WRITE(iunit,"(2X,A,1X,'=',1X,E22.14)") 'DENS0',dens0
      WRITE(iunit,"(2X,A,1X,'=',1X,E22.14)") 'TETI',teti
      WRITE(iunit,"(2X,A,1X,'=',1X,E22.14)") 'TEMPRES',tempres
      WRITE(iunit,"(2X,A,1X,'=',1X,E22.14)") 'DAMP',damp
      WRITE(iunit,"(2X,A,1X,'=',1X,E22.14)") 'DAMP_BS',damp_bs
      WRITE(iunit,"(2X,A,1X,'=',4(1X,E22.14))") 
     1             'ATE',(ate(n), n=0,11)
      WRITE(iunit,"(2X,A,1X,'=',4(1X,E22.14))") 
     1             'ATI',(ati(n), n=0,11)
      WRITE(iunit,'(A)') '/'
      END SUBROUTINE write_bootsj_input

      SUBROUTINE BCAST_BOOTSJ_INPUT(local_master,comm,istat)
!DEC$ IF DEFINED (MPI_OPT)
      USE mpi
!DEC$ ENDIF
      IMPLICIT NONE
      INTEGER, INTENT(inout) :: comm
      INTEGER, INTENT(in)    :: local_master
      INTEGER, INTENT(inout) :: istat
      IF (istat .ne. 0) RETURN
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BCAST(nrho,1,MPI_INTEGER,local_master,comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(mbuse,1,MPI_INTEGER,local_master,comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(nbuse,1,MPI_INTEGER,local_master,comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(isymm0,1,MPI_INTEGER,local_master,comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(zeff1,1,MPI_DOUBLE_PRECISION,local_master,
     1               comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(dens0,1,MPI_DOUBLE_PRECISION,local_master,
     1               comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(teti,1,MPI_DOUBLE_PRECISION,local_master,
     1               comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(tempres,1,MPI_DOUBLE_PRECISION,local_master,
     1               comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(damp,1,MPI_DOUBLE_PRECISION,local_master,
     1               comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(damp_bs,1,MPI_DOUBLE_PRECISION,local_master,
     1               comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(ate,SIZE(ate),MPI_DOUBLE_PRECISION,local_master,
     1               comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(ati,size(ati),MPI_DOUBLE_PRECISION,local_master,
     1               comm,istat)
      IF (istat .ne. 0) RETURN
!DEC$ ENDIF
      END SUBROUTINE BCAST_BOOTSJ_INPUT


      END MODULE bootsj_input
