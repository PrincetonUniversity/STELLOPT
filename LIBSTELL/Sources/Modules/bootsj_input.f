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

      nrho = 0; mbuse = 0; nbuse = 0; zeff1 = 1
      dens0 = 0; teti = 0; tempres = 0;
      damp = 0; damp_bs = 0; isymm0 = 0; ate = 0; ati = 0
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

      END MODULE bootsj_input
