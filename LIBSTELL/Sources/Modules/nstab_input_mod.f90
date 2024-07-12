!-----------------------------------------------------------------------
!     Module:        nstab_input_mod
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          08/15/13
!     Description:   Reads and NSTAB_INPUT namelist.
!-----------------------------------------------------------------------
      MODULE nstab_input_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE vparams, ONLY: rprec
      IMPLICIT NONE
   
!-----------------------------------------------------------------------
!     Module Variables 
!-----------------------------------------------------------------------
      LOGICAL :: stell, rays
      INTEGER :: numit, icont,numext,starth,stoph,mh,nh,ic
      REAL(rprec) :: p0,xpr,ypr, iota0, iota1, mode

      NAMELIST /nstab_input/ numit, icont, numext, starth, stoph, mh, nh,&
                             ic, p0, xpr, ypr, iota0, iota1, mode, stell, &
                             rays
      
!-----------------------------------------------------------------------
!     SUBROUTINES
!-----------------------------------------------------------------------
      CONTAINS

      SUBROUTINE read_nstab_namelist (iunit, istat)
      INTEGER :: iunit, istat
      
      stell = .FALSE.
      rays  = .FALSE.
      icont = 0
      numit = 1000
      p0    = -1
      xpr   = 0
      ypr   = 0
      iota0 = -1
      iota1 = -1
      numext = 0
      starth = 0
      stoph  = 0
      mh     = 0
      nh     = 0
      ic     = 200
      mode   = 0.0
      
      READ (iunit, nml=nstab_input, iostat=istat)

      END SUBROUTINE read_nstab_namelist
      
      SUBROUTINE write_nstab_namelist (iunit, istat)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: iunit
      INTEGER, INTENT(inout) :: istat
      INTEGER :: iftol,i,n,m
      INTEGER, DIMENSION(1) :: ins
      CHARACTER(LEN=*), PARAMETER :: outboo  = "(2X,A,1X,'=',1X,L1)"
      CHARACTER(LEN=*), PARAMETER :: outint  = "(2X,A,1X,'=',1X,I0)"
      CHARACTER(LEN=*), PARAMETER :: outint1 = "(2X,A,1X,'=',1X,I1.1)"
      CHARACTER(LEN=*), PARAMETER :: outint2 = "(2X,A,1X,'=',1X,I2.2)"
      CHARACTER(LEN=*), PARAMETER :: outint3 = "(2X,A,1X,'=',1X,I3.3)"
      CHARACTER(LEN=*), PARAMETER :: outint4 = "(2X,A,1X,'=',1X,I4.4)"
      CHARACTER(LEN=*), PARAMETER :: outint5 = "(2X,A,1X,'=',1X,I5.5)"
      CHARACTER(LEN=*), PARAMETER :: outint6 = "(2X,A,1X,'=',1X,I6.6)"
      CHARACTER(LEN=*), PARAMETER :: outflt="(2X,A,1X,'=',1X,ES22.12E3)"
      CHARACTER(LEN=*), PARAMETER :: outexp="(2X,A,1X,'=',1X,ES22.12E3)"
      IF (istat < 0) RETURN
      WRITE(iunit,'(A)') '!----- Runtime Parameters -----'
      WRITE(iunit,'(A)') '&NSTAB_INPUT'
      WRITE(iunit,outint) 'NUMIT',numit
      WRITE(iunit,outint) 'ICONT',icont
      WRITE(iunit,outint) 'NUMEXT',numext
      WRITE(iunit,outint) 'STARTH',starth
      WRITE(iunit,outint) 'STOPH',stoph
      WRITE(iunit,outint) 'MH',mh
      WRITE(iunit,outint) 'NH',nh
      WRITE(iunit,outflt) 'MODE',mode
      WRITE(iunit,outflt) 'P0',p0
      WRITE(iunit,outflt) 'XPR',xpr
      WRITE(iunit,outflt) 'YPR',ypr
      WRITE(iunit,outflt) 'IOTA0',iota0
      WRITE(iunit,outflt) 'IOTA1',iota1
      WRITE(iunit,outboo) 'STELL',stell
      WRITE(iunit,outboo) 'RAYS',rays
      WRITE(iunit,'(A)') '/'
      RETURN
      END SUBROUTINE write_nstab_namelist

      END MODULE nstab_input_mod

