!-----------------------------------------------------------------------
!     Module:        gist_mod
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/18/2013
!     Description:   This module interfaces to the GIST input namelists.
!-----------------------------------------------------------------------
      MODULE gist_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      
!-----------------------------------------------------------------------
!     Module Variables
!         
!-----------------------------------------------------------------------
      LOGICAL     ::  pest, boozer, global_y, benchmarks, initialize,&
                      gs2, gs2_press, vvbal, vmec_6_90, tracer, &
                      viewer, stellopt, terpsichore, nodrift, profw7
      INTEGER     ::  itrmax, mercie, lresis, modes, lmodel, nstart, &
                      nstop, surf, nz0, nalpha0
      INTEGER, DIMENSION(16) :: ldiag
      REAL(rprec) ::  ci, xmax, abserr, relerr, gmin, gmax, ginit,&
                      dg1, dg2, gscal1, gscal2, jstart, jstop, sresis,&
                      dlta, sk, x0, sloc_fac, pol_turns, theta0, &
                      rho_tr, z_tr, t_tr, alpha0, alpha0_start, &
                      alpha0_end, s0
      CHARACTER(30) :: name
      CHARACTER(256) :: vmec_file, table_tag,out_tag,vmec_dir
      
      

      NAMELIST /coordinates/ pest, boozer

      NAMELIST /in_out/  itrmax ,ci ,mercie, &                
         xmax   ,abserr ,relerr ,gmin   ,gmax   ,ginit,  &                
         dg1    ,dg2    ,gscal1 ,gscal2 ,jstart ,jstop,  &                
         lresis ,sresis ,dlta   ,sk     ,x0,    &                
         modes  ,lmodel ,ldiag  ,nstart   ,nstop, benchmarks,  &
         initialize, name, gs2 , gs2_press, vvbal, vmec_6_90, tracer, &
         vmec_file, table_tag, out_tag, vmec_dir, global_y, viewer, stellopt, &
         terpsichore, nodrift, sloc_fac
   
                   
      NAMELIST /setup/  pol_turns,theta0,surf,rho_tr,z_tr,t_tr,alpha0,nz0,profw7, &
         nalpha0,alpha0_start,alpha0_end,s0
      
!-----------------------------------------------------------------------
!     Subroutines
!         read_gist_namelist:    Return boozer s,theta,phi coordiante
!-----------------------------------------------------------------------
      CONTAINS
      
      SUBROUTINE read_gist_namelist (iunit, istat)
      INTEGER :: iunit, istat, istat1, istat2
      nodrift = .FALSE.
      initialize = .FALSE.
      vvbal   = .FALSE.
      gs2     = .FALSE.
      tracer  = .FALSE.
      boozer  = .FALSE.
      pest    = .FALSE.
      global_y = .FALSE.
      benchmarks = .FALSE.
      gs2_press = .FALSE.
      vmec_6_90 = .FALSE.
      viewer = .FALSE.
      stellopt = .FALSE.
      terpsichore = .FALSE.
      nodrift = .FALSE.
      profw7 = .FALSE.
      itrmax = 1
      mercier = 0
      lresis = 0
      modes = 0
      lmodel = 0
      nstart = 1
      nstop = 2
      surf = 1
      nz0 = 1
      nalpha0 = 1
      ldiag = 0
      ci = 0.0
      xmax = 0.0
      sloc_fac = 1.0
      abserr = 1.0E-9
      relerr = 1.0E-3
      gmin = 0
      gmax = 1
      name = ''
      vmec_file = ''
      table_tag = ''
      out_tag = ''
      vmec_dir = ''
      
      !READ (iunit, nml=coordinates, iostat=istat)
      !IF (istat /=0) RETURN
      boozer = .FALSE.
      pest   = .TRUE.
      istat1 = 0
      istat2 = 0
      READ (iunit, nml=in_out, iostat=istat1)
      IF (istat1 /=0) PRINT *,'ERROR READING GIST &IN_OUT NAMELIST'
      READ (iunit, nml=setup, iostat=istat2)
      IF (istat2 /=0) PRINT *,'ERROR READING GIST &SETUP NAMELIST'
      IF (istat1 /=0) istat = istat1
      IF (istat2 /=0) istat = istat2
      RETURN
      
      END SUBROUTINE read_gist_namelist
      
      SUBROUTINE write_gist_namelist (iunit, istat)
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
      WRITE(iunit,'(A)') '&COORDINATES'
      WRITE(iunit,outboo) 'BOOZER',boozer
      WRITE(iunit,outboo) 'PEST',pest
      WRITE(iunit,'(A)') '/'
      WRITE(iunit,'(A)') '&IN_OUT'
      WRITE(iunit,outboo) 'GLOBAL_Y',global_y
      !WRITE (iunit,'(2x,3a)') "OUT_TAG = '", TRIM(out_tag),"'" 
      WRITE(iunit,'(A)') '/'
      WRITE(iunit,'(A)') '&SETUP'
      WRITE(iunit,outflt) 'S0',s0
      WRITE(iunit,outflt) 'ALPHA0',alpha0
      WRITE(iunit,outflt) 'ALPHA0_START',alpha0_start
      WRITE(iunit,outflt) 'ALPHA0_END',alpha0_end
      WRITE(iunit,outint) 'NALPHA0',nalpha0
      WRITE(iunit,outflt) 'POL_TURNS',pol_turns
      WRITE(iunit,outint) 'NZ0',nz0
      WRITE(iunit,'(A)') '/'
      RETURN
      END SUBROUTINE write_gist_namelist
      
      END MODULE gist_mod
