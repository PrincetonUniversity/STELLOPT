!-----------------------------------------------------------------------
!     Module:        diagno_input_mod
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/02/2012
!     Description:   This module contains the DIAGNO input namelist and
!                    subroutine which initializes and reads the
!                    DIAGNO input namelist.
!-----------------------------------------------------------------------
      MODULE diagno_input_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE diagno_runtime
      USE safe_open_mod, ONLY: safe_open
      
!-----------------------------------------------------------------------
!     Module Variables
!         
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL :: lexist
!-----------------------------------------------------------------------
!     Input Namelists
!         &DIAGNO_IN
!           nu                    Number of poloidal gridpoints
!           nv                    Number of toroidal gridpoints
!           flux_diag_file        Fluxloop diagnostic specification
!           bprobes file          B-field probe specification
!           mirnov_file           Mirnov Array specification
!           seg_rog_file          Rogowski coils specification
!           bfield_points_file    B-field points specfication
!           flux_turns            Flux loop integer scale factors
!           units                 Units (Assume all quantities in [m])
!           int_type              Diagnostic integration method ('midpoint','simpson','bode')
!           int_step              Diagnostic integration substeps
!           lrphiz                R,PHI(rad),Z instead of X,Y,Z diagnostic specifciation
!           vc_adapt_tol          Adaptive integration absoulte tollerance
!           vc_adapt_rel          Adaptive integration relative tollerance
!           flux_mut_file         Mutual induction file (not fully implemented)
!           lvc_field             Use virtual casing instead of volume integral for free boundary
!-----------------------------------------------------------------------
      namelist /diagno_in/ nu, nv, &
           flux_diag_file, bprobes_file, mirnov_file, seg_rog_file,  &
           bfield_points_file, flux_turns, units, int_type, &
           int_step, lrphiz, vc_adapt_tol, vc_adapt_rel,&
           flux_mut_file, lvc_field, bprobe_turns, luse_extcur, &
           bprobes_mut_file, mir_mut_file, rog_mut_file
      
!-----------------------------------------------------------------------
!     Subroutines
!         read_diagno_input:   Reads diagno_in namelist
!-----------------------------------------------------------------------
      CONTAINS
      
      SUBROUTINE read_diagno_input(filename, istat)
      CHARACTER(*), INTENT(in) :: filename
      INTEGER, INTENT(out) :: istat
      INTEGER :: i, iunit
      ! Initializations
      nu     = 100
      nv     = 100
      flux_diag_file     = ''
      bprobes_file       = ''
      mirnov_file        = ''
      seg_rog_file       = ''
      bfield_points_file = ''
      bprobes_mut_file   = ''
      mir_mut_file      = ''
      rog_mut_file      = ''
      flux_mut_file      = ''
      flux_turns     = 1.0_rprec
      bprobe_turns   = 1.0_rprec
      units          = 1.0_rprec
      int_type       = 'simpson'
      int_step       = 2.0_rprec
      lrphiz         = .FALSE.
      luse_mut       = .FALSE.
      vc_adapt_tol   = 1.0E-06_rprec
      vc_adapt_rel   = 1.0E-04_rprec
      lvc_field      = .TRUE.
      luse_extcur(:) = .TRUE.  ! Do this so we default to using the whole coil if the user forgets
      ! Read namelist
      istat=0; iunit = 25
      INQUIRE(FILE='input.'//TRIM(filename),EXIST=lexist)
      IF (lexist) THEN
         CALL safe_open(iunit,istat,'input.'//TRIM(filename),'old','formatted')
      ELSE
         CALL safe_open(iunit,istat,'diagno.control','old','formatted')
         WRITE(6,*) '   -Reading diango_in namelist from diagno.control file'
      END IF
      IF (lverb .and. istat /=0 ) THEN
         WRITE(6,*) '   -Could not find input.'//TRIM(filename)//' file'
         WRITE(6,*) '    or DIAGNO control file'
         stop
      END IF
      READ(iunit,NML=diagno_in,IOSTAT=istat)
      IF (istat /= 0) THEN
         CLOSE(iunit)
         IF (lverb) THEN
            WRITE(6,*) '   -Could not find DIAGNO_IN namelist in input.'//TRIM(filename)//' file'
            WRITE(6,*) '   -Reading DIAGNO_IN namelist from diagno.control file'
         END IF
         iunit = 25
         CALL safe_open(iunit,istat,'diagno.control','old','formatted')
         IF (istat /= 0) THEN
            WRITE(6,*) '   -Could not open diagno.control file'
            WRITE(6,*) '   -Dumping DIAGNO_IN to screen'
            CALL write_diagno_input(6,istat)
            stop
         END IF
         READ(iunit,NML=diagno_in,IOSTAT=istat)
         IF (istat /= 0) THEN
            WRITE(6,*) '   -Could not read DIAGNO_IN from diagno.control file'
            WRITE(6,*) '   -Dumping DIAGNO_IN to screen'
            CALL write_diagno_input(6,istat)
            stop
         END IF
      END IF
      CALL FLUSH(iunit)
      CALL FLUSH(6)
      CLOSE(iunit)
      ! Clean up string variables
      flux_diag_file = TRIM(flux_diag_file)
      flux_diag_file = ADJUSTL(flux_diag_file)
      bprobes_file = TRIM(bprobes_file)
      bprobes_file = ADJUSTL(bprobes_file)
      mirnov_file = TRIM(mirnov_file)
      mirnov_file = ADJUSTL(mirnov_file)
      seg_rog_file = TRIM(seg_rog_file)
      seg_rog_file = ADJUSTL(seg_rog_file)
      bfield_points_file = TRIM(bfield_points_file)
      bfield_points_file = ADJUSTL(bfield_points_file)
      bprobes_mut_file = TRIM(bprobes_mut_file)
      bprobes_mut_file = ADJUSTL(bprobes_mut_file)
      mir_mut_file = TRIM(mir_mut_file)
      mir_mut_file = ADJUSTL(mir_mut_file)
      rog_mut_file = TRIM(rog_mut_file)
      rog_mut_file = ADJUSTL(rog_mut_file)
      flux_mut_file = TRIM(flux_mut_file)
      flux_mut_file = ADJUSTL(flux_mut_file)
      int_type = TRIM(int_type)
      int_type = ADJUSTL(int_type)

      IF (vc_adapt_tol .EQ. 0.0_rprec .and. lverb) THEN
         vc_adapt_tol = 1.E-6_rprec
         WRITE (6,*)'   !!!!! VC_ADAPT_TOL was reset to 1.0E-6 !!!!!'
      END IF


      ! Handle Mutual Induction
      IF (LEN_TRIM(flux_mut_file)>1) luse_mut = .TRUE.
      IF (LEN_TRIM(bprobes_mut_file)>1) luse_mut = .TRUE.
      IF (LEN_TRIM(mir_mut_file)>1) luse_mut = .TRUE.
      IF (LEN_TRIM(rog_mut_file)>1) luse_mut = .TRUE.
      IF (lmut) luse_mut = .FALSE.
      END SUBROUTINE read_diagno_input
      
      SUBROUTINE write_diagno_input(iunit,istat)
      INTEGER, INTENT(in)    :: iunit
      INTEGER, INTENT(inout) :: istat
      INTEGER ::  n
      CHARACTER(LEN=*), PARAMETER :: outstr  = "(2X,A,1X,'=',1X,'''',A,'''')"
      WRITE(iunit,'(A)') '&DIAGNO_IN'
      WRITE(iunit,"(2X,A,1X,'=',1X,I0)") 'NU',nu
      WRITE(iunit,"(2X,A,1X,'=',1X,I0)") 'NV',nv
      WRITE(iunit,"(2X,A,1X,'=',1X,E22.14)") 'UNITS',units
      WRITE(iunit,"(2X,A,1X,'=',1X,E22.14)") 'VC_ADAPT_TOL',vc_adapt_tol
      WRITE(iunit,"(2X,A,1X,'=',1X,E22.14)") 'VC_ADAPT_REL',vc_adapt_rel
      WRITE(iunit,outstr) 'INT_TYPE',TRIM(int_type)
      WRITE(iunit,"(2X,A,1X,'=',1X,I0)") 'INT_STEP',int_step
      WRITE(iunit,outstr) 'BFIELD_POINTS_FILE',TRIM(bfield_points_file)
      WRITE(iunit,outstr) 'BPROBES_FILE',TRIM(bprobes_file)
      WRITE(iunit,outstr) 'MIRNOV_FILE',TRIM(mirnov_file)
      WRITE(iunit,outstr) 'SEG_ROG_FILE',TRIM(seg_rog_file)
      WRITE(iunit,outstr) 'FLUX_DIAG_FILE',TRIM(flux_diag_file)
      WRITE(iunit,outstr) 'FLUX_MUT_FILE',TRIM(flux_mut_file)
      WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'BPROBE_TURNS',(bprobe_turns(n), n=1,256)
      WRITE(iunit,"(2X,A,1X,'=',10(1X,E22.14))") 'FLUX_TURNS',(flux_turns(n), n=1,256)
      WRITE(iunit,"(2X,A,1X,'=',1X,L1)") 'LRPHIZ',lrphiz
      WRITE(iunit,"(2X,A,1X,'=',1X,L1)") 'LVC_FIELD',lvc_field
      WRITE(iunit,'(A)') '/'
      CALL FLUSH(iunit)
      END SUBROUTINE write_diagno_input

      END MODULE diagno_input_mod
