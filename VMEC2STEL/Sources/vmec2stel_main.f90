!-----------------------------------------------------------------------
!     Program:       VMEC2STEL
!     Authors:       S. Lazerson
!     Date:          04/15/2014
!     Description:   The VMEC2STEL code reads a VMEC input namelist and
!                    returns a STELLOPTV2 OPTIMUM namelist.
!-----------------------------------------------------------------------
      PROGRAM VMEC2STEL
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds
      USE vmec_input, arg1_input => arg1
      USE safe_open_mod
!-----------------------------------------------------------------------
!     Local Variables
!          numargs      Number of input arguments
!          i            Index
!          arg_len      Length of input strings
!          arg1         Input file
!          args         Input arguments
!-----------------------------------------------------------------------
      IMPLICIT NONE
      ! Command line related
      integer                                      :: numargs,i,ier,iunit
      integer, parameter                           :: arg_len =256
      character*(arg_len)                          :: arg1
      character*(arg_len),allocatable,dimension(:) :: args
      ! Runtime Related
      INTEGER, PARAMETER     :: nu = 64
      INTEGER, PARAMETER     :: nv = 64
      LOGICAL                :: lvmec, lexist, lminmax, lrbc, lrhomn,ldeltamn, &
                                lfix_ntor, llmdif, lgade, lswarm, lmap, lbasic, &
                                lqas, lneed_booz, lqps, lhelical, lballoon, lneo, &
                                ldkes, lbootsj, ltxport, lmap_plane, ljdotb0, liota, &
                                lkink, lvaciota, ljcurv, loutput_harm, lmode, lorbit, &
                                lac, lam, lai, lphiedge, lpscale, lcurtor, lkappa, &
                                lwell, lcurvature, lfieldlines
      INTEGER                :: m,n,ns,j
      REAL(rprec)            :: bound_min, bound_max, var, var_min, var_max, &
                                temp, rho_exp,r1t,r2t,z1t, delta, filter_harm, pi2
      REAL(rprec), DIMENSION(-ntord:ntord,-mpol1d:mpol1d) ::  deltamn
      REAL(rprec), DIMENSION(-ntord:ntord,0:mpol1d)       ::  rhobc
      REAL(rprec), DIMENSION(-ntord:ntord,0:mpol1d)       :: rbc_temp,zbs_temp
      REAL(rprec), DIMENSION(nu,nv)                       :: rreal,zreal
      character(arg_len)     :: id_string, var_name
      
      ! Define output
      CHARACTER(LEN=*), PARAMETER :: onevar  = "(2X,A,1X,'= T ',3(2X,A,1X,'=',1X,ES22.12E3))"
      CHARACTER(LEN=*), PARAMETER :: vecvar  = "(2X,A,'(',I3.3,')',1X,'= T ',3(2X,A,'(',I3.3,')',1X,'=',1X,ES22.12E3))"
      CHARACTER(LEN=*), PARAMETER :: vecvar2 = "(21X,2(2X,A,'(',I3.3,')',1X,'=',1X,ES22.12E3))"
      CHARACTER(LEN=*), PARAMETER :: arrvar  = "(2X,A,'(',I4.3,',',I4.3,') = T',3(2X,A,'(',I4.3,',',I4.3,') = ',1X,ES22.12E3))"
      CHARACTER(LEN=*), PARAMETER :: arrvar2 = "(27x,3(2X,A,'(',I4.3,',',I4.3,') = ',1X,ES22.12E3))"
      CHARACTER(LEN=*), PARAMETER :: target2 = "(2(2X,A,'(',I3.3,')',1X,'=',1X,F8.4))"
      CHARACTER(LEN=*), PARAMETER :: target3 = "(3(2X,A,'(',I3.3,')',1X,'=',1X,F8.4))"
      
      REAL(rprec), PARAMETER :: VMEC2STEL_VERSION = 1.31_rprec

      REAL(rprec), EXTERNAL :: piota
!-----------------------------------------------------------------------
!     Begin Program
!-----------------------------------------------------------------------
      ! Set Defaults
      pi2 = 8.0 * ATAN(1.0)
      ier = 0
      lexist = .FALSE.
      lvmec = .FALSE.
      lminmax = .FALSE.
      lrhomn = .TRUE.
      ldeltamn = .FALSE.
      lrbc = .FALSE.
      lmode = .FALSE.
      lac = .FALSE.
      lam = .FALSE.
      lai = .FALSE.
      lphiedge= .FALSE.
      lpscale = .FALSE.
      lcurtor = .FALSE.
      lmap_plane = .FALSE.
      lfix_ntor = .FALSE.
      llmdif    = .TRUE.
      lgade     = .FALSE.
      lswarm    = .FALSE.
      lmap      = .FALSE.
      lbasic    = .FALSE.
      lneed_booz= .FALSE.
      lqas = .FALSE.
      lqps = .FALSE.
      lhelical = .FALSE.
      lballoon = .FALSE.
      lneo = .FALSE.
      ldkes = .FALSE.
      lbootsj = .FALSE.
      ltxport = .FALSE.
      ljdotb0 = .FALSE.
      liota = .FALSE.
      ljcurv = .FALSE.
      lkink = .FALSE.
      lvaciota = .FALSE.
      lorbit = .FALSE.
      lfieldlines = .FALSE.
      lkappa = .FALSE.
      lwell  = .FALSE.
      lcurvature = .FALSE.
      loutput_harm = .FALSE.
      bound_min = -1.0
      bound_max = 2.0
      filter_harm = 0
      iunit = 327
      rbc = 0; zbs = 0; rbc_temp=0; zbs_temp=0
      ! Hanlde in arguments
      numargs = 0
      i       = 1
      arg1    = ''
      CALL GETCARG(1,arg1,numargs)
      ALLOCATE(args(numargs))
      DO WHILE (i <= numargs)
         CALL GETCARG(i,args(i),numargs)
         SELECT CASE(args(i))
            CASE ("-bounds")
               lminmax = .TRUE.
               i = i + 1
               CALL GETCARG(i,args(i),numargs)
               READ(args(i),*) bound_min
               i = i + 1
               CALL GETCARG(i,args(i),numargs)
               READ(args(i),*) bound_max
            CASE ("-vmec")
               i = i + 1
               lvmec = .TRUE.
               CALL GETCARG(i,id_string,numargs)
            CASE ("-filter")
               i = i + 1
               CALL GETCARG(i,args(i),numargs)
               READ(args(i),*) filter_harm
            CASE ("-am")
               lam = .TRUE.
            CASE ("-ac")
               lac = .TRUE.
            CASE ("-ai")
               lai = .TRUE.
            CASE ("-pscale")
               lpscale = .TRUE.
            CASE ("-curtor")
               lcurtor = .TRUE.
            CASE ("-phiedge")
               lphiedge = .TRUE.
            CASE ("-rbc")
               lrbc     = .TRUE.
               lrhomn   = .FALSE.
               ldeltamn = .FALSE.
               lmode    = .FALSE.
            CASE ("-rhomn")
               lrbc     = .FALSE.
               lrhomn   = .TRUE.
               ldeltamn = .FALSE.
               lmode    = .FALSE.
            CASE ("-deltamn")
               lrbc     = .FALSE.
               lrhomn   = .FALSE.
               ldeltamn = .TRUE.
               lmode    = .FALSE.
            CASE ("-mode")
               lrbc     = .FALSE.
               lrhomn   = .FALSE.
               ldeltamn = .FALSE.
               lmode    = .TRUE.
            CASE ("-lmdif")
               llmdif = .TRUE.
               lgade  = .FALSE.
               lswarm = .FALSE.
               lmap   = .FALSE.
               lmap_plane = .FALSE.
            CASE ("-gade")
               llmdif = .FALSE.
               lgade  = .TRUE.
               lswarm = .FALSE.
               lmap   = .FALSE.
               lmap_plane = .FALSE.
            CASE ("-swarm")
               llmdif = .FALSE.
               lgade  = .FALSE.
               lswarm = .TRUE.
               lmap   = .FALSE.
               lmap_plane = .FALSE.
            CASE ("-map")
               llmdif = .FALSE.
               lgade  = .FALSE.
               lswarm = .FALSE.
               lmap   = .TRUE.
               lmap_plane = .FALSE.
            CASE ("-map_plane")
               llmdif = .FALSE.
               lgade  = .FALSE.
               lswarm = .FALSE.
               lmap   = .FALSE.
               lmap_plane = .TRUE.
            CASE ("-basic")
               lbasic = .TRUE.
            CASE ("-qas")
               lqas = .TRUE.
               lqps = .FALSE.
               lhelical = .FALSE.
               lneed_booz = .TRUE.
            CASE ("-qps")
               lqas = .FALSE.
               lqps = .TRUE.
               lhelical = .FALSE.
               lneed_booz = .TRUE.
            CASE ("-helical")
               lqas = .FALSE.
               lqps = .FALSE.
               lhelical = .TRUE.
               lneed_booz = .TRUE.
            CASE ("-balloon")
               lballoon = .TRUE.
            CASE ("-neo")
               lneo = .TRUE.
               lneed_booz = .TRUE.
            CASE ("-dkes")
               ldkes = .TRUE.
               lneed_booz = .TRUE.
            CASE ("-fix_ntor")
               lfix_ntor = .TRUE.
            CASE ("-bootstrap")
               lneed_booz = .TRUE.
               lbootsj = .TRUE.
            CASE ("-jdotb0")
               ljdotb0 = .TRUE.
            CASE ("-jcurv")
               ljdotb0 = .TRUE.
            CASE ("-vaciota")
               lvaciota = .TRUE.
            CASE ("-txport")
               ltxport = .TRUE.
            CASE ("-iota")
               liota = .TRUE.
            CASE ("-kink")
               lkink = .TRUE.
            CASE ("-orbit")
               lorbit = .TRUE.
            CASE ("-fieldlines")
               lfieldlines = .TRUE.
            CASE ("-harm")
               loutput_harm = .TRUE.
            CASE ("-kappa")
               lkappa= .TRUE.
            CASE ("-magwell")
               lwell= .TRUE.
            CASE ("-curvature")
               lcurvature= .TRUE.
            CASE ("-help","-h")
               WRITE(6,'(a,f5.2)') 'VMEC2STEL Version ',VMEC2STEL_VERSION
               WRITE(6,*) ' STELLOPTV2 Input Generation Utility'
               WRITE(6,*) ' Usage: xvmec2stel <options>'
               WRITE(6,*) '   <options>'
               WRITE(6,*) '   -vmec <ext>   VMEC input extension'
               WRITE(6,*) '   -bounds %min %max Min/Max scaling'
               WRITE(6,*) '   -phiedge          Vary PHIEDGE'
               WRITE(6,*) '   -pscale           Vary PSCALE'
               WRITE(6,*) '   -curtor           Vary CURTOR'
               WRITE(6,*) '   -am               Vary AM Array'
               WRITE(6,*) '   -ai               Vary AC Array'
               WRITE(6,*) '   -ac               Vary AI Array'
               WRITE(6,*) '   -filter %ratio    Boundary Harm. Filter'
               WRITE(6,*) '   -rbc              VMEC Boundary Representation'
               WRITE(6,*) '   -rhomn            H/B Boundary Representation (default)'
               WRITE(6,*) '   -deltamn          Garabedian Boundary Representation'
               WRITE(6,*) '   -mode             Mode pair targeting'
               WRITE(6,*) '   -harm             Output Harmonics (RBC/ZBS) 1% filter'
               WRITE(6,*) '   -fix_ntor         Fix m=0 modes (fixed boundary)'
               WRITE(6,*) '   -lmdif            Levenberg Optimization'
               WRITE(6,*) '   -gade             Differential Evolution'
               WRITE(6,*) '   -swarm            Particle Swarm'
               WRITE(6,*) '   -map              N-Dimensional Mapping'
               WRITE(6,*) '   -map_plane        Hyperplane (2D) Mapping'
               WRITE(6,*) '   -basic            Basic Targets'
               WRITE(6,*) '   -kappa            Plasma Elongation Target'
               WRITE(6,*) '   -qas              QAS Target'
               WRITE(6,*) '   -qps              QPS Target'
               WRITE(6,*) '   -helical          Helical Target'
               WRITE(6,*) '   -balloon          Ballooning Target'
               WRITE(6,*) '   -neo              Neoclassical (NEO) Target'
               WRITE(6,*) '   -dkes             Neoclassical (DKES) Target'
               WRITE(6,*) '   -bootstrap        Bootstrap (BOOTSJ) Target'
               WRITE(6,*) '   -txport           Turbulent Transport Target'
               WRITE(6,*) '   -iota             Iota Profile Target'
               WRITE(6,*) '   -kink             Kink Stability (TERPSICHORE) Target'
               WRITE(6,*) '   -jdotb0           <JdotB> Target'
               WRITE(6,*) '   -jcurv            <Jcurv> Target'
               WRITE(6,*) '   -vaciota          Vacuum Iota (-S12/S11)'
               WRITE(6,*) '   -magwell          Magnetic Well Target'
               WRITE(6,*) '   -curvature        Curvature Kurtosis Target'
               WRITE(6,*) '   -orbit            Orbit (BEAMS3D) Target'
               WRITE(6,*) '   -help             This help message'
               STOP
            CASE DEFAULT
               WRITE(6,'(a,f5.2)') 'VMEC2STEL Version ',VMEC2STEL_VERSION
               WRITE(6,*) ' STELLOPTV2 Input Generation Utility'
               WRITE(6,*) '    CALL with -h or -help for more info.'
               STOP
         END SELECT
         i = i + 1
      END DO
      DEALLOCATE(args)

      ! Generate the OPTIMUM NAMELIST
      WRITE(6,'(A)') '&OPTIMUM'
      IF (llmdif) THEN
         WRITE(6,'(A)')'!-----------------------------------------------------------------------'
         WRITE(6,'(A)')'!          OPTIMIZER RUN CONTROL PARAMETERS'
         WRITE(6,'(A)')'!-----------------------------------------------------------------------'
         WRITE(6,'(2X,A,I5.4)')       'NFUNC_MAX  = ',1000
         IF (lvmec) WRITE(6,'(2x,A)') 'EQUIL_TYPE = ''VMEC2000'''
         WRITE(6,'(2x,A)')            'OPT_TYPE   = ''LMDIF'''
         WRITE(6,'(2X,A,ES10.2)')     'FTOL       = ',1.0E-6
         WRITE(6,'(2X,A,ES10.2)')     'XTOL       = ',1.0E-30
         WRITE(6,'(2X,A,ES10.2)')     'GTOL       = ',1.0E-30
         WRITE(6,'(2X,A,F6.1)')       'FACTOR     = ',10.0
         WRITE(6,'(2X,A,ES10.2)')     'EPSFCN     = ',1.0E-6
         WRITE(6,'(2X,A,I1.1)')       'MODE       = ',1
         WRITE(6,'(2X,A,I1.1)')       'NOPTIMIZERS =',16
         WRITE(6,'(2X,A)')            'LKEEP_MINS = T'
      ELSEIF (lgade) THEN
         WRITE(6,'(A)')'!-----------------------------------------------------------------------'
         WRITE(6,'(A)')'!          OPTIMIZER RUN CONTROL PARAMETERS'
         WRITE(6,'(A)')'!-----------------------------------------------------------------------'
         WRITE(6,'(2X,A,I5.4)')       'NFUNC_MAX   = ',1000
         IF (lvmec) WRITE(6,'(2x,A)') 'EQUIL_TYPE  = ''VMEC2000'''
         WRITE(6,'(2x,A)')            'OPT_TYPE    = ''GADE'''
         WRITE(6,'(2X,A,F6.1)')       'FACTOR      = ',0.5
         WRITE(6,'(2X,A,F6.1)')       'EPSFCN      = ',0.3
         WRITE(6,'(2X,A,I5.4)')       'MODE        = ',2
         WRITE(6,'(2X,A,I5.4)')       'CR_STRATEGY = ',0
         WRITE(6,'(2X,A)')            'LKEEP_MINS  = T'
         WRITE(6,'(2X,A,I5.4)')       'NPOPULATION = ',1000
      ELSEIF (lswarm) THEN
         WRITE(6,'(A)')'!-----------------------------------------------------------------------'
         WRITE(6,'(A)')'!          OPTIMIZER RUN CONTROL PARAMETERS'
         WRITE(6,'(A)')'!-----------------------------------------------------------------------'
         WRITE(6,'(2X,A,I5.4)')       'NFUNC_MAX  = ',1000
         IF (lvmec) WRITE(6,'(2x,A)') 'EQUIL_TYPE = ''VMEC2000'''
         WRITE(6,'(2x,A)')            'OPT_TYPE   = ''PSO'''
         WRITE(6,'(2X,A,ES10.2)')     'FTOL       = ',1.0E-6
         WRITE(6,'(2X,A,ES10.2)')     'XTOL       = ',1.0E-6
         WRITE(6,'(2X,A,F6.1)')       'EPSFCN     = ',0.5
         WRITE(6,'(2X,A,F6.1)')       'FACTOR     = ',100.0
      ELSEIF (lmap) THEN
         WRITE(6,'(A)')'!-----------------------------------------------------------------------'
         WRITE(6,'(A)')'!          OPTIMIZER RUN CONTROL PARAMETERS'
         WRITE(6,'(A)')'!-----------------------------------------------------------------------'
         WRITE(6,'(2X,A,I5.4)')       'NFUNC_MAX  = ',1000
         IF (lvmec) WRITE(6,'(2x,A)') 'EQUIL_TYPE = ''VMEC2000'''
         WRITE(6,'(2x,A)')            'OPT_TYPE   = ''MAP'''
         WRITE(6,'(2X,A,I2.2)')       'MODE       = ',10
         WRITE(6,'(2X,A,I2.2)')       'NPOPULATION = ',16
      ELSEIF (lmap_plane) THEN
         WRITE(6,'(A)')'!-----------------------------------------------------------------------'
         WRITE(6,'(A)')'!          OPTIMIZER RUN CONTROL PARAMETERS'
         WRITE(6,'(A)')'!-----------------------------------------------------------------------'
         WRITE(6,'(2X,A,I5.4)')       'NFUNC_MAX  = ',1000
         IF (lvmec) WRITE(6,'(2x,A)') 'EQUIL_TYPE = ''VMEC2000'''
         WRITE(6,'(2x,A)')            'OPT_TYPE   = ''MAP_PLANE'''
         WRITE(6,'(2X,A,I2.2)')       'MODE       = ',10
         WRITE(6,'(2X,A,I2.2)')       'NPOPULATION = ',16
      ENDIF
      WRITE(6,'(A)')'!-----------------------------------------------------------------------'
      WRITE(6,'(A)')'!          OPTIMIZED QUANTITIES'
      WRITE(6,'(A)')'!-----------------------------------------------------------------------'
      IF (lvmec .and. .not. lmap_plane) THEN
         lfreeb = .false.
         ! Read the input file
         INQUIRE(FILE='input.'//TRIM(id_string),EXIST=lexist)
      	IF (.not.lexist) STOP '!!!!!COULD not file input file!!!!!'
         CALL safe_open(iunit,ier,'input.'//TRIM(id_string),'old','formatted')
         IF(ier /= 0) STOP '!!!!!ERROR opening input file!!!!!'
         CALL read_indata_namelist(iunit,ier)
         IF(ier /= 0) STOP '!!!!!ERROR reading &INDATA namelist!!!!!'
         CLOSE(iunit)
         ! Handle lfreeb=.true. default
         IF (ALL(extcur==0) .and. lfreeb) lfreeb = .FALSE.
         ! Fourier Transform rbc zbs
         rreal = 0; zreal=0;
         DO i = 1, nu
            DO j = 1, nv
               DO m = 0, mpol-1
                  DO n = -ntor,ntor
                     rreal(i,j) = rreal(i,j) + rbc(n,m)*cos(pi2*(m*DBLE(i-1)/(nu-1)-n*DBLE(j-1)/(nv-1)))
                     zreal(i,j) = zreal(i,j) + zbs(n,m)*sin(pi2*(m*DBLE(i-1)/(nu-1)-n*DBLE(j-1)/(nv-1)))
                  END DO
               END DO
            END DO
         END DO
         ! Handle the optimizer
         ! Handle PHIEDGE
         IF (lphiedge) THEN
            var = phiedge
            var_name = 'PHIEDGE'
            IF (var == 0) var = 1
            var_min = bound_min*var
            var_max = bound_max*var
            IF (var_max < var_min) THEN
               temp = var_max
               var_max = var_min
               var_min = temp
            END IF
            WRITE(6,onevar) 'L'//TRIM(var_name)//'_OPT',TRIM(var_name)//'_MIN',var_min,TRIM(var_name)//'_MAX',var_max,'D'//TRIM(var_name)//'_OPT',1.0
         END IF
         ! Handle CURTOR
         IF (lcurtor) THEN
            var = curtor
            var_name = 'CURTOR'
            IF (var == 0) var = 1
            var_min = bound_min*var
            var_max = bound_max*var
            IF (var_max < var_min) THEN
               temp = var_max
               var_max = var_min
               var_min = temp
            END IF
            WRITE(6,onevar) 'L'//TRIM(var_name)//'_OPT',TRIM(var_name)//'_MIN',var_min,TRIM(var_name)//'_MAX',var_max,'D'//TRIM(var_name)//'_OPT',1.0
         END IF
         ! Handle PSCALE
         IF (lpscale) THEN
            var = pres_scale
            var_name = 'PSCALE'
            IF (var == 0) var = 1
            var_min = bound_min*var
            var_max = bound_max*var
            IF (var_max < var_min) THEN
               temp = var_max
               var_max = var_min
               var_min = temp
            END IF
            IF (var_min < 0) var_min = 0  ! Positive Pressure
            WRITE(6,onevar) 'L'//TRIM(var_name)//'_OPT',TRIM(var_name)//'_MIN',var_min,TRIM(var_name)//'_MAX',var_max,'D'//TRIM(var_name)//'_OPT',1.0
         END IF
         IF (lfreeb) THEN
            ! Handle EXTCUR
            var_name = 'EXTCUR'
            DO i = LBOUND(extcur,DIM=1), UBOUND(extcur,DIM=1)
               IF (ABS(extcur(i)) > 0) THEN
                  var = extcur(i)
                  var_min = bound_min*var
                  var_max = bound_max*var
                  IF (var_max < var_min) THEN
                     temp = var_max
                     var_max = var_min
                     var_min = temp
                  END IF
                  WRITE(6,vecvar) 'L'//TRIM(var_name)//'_OPT',i,TRIM(var_name)//'_MIN',i,var_min,TRIM(var_name)//'_MAX',i,var_max,'D'//TRIM(var_name)//'_OPT',i,1.0
               END IF
            END DO
         END IF
         ! Pressure Profile
         IF (lam) THEN
            CALL tolower(pmass_type)
            SELECT CASE(pmass_type)
               CASE('cubic_spline','akima_spline')
                  var_name='AM_F'
                  DO i = LBOUND(am_aux_f,DIM=1),UBOUND(am_aux_f,DIM=1)
                     IF (am_aux_s(i) .ge. 0) THEN
                        var = am_aux_f(i)
                        var_min = bound_min*var
                        var_max = bound_max*var
                        IF (var_max < var_min) THEN
                           temp = var_max
                           var_max = var_min
                           var_min = temp
                        END IF
                        WRITE(6,vecvar) 'L'//TRIM(var_name)//'_OPT',i,TRIM(var_name)//'_MIN',i,var_min,TRIM(var_name)//'_MAX',i,var_max,'D'//TRIM(var_name)//'_OPT',i,1.0
                     END IF
                  END DO
               CASE DEFAULT
                  var_name='AM'
                  DO i = LBOUND(am,DIM=1),UBOUND(am,DIM=1)
                     IF (ABS(am(i)) > 0) THEN
                        var = am(i)
                        var_min = bound_min*var
                        var_max = bound_max*var
                        IF (var_max < var_min) THEN
                           temp = var_max
                           var_max = var_min
                           var_min = temp
                        END IF
                        WRITE(6,vecvar) 'L'//TRIM(var_name)//'_OPT',i,TRIM(var_name)//'_MIN',i,var_min,TRIM(var_name)//'_MAX',i,var_max,'D'//TRIM(var_name)//'_OPT',i,1.0
                     END IF
                  END DO
            END SELECT
         END IF 
         ! Current/Iota Profile
         IF (ncurr == 1 .and. lac) THEN
            CALL tolower(pcurr_type)
            SELECT CASE(pcurr_type)
               CASE('cubic_spline','akima_spline','cubic_spline_ip','akima_spline_ip')
                  var_name='AC_F'
                  DO i = LBOUND(ac_aux_f,DIM=1),UBOUND(ac_aux_f,DIM=1)
                     IF (ac_aux_s(i) .ge. 0) THEN
                        var = ac_aux_f(i)
                        var_min = bound_min*var
                        var_max = bound_max*var
                        IF (var_max < var_min) THEN
                           temp = var_max
                           var_max = var_min
                           var_min = temp
                        END IF
                        WRITE(6,vecvar) 'L'//TRIM(var_name)//'_OPT',i,TRIM(var_name)//'_MIN',i,var_min,TRIM(var_name)//'_MAX',i,var_max,'D'//TRIM(var_name)//'_OPT',i,1.0
                     END IF
                  END DO
               CASE DEFAULT
                  var_name='AC'
                  DO i = LBOUND(ac,DIM=1),UBOUND(ac,DIM=1)
                     IF (ABS(ac(i)) > 0) THEN
                        var = ac(i)
                        var_min = bound_min*var
                        var_max = bound_max*var
                        IF (var_max < var_min) THEN
                           temp = var_max
                           var_max = var_min
                           var_min = temp
                        END IF
                        WRITE(6,vecvar) 'L'//TRIM(var_name)//'_OPT',i,TRIM(var_name)//'_MIN',i,var_min,TRIM(var_name)//'_MAX',i,var_max,'D'//TRIM(var_name)//'_OPT',i,1.0
                     END IF
                  END DO
            END SELECT
         ELSEIF (ncurr==0 .and. lai) THEN
            CALL tolower(piota_type)
            SELECT CASE(piota_type)
               CASE('cubic_spline','akima_spline')
                  var_name='AI_F'
                  DO i = LBOUND(ai_aux_f,DIM=1),UBOUND(ai_aux_f,DIM=1)
                     IF (ai_aux_s(i) .ge. 0) THEN
                        var = ai_aux_f(i)
                        var_min = bound_min*var
                        var_max = bound_max*var
                        IF (var_max < var_min) THEN
                           temp = var_max
                           var_max = var_min
                           var_min = temp
                        END IF
                        WRITE(6,vecvar) 'L'//TRIM(var_name)//'_OPT',i,TRIM(var_name)//'_MIN',i,var_min,TRIM(var_name)//'_MAX',i,var_max,'D'//TRIM(var_name)//'_OPT',i,1.0
                     END IF
                  END DO
               CASE DEFAULT
                  var_name='AI'
                  DO i = LBOUND(ai,DIM=1),UBOUND(ai,DIM=1)
                     IF (ABS(ai(i)) > 0) THEN
                        var = ai(i)
                        var_min = bound_min*var
                        var_max = bound_max*var
                        IF (var_max < var_min) THEN
                           temp = var_max
                           var_max = var_min
                           var_min = temp
                        END IF
                        WRITE(6,vecvar) 'L'//TRIM(var_name)//'_OPT',i,TRIM(var_name)//'_MIN',i,var_min,TRIM(var_name)//'_MAX',i,var_max,'D'//TRIM(var_name)//'_OPT',i,1.0
                     END IF
                  END DO
            END SELECT
         END IF
         ! Magnetic Axis
         IF (.not.lfreeb) THEN
            var_name = 'AXIS'
            DO i = LBOUND(raxis_cc,DIM=1), UBOUND(raxis_cc,DIM=1)
               IF (ABS(raxis_cc(i)) > 0 .or. ABS(zaxis_cs(i)) > 0) THEN
                  var = raxis_cc(i)
                  var_min = bound_min*var
                  var_max = bound_max*var
                  IF (var_max < var_min) THEN
                     temp = var_max
                     var_max = var_min
                     var_min = temp
                  END IF
                  WRITE(6,vecvar) 'L'//TRIM(var_name)//'_OPT',i,'R'//TRIM(var_name)//'_MIN',i,var_min,'R'//TRIM(var_name)//'_MAX',i,var_max,'D'//TRIM(var_name)//'_OPT',i,1.0
                  var = zaxis_cs(i)
                  var_min = bound_min*var
                  var_max = bound_max*var
                  IF (var_max < var_min) THEN
                     temp = var_max
                     var_max = var_min
                     var_min = temp
                  END IF
                  WRITE(6,vecvar2) 'Z'//TRIM(var_name)//'_MIN',i,var_min,'Z'//TRIM(var_name)//'_MAX',i,var_max
               END IF
            END DO
            ! Now do boundary
            IF (lrbc) THEN
               var_name = 'BOUND'
               DO n = LBOUND(rbc,DIM=1), UBOUND(rbc,DIM=1)
                  DO m = LBOUND(rbc,DIM=2), UBOUND(rbc,DIM=2)
                     IF (ABS(rbc(n,m)) > 0 .or. ABS(zbs(n,m)) > 0) THEN
                        var = rbc(n,m)
                        var_min = bound_min*var
                        var_max = bound_max*var
                        IF (var_max < var_min) THEN
                           temp = var_max
                           var_max = var_min
                           var_min = temp
                        END IF
                        WRITE(6,arrvar) 'L'//TRIM(var_name)//'_OPT',n,m,'RBC_MIN',n,m,var_min,'RBC_MAX',n,m,var_max,'DBOUND_OPT',n,m,1.0
                        var = zbs(n,m)
                        var_min = bound_min*var
                        var_max = bound_max*var
                        IF (var_max < var_min) THEN
                           temp = var_max
                           var_max = var_min
                           var_min = temp
                        END IF
                        WRITE(6,arrvar2) 'ZBS_MIN',n,m,var_min,'ZBS_MAX',n,m,var_max
                     END IF
                  END DO
               END DO
            ELSEIF (lmode) THEN
               var_name = 'MODE'
               DO n = LBOUND(rbc,DIM=1), UBOUND(rbc,DIM=1)
                  DO m = LBOUND(rbc,DIM=2), UBOUND(rbc,DIM=2)
                     IF (ABS(rbc(n,m)) > 0 .or. ABS(zbs(n,m)) > 0) THEN
                        var = 0.5*(rbc(n,m) + zbs(n,m))
                        var_min = bound_min*var
                        var_max = bound_max*var
                        IF (var_max < var_min) THEN
                           temp = var_max
                           var_max = var_min
                           var_min = temp
                        END IF
                        WRITE(6,arrvar) 'L'//TRIM(var_name)//'_OPT',n,m,'BOUND_MIN',n,m,var_min,'BOUND_MAX',n,m,var_max,'DBOUND_OPT',n,m,1.0
                     END IF
                  END DO
               END DO
            ELSEIF (lrhomn) THEN
               rhobc = 0.0_rprec
               rho_exp = 4
               rbc_temp = rbc
               zbs_temp = zbs
               CALL convert_boundary(rbc_temp,zbs_temp,rhobc,mpol1d,ntord,rho_exp)
               rbc_temp = rbc
               zbs_temp = zbs
               WHERE(ABS(rhobc) < ABS(filter_harm*rhobc(0,1))) rhobc = 0
               IF (loutput_harm) THEN
                  DO n = -ntor,ntor
                     DO m = 0, mpol
                        IF (abs(rhobc(n,m)) > 0.0) &
                            WRITE(6,'(1(A,I2,A,I2,A,E22.12))') '!  RHOMN(',n,',',m,') = ',rhobc(n,m)
                     END DO
                  END DO
               END IF
               CALL unique_boundary(rbc_temp,zbs_temp,rhobc,mpol1d,ntord,mpol-1,ntor,mpol-1,rho_exp)
               IF (loutput_harm) THEN
                  DO n = -ntor,ntor
                     DO m = 0, mpol
                         IF (abs(rbc_temp(n,m)) > 0.0 .or. abs(zbs_temp(n,m)) > 0.0) &
                         WRITE(6,'(2(A,I2,A,I2,A,E22.12))') '  RBC(',n,',',m,') = ',rbc_temp(n,m),'  ZBS(',n,',',m,') = ',zbs_temp(n,m)
                     END DO
                  END DO
               END IF
               delta = 0
               DO m = 0, mpol
                  DO n = -ntor, ntor
                     IF (rbc(n,m) /= 0) delta = delta + (rbc(n,m)-rbc_temp(n,m))**2/rbc(n,m)**2
                     IF (zbs(n,m) /= 0) delta = delta + (zbs(n,m)-zbs_temp(n,m))**2/zbs(n,m)**2
                  END DO
               END DO
               delta = sqrt(delta)
               WRITE(6,'(A)')'!-----------------------------------------------------------------------'
               WRITE(6,'(A)')'!          HIRSHMAN-BRESLAU BOUNDARY REPRESENTATION'
               WRITE(6,'(A,F10.2,A)')'!            BOUNDARY CONVERSION ACCURACY: ',100*(1-delta),'%'
               WRITE(6,'(A)')'!-----------------------------------------------------------------------'
               ! Mimic LFIX_NTOR
               m=0
               var_name = 'RHO'
               DO n = 0,ntor
                  IF ((rbc(n,0) /= 0 .or. zbs(n,0) /= 0) .and. (n /= 0) .and. (.not. lfix_ntor)) THEN
                     var = rbc(n,0)
                     var_min = bound_min*var
                     var_max = bound_max*var
                     IF (var_max < var_min) THEN
                        temp = var_max
                        var_max = var_min
                        var_min = temp
                     END IF
                     WRITE(6,arrvar) 'LBOUND_OPT',n,m,'RBC_MIN',n,m,var_min,'RBC_MAX',n,m,var_max,'DBOUND_OPT',n,m,1.0
                     var = zbs(n,0)
                     var_min = bound_min*var
                     var_max = bound_max*var
                     IF (var_max < var_min) THEN
                        temp = var_max
                        var_max = var_min
                        var_min = temp
                     END IF
                     WRITE(6,arrvar2) 'ZBS_MIN',n,m,var_min,'ZBS_MAX',n,m,var_max
                  END IF
               END DO
               ! Now do rho
               !DO n = LBOUND(rhobc,DIM=1), UBOUND(rhobc,DIM=1)
               DO n = -ntor,ntor
                  DO m = LBOUND(rhobc,DIM=2), UBOUND(rhobc,DIM=2)
                     IF (ABS(rhobc(n,m)) > 0 .and. (m /= 0 .or. n >= 0)) THEN
                        var = rhobc(n,m)
                        var_min = bound_min*var
                        var_max = bound_max*var
                        IF (var_max < var_min) THEN
                           temp = var_max
                           var_max = var_min
                           var_min = temp
                        END IF
                        WRITE(6,arrvar) 'L'//TRIM(var_name)//'_OPT',n,m,'BOUND_MIN',n,m,var_min,'BOUND_MAX',n,m,var_max,'DRHO_OPT',n,m,1.0
                     END IF
                  END DO
               END DO
            ELSEIF (ldeltamn) THEN
               var_name = 'DELTA'
               deltamn = 0.0_rprec
               rbc_temp = rbc
               zbs_temp = zbs
               CALL convert_boundary_PG(rbc_temp,zbs_temp,deltamn,mpol1d,ntord)
               ! Quick and dirty LPKU Test
               !deltamn = 0
               !deltamn(-1,-1) = 0.087
               !deltamn( 0,-1) = 0.119
               !deltamn( 0, 0) = 1.000
               !deltamn( 1, 0) = 0.099! deltamn(0,0) = 5.178*2
               !deltamn( 0, 1) = 5.178
               !deltamn( 1, 1) = 0.107
               !deltamn( 2, 1) =-0.035
               !deltamn(-1, 2) = 0.049
               !deltamn( 0, 2) =-0.313
               !deltamn( 1, 2) =-0.350
               !deltamn( 2, 2) =-0.019
               !deltamn(-1, 3) =-0.020
               !deltamn( 0, 3) = 0.132
               !deltamn( 1, 3) = 0.069
               !deltamn( 2, 3) = 0.061
               !deltamn( 0, 4) =-0.002
               !deltamn( 1, 4) = 0.000
               !deltamn( 3, 4) =-0.017
               rbc_temp = rbc
               zbs_temp = zbs
               WHERE(ABS(deltamn) < ABS(filter_harm*deltamn(0,1))) deltamn = 0
               IF (loutput_harm) THEN
                  DO n = -ntor,ntor
                     DO m = -mpol, mpol
                        IF (abs(deltamn(n,m)) > 0.0) &
                            WRITE(6,'(1(A,I2,A,I2,A,E22.12))') '!  DELTAMN(',n,',',m,') = ',deltamn(n,m)
                     END DO
                  END DO
               END IF
               CALL unique_boundary_PG(rbc_temp,zbs_temp,deltamn,ntord,mpol1d,mpol,ntor)
               IF (loutput_harm) THEN
                  DO n = -ntor,ntor
                     DO m = 0, mpol
                         IF (abs(rbc_temp(n,m)) > 0.0 .or. abs(zbs_temp(n,m)) > 0.0) &
                         WRITE(6,'(2(A,I2,A,I2,A,E22.12))') '  RBC(',n,',',m,') = ',rbc_temp(n,m),'  ZBS(',n,',',m,') = ',zbs_temp(n,m)
                     END DO
                  END DO
               END IF
               delta = 0
               DO m = 0, mpol
                  DO n = -ntor, ntor
                     IF (rbc(n,m) /= 0) delta = delta + (rbc(n,m)-rbc_temp(n,m))**2/rbc(n,m)**2
                     IF (zbs(n,m) /= 0) delta = delta + (zbs(n,m)-zbs_temp(n,m))**2/zbs(n,m)**2
                  END DO
               END DO
               delta = sqrt(delta)
               WRITE(6,'(A)')'!-----------------------------------------------------------------------'
               WRITE(6,'(A)')'!          GARABEDIAN BOUNDARY REPRESENTATION'
               WRITE(6,'(A,F10.2,A)')'!            BOUNDARY CONVERSION ACCURACY: ',100*(1-delta),'%'
               WRITE(6,'(A)')'!-----------------------------------------------------------------------'
               !DO n = LBOUND(deltamn,DIM=1), UBOUND(deltamn,DIM=1)
               DO n = -ntor,ntor
                  DO m = LBOUND(deltamn,DIM=2), UBOUND(deltamn,DIM=2)
                     IF (ABS(deltamn(n,m)) > 0 &
                         .and. .not.(n == 0 .and. m == 0)) THEN
                        var = deltamn(n,m)
                        var_min = bound_min*var
                        var_max = bound_max*var
                        IF (var_max < var_min) THEN
                           temp = var_max
                           var_max = var_min
                           var_min = temp
                        END IF
                        WRITE(6,arrvar) 'L'//TRIM(var_name)//'MN_OPT',n,m,'DELTA_MIN',n,m,var_min,'DELTA_MAX',n,m,var_max,'DDELTAMN_OPT',n,m,1.0
                     END IF
                  END DO
               END DO
               ! Add in RBC(n,m=0) modes
               m=0
               DO n = 0, ntor
                  IF (rbc(n,0) .ne. 0 .or. zbs(n,0) .ne. 0) THEN
                  var = rbc(n,0)
                     var_min = bound_min*var
                     var_max = bound_max*var
                     IF (var_max < var_min) THEN
                        temp = var_max
                        var_max = var_min
                        var_min = temp
                     END IF
                     WRITE(6,arrvar) 'LBOUND_OPT',n,m,'RBC_MIN',n,m,var_min,'RBC_MAX',n,m,var_max,'DBOUND_OPT',n,m,1.0
                     var = zbs(n,0)
                     var_min = bound_min*var
                     var_max = bound_max*var
                     IF (var_max < var_min) THEN
                        temp = var_max
                        var_max = var_min
                        var_min = temp
                     END IF
                     WRITE(6,arrvar2) 'ZBS_MIN',n,m,var_min,'ZBS_MAX',n,m,var_max
                  END IF
               END DO
            END IF
         END IF
      ELSE IF (lvmec .and. lmap_plane) THEN
         ! Read the input file
         INQUIRE(FILE='input.'//TRIM(id_string),EXIST=lexist)
      	IF (.not.lexist) STOP '!!!!!COULD not file input file!!!!!'
         CALL safe_open(iunit,ier,'input.'//TRIM(id_string),'old','formatted')
         IF(ier /= 0) STOP '!!!!!ERROR opening input file!!!!!'
         CALL read_indata_namelist(iunit,ier)
         IF(ier /= 0) STOP '!!!!!ERROR reading &INDATA namelist!!!!!'
         CLOSE(iunit)
         WRITE(6,'(2X,A)')'LPHIEDGE_OPT = T'
         WRITE(6,'(2X,A)')'LPSCALE_OPT = T'
         WRITE(6,'(2X,A)')'LCURTOR_OPT = T'
         CALL tolower(pmass_type)
         SELECT CASE(pmass_type)
            CASE('cubic_spline','akima_spline')
               m = LBOUND(am_aux_f,DIM=1)
               n = UBOUND(am_aux_f,DIM=1)
               WRITE(6,'(2X,A,I5,A)')'LAM_F_OPT = ',n-m+1,'*T'
            CASE DEFAULT
               m = LBOUND(am,DIM=1)
               n = UBOUND(am,DIM=1)
               WRITE(6,'(2X,A,I5,A)')'LAM_OPT = ',n-m+1,'*T'
         END SELECT
         IF (ncurr == 1) THEN
            CALL tolower(pcurr_type)
            SELECT CASE(pcurr_type)
               CASE('cubic_spline','akima_spline','cubic_spline_ip','akima_spline_ip')
                  m = LBOUND(ac_aux_f,DIM=1)
                  n = UBOUND(ac_aux_f,DIM=1)
                  WRITE(6,'(2X,A,I5,A)')'LAC_F_OPT = ',n-m+1,'*T'
               CASE DEFAULT
                  m = LBOUND(ac,DIM=1)
                  n = UBOUND(ac,DIM=1)
                  WRITE(6,'(2X,A,I5,A)')'LAC_OPT = ',n-m+1,'*T'
            END SELECT
         ELSE
            CALL tolower(piota_type)
            SELECT CASE(piota_type)
               CASE('cubic_spline','akima_spline')
                  m = LBOUND(ai_aux_f,DIM=1)
                  n = UBOUND(ai_aux_f,DIM=1)
                  WRITE(6,'(2X,A,I5,A)')'LAI_F_OPT = ',n-m+1,'*T'
               CASE DEFAULT
                  m = LBOUND(ai,DIM=1)
                  n = UBOUND(ai,DIM=1)
                  WRITE(6,'(2X,A,I5,A)')'LAI_OPT = ',n-m+1,'*T'
            END SELECT
         END IF
         IF (lfreeb) THEN
            m = LBOUND(extcur,DIM=1)
            n = UBOUND(extcur,DIM=1)
            WRITE(6,'(2X,A,I5,A)')'LEXTCUR_OPT = ',n-m+1,'*T'
         END IF
         WRITE(6,'(2X,A,I5,A)')'LAXIS_OPT = ',ntor+1,'*T'
         IF (lrbc) THEN
            DO m = 0, mpol-1
               WRITE(6,"(2X,A,I4.3,':',I3.3,',',I4,A,I4,A)") 'LBOUND_OPT(',-ntor,ntor,m,') = ',2*ntor+1,'*T'
            END DO
         ELSEIF (lrhomn) THEN
            DO m = 0, mpol-2
               WRITE(6,"(2X,A,I4.3,':',I3.3,',',I4,A,I4,A)") 'LRHO_OPT(',-ntor,ntor,m,') = ',2*ntor+1,'*T'
            END DO
         ELSEIF (ldeltamn) THEN
            DO m = -mpol, mpol
               WRITE(6,"(2X,A,I4.3,':',I3.3,',',I4,A,I4,A)") 'LDELTAMN_OPT(',-ntor,ntor,m,') = ',2*ntor+1,'*T'
            END DO
         END IF
      END IF
      IF (lbasic) THEN
         WRITE(6,'(A)')'!-----------------------------------------------------------------------'
         WRITE(6,'(A)')'!          PLASMA PROPERTIES'
         WRITE(6,'(A)')'!-----------------------------------------------------------------------'
         WRITE(6,'(2X,A,F6.3,2X,A,ES10.1)') 'TARGET_ASPECT  = ',rbc(0,0)/rbc(0,1),'SIGMA_ASPECT = ',0.001
         WRITE(6,'(2X,A,F6.3,2X,A,ES10.1)') 'TARGET_BETA    = ',0.05,'SIGMA_BETA = ',0.001
         WRITE(6,'(2X,A,ES10.1,2X,A,ES10.1)') 'TARGET_CURTOR  = ',curtor,'SIGMA_CURTOR = ',0.01*curtor
         WRITE(6,'(2X,A,F6.3,2X,A,ES10.1)') 'TARGET_PHIEDGE = ',phiedge,'SIGMA_PHIEDGE = ',0.01*phiedge
         WRITE(6,'(2X,A,F6.3,2X,A,ES10.1)') 'TARGET_RBTOR   = ',1.0,'SIGMA_RBTOR = ',0.01
         WRITE(6,'(2X,A,F6.3,2X,A,ES10.1)') 'TARGET_R0      = ',raxis_cc(0),'SIGMA_R0 = ',0.01*raxis_cc(0)
         IF (lasym) &
         WRITE(6,'(2X,A,F6.3,2X,A,ES10.1)') 'TARGET_Z0      = ',zaxis_cs(0),'SIGMA_Z0 = ',0.01*raxis_cs(0)
         WRITE(6,'(2X,A,F6.3,2X,A,ES10.1)') 'TARGET_VOLUME  = ',1.0,'SIGMA_VOLUME = ',0.01
         WRITE(6,'(2X,A,F6.3,2X,A,ES10.1)') 'TARGET_WP      = ',1.0,'SIGMA_WP = ',0.01
      END IF
      IF (lkappa) THEN
         WRITE(6,'(A)')'!-----------------------------------------------------------------------'
         WRITE(6,'(A)')'!          PLASMA ELONGATION (kappa)'
         WRITE(6,'(A)')'!             PHI is toroidal angle [0,2pi] over whole device'
         WRITE(6,'(A)')'!             Values output are approximate'
         WRITE(6,'(A)')'!-----------------------------------------------------------------------'
         r1t = SUM(rbc(0,0:mpol-1))
         r2t = SUM(rbc(0,0:mpol-1:2)) - SUM(rbc(0,1:mpol-1:2))
         z1t = SUM(zbs(0,1:mpol-1:2))
         WRITE(6,'(2X,A,F6.3,2X,A,ES10.1,2X,A,F6.3)') 'TARGET_KAPPA     = ',2*z1t/(r1t-r2t),'SIGMA_KAPPA     = ',0.01*z1t/(r1t-r2t),' PHI_KAPPA = ',0.0
         WRITE(6,'(2X,A,F6.3,2X,A,ES10.1,2X,A,F6.3)') 'TARGET_KAPPA_BOX = ',2*z1t/(r1t-r2t),'SIGMA_KAPPA_BOX = ',0.01*z1t/(r1t-r2t),' PHI_KAPPA_BOX = ',0.0
         WRITE(6,'(2X,A,F6.3,2X,A,ES10.1)')           'TARGET_KAPPA_AVG = ',2*z1t/(r1t-r2t),'SIGMA_KAPPA_AVG = ',0.01*z1t/(r1t-r2t)
      END IF
      IF (lneed_booz) THEN
         WRITE(6,'(A)')'!-----------------------------------------------------------------------'
         WRITE(6,'(A)')'!          Boozer Coordinate Transformation'
         WRITE(6,'(A)')'!-----------------------------------------------------------------------'
         WRITE(6,'(2X,A,I3)') 'MBOZ = ',2 ** CEILING(log(REAL(MAX(6*mpol,2)))/log(2.0_rprec))
         WRITE(6,'(2X,A,I3)') 'NBOZ = ',2 ** CEILING(log(REAL(MAX(2*ntor-1,1)))/log(2.0_rprec))
      END IF
      ns = MAXVAL(ns_array)
      IF (liota) THEN
         WRITE(6,'(A)')'!------------------------------------------------------------------------'
         WRITE(6,'(A)')'!       IOTA PROFILE TARGETS'
         WRITE(6,'(A)')'!------------------------------------------------------------------------'
         DO i = 1, ns
            temp = REAL(i-1) / REAL(ns-1)
            WRITE(6,target3) 'TARGET_IOTA',i,piota(temp),'SIGMA_IOTA',i,0.01,'S_IOTA',i,temp
         END DO
      END IF
      IF (lvaciota) THEN
         WRITE(6,'(A)')'!------------------------------------------------------------------------'
         WRITE(6,'(A)')'!       VACUUM IOTA PROFILE TARGETS (-S12/S11)'
         WRITE(6,'(A)')'!------------------------------------------------------------------------'
         DO i = 1, ns
            temp = REAL(i-1) / REAL(ns-1)
            WRITE(6,target3) 'TARGET_VACIOTA',i,piota(temp),'SIGMA_VACIOTA',i,0.01,'S_VACIOTA',i,temp
         END DO
      END IF
      IF (lqas .or. lqps .or. lhelical) THEN
         WRITE(6,'(A)')'!-----------------------------------------------------------------------'
         WRITE(6,'(A)')'!       Boozer Coordinate Helicity'
         WRITE(6,'(A)')'!         Note that helicity targeting is by surface.  Axis (01) is ignored.'
         WRITE(6,'(A)')'!         (X,0): Quasi-Axisymetry'
         WRITE(6,'(A)')'!         (0,X): Quasi-Poloidal Symmetry'
         WRITE(6,'(A)')'!         (L,K): Quasi-Helical Symmetry (m *K + n*L)'
         WRITE(6,'(A)')'!-----------------------------------------------------------------------'
         IF (lqas) WRITE(6,'(2X,A)') 'HELICITY = (1,0)'
         IF (lqps) WRITE(6,'(2X,A)') 'HELICITY = (0,1)'
         IF (lhelical) WRITE(6,'(2X,A)') 'HELICITY = (2,1)'
         WRITE(6,'(2X,A,I3.3,A,I3.3,A,I3.3,A,I3.3,A)') 'TARGET_HELICITY(1:',ns,') = ',ns,'*0.0  SIGMA_HELICITY(1:',ns,') = ',ns,'*-1.0'
      END IF
      IF (lballoon) THEN
         WRITE(6,'(A)')'!------------------------------------------------------------------------'
         WRITE(6,'(A)')'!       Ballooning Stability (as calculated by COBRA_VMEC)'
         WRITE(6,'(A)')'!         Note that ballooning stability is by surface.  Axis (01) is ignored.'
         WRITE(6,'(A)')'!         THETA, ZETA: Ballooning angle perturbations'
         WRITE(6,'(A)')'!------------------------------------------------------------------------'
         WRITE(6,'(2X,A)') 'BALLOON_THETA = 0.0'
         WRITE(6,'(2X,A)') 'BALLOON_ZETA  = 0.0'
         WRITE(6,'(2X,A,I3.3,A,I3.3,A,I3.3,A,I3.3,A)') 'TARGET_BALLOON(1:',ns,') = ',ns,'*0.0  SIGMA_BALLOON(1:',ns,') = ',ns,'*1.0'
      END IF
      IF (lneo) THEN
         WRITE(6,'(A)')'!------------------------------------------------------------------------'
         WRITE(6,'(A)')'!       Neoclassical Transport Calculation (as calculated by NEO)'
         WRITE(6,'(A)')'!------------------------------------------------------------------------'
         WRITE(6,'(2X,A,I3.3,A,I3.3,A,I3.3,A,I3.3,A)') 'TARGET_NEO(1:',ns,') = ',ns,'*0.0  SIGMA_NEO(1:',ns,') = ',ns,'*1.0'
      END IF
      IF (ldkes) THEN
         WRITE(6,'(A)')'!------------------------------------------------------------------------'
         WRITE(6,'(A)')'!       Neoclassical Transport Calculation (as calculated by DKES)'
         WRITE(6,'(A)')'!------------------------------------------------------------------------'
         WRITE(6,'(2X,A,I3.3,A,I3.3,A,I3.3,A,I3.3,A)') 'TARGET_DKES(1:',ns,') = ',ns,'*0.0  SIGMA_DKES(1:',ns,') = ',ns,'*1.0'
      END IF
      IF (lbootsj) THEN
         ! Note boozer quantities are on half grid so nboot = ns-1
         WRITE(6,'(A)')'!------------------------------------------------------------------------'
         WRITE(6,'(A)')'!       Bootstrap Current Calculation (as calculated by BOOTSJ)'
         WRITE(6,'(A)')'!       NE/TE/TI profiles required!'
         WRITE(6,'(A)')'!------------------------------------------------------------------------'
         WRITE(6,'(2X,A,I3.3,A,I3.3,A,I3.3,A,I3.3,A)') 'TARGET_BOOTSTRAP(2:',ns,') = ',ns-1,'*0.0  SIGMA_BOOTSTRAP(1:',ns-1,') = ',ns-1,'*1.0'
      END IF
      IF (lkink) THEN
         WRITE(6,'(A)')'!------------------------------------------------------------------------'
         WRITE(6,'(A)')'!       Kink Stability Calculation (as calculated by TERPSICHORE)'
         WRITE(6,'(A)')'!           terpsichore_input file required in run directory!'
         WRITE(6,'(A)')'!------------------------------------------------------------------------'
         WRITE(6,'(2X,A)') 'MLMNB_KINK = 264 IVAC_KINK = 24'
         WRITE(6,'(2X,A)') 'TARGET_KINK(01) = 1.0E-3  SIGMA_KINK(01) = 1.0   ! terpsichore_input_00'
         WRITE(6,'(4X,A)') 'NJ_KINK(01) = 256 NK_KINK(01) = 256'
         WRITE(6,'(4X,A)') 'MLMNS_KINK(01) = 76 LSSD_KINK(01) = 4096 LSSL_KINK(01) = 4096'
         WRITE(6,'(2X,A)') 'TARGET_KINK(02) = 1.0E-3  SIGMA_KINK(02) = 1.0   ! terpsichore_input_01'
         WRITE(6,'(4X,A)') 'NJ_KINK(02) = 256 NK_KINK(02) = 256'
         WRITE(6,'(4X,A)') 'MLMNS_KINK(02) = 76 LSSD_KINK(02) = 4096 LSSL_KINK(02) = 4096'
      END IF
      IF (lorbit) THEN
         WRITE(6,'(A)')'!------------------------------------------------------------------------'
         WRITE(6,'(A)')'!       ORBIT (BEAMS3D) OPTIMIZATION'
         WRITE(6,'(A)')'!          ASSUMES ALPHAS'
         WRITE(6,'(A)')'!------------------------------------------------------------------------'
         WRITE(6,'(2X,A)') 'MASS_ORBIT = 6.64465675E-27'
         WRITE(6,'(2X,A)') 'Z_ORBIT = 2'
         WRITE(6,'(2X,A)') 'NU_ORBIT = 16'
         WRITE(6,'(2X,A)') 'NV_ORBIT = 16'
         WRITE(6,target2) 'TARGET_ORBIT',2,0.0,'SIGMA_ORBIT',2,1.0
         WRITE(6,'(2X,A)') 'NP_ORBIT = 8'
         DO i = 1, 8
            WRITE(6,target2) 'VLL_ORBIT',i,1.0,'VPERP_ORBIT',i,1.0
         END DO
      END IF
      IF (ltxport) THEN
         WRITE(6,'(A)')'!------------------------------------------------------------------------'
         WRITE(6,'(A)')'!       TURBULENT TRANSPORT'
         WRITE(6,'(A)')'!------------------------------------------------------------------------'
         WRITE(6,'(2X,A)') 'TXPORT_PROXY = ''prox1d'''
         WRITE(6,'(2X,A)') 'LGLOBAL_TXPORT = F'
         WRITE(6,'(2X,A)') 'NZ_TXPORT = 128'
         WRITE(6,'(2X,A)') 'NALPHA_TXPORT = 1'
         WRITE(6,'(2X,A)') 'ALPHA_START_TXPORT = 0.0'
         WRITE(6,'(2X,A,ES10.8)') 'ALPHA_END_TXPORT = ', pi2/(2*nfp)
         DO i = 2, ns
            WRITE(6,target3) 'TARGET_TXPORT',i,0.0,'SIGMA_TXPORT',i,1.0,'S_TXPORT',i,REAL(i)/REAL(ns)
         END DO
      END IF
      IF (ljdotb0) THEN
         WRITE(6,'(A)')'!------------------------------------------------------------------------'
         WRITE(6,'(A)')'!       Flux surface averaged parallel current density (<JDOTB>)'
         WRITE(6,'(A)')'!------------------------------------------------------------------------'
         WRITE(6,'(2X,A,I3.3,A,I3.3,A,I3.3,A,I3.3,A)') 'TARGET_JDOTB(1:',ns,') = ',ns,'*0.0  SIGMA_JDOTB(1:',ns,') = ',ns,'*1.0'
      END IF
      IF (ljcurv) THEN
         WRITE(6,'(A)')'!------------------------------------------------------------------------'
         WRITE(6,'(A)')'!       Flux surface averaged toroidal current density (<JCURV>)'
         WRITE(6,'(A)')'!------------------------------------------------------------------------'
         WRITE(6,'(2X,A,I3.3,A,I3.3,A,I3.3,A,I3.3,A)') 'TARGET_JCURV(1:',ns,') = ',ns,'*0.0  SIGMA_JCURV(1:',ns,') = ',ns,'*1.0'
      END IF
      IF (lwell) THEN
         WRITE(6,'(A)')'!------------------------------------------------------------------------'
         WRITE(6,'(A)')'!       Magnetic Well'
         WRITE(6,'(A)')'!       Stable:    +'
         WRITE(6,'(A)')'!       Un-stable: -'
         WRITE(6,'(A)')'!------------------------------------------------------------------------'
         WRITE(6,'(2X,A,I3.3,A,I3.3,A,I3.3,A,I3.3,A)') 'TARGET_MAGWELL(1:',ns,') = ',ns,'*1.0  SIGMA_MAGWELL(1:',ns,') = ',ns,'*1.0'
      END IF
      IF (lcurvature) THEN
         WRITE(6,'(A)')'!------------------------------------------------------------------------'
         WRITE(6,'(A)')'!       Curvature Kurtosis at Boundary'
         WRITE(6,'(A)')'!------------------------------------------------------------------------'
         WRITE(6,'(2X,A,F6.3,2X,A,ES10.1)') 'TARGET_CURVATURE = ',0.0,'SIGMA_CURVATURE = ',1.00
      END IF
      ! END OPTIMUM Namelist
      WRITE(6,'(A)') '/'
      ! Add namelists
      IF (lneo) THEN
         WRITE(6,'(A)') '&NEO_IN'
         WRITE(6,'(A)') '/'
      END IF
      IF (lbootsj) THEN
         WRITE(6,'(A)') '&BOOTIN'
         WRITE(6,'(A)')'!------------------------------------------------------------------------'
         WRITE(6,'(A)')'!       Notes'
         WRITE(6,'(A)')'!       ZEFF:  Effective Ion Charge'
         WRITE(6,'(A)')'!       DENS0: Maximum density (10^20 m^-3) IFF TEMPRES>=0'
         WRITE(6,'(A)')'!       TEMPRES:  Controls density profile shape'
         WRITE(6,'(A)')'!       TETI:  TETI Ratio'
         WRITE(6,'(A)')'!       Note defaults below use STELLOPT NE, TE, TI'
         WRITE(6,'(A)')'!------------------------------------------------------------------------'
         WRITE(6,'(A)') '  ZEFF1 = 1'
         WRITE(6,'(A)') '  TEMPRES = -1.0 ! Switch (-1: use STEL prof, >0: use alternative)'
         WRITE(6,'(A)') '  DENS0 = 0.5 ! Density (10^20 m^-3)'
         WRITE(6,'(A)') '  TETI = -1.0 ! TE/TI ratio (set negative to use STEL TI profile)'
         WRITE(6,'(A)') '/'
      END IF
      IF (lorbit) THEN
         WRITE(6,'(A)') '&BEASM3D'
         WRITE(6,'(2X,A)') 'NR = 128'
         WRITE(6,'(2X,A,I3)') 'NPHI = ',180/nfp
         WRITE(6,'(2X,A)') 'NZ = 128'
         WRITE(6,'(2X,A,F6.3)') 'RMAX = ',MAXVAL(MAXVAL(rreal,2),1)+0.1
         WRITE(6,'(2X,A,F6.3)') 'RMIN = ',MINVAL(MINVAL(rreal,2),1)-0.1
         WRITE(6,'(2X,A,F6.3)') 'ZMAX = ',MAXVAL(MAXVAL(zreal,2),1)+0.1
         WRITE(6,'(2X,A,F6.3)') 'ZMIN = ',MINVAL(MINVAL(zreal,2),1)-0.1
         WRITE(6,'(2X,A,F6.3)') 'PHIMIN = ',0.0
         WRITE(6,'(2X,A,F12.10)') 'PHIMAX = ',pi2/nfp
         WRITE(6,'(2X,A)') 'NPOINC = 512'
         WRITE(6,'(2X,A)') 'INT_TYPE = ''LSODE'''
         WRITE(6,'(2X,A)') 'FOLLOW_TOL = 1.0E-9'
         WRITE(6,'(2X,A)') 'T_END_IN = 2*1.0E-2'
         WRITE(6,'(2X,A)') 'R_START_IN(1) = 1.0'
         WRITE(6,'(A)') '/'
      END IF
      IF (lfieldlines) THEN
         WRITE(6,'(A)') '&FIELDLINES_INPUT'
         WRITE(6,'(2X,A)') 'NR = 128'
         WRITE(6,'(2X,A,I3)') 'NPHI = ',180/nfp
         WRITE(6,'(2X,A)') 'NZ = 128'
         WRITE(6,'(2X,A,F6.3)') 'RMAX = ',MAXVAL(MAXVAL(rreal,2),1)+0.5
         WRITE(6,'(2X,A,F6.3)') 'RMIN = ',MINVAL(MINVAL(rreal,2),1)-0.5
         WRITE(6,'(2X,A,F6.3)') 'ZMAX = ',MAXVAL(MAXVAL(zreal,2),1)+0.5
         WRITE(6,'(2X,A,F6.3)') 'ZMIN = ',MINVAL(MINVAL(zreal,2),1)-0.5
         WRITE(6,'(2X,A,F6.3)') 'PHIMIN = ',0.0
         WRITE(6,'(2X,A,F12.10)') 'PHIMAX = ',pi2/nfp
         WRITE(6,'(2X,A)') 'NPOINC = 16'
         WRITE(6,'(2X,A)') 'INT_TYPE = ''LSODE'''
         WRITE(6,'(2X,A)') 'FOLLOW_TOL = 1.0E-9'
         WRITE(6,'(2X,2(A,F6.3))') 'R_START = ',sum(raxis_cc),'  ',MAXVAL(rreal(:,1))
         WRITE(6,'(2X,2(A,F6.3))') 'Z_START = ',0.0,'  ',MAXVAL(zreal(:,1))
         WRITE(6,'(2X,A)') 'PHI_START = 2*0.0'
         WRITE(6,'(2X,A,E20.10)') 'PHI_END = 2*',1000*pi2/nfp
         WRITE(6,'(A)') '/'
      END IF
      CALL FLUSH(6)
!-----------------------------------------------------------------------
!     End Program
!-----------------------------------------------------------------------
      END PROGRAM VMEC2STEL
