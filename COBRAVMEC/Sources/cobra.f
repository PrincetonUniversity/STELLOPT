      PROGRAM cobra
c___________________________________________________________________________________
c                                                                                   |
c              COBRA (COde for Ballooning Rapid Analysis)                           |
c                     ==       =          =     =                                   |
c                                                                                   |
c      VERSION: 4.1  (FORTRAN 90: fixed FORMAT; standard I/O; VMEC coordinates)     |
c                    (Accepts ASYMMETRIC input from VMEC via LASYM)                 |
c                    (Accepts initial line location as geometrical or (ALHPA,ANGLE) |
c                    (Good for LRFP=T or F)                                         |
c      Last update: 03/27/13                                                        |
c                                                                                   |
c      AUTHORS:   R. Sanchez                                                        |
c             Universidad Carlos III de Madrid, Leganes 28911, SPAIN.               |
c                 S.P. Hirshman                                                     |
c             Oak Ridge National Laboratory, Oak Ridge, TN 37831-9701, USA          |
c                                                                                   |
c      REFERENCES:                                                                  |
c                                                                                   |
c        1."COBRA: and optimized code for fast analysis of ideal ballooning         |
c        stability of 3-D magnetic equilibria", R. Sanchez, S.P. Hirshman, J.C.     |
c        Whitson and A.S. Ware, submitted to Journal of Computational Physics (1999)|
c                                                                                   |
c        2."Improved magnetic coordinate representation for ideal ballooning        |
c        stability calculations with the COBRA code", R. Sanchez, S.P. Hirshman     |
c        H.V. Wong, to be submitted to Journal of Computational physics (2000)      |
c                                                                                   |
c      DISCLAIMER:  This code is under development by R.Sanchez at the Departamento |
c        de Fisica, Universidad Carlos III de Madrid, SPAIN. As a BETA version, the |
c        code is supplied on "as it is" basis, and the non-existence of "bugs" of   |
c        ANY kind is NOT guaranteed. Any problem or comment should be reported to   |
c        R. Sanchez at raul.sanchez@uc3m.es                                         |
c___________________________________________________________________________________|
c
      USE stel_kinds
      USE normalize_data
      USE ballooning_data
      USE readin_data
      USE safe_open_mod
      USE date_and_computer
      USE fmesh_quantities
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: unit_cobra = 12, input_cobra = 14
      CHARACTER*(*), PARAMETER :: banner =
     1   'THIS IS THE ASYMMETRIC COBRA BALLOONING CODE Version '
      CHARACTER*(5), PARAMETER :: cobra_version = '4.10'
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: grate
      REAL(rprec) :: t1,t2
      INTEGER :: j, nlis, i, iunit_out, numargs, iunit_in, ij, jk
      INTEGER :: istat, ns_surf, nini_zeta, nini_theta
      INTEGER,DIMENSION(:),allocatable:: in_surf, in_idx, bsurf
      INTEGER :: imon
      CHARACTER*120 :: extension
      CHARACTER*10 :: date0, time0, zone0
      CHARACTER*40 :: dateloc, arg1, arg2
!-----------------------------------------------
!
!     Read command line argument to get input file or sequence file
!
      lscreen = .true.

      CALL getcarg(1, arg1, numargs)
      IF (numargs .gt. 1) CALL getcarg(2, arg2, numargs)

      IF (numargs .lt. 1) THEN
         STOP 'Invalid command line'
      ELSE IF (arg1 .eq. '-h' .or. arg1 .eq. '/h') THEN
         PRINT *,
     1   ' ENTER INPUT FILE NAME ON COMMAND LINE'
         PRINT *,' For example: xcobravmec in_cobra.tftr'
         PRINT *
         PRINT *,' Optional command line argument (T/F):'
         PRINT *,' xcobravmec input_file (T or F)'
         PRINT *
         PRINT *,' where F will suppress screen output'
         STOP
      ELSE IF (numargs .gt. 1) THEN
         IF (arg2(1:1).eq.'f' .or. arg2(1:1).eq.'F') lscreen = .false.
      ENDIF

!      INPUT FILE
!      1st line:   extension of WOUT file
!      2nd line:   k_w, kth  (No. helical wells; mode number)
!      3rd line:   l_geom_input, l_tokamak_input
!      4th line:   nini_zeta (No. initial toroidal angles; or ALPHA labels if L_GEOM_INPUT=F .AND. L_TOKAMAK_INPUT=T)
!      5th line:   init_zeta_v (vector of initial toroidal angles (or labels), in degrees)
!      6th line:   nini_theta (No. initial polidal angles; or ALPHA labels if L_GEOM_INPUT=F .AND. L_TOKAMAK_INPUT=F)  
!      7th line:   init_theta_v (vector of initial poloidal angles (or labels), in degrees)
!      8th line:   ns_surf (No. surfaces where growth rate is to be computed)
!      9th line:   ns_v (vector of surfaces where growth rate is to be computed)


      iunit_in = input_cobra
      CALL safe_open(iunit_in, istat, TRIM(arg1), 'old', 'formatted')

      IF (istat .ne. 0) STOP ' Error opening input file in COBRA'

      READ (iunit_in, *, iostat=istat) extension                                ! Read extension of WOUT VMEC file
      READ (iunit_in, *, iostat=istat) k_w, kth                                 ! Number of helical wells; mode number (=1, most unstable)
      READ (iunit_in, *, iostat=istat) l_geom_input, l_tokamak_input
      READ (iunit_in, *, iostat=istat) nini_zeta                                ! Read number initial ZETA points (if STELLOPT or GEOM)  or labels (otherwise)
      ALLOCATE (init_zeta_v(nini_zeta), stat=istat)
      READ (iunit_in, *, iostat=istat) (init_zeta_v(j), j=1,nini_zeta)          ! Read vector of ballooning initial toroidal angles
      READ (iunit_in, *, iostat=istat) nini_theta                               ! Read number initial THETA points (if STELLOPT or GEOM)  or labels (otherwise)
      ALLOCATE (init_theta_v(nini_theta), stat=istat)
      READ (iunit_in, *, iostat=istat) (init_theta_v(j),j=1,nini_theta)         ! Read vector of ballooning initial poloidal angles
      READ (iunit_in, *, iostat=istat) ns_surf                                  ! READ no surfaces where growth rates are to be computed
      ALLOCATE (in_surf(ns_surf), in_idx(ns_surf))
      READ (iunit_in, *, iostat=istat) in_surf(1:ns_surf)                       ! READ vector of desried surfaces on the full grid
      CLOSE (iunit_in)

      iunit_out = unit_cobra
      CALL safe_open(iunit_out, istat, 'cobra_grate.'//extension,               ! Output file
     1     'replace', 'formatted')

      IF (istat .ne. 0) THEN
          WRITE (iunit_out, *) ' Error opening COBRA_GRATE'
          CLOSE (iunit_out)
          STOP ' Error opening COBRA_GRATE'
      END IF

      IF (lscreen) THEN
         PRINT 48
         PRINT *
         CALL date_and_time(date0,time0,zone0)
         READ (date0(5:6),'(i2)') imon
         WRITE (dateloc,100) months(imon),date0(7:8),date0(1:4),
     1      time0(1:2),time0(3:4),time0(5:6)
         WRITE (*,'(1x,2a,/,1x,2a)') banner, cobra_version,
     1      computer, dateloc
         PRINT *
      END IF
 100  FORMAT('DATE = ',a3,' ',a2,',',a4,' ',' TIME = ',2(a2,':'),a2)

      nlis = 0
      DO i = 1, ns_surf
         IF (in_surf(i) .gt. 1) THEN
            nlis = nlis+1
            in_idx(nlis) = in_surf(i)
         ENDIF
      ENDDO

      ALLOCATE (bsurf(nlis))
      bsurf(:nlis) = in_idx(:nlis)

!...  Read equilibrium data and form surface list

      CALL order_input(extension, nlis, bsurf, istat)
      IF (istat .ne. 0) THEN
         extension = ' Error reading WOUT.' // TRIM(extension)
     1         // ' in COBRAVMEC'
         WRITE (iunit_out, *) TRIM(extension)
         CLOSE (iunit_out)
         PRINT *, TRIM(extension)
         STOP
      END IF

!..   solve ballooning equation

      ALLOCATE (grate(ns_cob), radios(ns_cob))                                  ! ALLOCATE growth rate vector
      CALL second0(t1)

      theta: DO ij = 1, nini_theta                                             ! Loop over all poloidal, toroidal initial pairs
        zeta: DO jk = 1, nini_zeta

          IF (l_geom_input) THEN                                               ! Deal with different types of input                                    
            init_zeta = init_zeta_v(jk)
            init_zeta = MOD(init_zeta, 360._dp)
            init_theta = init_theta_v(ij)
            init_theta = MOD(init_theta, 360._dp)
          ELSE
            IF (l_tokamak_input) THEN
              init_alpha_tok = init_zeta_v(jk)
              init_alpha_tok = MOD(init_alpha_tok, 360._dp)
              init_thetak_tok = init_theta_v(ij)
              init_thetak_tok = MOD(init_thetak_tok, 360._dp)
            ELSE
              init_zetak_st = init_zeta_v(jk)
              init_zetak_st = MOD(init_zetak_st, 360._dp)
              init_alpha_st = init_theta_v(ij)
              init_alpha_st = MOD(init_alpha_st, 360._dp)
            ENDIF
          ENDIF
          lfail_balloon = .false.
          CALL get_ballooning_grate(grate)                                     ! get growth rates on desired surfaces
          CALL second0(t2)

!...   generate output

          IF (lscreen) THEN
            IF (l_geom_input) THEN
               WRITE(*,120) init_zeta, init_theta, t2-t1
            ELSE
              IF (l_tokamak_input) THEN
                WRITE(*,121) init_alpha_tok, init_thetak_tok, 
     1            t2-t1
              ELSE
                WRITE(*,122) init_alpha_st, init_zetak_st, t2-t1
              ENDIF
            ENDIF
               WRITE(*,48)
          END IF

          IF (l_geom_input) THEN
            WRITE (iunit_out, 12) init_zeta, init_theta, nlis
          ELSE
            IF (l_tokamak_input) THEN
              WRITE (iunit_out, 12) init_alpha_tok, init_thetak_tok, 
     1          nlis
            ELSE
              WRITE (iunit_out, 12) init_alpha_st, init_zetak_st, 
     1          nlis
            ENDIF
          ENDIF
          WRITE (iunit_out, 14) (bsurf(j), radios(bsurf(j)),
     1      grate(bsurf(j)),j=1,nlis)                                         ! WRITE output file

        ENDDO zeta
      ENDDO theta

12    FORMAT(1p,2e10.3,i5)
14    FORMAT(i4, 1p,e16.8, 1p,e16.8)
120   FORMAT(/,'ZETA0 = ', 1pe10.3,' THETA0 = ', 1pe10.3,
     1   ' TIME IN COBRA CODE:',1pe10.2,' SEC')
121   FORMAT(/,'ALPHA_TOK = ', 1pe10.3,' THETAK = ', 1pe10.3,
     1   ' TIME IN COBRA CODE:',1pe10.2,' SEC')
122   FORMAT(/,'ALPHA_ST = ', 1pe10.3,' ZETAK = ', 1pe10.3,
     1   ' TIME IN COBRA CODE:',1pe10.2,' SEC')

48    FORMAT('====================================================')
      CLOSE (unit=iunit_out)

      DEALLOCATE (grate, bsurf, in_surf, in_idx, init_theta_v,
     1       init_zeta_v, radios)     
      DEALLOCATE (hiota, hpres, hphip, rmncf, zmnsf, list,
     1   lmnsh, bmnch, bsupumnch, bsupvmnch, mercierf, xn_v, xm_v)
      IF (lasym_v) THEN                                                        ! 110909 RS = for ASYMMETRIC INPUT
        DEALLOCATE (rmnsf, zmncf, lmnch, bmnsh, bsupumnsh, bsupvmnsh)
      ENDIF

      END PROGRAM cobra
