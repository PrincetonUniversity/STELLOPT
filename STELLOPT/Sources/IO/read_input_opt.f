      SUBROUTINE read_input_opt (file, istat)
      USE optim_params, ONLY: lcoil_geom, lvac_opt
      USE mpi_params
      USE safe_open_mod
      USE vmec_input, ONLY: ntheta, nzeta, mpol, ntor
      USE vacfield_mod, ONLY: lsurf, lisize, lprt, ldisplace, 
     1                        read_vacopt_namelist
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(out) :: istat
      CHARACTER(LEN=*), INTENT(in) :: file
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: iunit=7
C-----------------------------------------------
!
!     OPENS INPUT "FILE" AND READS ALL THE VARIOUS NAMELISTS
!
      CALL safe_open (iunit, istat, file, 'old', 'formatted')
      IF (istat .ne. 0) GOTO 100

      CALL read_namelist (iunit, istat, 'indata')
      IF (istat .ne. 0) THEN
         IF (myid .eq. master) THEN                                      !MPI
         PRINT *,' indata namelist read error in read_input_opt,',
     1       ' istat = ', istat
         ENDIF                                                           !MPI
         GOTO 100
      END IF
      IF (nzeta .eq. 0) nzeta = 2*ntor + 4                               
      IF (ntor .eq. 0) nzeta = 1
      IF (ntheta .eq. 0) ntheta = 2*mpol + 6

      CALL read_namelist (iunit, istat, 'optimum')
      IF (istat .ne. 0) THEN
         IF (myid .eq. master) THEN                                      !MPI
            PRINT *,' optimum namelist read error in read_input_opt,',
     1      ' istat = ', istat
         ENDIF                                                           !MPI
         GOTO 100
      END IF

      CALL read_namelist (iunit, istat, 'bootin')

      CALL read_namelist (iunit, istat, 'ga_de')

      istat = 0
      IF (lcoil_geom) CALL read_namelist (iunit, istat, 'coilsin')
      IF (istat .ne. 0) THEN
         IF (myid .eq. master) THEN                                      !MPI
            PRINT *,' coilsin namelist read error in read_input_opt,',
     1      ' istat = ', istat
         ENDIF                                                           !MPI
         STOP
      END IF

      istat = 0
      IF (lvac_opt) THEN
         CALL read_vacopt_namelist (iunit, istat)
         IF (istat .ne. 0) THEN
            IF (myid .eq. master) THEN                                      !MPI
            PRINT *,' vacopt namelist read error in read_input_opt,',
     1      ' istat = ', istat
            ENDIF                                                           !MPI
            STOP
         END IF
    
         lsurf = .false.;  lisize = .false.; lprt = .false.
         ldisplace = .true.

      END IF

 100  CONTINUE
   
      IF (istat .eq. 0) THEN
         CLOSE (iunit, status='delete')
      ELSE
         PRINT *,' read_input_opt read error, istat= ',istat
         CLOSE (iunit)
      END IF

      END SUBROUTINE read_input_opt
