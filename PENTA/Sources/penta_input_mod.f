      module penta_input_mod
      use penta_kind_mod
      implicit none      
      INTEGER(iknd), parameter ::  num_ion_max = 20_iknd    !Maximum number of ion species to be read from namelist file
      INTEGER(iknd) :: num_ion_species ! Number of ion species
      REAL(rknd), DIMENSION(num_ion_max) :: Z_ion_init,miomp_init  

      namelist / ion_params / num_ion_species,Z_ion_init,miomp_init
      contains

      SUBROUTINE init_ion_params_nml
      IMPLICIT NONE
      num_ion_species = 1
      Z_ion_init = 1.0_rknd
      miomp_init = 1.0_rknd
      END SUBROUTINE init_ion_params_nml

      SUBROUTINE read_ion_params_nml(filename,istat)
      USE safe_open_mod
      IMPLICIT NONE
      ! Subroutine parameters
      CHARACTER(256), INTENT(INOUT) :: filename
      INTEGER(iknd), INTENT(INOUT) :: istat
      ! Local Variables
      LOGICAL :: lexist
      INTEGER(iknd) :: iunit
      CHARACTER(256) :: filename_local

      ! Handle filename
      filename_local = TRIM(filename)

      !Check which file is available
      INQUIRE(FILE=TRIM(filename_local),EXIST=lexist)
      IF (.not. lexist) THEN
         WRITE(6,*) '   -Could not open input file trying ion_params'
         filename_local = 'ion_params'
         INQUIRE(FILE=TRIM(filename_local),EXIST=lexist)
         IF (.not. lexist) THEN
            WRITE(6,*) '   -Could not open ion_params file'
            istat = -1
            RETURN
         ENDIF
      ENDIF

      ! Open and read the file
      CALL safe_open(iunit,istat,TRIM(filename_local),'old',
     1                  'formatted')
      IF (istat /= 0) THEN
         WRITE(6,*) 'ERROR: Could no open: ',TRIM(filename_local)
         RETURN
      ENDIF
      READ(iunit,NML=ion_params,IOSTAT=istat)
      IF (istat /= 0) THEN
         WRITE(6,*) 'ERROR: ion_params in ',TRIM(filename_local)
         !backspace(iu_nl)
         !read(iu_nl,fmt='(A)') line
         !write(6,'(A)') 'Invalid line in namelist: '//TRIM(line)
         RETURN
      ENDIF
      CLOSE(iunit)
      RETURN
      END SUBROUTINE read_ion_params_nml

      SUBROUTINE write_ion_params_nml(iunit)
      IMPLICIT NONE
      INTEGER(iknd), INTENT(inout) :: iunit
      CHARACTER(LEN=*), PARAMETER :: outint  = "(2X,A,1X,'=',1X,I0)"
      INTEGER(iknd) :: k
      WRITE(iunit,'(A)') '&ION_PARAMS'
      WRITE(iunit,'(A)') 
     1'----------------------------------------------------------------'
      WRITE(iunit,outint) 'NUM_ION_SPECIES',num_ion_species
      WRITE(iunit,"(2X,A,1X,'=',4(1X,ES22.12E3))") 
     1'Z_ION_INIT',(Z_ion_init(k), k=1,num_ion_species)
      WRITE(iunit,"(2X,A,1X,'=',4(1X,ES22.12E3))") 
     1'MIOMP_INIT',(miomp_init(k), k=1,num_ion_species)
      WRITE(iunit,'(A)') '/\n'
      END SUBROUTINE write_ion_params_nml

      SUBROUTINE write_ion_params_namelist_byfile(filename)
      USE safe_open_mod
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER(iknd) :: iunit, istat
      LOGICAL :: lexists
      
      iunit = 100
      istat = 0
      INQUIRE(FILE=TRIM(filename),exist=lexists)
      IF (lexists) THEN
         OPEN(unit=iunit, file=TRIM(filename), 
     1        iostat=istat, status="old", position="append")
      ELSE
         OPEN(unit=iunit, file=TRIM(filename), 
     1        iostat=istat, status="new")
      END IF
      IF (istat .ne. 0) RETURN
      CALL write_ion_params_nml(iunit)
      CLOSE(iunit)

      RETURN
      END SUBROUTINE write_ion_params_namelist_byfile
      end module penta_input_mod