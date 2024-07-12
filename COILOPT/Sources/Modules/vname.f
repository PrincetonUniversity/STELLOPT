      MODULE Vname
      USE stel_constants
#if defined(RISC) || defined(WIN32)
      INTEGER, PARAMETER :: num_processors    = 1
      INTEGER, PARAMETER :: num_levmar_params = 10
#else
      INTEGER, PARAMETER :: num_processors    = 8
      INTEGER, PARAMETER :: num_levmar_params = 10
#endif
      CHARACTER :: input_data_file*200
      CHARACTER :: for05_file*200, for05_file_new*206, path*180
      CHARACTER :: extension*150, wout_file*200
      REAL(rprec), PARAMETER :: min_chisq = 1000000
      REAL(rprec) :: chisq_min = min_chisq
      LOGICAL :: lgeom_only = .FALSE., lplasma
      END MODULE Vname
