      MODULE FFTPACK

      INTERFACE CFTFAX
          SUBROUTINE cftfax_g(n, ifax, trigs)
          USE stel_kinds
          INTEGER n
          INTEGER, DIMENSION(13) :: ifax
          REAL(rprec), DIMENSION(*) :: trigs
          END SUBROUTINE cftfax_g
      END INTERFACE

      INTERFACE CFTRIG
          SUBROUTINE cftrig_g (n, trigs)
          USE stel_kinds
          INTEGER n
          REAL(rprec), DIMENSION(*) :: trigs
          END SUBROUTINE CFTRIG_G
      END INTERFACE

      INTERFACE FACT
          SUBROUTINE fact_g(n, ifax)
          USE stel_kinds
          INTEGER n
          INTEGER, DIMENSION(13) :: ifax
          END SUBROUTINE fact_g
      END INTERFACE

      INTERFACE FFTFAX
          SUBROUTINE fftfax_g(n, ifax, trigs)
          USE stel_kinds
          INTEGER n
          INTEGER, DIMENSION(13) :: ifax
          REAL(rprec), DIMENSION(*) :: trigs
          END SUBROUTINE fftfax_g
      END INTERFACE

      INTERFACE FFTRIG
          SUBROUTINE fftrig_g (trigs, n, mode)
          USE stel_kinds
          INTEGER n, mode
          REAL(rprec), DIMENSION(*) :: trigs
          END SUBROUTINE FFTRIG_G
      END INTERFACE

      END MODULE FFTPACK
