      MODULE liprec

      INTEGER, PARAMETER :: SP=KIND(1.0), DP=KIND(1.0D0)

      INTERFACE li_gbfa

         SUBROUTINE sgbfa (ABD, LDA, N, ML, MU, IPVT, INFO)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER LDA, N, ML, MU, INFO
            INTEGER, DIMENSION(N) :: IPVT
            REAL(SP), DIMENSION(LDA,N) :: ABD
         END SUBROUTINE sgbfa

         SUBROUTINE dgbfa (ABD, LDA, N, ML, MU, IPVT, INFO)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER LDA, N, ML, MU, INFO
            INTEGER, DIMENSION(N) :: IPVT
            REAL(DP), DIMENSION(LDA,N) :: ABD
         END SUBROUTINE dgbfa

         SUBROUTINE sgbfa1 (ABD, LDA, N, ML, MU, IPVT, INFO)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER LDA, N, ML, MU, INFO
            INTEGER, DIMENSION(N) :: IPVT
            REAL(SP), DIMENSION(*) :: ABD
         END SUBROUTINE sgbfa1

         SUBROUTINE dgbfa1 (ABD, LDA, N, ML, MU, IPVT, INFO)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER LDA, N, ML, MU, INFO
            INTEGER, DIMENSION(N) :: IPVT
            REAL(DP), DIMENSION(*) :: ABD
         END SUBROUTINE dgbfa1

      END INTERFACE

      INTERFACE li_gbsl

         SUBROUTINE sgbsl (ABD, LDA, N, ML, MU, IPVT, B, JOB)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER LDA, N, ML, MU, JOB
            INTEGER, DIMENSION(N) :: IPVT
            REAL(SP), DIMENSION(LDA,N) :: ABD
            REAL(SP), DIMENSION(N) :: B
         END SUBROUTINE sgbsl

         SUBROUTINE sgbsl1 (ABD, LDA, N, ML, MU, IPVT, B, JOB)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER LDA, N, ML, MU, JOB
            INTEGER, DIMENSION(N) :: IPVT
            REAL(SP), DIMENSION(*) :: ABD
            REAL(SP), DIMENSION(N) :: B
         END SUBROUTINE sgbsl1

         SUBROUTINE dgbsl (ABD, LDA, N, ML, MU, IPVT, B, JOB)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER LDA, N, ML, MU, JOB
            INTEGER, DIMENSION(N) :: IPVT
            REAL(DP), DIMENSION(LDA,N) :: ABD
            REAL(DP), DIMENSION(N) :: B
         END SUBROUTINE dgbsl

         SUBROUTINE dgbsl1 (ABD, LDA, N, ML, MU, IPVT, B, JOB)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER LDA, N, ML, MU, JOB
            INTEGER, DIMENSION(N) :: IPVT
            REAL(DP), DIMENSION(*) :: ABD
            REAL(DP), DIMENSION(N) :: B
         END SUBROUTINE dgbsl1

      END INTERFACE

      INTERFACE li_gefa

         SUBROUTINE sgefa (A, LDA, N, IPVT, INFO)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER LDA, N, INFO
            INTEGER, DIMENSION(N) :: IPVT
            REAL(SP), DIMENSION(LDA,N) :: A
         END SUBROUTINE sgefa

         SUBROUTINE dgefa (A, LDA, N, IPVT, INFO)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER LDA, N, INFO
            INTEGER, DIMENSION(N) :: IPVT
            REAL(DP), DIMENSION(LDA,N) :: A
         END SUBROUTINE dgefa

         SUBROUTINE sgefa1 (A, LDA, N, IPVT, INFO)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER LDA, N, INFO
            INTEGER, DIMENSION(N) :: IPVT
            REAL(SP), DIMENSION(*) :: A
         END SUBROUTINE sgefa1

         SUBROUTINE dgefa1 (A, LDA, N, IPVT, INFO)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER LDA, N, INFO
            INTEGER, DIMENSION(N) :: IPVT
            REAL(DP), DIMENSION(*) :: A
         END SUBROUTINE dgefa1

      END INTERFACE

      INTERFACE li_gesl

         SUBROUTINE sgesl (A, LDA, N, IPVT, B, JOB)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER, INTENT(IN) :: LDA, N, JOB
            INTEGER, DIMENSION(N), INTENT(IN) :: IPVT
            REAL(SP), DIMENSION(LDA,N), INTENT(IN) :: A
            REAL(SP), DIMENSION(N), INTENT(INOUT) :: B
         END SUBROUTINE sgesl

         SUBROUTINE sgesl1 (A, LDA, N, IPVT, B, JOB)
            INTEGER, PARAMETER :: SP=KIND(1.0)
            INTEGER, INTENT(IN) :: LDA, N, JOB
            INTEGER, DIMENSION(N), INTENT(IN) :: IPVT
            REAL(SP), DIMENSION(*), INTENT(IN) :: A
            REAL(SP), DIMENSION(N), INTENT(INOUT) :: B
         END SUBROUTINE sgesl1

         SUBROUTINE dgesl (A, LDA, N, IPVT, B, JOB)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER, INTENT(IN) :: LDA, N, JOB
            INTEGER, DIMENSION(N), INTENT(IN) :: IPVT
            REAL(DP), DIMENSION(LDA,N), INTENT(IN) :: A
            REAL(DP), DIMENSION(N), INTENT(INOUT) :: B
         END SUBROUTINE dgesl

         SUBROUTINE dgesl1 (A, LDA, N, IPVT, B, JOB)
            INTEGER, PARAMETER :: DP=KIND(1.0D0)
            INTEGER, INTENT(IN) :: LDA, N, JOB
            INTEGER, DIMENSION(N), INTENT(IN) :: IPVT
            REAL(DP), DIMENSION(*), INTENT(IN) :: A
            REAL(DP), DIMENSION(N), INTENT(INOUT) :: B
         END SUBROUTINE dgesl1

      END INTERFACE

      END MODULE liprec
