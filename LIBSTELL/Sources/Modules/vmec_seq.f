      MODULE vmec_seq
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: nseqmax = 100
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!     NSEQ          Number of files to run
!     EXTENSION     List of extension to run
!     NSEQ_SELECT   Order to execute extensions
!     NSEQ_RESTART  Extension to use as restart
!                    Set to 0 to do full run.
!     NSEQ_RESTART(1) is always 0
!-----------------------------------------------
      INTEGER :: nseq
      INTEGER, DIMENSION(nseqmax) :: nseq_select
      CHARACTER(LEN=120), DIMENSION(nseqmax) :: extension
      INTEGER, DIMENSION(nseqmax) :: nseq_restart

      NAMELIST /vseq/ nseq, nseq_select, extension, nseq_restart

      END MODULE vmec_seq
