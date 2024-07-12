!-----------------------------------------------------------------------
!     Subroutine:    stellopt_bnorm
!     Authors:       Caoxiang Zhu
!     Date:          08/07/18
!     Description:   This subroutine is called to invoke the
!                    BNORM code.
!-----------------------------------------------------------------------
SUBROUTINE stellopt_bnorm(file_str,lscreen)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_input_mod
      USE stellopt_vars
      USE mpi_params     
      
!-----------------------------------------------------------------------
!     Subroutine Parameters
!----------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(256), INTENT(inout)    :: file_str
      LOGICAL, INTENT(inout)           :: lscreen
!-----------------------------------------------------------------------
!     Local Variables
!        istat         Error status
!        iunit         File unit number
!        bnfou_c      B-Normal Fourier coefficients
!----------------------------------------------------------------------
      INTEGER :: istat, nu, nv, mf, nf, md, nd, iunit, m, n, &
                 ivmec, ispline_file
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: bnfou_c

!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      ! Have master run bnorm
      nu = nu_bnorm
      nv = nv_bnorm
      mf = mf_bnorm
      nf = nf_bnorm
      md = md_bnorm
      nd = nd_bnorm
      ! Run BNORM code
      ALLOCATE(bnfou_c(0:mf_bnorm,-nf_bnorm:nf_bnorm),STAT=istat)
      IF (lscreen) WRITE(6,"(A)") '   - Calculating B-Normal File'
      CALL bnormal(nu,nv,mf,nf,md,nd,bnfou,bnfou_c,TRIM(file_str)//'.nc')
      IF (lscreen) WRITE(6,"(A,ES22.12E3)") '      Max. B-Normal: ',MAXVAL(MAXVAL(bnfou,DIM=2),DIM=1)
      IF (lscreen) WRITE(6,"(A,ES22.12E3)") '      MIN. B-Normal: ',MINVAL(MINVAL(bnfou,DIM=2),DIM=1)
      ! WRITE BNORMAL
      CALL safe_open(iunit, istat, 'bnorm.' // TRIM(file_str), 'replace','formatted')
      DO m = 0, mf
         DO n = -nf, nf
            WRITE(iunit,"(1x,2i5,ES22.12E3)") m,n,bnfou(m,n)
         END DO
      END DO
      CLOSE(iunit)
      DEALLOCATE(bnfou_c)
      IF (lscreen) WRITE(6,"(A)") '      Coefficients output to:   '//'bnorm.' // TRIM(file_str)
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
END SUBROUTINE stellopt_bnorm
