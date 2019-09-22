!-----------------------------------------------------------------------
!     Subroutine:    stellopt_renorm
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          09/20/2019
!     Description:   This routine renomalizes the sigmas so that each
!                    chi-squared functional returns a value = 1 over
!                    that diagnostic
!                    f = (xt-x)/sigmas
!                    1 = chisq*c
!                    chisq = 1/c
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_renorm(m1,fvec)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_input_mod
      USE stellopt_vars
      USE mpi_params
      USE mpi_inc
      
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      IMPLICIT NONE
      !LOGICAL ::  lrestart
      INTEGER, INTENT(in) :: m1
      REAL(rprec), INTENT(in) :: fvec(m1)

      ! For renomalization
      INTEGER :: nuniq, cid,i,j
      INTEGER, ALLOCATABLE :: iddex(:)
      REAL(rprec) :: temp
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      ! Count number of unique types
      cid = 0
      nuniq = 0
      DO i = 1, mtargets
         IF (target_dex(i).ne.cid) THEN
            nuniq=nuniq+1
            cid = target_dex(i)
         ENDIF
      END DO
      ! Construct array of unique IDs
      ALLOCATE(iddex(nuniq))
      cid = 0; j = 1
      DO i = 1, mtargets
         IF (target_dex(i).ne.cid) THEN
            cid = target_dex(i)
            iddex(j) = target_dex(i)
            j = j+1
         ENDIF
      END DO
      ! For each ID renormalize
      DO i = 1, nuniq
         temp = 0
         DO j = 1, mtargets
            IF (iddex(i)==target_dex(j)) THEN
               temp = temp + fvec(j)*fvec(j)
            END IF
         END DO
         IF (temp == 0) temp = 1
         temp = 1./sqrt(temp)
         ! now renormalize sigmas (BIG LIST)
         SELECT CASE(iddex(i))
            !CASE(jtarget_curtor)
            !   sigma_curtor = sigma_curtor/temp
            CASE(jtarget_separatrix)
               WHERE(sigma_separatrix<bigno) sigma_separatrix = sigma_separatrix/temp
            CASE(jtarget_ne)
               WHERE(sigma_ne<bigno_ne) sigma_ne = sigma_ne/temp
            CASE(jtarget_te)
               WHERE(sigma_te<bigno) sigma_te = sigma_te/temp
            CASE(jtarget_ti)
               WHERE(sigma_ti<bigno) sigma_ti = sigma_ti/temp
            CASE(jtarget_line_ne)
               WHERE(sigma_ne_line<bigno_ne) sigma_ne_line = sigma_ne_line/temp
            CASE(jtarget_line_te)
               WHERE(sigma_te_line<bigno) sigma_te_line = sigma_te_line/temp
            CASE(jtarget_line_ti)
               WHERE(sigma_ti_line<bigno) sigma_ti_line = sigma_ti_line/temp
            CASE(jtarget_xics)
               WHERE(sigma_xics<bigno) sigma_xics = sigma_xics/temp
            CASE(jtarget_xics_bright)
               WHERE(sigma_xics_bright<bigno) sigma_xics_bright = sigma_xics_bright/temp
            CASE(jtarget_xics_w3)
               WHERE(sigma_xics_w3<bigno) sigma_xics_w3 = sigma_xics_w3/temp
            CASE(jtarget_xics_v)
               WHERE(sigma_xics_v<bigno) sigma_xics_v = sigma_xics_v/temp
            CASE(jtarget_mse)
               WHERE(sigma_mse<bigno) sigma_mse = sigma_mse/temp
            CASE(jtarget_faraday)
               WHERE(sigma_faraday<bigno) sigma_faraday = sigma_faraday/temp
            CASE(jtarget_sxr)
               WHERE(sigma_sxr<bigno) sigma_sxr = sigma_sxr/temp
            CASE(jtarget_ece)
               WHERE(sigma_ece<bigno) sigma_ece = sigma_ece/temp
            CASE(jtarget_bprobe)
               WHERE(sigma_bprobe<bigno) sigma_bprobe = sigma_bprobe/temp
            CASE(jtarget_fluxloop)
               WHERE(sigma_fluxloop<bigno) sigma_fluxloop = sigma_fluxloop/temp
            CASE(jtarget_segrog)
               WHERE(sigma_segrog<bigno) sigma_segrog = sigma_segrog/temp
            CASE DEFAULT
               WRITE(6,'(A,I3.3,A)') '!!! JTARGET=',iddex(i),' not supported'
               CALL write_targets(6,iddex(i))
         END SELECT
      END DO
      DEALLOCATE(iddex)
        
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_renorm
