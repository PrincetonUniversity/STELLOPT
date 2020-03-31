!-----------------------------------------------------------------------
!     Subroutine:    chisq_coil_bnorm
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          1/23/2015
!     Description:   This subroutine calculates chisqared for coil
!                    optimized bnormal.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_coil_bnorm(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_vals
      USE safe_open_mod
      USE iso_c_binding
      USE mpi_params ! MPI
      
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  target
      REAL(rprec), INTENT(in)    ::  sigma
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
!DEC$ IF DEFINED (COILOPTPP)
      INCLUDE 'coilopt_f.h'
!DEC$ ENDIF
      
!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      INTEGER :: iunit, ik, i, j, nu, nv, ier
      REAL(rprec) :: val
      REAL(rprec), ALlOCATABLE :: uvec(:),vvec(:),bnormf(:),bnorm0(:)
      CHARACTER(16) :: file_str
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (niter >= 0) THEN
         ! Read the initial file
         iunit = 12
         CALL safe_open(iunit,ier,'b_norm_eq_'//TRIM(proc_string_old)//'.dat','old','formatted')
         IF (ier < 0) THEN; iflag=ier; RETURN; END IF
         READ(iunit,*) nu, nv
         IF (nu .ne. nu_bnorm) STOP ' NU /= NU_BNORM in b_norm_eq.dat from COILOPT++'
         IF (nv .ne. nv_bnorm) STOP ' NU /= NU_BNORM in b_norm_eq.dat from COILOPT++'
         ik = nu*nv
         ALLOCATE(uvec(ik),vvec(ik),bnorm0(ik))
         DO i = 1, nu*nv
            READ(iunit,*) uvec(i), vvec(i), bnorm0(i)
         END DO
         CLOSE(iunit)
         DEALLOCATE(uvec,vvec)
         ! Read the final file
         iunit = 12
         CALL safe_open(iunit,ier,'b_norm_final_'//TRIM(proc_string_old)//'.dat','old','formatted')
         IF (ier < 0) THEN; iflag=ier; RETURN; END IF
         READ(iunit,*) nu, nv
         IF (nu .ne. nu_bnorm) STOP ' NU /= NU_BNORM in b_norm_final.dat from COILOPT++'
         IF (nv .ne. nv_bnorm) STOP ' NU /= NU_BNORM in b_norm_final.dat from COILOPT++'
         ik = nu*nv
         ALLOCATE(uvec(ik),vvec(ik),bnormf(ik))
         DO i = 1, nu*nv
            READ(iunit,*) uvec(i), vvec(i), bnormf(i)
         END DO
         CLOSE(iunit)
         ! Write the header
         IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I7.7))') 'COIL_BNORM ',ik,7
         IF (iflag == 1) WRITE(iunit_out,'(A)') 'TARGET  SIGMA  VAL  UVEC  VVEC  BNEQ  BNF '
         DO ik = 1, nu*nv
            mtargets=mtargets+1
            targets(mtargets) = target
            sigmas(mtargets) = sigma
            vals(mtargets) = bnormf(ik) - bnorm0(ik)
            IF (iflag == 1) WRITE(iunit_out,'(7ES22.12E3)') target,sigma,vals(mtargets),uvec(ik),vvec(ik),bnormf(ik),bnorm0(ik)
         END DO
      ELSE
         IF (sigma < bigno) THEN
            DO i = 1, nv_bnorm
               DO j = 1, nu_bnorm
                  mtargets = mtargets + 1
                  IF (niter == -2) target_dex(mtargets)=jtarget_coil_bnorm
               END DO
            END DO
!DEC$ IF DEFINED (COILOPTPP)
            file_str = 'coilopt_params'//CHAR(0)
            CALL init_settings(MPI_COMM_STEL,file_str)
            CALL coilopt_get_numws(numws)
!DEC$ ENDIF
         END IF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_coil_bnorm
