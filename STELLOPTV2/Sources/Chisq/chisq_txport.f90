!-----------------------------------------------------------------------
!     Subroutine:    chisq_txport
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2012
!     Description:   This subroutine calculates chisqared for turbulent
!                    transport optimization.
!-----------------------------------------------------------------------
      SUBROUTINE chisq_txport(target,sigma,niter,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_targets
      USE equil_vals
      USE gist_mod, ONLY: read_gist_namelist
      USE safe_open_mod, ONLY: safe_open
      USE vmec_input, ONLY: nfp_vmec => nfp
!DEC$ IF DEFINED (GENE)
      USE parameters_IO, ONLY: file_extension, read_parameters
      USE gene_scan, ONLY: check_for_scan
      USE geometry, ONLY: geomdir, magn_geometry, geomfile
      USE par_in, ONLY: beta_gene => beta
      USE par_other, ONLY: print_ini_msg
!      USE check_parameters, ONLY: check_for_diagdir ! causes crash
!DEC$ ENDIF
      
      
!-----------------------------------------------------------------------
!     Input/Output Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    ::  target(nsd)
      REAL(rprec), INTENT(in)    ::  sigma(nsd)
      INTEGER,     INTENT(in)    ::  niter
      INTEGER,     INTENT(inout) ::  iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!
!-----------------------------------------------------------------------
      LOGICAL :: mult_par = .true.
      INTEGER :: iunit, ik, i, j, nu, nv, np, ip
      REAL(rprec) :: val
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      np = 1
      IF (TRIM(txport_proxy(1:3)) == 'all') np = 6
      ik = COUNT(sigma < bigno)*np
      IF (iflag == 1) WRITE(iunit_out,'(A,2(2X,I5.5))') 'TXPORT ',ik,4
      IF (iflag == 1) WRITE(iunit_out,'(A,A,A)') 'TARGET  SIGMA  Q  S      (',TRIM(txport_proxy),')'
      IF (niter >= 0) THEN
         DO ik = 1, nsd
            IF (sigma(ik) >= bigno) CYCLE
            DO ip = 1, np
               IF (ALLOCATED(txport_q_all)) THEN
                  val = 0.0
                  nu  = SIZE(txport_q_all,DIM=3)
                  nv  = SIZE(txport_q_all,DIM=4)
                  DO i = 1, nu
                     DO j = 1, nv
                       val = val + txport_q_all(ip,ik,i,j)
                     END DO
                  END DO
                  val = val / (nu*nv)
               ELSE IF (ALLOCATED(txport_q)) THEN
                  val = 0.0
                  nu  = SIZE(txport_q,DIM=2)
                  nv  = SIZE(txport_q,DIM=3)
                  DO i = 1, nu
                     DO j = 1, nv
                       val = val + txport_q(ik,i,j)
                     END DO
                  END DO
                  val = val / (nu*nv)
               ELSE
                  val = bigno
               END IF
               mtargets = mtargets + 1
               targets(mtargets) = target(ik)
               sigmas(mtargets)  = sigma(ik)
               vals(mtargets)    = val
               IF (iflag == 1) WRITE(iunit_out,'(4ES22.12E3)') target(ik),sigma(ik),vals(mtargets),s_txport(ik)
            END DO
         END DO
      ELSE
         IF (alpha_end_txport > pi2/(2*nfp_vmec)) alpha_end_txport = pi2/(2*nfp_vmec)
         DO ik = 1, nsd
            IF (sigma(ik) >= bigno) CYCLE
            DO ip = 1, np
               mtargets = mtargets + 1
               IF (niter == -2) target_dex(mtargets)=jtarget_txport
            END DO
         END DO
         !iunit=12
         !CALL safe_open(iunit,iflag,'input.'//TRIM(id_string),'old','formatted')
         !IF (iflag < 0) RETURN
         !CALL read_gist_namelist (iunit, iflag)
         !IF (iflag < 0) RETURN
         !CLOSE(iunit)
!DEC$ IF DEFINED (GENE)
         IF (TRIM(txport_proxy(1:4)) == 'gene')  THEN
            ! We do this so parameters sets our defaults
            file_extension = ''
            geomdir = '.'
            magn_geometry = 'gist'
            geomfile = '.'
            beta_gene = 0.0
            print_ini_msg = .false.
            CALL check_for_scan(mult_par)
            CALL read_parameters('')
         END IF
!DEC$ ENDIF
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE chisq_txport
