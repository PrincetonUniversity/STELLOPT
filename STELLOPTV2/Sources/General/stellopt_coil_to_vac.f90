!-----------------------------------------------------------------------
!     Subroutine:    stellopt_coil_to_vac
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/15/2016
!     Description:   This subroutine reads a coils file and generates
!                    the appropriate vacuum grid file.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_coil_to_vac(lscreen,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime, ONLY: proc_string, id_string
      USE stellopt_vars, ONLY: equil_type, lcreate_coils
      USE write_mgrid, only: mgrid_ext, lstell_sym, kp
      USE makegrid_global, only: task, lscreen_mgrid => lscreen
      USE vmec_input, ONLY:  lfreeb, nfp, ntor, mpol, rbc, zbs, &
                              INIT_AXIS_MIDPOINT, raxis_cc, zaxis_cs
      USE read_wout_mod, ONLY: mnmax, ns, xm, xn, rmnc, zmns, isigng
      IMPLICIT NONE
      
!-----------------------------------------------------------------------
!     Input Variables
!        lscreen   Terminal output
!----------------------------------------------------------------------
      LOGICAL, INTENT(inout) :: lscreen
      INTEGER, INTENT(inout) :: iflag

!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      LOGICAL, SAVE :: lfirst_pass = .TRUE.
      INTEGER :: n, m

!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      ! Note that this subroutine is only called in a VMEC loop
      ! First we USE RBC/ZBS to initialize the wout variables from
      ! the RBC/ZBS data
      mnmax = 0
      DO n = -ntor,ntor
         DO m = 0, mpol+1
            IF ((rbc(n,m)>0) .or. (zbs(n,m)>0)) mnmax = mnmax + 1
         END DO
      END DO
      IF (ALLOCATED(xm)) DEALLOCATE(xm)
      IF (ALLOCATED(xn)) DEALLOCATE(xn)
      IF (ALLOCATED(rmnc)) DEALLOCATE(rmnc)
      IF (ALLOCATED(zmns)) DEALLOCATE(zmns)
      ns = 2
      ALLOCATE(xm(mnmax),xn(mnmax),rmnc(mnmax,ns),zmns(mnmax,ns))
      mnmax = 0
      DO n = -ntor,ntor
         DO m = 0, mpol+1
            IF ((rbc(n,m)>0) .or. (zbs(n,m)>0)) THEN
               mnmax = mnmax + 1
               xm(mnmax) = m
               xn(mnmax) = n
               rmnc(mnmax,ns) = rbc(n,m)
               zmns(mnmax,ns) = zbs(n,m)
            END IF
         END DO
      END DO
      CALL INIT_AXIS_MIDPOINT()
      DO m = 1, mnmax
         n = xn(m)
         IF (xm(m)==0 .and. n>=0) THEN
            rmnc(m,1) = raxis_cc(n)
            zmns(m,1) = zaxis_cs(n)
         ENDIF
      END DO
      xn = xn * nfp
      ! Now generate the coil
      IF (lcreate_coils) CALL stellopt_generate_coils(lscreen,iflag)
      ! Read MGRID namelist from input.EXT if first time through
      IF (lfirst_pass) CALL namelist_input_makegrid('input.'//TRIM(id_string))
      ! First run VMEC in fixed boundary
      !lfreeb = .FALSE.
      !CALL stellopt_run_vmec(lscreen,iflag)
      ! Generate the MGRID
      task='MGRID'
      mgrid_ext=TRIM(proc_string)
      lscreen_mgrid = lscreen
      kp = nfp
      CALL task_mgrid()
      ! Reset free boundary
      lfreeb = .TRUE.
      lfirst_pass = .FALSE.
      ! Now re-read the 
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_coil_to_vac
