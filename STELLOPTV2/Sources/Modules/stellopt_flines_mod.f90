!-----------------------------------------------------------------------
!     Module:        stellopt_flines_mod
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          07/12/2025
!     Description:   This module contains routines for following
!                    field lines.
!-----------------------------------------------------------------------
      MODULE stellopt_flines_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE biotsavart, ONLY: cleanup_biotsavart, parse_coils_file, &
            bfield, bsc_coil
      
!-----------------------------------------------------------------------
!     Module Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, PARAMETER, PRIVATE :: pi2 = 6.283185482025146D+00
      DOUBLE PRECISION, PARAMETER, PRIVATE :: zero = 0.0D+00
      DOUBLE PRECISION, PARAMETER, PRIVATE :: one = 1.0D+00
      DOUBLE PRECISION, PARAMETER, PRIVATE :: REL_TOL = 1.0D-9
      DOUBLE PRECISION, PARAMETER, PRIVATE :: ABS_TOL = 1.0D-9
      
      
!-----------------------------------------------------------------------
!     Module SUBROUTINES/FUNCTIONS
!-----------------------------------------------------------------------
      CONTAINS

      SUBROUTINE stellopt_follow_init(coil_string)
      IMPLICIT NONE
      CHARACTER(256), INTENT(INOUT) :: coil_string
      CALL cleanup_biotsavart()
      CALL parse_coils_file(TRIM(coil_string))
      RETURN
      END SUBROUTINE stellopt_follow_init

      SUBROUTINE stellopt_follow_single(R0,PHI0,Z0,R1,PHI1,Z1,ier)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(INOUT) :: R0,PHI0,Z0,R1,PHI1,Z1
      DOUBLE PRECISION, INTENT(INOUT) :: ier
      INTEGER     :: neqs, iopt, mf, itol, itask, istate,&
                     lrw, liw
      INTEGER     :: iwork(20)
      DOUBLE PRECISION :: q(2)
      DOUBLE PRECISION :: w(52) !20*16*neqs
      IF (ier .ne. 0) RETURN
      neqs = 2; itol = 1; itask = 1; istate = 1
      w = zero; lrw = 52
      iwork = 0; liw = 20
      mf = 10
      q(1) = R0; q(2) = Z0
      iopt = 1; iwork(6) = 50000 ! Need this because 500 over a field period is not much
      CALL DLSODE(stellopt_fbline,neqs,q,PHI0,PHI1,itol,REL_TOL,ABS_TOL,&
                  itask,istate,iopt,w,lrw,iwork,liw,stellopt_fbline_jacobian,mf)
      R1 = q(1)
      Z1 = q(2)
      IF (istate < -1) ier = istate
      RETURN
      END SUBROUTINE stellopt_follow_single

      SUBROUTINE stellopt_fbline(neq,phi,q,qdot)
      IMPLICIT NONE
      DOUBLE PRECISION :: phi, q(2), qdot(2)
      INTEGER :: neq, ier
      DOUBLE PRECISION :: br,bphi,bz
      ier = 0
      qdot = 0
      CALL bfield(q(1),phi,q(2),br,bphi,bz,ier)
      qdot(1) = q(1)*br/bphi
      qdot(2) = q(1)*bz/bphi
      RETURN
      END SUBROUTINE stellopt_fbline

      SUBROUTINE stellopt_fbline_jacobian(neq,phi,q,ml,mp,pd,nrpd)
      IMPLICIT NONE
      INTEGER          :: neq, nrpd, ml, mp
      DOUBLE PRECISION :: phi, q(6), pd(nrpd,6)
      pd = zero
      RETURN
      END SUBROUTINE stellopt_fbline_jacobian

!-----------------------------------------------------------------------
!     End Module
!-----------------------------------------------------------------------
      END MODULE stellopt_flines_mod