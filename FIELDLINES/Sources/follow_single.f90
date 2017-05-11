!-----------------------------------------------------------------------
!     Function:      follow_single
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/21/2012
!     Description:   This subroutine follows a single field line and
!                    returns it's final position.
!-----------------------------------------------------------------------
      SUBROUTINE follow_single(r0,phi0,z0,phi1,tan_map)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE fieldlines_grid
      USE fieldlines_runtime
!-----------------------------------------------------------------------
!     Input Variables
!          r0         Initial R position
!          phi0       Initial phi position
!          z0         Initial z position
!          phi1       Final phi position
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in)    :: phi0, phi1
      REAL(rprec), INTENT(inout) :: r0, z0
      REAL(rprec), INTENT(out), OPTIONAL :: tan_map(4)
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      INTEGER     :: neqs, ier, iopt, mf, itol, itask, istate,&
                     lrw, liw
      INTEGER     :: iwork(20)
      DOUBLE PRECISION :: tol_nag
      DOUBLE PRECISION :: q(6), w(6*21+28), w2(20+16*6)
      DOUBLE PRECISION :: phi0_nag, phi1_nag, rtol
      DOUBLE PRECISION :: atol(6)
      DOUBLE PRECISION :: rkh_work(6,2)
      CHARACTER*1 :: relab
!-----------------------------------------------------------------------
!     External Functions
!          fblin_nag            RHS of ODE integrator (for NAG)
!          fblin_lsode          RHS of ODE integrator (for LSODE)
!          fblin_rkh68          RHS of ODE integrator (for RKH68)
!          out_fieldlines_nag   Fieldline output (for NAG)
!          D02CJF               NAG ODE Solver
!          D02CJW               NAG Dummy function
!          D02CJX               NAG Dummy function
!          jacobian_lsode       Jacobian function (for LSODE, not currently utilized)
!-----------------------------------------------------------------------
      EXTERNAL fblin_tanmap_nag, out_fieldlines_nag, D02CJF, D02CJW, D02CJX
      EXTERNAL fblin_tanmap_lsode, jacobian_tanmap_lsode
      EXTERNAL fblin_tanmap_rkh68
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      neqs = 6
      q(1) = r0
      q(2) = z0
      q(3) = 1
      q(4) = 0
      q(5) = 0
      q(6) = 1
      relab = "M"
      mf = 10
      tol_nag = follow_tol
      phi0_nag = phi0
      phi1_nag = phi1
      ier = 1
      SELECT CASE (TRIM(int_type))
         CASE ("NAG")
!DEC$ IF DEFINED (NAG)
            CALL D02CJF(phi0_nag,phi1_nag,neqs,q,fblin_tanmap_nag,tol_nag,relab,D02CJX,D02CJW,w,ier)
            IF (ier < 0) CALL handle_err(D02CJF_ERR,'follow_single',ier)
!DEC$ ELSE
            ier = -1
            CALL handle_err(NAG_ERR,'follow_single',ier)
!DEC$ ENDIF 
         CASE ("RKH68")
            ier = 0
            CALL drkhvg(phi0_nag,q,neqs,phi1_nag-phi0_nag,2,fblin_tanmap_rkh68,rkh_work,iopt,ier)
            IF (ier < 0) CALL handle_err(RKH68_ERR,'follow_single',ier)
            q(1)=rkh_work(1,2)
            q(2)=rkh_work(2,2)
            q(3)=rkh_work(3,2)
            q(4)=rkh_work(4,2)
            q(5)=rkh_work(5,2)
            q(6)=rkh_work(6,2)
         CASE ("LSODE","DLSODE")
            itol = 2; rtol = follow_tol; atol(:) = follow_tol
            iopt = 0 ! No optional output
            lrw = 20 + 16 * 6
            liw = 20
            ier = 0
            w2=0; iwork = 0; itask = 1; istate = 1;
            iopt = 1; iwork(6) = 50000 ! Need this because 500 over a field period is not much
            CALL DLSODE(fblin_tanmap_lsode,neqs,q,phi0_nag,phi1_nag,itol,rtol,atol,&
                        itask,istate,iopt,w2,lrw,iwork,liw,jacobian_tanmap_lsode,mf)
            IF (istate < -1) CALL handle_err(LSODE_ERR,'follow_single',istate)
      END SELECT
      r0 = q(1)
      z0 = q(2)
      IF (PRESENT(tan_map)) THEN
         tan_map(1) = q(3)
         tan_map(2) = q(4)
         tan_map(3) = q(5)
         tan_map(4) = q(6)
      END IF
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE follow_single
