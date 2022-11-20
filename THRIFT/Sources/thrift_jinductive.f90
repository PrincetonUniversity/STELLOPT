!-----------------------------------------------------------------------
!     Subroutine:    thrift_jinductive
!     Authors:       L. van Ham
!     Date:          11/14/2022
!     Description:   This subroutine updates the plasma response to
!                    source currents.  
!-----------------------------------------------------------------------
      SUBROUTINE thrift_jinductive
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE thrift_runtime
      USE thrift_vars
      USE thrift_equil
      USE thrift_profiles_mod
      USE stel_tools
      USE EZspline_obj
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: i, ier
      REAL(rprec) :: rho,s,s11,s12,s21,s22,iota,phip, etapara
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: f1, f2, f3
      TYPE(EZspline1_r8) :: f_spl
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
!     The general idea here is to use equation (22) in form to update
!     THRIFT_JPLASMA.  Where THRIFT_JSOURCE is the J in <J_s.B>.

      ! Allocations
      ALLOCATE(f1(nrho),f2(nrho),f3(nrho))
      f1=0; f2=0; f3=0;

      ! First we calculate helpers
      DO i = 1, nrho
         rho = THRIFT_RHO(i)
         s   = rho*rho
         ier = 0
         CALL get_equil_sus(s,s11,s12,s21,s22,ier)
         CALL EZspline_interp(iota_spl,rho,iota,ier)
         f1(i) = (s11*iota+s12)
         f2(i) = (s21*iota+s22)
      END DO

      ! Now we need to allocate a helper Sup/Sdn and take derivaives
      bcs1=(/ 0, 0/)
      CALL EZspline_init(f_spl,nrho,bcs1,ier)
      f_spl%isHermite   = 0
      CALL EZspline_setup(f_spl,f1/f2,ier,EXACT_DIM=.true.)
      DO i = 1, nrho
         CALL EZspline_derivative(f_spl,1,THRIFT_RHO(i),f3(i),ier)
         CALL EZspline_interp(phip_spl,THRIFT_RHO(i),phip,ier)
         f3(i)=f3(i)*phip
         CALL get_equil_Bav(s,s11,s22,ier)
         f1(i) = s11 ! Bsq <B>
      END DO
      CALL EZspline_free(f_spl,ier)
      f3 = f3*f2*f2
      f1 = f1*2*mu0

      ! Now we subtrac the driven current
      !   THRIFT_JSOURCE is dI/ds  but we need <j.B>
      !    <j.B> = a * rho * <B> * dI/ds / dV/ds [A*T/m^2]
      DO i = 1, nrho
         CALL get_prof_etapara(THRIFT_RHO(i),THRIFT_T(mytimestep),etapara)
         CALL EZspline_interp(vp_spl,THRIFT_RHO(i),s11,ier)
         f3(i) = (f3(i)/mu0 - f1(i)*THRIFT_JSOURCE(i,mytimestep)*s11)*etapara
      END DO

      ! Now we take the derivative d/dPhi = 1/(dPhi/drho) dShelp/drho


      
      ! Deallocate
      DEALLOCATE(f1,f2)

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_jinductive

