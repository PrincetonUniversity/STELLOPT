!-----------------------------------------------------------------------
!     Subroutine:    thrift_jinductive
!     Authors:       L. van Ham
!     Date:          01/23/2023
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
      LOGICAL :: lverbose
      INTEGER :: i, ier
      REAL(rprec) :: rho,s,s11,s12,Bav,Bsqav,iota,phip,etapara
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: f1, f2, f3, f4
      TYPE(EZspline1_r8) :: f_spl
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
!     The general idea here is to use equation (22) in form to update
!     THRIFT_JPLASMA.  Where THRIFT_JSOURCE is the J in <J_s.B>.

      lverbose = .true.
      ! Allocations
      ALLOCATE(f1(nrho),f2(nrho),f3(nrho),f4(nrho))
      f1=0; f2=0; f3=0; f4=0;

      DO i = 1, nrho
         rho = THRIFT_RHO(i)
         s   = rho*rho
         ier = 0
         CALL get_equil_sus(s,s11,s12,Bav,Bsqav,ier) 
         CALL EZspline_interp(iota_spl,rho,iota,ier)   
         CALL EZspline_interp(phip_spl,rho,phip,ier)
         f1(i) = invmu0*phip*(s11*iota+s12) ! I(rho)
         f2(i) = f1(i)/phip ! I/Phi'     
      END DO

      bcs1=(/ 0, 0/)
      CALL EZspline_init(f_spl,nrho,bcs1,ier)
      f_spl%isHermite   = 0
      CALL EZspline_setup(f_spl,f2,ier,EXACT_DIM=.true.)
     
      DO i = 1, nrho
         rho = THRIFT_RHO(i)
         s   = rho*rho
         ier = 0

         CALL EZspline_derivative(f_spl,1,rho,f3(i),ier)
         CALL get_equil_Bav(s,Bav,Bsqav,ier) ! Bav =<B>;Bsqav=<B^2>
         CALL get_prof_pprime(rho,THRIFT_T(mytimestep),s12) !s12 = p'
         CALL get_prof_etapara(rho,THRIFT_T(mytimestep),etapara) !eta = eta
         CALL EZspline_interp(vp_spl,rho,phip,ier) ! phip = V'
         f4(i)= etapara*phip*( mu0*s12*f2(i) + 0.5*Bsqav/rho*f3(i) &
               - THRIFT_JSOURCE(i,mytimestep)*Bav)      
      END DO

      CALL EZspline_free(f_spl,ier)

      bcs1=(/ 0, 0/)
      CALL EZspline_init(f_spl,nrho,bcs1,ier)
      f_spl%isHermite   = 0
      CALL EZspline_setup(f_spl,f4,ier,EXACT_DIM=.true.)
      !f4 = THRIFT_T(mytimestep+1)-THRIFT_T(mytimestep)

      DO i = 1, nrho
         rho = THRIFT_RHO(i)
         s   = rho*rho
         ier = 0
         CALL get_equil_sus(s,s11,s12,Bav,Bsqav,ier) 
         s12 = THRIFT_T(mytimestep+1)-THRIFT_T(mytimestep) !timestep
         CALL EZspline_derivative(f_spl,1,rho,f3(i),ier)
         CALL Ezspline_interp(phip_spl,rho,phip,ier)
         f2(i) = f1(i) + 0.5*invmu0*s12*phip*s11*f3(i)/(rho*eq_phiedge**2) !I next
      END DO

      CALL EZspline_free(f_spl,ier)
      bcs1=(/ 0, 0/)
      CALL EZspline_init(f_spl,nrho,bcs1,ier)
      f_spl%isHermite   = 0
      CALL EZspline_setup(f_spl,f2,ier,EXACT_DIM=.true.)
      DO i = 1, nrho
         rho = THRIFT_RHO(i)
         CALL EZspline_derivative(f_spl,1,rho,f1(i),ier)
         CALL Ezspline_interp(vp_spl,rho,phip,ier)
         THRIFT_JPLASMA(i,mytimestep) = 2*pi*eq_Rmajor*f1(i)/phip
      END DO     
      CALL EZspline_free(f_spl,ier)

      ! Deallocate
      DEALLOCATE(f1,f2,f3,f4)

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_jinductive

