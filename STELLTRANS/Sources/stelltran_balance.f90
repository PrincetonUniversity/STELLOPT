!-----------------------------------------------------------------------
!     Subroutine:    stelltran_balance
!     Authors:       J. Mittelstaedt (jmittelstaedt@uchicago.edu)
!     Date:          09/01/2016
!     Description:   This subroutine calculates transport coefficients
!                    from input profiles using transport equations.
!-----------------------------------------------------------------------
      SUBROUTINE stelltran_balance
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE stelltran_input_mod, ONLY: rad_param ! NOTE: This should be temporary. Just to get the rad_param thing
      USE stelltran_runtime
      USE EZspline_obj
      USE EZspline
      USE stelltran_equilutils
      USE stelltran_vars
!-----------------------------------------------------------------------
!     Local Variables
!          zeta           Location of boundary condition
!          z              Array of relevant derivatives
!          u              Solution vector
!          tol            Tolerace on solution and its derivatives (pretty sure is relative, constrained by tol + tol*z(u) where z(u) is the solution)
!          fixpnt         Points to always include in mesh while solving
!          fspace         Work array
!          ispace         Integer work array
!          m              Array of orders of differential equations
!          ipar           Array of parameters influencing how calculations are done
!          iflag          Mode of return of Colnew, diagnoses possible errors
!          ltol           Array specifying which tolerances correspond to which elements of u
!          n_max          Maximum number of subintervals
!          RHS_e          Right hand side of electron power balance equation
!          RHS_i          Right hand side of ion power balance equation
!          VpQe             Convolutoin of V' and Q_e
!          VpQi             Convolutoin of V' and Q_i
!          bcs1           Boundary conditions for spline
!          rho            Normalized radial coordinate, averaged around a flux surface.
!          gradrho        Geometric factor describing the average gradient of rho around flux surface
!          gradrhosq      Geometric factor describing the average gradient of rho around flux surface squared
!          DteDr          Derivative of te wrt rho
!          DtiDr          Derivative of ti wrt rho
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), DIMENSION(prof_length) :: Vps, gradrho, rho_coord
      REAL(rprec), DIMENSION(prof_length) :: VpQe, VpQi, VpGe
      REAL(rprec), DIMENSION(prof_length) :: RHS_powe, RHS_powi, RHS_parte, DteDr, DtiDr, DneDr
      REAL(rprec), DIMENSION(prof_length) :: tau_e, lambda, pow_density
      REAL(rprec) :: irrelevant, mi, s, perf_temp
      REAL(rprec) :: rad_pow(prof_length), tot_ECRH_pow, pt_pow, delta_rho
      INTEGER :: j, ier, dstatus, astatus


      ! COLNEW variables
      INTEGER :: m(3), n_max, ipar(11), ltol(3)
      REAL*8 :: zeta(3), z(3), u(3), tol(3), fixpnt(1)
      REAL*8, ALLOCATABLE :: fspace(:)
      INTEGER, ALLOCATABLE :: ispace(:)
      EXTERNAL :: dummy_guess
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      Er(:,2) = 0 ! Tentative while we figure things out.
      mi = 1860.*me ! tentative, idk where this comes from.
      DO j=1,prof_length
            s = real(j-1)/real(prof_length-1) ! sqrt to go from s to rho
            S_pe(j,2) = 1.D26*EXP(-s**2/1.0D-2)
      END DO
      ni(:,2) = ne(:,2)/zeff(:,2)

      ! set COLNEW parameters for electron power balance
      n_max = 1024
      ALLOCATE(fspace(1:286*n_max))
      ALLOCATE(ispace(1:18*n_max))
      zeta = 0.D0 ! Say boundary conditions are at r=0
      m = 1 ! orders of equations (all 1)
      ipar = 0 ! use mostly defaults
      ipar(3) = 2 ! number of points in initial mesh
      ipar(4) = 3 ! number of tolerances supplied
      ipar(5) = 286*n_max ! dimension of fspace
      ipar(6) = 18*n_max ! dimension of ispace
      ipar(7) = 1 ! printout control. 1= none, 0= selected, -1=full diagnostic
      ltol = (/ 1,2,3/) ! Location of tolerances in solution vector
      tol = (/1.D-3, 1.D-6, 1.D-6/) ! tolerance on the solution

      ! Find needed geometric factors and move to rho grid
      pow_density = 0
      DO j = 1, prof_length
         ! Given S this function returns rho,Vp,<|grad(rho)|>,<|grad(rho)|^2>
         CALL get_equil_rho(Vp(j,1),rho_coord(j),Vps(j),gradrho(j),gradrhosq(j,2),ier)
         CALL eval_prof_spline(prof_length,sqrt(dPdV_erf(:,1)),dPdV_erf(:,2),rho_coord(j),perf_temp,ier)
         pow_density(j) = perf_temp
      END DO

      ! Convert dV/dPhi to dV/drho
      Vp(:,2) = 2*rho_coord(:)*Vps(:)
      Vp(1,2) = Vp(2,2) ! TENTATIVE while we decide what to do about division by zero stuff

      ! Calculate colisional power
      CALL calc_pi_coll(ne(:,2),te(:,2),ti(:,2),zeff(:,2),mi,pi_col(:,2),lambda,tau_e)

      ! Set up arrays of RHS of balance equatoins to spline over
      RHS_powe = Vp(:,2) * (pow_density(:) - pi_col(:,2))
      RHS_parte = Vp(:,2) * S_pe(:,2)
      RHS_powi = Vp(:,2) * pi_col(:,2)

      bcs1 = (/ 0, 0 /)

      ! set up splines
      IF (EZspline_allocated(Pow_e_spl)) CALL EZspline_free(Pow_e_spl,ier)
      CALL EZspline_init(Pow_e_spl,prof_length,bcs1,ier)
      Pow_e_spl%isHermite = 1
      Pow_e_spl%x1 = rho_coord(:)
      CALL EZspline_setup(Pow_e_spl,RHS_powe,ier,EXACT_DIM=.true.)
      IF (EZspline_allocated(Pow_i_spl)) CALL EZspline_free(Pow_i_spl,ier)
      CALL EZspline_init(Pow_i_spl,prof_length,bcs1,ier)
      Pow_i_spl%isHermite = 1
      Pow_e_spl%x1 = rho_coord(:)
      CALL EZspline_setup(Pow_i_spl,RHS_powi,ier,EXACT_DIM=.true.)
      IF (EZspline_allocated(Er_spl)) CALL EZspline_free(Er_spl,ier)
      CALL EZspline_init(Er_spl,prof_length,bcs1,ier)
      Pow_i_spl%isHermite = 1
      Pow_e_spl%x1 = rho_coord(:)
      CALL EZspline_setup(Er_spl,Er(:,2),ier,EXACT_DIM=.true.)
      IF (EZspline_allocated(Part_sple)) CALL EZspline_free(Part_sple,ier)
      CALL EZspline_init(Part_sple,prof_length,bcs1,ier)
      Part_sple%isHermite = 1
      Part_sple%x1 = rho_coord(:)
      CALL EZspline_setup(Part_sple,RHS_parte,ier,EXACT_DIM=.true.)

      ! Should run COLNEW, increasing possible number of subintervals until we get a 'serious' error or we have normal output.
      ier = -1
      DO WHILE (ier /= 1)
            ! Integrate the transport equation to find QeV'
            CALL colnew(3, m, 0.D0, 1.D0, zeta, ipar, ltol, tol, fixpnt, ispace,& 
            fspace, ier, bal_rhs, dbal_rhs, bal_bound_cond, dbal_bound_cond, dummy_guess)
            write(*,'(A,I3.1,A,I5.1,A)'), '****** COLNEW Diagnostic flag is: ', ier, ' nmax is: ', n_max, ' ******'
            IF (ier == -1) THEN! In case it fails, increase max subintervals
                  n_max = 2*n_max
                  ipar(5) = 286*n_max
                  ipar(6) = 18*n_max
                  DEALLOCATE(fspace,STAT=dstatus)
                  DEALLOCATE(ispace,STAT=dstatus)
                  ALLOCATE(fspace(1:286*n_max),STAT=astatus)
                  ALLOCATE(ispace(1:18*n_max),STAT=astatus)
            ELSE IF (ier == 1) THEN
                  CONTINUE
            ELSE ! 'serious' error not fixed by increasing # of subint
                  print*, 'More serious COLNEW error encountered'
                  EXIT
            END IF
      END DO
      ! Apply found solution
      DO j=1,prof_length
            CALL appsln(rho_coord(j),z,fspace,ispace)
            VpGe(j) = z(1)
            VpQe(j) = z(2)
            VpQi(j) = z(3)
      END DO
      Ge(:,2) = VpGe/(gradrho(:)*Vp(:,2))
      Qe(:,2) = VpQe/(gradrho(:)*Vp(:,2))  ! gradrho from geometric factors ignored in Turkin
      Qi(:,2) = VpQi/(gradrho(:)*Vp(:,2))

      ! Take derivative of the spline
      DO j=1,prof_length
            CALL eval_prof_spline(prof_length,rho_coord(:),te(:,2)*ec,rho_coord(j),irrelevant,ier,DteDr(j))
      END DO
      DO j=1,prof_length
            CALL eval_prof_spline(prof_length,rho_coord(:),ti(:,2)*ec,rho_coord(j),irrelevant,ier,DtiDr(j))
      END DO
      DO j=1,prof_length
            CALL eval_prof_spline(prof_length,rho_coord(:),ne(:,2),rho_coord(j),irrelevant,ier,DneDr(j))
      END DO

      De(:,2) =  VpGe(:)/(Vp(:,2)*gradrhosq(:,2)*DneDr(:))
      Xe(:,2) = VpQe(:)/(gradrhosq(:,2)*Vp(:,2)*DteDr(:)*ne(:,2))
      Xi(:,2) = VpQi(:)/(Vp(:,2)*gradrhosq(:,2)*DtiDr(:)*ni(:,2))

      ! Just for testing, write stuff to an output file.
      CALL save_ST_param(8,'./Ge.txt',prof_length,Ge(:,2),Ge(:,1))
      CALL save_ST_param(8,'./Qe.txt',prof_length,Qe(:,2),Qe(:,1))
      CALL save_ST_param(8,'./Qi.txt',prof_length,Qi(:,2),Qi(:,1))
      CALL save_ST_param(8,'./De.txt',prof_length,De(:,2),De(:,1))
      CALL save_ST_param(8,'./Xe.txt',prof_length,Xe(:,2),Xe(:,1))
      CALL save_ST_param(8,'./Xi.txt',prof_length,Xi(:,2),Xi(:,1))
      CALL save_ST_param(8,'./Vp.txt',prof_length,Vp(:,2),Vp(:,1))
      CALL save_ST_param(8,'./te.txt',prof_length,te(:,2),te(:,1))
      CALL save_ST_param(8,'./ne.txt',prof_length,ne(:,2),ne(:,1))
      CALL save_ST_param(8,'./ti.txt',prof_length,ti(:,2),ti(:,1))
      CALL save_ST_param(10,'./S_pe.txt',prof_length,S_pe(:,2),S_pe(:,1))
      CALL save_ST_param(17,'./pow_density.txt',prof_length,pow_density(:),Qe(:,1))
      CALL save_ST_param(10,'./Ptot.txt',prof_length,Ptot(:,2),Ptot(:,1))
      CALL save_ST_param(10,'./dPdV_erf.txt',prof_length,dPdV_erf(:,2),dPdV_erf(:,1))
      CALL save_ST_param(12,'./pi_col.txt',prof_length,pi_col(:,2),pi_col(:,1))
      CALL save_ST_param(12,'./lambda.txt',prof_length,lambda(:),Qe(:,1))
      CALL save_ST_param(11,'./tau_e.txt',prof_length,tau_e(:),Qe(:,1))
      CALL save_ST_param(11,'./DneDr.txt',prof_length,DneDr(:),Qe(:,1))
      CALL save_ST_param(11,'./DteDr.txt',prof_length,DteDr(:),Qe(:,1))
      CALL save_ST_param(11,'./DtiDr.txt',prof_length,DtiDr(:),Qe(:,1))
      CALL save_ST_param(13,'./rad_pow.txt',prof_length,rad_pow(:),Qe(:,1))
      CALL save_ST_param(15,'./rho_coord.txt',prof_length,rho_coord(:),Qe(:,1))
      CALL save_ST_param(13,'./gradrho.txt',prof_length,gradrho(:),Qe(:,1))
      CALL save_ST_param(15,'./gradrhosq.txt',prof_length,gradrhosq(:,2),gradrhosq(:,1))

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE stelltran_balance

      SUBROUTINE dummy_guess ! would be the initial guess for the solution but we are not supplying one.
      END SUBROUTINE dummy_guess

      SUBROUTINE calc_pi_coll(n_e,t_e,t_i,z_eff,m_i,p_i,lambda,tau_e)
            USE stel_kinds, ONLY: rprec
            USE stelltran_vars, ONLY: prof_length, ec, me
            REAL(rprec), DIMENSION(prof_length), INTENT(IN) :: n_e, t_e, t_i, z_eff
            REAL(rprec), INTENT(IN) :: m_i
            REAL(rprec), DIMENSION(prof_length), INTENT(OUT) :: p_i
            REAL(rprec), DIMENSION(prof_length), INTENT(OUT) :: lambda, tau_e ! TENTATIVE change to purely internal later
            INTEGER :: j
            ! NOTE: Formulary eqns are in cgs so need to do converting
            ! Determine appropriate Coulomb Logarithm, NRL Formulary p. 34
            DO j=1,prof_length
                  IF (t_e(j) < 10.*z_eff(j)**2.) THEN
                        lambda(j) = 23. - LOG(SQRT(n_e(j)/1.D6)*z_eff(j)/t_e(j)**(1.5))
                  ELSE
                        lambda(j) = 24. - LOG(SQRT(n_e(j)/1.D6)/t_e(j))
                  END IF
            END DO
            ! Find electron -> ion energy from collisions. NRL Formulary p. 37
            tau_e = 3.44D5*t_e(:)**1.5/(lambda(:)*(n_e(:)/1.D6))
            p_i = 3.0*me*n_e(:)*ec*(t_e(:)-t_i(:))/(m_i*tau_e(:)) ! Should be good to not use cgs here since tau is in sec, ec is for eV -> J conversion
      END SUBROUTINE calc_pi_coll

      SUBROUTINE save_ST_param(str_len,fname_saving,npts,indep,dep)
            INTEGER, INTENT(IN) :: npts, str_len
            CHARACTER(str_len), INTENT(IN) :: fname_saving
            REAL*8, DIMENSION(npts), INTENT(IN) :: indep, dep
            INTEGER :: j
            IF (access(fname_saving,' ') == 0) THEN
                  OPEN(UNIT=12, FILE=fname_saving, ACTION="write", STATUS="old",POSITION="append")
                  do j=1,npts
                        write(12,*) dep(j), indep(j)
                  end do
                  CLOSE(UNIT=12)
            ELSE
                  OPEN(UNIT=12, FILE=fname_saving, ACTION="write", STATUS="new")
                  write(12,'(I2.2)') npts
                  do j=1,npts
                        write(12,*) dep(j), indep(j)
                  end do
                  CLOSE(UNIT=12)
            END IF
      END SUBROUTINE save_ST_param
