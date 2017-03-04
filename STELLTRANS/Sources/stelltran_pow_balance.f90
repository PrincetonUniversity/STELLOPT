!-----------------------------------------------------------------------
!     Subroutine:    stelltran_pow_balance
!     Authors:       J. Mittelstaedt (jmittelstaedt@uchicago.edu)
!     Date:          06/15/2016
!     Description:   This subroutine calculates the electron and ion heat
!                    flux coefficients by integrating the transport equation.
!-----------------------------------------------------------------------
      SUBROUTINE stelltran_pow_balance
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
      REAL(rprec), DIMENSION(prof_length) :: VpQe, VpQi, RHS_e
      REAL(rprec), DIMENSION(prof_length) :: RHS_i, DteDr, DtiDr
      REAL(rprec), DIMENSION(prof_length) :: tau_e, lambda
      REAL(rprec) :: irrelevant, mi
      REAL(rprec) :: rad_pow(prof_length), tot_ECRH_pow, pt_pow, delta_rho
      INTEGER :: j, ier, dstatus, astatus


      ! COLNEW variables
      INTEGER :: isize, fsize, n_max, m(1), ipar(11), ltol(1), iflag
      REAL*8 :: zeta(1), z(1), u(1), tol(1), fixpnt(1)
      REAL*8, ALLOCATABLE :: fspace(:)
      INTEGER, ALLOCATABLE :: ispace(:)
      EXTERNAL :: dummy_guess_powe, dummy_guess_powi
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      Er(:,2) = 0 ! Tentative while we figure things out.
      mi = 2000.*me ! tentative, idk where this comes from.
!       rad_param = 0.05 ! Fraction of total power to be lost as radiation

      ! set COLNEW parameters for electron power balance
      n_max = 25
      isize = 8*n_max
      fsize = 56*n_max
      ALLOCATE(fspace(1:fsize))
      ALLOCATE(ispace(1:isize))
      zeta(1) = 0.D0 ! Say boundary condition is at r=0
      m(1) = 1
      ipar = 0 ! use mostly defaults
      ipar(3) = 2 ! number of points in initial mesh
      ipar(4) = 1 ! number of tolerances supplied
      ipar(5) = fsize ! dimension of fspace
      ipar(6) = isize ! dimension of ispace
      ipar(7) = 1 ! printout control. 1= none, 0= selected, -1=full diagnostic
      ltol(1) = 1
      tol(1) = 1.D-6 ! tolerance on the solution


      ! Find needed geometric factors.
      DO j = 1, prof_length
         ! Given S this function returns rho,Vp,<|grad(rho)|>,<|grad(rho)|^2>
         CALL get_equil_rho(Vp(j,1),rho_coord(j),Vps(j),gradrho(j),gradrhosq(j,2),ier)
      END DO

      ! Convert dV/dPhi to dV/drho
      Vp(:,2) = 2*rho_coord(:)*Vps(:)
      Vp(1,2) = Vp(2,2) ! TENTATIVE while we decide what to do about division by zero stuff

      ! Make an estimate for radiated power
      IF (rad_pow /= 0.0) THEN
            tot_ECRH_pow = 0
            DO j=1,prof_length-1
                  pt_pow = (Vp(j+1,2)*pe_rf(j+1,2) + Vp(j,2)*pe_rf(j,2))/2.
                  delta_rho = (rho_coord(j+1) - rho_coord(j))
                  tot_ECRH_pow = tot_ECRH_pow + delta_rho*pt_pow
            END DO
            ! When integrated over volume, this will give a fraction (rad param) of the total input ECRH power
            rad_pow = rad_param*tot_ECRH_pow * (EXP(rho_coord(:,2))-1.)/((EXP(1.) - 2.)*Vp(:,2))
      ELSE
            rad_pow = 0.0
      END IF

      CALL calc_pi_coll(ne(:,2),te(:,2),ti(:,2),zeff(:,2),mi,pi_col(:,2),lambda,tau_e)

      RHS_e = Vp(:,2) * (pe_rf(:,2) - pi_col(:,2) - rad_pow(:) - Ge(:,2)*Er(:,2))
      RHS_i = Vp(:,2) * (pi_col(:,2) + zeff(:,2)*Gi(:,2)*Er(:,2))

      bcs1 = (/ 0, 0 /)

      ! set up spline for this particular line
      IF (EZspline_allocated(Pow_e_spl)) CALL EZspline_free(Pow_e_spl,ier)
      CALL EZspline_init(Pow_e_spl,prof_length,bcs1,ier)
      Pow_e_spl%isHermite = 1
      Pow_e_spl%x1 = rho_coord(:,2)
      CALL EZspline_setup(Pow_e_spl,RHS_e,ier,EXACT_DIM=.true.)

      IF (EZspline_allocated(Pow_i_spl)) CALL EZspline_free(Pow_i_spl,ier)
      CALL EZspline_init(Pow_i_spl,prof_length,bcs1,ier)
      Pow_i_spl%isHermite = 1
      Pow_e_spl%x1 = rho_coord(:,2)
      CALL EZspline_setup(Pow_i_spl,RHS_i,ier,EXACT_DIM=.true.)

      ! Should run COLNEW, increasing possible number of subintervals until we get a 'serious' error or we have normal output.
      iflag = -1
      DO WHILE (iflag /= 1)
            ! Integrate the transport equation to find QeV'
            CALL colnew(1, m, 0.D0, 1.D0, zeta, ipar, ltol, tol, fixpnt, ispace,& 
            fspace, iflag, pow_e_rhs, dpow_e_rhs, powe_bound_cond, dpowe_bound_cond, dummy_guess_powe)
            write(*,*) 'in electron pow balance'
            write(*,'(A,I3.1,A,I5.1,A)'), '****** COLNEW Diagnostic flag is: ', iflag, ' nmax is: ', n_max, ' ******'
            IF (iflag == -1) THEN! In case it fails, increase max subintervals
                  n_max = n_max + 5
                  isize = 8*n_max
                  fsize = 56*n_max
                  ipar(5) = fsize
                  ipar(6) = isize
                  DEALLOCATE(fspace,STAT=dstatus)
                  DEALLOCATE(ispace,STAT=dstatus)
                  ALLOCATE(fspace(1:fsize),STAT=astatus)
                  ALLOCATE(ispace(1:isize),STAT=astatus)
            ELSE IF (iflag == 1) THEN
                  CONTINUE
            ELSE ! 'serious' error not fixed by increasing # of subint
                  print*, 'More serious COLNEW error encountered'
                  EXIT
            END IF
      END DO
      ! Apply found solution
      DO j=1,prof_length
            CALL appsln(rho_coord(j),z,fspace,ispace)
            VpQe(j) = z(1)
      END DO
      Qe(:,2) = VpQe/(gradrho(:)*Vp(:,2))  ! gradrho from geometric factors ignored in Turkin


      ! Take derivative of the spline and put it in DteDr
      DO j=1,prof_length
            CALL eval_prof_spline(prof_length,rho_coord(:),te(:,2)*ec,rho_coord(j),irrelevant,ier,DteDr(j))
      END DO

      Xe(:,2) = VpQe(:)/(gradrhosq(:,2)*Vp(:,2)*DteDr(:)*ne(:,2))

      ! set COLNEW parameters for ion power balance
      n_max = 25
      isize = 8*n_max
      fsize = 56*n_max
      DEALLOCATE(fspace,STAT=dstatus)
      DEALLOCATE(ispace,STAT=dstatus)
      ALLOCATE(fspace(1:fsize),STAT=astatus)
      ALLOCATE(ispace(1:isize),STAT=astatus)
      zeta(1) = 0.D0 ! Say boundary condition is at r=0
      m(1) = 1
      ipar = 0 ! use mostly defaults
      ipar(3) = 2 ! number of points in initial mesh
      ipar(4) = 1 ! number of tolerances supplied
      ipar(5) = fsize ! dimension of fspace
      ipar(6) = isize ! dimension of ispace
      ipar(7) = 1 ! printout control. 1= none, 0= selected, -1=full diagnostic
      ltol(1) = 1
      tol(1) = 1.D-6 ! tolerance on the solution

      iflag = -1
      DO WHILE (iflag /= 1)
            ! Integrate the transport equation to find QiV'
            CALL colnew(1, m, 0.D0, 1.D0, zeta, ipar, ltol, tol, fixpnt, ispace,& 
            fspace, iflag, pow_i_rhs, dpow_i_rhs, powi_bound_cond, dpowi_bound_cond, dummy_guess_powi)
            write(*,*) 'in ion pow balance'
            write(*,'(A,I3.1,A,I5.1,A)'), '****** COLNEW Diagnostic flag is: ', iflag, ' nmax is: ', n_max, ' ******'
            IF (iflag == -1) THEN! In case it fails, increase max subintervals
                  n_max = n_max + 5 
                  isize = 8*n_max
                  fsize = 56*n_max
                  ipar(5) = fsize
                  ipar(6) = isize
                  DEALLOCATE(fspace,STAT=dstatus)
                  DEALLOCATE(ispace,STAT=dstatus)
                  ALLOCATE(fspace(1:fsize),STAT=astatus)
                  ALLOCATE(ispace(1:isize),STAT=astatus)
            ELSE IF (iflag == 1) THEN
                  CONTINUE
            ELSE ! 'serious' error not fixed by increasing # of subint
                  print*, 'More serious COLNEW error encountered'
                  EXIT
            END IF
      END DO
      ! Apply found solution
      DO j=1,prof_length
            CALL appsln(rho_coord(j),z,fspace,ispace)
            VpQi(j) = z(1)
      END DO
      Qi(:,2) = VpQi/(gradrho(:)*Vp(:,2)) ! gradrho from geometric factors ignored in Turkin

      ! Take derivative of the spline and put it in DtiDr
      DO j=1,prof_length
            CALL eval_prof_spline(prof_length,rho_coord(:),ti(:,2)*ec,rho_coord(j),irrelevant,ier,DtiDr(j))
      END DO

      Xi(:,2) = VpQi(:)/(Vp(:,2)*gradrhosq(:,2)*DtiDr(:)*ni(:,2))

      ! Just for testing, write stuff to an output file.
      IF (access("./Qe.txt",' ') == 0) THEN
            OPEN(UNIT=12, FILE="./Qe.txt", ACTION="write", STATUS="old",POSITION="append")
            do j=1,prof_length
                  write(12,*) Qe(j,:)
            end do
            CLOSE(UNIT=12)
      ELSE
            OPEN(UNIT=12, FILE="./Qe.txt", ACTION="write", STATUS="new")
            write(12,'(I2.2)') prof_length
            do j=1,prof_length
                  write(12,*) Qe(j,:)
            end do
            CLOSE(UNIT=12)
      END IF

      IF (access("./Xe.txt",' ') == 0) THEN
            OPEN(UNIT=12, FILE="./Xe.txt", ACTION="write", STATUS="old",POSITION="append")
            do j=1,prof_length
                  write(12,*) Xe(j,:)
            end do
            CLOSE(UNIT=12)
      ELSE
            OPEN(UNIT=12, FILE="./Xe.txt", ACTION="write", STATUS="new")
            write(12,'(I2.2)') prof_length
            do j=1,prof_length
                  write(12,*) Xe(j,:)
            end do
            CLOSE(UNIT=12)
      END IF

      IF (access("./Qi.txt",' ') == 0) THEN
            OPEN(UNIT=12, FILE="./Qi.txt", ACTION="write", STATUS="old",POSITION="append")
            do j=1,prof_length
                  write(12,*) Qi(j,:)
            end do
            CLOSE(UNIT=12)
      ELSE
            OPEN(UNIT=12, FILE="./Qi.txt", ACTION="write", STATUS="new")
            write(12,'(I2.2)') prof_length
            do j=1,prof_length
                  write(12,*) Qi(j,:)
            end do
            CLOSE(UNIT=12)
      END IF

      IF (access("./Xi.txt",' ') == 0) THEN
            OPEN(UNIT=12, FILE="./Xi.txt", ACTION="write", STATUS="old",POSITION="append")
            do j=1,prof_length
                  write(12,*) Xi(j,:)
            end do
            CLOSE(UNIT=12)
      ELSE
            OPEN(UNIT=12, FILE="./Xi.txt", ACTION="write", STATUS="new")
            write(12,'(I2.2)') prof_length
            do j=1,prof_length
                  write(12,*) Xi(j,:)
            end do
            CLOSE(UNIT=12)
      END IF

      IF (access("./Vp.txt",' ') == 0) THEN
            OPEN(UNIT=12, FILE="./Vp.txt", ACTION="write", STATUS="old",POSITION="append")
            do j=1,prof_length
                  write(12,*) Vp(j,:)
            end do
            CLOSE(UNIT=12)
      ELSE
            OPEN(UNIT=12, FILE="./Vp.txt", ACTION="write", STATUS="new")
            write(12,'(I2.2)') prof_length
            do j=1,prof_length
                  write(12,*) Vp(j,:)
            end do
            CLOSE(UNIT=12)
      END IF

      IF (access("./te.txt",' ') == 0) THEN
            OPEN(UNIT=12, FILE="./te.txt", ACTION="write", STATUS="old",POSITION="append")
            do j=1,prof_length
                  write(12,*) te(j,:)
            end do
            CLOSE(UNIT=12)
      ELSE
            OPEN(UNIT=12, FILE="./te.txt", ACTION="write", STATUS="new")
            write(12,'(I2.2)') prof_length
            do j=1,prof_length
                  write(12,*) te(j,:)
            end do
            CLOSE(UNIT=12)
      END IF

      IF (access("./ne.txt",' ') == 0) THEN
            OPEN(UNIT=12, FILE="./ne.txt", ACTION="write", STATUS="old",POSITION="append")
            do j=1,prof_length
                  write(12,*) ne(j,:)
            end do
            CLOSE(UNIT=12)
      ELSE
            OPEN(UNIT=12, FILE="./ne.txt", ACTION="write", STATUS="new")
            write(12,'(I2.2)') prof_length
            do j=1,prof_length
                  write(12,*) ne(j,:)
            end do
            CLOSE(UNIT=12)
      END IF

      IF (access("./ti.txt",' ') == 0) THEN
            OPEN(UNIT=12, FILE="./ti.txt", ACTION="write", STATUS="old",POSITION="append")
            do j=1,prof_length
                  write(12,*) ti(j,:)
            end do
            CLOSE(UNIT=12)
      ELSE
            OPEN(UNIT=12, FILE="./ti.txt", ACTION="write", STATUS="new")
            write(12,'(I2.2)') prof_length
            do j=1,prof_length
                  write(12,*) ti(j,:)
            end do
            CLOSE(UNIT=12)
      END IF

      IF (access("./pe_rf.txt",' ') == 0) THEN
            OPEN(UNIT=12, FILE="./pe_rf.txt", ACTION="write", STATUS="old",POSITION="append")
            do j=1,prof_length
                  write(12,*) pe_rf(j,:)
            end do
            CLOSE(UNIT=12)
      ELSE
            OPEN(UNIT=12, FILE="./pe_rf.txt", ACTION="write", STATUS="new")
            write(12,'(I2.2)') prof_length
            do j=1,prof_length
                  write(12,*) pe_rf(j,:)
            end do
            CLOSE(UNIT=12)
      END IF

      IF (access("./pi_col.txt",' ') == 0) THEN
            OPEN(UNIT=12, FILE="./pi_col.txt", ACTION="write", STATUS="old",POSITION="append")
            do j=1,prof_length
                  write(12,*) pi_col(j,:)
            end do
            CLOSE(UNIT=12)
      ELSE
            OPEN(UNIT=12, FILE="./pi_col.txt", ACTION="write", STATUS="new")
            write(12,'(I2.2)') prof_length
            do j=1,prof_length
                  write(12,*) pi_col(j,:)
            end do
            CLOSE(UNIT=12)
      END IF

      IF (access("./lambda.txt",' ') == 0) THEN
            OPEN(UNIT=12, FILE="./lambda.txt", ACTION="write", STATUS="old",POSITION="append")
            do j=1,prof_length
                  write(12,*) Vp(j,1), lambda(j)
            end do
            CLOSE(UNIT=12)
      ELSE
            OPEN(UNIT=12, FILE="./lambda.txt", ACTION="write", STATUS="new")
            write(12,'(I2.2)') prof_length
            do j=1,prof_length
                  write(12,*) Vp(j,1), lambda(j)
            end do
            CLOSE(UNIT=12)
      END IF

      IF (access("./tau_e.txt",' ') == 0) THEN
            OPEN(UNIT=12, FILE="./tau_e.txt", ACTION="write", STATUS="old",POSITION="append")
            do j=1,prof_length
                  write(12,*) Vp(j,1), tau_e(j)
            end do
            CLOSE(UNIT=12)
      ELSE
            OPEN(UNIT=12, FILE="./tau_e.txt", ACTION="write", STATUS="new")
            write(12,'(I2.2)') prof_length
            do j=1,prof_length
                  write(12,*) Vp(j,1), tau_e(j)
            end do
            CLOSE(UNIT=12)
      END IF

      IF (access("./DtiDr.txt",' ') == 0) THEN
            OPEN(UNIT=12, FILE="./DtiDr.txt", ACTION="write", STATUS="old",POSITION="append")
            do j=1,prof_length
                  write(12,*) Vp(j,1), DtiDr(j)
            end do
            CLOSE(UNIT=12)
      ELSE
            OPEN(UNIT=12, FILE="./DtiDr.txt", ACTION="write", STATUS="new")
            write(12,'(I2.2)') prof_length
            do j=1,prof_length
                  write(12,*) Vp(j,1), DtiDr(j)
            end do
            CLOSE(UNIT=12)
      END IF

      IF (access("./DteDr.txt",' ') == 0) THEN
            OPEN(UNIT=12, FILE="./DteDr.txt", ACTION="write", STATUS="old",POSITION="append")
            do j=1,prof_length
                  write(12,*) Vp(j,1), DteDr(j)
            end do
            CLOSE(UNIT=12)
      ELSE
            OPEN(UNIT=12, FILE="./DteDr.txt", ACTION="write", STATUS="new")
            write(12,'(I2.2)') prof_length
            do j=1,prof_length
                  write(12,*) Vp(j,1), DteDr(j)
            end do
            CLOSE(UNIT=12)
      END IF

      IF (access("./rad_pow.txt",' ') == 0) THEN
            OPEN(UNIT=12, FILE="./rad_pow.txt", ACTION="write", STATUS="old",POSITION="append")
            do j=1,prof_length
                  write(12,*) Vp(j,1), rad_pow(j)
            end do
            CLOSE(UNIT=12)
      ELSE
            OPEN(UNIT=12, FILE="./rad_pow.txt", ACTION="write", STATUS="new")
            write(12,'(I2.2)') prof_length
            do j=1,prof_length
                  write(12,*) Vp(j,1), rad_pow(j)
            end do
            CLOSE(UNIT=12)
      END IF

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE stelltran_pow_balance

      SUBROUTINE dummy_guess_powe ! would be the initial guess for the solution but we are not supplying one.
      END SUBROUTINE dummy_guess_powe

      SUBROUTINE dummy_guess_powi ! would be the initial guess for the solution but we are not supplying one.
      END SUBROUTINE dummy_guess_powi

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
