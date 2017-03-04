!-----------------------------------------------------------------------
!     Subroutine:    stelltran_part_balance
!     Authors:       J. Mittelstaedt (jmittelstaedt@uchicago.edu)
!     Date:          06/21/2016
!     Description:   This subroutine calculates particle flux coefficient
!                    by integrating the transport equation
!-----------------------------------------------------------------------
      SUBROUTINE stelltran_part_balance
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
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
!          tol            Tolerace on solution and its derivatives
!          fixpnt         Points to always include in mesh while solving
!          fspace         Work array
!          m              Array of orders of differential equations
!          ipar           Array of parameters influencing how calculations are done
!          iflag          Mode of return of Colnew, diagnoses possible errors
!          ltol           Array specifying which tolerances correspond to which elements of u
!          ispace         Integer work array
!          RHS            Right hand side of particle balance equation
!          VpGe           Convolutoin of V' and G_e
!          VpGi           Convolutoin of V' and G_i
!          bcs1           Boundary conditions for spline
!          n_max          Maximum number of subintervals
!          rho            Normalized radial coordinate, averaged around a flux surface.
!          gradrho        Geometric factor describing the average gradient of rho around flux surface
!          gradrhosq      Geometric factor describing the average gradient of rho around flux surface squared
!          DneDr          Derivative of ne wrt rho
!          DniDr          Derivative of ni wrt rho
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), DIMENSION(prof_length) :: Vps, VpGe, VpGi, RHS_e, DneDr, DniDr, RHS_i
      REAL(rprec) :: s, irrelevant
      REAL*8 :: zeta(1), z(1), u(1), tol(1), fixpnt(1)
      REAL*8, ALLOCATABLE :: fspace(:)
      INTEGER ::  m(1), ipar(11), ltol(1), iflag
      INTEGER, ALLOCATABLE :: ispace(:)
      INTEGER :: j, ier, n_max, isize, fsize, dstatus, astatus
      EXTERNAL :: dummy_guess_parte, dummy_guess_parti
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Make dummy S_p values while we're still making it for real.
      DO j=1,prof_length
            s = real(j-1)/real(prof_length-1) ! sqrt to go from s to rho
            S_pe(j,2) = 1.D5*EXP(-s**2/1.0D-1)
      END DO
      S_pi(:,2) = S_pe(:,2)/zeff(:,2) ! Assuming we maintain quasi-neutrality for now.

      ni(:,2) = ne(:,2)/zeff(:,2) ! Also tentative, assuming somewhat local neutrality

      n_max = 35
      isize = 8*n_max
      fsize = 56*n_max
      ALLOCATE(fspace(1:fsize))
      ALLOCATE(ispace(1:isize))
      zeta(1) = 0.D0 ! Say boundary condition is at r=0 NEED TO THINK ABOUT THIS NOW *****************
      m(1) = 1
      ipar = 0 ! use mostly defaults
      ipar(3) = 2 ! number of points in initial mesh
      ipar(4) = 1 ! number of tolerances supplied
      ipar(5) = 56*n_max ! dimension of fspace
      ipar(6) = 8*n_max ! dimension of ispace
      ipar(7) = 1 ! printout control. 1= none, 0= selected, -1=full diagnostic
      ltol(1) = 1
      tol(1) = 1.D-3 ! tolerance on the solution

      ! Find needed geometric factors.
      DO j = 1, prof_length
         ! Given S this function returns rho,Vp,<|grad(rho)|>,<|grad(rho)|^2>
         CALL get_equil_rho(Vp(j,1),rho_coord(j,2),Vps(j),gradrho(j,2),gradrhosq(j,2),ier)
      END DO

      ! Convert dV/dPhi to dV/drho
      Vp(:,2) = 2*rho_coord(:,2)*Vps(:)
      Vp(1,2) = Vp(2,2) ! TENTATIVE while we decide what to do about division by zero stuff

      RHS_e = Vp(:,2) * S_pe(:,2)
      RHS_i = Vp(:,2) * S_pi(:,2)

      bcs1 = (/ 0, 0 /)

      ! set up spline for this particular line
      IF (EZspline_allocated(Part_sple)) CALL EZspline_free(Part_sple,ier)
      CALL EZspline_init(Part_sple,prof_length,bcs1,ier)
      Part_sple%isHermite = 1
      Part_sple%x1 = rho_coord(:,2)
      CALL EZspline_setup(Part_sple,RHS_e,ier,EXACT_DIM=.true.)

      iflag = -1
      DO WHILE (iflag /= 1)
            ! Integrate the transport equation to find GeV'
            CALL colnew(1, m, 0.D0, 1.D0, zeta, ipar, ltol, tol, fixpnt, ispace,& 
            fspace, iflag, parte_rhs, dparte_rhs, parte_bound_cond, dparte_bound_cond, dummy_guess_parte)
            write(*,*) 'in electron particle balance'
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
            CALL appsln(rho_coord(j,2),z,fspace,ispace)
            VpGe(j) = z(1)
      END DO
      Ge(:,2) = VpGe(:)/(Vp(:,2)*gradrho(:,2))

      ! Take derivative of the spline and put it in DneDr
      DO j=1,prof_length
            CALL eval_prof_spline(prof_length,rho_coord(:,2),ne(:,2),rho_coord(j,2),irrelevant,ier,DneDr(j))
      END DO

      De(:,2) =  VpGe(:)/(Vp(:,2)*gradrhosq(:,2)*DneDr(:))

      n_max = 35
      isize = 8*n_max
      fsize = 56*n_max
      DEALLOCATE(fspace)
      DEALLOCATE(ispace)
      ALLOCATE(fspace(1:fsize))
      ALLOCATE(ispace(1:isize))
      zeta(1) = 0.D0 ! Say boundary condition is at r=0 NEED TO THINK ABOUT THIS NOW *****************
      m(1) = 1
      ipar = 0 ! use mostly defaults
      ipar(3) = 2 ! number of points in initial mesh
      ipar(4) = 1 ! number of tolerances supplied
      ipar(5) = 56*n_max ! dimension of fspace
      ipar(6) = 8*n_max ! dimension of ispace
      ipar(7) = 1 ! printout control. 1= none, 0= selected, -1=full diagnostic
      ltol(1) = 1
      tol(1) = 1.D-3 ! tolerance on the solution

      ! set up spline for this particular line
      IF (EZspline_allocated(Part_spli)) CALL EZspline_free(Part_spli,ier)
      CALL EZspline_init(Part_spli,prof_length,bcs1,ier)
      Part_spli%isHermite = 1
      Part_spli%x1 = rho_coord(:,2)
      CALL EZspline_setup(Part_spli,RHS_i,ier,EXACT_DIM=.true.)

      iflag = -1
      DO WHILE (iflag /= 1)
            ! Integrate the transport equation to find GiV'
            CALL colnew(1, m, 0.D0, 1.D0, zeta, ipar, ltol, tol, fixpnt, ispace,& 
            fspace, iflag, parti_rhs, dparti_rhs, parti_bound_cond, dparti_bound_cond, dummy_guess_parti)
            write(*,*) 'in ion particle balance'
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
            CALL appsln(rho_coord(j,2),z,fspace,ispace)
            VpGi(j) = z(1)
      END DO
      Gi(:,2) = VpGi/(gradrho(:,2)*Vp(:,2))

      ! Take derivative of the spline and put it in DniDr
      DO j=1,prof_length
            CALL eval_prof_spline(prof_length,rho_coord(:,2),ni(:,2),rho_coord(j,2),irrelevant,ier,DniDr(j))
      END DO

      Di(:,2) = VpGi(:)/(Vp(:,2)*gradrhosq(:,2)*DniDr(:))

      ! Just for testing, write the result to an output file.
      IF (access("./Ge.txt",' ') == 0) THEN
            OPEN(UNIT=12, FILE="./Ge.txt", ACTION="write", STATUS="old",POSITION="append")
            do j=1,prof_length
                  write(12,*) Ge(j,:)
            end do
            CLOSE(UNIT=12)
      ELSE
            OPEN(UNIT=12, FILE="./Ge.txt", ACTION="write", STATUS="new")
            write(12,'(I2.2)') prof_length
            do j=1,prof_length
                  write(12,*) Ge(j,:)
            end do
            CLOSE(UNIT=12)
      END IF

      IF (access("./De.txt",' ') == 0) THEN
            OPEN(UNIT=12, FILE="./De.txt", ACTION="write", STATUS="old",POSITION="append")
            do j=1,prof_length
                  write(12,*) De(j,:)
            end do
            CLOSE(UNIT=12)
      ELSE
            OPEN(UNIT=12, FILE="./De.txt", ACTION="write", STATUS="new")
            write(12,'(I2.2)') prof_length
            do j=1,prof_length
                  write(12,*) De(j,:)
            end do
            CLOSE(UNIT=12)
      END IF

      IF (access("./Gi.txt",' ') == 0) THEN
            OPEN(UNIT=12, FILE="./Gi.txt", ACTION="write", STATUS="old",POSITION="append")
            do j=1,prof_length
                  write(12,*) Gi(j,:)
            end do
            CLOSE(UNIT=12)
      ELSE
            OPEN(UNIT=12, FILE="./Gi.txt", ACTION="write", STATUS="new")
            write(12,'(I2.2)') prof_length
            do j=1,prof_length
                  write(12,*) Gi(j,:)
            end do
            CLOSE(UNIT=12)
      END IF

      IF (access("./Di.txt",' ') == 0) THEN
            OPEN(UNIT=12, FILE="./Di.txt", ACTION="write", STATUS="old",POSITION="append")
            do j=1,prof_length
                  write(12,*) Di(j,:)
            end do
            CLOSE(UNIT=12)
      ELSE
            OPEN(UNIT=12, FILE="./Di.txt", ACTION="write", STATUS="new")
            write(12,'(I2.2)') prof_length
            do j=1,prof_length
                  write(12,*) Di(j,:)
            end do
            CLOSE(UNIT=12)
      END IF

      IF (access("./Spe.txt",' ') == 0) THEN
            OPEN(UNIT=12, FILE="./Spe.txt", ACTION="write", STATUS="old",POSITION="append")
            do j=1,prof_length
                  write(12,*) S_pe(j,:)
            end do
            CLOSE(UNIT=12)
      ELSE
            OPEN(UNIT=12, FILE="./Spe.txt", ACTION="write", STATUS="new")
            write(12,'(I2.2)') prof_length
            do j=1,prof_length
                  write(12,*) S_pe(j,:)
            end do
            CLOSE(UNIT=12)
      END IF

      IF (access("./Spi.txt",' ') == 0) THEN
            OPEN(UNIT=12, FILE="./Spi.txt", ACTION="write", STATUS="old",POSITION="append")
            do j=1,prof_length
                  write(12,*) S_pi(j,:)
            end do
            CLOSE(UNIT=12)
      ELSE
            OPEN(UNIT=12, FILE="./Spi.txt", ACTION="write", STATUS="new")
            write(12,'(I2.2)') prof_length
            do j=1,prof_length
                  write(12,*) S_pi(j,:)
            end do
            CLOSE(UNIT=12)
      END IF

      IF (access("./DneDr.txt",' ') == 0) THEN
            OPEN(UNIT=12, FILE="./DneDr.txt", ACTION="write", STATUS="old",POSITION="append")
            do j=1,prof_length
                  write(12,*) Vp(j,1), DneDr(j)
            end do
            CLOSE(UNIT=12)
      ELSE
            OPEN(UNIT=12, FILE="./DneDr.txt", ACTION="write", STATUS="new")
            write(12,'(I2.2)') prof_length
            do j=1,prof_length
                  write(12,*) Vp(j,1), DneDr(j)
            end do
            CLOSE(UNIT=12)
      END IF

      IF (access("./DniDr.txt",' ') == 0) THEN
            OPEN(UNIT=12, FILE="./DniDr.txt", ACTION="write", STATUS="old",POSITION="append")
            do j=1,prof_length
                  write(12,*) Vp(j,1), DniDr(j)
            end do
            CLOSE(UNIT=12)
      ELSE
            OPEN(UNIT=12, FILE="./DniDr.txt", ACTION="write", STATUS="new")
            write(12,'(I2.2)') prof_length
            do j=1,prof_length
                  write(12,*) Vp(j,1), DniDr(j)
            end do
            CLOSE(UNIT=12)
      END IF

      IF (access("./rho_coord.txt",' ') == 0) THEN
            OPEN(UNIT=12, FILE="./rho_coord.txt", ACTION="write", STATUS="old",POSITION="append")
            do j=1,prof_length
                  write(12,*) rho_coord(j,:)
            end do
            CLOSE(UNIT=12)
      ELSE
            OPEN(UNIT=12, FILE="./rho_coord.txt", ACTION="write", STATUS="new")
            write(12,'(I2.2)') prof_length
            do j=1,prof_length
                  write(12,*) rho_coord(j,:)
            end do
            CLOSE(UNIT=12)
      END IF

      IF (access("./gradrho.txt",' ') == 0) THEN
            OPEN(UNIT=12, FILE="./gradrho.txt", ACTION="write", STATUS="old",POSITION="append")
            do j=1,prof_length
                  write(12,*) gradrho(j,:)
            end do
            CLOSE(UNIT=12)
      ELSE
            OPEN(UNIT=12, FILE="./gradrho.txt", ACTION="write", STATUS="new")
            write(12,'(I2.2)') prof_length
            do j=1,prof_length
                  write(12,*) gradrho(j,:)
            end do
            CLOSE(UNIT=12)
      END IF

      IF (access("./gradrhosq.txt",' ') == 0) THEN
            OPEN(UNIT=12, FILE="./gradrhosq.txt", ACTION="write", STATUS="old",POSITION="append")
            do j=1,prof_length
                  write(12,*) gradrhosq(j,:)
            end do
            CLOSE(UNIT=12)
      ELSE
            OPEN(UNIT=12, FILE="./gradrhosq.txt", ACTION="write", STATUS="new")
            write(12,'(I2.2)') prof_length
            do j=1,prof_length
                  write(12,*) gradrhosq(j,:)
            end do
            CLOSE(UNIT=12)
      END IF
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE stelltran_part_balance

      SUBROUTINE dummy_guess_parte ! would be the initial guess for the solution but we are not supplying one.
      END SUBROUTINE dummy_guess_parte

      SUBROUTINE dummy_guess_parti ! would be the initial guess for the solution but we are not supplying one.
      END SUBROUTINE dummy_guess_parti
