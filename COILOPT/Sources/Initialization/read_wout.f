      SUBROUTINE read_wout_coilopt
      USE read_wout_mod, nfp_w => nfp, ns_w => ns, mpol_w => mpol,
     1   ntor_w => ntor, mnmax_w => mnmax, xm_w => xm, xn_w => xn
      USE boundary
      USE mpi_params                                         !mpi stuff
      USE Vname, ONLY: extension
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: istat
!-----------------------------------------------
      CALL read_wout_file(extension, istat)
      IF (istat .ne. 0) THEN
         IF( myid .eq. master ) THEN                         !mpi stuff
           PRINT *,' Error reading wout file: ierr = ', istat
         END IF     !IF( myid .eq. master )                  !mpi stuff
         STOP
      END IF

!     Load values from read_wout_mod module into boundary arrays

      nfp = nfp_w
      ns  = ns_w
      mpol= mpol_w
      ntor= ntor_w
      mnmax = mnmax_w

!     ALLOCATE ARRAYS FIRST TIME THRU

      ALLOCATE (rmnc_b(mnmax), zmns_b(mnmax), rmnc_a(mnmax),
     1          zmns_a(mnmax), xm_b(mnmax), xn_b(mnmax),
     2          lmns_b(mnmax), stat=istat)
      IF (myid.EQ.master .AND. istat.NE.0)
     1  STOP 'allocation error in read_wout'

!     ONLY SAVE BOUNDARY COEFFICIENTS

      xm_b = zero
      xn_b = zero
      rmnc_b = zero
      zmns_b = zero
      rmnc_a = zero
      zmns_a = zero
      lmns_b = zero

      xm_b(:) = xm_w(:)
      xn_b(:) = xn_w(:)
      rmnc_b(:) = rmnc(:,ns)
      zmns_b(:) = zmns(:,ns)
      rmnc_a(:) = rmnc(:,1)
      zmns_a(:) = zmns(:,1)
      lmns_b(:) = 1.5_dp*lmns(:,ns) - 0.5_dp*lmns(:,ns-1)
      iota_b = 1.5_dp*iotas(ns) - 0.5_dp*iotas(ns-1)
      phip_b = 1.5_dp*phip(ns) - 0.5_dp*phip(ns-1)

      IF (myid .EQ. master) THEN
        WRITE(6,'(A)')   '----- Equilibrium Information -----'
        WRITE(6,'(3X,A,F10.5,A,F10.5,A)') 'R       = [',rmin_surf,
     1          ',',rmax_surf,'] [m]'
        WRITE(6,'(3X,A,F10.5,A)')       'Z       =',zmax_surf, '[m]'
        WRITE(6,'(3X,A,F10.5)')         'Beta    =',betatot
        IF (ABS(Itor) >= 1.0E6) THEN
           WRITE(6,'(3X,A,F10.5,A)')       'Current =',
     1            Itor/1.0E+06, '[MA]'
        ELSE IF (ABS(Itor) >= 1.0E3) THEN
           WRITE(6,'(3X,A,F10.5,A)')       'Current =',
     1            Itor/1.0E+03, '[kA]'
        ELSE
           WRITE(6,'(3X,A,F10.5,A)')   'Current =',Itor, '[A]'
        END IF
        WRITE(6,'(3X,A,F10.5,A)')  'Flux    =',phi(ns), '[Wb]'
        WRITE(6,'(3X,A,F10.5)')   'Iota(edge)=',iota_b
        WRITE(6,'(3X,A,F10.5)')  'Phip(edge)=',phip_b
        WRITE(6,'(3X,A,E14.5)')'Ipol     =',twopi*bvco(ns)/mu0/nfp
        WRITE(6,'(3X,A,F7.2)')   'VMEC v.',version_
        CALL FLUSH(6)
        !PRINT 100, iota_b, phip_b, twopi*bvco(ns)/mu0/nfp
      END IF
!  100 FORMAT (" iota = ",f7.3,", phip = ",f7.3,", ipol = ",1pe14.5)

      nu = 2*mpol_w + 6
      nv = 4*ntor_w

!      nu = 24
!      nv = 24

      nuv = nu * nv

!     DEALLOCATE memory

 1000 CALL read_wout_deallocate

      END SUBROUTINE read_wout_coilopt
