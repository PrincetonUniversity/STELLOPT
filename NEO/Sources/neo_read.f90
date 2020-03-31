
SUBROUTINE neo_read
! Read in Boozer Coordinate and Magnetic Data
!***********************************************************************
! Modules
!***********************************************************************
  USE neo_input
  USE neo_units
  USE neo_control
  USE neo_work
  USE neo_exchange
  use safe_open_mod                                   ! SPH
!***********************************************************************
! Local definitions
!***********************************************************************
  IMPLICIT NONE

  INTEGER :: i,j,j_m,j_n,istat
  INTEGER :: m,n,num_m,num_n,m_found,n_found
  INTEGER :: i_alloc
  INTEGER :: id1,id2,id3,id4,id5,id6,id7
  INTEGER :: k, ns_all                                      ! LPK (used to pack arrays)
  CHARACTER(5) :: dummy
  CHARACTER(45) :: cdum
!***********************************************************************
! Open input-unit and read first quantities
!***********************************************************************
  lasym = .false.
  IF( inp_swi .eq. 0 ) THEN
     CALL READ_BOOZ_IN

  ELSE IF (inp_swi .EQ. 1) THEN            !PPPL-style Bmns file
     call safe_open(r_u1, istat, trim(in_file)//"."//extension, 'old',    &
         'formatted')
     IF (istat .ne. 0) THEN
        write (w_us, *) 'IN_FILE.' // trim(extension) //' NOT FOUND IN NEO_READ'
        stop
     END IF
     READ (r_u1,*) dummy
     READ (r_u1,*) m0b,n0b,ns,nfp,flux
     m_max = m0b+1
     n_max = 2*n0b+1
     mnmax = m_max*n_max

     ns_all = ns                                         ! LPK
     if (no_fluxs>0 .and. no_fluxs<=ns ) ns = no_fluxs    ! LPK/SPH

     if (max_m_mode .le. 0) max_m_mode = m0b             ! SPH
     if (max_n_mode .le. 0) max_n_mode = n0b * nfp       ! SPH

! **********************************************************************
! Allocate storage arrays
! **********************************************************************
     ALLOCATE(ixm(mnmax), ixn(mnmax), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

     ALLOCATE(pixm(mnmax), pixn(mnmax), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays pointers failed!'

     ALLOCATE(i_m(m_max), i_n(n_max), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

     ALLOCATE(es(ns), iota(ns), curr_pol(ns), curr_tor(ns),               &
          pprime(ns), sqrtg00(ns), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for real arrays failed!'

     ALLOCATE(rmnc(ns,mnmax), zmns(ns,mnmax), lmns(ns,mnmax),             &
          bmnc(ns,mnmax),                                                 &
          stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for fourier arrays (1) failed!'
!***********************************************************************
! Read input arrays
! Only store REQUIRED arrays (ns), but need to scan all of them          ! LPK
!***********************************************************************
     DO k = 1, ns_all                                                    ! LPK
        IF (allocated (fluxs_arr)) THEN                                  ! SPH
           i = 0
           DO j = 1, no_fluxs
              IF (fluxs_arr(j) .eq. k) THEN
                  i = j
                  exit
              END IF
           END DO
        ELSE
           i = k
        END IF                                                           ! SPH END IF

        IF (i .gt. 0) THEN
           READ(r_u1,*) dummy
           READ(r_u1,*) es(i),iota(i),curr_pol(i),curr_tor(i),               &
                pprime(i),sqrtg00(i)
           READ(r_u1,*) dummy
           DO j=1,mnmax
              READ(r_u1,*) ixm(j),ixn(j),                                    &
                   rmnc(i,j),zmns(i,j),lmns(i,j),                            &
                   bmnc(i,j)
           END DO

        ELSE
           DO j=1,mnmax+3                             ! LPK
              READ(r_u1,*) dummy                      ! LPK
           END DO                                     ! LPK

        END IF

     END DO                                                             ! LPK NS_ALL LOOP

  ELSE IF (inp_swi .EQ. 2) THEN        !ORNL-style Boozer file
     call safe_open(r_u1, istat, trim(in_file), 'old', 'formatted')
     IF (istat .ne. 0) THEN
        write (w_us, *) trim(in_file) // ' NOT FOUND IN NEO_READ'
        stop
     END IF

     READ(r_u1,'(a)') cdum
     READ(r_u1,'(a)') cdum
     READ(r_u1,'(a45,10i5)') cdum,ns,id1,id2,id3,id4,m0b,n0b,id5,nfp,mnmax
     READ(r_u1,'(a44,10i5)') cdum,id1,id2,id3,id4,id5,id6,id7

     ns_all = ns                                         ! LPK
     if (no_fluxs>0 .and. no_fluxs<=ns ) ns = no_fluxs    ! LPK/SPH
     if (max_m_mode .le. 0) max_m_mode = m0b
     if (max_n_mode .le. 0) max_n_mode = n0b * nfp

     DO i=1,4
        READ(r_u1,'(a)')cdum
     END DO
     m_max = m0b+1
     n_max = 2*n0b+1
!     mnmax = m_max*n_max - n0b
! **********************************************************************
! Allocate storage arrays
! **********************************************************************
     ALLOCATE(ixm(mnmax), ixn(mnmax), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

     ALLOCATE(pixm(mnmax), pixn(mnmax), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays pointers failed!'

     ALLOCATE(i_m(m_max), i_n(n_max), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

     ALLOCATE(es(ns), iota(ns), curr_pol(ns), curr_tor(ns),               &
          pprime(ns), sqrtg00(ns), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for real arrays failed!'

     ALLOCATE(rmnc(ns,mnmax), zmnc(ns,mnmax), lmnc(ns,mnmax),             &
          bmnc(ns,mnmax),                                                 &
          stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for fourier arrays (1) failed!'
!***********************************************************************
! Read input arrays
!***********************************************************************
     DO  k = 1,ns_all
        IF (allocated (fluxs_arr)) THEN                                  ! SPH
           i = 0
           DO j = 1, no_fluxs
              IF (fluxs_arr(j) .eq. k) THEN
                  i = j
                  exit
              END IF
           END DO
        ELSE
           i = k
        END IF                                                           ! SPH END IF

        IF (i .gt. 0) THEN
           READ(r_u1,'(a)') cdum
           READ(r_u1,'(5e12.4)') es(i),iota(i),curr_pol(i),curr_tor(i),flux
           pprime(i) = 0; sqrtg00(i) = 0
           READ(r_u1,'(a)') cdum
           READ(r_u1,"(2i5,1p,4e16.8)") (ixm(j),ixn(j),rmnc(i,j),zmns(i,j),  &
              lmns(i,j),bmnc(i,j),j=1,mnmax)
           READ(r_u1,'(a)') cdum
        ELSE
           DO j = 1, mnmax+4
              READ(r_u1,'(a)') dummy                      ! SPH-EAT LINES ON UNWANTED SURFACES
           END DO
        END IF
     END DO
!   do i=1,ns                              !DAS test
!    curr_pol(i) = curr_pol(i)*TWOPI/nfp
!    curr_tor(i) = curr_tor(i)*TWOPI
!   end do
  ELSE
     WRITE (w_us,*) 'FATAL: There is yet no other input type defined'
     STOP
  END IF
!
! Filling of i_m and i_n
! and pointers pixm from ixm to i_m, and pixn from ixn to i_n
  DO j = 1,mnmax
     m = ixm(j)
     n = ixn(j)
     IF (j .EQ. 1) THEN
        num_m = 1
        i_m(num_m) = m
        pixm(j) = num_m
        num_n = 1
        i_n(num_n) = n
        pixn(j) = num_n
     ELSE
        m_found = 0
        DO j_m = 1, num_m
           IF (m .EQ. i_m(j_m)) THEN
              pixm(j) = j_m
              m_found = 1
           END IF
        END DO
        IF (m_found .EQ. 0) THEN
           num_m = num_m + 1
           i_m(num_m) = m
           pixm(j) = num_m
        END IF
        n_found = 0
        DO j_n = 1, num_n
           IF (n .EQ. i_n(j_n)) THEN
              pixn(j) = j_n
              n_found = 1
           END IF
        END DO
        IF (n_found .EQ. 0) THEN
           num_n = num_n + 1
           i_n(num_n) = n
           pixn(j) = num_n
        END IF
     END IF
  END DO
  IF (lab_swi .EQ. 1) THEN
!
! ATTENTION: Switch n TO -n
!            Toroidal mode numbers have to multiplied by number of field periods
!            Change iota to iota*nfp (PRINCETON)
!            it is named iota but actually it is iotabar
     ixn = - ixn * nfp
     i_n = - i_n * nfp
     max_n_mode = max_n_mode * nfp
     iota = iota*nfp
  ELSE IF  (lab_swi .EQ. 2) THEN         !ORNL Boozer file
     DO i=1,ns
        DO j=1,mnmax
           lmnc(i,j) = -nfp*lmns(i,j)/TWOPI
        END DO
     END DO
  ELSE IF (lab_swi .NE. 0) THEN
     WRITE(w_us,*) 'LAB_SWI = ', lab_swi,' NOT HANDLED YET!'
     STOP
  END IF
!
! ATTENTION THIS IS JUST FOR TESTING
!
! For scaling of B change the following three with the same factor
! eps_eff should then stay unchanged if the reference for B and R is the same!
!
! bmnc = bmnc * 2.0_dp
! curr_pol = curr_pol * 2.0_dp
! curr_tor = curr_tor * 2.0_dp
!
! For scaling of R change the following four with the same factor
! eps_eff should then stay unchanged if the reference for B and R is the same!
!
! rmnc = rmnc * 2.0_dp
! zmnc = zmnc * 2.0_dp
! curr_pol = curr_pol * 2.0_dp
! curr_tor = curr_tor * 2.0_dp
!
  CLOSE (unit=r_u1)
! **********************************************************************
! Write optional output for Plotting
! **********************************************************************
  IF (write_output_files .NE. 0) THEN
     IF (write_progress .NE. 0) WRITE (w_us,*) 'write dimension.dat'
     OPEN(unit=w_u1,file='dimension.dat',status='replace',form='formatted')
     WRITE (w_u1,*) ns
     WRITE (w_u1,*) mnmax
     WRITE (w_u1,*) nfp
     WRITE (w_u1,*) theta_n
     WRITE (w_u1,*) phi_n
     WRITE (w_u1,*) s_ind_in
     CLOSE(unit=w_u1)
  ENDIF
  IF (write_output_files .NE. 0) THEN
     IF (write_progress .NE. 0) WRITE (w_us,*) 'write es_arr.dat'
     OPEN(unit=w_u1,file='es_arr.dat',status='replace',form='formatted')
     DO j=1,ns
        WRITE(w_u1,format220) es(j),iota(j),curr_pol(j),curr_tor(j),       &
             pprime(j),sqrtg00(j)
     ENDDO
     CLOSE(unit=w_u1)

     IF (write_progress .NE. 0) WRITE (w_us,*) 'write mn_arr.dat'
     OPEN(unit=w_u1,file='mn_arr.dat',status='replace',form='formatted')
     DO j = 1,mnmax
        WRITE(w_u1,*)ixm(j),ixn(j)
     ENDDO
     CLOSE(unit=w_u1)

     IF (write_progress .NE. 0) WRITE (w_us,*) 'write rmnc_arr.dat'
     OPEN(unit=w_u1,file='rmnc_arr.dat',status='replace',form='formatted')
     DO i=1,ns
        DO j=1,mnmax
           WRITE(w_u1,*) rmnc(i,j)
        END DO
     END DO
     CLOSE(unit=w_u1)

     IF (write_progress .NE. 0) WRITE (w_us,*) 'write zmns_arr.dat'
     OPEN(unit=w_u1,file='zmns_arr.dat',status='replace',form='formatted')
     DO i=1,ns
        DO j=1,mnmax
           WRITE(w_u1,*) zmns(i,j)
        END DO
     END DO
     CLOSE(unit=w_u1)

     IF (write_progress .NE. 0) WRITE (w_us,*) 'write lmns_arr.dat'
     OPEN(unit=w_u1,file='lmns_arr.dat',status='replace',form='formatted')
     DO i=1,ns
        DO j=1,mnmax
           WRITE(w_u1,*) lmns(i,j)
        END DO
     END DO
     CLOSE(unit=w_u1)

     IF (write_progress .NE. 0) WRITE (w_us,*) 'write bmnc_arr.dat'
     OPEN(unit=w_u1,file='bmnc_arr.dat',status='replace',form='formatted')
     DO i=1,ns
        DO j=1,mnmax
           WRITE(w_u1,*) bmnc(i,j)
        END DO
     END DO
     CLOSE(unit=w_u1)
  ENDIF
! **********************************************************************
  RETURN
END SUBROUTINE neo_read
