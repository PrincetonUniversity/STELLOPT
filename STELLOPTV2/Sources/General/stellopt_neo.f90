!-----------------------------------------------------------------------
!     Subroutine:    stellopt_neo
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          06/06/2012
!     Description:   This subroutine calculates neoclassical transport
!                    via the NEO code.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_neo(lscreen,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime, ONLY:  proc_string, bigno
      USE stellopt_vars, ONLY: equil_type
      USE equil_vals, ONLY: eff_ripple
      USE stellopt_targets, ONLY: sigma_neo
      USE mpi_params
      ! NEO LIBRARIES
!DEC$ IF DEFINED (NEO_OPT)
      USE neo_precision
      USE neo_units
      USE neo_parameters
      USE neo_control
      USE neo_input
      USE neo_work
      USE neo_exchange
      USE neo_output
      USE sizey_bo
!DEC$ ENDIF
      USE safe_open_mod                          ! SPH
      ! BOOZER_XFORM
      USE read_boozer_mod, dp1 => dp
      
!-----------------------------------------------------------------------
!     Subroutine Parameters
!        iflag         Error flag
!----------------------------------------------------------------------
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      INCLUDE 'mpif.h'
!DEC$ ENDIF
      LOGICAL, INTENT(in)    :: lscreen
      INTEGER, INTENT(inout) :: iflag
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      INTEGER       ::                     i, j, k, i_surf, mn0, npsi, &
                                           istat, imn, m, n, num_m, num_n, &
                                           m_found, j_m, n_found, j_n
      INTEGER       ::                     fluxs_arr_i
      REAL(rprec)   ::                     hs
      REAL(rprec)   ::                     reff
      REAL(rprec)   ::                     psi,dpsi
      REAL(rprec)   ::                     b_ref,r_ref
      REAL(rprec), ALLOCATABLE ::          save_array(:,:)
      INTEGER, ALLOCATABLE :: mnum(:)
      INTEGER       :: mystart,myend, chunk, numprocs_local
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
!DEC$ IF DEFINED (NEO_OPT)
      IF (lscreen) WRITE(6,'(a)') ' ---------------------  NEOCLASSICAL TRANSPORT CALCULATION  ------------------'
      SELECT CASE(TRIM(equil_type))
         CASE('vmec2000','animec','flow','satire','parvmec','paravmec','vboot')
!DEC$ IF DEFINED (MPI_OPT)
            CALL MPI_COMM_SIZE( MPI_COMM_MYWORLD, numprocs_local, ierr_mpi )
            CALL MPI_BCAST(mnboz_b,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(mboz_b,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(nboz_b,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(nfp_b,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(ns_b,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(aspect_b,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(rmax_b,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(rmin_b,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(betaxis_b,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(lasym_b,1,MPI_LOGICAL,master,MPI_COMM_MYWORLD,ierr_mpi)
            ! Handle allocations
            IF (myworkid /= master) THEN
               !1D
               IF (ALLOCATED(idx_b))  DEALLOCATE(idx_b);  ALLOCATE(idx_b(ns_b))
               IF (ALLOCATED(iota_b)) DEALLOCATE(iota_b); ALLOCATE(iota_b(ns_b))
               IF (ALLOCATED(pres_b)) DEALLOCATE(pres_b); ALLOCATE(pres_b(ns_b))
               IF (ALLOCATED(phip_b)) DEALLOCATE(phip_b); ALLOCATE(phip_b(ns_b))
               IF (ALLOCATED(phi_b))  DEALLOCATE(phi_b);  ALLOCATE(phi_b(ns_b))
               IF (ALLOCATED(beta_b)) DEALLOCATE(beta_b); ALLOCATE(beta_b(ns_b))
               IF (ALLOCATED(buco_b)) DEALLOCATE(buco_b); ALLOCATE(buco_b(ns_b))
               IF (ALLOCATED(bvco_b)) DEALLOCATE(bvco_b); ALLOCATE(bvco_b(ns_b))
               IF (ALLOCATED(ixm_b))  DEALLOCATE(ixm_b);  ALLOCATE(ixm_b(mnboz_b))
               IF (ALLOCATED(ixn_b))  DEALLOCATE(ixn_b);  ALLOCATE(ixn_b(mnboz_b))
               !2D
               IF (ALLOCATED(rmnc_b))  DEALLOCATE(rmnc_b);  ALLOCATE(rmnc_b(mnboz_b,ns_b))
               IF (ALLOCATED(zmns_b))  DEALLOCATE(zmns_b);  ALLOCATE(zmns_b(mnboz_b,ns_b))
               IF (ALLOCATED(bmnc_b))  DEALLOCATE(bmnc_b);  ALLOCATE(bmnc_b(mnboz_b,ns_b))
               IF (ALLOCATED(pmns_b))  DEALLOCATE(pmns_b);  ALLOCATE(pmns_b(mnboz_b,ns_b))
               IF (ALLOCATED(gmnc_b))  DEALLOCATE(gmnc_b);  ALLOCATE(gmnc_b(mnboz_b,ns_b))
               IF (lasym_b) THEN
                  IF (ALLOCATED(rmns_b))  DEALLOCATE(rmns_b);  ALLOCATE(rmns_b(mnboz_b,ns_b))
                  IF (ALLOCATED(zmnc_b))  DEALLOCATE(zmnc_b);  ALLOCATE(zmnc_b(mnboz_b,ns_b))
                  IF (ALLOCATED(bmns_b))  DEALLOCATE(bmns_b);  ALLOCATE(bmns_b(mnboz_b,ns_b))
                  IF (ALLOCATED(pmnc_b))  DEALLOCATE(pmnc_b);  ALLOCATE(pmnc_b(mnboz_b,ns_b))
                  IF (ALLOCATED(gmns_b))  DEALLOCATE(gmns_b);  ALLOCATE(gmns_b(mnboz_b,ns_b))
               END IF
            END IF
            ! Broadcast the variables
            CALL MPI_BCAST(idx_b,ns_b,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(ixm_b,mnboz_b,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(ixn_b,mnboz_b,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(iota_b,ns_b,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(pres_b,ns_b,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(phip_b,ns_b,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(phi_b,ns_b,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(beta_b,ns_b,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(buco_b,ns_b,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(bvco_b,ns_b,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(rmnc_b,ns_b*mnboz_b,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(zmns_b,ns_b*mnboz_b,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(bmnc_b,ns_b*mnboz_b,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(gmnc_b,ns_b*mnboz_b,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(pmns_b,ns_b*mnboz_b,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            IF (lasym_b) THEN
               CALL MPI_BCAST(rmns_b,ns_b*mnboz_b,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(zmnc_b,ns_b*mnboz_b,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(bmns_b,ns_b*mnboz_b,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(gmns_b,ns_b*mnboz_b,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(pmnc_b,ns_b*mnboz_b,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            END IF
!DEC$ ENDIF
            ! CALL neo_read_control (Assume NEO namelist has been read)
            in_file  = 'dummy'
            extension = TRIM(proc_string)
            out_file = 'neo_out.'//TRIM(proc_string)
            cur_file = 'neo_cur.'//TRIM(proc_string)
            IF (no_fluxs == 0) THEN
               no_fluxs = COUNT(sigma_neo < bigno)
               IF (ALLOCATED(fluxs_arr)) DEALLOCATE(fluxs_arr)
               ALLOCATE(fluxs_arr(no_fluxs), STAT = iflag)
               IF (iflag /= 0) RETURN
               k = 1
               DO i = 1, UBOUND(sigma_neo, DIM=1)
                  IF (sigma_neo(i) < bigno) THEN
                     fluxs_arr(k) = i
                     k = k + 1
                  END IF
                  IF (k > no_fluxs) EXIT
               END DO
            END IF
            ! CALL neo_init(npsi)
            !   CALL neo_read (boozer data)
            lasym = lasym_b
            m0b = mboz_b-1
            n0b = nboz_b
            mnmax = mnboz_b   
            ns = ns_b
            nfp = nfp_b
            flux = abs( phi_b(ns_b) )
            !max_n_mode = max_n_mode * nfp
            k = 0
            if (max_m_mode .le. 0) max_m_mode = m0b
            if (max_n_mode .le. 0) max_n_mode = n0b * nfp
            do j = 1, mnboz_b
              if( ixm_b(j) .gt. max_m_mode ) cycle
              if( iabs(ixn_b(j) ) .gt. max_n_mode ) cycle
              k = k + 1
            end do
            mnmax = k
            m_max = min0(m0b, max_m_mode)
            m_max = m_max + 1
            n_max = min0(n0b, max_n_mode/nfp)
            n_max = 2*n_max+1
            if ((no_fluxs.gt.0) .and. (no_fluxs.le.ns_b)) ns = no_fluxs
            IF (ALLOCATED(ixm)) DEALLOCATE(ixm); ALLOCATE(ixm(mnmax))
            IF (ALLOCATED(ixn)) DEALLOCATE(ixn); ALLOCATE(ixn(mnmax))
            IF (ALLOCATED(pixm)) DEALLOCATE(pixm); ALLOCATE(pixm(mnmax))
            IF (ALLOCATED(pixn)) DEALLOCATE(pixn); ALLOCATE(pixn(mnmax))
            IF (ALLOCATED(i_m)) DEALLOCATE(i_m); ALLOCATE(i_m(m_max))
            IF (ALLOCATED(i_n)) DEALLOCATE(i_n); ALLOCATE(i_n(n_max))
            IF (ALLOCATED(es)) DEALLOCATE(es); ALLOCATE(es(ns))
            IF (ALLOCATED(iota)) DEALLOCATE(iota); ALLOCATE(iota(ns))
            IF (ALLOCATED(curr_pol)) DEALLOCATE(curr_pol); ALLOCATE(curr_pol(ns))
            IF (ALLOCATED(curr_tor)) DEALLOCATE(curr_tor); ALLOCATE(curr_tor(ns))
            IF (ALLOCATED(pprime)) DEALLOCATE(pprime); ALLOCATE(pprime(ns))
            IF (ALLOCATED(sqrtg00)) DEALLOCATE(sqrtg00); ALLOCATE(sqrtg00(ns))
            IF (ALLOCATED(rmnc)) DEALLOCATE(rmnc); ALLOCATE(rmnc(ns,mnmax))
            IF (ALLOCATED(zmns)) DEALLOCATE(zmns); ALLOCATE(zmns(ns,mnmax))
            IF (ALLOCATED(lmns)) DEALLOCATE(lmns); ALLOCATE(lmns(ns,mnmax))
            IF (ALLOCATED(bmnc)) DEALLOCATE(bmnc); ALLOCATE(bmnc(ns,mnmax))
            IF (lasym) THEN
               IF (ALLOCATED(rmns)) DEALLOCATE(rmns); ALLOCATE(rmns(ns,mnmax))
               IF (ALLOCATED(zmnc)) DEALLOCATE(zmnc); ALLOCATE(zmnc(ns,mnmax))
               IF (ALLOCATED(lmnc)) DEALLOCATE(lmnc); ALLOCATE(lmnc(ns,mnmax))
               IF (ALLOCATED(bmns)) DEALLOCATE(bmns); ALLOCATE(bmns(ns,mnmax))
            END IF
            k = 0;   mn0 = 0
            do j = 1, mnboz_b
               if( ixm_b(j) .gt. max_m_mode ) cycle
               if( iabs(ixn_b(j) ) .gt. max_n_mode ) cycle
               if (ixm_b(j).eq.0 .and. ixn_b(j).eq.0) mn0 = j
               k = k + 1
               ixm(k) = ixm_b(j)
               ixn(k) = ixn_b(j)
            enddo
            hs = ONE/(ns_b -1)
            do i = 1, ns                  !!NEED TO INCLUDE IN LOOP CASE FOR NO_FLUXS < 0 (SPH)
               if (no_fluxs > 0) then
                  i_surf = fluxs_arr(i)  !!NOT ALLOCATED YET IF NO_FLUXS<0
               else
                  i_surf = i
               end if
               k = 0
               do j = 1, mnboz_b
                 if( ixm_b(j) .gt. max_m_mode ) cycle
                 if( abs(ixn_b(j)) .gt. max_n_mode ) cycle
                 k = k + 1
                 rmnc(i,k) = rmnc_b(j,i_surf)
                 zmns(i,k) = zmns_b(j,i_surf)
                 lmns(i,k) = -pmns_b(j,i_surf)*nfp_b/twopi
                 bmnc(i,k) = bmnc_b(j,i_surf)
                 if (ixm_b(j).eq.0 .and. ixn_b(j).eq.0)            &
                     sqrtg00(i) = gmnc_b(j,i_surf)
                 if (lasym) then
                    rmns(i,k) = rmns_b(j,i_surf)
                    zmnc(i,k) = zmnc_b(j,i_surf)
                    lmnc(i,k) = -pmnc_b(j,i_surf)*nfp_b/twopi
                    bmns(i,k) = bmns_b(j,i_surf)
                    if (ixm_b(j).eq.0 .and. ixn_b(j).eq.0)            &
                     sqrtg00(i) = gmnc_b(j,i_surf) + gmns_b(j,i_surf)
                 end if
               enddo
               if (rmnc(i,mn0) .eq. ZERO) then
                  IF (lscreen) write (6, *)' The surface i = ', i_surf, ' is absent in the BOOZMN file'
                  iflag = -1
                  RETURN
               end if
               es(i) = (i_surf-1.5_dp)*hs
               iota(i) = iota_b(i_surf)
               if (i_surf .lt. ns_b) then
                  pprime(i) = (pres_b(i_surf+1) - pres_b(i_surf))/hs
               else
                  pprime(i) = 0
               end if
               curr_pol(i) = bvco_b(i_surf)
               curr_tor(i) = buco_b(i_surf)
            enddo
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
            npsi = ns
            !   CALL neo_prep (boozer data)
            CALL neo_prep
            rt0=0
            bmref=0
            DO imn=1,mnmax
              IF(ixm(imn).EQ.0 .AND. ixn(imn).EQ.0) THEN
                 rt0 = rmnc(1,imn)
                 bmref = bmnc(1,imn)
                 rt0_g = rt0
                 bmref_g = bmref
              ENDIF
            ENDDO
            IF(rt0.EQ.ZERO .OR. bmref.EQ.ZERO) THEN
               IF (lscreen) WRITE (6,*) ' NEO_INIT: Fatal problem setting rt0 or bmref'
               iflag = -1
               RETURN
            ENDIF
            nper = nfp
            !Open files
            w_u6_open = 0
            w_u3 = 222
            istat = 0
            IF (myworkid == master) call safe_open(w_u3, istat, out_file, 'replace', 'formatted')
!DEC$ IF DEFINED (MPI_OPT)
            CALL MPI_BCAST(istat,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
!DEC$ ENDIF
            !if (istat .ne. 0) stop 'Error opening NEO output file'
            IF (istat .ne. 0) THEN
               IF (myworkid == master) PRINT *,istat,out_file
               IF (myworkid == master) WRITE(6,*) 'Error opening NEO output file:',TRIM(out_file),istat
               iflag = -1
               RETURN
            END IF
            istat = 0
            IF (calc_cur .EQ. 1 .and. myworkid == master) THEN
              !OPEN(unit=w_u9,file=cur_file)
              call safe_open(w_u9, istat, cur_file, 'replace', 'formatted')
            END IF
!DEC$ IF DEFINED (MPI_OPT)
            CALL MPI_BCAST(istat,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
!DEC$ ENDIF
            IF (istat .ne. 0) THEN
               IF (myworkid == master) PRINT *,istat,out_file
               IF (myworkid == master) WRITE(6,*) 'Error opening NEO Current file:',TRIM(out_file),istat
               iflag = -1
               RETURN
            END IF
            ! Loop over magnetic surfaces
            reff = 0
            IF (ALLOCATED(eff_ripple)) DEALLOCATE(eff_ripple)
            ALLOCATE(eff_ripple(ns_b))
            eff_ripple = 0
            IF (lscreen) THEN
               WRITE(6,'(A)')           '==========================================='
               WRITE(6,'(A)')           '=================  N E O =================='
               WRITE(6,'(2X,2(A,I3.3))')'MAX_M_MODE = ',max_m_mode,';   M_BOOZER = ',mboz_b 
               WRITE(6,'(2X,2(A,I3.3))')'MAX_N_MODE = ',max_n_mode,';   N_BOOZER = ',nboz_b
               WRITE(6,'(2X,2(A,I3.3))')'THETA_N = ',theta_n,';   PHI_N = ',phi_n
               WRITE(6,'(2X,1(A,I3.3))')'NPARTICLES = ',npart
               WRITE(6,'(2X,A,E14.4)')  'ACCURACY   = ',acc_req
               WRITE(6,'(2X,2(A,I4.4))')'NSTEP_MIN = ',nstep_min,';   NSTEP_MAX = ',nstep_max
               WRITE(6,'(2X,A,I3.3)')   'NSTEP_PER = ',nstep_per
               WRITE(6,'(A)')           '------------------------------------------------------------------'
               WRITE(6,'(A)')           '  SURF  EPS_EFF        REFF       IOTA        B_REF       R_REF'
               WRITE(6,'(A)')           '------------------------------------------------------------------'
               CALL FLUSH(6)
            END IF
 
            ! Divide up work
            IF (ALLOCATED(mnum)) DEALLOCATE(mnum)
            ALLOCATE(mnum(numprocs_local))
            mnum=0
            i = 1
            DO
               IF (SUM(mnum,DIM=1) == no_fluxs) EXIT  ! Have to use ns_b because of logic
               IF (i > numprocs_local) i = 1
               mnum(i) = mnum(i) + 1
               i=i+1
            END DO
            mystart = 1
            DO i = 1, myworkid
               mystart = SUM(mnum(1:i))+1
            END DO
            myend = mystart + mnum(myworkid+1) - 1
            IF (myend < mystart) myend = mystart
            IF (mnum(myworkid+1) == 0) mystart = myend + 1
            DEALLOCATE(mnum)

            ALLOCATE(save_array(no_fluxs,17)); save_array=0
            DO fluxs_arr_i = mystart,myend
               psi_ind = fluxs_arr_i
               !IF (psi_ind .GE. 1 .AND. psi_ind .LE. npsi) THEN
               CALL neo_init_s(psi,dpsi)
               IF(psi_ind.EQ.1) dpsi=psi
               CALL flint_bo()
               reff=reff+drdpsi*dpsi
               IF (reff /= reff) reff = 0
               IF (calc_cur .EQ. 1) CALL flint_cur()
               IF (ref_swi .EQ. 1) THEN
                  b_ref = bmref_g
                  r_ref = rt0_g
               ELSE IF (ref_swi .EQ. 2) THEN
                  b_ref = bmref
                  r_ref = rt0
               ELSE
                  ! This is checked above and would break the code if caught here
                  !IF (lscreen) WRITE (6,*) 'FATAL: This ref_swi ',ref_swi,' is not implemented!'
                  !IF (myworkid == master) CLOSE(w_u3)
                  !iflag = -1
                  !RETURN
               END IF
               epstot = epstot * (b_ref/bmref)**2 * (r_ref/rt0)**2
               epspar = epspar * (b_ref/bmref)**2 * (r_ref/rt0)**2
               IF (epstot /= epstot) epstot = 0 
               save_array(fluxs_arr_i,1) = epstot
               save_array(fluxs_arr_i,2) = reff
               save_array(fluxs_arr_i,3) = iota(psi_ind)
               save_array(fluxs_arr_i,4) = b_ref
               save_array(fluxs_arr_i,5) = r_ref
               save_array(fluxs_arr_i,6) = epspar(1)
               IF (multra >1) save_array(fluxs_arr_i,7) = epspar(2)
               save_array(fluxs_arr_i,8) = ctrone
               save_array(fluxs_arr_i,9) = ctrtot
               save_array(fluxs_arr_i,10) = bareph
               save_array(fluxs_arr_i,11) = barept
               save_array(fluxs_arr_i,12) = yps
               save_array(fluxs_arr_i,13) = lambda_b
               save_array(fluxs_arr_i,14) = lambda_ps1
               save_array(fluxs_arr_i,15) = lambda_ps2
               save_array(fluxs_arr_i,16) = lambda_b1
               save_array(fluxs_arr_i,17) = lambda_b2
            END DO
!DEC$ IF DEFINED (MPI_OPT)
            IF (myworkid == master) THEN
               CALL MPI_REDUCE(MPI_IN_PLACE,save_array,no_fluxs*17,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            ELSE
               CALL MPI_REDUCE(save_array,save_array,no_fluxs*17,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            END IF
!DEC$ ENDIF
            IF (myworkid == master) THEN
               DO fluxs_arr_i = 1, no_fluxs
                  psi_ind = fluxs_arr_i
                  epstot = save_array(fluxs_arr_i,1)
                  reff = save_array(fluxs_arr_i,2)
                  iota(psi_ind) = save_array(fluxs_arr_i,3)
                  b_ref = save_array(fluxs_arr_i,4)
                  r_ref = save_array(fluxs_arr_i,5)
                  epspar(1) = save_array(fluxs_arr_i,6)
                  IF (multra >1) epspar(2) = save_array(fluxs_arr_i,7)
                  ctrone = save_array(fluxs_arr_i,8)
                  ctrtot = save_array(fluxs_arr_i,9)
                  bareph = save_array(fluxs_arr_i,10)
                  barept = save_array(fluxs_arr_i,11)
                  yps = save_array(fluxs_arr_i,12)
                  lambda_b = save_array(fluxs_arr_i,13)
                  lambda_ps1 = save_array(fluxs_arr_i,14)
                  lambda_ps2 = save_array(fluxs_arr_i,15)
                  lambda_b1 = save_array(fluxs_arr_i,16)
                  lambda_b2 = save_array(fluxs_arr_i,17)
                  IF (eout_swi .EQ. 1) THEN
                     WRITE(w_u3,'(1(1x,i8),5(1x,e17.10))')                    &
                           fluxs_arr(fluxs_arr_i),                             &
                           epstot,reff,iota(psi_ind),b_ref,r_ref
                  ELSEIF (eout_swi .EQ. 2) THEN
                     IF (multra >1)THEN
                         WRITE(w_u3,'(1(1x,i8),12(1x,e17.10))')                   &
                           fluxs_arr(fluxs_arr_i),                             &
                           epstot,reff,iota(psi_ind),b_ref,r_ref,              &
                           epspar(1),epspar(2),ctrone,ctrtot,bareph,barept,yps
                     ELSE 
                           WRITE(w_u3,'(1(1x,i8),12(1x,e17.10))')                   &
                           fluxs_arr(fluxs_arr_i),                             &
                           epstot,reff,iota(psi_ind),b_ref,r_ref,              &
                           epspar(1),0.0,ctrone,ctrtot,bareph,barept,yps
                     END IF
                  ELSEIF (eout_swi .EQ. 10) THEN                              !LPK
                     WRITE(w_u3,*) b_ref, r_ref, epstot                       !LPK
                  ELSE
                     IF (lscreen) WRITE(6,*) 'FATAL: This eout_swi ',eout_swi,' is not implemented!'
                     IF (myworkid == master) CLOSE(w_u3)
                     iflag = -1
                     RETURN
                  END IF
                  IF (calc_cur .EQ. 1) THEN
                     WRITE(w_u9,'(1(1x,i8),5(1x,e17.10))')                       &
                        psi_ind,                                               &
                        lambda_b,                                              &
                        lambda_ps1,lambda_ps2,                                 &
                        lambda_b1,lambda_b2
                  END IF
                  eff_ripple(fluxs_arr(fluxs_arr_i)) = epstot
                  IF (lscreen) WRITE(6,'(2X,I3,5(2X,E11.4))') fluxs_arr(fluxs_arr_i),epstot,reff,iota(psi_ind),b_ref,r_ref
                  CALL FLUSH(6)
               END DO
               CLOSE(w_u3)
               IF (calc_cur .EQ. 1) THEN
                 CLOSE(w_u9)
               END IF
            END IF
            CALL neo_dealloc
            no_fluxs = 0  ! For next pass through
            IF (lscreen) WRITE(6,'(A)')           '------------------------------------------------------------------'
         CASE('spec')
      END SELECT
      IF (lscreen) WRITE(6,'(a)') ' -----------------  NEOCLASSICAL TRANSPORT CALCULATION (DONE) ----------------'
!DEC$ ENDIF
      RETURN
  90  format(5e16.8)
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_neo
