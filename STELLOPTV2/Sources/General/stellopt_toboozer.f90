!-----------------------------------------------------------------------
!     Subroutine:    stellopt_toboozer
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          06/06/2012
!     Description:   This subroutine handles conversion to boozer
!                    coordinates.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_toboozer(lscreen,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_input_mod
      USE stellopt_vars
      USE vmec_input
      USE booz_params, mboz_xboozer => mboz, nboz_xboozer => nboz,&
                       lscreen_xboozer => lscreen, nfp_xboozer => nfp, &
                       ns_xboozer => ns, lasym_xboozer => lasym_b, &
                       mpol_b => mpol, ntor_b => ntor
      USE read_wout_mod, ONLY: input_extension_wout => input_extension,&
                               ns_vmec => ns, aspect_vmec => aspect,&
                               rmax_vmec => rmax_surf, rmin_vmec=> rmin_surf,&
                               betaxis_vmec => betaxis
      USE safe_open_mod, ONLY: safe_open
      USE booz_persistent
      USE read_boozer_mod
      USE equil_vals
      USE equil_utils, ONLY: mntouv
      USE mpi_params
      
!-----------------------------------------------------------------------
!     Subroutine Parameters
!        iflag         Error flag
!----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL, INTENT(in)    :: lscreen
      INTEGER, INTENT(inout) :: iflag
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      INTEGER ::  ier, ik, ij, iunit, nu, nv, u, v, mn, mystart,myend, chunk, numprocs_local, num_booz
      INTEGER,SAVE ::  irun_setup_booz = 0
      INTEGER, ALLOCATABLE :: im(:), in(:), mnum(:)
      REAL(rprec), ALLOCATABLE :: xu(:), xv(:)
      REAL(rprec), ALLOCATABLE :: r_temp(:,:,:), z_temp(:,:,:), p_temp(:,:,:), g_temp(:,:,:), b_temp(:,:,:)
      REAL(rprec), ALLOCATABLE :: ru(:,:,:), rv(:,:,:)
      REAL(rprec), ALLOCATABLE :: zu(:,:,:), zv(:,:,:)
      REAL(rprec), ALLOCATABLE :: pu(:,:,:), pv(:,:,:)
      REAL(rprec), ALLOCATABLE :: fmn_temp(:,:)
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (lscreen) WRITE(6,'(a)') ' ---------------------------  BOOZER TRANSFORMATION  -------------------------'
      SELECT CASE(TRIM(equil_type))
         CASE('vmec2000','animec','flow','satire','parvmec','paravmec','vboot')
            lscreen_xboozer = lscreen
            ! We need to pass mboz and nboz to the boozer routines
            mboz_xboozer = mboz
            nboz_xboozer = nboz
!DEC$ IF DEFINED (MPI_OPT)
            CALL MPI_COMM_SIZE( MPI_COMM_MYWORLD, numprocs_local, ierr_mpi )
            CALL MPI_BCAST(mboz_xboozer,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(mboz_xboozer,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(ns_vmec,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(lbooz,nsd,MPI_LOGICAL,master,MPI_COMM_MYWORLD,ierr_mpi)
!DEC$ ENDIF
            ! Now read the wout file
            IF (ALLOCATED(lsurf_boz)) DEALLOCATE(lsurf_boz)
            ALLOCATE(lsurf_boz(ns_vmec))
            lsurf_boz(1:ns_vmec) = lbooz(1:ns_vmec)
            num_booz = -COUNT(lbooz) ! This tricks the code into not reading the wout or input file
            ier = 0
            !CALL read_wout_booz(TRIM(input_extension),iunit,ier)
            !CLOSE(iunit)
            !IF (ier .ne. 0) THEN; iflag=ier; RETURN; END IF
!DEC$ IF DEFINED (MPI_OPT)
            IF (myworkid == master) THEN
               CALL read_wout_booz(TRIM(input_extension),num_booz,ier)
               CLOSE(iunit)
            END IF
            CALL MPI_BCAST(ier,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            IF (ier .ne. 0) THEN; iflag=ier; RETURN; END IF
            IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_toboozer: BARRIER1',ierr_mpi)
            CALL MPI_BCAST(lasym_xboozer,1,MPI_LOGICAL,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(ier,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(mpol_b,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(mpol1,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(mpol_nyq,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(ntor_nyq,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(nfp_xboozer,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(ns_xboozer,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(ntor_b,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(mnmax,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(mnmax_nyq,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(mboz_xboozer,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(nboz_xboozer,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(nu_boz,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(nv_boz,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(nunv,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(mnboz,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(nu2_b,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            iunit = -COUNT(lbooz) ! This tricks the code into not reading the wout or input file
            IF (myworkid /= master) CALL allocate_boozer(iunit)
            CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(xm,mnmax,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(xn,mnmax,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(xm_nyq,mnmax_nyq,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(xn_nyq,mnmax_nyq,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(rmnc,mnmax*ns_xboozer,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(zmns,mnmax*ns_xboozer,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(lmns,mnmax*ns_xboozer,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(bsubumnc,mnmax_nyq*ns_xboozer,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(bsubvmnc,mnmax_nyq*ns_xboozer,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(bmodmnc,mnmax_nyq*ns_xboozer,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            IF (lasym_xboozer) THEN
               CALL MPI_BCAST(rmns,mnmax*ns_xboozer,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(zmnc,mnmax*ns_xboozer,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(lmnc,mnmax*ns_xboozer,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(bsubumns,mnmax_nyq*ns_xboozer,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(bsubvmns,mnmax_nyq*ns_xboozer,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_BCAST(bmodmns,mnmax_nyq*ns_xboozer,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            END IF
            CALL MPI_BCAST(hiota,ns_xboozer,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(phip,ns_xboozer,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(pres,ns_xboozer,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(beta_vol,ns_xboozer,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(phi,ns_xboozer,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(buco,ns_xboozer,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(bvco,ns_xboozer,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)

            IF (ALLOCATED(mnum)) DEALLOCATE(mnum)
            ALLOCATE(mnum(numprocs_local))
            mnum=0
            ij = 1
            DO
               IF (SUM(mnum,DIM=1) == ns_xboozer) EXIT  ! Have to use ns_b because of logic
               IF (ij > numprocs_local) ij = 1
               mnum(ij) = mnum(ij) + 1
               ij=ij+1
            END DO
            mystart = 1
            DO ij = 1, myworkid
               mystart = SUM(mnum(1:ij))+1
            END DO
            myend = mystart + mnum(myworkid+1) - 1
            IF (myend < mystart) myend = mystart
            IF (mnum(myworkid+1) == 0) mystart = myend + 1
            DEALLOCATE(mnum)

            ! Becasue of jsurf we first have everyone call on the 0th surface
            ik = 1
            DO
               IF (lbooz(ik)) EXIT
               ik = ik + 1
            END DO
            irun_setup_booz = 0
            CALL boozer_coords(ik,irun_setup_booz)
            IF (myworkid /= master) THEN
               rmncb = 0
               zmnsb = 0
               pmnsb = 0
               gmncb = 0
               bmncb = 0
               IF (lasym_xboozer) THEN
                  rmnsb = 0
                  zmncb = 0
                  pmncb = 0
                  gmnsb = 0
                  bmnsb = 0
               END IF
               irun_setup_booz = COUNT(lbooz(1:mystart-1))
            ELSE
               mystart = ik+1
            END IF
!DEC$ ENDIF
            !DO ik = 0,numprocs_local-1
            !   IF (myworkid == ik) PRINT *,myworkid,mystart,chunk,myend,irun_setup_booz
            !END DO
            ! Calculate Boozer Coordinates
            DO ik = mystart, myend
               IF (lbooz(ik)) CALL boozer_coords(ik,irun_setup_booz)
            END DO
!DEC$ IF DEFINED (MPI_OPT)
            CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi)
            ik = SIZE(rmncb)
            IF (myworkid == master) THEN
               CALL MPI_REDUCE(MPI_IN_PLACE,bmncb,ik,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,rmncb,ik,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,zmnsb,ik,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,pmnsb,ik,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,gmncb,ik,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               IF (lasym_xboozer) THEN
                  CALL MPI_REDUCE(MPI_IN_PLACE,bmnsb,ik,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
                  CALL MPI_REDUCE(MPI_IN_PLACE,rmnsb,ik,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
                  CALL MPI_REDUCE(MPI_IN_PLACE,zmncb,ik,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
                  CALL MPI_REDUCE(MPI_IN_PLACE,pmncb,ik,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
                  CALL MPI_REDUCE(MPI_IN_PLACE,gmnsb,ik,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               END IF
            ELSE
               CALL MPI_REDUCE(bmncb,bmncb,ik,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(rmncb,rmncb,ik,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(zmnsb,zmnsb,ik,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(pmnsb,pmnsb,ik,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(gmncb,gmncb,ik,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               IF (lasym_xboozer) THEN
                  CALL MPI_REDUCE(bmnsb,bmnsb,ik,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
                  CALL MPI_REDUCE(rmnsb,rmnsb,ik,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
                  CALL MPI_REDUCE(zmncb,zmncb,ik,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
                  CALL MPI_REDUCE(pmncb,pmncb,ik,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
                  CALL MPI_REDUCE(gmnsb,gmnsb,ik,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               END IF
               CALL free_mem_boozer
               irun_setup_booz = 0
               RETURN
            END IF
!DEC$ ENDIF
            ! Output the data
            CALL write_boozmn(TRIM(input_extension))
            ! Load the boozer data from memory
            mnboz_b = mnboz
            mboz_b  = mboz_xboozer
            nboz_b  = nboz_xboozer
            nfp_b   = nfp_xboozer
            ns_b    = ns_xboozer
            aspect_b = aspect_vmec
            rmax_b   = rmax_vmec
            rmin_b   = rmin_vmec
            betaxis_b = betaxis_vmec
            IF (ALLOCATED(idx_b)) DEALLOCATE(idx_b); ALLOCATE(idx_b(ns_b)); idx_b=0
            IF (ALLOCATED(iota_b)) DEALLOCATE(iota_b); ALLOCATE(iota_b(ns_b)); iota_b=0
            IF (ALLOCATED(pres_b)) DEALLOCATE(pres_b); ALLOCATE(pres_b(ns_b)); pres_b=0
            IF (ALLOCATED(phip_b)) DEALLOCATE(phip_b); ALLOCATE(phip_b(ns_b)); phip_b=0
            IF (ALLOCATED(phi_b)) DEALLOCATE(phi_b); ALLOCATE(phi_b(ns_b)); phi_b=0
            IF (ALLOCATED(beta_b)) DEALLOCATE(beta_b); ALLOCATE(beta_b(ns_b)); beta_b=0
            IF (ALLOCATED(buco_b)) DEALLOCATE(buco_b); ALLOCATE(buco_b(ns_b)); buco_b=0
            IF (ALLOCATED(bvco_b)) DEALLOCATE(bvco_b); ALLOCATE(bvco_b(ns_b)); bvco_b=0
            idx_b = 0
            WHERE(lbooz(1:ns_b)) idx_b(1:ns_b) = 1
            iota_b = hiota
            pres_b = pres
            phip_b = phip
            phi_b  = phi
            beta_b = beta_vol
            buco_b = buco
            bvco_b = bvco
            IF (ALLOCATED(ixm_b)) DEALLOCATE(ixm_b); ALLOCATE(ixm_b(mnboz_b)); ixm_b=0
            IF (ALLOCATED(ixn_b)) DEALLOCATE(ixn_b); ALLOCATE(ixn_b(mnboz_b)); ixn_b=0
            ixm_b = NINT(xmb)
            ixn_b = NINT(xnb)
            IF (ALLOCATED(bmnc_b)) DEALLOCATE(bmnc_b); ALLOCATE(bmnc_b(mnboz_b,ns_b)); bmnc_b=0
            IF (ALLOCATED(rmnc_b)) DEALLOCATE(rmnc_b); ALLOCATE(rmnc_b(mnboz_b,ns_b)); rmnc_b=0  
            IF (ALLOCATED(zmns_b)) DEALLOCATE(zmns_b); ALLOCATE(zmns_b(mnboz_b,ns_b)); zmns_b=0  
            IF (ALLOCATED(pmns_b)) DEALLOCATE(pmns_b); ALLOCATE(pmns_b(mnboz_b,ns_b)); pmns_b=0  
            IF (ALLOCATED(gmnc_b)) DEALLOCATE(gmnc_b); ALLOCATE(gmnc_b(mnboz_b,ns_b)); gmnc_b=0
            lasym_b = .FALSE.
            IF (lasym_xboozer) THEN
               lasym_b = .TRUE.
               IF (ALLOCATED(bmns_b)) DEALLOCATE(bmns_b); ALLOCATE(bmns_b(mnboz_b,ns_b)); bmns_b=0
               IF (ALLOCATED(rmns_b)) DEALLOCATE(rmns_b); ALLOCATE(rmns_b(mnboz_b,ns_b)); rmns_b=0
               IF (ALLOCATED(zmnc_b)) DEALLOCATE(zmnc_b); ALLOCATE(zmnc_b(mnboz_b,ns_b)); zmnc_b=0
               IF (ALLOCATED(pmnc_b)) DEALLOCATE(pmnc_b); ALLOCATE(pmnc_b(mnboz_b,ns_b)); pmnc_b=0
               IF (ALLOCATED(gmns_b)) DEALLOCATE(gmns_b); ALLOCATE(gmns_b(mnboz_b,ns_b)); gmns_b=0  
            END IF
            ! The internal BOOZER variables are in packed form (meaning that they run
            ! over the number of surfaces calculated not the total number of
            ! equilibrium surfaces.  Thus we need to get them in the read_boozer_mod
            ! form where the radial index runs over the total number of ns_b
            ! surfaces.
            ij = 1
            DO ik = 1, ns_b
              IF (idx_b(ik) == 0) CYCLE
              bmnc_b(:,ik) = bmncb(:,ij)
              rmnc_b(:,ik) = rmncb(:,ij)
              zmns_b(:,ik) = zmnsb(:,ij)
              pmns_b(:,ik) = pmnsb(:,ij)
              gmnc_b(:,ik) = gmncb(:,ij)
              IF (lasym_xboozer) THEN
                 bmns_b(:,ik) = bmnsb(:,ij)
                 rmns_b(:,ik) = rmnsb(:,ij)
                 zmnc_b(:,ik) = zmncb(:,ij)
                 pmnc_b(:,ik) = pmncb(:,ij)
                 gmns_b(:,ik) = gmnsb(:,ij)
              END IF
              ij = ij + 1
            END DO
            ! Free the memory
            CALL free_mem_boozer
            irun_setup_booz = 0
            return
         CASE('spec')
      END SELECT
      IF (lscreen) WRITE(6,'(a)') ' -------------------------  BOOZER TRANSFORMATION DONE  ----------------------'
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_toboozer
