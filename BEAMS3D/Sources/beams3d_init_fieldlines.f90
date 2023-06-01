!-----------------------------------------------------------------------
!     Subroutine:    beams3d_init_fieldlines
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          05/31/2023
!     Description:   This subroutine reads the FIELDLINS files.
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_init_fieldlines
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE beams3d_runtime
      USE beams3d_grid, ONLY: raxis_g => raxis, phiaxis, &
                                 zaxis_g => zaxis, nr, nphi, nz, &
                                 rmin, rmax, zmin, zmax, phimin, &
                                 phimax, vc_adapt_tol, B_R, B_Z, B_PHI,&
                                 BR_spl, BZ_spl, TE_spl_s, NE_spl_s, TI_spl_s, &
                                 nte, nne, nti, TE, NE, TI, Vp_spl_s, S_ARR,&
                                 U_ARR, POT_ARR, POT_spl_s, nne, nte, nti, npot, &
                                 ZEFF_spl_s, nzeff, ZEFF_ARR, req_axis, zeq_axis, &
                                 phiedge_eq, reff_eq, NI_spl_s, NI
      USE beams3d_lines, ONLY: GFactor, ns_prof1
      USE read_fieldlines_mod, ONLY: get_fieldlines_grid, get_fieldlines_B, &
                               read_fieldlines_deallocate, setup_fieldlines_rhogrid, &
                               get_fieldlines_magaxis, get_fieldlines_gridB
      USE mpi_params
      USE mpi_inc
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, PARAMETER :: BYTE_8 = SELECTED_INT_KIND (8)
#if defined(MPI_OPT)
      INTEGER :: numprocs_local, mylocalid, mylocalmaster, mystart, myend
      INTEGER :: MPI_COMM_LOCAL
#endif
      LOGICAL :: lcreate_wall
      INTEGER :: ier, s, i, j, k, u
      REAL(rprec) :: brtemp, bptemp, bztemp, betatot, sflx, uflx, &
                     tetemp,netemp,titemp,zetemp,pottemp, rminor
      INTEGER :: nrh,nzh,nph
      REAL(rprec) :: rmin_hint, rmax_hint, zmin_hint, zmax_hint, &
                     pmax_hint, pres_max
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Use the EQDSK limiter if not wall supplied
      !lcreate_wall = .not. lvessel

      ! Divide up Work
      mylocalid = myworkid
      numprocs_local = 1
      ierr_mpi = 0
#if defined(MPI_OPT)
      CALL MPI_COMM_DUP( MPI_COMM_SHARMEM, MPI_COMM_LOCAL, ierr_mpi)
      CALL MPI_COMM_RANK( MPI_COMM_LOCAL, mylocalid, ierr_mpi )              ! MPI
      CALL MPI_COMM_SIZE( MPI_COMM_LOCAL, numprocs_local, ierr_mpi )          ! MPI
#endif
      mylocalmaster = master

      ! Write info to screen
      IF (lverb) THEN
         betatot = 0
         CALL get_fieldlines_grid(nrh,nzh,nph,rmin_hint,rmax_hint,zmin_hint,zmax_hint,pmax_hint)
         WRITE(6,'(A)')               '----- FIELDLINES Information -----'
         WRITE(6,'(A,F9.5,A,F9.5,A,I4)') '   R   = [',rmin_hint,',',rmax_hint,'];  NR:   ',nrh
         WRITE(6,'(A,F8.5,A,F8.5,A,I4)') '   PHI = [',0.0,',',pmax_hint,'];  NPHI: ',nph
         WRITE(6,'(A,F8.5,A,F8.5,A,I4)') '   Z   = [',zmin_hint,',',zmax_hint,'];  NZ:   ',nzh
         WRITE(6,'(A,F6.3,A)')           '   AMINOR_NORM = ',rminor_norm,' [m]'
      END IF

      CALL setup_fieldlines_rhogrid('fieldlines_'//TRIM(id_string)//'.h5',MPI_COMM_SHARMEM,ier)

      ! Calculate the axis values
      DO s = 1, nphi
         CALL get_fieldlines_magaxis(phiaxis(s),req_axis(s),zeq_axis(s))
      END DO

      
      IF (lverb) THEN
         WRITE(6,'(5X,A,I3.3,A)',ADVANCE='no') 'Plasma Field Calculation [',0,']%'
         CALL FLUSH(6)
      END IF
      
      ! Break up the Work
      CALL MPI_CALC_MYRANGE(MPI_COMM_LOCAL, 1, nr*nphi*nz, mystart, myend)

      IF (mylocalid == mylocalmaster) THEN
         TE = 0; NE = 0; TI=0; S_ARR=1.5; U_ARR=0; POT_ARR=0; ZEFF_ARR = 1;
      END IF
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
#endif
      DO s = mystart, myend
         i = MOD(s-1,nr)+1
         j = MOD(s-1,nr*nphi)
         j = FLOOR(REAL(j) / REAL(nr))+1
         k = CEILING(REAL(s) / REAL(nr*nphi))
         sflx = 0.0

         ! Bfield
         CALL get_fieldlines_gridB(i,j,k,brtemp,bptemp,bztemp,rminor,uflx)
         B_R(i,j,k) = brtemp
         B_PHI(i,j,k) = bptemp
         B_Z(i,j,k) = bztemp
         sflx = (rminor/rminor_norm)**2
         S_ARR(i,j,k) = sflx
         U_ARR(i,j,k) = uflx

         IF (sflx <= 1) THEN
            tetemp = 0; netemp = 0; titemp=0; pottemp=0; zetemp=0
            IF (nte > 0) CALL EZspline_interp(TE_spl_s,sflx,tetemp,ier)
            IF (nne > 0) CALL EZspline_interp(NE_spl_s,sflx,netemp,ier)
            IF (nti > 0) CALL EZspline_interp(TI_spl_s,sflx,titemp,ier)
            IF (npot > 0) CALL EZspline_interp(POT_spl_s,sflx,pottemp,ier)
            IF (nzeff > 0) THEN 
               CALL EZspline_interp(ZEFF_spl_s,sflx,ZEFF_ARR(i,j,k),ier)
               DO u=1, NION
                  CALL EZspline_interp(NI_spl_s(u),sflx,NI(u,i,j,k),ier)
               END DO
            END IF
            NE(i,j,k) = netemp; TE(i,j,k) = tetemp; TI(i,j,k) = titemp
            POT_ARR(i,j,k) = pottemp
         END IF
         IF (MOD(s,nr) == 0) THEN
            IF (lverb) THEN
               CALL backspace_out(6,6)
               WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*s)/(myend-mystart+1)),']%'
               CALL FLUSH(6)
            END IF
         END IF
      END DO
      
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
#endif
      

      ! Calculate Gfactor for mgrid
      IF (myworkid == master) THEN
         ! Only master has Gfactor
         ALLOCATE(Gfactor(ns_prof1))
         Gfactor = 1
!         DO s = 1, ns_prof1
!            sflx = REAL(s-1 + 0.5)/ns_prof1 ! Half grid rho
!            sflx = sflx*sflx
!            CALL EZspline_interp(ZEFF_spl_s,sflx,uflx,ier)
            ! sflx=s uflx=Zeff
!            CALL vmec_ohkawa(sflx,uflx,Gfactor(s)) ! (1-l31)/Zeff
!         END DO
      ENDIF

      ! Deallocations
!      IF (myworkid == master) THEN
         CALL read_fieldlines_deallocate
!      ELSE
!         DEALLOCATE(vp,phi)
!         DEALLOCATE(xm,xn,xm_nyq,xn_nyq)
!         DEALLOCATE(rmnc,zmns,bsupumnc,bsupvmnc)
!         IF (lasym) DEALLOCATE(rmns,zmnc,bsupumns,bsupvmns)
!      END IF
      
      IF (lverb) THEN
         CALL backspace_out(6,36)
         CALL FLUSH(6)
         WRITE(6,'(36X)',ADVANCE='no')
         CALL FLUSH(6)
         CALL backspace_out(6,36)
         WRITE(6,*)
         CALL FLUSH(6)
      END IF    

      ! Fix ZEFF
      IF (mylocalid == mylocalmaster) WHERE(ZEFF_ARR < 1) ZEFF_ARR = 1

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
      CALL MPI_COMM_FREE(MPI_COMM_LOCAL,ierr_mpi)
      CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'beams3d_init_fieldlines',ierr_mpi)
#endif

      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE beams3d_init_fieldlines
