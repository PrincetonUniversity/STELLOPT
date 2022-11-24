!-----------------------------------------------------------------------
!     Subroutine:    fieldlines_init_eqdsk
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          02/15/2021
!     Description:   This subroutine reads the EFIT files and sets up
!                    grids.
!-----------------------------------------------------------------------
      SUBROUTINE fieldlines_init_eqdsk
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
        USE stel_kinds, ONLY: rprec
        USE read_eqdsk_mod, ONLY: totcur, psisep, &
                                raxis_eqdsk => raxis, &
                                zaxis_eqdsk => zaxis, get_eqdsk_Bspl, &
                                get_eqdsk_fluxspl, read_eqdsk_deallocate, &
                                nlim, xlim, zlim, ntitle, btor, &
                                rcenter, sp, nbndry, xbndry, zbndry, COCOS_ID
        USE fieldlines_runtime
        USE fieldlines_grid, ONLY: raxis_g => raxis, phiaxis, &
                                    zaxis_g => zaxis, nr, nphi, nz, &
                                    rmin, rmax, zmin, zmax, phimin, &
                                    phimax, vc_adapt_tol, B_R, B_Z, B_PHI, &
                                    BR_spl, BZ_spl!, S_ARR,U_ARR
                                    !req_axis, zeq_axis, phiedge_eq, reff_eq
        !USE fieldlines_lines, ONLY: GFactor, ns_prof1
        USE wall_mod, ONLY: wall_load_seg
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
        LOGICAL :: lcreate_wall, lverb_wall
        INTEGER :: ier, s, i, j, k, u
        REAL(rprec) :: brtemp, bptemp, bztemp, betatot, sflx, uflx, &
                        tetemp,netemp,titemp,zetemp,pottemp, rsmax, rsmin,&
                        zsmax, zsmin,rhoflx, nitemp

!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
        ! Use the EQDSK limiter if not wall supplied
        lcreate_wall = .not. lvessel

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

        ! Calculate the minor radius
        IF (myworkid == master) THEN
            brtemp = 0.5*SQRT((raxis_g(nr)-raxis_g(1))**2 + (zaxis_g(nr)-zaxis_g(1))**2)
            k = FLOOR(brtemp/0.001)
            !reff_eq = 0
            s = 0
            DO i = 0, 359, 1
            INNERLOOP: DO j = 1, k
                brtemp = raxis_eqdsk + 0.001*j*COS(i*pi2/360.0)
                bztemp = zaxis_eqdsk + 0.001*j*SIN(i*pi2/360.0)
                CALL get_eqdsk_fluxspl(brtemp,bztemp,rhoflx,uflx)
                IF (rhoflx .ge. 1) THEN
                    !reff_eq = reff_eq + 0.001*(j-0.5) ! halfway between last gridpoint and here.
                    s = s + 1
                    EXIT INNERLOOP
                END IF
            END DO INNERLOOP
            END DO
            !reff_eq = reff_eq / s
        end IF

        ! Write info to screen
        IF (lverb) THEN
            betatot = 0
            WRITE(6,'(A)')               '----- EQDSK Information -----'
            WRITE(6,'(A,I2)') '  COCOS_ID : ',COCOS_ID
            WRITE(6,'(A)')    '  HEADER: '//ntitle(1)//ntitle(2)//ntitle(3)//ntitle(4)//ntitle(5)
            WRITE(6,'(A,F7.2,A,F7.2,A)') '   B  = ',Btor,' [T] @ R = ',rcenter,' [m]'
            IF (ABS(totcur) > 1E8) THEN
            WRITE(6,'(A,F7.2,A)') '   I  = ',totcur*1E-9,' [GA]'
            ELSEIF (ABS(totcur) > 1E5) THEN
            WRITE(6,'(A,F7.2,A)') '   I  = ',totcur*1E-6,' [MA]'
            ELSEIF (ABS(totcur) > 1E2) THEN
            WRITE(6,'(A,F7.2,A)') '   I  = ',totcur*1E-3,' [kA]'
            ELSE
            WRITE(6,'(A,F7.2,A)') '   I  = ',totcur,' [A]'
            END IF
            !WRITE(6,'(A,F7.2,A)')        '   AMINOR  = ',reff_eq,' [m]'
            WRITE(6,'(A,F7.2,A)')        '   P_MAX = ',MAXVAL(sp)*1E-3,' [kPa]'
            WRITE(6,'(A,F7.2,A)')        '   PSIEDGE = ',psisep,' [Wb]'
            IF (lcreate_wall) THEN
            IF (lplasma_only) THEN
                WRITE(6,'(A)')         '   Using EQDSK Separatrix as wall'
            ELSE
                WRITE(6,'(A)')         '   Using EQDSK Limiter as wall'
            END IF
            END IF
        END IF

        ! ! Calculate the axis values
        ! req_axis(:) = raxis_eqdsk
        ! zeq_axis(:) = zaxis_eqdsk

        ! Helpers
        !phiedge_eq = 1.0 ! normalized value

        ! Use the vessel as a mask for rho
        zsmax = MAXVAL(zbndry)
        zsmin = MINVAL(zbndry)
        rsmax = MAXVAL(xbndry)
        rsmin = MINVAL(xbndry)

        ! If we ask for a plasma-only run and do not provide a vessel then make one.
        IF (lcreate_wall) THEN
            lvessel = .TRUE.
            k = 128
            ier = 0
            lverb_wall=.TRUE.
            IF (lplasma_only) THEN
            CALL wall_load_seg(nbndry,xbndry,zbndry,k,ier,VERB=lverb_wall,COMM=MPI_COMM_LOCAL)
            IF (ier==-327) WRITE(6,'(A)') 'ERROR: EQDSK Boundary has repeated index.'
            ELSE
            CALL wall_load_seg(nlim,xlim,zlim,k,ier,VERB=lverb_wall,COMM=MPI_COMM_LOCAL)
            IF (ier==-327) WRITE(6,'(A)') 'ERROR: EQDSK Limiter has repeated index.'
            END IF
            ier = 0
        END IF
        
        IF (lverb) THEN
            WRITE(6,'(5X,A,I3.3,A)',ADVANCE='no') 'Plasma Field Calculation [',0,']%'
            CALL FLUSH(6)
        END IF
        
        ! Break up the Work
        CALL MPI_CALC_MYRANGE(MPI_COMM_LOCAL, 1, nr*nz, mystart, myend)

#if defined(MPI_OPT)
        CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
#endif
        DO s = mystart, myend
            i = MOD(s-1,nr)+1
            k = MOD(s-1,nr*nz)
            k = FLOOR(REAL(k) / REAL(nr))+1
            sflx = 0.0

            ! Bfield
            CALL get_eqdsk_Bspl(raxis_g(i),zaxis_g(k),brtemp,bptemp,bztemp)
            B_R(i,:,k) = B_R(i,:,k) + brtemp
            B_PHI(i,:,k) = B_PHI(i,:,k) + bptemp
            B_Z(i,:,k) = B_Z(i,:,k) + bztemp

            ! ! Flux
            ! CALL get_eqdsk_fluxspl(raxis_g(i),zaxis_g(k),rhoflx,uflx)
            ! sflx = rhoflx*rhoflx
            ! IF (zaxis_g(k)>zsmax .or. zaxis_g(k)<zsmin .or. &
            !     raxis_g(i)>rsmax .or. raxis_g(i)<rsmin) THEN
            ! IF (sflx<1) sflx = 2-sflx ! Handle flux issue
            ! END IF
            ! S_ARR(i,:,k)=sflx
            ! IF (uflx<0)  uflx = uflx+pi2
            ! U_ARR(i,:,k)=uflx
            ! IF (MOD(s,nr) == 0) THEN
            ! IF (lverb) THEN
            !     CALL backspace_out(6,6)
            !     WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*s)/(myend-mystart+1)),']%'
            !     CALL FLUSH(6)
            ! END IF
            ! END IF
        END DO
        
#if defined(MPI_OPT)
        CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
#endif
        

!       ! Calculate Gfactor for mgrid
!       !   Still needs to be done for EQDSK
!       IF (myworkid == master) THEN
!          ! Only master has Gfactor
!          ALLOCATE(Gfactor(ns_prof1))
!          Gfactor = 1
! !         DO s = 1, ns_prof1
! !            sflx = REAL(s-1 + 0.5)/ns_prof1 ! Half grid rho
! !            sflx = sflx*sflx
! !            CALL EZspline_interp(ZEFF_spl_s,sflx,uflx,ier)
!             ! sflx=s uflx=Zeff
! !            CALL vmec_ohkawa(sflx,uflx,Gfactor(s)) ! (1-l31)/Zeff
! !         END DO
!       ENDIF

        ! Deallocations
        CALL read_eqdsk_deallocate
        
        IF (lverb) THEN
            CALL backspace_out(6,36)
            CALL FLUSH(6)
            WRITE(6,'(36X)',ADVANCE='no')
            CALL FLUSH(6)
            CALL backspace_out(6,36)
            WRITE(6,*)
            CALL FLUSH(6)
        END IF    

#if defined(MPI_OPT)
        CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
        CALL MPI_COMM_FREE(MPI_COMM_LOCAL,ierr_mpi)
        CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
        IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init_eqdsk',ierr_mpi)
#endif

        RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
        END SUBROUTINE fieldlines_init_eqdsk
