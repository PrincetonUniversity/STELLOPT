!-----------------------------------------------------------------------
!     Subroutine:    torlines_init_external
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          04/30/2013
!-----------------------------------------------------------------------
      SUBROUTINE torlines_init_external
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE read_wout_mod, nfp_in => nfp, lasym_in => lasym, &
                      lfreeb_in => lfreeb, mnmax_in => mnmax, &
                      rmnc_in => rmnc, zmns_in => zmns, &
                      rmns_in => rmns, zmnc_in => zmnc, &
                      bsupumnc_in => bsupumnc, bsupumns_in => bsupumns,&
                      bsupvmnc_in => bsupvmnc, bsupvmns_in => bsupvmns,&
                      xn_in => xn, xm_in => xm, mgrid_file_in => mgrid_file,&
                      nextcur_in => nextcur, extcur_in => extcur, phi_in => phi
      USE torlines_background
      USE torlines_realspace
      USE torlines_runtime
      USE torlines_fieldlines
      USE virtual_casing_mod, pi2_vc => pi2
      USE mgrid_field_mod, pi2_mgrid => pi2
      USE biotsavart
      USE EZspline_obj
      USE EZspline
      USE mpi_params 
!-----------------------------------------------------------------------
!     Local Variables
!          iunit       File Unit Number
!          ierr        Error flag
!          mn          Dummy index over modes
!          im          Dummy index over modes
!          in          Dummy index over modes
!          uv          Dummy index over real-space
!          ik          Dummy index over radial surfaces
!          vsurf       VMEC surface in PIES background coordinates
!          rho_in      VMEC Rho
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, PARAMETER :: BYTE_8 = SELECTED_INT_KIND (8)
!DEC$ IF DEFINED (MPI_OPT)
!      INCLUDE 'mpif.h'   ! MPI - Don't need because of fieldlines_runtime
      INTEGER :: sender
      INTEGER :: status(MPI_STATUS_SIZE)
      INTEGER(KIND=BYTE_8),ALLOCATABLE :: mnum(:), moffsets(:)
      !REAL(rprec), ALLOCATABLE :: buffer_mast(:,:,:),buffer_slav(:,:,:)
!DEC$ ENDIF  
      INTEGER(KIND=BYTE_8) :: icount, chunk
      INTEGER :: ier, s, u, v, i, j, iunit, s1, nv_in
      REAL(rprec) :: x_vc, y_vc, z_vc, bx_vc, by_vc, bz_vc, br_vc, bphi_vc
      REAL(rprec) :: r, phi, z, br_temp, bphi_temp, bz_temp, br, bphi, bz,&
                    current, current_first, rmin, rmax, zmin, zmax
      LOGICAL :: lcoil_open
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Handle lvac
      ALLOCATE(brreal(1:k,1:nu,1:nv),bphireal(1:k,1:nu,1:nv),bzreal(1:k,1:nu,1:nv))
      s1 = vsurf+1
      IF (lvac .and. (lcoil .or. lmgrid)) s1 = 1
      brreal   = 0.0
      bphireal = 1.0
      bzreal   = 0.0
      
      IF (lverb) THEN
         WRITE(6,'(A)')   '----- Vacuum Grid Info. -----'
         WRITE(6,'(A,I6)')'       Eq. Surface: ',vsurf
         WRITE(6,'(A,I8)')'       Grid Points: ',nu*nv*(k-s1+1)
      END IF

      ! Break up the Work
      chunk = FLOOR(REAL(k*nu*nv) / REAL(numprocs))
      mystart = myid*chunk + 1
      myend = mystart + chunk - 1

      ! This section sets up the work so we can use ALLGATHERV
!DEC$ IF DEFINED (MPI_OPT)
      IF (ALLOCATED(mnum)) DEALLOCATE(mnum)
      IF (ALLOCATED(moffsets)) DEALLOCATE(moffsets)
      ALLOCATE(mnum(numprocs), moffsets(numprocs))
      CALL MPI_ALLGATHER(chunk,1,MPI_INTEGER,mnum,1,MPI_INTEGER,MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_ALLGATHER(mystart,1,MPI_INTEGER,moffsets,1,MPI_INTEGER,MPI_COMM_FIELDLINES,ierr_mpi)
      i = 1
      DO
         IF ((moffsets(numprocs)+mnum(numprocs)-1) == k*nu*nv) EXIT
         IF (i == numprocs) i = 1
         mnum(i) = mnum(i) + 1
         moffsets(i+1:numprocs) = moffsets(i+1:numprocs) + 1
         i=i+1
      END DO
      mystart = moffsets(myid+1)
      chunk  = mnum(myid+1)
      myend   = mystart + chunk - 1
!DEC$ ENDIF

      ! Get vacuum fields
      IF (lcoil) THEN
         ! Initalize the coil
         CALL cleanup_biotsavart
         INQUIRE(FILE=TRIM(coil_string),NUMBER=iunit,OPENED=lcoil_open)
         IF (lcoil_open) CLOSE(iunit)
         CALL parse_coils_file(TRIM(coil_string))
         nextcur = SIZE(coil_group) !SAL
         DO i = 1, nextcur
            DO j = 1, coil_group(i) % ncoil
               current = coil_group(i) % coils(j) % current
               IF (j .eq. 1) current_first = current
               IF (lraw) THEN
                  coil_group(i) % coils(j) % current = current*extcur_in(i)
               ELSE
                  IF (current_first .ne. zero) coil_group(i) % coils(j) % current = (current/current_first)*extcur_in(i)
               END IF
            END DO
         END DO
         IF (lverb) THEN
            WRITE(6,'(A)')   '----- COILS Information -----'
            WRITE(6,'(A,A)') '   FILE: ',TRIM(coil_string)
            WRITE(6,'(A,I3)')'   Coil Periodicity: ',nfp_bs
            WRITE(6,'(A,I3)')'   Current Systems: ',nextcur
            IF (lraw) THEN
               WRITE(6,'(A)')   '   Current Type:      RAW'
            ELSE
               WRITE(6,'(A)')   '   Current Type:      SCALED'
            END IF
            DO i = 1, nextcur
               IF (ABS(extcur(i)).ge. 1.0E9) THEN
                  WRITE(6,'(A,I4,A,F8.3,A)')'   Num Coils  = ',coil_group(i)%ncoil,'  EXTCUR = ',extcur(i)/1.0E6,' [GA]'
               ELSE IF (ABS(extcur(i)).ge. 1.0E6) THEN
                  WRITE(6,'(A,I4,A,F8.3,A)')'   Num Coils  = ',coil_group(i)%ncoil,'  EXTCUR = ',extcur(i)/1.0E6,' [MA]'
               ELSE IF (ABS(extcur(i)).ge. 1.0E3) THEN
                  WRITE(6,'(A,I4,A,F8.3,A)')'   Num Coils  = ',coil_group(i)%ncoil,'  EXTCUR = ',extcur(i)/1.0E3,' [kA]'
               ELSE
                  WRITE(6,'(A,I4,A,F8.3,A)')'   Num Coils  = ',coil_group(i)%ncoil,'  EXTCUR = ',extcur(i),' [A]'
               END IF
            END DO
            WRITE(6,'(5X,A,I3.3,A)',ADVANCE='no') 'Vacuum Field Calculation [',0,']%'
            CALL FLUSH(6)
         END IF
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init',ierr_mpi)
!DEC$ ENDIF
         DO i = mystart, myend
            s = MOD(i-1,k) + 1
            u = MOD(i-1,k*nu)
            u = FLOOR(REAL(u) / REAL(k))+1
            v = CEILING(REAL(i) / REAL(k*nu))
            IF (s < s1) CYCLE
            !u = MOD(s,k)
            !IF (u < s1 .and. u /= 0) CYCLE
            br = 0; bphi = 0; bz = 0
            br_temp   = 0
            bphi_temp = 0
            bz_temp   = 0
            r = rreal(s,u,v)
            phi = pi2*xv(v)/nfp
            z = zreal(s,u,v)
            DO j = 1, nextcur
               CALL bfield(r, phi, z, br, bphi, bz, IG = j)
               br_temp = br_temp + br
               bphi_temp = bphi_temp + bphi
               bz_temp = bz_temp + bz
            END DO
            brreal(s,u,v) = br_temp
            bphireal(s,u,v) = bphi_temp
            bzreal(s,u,v) = bz_temp
            breal(s,u,v)  = sqrt(br_temp*br_temp+bphi_temp*bphi_temp+bz_temp*bz_temp)
            IF (lverb .and. (MOD(i,k) == 0)) THEN
               CALL backspace_out(6,6)
               WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*i)/(myend-mystart+1)),']%'
               CALL FLUSH(6)
            END IF
         END DO
      ELSEIF (lmgrid) THEN
         nv_in = nv
         nfp_in = nfp
         ier = 0
         CALL mgrid_load(mgrid_string,extcur,nextcur,nv_in,nfp_in,ier,myid)
         rmin = MINVAL(rreal)
         rmax = MAXVAL(rreal)
         zmin = MINVAL(zreal)
         zmax = MAXVAL(zreal)
         IF ((rmin < rminb) .or. (rmax > rmaxb) .or. &
             (zmin < zminb) .or. (zmax > zmaxb)) THEN
            IF (lverb) THEN
               WRITE(6,'(A)') '!!!!!!!!!!!! WARNING !!!!!!!!!!!!'
               WRITE(6,'(A)') '!!  Desired Grid larger than   !!'
               WRITE(6,'(A)') '!!  mgrid extent.  Lower grid  !!'
               WRITE(6,'(A)') '!!  dimensions or use coils    !!'
               WRITE(6,'(A)') '!!  get vacuum field.          !!'
               WRITE(6,'(A)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            END IF
            stop
         END IF
         IF (lverb) THEN
            CALL mgrid_info(6)
            WRITE(6,'(5X,A,I3.3,A)',ADVANCE='no') 'Vacuum Field Calculation [',0,']%'
            CALL FLUSH(6)
         END IF
         DO i = mystart, myend
            s = MOD(i-1,k) + 1
            u = MOD(i-1,k*nu)
            u = FLOOR(REAL(u) / REAL(k))+1
            v = CEILING(REAL(i) / REAL(k*nu))
            IF (s < s1) CYCLE
            br = 0; bphi = 0; bz = 0
            r = rreal(s,u,v)
            phi = pi2*xv(v)/nfp
            z = zreal(s,u,v)
            ier = 0
            CALL mgrid_bcyl(r, phi, z, br, bphi, bz, ier)
            brreal(s,u,v) = br
            bphireal(s,u,v) = bphi
            bzreal(s,u,v) = bz
            breal(s,u,v)  = sqrt(br*br+bphi*bphi+bz*bz)
            IF (lverb .and. (MOD(i,k) == 0)) THEN
               CALL backspace_out(6,6)
               WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*i)/(myend-mystart+1)),']%'
               CALL FLUSH(6)
            END IF
         END DO
         CALL mgrid_free(ier)
      END IF
      
      
      ! Now we need the fields on the background grid
      IF (.not.lvac) THEN
         IF (lverb) THEN
            WRITE(6,*) ''
            CALL virtual_casing_info(6)
            WRITE(6,'(5X,A,I3.3,A)',ADVANCE='no') 'Plasma Field Calculation [',0,']%'
            CALL FLUSH(6)
         END IF
         DO i = mystart, myend
            s = MOD(i-1,k) + 1
            u = MOD(i-1,k*nu)
            u = FLOOR(REAL(u) / REAL(k))+1
            v = CEILING(REAL(i) / REAL(k*nu))
            IF (s < s1) CYCLE
            phi = pi2*xv(v)/nfp
            x_vc = rreal(s,u,v)*cos(phi)
            y_vc = rreal(s,u,v)*sin(phi)
            z_vc = zreal(s,u,v)
            ier = 1 ! Keep on truckin
            bx_vc = 0.0; by_vc = 0.0; bz_vc = 0.0
            adapt_rel = 0.0
            adapt_tol = vc_adapt_tol*sqrt(brreal(s,u,v)**2+bphireal(s,u,v)**2+bzreal(s,u,v)**2)
            CALL bfield_vc(x_vc,y_vc,z_vc,bx_vc,by_vc,bz_vc,ier)
            IF (ier == 0) THEN
               br_vc   = bx_vc * cos(phi) + by_vc * sin(phi)
               bphi_vc = by_vc * cos(phi) - bx_vc * sin(phi)
               brreal(s,u,v) = brreal(s,u,v) + br_vc
               bphireal(s,u,v) = bphireal(s,u,v) + bphi_vc
               bzreal(s,u,v) = bzreal(s,u,v) + bz_vc
               breal(s,u,v)  = sqrt(brreal(s,u,v)**2+bphireal(s,u,v)*2+bzreal(s,u,v)**2)
            END IF
            IF (lverb .and. (MOD(i,k) == 0)) THEN
               CALL backspace_out(6,6)
               WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*i)/(myend-mystart+1)),']%'
               CALL FLUSH(6)
            END IF
         END DO
         CALL free_virtual_casing
      END IF

      CALL read_wout_deallocate


!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'torlines_init_external',ierr_mpi)
!       ! Adjust indexing to send 2D arrays
       CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
                        brreal,mnum,moffsets-1,MPI_DOUBLE_PRECISION,&
                        MPI_COMM_FIELDLINES,ierr_mpi)
       CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
                        bphireal,mnum,moffsets-1,MPI_DOUBLE_PRECISION,&
                        MPI_COMM_FIELDLINES,ierr_mpi)
       CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
                        bzreal,mnum,moffsets-1,MPI_DOUBLE_PRECISION,&
                        MPI_COMM_FIELDLINES,ierr_mpi)
       CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
                        breal,mnum,moffsets-1,MPI_DOUBLE_PRECISION,&
                        MPI_COMM_FIELDLINES,ierr_mpi)
       DEALLOCATE(mnum)
       DEALLOCATE(moffsets)
!DEC$ ENDIF

      ! Now we need to calculate the Bs,Bu,Bv compoments of the field
      z = 1
      CALL cyl2suv(s1,k,nu,nv,rreal(s1:k,:,:),brreal(s1:k,:,:),bphireal(s1:k,:,:),bzreal(s1:k,:,:),bsreal(s1:k,:,:),bureal(s1:k,:,:),bvreal(s1:k,:,:),z)
     
      IF (lvac) THEN
         bsreal(1,:,:) = 0.0
         bureal(1,:,:) = 0.0
      END IF 
      
      
      
      DEALLOCATE(brreal,bphireal,bzreal)
      
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE torlines_init_external
