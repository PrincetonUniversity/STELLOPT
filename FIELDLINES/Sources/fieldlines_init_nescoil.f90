!-----------------------------------------------------------------------
!     Subroutine:    fieldlines_init_nescoil
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          04/28/2016
!     Description:   This subroutine reads the NESCOIL output file
!                    and calculates the field from that file.
!-----------------------------------------------------------------------
      SUBROUTINE fieldlines_init_nescoil
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE safe_open_mod
      USE fieldlines_runtime
      USE fieldlines_grid, ONLY: raxis_g => raxis, phiaxis, &
                                 zaxis_g => zaxis, nr, nphi, nz, &
                                 rmin, rmax, zmin, zmax, phimin, &
                                 phimax, B_R, B_Z, B_PHI,&
                                 BR_spl, BZ_spl
      USE fieldlines_lines, ONLY: nlines
      USE mpi_params
      USE mpi_inc
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, PARAMETER :: BYTE_8 = SELECTED_INT_KIND (8)
!DEC$ IF DEFINED (MPI_OPT)
      INTEGER :: sender, status(MPI_STATUS_SIZE)                     !mpi stuff
      INTEGER(KIND=BYTE_8),ALLOCATABLE :: mnum(:), moffsets(:)
!DEC$ ENDIF  
      INTEGER(KIND=BYTE_8) :: icount, chunk
      INTEGER :: s, i, j, k, nu, nv, mystart,myend

      INTEGER :: iunit, ier, mn_surf, npos, temp_int, mn_pot, &
                 nuv, n, iper, nfp, ibex_nes
      INTEGER, ALLOCATABLE :: xm_surf(:), xn_surf(:), xm_pot(:), xn_pot(:)
      REAL(rprec) :: temp_real, alp, fnuv, x, y, z, bx, by, bz, &
                     bx_temp, by_temp, rmax_nes, rmin_nes, zmax_nes, &
                     zmin_nes, cut_nes, cup_nes
      REAL(rprec), ALLOCATABLE :: rmnc_surf(:), zmns_surf(:), pmns(:),&
                                  x_surf(:), y_surf(:), z_surf(:), &
                                  xu_surf(:),yu_surf(:),zu_surf(:),&
                                  xv_surf(:),yv_surf(:),zv_surf(:),&
                                  snx(:),sny(:),snz(:), &
                                  pot_surf(:), xcur(:), ycur(:), zcur(:),&
                                  dbx(:), dby(:), dbz(:), dn(:), sq2(:),&
                                  sq3(:), dx(:), dy(:), dz(:)
      DOUBLE PRECISION, ALLOCATABLE :: temp_3D(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: xu(:), xv(:)
      DOUBLE PRECISION, ALLOCATABLE :: fmn_temp(:,:)
      CHARACTER(256) :: temp_str
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Set some defaults
      nu = 128; nv = 128; nuv = nu * nv; fnuv = REAL(1)/REAL(nuv)

      ! Read the NESCOIL FILE
      IF (myworkid == master) THEN
          npos = 1
          CALL safe_open(iunit,ier,'nescout.'//TRIM(nescoil_string),'unknown','formatted')
          DO
          READ(iunit,'(A256)') temp_str
          npos = npos + 1
          IF (temp_str(1:40) == '----- Plasma information from VMEC -----') THEN
               READ(iunit,'(A256)') temp_str ! 'np, iota_edge, phip_edge, curpol'
               npos = npos + 1
               READ(iunit,*)  nfp, temp_real, temp_real, temp_real
               npos = npos + 1
          END IF
          IF (temp_str(1:28) == '----- Current Controls -----') THEN
               READ(iunit,'(A256)') temp_str ! 'cut, cup, ibex'
               npos = npos + 1
               READ(iunit,*)  cut_nes, cup_nes, ibex_nes
               npos = npos + 1
          END IF
          IF (temp_str(1:24) == '----- Coil Surface -----') THEN
               READ(iunit,'(A256)') temp_str ! 'Number of fourier modes in table'
               npos = npos + 1
               READ(iunit,*) mn_surf
               npos = npos + 1
               ALLOCATE(xm_surf(mn_surf),xn_surf(mn_surf),rmnc_surf(mn_surf),zmns_surf(mn_surf))
               READ(iunit,'(A256)') temp_str ! '----- Coil surface fourier coefficients -----'
               npos = npos + 1
               READ(iunit,'(A256)') temp_str ! '    m    n         R(m,n)         Z(m,n)'
               npos = npos + 1
               DO i = 1, mn_surf
                    READ(iunit,*)  xm_surf(i), xn_surf(i), rmnc_surf(i),zmns_surf(i), temp_real, temp_real
                    npos = npos + 1
               END DO
          END IF
          IF (temp_str(1:35) == '---- Phi(m,n) for least squares ---') THEN
               mn_pot = 0
               DO
                    ier = 0
                    READ(UNIT=iunit,FMT=*,IOSTAT=ier) temp_int, temp_int, temp_real
                    IF (ier /=0) EXIT
                    mn_pot = mn_pot + 1
               END DO
               ALLOCATE(xm_pot(mn_pot),xn_pot(mn_pot),pmns(mn_pot))
               REWIND(iunit)
               DO i = 1, npos
                    READ(iunit,'(A256)') temp_str
               END DO
               DO i = 1, mn_pot
                    ier = 0
                    READ(UNIT=iunit,FMT=*,IOSTAT=ier) xm_pot(i), xn_pot(i), pmns(i)
               END DO
               EXIT
          END IF
          END DO
          CLOSE(iunit)
     
          ! Transform to real space
          alp  = pi2/nfp
          ALLOCATE(xu(nu),xv(nv))
          FORALL(i = 1:nu) xu(i) = REAL(i-1)/REAL(nu)
          FORALL(i = 1:nv) xv(i) = REAL(i-1)/REAL(nv)
     
          ! Transform Boundary (Note NESCOIL convention is mu+nv)
          ALLOCATE(temp_3D(1,nu,nv),fmn_temp(mn_surf,1))
          ALLOCATE(x_surf(nuv), y_surf(nuv), z_surf(nuv), pot_surf(nuv))
          FORALL(i = 1:mn_surf) fmn_temp(i,1) = rmnc_surf(i)
          STOP "NEED TO FIX FIELDLINES_INIT_NESCOIL TO SUPPORT mntouv_local NOW PRIVATE"
          !CALL mntouv(1,1,mn_surf,nu,nv,xu,xv,fmn_temp,xm_surf,xn_surf,temp_3D,0,1)
          rmax_nes = maxval(maxval(maxval(temp_3D,DIM=3),DIM=2),DIM=1)
          rmin_nes = minval(minval(minval(temp_3D,DIM=3),DIM=2),DIM=1)
          n = 1
          DO i = 1, nu
             DO j = 1, nv
                x_surf(n) = temp_3D(1,i,j)*cos(pi2*xv(j)/nfp)
                y_surf(n) = temp_3D(1,i,j)*sin(pi2*xv(j)/nfp)
                n = n + 1
             END DO
          END DO
          temp_3D = 0
          FORALL(i = 1:mn_surf) fmn_temp(i,1) = zmns_surf(i)
          !CALL mntouv(1,1,mn_surf,nu,nv,xu,xv,fmn_temp,xm_surf,xn_surf,temp_3D,1,0)
          zmax_nes = maxval(maxval(maxval(temp_3D,DIM=3),DIM=2),DIM=1)
          zmin_nes = minval(minval(minval(temp_3D,DIM=3),DIM=2),DIM=1)
          n = 1
          DO i = 1, nu
             DO j = 1, nv
                z_surf(n) = temp_3D(1,i,j)
                n = n + 1
             END DO
          END DO
     
          ! Calculate Surface normals
          ALLOCATE(xu_surf(nuv), yu_surf(nuv))
          FORALL(i = 1:mn_surf) fmn_temp(i,1) = -rmnc_surf(i)*xm_surf(i)
          !CALL mntouv(1,1,mn_surf,nu,nv,xu,xv,fmn_temp,xm_surf,xn_surf,temp_3D,1,0)
          n = 1
          DO i = 1, nu
             DO j = 1, nv
                xu_surf(n) = temp_3D(1,i,j)*cos(pi2*xv(j))
                yu_surf(n) = temp_3D(1,i,j)*sin(pi2*xv(j))
                n = n + 1
             END DO
          END DO
          ALLOCATE(xv_surf(nuv), yv_surf(nuv))
          FORALL(i = 1:mn_surf) fmn_temp(i,1) = -rmnc_surf(i)*xn_surf(i)
          !CALL mntouv(1,1,mn_surf,nu,nv,xu,xv,fmn_temp,xm_surf,xn_surf,temp_3D,1,0)
          n = 1
          DO i = 1, nu
             DO j = 1, nv
                xv_surf(n) = temp_3D(1,i,j)*cos(pi2*xv(j))
                yv_surf(n) = temp_3D(1,i,j)*sin(pi2*xv(j))
                n = n + 1
             END DO
          END DO
          xv_surf = xv_surf - y_surf
          yv_surf = yv_surf + x_surf
          ALLOCATE(zu_surf(nuv))
          FORALL(i = 1:mn_surf) fmn_temp(i,1) = zmns_surf(i)*xm_surf(i)
          !CALL mntouv(1,1,mn_surf,nu,nv,xu,xv,fmn_temp,xm_surf,xn_surf,temp_3D,0,0)
          n = 1
          DO i = 1, nu
             DO j = 1, nv
                zu_surf(n) = temp_3D(1,i,j)
                n = n + 1
             END DO
          END DO
          ALLOCATE(zv_surf(nuv))
          FORALL(i = 1:mn_surf) fmn_temp(i,1) = zmns_surf(i)*xn_surf(i)
          !CALL mntouv(1,1,mn_surf,nu,nv,xu,xv,fmn_temp,xm_surf,xn_surf,temp_3D,0,0)
          n = 1
          DO i = 1, nu
             DO j = 1, nv
                zv_surf(n) = temp_3D(1,i,j)
                n = n + 1
             END DO
          END DO
          DEALLOCATE(rmnc_surf, zmns_surf, xm_surf, xn_surf)
          DEALLOCATE(temp_3D, fmn_temp)
          ALLOCATE(snx(nuv),sny(nuv),snz(nuv),xcur(nuv),ycur(nuv),zcur(nuv))
          snx = yu_surf * zv_surf - zu_surf * yv_surf
          sny = zu_surf * xv_surf - xu_surf * zv_surf
          snz = xu_surf * yv_surf - yu_surf * xv_surf
          xcur = cut_nes*xv_surf-cup_nes*xu_surf
          ycur = cut_nes*yv_surf-cup_nes*yu_surf
          zcur = cut_nes*zv_surf-cup_nes*zu_surf
          DEALLOCATE(xu_surf, yu_surf, zu_surf, xv_surf, yv_surf, zv_surf)
          ALLOCATE(dbx(nuv),dby(nuv),dbz(nuv),dn(nuv),sq2(nuv),sq3(nuv),dx(nuv),dy(nuv),dz(nuv))
     
          ALLOCATE(temp_3D(1,nu,nv),fmn_temp(mn_pot,1))
          temp_3D = 0
          FORALL(i = 1:mn_pot) fmn_temp(i,1) = pmns(i)
          !CALL mntouv(1,1,mn_pot,nu,nv,xu,xv,fmn_temp,xm_pot,xn_pot,temp_3D,1,1)
          n = 1
          DO i = 1, nu
             DO j = 1, nv
                pot_surf(n) = temp_3D(1,i,j)
                n = n + 1
             END DO
          END DO
          DEALLOCATE(pmns, xn_pot, xm_pot)
          DEALLOCATE(temp_3D, fmn_temp, xu, xv)
      
      ELSE
         ALLOCATE(x_surf(nuv), y_surf(nuv), z_surf(nuv), pot_surf(nuv), &
                  snx(nuv), sny(nuv), snz(nuv), xcur(nuv), ycur(nuv), zcur(nuv))
         ALLOCATE(dbx(nuv),dby(nuv),dbz(nuv),dn(nuv),sq2(nuv),sq3(nuv),dx(nuv),dy(nuv),dz(nuv))
      END IF
      
      !


!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_nescoil 1',ierr_mpi)
      CALL MPI_BCAST(nfp,1,MPI_INTEGER, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_nescoil 2',ierr_mpi)
      CALL MPI_BCAST(alp,1,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_nescoil 3',ierr_mpi)
      CALL MPI_BCAST(x_surf,nuv,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_nescoil 5',ierr_mpi)
      CALL MPI_BCAST(y_surf,nuv,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_nescoil 6',ierr_mpi)
      CALL MPI_BCAST(z_surf,nuv,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_nescoil 7',ierr_mpi)
      CALL MPI_BCAST(pot_surf,nuv,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_nescoil 8',ierr_mpi)
      CALL MPI_BCAST(snx,nuv,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_nescoil 9',ierr_mpi)
      CALL MPI_BCAST(sny,nuv,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_nescoil 10',ierr_mpi)
      CALL MPI_BCAST(snz,nuv,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_nescoil 11',ierr_mpi)
      CALL MPI_BCAST(xcur,nuv,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_nescoil 12',ierr_mpi)
      CALL MPI_BCAST(ycur,nuv,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_nescoil 13',ierr_mpi)
      CALL MPI_BCAST(zcur,nuv,MPI_DOUBLE_PRECISION, master, MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_nescoil 14',ierr_mpi)
!DEC$ ENDIF
      
      IF (lverb) THEN
         WRITE(6,'(A)') '----- NESCOIL Current Surface -----'
         WRITE(6,'(A,F8.5,A,F8.5,A,I4)') '   R   = [',rmin_nes,',',rmax_nes,'];  NU:   ',nu
         WRITE(6,'(A,F8.5,A,F8.5,A,I4)') '   Z   = [',zmin_nes,',',zmax_nes,'];  NV:   ',nv
         !IF (.not.lplasma_only) CALL virtual_casing_info(6)
         WRITE(6,'(5X,A,I3.3,A)',ADVANCE='no') 'Vacuum Field Calculation [',0,']%'
         CALL FLUSH(6)
      END IF
      
      ! Break up the Work
      chunk = FLOOR(REAL(nr*nphi*nz) / REAL(nprocs_fieldlines))
      mystart = myworkid*chunk + 1
      myend = mystart + chunk - 1

      ! This section sets up the work so we can use ALLGATHERV
!DEC$ IF DEFINED (MPI_OPT)
      IF (ALLOCATED(mnum)) DEALLOCATE(mnum)
      IF (ALLOCATED(moffsets)) DEALLOCATE(moffsets)
      ALLOCATE(mnum(nprocs_fieldlines), moffsets(nprocs_fieldlines))
      CALL MPI_ALLGATHER(chunk,1,MPI_INTEGER,mnum,1,MPI_INTEGER,MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_ALLGATHER(mystart,1,MPI_INTEGER,moffsets,1,MPI_INTEGER,MPI_COMM_FIELDLINES,ierr_mpi)
      i = 1
      DO
         IF ((moffsets(nprocs_fieldlines)+mnum(nprocs_fieldlines)-1) == nr*nphi*nz) EXIT
         IF (i == nprocs_fieldlines) i = 1
         mnum(i) = mnum(i) + 1
         moffsets(i+1:nprocs_fieldlines) = moffsets(i+1:nprocs_fieldlines) + 1
         i=i+1
      END DO
      mystart = moffsets(myworkid+1)
      chunk  = mnum(myworkid+1)
      myend   = mystart + chunk - 1
!DEC$ ENDIF
      IF (lafield_only) THEN
            ! This is not supported
      ELSE
         DO s = mystart, myend
            i = MOD(s-1,nr)+1
            j = MOD(s-1,nr*nphi)
            j = FLOOR(REAL(j) / REAL(nr))+1
            k = CEILING(REAL(s) / REAL(nr*nphi))
            x = raxis_g(i)*COS(phiaxis(j))
            y = raxis_g(i)*SIN(phiaxis(j))
            z = zaxis_g(k)
            bx = 0; by = 0; bz = 0;
            DO iper = 0, nfp-1  ! Take advantage of symmetry
               dx = x_surf * COS(alp * iper) + y_surf * SIN(alp * iper) - x
               dy = y_surf * COS(alp * iper) - x_surf * SIN(alp * iper) - y
               dz = z_surf - z
               sq2 = 1/(dx*dx + dy*dy + dz*dz)
               sq3 = sq2*SQRT(sq2)
               dn  = 3*sq2*(snx*dx + sny*dy + snz*dz)*pot_surf
               dbx = (snx*pot_surf - dn*dx + ycur*dz - zcur*dy)*sq3
               dby = (sny*pot_surf - dn*dy + zcur*dx - xcur*dz)*sq3
               dbz = (snz*pot_surf - dn*dz + xcur*dy - ycur*dx)*sq3
               bx_temp = SUM(dbx)*fnuv
               by_temp = SUM(dby)*fnuv
               bx = bx + (bx_temp * COS(alp * iper) - by_temp * SIN(alp * iper))
               by = by + (by_temp * COS(alp * iper) + bx_temp * SIN(alp * iper))
               !IF (myworkid == 10) THEN; WRITE(6,*) fnuv; CALL FLUSH(6); STOP; END IF
               B_Z(i,j,k) = B_Z(i,j,k) + SUM(dbz)*fnuv
            END DO
            B_R(i,j,k) = B_R(i,j,k) + bx * COS(phiaxis(j)) + by * sin(phiaxis(j))
            B_PHI(i,j,k) = B_PHI(i,j,k) + by * COS(phiaxis(j)) - bx * sin(phiaxis(j))
            IF (lverb .and. (MOD(s,nr) == 0)) THEN
               CALL backspace_out(6,6)
               WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*s)/(myend-mystart+1)),']%'
               CALL FLUSH(6)
            END IF
         END DO
      END IF
      
      IF (lverb) THEN
         CALL backspace_out(6,36)
         CALL FLUSH(6)
         WRITE(6,'(36X)',ADVANCE='no')
         CALL FLUSH(6)
         CALL backspace_out(6,36)
         CALL FLUSH(6)
      END IF    
      
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'torlines_init_external',ierr_mpi)
!       ! Adjust indexing to send 2D arrays
       CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
                        B_R,mnum,moffsets-1,MPI_DOUBLE_PRECISION,&
                        MPI_COMM_FIELDLINES,ierr_mpi)
       CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
                        B_PHI,mnum,moffsets-1,MPI_DOUBLE_PRECISION,&
                        MPI_COMM_FIELDLINES,ierr_mpi)
       CALL MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&
                        B_Z,mnum,moffsets-1,MPI_DOUBLE_PRECISION,&
                        MPI_COMM_FIELDLINES,ierr_mpi)
       DEALLOCATE(mnum)
       DEALLOCATE(moffsets)
!DEC$ ENDIF
      WHERE(B_PHI==0) B_PHI=1

!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init_vmec',ierr_mpi)
!DEC$ ENDIF
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE fieldlines_init_nescoil
