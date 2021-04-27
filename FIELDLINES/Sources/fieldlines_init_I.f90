!-----------------------------------------------------------------------
!     Module:        fieldlines_init_I
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          04/27/2012
!     Description:   This subroutine places CURTOR current on the
!                    magnetic axis.
!-----------------------------------------------------------------------
      SUBROUTINE fieldlines_init_I
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE fieldlines_runtime
      USE fieldlines_grid, ONLY: raxis_g=>raxis,phiaxis,zaxis_g=>zaxis, nr, nphi, nz, &
                                 rmin, rmax, zmin, zmax, phimin, &
                                 phimax, B_R, B_Z, B_PHI
      USE vmec_input,  ONLY: curtor, raxis_cc, raxis_cs, zaxis_cc, &
                             zaxis_cs, ntor, nfp, read_indata_namelist, ntord
      USE mpi_params
      USE mpi_inc
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, PARAMETER :: BYTE_8 = SELECTED_INT_KIND (8)
      INTEGER(KIND=BYTE_8),ALLOCATABLE :: mnum(:), moffsets(:)
      INTEGER :: numprocs_local, mylocalid, mylocalmaster
      INTEGER :: MPI_COMM_LOCAL
      INTEGER(KIND=BYTE_8) :: chunk
      INTEGER :: ier, iunit, s, i, j, mystart, myend, k, mn, nv
      REAL(rprec)  :: fac, xp, yp, zp, ax, ay, az, bx, by, bz, br, bphi
      REAL(rprec), ALLOCATABLE :: xv(:), raxis(:), xaxis(:), yaxis(:),& 
                                  zaxis(:), dx(:), dy(:), dz(:), vx(:),&
                                  vy(:), vz(:), &
                                  x1(:), y1(:), z1(:), fa(:), r(:),&
                                  r12(:)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------

      ! Divide up Work
      numprocs_local = 1; mylocalid = master
#if defined(MPI_OPT)
      CALL MPI_COMM_DUP( MPI_COMM_SHARMEM, MPI_COMM_LOCAL, ierr_mpi)
      CALL MPI_COMM_RANK( MPI_COMM_LOCAL, mylocalid, ierr_mpi )              ! MPI
      CALL MPI_COMM_SIZE( MPI_COMM_LOCAL, numprocs_local, ierr_mpi )          ! MPI
#endif
      mylocalmaster = master


      ! Reread INDATA
      IF (mylocalid == mylocalmaster) THEN
         iunit = 11
         OPEN(UNIT=iunit, FILE='input.' // TRIM(id_string), STATUS='OLD', IOSTAT=ier)
         IF (ier /= 0) CALL handle_err(FILE_OPEN_ERR,id_string,ier)
         CALL read_indata_namelist(iunit,ier)
         IF (ier /= 0) CALL handle_err(VMEC_INPUT_ERR,id_string,ier)
         CLOSE(iunit)
      ENDIF
#if defined(MPI_OPT)
      CALL MPI_BCAST(nfp,1,MPI_INTEGER, mylocalmaster, MPI_COMM_LOCAL,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_coil',ierr_mpi)
      CALL MPI_BCAST(ntor,1,MPI_INTEGER, mylocalmaster, MPI_COMM_LOCAL,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_coil',ierr_mpi)
      CALL MPI_BCAST(curtor,1,MPI_REAL, mylocalmaster, MPI_COMM_LOCAL,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_coil',ierr_mpi)
      CALL MPI_BCAST(raxis_cc(0:ntord),ntord+1,MPI_REAL, mylocalmaster, MPI_COMM_LOCAL,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_coil',ierr_mpi)
      CALL MPI_BCAST(raxis_cs(0:ntord),ntord+1,MPI_REAL, mylocalmaster, MPI_COMM_LOCAL,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_coil',ierr_mpi)
      CALL MPI_BCAST(zaxis_cc(0:ntord),ntord+1,MPI_REAL, mylocalmaster, MPI_COMM_LOCAL,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_coil',ierr_mpi)
      CALL MPI_BCAST(zaxis_cs(0:ntord),ntord+1,MPI_REAL, mylocalmaster, MPI_COMM_LOCAL,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_coil',ierr_mpi)
#endif

      ! Construct the axis filament
      fac = curtor*1.0E-7
      nv = nphi*nfp+1
      ALLOCATE(xv(1:nv),STAT=ier)
      FORALL(i=1:nv) xv(i) = REAL(i-1)/REAL(nv)
      ALLOCATE(raxis(1:nv+1),zaxis(1:nv+1),STAT=ier)
      ALLOCATE(xaxis(nv+1),yaxis(nv+1),STAT=ier)
      raxis = 0.0; zaxis = 0.0
      DO i = 1, nv
         DO mn = 0, ntor
            raxis(i) = raxis(i) + raxis_cc(mn)*cos(-mn*nfp*pi2*xv(i))
            raxis(i) = raxis(i) + raxis_cs(mn)*sin(-mn*nfp*pi2*xv(i))
            zaxis(i) = zaxis(i) + zaxis_cc(mn)*cos(-mn*nfp*pi2*xv(i))
            zaxis(i) = zaxis(i) + zaxis_cs(mn)*sin(-mn*nfp*pi2*xv(i))
         END DO
         xaxis(i) = raxis(i) * cos(pi2*xv(i))
         yaxis(i) = raxis(i) * sin(pi2*xv(i))
         !IF (lverb) PRINT *,'axis',raxis(i),xaxis(i),yaxis(i),zaxis(i)
      END DO
      raxis(nv+1) = raxis(1)
      zaxis(nv+1) = zaxis(1)
      xaxis(nv+1) = xaxis(1)
      yaxis(nv+1) = yaxis(1)
      

      ! Output the Axis Information
      IF (lverb) THEN
         WRITE(6,'(A)')   '----- Axis Current Information -----'
         WRITE(6,'(A,I3,A,I3,A)')    '       N:  [',-ntor,',',ntor,']'
         WRITE(6,'(A,F8.3,A,F8.3,A)')'   R0,Z0: ',raxis(1),' ',zaxis(1),' [m]'
         IF (ABS(curtor) .ge. 1.0E6) THEN
            WRITE(6,'(A,F8.3,A)')    '      I0: ',curtor/1.0E6,' [MA]'
         ELSE IF (ABS(curtor) .ge. 1.0E3) THEN
            WRITE(6,'(A,F8.3,A)')    '      I0: ',curtor/1.0E3,' [kA]'
         ELSE
            WRITE(6,'(A,F8.3,A)')    '      I0: ',curtor,' [A]'
         END IF
         WRITE(6,'(5X,A,I3.3,A)',ADVANCE='no') 'Axis Current Calculation [',0,']%'
         CALL FLUSH(6)
      END IF
      
      ! Precalculate various arrays
      ALLOCATE(dx(nv),dy(nv),dz(nv),STAT=ier)
      ALLOCATE(vx(nv),vy(nv),vz(nv),STAT=ier)
      ALLOCATE(x1(nv+1),y1(nv+1),z1(nv+1),STAT=ier)
      ALLOCATE(r(nv+1),r12(nv),fa(nv),STAT=ier)
      dx = 0.0; dy = 0.0; dz = 0.0; vx = 0.0; vy = 0.0; vz = 0.0
      dx(1:nv) = (xaxis(2:nv+1)-xaxis(1:nv))
      dy(1:nv) = (yaxis(2:nv+1)-yaxis(1:nv))
      dz(1:nv) = (zaxis(2:nv+1)-zaxis(1:nv))
      vx(1:nv) = yaxis(1:nv)*dz(1:nv) - zaxis(1:nv)*dy(1:nv)
      vy(1:nv) = zaxis(1:nv)*dx(1:nv) - xaxis(1:nv)*dz(1:nv)
      vz(1:nv) = xaxis(1:nv)*dy(1:nv) - yaxis(1:nv)*dx(1:nv)


      ! Break up the Work
      CALL MPI_CALC_MYRANGE(MPI_COMM_LOCAL,1, nr*nphi*nz, mystart, myend)

      ! Put the field on the grid
      DO s = mystart, myend
         i = MOD(s-1,nr)+1
         j = MOD(s-1,nr*nphi)
         j = FLOOR(REAL(j) / REAL(nr))+1
         k = CEILING(REAL(s) / REAL(nr*nphi))
         xp  = raxis_g(i)*cos(phiaxis(j))
         yp  = raxis_g(i)*sin(phiaxis(j))
         zp  = zaxis_g(k)
         x1  = xp - xaxis
         y1  = yp - yaxis
         z1  = zp - zaxis
         r   = SQRT(x1*x1+y1*y1+z1*z1)
         IF (ANY(r .eq. 0)) CYCLE
         r12(1:nv) = r(2:nv+1)*r(2:nv+1)
         fa(1:nv)  = (r(2:nv+1)+r(1:nv))/&
                     (r12(1:nv) * ( r12(1:nv) &
                                   + x1(2:nv+1)*x1(1:nv)&
                                   + y1(2:nv+1)*y1(1:nv)&
                                   + z1(2:nv+1)*z1(1:nv)) )
         ax   = SUM(fa(1:nv)*dx(1:nv))
         ay   = SUM(fa(1:nv)*dy(1:nv))
         az   = SUM(fa(1:nv)*dz(1:nv))
         bx   = fac*(SUM(fa(1:nv)*vx(1:nv) - yp*az+zp*ay))
         by   = fac*(SUM(fa(1:nv)*vy(1:nv) - zp*ax+xp*az))
         bz   = fac*(SUM(fa(1:nv)*vz(1:nv) - xp*ay+yp*ax))
         br   = bx * cos(phiaxis(j)) + by * sin(phiaxis(j))
         bphi = by * cos(phiaxis(j)) - bx * sin(phiaxis(j))
         IF (lverb) PRINT *,'x,y,z,bx,by,bz',xp,yp,zp,bx,by,bz
         B_R(i,j,k) = B_R(i,j,k) + br
         B_PHI(i,j,k) = B_PHI(i,j,k) + bphi
         B_Z(i,j,k) = B_Z(i,j,k) + bz
         IF (lverb .and. (MOD(s,nr) == 0)) THEN
            CALL backspace_out(6,6)
            WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*s)/(myend-mystart+1)),']%'
            CALL FLUSH(6)
         END IF
      END DO
      
      ! Clean up the progress bar
      IF (lverb) THEN
         CALL backspace_out(6,38)
         WRITE(6,'(38X)',ADVANCE='no')
         CALL backspace_out(6,38)
         CALL FLUSH(6)
      END IF    
      
      ! Free Variables
      DEALLOCATE(xv, xaxis, yaxis, zaxis, raxis)
      DEALlOCATE(dx, dy, dz, vx, vy, vz)
      DEALLOCATE(x1, y1, z1, r, r12, fa)

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init_coil1',ierr_mpi)
      CALL MPI_COMM_FREE(MPI_COMM_LOCAL,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'fieldlines_init_coil: MPI_COMM_LOCAL',ierr_mpi)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init_coil',ierr_mpi)
#endif
      
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE fieldlines_init_I
