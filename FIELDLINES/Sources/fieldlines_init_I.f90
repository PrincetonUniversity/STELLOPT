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
                             zaxis_cs, ntor, nfp, read_indata_namelist
      USE mpi_params                       
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      INCLUDE 'mpif.h'   ! MPI - Don't need because of fieldlines_runtime
      INTEGER :: sender
      INTEGER :: status(MPI_STATUS_SIZE)
      REAL(rprec), ALLOCATABLE :: buffer_mast(:,:,:),buffer_slav(:,:,:)
!DEC$ ENDIF  
      INTEGER :: ier, iunit, i, j, myzs, myze, ik, chunk,&
                 count, nv, mn
      REAL(rprec)  :: fac, xp, yp, zp, ax, ay, az, bx, by, bz, br, bphi
      REAL(rprec), ALLOCATABLE :: xv(:), raxis(:), xaxis(:), yaxis(:),& 
                                  zaxis(:), dx(:), dy(:), dz(:), vx(:),&
                                  vy(:), vz(:), &
                                  x1(:), y1(:), z1(:), fa(:), r(:),&
                                  r12(:)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Reread INDATA
      iunit = 11
      OPEN(UNIT=iunit, FILE='input.' // TRIM(id_string), STATUS='OLD', IOSTAT=ier)
      IF (ier /= 0) CALL handle_err(FILE_OPEN_ERR,id_string,ier)
      CALL read_indata_namelist(iunit,ier)
      IF (ier /= 0) CALL handle_err(VMEC_INPUT_ERR,id_string,ier)
      CLOSE(iunit)
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
      
      
      ! Get the fields
      chunk  = nz/numprocs
      IF (MOD(nz,numprocs) /= 0) chunk = chunk + 1
      myzs = myid*chunk + 1
      myze = myzs + chunk - 1
      if (myze > nz) myze = nz
      count = 0
      IF (myzs <= nz) THEN
         DO ik = myzs, myze
            DO j = 1, nphi
               DO i = 1, nr
                  xp  = raxis_g(i)*cos(phiaxis(j))
                  yp  = raxis_g(i)*sin(phiaxis(j))
                  zp  = zaxis_g(ik)
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
                  B_R(i,j,ik) = B_R(i,j,ik) + br
                  B_PHI(i,j,ik) = B_PHI(i,j,ik) + bphi
                  B_Z(i,j,ik) = B_Z(i,j,ik) + bz
                  count = count + 1
               END DO
               IF (lverb) THEN
                  CALL backspace_out(6,6)
                  WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*count)/(nr*nphi*myze)),']%'
                  CALL FLUSH(6)
               END IF
            END DO
         END DO
      END IF
      
      ! Clean up the progress bar
      IF (lverb) THEN
         CALL backspace_out(6,38)
         WRITE(6,'(38X)',ADVANCE='no')
         CALL backspace_out(6,38)
         CALL FLUSH(6)
      END IF    
      
      ! Now recompose the array and send to everyone.
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init_I',ierr_mpi)
      IF (myid == master) THEN
         ALLOCATE(buffer_mast(nr,nphi,3))
         DO i = myze+1, nz
            CALL MPI_RECV(buffer_mast,nr*nphi*3,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_FIELDLINES,status,ierr_mpi)
            IF (ierr_mpi /=0) CALL handle_err(MPI_RECV_ERR,'fieldlines_init_I',ierr_mpi)
            sender = status(MPI_SOURCE)
            j      = status(MPI_TAG)
            B_R(:,:,j) = buffer_mast(:,:,1)
            B_PHI(:,:,j) = buffer_mast(:,:,2)
            B_Z(:,:,j) = buffer_mast(:,:,3)
         END DO
         DEALLOCATE(buffer_mast)
      ELSE
         IF (myzs <= nz) THEN
            ALLOCATE(buffer_slav(nr,nphi,3))
            DO j = myzs, myze
               buffer_slav(:,:,1) = B_R(:,:,j)
               buffer_slav(:,:,2) = B_PHI(:,:,j)
               buffer_slav(:,:,3) = B_Z(:,:,j)
               CALL MPI_SEND(buffer_slav,nr*nphi*3,MPI_DOUBLE_PRECISION,master,j,MPI_COMM_FIELDLINES,ierr_mpi)
               IF (ierr_mpi /=0) CALL handle_err(MPI_SEND_ERR,'fieldlines_init_I',ierr_mpi)
            END DO
            DEALLOCATE(buffer_slav)
         END IF
      END IF
      
      ! Have master send info to slaves, in parts
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init_I',ierr_mpi)
      IF (myid == master) THEN
         ALLOCATE(buffer_mast(nr,nphi,3))
         DO i = 1, nz
            buffer_mast(:,:,1) = B_R(:,:,i)
            buffer_mast(:,:,2) = B_PHI(:,:,i)
            buffer_mast(:,:,3) = B_Z(:,:,i)
            DO j = 1, numprocs - 1
               CALL MPI_SEND(buffer_mast,nr*nphi*3,MPI_DOUBLE_PRECISION,j,i,MPI_COMM_FIELDLINES,ierr_mpi)
               IF (ierr_mpi /=0) CALL handle_err(MPI_SEND_ERR,'fieldlines_init_I',ierr_mpi)
            END DO
         END DO
         DEALLOCATE(buffer_mast)
      ELSE
         ALLOCATE(buffer_slav(nr,nphi,3))
         DO i = 1, nz
            CALL MPI_RECV(buffer_slav,nr*nphi*3,MPI_DOUBLE_PRECISION,master,MPI_ANY_TAG,MPI_COMM_FIELDLINES,status,ierr_mpi)
            IF (ierr_mpi /=0) CALL handle_err(MPI_RECV_ERR,'fieldlines_init_I',ierr_mpi)
            j      = status(MPI_TAG)
            B_R(:,:,j) = buffer_slav(:,:,1)
            B_PHI(:,:,j) = buffer_slav(:,:,2)
            B_Z(:,:,j) = buffer_slav(:,:,3)
         END DO
         DEALLOCATE(buffer_slav)
      END IF
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      !Old dumb way
      !IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init_mgrid',ierr_mpi)
      !CALL MPI_BCAST(B_R,nr*nphi*nz,MPI_DOUBLE_PRECISION,master,MPI_COMM_FIELDLINES,ierr_mpi)
      !IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_mgrid',ierr_mpi)
      !CALL MPI_BCAST(B_Z,nr*nphi*nz,MPI_DOUBLE_PRECISION,master,MPI_COMM_FIELDLINES,ierr_mpi)
      !IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_mgrid',ierr_mpi)
!DEC$ ENDIF
      
      ! Free Variables
      DEALLOCATE(xv, xaxis, yaxis, zaxis, raxis)
      DEALlOCATE(dx, dy, dz, vx, vy, vz)
      DEALLOCATE(x1, y1, z1, r, r12, fa)
      
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init_I',ierr_mpi)
!DEC$ ENDIF
      
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE fieldlines_init_I
