!-----------------------------------------------------------------------
!     Module:        diagno_bench
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          07/30/2013
!     Description:   This subroutine is utilized for debugging
!-----------------------------------------------------------------------
      SUBROUTINE diagno_bench
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE diagno_runtime, pi2_runtime => pi2
      USE virtual_casing_mod
      USE read_wout_mod, ONLY: nfp, zmax_surf
      USE biotsavart
      USE safe_open_mod
      USE mpi_params
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, PARAMETER :: BYTE_8 = SELECTED_INT_KIND (8)
#if defined(MPI_OPT)
      INCLUDE 'mpif.h'
      INTEGER(KIND=BYTE_8),ALLOCATABLE :: mnum(:), moffsets(:)
      INTEGER :: numprocs_local, mylocalid, mylocalmaster
      INTEGER :: MPI_COMM_LOCAL
#endif
      INTEGER :: ier, iunit, ncoils, i, ig, iunit_out, ier1
      REAL(rprec) :: xp, yp, rp, phip, zp, bx, by, br, bphi, bz, modb,&
                     bxp, byp, bzp, c_loop, r_loop
      REAL(rprec) :: xvec(3),bvec(3)
      
      INTEGER :: nsect, npts, j, k, nrad, nzed, chunk
      REAL(rprec) :: rs,rm, zs, axp, ayp, azp, r1, z1, dr, dz, r2, z2,&
                     my_rtol
      REAL(rprec), ALLOCATABLE :: bfield_data(:,:), bfield_data2(:,:)
      CHARACTER(256)  ::  id_string_temp
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Basic copy of MPI_COMM_DIANGO
      mylocalid = myworkid
      mylocalmaster = master
      MPI_COMM_LOCAL = MPI_COMM_DIAGNO
      numprocs_local = nprocs_diagno

      if(lverb) write(6,*)' -BENCHMARK CALCULATION'
      !CALL free_virtual_casing
      npts=10
      nsect=3
      SELECT CASE (TRIM(id_string))
         CASE ('bigtok','bigtok_phi1') ! BIGTOK
            r1 = 101.01; r2 = 102.01; z1 = 0.000; z2 = 1.000; nrad = 5 ; nzed = 5; c_loop = 100; r_loop = 1.5; my_rtol = 1.0E-3
         CASE ('a3tok') ! A3TOK
            r1 = 4.01;   r2 = 5.01;   z1 = 0.000; z2 = 1.000; nrad = 5 ; nzed = 5; c_loop = 3;   r_loop = 1.25; my_rtol = 1.0E-3
         CASE ('DIIID_m24n0s99_nfp1','DIIID_m24n0s99_nfp10','DIIID_m20n0s128_nfp1_lasym')
            r1 = 2.3;    r2 = 3.300;  z1 = 0.000; z2 = 1.000; nrad = 5 ; nzed = 5; c_loop = 1.6; r_loop = 1.25; my_rtol = 1.0E-3
         CASE ('ncsx') ! NCSX
            r1 = 1.79;   r2 = 2.79;   z1 = 0.000; z2 = 0.500; nrad = 5 ; nzed = 5; c_loop = 1.4; r_loop = 0.80; my_rtol = 1.0E-3
         CASE ('w7x') ! W7X
            r1 = 6.50;   r2 = 7.50;   z1 = 0.000; z2 = 0.500; nrad = 5 ; nzed = 5; c_loop = 5.5; r_loop = 1.0; my_rtol = 1.0E-2
         CASE('iter') ! ITER
            r1 = 8.20;   r2 = 9.200;  z1 = 0.000; z2 = 4.000; nrad = 5 ; nzed = 5; c_loop = 6.0; r_loop = 4.5; my_rtol = 1.0E-3
         CASE('lhd')  ! LHD
            r1 = 4.75;   r2 = 5.50;   z1 = 0.000; z2 = 1.000; nrad = 5 ; nzed = 5; c_loop = 3.6; r_loop = 1.0; my_rtol = 1.0E-3
         CASE('hsx') !HSX
            r1 = 1.55;   r2 = 2.0;   z1 = 0.000; z2 = 0.300; nrad = 5 ; nzed = 5; c_loop = 1.35; r_loop = 0.3; my_rtol = 1.0E-2
      END SELECT
      ! Setup Grid
      ig = 1
      dr = r2-r1
      dz = z2-z1
      ALLOCATE(bfield_data(nsect*nrad*nzed,15),bfield_data2(nsect*nrad*nzed,3))
      bfield_data(:,:) = 0.0
      DO i = 1, nsect
         DO j = 1, nrad
            DO k = 1, nzed
               rp = r1+dr*REAL(j-1)/REAL(nrad)
               xp = rp*COS(pi2*REAL(i-1)/REAL(nsect)/10)
               yp = rp*SIN(pi2*REAL(i-1)/REAL(nsect)/10)
               zp = z1+dz*REAL(k-1)/REAL(nzed)
               bfield_data(ig,1)=xp
               bfield_data(ig,2)=yp
               bfield_data(ig,3)=zp
               ig = ig + 1
            END DO
         END DO
      END DO
      bfield_data2(:,1:3)=bfield_data(:,1:3)
      ig = ig - 1
      ! First do volint
      lvc_field = .false.
      vc_adapt_tol = 0.0E-00
      vc_adapt_rel = my_rtol
      CALL diagno_init_vmec
      MIN_CLS = 0
      ! Create a fluxloop files
      ! Note that for DIA_LOOP we need to subtract off PHIEDGE if using VC but not J
      IF (myworkid == master) THEN
         CALL safe_open(iunit_out,ier,'test_loops_j.'//TRIM(id_string),'replace','formatted')
         WRITE(iunit_out,'(I6)') 3
         WRITE(iunit_out,'(3I6,A48)') 36,0,0,'DIA_LOOP'
         DO i = 1, 36
            xp = c_loop+r_loop*DCOS(pi2*(i-1)/36)
            yp = 0.0
            zp = r_loop*DSIN(pi2*(i-1)/36)
            WRITE(iunit_out,'(3ES22.12E3)') xp,yp,zp
         END DO
         WRITE(iunit_out,'(3I6,A48)') 36,0,0,'TOR_LOOP'
         DO i = 1, 36
            xp = (c_loop+r_loop)*DCOS(pi2*(i-1)/36)
            yp = (c_loop+r_loop)*DSIN(pi2*(i-1)/36)
            zp = 0.0
            WRITE(iunit_out,'(3ES22.12E3)') xp,yp,zp
         END DO
         WRITE(iunit_out,'(3I6,A48)') 36,0,0,'RAD_LOOP'
         DO i = 1, 36
            xp = (c_loop+r_loop)
            yp = 0.5*zmax_surf*DCOS(pi2*(i-1)/36)
            zp = 0.5*zmax_surf*DSIN(pi2*(i-1)/36)
            WRITE(iunit_out,'(3ES22.12E3)') xp,yp,zp
         END DO
         CLOSE(iunit_out)
         CALL safe_open(iunit_out,ier,'test_loops_b.'//TRIM(id_string),'replace','formatted')
         WRITE(iunit_out,'(I6)') 3
         WRITE(iunit_out,'(3I6,A48)') 36,0,1,'DIA_LOOP'
         DO i = 1, 36
            xp = c_loop+r_loop*DCOS(pi2*(i-1)/36)
            yp = 0.0
            zp = r_loop*DSIN(pi2*(i-1)/36)
            WRITE(iunit_out,'(3ES22.12E3)') xp,yp,zp
         END DO
         WRITE(iunit_out,'(3I6,A48)') 36,0,0,'TOR_LOOP'
         DO i = 1, 36
            xp = (c_loop+r_loop)*DCOS(pi2*(i-1)/36)
            yp = (c_loop+r_loop)*DSIN(pi2*(i-1)/36)
            zp = 0.0
            WRITE(iunit_out,'(3ES22.12E3)') xp,yp,zp
         END DO
         WRITE(iunit_out,'(3I6,A48)') 36,0,0,'RAD_LOOP'
         DO i = 1, 36
            xp = (c_loop+r_loop)
            yp = 0.5*zmax_surf*DCOS(pi2*(i-1)/36)
            zp = 0.5*zmax_surf*DSIN(pi2*(i-1)/36)
            WRITE(iunit_out,'(3ES22.12E3)') xp,yp,zp
         END DO
         CLOSE(iunit_out)
      END IF

#if defined(MPI_OPT)
      ! Broadcast Array sizes
      CALL MPI_BARRIER(MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'diagno_flux1',ierr_mpi)
#endif

      ! Divide up the work
      chunk = FLOOR(REAL(ig) / REAL(numprocs_local))
      mystart = myworkid*chunk + 1
      myend = mystart + chunk - 1
#if defined(MPI_OPT)
      IF (ALLOCATED(mnum)) DEALLOCATE(mnum)
      IF (ALLOCATED(moffsets)) DEALLOCATE(moffsets)
      ALLOCATE(mnum(numprocs_local), moffsets(numprocs_local))
      CALL MPI_ALLGATHER(chunk,1,MPI_INTEGER,mnum,1,MPI_INTEGER,MPI_COMM_LOCAL,ierr_mpi)
      CALL MPI_ALLGATHER(mystart,1,MPI_INTEGER,moffsets,1,MPI_INTEGER,MPI_COMM_LOCAL,ierr_mpi)
      i = 1
      DO
         IF ((moffsets(numprocs_local)+mnum(numprocs_local)-1) == ig) EXIT
         IF (i == numprocs_local) i = 1
         mnum(i) = mnum(i) + 1
         moffsets(i+1:numprocs_local) = moffsets(i+1:numprocs_local) + 1
         i=i+1
      END DO
      mystart = moffsets(mylocalid+1)
      chunk  = mnum(mylocalid+1)
      myend   = mystart + chunk - 1
#endif
      PRINT *,myworkid, mystart,myend,ig

      if(lverb) write(6,*)' ---VOLINT'
      bfield_data(:,4:15) = 0
!      IF (myworkid == master) THEN
         DO i = mystart, myend
         !DO i = 1, ig
            bxp = 0.0; byp = 0.0; bzp = 0.0;
            axp = 0.0; ayp = 0.0; azp = 0.0;
            xp  = bfield_data(i,1)
            yp  = bfield_data(i,2)
            zp  = bfield_data(i,3)
            ier = 1
            !CALL bfield_vc(xp,yp,zp,bxp,byp,bzp,ier)
            ier1 = 1
            !CALL vecpot_vc(xp,yp,zp,axp,ayp,azp,ier1)
            !PRINT *,myworkid,xp,yp,zp,ier,ier1,nlastcall
            bfield_data(i,4)=axp
            bfield_data(i,5)=ayp
            bfield_data(i,6)=azp
            bfield_data(i,7)=bxp
            bfield_data(i,8)=byp
            bfield_data(i,9)=bzp
         END DO
!      END IF
#if defined(MPI_OPT)
      ! Broadcast Array sizes
      CALL MPI_BARRIER(MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'diagno_flux1',ierr_mpi)
      DEALLOCATE(mnum, moffsets)
#endif
      flux_diag_file = 'test_loops_j.'//TRIM(id_string)
      id_string_temp = id_string
      id_string = TRIM(id_string) // '_j'
!      CALL diagno_flux
      id_string = id_string_temp
      CALL free_virtual_casing
      ! Now do VC
      if(lverb) write(6,*)' ---VIRTUAL CASING'
      vc_adapt_tol = 0.0
      vc_adapt_rel = 1.0E-03
      lvc_field = .true.
      CALL diagno_init_vmec
      MIN_CLS = 0

      ! Divide up the work
      chunk = FLOOR(REAL(ig) / REAL(numprocs_local))
      mystart = myworkid*chunk + 1
      myend = mystart + chunk - 1
#if defined(MPI_OPT)
      IF (ALLOCATED(mnum)) DEALLOCATE(mnum)
      IF (ALLOCATED(moffsets)) DEALLOCATE(moffsets)
      ALLOCATE(mnum(numprocs_local), moffsets(numprocs_local))
      CALL MPI_ALLGATHER(chunk,1,MPI_INTEGER,mnum,1,MPI_INTEGER,MPI_COMM_LOCAL,ierr_mpi)
      CALL MPI_ALLGATHER(mystart,1,MPI_INTEGER,moffsets,1,MPI_INTEGER,MPI_COMM_LOCAL,ierr_mpi)
      i = 1
      DO
         IF ((moffsets(numprocs_local)+mnum(numprocs_local)-1) == ig) EXIT
         IF (i == numprocs_local) i = 1
         mnum(i) = mnum(i) + 1
         moffsets(i+1:numprocs_local) = moffsets(i+1:numprocs_local) + 1
         i=i+1
      END DO
      mystart = moffsets(mylocalid+1)
      chunk  = mnum(mylocalid+1)
      myend   = mystart + chunk - 1
#endif
      PRINT *,myworkid, mystart,myend,ig
!      IF (myworkid == master) THEN
         DO i = mystart, myend
         !DO i = 1, ig
            bxp = 0.0; byp = 0.0; bzp = 0.0;
            axp = 0.0; ayp = 0.0; azp = 0.0;
            xp  = bfield_data(i,1)
            yp  = bfield_data(i,2)
            zp  = bfield_data(i,3)
            ier = 1
            CALL bfield_vc(xp,yp,zp,bxp,byp,bzp,ier)
            ier1 = 1
            CALL vecpot_vc(xp,yp,zp,axp,ayp,azp,ier1)
            !PRINT *,'*',xp,yp,zp,ier,ier1,nlastcall
            bfield_data(i,10)=axp
            bfield_data(i,11)=ayp
            bfield_data(i,12)=azp
            bfield_data(i,13)=bxp
            bfield_data(i,14)=byp
            bfield_data(i,15)=bzp
         END DO
!      END IF

#if defined(MPI_OPT)
      ! Broadcast Array sizes
      CALL MPI_BARRIER(MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'diagno_flux1',ierr_mpi)
      IF (myworkid == master) THEN
         CALL MPI_REDUCE(MPI_IN_PLACE,bfield_data,15*ig,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
      ELSE
         CALL MPI_REDUCE(bfield_data,bfield_data,15*ig,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
      END IF
      bfield_data(:,1:3)=bfield_data2(:,1:3)
      DEALLOCATE(mnum, moffsets)
#endif

      flux_diag_file = 'test_loops_b.'//TRIM(id_string)
      id_string_temp = id_string
      id_string = TRIM(id_string) // '_b'
      CALL diagno_flux
      id_string = id_string_temp
      IF (myworkid == master) THEN
         CALL safe_open(iunit_out,ier,'diagno_bench.'//TRIM(id_string),'replace','formatted')
         DO i = 1, ig
               WRITE(iunit_out,'(15E20.10)') bfield_data(i,1:15)
         END DO
         CLOSE(iunit_out)
         CALL safe_open(iunit_out,ier,'diagno_surf.'//TRIM(id_string),'replace','formatted')
         IF (.true.) CALL virtual_casing_surf_dump(iunit_out)
         CLOSE(iunit_out)
      END IF

      DEALLOCATE(bfield_data,bfield_data2)

#if defined(MPI_OPT)
      ! Broadcast Array sizes
      CALL MPI_BARRIER(MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'diagno_flux1',ierr_mpi)
#endif
      
!      STOP
 
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE diagno_bench
