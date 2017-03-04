!-----------------------------------------------------------------------
!     Subroutine:    torlines_emc3
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2015
!     Description:   This subroutine generates a grid file for EMC3.
!-----------------------------------------------------------------------
      SUBROUTINE torlines_emc3
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE safe_open_mod, ONLY: safe_open
      USE torlines_fieldlines
      USE torlines_background
      USE torlines_realspace
      USE torlines_runtime
      USE ez_hdf5
      USE mpi_params
!-----------------------------------------------------------------------
!     Local Variables
!          ierr        Error flag
!          iz          Index over zones
!          NZONET      Number of toroidal zones
!          
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, iunit_geo, iunit_grid, iunit_field
      INTEGER :: iz, ik, iu, ir, i1 ,i2, SRF_TOTAL
      REAL(rprec) :: delta_s, phi_offset
      REAL(rprec) :: rho_local, theta_local, phi_local
      REAL(rprec) :: r_inner, z_inner, b_inner
      REAL(rprec) :: r_outer, z_outer, b_outer
      REAL(rprec), ALLOCATABLE :: s_arr(:), phi_arr(:,:)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Handle Inner Grid point
      i1=COUNT(rho <= s_inner)
      IF (i1 > 12) THEN
         s_inbg  = rho(i1-6)
      ELSE
         i1 = 12
         s_inbg = rho(6)
      END IF
      ! Outter grid point
      i2=COUNT(rho <= s_outer)
      SRF_RADI = i2-i1+1  !Actual number of field lines followed
      SRF_TOTAL = SRF_RADI + 2
      ! Helper Arrays
      ALLOCATE(s_arr(SRF_RADI),phi_arr(2,NZONET))
      FORALL(iz = 1:SRF_RADI) s_arr(iz) = rho(i1+iz-1) ! S values to follow
      ! Screen output
      IF (lverb) THEN
         WRITE(6,'(a)')'----     WRITING: input.geo     ----'
         WRITE(6,'(A,I5)')'   TOROIDAL ZONES: ',NZONET
         WRITE(6,'(A,F6.3,A,F6.3,A,I4)')'     RHO=[',s_inner,',',s_outer,']; NR: ',SRF_RADI
         WRITE(6,'(A,F6.3,A,F6.3,A,I4)')'   THETA=[',0.0,',',pi2,']; NTHETA: ',SRF_POLO
         WRITE(6,'(A,F6.3,A,F6.3,A,I4)')'    ZETA=[',0.0,',',pi2/nfp/NZONET,']; NZETA: ',SRF_TORO
         CALL FLUSH(6)
      END IF
      ! Write GEO File
      IF (myid == master) THEN
         CALL safe_open(iunit_geo,ier,'input.geo','replace','formatted')
         WRITE(iunit_geo,'(A)') '* geometry information for EMC3'
         WRITE(iunit_geo,'(A)') '*--------------------------------'
         WRITE(iunit_geo,'(A)') '*** 1. grid resolution'
         WRITE(iunit_geo,'(A)') '*--------------------------------'
         WRITE(iunit_geo,'(A)') '* number of zones/blocks'
         WRITE(iunit_geo,'(I5)') NZONET
         WRITE(iunit_geo,'(A)') '* number of radial, poloidal and toroidal grid points'
         DO iz = 1, NZONET
            WRITE(iunit_geo,'(3I12)') SRF_TOTAL,SRF_POLO,SRF_TORO
         END DO
         WRITE(iunit_geo,'(A)') '*--------------------------------'
         WRITE(iunit_geo,'(A)') '*** 2. surface deffinition'
         WRITE(iunit_geo,'(A)') '*--------------------------------'
         WRITE(iunit_geo,'(A)') '*** 2.1 non default surface'
         WRITE(iunit_geo,'(A)') '* radial'
         WRITE(iunit_geo,'(I5)') 0
         WRITE(iunit_geo,'(A)') '* poloidal'
         WRITE(iunit_geo,'(I5)') 2*NZONET
         DO iz = 1, NZONET
            WRITE(iunit_geo,'(3I12)') 0,iz-1,1
            WRITE(iunit_geo,'(4I12)') 0,SRF_TOTAL-2,0,SRF_TORO-2
            WRITE(iunit_geo,'(3I12)') SRF_POLO-1,iz-1,1
            WRITE(iunit_geo,'(4I12)') 0,SRF_TOTAL-2,0,SRF_TORO-2
         END DO
         WRITE(iunit_geo,'(A)') '* toroidal'
         WRITE(iunit_geo,'(I5)') 2*NZONET
         DO iz = 1, NZONET
            WRITE(iunit_geo,'(3I12)') 0,iz-1,3
            WRITE(iunit_geo,'(4I12)') 0,SRF_TOTAL-2,0,SRF_POLO-2
            WRITE(iunit_geo,'(3I12)') SRF_TORO-1,iz-1,3
            WRITE(iunit_geo,'(4I12)') 0,SRF_TOTAL-2,0,SRF_POLO-2
         END DO
         WRITE(iunit_geo,'(A)') '***2.2 non transparent surface (Boundary condition must be defined)'
         WRITE(iunit_geo,'(A)') '* radial'
         WRITE(iunit_geo,'(I5)') 2*NZONET
         DO iz = 1, NZONET
            WRITE(iunit_geo,'(3I12)') 1,iz-1,1
            WRITE(iunit_geo,'(4I12)') 0,SRF_POLO-2,0,SRF_TORO-2
            WRITE(iunit_geo,'(3I12)') SRF_RADI-2,iz-1,-1
            WRITE(iunit_geo,'(4I12)') 0,SRF_POLO-2,0,SRF_TORO-2
         END DO
         WRITE(iunit_geo,'(A)') '* poloidal'
         WRITE(iunit_geo,'(I5)') 0
         WRITE(iunit_geo,'(A)') '* toroidal'
         WRITE(iunit_geo,'(I5)') 0
         WRITE(iunit_geo,'(A)') '*** 2.3 plate surface (Bohm Boundary condition)'
         WRITE(iunit_geo,'(A)') '* radial'
         WRITE(iunit_geo,'(I5)') -3
         WRITE(iunit_geo,'(A)') '* poloidal'
         WRITE(iunit_geo,'(I5)') -3
         WRITE(iunit_geo,'(A)') '* toroidal'
         WRITE(iunit_geo,'(I5)') -3
         WRITE(iunit_geo,'(A)') '*--------------------------------'
         WRITE(iunit_geo,'(A)') '*** 3. physical cell definition'
         WRITE(iunit_geo,'(A)') '*--------------------------------'
         WRITE(iunit_geo,'(I5)') -1
         WRITE(iunit_geo,'(A)') '* run cell check?'
         WRITE(iunit_geo,'(A)') 'F'
         CLOSE(iunit_geo)
      END IF

      ! WRITE GRID File
      IF (myid == master) CALL safe_open(iunit_grid,ier,'grid3d.dat','replace','formatted')
      IF (myid == master) CALL safe_open(iunit_field,ier,'field.dat','replace','formatted')
      ! Loop over zones
      phi_offset = phmx/(NZONET*2)
      NPOINC = (SRF_TORO-1)/2 !Index from 0 to NPOINC
      phi_arr(1,1) = 0.0
      phi_arr(2,1) = 2*phi_offset
      DO iz = 2, NZONET
         phi_arr(1,iz) = phi_arr(2,iz-1)
         phi_arr(2,iz) = phi_arr(1,iz) + 2*phi_offset
      END DO
      DO iz = 1, NZONET
         IF (myid == master) WRITE(iunit_grid,'(3I12)') SRF_TOTAL,SRF_POLO,SRF_TORO
         ! Loop over Radial grid
         R_START = -1; Z_START = -1; PHI_START = -1; PHI_END = -1; i1 = 0; i2 = 0
         ! Backwards
         DO iu = 1, SRF_POLO
            i1 = (iu-1)*SRF_RADI+1
            i2 = i1 + SRF_RADI - 1
            R_START(i1:i2) = s_arr(1:SRF_RADI)
            Z_START(i1:i2) = (iu-1)*pi2/(SRF_POLO-1)
            PHI_START(i1:i2) = phi_arr(1,iz)+phi_offset
            PHI_END(i1:i2)   = phi_arr(1,iz)
         END DO
         NLINES = i2
         IF (myid == master .AND. iz == 2) WRITE(6,*) '' 
         IF (myid == master) WRITE(6,'(A,I5,A,I5)')'---- Working on SECTOR:',iz,' of ',NZONET
         CALL FLUSH(6)
         CALL torlines_follow
         IF (myid == master) THEN
            !PHI_lines(1,NPOINC) = phi_arr(1,iz)
            DO ik = NPOINC, 0, -1
               phi_local   = PHI_lines(1,ik)
               WRITE(iunit_grid,'(F21.16)') phi_local*360/pi2
               DO iu = 1, SRF_POLO
                  i1 = (iu-1)*SRF_RADI+1
                  i2 = i1 + SRF_RADI - 1
                  ! Add inner and outer gridpoint
                  rho_local   = s_inbg
                  theta_local = U_lines(i1,ik)
                  CALL EZspline_interp(R_spl,rho_local,theta_local,phi_local,r_inner,ier)
                  CALL EZspline_interp(B_spl,rho_local,theta_local,phi_local,b_inner,ier)
                  rho_local = 1
                  theta_local = U_lines(i2,ik)
                  CALL EZspline_interp(R_spl,rho_local,theta_local,phi_local,r_outer,ier)
                  CALL EZspline_interp(B_spl,rho_local,theta_local,phi_local,b_outer,ier)
                  WRITE(iunit_grid,'(6F12.6)') r_inner*100,R_lines(i1:i2,ik)*100,r_outer*100 ! EMC3 in cm
                  WRITE(iunit_field,'(6F12.6)') b_inner,B_lines(i1:i2,ik),b_outer ! EMC3 in T
               END DO
               DO iu = 1, SRF_POLO
                  i1 = (iu-1)*SRF_RADI+1
                  i2 = i1 + SRF_RADI - 1
                  ! Add inner and outer gridpoint
                  rho_local   = s_inbg
                  theta_local = U_lines(i1,ik)
                  phi_local   = PHI_lines(1,ik)
                  CALL EZspline_interp(Z_spl,rho_local,theta_local,phi_local,z_inner,ier)
                  rho_local = 1
                  theta_local = U_lines(i2,ik)
                  CALL EZspline_interp(Z_spl,rho_local,theta_local,phi_local,z_outer,ier)
                  WRITE(iunit_grid,'(6F12.6)') z_inner*100,Z_lines(i1:i2,ik)*100,z_outer*100 ! EMC3 in cm
               END DO
            END DO
            CALL FLUSH(iunit_grid)
         END IF
         IF (iz >= 1) lverb = .false.
         ! Forwards
         R_START = -1; Z_START = -1; PHI_START = -1; PHI_END = -1; i1 = 0; i2 = 0
         DO iu = 1, SRF_POLO
            i1 = (iu-1)*SRF_RADI+1
            i2 = i1 + SRF_RADI - 1
            R_START(i1:i2) = s_arr(1:SRF_RADI)
            Z_START(i1:i2) = (iu-1)*pi2/(SRF_POLO-1)
            PHI_START(i1:i2) = phi_arr(2,iz)-phi_offset
            PHI_END(i1:i2)   = phi_arr(2,iz)
         END DO
         NLINES = i2
         CALL torlines_follow
         IF (myid == master) THEN
            DO ik = 1, NPOINC
               phi_local   = PHI_lines(1,ik)
               WRITE(iunit_grid,'(F21.16)') phi_local*360/pi2
               DO iu = 1, SRF_POLO
                  i1 = (iu-1)*SRF_RADI+1
                  i2 = i1 + SRF_RADI - 1
                  ! Add inner and outer gridpoint
                  rho_local   = s_inbg
                  theta_local = U_lines(i1,ik)
                  CALL EZspline_interp(R_spl,rho_local,theta_local,phi_local,r_inner,ier)
                  CALL EZspline_interp(B_spl,rho_local,theta_local,phi_local,b_inner,ier)
                  rho_local = 1
                  theta_local = U_lines(i2,ik)
                  CALL EZspline_interp(R_spl,rho_local,theta_local,phi_local,r_outer,ier)
                  CALL EZspline_interp(B_spl,rho_local,theta_local,phi_local,b_outer,ier)
                  WRITE(iunit_grid,'(6F12.6)') r_inner*100,R_lines(i1:i2,ik)*100,r_outer*100 ! EMC3 in cm
                  WRITE(iunit_field,'(6F12.6)') b_inner,B_lines(i1:i2,ik),b_outer ! EMC3 in T
               END DO
               DO iu = 1, SRF_POLO
                  i1 = (iu-1)*SRF_RADI+1
                  i2 = i1 + SRF_RADI - 1
                  ! Add inner and outer gridpoint
                  rho_local   = s_inbg
                  theta_local = U_lines(i1,ik)
                  phi_local   = PHI_lines(1,ik)
                  CALL EZspline_interp(Z_spl,rho_local,theta_local,phi_local,z_inner,ier)
                  rho_local = 1
                  theta_local = U_lines(i2,ik)
                  CALL EZspline_interp(Z_spl,rho_local,theta_local,phi_local,z_outer,ier)
                  WRITE(iunit_grid,'(6F12.6)') z_inner*100,Z_lines(i1:i2,ik)*100,z_outer*100 ! EMC3 in cm
               END DO
            END DO
            CALL FLUSH(iunit_grid)
         END IF
      END DO
      IF (myid == master) lverb = .true.
      ! Close files
      IF (myid == master) CLOSE(iunit_grid)
      IF (myid == master) CLOSE(iunit_field)
      ! Deallocate Arrays
      DEALLOCATE(s_arr,phi_arr)
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE torlines_emc3
