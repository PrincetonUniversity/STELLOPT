!-----------------------------------------------------------------------
!     Module:        beams3d_write_ascoth4
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          08/21/2019
!     Description:   This subroutine outputs simulation data for the
!                    ASCOT4 code. (FOR particleGenerator ONLY)
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_write_ascoth4(write_type)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
#ifdef LHDF5
      USE hdf5
      USE ez_hdf5
#endif 
      USE beams3d_lines
      USE beams3d_grid, ONLY: nr, nphi, nz, B_R, B_PHI, B_Z, raxis, &
                                 zaxis, phiaxis, S_ARR, U_ARR, POT_ARR, &
                                 ZEFF_ARR, TE, TI, NE, req_axis, zeq_axis, npot, &
                                 POT_SPL_S, ezspline_interp, phiedge_eq, TE_spl_s, &
                                 NE_spl_s, TI_spl_s, ZEFF_spl_s, POT_spl_s, vp_spl_s, &
                                 nne, nte, nti, nzeff, npot
      USE beams3d_runtime, ONLY: id_string, npoinc, nbeams, beam, t_end, lverb, &
                                    lvmec, lpies, lspec, lcoil, lmgrid, lbeam, &
                                    lvessel, lvac, lbeam_simple, handle_err, nparticles_start, &
                                    HDF5_OPEN_ERR,HDF5_WRITE_ERR,&
                                    HDF5_CLOSE_ERR, BEAMS3D_VERSION, weight, e_beams, p_beams,&
                                    charge, Zatom, mass, ldepo, v_neut,lcollision, pi, pi2, &
                                    t_end_in, nprocs_beams
      USE safe_open_mod, ONLY: safe_open
      USE wall_mod, ONLY: nface,nvertex,face,vertex,ihit_array, wall_free
      USE beams3d_write_par
      USE mpi_params
      USE mpi
!-----------------------------------------------------------------------
!     Input Variables
!          write_type  Type of write to preform
!-----------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN):: write_type
!-----------------------------------------------------------------------
!     Local Variables
!          ier          Error Flag
!          iunit        File ID
!-----------------------------------------------------------------------
      INTEGER :: ier, iunit, i, d1, d2, d3, k, k1, k2, kmax
      INTEGER(HID_T) :: options_gid, bfield_gid, efield_gid, plasma_gid, &
                        neutral_gid, wall_gid, marker_gid, qid_gid
      INTEGER, ALLOCATABLE :: i1dtemp(:)
      DOUBLE PRECISION :: dbl_temp
      DOUBLE PRECISION, ALLOCATABLE :: rtemp(:,:,:), r1dtemp(:)
      CHARACTER(LEN=10) ::  qid_str
      CHARACTER(LEN=8) :: temp_str8


      DOUBLE PRECISION, PARAMETER :: e_charge      = 1.60217662E-19 !e_c
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      SELECT CASE (TRIM(write_type))
         CASE('INIT')
            IF (myworkid == master) THEN
               CALL DATE_AND_TIME(DATE=temp_str8)
               CALL open_hdf5('ascot4_'//TRIM(id_string)//'.h5',fid,ier,LCREATE=.true.)
               IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,'beams3d_ascot5_'//TRIM(id_string)//'.h5',ier)
               ! Write the ID
               !CALL write_var_hdf5(fid,'/options/',nz,ier,DBLVAR=zaxis,ATT='Vertical Axis [rad]',ATT_NAME='description')

               !--------------------------------------------------------------
               !           B-FIELD
               !--------------------------------------------------------------
               CALL h5gcreate_f(fid,'bfield', bfield_gid, ier)
               CALL h5gcreate_f(bfield_gid,'stellarator', qid_gid, ier)
               CALL h5gcreate_f(qid_gid,'profiles', plasma_gid, ier)
               CALL write_att_hdf5(qid_gid,'date',temp_str8,ier)
               CALL write_att_int_hdf5(qid_gid,'version',2,ier) ! Otherwise stellarator symmetry is assumed
               CALL write_att_hdf5(qid_gid,'description','Data initialized from BEAMS3D',ier)
               ALLOCATE(rtemp(nr,nphi,nz))
               rtemp = RESHAPE(B_R,(/nr,nphi,nz/),ORDER=(/1,2,3/))
               CALL write_var_hdf5(qid_gid,'br',nr,nphi,nz,ier,DBLVAR=rtemp)
               rtemp = RESHAPE(B_PHI,(/nr,nphi,nz/),ORDER=(/1,2,3/))
               CALL write_var_hdf5(qid_gid,'bphi',nr,nphi,nz,ier,DBLVAR=rtemp)
               rtemp = RESHAPE(B_Z,(/nr,nphi,nz/),ORDER=(/1,2,3/))
               CALL write_var_hdf5(qid_gid,'bz',nr,nphi,nz,ier,DBLVAR=rtemp)
               rtemp = RESHAPE(S_ARR,(/nr,nphi,nz/),ORDER=(/1,2,3/))
               WHERE(rtemp < 1.0E-3) rtemp = 1.0E-3
               CALL write_var_hdf5(qid_gid,'s',nr,nphi,nz,ier,DBLVAR=rtemp)
               DEALLOCATE(rtemp)
               CALL write_var_hdf5(qid_gid,'phi',nphi,ier,DBLVAR=360*phiaxis/pi2)
               CALL write_var_hdf5(qid_gid,'r',nr,ier,DBLVAR=raxis)
               CALL write_var_hdf5(qid_gid,'z',nz,ier,DBLVAR=zaxis)
               CALL write_var_hdf5(qid_gid,'toroidalPeriods',ier,INTVAR=FLOOR(pi2/phiaxis(nphi)))
               CALL write_var_hdf5(qid_gid,'symmetrymode',ier,INTVAR=1)
               d3 = FLOOR(pi2/phiaxis(nphi))
               ALLOCATE(rtemp((nphi-1)*d3+1,3,1))
               DO i = 1, d3
                  d1 = (i-1)*(nphi-1)+1
                  d2 = (i)*(nphi-1)
                  rtemp(d1:d2,1,1)=req_axis(1:(nphi-1))
                  rtemp(d1:d2,2,1)=zeq_axis(1:(nphi-1))
                  rtemp(d1:d2,3,1)=360*(phiaxis(1:(nphi-1))+(i-1)*pi2/d3)/pi2
               END DO
               rtemp((nphi-1)*d3+1,1:2,1) = rtemp(1,1:2,1)
               rtemp((nphi-1)*d3+1,3,1) = 360
               CALL write_var_hdf5(qid_gid,'axis_R',(nphi-1)*d3+1,ier,DBLVAR=rtemp(:,1,1))
               CALL write_var_hdf5(qid_gid,'axis_z',(nphi-1)*d3+1,ier,DBLVAR=rtemp(:,2,1))
               CALL write_var_hdf5(qid_gid,'axis_phi',(nphi-1)*d3+1,ier,DBLVAR=rtemp(:,3,1))
               DEALLOCATE(rtemp)

               ALLOCATE(rtemp(nr,2,1))
               rtemp = 0
               rtemp(:,2,1) = 0
               DO i = 1, nr
                  rtemp(i,1,1)=DBLE(i-1)/DBLE(nr-1)
               END DO
               CALL EZspline_interp( vp_spl_s, nr, rtemp(:,1,1), rtemp(:,2,1), ier)
               DO i = 2, nr
                  rtemp(i,2,1)=rtemp(i,2,1)+(rtemp(i-1,2,1)) ! Integral So volume
               END DO
               CALL write_var_hdf5(plasma_gid,'s',nr,ier,DBLVAR=sqrt(rtemp(:,1,1)))
               CALL write_var_hdf5(plasma_gid,'volume',nr,ier,DBLVAR=rtemp(:,2,1))
               DEALLOCATE(rtemp)
               
               CALL h5gclose_f(plasma_gid, ier)
               CALL h5gclose_f(qid_gid, ier)
               CALL h5gclose_f(bfield_gid, ier)

               !--------------------------------------------------------------
               !           WALL
               !--------------------------------------------------------------
               CALL h5gcreate_f(fid,'wall', wall_gid, ier)
               CALL h5gcreate_f(wall_gid,'3d', qid_gid, ier)
               CALL write_att_hdf5(qid_gid,'date',temp_str8,ier)
               CALL write_att_hdf5(qid_gid,'description','Data initialized from BEAMS3D',ier)
               CALL write_var_hdf5(qid_gid,'nelements',ier,INTVAR=nface)
               ALLOCATE(rtemp(3,nface,3))
               ALLOCATE(i1dtemp(nface))
               i1dtemp = 0
               DO i = 1, nface
                  d1 = face(i,1)
                  d2 = face(i,2)
                  d3 = face(i,3)
                  ! X
                  rtemp(1,i,1) = vertex(d1,1)
                  rtemp(2,i,1) = vertex(d2,1)
                  rtemp(3,i,1) = vertex(d3,1)
                  ! Y
                  rtemp(1,i,2) = vertex(d1,2)
                  rtemp(2,i,2) = vertex(d2,2)
                  rtemp(3,i,2) = vertex(d3,2)
                  ! Z
                  rtemp(1,i,3) = vertex(d1,3)
                  rtemp(2,i,3) = vertex(d2,3)
                  rtemp(3,i,3) = vertex(d3,3)
               END DO
               CALL write_var_hdf5( qid_gid, 'triangles_flag', nface, ier, INTVAR=i1dtemp)
               CALL write_var_hdf5( qid_gid, 'triangles_x1x2x3', 3, nface, ier, DBLVAR=rtemp(:,:,1))
               CALL write_var_hdf5( qid_gid, 'triangles_y1y2y3', 3, nface, ier, DBLVAR=rtemp(:,:,2))
               CALL write_var_hdf5( qid_gid, 'triangles_z1z2z3', 3, nface, ier, DBLVAR=rtemp(:,:,3))
               DEALLOCATE(rtemp,i1dtemp)
               CALL h5gclose_f(qid_gid, ier)
               CALL h5gclose_f(wall_gid, ier)

               ! Close file
               CALL close_hdf5(fid,ier)
               IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,'beams3d_ascot4_'//TRIM(id_string)//'.h5',ier)

               ! Write the input text file
               ier = 0; iunit = 411
               CALL safe_open(iunit,ier,'input.plasma_1d','replace','formatted')
               WRITE(iunit,'(A)') 'Created by BEAMS3D from run: '//TRIM(id_string)
               WRITE(iunit,'(A)') 'Profiles from input file.  Assuming Carbon Impurity.'
               WRITE(iunit,'(A)') temp_str8
               ALLOCATE(rtemp(nr,5,1))
               rtemp = 0
               rtemp(:,5,1) = 1
               DO i = 1, nr
                  rtemp(i,1,1)=DBLE(i-1)/DBLE(nr-1)
               END DO ! Treat rtemp(:,1,1) as rho not s.
               IF (nte > 0)   CALL EZspline_interp( TE_spl_s,   nr, rtemp(:,1,1)**2, rtemp(:,2,1), ier)
               IF (nne > 0)   CALL EZspline_interp( NE_spl_s,   nr, rtemp(:,1,1)**2, rtemp(:,3,1), ier)
               IF (nti > 0)   CALL EZspline_interp( TI_spl_s,   nr, rtemp(:,1,1)**2, rtemp(:,4,1), ier)
               IF (nzeff > 0) CALL EZspline_interp( ZEFF_spl_s, nr, rtemp(:,1,1)**2, rtemp(:,5,1), ier)
               rtemp(nr,1,1) = 1.0; rtemp(nr,2,1) = 1.0; rtemp(nr,4,1) = 1.0 ! Default Te and Ti to almost zero
               WRITE(iunit,'(2X,I4,2X,I4,2X,A)') nr,2,'# Nrad,Nion'
               WRITE(iunit,'(2X,I4,2X,I4,2X,A)') 1,6,'# ion Znum'
               WRITE(iunit,'(2X,I4,2X,I4,2X,A)') 1,12,'# ion Anum'
               WRITE(iunit,'(2X,I4,2X,I4,2X,I4,2X,A)') 1,1,1,'# OBSOLETE VALUES. PUT 1'
               WRITE(iunit,'(A)') 'RHO (pol)     Te (eV)         Ne (1/m3)       Vtor_I (rad/s)  Ti1 (eV)        Ni1 (1/m3)      Ni2 (1/m3) ...'
               DO i = 1, nr
                  dbl_temp = (rtemp(i,5,1)-1)/6.0 ! Frac nH=frac*nC;   Zeff*ni=nH+6*nC=(1+6*frac)*nH = ne; nH = ne/(1*6*frac)
                  WRITE(iunit,'(7(2X,ES18.10))') rtemp(i,1,1),rtemp(i,2,1),rtemp(i,3,1),-999.0,rtemp(i,4,1),&
                                              rtemp(i,3,1)/(1+6*dbl_temp),MAX(dbl_temp*rtemp(i,3,1)/(1+6*dbl_temp),1.0E16)
               END DO
               DEALLOCATE(rtemp)
               CLOSE(iunit)

               ! Write the input erad file
               ier = 0; iunit = 411
               CALL safe_open(iunit,ier,'input.erad','replace','formatted')
               WRITE(iunit,'(A)') 'Created by BEAMS3D from run: '//TRIM(id_string)
               ALLOCATE(rtemp(nr,2,1))
               rtemp = 0
               rtemp(:,2,1) = 0
               DO i = 1, nr
                  rtemp(i,1,1)=DBLE(i-1)/DBLE(nr-1)
               END DO
               IF (npot > 0)   CALL EZspline_interp( POT_spl_s,   nr, rtemp(:,1,1), rtemp(:,2,1), ier)
               WRITE(iunit,'(2X,I4,2X,A)') nr,'# Nradial points'
               DO i = 1, nr
                  WRITE(iunit,'(2(2X,ES18.10))') rtemp(i,1,1),rtemp(i,2,1)
               END DO
               DEALLOCATE(rtemp)

               CLOSE(iunit)
            END IF
         CASE('MARKER')
            ! Downselect the markers
            d1 = LBOUND(R_lines,DIM=2)
            d2 = UBOUND(R_lines,DIM=2)
            d3 = 0
            IF (lbeam) d3 = 1 ! ASCOT4 wants the neutral information
            ALLOCATE(i1dtemp(nprocs_beams))
            i1dtemp = 0
            i1dtemp(myworkid+1) = COUNT(end_state(d1:d2)==0)
            CALL MPI_ALLREDUCE(MPI_IN_PLACE, i1dtemp, nprocs_beams, MPI_INTEGER, MPI_SUM, MPI_COMM_BEAMS, ierr_mpi)
            k2 = SUM(i1dtemp(1:myworkid+1))
            k1 = k2 - i1dtemp(myworkid+1) + 1
            kmax = SUM(i1dtemp)
            DEALLOCATE(i1dtemp)
            ALLOCATE(rtemp(k1:k2,16,1))
            k = k1
            DO i = d1, d2
               IF (end_state(i)>0) CYCLE
               rtemp(k,1,1) = NINT(mass(i)*5.97863320194E26) ! Anum
               rtemp(k,2,1) = mass(i)*5.97863320194E26 ! mass
               rtemp(k,3,1) = Zatom(i)
               rtemp(k,4,1) = Zatom(i)
               rtemp(k,5,1) = 0.5*mass(i)*SUM(v_neut(:,i)*v_neut(:,i),DIM=1)*e_charge
               rtemp(k,6,1) = SQRT(S_lines(d3,i))
               dbl_temp = PHI_lines(d3,i)
               rtemp(k,7,1) = dbl_temp*180/pi
               rtemp(k,8,1) = R_lines(d3,i)
               rtemp(k,9,1) = Z_lines(d3,i)
               ! Now we get a little out of order since we need the components of the velocity
               rtemp(k,10,1) = v_neut(2,i)*cos(dbl_temp)-v_neut(1,i)*sin(dbl_temp) ! Vphi
               rtemp(k,11,1) = v_neut(1,i)*cos(dbl_temp)+v_neut(2,i)*sin(dbl_temp) ! Vr
               rtemp(k,12,1) = v_neut(3,i) !Vz
               rtemp(k,13,1) = Beam(i)
               rtemp(k,14,1) = weight(i) ! weight
               rtemp(k,15,1) = k
               rtemp(k,16,1) = t_end(k)
               ! Now 17-19 are Bphi, Br, Bz
               ! Set to zero for now.
               k=k+1
            END DO

            ! First write the header
            IF (myworkid==master) THEN
               ! Write the input text file
               ier = 0; iunit = 411
               CALL safe_open(iunit,ier,'input.particles','replace','formatted')
               WRITE(iunit,'(A)') ' PARTICLE INITIAL DATA FOR ASCOT'
               WRITE(iunit,'(A)') ' 4 VERSION ====================='
               WRITE(iunit,'(A)') ' '
               WRITE(iunit,'(A)') ' 2  # Number of comment lines, max length 256 char, (defined in prtRead_lineLength)'
               WRITE(iunit,'(A,F5.2)') 'Particles Generated by BEAMS3D v.',BEAMS3D_VERSION
               WRITE(iunit,'(A)') 'Created from run: '//TRIM(id_string)
               WRITE(iunit,'(A)') ' '
               WRITE(iunit,'(I8,A)') kmax,' # Number of particles (-1 means unknown number)'
               WRITE(iunit,'(A)') ' '
               WRITE(iunit,'(A)') '16 # Number of different fields for each particle [10 first letters are significant]'
               WRITE(iunit,'(A)') 'Anum      - mass number of particle        (integer)'
               WRITE(iunit,'(A)') 'mass      - mass of the particle           (amu)'
               WRITE(iunit,'(A)') 'Znum      - charge number of particle      (integer)'
               WRITE(iunit,'(A)') 'charge    - charge of the particle         (elemental charge)'
               WRITE(iunit,'(A)') 'energy    - kinetic energy of particle     (eV)'
               WRITE(iunit,'(A)') 'rho       - normalized poloidal flux       (real)'
               WRITE(iunit,'(A)') 'phiprt    - toroidal angle of particle     (deg)'
               WRITE(iunit,'(A)') 'Rprt      - R of particle                  (m)'
               WRITE(iunit,'(A)') 'zprt      - z of particle                  (m)'
               WRITE(iunit,'(A)') 'vphi      - toroidal velocity of particle  (m/s)'
               WRITE(iunit,'(A)') 'vR        - radial velocity of particle    (m/s)'
               WRITE(iunit,'(A)') 'vz        - vertical velocity of particle  (m/s)'
               WRITE(iunit,'(A)') 'origin    - origin of the particle         ()'
               WRITE(iunit,'(A)') 'weight    - weight factor of particle      (particle/second)'
               WRITE(iunit,'(A)') 'id        - unique identifier of particle  (integer)'
               WRITE(iunit,'(A)') 'Tmax      - maximum time to follow the prt (s)'
               WRITE(iunit,'(A)') ' '
               CALL FLUSH(iunit)
               CLOSE(iunit)
            END IF
            ! Now write the file
            DO i = 0, nprocs_beams-1
               IF (myworkid == i) THEN
                  CALL safe_open(iunit,ier,'input.particles','old','formatted',ACCESS_IN='APPEND')
                  DO k = k1,k2
                     WRITE(iunit,'(2(I8,E14.5E3),8(E18.9E3),2(I8,E18.9E3))') &
                        NINT(rtemp(k,1,1)),rtemp(k,2,1),NINT(rtemp(k,3,1)),rtemp(k,4:12,1),&
                        NINT(rtemp(k,13,1)),rtemp(k,14,1),NINT(rtemp(k,15,1)),rtemp(k,16,1)
                  END DO
                  CALL FLUSH(iunit)
                  CLOSE(iunit)
                  DEALLOCATE(rtemp)
               END IF
               CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
            END DO
            IF (myworkid == master) THEN
               CALL safe_open(iunit,ier,'input.particles','old','formatted',ACCESS_IN='APPEND')
               WRITE(iunit,'(A)') '#EOF'
               CALL FLUSH(iunit)
               CLOSE(iunit)
            END IF

      END SELECT

      RETURN

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE beams3d_write_ascoth4
