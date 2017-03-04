!-----------------------------------------------------------------------
!     Module:        fieldlines_init
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/21/2012
!     Description:   This subroutine initializes the fields on the
!                    R, phi, Z grid.  
!-----------------------------------------------------------------------
      SUBROUTINE fieldlines_init
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE vmec_input,  ONLY: extcur_in => extcur, read_indata_namelist,&
                             nv_in => nzeta, nfp_in => nfp, nigroup
      USE read_wout_mod, ONLY:  read_wout_file, rmin_surf, rmax_surf, zmax_surf
      USE fieldlines_runtime
      USE fieldlines_grid, ONLY: raxis,phiaxis,zaxis, nr, nphi, nz, &
                                 rmin, rmax, zmin, zmax, phimin, &
                                 phimax, B_R, B_Z, B_PHI, MU3D, &
                                 BR_spl, BZ_spl, MU_spl, delta_phi, MODB_spl
      USE fieldlines_input_mod, ONLY: read_fieldlines_input
      USE fieldlines_lines, ONLY: nlines
      USE wall_mod
      USE random, ONLY: random_normal
      USE mpi_params                                                    ! MPI
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
!DEC$ IF DEFINED (NAG)
!      LOGICAL        :: licval
!      INTEGER        :: mkmaj, mkmin
!      CHARACTER(128) :: impl,prec,pcode,hdware,opsys,fcomp,vend
!DEC$ ENDIF
!DEC$ IF DEFINED (MPI_OPT)
      INCLUDE 'mpif.h'                                                          ! MPI
!DEC$ ENDIF  
      INTEGER :: i,j,k,ier, iunit
      INTEGER :: bcs1(2), bcs2(2), bcs3(2)
      REAL(rprec) :: rmin_temp, rmax_temp, zmin_temp, zmax_temp,&
                     phimax_temp, bx_err, by_err, bz_err, br_err, bphi_err, b_err,&
                     ang_err, phimin_temp
      REAL(rprec) :: xw,yw,zw
      LOGICAL     :: lhit
      DOUBLE PRECISION :: phi_q, q(2), qdot(2)
!-----------------------------------------------------------------------
!     External Functions
!          A00ADF               NAG Detection
!-----------------------------------------------------------------------
!      EXTERNAL A00ADF
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! TESTING
      IF (lvessel .and. lverb .and. .false.) THEN
         CALL wall_load_txt(TRIM(vessel_string),ier)
         IF (lverb) CALL wall_info(6)
         !CALL collide(5.0_rprec,5.0_rprec,0.0_rprec,5.0_rprec,5.0_rprec,1.5_rprec,xw,yw,zw,lhit)
         CALL collide(6*1/100._rprec,6*1/100._rprec,0.0_rprec,6*1/100._rprec,6*1/100._rprec,1.5_rprec,xw,yw,zw,lhit)
         WRITE(*,*) xw,yw,zw,lhit
         DO i = 1, 100
            DO j = 1, 100
               CALL collide(6*i/100._rprec,6*j/100._rprec,0.0_rprec,6*i/100._rprec,6*j/100._rprec,1.5_rprec,xw,yw,zw,lhit)
               IF (lhit) WRITE(328,*) xw,yw,zw
               CALL FLUSH(328)
            END DO
         END DO
         CALL wall_dump('test',ier)
      END IF


      ! Limit what the user can do
      IF (laxis_i .and. .not.lvac) lvac = .true.
      
      
      ! First Read The Input Namelist
      iunit = 11
      IF (lverb) WRITE(6,'(A)') '----- Input Parameters -----'
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init',ierr_mpi)
!DEC$ ENDIF
      IF (lvmec) THEN
         CALL read_fieldlines_input('input.' // TRIM(id_string),ier,myid)
         IF (lverb) WRITE(6,'(A)') '   FILE: input.' // TRIM(id_string)
         IF (.not. lvac) CALL read_wout_file(TRIM(id_string),ier)
      ELSE IF (lpies) THEN
         CALL read_fieldlines_input(TRIM(id_string) // '.in',ier,myid)
         IF (lverb) WRITE(6,'(A)') '   FILE: ' // TRIM(id_string) // '.in'
      ELSE IF (lspec) THEN
         CALL read_fieldlines_input('input.' // TRIM(id_string),ier,myid)
         IF (lverb) WRITE(6,'(A)') '   FILE: input.' // TRIM(id_string)
      END IF

      ! TESTING LMU
      IF (lmu .and. .false.) THEN
         q = 0
         DO i = 1,100000
            phi_q = 0
            WRITE(327,*) q(1),q(2)
            q(1) = q(1) + mu * random_normal()
            q(2) = q(2) + mu * random_normal()
         END DO
         STOP
      END IF
      
      IF (lplasma_only) THEN
         IF (lvmec) THEN
            rmin = 1.2*rmin_surf
            rmax = 1.2*rmax_surf
            zmin = -1.2*zmax_surf
            zmax = 1.2*zmax_surf
         ELSE IF (lpies) THEN
         ELSE IF (lspec) THEN
         END IF
      END IF
      
      IF (lauto) THEN
         nlines = nr
         IF (nruntype == runtype_full) THEN
            rmin_temp = r_start(1)
            zmin_temp = z_start(1)
            phimin_temp = phi_start(1)
         ELSE
            !rmin_temp = MINVAL(r_start,MASK = r_start > 0)
            !zmin_temp = MINVAL(z_start,MASK = r_start > 0)
            
            rmin_temp = r_start(1)
            zmin_temp = z_start(1)
            phimin_temp = phi_start(1)
         END IF
         rmax_temp = MAXVAL(r_start,MASK = r_start > 0)
         zmax_temp = MAXVAL(z_start,MASK = r_start > 0)
         phimax_temp = MAXVAL(phi_end,MASK = r_start > 0)
         r_start = -1
         z_start = -1
         phi_end = 0.0
         IF (nlines > MAXLINES) nlines = MAXLINES-2
         DO i = 1, nlines
            r_start(i)   = (rmax_temp-rmin_temp)*REAL(i-1)/REAL(nlines-1) + rmin_temp
            z_start(i)   = (zmax_temp-zmin_temp)*REAL(i-1)/REAL(nlines-1) + zmin_temp
            phi_start(i) = phimin_temp
            phi_end(i)   = phimax_temp
         END DO
      END IF
      
      ! Output some information
      IF (lverb .and. .not.lrestart) THEN
!DEC$ IF DEFINED (NAG)
!         This doesn't work because A00ADF not supported in all NAG Marks
!         CALL A00ADF(impl,prec,pcode,mkmaj,mkmin,hdware,opsys,fcomp,vend,licval)
!         IF (.not.licval) THEN
!            WRITE(6,'(A)') '!!!!!!  ERROR: Invalid NAG License  !!!!!!'
!            ier = -1
!            CALL handle_err(NAG_ERR,'fieldlines_init',ier)
!         END IF
!         WRITE(6,'(A,I3,A1,I2,A)') '   NAG Mark ',mkmaj,mkmin,' Detected!'
!DEC$ ELSE
!         CALL handle_err(NAG_ERR,'fieldlines_init',ier)
!DEC$ ENDIF  
         WRITE(6,'(A,F8.5,A,F8.5,A,I4)') '   R   = [',rmin,',',rmax,'];  NR:   ',nr
         WRITE(6,'(A,F8.5,A,F8.5,A,I4)') '   PHI = [',phimin,',',phimax,'];  NPHI: ',nphi
         WRITE(6,'(A,F8.5,A,F8.5,A,I4)') '   Z   = [',zmin,',',zmax,'];  NZ:   ',nz
         IF (lauto) WRITE(6,'(A)') '   AUTO CALCULATED STARTING POINTS!'
         WRITE(6,'(A,I6)')               '   # of Fieldlines: ',nlines
         IF (lvac) WRITE(6,'(A)') '   VACUUM FIELDS ONLY!'
         IF (lplasma_only) WRITE(6,'(A)') '   PLASMA FIELDS ONLY!'
         IF (lbfield_only) WRITE(6,'(A)') '   ONLY MAGNETIC FIELDS WILL BE OUTPUT! (NO FIELDLINES)'
         IF (lafield_only) WRITE(6,'(A)') '   ONLY VECTOR POTENTIAL WILL BE OUTPUT! (NO FIELDLINES)'
         IF (lemc3) WRITE(6,'(A)') '   EMC3-EIRENE OUTPUT WILL BE GENERATED! (NO FIELDLINES)'
         IF (lmodb) WRITE(6,'(A)') '   SAVING |B| ALONG FIELDLINE!'
         IF (lmu) WRITE(6,'(A)')   '   DIFFUSION OPERATOR TURNED ON!'
         CALL FLUSH(6)
      END IF
       
      ! Create the background grid and initialize if necessary
      IF (lrestart) THEN
         CALL fieldlines_init_restart
      ELSE
         ALLOCATE(raxis(nr),zaxis(nz),phiaxis(nphi),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'RAXIS ZAXIS PHIAXIS',ier)
         ALLOCATE(B_R(nr,nphi,nz),B_PHI(nr,nphi,nz),B_Z(nr,nphi,nz),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'BR BZ',ier)
         FORALL(i = 1:nr) raxis(i) = (i-1)*(rmax-rmin)/(nr-1) + rmin
         FORALL(i = 1:nz) zaxis(i) = (i-1)*(zmax-zmin)/(nz-1) + zmin
         FORALL(i = 1:nphi) phiaxis(i) = (i-1)*(phimax-phimin)/(nphi-1) + phimin
         B_R = 0; B_PHI = 0; B_Z = 0
         ! Put the vacuum field on the background grid
         IF (lmgrid) THEN
            CALL fieldlines_init_mgrid
         ELSE IF (lcoil) THEN
            CALL fieldlines_init_coil
         ELSE IF (lnescoil) THEN
            CALL fieldlines_init_nescoil
         END IF
      END IF

      ! Put the plasma field on the background grid
      IF (lrestart) THEN
         ! Do NOTHING
      ELSE IF (lvmec .and. .not.lvac) THEN
         CALL fieldlines_init_vmec
         IF (ledge_start) CALL fieldlines_init_vmec_edgestart
      ELSE IF (lpies .and. .not.lvac) THEN
         !CALL fieldlines_init_pies
      ELSE IF (lspec .and. .not.lvac) THEN
         !CALL fieldlines_init_spec
      END IF

      !ERROR FIELD code section
      ! Note that for this to make sense the code must be run with NFP=1 and PHIMAX=2*pi
      IF (lerror_field) THEN
         IF (lverb) WRITE(6,'(A)') '!!!!!ADDING STATIC ERROR FIELD!!!!!'
         b_err   = 2.5*2.0E-5
         ang_err = pi2*0./360.
         b_err  = SQRT((b_err*b_err)/2)
         bx_err = b_err * COS(ang_err)
         by_err = b_err * SIN(ang_err)
         bz_err = 0
         DO k = 1, nz
            DO j = 1, nphi
               DO i = 1, nr
                  br_err = bx_err*cos(phiaxis(j)) + by_err*sin(phiaxis(j))
                  bphi_err = by_err*cos(phiaxis(j)) - bx_err*sin(phiaxis(j))
                  B_R(i,j,k) = B_R(i,j,k) + br_err
                  B_PHI(i,j,k) = B_PHI(i,j,k) + bphi_err
                  B_Z(i,j,k) = B_Z(i,j,k) + bz_err
               END DO
            END DO
         END DO
      END IF

      ! Put curtor on axis and calculate the field
      IF (laxis_i)  CALL fieldlines_init_I
      
      IF (ANY(B_PHI .eq. 0)) THEN
!DEC$ IF DEFINED (MPI_OPT)
         CALL MPI_FINALIZE(ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_FINE_ERR,'fieldlines_init',ierr_mpi)
!DEC$ ENDIF
         stop 'ERROR: B_PHI = 0 Found'
      END IF
      
      ! Handle outputting the B-FIELD
      IF (lemc3 .or. lbfield_only .or. lafield_only) THEN
         IF (lemc3 .and. myid==master) CALL fieldlines_write_emc3
!DEC$ IF DEFINED (MPI_OPT)
         CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init',ierr_mpi)
!DEC$ ENDIF
         RETURN
      END IF

      ! Handle reversing the field
      IF (lreverse) phi_end = -phi_end

      ! Handle MOD_B
      IF (lmodb) THEN
         ALLOCATE(MU3D(nr,nphi,nz),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'MU3D1',ier)
         MU3D = sqrt(B_R**2+B_PHI**2+B_Z**2)
         bcs1=(/ 0, 0/)
         bcs2=(/-1,-1/)
         bcs3=(/ 0, 0/)
         CALL EZspline_init(MODB_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'fieldlines_init:MODB_spl',ier)
         MODB_spl%isHermite = 1
         MODB_spl%x1 = raxis
         MODB_spl%x2 = phiaxis
         MODB_spl%x3 = zaxis
         CALL EZspline_setup(MODB_spl,MU3D,ier,EXACT_DIM=.true.)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'fieldlines_init:MODB_spl',ier)
         DEALLOCATE(MU3D)
      END IF

      ! Handle mu
      IF (lmu) THEN
         CALL init_random_seed()
         ALLOCATE(MU3D(nr,nphi,nz),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'MU3D2',ier)
         MU3D = ABS(2*mu*sqrt(B_R**2+B_PHI**2+B_Z**2)/B_PHI)
         DO i = 1, nr
            MU3D(i,:,:) = raxis(i) * MU3D(i,:,:)
         END DO
         bcs1=(/ 0, 0/)
         bcs2=(/-1,-1/)
         bcs3=(/ 0, 0/)
         CALL EZspline_init(MU_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'fieldlines_init:MU_SPL',ier)
         MU_spl%isHermite = 1
         MU_spl%x1 = raxis
         MU_spl%x2 = phiaxis
         MU_spl%x3 = zaxis
         CALL EZspline_setup(MU_spl,MU3D,ier,EXACT_DIM=.true.)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'fieldlines_init:MU_SPL',ier)
         DEALLOCATE(MU3D)
         CALL RANDOM_SEED()
      END IF

      ! Get setup vessel
      IF (lvessel) THEN
         CALL wall_load_txt(TRIM(vessel_string),ier)
         IF (ier /= 0) CALL handle_err(WALL_ERR,'fieldlines_init',ier)
         IF (myid /= master) DEALLOCATE(vertex,face) ! Do this to save memory
         IF (lverb) THEN
            CALL wall_info(6)
            IF (lwall_trans) WRITE(6,'(A)')'   !!!!! Poincare Screen  !!!!!'
         END IF
      END IF

      ! Now we need to reformulate B_R and B_Z as functions of phi
      DO k = 1, nz
         DO j = 1, nphi
            DO i = 1, nr
               B_R(i,j,k) = raxis(i) * B_R(i,j,k) / B_PHI(i,j,k)
               B_Z(i,j,k) = raxis(i) * B_Z(i,j,k) / B_PHI(i,j,k)
            END DO
         END DO
      END DO

      ! Construct Splines
      bcs1=(/ 0, 0/)
      bcs2=(/-1,-1/)
      bcs3=(/ 0, 0/)
      CALL EZspline_init(BR_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
      IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'fieldlines_init',ier)
      CALL EZspline_init(BZ_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
      IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'fieldlines_init',ier)
      BR_spl%isHermite = 0
      BZ_spl%isHermite = 0
      BR_spl%x1 = raxis
      BZ_spl%x1 = raxis
      BR_spl%x2 = phiaxis
      BZ_spl%x2 = phiaxis
      BR_spl%x3 = zaxis
      BZ_spl%x3 = zaxis
      IF (lverb) THEN
         WRITE(6,'(A)')'----- Constructing Splines -----'
         WRITE(6,'(A,F8.5,A,F8.5,A,I4)') '   R   = [',MINVAL(raxis),',',MAXVAL(raxis),'];  NR:   ',BR_spl%n1
         WRITE(6,'(A,F8.5,A,F8.5,A,I4)') '   PHI = [',MINVAL(phiaxis),',',MAXVAL(phiaxis),'];  NPHI: ',BR_spl%n2
         WRITE(6,'(A,F8.5,A,F8.5,A,I4)') '   Z   = [',MINVAL(zaxis),',',MAXVAL(zaxis),'];  NZ:   ',BR_spl%n3
         IF (lmu) WRITE(6,'(A,E15.5)')   '   MU  = ',mu
         WRITE(6,'(A,I1)')               '   HERMITE FORM: ',BR_spl%isHermite
      END IF
      CALL EZspline_setup(BR_spl,B_R,ier,EXACT_DIM=.true.)
      IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'fieldlines_init',ier)
      CALL EZspline_setup(BZ_spl,B_Z,ier,EXACT_DIM=.true.)
      IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'fieldlines_init',ier)

      IF (.FALSE.) THEN
         qdot = 0
         q(1) = 0.7
         q(2) = 0.0
         phi_q = 0.0
         delta_phi = 0.01
         CALL fblin_nag(phi_q,q,qdot)
         PRINT *,phi_q,q,qdot
      END IF
      
      ! DEALLOCATE Variables
      ! Only master retains a copy for output
      IF (myid /= master) THEN
         DEALLOCATE(raxis,zaxis,phiaxis)
         DEALLOCATE(B_R,B_Z,B_PHI)
      END IF

!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init',ierr_mpi)
!DEC$ ENDIF
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE fieldlines_init
