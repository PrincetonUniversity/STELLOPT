!-----------------------------------------------------------------------
!     Module:        beams3d_init
!     Authors:       S. Lazerson (lazerson@pppl.gov), M. McMillan (matthew.mcmillan@my.wheaton.edu)
!     Date:          06/20/2012
!     Description:   This subroutine initializes the fields on the
!                    R, phi, Z grid.  
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_init
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE vmec_input,  ONLY: extcur_in => extcur, read_indata_namelist,&
                             nv_in => nzeta, nfp_in => nfp, nigroup
      USE beams3d_runtime
      USE beams3d_grid, ONLY: raxis,phiaxis,zaxis, nr, nphi, nz, nte, &
                                 nne, nti, npot, &
                                 rmin, rmax, zmin, zmax, phimin, &
                                 phimax, B_R, B_Z, B_PHI, MODB, &
                                 BR_spl, BZ_spl, BPHI_spl, MODB_spl, &
                                 TE_spl, NE_spl, TI_spl, TE, NE, TI, &
                                 TE_spl_s, NE_spl_s, TI_spl_s, S_ARR,&
                                 S_spl, U_ARR, U_spl, &
                                 POT_ARR, POT_spl_s, POT_spl
      USE beams3d_input_mod, ONLY: read_beams3d_input
      USE beams3d_lines, ONLY: nparticles
      USE wall_mod
      USE mpi_params                                                    ! MPI
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
!DEC$ IF DEFINED (NAG)
      LOGICAL        :: licval
      INTEGER        :: mkmaj, mkmin
      CHARACTER(128) :: impl,prec,pcode,hdware,opsys,fcomp,vend
!DEC$ ENDIF
!DEC$ IF DEFINED (MPI_OPT)
      INCLUDE 'mpif.h'                                                          ! MPI
!DEC$ ENDIF  
      INTEGER :: i,j,k,ier, iunit, nextcur_in
      INTEGER :: bcs1(2), bcs2(2), bcs3(2), bcs1_s(2)
      REAL(rprec) :: br, bphi, bz, ti_temp
!-----------------------------------------------------------------------
!     External Functions
!          A00ADF               NAG Detection
!-----------------------------------------------------------------------
!      EXTERNAL A00ADF
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      bcs1=(/ 0, 0/)
      bcs2=(/-1,-1/)
      bcs3=(/ 0, 0/)
      
      ! Handle Consistency checks
      !IF (ldepo .and. lvessel) lplasma_only = .FALSE.
      
      ! First Read The Input Namelist
      iunit = 11
      IF (lverb) WRITE(6,'(A)') '----- Input Parameters -----'
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'beams3d_init',ierr_mpi)
!DEC$ ENDIF

      IF (lvmec .and. lread_input) THEN
         CALL read_beams3d_input('input.' // TRIM(id_string),ier)
         IF (lverb) WRITE(6,'(A)') '   FILE: input.' // TRIM(id_string)
      ELSE IF (lpies .and. lread_input) THEN
         CALL read_beams3d_input(TRIM(id_string) // '.in',ier)
         IF (lverb) WRITE(6,'(A)') '   FILE: ' // TRIM(id_string) // '.in'
      ELSE IF (lspec .and. lread_input) THEN
         CALL read_beams3d_input('input.' // TRIM(id_string),ier)
         IF (lverb) WRITE(6,'(A)') '   FILE: input.' // TRIM(id_string)
      END IF

      ! Output some information
      IF (lverb .and. .not.lrestart) THEN
         WRITE(6,'(A,F9.5,A,F9.5,A,I4)') '   R   = [',rmin,',',rmax,'];  NR:   ',nr
         WRITE(6,'(A,F8.5,A,F8.5,A,I4)') '   PHI = [',phimin,',',phimax,'];  NPHI: ',nphi
         WRITE(6,'(A,F8.5,A,F8.5,A,I4)') '   Z   = [',zmin,',',zmax,'];  NZ:   ',nz
         IF (lbeam) THEN
            WRITE(6,'(A,I6)')               '   # of Particles to Start: ', nparticles_start
            WRITE(6,'(A,I6)')                          '   # of Beams: ', nbeams
         ELSE
            WRITE(6,'(A,I6)')               '   # of Particles to Start: ', nparticles
         END IF
         IF (lcollision) WRITE(6,'(A)') '   COLLISION OPERATOR ON!'
         IF (lvac)  WRITE(6,'(A)') '   VACUUM FIELDS ONLY!'
         IF (ldepo) WRITE(6,'(A)') '   DEPOSITION ONLY!'
         IF (lw7x) WRITE(6,'(A)') '   W7-X BEAM Model!'
         CALL FLUSH(6)
      END IF

      ! Construct 1D splines
      bcs1_s=(/ 0, 0 /)
      IF (lvmec .and. .not.lvac) THEN
          ! TE
          IF (nte>0) THEN
             CALL EZspline_init(TE_spl_s,nte,bcs1_s,ier)
             IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init1',ier)
             TE_spl_s%isHermite   = 1
             TE_spl_s%x1          = TE_AUX_S(1:nte)
             CALL EZspline_setup(TE_spl_s,TE_AUX_F(1:nte),ier,EXACT_DIM=.true.)
             IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init2',ier)
          END IF
          ! TI
          IF (nti>0) THEN
             CALL EZspline_init(TI_spl_s,nti,bcs1_s,ier)
             IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init3',ier)
             TI_spl_s%isHermite   = 1
             TI_spl_s%x1          = TI_AUX_S(1:nti)
             CALL EZspline_setup(TI_spl_s,TI_AUX_F(1:nti),ier,EXACT_DIM=.true.)
             IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init4',ier)
          END IF
          ! NE
          IF (nne>0) THEN
             CALL EZspline_init(NE_spl_s,nne,bcs1_s,ier)
             IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init5',ier)
             NE_spl_s%x1          = NE_AUX_S(1:nne)
             NE_spl_s%isHermite   = 1
             CALL EZspline_setup(NE_spl_s,NE_AUX_F(1:nne),ier,EXACT_DIM=.true.)
             IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init6',ier)
          END IF
          ! POTENTIAL
          IF (npot>0) THEN
             CALL EZspline_init(POT_spl_s,npot,bcs1_s,ier)
             IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init7',ier)
             POT_spl_s%x1          = POT_AUX_S(1:npot)
             POT_spl_s%isHermite   = 1
             CALL EZspline_setup(POT_spl_s,POT_AUX_F(1:npot),ier,EXACT_DIM=.true.)
             IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init8',ier)
          END IF
         
      END IF


      IF (lrestart) THEN
         CALL beams3d_init_restart
      ELSE
         ! Create the background grid
         ALLOCATE(raxis(nr),zaxis(nz),phiaxis(nphi),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'RAXIS PHIAXIS ZAXIS',ier)
         ALLOCATE(B_R(nr,nphi,nz),B_PHI(nr,nphi,nz),B_Z(nr,nphi,nz),MODB(nr,nphi,nz),&
                  TE(nr,nphi,nz), NE(nr,nphi,nz), TI(nr,nphi,nz), POT_ARR(nr,nphi,nz),&
                  S_ARR(nr,nphi,nz), U_ARR(nr,nphi,nz), STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'BR BZ',ier)
         FORALL(i = 1:nr) raxis(i) = (i-1)*(rmax-rmin)/(nr-1) + rmin
         FORALL(i = 1:nz) zaxis(i) = (i-1)*(zmax-zmin)/(nz-1) + zmin
         FORALL(i = 1:nphi) phiaxis(i) = (i-1)*(phimax-phimin)/(nphi-1) + phimin
         S_ARR = 1.5
         POT_ARR = 0
         ! Put the vacuum field on the background grid
         IF (lmgrid) THEN
            CALL beams3d_init_mgrid
         ELSE IF (lcoil) THEN
            CALL beams3d_init_coil
         END IF
      END IF


  
   
      ! Put the plasma field on the background grid
      IF (lvmec .and. .not.lvac .and. nte > 0) THEN
         CALL beams3d_init_vmec
      ELSE IF (lpies .and. .not.lvac) THEN
         !CALL beams3d_init_pies
      ELSE IF (lspec .and. .not.lvac) THEN
         !CALL beams3d_init_spec
      END IF

      ! Construct 3D Profile Splines
      IF (.not. lvac) THEN
         CALL EZspline_init(TE_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init: TE',ier)
         CALL EZspline_init(NE_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init: NE',ier)
         CALL EZspline_init(TI_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init: TI',ier)
         CALL EZspline_init(POT_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init: POT',ier)
         TE_spl%isHermite   = 1
         NE_spl%isHermite   = 1
         TI_spl%isHermite   = 1
         POT_spl%isHermite   = 1
         TE_spl%x1   = raxis
         NE_spl%x1   = raxis
         TI_spl%x1   = raxis
         POT_spl%x1   = raxis
         TE_spl%x2   = phiaxis
         NE_spl%x2   = phiaxis
         TI_spl%x2   = phiaxis
         POT_spl%x2   = phiaxis
         TE_spl%x3   = zaxis
         NE_spl%x3   = zaxis
         TI_spl%x3   = zaxis
         POT_spl%x3   = zaxis
         CALL EZspline_setup(TE_spl,TE,ier,EXACT_DIM=.true.)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init: TE',ier)
         CALL EZspline_setup(NE_spl,NE,ier,EXACT_DIM=.true.)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init: NE',ier)
         CALL EZspline_setup(TI_spl,TI,ier,EXACT_DIM=.true.)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init: TI',ier)
         CALL EZspline_setup(POT_spl,POT_ARR,ier,EXACT_DIM=.true.)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init: POT',ier)
         IF (myworkid /= master) DEALLOCATE(TE, NE, TI, POT_ARR)
      END IF
         
      ! Construct MODB
      MODB = SQRT(B_R*B_R+B_PHI*B_PHI+B_Z*B_Z)
      
      ! Get setup vessel
      IF (lvessel .and. (.not. lplasma_only .or. ldepo)) THEN
         IF (myworkid == master) PRINT *,'Loading vessel:',TRIM(vessel_string)
         CALL wall_load_txt(TRIM(vessel_string),ier)
         IF (myworkid /= master) DEALLOCATE(vertex,face) ! Do this to save memory
         IF (lverb) CALL wall_info(6)
         CALL FLUSH(6)
      END IF

      ! Initialize Random Number generator
      CALL RANDOM_SEED
      
!      Initialize beams (define a distribution of directions and weights)
      IF (lbeam) THEN
          IF (lw7x) THEN
             CALL beams3d_init_beams_w7x
          ELSE
             CALL beams3d_init_beams
          END IF
      ELSE
          ALLOCATE(  R_start(nparticles), phi_start(nparticles), Z_start(nparticles), &
          & v_neut(3,nparticles), mass(nparticles), charge(nparticles), &
          & mu_start(nparticles), Zatom(nparticles), t_end(nparticles), vll_start(nparticles), &
          & beam(nparticles), weight(nparticles_start, nbeams)  )

          R_start = r_start_in(1:nparticles)
          phi_start = phi_start_in(1:nparticles)
          Z_start = z_start_in(1:nparticles)
          vll_start = vll_start_in(1:nparticles)
          v_neut = 0.0
          weight = 1.0
          Zatom = Zatom_in(1:nparticles)
          mass = mass_in(1:nparticles)
          charge = charge_in(1:nparticles)
          mu_start = mu_start_in(1:nparticles)
          t_end = t_end_in(1:nparticles)
          beam  = 0
      END IF; CALL FLUSH(6)

      ! Do a reality check
      IF (ANY(ABS(vll_start)>3E8) .and. lverb) THEN
            ! This is an error code check
            PRINT *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            PRINT *,'!!!!!  Super-luminal particle velocity detected  !!!!!'
            PRINT *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            STOP
      END IF


      ! Construct Splines
      bcs1=(/ 0, 0/)
      bcs2=(/-1,-1/)
      bcs3=(/ 0, 0/)
      CALL EZspline_init(BR_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
      IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init:BR_spl',ier)
      CALL EZspline_init(BZ_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
      IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init:BZ_spl',ier)
      CALL EZspline_init(BPHI_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
      IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init:BPHI_spl',ier)
      CALL EZspline_init(MODB_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
      IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init:MODB_spl',ier)
      CALL EZspline_init(S_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
      IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init:S_spl',ier)
      CALL EZspline_init(U_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
      IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init:S_spl',ier)
      BR_spl%isHermite   = 0
      BZ_spl%isHermite   = 0
      BPHI_spl%isHermite = 0
      MODB_spl%isHermite = 0
      S_spl%isHermite = 0
      U_spl%isHermite = 0
      BR_spl%x1   = raxis
      BZ_spl%x1   = raxis
      BPHI_spl%x1 = raxis
      MODB_spl%x1 = raxis
      S_spl%x1 = raxis
      U_spl%x1 = raxis
      BR_spl%x2   = phiaxis
      BZ_spl%x2   = phiaxis
      BPHI_spl%x2 = phiaxis
      MODB_spl%x2 = phiaxis
      S_spl%x2 = phiaxis
      U_spl%x2 = phiaxis
      BR_spl%x3   = zaxis
      BZ_spl%x3   = zaxis
      BPHI_spl%x3 = zaxis
      MODB_spl%x3 = zaxis
      S_spl%x3 = zaxis
      U_spl%x3 = zaxis
      IF (lverb) THEN
         WRITE(6,'(A)')'----- Constructing Splines -----'
         WRITE(6,'(A,F9.5,A,F9.5,A,I4)') '   R   = [',MINVAL(raxis),',',MAXVAL(raxis),'];  NR:   ',BR_spl%n1
         WRITE(6,'(A,F8.5,A,F8.5,A,I4)') '   PHI = [',MINVAL(phiaxis),',',MAXVAL(phiaxis),'];  NPHI: ',BR_spl%n2
         WRITE(6,'(A,F8.5,A,F8.5,A,I4)') '   Z   = [',MINVAL(zaxis),',',MAXVAL(zaxis),'];  NZ:   ',BR_spl%n3
         WRITE(6,'(A,I1)')               '   HERMITE FORM: ',BR_spl%isHermite
         CALL FLUSH(6)
      END IF


      CALL EZspline_setup(BR_spl,B_R,ier,EXACT_DIM=.true.)
      IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init:BR_spl',ier)
      CALL EZspline_setup(BZ_spl,B_Z,ier,EXACT_DIM=.true.)
      IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init:BZ_spl',ier)
      CALL EZspline_setup(BPHI_spl,B_PHI,ier,EXACT_DIM=.true.)
      IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init:BPHI_spl',ier)
      CALL EZspline_setup(MODB_spl,MODB,ier,EXACT_DIM=.true.)
      IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init:MODB_spl',ier)
      CALL EZspline_setup(S_spl,S_ARR,ier,EXACT_DIM=.true.)
      IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init:S_spl',ier)
      CALL EZspline_setup(U_spl,U_ARR,ier,EXACT_DIM=.true.)
      IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init:U_spl',ier)

      ! Output Grid
      CALL beams3d_write('GRID_INIT')

      ! DEALLOCATE Variables
      IF (lvmec .and. .not.lvac) THEN
         IF (nte > 0) CALL EZspline_free(TE_spl_s,ier)
         IF (nne > 0) CALL EZspline_free(NE_spl_s,ier)
         IF (nti > 0) CALL EZspline_free(TI_spl_s,ier)
         IF (npot > 0) CALL EZspline_free(POT_spl_s,ier)
      END IF

      IF (myworkid /= master) THEN
         DEALLOCATE(raxis,zaxis,phiaxis)
         DEALLOCATE(B_R,B_Z,B_PHI,MODB,S_ARR,U_ARR)
      END IF
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'beams3d_init',ierr_mpi)
!DEC$ ENDIF
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE beams3d_init
