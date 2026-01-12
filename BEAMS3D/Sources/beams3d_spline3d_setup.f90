!-----------------------------------------------------------------------
!     Module:        beams3d_spline3d_setup
!     Authors:       S. Lazerson
!     Date:          10/22/2025
!     Description:   Automate allocation and deallocation of 3D splines
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_spline3d_setup
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE beams3d_globals, ONLY: nr, nphi, nz
      USE beams3d_runtime, ONLY: handle_err, lvac, EZSPLINE_ERR, &
         MPI_BARRIER_ERR
      USE mpi_params, ONLY: master, myid_sharmem, MPI_COMM_SHARMEM
      USE mpi_sharmem, ONLY: mpialloc, mpidealloc
      USE mpi_inc
      USE EZspline_obj, ONLY: EZspline3_r8
      USE EZspline, ONLY: EZspline_init, EZspline_setup, EZspline_free
      USE stel_kinds, ONLY: rprec
      USE beams3d_grid

!-----------------------------------------------------------------------
!     Input parameters
!-----------------------------------------------------------------------
      IMPLICIT NONE


!-----------------------------------------------------------------------
!     Local variables
!-----------------------------------------------------------------------
      INTEGER :: i, ier
      INTEGER :: bcs1(2), bcs2(2), bcs3(2)
      TYPE(EZspline3_r8) :: F_spl

!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------  
      ier = 0 ; bcs1=(/ 0, 0/); bcs2=(/-1,-1/); bcs3=(/ 0, 0/)

      IF (.not. lvac) THEN
!-----------------------------------------------------------------------
!        TE
!----------------------------------------------------------------------- 
         CALL mpialloc(TE4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_TE4D)
         IF (myid_sharmem == master) THEN
            CALL EZspline_init(F_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
            IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZsplineinit',ier)        
            F_spl%isHermite   = 1
            F_spl%x1   = raxis
            F_spl%x2   = phiaxis
            F_spl%x3   = zaxis
            CALL EZspline_setup(F_spl,TE,ier,EXACT_DIM=.true.)
            IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZspline_setup',ier)
            TE4D = F_spl%fspl
            CALL EZspline_free(F_spl,ier)
         END IF
         CALL mpidealloc(TE,win_TE)
         IF (nte > 0) CALL EZspline_free(TE_spl_s,ier)

!-----------------------------------------------------------------------
!        NE
!----------------------------------------------------------------------- 
         CALL mpialloc(NE4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_NE4D)
         IF (myid_sharmem == master) THEN
            CALL EZspline_init(F_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
            IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZsplineinit',ier)        
            F_spl%isHermite   = 1
            F_spl%x1   = raxis
            F_spl%x2   = phiaxis
            F_spl%x3   = zaxis
            CALL EZspline_setup(F_spl,NE,ier,EXACT_DIM=.true.)
            IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZspline_setup',ier)
            NE4D = F_spl%fspl
            CALL EZspline_free(F_spl,ier)
         END IF
         CALL mpidealloc(NE,win_NE)
         IF (nne > 0) CALL EZspline_free(NE_spl_s,ier)

!-----------------------------------------------------------------------
!        TI
!----------------------------------------------------------------------- 
         CALL mpialloc(TI4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_TI4D)
         IF (myid_sharmem == master) THEN
            CALL EZspline_init(F_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
            IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZsplineinit',ier)        
            F_spl%isHermite   = 1
            F_spl%x1   = raxis
            F_spl%x2   = phiaxis
            F_spl%x3   = zaxis
            CALL EZspline_setup(F_spl,TI,ier,EXACT_DIM=.true.)
            IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZspline_setup',ier)
            TI4D = F_spl%fspl
            CALL EZspline_free(F_spl,ier)
         END IF
         CALL mpidealloc(TI,win_TI)
         IF (nti > 0) CALL EZspline_free(TI_spl_s,ier)

!-----------------------------------------------------------------------
!        ZEFF
!----------------------------------------------------------------------- 
         CALL mpialloc(ZEFF4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_ZEFF4D)
         IF (myid_sharmem == master) THEN
            CALL EZspline_init(F_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
            IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZsplineinit',ier)        
            F_spl%isHermite   = 1
            F_spl%x1   = raxis
            F_spl%x2   = phiaxis
            F_spl%x3   = zaxis
            CALL EZspline_setup(F_spl,ZEFF_ARR,ier,EXACT_DIM=.true.)
            IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZspline_setup',ier)
            ZEFF4D = F_spl%fspl
            CALL EZspline_free(F_spl,ier)
         END IF
         CALL mpidealloc(ZEFF_ARR,win_ZEFF_ARR)
         IF (nzeff > 0) CALL EZspline_free(ZEFF_spl_s,ier)

!-----------------------------------------------------------------------
!        OMEG
!----------------------------------------------------------------------- 
         CALL mpialloc(OMEG4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_OMEG4D)
         IF (myid_sharmem == master) THEN
            CALL EZspline_init(F_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
            IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZsplineinit',ier)        
            F_spl%isHermite   = 1
            F_spl%x1   = raxis
            F_spl%x2   = phiaxis
            F_spl%x3   = zaxis
            CALL EZspline_setup(F_spl,OMEG_ARR,ier,EXACT_DIM=.true.)
            IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZspline_setup',ier)
            OMEG4D = F_spl%fspl
            CALL EZspline_free(F_spl,ier)
         END IF
         CALL mpidealloc(OMEG_ARR,win_OMEG_ARR)
         IF (nomeg > 0) CALL EZspline_free(OMEG_spl_s,ier)

!-----------------------------------------------------------------------
!        NI 5D
!----------------------------------------------------------------------- 
         CALL mpialloc(NI5D, 8, nr, nphi, nz, NION, myid_sharmem, 0, MPI_COMM_SHARMEM, win_NI5D)
         IF (myid_sharmem == master) THEN
            DO i = 1, NION
               CALL EZspline_init(F_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
               IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZsplineinit',ier)        
               F_spl%isHermite   = 1
               F_spl%x1   = raxis
               F_spl%x2   = phiaxis
               F_spl%x3   = zaxis
               CALL EZspline_setup(F_spl,NI(i,:,:,:),ier,EXACT_DIM=.true.)
               IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZspline_setup',ier)
               NI5D(:,:,:,:,i) = F_spl%fspl
               CALL EZspline_free(F_spl,ier)
            END DO
         END IF
         CALL mpidealloc(NI,win_NI)
         IF (nzeff > 0) THEN
            DO i = 1, NION
               CALL EZspline_free(NI_spl_s(i),ier)
            END DO
         END IF

      END IF !.not.lvac

!-----------------------------------------------------------------------
!     B_R
!----------------------------------------------------------------------- 
      CALL mpialloc(BR4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_BR4D)
      IF (myid_sharmem == master) THEN
         CALL EZspline_init(F_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZsplineinit',ier)        
         F_spl%isHermite   = 1
         F_spl%x1   = raxis
         F_spl%x2   = phiaxis
         F_spl%x3   = zaxis
         CALL EZspline_setup(F_spl,B_R,ier,EXACT_DIM=.true.)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZspline_setup',ier)
         BR4D = F_spl%fspl
         CALL EZspline_free(F_spl,ier)
      END IF
      CALL mpidealloc(B_R,win_B_R)

!-----------------------------------------------------------------------
!     B_PHI
!----------------------------------------------------------------------- 
      CALL mpialloc(BPHI4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_BPHI4D)
      IF (myid_sharmem == master) THEN
         CALL EZspline_init(F_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZsplineinit',ier)        
         F_spl%isHermite   = 1
         F_spl%x1   = raxis
         F_spl%x2   = phiaxis
         F_spl%x3   = zaxis
         CALL EZspline_setup(F_spl,B_PHI,ier,EXACT_DIM=.true.)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZspline_setup',ier)
         BPHI4D = F_spl%fspl
         CALL EZspline_free(F_spl,ier)
      END IF
      CALL mpidealloc(B_PHI,win_B_PHI)

!-----------------------------------------------------------------------
!     B_Z
!----------------------------------------------------------------------- 
      CALL mpialloc(BZ4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_BZ4D)
      IF (myid_sharmem == master) THEN
         CALL EZspline_init(F_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZsplineinit',ier)        
         F_spl%isHermite   = 1
         F_spl%x1   = raxis
         F_spl%x2   = phiaxis
         F_spl%x3   = zaxis
         CALL EZspline_setup(F_spl,B_Z,ier,EXACT_DIM=.true.)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZspline_setup',ier)
         BZ4D = F_spl%fspl
         CALL EZspline_free(F_spl,ier)
      END IF
      CALL mpidealloc(B_Z,win_B_Z)

!-----------------------------------------------------------------------
!     MODB
!----------------------------------------------------------------------- 
      CALL mpialloc(MODB4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_MODB4D)
      IF (myid_sharmem == master) THEN
         CALL EZspline_init(F_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZsplineinit',ier)        
         F_spl%isHermite   = 1
         F_spl%x1   = raxis
         F_spl%x2   = phiaxis
         F_spl%x3   = zaxis
         CALL EZspline_setup(F_spl,MODB,ier,EXACT_DIM=.true.)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZspline_setup',ier)
         MODB4D = F_spl%fspl
         CALL EZspline_free(F_spl,ier)
      END IF
      CALL mpidealloc(MODB,win_MODB)

!-----------------------------------------------------------------------
!     U_ARR
!----------------------------------------------------------------------- 
      CALL mpialloc(U4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_U4D)
      IF (myid_sharmem == master) THEN
         CALL EZspline_init(F_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZsplineinit',ier)        
         F_spl%isHermite   = 1
         F_spl%x1   = raxis
         F_spl%x2   = phiaxis
         F_spl%x3   = zaxis
         CALL EZspline_setup(F_spl,U_ARR,ier,EXACT_DIM=.true.)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZspline_setup',ier)
         U4D = F_spl%fspl
         CALL EZspline_free(F_spl,ier)
      END IF
      CALL mpidealloc(U_ARR,win_U_ARR)

!-----------------------------------------------------------------------
!     POT_ARR
!----------------------------------------------------------------------- 
      CALL mpialloc(POT4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_POT4D)
      IF (myid_sharmem == master) THEN
         CALL EZspline_init(F_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZsplineinit',ier)        
         F_spl%isHermite   = 1
         F_spl%x1   = raxis
         F_spl%x2   = phiaxis
         F_spl%x3   = zaxis
         CALL EZspline_setup(F_spl,POT_ARR,ier,EXACT_DIM=.true.)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZspline_setup',ier)
         POT4D = F_spl%fspl
         CALL EZspline_free(F_spl,ier)
      END IF
      CALL mpidealloc(POT_ARR,win_POT_ARR)

!-----------------------------------------------------------------------
!     RHO_ARR
!----------------------------------------------------------------------- 
      CALL mpialloc(RHO4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_RHO4D)
      IF (myid_sharmem == master) THEN
         CALL EZspline_init(F_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZsplineinit',ier)        
         F_spl%isHermite   = 1
         F_spl%x1   = raxis
         F_spl%x2   = phiaxis
         F_spl%x3   = zaxis
         CALL EZspline_setup(F_spl,RHO_ARR,ier,EXACT_DIM=.true.)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZspline_setup',ier)
         RHO4D = F_spl%fspl
         CALL EZspline_free(F_spl,ier)
      END IF
      CALL mpidealloc(RHO_ARR,win_RHO_ARR)

!-----------------------------------------------------------------------
!     XRHO_ARR
!----------------------------------------------------------------------- 
      CALL mpialloc(XRHO4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_XRHO4D)
      IF (myid_sharmem == master) THEN
         CALL EZspline_init(F_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZsplineinit',ier)        
         F_spl%isHermite   = 1
         F_spl%x1   = raxis
         F_spl%x2   = phiaxis
         F_spl%x3   = zaxis
         CALL EZspline_setup(F_spl,XRHO_ARR,ier,EXACT_DIM=.true.)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZspline_setup',ier)
         XRHO4D = F_spl%fspl
         CALL EZspline_free(F_spl,ier)
      END IF
      CALL mpidealloc(XRHO_ARR,win_XRHO_ARR)

!-----------------------------------------------------------------------
!     YRHO_ARR
!----------------------------------------------------------------------- 
      CALL mpialloc(YRHO4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_YRHO4D)
      IF (myid_sharmem == master) THEN
         CALL EZspline_init(F_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZsplineinit',ier)        
         F_spl%isHermite   = 1
         F_spl%x1   = raxis
         F_spl%x2   = phiaxis
         F_spl%x3   = zaxis
         CALL EZspline_setup(F_spl,YRHO_ARR,ier,EXACT_DIM=.true.)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZspline_setup',ier)
         YRHO4D = F_spl%fspl
         CALL EZspline_free(F_spl,ier)
      END IF
      CALL mpidealloc(YRHO_ARR,win_YRHO_ARR)

      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE beams3d_spline3d_setup
