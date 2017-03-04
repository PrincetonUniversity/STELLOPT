!-----------------------------------------------------------------------
!     Subroutine:    stelltran_ecrh
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/14/2016
!     Description:   Calculate ECRH heating and current drive
!-----------------------------------------------------------------------
      SUBROUTINE stelltran_ecrh(itime)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE stelltran_runtime
      USE stelltran_vars
      USE stelltran_data, ONLY: P_ecrh
      USE stelltran_equilutils
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER,INTENT(in):: itime
      INTEGER        :: temp_len,i, iflag,j,n, ier
      REAL(rprec)    :: Pabs,jll,It
      
      ! FOR TRAVIS
      INTEGER(4)     :: nra,nphi,QOemul, nbeams
      INTEGER(4)     :: travisoutTEMP,travisScrOut
      INTEGER(4)     :: maxSteps,stopray,npass,maxHarm,nrho,nu
      INTEGER(8)     :: mconf8=0,mVessel=0
      REAL(8)        :: hgrid,dphi,B_scale,B0_ref,phi_ref,phibx
      REAL(8)        :: maxlength,minStepSize,maxStepSize,odetolerance,umax
      REAL(8)        :: antennaPosition(3),targetPosition(3),rbeam(2),rfocus(2)
      REAL(8), ALLOCATABLE :: te_prof(:),ne_prof(:),z_prof(:),s_prof(:)
      CHARACTER,DIMENSION(3)   :: wmode
      CHARACTER(4)   :: antennaCoordType,targetPositionType
      CHARACTER(256) :: B0type,labelType,equiname
      CHARACTER(256) :: hamilt,dieltensor

      INTERFACE
         INTEGER(4) FUNCTION mcload(mconf,name) 
            INTEGER(8)   :: mconf
            CHARACTER(*) :: name
         END FUNCTION mcLoad
      END INTERFACE
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      SELECT CASE(TRIM(equil_type))
         CASE('vmec2000','animec','flow','satire')
!DEC$ IF DEFINED (TRAVIS)
            IF (lveryfirst_pass) THEN
               WRITE(6,'(A)')             '==========================================='
               WRITE(6,'(A,A,A)')         '==============  T R A V I S  =============='
            END IF
            ! Load travis input file
            labelType = 'ech'
            CALL set_hedi_f77(labelType)
 
            ! Load equilibrium
            equiname = 'wout_'//TRIM(proc_string)//'.nc'
            ier = mcLoad(mconf8,equiname)
            IF (ier /=1) THEN
               iflag = -1
               RETURN
            END IF
 
            ! Setup profiles (on uniform internal grid)
            ALLOCATE(s_prof(prof_length),te_prof(prof_length),&
                     ne_prof(prof_length),z_prof(prof_length))
            s_prof = te(:,1)
            te_prof  = te(:,2)*1.0E-3 ! in keV
            ne_prof  = ne(:,2)
            z_prof   = zeff(:,2)

            ! Handle Output
            travisOutTemp = 0; travisScrOut=0;
            !IF (lscreen) travisScrOut = 0;
            CALL set_TracingOutput_f77(travisOutTemp,travisScrOut)

            ! Handle Vessel and Mirrors
            vessel_ecrh = TRIM(vessel_ecrh)
            mirror_ecrh = TRIM(mirror_ecrh)
            CALL set_VesselMirrors_f77(vessel_ecrh,mirror_ecrh)
            ! Control Output
            CALL Disable_Output_f77

            ! Init Equilibrium
            hgrid   = Aminor/30; dphi = 2; B0type  = 'at angle on magn.axis'; phi_ref = 0; B_scale = 1; B0_ref = Baxis
            CALL initEquilibr_f77(hgrid, dphi, B0_ref, phi_ref, B_scale, B0type)

            ! Set Configuration
            CALL SET_MAGCONFIG_F77(mConf8,mVessel)

            ! Set Profiles
            labelType = 'tor_rho'
            CALL set_nTZ_profs_f77(prof_length, s_prof, ne_prof, te_prof, z_prof, labelType)

            ! Setup Beams
            antennaCoordType     = antennatype_ecrh(1:4)
            targetPositionType     = targettype_ecrh(1:4)
            nra = nra_ecrh
            nphi = nphi_ecrh
            nbeams = 0
            DO i = 1,nsys
               IF (ANY(antennaposition_ecrh(i,:) .ne. 0)) nbeams = nbeams + 1
            END DO
            CALL setup_ECRH_beams_f77(nbeams)
            DO i = 1,nbeams
               wmode(1)      = wmode_ecrh(i,1)
               wmode(2)      = wmode_ecrh(i,2)
               wmode(3)      = wmode_ecrh(i,3)
               antennaPosition(1:3) = antennaPosition_ecrh(i,1:3)
               targetPosition(1:3) = targetPosition_ecrh(i,1:3)
               rbeam(1:2) = rbeam_ecrh(i,1:2)
               QOemul = 0
               IF (rbeam_ecrh(i,3) > 0) QOemul = 1
               Rfocus(1:2) = rfocus_ecrh(i,1:2)
               phibx       = rfocus_ecrh(i,3)
               CALL setup_ECRH_beam_f77(i,freq_ecrh(i),P_ecrh(i,itime),&
                                        rbeam,Rfocus,phibx,QOemul,nra,nphi,&
                                        antennaPosition,targetPosition,&
                                        antennaCoordType,targetPositionType,wmode)
            END DO
            CALL BeamInfoAllocate_f77()

            ! Run the beams
            CALL run_beams_f77

            IF (lveryfirst_pass .and. lverb) WRITE(6,'(7X,A,6(9X,A))') 'i', 'rho','dPdV [MW/m^3]','   jll','    jrf','    It','    Ptot'
            DO i = 1, prof_length
!                CALL get_ECRH_deposition_f77(s_prof(i),dPdV,Pabs,jll,jrf(i,2),It)
!                IF (lveryfirst_pass .and. lverb) WRITE(6,'(5X,i3,3x,f8.2,4(5x,f10.3))') i,s_prof(i),&
!                                                                    dPdV,jll,jrf(i,2),It
               CALL get_ECRH_deposition_f77(s_prof(i),dPdV_erf(i,2),Ptot(i,2),jll,jrf(i,2),It)
               IF (lveryfirst_pass .and. lverb) WRITE(6,'(5X,i3,3x,f8.2,5(5x,f10.3))') i,s_prof(i),&
                                                                   dPdV_erf(i,2)/1E6,jll,jrf(i,2),It,Ptot(i,2)/1E6
            END DO

!DEC$ ENDIF
      END SELECT
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE stelltran_ecrh
