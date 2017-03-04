      SUBROUTINE chisq_coilgeom (num, nopt, iflag, extension, lscreen)
      USE stel_kinds
      USE real_ptr_type
      USE chisq_mod
      USE coilsnamin
      USE safe_open_mod
      USE optim, ONLY: home_dir, exe_suffix, remove, sigma_berr_avg, 
     1                 sigma_berr_max, bigno
      USE system_mod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nopt
      INTEGER :: num, iflag
      CHARACTER(LEN=*) :: extension
      LOGICAL :: lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: n_list=9
      INTEGER, PARAMETER :: coilgeom0 = 28
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: j, k, iunit, istat, ierr, ndim_cc
      TYPE (real_ptr), DIMENSION(n_list) :: weights, weightm
      REAL(rprec), DIMENSION(1), TARGET :: scalar1
      real(rprec), dimension(2) :: no_weights
      LOGICAL :: access, vf_reg, polcur, mplcoil, splcoil, vf_rmax,
     1           vac_match, lccflag, berrflag
      LOGICAL, DIMENSION(n_list) :: mflag, sflag
      INTEGER, DIMENSION(n_list) :: ivars
      CHARACTER(LEN=30) :: header
      CHARACTER(LEN=200) :: output_file, cmd

!****************************************************************************
!
!   **** NOTE OF IMPORTANCE ****
!
!  This routine must be maintained in concert with COILGEOM/LSFUN1
!
!****************************************************************************
!
!   Setup the flags indicating what is to be done
!   initialize

      mflag = .false.
      sflag = .false.
      lccflag = .false.
      access = .false.
      vf_reg = .false.
      vf_rmax = .false.
      polcur = .false.
      mplcoil = .false.
      splcoil = .false.
      vac_match = vacfld_wgt .ne. zero
      ndim_cc = 0

      no_weights = 1
!   setup ivars for the initial 9 types of penalties

      ivars(1) = ivar_coil_len
      ivars(2) = ivar_coil_ymin
      ivars(3) = ivar_coil_sep
      ivars(4) = ivar_coil_plcross
      ivars(5) = ivar_coil_curv
      ivars(6) = ivar_coil_icurv               !!integrated curvature along coil
      ivars(7) = ivar_coil_reg
      ivars(8) = ivar_coil_curdens
      ivars(9) = ivar_coil_rmax                !!largest R of coil (keeps bounded in bell jar)


!   here we set up flags indicating whether there are any values we want
!   to include for the initial 5 types of penalties, for both modular and
!   saddle representations.   Note that we must carefully parallel the logic
!   used in COILGEOM in writing these out.

!   'modular' representation

      IF (lmodular .and. nmid >= 1) THEN
         mflag(1) = (.not. lsaddle) .and. ANY(lmod_wgt(1:nmid).ne.zero)
         mflag(2) = (.not. lsaddle) .and. ANY(ymin_wgt(1:nmid).ne.zero)
         mflag(3) = ANY( dcc_wgt(1:nmid) .ne. zero)
         mflag(4) = .false.                                      ! someday
         mflag(5) = ANY( rc_wgt(1:nmid) .ne. zero)
         mflag(6) = ANY( cu_wgt(1:nmid) .ne. zero)
         mflag(7) = .false.
         mflag(8) = .false.
         mflag(9) = .false.

         weightm(1)%x => lmod_wgt
         weightm(2)%x => ymin_wgt
         weightm(3)%x => dcc_wgt
         weightm(5)%x => rc_wgt
         weightm(6)%x => cu_wgt

         mplcoil = dcp_wgt .ne. zero
      ENDIF

!   'saddle' representation

      IF (lsaddle .and. nsmid >= 1) THEN
         sflag(1) = ANY(lsad_wgt(1:nsmid) .ne. zero)
         sflag(2) = ANY(ymin_wgt(1:nsmid) .ne. zero)
         sflag(3) = ANY(dsc_wgt(1:nsmid) .ne. zero)
         sflag(4) = ANY(dscxp_wgt(1:nsmid) .ne. zero)
         sflag(5) = ANY(rs_wgt(1:nsmid) .ne. zero)
         sflag(6) = ANY(cs_wgt(1:nsmid) .ne. zero)
         sflag(7) = ANY(csc_wgt(1:nsmid) .ne. zero)
         sflag(8) = ANY(scd_wgt(1:nsmid) .ne. zero)
         sflag(9) = ANY(rmax_wgt(1:nsmid) .ne. zero)
         
         weights(1)%x => lsad_wgt
         weights(2)%x => ymin_wgt
         weights(3)%x => dsc_wgt
         weights(4)%x => dscxp_wgt
         weights(5)%x => rs_wgt
         weights(6)%x => cs_wgt
         weights(7)%x => csc_wgt
         weights(8)%x => scd_wgt
         weights(9)%x => rmax_wgt
        
         lccflag = ANY(sc_dmin_wgt .ne. zero)
         DO k = 1, nsad_unique_coils
            sc_dmin_wgt(k,k) = 0
         END DO
         ndim_cc = nsad_unique_coils*nsad_coils

         splcoil = dcp_wgt .ne. zero
      ENDIF

!   set flags for the penalties that are independent of representation

!   Berr
      berrflag = lbnorm .and. (sigma_berr_avg < bigno .or.
     1                       sigma_berr_max < bigno)
     
!   access

      IF (laccess .and. n_access > 0)
     1   access = ANY (dac_wgt(1:n_access) .ne. zero)

!   PF specific: regularization and rmax

      IF (lvf .and. num_vf > 0 ) THEN
         vf_reg = ANY(cvf_wgt(1:num_vf) .ne. zero)
         vf_rmax = ANY(rvf_wgt(1:num_vf) .ne. zero) .and. lvfvar
      ENDIF

!   total poloidal current

      IF( lpolcur) polcur = dpc_wgt .ne. zero

! ------------------------------------------------

      IF (nopt > 0) THEN

!   if necessary, calculate Bnorm

         if( lbnorm) then
            cmd = trim(home_dir) // '/xbnorm wout.' // trim(extension)
            write (cmd,'(a,1pe12.3)') trim(cmd), 0.20
            if (lscreen) then
               print '(/,a)',' Running BNORM code...'
            end if
            call system(cmd)                     !!Produce bnorm file
         endif

!   execute xcoilgeom again (called previously from generate_mgrid) to generate plasmas-dependent penalties

         output_file = 'coil_targets.' // TRIM(extension)
         cmd = TRIM(home_dir) // 'xcoilgeom' // exe_suffix // 
     1         TRIM(extension) // ' -P > coildump.' // extension
         CALL system(cmd, ierr)
         IF (ierr.lt.127 .and. ierr.ne.0) THEN
            IF (lscreen) PRINT *,
     1                   'XCOILGEOM failed in chisq_coilgeom call: ',
     2                           'ierr = ',ierr
            iflag = -22
            RETURN
         END IF
         
         cmd = remove // 'coildump.' // TRIM(extension)
         CALL system(cmd)

!
!  Open and read coil_targets file generated when coils file written
!
         iunit = coilgeom0

         CALL safe_open(iunit, k, output_file, 'old', 'formatted')
         IF (k .ne. 0) THEN
            iflag = -25
            RETURN
         ENDIF

         header = ' '
         DO WHILE (INDEX(header,'coilgeom penalties') == 0)
            READ(iunit,'(a)',iostat=istat) header
            IF (istat .ne. 0) STOP ' ISTAT != 0 in chisq_coilgeom(1)'
         ENDDO

!   the first penalties are for Berr avg. & max

         call read_coilgeom_elem(iunit, num, 2, no_weights,
     1                              berrflag, ivar_coil_aBerr)
         if( berrflag) then
            index_array(num) = ivar_coil_mBerr
            chisq_descript(num) = descript(ivar_coil_mBerr)
            wegt(num) = sigma_berr_max
            wegt(num-1) = sigma_berr_avg
         endif

!   The second penalty is the chi-sq due to matching the VACUUM bnorm field on 
!   the FINITE BETA plasma boundary (only done if vac_match = .true.)

         scalar1(1) = vacfld_wgt
         CALL read_coilgeom_elem(iunit, num, 1, scalar1, 
     1                           vac_match, ivar_vacuum_field)

!   The next n_list  penalties occur in pairs, first for the modular representation,
!   then for the saddle representation.  This pairing is broken due to some
!   penalties not being implemented for one representation or the other in
!   COILGEOM/COILOPT, causing the extra tests

         DO j=1, n_list
            IF (lmodular .and. (j .ne. 2 .or. (.not. lsaddle))) THEN
               CALL read_coilgeom_elem(iunit, num, nmid, weightm(j)%x,
     1                                 mflag(j), ivars(j))
            END IF

            IF (lsaddle) THEN
               CALL read_coilgeom_elem(iunit, num, nsmid, weights(j)%x,
     1                                 sflag(j), ivars(j))
            END IF
         END DO

!   Coil-Coil distances (matrix form)
         DO k = 1, nsad_coils
            weights(1)%x => sc_dmin_wgt(1:nsad_unique_coils,k)
            IF (lsaddle) CALL read_coilgeom_elem(iunit, num, 
     1                        nsad_unique_coils, weights(1)%x, 
     2                        lccflag, ivar_coil_sep_mat)
         END DO

!   Coil-Plasma distances
         scalar1(1) = dcp_wgt;  weights(1)%x => scalar1
         IF (lmodular) CALL read_coilgeom_elem(iunit, num, 1,
     1                    weights(1)%x, mplcoil, ivar_coil_pldist)

         IF (lsaddle) CALL read_coilgeom_elem(iunit, num, 1,
     1                          weights(1)%x, splcoil, ivar_coil_pldist)

!   Now we have penalties that do not depend on representations used

         IF (laccess) THEN
            weights(1)%x => dac_wgt
            CALL read_coilgeom_elem(iunit, num, n_access, weights(1)%x,
     1                              access, ivar_coil_access)
         ENDIF

         IF (lvf) THEN
            weights(1)%x => cvf_wgt
            CALL read_coilgeom_elem(iunit, num, num_vf, weights(1)%x,
     1                              vf_reg, ivar_coil_vfreg)
         ENDIF

         IF (lvfvar) THEN
            weights(1)%x => rvf_wgt
            CALL read_coilgeom_elem(iunit, num, num_vf, weights(1)%x,
     1                              vf_rmax, ivar_coil_vfrmax)
         ENDIF

         IF (lpolcur) THEN
            scalar1(1) = dpc_wgt;    weights(1)%x => scalar1
            CALL read_coilgeom_elem(iunit, num, 1, weights(1)%x, polcur,
     1                              ivar_coil_polcur)
         ENDIF

         CLOSE (iunit)

      ELSE
         IF (vac_match) THEN 
            num = num + 1
            IF (nopt .eq. -2) chisq_descript(num) = 
     1                        descript(ivar_vacuum_field)
         END IF

         DO j=1, n_list
            IF (mflag(j)) THEN
               IF (nopt .eq. -2) chisq_descript(num+1:num+nmid) =
     1            descript(ivars(j))
               num = num + nmid
            ENDIF

            IF (sflag(j)) THEN
               IF (nopt .eq. -2) chisq_descript(num+1:num+nsmid) =
     1            descript(ivars(j))
               num = num + nsmid
             ENDIF
         ENDDO
         
         IF (lccflag) THEN
            DO k = 1, nsad_coils
               IF (nopt .eq. -2) 
     1            chisq_descript(num+1:num+nsad_unique_coils) =
     2            descript(ivar_coil_sep_mat)
               num = num+nsad_unique_coils
            ENDDO
         ENDIF

         IF (lmodular .and. mplcoil) THEN
            num = num+1
            IF (nopt .eq. -2) chisq_descript(num) =
     1                        descript(ivar_coil_pldist)
         ENDIF

         IF (lsaddle .and. splcoil) THEN
            num = num+1
            IF (nopt .eq. -2) chisq_descript(num) =
     1                        descript(ivar_coil_pldist)
         ENDIF

         IF (berrflag) THEN
            IF (nopt .eq. -2) chisq_descript(num+1:num+2) =
     1                        descript(ivar_coil_access)
            num = num + 2
         ENDIF

         IF (access) THEN
            IF (nopt .eq. -2) chisq_descript(num+1:num+n_access) =
     1                        descript(ivar_coil_access)
            num = num+n_access
         ENDIF

         IF (vf_reg) THEN
            IF (nopt .eq. -2) chisq_descript(num+1:num+num_vf) =
     1                        descript(ivar_coil_vfreg)
            num = num + num_vf
         ENDIF

         IF (vf_rmax) THEN
            IF (nopt .eq. -2) chisq_descript(num+1:num+num_vf) =
     1                        descript(ivar_coil_vfrmax)
            num = num + num_vf
         ENDIF

         IF (polcur) THEN
            num = num + 1
            IF (nopt .eq. -2) chisq_descript(num) =
     1                        descript(ivar_coil_polcur)
         ENDIF

         iflag = 0
      END IF

      END SUBROUTINE chisq_coilgeom
