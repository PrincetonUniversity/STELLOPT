!-----------------------------------------------------------------------
!     Subroutine:    stellopt_bcast_vmec
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2012
!     Description:   This subroutine broadcasts the entire VMEC INDATA
!                    namelist because the processors in the local_comm
!                    never read the INDATA namelist.  Essentially it
!                    mimics the VMEC initialization routines.  This was
!                    probably overkill but at least the whole namelist
!                    gets updated.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_bcast_vmec(local_master,local_comm,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE vmec_input
      USE mpi_params
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      INCLUDE 'mpif.h'                                                          ! MPI
!DEC$ ENDIF  
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      INTEGER,INTENT(INOUT) :: local_master
      INTEGER,INTENT(INOUT) :: local_comm
      INTEGER,INTENT(INOUT) :: iflag
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
!DEC$ IF DEFINED (MPI_OPT)
      iflag = 0
      ! Logicals
      CALL MPI_BCAST(lpofr,1,MPI_LOGICAL,local_master,local_comm,iflag)
      CALL MPI_BCAST(lmac,1,MPI_LOGICAL,local_master,local_comm,iflag)
      CALL MPI_BCAST(lfreeb,1,MPI_LOGICAL,local_master,local_comm,iflag)
      CALL MPI_BCAST(lrecon,1,MPI_LOGICAL,local_master,local_comm,iflag)
      CALL MPI_BCAST(loldout,1,MPI_LOGICAL,local_master,local_comm,iflag)
      CALL MPI_BCAST(ledge_dump,1,MPI_LOGICAL,local_master,local_comm,iflag)
      CALL MPI_BCAST(lasym,1,MPI_LOGICAL,local_master,local_comm,iflag)
      CALL MPI_BCAST(lforbal,1,MPI_LOGICAL,local_master,local_comm,iflag)
      CALL MPI_BCAST(lrfp,1,MPI_LOGICAL,local_master,local_comm,iflag)
      CALL MPI_BCAST(lmovie,1,MPI_LOGICAL,local_master,local_comm,iflag)
      CALL MPI_BCAST(lmove_axis,1,MPI_LOGICAL,local_master,local_comm,iflag)
      CALL MPI_BCAST(lwouttxt,1,MPI_LOGICAL,local_master,local_comm,iflag)
      CALL MPI_BCAST(ldiagno,1,MPI_LOGICAL,local_master,local_comm,iflag)
      CALL MPI_BCAST(lmoreiter,1,MPI_LOGICAL,local_master,local_comm,iflag)
      CALL MPI_BCAST(lfull3d1out,1,MPI_LOGICAL,local_master,local_comm,iflag)
      CALL MPI_BCAST(l_v3fit,1,MPI_LOGICAL,local_master,local_comm,iflag)
      CALL MPI_BCAST(lspectrum_dump,1,MPI_LOGICAL,local_master,local_comm,iflag)
      CALL MPI_BCAST(loptim,1,MPI_LOGICAL,local_master,local_comm,iflag)
      CALL MPI_BCAST(lgiveup,1,MPI_LOGICAL,local_master,local_comm,iflag)
      CALL MPI_BCAST(lbsubs,1,MPI_LOGICAL,local_master,local_comm,iflag)
      CALL MPI_BARRIER(local_comm,iflag)
      ! Integers
      CALL MPI_BCAST(nfp,1,MPI_INTEGER,local_master,local_comm,iflag)
      CALL MPI_BCAST(ncurr,1,MPI_INTEGER,local_master,local_comm,iflag)
      CALL MPI_BCAST(nsin,1,MPI_INTEGER,local_master,local_comm,iflag)
      CALL MPI_BCAST(niter,1,MPI_INTEGER,local_master,local_comm,iflag)
      CALL MPI_BCAST(nstep,1,MPI_INTEGER,local_master,local_comm,iflag)
      CALL MPI_BCAST(nvacskip,1,MPI_INTEGER,local_master,local_comm,iflag)
      CALL MPI_BCAST(mpol,1,MPI_INTEGER,local_master,local_comm,iflag)
      CALL MPI_BCAST(ntor,1,MPI_INTEGER,local_master,local_comm,iflag)
      CALL MPI_BCAST(ntheta,1,MPI_INTEGER,local_master,local_comm,iflag)
      CALL MPI_BCAST(nzeta,1,MPI_INTEGER,local_master,local_comm,iflag)
      CALL MPI_BCAST(mfilter_fbdy,1,MPI_INTEGER,local_master,local_comm,iflag)
      CALL MPI_BCAST(nfilter_fbdy,1,MPI_INTEGER,local_master,local_comm,iflag)
      CALL MPI_BCAST(max_main_iterations,1,MPI_INTEGER,local_master,local_comm,iflag)
      CALL MPI_BCAST(imse,1,MPI_INTEGER,local_master,local_comm,iflag)
      CALL MPI_BCAST(isnodes,1,MPI_INTEGER,local_master,local_comm,iflag)
      CALL MPI_BCAST(itse,1,MPI_INTEGER,local_master,local_comm,iflag)
      CALL MPI_BCAST(ipnodes,1,MPI_INTEGER,local_master,local_comm,iflag)
      CALL MPI_BCAST(iopt_raxis,1,MPI_INTEGER,local_master,local_comm,iflag)
      CALL MPI_BCAST(imatch_phiedge,1,MPI_INTEGER,local_master,local_comm,iflag)
      CALL MPI_BCAST(nflxs,1,MPI_INTEGER,local_master,local_comm,iflag)
      CALL MPI_BARRIER(local_comm,iflag)
      ! Integer Arrays
      CALL MPI_BCAST(ns_array,100,MPI_INTEGER,local_master,local_comm,iflag)
      CALL MPI_BCAST(niter_array,100,MPI_INTEGER,local_master,local_comm,iflag)
      CALL MPI_BCAST(nbfld,nbsetsp,MPI_INTEGER,local_master,local_comm,iflag)
      CALL MPI_BCAST(indxflx,nfloops,MPI_INTEGER,local_master,local_comm,iflag)
      CALL MPI_BCAST(indxbfld,nbcoilsp*nbsetsp,MPI_INTEGER,local_master,local_comm,iflag)
      CALL MPI_BARRIER(local_comm,iflag)
      ! Reals
      CALL MPI_BCAST(time_slice,1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(curtor,1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(delt,1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(ftol,1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(tcon0,1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(gamma,1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(phiedge,1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(phidiam,1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(sigma_current,1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(sigma_delphid,1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(tensi,1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(tensp,1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(tensi2,1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(fpolyi,1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(presfac,1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(mseangle_offset,1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(pres_offset,1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(mseangle_offsetm,1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(spres_ped,1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(bloat,1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(pres_scale,1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(prec2d_threshold,1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(bcrit,1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(fgiveup,1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BARRIER(local_comm,iflag)
      ! Real Arrays
      CALL MPI_BCAST(rbc,(2*ntord+1)*(mpol1d+1),MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(zbs,(2*ntord+1)*(mpol1d+1),MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(rbs,(2*ntord+1)*(mpol1d+1),MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(zbc,(2*ntord+1)*(mpol1d+1),MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(am,21,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(ai,21,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(ac,21,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(aphi,20,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(ah,21,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(at,21,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(am_aux_s,ndatafmax,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(am_aux_f,ndatafmax,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(ac_aux_s,ndatafmax,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(ac_aux_f,ndatafmax,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(ai_aux_s,ndatafmax,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(ai_aux_f,ndatafmax,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(ah_aux_s,ndatafmax,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(ah_aux_f,ndatafmax,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(at_aux_s,ndatafmax,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(at_aux_f,ndatafmax,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(raxis,ntord+1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(zaxis,ntord+1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(raxis_cc,ntord+1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(raxis_cs,ntord+1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(zaxis_cc,ntord+1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(zaxis_cs,ntord+1,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(ftol_array,100,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(extcur,nigroup,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(mseprof,nmse,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(rthom,ntse,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(datathom,ntse,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(sigma_thom,ntse,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(rstark,nmse,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(datastark,nmse,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(sigma_stark,nmse,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(dsiobt,nfloops,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(sigma_flux,nfloops,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(bbc,nbcoilsp*nbsetsp,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(sigma_b,nbcoilsp*nbsetsp,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(psa,ndatafmax,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(pfa,ndatafmax,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(isa,ndatafmax,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BCAST(ifa,ndatafmax,MPI_DOUBLE_PRECISION,local_master,local_comm,iflag)
      CALL MPI_BARRIER(local_comm,iflag)
      ! Strings
      CALL MPI_BCAST(pcurr_type,20,MPI_CHARACTER,local_master,local_comm,iflag)
      CALL MPI_BCAST(piota_type,20,MPI_CHARACTER,local_master,local_comm,iflag)
      CALL MPI_BCAST(pmass_type,20,MPI_CHARACTER,local_master,local_comm,iflag)
      CALL MPI_BCAST(pt_type,20,MPI_CHARACTER,local_master,local_comm,iflag)
      CALL MPI_BCAST(ph_type,20,MPI_CHARACTER,local_master,local_comm,iflag)
      CALL MPI_BCAST(mgrid_file,200,MPI_CHARACTER,local_master,local_comm,iflag)
      CALL MPI_BCAST(trip3d_file,200,MPI_CHARACTER,local_master,local_comm,iflag)
      CALL MPI_BCAST(precon_type,10,MPI_CHARACTER,local_master,local_comm,iflag)
      CALL MPI_BCAST(arg1,120,MPI_CHARACTER,local_master,local_comm,iflag)
      CALL MPI_BCAST(input_extension,100,MPI_CHARACTER,local_master,local_comm,iflag)
      CALL MPI_BARRIER(local_comm,iflag)
      IF (iflag /= 0) RETURN
!DEC$ ENDIF

      ! DO A CHECK HERE
      iflag = 0
      !IF (rbc(0,1) < 0) iflag = -1
      !IF (rbc(0,0) < (SUM(SUM(rbc,DIM=2),DIM=1)-rbc(0,0))) iflag = -1
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_bcast_vmec
