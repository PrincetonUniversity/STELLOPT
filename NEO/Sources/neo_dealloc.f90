
SUBROUTINE neo_dealloc
! **********************************************************************
! Modules
! **********************************************************************
  USE neo_precision
  USE neo_input
  USE neo_work
  USE neo_exchange
  USE neo_units
  USE neo_parameters
  USE neo_control
  USE neo_spline
! **********************************************************************
! Local Definitions
! **********************************************************************
  IMPLICIT NONE
! *******************************************************************
! DeAllocate Storage Arrays
! *******************************************************************
  !DEALLOCATE (es,iota,curr_pol,curr_tor,pprime,sqrtg00)
  !DEALLOCATE (theta_arr,phi_arr)
  !DEALLOCATE (rmnc,zmns,lmns,bmnc)
  !DEALLOCATE (ixm, ixn)
  !DEALLOCATE (pixm, pixn)
  !DEALLOCATE (i_m, i_n)
  !DEALLOCATE (b_spl,k_spl,g_spl,p_spl)
  !DEALLOCATE (b,sqrg11,kg,pard)
  !IF (calc_cur .EQ. 1) THEN
  !   DEALLOCATE (bqtphi)
  !   DEALLOCATE (q_spl)
  !END IF
  !DEALLOCATE (fluxs_arr)
  ! For protection SAL 05/29/14
  IF (ALLOCATED(es)) DEALLOCATE(es)
  IF (ALLOCATED(iota)) DEALLOCATE(iota)
  IF (ALLOCATED(curr_pol)) DEALLOCATE(curr_pol)
  IF (ALLOCATED(curr_tor)) DEALLOCATE(curr_tor)
  IF (ALLOCATED(pprime)) DEALLOCATE(pprime)
  IF (ALLOCATED(sqrtg00)) DEALLOCATE(sqrtg00)
  IF (ALLOCATED(theta_arr)) DEALLOCATE(theta_arr)
  IF (ALLOCATED(phi_arr)) DEALLOCATE(phi_arr)
  IF (ALLOCATED(rmnc)) DEALLOCATE(rmnc)
  IF (ALLOCATED(zmns)) DEALLOCATE(zmns)
  IF (ALLOCATED(lmns)) DEALLOCATE(lmns)
  IF (ALLOCATED(bmnc)) DEALLOCATE(bmnc)
  IF (ALLOCATED(ixm)) DEALLOCATE(ixm)
  IF (ALLOCATED(ixn)) DEALLOCATE(ixn)
  IF (ALLOCATED(pixm)) DEALLOCATE(pixn)
  IF (ALLOCATED(i_m)) DEALLOCATE(i_n)
  IF (ALLOCATED(b_spl)) DEALLOCATE(b_spl)
  IF (ALLOCATED(k_spl)) DEALLOCATE(k_spl)
  IF (ALLOCATED(g_spl)) DEALLOCATE(g_spl)
  IF (ALLOCATED(p_spl)) DEALLOCATE(p_spl)
  IF (ALLOCATED(b)) DEALLOCATE(b)
  IF (ALLOCATED(sqrg11)) DEALLOCATE(sqrg11)
  IF (ALLOCATED(bqtphi)) DEALLOCATE(bqtphi)
  IF (ALLOCATED(q_spl)) DEALLOCATE(q_spl)
  IF (ALLOCATED(fluxs_arr)) DEALLOCATE(fluxs_arr)
! *******************************************************************
  RETURN
END SUBROUTINE neo_dealloc
