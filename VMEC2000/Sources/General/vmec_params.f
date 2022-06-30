      MODULE vmec_params
      USE stel_kinds, ONLY: rprec, dp
      USE vparams, ONLY: mpold
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: meven = 0, modd = 1
      INTEGER, PARAMETER :: ndamp = 10
      INTEGER, PARAMETER :: ns4 = 25

      INTEGER, PRIVATE :: ink
      INTEGER, PARAMETER, DIMENSION(0:mpold) ::
     1  jmin1 = (/ 1,1,(2,ink=2,mpold) /),        !starting js(m) values where R,Z are non-zero
     2  jmin2 = (/ 1,2,(2,ink=2,mpold) /),        !starting js(m) values for which R,Z are evolved
     3  jlam  = (/ 2,2,(2,ink=2,mpold) /)         !starting js(m) values for which Lambda is evolved

!  Besure to update werror in fileout.f when adding more error flags.
      INTEGER, PARAMETER :: norm_term_flag=0, bad_jacobian_flag=1, 
     1                      more_iter_flag=2, 
     2                      jac75_flag=4, input_error_flag=5,
     3                      phiedge_error_flag=7,
     4                      ns_error_flag=8,
     5                      misc_error_flag=9,
     6                      successful_term_flag=11, !ftol force criterion has been met
     7                      bsub_bad_js1_flag=12,
     8                      r01_bad_value_flag=13,
     9                      arz_bad_value_flag=14,
     1                      imas_read_flag=15
      INTEGER, PARAMETER :: restart_flag=1, readin_flag=2,
     1                      timestep_flag=4,output_flag=8, 
     2                      cleanup_flag=16, reset_jacdt_flag=32,
     3                      imasrun_flag = 64
    
      REAL(rprec), PARAMETER :: pdamp = 0.05_dp  
      CHARACTER(LEN=*), PARAMETER :: version_ = '9.0'
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ntmax, rcc, rss, rsc, rcs, zsc, zcs, zcc, zss
      INTEGER :: mnyq, nnyq
      INTEGER, ALLOCATABLE :: uminus(:)
      REAL(rprec), ALLOCATABLE :: mscale(:), nscale(:)
      REAL(rprec) :: signgs, lamscale
!-----------------------------------------------
!
!     VERSION INFORMATION HISTORY
!
!     8.00 
!          a) added lforbal logical to fbal module to control whether to compute the flux-averaged
!             force balance equation printed in the threed1 file. This requires a modification of
!             the m=1,n=0 forces for R,Z in tomnsps subroutine. This works well, generally, and
!             yields an improved <EQUIF> in threed1 file. However, there have been some cases where
!             this non-variational departure fails to converge.
!          b) added "bias" to iotas in bcovar when determining preconditioner for very small iota
!             values. Seems to need this for improving convergence of preconditioned current-hole cases.
!             Eliminated in v8.20.
!          c) added "totzsps,a_hess" to recompute r,z,l rapidly for only those modes that are jogged during
!             Hessian calculation. NOTE: need to do this for lasym=true case, as well, eventually
!     8.20  (August, 2004)
!          a) removed 2-pt tie at origin for m=1 modes, to preserve tri-diagonal structure of Hessian.
!             This is needed for preconditioning, which assumes block-tridi structure of equations
!          b) fixed problem with free-boundary preconditioner, namely, ctor can not be extrapolated
!             at edge when computing preconditioner, because this breaks tri-diagonal structure
!          c) added new variables to input file to control preconditioning:
!             1) PRECON_TYPE: = 'default', default tri-di (block size = 1)
!                             = 'cg',      block tri-di, conjugate-gradient time-stepper
!                             = 'gmres',   "          ", gmres time-stepper
!                             = 'tfqmr',   "          ", transpose free qmr
!             2) PREC2D_THRESHOLD: value of (unpreconditioned) forces at which block (2D) preconditioner
!                                  is turned on (=0 block preconditioner never turned on); recommended 
!                                  (default) value ~ 1.E-10, or smaller, if convergence is poor
!             3) LFORBAL: logical variable (default = .true.); when true, the force balance
!                         used in the threed1 file is used to evolve the m=1 R,Z components. This
!                         is a non-variational departure from the force equations for these modes,
!                         but generally does not have an unfavorable impact on convergence.
!          d) added new internal variable, ICTRL_PREC2D, to precon2d module. Replaces previous lprec2d
!             and iequi>1 variables.
!          e) removed lsweep_fast option from module precon2d. This slows the computation of the Hessian
!             by about 2/3, but is more accurate (includes pdamp, liota, lforbal correctly)
!          f) removed lflam logicals from bcovar and tomnsps, since now we must compute dFlam/dR,Z by
!             jogging
!          g) removed Compute_Hess_Flam_RZ from lamblks; this is now computed via jogging
!             (also removed Get_dFlam_dRZ, FFT2Hessian, Forbal_avg, GetGsqrtVar supporting routines)
!          h) removed internal liota logic, used to push iota profile rather than solving for it. Had
!             needed this for symmetric Hessian (lsweep_fast=true option), but no longer required. Also,
!             it was not implemented entirely correctly for lforbal=true case
!          i) for lasym m=1 constraint rsc = zcc, changed xc array so that R+ = .5*(rsc + zcc) is stored at
!             xc(rsc,m=1) and R- = .5*(rsc - zcc) is stored at xc(zcc,m=1). In residue, gcz(R-) == gcz(zcc)
!             is zeroed by "evolving" gcr(zcc) = azd*[xc(zcc)-xcint], and gcr(rsc) => .5*[gcr(rsc) + gcz(zcc)] 
!             is evolved. In totzspa, the real rsc and zcc are computed from the internal representations 
!             (check convert call, too) by calling a new routine convert_asym (also called from wrout before 
!             writing out xc info). In readin, the original R+,R- are stored, so that for external "jogs", 
!             there will be no change in forces. All these changes are needed to obtain an invertible Hessian.
!          j) added m=1 constraint for 3D case (similar to (i)), rss(n) = zcs(n), for n != 0. Imposed this
!             on forces by adding routine constrain_m1 in residue. Added convert_sym routine to totzsp to convert
!             from internal xc representation TO internal one. 
!          k) Decreased exponent on pdamp factor r2 (in bcovar) from 2 to 1, to give better conditioning
!             especially for current hole cases
!          l) Eliminated iotas bias for determining preconditioner, previously added in v8.00 for stabilizing
!             current hole cases (not needed with corrected preconditioner)
!     8.30  (October, 2004)
!          a) Implemented flags for "reverse-communication" mode of vmec
!     8.40 a) Converted the m=1 constraints for 3D and asym back to old way; did not always
!             converge well with the new constraints introduced in 8.20 (i-j)
!     8.45  (December, 2005)
!          a) Added the lconm1 logical. If = True, new constraint; if = False, old m=1 constraint used
!          b) Added "perturbation" computation for lasym=TRUE case (totzspa_hess)
!     8.46  (June, 2009)
!          a) Added LRFP logical to allow easy switching of profiles between Stellarator/tokamak (PHIP=1, LRFP=F)
!             and RFP (CHIP=1, LRFP=T). When LRFP=T, AI coefficients are expansion of q = 1/iota. Added lrfp to
!             LIBSTELL/vmec_input module.
!     8.47  (July, 2010)
!          a) Rewrote magnetic field representation so that phip*lambda = new internal lambda. This greatly improves
!             the conditioning of the lambda equations which otherwise become singular at the RFP reversal point
!     8.48  (March 2012 - JDH)
!          a) Accumulated small changes from SPH & JDH
!          b) Modifications from J Geiger, March 2012
!             - to be able to get additional main iterations if the force tolerance is
!                  not met. Parameter MAX_MAIN_ITERATIONS specifies how many main iteration
!                  cycles should be run.
!             - to get a full output in the threed1-file if the force tolerance is not
!                  met. Specify the logical LFULL3D1OUT to true for this.
!             - if vmec2000 is compiled with netcdf, you can still get the ascii-output
!                  if you specify the logical LWOUTTXT as true.
!             - you get the output for diagno 1.0 and 1.5 if the logical LDIAGNO set true.
!             - you get a rather old fort.8-output if you specify LOLDOUT as true.
!             
!             If none of these new variables is set, the behavior of vmec2000 is as
!             expected from the version without the changes.
!     8.49  (June, 2012)
!          a) Fixed bug in bcovar when averaging half-grid covariant components onto full grid: did not
!             zero components at (phantom) js=1 point, so edge force averages were incorrect
!          b) Added lamscale factor to scale lambda in contravariant B-components. Modify wrout to maintain
!             old-style lambda output
!          c) Set lbsubs=F in jxbforce by default to capture current sheets
!          d) Added lmove_axis INPUT logical (=T by default) so user can control whether or not the magnetic
!             axis can be initially shifted to improve the initial force residuals. It is important NOT to move
!             the helical axis for RFP equilibria requiring a helical seed (set l_moveaxis=F for this case!)
!
!     8.50   (Jan, 2013)
!          a) Improved scaling of lambda forces with respect to lamscale
!          b) Fixed fnorm1 scaling (remove hs dependence)
!          c) Added lgiveup logical (M. Drevlak/J. Geiger)
!-----------------------------------------------
      END MODULE vmec_params
