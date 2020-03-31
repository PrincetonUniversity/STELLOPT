
      subroutine bootsj(aibstot, extension, iunit_in)
c
!  The code BOOTSJ calculates the bootstrap current for 3D configurations.
!  The main part of the code, calculation of the geometrical factor GBSNORM,
!  was written by Johnny S. Tolliver of ORNL on the basis of
!  Ker-Chung Shaing's formulas. Other parts of the code, as well as modifications
!  to the Tolliver's part, were written by Paul Moroz of UW-Madison.
c
!  References:
c
!  1. K.C. Shaing, B.A. Carreras, N. Dominguez, V.E. Lynch, J.S. Tolliver
!    "Bootstrap current control in stellarators", Phys. Fluids B1, 1663 (1989).
!  2. K.C. Shaing, E.C. Crume, Jr., J.S. Tolliver, S.P. Hirshman, W.I. van Rij
!     "Bootstrap current and parallel viscosity in the low collisionality
!     regime in toroidal plasmas", Phys. Fluids B1, 148 (1989).
!  3. K.C. Shaing, S.P. Hirshman, J.S. Tolliver "Parallel viscosity-driven
!     neoclassical fluxes in the banana regime in nonsymmetri! toroidal
!     plasmas", Phys. Fluids 29, 2548 (1986).
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      use parambs
      use safe_open_mod
      use trig
      use vmec0
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer     :: iunit_in
      real(rprec) :: aibstot
      character*(*) :: extension
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: nfax = 13
      integer, parameter :: indata0 = 7
      integer, parameter :: jbs_file=59, ans_file=18, ans_dat_file=19
      real(rprec) :: one = 1, p5 = .5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C----------------------------------------------- 
      integer, dimension(nfax) :: ifaxu, ifaxv
      integer :: ntrigu, ntrigv
      integer :: irho, irho1, ierr, iunit, ijbs, ians, ians_plot
      real(rprec), dimension(:), allocatable :: cputimes
      real(rprec) :: time1, timecpu, unit, file, status, err,
     1   time2, r, x, al31t, gradbs1, gradbs2,
     2   gradbs3, gradbs4,  al31s
      integer :: ihere = 0
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      real(rprec) , EXTERNAL :: al31
C-----------------------------------------------

c
c  define the starting CPU time
c
      call second0 (time1)
      timecpu = time1
      if (lscreen) write (*, 4) version_
    4 format(/,' Start BOOTSJ: Version ', a)

c  open files

      iunit = indata0
      call safe_open(iunit, ierr, 'input.' // trim(extension), 'old',
     1     'formatted')
      if (ierr .ne. 0) then
         print *,' Error opening input file: input.', trim(extension)
         return
      end if

      ijbs = jbs_file
      call safe_open(ijbs, ierr, 'jBbs.'//trim(extension), 'replace',
     1     'formatted')

      ians = ans_file
      call safe_open(ians, ierr, 'answers.'//trim(extension), 'replace',
     1     'formatted')

      ians_plot = ans_dat_file
      call safe_open(ians_plot, ierr, 'answers_plot.'//trim(extension),
     1     'replace', 'formatted')

c  read and initialize data

      call datain(trim(extension), iunit_in, iunit, ijbs, ians)
      close (iunit)

      ntrigu = 3*nthetah/2 + 1
      ntrigv = 2*nzetah
      allocate (cputimes(irup))
      allocate (
     1  dmn(-mbuse:mbuse,0:nbuse), fmn(-mbuse:mbuse,0:nbuse),
     2  rfmn(-mbuse:mbuse,0:nbuse),alpha1mn(-mbuse:mbuse,0:nbuse),
     3  trigsv(ntrigv),trigsu(ntrigu), stat=irho)

      if (irho .ne. 0) stop 'allocation error in bootsj main'
c     convert jlist to idx form, then "and" the two lists
c     SPH: Corrected error here, writing off end of jlist error when jlist <= 1
c


      jlist_idx = 0
      do irho = 1, irup
        if (jlist(irho) .gt. 1) jlist_idx(jlist(irho)-1) = 1
        idx(irho) = idx(irho)*jlist_idx(irho)
      enddo


      l_boot_all = .true.

!  if any of the boozer evaluation surfaces are missing, or not requested
!  then set l_boot_all=false, total current can not be calculated

      do irho=1, irup
         if(idx(irho) .eq. 0) l_boot_all = .false.
      enddo
      if(.not.l_boot_all) then
         if (lscreen) write (*,*) 'partial surface evaluation'
         write (ians,*) 'partial surface evaluation'
      endif


      call fftfax_g (nthetah, ifaxu, trigsu)
      call cftfax_g (nzetah, ifaxv, trigsv)


c start main radial loop

      do irho = 1, irup

c  initialize timing for radius step

         call second0 (time2)
         timecpu = time2 - time1
         cputimes(irho) = timecpu

c  if there is no boozer information available, skip radial point

         if(idx(irho) .eq. 0)  cycle
         irho1 = irho - 1
         r = sqrt(rhoar(irho) + 1.E-36_dp)

c  initialize  angle grids grid for first radial evaluation.  For this
c  and all subsequent radial evaluate B and related and quantites as well
c  as plasma derivatives and the tokamak trapped fraction.

         call bongrid(irho, ians, ihere)

c  calculate bootstrap current for equivalent tokamak

         x = fttok(irho)/(fptok(irho)+1.E-36_dp)

         al31t = al31(x,zeff1,alphae,alphai)

c  calculate gradient factors overall normalization inlcuding q and
c  boozer g.

         call grad (gradbs1, gradbs2, gradbs3, gradbs4, irho)

         bsdenste(irho) = gradbs1*al31t          !due to dens gradient
         bsdensti(irho) = gradbs2*al31t          !due to dens gradient
         bstempte(irho) = gradbs3*al31t          !due to temp gradient
         bstempti(irho) = gradbs4*al31t          !due to temp gradient

         dibst(irho) = bsdenste(irho) + bsdensti(irho) + bstempte(irho)
     1       + bstempti(irho)                    !total Jbst

         if (l_boot_all) then
            if (irho .eq. 1) then
               aibst(1) = dibst(1)*d_rho(1)
            else
               aibst(irho) = aibst(irho1)+dibst(irho)*d_rho(irho)
            endif
         end if

c  Now start the general evaluation.
c  Find coefficients d(n,m) and evaluate the fraction trapped.

         call denmf (trigsu, trigsv, ifaxu, ifaxv, irho)

c  Evaluate R, S, and H2.

         call caprsh2(irho)

c  Evaluate the integral term.

         call woflam (trigsu, trigsv, ifaxu, ifaxv, irho)

c  Evaluate the summation term in W(lambda) that does not depend on lambda.
c  note that while the paper shows an itegral, it is the fpassing integral
c  that cancels the fpassing in the over all multiplier

         call othersums(irho)

c  Calculate the final answer

         amain(irho) = p5*(one - aiogar(irho)/qsafety(irho)) + p5*(one +
     1      aiogar(irho)/qsafety(irho))*h2(irho)

         gbsnorm(irho) = amain(irho) + other1(irho) +  aiterm1(irho)

c- derivative of the enclosed bootstrap current over the normalized
c  toroidal flux, dIbs/ds, and the enclosed Ibs (in MA)

         x = ftrapped(irho)/(fpassing(irho)+1.E-36_dp)
         al31s = al31(x,zeff1,alphae,alphai)
         call grad (gradbs1, gradbs2, gradbs3, gradbs4, irho)
c
c                                                !due to dens gradient
         bsdense(irho) = gbsnorm(irho)*gradbs1*al31s
c                                                !due to dens gradient
         bsdensi(irho) = gbsnorm(irho)*gradbs2*al31s
c                                                !due to temp gradient
         bstempe(irho) = gbsnorm(irho)*gradbs3*al31s
c                                                !due to temp gradient
         bstempi(irho) = gbsnorm(irho)*gradbs4*al31s

         dibs(irho) = bsdense(irho) + bsdensi(irho) + bstempe(irho) +
     1      bstempi(irho)                        !total Jbst (dI/ds)

c   convert to j dot B.  2*dmu0 from beta, 10**6 from MA, psimax is 1/dpsi/ds
c   and sign_jacobian takes out the sign previously used for dpsi to dA
CCCCC  flux was changed to real flux in boot vmec so that the psimax now
CCCCC  needs an additional sign_jacobian so that the two will cancel
             ajBbs(irho) = (2.0e6_dp)*dmu0*dibs(irho)*
     1      (pres1(irho)/betar(irho))/psimax

         if (l_boot_all) then
            if (irho .eq. 1) then
               aibs(1) = dibs(1)*d_rho(1)
            else
               aibs(irho) = aibs(irho1) + dibs(irho)*d_rho(irho)
            endif
         end if

c- the ratio of bootstrap current to that in ESD (equivalent symmetric device):

         bsnorm(irho) = dibs(irho)/(dibst(irho)+1.E-36_dp)

c  get time at end of loop if completed

         call second0 (time2)
         timecpu = time2 - time1
         cputimes(irho) = timecpu
c

      end do
c
c- Output answers for BOOTSJ
c
      call output (cputimes, aibstot, ijbs, ians, ians_plot)
      close(ians)
      close(ians_plot)
      close(ijbs)

      call deallocate_all
      if (lscreen) write (*,400) (cputimes(irup)-cputimes(1))
  400 format(1x,'Finished BOOTSJ, time =  ', f8.3, '  sec')

      deallocate (cputimes, trigsu, trigsv)
      deallocate (dmn, fmn, rfmn, alpha1mn)

      end subroutine bootsj
