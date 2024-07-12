      PROGRAM dkes2
c-----------------------------------------------------------------------
!!!!
!                              D   I   S   C   L   A   I   M   E   R
!
!       You are using a BETA version of the PROGRAM DKES, which is currently
!       under development by S. P. Hirshman at the Fusion Energy Division,
!       Oak Ridge National Laboratory.  Please report ANY problems or comments
!       to him.  As a BETA version, this PROGRAM is subject to change
!       and improvement without notice.
!
!!!!
c-----------------------------------------------------------------------
c     SPH:         implemented in UNIX, F77                 February 21, 1994
c     SPH:         upgraded to F90, dynamic memory allocation,
c                  matmul routines                          May 1999
c-----------------------------------------------------------------------
c
c
c   Example usage:
c
c   OR   xdkes INPUT_FILE_NAME
c
c        xdkes BOOZ_FILE_NAME NSURF CMUL EFIELD LSCREEN
c
c
c  WHERE
c
c        INPUT_FILE_NAME:   name of the input file prepared by dkes calling
c                           subroutine dkes_input_prepare (the old way of running DKES)
c
c        BOOZ_FILE_NAME:    name of Boozer coordinate file; see dkes_input_prepare SUBROUTINE for
c                           description of this and nsurf, cmul, efield, lscreen parameters
c
c
c  DKES.VAR (Drift Kinetic Equation Solver, Variational) solves a set
c  of 3-D drift kinetic equations to obtain upper and lower bounds for
c  the diffusion coefficients of a prescribed toroidal plasma
c  equilibrium. the 3 dimensions are theta (poloidal angle),
c  zeta (toroidal angle), and alpha (pitch angle). Straight-line flux
c  coordinates are used to describe the equilibrium, which satisfies
c  the stellarator symmetry conditions R(theta,zeta) = R(-theta,-zeta),
c  Z(theta,zeta) = - Z(-theta,-zeta).
c
c  Reference: W. I. van Rij and S. P. Hirshman, Variational Bounds for
c  Transport Coefficients in Three-Dimensional Plasmas,
c  Phys. Fluids B 1,563(1989).
c
c  Boozer Coordinate Version
c
c-----------------------------------------------------------------------
c
c  mpnt      : The total number of distribution modes is computed from
c              the mmnn matrix in the MNSET subroutine
c  iswpm = 1 : "Plus" sources, distributions used  (Maximizing)
c  iswpm = 2 : "Minus" sources,distributions used  (Minimizing)
c
c-----------------------------------------------------------------------
c
c  SUBROUTINES required:
c
c       name:     purpose:
c
c  DKES2 essential routines:
c       blk5d     solves the block-pentadiagonal system of
c                 equations
c       blox      forms the l-row block matrices and rows
c       cescale   cmul and efield scaling and
c                 the dominant-diagonal scaling arrays
c       ftconv    computes magnetic field spectral arrays and Fourier transforms
c       lcalc     degree of Legendre polynomial arrays
c       printout  calculate diffusion coefficients and output
c       residue   calculates the solution residuals
c       reverse   returns solution arrays to their original
c                 order with respect to l
c       scalel    dominant-diagonal l matrix scaling
c       wrout     output of magnetic field and of final
c                 distribution functions
c
c  DKES2 auxiliary routines:
c       rddisk    auxilary routine for fast internal disk i/o
c       wrdisk    auxilary routine for fast internal disk i/o
c
c-----------------------------------------------------------------------
c
c     MAIN PROGRAM
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE Vimatrix
      USE Vnamecl2
      USE Vcmain
      USE dkes_input
      USE dkes_realspace
      USE safe_open_mod                              !record file addition
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ir, irun, istat, neqs
      REAL(rprec), DIMENSION(:), POINTER :: f0p1, f0p2, f0m1, f0m2
      REAL(rprec) :: tcpu0, tcpu1, tcpui, tcput, tcpu, tcpua
      CHARACTER :: tb*1                              !record file addition
C-----------------------------------------------------------------------
      DKES_rad_dex = 0                               ! SAL for STELLOPT
      tb = char(9)                                   !record file addition
      CALL second0 (tcpu0)
!
!     OPEN FILES FOR IO AND READ INPUT FILE DATA
!
      CALL read_dkes_input

!
!     set-up : calculate spectral arrays, matrix elements, and l arrays
!
      CALL set_mndim

      neqs = mpnt*(lalpha + 1)                                           !add l = -1 "constraint"
      ALLOCATE (fzerop(2*neqs), fzerom(2*neqs), stat=istat)
      IF (istat .ne. 0) STOP 'allocation error(1) in dkes2!'


      f0p1 => fzerop(1:neqs);   f0p2 => fzerop(neqs+1:)                  !+ functions
      f0m1 => fzerom(1:neqs);   f0m2 => fzerom(neqs+1:)                  !- functions

      CALL ftconv

      CALL lcalc

      CALL second0 (tcpu1)
      tcpui = tcpu1 - tcpu0

!  header for file DKESOUT; output magnetic field

      CALL header

      CALL second0 (tcpu0)
      tcput = zero

!  cmul, efield loops
      nrun = MAX(nrun, 1)
      nrun = MIN(nrun, krun)

      DO ir = 1, nrun
         irun = ir
         efield1 = efield(irun)
         cmul1 = cmul(irun)

c     record file addition
      if(ir .eq. 1) then
        call safe_open(itab_out, istat, summary_file,'unknown',
     >  'formatted')
        write(itab_out,'("*",/,"cmul",a1,"efield",a1,"weov",a1,"wtov",
     >   a1,"L11m",a1,"L11p",a1,"L31m",a1,"L31p",a1,"L33m",a1,"L33p",
     >   a1,"scal11",a1,"scal13",a1,"scal33",a1,"max_residual",
     >   a1,"chip",a1,"psip",a1,"btheta",a1,"bzeta",a1,"vp")')
     >   tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb
      else if(ir .gt. 1) then
        open(itab_out,file=summary_file,status='unknown',
     >   position='append',form='formatted')
      endif
c     end record file addition

!  cmul and efield scaling; scaling arrays

         CALL cescale (srces0)

         WRITE (ioout, 950) dashes, cmul1, efield1, weov, wtov,
     1          wcyclo, vthermi
  950    FORMAT(/9x,'CMUL',7x,'EFIELD',4x,'OMEGA-E/v',4x,'OMEGA-T/v',
     1         5x,'OMEGA-ci',3x,'VI-THERMAL'/
     2         7x,'(nu/v)',7x,'(Es/v)',2x,'(ExB drift)',4x,'(transit)',
     3         3x, '(H+, B=B00)',2x,'(H+, Ti=1keV)'/,131a,/
     4         1p,6e13.4/)

!  block-pentadiagonal solutions and residuals

         IF (ipmb < 2) THEN
            iswpm = 1
            CALL blk5d (blk1, blk2, blk3, blk4, blk5, blk6,
     1               blk7, f0p1, f0p2, srces0)

            IF (ier .ne. 0) THEN
               WRITE (ioout, 1000) ier
 1000          FORMAT(/' blk5d error in block = ',i5)
               STOP
            ENDIF
            CALL residue_dkes (blk1, blk2, blk3, blk4, blk5, blk6,
     1           f0p1, f0p2, srces0, rsd1p, rsd3p, g11p, g33p,
     2           g31p, g13p, crs1p, crs3p)
         ENDIF

         IF (ipmb .ne. 1) THEN
            iswpm = 2
            CALL blk5d (blk1, blk2, blk3, blk4, blk5, blk6,
     1         blk7, f0m1, f0m2, srces0)

            IF (ier .ne. 0) THEN
               WRITE (ioout, 1050) ier
 1050          FORMAT(/' blk5d error : ierm = ',i5)
               STOP
            ENDIF
            CALL residue_dkes (blk1, blk2, blk3, blk4, blk5, blk6,
     1           f0m1, f0m2, srces0, rsd1m, rsd3m, g11m, g33m,
     2           g31m, g13m, crs1m, crs3m)
         ENDIF

!  calculate and output diffusion coefficients

         CALL dkes_printout (f0p1, f0m1, f0p2, f0m2, srces0)

!  timing and check remaining run time

         CALL second0 (tcpu1)
         tcpu = tcpu1 - tcpu0
         tcpu0 = tcpu1
         tcput = tcput + tcpu
         tcpua = tcput/irun
         WRITE (ioout, 1100) tcpu
!        IF (irun<nrun .and. tcpu1<1.05*tcpua+3.0) EXIT
 1100    FORMAT(/' time used:    tcpu =',1p,e10.2,'  sec'/)
         close(unit=itab_out)                         !record file addition
      END DO
      CLOSE(unit=ioout_opt)

!  output magnetic field and final distributions

      IF (lfout .ne. 0) CALL wrout (f0p1, f0m1, f0p2, f0m2,
     1   srces0)



!  clean-up memory
      CALL free_mndim
      DEALLOCATE (cols, al1, al2, al3, al4, bl1, bl2, bl3, bl4, cl1,
     1   cl2, cl3, cl4, cols0, omgl, al01, al02, al03, al04, bl01,
     2   bl02, bl03, bl04, cl01, cl02, cl03, cl04, fzerop, fzerom)

!  time and date

      tcput = tcput + tcpui
      WRITE (ioout, 1200) dashes, tcpui, tcpua, tcput
      IF (irun < nrun) THEN
         WRITE (ioout, '(a,i3,a,i3)')
     1      'DKES2 run not completed (time limit):   irun = ',
     2      irun, '   nrun =', nrun
      ELSE
         WRITE (ioout, '(a,i3)') ' DKES2 run completed: nrun = ', nrun
      ENDIF

      WRITE (ioout, '(1x,131a)') dashes
 1200 FORMAT(1x,131a,/,1p,' tcpui =',e9.2,' s',5x,'tcpua =',e9.2,
     1    ' s', 5x,'tcput =',e9.2,' s')

      END PROGRAM dkes2
