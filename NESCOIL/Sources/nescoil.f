! ----------------------------------------------------------------------
      subroutine nescoil (loopcount)
!................................................................
      USE stel_kinds
      USE Vvacuum2, ONLY: iota_edge, phip_edge, curpol
      USE Vvacuum3
      USE Vprecal1, ONLY: np
      USE Vbiopre
      USE Vculine1
      USE Vprecal4
      use SvdCtrl, ONLY : noaccuracy
      use LoopCtrl
      use Vbnorm, ONLY : extension
      use NumParams, ONLY : inesc
      use safe_open_mod
      implicit none
!................................................................
C   D u m m y   A r g u m e n t s
c................................................................
      integer :: loopcount
c................................................................
C   L o c a l   V a r i a b l e s
c................................................................
      integer :: istat
      real(rprec) :: unit, file, t1, t2, tbeg, tend
      character*10 :: date, time

c................................................................
c
c      The NESCOIL code is based on a method described in
c         [1]  P.Merkel,Nucl.Fusion 27 (1987) 867.
c         [2]  P.Merkel, In: Theory of Fusion Plasmas,Eds. A.Bondeson
c              E.Sindoni,and F.Troyon, Varenna,Italy,EUR 11336 EN,25-46.
c
c      Introducing two nested toroidally closed surfaces a current
c      on the outer surface is determined such that the normal
c      component of the magnetic field on the inner surface is
c      minimized.
c
c      The current carrying surface is mapped onto the unit square:
c      0< u <1, 0< v <1.
c
c
c         ^
c       1.|
c         |                  cup - poloidal current per period
c         |                  cut - toroidal current
c         |
c   ^     |
c   |   u |
c   |     |
c         |
c  cup    |
c         |
c         |
c         |
c       0 ----------------------------->
c         0                v           1.
c                     <----- cut
c
c
c
c................................................................
c     PARAMETERS READ IN FROM INPUT FILE:
c
c     1. DIMENSIONS:
c        nu    - number of poloidal meshpoints for coil surface
c        nv    - number of toroidal meshpoints for coil surface
c
c        nu1   - number of poloidal meshpoints for plasma surface
c        nv1   - number of toroidal meshpoints for plasma surface
c
c        mf    - number of poloidal modes for current potential
c        nf    - number of toroidal modes for current potential
c
c        md    - number of poloidal modes for plasma and coil surface shape
c        nd    - number of toroidal modes for plasma and coil surface shape
c
c        npol  - number of segments of a modular or helical filament
c        ntor  - number of filaments per period
c
c        np    - number of field periods
c
c     2. NESCOIL CONTROL PARAMETERS:
c        cut   - Net toroidal current +-1 or 0
c        cup   - Net poloidal current +-1 or 0
c
c        ibex  - An  external field can be added by supplying
c                a subroutine bexter with IBEX = 1
c
c     3. INFORMATION FROM VMEC PLASMA:
c        iota_edge - iota at plasma edge
c        phip_edge - phi_prime at plasma edge
c        curpol - Total poloidal current Amps per field period from VMEC
c
c     4. SVD CONTROL PARAMETERS:
c        1. MSTRT = Method + svdscan start if >1
c                 >=0:Berr, <0:Xerr, unless MSTEP <=0: Least square
c        2. MSTEP = Method + svdscan stepsize:
c                 <=0: LeastSquare, =0: use old f04abe (now SOLVER), no svd
c        3. MKEEP = svd/scan control: =0: svdscan, else keep |nkeep| wgts
c                 <0: write all weights to output
c        4. MDSPW = 2+exponent of dsur multiplying bfn,ben:
c                 <0: post-calculate Xerr svdscan and write to output
c        5. CURWT = Weight for surface current minimization
c                   Works ONLY in LSQ branch
c        6. TRGWT : Not implemented yet (PMV)

c................................................................
c        Secondary parameters:
c
c        md,nd - dimension of arrays: md>=mf,nd>=nf
c        nmax  - dimension of array:  nmax> 3*(nu+2)*nv
c
c................................................................
c      Set loop counter iloop for internal use
      iloop = loopcount

c     Strategy based on iloop:
c      1. iloop = 0 or  1 : single run or first call, allocate
c      2. iloop = 0 or -1 : single run or last  call, deallocate
c................................................................
c     Read the input file only in the first loop
c     Read input data from file specified on command line:
c     xnesopt input
c     NOTE: SINCE ALL DIMENSIONS ARE ALLOCATABLE NOW,
c           YOU MUST FIRST READ INPUT FILE (TO GET DIMENSIONS)
c           BEFORE DOING ANYTHING ELSE
      if (iloop == 0 .or. iloop == 1) call read_nescoil_input
c     This ends first-loop tasks, rest is done in each loop
c................................................................
c
      write (inesc, '(a,i3,a)') '-----  Begin Nescoil run num ', iloop,
     1    ' -----'

c     Now that all inputs are here, begin nescoil run
      call date_and_time(date,time)
      write (inesc, 100) date(5:6),date(7:8),date(1:4),
     1  time(1:2),time(3:4),time(5:6)
 100  format('DATE = ',a2,'-',a2,'-',a4,' ',' TIME = ',2(a2,':'),a2)
      call second0(tbeg)   !This one for timing the whole run
c
c     setup constants, arrays, and allocations based on nwire
c     uncomment for NT      USE NUMERICAL_LIBRARIES
c
      nw = ntor*np*(npol+1)
      istat=0
      if (.not.allocated(vx))
     1 allocate (vx(nw), vy(nw), vz(nw), dx(nw), dy(nw), dz(nw),
     1           xw(nw), yw(nw), zw(nw), curre(nw),
     2           cok(0:np - 1), sik(0:np - 1), stat=istat)
      curre = 0; xw = 0; yw = 0; zw = 0                                !!MUST be initialized (SPH)
      if (istat .ne. 0) stop 'allocation error nw in NESCOIL'

c     Do all precalculations
      call precal
      call second0(t2)
      write (inesc,"('PRECAL took ',g12.3,' sec')") t2-tbeg

c     Compute quantities on plasma surface
      write (inesc, 10) '----- Calling Surface_Plasma -----'
      call second0(t1)
      call surface_plas
      call second0(t2)
      write (inesc,"('Time in Surface_Plasma: ',g12.3,' sec')") t2-t1

c     Compute quantities on coil surface
      write (inesc, 10) '----- Calling Surface_Coil -----'
      call second0(t1)
      call surface_coil
      call second0(t2)
      write (inesc,"('Time in Surface_Coil: ',g12.3,' sec')") t2-t1

c................................................................
c     Solve boundary value problem
      write (inesc, 10) '----- Calling Solver -----'
      call second0(t1)
      call solver_nescoil
      call second0(t2)
      write (inesc,"('Time in Solver: ',g12.3,' sec')") t2-t1

c................................................................
c     Post-process solution for various answers
      write (inesc, 10) '----- Calling Surfcur_Diag -----'
      call second0(t1)
      call surfcur_diag
      call second0(t2)
      write (inesc,"('Time in Surfcur_Diag: ',g12.3,' sec')") t2-t1

      if(noaccuracy .eq. 0) then
         write (inesc, 10) '----- Calling Accuracy -----'
         call second0(t1)
         call accuracy
         call second0(t2)
         write (inesc,"('ACCURACY took ',g12.3,' sec')") t2-t1
      endif

      call second0(tend)
      write (inesc,"('ONE NESCOIL RUN took ',g12.3,' sec')")
     1    tend-tbeg

 10   format(a)
c................................................................
c     Tasks to be done in last loop only
      if (iloop .le. 0) call nescoil_cleanup
c................................................................

      end subroutine nescoil
