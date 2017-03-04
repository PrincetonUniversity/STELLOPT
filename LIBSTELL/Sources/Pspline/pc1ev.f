      subroutine pc1ev(xget,x,nx,ilinx,f,ict,fval,ier)
C
C  evaluate a 1d piecewise linear interpolant -- this is C0.
C    ...a derivative can be evaluated but it is not continuous.
C
C  this subroutine calls two subroutines:
C     herm1x   -- find cell containing (xget,yget)
C     pc1fcn -- evaluate interpolant function and (optionally) derivatives
C
C  input arguments:
C  ================
C
      real xget                         ! target of this interpolation
      real x(nx)                        ! ordered x grid
      integer ilinx                     ! ilinx=1 => assume x evenly spaced
C
      real f(nx)                        ! function data
C
C       contents of f:
C
C  f(i) = f @ x(i)
C
      integer ict(2)                    ! code specifying output desired
C
C  ict(1)=1 -- return f  (0, don't)
C  ict(2)=1 -- return df/dx  (0, don't)
C
C output arguments:
C =================
C
      real fval(*)                      ! output data
      integer ier                       ! error code =0 ==> no error
C
C  fval(1) receives the first output (depends on ict(...) spec)
C  fval(2) receives the second output (depends on ict(...) spec)
C
C  examples:
C    on input ict = [1,1]
C   on output fval= [f,df/dx]
C
C    on input ict = [1,0]
C   on output fval= [f] ... element 2 never referenced
C
C    on input ict = [0,1]
C   on output fval= [df/dx] ... element 2 never referenced
C
C  ier -- completion code:  0 means OK
C-------------------
C  local:
C
      integer :: i=0                         ! cell index
C
C  normalized displacement from (x(i)) corner of cell.
C    xparam=0 @x(i)  xparam=1 @x(i+1)
C
      real xparam
C
C  cell dimensions and
C  inverse cell dimensions hxi = 1/(x(i+1)-x(i))
C
      real hx
      real hxi
C
C  0 .le. xparam .le. 1
C
C---------------------------------------------------------------------
C
      call herm1x(xget,x,nx,ilinx,i,xparam,hx,hxi,ier)
      if(ier.ne.0) return
c
      call pc1fcn(ict,1,1,fval,i,xparam,hx,hxi,f,nx)
C
      return
      end
C---------------------------------------------------------------------
C  evaluate C0 piecewise linear function interpolation -- 1d fcn
C   --vectorized-- dmc 10 Feb 1999
C
      subroutine pc1fcn(ict,ivec,ivecd,
     >   fval,ii,xparam,hx,hxi,fin,nx)
C
C  input:
C
      integer ict(2)                    ! requested output control
      integer ivec                      ! vector length
      integer ivecd                     ! vector dimension (1st dim of fval)
C
      integer ii(ivec)                  ! target cells
      real xparam(ivec)
                                        ! normalized displacements in cells
C
      real hx(ivec)                     ! grid spacing, and
      real hxi(ivec)                    ! inverse grid spacing 1/(x(i+1)-x(i))
C
      real fin(nx)                      ! the data
C
C  output:
C
      real fval(ivecd,*)                ! interpolation results
C
C  for detailed description of fin, ict and fval see subroutine
C  pc1ev comments.  Note ict is not vectorized -- the same output
C  is expected to be returned for all input vector data points.
C
C  note that the index inputs ii, and parameter inputs
C     xparam,hx,hxi,are vectorized, and the
C
C     output array fval has a vector ** 1st dimension ** whose
C     size must be given as a separate argument; ivecd.ge.ivec
C     expected!
C
C  to use this routine in scalar mode, pass in ivec=ivecd=1
C
C---------------
C
      integer v
C
      do v=1,ivec
         i=ii(v)
         xp=xparam(v)
         xpi=1.0-xp
C
         iadr=0
C
C  get desired values:
C
         if(ict(1).eq.1) then
C
C  function value:
C
            iadr=iadr+1
            fval(v,iadr)=xpi*fin(i)+xp*fin(i+1)
C
         endif
C
         if(ict(2).eq.1) then
C
C  df/dx:
C
            iadr=iadr+1
C
            fval(v,iadr)=(fin(i+1)-fin(i))*hxi(v)
         endif
C
      enddo                             ! vector loop
C
      return
      end
