c                                                                                      
c  L-BFGS-B is released under the “New BSD License” (aka “Modified BSD License”        
c  or “3-clause license”)                                                              
c  Please read attached file License.txt                                               
c                                        

      double precision function dnrm2(n,x,incx)
      integer n,incx
      double precision x(n)
c     **********
c
c     Function dnrm2
c
c     Given a vector x of length n, this function calculates the
c     Euclidean norm of x with stride incx.
c
c     The function statement is
c
c       double precision function dnrm2(n,x,incx)
c
c     where
c
c       n is a positive integer input variable.
c
c       x is an input array of length n.
c
c       incx is a positive integer variable that specifies the 
c         stride of the vector.
c
c     Subprograms called
c
c       FORTRAN-supplied ... abs, max, sqrt
c
c     MINPACK-2 Project. February 1991.
c     Argonne National Laboratory.
c     Brett M. Averick.
c
c     **********
      integer i
      double precision scale

      dnrm2 = 0.0d0
      scale = 0.0d0

      do 10 i = 1, n, incx
         scale = max(scale, abs(x(i)))
   10 continue

      if (scale .eq. 0.0d0) return

      do 20 i = 1, n, incx
         dnrm2 = dnrm2 + (x(i)/scale)**2
   20 continue

      dnrm2 = scale*sqrt(dnrm2)

 
      return

      end
      
c====================== The end of dnrm2 ===============================


