        subroutine init_dgmres(icntl,cntl)
        implicit none
*
*  Purpose
*  =======
*    Set default values for the parameters defining the characteristics
* of the Gmres algorithm.
*  See the User's Guide for an example of use.
*
*
* Written : April 1997
* Authors :  Valerie Fraysse, Luc Giraud, Serge Gratton
*             Parallel Algorithms - CERFACS
*
*
*  Arguments
*  =========
*
* icntl    (input) INTEGER array. length 6
*            icntl(1) : stdout for error messages
*            icntl(2) : stdout for warnings
*            icntl(3) : stdout for convergence history
*            icntl(4) : 0 - no preconditioning
*                       1 - left preconditioning
*                       2 - right preconditioning
*                       3 - double side preconditioning
*                       4 - error, default set in Init
*            icntl(5) : 0 - modified Gram-Schmidt
*                       1 - iterative modified Gram-Schmidt
*                       2 - classical Gram-Schmidt
*                       3 - iterative classical Gram-Schmidt
*            icntl(6) : 0 - default initial guess x_0 = 0 (to be set)
*                       1 - user supplied initial guess
*            icntl(7) : maximum number of iterations
*            icntl(8) : 1 - default compute the true residual at each restart
*                       0 - use recurrence formula at restart
*
* cntl     (input) real*8 array, length 5
*            cntl(1) : tolerance for convergence
*            cntl(2) : scaling factor for normwise perturbation on A
*            cntl(3) : scaling factor for normwise perturbation on b
*            cntl(4) : scaling factor for normwise perturbation on the
*                      preconditioned matrix
*            cntl(5) : scaling factor for normwise perturbation on 
*                      preconditioned right hand side
*
* Output variables
* ----------------
       integer icntl(*)
       real*8   cntl(*)
*
       icntl(1) = 6
       icntl(2) = 6
       icntl(3) = 0
       icntl(4) = 4
       icntl(5) = 0
       icntl(6) = 0
       icntl(7) = -1
       icntl(8) = 1
* 
       cntl(1) = 1.0 d -5
       cntl(2) = 0.0 d 0
       cntl(3) = 0.0 d 0
       cntl(4) = 0.0 d 0
       cntl(5) = 0.0 d 0
*
       return
       end

