      MODULE coiltypes
      USE stel_constants

         TYPE modularcoil
            REAL(rprec), DIMENSION(:), POINTER :: rhoc, rhos,
     1         phic, phis, rho, phi
            REAL(rprec) :: current
         END TYPE modularcoil

         TYPE saddlecoil
            REAL(rprec), DIMENSION(:), POINTER ::
     1         v_c, v_s, u_c, u_s, phi, theta
            REAL(rprec) :: current
         END TYPE saddlecoil

         TYPE vfcoil
            REAL(rprec) :: radius, height, current
         END TYPE vfcoil

      END MODULE coiltypes
